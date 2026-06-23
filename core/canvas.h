/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// platform.h first: on device it defines NDEBUG, which must be set before
// <cassert> expands the assert macro — otherwise assert-stripping would depend
// on a prior TU having pulled in platform.h, making this header non-self-sufficient.
#include "platform.h"
#include <cstring>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <variant>
#include <atomic>
#include "constants.h"
#include "color.h"
#include <array>

class Canvas;

/**
 * @brief Base class for all visual effects.
 * @details Manages double buffering, persistence, and provides an interface for
 * drawing a frame.
 */
class Effect {
  friend class Canvas;

public:
  bool debug_visuals = false; /**< Flag to enable visual debugging overlays. */

  /**
   * @brief Constructs an Effect instance.
   * @param W The width (resolution) of the effect.
   * @param H The height (resolution) of the effect.
   */
  Effect(int W, int H) : persist_pixels(false), width_(W), height_(H) {
    // Single-live-Effect precondition (see buffer_a/buffer_b): every Effect
    // aliases the same two static buffers, so two simultaneously-live instances
    // would scribble over each other's frames. The contract holds everywhere
    // today — WASM resets currentEffect before constructing the next, the device
    // driver deletes the old effect before building its successor, and the test
    // harness constructs effects one at a time — so this traps a future overlap
    // loudly at the cold construction site instead of silently tearing frames.
    HS_CHECK(!s_alive,
             "Effect: a second Effect was constructed while one is still alive; "
             "buffer_a/buffer_b are shared static storage (one live Effect only)");
    s_alive = true;
    // Point bufs_ at the shared static storage and clear both buffers. Factored
    // into a noinline helper so the two MAX_W*MAX_H fills are emitted once rather
    // than duplicated into both GCC constructor variants (C1 complete + C2 base)
    // — the same code-size concern init() exists to address. They stay in
    // construction rather than moving to init() because the clear must run for
    // every Effect and derived init() overrides do not chain to a base init().
    clear_buffers();
    clip.w = W;
    clip.h = H;
    clip.y_end = H;
    clip.x_end = W;
  }

  /**
   * @brief Destroys the Effect instance.
   * @details Clears the single-live-Effect guard so the next construction
   * (effect switch) is admitted. Ordering is safe because every effect-swap
   * path destroys the outgoing instance before building its replacement.
   */
  virtual ~Effect() { s_alive = false; };

  /**
   * @brief Post-construction initialization. Override to move heavy
   * setup logic here (avoids GCC C1/C2 constructor duplication).
   */
  virtual void __attribute__((noinline)) init() {}

  /**
   * @brief Abstract method to be implemented by derived classes to generate a
   * frame of graphics data.
   */
  virtual void draw_frame() = 0;
  /**
   * @brief POV-display strobe control: whether each LED column is blanked to
   *        black immediately after it is shown, instead of persisting on the
   *        strip until the next column is drawn.
   *
   * Read ONLY on the hardware column-draw path (hardware/pov_single.h,
   * hardware/pov_segmented.h), once per swept column. It governs inter-column
   * behavior on the spinning strip, NOT framebuffer contents: it never clears
   * or persists the canvas — that is `persist_pixels` (see advance_buffer),
   * an orthogonal mechanism the two were historically conflated with.
   *
   * @return true  to strobe: after a column is lit the strip is re-shown black
   *               (FastLED.showColor / a trailing all-black DMA frame), so each
   *               column reads as a sharp slice with dark inter-column gaps.
   * @return false to persist: the lit column is held on the strip until the
   *               next column overwrites it, filling its full angular cell.
   */
  virtual bool strobe_columns() const = 0;

  /**
   * @brief Whether this effect must render the FULL canvas per simulator worker
   *        rather than be clipped to a segment band.
   * @return True for a cross-segment stateful effect (its per-frame state reads
   *         pixels outside the band — e.g. MeshFeedback's unbounded warp);
   *         false (the default) for every effect whose output in a band depends
   *         only on that band, which keeps segmented rendering's clipping win.
   * @details Read by the WASM segment driver (targets/wasm/wasm.cpp setClip) to
   *          leave the clip at full canvas for stateful effects. The device path
   *          already renders full-frame per board, so this is simulator-only.
   *          An effect derives the answer from its filter pipeline's compile-time
   *          `any_crosses_segments` fold (see core/filter.h and
   *          docs/segmented_stateful_effects_spec.md). The base default is false;
   *          only effects with a cross-segment filter override it.
   */
  [[nodiscard]] virtual bool needs_full_frame() const { return false; }

  /** @brief Segment clip region (display + render margin). */
  ClipRegion clip;

  /** @brief Accumulated rasterization time for the current frame (µs).
   *  WASM-only telemetry: only emscripten builds stamp and read it (the
   *  ScopedRenderTimer and wasm bridge), so the field and its accessors below
   *  are compiled out on the device rather than carried as dead storage. */
#ifdef __EMSCRIPTEN__
  double render_us = 0.0;
#endif

  /**
   * @brief Driver sets display bounds (which segment this Teensy owns).
   * @param y0 Inclusive start row of the owned segment.
   * @param y1 Exclusive end row of the owned segment.
   * @param x0 Inclusive start column of the owned segment.
   * @param x1 Exclusive end column of the owned segment.
   */
  void set_clip(int y0, int y1, int x0, int x1) {
    clip.y_start = y0;
    clip.y_end = y1;
    clip.x_start = x0;
    clip.x_end = x1;
  }
  /**
   * @brief Update only the horizontal clip band, leaving the y bounds intact.
   *
   * No device driver calls this — the firmware sets the full segment once via
   * set_clip(). It exists for callers that retune just the x-band (the web
   * simulator's segment-mode bridge writes the clip through set_clip(); the
   * canvas unit tests exercise this narrowing path directly).
   * @param x0 Inclusive start column of the horizontal clip band.
   * @param x1 Exclusive end column of the horizontal clip band.
   */
  void set_clip_x(int x0, int x1) {
    clip.x_start = x0;
    clip.x_end = x1;
  }
  /**
   * @brief Effect sets render margin for stateful filters.
   * @param m Render margin width in pixels.
   * @details Cold setup-time guard: ClipRegion's cylindrical wrap
   *          (render_x_start/render_x_end) only corrects a single period of
   *          underflow, so its documented [0, w) contract holds only while
   *          margin < w. Trap a margin that would wrap past the seam here, at
   *          configuration time, rather than letting a negative column leak into
   *          the per-fragment clip predicates on the hot path.
   */
  void set_margin(int m) {
    HS_CHECK(m < clip.w, "render margin must be < canvas width");
    clip.margin = m;
  }

  /**
   * @brief Retrieves the color of a pixel from the currently displayed buffer.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return The Pixel color reference.
   */
  virtual const Pixel &get_pixel(int x, int y) const {
    // Debug-only bounds guard, matching the write-path accessors
    // (operator()/prev): stripped on the device (NDEBUG, "No bounds checking"),
    // so it adds nothing to the ISR/WASM readback in release, but catches an
    // out-of-range display read in the test/sim build the same way the sibling
    // paths do.
    assert(x >= 0 && x < width_ && y >= 0 && y < height_);
    return bufs_[prev_.load(std::memory_order_relaxed)][y * width_ + x];
  }

  /**
   * @brief Raw base pointer to the currently displayed buffer (row-major,
   *        `width()` stride): `display_buffer()[y * width() + x]` equals
   *        `get_pixel(x, y)` for any effect that does not override `get_pixel`.
   *
   * ISR fast path: lets a column loop index pixels directly and skip the
   * per-pixel virtual `get_pixel` dispatch. Valid only until the next
   * `advance_display()` flip — re-fetch after one. Effects that override
   * `get_pixel` (e.g. the RingTwist scroller) must NOT use this, as it bypasses
   * their transform; check `overrides_get_pixel()` first.
   */
  [[nodiscard]] const Pixel *display_buffer() const {
    return bufs_[prev_.load(std::memory_order_relaxed)];
  }

  /**
   * @brief Whether this effect overrides get_pixel with a per-pixel transform
   *        that display_buffer() does NOT reflect.
   *
   * Base effects return false: `display_buffer()[y * width() + x]` equals
   * `get_pixel(x, y)`, so ISR fast paths may index the buffer directly. The
   * RingTwist scroller reads through a per-row offset in its get_pixel override
   * and returns true here, forcing those fast paths to fall back to virtual
   * dispatch. One cheap virtual call per column gates the bulk loop.
   */
  [[nodiscard]] virtual bool overrides_get_pixel() const { return false; }

  /**
   * @brief Gets the width of the effect.
   * @return The width.
   */
  [[nodiscard]] inline int width() const { return width_; }
  /**
   * @brief Gets the height of the effect.
   * @return The height.
   */
  [[nodiscard]] inline int height() const { return height_; }
  /**
   * @brief Checks whether the queued frame has been picked up for display.
   * @return True when `prev_ == next_` (no frame still waiting to be shown), so
   *         the writer is free to claim the other buffer.
   */
  [[nodiscard]] inline bool buffer_free() const {
    return prev_.load(std::memory_order_relaxed) ==
           next_.load(std::memory_order_relaxed);
  }
  /**
   * @brief Advances the display buffer pointer to the next queued frame.
   */
  inline void advance_display() {
    prev_.store(next_.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
  }
  /**
   * @brief Advances the drawing buffer pointer to the next available buffer.
   * @details If `persist_pixels` is true, copies the previous frame's content
   * to the new buffer.
   */
  inline void advance_buffer() {
    int c = cur_.load(std::memory_order_relaxed) ? 0 : 1;
    // The new write buffer must not be the one the ISR is currently scanning out
    // (prev_). With two physical buffers this holds only if buffer_free() gated
    // the advance; a caller that skips that gate would aim writes at the live
    // display buffer. Trap it here (once per frame, cold) instead of tearing.
    HS_CHECK(c != prev_.load(std::memory_order_relaxed));
    cur_.store(c, std::memory_order_relaxed);
    if (persist_pixels) {
      memcpy(bufs_[c], bufs_[prev_.load(std::memory_order_relaxed)],
             sizeof(Pixel) * width_ * height_);
    }
  }

  /**
   * @brief Queues the newly drawn frame to be displayed.
   */
  inline void queue_frame() {
    hs::disable_interrupts();
    next_.store(cur_.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
    hs::enable_interrupts();
  }

  /**
   * @brief Defines a runtime-adjustable parameter.
   */
  struct ParamDef {
    const char *name; /**< Parameter name. */
    std::variant<float *, bool *>
        target;             /**< Type-safe pointer to the variable. */
    float min = 0;          /**< Minimum value (for floats). */
    float max = 1;          /**< Maximum value (for floats). */
    float defaultValue = 0; /**< Default value. */
    bool animated = false;  /**< True if an animation drives this member; the GUI
                               surfaces these as auto-pausing sliders. */
    bool readonly = false;  /**< True if this is engine-written telemetry; the
                               GUI shows it live but disables editing. */

    /**
     * @brief Read the current value as float (bool maps to 0/1).
     * @return The target's value as a float; a bool target yields 1.0 or 0.0.
     */
    float get() const {
      return std::visit(
          [](auto *p) -> float {
            using T = std::remove_pointer_t<decltype(p)>;
            if constexpr (std::is_same_v<T, bool>)
              return *p ? 1.0f : 0.0f;
            else
              return *p;
          },
          target);
    }

    /**
     * @brief Write a float value (bool threshold at 0.5).
     * @param v Value to store; a bool target is set true when v > 0.5.
     * @warning Raw write: applies no readonly/finite/[min,max] gate. The
     * readonly-reject, non-finite-reject and range-clamp contract lives solely
     * in Effect::updateParameter — any value write originating outside trusted
     * engine code (UI bridge, animation driver) must route through there, not
     * call set() on a handle from ParamList::find().
     */
    void set(float v) {
      std::visit(
          [v](auto *p) {
            using T = std::remove_pointer_t<decltype(p)>;
            if constexpr (std::is_same_v<T, bool>)
              *p = (v > 0.5f);
            else
              *p = v;
          },
          target);
    }

    /**
     * @brief Check if this parameter targets a bool.
     * @return True if the target is a bool pointer, false if a float pointer.
     */
    bool is_bool() const { return std::holds_alternative<bool *>(target); }
  };

  /**
   * @brief Fixed-capacity registry of an effect's runtime parameters.
   * @details Stack-allocated array (no heap) to uphold the WASM no-realloc
   * memory-view invariant; capacity 32 is enforced at registration time.
   */
  struct ParamList {
    // Effect is the sole trusted mutator. registerParam fills the array and
    // markAnimated/markReadonly/updateParameter reach the writable handles in the
    // private section below; befriending Effect keeps those engine methods able
    // to mutate. Every OTHER caller — effect subclasses, and the WASM bridge via
    // the const getParameters() — sees only the public const overloads and so
    // gets a read-only ParamDef. A const handle cannot call the non-const
    // ParamDef::set(), so updateParameter (with its readonly/finite/[min,max]
    // gate) stays the single value-write path: the mutable bypass handle is no
    // longer reachable from outside the engine's own trusted methods.
    friend class Effect;

    std::array<ParamDef, 32> elements; /**< Fixed-capacity backing storage. */
    size_t count = 0;                  /**< Number of registered parameters. */

    /**
     * @brief Const iterator to the first registered parameter.
     * @return Pointer to the first element.
     */
    const ParamDef *begin() const { return elements.data(); }
    /**
     * @brief Const one-past-the-end iterator over registered parameters.
     * @return Pointer just past the last registered element.
     */
    const ParamDef *end() const { return elements.data() + count; }
    /**
     * @brief Looks up a registered parameter by name (the public, read-only
     * lookup).
     * @param name Parameter name to match (exact string compare).
     * @return Const pointer to the matching parameter, or nullptr if not found.
     */
    const ParamDef *find(const char *name) const {
      for (size_t i = 0; i < count; ++i) {
        if (std::strcmp(elements[i].name, name) == 0)
          return &elements[i];
      }
      return nullptr;
    }
    /**
     * @brief Number of registered parameters.
     * @return The count of registered parameters.
     */
    size_t size() const { return count; }

  private:
    // Writable accessors — reachable only by the friended Effect (see the note at
    // the top of the struct). They bypass updateParameter's write gate by design:
    // registerParam builds the list through `elements`/`count` and
    // markAnimated/markReadonly flip flags via the mutable find(). Kept private so
    // a UI/animation value write cannot reach a writable handle here and must
    // route through updateParameter, keeping the value contract single-sourced.
    // See ParamDef::set().
    ParamDef *begin() { return elements.data(); }
    ParamDef *end() { return elements.data() + count; }
    ParamDef *find(const char *name) {
      for (size_t i = 0; i < count; ++i) {
        if (std::strcmp(elements[i].name, name) == 0)
          return &elements[i];
      }
      return nullptr;
    }
  };

  // Parameter System
  /**
   * @brief Updates a parameter's value by name.
   * @param name The name of the parameter.
   * @param value The new value (mapped to bool if necessary).
   * @return true if the value was applied; false if the name is unknown (a
   *         stale/typo'd UI string), the parameter is readonly (engine-written
   *         telemetry), or the value was rejected as non-finite. The WASM bridge
   *         propagates this so the frontend can detect a no-op rather than
   *         silently dropping the write.
   */
  bool updateParameter(const char *name, float value) {
    auto *def = parameters.find(name);
    if (def == nullptr)
      return false;
    // This is the (untrusted) JS boundary — the only caller is WASM
    // setParameter. Reject writes the GUI is not allowed to make: readonly
    // params are engine-written telemetry the GUI shows live but disables
    // editing, so a write here is a stale/malicious frontend — drop it rather
    // than trust the client. Then reject non-finite input outright (a NaN/Inf
    // would silently poison render math), and clamp floats to the registered
    // [min,max] the GUI advertises. Bools are not range-clamped (set()
    // thresholds them at 0.5).
    if (def->readonly)
      return false;
    if (!std::isfinite(value))
      return false;
    if (!def->is_bool())
      value = hs::clamp(value, def->min, def->max);
    def->set(value);
    return true;
  }

  /**
   * @brief Retrieves the list of registered parameters.
   * @return Const reference to the parameter list.
   */
  const ParamList &getParameters() const { return parameters; }

  /**
   * @brief Pause/resume the effect's parameter-driving animations.
   * @details Wired to the `Mutation`/`Driver` gate via `anims_paused_`. Paused,
   * those animations freeze and the GUI slider bound to the same member is the
   * sole writer, so a user edit holds; resuming hands the member back to the
   * animation. Ambient motion (rotation/camera/palette) is not gated and keeps
   * running. The GUI surfaces this as the standard "Pause Animation" toggle.
   * @param paused True to freeze parameter-driving animations, false to resume.
   */
  void setAnimationsPaused(bool paused) { anims_paused_ = paused; }
  /**
   * @brief Reports whether parameter-driving animations are paused.
   * @return True if those animations are currently frozen.
   */
  bool animationsPaused() const { return anims_paused_; }

protected:
  /**
   * @brief Flag indicating if the previous frame's pixels should be copied to
   * the new buffer (for trails/decay).
   */
  bool persist_pixels;
  ParamList parameters; /**< List of parameters. */
  bool anims_paused_ = false; /**< Pause gate for parameter-driving animations;
                                 pass `&anims_paused_` to Mutation/Driver. */

  /**
   * @brief Flag a registered param as animation-driven.
   * @details The GUI renders flagged params as auto-pausing sliders (touching
   * one engages "Pause Animation"). Call after the matching registerParam.
   */
  void markAnimated(const char *name) {
    auto *def = parameters.find(name);
    // A misspelled name silently no-ops, leaving the param un-flagged and the
    // "Pause Animation" gate broken with no diagnostic. Trap instead — this is
    // cold setup code, the fail-fast doctrine's intended home.
    HS_CHECK(def, "markAnimated: unknown parameter name");
    def->animated = true;
  }

  /**
   * @brief Flag a registered param as engine-written telemetry (read-only).
   * @details The GUI keeps showing its live value but disables editing. Use for
   * output-only values clobbered every frame (e.g. an active-particle count).
   */
  void markReadonly(const char *name) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "markReadonly: unknown parameter name");
    def->readonly = true;
  }

  /**
   * @brief Registers a floating-point parameter.
   * @param name The name to expose.
   * @param ptr Pointer to the float variable.
   * @param min Minimum value.
   * @param max Maximum value.
   */
  void registerParam(const char *name, float *ptr, float min = 0.0f,
                     float max = 1.0f) {
    // Overflowing the fixed ParamList is an effect-authoring bug (too many
    // params); silently dropping the registration hides it and desyncs the GUI.
    // Trap instead. (Also upholds the WASM no-realloc memory-view invariant.)
    HS_CHECK(parameters.count < parameters.elements.size(),
             "registerParam: exceeded ParamList capacity");
    // A duplicate name is an authoring bug: find() returns the FIRST match, so a
    // second registration silently shadows (its slot is unreachable by name and
    // any UI/animation write lands on the first). Trap at this cold init seam.
    HS_CHECK(parameters.find(name) == nullptr,
             "registerParam: duplicate parameter name");
    // An inverted range feeds hs::clamp(value, min, max) with lo > hi
    // (implementation-defined), pinning the slider to garbage. Trap the
    // authoring bug at this cold init seam, like the capacity/duplicate guards.
    HS_CHECK(min <= max, "registerParam: min must be <= max");
    // The captured default is *ptr verbatim — registration never clamps it. Every
    // later updateParameter clamps into [min,max], so a default authored outside
    // the range would advertise an out-of-range value that snaps on the first GUI
    // edit. Trap that authoring mismatch at this cold init seam rather than ship a
    // self-inconsistent default. (NaN fails the comparison and traps too.)
    HS_CHECK(*ptr >= min && *ptr <= max,
             "registerParam: default *ptr outside [min,max]");
    parameters.elements[parameters.count++] = {name, ptr, min, max, *ptr};
  }

  /**
   * @brief Registers a boolean parameter.
   * @param name The name to expose.
   * @param ptr Pointer to the bool variable; its current value becomes the GUI
   *   default. Registration never mutates the target — symmetric with the float
   *   overload, which likewise captures `*ptr` and leaves it untouched.
   */
  void registerParam(const char *name, bool *ptr) {
    HS_CHECK(parameters.count < parameters.elements.size(),
             "registerParam: exceeded ParamList capacity");
    // Duplicate name guard, see the float overload.
    HS_CHECK(parameters.find(name) == nullptr,
             "registerParam: duplicate parameter name");
    parameters.elements[parameters.count++] = {name, ptr, 0.0f, 1.0f,
                                               (float)*ptr};
  }

  /**
   * @brief Registers a float param and flags it animation-driven in one call.
   * @details Convenience for the registerParam + markAnimated pair, so the name
   * literal is written once instead of twice (a typo in the second copy would
   * silently leave the param un-flagged — a dead-slider lint failure rather than
   * a compile error). Flags the just-registered element directly, skipping the
   * name lookup markAnimated would redo.
   */
  void registerAnimatedParam(const char *name, float *ptr, float min = 0.0f,
                             float max = 1.0f) {
    registerParam(name, ptr, min, max);
    parameters.elements[parameters.count - 1].animated = true;
  }

  /**
   * @brief Registers a float param and flags it engine-written telemetry in one
   * call.
   * @details Convenience for the registerParam + markReadonly pair; see
   * registerAnimatedParam for the single-source-the-literal rationale.
   */
  void registerReadonlyParam(const char *name, float *ptr, float min = 0.0f,
                             float max = 1.0f) {
    registerParam(name, ptr, min, max);
    parameters.elements[parameters.count - 1].readonly = true;
  }

private:
  /**
   * @brief Points bufs_ at the shared static storage and zeroes both buffers.
   * @details noinline so the two full-frame fills are emitted once instead of
   * being inlined into both GCC constructor variants (C1/C2). Invoked from the
   * ctor rather than init() because the clear must run for every Effect and
   * derived init() overrides do not chain to a base Effect::init().
   */
  void __attribute__((noinline)) clear_buffers() {
    bufs_[0] = buffer_a;
    std::fill_n(bufs_[0], MAX_W * MAX_H, Pixel(0, 0, 0));
    bufs_[1] = buffer_b;
    std::fill_n(bufs_[1], MAX_W * MAX_H, Pixel(0, 0, 0));
  }

  std::atomic<int> prev_{0}; /**< Buffer the ISR is currently reading. */
  std::atomic<int> cur_{0};  /**< Buffer the main loop is currently writing. */
  std::atomic<int> next_{0}; /**< Last completed frame, queued for display. */
  int width_;                /**< The width of the effect. */
  int height_;               /**< The height of the effect. */
  // Shared static storage for the double buffer. PRECONDITION: at most one
  // Effect may be live at a time — every instance points bufs_ at these same two
  // arrays, so a second concurrent Effect would corrupt both frames and the
  // prev_/cur_/next_ indices would alias a shared buffer. The Effect ctor/dtor
  // enforce this with the s_alive guard below; there is no RAM for per-instance
  // buffers (or a third buffer — see Canvas).
  static DMAMEM Pixel
      buffer_a[MAX_W * MAX_H]; /**< Static storage for buffer A (shared). */
  static DMAMEM Pixel
      buffer_b[MAX_W * MAX_H]; /**< Static storage for buffer B (shared). */
  Pixel *bufs_[2]; /**< Pointers to the two buffer storage locations. */
  // True while an Effect is constructed-but-not-destroyed. Guards the
  // single-live-Effect precondition on the shared buffer_a/buffer_b: the ctor
  // traps if it is already set, the dtor clears it. Cold path (effect switch
  // only), so the load/store is free.
  static bool s_alive;
};

/**
 * @brief Context class providing a safe, scoped interface to the current
 * drawing buffer.
 * @details Ensures the buffer advances correctly when the object is constructed
 * and queued when destroyed.
 */
class Canvas {
public:
  /**
   * @brief Constructs the Canvas, advancing the effect buffer and optionally
   * clearing it.
   * @param effect The effect instance owning the buffer.
   */
  Canvas(Effect &effect) : effect_(effect) {
    // NOT a TOCTOU race, despite the check-then-advance running outside an
    // interrupts-disabled window. Single-core Teensy, strict index ownership:
    // the main loop solely writes cur_/next_ (advance_buffer/queue_frame); the
    // ISR solely writes prev_ (advance_display: prev_ = next_). buffer_free()
    // (prev_ == next_) can only be *falsified* by the main loop's next
    // queue_frame() — this same thread, after this Canvas is destroyed. The
    // ISR's advance_display is idempotent once prev_ == next_, so it cannot
    // invalidate the gate between the check and the flip. advance_buffer() then
    // moves cur_ to the OTHER of the two buffers, which — since prev_ == next_
    // pins one — is guaranteed != prev_: the write buffer is never the ISR's
    // read buffer. Relaxed atomics suffice (the ISR preempts the main loop =
    // same observer; no multi-core reordering to fence). The 2-buffer + this
    // gate is correct by design — there is no RAM for a third buffer, and the
    // gate (not a third buffer) is what prevents tearing.
    //
    // Watchdog: buffer_free() is re-satisfied only when the display ISR runs
    // advance_display() at a frame boundary — at most once per revolution. An
    // unbounded wait here means that ISR stopped advancing (stalled rotor,
    // detached timer, priority inversion): a silent headless freeze, which the
    // fail-fast doctrine says must trap, not hang. The bound is far above one
    // revolution at any sane RPM, so a slow spin-up never false-trips (the
    // first frame doesn't wait at all — prev_==next_==0 at construction). On
    // timeout, route through HS_CHECK so the failure flows through the project's
    // single check_fail() trap path (log + flush + file:line breadcrumb +
    // __builtin_trap) rather than open-coding it; check_fail's breadcrumb names
    // the site, so the trap is no longer a bare illegal instruction with no
    // message. The outer buffer_free() guard skips the whole block — and thus
    // every micros() read — on the common no-wait path, keeping this zero-cost
    // when the buffer is already free. On the slow path, wait_start is sampled
    // once at wait entry; the deadline is then re-evaluated via micros() inside
    // the spin on every iteration, so the bound measures from wait entry.
    if (!effect_.buffer_free()) {
      const unsigned long wait_start = micros();
      while (!effect_.buffer_free()) {
        HS_CHECK(micros() - wait_start < kBufferFreeWatchdogUs,
                 "buffer_free watchdog timeout — display ISR stalled");
      }
    }
    effect_.advance_buffer();
    if (!effect_.persist_pixels) {
      clear_buffer();
    }
  }

  /**
   * @brief Destructor. Queues the finished frame to be displayed.
   */
  ~Canvas() { effect_.queue_frame(); }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Reference to the Pixel.
   */
  inline Pixel &operator()(int x, int y) {
    // Use assert here, NOT HS_CHECK — this is the one hot-loop exception. The
    // bounds guard is debug-only by design (stripped on the device, as the
    // README documents — "No bounds checking"); an always-on HS_CHECK branch on
    // every pixel access is the single place HS_CHECK's own contract forbids it.
    assert(x >= 0 && x < effect_.width_ && y >= 0 && y < effect_.height_);
    return effect_.bufs_[effect_.cur_.load(std::memory_order_relaxed)]
                        [y * effect_.width_ + x];
  }

  /**
   * @brief Accesses a pixel in the previous drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Copy of the Pixel from the previous frame.
   * @details Returns by value, unlike `operator()`'s by-reference access to the
   *          current buffer: the previous frame is read-only (the next frame is
   *          composed in the current buffer), so handing out a reference into it
   *          would invite an accidental write to a buffer that is about to be
   *          recycled. A Pixel is small enough that the copy is free.
   */
  inline Pixel prev(int x, int y) const {
    assert(x >= 0 && x < effect_.width_ && y >= 0 && y < effect_.height_);
    return effect_.bufs_[effect_.prev_.load(std::memory_order_relaxed)]
                        [y * effect_.width_ + x];
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 1D index.
   * @param xy The 1D index.
   * @return Reference to the Pixel.
   */
  inline Pixel &operator()(int xy) {
    assert(xy >= 0 && xy < effect_.width_ * effect_.height_);
    return effect_.bufs_[effect_.cur_.load(std::memory_order_relaxed)][xy];
  }

  /**
   * @brief Clears the entire current drawing buffer to black.
   */
  void clear_buffer() {
    int c = effect_.cur_.load(std::memory_order_relaxed);
    std::fill_n(effect_.bufs_[c], effect_.width_ * effect_.height_,
                Pixel(0, 0, 0));
  }

  /**
   * @brief Gets the width of the underlying effect.
   * @return The width in pixels.
   */
  [[nodiscard]] inline int width() const { return effect_.width(); }
  /**
   * @brief Gets the height of the underlying effect.
   * @return The height in pixels.
   */
  [[nodiscard]] inline int height() const { return effect_.height(); }
  /**
   * @brief Gets the effect's current clip region.
   * @return Const reference to the clip region (display + render margin).
   */
  [[nodiscard]] inline const ClipRegion &clip() const { return effect_.clip; }

#ifdef __EMSCRIPTEN__
  // WASM-only render-time telemetry (see Effect::render_us). Compiled out on the
  // device, where there is no JS perf clock and nothing reads the counter.
  /**
   * @brief Resets the accumulated rasterization time for the current frame.
   */
  inline void reset_render_us() { effect_.render_us = 0.0; }
  /**
   * @brief Adds to the accumulated rasterization time for the current frame.
   * @param us Elapsed time to add, in microseconds.
   */
  inline void add_render_us(double us) { effect_.render_us += us; }
  /**
   * @brief Gets the accumulated rasterization time for the current frame.
   * @return The accumulated render time, in microseconds.
   */
  [[nodiscard]] inline double get_render_us() const { return effect_.render_us; }
#endif
  /**
   * @brief Checks if debug visuals are enabled.
   * @return True if debugging is active.
   */
  inline bool debug() const { return effect_.debug_visuals; }

private:
  /** Watchdog bound for the ctor buffer_free() spin (µs). One display
   *  revolution is tens-to-hundreds of ms even at low RPM; 2 s is well above
   *  that, so only a genuinely stalled display ISR trips it. */
  static constexpr unsigned long kBufferFreeWatchdogUs = 2000000UL;
  Effect &effect_; /**< Reference to the owning Effect instance. */
};

