/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// platform.h first: on device it defines NDEBUG, which must be set before
// <cassert> expands the assert macro — otherwise assert-stripping would depend
// on a prior TU having pulled in platform.h, making this header non-self-sufficient.
#include "engine/platform.h"
#include <cstring>
#include <cassert>
#include <algorithm>
#include <variant>
#include <atomic>
#include <utility>
#include "engine/constants.h"
#include "color/color.h"
#include <array>

class Canvas;

/**
 * @brief Construction-time flags for an Effect (see the accessors of the same
 *        name). Defaults suit a plain, non-strobing, band-clippable effect.
 */
struct EffectConfig {
  bool strobe = false;     /**< POV column strobe (Effect::strobe_columns). */
  bool persist = false;    /**< Copy previous frame forward (persists_pixels). */
  bool full_frame = false; /**< Force full-canvas render (needs_full_frame). */
};

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
   * @param cfg Construction-time flags (strobe / persist / full-frame); see
   *        EffectConfig. Defaults to a plain band-clippable effect.
   */
  Effect(int W, int H, EffectConfig cfg = {})
      : persist_pixels(cfg.persist), full_frame(cfg.full_frame),
        strobe(cfg.strobe), width_(W), height_(H) {
    // Single-live-Effect precondition: every Effect aliases the same two static
    // buffers (buffer_a/buffer_b), so a second live instance corrupts both frames.
    HS_CHECK(!s_alive,
             "Effect: a second Effect was constructed while one is still alive; "
             "buffer_a/buffer_b are shared static storage (one live Effect only)");
    s_alive = true;
    // Point bufs_ at the shared static storage and clear both buffers.
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
  virtual ~Effect() { s_alive = false; }

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
   * Read on the hardware column-draw path once per swept column. Governs
   * inter-column strip behavior, not framebuffer contents (that is
   * `persist_pixels`).
   *
   * @return true  to strobe: after a column is lit the strip is re-shown black
   *               (FastLED.showColor / a trailing all-black DMA frame), so each
   *               column reads as a sharp slice with dark inter-column gaps.
   * @return false to persist: the lit column is held on the strip until the
   *               next column overwrites it, filling its full angular cell.
   */
  [[nodiscard]] bool strobe_columns() const { return strobe; }

  /**
   * @brief Whether this effect must render the FULL canvas per simulator worker
   *        rather than be clipped to a segment band.
   * @return True for a cross-segment stateful effect (its per-frame state reads
   *         pixels outside the band — e.g. MeshFeedback's unbounded warp);
   *         false (the default) for every effect whose output in a band depends
   *         only on that band, which keeps segmented rendering's clipping win.
   * @details Read by the WASM segment driver (targets/wasm/wasm.cpp setClip) and
   *          the device driver (hardware/pov_segmented.h clip_to_segment) to leave
   *          the clip at full canvas for stateful effects. Set once at
   *          construction from the effect's filter pipeline compile-time
   *          `any_crosses_segments` fold (see core/render/filter.h and
   *          docs/segmented_stateful_effects_spec.md); defaults to false.
   */
  [[nodiscard]] bool needs_full_frame() const { return full_frame; }

  /**
   * @brief Whether this effect copies its previous frame forward (trails/decay).
   * @return True when persist_pixels is set.
   * @details The segmented device driver leaves a persisting effect at full
   *          canvas: its per-frame arm-half alternation would break trail
   *          continuity if the render were clipped to one quadrant.
   */
  [[nodiscard]] bool persists_pixels() const { return persist_pixels; }

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
    // Guard the render-bound helpers' invariant: non-inverted band, non-negative
    // origins. The upper bound (x1 <= w) is checked downstream at the LUT-domain
    // HS_CHECK in Scan::Shader::draw, not here.
    HS_CHECK(y0 >= 0 && y0 <= y1 && x0 >= 0 && x0 <= x1,
             "set_clip band must be non-inverted with non-negative origins");
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
    // Same non-inverted/non-negative invariant as set_clip (x-band only); the
    // upper bound is likewise deferred to the downstream LUT-domain guard.
    HS_CHECK(x0 >= 0 && x0 <= x1,
             "set_clip_x band must be non-inverted with a non-negative origin");
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
    HS_CHECK(m >= 0 && m < clip.w,
             "render margin must be in [0, canvas width)");
    clip.margin = m;
  }

  /**
   * @brief Retrieves the color of a pixel from the currently displayed buffer.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return The Pixel color reference.
   * @note An override that applies a per-pixel transform (rather than a plain
   *       buffer read) must also override overrides_get_pixel() to return true;
   *       otherwise ISR/readback fast paths that index display_buffer() directly
   *       bypass the transform.
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
      // The trail base is the last COMPLETED frame (next_), not the one the ISR
      // is currently displaying (prev_). The buffer_free() gate that must precede
      // every advance forces prev_ == next_, so the two coincide here — copy from
      // next_ and assert the equality, rather than read prev_ and depend on the
      // gate silently across methods.
      int last = next_.load(std::memory_order_relaxed);
      HS_CHECK(last == prev_.load(std::memory_order_relaxed));
      memcpy(bufs_[c], bufs_[last], sizeof(Pixel) * width_ * height_);
    }
  }

  /**
   * @brief Queues the newly drawn frame to be displayed.
   * @details Publishes `cur_` as the new `next_`. Correctness relies on
   * `hs::disable_interrupts()` being a compiler barrier (single-core target), not
   * on the relaxed atomics: the drawn buffer is fully written before the flip is
   * observable to the ISR.
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
    // Effect is the sole trusted mutator; the writable accessors below are
    // private so every other caller sees only the const overloads and must route
    // value writes through updateParameter.
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
    // Writable accessors, reachable only by the friended Effect (see the note at
    // the top of the struct). Kept private so value writes route through
    // updateParameter.
    ParamDef *begin() { return elements.data(); }
    ParamDef *end() { return elements.data() + count; }
    ParamDef *find(const char *name) {
      return const_cast<ParamDef *>(std::as_const(*this).find(name));
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
    // Untrusted JS boundary (WASM setParameter): reject readonly-param writes and
    // non-finite input, and clamp floats to [min,max]. Bools are thresholded at
    // 0.5 by set(), not range-clamped.
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
  /**
   * @brief Full-canvas render gate (see needs_full_frame()); set once at
   * construction from the filter pipeline `any_crosses_segments` trait.
   */
  bool full_frame;
  /**
   * @brief POV column-strobe flag (see strobe_columns()); set at construction.
   */
  bool strobe;
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
    // Trap a default *ptr outside [min,max]: registration captures it verbatim but
    // every later updateParameter clamps, so it would snap on the first GUI edit.
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
  // Shared static storage for the double buffer. PRECONDITION: at most one Effect
  // live at a time (enforced by the s_alive guard); a second would alias these
  // arrays and the prev_/cur_/next_ indices.
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
    // Not a TOCTOU race: single-core, strict index ownership — the main loop
    // writes cur_/next_, the ISR only sets prev_ = next_ (idempotent once
    // prev_==next_), so the gate can't be falsified between check and flip.
    // advance_buffer() then picks the other buffer, guaranteed != prev_.
    //
    // Watchdog: buffer_free() is re-satisfied only when the display ISR advances
    // at a frame boundary; an unbounded wait means that ISR stalled, so trap
    // (fail-fast) rather than hang. The outer buffer_free() guard skips the
    // micros() reads on the common no-wait path.
    if (!effect_.buffer_free()) {
      const unsigned long wait_start = micros();
      while (!effect_.buffer_free()) {
        HS_CHECK(micros() - wait_start < BUFFER_FREE_WATCHDOG_US,
                 "buffer_free watchdog timeout — display ISR stalled");
#ifdef HS_TEST_BUILD
        s_buffer_free_spins.fetch_add(1, std::memory_order_relaxed);
#endif
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

  Canvas(const Canvas &) = delete;
  Canvas(Canvas &&) = delete;
  Canvas &operator=(const Canvas &) = delete;
  Canvas &operator=(Canvas &&) = delete;

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

#ifdef HS_TEST_BUILD
  /**
   * @brief Test-only count of buffer_free() spin iterations across all Canvas
   *        ctors in this process.
   * @return The running spin-iteration count.
   * @details Lets a test detect that a ctor has actually entered the wait loop
   *          and release it on observed progress, rather than racing a fixed
   *          sleep. Compiled out of the device/sim image.
   */
  static unsigned long buffer_free_spin_count() {
    return s_buffer_free_spins.load(std::memory_order_relaxed);
  }
#endif

private:
  /** Watchdog bound for the ctor buffer_free() spin (µs). One display
   *  revolution is tens-to-hundreds of ms even at low RPM; 2 s is well above
   *  that, so only a genuinely stalled display ISR trips it. */
  static constexpr unsigned long BUFFER_FREE_WATCHDOG_US = 2000000UL;
  Effect &effect_; /**< Reference to the owning Effect instance. */
#ifdef HS_TEST_BUILD
  inline static std::atomic<unsigned long> s_buffer_free_spins{0};
#endif
};

