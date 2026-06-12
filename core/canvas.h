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
    bufs_[0] = buffer_a;
    std::fill_n(bufs_[0], MAX_W * MAX_H, Pixel(0, 0, 0));
    bufs_[1] = buffer_b;
    std::fill_n(bufs_[1], MAX_W * MAX_H, Pixel(0, 0, 0));
    clip.w = W;
    clip.h = H;
    clip.y_end = H;
    clip.x_end = W;
  }

  virtual ~Effect() {};

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
   * @brief Indicates whether the display should show black between POV sweeps.
   * @return True if the background should be shown (often black).
   */
  virtual bool show_bg() const = 0;

  /** @brief Segment clip region (display + render margin). */
  ClipRegion clip;

  /** @brief Accumulated rasterization time for the current frame (µs, WASM only). */
  double render_us = 0.0;

  /** @brief Driver sets display bounds (which segment this Teensy owns). */
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
   */
  void set_clip_x(int x0, int x1) {
    clip.x_start = x0;
    clip.x_end = x1;
  }
  /** @brief Effect sets render margin for stateful filters. */
  void set_margin(int m) { clip.margin = m; }

  /**
   * @brief Retrieves the color of a pixel from the currently displayed buffer.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return The Pixel color reference.
   */
  virtual const Pixel &get_pixel(int x, int y) const {
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

    /** @brief Read the current value as float (bool maps to 0/1). */
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

    /** @brief Write a float value (bool threshold at 0.5). */
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

    /** @brief Check if this parameter targets a bool. */
    bool is_bool() const { return std::holds_alternative<bool *>(target); }
  };

  /**
   * @brief Fixed-capacity registry of an effect's runtime parameters.
   * @details Stack-allocated array (no heap) to uphold the WASM no-realloc
   * memory-view invariant; capacity 32 is enforced at registration time.
   */
  struct ParamList {
    std::array<ParamDef, 32> elements;
    size_t count = 0;

    const ParamDef *begin() const { return elements.data(); }
    const ParamDef *end() const { return elements.data() + count; }
    ParamDef *begin() { return elements.data(); }
    ParamDef *end() { return elements.data() + count; }

    ParamDef *find(const char *name) {
      for (size_t i = 0; i < count; ++i) {
        if (std::strcmp(elements[i].name, name) == 0)
          return &elements[i];
      }
      return nullptr;
    }
    const ParamDef *find(const char *name) const {
      for (size_t i = 0; i < count; ++i) {
        if (std::strcmp(elements[i].name, name) == 0)
          return &elements[i];
      }
      return nullptr;
    }
    size_t size() const { return count; }
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
   */
  void setAnimationsPaused(bool paused) { anims_paused_ = paused; }
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
    parameters.elements[parameters.count++] = {name, ptr, 0.0f, 1.0f,
                                               (float)*ptr};
  }

private:
  std::atomic<int> prev_{0}; /**< Buffer the ISR is currently reading. */
  std::atomic<int> cur_{0};  /**< Buffer the main loop is currently writing. */
  std::atomic<int> next_{0}; /**< Last completed frame, queued for display. */
  int width_;                /**< The width of the effect. */
  int height_;               /**< The height of the effect. */
  static DMAMEM Pixel
      buffer_a[MAX_W * MAX_H]; /**< Static storage for buffer A. */
  static DMAMEM Pixel
      buffer_b[MAX_W * MAX_H]; /**< Static storage for buffer B. */
  Pixel *bufs_[2]; /**< Pointers to the two buffer storage locations. */
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
    // timeout, name the site on Serial (the trap is a bare illegal instruction
    // with no message) then trap; a truncated Serial write is harmless since we
    // halt immediately. The outer buffer_free() guard skips the whole block —
    // and thus every micros() read — on the common no-wait path, keeping this
    // zero-cost when the buffer is already free. On the slow path, wait_start is
    // sampled once at wait entry; the deadline is then re-evaluated via micros()
    // inside the spin on every iteration, so the bound measures from wait entry.
    if (!effect_.buffer_free()) {
      const unsigned long wait_start = micros();
      while (!effect_.buffer_free()) {
        if (micros() - wait_start >= kBufferFreeWatchdogUs) {
          hs::log("FATAL: buffer_free watchdog timeout — display ISR stalled");
          __builtin_trap();
        }
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

  [[nodiscard]] inline int width() const { return effect_.width(); }
  [[nodiscard]] inline int height() const { return effect_.height(); }
  [[nodiscard]] inline const ClipRegion &clip() const { return effect_.clip; }

  inline void reset_render_us() { effect_.render_us = 0.0; }
  inline void add_render_us(double us) { effect_.render_us += us; }
  [[nodiscard]] inline double get_render_us() const { return effect_.render_us; }
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

