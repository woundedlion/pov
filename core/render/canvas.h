/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

// platform.h defines NDEBUG on device; include before <cassert> so assert
// stripping does not depend on include order.
#include "engine/platform.h"
#include <cstring>
#include <cassert>
#include <algorithm>
#include <atomic>
#include <type_traits>
#include <utility>
#include "engine/constants.h"
#include "color/color.h"
#include <array>

/**
 * @file canvas.h
 * @brief EffectConfig, the Effect base class, and the Canvas pixel buffer.
 */

class Canvas;

/**
 * @brief Outcome of a named parameter write (Effect::updateParameter, surfaced
 *        to JS by the WASM bridge's setParameter).
 * @details NO_EFFECT is produced only by the bridge, which may have no effect
 *          installed to forward the write to; Effect::updateParameter itself
 *          never returns it.
 */
enum class ParamSetResult {
  APPLIED,       /**< Value written (floats clamped to [min,max] first). */
  NO_EFFECT,     /**< No effect is installed to receive the write. */
  UNKNOWN_PARAM, /**< No registered parameter has this name. */
  READONLY,      /**< Parameter is engine-written telemetry. */
  NON_FINITE,    /**< Value is NaN or infinite. */
};

/**
 * @brief Construction-time flags for an Effect (see the accessors of the same
 *        name). Defaults suit a plain, non-strobing, band-clippable effect.
 */
struct EffectConfig {
  bool strobe = false;  /**< POV column strobe (Effect::strobe_columns). */
  bool persist = false; /**< Copy previous frame forward (persists_pixels). */
  bool full_frame = false; /**< Force full-canvas render (needs_full_frame). */
  bool reads_outside_band =
      false; /**< Clear pixels outside the display band. */
  /**
   * @brief Render-bound expansion the effect's filters need, in pixels
   *        (ClipRegion::margin). Raised to the ClipRegion default when lower.
   */
  int margin = ClipRegion{}.margin;
};

/**
 * @brief Base class for all visual effects.
 * @details Manages double buffering, persistence, and provides an interface for
 * drawing a frame.
 *
 * Naming split: methods on the JS/embind boundary (updateParameter,
 * getParameters, setAnimationsPaused) are camelCase to match the WASM bridge;
 * the internal C++ API (register_param, needs_full_frame, set_clip) is
 * snake_case.
 */
class Effect {
  friend class Canvas;

public:
  using BufferReadyHook = void (*)(Effect &);

  bool debug_visuals = false; /**< Flag to enable visual debugging overlays. */
#ifdef HS_TEST_BUILD
  /** @brief Forces the whole-buffer clear, for clip-clear parity tests. */
  bool force_full_buffer_clear = false;
#endif

  /**
   * @brief Constructs an Effect instance.
   * @param W The width (resolution) of the effect, in [1, MAX_W].
   * @param H The height (resolution) of the effect, in [1, MAX_H].
   * @param cfg Construction-time behavior flags; see EffectConfig. Defaults to
   *        a plain band-clippable effect.
   */
  HS_COLD_MEMBER Effect(int W, int H, EffectConfig cfg = {})
      : persist_pixels(cfg.persist), full_frame(cfg.full_frame),
        reads_outside_band(cfg.reads_outside_band), strobe(cfg.strobe),
        frame_width(W), frame_height(H) {
    HS_CHECK(W > 0 && W <= MAX_W && H > 0 && H <= MAX_H,
             "Effect dimensions %d x %d are outside 1..%d x 1..%d", W, H, MAX_W,
             MAX_H);
    // Single-live-Effect precondition: every Effect aliases the same two static
    // buffers, so a second live instance corrupts both frames.
    HS_CHECK(
        !s_alive,
        "Effect: a second Effect was constructed while one is still alive; "
        "buffer_a/buffer_b are shared static storage (one live Effect only)");
    s_alive = true;
    // Point bufs at the shared static storage and clear both buffers.
    clear_buffers();
    clip_region.w = W;
    clip_region.h = H;
    clip_region.y_end = H;
    clip_region.x_end = W;
    // Only widen: a filter fold of 0 must not shrink the coverage every effect
    // already gets from the ClipRegion default.
    if (cfg.margin > clip_region.margin)
      set_margin(cfg.margin);
  }

  /**
   * @brief Destroys the Effect instance.
   * @details Clears the single-live-Effect guard so the next construction is
   * admitted; every effect-swap path destroys the outgoing instance first.
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
   *        black immediately after it is shown. Governs inter-column strip
   *        behavior, not framebuffer contents (that is `persist_pixels`).
   * @return true to strobe (each column a sharp slice with dark gaps); false to
   *         persist the lit column until the next overwrites it.
   */
  [[nodiscard]] bool strobe_columns() const { return strobe; }

  /**
   * @brief Whether this effect must render the FULL canvas per simulator worker
   *        rather than be clipped to a segment band.
   * @return True for an effect whose output can cross segment boundaries (for
   *         example MeshFeedback's unbounded warp or World::Trails);
   *         false (the default) when each segment can render independently.
   * @details Read by the segment drivers to leave the clip at full canvas for
   *          stateful effects. Set once at construction from the filter
   *          pipeline's `any_crosses_segments` fold; defaults to false.
   */
  [[nodiscard]] bool needs_full_frame() const { return full_frame; }

  /**
   * @brief Whether this effect copies its previous frame forward (trails/decay).
   * @return True when persist_pixels is set.
   * @details The segmented device driver leaves a persisting effect at full
   *          canvas, else its per-frame arm-half alternation breaks trail
   *          continuity.
   */
  [[nodiscard]] bool persists_pixels() const { return persist_pixels; }

  /**
   * @brief The effect's current clip region (display band + render margin).
   * @return Const reference to the clip region; mutate only through set_clip(),
   *         set_clip_x(), and set_margin(), which enforce its invariants.
   */
  [[nodiscard]] const ClipRegion &clip() const { return clip_region; }

  /**
   * @brief Driver sets display bounds (which segment this Teensy owns).
   * @param y0 Inclusive start row of the owned segment.
   * @param y1 Exclusive end row of the owned segment.
   * @param x0 Inclusive start column of the owned segment.
   * @param x1 Exclusive end column of the owned segment.
   */
  void set_clip(int y0, int y1, int x0, int x1) {
    HS_CHECK(y0 >= 0 && y0 <= y1 && y1 <= clip_region.h && x0 >= 0 &&
                 x0 <= x1 && x1 <= clip_region.w,
             "set_clip band must be non-inverted and within canvas bounds");
    clip_region.y_start = y0;
    clip_region.y_end = y1;
    clip_region.x_start = x0;
    clip_region.x_end = x1;
  }
  /**
   * @brief Update only the horizontal clip band, leaving the y bounds intact.
   *
   * For callers that retune just the x-band (canvas unit tests); no device
   * driver calls this.
   * @param x0 Inclusive start column of the horizontal clip band.
   * @param x1 Exclusive end column of the horizontal clip band.
   */
  void set_clip_x(int x0, int x1) {
    HS_CHECK(x0 >= 0 && x0 <= x1 && x1 <= clip_region.w,
             "set_clip_x band must be non-inverted and within canvas width");
    clip_region.x_start = x0;
    clip_region.x_end = x1;
  }
  /**
   * @brief Effect sets render margin for stateful filters.
   * @param m Render margin width in pixels.
   * @details ClipRegion's cylindrical wrap only corrects a single period of
   *          underflow, so its [0, w) contract holds only while margin < w. Trap
   *          a wrapping margin here rather than let a negative column leak into
   *          the per-fragment clip predicates.
   */
  void set_margin(int m) {
    HS_CHECK(m >= 0 && m < clip_region.w,
             "render margin must be in [0, canvas width)");
    clip_region.margin = m;
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
    // Debug-only bounds guard, matching the write-path accessors (stripped on
    // device, catches an out-of-range display read in test/sim).
    assert(x >= 0 && x < frame_width && y >= 0 && y < frame_height);
    return bufs[prev.load(std::memory_order_relaxed)][y * frame_width + x];
  }

  /**
   * @brief Raw base pointer to the currently displayed buffer (row-major,
   *        `width()` stride): `display_buffer()[y * width() + x]` equals
   *        `get_pixel(x, y)` for any effect that does not override `get_pixel`.
   *
   * ISR fast path: index pixels directly, skipping the virtual `get_pixel`
   * dispatch. Valid only until the next `advance_display()` flip. Effects that
   * override `get_pixel` must NOT use this (it bypasses their transform); check
   * `overrides_get_pixel()` first.
   */
  [[nodiscard]] const Pixel *display_buffer() const {
    return bufs[prev.load(std::memory_order_relaxed)];
  }

  /**
   * @brief Whether this effect overrides get_pixel with a per-pixel transform
   *        that display_buffer() does NOT reflect.
   *
   * Base effects return false: `display_buffer()[y * width() + x]` equals
   * `get_pixel(x, y)`, so ISR fast paths may index the buffer directly. An
   * override that reads through a per-pixel transform (e.g. the RingTwist
   * scroller) returns true, forcing those paths back to virtual dispatch.
   */
  [[nodiscard]] virtual bool overrides_get_pixel() const { return false; }

  /**
   * @brief Gets the width of the effect.
   * @return The width.
   */
  [[nodiscard]] inline int width() const { return frame_width; }
  /**
   * @brief Gets the height of the effect.
   * @return The height.
   */
  [[nodiscard]] inline int height() const { return frame_height; }
  /**
   * @brief Checks whether the queued frame has been picked up for display.
   * @return True when `prev == next` (no frame still waiting to be shown), so
   *         the writer is free to claim the other buffer.
   */
  [[nodiscard]] inline bool buffer_free() const {
    return prev.load(std::memory_order_relaxed) ==
           next.load(std::memory_order_relaxed);
  }
  /**
   * @brief Installs a callback run once per frame from the Canvas constructor.
   * @param hook Callback to install, or null to disable it.
   * @details Runs after the buffer_free() wait returns and before
   * `advance_buffer()` and the stale-pixel clear. That slot is the only one that
   * sees the segment window the ISR settled on during the wait *and* still
   * precedes the clear, so it is where the segmented driver sets the clip that
   * governs which band gets cleared.
   */
  void set_buffer_ready_hook(BufferReadyHook hook) { buffer_ready_hook = hook; }
  /**
   * @brief Advances the display buffer pointer to the next queued frame.
   * @details The acquire fence pairs with `queue_frame()`'s release fence, so
   * the queued frame's pixel writes are visible before any read through
   * `prev`; the per-pixel loads themselves stay relaxed.
   */
  inline void advance_display() {
    int n = next.load(std::memory_order_relaxed);
    std::atomic_thread_fence(std::memory_order_acquire);
    prev.store(n, std::memory_order_relaxed);
  }
  /**
   * @brief Defines a runtime-adjustable parameter.
   */
  struct ParamDef {
    /** @brief Runtime representation of the parameter target. */
    enum class TargetType : uint8_t {
      FLOAT,
      BOOL,
      ENUM_I8,
      ENUM_U8,
      ENUM_I16,
      ENUM_U16,
      ENUM_I32,
      ENUM_U32,
    };

    const char *name = nullptr; /**< Parameter name. */
    void *target = nullptr;     /**< Pointer to the target variable. */
    const char *const *options = nullptr; /**< Option labels for an enumerated
                               param (GUI dropdown), or null for a plain param.
                               Must outlive the effect (string literals). */
    const char *const *export_options =
        nullptr;   /**< C++ enum literals indexed like options, or null. */
    float min = 0; /**< Minimum value (for floats). */
    float max = 1; /**< Maximum value (for floats). */
    int option_count = 0; /**< Number of labels; > 0 marks an enum target. */
    TargetType target_type = TargetType::FLOAT; /**< Target storage format. */
    bool animated = false; /**< True if an animation drives this member; the GUI
                               surfaces these as auto-pausing sliders. */
    bool readonly = false; /**< True if this is engine-written telemetry; the
                               GUI shows it live but disables editing. */
    bool preset = true; /**< Whether preset exports include this parameter. */

    template <typename Integer> float get_integer() const {
      Integer value;
      std::memcpy(&value, target, sizeof(value));
      return static_cast<float>(value);
    }

    template <typename Integer> void set_integer(float value) {
      const Integer stored = static_cast<Integer>(value);
      std::memcpy(target, &stored, sizeof(stored));
    }

    /**
     * @brief Read the current value as float (bool maps to 0/1).
     * @return The target's value as a float; a bool target yields 1.0 or 0.0.
     */
    float get() const {
      switch (target_type) {
      case TargetType::FLOAT:
        return *static_cast<const float *>(target);
      case TargetType::BOOL:
        return *static_cast<const bool *>(target) ? 1.0f : 0.0f;
      case TargetType::ENUM_I8:
        return get_integer<int8_t>();
      case TargetType::ENUM_U8:
        return get_integer<uint8_t>();
      case TargetType::ENUM_I16:
        return get_integer<int16_t>();
      case TargetType::ENUM_U16:
        return get_integer<uint16_t>();
      case TargetType::ENUM_I32:
        return get_integer<int32_t>();
      case TargetType::ENUM_U32:
        return get_integer<uint32_t>();
      }
      __builtin_unreachable();
    }

    /**
     * @brief Write a float value (bool threshold at 0.5).
     * @param v Value to store; a bool target is set true when v > 0.5.
     * @warning Raw write: applies no readonly/finite/[min,max] gate. That
     * contract lives solely in Effect::updateParameter — any write from outside
     * trusted engine code must route through there, not call set() on a handle
     * from ParamList::find().
     */
    void set(float v) {
      switch (target_type) {
      case TargetType::FLOAT:
        *static_cast<float *>(target) = v;
        return;
      case TargetType::BOOL:
        *static_cast<bool *>(target) = v > 0.5f;
        return;
      case TargetType::ENUM_I8:
        return set_integer<int8_t>(v);
      case TargetType::ENUM_U8:
        return set_integer<uint8_t>(v);
      case TargetType::ENUM_I16:
        return set_integer<int16_t>(v);
      case TargetType::ENUM_U16:
        return set_integer<uint16_t>(v);
      case TargetType::ENUM_I32:
        return set_integer<int32_t>(v);
      case TargetType::ENUM_U32:
        return set_integer<uint32_t>(v);
      }
      __builtin_unreachable();
    }

    /**
     * @brief Check if this parameter targets a bool.
     * @return True if the target is a bool pointer, false if a float pointer.
     */
    bool is_bool() const { return target_type == TargetType::BOOL; }

    /**
     * @brief Check if this parameter is enumerated (renders as a dropdown).
     * @return True if option labels are attached.
     */
    bool is_enum() const { return option_count > 0; }
  };
  static_assert(sizeof(void *) != 4 || sizeof(ParamDef) == 32,
                "ParamDef must keep its 32-bit device footprint");

  /**
   * @brief Fixed-capacity registry of an effect's runtime parameters.
   * @details Stack-allocated array (no heap) to uphold the WASM no-realloc
   * memory-view invariant; capacity 32 enforced at registration time.
   */
  struct ParamList {
    // Effect is the sole trusted mutator; the writable accessors below are
    // private so other callers see only the const overloads and route value
    // writes through updateParameter.
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
   * @return APPLIED if the value was written; otherwise the rejection reason
   *         (UNKNOWN_PARAM, READONLY, or NON_FINITE). The WASM bridge forwards
   *         this so the frontend can report why a write was dropped.
   * @details An accepted write to an animated parameter engages the effect's
   *          animation pause before storing the manual value.
   */
  ParamSetResult updateParameter(const char *name, float value) {
    auto *def = parameters.find(name);
    if (def == nullptr)
      return ParamSetResult::UNKNOWN_PARAM;
    // Untrusted JS boundary: reject readonly-param writes and non-finite input,
    // and clamp floats to [min,max]. Bools are thresholded at 0.5 by set().
    if (def->readonly)
      return ParamSetResult::READONLY;
    if (!std::isfinite(value))
      return ParamSetResult::NON_FINITE;
    // Enum targets hold an option index: snap a fractional write (e.g. a stale
    // deep link) to the nearest option before the range clamp.
    if (def->is_enum())
      value = roundf(value);
    if (!def->is_bool())
      value = hs::clamp(value, def->min, def->max);
    if (def->animated)
      setAnimationsPaused(true);
    def->set(value);
    return ParamSetResult::APPLIED;
  }

  /**
   * @brief Retrieves the list of registered parameters.
   * @return Const reference to the parameter list.
   */
  const ParamList &getParameters() const { return parameters; }

  /**
   * @brief Pause/resume the effect's parameter-driving animations.
   * @details Parameter animations and preset transitions wire this flag into
   * their timeline events. Paused, their active-time clocks and callbacks
   * freeze while ambient motion keeps running.
   * @param paused True to freeze parameter-driving animations, false to resume.
   */
  void setAnimationsPaused(bool paused) { anims_paused = paused; }
  /**
   * @brief Reports whether parameter-driving animations are paused.
   * @return True if those animations are currently frozen.
   */
  bool animations_paused() const { return anims_paused; }

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
  /** @brief Whether frame generation samples pixels outside the display band. */
  bool reads_outside_band;
  /**
   * @brief POV column-strobe flag (see strobe_columns()); set at construction.
   */
  bool strobe;
  ParamList parameters; /**< List of parameters. */
  /**
   * @brief Pause gate for parameter-driving animations.
   * @details Pass `&anims_paused` to Timeline::add_pausable, which freezes the
   * whole event: a not-yet-started event's delay counts active frames only.
   * Mutation/Driver/Lerp/Sprite also take an animation-level `paused` pointer,
   * which freezes stepping alone — a pending start delay keeps elapsing.
   */
  bool anims_paused = false;

  /**
   * @brief Flag a registered param as engine-written telemetry (read-only).
   * @details The GUI keeps showing its live value but disables editing. Use for
   * output-only values clobbered every frame (e.g. an active-particle count).
   */
  void mark_readonly(const char *name) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_readonly: unknown parameter name");
    def->readonly = true;
  }

  /** @brief Excludes a global parameter from preset exports. */
  void mark_global(const char *name) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_global: unknown parameter name");
    def->preset = false;
  }

  /**
   * @brief Registers a floating-point parameter.
   * @param name The name to expose.
   * @param ptr Pointer to the float variable.
   * @param min Minimum value.
   * @param max Maximum value.
   */
  HS_COLD_MEMBER void register_param(const char *name, float *ptr,
                                     float min = 0.0f, float max = 1.0f) {
    // Overflowing the fixed ParamList is an authoring bug (also upholds the WASM
    // no-realloc memory-view invariant).
    HS_CHECK(parameters.count < parameters.elements.size(),
             "register_param: exceeded ParamList capacity");
    // A duplicate name shadows: find() returns the FIRST match, so a second
    // registration's slot is unreachable by name.
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name");
    // An inverted range feeds hs::clamp() lo > hi (implementation-defined).
    HS_CHECK(min <= max, "register_param: min must be <= max");
    // A starting *ptr outside [min,max] would snap on the first GUI edit (every
    // updateParameter clamps).
    HS_CHECK(*ptr >= min && *ptr <= max,
             "register_param: default *ptr outside [min,max]");
    auto &def = parameters.elements[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.min = min;
    def.max = max;
  }

  /**
   * @brief Registers an enumerated parameter, rendered by the GUI as a dropdown.
   * @param name The name to expose.
   * @param ptr Pointer to the float variable holding the selected option index.
   * @param options Array of option labels indexed by the target's value; must
   *   outlive the effect (string literals).
   * @param option_count Number of labels; the value range is [0, option_count-1].
   */
  void register_param(const char *name, float *ptr, const char *const *options,
                      int option_count) {
    HS_CHECK(options != nullptr && option_count > 0,
             "register_param: enum needs at least one option");
    register_param(name, ptr, 0.0f, static_cast<float>(option_count - 1));
    parameters.elements[parameters.count - 1].options = options;
    parameters.elements[parameters.count - 1].option_count = option_count;
  }

  /**
   * @brief Registers a typed enum parameter, rendered by the GUI as a dropdown.
   * @tparam Enum Enum type stored by the target.
   * @param name The name to expose.
   * @param ptr Pointer to the enum variable.
   * @param options GUI labels indexed by the enum's underlying value.
   * @param export_options C++ enum literals indexed like @p options.
   * @param option_count Number of labels and literals.
   */
  template <typename Enum>
    requires std::is_enum_v<Enum>
  void register_param(const char *name, Enum *ptr, const char *const *options,
                      const char *const *export_options, int option_count) {
    HS_CHECK(options != nullptr && option_count > 0,
             "register_param: enum needs at least one option");
    HS_CHECK(parameters.count < parameters.elements.size(),
             "register_param: exceeded ParamList capacity");
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name");
    const float value =
        static_cast<float>(static_cast<std::underlying_type_t<Enum>>(*ptr));
    HS_CHECK(value >= 0.0f && value < static_cast<float>(option_count),
             "register_param: default enum outside option range");
    using Underlying = std::underlying_type_t<Enum>;
    static_assert(sizeof(Underlying) <= sizeof(uint32_t),
                  "register_param: enum underlying type exceeds 32 bits");
    constexpr auto TARGET_TYPE = [] {
      if constexpr (sizeof(Underlying) == sizeof(uint8_t))
        return std::is_signed_v<Underlying> ? ParamDef::TargetType::ENUM_I8
                                            : ParamDef::TargetType::ENUM_U8;
      if constexpr (sizeof(Underlying) == sizeof(uint16_t))
        return std::is_signed_v<Underlying> ? ParamDef::TargetType::ENUM_I16
                                            : ParamDef::TargetType::ENUM_U16;
      return std::is_signed_v<Underlying> ? ParamDef::TargetType::ENUM_I32
                                          : ParamDef::TargetType::ENUM_U32;
    }();
    auto &def = parameters.elements[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.min = 0.0f;
    def.max = static_cast<float>(option_count - 1);
    def.options = options;
    def.option_count = option_count;
    def.export_options = export_options;
    def.target_type = TARGET_TYPE;
  }

  /**
   * @brief Registers a boolean parameter.
   * @param name The name to expose.
   * @param ptr Pointer to the bool variable; registration never mutates the
   *   target, symmetric with the float overload.
   */
  HS_COLD_MEMBER void register_param(const char *name, bool *ptr) {
    HS_CHECK(parameters.count < parameters.elements.size(),
             "register_param: exceeded ParamList capacity");
    // Duplicate name guard, see the float overload.
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name");
    auto &def = parameters.elements[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.max = 1.0f;
    def.target_type = ParamDef::TargetType::BOOL;
  }

  /**
   * @brief Registers a float param and flags it animation-driven in one call.
   * @details Accepted writes engage the effect pause; the GUI also renders the
   * flagged param as an auto-pausing slider.
   */
  HS_COLD_MEMBER void register_animated_param(const char *name, float *ptr,
                                              float min = 0.0f,
                                              float max = 1.0f) {
    register_param(name, ptr, min, max);
    parameters.elements[parameters.count - 1].animated = true;
  }

  /**
   * @brief Registers a boolean param and flags it animation-driven.
   * @param name The name to expose.
   * @param ptr Pointer to the bool variable.
   */
  HS_COLD_MEMBER void register_animated_param(const char *name, bool *ptr) {
    register_param(name, ptr);
    parameters.elements[parameters.count - 1].animated = true;
  }

  /** @brief Registers a typed enum param and flags it animation-driven. */
  template <typename Enum>
    requires std::is_enum_v<Enum>
  HS_COLD_MEMBER void register_animated_param(const char *name, Enum *ptr,
                                              const char *const *options,
                                              const char *const *export_options,
                                              int option_count) {
    register_param(name, ptr, options, export_options, option_count);
    parameters.elements[parameters.count - 1].animated = true;
  }

  /**
   * @brief Registers a float param and flags it engine-written telemetry in one
   * call.
   */
  void register_readonly_param(const char *name, float *ptr, float min = 0.0f,
                               float max = 1.0f) {
    register_param(name, ptr, min, max);
    parameters.elements[parameters.count - 1].readonly = true;
  }

private:
  /**
   * @brief Advances the drawing buffer pointer to the next available buffer.
   * @details If `persist_pixels` is true, copies the previous frame's content
   * to the new buffer.
   */
  inline void advance_buffer() {
    int c = cur.load(std::memory_order_relaxed) ? 0 : 1;
    // The new write buffer must not be the one the ISR is scanning out (prev);
    // with two buffers this holds only if buffer_free() gated the advance. Trap
    // it here (once per frame, cold) instead of tearing.
    HS_CHECK(c != prev.load(std::memory_order_relaxed));
    cur.store(c, std::memory_order_relaxed);
    if (persist_pixels) {
      // The trail base is the last COMPLETED frame (next). The buffer_free()
      // gate forces prev == next, so copy from next and assert the equality
      // rather than depend on the gate silently across methods.
      int last = next.load(std::memory_order_relaxed);
      HS_CHECK(last == prev.load(std::memory_order_relaxed));
      memcpy(bufs[c], bufs[last], sizeof(Pixel) * frame_width * frame_height);
    }
  }

  /**
   * @brief Queues the newly drawn frame to be displayed.
   * @details Publishes `cur` as the new `next`. The release fence orders the
   * frame's pixel writes before the publish, pairing with the acquire fence in
   * `advance_display()`; the IRQ-off bracket keeps the publish atomic against
   * the on-device display ISR.
   */
  inline void queue_frame() {
    hs::disable_interrupts();
    std::atomic_thread_fence(std::memory_order_release);
    next.store(cur.load(std::memory_order_relaxed), std::memory_order_relaxed);
    hs::enable_interrupts();
  }

  /**
   * @brief Points bufs at the shared static storage and zeroes both buffers.
   * @details noinline so the two full-frame fills are emitted once, not inlined
   * into both GCC constructor variants (C1/C2). Invoked from the ctor (not init())
   * because derived init() overrides do not chain to Effect::init().
   */
  void __attribute__((noinline)) clear_buffers() {
    bufs[0] = buffer_a;
    std::fill_n(bufs[0], MAX_W * MAX_H, Pixel(0, 0, 0));
    bufs[1] = buffer_b;
    std::fill_n(bufs[1], MAX_W * MAX_H, Pixel(0, 0, 0));
  }

  inline void notify_buffer_ready() {
    if (buffer_ready_hook)
      buffer_ready_hook(*this);
  }

  std::atomic<int> prev{0}; /**< Buffer the ISR is currently reading. */
  std::atomic<int> cur{0};  /**< Buffer the main loop is currently writing. */
  std::atomic<int> next{0}; /**< Last completed frame, queued for display. */
  int frame_width;          /**< The width of the effect. */
  int frame_height;         /**< The height of the effect. */
  ClipRegion clip_region; /**< Segment clip region (display + render margin). */
  // Shared static storage for the double buffer. PRECONDITION: at most one Effect
  // live at a time (s_alive guard); a second would alias these arrays and the
  // prev/cur/next indices.
  static DMAMEM Pixel
      buffer_a[MAX_W * MAX_H]; /**< Static storage for buffer A (shared). */
  static DMAMEM Pixel
      buffer_b[MAX_W * MAX_H]; /**< Static storage for buffer B (shared). */
  Pixel *bufs[2]; /**< Pointers to the two buffer storage locations. */
  // True while an Effect is constructed-but-not-destroyed. Guards the
  // single-live-Effect precondition on the shared buffer_a/buffer_b: the ctor
  // traps if already set, the dtor clears it.
  static bool s_alive;
  BufferReadyHook buffer_ready_hook = nullptr;
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
   * @brief Constructs the Canvas, advancing the effect buffer and clearing
   * whatever of it can still hold stale pixels.
   * @param owner The effect instance owning the buffer.
   */
  Canvas(Effect &owner) : effect(owner) {
    wait_for_free_buffer();
    // Ordering is load-bearing: the hook runs post-wait but pre-clear, so the
    // clip it sets is the one clear_stale_pixels() honours.
    effect.notify_buffer_ready();
    effect.advance_buffer();
    if (!effect.persist_pixels) {
      HS_PROFILE(canvas_clear);
      clear_stale_pixels();
    }
  }

  /**
   * @brief Destructor. Queues the finished frame to be displayed.
   */
  ~Canvas() { effect.queue_frame(); }

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
    // assert, NOT HS_CHECK — the one hot-loop exception: a debug-only bounds
    // guard (stripped on device), since an always-on branch on every pixel access
    // is the single place HS_CHECK's contract forbids it.
    assert(x >= 0 && x < effect.frame_width && y >= 0 &&
           y < effect.frame_height);
    return effect.bufs[effect.cur.load(std::memory_order_relaxed)]
                      [y * effect.frame_width + x];
  }

  /**
   * @brief Accesses a pixel in the previous drawing buffer by 2D coordinates.
   * @param x The horizontal coordinate.
   * @param y The vertical coordinate.
   * @return Copy of the Pixel from the previous frame.
   * @details Returns by value, unlike `operator()`: the previous frame is
   *          read-only, so a reference into it would invite an accidental write
   *          to a buffer about to be recycled. A Pixel is small enough the copy
   *          is free.
   */
  inline Pixel prev(int x, int y) const {
    assert(x >= 0 && x < effect.frame_width && y >= 0 &&
           y < effect.frame_height);
    return effect.bufs[effect.prev.load(std::memory_order_relaxed)]
                      [y * effect.frame_width + x];
  }

  /**
   * @brief Raw base pointer to the previous-frame buffer (row-major,
   *        `width()` stride): `prev_data()[y * width() + x]` equals
   *        `prev(x, y)`.
   * @return Const base pointer, valid until this Canvas is destroyed.
   * @details Const because `prev()` returns by value specifically so a
   *          reference cannot be used to write into a buffer about to be
   *          recycled; the const pointer preserves that guarantee.
   */
  [[nodiscard]] inline const Pixel *prev_data() const {
    return effect.bufs[effect.prev.load(std::memory_order_relaxed)];
  }

  /**
   * @brief Raw base pointer to the current drawing buffer (row-major,
   *        `width()` stride): `data()[y * width() + x]` equals
   *        `(*this)(x, y)`.
   * @return Mutable base pointer, valid until this Canvas is destroyed.
   */
  [[nodiscard]] inline Pixel *data() {
    return effect.bufs[effect.cur.load(std::memory_order_relaxed)];
  }

  /**
   * @brief Accesses a pixel in the current drawing buffer by 1D index.
   * @param xy The 1D index.
   * @return Reference to the Pixel.
   */
  inline Pixel &operator()(int xy) {
    assert(xy >= 0 && xy < effect.frame_width * effect.frame_height);
    return effect.bufs[effect.cur.load(std::memory_order_relaxed)][xy];
  }

  /**
   * @brief Clears the entire current drawing buffer to black.
   */
  void clear_buffer() {
    int c = effect.cur.load(std::memory_order_relaxed);
    std::fill_n(effect.bufs[c], effect.frame_width * effect.frame_height,
                Pixel(0, 0, 0));
  }

  /**
   * @brief Gets the width of the underlying effect.
   * @return The width in pixels.
   */
  [[nodiscard]] inline int width() const { return effect.width(); }
  /**
   * @brief Gets the height of the underlying effect.
   * @return The height in pixels.
   */
  [[nodiscard]] inline int height() const { return effect.height(); }
  /**
   * @brief Gets the effect's current clip region.
   * @return Const reference to the clip region (display + render margin).
   */
  [[nodiscard]] inline const ClipRegion &clip() const {
    return effect.clip_region;
  }

  /**
   * @brief Checks if debug visuals are enabled.
   * @return True if debugging is active.
   */
  inline bool debug() const { return effect.debug_visuals; }

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
  /**
   * @brief Spins until the effect has a free back buffer.
   * @details Not a TOCTOU race: single-core, strict index ownership — the main
   * loop writes cur/next, the ISR only sets prev = next, so the gate can't be
   * falsified between check and flip. buffer_free() is re-satisfied only when the
   * display ISR advances at a frame boundary, so an unbounded wait means that ISR
   * stalled and the watchdog traps rather than hanging. The outer buffer_free()
   * guard skips the micros() reads on the no-wait path. The counter name is
   * load-bearing: Profile.ino derives each frame's render as wall minus the
   * counter whose name ends in _buffer_wait.
   */
  void wait_for_free_buffer() {
    HS_PROFILE(canvas_buffer_wait);
    if (effect.buffer_free())
      return;
    const unsigned long wait_start = micros();
    while (!effect.buffer_free()) {
      HS_CHECK(micros() - wait_start < BUFFER_FREE_WATCHDOG_US,
               "buffer_free watchdog timeout — display ISR stalled");
#ifdef HS_TEST_BUILD
      s_buffer_free_spins.fetch_add(1, std::memory_order_relaxed);
#endif
    }
  }

  /**
   * @brief Clears whatever of the freshly acquired buffer can still show stale
   *        pixels from the frame that last wrote it.
   * @details For an effect that does not read outside its display clip, the
   *          leftovers there are invisible and only the display band has to
   *          be cleared. It does draw into the margin-expanded render bounds,
   *          but that band is write-only scratch: nothing samples or displays
   *          it. A filter that samples across the band edge needs the whole
   *          buffer cleared, independently of whether its output crosses
   *          segment boundaries. Unsegmented targets clip to the full canvas,
   *          where the two are the same fill.
   */
  void clear_stale_pixels() {
#ifdef HS_TEST_BUILD
    if (effect.force_full_buffer_clear) {
      clear_buffer();
      return;
    }
#endif
    if (effect.reads_outside_band)
      clear_buffer();
    else
      clear_display_clip_buffer();
  }

  /**
   * @brief Clears the current display clip, excluding its render margin.
   * @details Device builds keep this helper out of line and outside ITCM.
   */
  HS_COLD_MEMBER void clear_display_clip_buffer() {
    const int c = effect.cur.load(std::memory_order_relaxed);
    const ClipRegion &clip = effect.clip_region;
    const int span = clip.x_end - clip.x_start;
    Pixel *const buffer = effect.bufs[c];

    if (span == effect.frame_width) {
      std::fill_n(buffer + clip.y_start * effect.frame_width,
                  span * (clip.y_end - clip.y_start), Pixel(0, 0, 0));
      return;
    }

    for (int y = clip.y_start; y < clip.y_end; ++y) {
      std::fill_n(buffer + y * effect.frame_width + clip.x_start, span,
                  Pixel(0, 0, 0));
    }
  }

  /** Watchdog bound for the ctor buffer_free() spin (µs). One display
   *  revolution is tens-to-hundreds of ms even at low RPM; 2 s is well above
   *  that, so only a genuinely stalled display ISR trips it. */
  static constexpr unsigned long BUFFER_FREE_WATCHDOG_US = 2000000UL;
  Effect &effect; /**< Reference to the owning Effect instance. */
#ifdef HS_TEST_BUILD
  inline static std::atomic<unsigned long> s_buffer_free_spins{0};
#endif
};
