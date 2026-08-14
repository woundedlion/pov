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
#include <cstdint>
#include <type_traits>
#include <utility>
#include "engine/constants.h"
#include "engine/effect_params.h"
#include "render/clip.h"
#include "color/color.h"
#include <array>
#if HS_ENABLE_TEST_HOOKS
#include <cstdlib>
#endif

/**
 * @file canvas.h
 * @brief EffectConfig, the Effect base class, and the Canvas pixel buffer.
 */

class Canvas;

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
 * @brief Folds a filter pipeline's compile-time segment traits into an
 *        EffectConfig.
 * @tparam PipelineT Filter pipeline or sink exposing `any_crosses_segments`,
 *         `any_reads_outside_band`, and `total_segment_margin`.
 * @param base Config carrying the effect's own flags, including any full_frame,
 *        reads_outside_band or margin the effect needs for its own reasons.
 * @return @p base with the pipeline's full_frame, reads_outside_band and margin
 *         requirements combined in; the other fields are untouched.
 * @details The single definition of the fold: an effect that stacks a filter
 *          crossing segment boundaries gets the full-frame render without
 *          restating the three traits at its base initializer. All three are
 *          "at least this much" requirements, so the fold widens and never
 *          clears what the caller asked for.
 */
template <typename PipelineT>
constexpr EffectConfig pipeline_config(EffectConfig base = {}) {
  base.full_frame = base.full_frame || PipelineT::any_crosses_segments;
  base.reads_outside_band =
      base.reads_outside_band || PipelineT::any_reads_outside_band;
  base.margin = std::max(base.margin, PipelineT::total_segment_margin);
  return base;
}

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
#if HS_ENABLE_TEST_HOOKS
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

  /** @brief Number of presets exposed for manual navigation. */
  size_t getPresetCount() const { return preset_count; }
  /** @brief Index of the currently selected preset. */
  size_t getPresetIndex() const { return preset_index; }
  /** @brief Selects and pauses one preset, or reports that it is unavailable. */
  bool selectPreset(size_t index) {
    if (!change_preset(index, PresetChangeOrigin::MANUAL))
      return false;
    setAnimationsPaused(true);
    return true;
  }
  /** @brief Selects one preset without changing the animation pause state. */
  bool synchronizePreset(size_t index) {
    return preset_count > 0 &&
           (index == preset_index ||
            change_preset(index, PresetChangeOrigin::SYNCHRONIZED));
  }
  /** @brief Selects and pauses the next preset. */
  bool nextPreset() {
    return preset_count > 0 && selectPreset((preset_index + 1) % preset_count);
  }
  /** @brief Selects and pauses the previous preset. */
  bool previousPreset() {
    return preset_count > 0 &&
           selectPreset((preset_index + preset_count - 1) % preset_count);
  }

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
   * @details The acquire load pairs with the release half of
   * `advance_display()`'s fence: everything the display ISR stored before the
   * flip — the segmented driver's display-window half — is visible to the
   * writer this gate releases.
   */
  [[nodiscard]] inline bool buffer_free() const {
    return prev.load(std::memory_order_acquire) ==
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
   * @details The fence's acquire half pairs with `queue_frame()`'s release
   * fence, so the queued frame's pixel writes are visible before any read
   * through `prev`; its release half orders whatever the caller stored before
   * the flip — the segmented driver's display-window publish — ahead of the
   * `prev` store, pairing with `buffer_free()`'s acquire load. An acq_rel fence
   * carries both halves in the one barrier an acquire fence already cost. The
   * per-pixel loads themselves stay relaxed.
   */
  inline void advance_display() {
    int n = next.load(std::memory_order_relaxed);
    std::atomic_thread_fence(std::memory_order_acq_rel);
    prev.store(n, std::memory_order_relaxed);
  }
  /** @brief Runtime parameter descriptor (see engine/effect_params.h). */
  using ParamDef = ::ParamDef;
  /** @brief Fixed-capacity parameter registry (see engine/effect_params.h). */
  using ParamList = ::ParamList;

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
    // Enum and integer targets hold whole numbers: snap a fractional write
    // (e.g. a stale deep link) to the nearest before the range clamp. An enum
    // may be float-backed, so neither predicate implies the other.
    if (def->is_enum() || def->is_integer())
      value = roundf(value);
    if (!def->is_bool())
      value = hs::clamp(value, def->min, def->max);
    if (def->animated)
      setAnimationsPaused(true);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const char *updated_name = def->name;
    const bool updated_enum = def->is_enum();
    def->set(value);
    if (parameter_updated_hook != nullptr)
      parameter_updated_hook(this, updated_name, updated_enum);
#else
    def->set(value);
#endif
    return ParamSetResult::APPLIED;
  }

  /**
   * @brief Retrieves the list of registered parameters.
   * @return Const reference to the parameter list.
   */
  const ParamList &getParameters() const { return parameters; }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Refreshes values exposed through parameter display mirrors. */
  virtual void refresh_parameter_display() {}

  /**
   * @brief Reads the value accepted for rendering for one parameter.
   * @param parameter Registered parameter descriptor.
   * @return Accepted value in the parameter's native numeric representation.
   */
  virtual float accepted_parameter_value(const ParamDef &parameter) const {
    return parameter.get_requested();
  }

  /**
   * @brief Reports an actionable GUI warning for one parameter.
   * @param name Registered parameter name.
   * @return Warning text, or null when the parameter is valid.
   */
  virtual const char *parameter_warning(const char *name) const {
    (void)name;
    return nullptr;
  }
#endif

  /** @brief Token identifying the current ordered parameter descriptor schema. */
  uint32_t getParameterSchemaGeneration() const {
    return parameters.schema_generation();
  }

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
  /** @brief Whether a preset move came from choreography or a user control. */
  enum class PresetChangeOrigin : uint8_t { AUTOMATIC, MANUAL, SYNCHRONIZED };

  /** @brief One validated preset transition handed to the effect hook. */
  struct PresetChange {
    size_t from;               /**< Previously committed preset index. */
    size_t to;                 /**< Candidate preset index. */
    PresetChangeOrigin origin; /**< Source of the preset change. */
  };

  /** @brief Enables the shared preset controller for this effect. */
  HS_FLASH_MEMBER void configure_presets(size_t count) {
    HS_CHECK(count > 0, "preset count must be positive");
    HS_CHECK(preset_count == 0, "presets already configured");
    preset_count = count;
  }

  /** @brief Advances choreography through the shared preset controller. */
  HS_FLASH_MEMBER bool advancePreset() {
    return preset_count > 0 && change_preset((preset_index + 1) % preset_count,
                                             PresetChangeOrigin::AUTOMATIC);
  }

  /** @brief Applies a candidate preset before its index is committed. */
  virtual bool apply_preset(const PresetChange &) { return false; }
  /** @brief Runs after a successful preset change has been committed. */
  virtual void preset_changed(const PresetChange &) {}

  bool change_preset(size_t index, PresetChangeOrigin origin) {
    if (index >= preset_count)
      return false;
    const PresetChange change{preset_index, index, origin};
    if (!apply_preset(change))
      return false;
    preset_index = index;
    preset_changed(change);
    return true;
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  using ParameterUpdatedHook = void (*)(Effect *, const char *, bool);

  /** @brief Installs an opt-in reaction to accepted GUI parameter writes. */
  void set_parameter_updated_hook(ParameterUpdatedHook hook) {
    parameter_updated_hook = hook;
  }
#endif

  void use_parameter_storage(ParamDef *storage, size_t capacity) {
    HS_CHECK(parameters.count == 0,
             "use_parameter_storage: parameters already registered");
    HS_CHECK(storage != nullptr && capacity > 0,
             "use_parameter_storage: invalid external storage");
#if HS_PARAM_EXTERNAL_STORAGE
    parameters.external_elements = storage;
    parameters.external_capacity = capacity;
#else
    (void)storage;
    (void)capacity;
    HS_CHECK(false, "use_parameter_storage: external storage is disabled");
#endif
  }

  template <size_t CAPACITY>
  void use_parameter_storage(std::array<ParamDef, CAPACITY> &storage) {
    use_parameter_storage(storage.data(), storage.size());
  }

  void reset_parameters() {
    parameters.count = 0;
    parameters.bump_schema_generation();
  }

  /**
   * @brief Reads registered values from a live mirror while writes target the
   *        corresponding requested state.
   * @details Targets outside @p requested retain their own value source.
   */
  template <typename State>
  void mirror_parameter_display_state(const State &requested,
                                      const State &displayed) {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const std::uintptr_t requested_begin =
        reinterpret_cast<std::uintptr_t>(&requested);
    const std::uintptr_t requested_end = requested_begin + sizeof(State);
    const std::uintptr_t displayed_begin =
        reinterpret_cast<std::uintptr_t>(&displayed);
    for (ParamDef &parameter : parameters) {
      const std::uintptr_t target =
          reinterpret_cast<std::uintptr_t>(parameter.target);
      if (target < requested_begin || target >= requested_end) {
        parameter.display_target = nullptr;
        continue;
      }
      parameter.display_target = reinterpret_cast<const void *>(
          displayed_begin + (target - requested_begin));
    }
#else
    (void)requested;
    (void)displayed;
#endif
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Shows a parameter's requested value instead of its render mirror. */
  void show_requested_parameter_value(const char *name) {
    ParamDef *parameter = parameters.find(name);
    if (parameter != nullptr)
      parameter->display_target = nullptr;
  }
#endif

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
#if HS_ENABLE_PARAM_GUI_BRIDGE
  ParameterUpdatedHook parameter_updated_hook = nullptr;
#endif
  /**
   * @brief Pause gate for parameter-driving animations.
   * @details Pass `&anims_paused` to Timeline::add_pausable, which freezes the
   * whole event: a not-yet-started event's delay counts active frames only.
   * Mutation/Driver/Lerp/Sprite also take an animation-level `paused` pointer,
   * which freezes stepping alone — a pending start delay keeps elapsing.
   */
  bool anims_paused = false;
  size_t preset_count = 0;
  size_t preset_index = 0;

  /**
   * @brief Flag a registered param as engine-written telemetry (read-only).
   * @details The GUI keeps showing its live value but disables editing. Use for
   * output-only values clobbered every frame (e.g. an active-particle count).
   */
  void mark_readonly(const char *name) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_readonly: unknown parameter name");
    def->readonly = true;
    parameters.bump_schema_generation();
  }

  /** @brief Excludes a global parameter from preset exports. */
  void mark_global(const char *name) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_global: unknown parameter name");
    def->preset = false;
    parameters.bump_schema_generation();
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
    HS_CHECK(parameters.count < parameters.capacity(),
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
    auto &def = parameters.data()[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.min = min;
    def.max = max;
    parameters.bump_schema_generation();
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /**
   * @brief Registers an animated float while preserving its requested value.
   * @details The requested value may lie outside the published range.
   */
  HS_COLD_MEMBER void register_animated_param_preserving_value(const char *name,
                                                               float *ptr,
                                                               float min,
                                                               float max) {
    HS_CHECK(parameters.count < parameters.capacity(),
             "register_param: exceeded ParamList capacity");
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name");
    HS_CHECK(min <= max, "register_param: min must be <= max");
    auto &def = parameters.data()[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.min = min;
    def.max = max;
    def.animated = true;
    parameters.bump_schema_generation();
  }
#endif

  /**
   * @brief Registers an enumerated parameter, rendered by the GUI as a dropdown.
   * @param name The name to expose.
   * @param ptr Pointer to the float variable holding the selected option index.
   * @param options Array of option labels indexed by the target's value; must
   *   outlive the effect (string literals).
   * @param option_count Number of labels; the value range is [0, option_count-1].
   * @details Preset exports write the selected index as a numeric literal — a
   * float target names no enum type. Use the Enum* overload to export C++
   * enumerators instead.
   */
  HS_COLD_MEMBER void register_param(const char *name, float *ptr,
                                     const char *const *options,
                                     int option_count) {
    HS_CHECK(options != nullptr && option_count > 0,
             "register_param: enum needs at least one option");
    register_param(name, ptr, 0.0f, static_cast<float>(option_count - 1));
    parameters.data()[parameters.count - 1].options = options;
    parameters.data()[parameters.count - 1].option_count = option_count;
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
  HS_COLD_MEMBER void
  register_param(const char *name, Enum *ptr, const char *const *options,
                 const char *const *export_options, int option_count) {
    HS_CHECK(options != nullptr && option_count > 0,
             "register_param: enum needs at least one option");
    HS_CHECK(parameters.count < parameters.capacity(),
             "register_param: exceeded ParamList capacity");
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name");
    const float value =
        static_cast<float>(static_cast<std::underlying_type_t<Enum>>(*ptr));
    HS_CHECK(value >= 0.0f && value < static_cast<float>(option_count),
             "register_param: default enum outside option range");
    constexpr auto TARGET_TYPE =
        integer_target_type<std::underlying_type_t<Enum>>();
    auto &def = parameters.data()[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.min = 0.0f;
    def.max = static_cast<float>(option_count - 1);
    def.options = options;
    def.option_count = option_count;
    def.export_options = export_options;
    def.target_type = TARGET_TYPE;
    parameters.bump_schema_generation();
  }

  /**
   * @brief Registers an integer parameter, rendered by the GUI as a unit-step
   *   slider.
   * @tparam Integer Integral type stored by the target, at most 32 bits.
   * @param name The name to expose.
   * @param ptr Pointer to the integer variable.
   * @param min Minimum value, inclusive.
   * @param max Maximum value, inclusive.
   * @details For a quantity whose target is a count rather than a choice: the
   * range carries the bound, so no label array is needed and preset exports
   * write a plain integer literal. A target with one distinct meaning per value
   * belongs on the Enum* overload, which exports enumerators instead.
   */
  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool>)
  HS_COLD_MEMBER void register_int_param(const char *name, Integer *ptr,
                                         int min, int max) {
    HS_CHECK(parameters.count < parameters.capacity(),
             "register_int_param: exceeded ParamList capacity");
    HS_CHECK(parameters.find(name) == nullptr,
             "register_int_param: duplicate parameter name");
    HS_CHECK(min <= max, "register_int_param: min must be <= max");
    const int value = static_cast<int>(*ptr);
    HS_CHECK(value >= min && value <= max,
             "register_int_param: default *ptr outside [min,max]");
    auto &def = parameters.data()[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.min = static_cast<float>(min);
    def.max = static_cast<float>(max);
    def.target_type = integer_target_type<Integer>();
    parameters.bump_schema_generation();
  }

  /** @brief Registers an integer param and flags it animation-driven. */
  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool>)
  HS_COLD_MEMBER void register_animated_int_param(const char *name,
                                                  Integer *ptr, int min,
                                                  int max) {
    register_int_param(name, ptr, min, max);
    parameters.data()[parameters.count - 1].animated = true;
  }

  /**
   * @brief Registers a boolean parameter.
   * @param name The name to expose.
   * @param ptr Pointer to the bool variable; registration never mutates the
   *   target, symmetric with the float overload.
   */
  HS_COLD_MEMBER void register_param(const char *name, bool *ptr) {
    HS_CHECK(parameters.count < parameters.capacity(),
             "register_param: exceeded ParamList capacity");
    // Duplicate name guard, see the float overload.
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name");
    auto &def = parameters.data()[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.max = 1.0f;
    def.target_type = ParamDef::TargetType::BOOL;
    parameters.bump_schema_generation();
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
    parameters.data()[parameters.count - 1].animated = true;
  }

  /**
   * @brief Registers a boolean param and flags it animation-driven.
   * @param name The name to expose.
   * @param ptr Pointer to the bool variable.
   */
  HS_COLD_MEMBER void register_animated_param(const char *name, bool *ptr) {
    register_param(name, ptr);
    parameters.data()[parameters.count - 1].animated = true;
  }

  /** @brief Registers a typed enum param and flags it animation-driven. */
  template <typename Enum>
    requires std::is_enum_v<Enum>
  HS_COLD_MEMBER void register_animated_param(const char *name, Enum *ptr,
                                              const char *const *options,
                                              const char *const *export_options,
                                              int option_count) {
    register_param(name, ptr, options, export_options, option_count);
    parameters.data()[parameters.count - 1].animated = true;
  }

  /**
   * @brief Registers a float param and flags it engine-written telemetry in one
   * call.
   */
  HS_COLD_MEMBER void register_readonly_param(const char *name, float *ptr,
                                              float min = 0.0f,
                                              float max = 1.0f) {
    register_param(name, ptr, min, max);
    parameters.data()[parameters.count - 1].readonly = true;
  }

private:
  /** @brief TargetType matching an integral storage type's width and sign. */
  template <typename Integer>
  static constexpr ParamDef::TargetType integer_target_type() {
    static_assert(sizeof(Integer) <= sizeof(uint32_t),
                  "parameter integer type exceeds 32 bits");
    if constexpr (sizeof(Integer) == sizeof(uint8_t))
      return std::is_signed_v<Integer> ? ParamDef::TargetType::INT_I8
                                       : ParamDef::TargetType::INT_U8;
    else if constexpr (sizeof(Integer) == sizeof(uint16_t))
      return std::is_signed_v<Integer> ? ParamDef::TargetType::INT_I16
                                       : ParamDef::TargetType::INT_U16;
    else
      return std::is_signed_v<Integer> ? ParamDef::TargetType::INT_I32
                                       : ParamDef::TargetType::INT_U32;
  }

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
  HS_FLASH_MEMBER void __attribute__((noinline)) clear_buffers() {
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

#if HS_ENABLE_TEST_HOOKS
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
#if HS_ENABLE_TEST_HOOKS
    const unsigned long WATCHDOG_US = buffer_free_watchdog_us();
#else
    constexpr unsigned long WATCHDOG_US = BUFFER_FREE_WATCHDOG_US;
#endif
    const unsigned long wait_start = micros();
    while (!effect.buffer_free()) {
      HS_CHECK(micros() - wait_start < WATCHDOG_US,
               "buffer_free watchdog timeout — display ISR stalled");
#if HS_ENABLE_TEST_HOOKS
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
#if HS_ENABLE_TEST_HOOKS
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

#if HS_ENABLE_TEST_HOOKS
  /**
   * @brief Test-only watchdog bound, overridable via HS_BUFFER_FREE_WATCHDOG_US.
   * @return The bound in µs; BUFFER_FREE_WATCHDOG_US unless the environment
   *         supplies a positive override.
   * @details The watchdog is a trap, so tripping it kills the whole test
   *          process rather than failing one test; a sanitizer job whose
   *          threads deschedule past the shipping bound raises it here. Read
   *          once per process, outside the spin loop.
   */
  static unsigned long buffer_free_watchdog_us() {
    static const unsigned long RESOLVED = [] {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
      const char *e = std::getenv("HS_BUFFER_FREE_WATCHDOG_US");
#pragma clang diagnostic pop
      if (e) {
        const unsigned long v = std::strtoul(e, nullptr, 10);
        if (v > 0)
          return v;
      }
      return BUFFER_FREE_WATCHDOG_US;
    }();
    return RESOLVED;
  }
#endif

  Effect &effect; /**< Reference to the owning Effect instance. */
#if HS_ENABLE_TEST_HOOKS
  inline static std::atomic<unsigned long> s_buffer_free_spins{0};
#endif
};
