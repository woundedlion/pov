/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file param_host.h
 * @brief ParamHost: the runtime parameter registry an effect exposes, the
 *        validated write gate the WASM bridge calls, and the animation pause
 *        an accepted animated write engages.
 */

#include "control/params.h"
#include "platform/platform.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

/**
 * @brief Owns an effect's registered parameters and the pause gate their
 *        writes engage.
 * @details Naming split: methods on the JS/embind boundary (updateParameter,
 * getParameters, setAnimationsPaused) are camelCase to match the WASM bridge;
 * the internal C++ API (register_param, reset_parameters) is snake_case.
 */
class ParamHost {
public:
  /** @brief Runtime parameter descriptor (see control/params.h). */
  using ParamDef = ::ParamDef;
  /** @brief Fixed-capacity parameter registry (see control/params.h). */
  using ParamList = ::ParamList;

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
    const bool animated = def->animated;
    if (animated)
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
    if (animated)
      animated_parameter_written();
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
  /**
   * @brief Runs after an accepted write to an animated parameter.
   * @details The write has already engaged the animation pause. An effect whose
   * choreography rewrites the whole parameter set — a preset crossfade — must
   * stop doing so here, or the next frame overwrites the value just written.
   */
  virtual void animated_parameter_written() {}

#if HS_ENABLE_PARAM_GUI_BRIDGE
  using ParameterUpdatedHook = void (*)(ParamHost *, const char *, bool);

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
    // set() narrows through static_cast<Integer>(float), which is UB outside
    // the target's range.
    const bool range_fits =
        static_cast<int64_t>(min) >=
            static_cast<int64_t>(std::numeric_limits<Integer>::min()) &&
        static_cast<int64_t>(max) <=
            static_cast<int64_t>(std::numeric_limits<Integer>::max());
    HS_CHECK(range_fits,
             "register_int_param: [min,max] must fit the target integer type");
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
   * @brief Registers a uint8_t-backed enumerated param (GUI dropdown) and
   *        flags it animation-driven.
   * @param name The name to expose.
   * @param ptr Pointer to the uint8_t variable holding the selected index.
   * @param options GUI labels indexed by the target's value; must outlive the
   *   effect.
   * @param option_count Number of labels; the value range is
   *   [0, option_count - 1].
   * @details For enum-shaped storage that is not a C++ enum type — the chain
   * interpreter's topology enum8s, whose option lists exist only as runtime
   * table data. No export literals: preset export names no enum type.
   */
  HS_COLD_MEMBER void register_animated_enum8_param(const char *name,
                                                    uint8_t *ptr,
                                                    const char *const *options,
                                                    int option_count) {
    HS_CHECK(options != nullptr && option_count > 0,
             "register_param: enum needs at least one option");
    HS_CHECK(parameters.count < parameters.capacity(),
             "register_param: exceeded ParamList capacity");
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name");
    HS_CHECK(*ptr < option_count,
             "register_param: default enum outside option range");
    auto &def = parameters.data()[parameters.count++];
    def = {};
    def.name = name;
    def.target = ptr;
    def.min = 0.0f;
    def.max = static_cast<float>(option_count - 1);
    def.options = options;
    def.option_count = option_count;
    def.target_type = ParamDef::TargetType::INT_U8;
    def.animated = true;
    parameters.bump_schema_generation();
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
};
