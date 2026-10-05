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
#include "control/param_spec.h"
#include "engine/memory.h"
#include "platform/platform.h"
#include <array>
#include <bit>
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
   *         (UNKNOWN_PARAM, READONLY, NON_FINITE, or INADMISSIBLE). The WASM bridge forwards
   *         this so the frontend can report why a write was dropped.
   * @details An accepted write to an animated parameter engages the effect's
   *          animation pause before storing the manual value.
   */
  ParamSetResult updateParameter(const char *name, float value) {
    check_parameter_storage();
    auto *def = parameters.find(name);
    if (def == nullptr)
      return ParamSetResult::UNKNOWN_PARAM;
    const ParamSetResult result = def->normalize(value);
    if (result != ParamSetResult::APPLIED)
      return result;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if (!parameter_write_admitted(*def, value))
      return ParamSetResult::INADMISSIBLE;
#endif
    apply_parameter(*def, value);
    return ParamSetResult::APPLIED;
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Restores trusted captured writes in order, then validates the state. */
  template <typename Values> void restore_parameters(const Values &values) {
    check_parameter_storage();
    for (const auto &[name, value] : values) {
      auto *def = parameters.find(name.c_str());
      HS_CHECK(def != nullptr && !def->readonly,
               "restore_parameters: unknown or readonly parameter");
      apply_parameter(*def, value);
    }
    for (const auto &def : parameters)
      if (!def.readonly)
        HS_CHECK(parameter_write_admitted(def, def.get_requested()),
                 "restore_parameters: inadmissible final state");
  }
#endif

  /**
   * @brief Retrieves the list of registered parameters.
   * @return Const reference to the parameter list.
   */
  const ParamList &getParameters() const {
    check_parameter_storage();
    return parameters;
  }

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
   * @return Borrowed warning text, or null when the parameter is valid.
   * Copy the text before the next warning query or mutation of this host.
   */
  virtual const char *parameter_warning(const char *name) const {
    (void)name;
    return nullptr;
  }
#endif

  /** @brief Ordered parameter schema-change token; 0 without the GUI bridge. */
  uint32_t getParameterSchemaGeneration() const {
    return parameters.schema_generation();
  }

  /**
   * @brief Pause/resume the effect's parameter-driving animations.
   * @details Events wired to this flag freeze their active-time clocks and
   * callbacks while paused. Ambient motion and unpausable preset blends keep
   * running.
   * @param paused True to freeze parameter-driving animations, false to resume.
   */
  void setAnimationsPaused(bool paused) { anims_paused = paused; }
  /**
   * @brief Reports whether parameter-driving animations are paused.
   * @return True if those animations are currently frozen.
   */
  bool animations_paused() const { return anims_paused; }

protected:
  ~ParamHost() = default;

  /** @brief Stores a trusted internal value without edit policy or callbacks. */
  static void write_parameter_unchecked(ParamDef &parameter, float value) {
    parameter.write_unchecked(value);
  }

  /**
   * @brief Runs after any accepted write to a preset (non-global) parameter.
   * @details Preset crossfades must stop rewriting the manually edited state.
   */
  virtual void parameter_written() {}

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Validates a candidate value before any parameter state changes. */
  virtual bool parameter_write_admitted(const ParamDef &, float) {
    return true;
  }

  using ParameterUpdatedHook = void (*)(ParamHost *, const char *, bool);

  /** @brief Installs an opt-in reaction to accepted GUI parameter writes. */
  void set_parameter_updated_hook(ParameterUpdatedHook hook) {
    parameter_updated_hook = hook;
  }
#endif

  /**
   * @brief Repoints the registry at external descriptor storage.
   * @param storage Descriptor array the registry registers into.
   * @param capacity Slots @p storage holds.
   * @details Relocates the array getParameters() hands out, so it bumps the
   * schema generation a cached descriptor view is keyed on.
   */
  void use_parameter_storage(ParamDef *storage, size_t capacity) {
    HS_CHECK(parameters.count == 0,
             "use_parameter_storage: parameters already registered");
    HS_CHECK(storage != nullptr && capacity > 0,
             "use_parameter_storage: invalid external storage");
#if HS_PARAM_EXTERNAL_STORAGE
#ifndef NDEBUG
    parameter_storage_stamp.clear();
#endif
    parameters.external_elements = storage;
    parameters.external_capacity = capacity;
    parameters.bump_schema_generation();
#else
    (void)storage;
    (void)capacity;
    HS_CHECK(false, "use_parameter_storage: external storage is disabled");
#endif
  }

  /** @brief Binds arena-owned descriptor storage and records its lifetime. */
  void use_parameter_storage(Arena &arena, ParamDef *storage, size_t capacity) {
    use_parameter_storage(storage, capacity);
#ifndef NDEBUG
    parameter_storage_stamp.record(arena);
#else
    (void)arena;
#endif
  }

  /**
   * @brief Repoints the registry at a fixed-size descriptor array.
   * @tparam CAPACITY Slots the array holds.
   * @param storage Descriptor array the registry registers into.
   */
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
   * @details Targets outside @p requested read their requested value directly.
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
  void mark_readonly(const char *name, bool readonly = true) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_readonly: unknown parameter name name=%s", name);
    if (def->readonly != readonly) {
      def->readonly = readonly;
      parameters.bump_schema_generation();
    }
  }

  /** @brief Excludes a global parameter from preset exports. */
  void mark_global(const char *name) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_global: unknown parameter name name=%s", name);
    def->preset = false;
    parameters.bump_schema_generation();
  }

  /**
   * @brief Registers one typed description without changing its target.
   * @details Registration appends one descriptor and advances its schema token
   * once. PRESERVE_REQUESTED_FLOAT admits finite out-of-range GUI requests only
   * for ordinary float targets; subsequent edits still obey the published bounds.
   */
  template <typename T>
  HS_COLD_MEMBER void register_param(const char *name, T *ptr,
                                     const ParamSpec<T> &spec) {
    HS_CHECK(name != nullptr, "register_param: null parameter name");
    HS_CHECK(ptr != nullptr, "register_param: null target name=%s", name);
    const auto min = spec.min;
    const auto max = spec.max;
    const char *const *options = spec.options;
    const int option_count = spec.option_count;
    HS_CHECK((options == nullptr) == (option_count == 0),
             "register_param: inconsistent options and count");
    HS_CHECK(spec.export_options == nullptr || options != nullptr,
             "register_param: export labels require options name=%s", name);
    HS_CHECK(spec.initial_value == ParamInitialValue::REQUIRE_IN_RANGE ||
                 spec.initial_value ==
                     ParamInitialValue::PRESERVE_REQUESTED_FLOAT,
             "register_param: invalid initial-value policy name=%s", name);
    if (options != nullptr)
      HS_CHECK(option_count > 0 && min == 0 && max == option_count - 1,
               "register_param: option range does not match labels");

    ParamDef::TargetType target_type;
    if constexpr (std::is_same_v<T, float>) {
      HS_CHECK(
          min <= max,
          "register_param: min must be <= max name=%s min_bits=%08lx max_bits=%08lx",
          name, static_cast<unsigned long>(std::bit_cast<uint32_t>(min)),
          static_cast<unsigned long>(std::bit_cast<uint32_t>(max)));
      HS_CHECK(std::isfinite(min) && std::isfinite(max),
               "register_param: bounds must be finite name=%s", name);
      const bool preserve =
          spec.initial_value == ParamInitialValue::PRESERVE_REQUESTED_FLOAT;
#if HS_ENABLE_PARAM_GUI_BRIDGE
      HS_CHECK(
          !preserve || options == nullptr,
          "register_param: preserving requested values requires a non-enum float name=%s",
          name);
#else
      HS_CHECK(
          !preserve,
          "register_param: preserving requested values requires the GUI bridge name=%s",
          name);
#endif
      if (preserve)
        HS_CHECK(std::isfinite(*ptr),
                 "register_param: requested value must be finite name=%s",
                 name);
      else
        HS_CHECK(
            *ptr >= min && *ptr <= max,
            "register_param: default *ptr outside [min,max] name=%s value_bits=%08lx min_bits=%08lx max_bits=%08lx",
            name, static_cast<unsigned long>(std::bit_cast<uint32_t>(*ptr)),
            static_cast<unsigned long>(std::bit_cast<uint32_t>(min)),
            static_cast<unsigned long>(std::bit_cast<uint32_t>(max)));
      target_type = ParamDef::TargetType::FLOAT;
    } else {
      HS_CHECK(
          spec.initial_value == ParamInitialValue::REQUIRE_IN_RANGE,
          "register_param: preserving requested values requires a non-enum float name=%s",
          name);
      if constexpr (std::is_same_v<T, bool>) {
        HS_CHECK(
            options == nullptr && min == 0 && max == 1,
            "register_param: bool range must be [0,1] without options name=%s",
            name);
        target_type = ParamDef::TargetType::BOOL;
      } else {
        using Integer = typename ParamInteger<T>::Type;
        if constexpr (std::is_enum_v<T>) {
          HS_CHECK(
              options != nullptr && option_count > 0,
              "register_param: enum needs at least one option name=%s count=%d",
              name, option_count);
          HS_CHECK(
              static_cast<int64_t>(option_count - 1) <=
                  static_cast<int64_t>(std::numeric_limits<Integer>::max()),
              "register_param: options must fit the target enum type name=%s count=%d",
              name, option_count);
          HS_CHECK(
              static_cast<int64_t>(static_cast<float>(option_count - 1)) ==
                  option_count - 1,
              "register_param: enum bound must be exactly representable as float name=%s count=%d",
              name, option_count);
        }
        HS_CHECK(
            min <= max,
            "register_int_param: min must be <= max name=%s min=%lld max=%lld",
            name, static_cast<long long>(min), static_cast<long long>(max));
        const bool range_fits =
            min >= static_cast<int64_t>(std::numeric_limits<Integer>::min()) &&
            max <= static_cast<int64_t>(std::numeric_limits<Integer>::max());
        HS_CHECK(
            range_fits,
            "register_int_param: [min,max] must fit the target integer type name=%s min=%lld max=%lld",
            name, static_cast<long long>(min), static_cast<long long>(max));
        const bool bounds_exact =
            static_cast<int64_t>(static_cast<float>(min)) == min &&
            static_cast<int64_t>(static_cast<float>(max)) == max;
        HS_CHECK(
            bounds_exact,
            "register_int_param: bounds must be exactly representable as float name=%s min=%lld max=%lld",
            name, static_cast<long long>(min), static_cast<long long>(max));
        const int64_t value = static_cast<int64_t>(*ptr);
        HS_CHECK(
            value >= min && value <= max,
            "register_int_param: default *ptr outside [min,max] name=%s value=%lld min=%lld max=%lld",
            name, static_cast<long long>(value), static_cast<long long>(min),
            static_cast<long long>(max));
        target_type = integer_target_type<Integer>();
      }
    }
    for (int i = 0; i < option_count; ++i) {
      HS_CHECK(options[i] != nullptr,
               "register_param: null option label name=%s index=%d", name, i);
      HS_CHECK(spec.export_options == nullptr ||
                   spec.export_options[i] != nullptr,
               "register_param: null export label name=%s index=%d", name, i);
    }
    auto &def = append_parameter(name);
    def.target = ptr;
    def.min = static_cast<float>(min);
    def.max = static_cast<float>(max);
    def.target_type = target_type;
    def.options = options;
    def.export_options = spec.export_options;
    def.option_count = option_count;
    def.animated = spec.animated;
    def.readonly = spec.readonly;
    def.preset = spec.preset;
    parameters.bump_schema_generation();
  }

  void register_param(const char *, float *, int, int) = delete;

  /** @brief Registers a float slider, optionally with dropdown labels. */
  HS_COLD_MEMBER void
  register_param(const char *name, float *ptr, float min = 0.0f,
                 float max = 1.0f, bool animated = false, bool readonly = false,
                 const char *const *options = nullptr, int option_count = 0) {
    register_param(name, ptr,
                   ParamSpec<float>{.min = min,
                                    .max = max,
                                    .animated = animated,
                                    .readonly = readonly,
                                    .options = options,
                                    .option_count = option_count});
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Registers an animated float with a finite requested value. */
  HS_COLD_MEMBER void register_animated_param_preserving_value(const char *name,
                                                               float *ptr,
                                                               float min,
                                                               float max) {
    register_param(
        name, ptr,
        ParamSpec<float>{.min = min,
                         .max = max,
                         .animated = true,
                         .initial_value =
                             ParamInitialValue::PRESERVE_REQUESTED_FLOAT});
  }
#endif

  /** @brief Registers a float-backed dropdown with numeric preset exports. */
  HS_COLD_MEMBER void register_param(const char *name, float *ptr,
                                     const char *const *options,
                                     int option_count) {
    register_param(name, ptr,
                   ParamSpec<float>::enumerated(options, option_count));
  }

  /** @brief Registers a typed dropdown with optional C++ export literals. */
  template <typename Enum>
    requires std::is_enum_v<Enum>
  HS_COLD_MEMBER void register_param(const char *name, Enum *ptr,
                                     const char *const *options,
                                     const char *const *export_options,
                                     int option_count, bool animated = false) {
    auto spec =
        ParamSpec<Enum>::enumerated(options, option_count, export_options);
    spec.animated = animated;
    register_param(name, ptr, spec);
  }

  /** @brief Registers an integer slider with exactly float-representable bounds. */
  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool>)
  HS_COLD_MEMBER void register_int_param(const char *name, Integer *ptr,
                                         int min, int max,
                                         bool animated = false) {
    register_param(
        name, ptr,
        ParamSpec<Integer>{.min = min, .max = max, .animated = animated});
  }

  /** @brief Registers an animation-driven integer slider. */
  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool>)
  HS_COLD_MEMBER void register_animated_int_param(const char *name,
                                                  Integer *ptr, int min,
                                                  int max) {
    register_param(
        name, ptr,
        ParamSpec<Integer>{.min = min, .max = max, .animated = true});
  }

  /** @brief Registers a bool without changing its initial value. */
  HS_COLD_MEMBER void register_param(const char *name, bool *ptr,
                                     bool animated = false) {
    register_param(name, ptr, ParamSpec<bool>{.animated = animated});
  }

  /** @brief Registers an animation-driven float slider. */
  HS_COLD_MEMBER void register_animated_param(const char *name, float *ptr,
                                              float min = 0.0f,
                                              float max = 1.0f) {
    register_param(name, ptr,
                   ParamSpec<float>{.min = min, .max = max, .animated = true});
  }

  /** @brief Registers an animation-driven bool. */
  HS_COLD_MEMBER void register_animated_param(const char *name, bool *ptr) {
    register_param(name, ptr, ParamSpec<bool>{.animated = true});
  }

  /** @brief Registers an animation-driven typed dropdown. */
  template <typename Enum>
    requires std::is_enum_v<Enum>
  HS_COLD_MEMBER void register_animated_param(const char *name, Enum *ptr,
                                              const char *const *options,
                                              const char *const *export_options,
                                              int option_count) {
    auto spec =
        ParamSpec<Enum>::enumerated(options, option_count, export_options);
    spec.animated = true;
    register_param(name, ptr, spec);
  }

  /** @brief Registers an animation-driven uint8_t-backed dropdown. */
  HS_COLD_MEMBER void register_animated_enum8_param(const char *name,
                                                    uint8_t *ptr,
                                                    const char *const *options,
                                                    int option_count) {
    auto spec = ParamSpec<uint8_t>::enumerated(options, option_count);
    spec.animated = true;
    register_param(name, ptr, spec);
  }

  /** @brief Registers a readonly float slider. */
  HS_COLD_MEMBER void register_readonly_param(const char *name, float *ptr,
                                              float min = 0.0f,
                                              float max = 1.0f) {
    register_param(name, ptr,
                   ParamSpec<float>{.min = min, .max = max, .readonly = true});
  }

private:
  void apply_parameter(ParamDef &parameter, float value) {
    auto *def = &parameter;
    const bool animated = def->animated;
    const bool preset = def->preset;
    if (animated)
      setAnimationsPaused(true);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const char *updated_name = def->name;
    const bool updated_enum = def->is_enum();
    def->write_unchecked(value);
    if (parameter_updated_hook != nullptr)
      parameter_updated_hook(this, updated_name, updated_enum);
#else
    def->write_unchecked(value);
#endif
    if (preset)
      parameter_written();
  }

  HS_COLD_MEMBER ParamDef &append_parameter(const char *name) {
    HS_CHECK(parameters.count < parameters.capacity(),
             "register_param: exceeded ParamList capacity name=%s", name);
    HS_CHECK(parameters.find(name) == nullptr,
             "register_param: duplicate parameter name name=%s", name);
    auto &def = parameters.data()[parameters.count++];
    def = {};
    def.name = name;
    return def;
  }

#ifndef NDEBUG
  ArenaBlockStamp parameter_storage_stamp;
#endif
  void check_parameter_storage() const {
#if HS_PARAM_EXTERNAL_STORAGE
    HS_ASSERT_BLOCK_ALIVE(parameter_storage_stamp, parameters.external_elements,
                          parameters.external_capacity * sizeof(ParamDef),
                          "ParamHost");
#endif
  }

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
