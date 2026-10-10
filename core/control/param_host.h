/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file param_host.h
 * @brief ParamHost: an effect's named, GUI-editable parameters, validated
 *        writes to them, and the animation pause flag.
 */

#include "control/params.h"
#include "control/param_spec.h"
#include "memory.h"
#include "platform/platform.h"
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

/**
 * @brief Holds an effect's registered parameters and validates writes to
 *        them by name.
 * @details JS/embind boundary methods are camelCase; the internal C++ API is
 * snake_case.
 */
class ParamHost {
public:
  /** @brief Runtime parameter descriptor. */
  using ParamDef = ::ParamDef;
  /** @brief Fixed-capacity parameter registry. */
  using ParamList = ::ParamList;

  /**
   * @brief Updates a parameter's value by name.
   * @param name The name of the parameter.
   * @param value The new value (mapped to bool if necessary).
   * @return APPLIED if the value was written; otherwise the rejection reason
   *         (UNKNOWN_PARAM, READONLY, NON_FINITE, or INADMISSIBLE).
   * @details Writing a parameter registered as animated first pauses the
   *          effect's animations.
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
  /** @brief Reapplies recorded writes in order without per-write checks, then
      fails fatally if the final values are not admissible.
   * @param values Ordered (name, value) pairs.
   */
  template <typename Values>
  void replay_parameter_writes(const Values &values) {
    check_parameter_storage();
    for (const auto &[name, value] : values) {
      auto *def = parameters.find(name.c_str());
      HS_CHECK(def != nullptr && !def->readonly,
               "replay_parameter_writes: unknown or readonly parameter");
      apply_parameter(*def, value);
    }
    for (const auto &def : parameters)
      if (!def.readonly)
        HS_CHECK(parameter_write_admitted(def, def.get_requested()),
                 "replay_parameter_writes: inadmissible final state");
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
  /** @brief Updates the displayed values of mirrored parameters. */
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

  /** @brief Counter bumped whenever the parameter list or a descriptor's
      flags change; always 0 without the GUI bridge.
   * @return The schema generation.
   */
  uint32_t getParameterSchemaGeneration() const {
    return parameters.schema_generation();
  }

  /**
   * @brief Pause/resume the effect's parameter-driving animations.
   * @details Animations bound to this flag stop advancing while it is set;
   * unbound animations and unpausable preset blends keep running.
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

  /**
   * @brief Stores a trusted internal value without edit policy or callbacks.
   * @param parameter Descriptor whose target is written.
   * @param value Value to store.
   */
  static void write_parameter_unchecked(ParamDef &parameter, float value) {
    parameter.write_unchecked(value);
  }

  /**
   * @brief Called after an accepted write to a parameter that presets
   *        include (not one marked global).
   */
  virtual void parameter_written() {}

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Returns false to reject a write; called before anything is
      stored.
   * @return False to reject the write.
   */
  virtual bool parameter_write_admitted(const ParamDef &, float) {
    return true;
  }

  /// Post-write callback: (host, parameter name, parameter is an enum).
  using ParameterUpdatedHook = void (*)(ParamHost *, const char *, bool);

  /**
   * @brief Sets a callback run after each accepted parameter write.
   * @param hook Callback; null clears it.
   */
  void set_parameter_updated_hook(ParameterUpdatedHook hook) {
    parameter_updated_hook = hook;
  }
#endif

  /**
   * @brief Repoints the registry at external descriptor storage.
   * @param storage Descriptor array the registry registers into.
   * @param capacity Slots @p storage holds.
   * @details Bumps the schema generation.
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

  /** @brief Repoints the registry at arena-owned descriptor storage; debug
      builds check the arena block is still alive on each access.
   * @param arena Arena that owns `storage`.
   * @param storage Descriptor array the registry registers into.
   * @param capacity Slots in `storage`.
   */
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

  /** @brief Clears every registration and bumps the schema generation. */
  void reset_parameters() {
    parameters.count = 0;
    parameters.bump_schema_generation();
  }

  /**
   * @brief Displays parameter values from @p displayed while writes still go
   *        to @p requested.
   * @details Each parameter whose target lies inside @p requested displays the
   * member at the same offset in @p displayed; others display their target.
   * @param requested State the parameter targets point into.
   * @param displayed State whose members are displayed.
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
  /**
   * @brief Shows a parameter's requested value instead of its render mirror.
   * @param name Parameter name; unknown names are ignored.
   */
  void show_requested_parameter_value(const char *name) {
    ParamDef *parameter = parameters.find(name);
    if (parameter != nullptr)
      parameter->display_target = nullptr;
  }
#endif
  ParamList parameters; /**< Registered parameters. */
#if HS_ENABLE_PARAM_GUI_BRIDGE
  ParameterUpdatedHook parameter_updated_hook = nullptr; ///< Null for none.
#endif
  /** @brief True while parameter-driving animations are paused. */
  bool anims_paused = false;

  /**
   * @brief Flag a registered param as engine-written telemetry (read-only).
   * @details The GUI shows its value but disables editing.
   * @param name Registered parameter name.
   * @param readonly New read-only flag.
   */
  void mark_readonly(const char *name, bool readonly = true) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_readonly: unknown parameter name name=%s", name);
    if (def->readonly != readonly) {
      def->readonly = readonly;
      parameters.bump_schema_generation();
    }
  }

  /**
   * @brief Excludes a global parameter from preset exports.
   * @param name Registered parameter name.
   */
  void mark_global(const char *name) {
    auto *def = parameters.find(name);
    HS_CHECK(def, "mark_global: unknown parameter name name=%s", name);
    def->preset = false;
    parameters.bump_schema_generation();
  }

  /**
   * @brief Registers a parameter described by @p spec; the target's current
   *        value is its initial value.
   * @details PRESERVE_REQUESTED_FLOAT allows a finite out-of-range initial
   * value for a non-enum float; later writes are clamped to [min, max].
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param spec Range, options and flags.
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
    if (options != nullptr && spec.option_values.empty())
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
          if (spec.option_values.empty()) {
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
        }
        HS_CHECK(min <= max,
                 "register_param: min must be <= max name=%s min=%lld max=%lld",
                 name, static_cast<long long>(min),
                 static_cast<long long>(max));
        const bool range_fits =
            min >= static_cast<int64_t>(std::numeric_limits<Integer>::min()) &&
            max <= static_cast<int64_t>(std::numeric_limits<Integer>::max());
        HS_CHECK(
            range_fits,
            "register_param: [min,max] must fit the target integer type name=%s min=%lld max=%lld",
            name, static_cast<long long>(min), static_cast<long long>(max));
        const bool bounds_exact =
            static_cast<int64_t>(static_cast<float>(min)) == min &&
            static_cast<int64_t>(static_cast<float>(max)) == max;
        HS_CHECK(
            bounds_exact,
            "register_param: bounds must be exactly representable as float name=%s min=%lld max=%lld",
            name, static_cast<long long>(min), static_cast<long long>(max));
        const int64_t value = static_cast<int64_t>(*ptr);
        HS_CHECK(
            value >= min && value <= max,
            "register_param: default *ptr outside [min,max] name=%s value=%lld min=%lld max=%lld",
            name, static_cast<long long>(value), static_cast<long long>(min),
            static_cast<long long>(max));
        target_type = integer_target_type<Integer>();
      }
    }
    HS_CHECK(spec.valid_option_values(*ptr),
             "register_param: invalid explicit option values name=%s", name);
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
#if HS_ENABLE_PARAM_GUI_BRIDGE
    def.export_options = spec.export_options;
#endif
    def.option_values =
        spec.option_values.empty() ? nullptr : spec.option_values.data();
    def.option_count = option_count;
    def.animated = spec.animated;
    def.readonly = spec.readonly;
    def.preset = spec.preset;
    parameters.bump_schema_generation();
  }

  void register_param(const char *, float *, int, int) = delete;

  /**
   * @brief Registers a float slider, optionally with dropdown labels.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member; its current value is the initial value.
   * @param min Lower bound.
   * @param max Upper bound.
   * @param animated True if an animation drives the target.
   * @param readonly True for engine-written telemetry.
   * @param options Dropdown labels, or null.
   * @param option_count Label count.
   */
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
  /** @brief Registers an animated float whose initial value may lie outside
      [min, max].
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member; its current value is the initial value.
   * @param min Lower bound for later writes.
   * @param max Upper bound for later writes.
   */
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

  /**
   * @brief Registers a float-backed dropdown over indices 0..option_count-1.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member holding the selected index.
   * @param options Dropdown labels.
   * @param option_count Label count.
   */
  HS_COLD_MEMBER void register_param(const char *name, float *ptr,
                                     const char *const *options,
                                     int option_count) {
    register_param(name, ptr,
                   ParamSpec<float>::enumerated(options, option_count));
  }

  /** @brief Registers an enum dropdown; @p export_options optionally names
      each option's C++ enumerator for preset export.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param options Dropdown labels, indexed by enumerator value.
   * @param export_options C++ enumerator names indexed like `options`, or null.
   * @param option_count Label count.
   * @param animated True if an animation drives the target.
   */
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

  /**
   * @brief Registers an integer slider with exactly float-representable bounds.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param min Lower bound.
   * @param max Upper bound.
   * @param animated True if an animation drives the target.
   */
  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool>)
  HS_COLD_MEMBER void register_int_param(const char *name, Integer *ptr,
                                         int min, int max,
                                         bool animated = false) {
    register_param(
        name, ptr,
        ParamSpec<Integer>{.min = min, .max = max, .animated = animated});
  }

  /**
   * @brief Registers an animation-driven integer slider.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param min Lower bound.
   * @param max Upper bound.
   */
  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool>)
  HS_COLD_MEMBER void register_animated_int_param(const char *name,
                                                  Integer *ptr, int min,
                                                  int max) {
    register_param(
        name, ptr,
        ParamSpec<Integer>{.min = min, .max = max, .animated = true});
  }

  /**
   * @brief Registers a bool without changing its initial value.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param animated True if an animation drives the target.
   */
  HS_COLD_MEMBER void register_param(const char *name, bool *ptr,
                                     bool animated = false) {
    register_param(name, ptr, ParamSpec<bool>{.animated = animated});
  }

  /**
   * @brief Registers an animation-driven float slider.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param min Lower bound.
   * @param max Upper bound.
   */
  HS_COLD_MEMBER void register_animated_param(const char *name, float *ptr,
                                              float min = 0.0f,
                                              float max = 1.0f) {
    register_param(name, ptr,
                   ParamSpec<float>{.min = min, .max = max, .animated = true});
  }

  /**
   * @brief Registers an animation-driven bool.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   */
  HS_COLD_MEMBER void register_animated_param(const char *name, bool *ptr) {
    register_param(name, ptr, ParamSpec<bool>{.animated = true});
  }

  /**
   * @brief Registers an animation-driven typed dropdown.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param options Dropdown labels, indexed by enumerator value.
   * @param export_options C++ enumerator names indexed like `options`, or null.
   * @param option_count Label count.
   */
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

  /**
   * @brief Registers an animation-driven uint8_t-backed dropdown.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member holding the selected index.
   * @param options Dropdown labels.
   * @param option_count Label count.
   */
  HS_COLD_MEMBER void register_animated_enum8_param(const char *name,
                                                    uint8_t *ptr,
                                                    const char *const *options,
                                                    int option_count) {
    auto spec = ParamSpec<uint8_t>::enumerated(options, option_count);
    spec.animated = true;
    register_param(name, ptr, spec);
  }

  /**
   * @brief Registers a readonly float slider.
   * @param name Parameter name; must outlive the host.
   * @param ptr Target member.
   * @param min Lower bound.
   * @param max Upper bound.
   */
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
