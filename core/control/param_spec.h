/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file param_spec.h
 * @brief ParamSpec, the typed construction-time description of a
 *        runtime parameter.
 */

#include <cstdint>
#include <span>
#include <type_traits>
#include "platform/attributes.h"

/** @brief Registration policy for a target's initial value. */
enum class ParamInitialValue {
  REQUIRE_IN_RANGE,
  PRESERVE_REQUESTED_FLOAT,
};

/**
 * @brief Integer representation of a parameter target type.
 * @tparam T Integral or enum target type.
 */
template <typename T, bool = std::is_enum_v<T>> struct ParamInteger {
  using Type = T; ///< T itself for a non-enum.
};
/**
 * @brief ParamInteger for an enum: its underlying type.
 * @tparam T Enum target type.
 */
template <typename T> struct ParamInteger<T, true> {
  using Type = std::underlying_type_t<T>;
};

/**
 * @brief Typed construction-time parameter description.
 * @details Label arrays, value maps and their strings must outlive the host.
 * Integer bounds retain their exact values until registration validates them.
 */
template <typename T> struct ParamSpec {
  static_assert(
      std::is_same_v<T, float> ||
          ((std::is_integral_v<T> || std::is_enum_v<T>) &&
           sizeof(T) <= sizeof(uint32_t)),
      "parameter target must be float, bool, or an integer/enum up to 32 bits");
  /// Bound type: float for a float target, else int64_t.
  using Bound = std::conditional_t<std::is_same_v<T, float>, float, int64_t>;

  Bound min = 0;         ///< Lower bound, inclusive.
  Bound max = 1;         ///< Upper bound, inclusive.
  bool animated = false; ///< True if an animation drives the target.
  bool readonly = false; ///< True for engine-written telemetry.
  bool preset = true;    ///< Whether preset exports include the parameter.
  /// Option labels (GUI dropdown), or null for a plain parameter.
  const char *const *options = nullptr;
  /// C++ enum literals indexed like `options`, or null.
  const char *const *export_options = nullptr;
  int option_count = 0; ///< Number of labels; > 0 marks an enum target.
  std::span<const int64_t> option_values{}; /**< Empty selects dense indices. */
  /// Registration policy for the target's initial value.
  ParamInitialValue initial_value = ParamInitialValue::REQUIRE_IN_RANGE;

  /**
   * @brief Checks explicit IDs, their bounds, and the selected initial ID.
   * @param initial Target's initial value.
   * @return True when `option_values` is empty or consistent with the labels,
   *         includes `min` and `max`, and includes `initial`.
   */
  HS_COLD_MEMBER constexpr bool valid_option_values(T initial) const {
    if (option_values.empty())
      return true;
    if (options == nullptr || option_count <= 0 ||
        option_values.size() != static_cast<size_t>(option_count))
      return false;
    bool includes_min = false;
    bool includes_max = false;
    bool includes_initial = false;
    for (size_t i = 0; i < option_values.size(); ++i) {
      const int64_t value = option_values[i];
      if (value < INT32_MIN || value > UINT32_MAX ||
          static_cast<int64_t>(static_cast<float>(value)) != value ||
          value < min || value > max)
        return false;
      includes_min |= value == min;
      includes_max |= value == max;
      includes_initial |= static_cast<Bound>(initial) == value;
      for (size_t j = 0; j < i; ++j)
        if (value == option_values[j])
          return false;
    }
    return includes_min && includes_max && includes_initial;
  }

  /**
   * @brief Describes a dropdown over the contiguous indices [0,count-1].
   * @param labels Dropdown labels.
   * @param count Label count.
   * @param exports C++ enumerator names indexed like `labels`, or null.
   * @return The dropdown description.
   */
  static constexpr ParamSpec enumerated(const char *const *labels, int count,
                                        const char *const *exports = nullptr) {
    return {.min = 0,
            .max = static_cast<Bound>(static_cast<int64_t>(count) - 1),
            .options = labels,
            .export_options = exports,
            .option_count = count};
  }
};
