/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <span>
#include <type_traits>
#include "platform/attributes.h"

/** @brief Registration policy for a target's initial value. */
enum class ParamInitialValue {
  REQUIRE_IN_RANGE,
  PRESERVE_REQUESTED_FLOAT,
};

template <typename T, bool = std::is_enum_v<T>> struct ParamInteger {
  using Type = T;
};
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
  using Bound = std::conditional_t<std::is_same_v<T, float>, float, int64_t>;

  Bound min = 0;
  Bound max = 1;
  bool animated = false;
  bool readonly = false;
  bool preset = true;
  const char *const *options = nullptr;
  const char *const *export_options = nullptr;
  int option_count = 0;
  std::span<const int64_t> option_values{}; /**< Empty selects dense indices. */
  ParamInitialValue initial_value = ParamInitialValue::REQUIRE_IN_RANGE;

  /** @brief Checks explicit IDs, their bounds, and the selected initial ID. */
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

  /** @brief Describes a dropdown over the contiguous indices [0,count-1]. */
  static constexpr ParamSpec enumerated(const char *const *labels, int count,
                                        const char *const *exports = nullptr) {
    return {.min = 0,
            .max = static_cast<Bound>(static_cast<int64_t>(count) - 1),
            .options = labels,
            .export_options = exports,
            .option_count = count};
  }
};
