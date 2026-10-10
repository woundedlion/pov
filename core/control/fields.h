/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file fields.h
 * @brief Typed parameter descriptions (Field, FieldGroup): GUI
 *        registration, validation and preset interpolation.
 */

#include <tuple>
#include <type_traits>
#include <limits>

#include "control/param_spec.h"
#include "math/interpolate.h"

namespace Control {

/** @brief Interpolation domain of a described parameter. */
enum class FieldCurve : uint8_t {
  LERP,
  LOG_POSITIVE,
  SHORTEST_PERIODIC,
  SHORTEST_TURN,
  SNAP,
  MIDPOINT,
  RAW_LINEAR
};

/**
 * @brief Blends two values along an interpolation curve.
 * @param curve Interpolation domain.
 * @param from Value at t = 0.
 * @param to Value at t = 1.
 * @param t Progress in [0, 1].
 * @return The blended value.
 */
HS_FLASH_INLINE inline float apply_curve(FieldCurve curve, float from, float to,
                                         float t) {
  switch (curve) {
  case FieldCurve::LERP:
    return interp::linear(from, to, t);
  case FieldCurve::LOG_POSITIVE:
    return interp::log_positive(from, to, t);
  case FieldCurve::SHORTEST_PERIODIC:
    return interp::shortest_periodic(from, to, t, math::TWO_PI_F);
  case FieldCurve::SHORTEST_TURN:
    return interp::shortest_periodic(from, to, t, 1.0f);
  case FieldCurve::MIDPOINT:
    return t < .5f ? from : to;
  case FieldCurve::RAW_LINEAR:
    return from + (to - from) * t;
  case FieldCurve::SNAP:
    return t < 1.0f ? from : to;
  }
  return from;
}

/** @brief Describes one parameter member: its GUI registration, valid range
    and interpolation curve. */
template <typename Owner, typename Value = float> struct Field {
  const char *id;        ///< Stable identifier.
  Value Owner::*member;  ///< Described member of Owner.
  const char *name;      ///< GUI label; null skips registration.
  ParamSpec<Value> spec; ///< Registration range and options.
  /// Interpolation curve; non-float fields support MIDPOINT and SNAP only.
  FieldCurve curve = std::is_same_v<Value, float> ? FieldCurve::RAW_LINEAR
                                                  : FieldCurve::MIDPOINT;
  bool interpolated = true; ///< False leaves the member out of interpolate().
  bool validated = true;    ///< False makes valid() always true.
  /// Lowest value valid() accepts.
  typename ParamSpec<Value>::Bound validation_min = spec.min;
  /// Highest value valid() accepts.
  typename ParamSpec<Value>::Bound validation_max = spec.max;

  /** @return Whether `curve` is supported for Value. */
  constexpr bool curve_supported() const {
    if constexpr (std::is_same_v<Value, float>)
      return true;
    else
      return curve == FieldCurve::MIDPOINT || curve == FieldCurve::SNAP;
  }

  /**
   * @brief Checks the member against the validation range and option IDs.
   * @param owner Aggregate to check.
   * @return True when unvalidated or in range; false on a malformed range.
   */
  constexpr bool valid(const Owner &owner) const {
    if (!validated)
      return true;
    const auto sample = owner.*member;
    if constexpr (std::is_same_v<Value, float>) {
      constexpr float LIMIT = std::numeric_limits<float>::max();
      if (!(validation_min >= -LIMIT && validation_max <= LIMIT &&
            validation_min <= validation_max))
        return false;
      return sample >= validation_min && sample <= validation_max;
    } else {
      using Integer = typename ParamInteger<Value>::Type;
      if (validation_min < std::numeric_limits<Integer>::lowest() ||
          validation_max > std::numeric_limits<Integer>::max() ||
          validation_min > validation_max)
        return false;
      const auto ID = static_cast<int64_t>(sample);
      if (ID < validation_min || ID > validation_max)
        return false;
      if (!spec.option_values.empty()) {
        for (const auto option : spec.option_values)
          if (ID == option)
            return true;
        return false;
      }
      return true;
    }
  }

  /**
   * @brief Writes the member blended along `curve`; no-op when not
   *        interpolated.
   * @param out Aggregate receiving the blended member.
   * @param from Start aggregate.
   * @param to End aggregate.
   * @param progress Blend progress in [0, 1].
   */
  HS_FLASH_INLINE void interpolate(Owner &out, const Owner &from,
                                   const Owner &to, float progress) const {
    if (!interpolated)
      return;
    if constexpr (std::is_same_v<Value, float>)
      out.*member = apply_curve(curve, from.*member, to.*member, progress);
    else
      out.*member = progress < (curve == FieldCurve::SNAP ? 1.0f : .5f)
                        ? from.*member
                        : to.*member;
  }

  /**
   * @brief Registers the member as a parameter when `name` is set.
   * @tparam Register Callable `(const char *, Value *, const ParamSpec &)`.
   * @param owner Aggregate holding the registered member.
   * @param add Registration callback.
   */
  template <typename Register>
  HS_FLASH_INLINE void register_to(Owner &owner, Register &add) const {
    if (name)
      add(name, &(owner.*member), spec);
  }
};

/** @brief Deduces Owner and Value from the member pointer. */
template <typename Owner, typename Value>
Field(const char *, Value Owner::*, const char *, ParamSpec<Value>)
    -> Field<Owner, Value>;

/** @brief Describes the fields of a nested struct member. */
template <typename Owner, typename Value, typename Fields> struct FieldGroup {
  Value Owner::*member; ///< Nested struct member of Owner.
  Fields fields;        ///< Tuple of Field descriptions over Value.

  /** @return Whether every nested field's curve is supported. */
  constexpr bool curve_supported() const {
    return std::apply(
        [](const auto &...field) { return (field.curve_supported() && ...); },
        fields);
  }

  /**
   * @brief Validates every nested field.
   * @param owner Aggregate holding the nested struct.
   * @return True when every nested field is valid.
   */
  constexpr bool valid(const Owner &owner) const {
    return std::apply(
        [&](const auto &...field) {
          return (field.valid(owner.*member) && ...);
        },
        fields);
  }
  /**
   * @brief Interpolates every nested field.
   * @param out Aggregate receiving the blended members.
   * @param from Start aggregate.
   * @param to End aggregate.
   * @param progress Blend progress in [0, 1].
   */
  HS_FLASH_INLINE void interpolate(Owner &out, const Owner &from,
                                   const Owner &to, float progress) const {
    std::apply(
        [&](const auto &...field) {
          (field.interpolate(out.*member, from.*member, to.*member, progress),
           ...);
        },
        fields);
  }
  /**
   * @brief Registers every nested field.
   * @tparam Register Registration callable, as for Field::register_to().
   * @param owner Aggregate holding the nested struct.
   * @param add Registration callback.
   */
  template <typename Register>
  HS_FLASH_INLINE void register_to(Owner &owner, Register &add) const {
    std::apply(
        [&](const auto &...field) {
          (field.register_to(owner.*member, add), ...);
        },
        fields);
  }
};

/** @brief Deduces Owner, Value and Fields from the constructor arguments. */
template <typename Owner, typename Value, typename Fields>
FieldGroup(Value Owner::*, Fields) -> FieldGroup<Owner, Value, Fields>;

/**
 * @brief Whether every non-float field uses MIDPOINT or SNAP.
 * @param fields Tuple of Field or FieldGroup descriptions.
 * @return True when every field's curve is supported.
 */
template <typename Fields>
constexpr bool curves_supported(const Fields &fields) {
  return std::apply(
      [](const auto &...field) { return (field.curve_supported() && ...); },
      fields);
}

/**
 * @brief Validates every described field of an aggregate.
 * @tparam Owner Parameter aggregate.
 * @tparam Fields Tuple of Field or FieldGroup descriptions.
 * @param owner Aggregate to check.
 * @param fields Descriptions to check against.
 * @return True when every field is valid.
 */
template <typename Owner, typename Fields>
constexpr bool valid_fields(const Owner &owner, const Fields &fields) {
  return std::apply(
      [&](const auto &...field) { return (field.valid(owner) && ...); },
      fields);
}

/** @brief Interpolates each described field of @p out between @p from and
    @p to; undescribed members and fields with interpolated=false are left
    unchanged.
 * @param out Aggregate receiving the blended members.
 * @param from Start aggregate.
 * @param to End aggregate.
 * @param progress Blend progress in [0, 1].
 * @param fields Tuple of Field or FieldGroup descriptions.
 */
template <typename Owner, typename Fields>
HS_FLASH_INLINE void interpolate_fields(Owner &out, const Owner &from,
                                        const Owner &to, float progress,
                                        const Fields &fields) {
  std::apply(
      [&](const auto &...field) {
        (field.interpolate(out, from, to, progress), ...);
      },
      fields);
}

/**
 * @brief Registers every named described field as a parameter.
 * @tparam Owner Parameter aggregate.
 * @tparam Fields Tuple of Field or FieldGroup descriptions.
 * @tparam Register Registration callable, as for Field::register_to().
 * @param owner Aggregate holding the registered members.
 * @param fields Descriptions to register.
 * @param add Registration callback.
 */
template <typename Owner, typename Fields, typename Register>
HS_FLASH_INLINE void register_fields(Owner &owner, const Fields &fields,
                                     Register add) {
  std::apply(
      [&](const auto &...field) { (field.register_to(owner, add), ...); },
      fields);
}

} // namespace Control
