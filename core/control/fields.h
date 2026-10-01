/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

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

/** @brief Typed registration, validation and interpolation of one member. */
template <typename Owner, typename Value = float> struct Field {
  const char *id;
  Value Owner::*member;
  const char *name;
  ParamSpec<Value> spec;
  FieldCurve curve = std::is_same_v<Value, float> ? FieldCurve::RAW_LINEAR
                                                  : FieldCurve::MIDPOINT;
  bool interpolated = true;
  bool validated = true;
  typename ParamSpec<Value>::Bound validation_min = spec.min;
  typename ParamSpec<Value>::Bound validation_max = spec.max;

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
      return static_cast<int64_t>(sample) >= validation_min &&
             static_cast<int64_t>(sample) <= validation_max;
    }
  }

  void interpolate(Owner &out, const Owner &from, const Owner &to,
                   float progress) const {
    if (!interpolated)
      return;
    if constexpr (std::is_same_v<Value, float>)
      out.*member = apply_curve(curve, from.*member, to.*member, progress);
    else
      out.*member = progress < (curve == FieldCurve::SNAP ? 1.0f : .5f)
                        ? from.*member
                        : to.*member;
  }

  template <typename Register>
  void register_to(Owner &owner, Register &add) const {
    if (name)
      add(name, &(owner.*member), spec);
  }
};

template <typename Owner, typename Value>
Field(const char *, Value Owner::*, const char *, ParamSpec<Value>)
    -> Field<Owner, Value>;

/** @brief Descriptions of a nested parameter aggregate. */
template <typename Owner, typename Value, typename Fields> struct FieldGroup {
  Value Owner::*member;
  Fields fields;

  constexpr bool valid(const Owner &owner) const {
    return std::apply(
        [&](const auto &...field) {
          return (field.valid(owner.*member) && ...);
        },
        fields);
  }
  void interpolate(Owner &out, const Owner &from, const Owner &to,
                   float progress) const {
    std::apply(
        [&](const auto &...field) {
          (field.interpolate(out.*member, from.*member, to.*member, progress),
           ...);
        },
        fields);
  }
  template <typename Register>
  void register_to(Owner &owner, Register &add) const {
    std::apply(
        [&](const auto &...field) {
          (field.register_to(owner.*member, add), ...);
        },
        fields);
  }
};

template <typename Owner, typename Value, typename Fields>
FieldGroup(Value Owner::*, Fields) -> FieldGroup<Owner, Value, Fields>;

template <typename Owner, typename Fields>
constexpr bool valid_fields(const Owner &owner, const Fields &fields) {
  return std::apply(
      [&](const auto &...field) { return (field.valid(owner) && ...); },
      fields);
}

/** @brief Writes described fields, preserving excluded and untabled state. */
template <typename Owner, typename Fields>
void interpolate_fields(Owner &out, const Owner &from, const Owner &to,
                        float progress, const Fields &fields) {
  std::apply(
      [&](const auto &...field) {
        (field.interpolate(out, from, to, progress), ...);
      },
      fields);
}

template <typename Owner, typename Fields, typename Register>
void register_fields(Owner &owner, const Fields &fields, Register add) {
  std::apply(
      [&](const auto &...field) { (field.register_to(owner, add), ...); },
      fields);
}

} // namespace Control
