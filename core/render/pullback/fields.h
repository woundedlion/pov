/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <string_view>

#include "math/interpolate.h"
#include "render/pullback/contract.h"

/**
 * @file fields.h
 * @brief Per-family field descriptors: one table per parameter family drives
 *        slider registration, preset interpolation and snapshot validation.
 */

namespace Pullback {

/** @brief How a field moves across a preset transition. */
enum class FieldCurve : uint8_t {
  LERP,              /**< Linear, exact at both endpoints. */
  LOG_POSITIVE,      /**< Geometric, for positive scales. */
  SHORTEST_PERIODIC, /**< Shortest arc on a 2*pi period. */
  SNAP               /**< Holds the start value until progress reaches 1. */
};

/** @brief When a field's slider is registered. */
enum class FieldGate : uint8_t {
  ALWAYS,
  ANIMATED_PROJECTION, /**< Only when the effect animates its projection. */
  CENTRAL_MERIDIAN /**< Only when the effect declares USES_CENTRAL_MERIDIAN. */
};

/**
 * @brief One scalar field of a parameter family.
 * @details `name == nullptr` marks a field that is interpolated and validated
 * but exposes no slider of its own: a warp slot's `speed`, which the runtime
 * registers under the slot's name, or a field only an envelope reads.
 * @tparam Owner The family struct the field belongs to.
 */
template <typename Owner> struct Field {
  const char *id; /**< Stable machine id; kebab-case, unique in the family. */
  float Owner::*member;
  const char *name;
  float min;
  float max;
  FieldCurve curve = FieldCurve::LERP;
  FieldGate gate = FieldGate::ALWAYS;
};

/** @brief Whether @p T carries a field-descriptor table. */
template <typename T>
concept HasFields = requires { T::FIELDS; };

/** @brief Whether every tabled field id is non-null and unique in the family. */
template <HasFields Family> consteval bool field_ids_unique() {
  for (size_t i = 0; i < Family::FIELDS.size(); ++i) {
    if (Family::FIELDS[i].id == nullptr)
      return false;
    for (size_t j = 0; j < i; ++j)
      if (std::string_view(Family::FIELDS[i].id) == Family::FIELDS[j].id)
        return false;
  }
  return true;
}

namespace Fields {

HS_FLASH_INLINE inline float apply_curve(FieldCurve curve, float from, float to,
                                         float t) {
  switch (curve) {
  case FieldCurve::LERP:
    return interp::linear(from, to, t);
  case FieldCurve::LOG_POSITIVE:
    return interp::log_positive(from, to, t);
  case FieldCurve::SHORTEST_PERIODIC:
    return interp::shortest_periodic(from, to, t, TWO_PI_F);
  case FieldCurve::SNAP:
    break;
  }
  return t < 1.0f ? from : to;
}

/**
 * @brief Interpolates every tabled field of a family.
 * @details Fields the table does not cover are left value-initialized; a
 * family with such members wraps this with its own overload.
 */
template <HasFields T>
HS_FLASH_INLINE inline T interpolate(const T &a, const T &b, float t) {
  T out{};
  for (const auto &field : T::FIELDS)
    out.*(field.member) =
        apply_curve(field.curve, a.*(field.member), b.*(field.member), t);
  return out;
}

/** @brief Whether every tabled field is finite and inside its range. */
template <HasFields T> HS_FLASH_INLINE inline bool valid(const T &value) {
  for (const auto &field : T::FIELDS) {
    const float sample = value.*(field.member);
    if (!std::isfinite(sample) || sample < field.min || sample > field.max)
      return false;
  }
  return true;
}

} // namespace Fields

} // namespace Pullback
