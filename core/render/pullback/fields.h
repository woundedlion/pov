/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <string_view>

#include "control/fields.h"
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
  SHORTEST_TURN,     /**< Shortest arc on a unit period, for fields in
                          turns. */
  SNAP               /**< Holds the start value until progress reaches 1. */
};

constexpr Control::FieldCurve control_curve(FieldCurve curve) {
  switch (curve) {
  case FieldCurve::LERP:
    return Control::FieldCurve::LERP;
  case FieldCurve::LOG_POSITIVE:
    return Control::FieldCurve::LOG_POSITIVE;
  case FieldCurve::SHORTEST_PERIODIC:
    return Control::FieldCurve::SHORTEST_PERIODIC;
  case FieldCurve::SHORTEST_TURN:
    return Control::FieldCurve::SHORTEST_TURN;
  case FieldCurve::SNAP:
    return Control::FieldCurve::SNAP;
  }
  return Control::FieldCurve::SNAP;
}

/** @brief When a field's slider is registered. */
enum class FieldGate : uint8_t {
  ALWAYS,
  ANIMATED_PROJECTION, /**< Only when the effect animates its projection. */
  CENTRAL_MERIDIAN,    /**< Only when the effect's projection reads it. */
  SINGULARITY_FADE     /**< Only when the effect's projection has a singular
                            locus to attenuate. */
};

/**
 * @brief Which topology enum8 of the same family decides whether a field is
 *        read, and the value indices that keep it live.
 * @details A null `field` marks a field every variant reads.
 */
struct TopologyGate {
  const char *field = nullptr;
  uint16_t values = 0; /**< Bit per live topology value index. */
};

/** Declared only: a call makes the enclosing constant evaluation ill-formed. */
void live_values_index_exceeds_gate_width();

/** @brief TopologyGate::values over the enumerators in @p values. */
template <typename... Values> consteval uint16_t live_values(Values... values) {
  constexpr unsigned GATE_BITS =
      std::numeric_limits<decltype(TopologyGate::values)>::digits;
  uint16_t mask = 0;
  for (const unsigned index : {static_cast<unsigned>(values)...}) {
    if (index >= GATE_BITS)
      live_values_index_exceeds_gate_width();
    mask = static_cast<uint16_t>(mask | (1U << index));
  }
  return mask;
}

/**
 * @brief One scalar field of a parameter family.
 * @details `name == nullptr` means no table-supplied display name. Table-driven
 * registration skips it; a consumer may register it under its own name.
 * It is still interpolated and validated.
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
  TopologyGate topology_gate{};

  constexpr Control::Field<Owner, float> description() const {
    return {.id = id,
            .member = member,
            .name = name,
            .spec = {.min = min, .max = max, .animated = true},
            .curve = control_curve(curve)};
  }
};

/**
 * @brief Rebuilds a base family's field table over a derived family.
 * @details A `float Base::*` converts implicitly to `float Combined::*`, so
 * the concatenation stays a homogeneous Field<Combined> array.
 */
template <typename Combined, typename BaseOwner, size_t N, size_t M>
consteval std::array<Field<Combined>, N + M>
concat_fields(const std::array<Field<BaseOwner>, N> &base,
              const std::array<Field<Combined>, M> &extra) {
  std::array<Field<Combined>, N + M> out{};
  for (size_t index = 0; index < N; ++index)
    out[index] = Field<Combined>{base[index].id,   base[index].member,
                                 base[index].name, base[index].min,
                                 base[index].max,  base[index].curve,
                                 base[index].gate, base[index].topology_gate};
  for (size_t index = 0; index < M; ++index)
    out[N + index] = extra[index];
  return out;
}

template <typename Owner>
constexpr Field<Owner> edge_width_field(float Owner::*member,
                                        const char *name = "Edge Width",
                                        TopologyGate topology_gate = {}) {
  return {"edge-width",      member,       name, 0.0f, 1.0f, FieldCurve::LERP,
          FieldGate::ALWAYS, topology_gate};
}

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

/** @brief Whether every tabled field's default lies within its declared
    range. */
template <HasFields Family> consteval bool field_defaults_in_range() {
  constexpr Family DEFAULTS{};
  for (const auto &field : Family::FIELDS) {
    const float value = DEFAULTS.*(field.member);
    if (!(value >= field.min && value <= field.max))
      return false;
  }
  return true;
}

namespace Fields {

HS_FLASH_INLINE inline float apply_curve(FieldCurve curve, float from, float to,
                                         float t) {
  return Control::apply_curve(control_curve(curve), from, to, t);
}

/**
 * @brief Interpolates every tabled field of a family.
 * @details Members the table does not cover — topology enum8s, non-float
 * carriers — snap: they hold @p a's value until progress reaches 1.
 */
template <HasFields T>
HS_FLASH_INLINE inline T interpolate(const T &a, const T &b, float t) {
  T out = t < 1.0f ? a : b;
  for (const auto &field : T::FIELDS)
    field.description().interpolate(out, a, b, t);
  return out;
}

/** @brief Whether every tabled field is finite and inside its range. */
template <HasFields T> HS_FLASH_INLINE inline bool valid(const T &value) {
  for (const auto &field : T::FIELDS) {
    if (!field.description().valid(value))
      return false;
  }
  return true;
}

} // namespace Fields

} // namespace Pullback
