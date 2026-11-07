/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <string_view>
#include <type_traits>

#include "core/control/fields.h"
#include "core/render/pullback/contract.h"

/**
 * @file fields.h
 * @brief Per-family field descriptors: one table per parameter family drives
 *        slider registration, preset interpolation and snapshot validation.
 */

namespace Pullback {

/**
 * @brief One topology enum8 of a parameter family: a structural switch that
 *        selects among eagerly constructed variants at run time.
 * @tparam Owner The family struct the enum8 lives in.
 */
template <typename Owner> struct TopologyField {
  const char *id;               /**< Stable machine id, kebab-case. */
  uint8_t Owner::*member;       /**< The enum8 member of @p Owner. */
  const char *const *value_ids; /**< Kebab-case value spellings. */
  uint8_t value_count;          /**< Entries in `value_ids`. */
  uint8_t def;                  /**< Default value index. */

  /**
   * @brief Binds an enum8 member to its value spellings.
   * @tparam N Number of values.
   * @param id Stable machine id, kebab-case.
   * @param member The enum8 member of @p Owner.
   * @param value_ids Kebab-case value spellings, one per value.
   * @param def Default value index.
   */
  template <size_t N>
  constexpr TopologyField(const char *id, uint8_t Owner::*member,
                          const char *const (&value_ids)[N], uint8_t def)
      : id(id), member(member), value_ids(value_ids), value_count(N), def(def) {
    static_assert(N > 0 && N <= UINT8_MAX);
  }

  /** @brief Binds an inherited selector to the derived owner.
   * @tparam Base Original owner.
   * @param field Selector of the base family. */
  template <typename Base>
    requires std::is_convertible_v<uint8_t Base::*, uint8_t Owner::*>
  constexpr explicit TopologyField(const TopologyField<Base> &field)
      : id(field.id), member(field.member), value_ids(field.value_ids),
        value_count(field.value_count), def(field.def) {}
};

/** @brief Interpolation curves shared with described controls. */
using FieldCurve = Control::FieldCurve;

/**
 * @brief One entry of an operator's registered parameter schema — the full
 *        union across topology variants, topology enum8s included.
 */
struct ParamFieldInfo {
  const char *id;              /**< Stable machine id, kebab-case. */
  const char *display_name;    /**< Nullable; cosmetic. */
  float min;                   /**< Float range minimum. */
  float max;                   /**< Float range maximum. */
  float def;                   /**< Float default. */
  FieldCurve curve;            /**< Preset-transition curve. */
  bool topology;               /**< Whether this is a topology enum8. */
  uint8_t enum_count;          /**< 0 marks a float field. */
  uint8_t enum_def;            /**< Default enum value index. */
  const char *const *enum_ids; /**< Enum value spellings; null for floats. */
  /** Id of the topology enum8 selecting this field; null when every variant
      reads it. */
  const char *gated_by;
  uint16_t gate_values; /**< Bit per live gating value index. */
};

/** @brief When a field's slider is registered. */
enum class FieldGate : uint8_t {
  ALWAYS,
  ANIMATED_PROJECTION, /**< Only when the effect animates its projection. */
  CENTRAL_MERIDIAN,    /**< Only when the effect's projection reads it. */
  SINGULARITY_FADE,    /**< Only when the effect's projection has a singular
                            locus to attenuate. */
  COORDINATE_SCALE,    ///< Atlas projections with an authored coordinate scale.
  LAYOUT_SCROLL,       ///< Peirce strip layouts.
  STANDARD_PARALLEL,   ///< Bonne reference latitude.
  AREA_CORRECTION,     ///< Polyhedral atlas blend toward equal area.
  OUTER_CAMERA,        ///< Shared outer-camera stage.
  AFFINE_ROTATES,      ///< Composed affine stages that rotate.
  AFFINE_TRANSLATES,   ///< Composed affine stages that translate.
  AFFINE_SCALES_SHEARS, ///< Composed affine stages that scale and shear.
  AFFINE_REPEAT, /**< Composed affine stages that translate, or that measure
                      path length over the next stage's lattice. */
  AFFINE_FALLBACK_LATTICE /**< Composed affine stages that translate where
                               the next stage reports no lattice. */
};

/**
 * @brief Which topology enum8 of the same family decides whether a field is
 *        read, and the value indices that keep it live.
 * @details A null `field` marks a field every variant reads.
 */
struct TopologyGate {
  const char *field = nullptr; /**< Topology field id, or null. */
  uint16_t values = 0;         /**< Bit per live topology value index. */
};

/** Declared only: a call makes the enclosing constant evaluation ill-formed. */
void live_values_index_exceeds_gate_width();

/** @brief TopologyGate::values over the enumerators in @p values.
 *  @tparam Values Enumerator types.
 *  @param values Live topology values; each index must fit the gate width.
 *  @return Mask with one bit per value index. */
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
  float Owner::*member; /**< The float member this field drives. */
  const char *name;     /**< Slider display name, or null. */
  float min;            /**< Inclusive lower bound. */
  float max;            /**< Inclusive upper bound. */
  FieldCurve curve = FieldCurve::LERP; /**< Transition curve. */
  FieldGate gate = FieldGate::ALWAYS;  /**< Slider registration gate. */
  TopologyGate topology_gate{};        /**< Topology value(s) that read it. */

  /** @brief The equivalent animated `Control::Field` descriptor.
   * @return Descriptor carrying id, member, name, range and curve. */
  HS_HOT_INLINE constexpr Control::Field<Owner, float> description() const {
    return {.id = id,
            .member = member,
            .name = name,
            .spec = {.min = min, .max = max, .animated = true},
            .curve = curve};
  }
};

/**
 * @brief Rebuilds a base family's field table over a derived family.
 * @details A `float Base::*` converts implicitly to `float Combined::*`, so
 * the concatenation stays a homogeneous Field<Combined> array.
 * @tparam Combined Derived family type.
 * @tparam BaseOwner Base family type.
 * @tparam N Base field count.
 * @tparam M Extra field count.
 * @param base Base family's field table.
 * @param extra Fields the derived family adds.
 * @return @p base followed by @p extra.
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

/**
 * @brief The shared `edge-width` field: range [0, 16], linear.
 * @details Uses `FieldGate::ALWAYS`; @p topology_gate limits applicable
 * topology values, and a null @p name leaves the field unregistered.
 * @tparam Owner The family struct the field belongs to.
 * @param member The float member it drives.
 * @param name Slider display name; null leaves the field unregistered.
 * @param topology_gate Topology value(s) that read it.
 * @return The field descriptor.
 */
template <typename Owner>
constexpr Field<Owner> edge_width_field(float Owner::*member,
                                        const char *name = "Edge Width",
                                        TopologyGate topology_gate = {}) {
  return {"edge-width",      member,       name, 0.0f, 16.0f, FieldCurve::LERP,
          FieldGate::ALWAYS, topology_gate};
}

/** @brief Whether @p T carries a field-descriptor table. */
template <typename T>
concept HasFields = requires { T::FIELDS; };

/** @brief Whether numeric and extra schema rows use registered pullback curves.
 * @tparam Family Parameter family with a field table.
 * @return False for midpoint, extrapolating or unknown curves. */
template <HasFields Family> consteval bool field_curves_supported() {
  const auto supported = [](const auto &fields) {
    for (const auto &field : fields)
      switch (field.curve) {
      case FieldCurve::LERP:
      case FieldCurve::LOG_POSITIVE:
      case FieldCurve::SHORTEST_PERIODIC:
      case FieldCurve::SHORTEST_TURN:
      case FieldCurve::SNAP:
        break;
      default:
        return false;
      }
    return true;
  };
  if (!supported(Family::FIELDS))
    return false;
  if constexpr (requires { Family::EXTRA_SCHEMA; })
    if (!supported(Family::EXTRA_SCHEMA))
      return false;
  return true;
}

/** @brief Whether every tabled field id is non-null and unique in the family.
 *  @tparam Family Parameter family.
 *  @return True when the ids are valid. */
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
    range.
    @tparam Family Parameter family.
    @return True when every default is in range. */
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

/**
 * @brief Interpolates every tabled field of a family.
 * @details Members the table does not cover — topology enum8s, non-float
 * carriers — snap: they hold @p a's value until progress reaches 1.
 * @tparam T Parameter family.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param t Progress fraction.
 * @return The interpolated family.
 */
template <HasFields T>
HS_FLASH_INLINE inline T interpolate(const T &a, const T &b, float t) {
  T out = t < 1.0f ? a : b;
  for (const auto &field : T::FIELDS)
    field.description().interpolate(out, a, b, t);
  return out;
}

/** @brief Whether every tabled field is finite and inside its range.
 *  @tparam T Parameter family.
 *  @param value Family to check.
 *  @return True when every field passes. */
template <HasFields T> HS_FLASH_INLINE inline bool valid(const T &value) {
  for (const auto &field : T::FIELDS) {
    if (!field.description().valid(value))
      return false;
  }
  if constexpr (requires { T::TOPOLOGY; })
    for (const auto &field : T::TOPOLOGY)
      if (value.*field.member >= field.value_count)
        return false;
  return true;
}

} // namespace Fields

} // namespace Pullback
