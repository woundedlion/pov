/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/math/lenses.h"

#include "core/math/mobius.h"
#include "core/math/tangent_field.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"

/**
 * @file lens.h
 * @brief Sphere-to-sphere lens policies.
 */

namespace Pullback {

namespace Lens {

/** @brief Empty parameter family of the parameterless lens operators. */
struct NoLensParams {
  /// Empty parameter registry.
  static constexpr std::array<Field<NoLensParams>, 0> FIELDS{};
};
static_assert(field_ids_unique<NoLensParams>());
static_assert(field_defaults_in_range<NoLensParams>());
/** @brief Authored rotation and boost of a Mobius lens. */
struct MobiusGeneratorParams {
  float spin_x = 0.0f; ///< Rotation about x, radians.
  float spin_y = 0.0f; ///< Rotation about y, radians.
  float spin_z = 0.0f; ///< Rotation about z, radians.
  float push_x = 0.0f; ///< Boost toward x.
  float push_y = 0.0f; ///< Boost toward y.
  float push_z = 0.0f; ///< Boost toward z.
  /// Refusal message for an authored generator outside its range.
  static constexpr char RANGE_WARNING[] =
      "Mobius generator is outside its authored range.";
  /// Authored generator controls.
  static constexpr auto FIELDS = std::array{
      Field<MobiusGeneratorParams>{"spin-x", &MobiusGeneratorParams::spin_x,
                                   "Spin X", -math::TWO_PI_F, math::TWO_PI_F,
                                   FieldCurve::LERP},
      Field<MobiusGeneratorParams>{"spin-y", &MobiusGeneratorParams::spin_y,
                                   "Spin Y", -math::TWO_PI_F, math::TWO_PI_F,
                                   FieldCurve::LERP},
      Field<MobiusGeneratorParams>{"spin-z", &MobiusGeneratorParams::spin_z,
                                   "Spin Z", -math::TWO_PI_F, math::TWO_PI_F,
                                   FieldCurve::LERP},
      Field<MobiusGeneratorParams>{"push-x", &MobiusGeneratorParams::push_x,
                                   "Push X", -6.0f, 6.0f, FieldCurve::LERP},
      Field<MobiusGeneratorParams>{"push-y", &MobiusGeneratorParams::push_y,
                                   "Push Y", -6.0f, 6.0f, FieldCurve::LERP},
      Field<MobiusGeneratorParams>{"push-z", &MobiusGeneratorParams::push_z,
                                   "Push Z", -6.0f, 6.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<MobiusGeneratorParams>());
static_assert(field_defaults_in_range<MobiusGeneratorParams>());

/** @brief Exponentiates an authored generator at its animation phase.
 * @param params Authored spin and push.
 * @param angle Right-handed rotation about y, radians.
 * @param magnitude Horizontal generator multiplier.
 * @return Prepared Mobius coefficients. */
HS_FLASH_MEMBER inline math::MobiusParams
prepare_mobius(const MobiusGeneratorParams &params, float angle = 0.0f,
               float magnitude = 1.0f) {
  const float c = cosf(angle) * magnitude;
  const float s = sinf(angle) * magnitude;
  const math::Vector rotation{c * params.spin_x + s * params.spin_z,
                              params.spin_y,
                              c * params.spin_z - s * params.spin_x};
  const math::Vector boost{c * params.push_x + s * params.push_z, params.push_y,
                           c * params.push_z - s * params.push_x};
  return math::mobius_exp(
      math::tangent_field::to_mobius_generator(rotation, boost));
}

/** @brief Rotates an exponentiated generator about the sphere's y axis.
 * @param base Coefficients prepared at angle zero.
 * @param angle Right-handed y rotation, radians.
 * @return Coefficients of the rotated generator. */
inline math::MobiusParams rotate_mobius_y(const math::MobiusParams &base,
                                          float angle) {
  const math::Complex phase(cosf(angle), -sinf(angle));
  return {base.a, base.b * phase, base.c * phase.conjugate(), base.d};
}

/**
 * @brief Applies a Mobius transform to a sphere direction.
 * @param input Unit direction on the sphere.
 * @param params Mobius coefficients.
 * @return `math::mobius_transform(input, params)`.
 */
__attribute__((always_inline)) inline math::Vector
mobius(const math::Vector &input, const math::MobiusParams &params) {
  return math::mobius_transform(input, params);
}

/** @brief Glitch lens: `lenses::glitch_lens`. */
struct Glitch : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return `lenses::glitch_lens(input)`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::glitch_lens(input);
  }
};

/** @brief Parameters of an animated axial twist lens. */
struct TwistLensParams {
  float twist_angle = lenses::TWIST_ANGLE; ///< Total pole-to-pole radians.
  float distribution = lenses::TWIST_DISTRIBUTION; ///< Polar redistribution.
  float axis_phi = 0.0f;    ///< Axis polar offset, in radians.
  float axis_lambda = 0.0f; ///< Axis azimuth offset, in radians.
  float speed = 0.02f;      ///< Polar rotation per frame, in radians.
  float drift = 0.002f;     ///< Azimuth rotation per frame, in radians.
  /// Authored twist controls.
  static constexpr auto FIELDS = std::array{
      Field<TwistLensParams>{"twist-angle", &TwistLensParams::twist_angle,
                             "Twist Angle (rad)", -96.0f, 96.0f,
                             FieldCurve::LERP},
      Field<TwistLensParams>{"twist-distribution",
                             &TwistLensParams::distribution,
                             "Twist Distribution", 0.0f,
                             lenses::TWIST_DISTRIBUTION_MAX, FieldCurve::LERP},
      Field<TwistLensParams>{"axis-phi", &TwistLensParams::axis_phi,
                             "Axis Phi (rad)", 0.0f, math::PI_F,
                             FieldCurve::LERP},
      Field<TwistLensParams>{"axis-lambda", &TwistLensParams::axis_lambda,
                             "Axis Lambda (rad)", -math::PI_F, math::PI_F,
                             FieldCurve::SHORTEST_PERIODIC},
      Field<TwistLensParams>{"speed", &TwistLensParams::speed, "Speed", 0.0f,
                             0.05f, FieldCurve::LERP},
      Field<TwistLensParams>{"drift", &TwistLensParams::drift, "Drift", 0.0f,
                             0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<TwistLensParams>());
static_assert(field_defaults_in_range<TwistLensParams>());

/** @brief Accumulated twist-axis rotation, in radians modulo 2pi. */
struct TwistPhase {
  float phi = 0.0f;    ///< Polar angle.
  float lambda = 0.0f; ///< Azimuth angle.
};

/** @brief Prepared axis and angular profile. */
struct TwistPrepared {
  math::Vector axis{0.0f, 1.0f, 0.0f}; ///< Unit twist axis.
  /// Normalized angular profile.
  lenses::TwistProfile profile = lenses::prepare_twist_profile();
};

/** @brief Advances both twist-axis clocks.
 * @param phase Accumulated axis angles.
 * @param params Lens controls. */
HS_FLASH_INLINE inline void advance_twist(TwistPhase &phase,
                                          const TwistLensParams &params) {
  phase.phi = fmodf(phase.phi + params.speed, math::TWO_PI_F);
  phase.lambda = fmodf(phase.lambda + params.drift, math::TWO_PI_F);
}

/** @brief Prepares the animated twist axis and angular profile.
 * @param params Lens controls.
 * @param phase Accumulated axis rotation, added to the authored axis.
 * @return Per-frame lens state. */
HS_FLASH_INLINE inline TwistPrepared
prepare_twist(const TwistLensParams &params, const TwistPhase &phase) {
  return {
      math::Vector::from_spherical(params.axis_lambda + phase.lambda,
                                   params.axis_phi + phase.phi),
      lenses::prepare_twist_profile(params.twist_angle, params.distribution)};
}

/** @brief Applies a prepared animated twist.
 * @param input Unit direction.
 * @param prepared Axis and angular profile.
 * @return Twisted direction. */
HS_HOT_FLASH_MEMBER inline math::Vector twist(const math::Vector &input,
                                              const TwistPrepared &prepared) {
  if (prepared.axis.y == 1.0f)
    return lenses::twist_lens(input, prepared.profile);
  const float height = math::dot(input, prepared.axis);
  const float angle = lenses::twist_angle(height, prepared.profile);
  const float c = math::fast_cosf(angle);
  const float sine = math::fast_sinf(angle);
  return input * c + math::cross(input, prepared.axis) * sine +
         prepared.axis * (height * (1.0f - c));
}

/** @brief Animated twist policy supplied by a typed provider.
 * @tparam State Parameter and phase provider. */
template <typename State> struct AnimatedTwist : ApproximationDefaults {
  /// Deforming lenses contribute to displacement modulation.
  static constexpr bool CONTRIBUTES_DISPLACEMENT = true;
  using Binding = typename State::Binding;       ///< Provider binding.
  using FrameState = typename State::FrameState; ///< Bound frame.
  using Prepared = TwistPrepared;                ///< Per-frame twist.
  /// Whether the provider uses the candidate binding.
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, CandidateBinding>;
  /** @brief Prepares one frame or column.
   * @param frame Bound frame.
   * @return Prepared axis and profile. */
  static Prepared prepare(const FrameState &frame) {
    return prepare_twist(State::params(frame), State::phase(frame));
  }
  /** @brief Applies the prepared twist.
   * @param input Unit direction.
   * @param prepared Axis and profile.
   * @return Twisted direction. */
  static math::Vector apply(const math::Vector &input, const FrameState &,
                            const Prepared &prepared) {
    return twist(input, prepared);
  }
};

/** @brief Twist lens: `lenses::twist_lens` at its default angle and distribution. */
struct Twist : ApproximationDefaults {
  /// Deforming lenses contribute to displacement modulation.
  static constexpr bool CONTRIBUTES_DISPLACEMENT = true;
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return `lenses::twist_lens(input)`.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::twist_lens(input);
  }
};

/** @brief Mobius generator lens policy.
 * @tparam State Authored controls and clock provider. */
template <typename State> struct MobiusGenerator : ApproximationDefaults {
  /// Deforming lenses contribute to displacement modulation.
  static constexpr bool CONTRIBUTES_DISPLACEMENT = true;
  using Binding = typename State::Binding;       ///< Provider binding.
  using FrameState = typename State::FrameState; ///< Bound frame.
  using Prepared = math::MobiusParams;           ///< Exponentiated generator.
  /// Whether the provider uses the candidate frame binding.
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, CandidateBinding>;
  /** @brief Prepares one frame or column.
   * @param frame Bound frame.
   * @return Exponentiated generator. */
  static Prepared prepare(const FrameState &frame) {
    const auto &phase = State::phase(frame);
    if constexpr (requires { State::mobius_base(frame); })
      if (const auto *base = State::mobius_base(frame))
        return rotate_mobius_y(*base, phase.angle);
    return prepare_mobius(State::params(frame), phase.angle, phase.magnitude);
  }
  /** @brief Applies the prepared map.
   * @param input Unit sphere direction.
   * @param prepared Mobius coefficients.
   * @return Transformed direction. */
  static math::Vector apply(const math::Vector &input, const FrameState &,
                            const Prepared &prepared) {
    return mobius(input, prepared);
  }
};

/** @brief Azimuthal kaleidoscope: `lenses::kaleidoscope_lens`. */
struct Kaleidoscope : ApproximationDefaults {
  /// Whether displayed columns are followed by black.
  static constexpr bool STROBE_COLUMNS = true;
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the wedge.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::kaleidoscope_lens(input);
  }
};

/** @brief Kaleidoscope folding using a polyhedral mirror set.
 *  @tparam MIRRORS Chamber mirror planes. */
template <const auto &MIRRORS> struct PolyhedralKaleidoscope : Kaleidoscope {
  /**
   * @brief Folds a sphere direction into the mirror chamber.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit sphere direction.
   * @return Symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input, MIRRORS);
  }
};

/** @brief Tetrahedral chamber folding. */
using TetrahedralKaleidoscope =
    PolyhedralKaleidoscope<lenses::TETRAHEDRAL_MIRRORS>;
/** @brief Octahedral chamber folding. */
using OctahedralKaleidoscope =
    PolyhedralKaleidoscope<lenses::OCTAHEDRAL_MIRRORS>;

/** @brief Kaleidoscope folding into the dodecahedral chamber. */
struct DodecahedralKaleidoscope : Kaleidoscope {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::dodecahedral_kaleidoscope_lens(input);
  }
};

/** @brief Canonical dodecahedral lens with validated interior matrix lookup. */
struct DodecahedralTableKaleidoscope : DodecahedralKaleidoscope {
  /** @brief Preserves literal endpoints and validates interior matrix folds.
   * @tparam FrameState Prepared column context.
   * @param input Unit view direction.
   * @param frame Column context.
   * @return Canonical chamber direction. */
  template <typename FrameState>
  HS_HOT_INLINE static math::Vector apply(const math::Vector &input,
                                          const FrameState &frame) {
    if constexpr (requires { frame.temporal.fraction; })
      if (frame.temporal.fraction > 0.0f && frame.temporal.fraction < 1.0f) {
        const auto *matrices = lenses::DodecahedralTable::GROUP.values.data();
        const uint8_t *faces = nullptr;
        if constexpr (requires { frame.dodecahedral_matrices; })
          if (frame.dodecahedral_matrices)
            matrices = frame.dodecahedral_matrices;
        if constexpr (requires { frame.dodecahedral_faces; })
          faces = frame.dodecahedral_faces;
        return lenses::dodecahedral_table_lens(input, matrices, faces);
      }
    return lenses::dodecahedral_kaleidoscope_lens(input);
  }
};

/** @brief Tetrahedral folding of a direction already in the dodecahedral
 * chamber: `lenses::tetrahedral_from_dodecahedral_chamber`. */
struct TetrahedralFromDodecahedralKaleidoscope : TetrahedralKaleidoscope {
  /**
   * @brief Folds a dodecahedral chamber direction into the tetrahedral chamber.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction in the dodecahedral chamber.
   * @return Direction in the tetrahedral chamber.
   */
  template <typename FrameState>
  HS_HOT_INLINE static math::Vector apply(const math::Vector &input,
                                          const FrameState &) {
    return lenses::tetrahedral_from_dodecahedral_chamber(input);
  }
};

/** @brief Triangular-prism chamber folding. */
using TriangularPrismKaleidoscope =
    PolyhedralKaleidoscope<lenses::TRIANGULAR_PRISM_MIRRORS>;
/** @brief Square-prism chamber folding. */
using SquarePrismKaleidoscope =
    PolyhedralKaleidoscope<lenses::SQUARE_PRISM_MIRRORS>;
/** @brief Pentagonal-prism chamber folding. */
using PentagonalPrismKaleidoscope =
    PolyhedralKaleidoscope<lenses::PENTAGONAL_PRISM_MIRRORS>;
/** @brief Hexagonal-prism chamber folding. */
using HexagonalPrismKaleidoscope =
    PolyhedralKaleidoscope<lenses::HEXAGONAL_PRISM_MIRRORS>;
/** @brief Octagonal-prism chamber folding. */
using OctagonalPrismKaleidoscope =
    PolyhedralKaleidoscope<lenses::OCTAGONAL_PRISM_MIRRORS>;

} // namespace Lens

} // namespace Pullback
