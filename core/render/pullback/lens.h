/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/math/lenses.h"

#include "core/math/mobius.h"
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
/** @brief Lens parameters for the Mobius map (Pullback::Lens::Mobius). */
struct MobiusLensParams {
  /// Message for coefficients that fail nondegenerate(); one object, so its
  /// address identifies the refusal in every translation unit.
  static constexpr char DEGENERATE_WARNING[] =
      "Mobius coefficients must have |ad - bc| >= 0.001.";
  static constexpr float COEFFICIENT_LIMIT = 4.0f;  ///< Max abs re/im part.
  static constexpr float MOBIUS_MIN_DET_SQ = 1e-6f; ///< Minimum |ad - bc|^2.

  /**
   * @brief Squared magnitude of the complex determinant ad - bc.
   * @param params Mobius coefficients.
   * @return |ad - bc|^2.
   */
  static constexpr float
  determinant_magnitude_squared(const math::MobiusParams &params) {
    const float re = params.a.re * params.d.re - params.a.im * params.d.im -
                     params.b.re * params.c.re + params.b.im * params.c.im;
    const float im = params.a.re * params.d.im + params.a.im * params.d.re -
                     params.b.re * params.c.im - params.b.im * params.c.re;
    return re * re + im * im;
  }

  /**
   * @brief Whether the coefficients define an invertible map.
   * @param params Mobius coefficients.
   * @return True when |ad - bc|^2 >= `MOBIUS_MIN_DET_SQ`.
   */
  static constexpr bool nondegenerate(const math::MobiusParams &params) {
    return determinant_magnitude_squared(params) >= MOBIUS_MIN_DET_SQ;
  }

  /** Mobius coefficients; the default is the identity map. */
  math::MobiusParams mobius{0.7071067811865475f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                            0.7071067811865475f, 0.0f};
};

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
  float speed = 0.02f;  ///< Polar rotation per frame, in radians.
  float drift = 0.002f; ///< Azimuth rotation per frame, in radians.
  /// Authored twist controls.
  static constexpr auto FIELDS = std::array{
      Field<TwistLensParams>{"twist-angle", &TwistLensParams::twist_angle,
                             "Twist Angle (rad)", -96.0f, 96.0f,
                             FieldCurve::LERP},
      Field<TwistLensParams>{"twist-distribution",
                             &TwistLensParams::distribution,
                             "Twist Distribution", 0.0f,
                             lenses::TWIST_DISTRIBUTION_MAX, FieldCurve::LERP},
      Field<TwistLensParams>{"speed", &TwistLensParams::speed, "Speed", 0.0f,
                             0.05f, FieldCurve::LERP},
      Field<TwistLensParams>{"drift", &TwistLensParams::drift, "Drift", 0.0f,
                             0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<TwistLensParams>());
static_assert(field_defaults_in_range<TwistLensParams>());

/** @brief Twist-axis angles, in radians modulo 2pi. */
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
inline void advance_twist(TwistPhase &phase, const TwistLensParams &params) {
  phase.phi = fmodf(phase.phi + params.speed, math::TWO_PI_F);
  phase.lambda = fmodf(phase.lambda + params.drift, math::TWO_PI_F);
}

/** @brief Prepares the animated twist axis and angular profile.
 * @param params Lens controls.
 * @param phase Accumulated axis angles.
 * @return Per-frame lens state. */
inline TwistPrepared prepare_twist(const TwistLensParams &params,
                                   const TwistPhase &phase) {
  return {
      math::Vector::from_spherical(phase.lambda, phase.phi),
      lenses::prepare_twist_profile(params.twist_angle, params.distribution)};
}

/** @brief Applies a prepared animated twist.
 * @param input Unit direction.
 * @param prepared Axis and angular profile.
 * @return Twisted direction. */
inline math::Vector twist(const math::Vector &input,
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
   * @param frame Bound frame.
   * @param prepared Axis and profile.
   * @return Twisted direction. */
  static math::Vector apply(const math::Vector &input, const FrameState &,
                            const Prepared &prepared) {
    return twist(input, prepared);
  }
};

/** @brief Twist lens: `lenses::twist_lens` at its default angle and distribution. */
struct Twist : ApproximationDefaults {
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

/**
 * @brief Mobius lens policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame) accessors.
 */
template <typename State> struct Mobius : ApproximationDefaults {
  using Binding = typename State::Binding;       ///< The provider's binding.
  using FrameState = typename State::FrameState; ///< The provider's frame.

  /**
   * @brief Whether `State` is a provider for `CandidateBinding` with every
   * accessor this policy reads.
   * @tparam CandidateBinding Chain binding to check against.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        State::params(frame).a;
        State::params(frame).b;
        State::params(frame).c;
        State::params(frame).d;
      } &&
      std::is_lvalue_reference_v<decltype(State::params(
          std::declval<const typename CandidateBinding::FrameState &>()))>;

  /**
   * @brief Applies this frame's Mobius transform.
   * @param input Unit direction on the sphere.
   * @param frame Frame state the coefficients are read from.
   * @return The transformed direction.
   */
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &frame) {
    return mobius(input, State::params(frame));
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
