/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "math/lenses.h"

#include "math/mobius.h"
#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

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
  /// Message for coefficients that fail nondegenerate().
  static constexpr const char *DEGENERATE_WARNING =
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

/** @brief Twist lens: `lenses::twist_lens` at its default rate. */
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

/** @brief Kaleidoscope folding into the tetrahedral chamber. */
struct TetrahedralKaleidoscope : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::TETRAHEDRAL_MIRRORS);
  }
};

/** @brief Kaleidoscope folding into the octahedral chamber. */
struct OctahedralKaleidoscope : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::OCTAHEDRAL_MIRRORS);
  }
};

/** @brief Kaleidoscope folding into the dodecahedral chamber. */
struct DodecahedralKaleidoscope : ApproximationDefaults {
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

/** @brief Kaleidoscope folding into the triangular-prism chamber. */
struct TriangularPrismKaleidoscope : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::TRIANGULAR_PRISM_MIRRORS);
  }
};

/** @brief Kaleidoscope folding into the square-prism chamber. */
struct SquarePrismKaleidoscope : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::SQUARE_PRISM_MIRRORS);
  }
};

/** @brief Kaleidoscope folding into the pentagonal-prism chamber. */
struct PentagonalPrismKaleidoscope : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::PENTAGONAL_PRISM_MIRRORS);
  }
};

/** @brief Kaleidoscope folding into the hexagonal-prism chamber. */
struct HexagonalPrismKaleidoscope : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::HEXAGONAL_PRISM_MIRRORS);
  }
};

/** @brief Kaleidoscope folding into the octagonal-prism chamber. */
struct OctagonalPrismKaleidoscope : ApproximationDefaults {
  /**
   * @brief Applies the lens to a sphere direction.
   * @tparam FrameState Frame state type; unused.
   * @param input Unit direction on the sphere.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <typename FrameState>
  __attribute__((always_inline)) static math::Vector
  apply(const math::Vector &input, const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::OCTAGONAL_PRISM_MIRRORS);
  }
};

} // namespace Lens

} // namespace Pullback
