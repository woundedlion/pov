/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

/**
 * @file lens.h
 * @brief Sphere-to-sphere lens policies.
 */

namespace Pullback {

namespace Lens {

/** @brief Lens placeholder for an effect with no parameterized lens. */
struct NoLensParams {
  static constexpr std::array<Field<NoLensParams>, 0> FIELDS{};
};
static_assert(field_ids_unique<NoLensParams>());
/** @brief Lens parameters for the Mobius map (Pullback::Lens::Mobius). */
struct MobiusLensParams {
  /** Mobius coefficients; the default is the identity map. */
  MobiusParams mobius{0.7071067811865475f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                      0.7071067811865475f, 0.0f};
};

template <typename Params>
__attribute__((always_inline)) inline Vector mobius(const Vector &input,
                                                    const Params &params) {
  float px = input.x;
  float pz = input.z;
  float scale = 1.0f - input.y;
  if (px * px + pz * pz < STEREO_DIV_NUM_EPS_SQ &&
      scale < STEREO_DIV_NUM_EPS_SQ) {
    px = 1.0f;
    pz = 0.0f;
    scale = 0.0f;
  }
  const float n_re = params.a.re * px - params.a.im * pz + params.b.re * scale;
  const float n_im = params.a.re * pz + params.a.im * px + params.b.im * scale;
  const float m_re = params.c.re * px - params.c.im * pz + params.d.re * scale;
  const float m_im = params.c.re * pz + params.c.im * px + params.d.im * scale;
  const float n_sq = n_re * n_re + n_im * n_im;
  const float m_sq = m_re * m_re + m_im * m_im;
  const float denominator = n_sq + m_sq;
  if (denominator < STEREO_DIV_NUM_EPS_SQ)
    return Vector(0.0f, 1.0f, 0.0f);
  const float inverse = 1.0f / denominator;
  return Vector(2.0f * (n_re * m_re + n_im * m_im) * inverse,
                (n_sq - m_sq) * inverse,
                2.0f * (n_im * m_re - n_re * m_im) * inverse);
}

struct Glitch : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::glitch_lens(input);
  }
};

struct Twist : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::twist_lens(input);
  }
};

template <typename State> struct Mobius : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

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

  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &frame) {
    return mobius(input, State::params(frame));
  }
};

struct Kaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::kaleidoscope_lens(input);
  }
};

struct TetrahedralKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::TETRAHEDRAL_MIRRORS);
  }
};

struct OctahedralKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::OCTAHEDRAL_MIRRORS);
  }
};

struct DodecahedralKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::dodecahedral_kaleidoscope_lens(input);
  }
};

struct TriangularPrismKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::TRIANGULAR_PRISM_MIRRORS);
  }
};

struct SquarePrismKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(input,
                                                lenses::SQUARE_PRISM_MIRRORS);
  }
};

struct PentagonalPrismKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::PENTAGONAL_PRISM_MIRRORS);
  }
};

struct HexagonalPrismKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::HEXAGONAL_PRISM_MIRRORS);
  }
};

struct OctagonalPrismKaleidoscope : ApproximationDefaults {
  template <typename FrameState>
  __attribute__((always_inline)) static Vector apply(const Vector &input,
                                                     const FrameState &) {
    return lenses::polyhedral_kaleidoscope_lens(
        input, lenses::OCTAGONAL_PRISM_MIRRORS);
  }
};

} // namespace Lens

} // namespace Pullback
