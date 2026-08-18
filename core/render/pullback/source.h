/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"
#include "math/3dmath.h"

/**
 * @file source.h
 * @brief Scalar source-field policies.
 */

namespace Pullback {

namespace Source {

/**
 * @brief Source parameters for the coupled sine grid
 *        (Pullback::Source::Grid).
 */
struct GridSourceParams {
  float pattern_freq = 1.0f; /**< Scale applied to the warped plane
                                   coordinates before the grid is sampled. */
  float speed = 0.0f;        /**< Per-frame advance of the primary phase. */
  float complexity = 0.0f;   /**< Amount of cross-axis coupling folded into the
                                   grid coordinates. */
  float pattern_mix = 0.0f;  /**< Blend from the coupled pattern at 0 to the
                                   direct sine product at 1. */
  float secondary_rate = 0.0f; /**< Secondary phase rate, as a multiple of
                                    `speed`. */
  float angle_rate = 0.0f;     /**< Per-frame advance of the source rotation. */

  static constexpr auto FIELDS = std::array{
      Field<GridSourceParams>{&GridSourceParams::pattern_freq, "Pattern Freq",
                              0.1f, 20.0f, FieldCurve::LERP},
      Field<GridSourceParams>{&GridSourceParams::speed, "Speed", 0.0f, 0.5f,
                              FieldCurve::LERP},
      Field<GridSourceParams>{&GridSourceParams::complexity, "Complexity", 0.0f,
                              3.0f, FieldCurve::LERP},
      Field<GridSourceParams>{&GridSourceParams::pattern_mix, "Pattern Mix",
                              0.0f, 1.0f, FieldCurve::LERP},
      Field<GridSourceParams>{&GridSourceParams::secondary_rate, "Drift", 0.0f,
                              1.25f, FieldCurve::LERP},
      Field<GridSourceParams>{&GridSourceParams::angle_rate,
                              "Source Angle Speed", 0.0f, 0.05f,
                              FieldCurve::LERP},
  };
};

/**
 * @brief Source parameters for the two-wave interference field
 *        (Pullback::Source::TwinWave).
 */
struct TwinWaveSourceParams {
  float pattern_freq = 1.0f;   /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;          /**< Per-frame advance of the primary phase. */
  float secondary_rate = 0.0f; /**< Secondary phase rate, as a multiple of
                                    `speed`. */
  float angle_rate = 0.0f; /**< Per-frame advance of the angle between the two
                                waves. */

  static constexpr auto FIELDS = std::array{
      Field<TwinWaveSourceParams>{&TwinWaveSourceParams::pattern_freq,
                                  "Pattern Freq", 0.1f, 20.0f,
                                  FieldCurve::LERP},
      Field<TwinWaveSourceParams>{&TwinWaveSourceParams::speed, "Speed", 0.0f,
                                  0.5f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{&TwinWaveSourceParams::secondary_rate,
                                  "Drift", 0.0f, 1.25f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{&TwinWaveSourceParams::angle_rate,
                                  "Source Angle Speed", 0.0f, 0.05f,
                                  FieldCurve::LERP},
  };
};

/**
 * @brief Source parameters for the rotating spiral field
 *        (Pullback::Source::Spiral).
 */
struct SpiralSourceParams {
  float pattern_freq = 1.0f; /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;        /**< Per-frame advance of the primary phase. */
  float angle_rate = 0.0f;   /**< Per-frame advance of the spiral rotation. */

  static constexpr auto FIELDS = std::array{
      Field<SpiralSourceParams>{&SpiralSourceParams::pattern_freq,
                                "Pattern Freq", 0.1f, 20.0f, FieldCurve::LERP},
      Field<SpiralSourceParams>{&SpiralSourceParams::speed, "Speed", 0.0f, 0.5f,
                                FieldCurve::LERP},
      Field<SpiralSourceParams>{&SpiralSourceParams::angle_rate,
                                "Source Angle Speed", 0.0f, 0.05f,
                                FieldCurve::LERP},
  };
};

/**
 * @brief Source parameters for the noise-contour sources
 *        (Pullback::Source::ProjectedNoise and SphericalNoise).
 */
struct NoiseSourceParams {
  float noise_scale = 1.0f;    /**< Spatial scale of the sampled field. */
  float noise_contrast = 0.0f; /**< Contour sharpening applied to the sample. */
  float noise_time_rate = 0.0f; /**< Per-frame advance of the noise time
                                     coordinate. */

  static constexpr auto FIELDS = std::array{
      Field<NoiseSourceParams>{&NoiseSourceParams::noise_scale,
                               "Source Noise Scale", 1.0f / 64.0f, 64.0f,
                               FieldCurve::LERP},
      Field<NoiseSourceParams>{&NoiseSourceParams::noise_contrast,
                               "Source Noise Contrast", 0.0f, 8.0f,
                               FieldCurve::LERP},
      Field<NoiseSourceParams>{&NoiseSourceParams::noise_time_rate,
                               "Source Noise Speed", -1.0f / 64.0f,
                               1.0f / 64.0f, FieldCurve::LERP},
  };
};

/**
 * @brief Source parameters for the per-cell primitive lattice
 *        (Pullback::Source::PrimitiveLattice).
 */
struct LatticeSourceParams {
  float lattice_cell_scale = 1.0f;  /**< Lattice cells per plane unit. */
  float lattice_shape_blend = 0.0f; /**< Cell primitive, from a circle at 0 to a
                                         rounded square at 1. */
  float lattice_softness = 0.05f;   /**< Half-width of the ramp across the
                                         primitive's boundary. */
  float lattice_radius = 0.25f;     /**< Primitive radius in cell units. */

  static constexpr auto FIELDS = std::array{
      Field<LatticeSourceParams>{&LatticeSourceParams::lattice_cell_scale,
                                 "Lattice Cell Scale", 1.0f / 64.0f, 8.0f,
                                 FieldCurve::LOG_POSITIVE},
      Field<LatticeSourceParams>{&LatticeSourceParams::lattice_shape_blend,
                                 "Lattice Shape", 0.0f, 1.0f, FieldCurve::LERP},
      Field<LatticeSourceParams>{&LatticeSourceParams::lattice_softness,
                                 "Lattice Softness", 1.0f / 1024.0f, 1.0f,
                                 FieldCurve::LOG_POSITIVE},
      Field<LatticeSourceParams>{&LatticeSourceParams::lattice_radius,
                                 "Lattice Radius", 1.0f / 64.0f, 0.49f,
                                 FieldCurve::LERP},
  };
};

/** @brief The source stage's phases, resolved once per frame. */
struct PreparedSource {
  float primary;   /**< Primary phase, wrapped into [0,2pi). */
  float secondary; /**< Secondary phase, wrapped into [0,2pi). */
  float angle;     /**< Source rotation, in radians. */
  float angle_cos; /**< Cosine of `angle`. */
  float angle_sin; /**< Sine of `angle`. */
};

/** @brief Wraps this frame's source phases with the rotation's cosine pair. */
HS_FLASH_INLINE inline PreparedSource prepare(float primary, float secondary,
                                              float angle) {
  return {primary, secondary, angle, fast_cosf(angle), fast_sinf(angle)};
}

template <typename State, typename Binding>
concept ParamsProvider =
    Detail::ProviderFor<State, Binding> &&
    requires(const typename Binding::FrameState &frame) {
      State::params(frame);
    } &&
    std::is_lvalue_reference_v<decltype(State::params(
        std::declval<const typename Binding::FrameState &>()))> &&
    std::is_const_v<std::remove_reference_t<decltype(State::params(
        std::declval<const typename Binding::FrameState &>()))>>;

template <typename State, typename Binding>
concept StateProvider = ParamsProvider<State, Binding> &&
                        requires(const typename Binding::FrameState &frame) {
                          State::prepare(frame);
                        };

template <typename Prepared>
HS_FLASH_MEMBER inline float twin_wave(const Complex &input,
                                       const Prepared &prepared) {
  const float rotated =
      input.re * prepared.angle_cos + input.im * prepared.angle_sin;
  return 0.5f * (fast_sinf(input.re + prepared.primary) +
                 fast_sinf(rotated + prepared.primary));
}

template <typename Prepared>
HS_FLASH_MEMBER inline float rings(const Complex &input,
                                   const Prepared &prepared) {
  return fast_sinf(sqrtf(input.re * input.re + input.im * input.im) -
                   prepared.primary);
}

template <typename Prepared>
HS_FLASH_MEMBER inline float spiral(const Complex &input,
                                    const Prepared &prepared) {
  const float radius = sqrtf(input.re * input.re + input.im * input.im);
  const float azimuth = fast_atan2(input.im, input.re);
  return fast_sinf(radius - 3.0f * (azimuth + prepared.angle) -
                   prepared.primary);
}

template <typename Params, typename Prepared>
HS_FLASH_MEMBER inline float grid(const Complex &input, const Params &params,
                                  const Prepared &prepared) {
  const float x = input.re * prepared.angle_cos + input.im * prepared.angle_sin;
  const float y =
      -input.re * prepared.angle_sin + input.im * prepared.angle_cos;
  if (params.pattern_mix == 1.0f)
    return fast_sinf(x + prepared.primary) * fast_cosf(y - prepared.secondary);
  float re = x;
  float im = y;
  if (params.complexity != 0.0f) {
    re += params.complexity * fast_sinf(y + prepared.primary);
    im += params.complexity * fast_cosf(x - prepared.secondary);
  }
  const float coupled = fast_sinf(re) * fast_cosf(im);
  if (params.pattern_mix == 0.0f)
    return coupled;
  const float direct =
      fast_sinf(x + prepared.primary) * fast_cosf(y - prepared.secondary);
  return hs::lerp(coupled, direct, params.pattern_mix);
}

template <typename Params>
HS_FLASH_MEMBER inline float primitive_lattice(const Complex &input,
                                               const Params &params) {
  const float x = wrap_t(params.lattice_cell_scale * input.re + 0.5f) - 0.5f;
  const float y = wrap_t(params.lattice_cell_scale * input.im + 0.5f) - 0.5f;
  const float circle = sqrtf(x * x + y * y) - params.lattice_radius;
  const float bx = fabsf(x) - params.lattice_radius;
  const float by = fabsf(y) - params.lattice_radius;
  const float square = sqrtf(std::max(bx, 0.0f) * std::max(bx, 0.0f) +
                             std::max(by, 0.0f) * std::max(by, 0.0f)) +
                       std::min(std::max(bx, by), 0.0f);
  const float distance = hs::lerp(circle, square, params.lattice_shape_blend);
  return 1.0f - 2.0f * ::smooth_ramp(-params.lattice_softness,
                                     params.lattice_softness, distance);
}

HS_FLASH_INLINE inline float noise_contour(const FastNoiseLite &noise,
                                           ::NoiseBasis basis,
                                           const Vector &coordinate,
                                           float contrast) {
  const float sample =
      hs::clamp(sample_noise_octaves(noise, basis, coordinate), -1.0f, 1.0f);
  return sample * (1.0f + contrast) / (1.0f + contrast * fabsf(sample));
}

template <typename State> struct TwinWave : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static float sample(const SourceInput &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return twin_wave(
        stereo_pattern_args(input.warped.coords, params.pattern_freq),
        prepared);
  }
};

template <typename State> struct Rings : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static float sample(const SourceInput &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return rings(stereo_pattern_args(input.warped.coords, params.pattern_freq),
                 prepared);
  }
};

template <typename State> struct Spiral : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).angle } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static float sample(const SourceInput &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return spiral(stereo_pattern_args(input.warped.coords, params.pattern_freq),
                  prepared);
  }
};

template <typename State> struct Grid : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::params(frame).pattern_mix } -> std::convertible_to<float>;
        { State::params(frame).complexity } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
        { State::prepare(frame).secondary } -> std::convertible_to<float>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static float sample(const SourceInput &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return grid(stereo_pattern_args(input.warped.coords, params.pattern_freq),
                params, prepared);
  }
};

template <typename State> struct PrimitiveLattice : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        {
          State::params(frame).lattice_cell_scale
        } -> std::convertible_to<float>;
        {
          State::params(frame).lattice_shape_blend
        } -> std::convertible_to<float>;
        { State::params(frame).lattice_softness } -> std::convertible_to<float>;
        { State::params(frame).lattice_radius } -> std::convertible_to<float>;
      };

  __attribute__((always_inline)) static float sample(const SourceInput &input,
                                                     const FrameState &frame) {
    return primitive_lattice(input.warped.coords, State::params(frame));
  }
};

template <typename State, ::NoiseBasis BasisV>
struct ProjectedNoise : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::noise_scale(frame) } -> std::same_as<float>;
        { State::noise_time(frame) } -> std::same_as<float>;
        { State::noise_contrast(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static float sample(const SourceInput &input,
                                                     const FrameState &frame) {
    return noise_contour(State::noise(frame), BasisV,
                         noise_projected_coordinate(input.warped.coords,
                                                    State::noise_scale(frame),
                                                    State::noise_time(frame)),
                         State::noise_contrast(frame));
  }
};

template <typename State, ::NoiseBasis BasisV>
struct SphericalNoise : ExactPolicy {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ProjectedNoise<State, BasisV>::template PROVIDER_VALID<CandidateBinding>;

  __attribute__((always_inline)) static float sample(const SourceInput &input,
                                                     const FrameState &frame) {
    return noise_contour(State::noise(frame), BasisV,
                         noise_sphere_coordinate(input.projected.sphere,
                                                 State::noise_scale(frame),
                                                 State::noise_time(frame)),
                         State::noise_contrast(frame));
  }
};

} // namespace Source

} // namespace Pullback
