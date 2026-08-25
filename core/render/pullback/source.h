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
      Field<GridSourceParams>{"pattern-freq", &GridSourceParams::pattern_freq,
                              "Pattern Freq", 0.1f, 64.0f, FieldCurve::LERP},
      Field<GridSourceParams>{"speed", &GridSourceParams::speed, "Speed", 0.0f,
                              0.5f, FieldCurve::LERP},
      Field<GridSourceParams>{"complexity", &GridSourceParams::complexity,
                              "Complexity", 0.0f, 3.0f, FieldCurve::LERP},
      Field<GridSourceParams>{"pattern-mix", &GridSourceParams::pattern_mix,
                              "Pattern Mix", 0.0f, 1.0f, FieldCurve::LERP},
      Field<GridSourceParams>{"drift", &GridSourceParams::secondary_rate,
                              "Drift", 0.0f, 1.25f, FieldCurve::LERP},
      Field<GridSourceParams>{"angle-speed", &GridSourceParams::angle_rate,
                              "Source Angle Speed", 0.0f, 0.05f,
                              FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<GridSourceParams>());

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
      Field<TwinWaveSourceParams>{
          "pattern-freq", &TwinWaveSourceParams::pattern_freq, "Pattern Freq",
          0.1f, 20.0f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{"speed", &TwinWaveSourceParams::speed,
                                  "Speed", 0.0f, 0.5f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{"drift",
                                  &TwinWaveSourceParams::secondary_rate,
                                  "Drift", 0.0f, 1.25f, FieldCurve::LERP},
      Field<TwinWaveSourceParams>{
          "angle-speed", &TwinWaveSourceParams::angle_rate,
          "Source Angle Speed", 0.0f, 0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<TwinWaveSourceParams>());

/**
 * @brief Source parameters for the rotating spiral field
 *        (Pullback::Source::Spiral).
 */
struct SpiralSourceParams {
  float pattern_freq = 1.0f; /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;        /**< Per-frame advance of the primary phase. */
  float angle_rate = 0.0f;   /**< Per-frame advance of the spiral rotation. */

  static constexpr auto FIELDS = std::array{
      Field<SpiralSourceParams>{"pattern-freq",
                                &SpiralSourceParams::pattern_freq,
                                "Pattern Freq", 0.1f, 20.0f, FieldCurve::LERP},
      Field<SpiralSourceParams>{"speed", &SpiralSourceParams::speed, "Speed",
                                0.0f, 0.5f, FieldCurve::LERP},
      Field<SpiralSourceParams>{"angle-speed", &SpiralSourceParams::angle_rate,
                                "Source Angle Speed", 0.0f, 0.05f,
                                FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<SpiralSourceParams>());

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
      Field<NoiseSourceParams>{"noise-scale", &NoiseSourceParams::noise_scale,
                               "Source Noise Scale", 1.0f / 64.0f, 64.0f,
                               FieldCurve::LERP},
      Field<NoiseSourceParams>{
          "noise-contrast", &NoiseSourceParams::noise_contrast,
          "Source Noise Contrast", 0.0f, 8.0f, FieldCurve::LERP},
      Field<NoiseSourceParams>{
          "noise-speed", &NoiseSourceParams::noise_time_rate,
          "Source Noise Speed", -1.0f / 64.0f, 1.0f / 64.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<NoiseSourceParams>());

/**
 * @brief Source family selecting the plane-domain noise contour
 *        (Pullback::Source::ProjectedNoise): the field is sampled in the
 *        projected chart, so it carries the projection's distortion.
 */
struct ProjectedNoiseSourceParams : NoiseSourceParams {};

/**
 * @brief Source family selecting the sphere-domain noise contour
 *        (Pullback::Source::SphericalNoise): the field is sampled on the
 *        pre-projection direction, so it is seamless and moves only with
 *        the projection frame.
 */
struct SphericalNoiseSourceParams : NoiseSourceParams {};

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
      Field<LatticeSourceParams>{
          "lattice-cell-scale", &LatticeSourceParams::lattice_cell_scale,
          "Lattice Cell Scale", 1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<LatticeSourceParams>{"lattice-shape",
                                 &LatticeSourceParams::lattice_shape_blend,
                                 "Lattice Shape", 0.0f, 1.0f, FieldCurve::LERP},
      Field<LatticeSourceParams>{
          "lattice-softness", &LatticeSourceParams::lattice_softness,
          "Lattice Softness", 1.0f / 1024.0f, 1.0f, FieldCurve::LOG_POSITIVE},
      Field<LatticeSourceParams>{
          "lattice-radius", &LatticeSourceParams::lattice_radius,
          "Lattice Radius", 1.0f / 64.0f, 0.49f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<LatticeSourceParams>());

/** @brief Source parameters for latitude bands on a moving sphere. */
struct SphericalRingsSourceParams {
  float ring_count = 6.0f;      /**< Number of bands from pole to pole. */
  float ring_thickness = 0.08f; /**< Band half-width, in radians. */
  float ring_softness = 0.02f;  /**< Angular width of the antialiased edge. */
  float speed = 0.0f;           /**< Per-frame advance of the band phase. */
  float spin_rate = 0.0f;       /**< Per-frame rotation of the band axis. */
  float wander = 0.0f;          /**< Fraction of the random walk applied. */

  static constexpr auto FIELDS = std::array{
      Field<SphericalRingsSourceParams>{
          "ring-count", &SphericalRingsSourceParams::ring_count, "Ring Count",
          1.0f, 32.0f, FieldCurve::SNAP},
      Field<SphericalRingsSourceParams>{
          "ring-thickness", &SphericalRingsSourceParams::ring_thickness,
          "Ring Thickness", 1.0f / 512.0f, 0.5f, FieldCurve::LOG_POSITIVE},
      Field<SphericalRingsSourceParams>{
          "ring-softness", &SphericalRingsSourceParams::ring_softness,
          "Ring Softness", 1.0f / 1024.0f, 0.25f, FieldCurve::LOG_POSITIVE},
      Field<SphericalRingsSourceParams>{
          "speed", &SphericalRingsSourceParams::speed, "Ring Speed", -0.5f,
          0.5f, FieldCurve::LERP},
      Field<SphericalRingsSourceParams>{
          "spin-speed", &SphericalRingsSourceParams::spin_rate,
          "Ring Spin Speed", -0.05f, 0.05f, FieldCurve::LERP},
      Field<SphericalRingsSourceParams>{
          "wander", &SphericalRingsSourceParams::wander, "Ring Wander", 0.0f,
          1.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<SphericalRingsSourceParams>());

/** @brief Source parameters for the quadratic escape-time fractal. */
struct FractalSourceParams {
  float scale = 0.5f;      /**< Plane-coordinate scale before iteration. */
  float iterations = 8.0f; /**< Escape iterations in [2,16]. */
  float julia_mix = 0.0f;  /**< Blend from Mandelbrot at 0 to Julia at 1. */
  float julia_re = -0.8f;  /**< Julia seed real component. */
  float julia_im = 0.156f; /**< Julia seed imaginary component. */
  float contours = 4.0f;   /**< Exterior contour cycles across the orbit. */
  float speed = 0.0f;      /**< Per-frame rotation of the Julia seed. */
  float angle_rate = 0.0f; /**< Per-frame rotation of the source plane. */

  static constexpr auto FIELDS = std::array{
      Field<FractalSourceParams>{"fractal-scale", &FractalSourceParams::scale,
                                 "Fractal Scale", 1.0f / 64.0f, 8.0f,
                                 FieldCurve::LOG_POSITIVE},
      Field<FractalSourceParams>{
          "fractal-iterations", &FractalSourceParams::iterations,
          "Fractal Iterations", 2.0f, 16.0f, FieldCurve::SNAP},
      Field<FractalSourceParams>{"julia-mix", &FractalSourceParams::julia_mix,
                                 "Julia Mix", 0.0f, 1.0f, FieldCurve::LERP},
      Field<FractalSourceParams>{"julia-real", &FractalSourceParams::julia_re,
                                 "Julia Real", -1.5f, 1.5f, FieldCurve::LERP},
      Field<FractalSourceParams>{
          "julia-imaginary", &FractalSourceParams::julia_im, "Julia Imaginary",
          -1.5f, 1.5f, FieldCurve::LERP},
      Field<FractalSourceParams>{
          "fractal-contours", &FractalSourceParams::contours,
          "Fractal Contours", 0.0f, 16.0f, FieldCurve::LERP},
      Field<FractalSourceParams>{"speed", &FractalSourceParams::speed,
                                 "Fractal Speed", -0.05f, 0.05f,
                                 FieldCurve::LERP},
      Field<FractalSourceParams>{
          "angle-speed", &FractalSourceParams::angle_rate, "Fractal Spin Speed",
          -0.05f, 0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<FractalSourceParams>());

enum class TessellationKind : uint8_t {
  TRIANGULAR = 0,
  SQUARE = 1,
  HEXAGONAL = 2
};

/** @brief Source parameters for periodic polygon edge tessellations. */
struct TessellationSourceParams {
  float cell_scale = 1.0f;      /**< Tessellation cells per plane unit. */
  float line_thickness = 0.04f; /**< Edge half-width in cell units. */
  float line_softness = 0.02f;  /**< Width of the antialiased edge. */
  float angle_rate = 0.0f;      /**< Per-frame rotation of the tessellation. */

  static constexpr auto FIELDS = std::array{
      Field<TessellationSourceParams>{
          "cell-scale", &TessellationSourceParams::cell_scale, "Cell Scale",
          1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<TessellationSourceParams>{
          "line-thickness", &TessellationSourceParams::line_thickness,
          "Line Thickness", 1.0f / 1024.0f, 0.25f, FieldCurve::LOG_POSITIVE},
      Field<TessellationSourceParams>{
          "line-softness", &TessellationSourceParams::line_softness,
          "Line Softness", 1.0f / 1024.0f, 0.25f, FieldCurve::LOG_POSITIVE},
      Field<TessellationSourceParams>{
          "angle-speed", &TessellationSourceParams::angle_rate,
          "Tessellation Spin Speed", -0.05f, 0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<TessellationSourceParams>());

/** @brief The source stage's phases, resolved once per frame. */
struct PreparedSource {
  float primary;   /**< Primary phase, wrapped into [0,2pi). */
  float secondary; /**< Secondary phase, wrapped into [0,2pi). */
  float angle;     /**< Source rotation, in radians. */
  float angle_cos; /**< Cosine of `angle`. */
  float angle_sin; /**< Sine of `angle`. */
};

/** @brief Per-frame axis and phase of the spherical ring source. */
struct PreparedSphericalRings {
  Vector axis; /**< Unit normal of the rings' equatorial plane. */
  float phase; /**< Angular band offset, in radians. */
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

template <typename Params>
HS_FLASH_MEMBER inline float
spherical_rings(const Vector &input, const Params &params,
                const PreparedSphericalRings &prepared) {
  const float axis_height = hs::clamp(dot(input, prepared.axis), -1.0f, 1.0f);
  const float latitude = fast_atan2(
      axis_height, sqrtf(std::max(0.0f, 1.0f - axis_height * axis_height)));
  const float count = std::max(params.ring_count, 1.0f);
  const float cycle =
      wrap_t((count * latitude - prepared.phase) / PI_F + 0.5f) - 0.5f;
  const float distance = fabsf(cycle) * PI_F / count;
  const float edge =
      ::smooth_ramp(params.ring_thickness,
                    params.ring_thickness + params.ring_softness, distance);
  return 1.0f - 2.0f * edge;
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
  float re = x + prepared.primary;
  float im = y - prepared.secondary;
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

template <typename Params, typename Prepared>
HS_FLASH_MEMBER inline float escape_fractal(const Complex &input,
                                            const Params &params,
                                            const Prepared &prepared) {
  const float x = params.scale * (input.re * prepared.angle_cos +
                                  input.im * prepared.angle_sin);
  const float y = params.scale * (-input.re * prepared.angle_sin +
                                  input.im * prepared.angle_cos);
  const float seed_cos = fast_cosf(prepared.primary);
  const float seed_sin = fast_sinf(prepared.primary);
  const float seed_re = params.julia_re * seed_cos - params.julia_im * seed_sin;
  const float seed_im = params.julia_re * seed_sin + params.julia_im * seed_cos;
  const float mix = hs::clamp(params.julia_mix, 0.0f, 1.0f);
  float z_re = x * mix;
  float z_im = y * mix;
  const float c_re = hs::lerp(x, seed_re, mix);
  const float c_im = hs::lerp(y, seed_im, mix);
  const int iterations =
      static_cast<int>(hs::clamp(params.iterations, 2.0f, 16.0f));
  for (int iteration = 0; iteration < iterations; ++iteration) {
    const float next_re = z_re * z_re - z_im * z_im + c_re;
    z_im = 2.0f * z_re * z_im + c_im;
    z_re = next_re;
    const float magnitude_squared = z_re * z_re + z_im * z_im;
    if (magnitude_squared > 4.0f) {
      const float orbit =
          (static_cast<float>(iteration) +
           hs::clamp((magnitude_squared - 4.0f) / 12.0f, 0.0f, 1.0f)) /
          static_cast<float>(iterations);
      return fast_cosf(TWO_PI_F * params.contours * orbit);
    }
  }
  return 1.0f;
}

HS_FLASH_INLINE inline float distance_to_lattice_line(float coordinate) {
  return fabsf(wrap_t(coordinate + 0.5f) - 0.5f);
}

template <typename Params, typename Prepared>
HS_FLASH_MEMBER inline float
tessellation(const Complex &input, const Params &params, TessellationKind kind,
             const Prepared &prepared) {
  constexpr float SQRT_3 = 1.7320508075688772f;
  const float x = params.cell_scale * (input.re * prepared.angle_cos +
                                       input.im * prepared.angle_sin);
  const float y = params.cell_scale * (-input.re * prepared.angle_sin +
                                       input.im * prepared.angle_cos);
  float distance;
  switch (kind) {
  case TessellationKind::TRIANGULAR:
    distance = std::min(
        distance_to_lattice_line(x),
        std::min(distance_to_lattice_line(0.5f * x + 0.5f * SQRT_3 * y),
                 distance_to_lattice_line(-0.5f * x + 0.5f * SQRT_3 * y)));
    break;
  case TessellationKind::SQUARE: {
    const float cell_x = wrap_t(x + 0.5f) - 0.5f;
    const float cell_y = wrap_t(y + 0.5f) - 0.5f;
    distance = 0.5f - std::max(fabsf(cell_x), fabsf(cell_y));
    break;
  }
  case TessellationKind::HEXAGONAL:
  default: {
    const float axial_x = (2.0f / 3.0f) * x;
    const float axial_z = y / SQRT_3 - 0.5f * axial_x;
    const float axial_y = -axial_x - axial_z;
    float cell_x = roundf(axial_x);
    float cell_y = roundf(axial_y);
    float cell_z = roundf(axial_z);
    const float error_x = fabsf(cell_x - axial_x);
    const float error_y = fabsf(cell_y - axial_y);
    const float error_z = fabsf(cell_z - axial_z);
    if (error_x > error_y && error_x > error_z)
      cell_x = -cell_y - cell_z;
    else if (error_y > error_z)
      cell_y = -cell_x - cell_z;
    else
      cell_z = -cell_x - cell_y;
    const float local_x = x - 1.5f * cell_x;
    const float local_y = y - SQRT_3 * (cell_z + 0.5f * cell_x);
    distance =
        1.0f -
        std::max(fabsf(local_x),
                 std::max(fabsf(0.5f * local_x + 0.5f * SQRT_3 * local_y),
                          fabsf(-0.5f * local_x + 0.5f * SQRT_3 * local_y)));
    break;
  }
  }
  const float edge =
      ::smooth_ramp(params.line_thickness,
                    params.line_thickness + params.line_softness, distance);
  return 1.0f - 2.0f * edge;
}

HS_FLASH_INLINE inline float noise_contour(const FastNoiseLite &noise,
                                           ::NoiseBasis basis,
                                           const Vector &coordinate,
                                           float contrast) {
  const float sample =
      hs::clamp(sample_noise_octaves(noise, basis, coordinate), -1.0f, 1.0f);
  return sample * (1.0f + contrast) / (1.0f + contrast * fabsf(sample));
}

template <typename State> struct TwinWave : ApproximationDefaults {
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

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return twin_wave(stereo_pattern_args(input.coords, params.pattern_freq),
                     prepared);
  }
};

template <typename State> struct Rings : ApproximationDefaults {
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

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return rings(stereo_pattern_args(input.coords, params.pattern_freq),
                 prepared);
  }
};

template <typename State> struct SphericalRings : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).ring_count } -> std::convertible_to<float>;
        { State::params(frame).ring_thickness } -> std::convertible_to<float>;
        { State::params(frame).ring_softness } -> std::convertible_to<float>;
        { State::prepare(frame).axis } -> std::convertible_to<Vector>;
        { State::prepare(frame).phase } -> std::convertible_to<float>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static float sample(const SphereSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return spherical_rings(input.dir, State::params(frame), prepared);
  }
};

template <typename State> struct Spiral : ApproximationDefaults {
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

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return spiral(stereo_pattern_args(input.coords, params.pattern_freq),
                  prepared);
  }
};

template <typename State> struct Grid : ApproximationDefaults {
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

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return grid(stereo_pattern_args(input.coords, params.pattern_freq), params,
                prepared);
  }
};

template <typename State> struct PrimitiveLattice : ApproximationDefaults {
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

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame) {
    return primitive_lattice(input.coords, State::params(frame));
  }
};

template <typename State> struct EscapeFractal : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).scale } -> std::convertible_to<float>;
        { State::params(frame).iterations } -> std::convertible_to<float>;
        { State::params(frame).julia_mix } -> std::convertible_to<float>;
        { State::params(frame).julia_re } -> std::convertible_to<float>;
        { State::params(frame).julia_im } -> std::convertible_to<float>;
        { State::params(frame).contours } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return escape_fractal(input.coords, State::params(frame), prepared);
  }
};

template <typename State, TessellationKind KindV>
struct Tessellation : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).cell_scale } -> std::convertible_to<float>;
        { State::params(frame).line_thickness } -> std::convertible_to<float>;
        { State::params(frame).line_softness } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return tessellation(input.coords, State::params(frame), KindV, prepared);
  }
};

template <typename State, ::NoiseBasis BasisV>
struct ProjectedNoise : ApproximationDefaults {
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

  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame) {
    return noise_contour(State::noise(frame), BasisV,
                         noise_projected_coordinate(input.coords,
                                                    State::noise_scale(frame),
                                                    State::noise_time(frame)),
                         State::noise_contrast(frame));
  }
};

template <typename State, ::NoiseBasis BasisV>
struct SphericalNoise : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ProjectedNoise<State, BasisV>::template PROVIDER_VALID<CandidateBinding>;

  /** @brief Post-projection form: samples the plane carrier's retained
      pre-projection point. */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame) {
    return noise_contour(State::noise(frame), BasisV,
                         noise_sphere_coordinate(input.sphere,
                                                 State::noise_scale(frame),
                                                 State::noise_time(frame)),
                         State::noise_contrast(frame));
  }
  __attribute__((always_inline)) static float sample(const SphereSample &input,
                                                     const FrameState &frame) {
    return noise_contour(State::noise(frame), BasisV,
                         noise_sphere_coordinate(input.dir,
                                                 State::noise_scale(frame),
                                                 State::noise_time(frame)),
                         State::noise_contrast(frame));
  }
};

} // namespace Source

} // namespace Pullback
