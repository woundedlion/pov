/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>

#include "math/noise_field.h"

#include "math/projection_patterns.h"
#include "render/pullback/contract.h"
#include "render/pullback/fields.h"
#include "math/3dmath.h"

/**
 * @file source.h
 * @brief Scalar source-field policies.
 */

namespace Pullback {

namespace Source {

/** @brief Parameters of the concentric ring source. */
struct RingsSourceParams {
  float pattern_freq = 1.0f; /**< Radial pattern frequency. */
  float speed = 0.0f;        /**< Phase advance per frame, in radians. */
  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<RingsSourceParams>{"pattern-freq", &RingsSourceParams::pattern_freq,
                               "Pattern Freq", 0.1f, 20.0f, FieldCurve::LERP},
      Field<RingsSourceParams>{"speed", &RingsSourceParams::speed, "Speed",
                               0.0f, 0.5f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<RingsSourceParams>());
static_assert(field_defaults_in_range<RingsSourceParams>());

/** @brief Advances the phase clocks declared by one scalar source family.
    @tparam Params Parameter block; `speed`, `secondary_rate` and
            `angle_rate` each drive a clock only when declared.
    @param params The family's parameters.
    @param primary Primary phase, wrapped into (-2pi, 2pi).
    @param secondary Secondary phase, wrapped into (-2pi, 2pi).
    @param angle Source rotation, wrapped into (-2pi, 2pi). */
template <typename Params>
inline void advance_clocks(const Params &params, float &primary,
                           float &secondary, float &angle) {
  if constexpr (requires { params.speed; })
    primary = fmodf(primary + params.speed, math::TWO_PI_F);
  if constexpr (requires { params.secondary_rate; })
    secondary =
        fmodf(secondary + params.speed * params.secondary_rate, math::TWO_PI_F);
  if constexpr (requires { params.angle_rate; })
    angle = fmodf(angle + params.angle_rate, math::TWO_PI_F);
}

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

  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<GridSourceParams>{"pattern-freq", &GridSourceParams::pattern_freq,
                              "Pattern Freq", 0.01f, 64.0f, FieldCurve::LERP},
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
static_assert(field_defaults_in_range<GridSourceParams>());

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

  /** @brief Per-parameter id, member, label, range and curve. */
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
static_assert(field_defaults_in_range<TwinWaveSourceParams>());

/**
 * @brief Source parameters for the rotating spiral field
 *        (Pullback::Source::Spiral).
 */
struct SpiralSourceParams {
  float pattern_freq = 1.0f; /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;        /**< Per-frame advance of the primary phase. */
  float angle_rate = 0.0f;   /**< Per-frame advance of the spiral rotation. */

  /** @brief Per-parameter id, member, label, range and curve. */
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
static_assert(field_defaults_in_range<SpiralSourceParams>());

/**
 * @brief Source parameters for the noise-contour sources
 *        (Pullback::Source::ProjectedNoise and SphericalNoise).
 */
struct NoiseSourceParams {
  float noise_scale = 1.0f;    /**< Spatial scale of the sampled field. */
  float noise_contrast = 0.0f; /**< Contour sharpening applied to the sample. */
  float noise_time_rate = 0.0f; /**< Per-frame advance of the noise time
                                     coordinate. */

  /** @brief Per-parameter id, member, label, range and curve. */
  static constexpr auto FIELDS = std::array{
      Field<NoiseSourceParams>{"noise-scale", &NoiseSourceParams::noise_scale,
                               "Source Noise Scale", 1.0f / 64.0f, 64.0f,
                               FieldCurve::LOG_POSITIVE},
      Field<NoiseSourceParams>{
          "noise-contrast", &NoiseSourceParams::noise_contrast,
          "Source Noise Contrast", 0.0f, 8.0f, FieldCurve::LERP},
      Field<NoiseSourceParams>{
          "noise-speed", &NoiseSourceParams::noise_time_rate,
          "Source Noise Speed", -1.0f / 64.0f, 1.0f / 64.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<NoiseSourceParams>());
static_assert(field_defaults_in_range<NoiseSourceParams>());

/**
 * @brief Source family selecting the plane-domain noise contour
 *        (Pullback::Source::ProjectedNoise): the field is sampled in the
 *        projected chart, so it carries the projection's distortion.
 */
struct ProjectedNoiseSourceParams : NoiseSourceParams {};

/**
 * @brief Source family selecting the sphere-domain noise contour
 *        (Pullback::Source::SphericalNoise): the field is sampled on the
 *        pre-projection direction, so it is seamless and does not carry
 *        the planar projection's distortion.
 */
struct SphericalNoiseSourceParams : NoiseSourceParams {};

/**
 * @brief Source parameters for the per-cell primitive lattice
 *        (Pullback::Source::PrimitiveLattice).
 */
struct LatticeSourceParams {
  float lattice_cell_scale = 1.0f;  /**< Lattice cells per plane unit. */
  float lattice_shape_blend = 0.0f; /**< Cell primitive, from a circle at 0 to a
                                         square at 1. */
  float lattice_softness = 0.05f;   /**< Half-width of the ramp across the
                                         primitive's boundary. */
  float lattice_radius = 0.25f;     /**< Primitive radius in cell units. */

  /** @brief Per-parameter id, member, label, range and curve. */
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
static_assert(field_defaults_in_range<LatticeSourceParams>());

/** @brief Source parameters for latitude bands on a moving sphere. */
struct SphericalRingsSourceParams {
  float ring_count = 6.0f;      /**< Number of bands from pole to pole. */
  float ring_thickness = 0.08f; /**< Band half-width, in radians. */
  float ring_softness = 0.02f;  /**< Angular width of the antialiased edge. */
  float speed = 0.0f;           /**< Per-frame advance of the band phase. */
  float spin_rate = 0.0f;       /**< Per-frame rotation of the band axis. */
  float wander = 0.0f;          /**< Fraction of the random walk applied. */

  /** @brief Per-parameter id, member, label, range and curve. */
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
static_assert(field_defaults_in_range<SphericalRingsSourceParams>());

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

  /** @brief Per-parameter id, member, label, range and curve. */
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
static_assert(field_defaults_in_range<FractalSourceParams>());

/** @brief Cell shape of a `Tessellation` source. */
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

  /** @brief Per-parameter id, member, label, range and curve. */
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
static_assert(field_defaults_in_range<TessellationSourceParams>());

/** @brief The source stage's phases, resolved once per frame. */
struct PreparedSource {
  float primary;   /**< Primary phase, wrapped into (-2pi,2pi). */
  float secondary; /**< Secondary phase, wrapped into (-2pi,2pi). */
  float angle;     /**< Source rotation, wrapped into (-2pi,2pi) radians. */
  float angle_cos; /**< Cosine of `angle`. */
  float angle_sin; /**< Sine of `angle`. */
};

/** @brief Frame-constant quadratic-fractal seed and iteration controls. */
struct PreparedFractal : PreparedSource {
  float seed_re;  ///< Julia seed, real part, rotated by `primary`.
  float seed_im;  ///< Julia seed, imaginary part, rotated by `primary`.
  float mix;      ///< Julia mix clamped to [0, 1].
  int iterations; ///< Escape iterations clamped to [2, 16].
};

/**
 * @brief Resolves the fractal's frame-constant seed and iteration controls.
 * @tparam Params FractalSourceParams-shaped parameter block.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param params Fractal parameters.
 * @param source This frame's phases; `primary` rotates the Julia seed.
 * @return @p source extended with the seed, mix and iteration count.
 */
template <typename Params, typename Prepared>
HS_FLASH_INLINE inline PreparedFractal prepare_fractal(const Params &params,
                                                       const Prepared &source) {
  const float seed_cos = math::fast_cosf(source.primary);
  const float seed_sin = math::fast_sinf(source.primary);
  return {{source.primary, source.secondary, source.angle, source.angle_cos,
           source.angle_sin},
          params.julia_re * seed_cos - params.julia_im * seed_sin,
          params.julia_re * seed_sin + params.julia_im * seed_cos,
          hs::clamp(params.julia_mix, 0.0f, 1.0f),
          static_cast<int>(hs::clamp(params.iterations, 2.0f, 16.0f))};
}

/** @brief Per-frame axis and phase of the spherical ring source. */
struct PreparedSphericalRings {
  math::Vector axis; /**< Unit normal of the rings' equatorial plane. */
  float phase;       /**< Angular band offset, in radians. */
};

/** @brief Wraps this frame's source phases with the rotation's cosine pair.
    @param primary Primary phase.
    @param secondary Secondary phase.
    @param angle Source rotation, in radians.
    @return The phases with the cosine and sine of @p angle. */
HS_FLASH_INLINE inline PreparedSource prepare(float primary, float secondary,
                                              float angle) {
  return {primary, secondary, angle, cosf(angle), sinf(angle)};
}

template <typename State, typename Binding>
concept StateProvider = Detail::ParamsProvider<State, Binding> &&
                        requires(const typename Binding::FrameState &frame) {
                          State::prepare(frame);
                        };

/**
 * @brief Two-wave interference: waves along x and along
 *        the axis at the source angle.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param input Pattern-space coordinate.
 * @param prepared Per-frame phases and rotation.
 * @return Field value in [-1, 1].
 */
template <typename Prepared>
HS_HOT_FLASH_MEMBER inline float twin_wave(const math::Complex &input,
                                           const Prepared &prepared) {
  const float rotated =
      input.re * prepared.angle_cos + input.im * prepared.angle_sin;
  return 0.5f * (math::fast_sinf(input.re + prepared.primary) +
                 math::fast_sinf(rotated + prepared.secondary));
}

/**
 * @brief Concentric rings about the origin, moving outward with `primary`.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param input Pattern-space coordinate.
 * @param prepared Per-frame phases.
 * @return Field value in [-1, 1].
 */
template <typename Prepared>
HS_HOT_FLASH_MEMBER inline float rings(const math::Complex &input,
                                       const Prepared &prepared) {
  return math::fast_sinf(input.magnitude() - prepared.primary);
}

/**
 * @brief Latitude bands about the prepared axis.
 * @tparam Params SphericalRingsSourceParams-shaped parameter block.
 * @param input Unit direction.
 * @param params Band count, thickness and softness.
 * @param prepared Band axis and phase.
 * @return 1 inside a band, -1 outside, ramped across the edge.
 */
template <typename Params>
HS_HOT_FLASH_MEMBER inline float
spherical_rings(const math::Vector &input, const Params &params,
                const PreparedSphericalRings &prepared) {
  const float axis_height =
      hs::clamp(math::dot(input, prepared.axis), -1.0f, 1.0f);
  const float latitude = math::fast_atan2(
      axis_height, sqrtf(fmaxf(0.0f, 1.0f - axis_height * axis_height)));
  const float count = fmaxf(params.ring_count, 1.0f);
  const float cycle =
      math::wrap_t((count * latitude - prepared.phase) / math::PI_F + 0.5f) -
      0.5f;
  const float distance = fabsf(cycle) * math::PI_F / count;
  const float edge = ::math::smooth_ramp(
      params.ring_thickness, params.ring_thickness + params.ring_softness,
      distance);
  return 1.0f - 2.0f * edge;
}

/**
 * @brief Three-armed Archimedean spiral turned by the source angle.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param input Pattern-space coordinate.
 * @param prepared Per-frame phase and rotation.
 * @return Field value in [-1, 1].
 */
template <typename Prepared>
HS_HOT_FLASH_MEMBER inline float spiral(const math::Complex &input,
                                        const Prepared &prepared) {
  const float radius = input.magnitude();
  const float azimuth = math::fast_atan2(input.im, input.re);
  return math::fast_sinf(radius - 3.0f * (azimuth + prepared.angle) -
                         prepared.primary);
}

/**
 * @brief Coupled sine grid in the rotated frame.
 * @tparam Params Parameter block of the source family.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param input Pattern-space coordinate.
 * @param params Complexity and pattern mix.
 * @param prepared Per-frame phases and rotation.
 * @return Field value in [-1, 1].
 */
template <typename Params, typename Prepared>
HS_HOT_FLASH_MEMBER inline float grid(const math::Complex &input,
                                      const Params &params,
                                      const Prepared &prepared) {
  const float x = input.re * prepared.angle_cos + input.im * prepared.angle_sin;
  const float y =
      -input.re * prepared.angle_sin + input.im * prepared.angle_cos;
  if (params.pattern_mix == 1.0f)
    return math::fast_sinf(x + prepared.primary) *
           math::fast_cosf(y - prepared.secondary);
  float re = x + prepared.primary;
  float im = y - prepared.secondary;
  if (params.complexity != 0.0f) {
    re += params.complexity * math::fast_sinf(y + prepared.primary);
    im += params.complexity * math::fast_cosf(x - prepared.secondary);
  }
  const float coupled = math::fast_sinf(re) * math::fast_cosf(im);
  if (params.pattern_mix == 0.0f)
    return coupled;
  const float direct = math::fast_sinf(x + prepared.primary) *
                       math::fast_cosf(y - prepared.secondary);
  return hs::lerp(coupled, direct, params.pattern_mix);
}

/**
 * @brief Repeating cell primitive blended between circle and square.
 * @tparam Params LatticeSourceParams-shaped parameter block.
 * @param input Plane coordinate.
 * @param params Cell scale, shape blend, radius and softness.
 * @return 1 inside the primitive, -1 outside, ramped across the edge.
 */
template <typename Params>
HS_HOT_FLASH_MEMBER inline float primitive_lattice(const math::Complex &input,
                                                   const Params &params) {
  const float x =
      math::wrap_t(params.lattice_cell_scale * input.re + 0.5f) - 0.5f;
  const float y =
      math::wrap_t(params.lattice_cell_scale * input.im + 0.5f) - 0.5f;
  const float circle = sqrtf(x * x + y * y) - params.lattice_radius;
  const float bx = fabsf(x) - params.lattice_radius;
  const float by = fabsf(y) - params.lattice_radius;
  const float square = sqrtf(fmaxf(bx, 0.0f) * fmaxf(bx, 0.0f) +
                             fmaxf(by, 0.0f) * fmaxf(by, 0.0f)) +
                       fminf(fmaxf(bx, by), 0.0f);
  const float distance = hs::lerp(circle, square, params.lattice_shape_blend);
  return 1.0f - 2.0f * ::math::smooth_ramp(-params.lattice_softness,
                                           params.lattice_softness, distance);
}

/**
 * @brief Quadratic escape-time fractal, blended from Mandelbrot to Julia.
 * @tparam Params FractalSourceParams-shaped parameter block.
 * @tparam Prepared PreparedFractal-shaped per-frame state.
 * @param input Plane coordinate.
 * @param params Scale and contour count.
 * @param prepared Rotation, seed, mix and iteration count.
 * @return Cosine contour of the escape orbit; 1 for points that never escape.
 */
template <typename Params, typename Prepared>
HS_HOT_FLASH_MEMBER inline float escape_fractal(const math::Complex &input,
                                                const Params &params,
                                                const Prepared &prepared) {
  const float x = params.scale * (input.re * prepared.angle_cos +
                                  input.im * prepared.angle_sin);
  const float y = params.scale * (-input.re * prepared.angle_sin +
                                  input.im * prepared.angle_cos);
  float z_re = x * prepared.mix;
  float z_im = y * prepared.mix;
  const float c_re = hs::lerp(x, prepared.seed_re, prepared.mix);
  const float c_im = hs::lerp(y, prepared.seed_im, prepared.mix);
  const int iterations = prepared.iterations;
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
      return math::fast_cosf(math::TWO_PI_F * params.contours * orbit);
    }
  }
  return 1.0f;
}

/**
 * @brief Distance to the nearest integer.
 * @param coordinate Lattice coordinate.
 * @return Distance in [0, 0.5].
 */
HS_O3_FN inline float distance_to_lattice_line(float coordinate) {
  return fabsf(math::wrap_t(coordinate + 0.5f) - 0.5f);
}

/** @brief Distance from the scaled, rotated plane point to the nearest cell
    edge of @p kind.
    @param x Cell-space x.
    @param y Cell-space y.
    @param kind Cell shape.
    @return Nonnegative distance in cell units. */
__attribute__((always_inline)) inline float
tessellation_distance(float x, float y, TessellationKind kind) {
  constexpr float SQRT_3 = 1.7320508075688772f;
  switch (kind) {
  case TessellationKind::TRIANGULAR:
    return fminf(
        distance_to_lattice_line(x),
        fminf(distance_to_lattice_line(0.5f * x + 0.5f * SQRT_3 * y),
              distance_to_lattice_line(-0.5f * x + 0.5f * SQRT_3 * y)));
  case TessellationKind::SQUARE: {
    const float cell_x = math::wrap_t(x + 0.5f) - 0.5f;
    const float cell_y = math::wrap_t(y + 0.5f) - 0.5f;
    return 0.5f - fmaxf(fabsf(cell_x), fabsf(cell_y));
  }
  case TessellationKind::HEXAGONAL:
    break;
  }
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
  return 0.5f * SQRT_3 -
         fmaxf(fabsf(local_y),
               fmaxf(fabsf(0.5f * SQRT_3 * local_x + 0.5f * local_y),
                     fabsf(0.5f * SQRT_3 * local_x - 0.5f * local_y)));
}

/**
 * @brief Antialiased edge lines of a periodic polygon tessellation.
 * @tparam Params TessellationSourceParams-shaped parameter block.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param input Plane coordinate.
 * @param params Cell scale, line thickness and softness.
 * @param kind Cell shape.
 * @param prepared Per-frame rotation.
 * @return 1 on an edge, -1 inside a cell, ramped between.
 */
template <typename Params, typename Prepared>
HS_HOT_FLASH_MEMBER inline float
tessellation(const math::Complex &input, const Params &params,
             TessellationKind kind, const Prepared &prepared) {
  const float x = params.cell_scale * (input.re * prepared.angle_cos +
                                       input.im * prepared.angle_sin);
  const float y = params.cell_scale * (-input.re * prepared.angle_sin +
                                       input.im * prepared.angle_cos);
  const float distance = tessellation_distance(x, y, kind);
  const float edge = ::math::smooth_ramp(
      params.line_thickness, params.line_thickness + params.line_softness,
      distance);
  return 1.0f - 2.0f * edge;
}

/**
 * @brief Octave noise sample with contrast sharpening.
 * @param noise Noise generator.
 * @param basis Octave basis.
 * @param coordinate Noise-space sample point.
 * @param contrast Nonnegative sharpening; 0 leaves the sample unchanged.
 * @return Sharpened sample in [-1, 1].
 */
HS_O3_FN inline float noise_contour(const FastNoiseLite &noise,
                                    math::NoiseBasis basis,
                                    const math::Vector &coordinate,
                                    float contrast) {
  const float sample = hs::clamp(
      math::sample_noise_octaves(noise, basis, coordinate), -1.0f, 1.0f);
  return sample * (1.0f + contrast) / (1.0f + contrast * fabsf(sample));
}

/**
 * @brief TwinWave source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct TwinWave : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the pattern
   *        frequency, rotation and both phases.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).angle_cos } -> std::convertible_to<float>;
        { State::prepare(frame).angle_sin } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
        { State::prepare(frame).secondary } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  /**
   * @brief Resolves the frame's phases from `State`.
   * @param frame Frame state.
   * @return `State::prepare(frame)`.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return twin_wave(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        prepared);
  }
};

/**
 * @brief Rings source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct Rings : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the pattern
   *        frequency and primary phase.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  /**
   * @brief Resolves the frame's phases from `State`.
   * @param frame Frame state.
   * @return `State::prepare(frame)`.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return rings(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        prepared);
  }
};

/**
 * @brief SphericalRings source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct SphericalRings : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the band
   *        parameters, axis and phase.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).ring_count } -> std::convertible_to<float>;
        { State::params(frame).ring_thickness } -> std::convertible_to<float>;
        { State::params(frame).ring_softness } -> std::convertible_to<float>;
        { State::prepare(frame).axis } -> std::convertible_to<math::Vector>;
        { State::prepare(frame).phase } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  /**
   * @brief Resolves the frame's phases from `State`.
   * @param frame Frame state.
   * @return `State::prepare(frame)`.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  /**
   * @brief Samples the signed field.
   * @param input Sphere carrier; its direction is sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const SphereSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return spherical_rings(input.dir, State::params(frame), prepared);
  }
};

/**
 * @brief Spiral source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct Spiral : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the pattern
   *        frequency, angle and primary phase.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      StateProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).pattern_freq } -> std::convertible_to<float>;
        { State::prepare(frame).angle } -> std::convertible_to<float>;
        { State::prepare(frame).primary } -> std::convertible_to<float>;
      };

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  /**
   * @brief Resolves the frame's phases from `State`.
   * @param frame Frame state.
   * @return `State::prepare(frame)`.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return spiral(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        prepared);
  }
};

/**
 * @brief Grid source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct Grid : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the grid
   *        parameters, rotation and both phases.
   * @tparam CandidateBinding Binding being checked.
   */
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

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  /**
   * @brief Resolves the frame's phases from `State`.
   * @param frame Frame state.
   * @return `State::prepare(frame)`.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    const auto &params = State::params(frame);
    return grid(
        projections::stereo_pattern_args(input.coords, params.pattern_freq),
        params, prepared);
  }
};

/**
 * @brief PrimitiveLattice source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame) accessors.
 */
template <typename State> struct PrimitiveLattice : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the lattice
   *        parameters.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
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

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame) {
    return primitive_lattice(input.coords, State::params(frame));
  }
};

/**
 * @brief EscapeFractal source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State> struct EscapeFractal : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the fractal
   *        parameters and phases.
   * @tparam CandidateBinding Binding being checked.
   */
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
        { State::prepare(frame).secondary } -> std::convertible_to<float>;
        { State::prepare(frame).angle } -> std::convertible_to<float>;
      };

  /// Frame-constant seed and iteration controls.
  using Prepared = PreparedFractal;

  /**
   * @brief Resolves the fractal seed and controls for the frame.
   * @param frame Frame state.
   * @return `prepare_fractal` of `State`'s parameters and phases.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return prepare_fractal(State::params(frame), State::prepare(frame));
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return escape_fractal(input.coords, State::params(frame), prepared);
  }
};

/**
 * @brief Tessellation source policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame) accessors.
 */
template <typename State, TessellationKind KindV>
struct Tessellation : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the line
   *        parameters and rotation.
   * @tparam CandidateBinding Binding being checked.
   */
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

  /// Phases `State::prepare` resolves for the frame.
  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  /**
   * @brief Resolves the frame's phases from `State`.
   * @param frame Frame state.
   * @return `State::prepare(frame)`.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return tessellation(input.coords, State::params(frame), KindV, prepared);
  }
};

/**
 * @brief ProjectedNoise source policy.
 * @tparam State Provider with Binding and FrameState types and
 * noise(frame), noise_scale(frame), noise_time(frame),
 * noise_contrast(frame) accessors.
 */
template <typename State, math::NoiseBasis BasisV>
struct ProjectedNoise : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the noise
   *        generator, scale, time and contrast.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::noise_scale(frame) } -> std::same_as<float>;
        { State::noise_time(frame) } -> std::same_as<float>;
        { State::noise_contrast(frame) } -> std::same_as<float>;
      };

  /// Noise-space loop offset for the frame's noise time.
  using Prepared = math::Vector;

  /**
   * @brief Resolves the loop offset for the frame's noise time.
   * @param frame Frame state.
   * @return Offset added to the noise-space coordinate.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return math::noise_projected_loop_offset(State::noise_time(frame));
  }

  /**
   * @brief Samples the signed field.
   * @param input Plane carrier; its coords are sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return noise_contour(State::noise(frame), BasisV,
                         math::noise_projected_coordinate(
                             input.coords, State::noise_scale(frame), prepared),
                         State::noise_contrast(frame));
  }
};

/**
 * @brief SphericalNoise source policy.
 * @tparam State Provider with Binding and FrameState types and
 * noise_time(frame), noise(frame), noise_scale(frame),
 * noise_contrast(frame) accessors.
 */
template <typename State, math::NoiseBasis BasisV>
struct SphericalNoise : ApproximationDefaults {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies what
   *        `ProjectedNoise` requires.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ProjectedNoise<State, BasisV>::template PROVIDER_VALID<CandidateBinding>;

  /// Noise-space loop offset for the frame's noise time.
  using Prepared = math::Vector;

  /**
   * @brief Resolves the loop offset for the frame's noise time.
   * @param frame Frame state.
   * @return Offset added to the noise-space coordinate.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return math::noise_sphere_loop_offset(State::noise_time(frame));
  }

  /** @brief Post-projection form: samples the plane carrier's retained
      pre-projection point.
      @param input Plane carrier; its `sphere` point is sampled.
      @param frame Frame state.
      @param prepared This frame's `prepare` result.
      @return Field value in [-1, 1]. */
  __attribute__((always_inline)) static float sample(const PlaneSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return noise_contour(State::noise(frame), BasisV,
                         math::noise_sphere_coordinate(
                             input.sphere, State::noise_scale(frame), prepared),
                         State::noise_contrast(frame));
  }
  /**
   * @brief Samples the signed field.
   * @param input Sphere carrier; its direction is sampled.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Field value in [-1, 1].
   */
  __attribute__((always_inline)) static float sample(const SphereSample &input,
                                                     const FrameState &frame,
                                                     const Prepared &prepared) {
    return noise_contour(State::noise(frame), BasisV,
                         math::noise_sphere_coordinate(
                             input.dir, State::noise_scale(frame), prepared),
                         State::noise_contrast(frame));
  }
};

} // namespace Source

} // namespace Pullback
