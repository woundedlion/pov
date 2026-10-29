/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sources.h
 * @brief Shared analytic source kernels and pure preparation. */

#include "core/pattern/field.h"
#include "core/math/geometry.h"

namespace Pattern {

/** @brief Maximum admitted axial cycles per unit coordinate. */
inline constexpr float MAX_AXIAL_FREQUENCY = 64.0f;
/** @brief Orthogonality tolerance for normalized grid axes. */
inline constexpr float GRID_AXIS_TOLERANCE = 1.0e-5f;

/** @brief Checks finite sine parameters and the qualified axial range.
 * @param params Axial wave controls.
 * @return Whether preparation can admit the recipe. */
inline bool valid_sine(const SineParams &params) {
  return SphereDirection::valid(params.axis) && finite(params.frequency) &&
         params.frequency >= 0.0f && params.frequency <= MAX_AXIAL_FREQUENCY &&
         finite(params.offset) && finite(params.speed);
}

/** @brief Directional sine with an explicit turn clock. */
struct DirectionalSine {
  using Coordinate = SphereDirection;
  using Params = SineParams;
  using State = PhaseState;
  /** @brief Complete immutable wave inputs. */
  struct Prepared {
    math::Vector axis;
    float frequency;
    float phase_turns;
  };
  /** @brief Checks a source recipe.
   * @param params Source controls.
   * @return Whether the recipe is admitted. */
  static bool valid(const Params &params) { return valid_sine(params); }
  /** @brief Initializes the source clock.
   * @param params Source controls.
   * @return Zero phase in turns. */
  static State init_state(const Params &params) {
    (void)params;
    return {};
  }
  /** @brief Advances the source clock by one engine update.
   * @param state Previous phase.
   * @param params Source controls; speed is turns per update.
   * @return Wrapped phase. */
  static State advance(const State &state, const Params &params) {
    HS_CHECK(valid(params) && finite(state.phase_turns),
             "pattern sine: invalid state");
    return {math::wrap_t(math::wrap_t(state.phase_turns) +
                         math::wrap_t(params.speed))};
  }
  /** @brief Resolves unit axis and wrapped phase.
   * @param params Source controls.
   * @param state Explicit source clock.
   * @return Complete immutable source inputs. */
  static Prepared prepare(const Params &params, const State &state) {
    HS_CHECK(valid(params) && finite(state.phase_turns),
             "pattern sine: invalid parameters");
    return {SphereDirection::from(params.axis).value(), params.frequency,
            math::wrap_t(math::wrap_t(state.phase_turns) +
                         math::wrap_t(params.offset))};
  }
  /** @brief Evaluates an axial wave with a range-reduced argument.
   * @param coordinate Admitted unit direction.
   * @param prepared Immutable wave inputs.
   * @return Signed wave sample. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    const float turns =
        prepared.frequency * math::dot(prepared.axis, coordinate.value()) +
        prepared.phase_turns;
    return math::sinf_0_2pi(math::TWO_PI_F * math::wrap_t(turns));
  }
};

/** @brief Half-sum of two directional waves sharing one turn clock. */
struct SphereInterference {
  using Coordinate = SphereDirection;
  using Params = InterferenceParams;
  using State = PhaseState;
  /** @brief Complete immutable pair of waves. */
  struct Prepared {
    DirectionalSine::Prepared a;
    DirectionalSine::Prepared b;
  };
  /** @brief Checks both axes and wave controls.
   * @param params Source controls.
   * @return Whether both waves are admitted. */
  static bool valid(const Params &params) {
    return valid_sine(params) && valid_sine({params.axis_b, params.frequency_b,
                                             params.offset_b, params.speed});
  }
  /** @brief Initializes the shared clock.
   * @param params Source controls.
   * @return Zero phase. */
  static State init_state(const Params &params) {
    (void)params;
    return {};
  }
  /** @brief Advances the shared turn clock.
   * @param state Previous phase.
   * @param params Source controls.
   * @return Wrapped phase. */
  static State advance(const State &state, const Params &params) {
    HS_CHECK(valid(params), "pattern interference: invalid advance parameters");
    return DirectionalSine::advance(state, params);
  }
  /** @brief Prepares both waves with their shared clock.
   * @param params Source controls.
   * @param state Explicit source clock.
   * @return Immutable wave pair. */
  static Prepared prepare(const Params &params, const State &state) {
    HS_CHECK(valid(params),
             "pattern interference: invalid preparation parameters");
    return {DirectionalSine::prepare(params, state),
            DirectionalSine::prepare({params.axis_b, params.frequency_b,
                                      params.offset_b, params.speed},
                                     state)};
  }
  /** @brief Samples the fixed half-gain interference.
   * @param coordinate Admitted unit direction.
   * @param prepared Immutable wave pair.
   * @return Scalar in [-1, 1]. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return 0.5f * (DirectionalSine::sample(coordinate, prepared.a) +
                   DirectionalSine::sample(coordinate, prepared.b));
  }
};

/** @brief Product of two orthogonal directional waves. */
struct SphereSmoothGrid {
  using Coordinate = SphereDirection;
  using Params = GridParams;
  using State = PhaseState;
  using Prepared = SphereInterference::Prepared;
  /** @brief Checks the qualified waves and unit-axis orthogonality.
   * @param params Source controls.
   * @return Whether the grid recipe is admitted. */
  static bool valid(const Params &params) {
    return SphereInterference::valid(params) &&
           fabsf(math::dot(SphereDirection::from(params.axis).value(),
                           SphereDirection::from(params.axis_b).value())) <=
               GRID_AXIS_TOLERANCE;
  }
  /** @brief Initializes the shared turn clock.
   * @param params Source controls.
   * @return Zero phase. */
  static State init_state(const Params &params) {
    (void)params;
    return {};
  }
  /** @brief Advances the shared turn clock.
   * @param state Previous phase.
   * @param params Source controls.
   * @return Wrapped phase. */
  static State advance(const State &state, const Params &params) {
    HS_CHECK(valid(params), "pattern grid: invalid advance parameters");
    return SphereInterference::advance(state, params);
  }
  /** @brief Prepares the admitted orthogonal waves.
   * @param params Source controls.
   * @param state Explicit source clock.
   * @return Immutable wave pair. */
  static Prepared prepare(const Params &params, const State &state) {
    HS_CHECK(valid(params), "pattern grid: invalid preparation parameters");
    return SphereInterference::prepare(params, state);
  }
  /** @brief Samples the smooth sine product.
   * @param coordinate Admitted unit direction.
   * @param prepared Immutable wave pair.
   * @return Scalar in [-1, 1]. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return DirectionalSine::sample(coordinate, prepared.a) *
           DirectionalSine::sample(coordinate, prepared.b);
  }
};

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

/** @brief The source stage's phases, resolved once per frame. */
struct PreparedSource {
  float primary;   /**< Primary phase, wrapped into (-2pi,2pi). */
  float secondary; /**< Secondary phase, wrapped into (-2pi,2pi). */
  float angle;     /**< Source rotation, wrapped into (-2pi,2pi) radians. */
  float angle_cos; /**< Cosine of `angle`. */
  float angle_sin; /**< Sine of `angle`. */
};

/** @brief Grid phases and the stationary uncoupled kernel selection. */
struct PreparedGrid : PreparedSource {
  bool stationary_uncoupled; ///< Identity rotation, zero phases and coupling.
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

/**
 * @brief Two-wave interference: waves along x and along
 *        the axis at the source angle.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param input Pattern-space coordinate.
 * @param prepared Per-frame phases and rotation.
 * @return Field value in [-1, 1].
 */
template <typename Prepared>
HS_HOT_FLASH_MEMBER inline float twin_wave(const PlanePosition &coordinate,
                                           const Prepared &prepared) {
  const auto &input = coordinate.value;
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
HS_HOT_FLASH_MEMBER inline float rings(const PlanePosition &coordinate,
                                       const Prepared &prepared) {
  const auto &input = coordinate.value;
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
spherical_rings(const SphereDirection &coordinate, const Params &params,
                const PreparedSphericalRings &prepared) {
  const auto &input = coordinate.value();
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
HS_HOT_FLASH_MEMBER inline float spiral(const PlanePosition &coordinate,
                                        const Prepared &prepared) {
  const auto &input = coordinate.value;
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
HS_HOT_ITCM_MEMBER inline float grid(const PlanePosition &coordinate,
                                     const Params &params,
                                     const Prepared &prepared) {
  const auto &input = coordinate.value;
  if constexpr (requires { prepared.stationary_uncoupled; })
    if (prepared.stationary_uncoupled)
      return math::fast_sinf(input.re) * math::fast_cosf(input.im);
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
HS_HOT_FLASH_MEMBER inline float
primitive_lattice(const PlanePosition &coordinate, const Params &params) {
  const auto &input = coordinate.value;
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
HS_HOT_FLASH_MEMBER inline float escape_fractal(const PlanePosition &coordinate,
                                                const Params &params,
                                                const Prepared &prepared) {
  const auto &input = coordinate.value;
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
tessellation(const PlanePosition &coordinate, const Params &params,
             TessellationKind kind, const Prepared &prepared) {
  const auto &input = coordinate.value;
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

} // namespace Pattern
