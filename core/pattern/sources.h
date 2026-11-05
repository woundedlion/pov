/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sources.h
 * @brief Shared analytic source kernels and pure preparation. */

#include "core/pattern/field.h"
#include "core/math/geometry.h"
#include "core/math/projection_patterns.h"

namespace Pattern {

/** @brief Maximum admitted axial cycles per unit coordinate. */
inline constexpr float MAX_AXIAL_FREQUENCY = 64.0f;

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
  using Coordinate = SphereDirection; ///< Admitted unit direction.
  using Params = SineParams;          ///< Axial wave controls.
  using State = PhaseState;           ///< Shared loop phase in turns.
  /** @brief Complete immutable wave inputs. */
  struct Prepared {
    math::Vector axis; ///< Unit wave axis.
    float frequency;   ///< Cycles per unit coordinate.
    float phase_turns; ///< Wrapped wave phase in turns.
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
                         std::fmod(params.speed, 1.0f))};
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
                         std::fmod(params.offset, 1.0f))};
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
  using Coordinate = SphereDirection; ///< Admitted unit direction.
  using Params = InterferenceParams;  ///< Paired wave controls.
  using State = PhaseState;           ///< Shared loop phase in turns.
  /** @brief Complete immutable pair of waves. */
  struct Prepared {
    DirectionalSine::Prepared a; ///< First prepared wave.
    DirectionalSine::Prepared b; ///< Second prepared wave.
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

/** @brief Product of two independently oriented directional waves. */
struct SphereSmoothGrid {
  using Coordinate = SphereDirection; ///< Admitted unit direction.
  using Params = GridParams;          ///< Smooth grid controls.
  using State = PhaseState;           ///< Shared loop phase in turns.
  using Prepared = SphereInterference::Prepared; ///< Prepared wave pair.
  /** @brief Checks both axes and wave controls.
   * @param params Source controls.
   * @return Whether the grid recipe is admitted. */
  static bool valid(const Params &params) {
    return SphereInterference::valid(params);
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
  /** @brief Prepares the independently normalized waves.
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

/** @brief Bundles this frame's source phases with the rotation's cosine pair.
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
 * @param coordinate Pattern-space coordinate.
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
 * @brief Stationary Archimedean spiral sin(r - arms * azimuth).
 * @param coordinate Pattern-space coordinate.
 * @param arms Whole arm count; 0 draws concentric rings.
 * @return Field value in [-1, 1].
 */
HS_HOT_FLASH_MEMBER inline float spiral(const PlanePosition &coordinate,
                                        float arms) {
  const auto &input = coordinate.value;
  const float radius = input.magnitude();
  if (arms == 0.0f)
    return math::fast_sinf(radius);
  return math::fast_sinf(radius - arms * math::fast_atan2(input.im, input.re));
}

/**
 * @brief Coupled sine grid sin(x + c sin(y + phase)) cos(y + c cos(x -
 *        phase)), blended toward the direct product sin(x) cos(y).
 * @tparam Params Parameter block of the source family.
 * @tparam Prepared PreparedSource-shaped per-frame phases.
 * @param coordinate Pattern-space coordinate.
 * @param params Complexity c and pattern mix.
 * @param prepared Per-frame coupling phase.
 * @return Field value in [-1, 1].
 */
template <typename Params, typename Prepared>
HS_HOT_ITCM_MEMBER inline float grid(const PlanePosition &coordinate,
                                     const Params &params,
                                     const Prepared &prepared) {
  const float x = coordinate.value.re;
  const float y = coordinate.value.im;
  if (params.pattern_mix == 1.0f || params.complexity == 0.0f)
    return math::fast_sinf(x) * math::fast_cosf(y);
  const float coupled =
      math::fast_sinf(x + params.complexity *
                              math::fast_sinf(y + prepared.primary)) *
      math::fast_cosf(y + params.complexity *
                              math::fast_cosf(x - prepared.primary));
  if (params.pattern_mix == 0.0f)
    return coupled;
  return hs::lerp(coupled, math::fast_sinf(x) * math::fast_cosf(y),
                  params.pattern_mix);
}

/** @brief Signed distance to a centred square, rounded past its corners.
 * @param x Cell-local x.
 * @param y Cell-local y.
 * @param radius Square inradius.
 * @return Negative inside, positive outside. */
__attribute__((always_inline)) inline float square_distance(float x, float y,
                                                            float radius) {
  const float bx = fabsf(x) - radius;
  const float by = fabsf(y) - radius;
  const float EDGE = fmaxf(bx, by);
  const float CORNER = EDGE < 0.0f ? 0.0f
                                   : sqrtf(fmaxf(bx, 0.0f) * fmaxf(bx, 0.0f) +
                                           fmaxf(by, 0.0f) * fmaxf(by, 0.0f));
  return CORNER + fminf(EDGE, 0.0f);
}

/** @brief Signed distance to a centred hexagon whose edges face the x axis
    and its 60-degree rotations, rounded past its corners.
 * @param x Cell-local x.
 * @param y Cell-local y.
 * @param radius Hexagon inradius.
 * @return Negative inside, positive outside. */
__attribute__((always_inline)) inline float hexagon_distance(float x, float y,
                                                             float radius) {
  constexpr float HALF_SQRT_3 = 0.8660254037844386f;
  constexpr float INV_SQRT_3 = 0.5773502691896258f;
  float across = fabsf(y);
  float along = fabsf(x);
  const float fold = fminf(0.5f * along - HALF_SQRT_3 * across, 0.0f);
  across += 2.0f * HALF_SQRT_3 * fold;
  along -= fold;
  const float half_edge = INV_SQRT_3 * radius;
  across -= fminf(fmaxf(across, -half_edge), half_edge);
  along -= radius;
  return copysignf(sqrtf(across * across + along * along), along);
}

/** @brief Signed distance to a centred upward equilateral triangle, rounded
    past its corners.
 * @param x Cell-local x.
 * @param y Cell-local y; the apex lies on +y.
 * @param radius Triangle inradius.
 * @return Negative inside, positive outside. */
__attribute__((always_inline)) inline float triangle_distance(float x, float y,
                                                              float radius) {
  constexpr float SQRT_3 = 1.7320508075688772f;
  const float half_side = SQRT_3 * radius;
  float px = fabsf(x) - half_side;
  float py = y + radius;
  if (px + SQRT_3 * py > 0.0f) {
    const float folded = 0.5f * (px - SQRT_3 * py);
    py = 0.5f * (-SQRT_3 * px - py);
    px = folded;
  }
  px -= fminf(fmaxf(px, -2.0f * half_side), 0.0f);
  return -copysignf(sqrtf(px * px + py * py), py);
}

/** @brief Field of a primitive blended from a circle toward a cell polygon.
 * @tparam Kind Cell shape; its polygon is the primitive at full blend.
 * @tparam Params TessellationCellParams-shaped parameter block.
 * @param x Cell-local x, in cell units.
 * @param y Cell-local y, in cell units.
 * @param params Shape blend, radius and softness.
 * @return 1 inside the primitive, -1 outside, ramped across the edge. */
template <TessellationKind Kind, typename Params>
__attribute__((always_inline)) inline float
cell_primitive(float x, float y, const Params &params) {
  const auto polygon = [](float px, float py,
                          float radius) __attribute__((always_inline)) {
    if constexpr (Kind == TessellationKind::SQUARE)
      return square_distance(px, py, radius);
    else if constexpr (Kind == TessellationKind::HEXAGONAL)
      return hexagon_distance(px, py, radius);
    else
      return triangle_distance(px, py, radius);
  };
  const auto evaluate = [&](float blend) __attribute__((always_inline)) {
    const float circle = sqrtf(x * x + y * y) - params.cell_radius;
    if (blend == 0.0f)
      return 1.0f - 2.0f * ::math::smooth_ramp(-params.cell_softness,
                                               params.cell_softness, circle);
    const float distance =
        hs::lerp(circle, polygon(x, y, params.cell_radius), blend);
    return 1.0f - 2.0f * ::math::smooth_ramp(-params.cell_softness,
                                             params.cell_softness, distance);
  };
  if (params.cell_shape == 1.0f)
    return evaluate(1.0f);
  return evaluate(params.cell_shape);
}

/** @brief Offset from the nearest point of the hexagonal lattice with basis
    (1, 0), (1/2, sqrt(3)/2).
 * @param u Plane x in lattice periods.
 * @param v Plane y in lattice periods.
 * @param x Receives the offset's x, in periods.
 * @param y Receives the offset's y, in periods. */
__attribute__((always_inline)) inline void
nearest_hexagonal_offset(float u, float v, float &x, float &y) {
  constexpr float SQRT_3 = 1.7320508075688772f;
  constexpr float INV_SQRT_3 = 0.5773502691896258f;
  const float row = v * INV_SQRT_3;
  const float ax = math::wrap_t(u + 0.5f) - 0.5f;
  const float ay = SQRT_3 * (math::wrap_t(row + 0.5f) - 0.5f);
  const float bx = math::wrap_t(u) - 0.5f;
  const float by = SQRT_3 * (math::wrap_t(row) - 0.5f);
  const bool first = ax * ax + ay * ay <= bx * bx + by * by;
  x = first ? ax : bx;
  y = first ? ay : by;
}

/** @brief Offset from the centroid of the triangle containing a point, in
    that triangle's apex-up frame.
 * @details Triangle vertices lie on the lattice of `nearest_hexagonal_offset`.
 * @param u Plane x in lattice periods.
 * @param v Plane y in lattice periods.
 * @param x Receives the offset's x, in cell units.
 * @param y Receives the offset's y, in cell units. */
__attribute__((always_inline)) inline void
triangular_cell_offset(float u, float v, float &x, float &y) {
  constexpr float SQRT_3 = 1.7320508075688772f;
  constexpr float INV_SQRT_3 = 0.5773502691896258f;
  const float b = 2.0f * INV_SQRT_3 * v;
  const float fa = math::wrap_t(u - 0.5f * b);
  const float fb = math::wrap_t(b);
  const bool up = fa + fb < 1.0f;
  const float centroid = up ? 1.0f / 3.0f : 2.0f / 3.0f;
  const float da = fa - centroid;
  const float db = fb - centroid;
  x = SQRT_3 * (da + 0.5f * db);
  y = (up ? 1.5f : -1.5f) * db;
}

/**
 * @brief Primitive centred on each cell of the square lattice.
 * @tparam Params TessellationCellParams-shaped parameter block.
 * @param coordinate Plane coordinate.
 * @param params Cell scale, shape blend, radius and softness.
 * @return 1 inside the primitive, -1 outside, ramped across the edge.
 */
template <typename Params>
HS_HOT_ITCM_MEMBER inline float square_cells(const PlanePosition &coordinate,
                                             const Params &params) {
  const auto &input = coordinate.value;
  const float x = math::wrap_t(params.cell_scale * input.re + 0.5f) - 0.5f;
  const float y = math::wrap_t(params.cell_scale * input.im + 0.5f) - 0.5f;
  return cell_primitive<TessellationKind::SQUARE>(x, y, params);
}

/**
 * @brief Primitive centred on each cell of the hexagonal lattice.
 * @tparam Params TessellationCellParams-shaped parameter block.
 * @param coordinate Plane coordinate.
 * @param params Cell scale, shape blend, radius and softness.
 * @return 1 inside the primitive, -1 outside, ramped across the edge.
 */
template <typename Params>
HS_HOT_FLASH_MEMBER inline float
hexagonal_cells(const PlanePosition &coordinate, const Params &params) {
  float x, y;
  nearest_hexagonal_offset(params.cell_scale * coordinate.value.re,
                           params.cell_scale * coordinate.value.im, x, y);
  return cell_primitive<TessellationKind::HEXAGONAL>(x, y, params);
}

/**
 * @brief Primitive centred on each triangle of the triangular tessellation.
 * @tparam Params TessellationCellParams-shaped parameter block.
 * @param coordinate Plane coordinate.
 * @param params Cell scale, shape blend, radius and softness.
 * @return 1 inside the primitive, -1 outside, ramped across the edge.
 */
template <typename Params>
HS_HOT_FLASH_MEMBER inline float
triangular_cells(const PlanePosition &coordinate, const Params &params) {
  float x, y;
  triangular_cell_offset(params.cell_scale * coordinate.value.re,
                         params.cell_scale * coordinate.value.im, x, y);
  return cell_primitive<TessellationKind::TRIANGULAR>(x, y, params);
}

/**
 * @brief Distance to the nearest integer.
 * @param coordinate Lattice coordinate.
 * @return Distance in [0, 0.5].
 */
HS_O3_FN inline float distance_to_lattice_line(float coordinate) {
  return fabsf(math::wrap_t(coordinate + 0.5f) - 0.5f);
}

/** @brief Distance from a plane point to the nearest cell edge of @p Kind.
    @tparam Kind Cell shape.
    @param u Plane x in lattice periods.
    @param v Plane y in lattice periods.
    @return Distance in cell units; one half at a cell centre. */
template <TessellationKind Kind>
__attribute__((always_inline)) inline float
tessellation_edge_distance(float u, float v) {
  constexpr float HALF_SQRT_3 = 0.8660254037844386f;
  if constexpr (Kind == TessellationKind::SQUARE) {
    const float cell_x = math::wrap_t(u + 0.5f) - 0.5f;
    const float cell_y = math::wrap_t(v + 0.5f) - 0.5f;
    return 0.5f - fmaxf(fabsf(cell_x), fabsf(cell_y));
  } else if constexpr (Kind == TessellationKind::HEXAGONAL) {
    float x, y;
    nearest_hexagonal_offset(u, v, x, y);
    return 0.5f - fmaxf(fabsf(x), fmaxf(fabsf(0.5f * x + HALF_SQRT_3 * y),
                                        fabsf(0.5f * x - HALF_SQRT_3 * y)));
  } else {
    const float b = v / HALF_SQRT_3;
    const float a = u - 0.5f * b;
    return 1.5f * fminf(distance_to_lattice_line(a),
                        fminf(distance_to_lattice_line(b),
                              distance_to_lattice_line(a + b)));
  }
}

/**
 * @brief Antialiased edge band of a periodic polygon tessellation.
 * @tparam Kind Cell shape.
 * @tparam Params TessellationEdgeParams-shaped parameter block.
 * @param coordinate Plane coordinate.
 * @param params Cell scale, line thickness and softness.
 * @return 1 on an edge, -1 inside a cell, ramped between.
 */
template <TessellationKind Kind, typename Params>
HS_HOT_FLASH_MEMBER inline float
tessellation_edges(const PlanePosition &coordinate, const Params &params) {
  const float distance =
      tessellation_edge_distance<Kind>(params.cell_scale * coordinate.value.re,
                                       params.cell_scale * coordinate.value.im);
  const float edge = ::math::smooth_ramp(
      params.line_thickness, params.line_thickness + params.line_softness,
      distance);
  return 1.0f - 2.0f * edge;
}

/**
 * @brief Static periodic tessellation field.
 * @tparam Kind Cell shape.
 * @tparam Fill Cell primitives or edge band.
 * @tparam Params Parameter block with the controls @p Fill reads.
 * @param coordinate Plane coordinate.
 * @param params Tessellation controls.
 * @return Field value in [-1, 1].
 */
template <TessellationKind Kind, TessellationFill Fill, typename Params>
__attribute__((always_inline)) inline float
tessellation(const PlanePosition &coordinate, const Params &params) {
  if constexpr (Fill == TessellationFill::EDGES)
    return tessellation_edges<Kind>(coordinate, params);
  else if constexpr (Kind == TessellationKind::SQUARE)
    return square_cells(coordinate, params);
  else if constexpr (Kind == TessellationKind::HEXAGONAL)
    return hexagonal_cells(coordinate, params);
  else
    return triangular_cells(coordinate, params);
}

/**
 * @brief Static periodic tessellation field with a run-time kind and fill.
 * @tparam Params TessellationSourceParams-shaped parameter block.
 * @param coordinate Plane coordinate.
 * @param params Tessellation controls.
 * @param kind Checked cell shape.
 * @param fill Checked fill.
 * @return Field value in [-1, 1].
 */
template <typename Params>
inline float tessellation(const PlanePosition &coordinate, const Params &params,
                          TessellationKind kind, TessellationFill fill) {
  using K = TessellationKind;
  using F = TessellationFill;
  if (fill == F::EDGES) {
    if (kind == K::TRIANGULAR)
      return tessellation<K::TRIANGULAR, F::EDGES>(coordinate, params);
    if (kind == K::SQUARE)
      return tessellation<K::SQUARE, F::EDGES>(coordinate, params);
    return tessellation<K::HEXAGONAL, F::EDGES>(coordinate, params);
  }
  if (kind == K::TRIANGULAR)
    return tessellation<K::TRIANGULAR, F::CELLS>(coordinate, params);
  if (kind == K::SQUARE)
    return tessellation<K::SQUARE, F::CELLS>(coordinate, params);
  return tessellation<K::HEXAGONAL, F::CELLS>(coordinate, params);
}

/** @brief Checks a finite source parameter against an admitted interval.
 * @param value Parameter value.
 * @param low Inclusive lower bound.
 * @param high Inclusive upper bound.
 * @return Whether the parameter is admitted. */
inline bool parameter_in_range(float value, float low, float high) {
  return finite(value) && value >= low && value <= high;
}

/** @brief Validates a source parameter record's declared controls.
 * @tparam Params Source parameter record.
 * @param params Source controls.
 * @return Whether all declared controls are admitted. */
template <class Params> bool valid_source_params(const Params &params) {
  bool valid = true;
  if constexpr (requires { params.pattern_freq; })
    valid &= parameter_in_range(params.pattern_freq, 0.01f, 64.0f);
  if constexpr (requires { params.speed; })
    valid &= parameter_in_range(params.speed, -0.5f, 0.5f);
  if constexpr (requires { params.secondary_rate; })
    valid &= parameter_in_range(params.secondary_rate, 0.0f, 2.0f);
  if constexpr (requires { params.angle_rate; })
    valid &= parameter_in_range(params.angle_rate, -0.05f, 0.05f);
  if constexpr (requires { params.angle; })
    valid &= parameter_in_range(params.angle, 0.0f, math::TWO_PI_F);
  if constexpr (requires { params.arms; })
    valid &= parameter_in_range(params.arms, 0.0f, 12.0f);
  if constexpr (requires { params.complexity; })
    valid &= parameter_in_range(params.complexity, 0.0f, 3.0f) &&
             parameter_in_range(params.pattern_mix, 0.0f, 1.0f);
  if constexpr (requires { params.cell_scale; })
    valid &= parameter_in_range(params.cell_scale, 1.0f / 64.0f, 8.0f);
  if constexpr (requires { params.cell_radius; })
    valid &= parameter_in_range(params.cell_shape, 0.0f, 1.0f) &&
             parameter_in_range(params.cell_softness, 1.0f / 1024.0f, 1.0f) &&
             parameter_in_range(params.cell_radius, 1.0f / 64.0f, 0.49f);
  if constexpr (requires { params.line_thickness; })
    valid &= parameter_in_range(params.line_thickness, 1.0f / 1024.0f, 0.25f) &&
             parameter_in_range(params.line_softness, 1.0f / 1024.0f, 0.25f);
  return valid;
}

/** @brief Pure radian-clock protocol shared by analytic sources.
 * @tparam Parameters The source's passive parameter record. */
template <class Parameters> struct RadianSource {
  using Coordinate =
      PlanePosition;         ///< Plane coordinate supplied by the consumer.
  using Params = Parameters; ///< Canonical passive source parameters.
  using State = RadianState; ///< Primary, secondary and angle clocks.
  /** @brief Complete immutable controls and frame-constant phase math. */
  struct Prepared {
    Params params;         ///< Owned source controls.
    PreparedSource phases; ///< Prepared source clocks and rotation.
  };
  /** @brief Validates source controls.
   * @param params Source controls.
   * @return Whether preparation admits the controls. */
  static bool valid(const Params &params) {
    return valid_source_params(params);
  }
  /** @brief Initializes source-specific radian clocks.
   * @param params Source controls.
   * @return Zero phases. */
  static State init_state(const Params &params) {
    (void)params;
    return {};
  }
  /** @brief Advances the declared clocks; the secondary rate scales with speed.
   * @param state Previous clocks.
   * @param params Source controls.
   * @return Wrapped next clocks. */
  static State advance(const State &state, const Params &params) {
    HS_CHECK(valid(params) && valid_state(state),
             "pattern radian source: invalid advance inputs");
    State next = state;
    advance_clocks(params, next.primary, next.secondary, next.angle);
    return next;
  }
  /** @brief Owns controls and resolves frame-constant rotation math.
   * @param params Source controls.
   * @param state Explicit source clocks.
   * @return Immutable source inputs. */
  static Prepared prepare(const Params &params, const State &state) {
    HS_CHECK(valid(params) && valid_state(state),
             "pattern radian source: invalid preparation inputs");
    return {params,
            Pattern::prepare(state.primary, state.secondary, state.angle)};
  }
  /** @brief Checks finite bounded radian clocks.
   * @param state Source clocks.
   * @return Whether each clock is within its declared period. */
  static bool valid_state(const State &state) {
    return parameter_in_range(state.primary, -math::TWO_PI_F, math::TWO_PI_F) &&
           parameter_in_range(state.secondary, -math::TWO_PI_F,
                              math::TWO_PI_F) &&
           parameter_in_range(state.angle, -math::TWO_PI_F, math::TWO_PI_F);
  }
};

/** @brief Pure protocol shared by stationary sources: empty
    passive state and immutable prepared controls.
 * @tparam Parameters The source's passive parameter record. */
template <class Parameters> struct StationarySource {
  using Coordinate = PlanePosition; ///< Consumer-supplied plane coordinate.
  using Params = Parameters;        ///< Canonical passive source parameters.
  using State = EmptyState;         ///< Stationary source state.
  using Prepared = Params;          ///< Complete immutable source inputs.
  /** @brief Checks source controls.
   * @param params Source controls.
   * @return Whether the controls are admitted. */
  static bool valid(const Params &params) {
    return valid_source_params(params);
  }
  /** @brief Initializes empty source state.
   * @param params Source controls.
   * @return Empty state. */
  static State init_state(const Params &params) {
    (void)params;
    return {};
  }
  /** @brief Preserves empty source state.
   * @param state Previous state.
   * @param params Source controls.
   * @return Empty state. */
  static State advance(const State &state, const Params &params) {
    (void)params;
    return state;
  }
};

/** @brief Parallel sine bands in the consumer's projection plane. */
struct PlaneBands : StationarySource<BandsSourceParams> {
  /** @brief Immutable frequency and band orientation. */
  struct Prepared {
    float pattern_freq; ///< Radians per plane-coordinate unit.
    float angle_cos;    ///< Cosine of the band angle.
    float angle_sin;    ///< Sine of the band angle.
  };
  /** @brief Prepares the band's frequency and orientation.
   * @param params Source controls.
   * @param state Empty source state.
   * @return Immutable source inputs. */
  static Prepared prepare(const Params &params, const State &state) {
    (void)state;
    HS_CHECK(valid(params), "pattern bands: invalid parameters");
    return {params.pattern_freq, std::cos(params.angle),
            std::sin(params.angle)};
  }
  /** @brief Samples parallel sine bands; angle zero varies along Y.
   * @param coordinate Plane coordinate.
   * @param prepared Immutable controls.
   * @return Signed field value. */
  HS_HOT_FLASH_MEMBER static float sample(const Coordinate &coordinate,
                                          const Prepared &prepared) {
    const auto point = projections::stereo_pattern_args(coordinate.value,
                                                        prepared.pattern_freq);
    return math::fast_sinf(point.im * prepared.angle_cos -
                           point.re * prepared.angle_sin);
  }
};

/** @brief Pure plane source policy for the stationary spiral kernel. */
struct PlaneSpiral : StationarySource<SpiralSourceParams> {
  /** @brief Owns admitted controls with the arm count rounded to whole.
   * @param params Source controls.
   * @param state Empty source state.
   * @return Immutable source inputs. */
  static Prepared prepare(const Params &params, const State &state) {
    (void)state;
    HS_CHECK(valid(params), "pattern spiral: invalid parameters");
    Prepared prepared = params;
    prepared.arms = std::round(params.arms);
    return prepared;
  }
  /** @brief Evaluates the prepared source at a consumer-supplied plane point.
   * @param coordinate Plane coordinate.
   * @param prepared Immutable controls.
   * @return Signed field value. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return spiral(PlanePosition{projections::stereo_pattern_args(
                      coordinate.value, prepared.pattern_freq)},
                  prepared.arms);
  }
};

/** @brief Pure plane source policy for the twin_wave kernel. */
struct PlaneInterference : RadianSource<TwinWaveSourceParams> {
  /** @brief Evaluates the prepared source at a consumer-supplied plane point.
   * @param coordinate Plane coordinate.
   * @param prepared Immutable controls and clocks.
   * @return Signed field value. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return twin_wave(PlanePosition{projections::stereo_pattern_args(
                         coordinate.value, prepared.params.pattern_freq)},
                     prepared.phases);
  }
};

/** @brief Pure plane source policy for the grid kernel. */
struct PlaneCoupledGrid : RadianSource<GridSourceParams> {
  /** @brief Evaluates the prepared source at a consumer-supplied plane point.
   * @param coordinate Plane coordinate.
   * @param prepared Immutable controls and clocks.
   * @return Signed field value. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return grid(PlanePosition{projections::stereo_pattern_args(
                    coordinate.value, prepared.params.pattern_freq)},
                prepared.params, prepared.phases);
  }
};

/** @brief Pure static tessellation policy with empty passive state.
 * @tparam Kind Cell shape.
 * @tparam Fill Cell primitives or edge band. */
template <TessellationKind Kind, TessellationFill Fill>
struct PlaneTessellation : StationarySource<TessellationSourceParams> {
  static_assert(Kind <= TessellationKind::HEXAGONAL &&
                Fill <= TessellationFill::EDGES);
  /** @brief Owns admitted tessellation controls.
   * @param params Source controls.
   * @param state Empty source state.
   * @return Immutable source inputs. */
  static Prepared prepare(const Params &params, const State &state) {
    (void)state;
    HS_CHECK(valid(params), "pattern tessellation: invalid parameters");
    return params;
  }
  /** @brief Evaluates the canonical tessellation field.
   * @param coordinate Plane position.
   * @param prepared Immutable tessellation controls.
   * @return Signed tessellation field. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return tessellation<Kind, Fill>(coordinate, prepared);
  }
};

} // namespace Pattern
