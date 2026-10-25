/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file surface.h
 * @brief The two-arm sampling surface and simulation time windows.
 */

#include "core/animation/orientation.h"
#include "core/math/3dmath.h"
#include "core/platform/attributes.h"
#include "core/render/canvas.h"
#include <cmath>

namespace Spacetime {

#if defined(__EMSCRIPTEN__) || HS_ENABLE_TEST_HOOKS || HS_ENABLE_TEMPORAL_SWITCH
namespace detail {
inline bool flat = false; ///< Simulator sweep-time sampling selection.
}

/** @brief Selects sweep-time sampling for simulator presentation.
 * @param flat Whether every column samples the frame end. */
inline void set_flat_surface(bool flat) { detail::flat = flat; }

/** @brief Whether every column samples the frame end.
 * @return The sampling surface selection. */
inline bool flat_surface() { return detail::flat; }
#else
/** @brief Whether every column samples the frame end.
 * @return False for physical display rendering. */
constexpr bool flat_surface() { return false; }
#endif

/** @brief Column time within a two-arm half-turn.
 * @param column Absolute canvas column.
 * @param width Even canvas width.
 * @return Fraction of one effect frame, in [0, 1). */
constexpr float surface_time(int column, int width) {
  return static_cast<float>(column % (width / 2)) / (width / 2);
}

/// Longest rotation rotate_by() takes by its series coefficients, radians.
inline constexpr float SHORT_TURN_ANGLE = 0.25f;

/** @brief Rodrigues coefficients of a rotation by angle `a`. */
struct TurnCoefficients {
  float sine;    ///< sin(a) / a.
  float versine; ///< (1 - cos(a)) / a^2.
};

/**
 * @brief Rodrigues coefficients of a rotation past SHORT_TURN_ANGLE.
 * @param angle_sq Squared rotation angle, below pi squared.
 * @return The coefficients, from the series at the angle halved until it is
 *         short, doubled back by the double-angle formulas.
 */
HS_FLASH_MEMBER inline TurnCoefficients long_turn_coefficients(float angle_sq) {
  const float angle = sqrtf(angle_sq);
  float half = angle;
  int doublings = 0;
  while (half > SHORT_TURN_ANGLE) {
    half *= 0.5f;
    ++doublings;
  }
  const float h2 = half * half;
  float s = half * (1.0f - h2 * (1.0f / 6.0f - h2 / 120.0f));
  float c = 1.0f - h2 * (0.5f - h2 * (1.0f / 24.0f - h2 / 720.0f));
  float versine = 1.0f - c;
  for (int i = 0; i < doublings; ++i) {
    versine = 2.0f * s * s;
    const float next = c * c - s * s;
    s = 2.0f * s * c;
    c = next;
  }
  return {s / angle, versine / angle_sq};
}

/**
 * @brief Rodrigues coefficients of a rotation shorter than a half turn.
 * @param angle_sq Squared rotation angle.
 * @return Series coefficients up to SHORT_TURN_ANGLE, long_turn_coefficients()
 *         past it.
 */
__attribute__((always_inline)) inline TurnCoefficients
turn_coefficients(float angle_sq) {
  if (angle_sq <= SHORT_TURN_ANGLE * SHORT_TURN_ANGLE)
    // Series for sin(a) / a and (1 - cos(a)) / a^2.
    return {1.0f - angle_sq * (1.0f / 6.0f - angle_sq / 120.0f),
            0.5f - angle_sq * (1.0f / 24.0f - angle_sq / 720.0f)};
  return long_turn_coefficients(angle_sq);
}

/**
 * @brief Rotates a direction by a rotation vector with its coefficients.
 * @param p Direction to rotate.
 * @param rotation Axis times angle.
 * @param k The rotation's coefficients, as turn_coefficients() gives them.
 * @return @p p rotated by Rodrigues' formula.
 */
__attribute__((always_inline)) inline math::Vector
rotate_by(const math::Vector &p, const math::Vector &rotation,
          const TurnCoefficients &k) {
  const math::Vector across = math::cross(rotation, p);
  return p + across * k.sine + math::cross(rotation, across) * k.versine;
}

/**
 * @brief Rotates a direction by a rotation vector.
 * @param p Direction to rotate.
 * @param rotation Axis times angle, shorter than a half turn.
 * @return @p p rotated by Rodrigues' formula.
 */
__attribute__((always_inline)) inline math::Vector
rotate_by(const math::Vector &p, const math::Vector &rotation) {
  return rotate_by(p, rotation,
                   turn_coefficients(math::dot(rotation, rotation)));
}

/** @brief First-arm column lit at a sweep fraction.
 * @param fraction A column's sweep fraction, as surface_time() returns it.
 * @param width Even canvas width.
 * @return The column in [0, width / 2) lit at @p fraction. */
constexpr int surface_column(float fraction, int width) {
  return static_cast<int>(fraction * static_cast<float>(width / 2) + 0.5f);
}

/**
 * @brief First half-turn column at or after a frame fraction.
 * @param fraction Frame fraction in [0, 1].
 * @param width Even canvas width.
 * @return Smallest column `c` in [0, width / 2] with
 *         `surface_time(c, width) >= fraction`, or width / 2 when none is.
 */
HS_FLASH_MEMBER inline int first_column_at(float fraction, int width) {
  const int half = width / 2;
  int column = static_cast<int>(
      std::ceil(hs::clamp(fraction, 0.0f, 1.0f) * static_cast<float>(half)));
  while (column > 0 && surface_time(column - 1, width) >= fraction)
    --column;
  while (column < half && surface_time(column, width) < fraction)
    ++column;
  return column;
}

/** @brief Half-turn columns of one sweep slab. */
struct TimeSlab {
  int begin; ///< First half-turn column.
  int end;   ///< One past the last half-turn column.
};

/** @brief Columns of one of `count` slabs that tile the half-turn in sweep
 * order, for a render-x window.
 * @param slab Slab index, in [0, count).
 * @param count Slab count, in [1, width / 2].
 * @param width Even canvas width.
 * @return Slab `slab`, whose columns' sweep fractions lie in
 * `[slab / count, (slab + 1) / count)` up to column rounding. */
constexpr TimeSlab time_slab(int slab, int count, int width) {
  return {slab * (width / 2) / count, (slab + 1) * (width / 2) / count};
}

/**
 * @brief Advances a simulation and shades it in sweep slabs, each slab after
 *        the substeps that end its share of the frame.
 * @tparam SUBSTEPS Substeps per frame.
 * @tparam SLABS Slab count, a divisor of `SUBSTEPS`.
 * @param canvas Target canvas; each slab shades under a TimeWindow.
 * @param advance Callable (int substeps) advancing the simulation.
 * @param shade Callable () shading the current state into the open window.
 * @details With temporal sampling disabled, runs every substep and then one
 * frame-wide shade.
 */
template <int SUBSTEPS, int SLABS, typename Advance, typename Shade>
HS_HOT_INLINE void simulate_in_slabs(Canvas &canvas, Advance &&advance,
                                     Shade &&shade) {
  static_assert(SLABS > 0 && SUBSTEPS % SLABS == 0);
  if (flat_surface()) {
    advance(SUBSTEPS);
    shade();
    return;
  }
  for (int slab = 0; slab < SLABS; ++slab) {
    advance(SUBSTEPS / SLABS);
    const TimeSlab columns = time_slab(slab, SLABS, canvas.width());
    ::TimeWindow window(canvas, columns.begin, columns.end);
    shade();
  }
}

/** @brief Two-arm map between columns and frame time.
 * @tparam W Even canvas width. */
template <int W> struct SamplingSurface {
  static constexpr int TIMES = W / 2; ///< Distinct surface times.
  /** @brief Time index of a column.
   * @param x Wrapped canvas column.
   * @return Shared half-turn time index. */
  static constexpr int index(int x) { return x % TIMES; }
  /** @brief Arm containing a column.
   * @param x Wrapped canvas column.
   * @return Arm index. */
  static constexpr int arm(int x) { return x >= TIMES; }
  /** @brief Column of a time index on an arm.
   * @param k Time index.
   * @param arm Arm index.
   * @return Canvas column. */
  static constexpr int column(int k, int arm) { return k + arm * TIMES; }
  /** @brief Surface time of a wrapped column.
   * @param x Wrapped canvas column.
   * @return Frame fraction. */
  static constexpr float time(int x) { return surface_time(x, W); }
  /** @brief Continuous sweep fraction of an unwrapped column.
 * @param column Column coordinate.
 * @return Fraction in [0, 1). */
  HS_HOT_INLINE static float time(float column) {
    const float half_turn = column / (W / 2.0f);
    const int whole = static_cast<int>(half_turn);
    const float fraction = half_turn - static_cast<float>(whole);
    return fraction < 0.0f ? fraction + 1.0f : fraction;
  }

  /** @brief Column nearest a surface time.
   * @param t Frame fraction.
   * @return First-arm column. */
  static constexpr int column_at(float t) { return surface_column(t, W); }
  /** @brief First column at or after a surface time.
   * @param t Frame fraction.
   * @return First-arm column, or TIMES at the frame end. */
  static int first_column_at(float t) {
    return Spacetime::first_column_at(t, W);
  }
};

} // namespace Spacetime
