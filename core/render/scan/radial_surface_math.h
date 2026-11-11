/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file radial_surface_math.h
 * @brief Chart arithmetic and accurate bounded angle evaluation. */

#include <cmath>
#include "core/math/3dmath.h"

namespace Scan {

/** @brief Signed determinant of chart vectors.
 * @tparam Point Chart point type.
 * @param a,b Chart vectors.
 * @return Signed determinant. */
template <class Point>
HS_HOT_INLINE inline float radial_cross(Point a, Point b) {
  return a.x * b.y - a.y * b.x;
}

/** @brief Complex product of chart vectors.
 * @tparam Point Chart point type.
 * @param a,b Complex chart vectors.
 * @return Complex product. */
template <class Point>
HS_HOT_INLINE inline Point radial_multiply(Point a, Point b) {
  return {a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y};
}

/** @brief Floors a chart coordinate to an integer.
 * @param value Coordinate within the integer range.
 * @return Floor of the coordinate. */
HS_HOT_INLINE inline int radial_floor_int(float value) {
  const int truncated = int(value);
  return truncated - (float(truncated) > value);
}

/** @brief Wraps a chart index by a positive count.
 * @param value Unwrapped index.
 * @param count Positive index count.
 * @return Index in [0, count). */
HS_HOT_INLINE inline int radial_wrap_index(int value, int count) {
  return math::wrap_positive(value, count);
}

/** @brief Evaluates atan2 with an octant-reduced degree-9 polynomial.
 * @param y,x Finite Cartesian coordinates.
 * @return Argument in [-pi,pi], preserving signed-zero axis conventions.
 */
HS_HOT_INLINE inline float radial_atan2(float y, float x) {
  const float ax = fabsf(x), ay = fabsf(y);
  const float maximum = fmaxf(ax, ay);
  if (maximum == 0.0f)
    return copysignf(std::signbit(x) ? math::PI_F : 0.0f, y);
  const float ratio = fminf(ax, ay) / maximum;
  const bool shifted = ratio > .4142135623730950488f;
  const float t = shifted ? (ratio - 1.0f) / (ratio + 1.0f) : ratio;
  const float square = t * t;
  const float polynomial =
      .99999990303436148f +
      square * (-.33332184265314274f +
                square * (.19961549539561524f +
                          square * (-.13751577271705973f +
                                    square * .077261777546314622f)));
  float angle = t * polynomial;
  if (shifted)
    angle += .25f * math::PI_F;
  if (ax < ay)
    angle = .5f * math::PI_F - angle;
  if (std::signbit(x))
    angle = math::PI_F - angle;
  return copysignf(angle, y);
}

} // namespace Scan
