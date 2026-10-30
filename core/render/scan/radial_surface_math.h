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
 * @tparam X_FIRST Add the x-led imaginary product first.
 * @tparam Point Chart point type.
 * @param a,b Complex chart vectors.
 * @return Complex product. */
template <bool X_FIRST = false, class Point>
HS_HOT_INLINE inline Point radial_multiply(Point a, Point b) {
  if constexpr (X_FIRST)
    return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
  else
    return {a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y};
}

/** @brief Floors a chart coordinate to an integer.
 * @param value Coordinate within the integer range.
 * @return Floor of the coordinate. */
HS_HOT_INLINE inline int radial_floor_int(float value) {
  const int TRUNCATED = int(value);
  return TRUNCATED - (float(TRUNCATED) > value);
}

/** @brief Wraps a chart index by a positive count.
 * @param value Unwrapped index.
 * @param count Positive index count.
 * @return Index in [0, count). */
HS_HOT_INLINE inline int radial_wrap_index(int value, int count) {
  const int REMAINDER = value % count;
  return REMAINDER < 0 ? REMAINDER + count : REMAINDER;
}

/** @brief Evaluates atan2 with an octant-reduced degree-9 polynomial.
 * @param y,x Finite Cartesian coordinates.
 * @return Argument in [-pi,pi], preserving signed-zero axis conventions.
 * @details Chebyshev truncation contributes less than 3.66e-9 radians;
 * floating-point range reduction and reconstruction add rounding error. */
HS_HOT_INLINE inline float radial_atan2(float y, float x) {
  const float AX = fabsf(x), AY = fabsf(y);
  const float MAXIMUM = fmaxf(AX, AY);
  if (MAXIMUM == 0.0f)
    return copysignf(std::signbit(x) ? math::PI_F : 0.0f, y);
  const float RATIO = fminf(AX, AY) / MAXIMUM;
  const bool SHIFTED = RATIO > .4142135623730950488f;
  const float T = SHIFTED ? (RATIO - 1.0f) / (RATIO + 1.0f) : RATIO;
  const float SQUARE = T * T;
  const float POLYNOMIAL =
      .99999990303436148f +
      SQUARE * (-.33332184265314274f +
                SQUARE * (.19961549539561524f +
                          SQUARE * (-.13751577271705973f +
                                    SQUARE * .077261777546314622f)));
  float angle = T * POLYNOMIAL;
  if (SHIFTED)
    angle += .25f * math::PI_F;
  if (AX < AY)
    angle = .5f * math::PI_F - angle;
  if (std::signbit(x))
    angle = math::PI_F - angle;
  return copysignf(angle, y);
}

} // namespace Scan
