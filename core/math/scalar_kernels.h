/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file scalar_kernels.h @brief Scalar curves and interpolation kernels. */

#include <cmath>
#include "core/math/math_constants.h"

namespace math {

/**
 * @brief Unit cup curve over [0, 1]: 1 at both ends, 0 at the midpoint.
 * @param t Parameter in [0, 1].
 * @return |2t - 1|.
 */
__attribute__((always_inline)) inline float unit_cup(float t) {
  return fabsf(2.0f * t - 1.0f);
}

/**
 * @brief Unit bell curve over [0, 1]: 0 at both ends, 1 at the midpoint.
 * @param t Parameter in [0, 1].
 * @return 1 - |2t - 1|.
 */
__attribute__((always_inline)) inline float unit_bell(float t) {
  return 1.0f - unit_cup(t);
}

/**
 * @brief Quintic interpolation kernel (smootherstep).
 * @param t Interpolation factor (clamped to 0.0 - 1.0).
 * @return The interpolated value: t^3 * (t * (t * 6 - 15) + 10).
 */
__attribute__((always_inline)) inline float quintic_kernel(float t) {
  t = hs::clamp(t, 0.0f, 1.0f);
  float t2 = t * t;
#if defined(CORE_TEENSY) && defined(__GNUC__)
  // Keep GCC -Os from replacing the cube with __powisf2.
  __asm__("" : "+t"(t2));
#endif
  return (t * t2) * (t * (t * 6.0f - 15.0f) + 10.0f);
}

/**
 * @brief Cubic interpolation kernel (smoothstep).
 * @param t Interpolation factor (clamped to 0.0 - 1.0).
 * @return The interpolated value: t^2 * (3 - 2t).
 */
__attribute__((always_inline)) inline float cubic_kernel(float t) {
  t = hs::clamp(t, 0.0f, 1.0f);
  return t * t * (3.0f - 2.0f * t);
}

/**
 * @brief Smoothstep between two edges.
 * @param edge0 Value mapped to 0.
 * @param edge1 Value mapped to 1; must differ from edge0.
 * @param value Argument to ramp.
 * @return cubic_kernel of the normalized position, clamped to 0.0 - 1.0.
 */
HS_O3_FN HS_HOT_INLINE inline float smooth_ramp(float edge0, float edge1,
                                                float value) {
  return cubic_kernel((value - edge0) / (edge1 - edge0));
}

/**
 * @brief Taylor sine through degree nine, folded about pi/2.
 * @param x Angle in [0, pi], in radians.
 * @return Approximate sine of x.
 */
constexpr float constexpr_sin(float x) {
  if (x > PI_F * 0.5f)
    x = PI_F - x;
  const float x2 = x * x;
  return x * (1.0f + x2 * (-1.0f / 6.0f +
                           x2 * (1.0f / 120.0f +
                                 x2 * (-1.0f / 5040.0f + x2 / 362880.0f))));
}

/**
 * @brief Compile-time square root by Newton-Raphson.
 * @param x Radicand.
 * @return sqrt(x) to within one ulp, or 0 for a non-positive @p x.
 * @details Seeded at or above sqrt(x), the iteration descends monotonically
 *          to its fixed point. For run-time use call sqrtf.
 */
constexpr float constexpr_sqrt(float x) {
  if (x <= 0.0f)
    return 0.0f;
  float root = 0.5f * (1.0f + x);
  for (int i = 0; i < 128; ++i) {
    const float next = 0.5f * (root + x / root);
    if (next >= root)
      break;
    root = next;
  }
  return root;
}

} // namespace math
