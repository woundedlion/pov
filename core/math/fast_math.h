/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file fast_math.h @brief Scalar approximations and Cartesian-to-spherical conversion. */

#include "core/math/vector.h"

namespace math {

/** @brief Absolute angular error bound in radians for unit-length inputs. */
inline constexpr double FAST_ATAN2_MAX_ERROR = 5e-3;

/**
 * @brief Fast atan2 using the 0.273-Hastings polynomial.
 * @param y Y (numerator) coordinate.
 * @param x X (denominator) coordinate.
 * @return The angle in radians in [-π, π]. Matches std::atan2 except at
 * negative zero: `-0.0f < 0.0f` is false, so `x<0, y==-0.0f` returns +π where
 * std::atan2 gives -π; -π is reached only for strictly negative y.
 * At the origin it returns π/2.
 * @details Peak abs error ~0.0038 rad (~0.22°). The origin nudge is added to
 * |y| unconditionally, so at tiny |(x, y)| the result is not scale-invariant
 * and the error grows toward a full quadrant.
 */
__attribute__((always_inline)) inline float fast_atan2(float y, float x) {
  // +1e-10f keeps abs_y strictly positive so the (0,0) origin stays finite.
  float abs_y = std::abs(y) + 1e-10f;
  float abs_x = std::abs(x);
  float r, angle;

  if (abs_x < abs_y) {
    r = abs_x / abs_y;
    angle = 1.57079633f - r * (0.78539816f + 0.273f * (1.0f - r));
  } else {
    r = abs_y / abs_x;
    angle = r * (0.78539816f + 0.273f * (1.0f - r));
  }

  if (x < 0.0f)
    angle = 3.14159265f - angle;
  if (y < 0.0f)
    angle = -angle;

  return angle;
}

/** @brief Azimuth in half-turns.
 * @param p Point on the sphere.
 * @return Signed azimuth divided by pi. */
HS_HOT_INLINE inline float half_turns(const Vector &p) {
  return fast_atan2(p.z, p.x) / PI_F;
}

/**
 * @brief Arctangent of a ratio in [0, 1], accurate to about 1e-5 rad.
 * @param r Ratio in [0, 1].
 * @return atan(r) in [0, pi/4].
 * @details Odd degree-11 minimax polynomial; precise_atan2 folds the other
 * octants onto it.
 */
__attribute__((always_inline)) inline float atan_unit(float r) {
  const float r2 = r * r;
  return r * (0.99997726f +
              r2 * (-0.33262347f +
                    r2 * (0.19354346f +
                          r2 * (-0.11643287f +
                                r2 * (0.05265332f + r2 * -0.01172120f)))));
}

/**
 * @brief atan2 accurate to about 1e-5 rad.
 * @param y Y (numerator) coordinate.
 * @param x X (denominator) coordinate.
 * @return The angle of (x, y) in [-pi, pi]; 0 at the origin.
 * @details One divide and an odd degree-11 minimax arctangent on [0, 1].
 */
__attribute__((always_inline)) inline float precise_atan2(float y, float x) {
  const float abs_y = std::abs(y);
  const float abs_x = std::abs(x);
  const float most = fmaxf(abs_x, abs_y);
  if (!(most > 0.0f))
    return 0.0f;
  float angle = atan_unit(fminf(abs_x, abs_y) / most);
  if (abs_y > abs_x)
    angle = 1.57079633f - angle;
  if (x < 0.0f)
    angle = 3.14159265f - angle;
  return y < 0.0f ? -angle : angle;
}

/**
 * @brief Diamond pseudo-angle of (x, y), monotonic with atan2.
 * @tparam FoldSeam Fold a rounded 4.0 endpoint onto zero when true.
 * @param y Y (numerator) coordinate.
 * @param x X (denominator) coordinate.
 * @return The diamond angle in [0, 4), or [0, 4] when FoldSeam is false;
 *         0 when |x| + |y| < 1e-20.
 * @details Arc position along the diamond |x| + |y| = 1, walked
 * counter-clockwise from +x. Strictly monotonic with atan2, so it bins a
 * direction without a trig call. The near-origin fallback is not scale invariant.
 */
template <bool FoldSeam = true>
HS_HOT_INLINE constexpr float diamond_angle(float y, float x) {
  float d = __builtin_fabsf(x) + __builtin_fabsf(y);
  if (d < 1e-20f)
    return 0.0f;
  float r = y / d;
  if (y >= 0.0f)
    return (x >= 0.0f) ? r : (2.0f - r);
  if (x < 0.0f)
    return 2.0f - r;
  // A tiny negative r can round 4 + r to exactly 4.
  const float a = 4.0f + r;
  if constexpr (FoldSeam)
    return a >= 4.0f ? 0.0f : a;
  return a;
}

/**
 * @brief Fast reciprocal for positive normal inputs up to 2^125.
 * @param x Positive normal float no greater than 2^125 (about 4.25e37).
 * @return An approximation of 1 / x.
 * @details Float-bit seed refined by two Newton steps; peak relative error is
 * below 7e-6 over the stated domain, where the seed estimate remains normal.
 */
HS_O3_FN __attribute__((always_inline)) inline float fast_reciprocal(float x) {
  uint32_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  bits = 0x7ef311c3u - bits;
  float estimate;
  std::memcpy(&estimate, &bits, sizeof(estimate));
  estimate *= 2.0f - x * estimate;
  estimate *= 2.0f - x * estimate;
  return estimate;
}

/**
 * @brief Fast reciprocal square root for positive normal floats.
 * @tparam NEWTON_STEPS Newton refinements, one or two.
 * @param x Positive normal float; subnormal inputs lose accuracy.
 * @return An approximation of 1 / sqrt(x).
 * @details Bit-hack initial guess refined by the selected Newton steps;
 * the default two-step result has peak relative error ~5e-6 and is one-sided
 * low up to the rounding of the final multiply.
 */
template <int NEWTON_STEPS = 2> HS_O3_FN inline float fast_rsqrt(float x) {
  static_assert(NEWTON_STEPS == 1 || NEWTON_STEPS == 2);
  uint32_t i;
  std::memcpy(&i, &x, sizeof(i));
  i = 0x5f3759dfu - (i >> 1);
  float y;
  std::memcpy(&y, &i, sizeof(y));
  const float half = 0.5f * x;
  y = y * (1.5f - half * y * y);
  if constexpr (NEWTON_STEPS == 2)
    y = y * (1.5f - half * y * y);
  return y;
}

/**
 * @brief One Halley cube-root step split into numerator and denominator.
 * @param x Input value; x <= 0 yields num = 0, den = 1, leaving a shared
 * product unchanged.
 * @param num Receives the Halley numerator; cbrt(x) ~= num / den.
 * @param den Receives the Halley denominator.
 */
HS_O3_FN inline void cbrt_halley_terms(float x, float &num, float &den) {
  if (x <= 0.0f) {
    num = 0.0f;
    den = 1.0f;
    return;
  }
  uint32_t i;
  std::memcpy(&i, &x, sizeof(i));
  i = i / 3u + 0x2a514067u;
  float y;
  std::memcpy(&y, &i, sizeof(y));
  float c = y * y * y;
  num = y * (c + 2.0f * x);
  den = 2.0f * c + x;
}

/**
 * @brief Fast cube root for x in [~1e-29, ~3e28].
 * @param x Input value; x <= 0 returns 0.
 * @return An approximation of the cube root of `x`.
 * @details Bit-hack initial guess refined by one Halley step. Relative error
 * is about 2.3e-5 over that range; results underflow to 0 below ~1e-34.
 */
HS_O3_FN inline float fast_cbrt(float x) {
  float numerator, denominator;
  cbrt_halley_terms(x, numerator, denominator);
  return numerator / denominator;
}

/**
 * @brief Fast cube roots of three values sharing a single division.
 * @param x1 First input; same domain as fast_cbrt (x <= 0 yields 0).
 * @param x2 Second input.
 * @param x3 Third input.
 * @param o1 Cube root of `x1`.
 * @param o2 Cube root of `x2`.
 * @param o3 Cube root of `x3`.
 * @details Uses fast_cbrt's seed and Halley step. The shared denominator
 * product can overflow or underflow, including for unequal inputs. For equal
 * inputs, overflow begins near 1.3e11 and accuracy degrades below about 3e-13.
 */
HS_O3_FN inline void fast_cbrt3(float x1, float x2, float x3, float &o1,
                                float &o2, float &o3) {
  float n1, d1, n2, d2, n3, d3;
  cbrt_halley_terms(x1, n1, d1);
  cbrt_halley_terms(x2, n2, d2);
  cbrt_halley_terms(x3, n3, d3);
  const float inv = 1.0f / (d1 * d2 * d3);
  o1 = n1 * (d2 * d3) * inv;
  o2 = n2 * (d1 * d3) * inv;
  o3 = n3 * (d1 * d2) * inv;
}

/**
 * @brief Fast cube roots of six values sharing a single division.
 * @param x Six inputs; same domain as fast_cbrt (x <= 0 yields 0).
 * @param o Six outputs, the cube roots of the corresponding inputs.
 * @details Uses fast_cbrt's seed and Halley step. The shared denominator
 * product can overflow or underflow, including for unequal inputs. For equal
 * inputs, overflow begins near 4.2e5.
 */
HS_O3_FN inline void fast_cbrt6(const float x[6], float o[6]) {
  float n[6], d[6];
  for (int i = 0; i < 6; ++i)
    cbrt_halley_terms(x[i], n[i], d[i]);
  // pre[i] * suf[i] is the product of every denominator but d[i].
  float pre[6], suf[6];
  pre[0] = 1.0f;
  for (int i = 1; i < 6; ++i)
    pre[i] = pre[i - 1] * d[i - 1];
  suf[5] = 1.0f;
  for (int i = 4; i >= 0; --i)
    suf[i] = suf[i + 1] * d[i + 1];
  const float inv = 1.0f / (pre[5] * d[5]);
  for (int i = 0; i < 6; ++i)
    o[i] = n[i] * (pre[i] * suf[i]) * inv;
}

// Forward declarations
inline float fast_acos(float x);

inline Spherical::Spherical(const Vector &v) {
  Vector n(v);
  n.normalize();
  // theta is in [-π, π], not [0, 2π), and arbitrary at the poles.
  theta = fast_atan2(n.z, n.x);
  phi = fast_acos(hs::clamp(n.y, -1.0f, 1.0f));
}

} // namespace math
