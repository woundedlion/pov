/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file scalar_noise.h @brief Scalar constants, interpolation kernels, and value noise. */

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <algorithm>
#include "core/platform/platform.h"

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

/** Exact 2^96 lift for divisors whose squared magnitude underflows the normal range. */
inline constexpr float COMPLEX_UNDERFLOW_LIFT = 0x1p96f;

/**
 * @brief The Golden Ratio constant (Phi).
 */
inline constexpr float PHI = 1.6180339887498948482045868343656f;
/**
 * @brief The inverse of the Golden Ratio (1/Phi).
 */
inline constexpr float INV_PHI = 1 / PHI;
/** @brief Generic floating-point comparison tolerance. */
inline constexpr float TOLERANCE = 1e-4f;
/// Geometric near-equality (positions, angles).
inline constexpr float EPS_GEOMETRIC = 1e-5f;
/// Squared edge length below which an edge is degenerate.
inline constexpr float EPS_LEN_SQ = 1e-6f;
/**
 * @brief Squared |cross(from, to)| below which make_rotation synthesizes a π
 *        turn instead of using the half-vector form.
 * @details Placed where the half-vector form's angular error
 * (~2*ULP(1)/|cross|) meets the synthesized turn's (|cross|).
 */
inline constexpr float EPS_ANTIPARALLEL_SQ = 4e-7f;
/// Squared cross-product magnitude below which it is degenerate.
inline constexpr float EPS_CROSS_SQ = 1e-8f;
/// Squared face-normal magnitude below which the face is degenerate.
inline constexpr float EPS_NORMAL_SQ = 1e-9f;
/// Squared length below which normalize() has no reliable direction.
inline constexpr float EPS_NORMALIZE_SQ = 1e-12f;
/** @brief Squared length below which a blend of two directions has cancelled. */
inline constexpr float EPS_BLEND_LEN_SQ = 1e-8f;
/// Slack on |q|^2 - 1 for is-unit quaternion assertions.
inline constexpr float EPS_UNIT_QUAT_SQ = 0.01f;
/// Slack on |v|^2 - 1 for is-unit vector assertions.
inline constexpr float EPS_UNIT_VEC_SQ = 0.02f;
/**
 * @brief Cosine above which a vector is treated as parallel to a reference axis.
 */
inline constexpr float COS_AXIS_PARALLEL = 1.0f - TOLERANCE;
/**
 * @brief Floating-point representation of PI.
 */
inline constexpr float PI_F = static_cast<float>(PI);
/**
 * @brief A full turn in radians (2π).
 * @details Suffixed _F because Arduino's wiring.h defines TWO_PI as a macro.
 */
inline constexpr float TWO_PI_F = 2.0f * PI_F;

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

/**
 * @brief Stable hash of an index and seed to a float in [0, 1).
 * @param i Index to hash.
 * @param seed Stream selector; different seeds give independent hashes.
 * @return The hashed value in [0, 1).
 * @details Stateless PCG output hash (RXS-M-XS). The seed enters both the
 * index and the output state; folding it into the index alone would only
 * permute the lattice.
 */
inline float hash01(uint32_t i, uint32_t seed) {
  uint32_t h = (i ^ seed) * 747796405u + 2891336453u;
  h = ((h >> ((h >> 28) + 4u)) ^ h ^ seed) * 277803737u;
  h = (h >> 22) ^ h;
  return static_cast<float>(h >> 8) * (1.0f / 16777216.0f);
}

/**
 * @brief 1D value noise: quintic-smoothed hash lattice.
 * @param x Sample coordinate; lattice points sit on the integers.
 * @pre floorf(x) is strictly between -2^31 and 2^31, checked only in audit builds.
 * @param seed Stream selector forwarded to the lattice hash.
 * @return Noise value in [0, 1), C1-continuous in x.
 */
inline float value_noise_1d(float x, uint32_t seed = 0) {
  float fx = floorf(x);
  // int32_t first: float -> unsigned is UB for negatives. Direct compare so a
  // NaN coordinate fails it too.
  HS_AUDIT_CHECK(std::fabs(fx) < 2147483648.0f);
  uint32_t ix = static_cast<uint32_t>(static_cast<int32_t>(fx));
  float f = quintic_kernel(x - fx);
  float a = hash01(ix, seed);
  float b = hash01(ix + 1u, seed);
  return a + (b - a) * f;
}

/**
 * @brief 2D value noise: bilinear quintic-smoothed hash lattice.
 * @param x First sample coordinate; lattice points sit on the integers.
 * @param y Second sample coordinate.
 * @pre floorf(x) and floorf(y) are strictly between -2^31 and 2^31, checked only
 *      in audit builds.
 * @param seed Stream selector forwarded to the lattice hash.
 * @return Noise value in [0, 1), C1-continuous in x and y.
 */
inline float value_noise_2d(float x, float y, uint32_t seed = 0) {
  float fx = floorf(x), fy = floorf(y);
  HS_AUDIT_CHECK(std::fabs(fx) < 2147483648.0f);
  HS_AUDIT_CHECK(std::fabs(fy) < 2147483648.0f);
  uint32_t ix = static_cast<uint32_t>(static_cast<int32_t>(fx));
  uint32_t iy = static_cast<uint32_t>(static_cast<int32_t>(fy));
  float u = quintic_kernel(x - fx);
  float v = quintic_kernel(y - fy);
  // Row stride is a large odd constant so adjacent rows decorrelate.
  uint32_t row0 = iy * 0x9E3779B9u, row1 = (iy + 1u) * 0x9E3779B9u;
  float a = hash01(ix + row0, seed);
  float b = hash01(ix + 1u + row0, seed);
  float c = hash01(ix + row1, seed);
  float d = hash01(ix + 1u + row1, seed);
  float ab = a + (b - a) * u;
  float cd = c + (d - c) * u;
  return ab + (cd - ab) * v;
}

} // namespace math
