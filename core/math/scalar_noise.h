/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#ifndef HS_MATH_SCALAR_NOISE_H
#define HS_MATH_SCALAR_NOISE_H

/** @file scalar_noise.h @brief Scalar hash and value noise. */

#include <cmath>
#include <cstdint>
#include "core/math/math_constants.h"
#include "core/math/scalar_kernels.h"

namespace math {

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

#endif
