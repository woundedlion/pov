/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "3dmath.h"

/**
 * @brief Namespace for custom utility and random number generation helpers.
 */
namespace hs {
  /**
   * @brief Generates a pseudo-random floating-point number between 0.0 and 1.0.
   * @return A random float in the range [0.0, 1.0].
   */
  float rand_f() {
    return static_cast<float>(::random(0, std::numeric_limits<int32_t>::max()))
      / std::numeric_limits<int32_t>::max();
  }

  /**
   * @brief Generates a pseudo-random integer within a specified range.
   * @param min The minimum value (inclusive).
   * @param max The maximum value (exclusive).
   * @return A random integer in the range [min, max).
   */
  int rand_int(int min, int max) {
    return ::random(min, max);
  }
}

/**
 * @brief Wraps a floating-point value around a modulo base (m).
 * @details Ensures the result is always non-negative and correctly handles small floating-point zero values.
 * @tparam T The type of the value being wrapped (e.g., float).
 * @tparam U The type of the modulo base.
 * @param x The value to wrap.
 * @param m The modulo base.
 * @return The wrapped value in the range [0, m).
 */
template <typename T, typename U>
inline T wrap(T x, U m) {
  if (std::abs(x) < TOLERANCE) {
    return 0;
  }

  T r = std::fmod(x, m);
  if (r < 0) {
    r += m;
  }
  return (r >= m) ? 0 : r;
}

/**
 * @brief Calculates the shortest distance (either forwards or backwards) between two points on a circular domain.
 * @param a The first position.
 * @param b The second position.
 * @param m The modulo base (length of the domain).
 * @return The shortest distance in the range [0, m/2].
 */
float shortest_distance(float a, float b, float m) {
  float d = std::fmod(std::fmod(a - b, m) + m, m);
  return std::min(d, m - d);
}

/**
 * @brief Calculates the forward (positive direction) distance from point a to point b on a circular domain.
 * @param a The starting position.
 * @param b The ending position.
 * @param m The modulo base (length of the domain).
 * @return The forward distance in the range [0, m).
 */
constexpr float fwd_distance(float a, float b, float m) {
  auto d = b - a;
  if (d < 0) {
    d += m;
  }
  return d;
}
