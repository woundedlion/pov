/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file periodic.h
 * @brief Scalar periodic wrapping and circular distances.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include "platform/platform.h"
#include <cassert>

namespace math {

/** @brief Wraps into [0, period); requires a positive period. */
template <typename T> inline T wrap_positive(T value, T period) {
  T result = std::fmod(value, period);
  if (result < 0)
    result += period;
  return result >= period ? T{0} : result;
}

/**
 * @brief Wraps a floating-point value around a modulo base (m).
 * @details Requires m > 0 (checked by a debug assert). The effective base is
 *   max(m, numeric_limits<R>::min()), where R is the common result type, so
 *   positive subnormal bases also use the smallest positive normal value.
 *   Under NDEBUG, non-positive or NaN bases use that floor too.
 *   An all-integral call is rejected by static_assert; use the
 *   exact `wrap(int, int)` overload.
 * @tparam T The type of the value being wrapped (e.g., float).
 * @tparam U The type of the modulo base.
 * @param x Finite value to wrap.
 * @param m Positive modulo base, floored to the result type's minimum normal.
 * @return The wrapped value in [0, effective base), as `std::common_type_t<T, U>`
 *   (preserves the float math, so `wrap(3, 2.5f)` yields 0.5f).
 */
template <typename T, typename U>
inline std::common_type_t<T, U> wrap(T x, U m) {
  static_assert(!(std::is_integral_v<T> && std::is_integral_v<U>),
                "wrap(integral, integral) loses precision past 2^53 via the "
                "double fmod path; use the exact wrap(int, int) overload");
  using R = std::common_type_t<T, U>;
  assert(m > 0);
  // fmax(NaN, y) == y also floors a NaN base.
  const R mm = std::fmax(static_cast<R>(m), std::numeric_limits<R>::min());
  return wrap_positive(static_cast<R>(x), mm);
}

/**
 * @brief Fast floating point modulo for 1.0.
 * @param t Finite value to wrap.
 * @return The wrapped value in the range [0.0, 1.0).
 * @details Wraps finite floats into [0.0, 1.0). For tiny negative t in (-2.98e-8, 0),
 *   `t - floorf(t)` rounds up to exactly 1.0f, violating the half-open contract;
 *   the guard folds that boundary back to 0.
 */
__attribute__((always_inline)) inline float wrap_t(float t) {
  float r = t - floorf(t);
  return r >= 1.0f ? 0.0f : r;
}

/**
 * @brief Wraps an integer value around a modulo base (m).
 * @param x The value to wrap.
 * @param m The modulo base.
 * @return The wrapped value in the range [0, m), or x when m <= 0.
 */
inline int wrap(int x, int m) {
  if (m <= 0)
    return x;
  int r = x % m;
  return r < 0 ? r + m : r;
}

/**
 * @brief Single-step integer wrap into [0, W), assuming x is at most one period
 * out of range (i.e. x is in [-W, 2*W)).
 * @param x The value to wrap; must lie in [-W, 2*W).
 * @param W The modulo base (period length).
 * @return The wrapped value in the range [0, W).
 */
inline int fast_wrap(int x, int W) {
  assert(x >= -W && x < 2 * W);
  if (x >= W)
    return x - W;
  if (x < 0)
    return x + W;
  return x;
}

/**
 * @brief Single-step wrap of a sub-pixel column into [0, W), assuming x is at
 * most one period out of range (i.e. x is in [-W, 2*W)).
 * @param x The value to wrap; must lie in [-W, 2*W).
 * @param W The modulo base (period length).
 * @return The wrapped value in the range [0, W).
 * @details For tiny negative x, `x + period` rounds up to exactly `period`,
 *   violating the half-open contract; the guard folds that boundary back to 0.
 */
inline float fast_wrap(float x, int W) {
  const float period = static_cast<float>(W);
  assert(x >= -period && x < 2.0f * period);
  if (x >= period)
    return x - period;
  if (x < 0.0f) {
    const float r = x + period;
    return r >= period ? 0.0f : r;
  }
  return x;
}

/**
 * @brief Calculates the shortest distance (either forwards or backwards)
 * between two points on a circular domain.
 * @param a The first position.
 * @param b The second position.
 * @param m Positive domain length, floored to the minimum normal float.
 * @return The shortest distance in [0, effective domain length / 2].
 */
inline float shortest_distance(float a, float b, float m) {
  assert(m > 0.0f);
  m = std::fmax(m, std::numeric_limits<float>::min());
  // Double-fmod for full range reduction of an arbitrary a - b into [0, m).
  float d = std::fmod(std::fmod(a - b, m) + m, m);
  return std::min(d, m - d);
}

} // namespace math
