/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// platform.h first: on device it defines NDEBUG, which must be set before
// <cassert> expands the assert macro — otherwise assert-stripping would depend
// on a prior TU having pulled in platform.h, making this header non-self-sufficient.
#include "platform.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>
#include "3dmath.h"
#include <memory>
#include <utility>

/**
 * @brief Wraps a floating-point value around a modulo base (m).
 * @details Ensures the result is always non-negative and in [0, m).
 *   Overload resolution: this `fmod`-based template is selected for any call
 *   that is not an exact `(int, int)` match; the non-template `wrap(int, int)`
 *   below is preferred only for two `int` arguments. A mixed call such as
 *   `wrap(int_x, float_m)` therefore resolves here (float math), not to the
 *   integer overload. Precondition: m > 0; the debug `assert(m > 0)` enforces it.
 *   Under NDEBUG (assert stripped) the `fmax(m, min())` floor below clamps a
 *   non-positive m to the smallest positive value, so the result collapses to a
 *   near-zero in [0, min()) rather than NaN/garbage.
 *   Exactness: the body wraps via `std::fmod` (floating-point), so the result is
 *   exact only within the working mantissa — ~2^24 for a `float` common type,
 *   ~2^53 for the `double` path. An all-integral call (e.g. `wrap(long, long)`)
 *   would silently lose precision past 2^53, so a `static_assert` rejects it at
 *   compile time; use the dedicated `wrap(int, int)` overload below, which uses
 *   integer `%` and is exact across its full range.
 * @tparam T The type of the value being wrapped (e.g., float).
 * @tparam U The type of the modulo base.
 * @param x The value to wrap.
 * @param m The modulo base.
 * @return The wrapped value in the range [0, m), as `std::common_type_t<T, U>`.
 *   Returning the common type rather than T preserves the float math, so
 *   `wrap(3, 2.5f)` yields 0.5f rather than a truncated 0.
 */
template <typename T, typename U>
inline std::common_type_t<T, U> wrap(T x, U m) {
  static_assert(!(std::is_integral_v<T> && std::is_integral_v<U>),
                "wrap(integral, integral) loses precision past 2^53 via the "
                "double fmod path; use the exact wrap(int, int) overload");
  using R = std::common_type_t<T, U>;
  assert(m > 0);
  // Branchless floor to the smallest positive value (the hot SDF angular-repeat
  // path forbids a branch); fmax(NaN, y) == y also blocks a NaN base.
  const R mm = std::fmax(static_cast<R>(m), std::numeric_limits<R>::min());
  R r = std::fmod(static_cast<R>(x), mm);
  if (r < 0) {
    r += mm;
  }
  return (r >= mm) ? R{0} : r;
}

/**
 * @brief Fast floating point modulo for 1.0.
 * @param t The value to wrap.
 * @return The wrapped value in the range [0.0, 1.0).
 * @details Safely wraps any float (positive or negative) into the [0.0, 1.0)
 *   range. For tiny negative t in (-2.98e-8, 0), `t - floorf(t)` rounds up to
 *   exactly 1.0f, violating the half-open [0,1) contract; the guard folds that
 *   boundary back to 0, mirroring the sibling `wrap()` template above. The
 *   compare lowers to a branchless select, so the per-call cost is negligible.
 */
inline float wrap_t(float t) {
  float r = t - floorf(t);
  return r >= 1.0f ? 0.0f : r;
}

/**
 * @brief Wraps an integer value around a modulo base (m).
 * @details Ensures the result is always non-negative.
 * @param x The value to wrap.
 * @param m The modulo base.
 * @return The wrapped value in the range [0, m).
 */
inline int wrap(int x, int m) {
  assert(m > 0);
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
 * @brief Calculates the shortest distance (either forwards or backwards)
 * between two points on a circular domain.
 * @param a The first position.
 * @param b The second position.
 * @param m The modulo base (length of the domain).
 * @return The shortest distance in the range [0, m/2].
 */
inline float shortest_distance(float a, float b, float m) {
  assert(m > 0.0f);
  // Double-fmod for full range reduction of an arbitrary a - b into [0, m).
  float d = std::fmod(std::fmod(a - b, m) + m, m);
  return std::min(d, m - d);
}

/**
 * @brief Invokes `apply(current)` only when `current` differs from `last`,
 * then latches `last = current`.
 * @details Encodes the "live-apply a slider value, but only on change" idiom:
 *   several effects retain an animation handle whose setter reschedules from
 *   "now" (e.g. `set_duration`/`set_period`), so calling it every frame would
 *   perpetually defer the trigger — it must fire only when the parameter
 *   actually moves. The comparison/branch is identical to the hand-rolled form
 *   and `apply` is invoked at most once per change, so this is a zero-overhead
 *   inline.
 * @param current The latest parameter value to test against the cached one.
 * @param last The cached value; updated to `current` when they differ.
 * @param apply Callable receiving the new value; run only on change.
 */
template <typename T, typename F>
inline void apply_if_changed(const T &current, T &last, F &&apply) {
  // T must compare with exact equality: a tolerance comparator (e.g. Vector's
  // operator!=) is non-transitive and could re-fire or latch incorrectly.
  static_assert(std::is_arithmetic_v<T>,
                "apply_if_changed requires an exact-equality T (scalar/int)");
  if (current != last) {
    last = current;
    apply(current);
  }
}

