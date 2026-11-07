/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file interpolate.h
 * @brief Interpolators per parameter domain — unconstrained scalar, positive
 *        scale, periodic angle, unit vector — plus a progress clamp. Each is
 *        exact at both endpoints and snaps to the nearer one outside [0,1].
 */

#include "core/math/periodic.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace interp {

/**
 * @brief Clamps a progress fraction into [0,1].
 * @param value Unclamped fraction.
 * @return @p value confined to [0,1]. Neither comparison holds for NaN, which
 *         therefore passes through unchanged.
 */
constexpr float clamp_progress(float value) {
  if (value <= 0.0f)
    return 0.0f;
  if (value >= 1.0f)
    return 1.0f;
  return value;
}

/**
 * @brief Linear interpolation with exact endpoints.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param progress Fraction; outside [0,1] it snaps to the nearer endpoint.
 * @return The interpolated value.
 */
constexpr float linear(float a, float b, float progress) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
  return a + (b - a) * progress;
}

/**
 * @brief Geometric (constant-ratio) interpolation for a positive scale.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param progress Fraction; outside [0,1] it snaps to the nearer endpoint.
 * @return The interpolated value, or the linear interpolation when either
 *         endpoint is non-positive and the logarithm is undefined.
 */
inline float log_positive(float a, float b, float progress) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
  if (a <= 0.0f || b <= 0.0f)
    return linear(a, b, progress);
  return std::exp(linear(std::log(a), std::log(b), progress));
}

namespace detail {

/** @brief Folds a displacement into [-period / 2, period / 2).
 * @param raw_delta Unwrapped displacement.
 * @param period Positive domain period.
 * @return Folded displacement. */
inline float shortest_periodic_delta(float raw_delta, float period) {
  const float half = period * 0.5f;
  float delta =
      std::isfinite(period) && raw_delta > -period && raw_delta < period
          ? raw_delta
          : std::fmod(raw_delta, period);
  if (delta < -half)
    delta += period;
  if (delta >= half)
    delta -= period;
  return delta;
}

/** @brief Interior interpolation with binary32 rounding at each step.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param progress Interior fraction in (0, 1).
 * @param period Positive domain period.
 * @param minimum Domain origin.
 * @return Interpolated value in [minimum, minimum + period]. */
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-associative-math", "fp-contract=off")))
#endif
inline float shortest_periodic_in_range(float a, float b, float progress,
                                        float period, float minimum) {
#if defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#endif
  const float delta = shortest_periodic_delta(b - a, period);
  const float VALUE = a + delta * progress;
  if (VALUE >= minimum && VALUE <= minimum + period)
    return VALUE;
  return minimum + math::wrap_positive(VALUE - minimum, period);
}

} // namespace detail

/**
 * @brief Interpolation along the shorter arc of a periodic domain.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param progress Fraction; outside [0,1] it snaps to the nearer endpoint.
 * @param period Domain period (e.g. TWO_PI_F for an angle).
 * @return An interior progress gives a value in [0, period); an endpoint
 *         gives that endpoint as supplied, unwrapped. Falls back to the
 *         linear interpolation for a non-positive period.
 */
inline float shortest_periodic(float a, float b, float progress, float period) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
  if (period <= 0.0f)
    return linear(a, b, progress);
  const float delta = detail::shortest_periodic_delta(b - a, period);
  return math::wrap_positive(a + delta * progress, period);
}

/**
 * @brief Shorter-arc interpolation in a domain beginning at @p minimum.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param progress Fraction; outside [0, 1] it snaps to the nearer endpoint.
 * @param period Domain period; non-positive values use linear interpolation.
 * @param minimum Domain origin.
 * @return An interior value in [minimum, minimum + period] for a nonzero
 *         minimum, with the same endpoint and fallback rules as the
 *         zero-origin overload.
 */
inline float shortest_periodic(float a, float b, float progress, float period,
                               float minimum) {
  if (minimum == 0.0f || progress <= 0.0f || progress >= 1.0f || period <= 0.0f)
    return shortest_periodic(a, b, progress, period);
  return detail::shortest_periodic_in_range(a, b, progress, period, minimum);
}

/**
 * @brief Outcome of normalized_linear().
 * @tparam Size Component count of the interpolated vector.
 */
template <size_t Size> struct NormalizedLinearResult {
  std::array<float, Size> values; /**< Interpolated vector; unit-length only
                                       for 0 < progress < 1 with valid true —
                                       a snapped endpoint is returned as
                                       supplied. */
  bool valid; /**< False when an interior interpolant's norm is at or below
                  epsilon. */
};

/**
 * @brief Linear interpolation of a vector, renormalized to unit length.
 * @tparam Size Component count.
 * @param a Vector at progress 0.
 * @param b Vector at progress 1.
 * @param progress Fraction; outside [0,1] it snaps to the nearer endpoint
 *        (returned as supplied, without renormalizing).
 * @param epsilon Norm below which the interpolant has no direction.
 * @pre The interior interpolant's squared norm is finite; @p epsilon is finite
 *      and non-negative.
 * @return The unit-length interpolant with valid true, or the un-normalized
 *         interpolant with valid false when its norm is at or below @p epsilon.
 */
template <size_t Size>
inline NormalizedLinearResult<Size>
normalized_linear(const std::array<float, Size> &a,
                  const std::array<float, Size> &b, float progress,
                  float epsilon) {
  if (progress <= 0.0f)
    return {a, true};
  if (progress >= 1.0f)
    return {b, true};
  std::array<float, Size> values{};
  float norm_squared = 0.0f;
  for (size_t index = 0; index < Size; ++index) {
    values[index] = linear(a[index], b[index], progress);
    norm_squared += values[index] * values[index];
  }
  if (norm_squared <= epsilon * epsilon)
    return {values, false};
  const float inverse_norm = 1.0f / std::sqrt(norm_squared);
  for (float &value : values)
    value *= inverse_norm;
  return {values, true};
}

} // namespace interp
