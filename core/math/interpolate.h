/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file interpolate.h
 * @brief One interpolator per parameter domain — unconstrained scalar,
 *        positive scale, periodic angle, unit vector — plus the progress clamp
 *        they share. Every interpolator is exact at both endpoints and snaps to
 *        the nearer one outside [0,1], so a caller that has already clamped its
 *        progress pays nothing for the guard.
 */

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

/**
 * @brief Interpolation along the shorter arc of a periodic domain.
 * @param a Value at progress 0.
 * @param b Value at progress 1.
 * @param progress Fraction; outside [0,1] it snaps to the nearer endpoint.
 * @param period Domain period (e.g. TWO_PI_F for an angle).
 * @return The interpolated value wrapped into [0,period), or the linear
 *         interpolation when @p period is non-positive and no arc is defined.
 */
inline float shortest_periodic(float a, float b, float progress, float period) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
  if (period <= 0.0f)
    return linear(a, b, progress);
  const float half = period * 0.5f;
  float delta = std::fmod(b - a, period);
  if (delta < -half)
    delta += period;
  if (delta >= half)
    delta -= period;
  float value = std::fmod(a + delta * progress, period);
  if (value < 0.0f)
    value += period;
  return value >= period ? 0.0f : value;
}

/**
 * @brief Outcome of normalized_linear().
 * @tparam Size Component count of the interpolated vector.
 */
template <size_t Size> struct NormalizedLinearResult {
  std::array<float, Size> values; /**< Interpolated vector; unit-length when
                                       valid, un-normalized otherwise. */
  bool valid; /**< False when the interpolant passed through the origin. */
};

/**
 * @brief Linear interpolation of a vector, renormalized to unit length.
 * @tparam Size Component count.
 * @param a Vector at progress 0.
 * @param b Vector at progress 1.
 * @param progress Fraction; outside [0,1] it snaps to the nearer endpoint
 *        (returned as supplied, without renormalizing).
 * @param epsilon Norm below which the interpolant has no direction.
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
