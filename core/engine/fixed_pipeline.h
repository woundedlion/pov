/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file fixed_pipeline.h
 * @brief Interpolation primitives for a fixed-look preset-to-preset edge:
 *        progress shaping plus one interpolator per parameter domain.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "math/geometry.h"

namespace FixedPipeline {

/** @brief Which interpolator a parameter's domain calls for. */
enum class InterpolationTrait : uint8_t {
  LINEAR,            /**< Unconstrained scalar; linear(). */
  LOG_POSITIVE,      /**< Strictly positive scale; log_positive(). */
  SHORTEST_PERIODIC, /**< Wrapping angle/phase; shortest_periodic(). */
  NORMALIZED_LINEAR  /**< Unit-length vector; normalized_linear(). */
};

/** @brief Whether the groups on an edge advance together or in sequence. */
enum class PathPolicy : uint8_t {
  PARALLEL,          /**< Every group sees the edge progress unchanged. */
  STAGGERED_ORDERED, /**< Groups run back to back; staggered_group_progress(). */
};

/** @brief Progress shaping applied on top of the raw evaluation fraction. */
enum class Easing : uint8_t { LINEAR, EASE_IN_OUT_SIN };

/** @brief Disposition of a requested edge the look table does not define. */
enum class AbsentEdgeFallback : uint8_t {
  SNAP,  /**< Adopt the destination with no interpolation. */
  REJECT /**< Refuse the edge and hold the source. */
};

/** @brief One edge's progress in both its raw and its eased form. */
struct EdgeProgress {
  float raw;   /**< Evaluation fraction, clamped to [0,1]. */
  float eased; /**< raw after the edge's Easing. */
};

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
 * @brief Shapes a progress fraction by an easing curve.
 * @param easing Curve to apply.
 * @param progress Fraction, clamped to [0,1] before shaping.
 * @return The shaped fraction; the endpoints 0 and 1 are exact under every
 *         curve.
 */
inline float ease(Easing easing, float progress) {
  const float t = clamp_progress(progress);
  if (t == 0.0f || t == 1.0f || easing == Easing::LINEAR)
    return t;
  constexpr float HALF_TURN_RADIANS = 3.14159265358979323846f;
  return (1.0f - std::cos(HALF_TURN_RADIANS * t)) * 0.5f;
}

/**
 * @brief Converts an edge's evaluation counter into raw and eased progress.
 * @param evaluation Elapsed evaluations on the edge.
 * @param duration Total evaluations the edge spans.
 * @param easing Curve applied to the raw fraction.
 * @return Both forms of the progress; a zero @p duration yields 0 progress.
 */
inline EdgeProgress edge_progress(uint32_t evaluation, uint32_t duration,
                                  Easing easing) {
  const float raw = duration == 0U
                        ? 0.0f
                        : clamp_progress(static_cast<float>(evaluation) /
                                         static_cast<float>(duration));
  return {raw, ease(easing, raw)};
}

/**
 * @brief Scales a rotation delta down to a fraction of its arc.
 * @param delta Full rotation delta; must be unit-length for the slerp path.
 * @param amount Fraction of the arc to keep.
 * @return @p delta at @p amount of its rotation; identity at 0, @p delta at 1.
 */
HS_HOT_FLASH_MEMBER inline Quaternion
scaled_rotation_delta(const Quaternion &delta, float amount) {
  if (amount == 1.0f)
    return delta;
  if (amount == 0.0f)
    return Quaternion();
  return slerp(Quaternion(), delta, amount);
}

/**
 * @brief Progress of one group under PathPolicy::STAGGERED_ORDERED.
 * @param progress Whole-edge progress.
 * @param group Zero-based group index.
 * @param group_count Number of groups sharing the edge.
 * @return The group's own [0,1] progress: it runs over the 1/group_count slice
 *         of the edge starting at @p group, and reads 0 for an out-of-range
 *         @p group or an empty @p group_count.
 */
constexpr float staggered_group_progress(float progress, size_t group,
                                         size_t group_count) {
  if (group_count == 0U || group >= group_count)
    return 0.0f;
  return clamp_progress(progress * static_cast<float>(group_count) -
                        static_cast<float>(group));
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

} // namespace FixedPipeline
