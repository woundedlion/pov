/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace FixedPipeline {

enum class InterpolationTrait : uint8_t {
  LINEAR,
  LOG_POSITIVE,
  SHORTEST_PERIODIC,
  NORMALIZED_LINEAR
};

enum class PathPolicy : uint8_t { PARALLEL, STAGGERED_ORDERED };

enum class Easing : uint8_t { LINEAR, EASE_IN_OUT_SIN };

enum class AbsentEdgeFallback : uint8_t { SNAP, REJECT };

struct EdgeProgress {
  float raw;
  float eased;
};

constexpr float clamp_progress(float value) {
  if (value <= 0.0f)
    return 0.0f;
  if (value >= 1.0f)
    return 1.0f;
  return value;
}

inline float ease(Easing easing, float progress) {
  const float t = clamp_progress(progress);
  if (t == 0.0f || t == 1.0f || easing == Easing::LINEAR)
    return t;
  constexpr float HALF_TURN_RADIANS = 3.14159265358979323846f;
  return (1.0f - std::cos(HALF_TURN_RADIANS * t)) * 0.5f;
}

inline EdgeProgress edge_progress(uint32_t evaluation, uint32_t duration,
                                  Easing easing) {
  const float raw = duration == 0U
                        ? 0.0f
                        : clamp_progress(static_cast<float>(evaluation) /
                                         static_cast<float>(duration));
  return {raw, ease(easing, raw)};
}

constexpr float staggered_group_progress(float progress, size_t group,
                                         size_t group_count) {
  if (group_count == 0U || group >= group_count)
    return 0.0f;
  return clamp_progress(progress * static_cast<float>(group_count) -
                        static_cast<float>(group));
}

constexpr float linear(float a, float b, float progress) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
  return a + (b - a) * progress;
}

inline float log_positive(float a, float b, float progress) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
  return std::exp(linear(std::log(a), std::log(b), progress));
}

inline float shortest_periodic(float a, float b, float progress, float period) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
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

template <size_t Size> struct NormalizedLinearResult {
  std::array<float, Size> values;
  bool valid;
};

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
