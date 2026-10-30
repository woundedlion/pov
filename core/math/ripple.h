/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file ripple.h
 * @brief Spherical Ricker-wavelet geometry and ordered wave composition. */

#include "core/math/3dmath.h"
#include <array>

namespace RippleKernel {

/** @brief Instantaneous wave parameters and cached cosine support. */
struct Wave {
  math::Vector center;             ///< Unit source direction.
  float amplitude = 0.0f;          ///< Signed angular amplitude.
  float phase = 0.0f;              ///< Wavefront distance in radians.
  float decay{5.0};                ///< Spatial decay rate.
  float thickness{1.0f};           ///< Wavelet thickness in radians.
  float cos_threshold_min = 1.0f;  ///< Upper cosine support bound.
  float cos_threshold_max = -1.0f; ///< Lower cosine support bound.

  /** @brief Returns the wavelet half-width, floored at 0.001 radians. */
  float half_width() const {
    float hw = thickness * 0.5f;
    return hw < 0.001f ? 0.001f : hw;
  }

  /** @brief Refreshes the cosine support bounds from phase and thickness. */
  void sync() {
    float hw = half_width();
    float d_min = hs::clamp(phase - hw * 2.0f, 0.0f, math::PI_F);
    float d_max = hs::clamp(phase + hw * 2.0f, 0.0f, math::PI_F);
    cos_threshold_min = cosf(d_min);
    cos_threshold_max = cosf(d_max);
  }
};

/** @brief Largest angle evaluated with the truncated quaternion series. */
inline constexpr float SMALL_ANGLE_MAX = 0.15f;

/** @brief Tests exact center identity. */
HS_HOT_INLINE inline bool same_center(const math::Vector &a,
                                      const math::Vector &b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

/** @brief Includes a wave's support in a cosine-band union. */
HS_HOT_INLINE inline void include_support(float &near, float &far,
                                          const Wave &wave) {
  near = fmaxf(near, wave.cos_threshold_min);
  far = fminf(far, wave.cos_threshold_max);
}

/** @brief Tests the closed cosine support band. */
HS_HOT_INLINE inline bool contains(const Wave &wave, float cosine) {
  return !(cosine > wave.cos_threshold_min || cosine < wave.cos_threshold_max);
}

/** @brief Evaluates the wavelet at an angular distance and phase. */
HS_HOT_INLINE inline float angle_at_distance(const Wave &wave, float distance,
                                             float phase) {
  float dist_from_peak = distance - phase;
  float half_width = wave.half_width();
  float t = (dist_from_peak / half_width) * 2.0f;
  return wave.amplitude * (1.0f - t * t) *
         math::fast_expf(-0.5f * t * t - wave.decay * distance);
}

/** @brief Rotates around a normalized axis with the ripple angle policy. */
HS_HOT_INLINE inline math::Vector
rotate(const math::Vector &v, const math::Vector &axis, float theta) {
  if (fabsf(theta) <= SMALL_ANGLE_MAX) {
    float h = 0.5f * theta;
    float h2 = h * h;
    float s = h * (1.0f - h2 * (1.0f / 6.0f));
    float c = 1.0f - h2 * (0.5f - h2 * (1.0f / 24.0f));
    return math::rotate(v, math::Quaternion(c, s * axis));
  }
  math::Quaternion q = math::make_rotation(axis, theta);
  return math::rotate(v, q);
}

/** @brief Applies one wave at a known distance without testing its support. */
HS_HOT_INLINE inline math::Vector transform_at_distance(const math::Vector &v,
                                                        const Wave &wave,
                                                        float distance,
                                                        float phase) {
  const float THETA = angle_at_distance(wave, distance, phase);
  math::Vector axis = math::cross(wave.center, v);
  float len_sq = math::dot(axis, axis);
  if (len_sq > 1e-6f) {
    axis = axis * (1.0f / sqrtf(len_sq));
    return rotate(v, axis, THETA);
  }
  return v;
}

/** @brief Applies one wave after amplitude and cosine support rejection. */
HS_HOT_INLINE inline math::Vector transform(const math::Vector &v,
                                            const Wave &wave) {
  if (wave.amplitude <= 0.001f)
    return v;
  float cos_d = math::dot(v, wave.center);
  if (!contains(wave, cos_d))
    return v;
  float d = math::fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
  return transform_at_distance(v, wave, d, wave.phase);
}

/** @brief Net radial rotation and travel; negative travel marks a pole crossing. */
struct RadialTransfer {
  float angle = 0.0f; ///< Signed net rotation in radians.
  float path = 0.0f;  ///< Absolute angular travel.
};

/** @brief Advances one angularly supported wave, rejecting pole crossings. */
HS_HOT_INLINE inline bool radial_step(const Wave &wave, float &distance,
                                      float &angle, float &path,
                                      bool track_path) {
  const float HALF_WIDTH = wave.half_width();
  if (wave.amplitude <= 0.001f || distance < wave.phase - 2.0f * HALF_WIDTH ||
      distance > wave.phase + 2.0f * HALF_WIDTH)
    return true;
  const float THETA = angle_at_distance(wave, distance, wave.phase);
  if (THETA == 0.0f)
    return true;
  distance += THETA;
  if (!(distance > 0.01f && distance < math::PI_F - 0.01f))
    return false;
  angle += THETA;
  if (track_path)
    path += fabsf(THETA);
  return true;
}

/** @brief Resolves a same-center wave sequence in radial coordinates.
 * @tparam WaveAt Inline accessor returning a wave reference. */
template <typename WaveAt>
HS_HOT_INLINE inline RadialTransfer
radial_transfer(float distance, const WaveAt &wave_at, uint8_t count) {
  RadialTransfer result;
  for (uint8_t i = 0; i < count; ++i)
    if (!radial_step(wave_at(i), distance, result.angle, result.path, true))
      return {0.0f, -1.0f};
  return result;
}

/** @brief Applies a grouped rotation when its axis and unit-vector bounds hold. */
HS_HOT_INLINE inline bool rotate_group(const math::Vector &start,
                                       const math::Vector &center, float angle,
                                       math::Vector &output) {
  const math::Vector CROSS = math::cross(center, start);
  const float CROSS_SQ = math::dot(CROSS, CROSS);
  if (!(CROSS_SQ >= 1e-4f) || !(fabsf(math::dot(start, start) - 1.0f) <= 1e-4f))
    return false;
  const math::Vector AXIS = CROSS * (1.0f / sqrtf(CROSS_SQ));
  output = rotate(start, AXIS, angle);
  return true;
}

/** @brief Applies waves sequentially and measures their actual angular travel.
 * @tparam Result Carrier with sphere and path_length members.
 * @tparam WaveAt Inline accessor returning a wave reference. */
template <typename Result, typename WaveAt>
HS_HOT_INLINE inline Result sequential(const math::Vector &input,
                                       const WaveAt &wave_at, uint8_t count,
                                       bool track_path) {
  Result result{input, 0.0f};
  for (uint8_t index = 0; index < count; ++index) {
    const math::Vector START = result.sphere;
    result.sphere = transform(START, wave_at(index));
    if (track_path) {
      const math::Vector CROSS = math::cross(START, result.sphere);
      result.path_length += atan2f(sqrtf(math::dot(CROSS, CROSS)),
                                   math::dot(START, result.sphere));
    }
  }
  return result;
}

/** @brief Composes consecutive same-center waves in one radial basis.
 * @tparam Result Carrier with sphere and path_length members.
 * @tparam WaveAt Inline accessor returning a wave reference.
 * @tparam Fallback Inline callable resolving the original input sequentially. */
template <typename Result, typename WaveAt, typename Fallback>
HS_HOT_INLINE inline Result grouped(const math::Vector &input,
                                    const WaveAt &wave_at, uint8_t count,
                                    bool track_path, const Fallback &fallback) {
  Result result{input, 0.0f};
  uint8_t index = 0;
  while (index < count) {
    const math::Vector START = result.sphere;
    const math::Vector CENTER = wave_at(index).center;
    const float DISTANCE_COS = math::dot(CENTER, START);
    float distance = math::fast_acos(hs::clamp(DISTANCE_COS, -1.0f, 1.0f));
    float total_angle = 0.0f;
    do {
      if (!radial_step(wave_at(index++), distance, total_angle,
                       result.path_length, track_path))
        return fallback();
    } while (index < count && same_center(wave_at(index).center, CENTER));
    if (total_angle != 0.0f &&
        !rotate_group(START, CENTER, total_angle, result.sphere))
      return fallback();
  }
  return result;
}

/** @brief Samples neighboring nodes from two radial-transfer curves. */
template <size_t Samples>
HS_HOT_INLINE inline bool
sample_transfer(const std::array<RadialTransfer, Samples> &previous,
                const std::array<RadialTransfer, Samples> &current,
                float distance, float fraction, bool track_path,
                RadialTransfer &output) {
  const float POSITION = distance * ((Samples - 1) / math::PI_F);
  if (!(POSITION >= 0.0f && POSITION <= Samples - 1))
    return false;
  const size_t LOW_INDEX = std::min(static_cast<size_t>(POSITION), Samples - 2);
  const float WEIGHT = POSITION - LOW_INDEX;
  const auto &a = previous[LOW_INDEX];
  const auto &b = previous[LOW_INDEX + 1];
  const auto &c = current[LOW_INDEX];
  const auto &d = current[LOW_INDEX + 1];
  if (a.path < 0.0f || b.path < 0.0f || c.path < 0.0f || d.path < 0.0f)
    return false;
  output.angle = hs::lerp(hs::lerp(a.angle, b.angle, WEIGHT),
                          hs::lerp(c.angle, d.angle, WEIGHT), fraction);
  if (track_path)
    output.path = hs::lerp(hs::lerp(a.path, b.path, WEIGHT),
                           hs::lerp(c.path, d.path, WEIGHT), fraction);
  return true;
}

/** @brief A wave pointer and bounds of its consecutive same-center group. */
template <typename Params> struct EndpointGroup {
  const Params *params = nullptr; ///< Wave in spawn order.
  float near = 0.0f;              ///< Upper cosine bound of the group.
  float far = 0.0f;               ///< Lower cosine bound of the group.
  int end = 0;                    ///< Exclusive end of the group.
};

/** @brief Appends one active wave to caller-owned endpoint group storage. */
template <typename Params>
HS_HOT_INLINE inline void append_endpoint(EndpointGroup<Params> *waves,
                                          int &count, int &first,
                                          const Params &params) {
  auto &wave = waves[count];
  wave.params = &params;
  if (count == 0 || !same_center(params.center, waves[first].params->center)) {
    first = count;
    wave.near = params.cos_threshold_min;
    wave.far = params.cos_threshold_max;
  } else
    include_support(waves[first].near, waves[first].far, params);
  waves[first].end = ++count;
}

/** @brief Rejects endpoint groups, then applies surviving waves sequentially.
 * @tparam Params Wave parameter type.
 * @tparam OnGroup Inline group-visit observer.
 * @tparam OnWave Inline individual-band observer. */
template <typename Params, typename OnGroup, typename OnWave>
HS_HOT_INLINE inline math::Vector
sequential_groups(math::Vector v, const EndpointGroup<Params> *waves, int count,
                  const OnGroup &on_group, const OnWave &on_wave) {
  for (int first = 0; first < count;) {
    const auto &group = waves[first];
    const float cosine = math::dot(v, group.params->center);
    on_group();
    if (cosine <= group.near && cosine >= group.far) {
      for (int k = first; k < group.end; ++k) {
        const auto &p = *waves[k].params;
        const float c = k == first ? cosine : math::dot(v, p.center);
        on_wave();
        if (!contains(p, c))
          continue;
        v = transform_at_distance(
            v, p, math::fast_acos(hs::clamp(c, -1.0f, 1.0f)), p.phase);
      }
    }
    first = group.end;
  }
  return v;
}

} // namespace RippleKernel
