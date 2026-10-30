/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_animation.h
 * @brief Ripple burst scheduling. */

#include "core/render/pullback/surface.h"
#include "core/render/pullback/runtime_seeds.h"

namespace Pullback::SphereAnimation {

/** @brief Four-wave ripple burst controls. */
struct BurstRippleParams : Surface::PeriodicRippleParams {
  float burst_rate = 1.0f / 128.0f; ///< Bursts per frame.
  uint8_t burst_size = 3;           ///< Ripples per burst minus one.
  uint8_t random_center = 1; ///< One random center per burst when enabled.
  static constexpr auto FIELDS = concat_fields<BurstRippleParams>(
      Surface::PeriodicRippleParams::FIELDS,
      std::array{Field<BurstRippleParams>{
          "burst-rate", &BurstRippleParams::burst_rate, "Burst Rate", 0.0f,
          1.0f / 64.0f, FieldCurve::LERP}});
};

/** @brief One scheduled wave in a ripple burst. */
struct RippleWave {
  float age = 1000.0f; ///< Frames since launch; negative before launch.
  uint32_t burst = 0;  ///< Burst identity shared by its waves.
};
/** @brief Ripple burst scheduler state. */
struct RippleState {
  float phase = 0;          ///< Fraction of the interval until the next burst.
  uint32_t seed = 0;        ///< Instance seed for burst origins.
  uint32_t burst_count = 0; ///< Bursts launched.
  uint32_t frames = 0;      ///< Frames advanced.
  std::array<RippleWave, 16> waves{}; ///< Scheduled ripple waves.
};

/** @brief Advances the scheduler by one frame.
 * @tparam Params Ripple controls.
 * @param state Scheduler state.
 * @param params Current controls. */
template <typename Params>
inline void advance_ripple(RippleState &state, const Params &params) {
  ++state.frames;
  for (auto &wave : state.waves)
    wave.age = fminf(wave.age + 1.0f, 1000.0f);
  state.phase += params.burst_rate;
  const bool FIRST = state.burst_count == 0;
  if (!FIRST && state.phase < 1.0f)
    return;
  if (!FIRST)
    state.phase -= 1.0f;
  const float AGE = FIRST ? 1.0f : state.phase / params.burst_rate;
  const size_t GROUP = (state.burst_count % 4) * 4;
  for (uint8_t i = 0; i <= params.burst_size; ++i)
    state.waves[GROUP + i] = {AGE - 16.0f * i, state.burst_count};
  ++state.burst_count;
}

/** @brief Deterministic uniform center for a burst.
 * @param seed Instance seed.
 * @param burst Burst identity.
 * @return Unit center direction. */
inline math::Vector ripple_center(uint32_t seed, uint32_t burst) {
  uint32_t hash = seed ^ (burst * 0x9e3779b9U);
  hash = (hash ^ (hash >> 16)) * 0x7feb352dU;
  hash = (hash ^ (hash >> 15)) * 0x846ca68bU;
  hash ^= hash >> 16;
  const float Y = 2.0f * static_cast<float>(hash & 0xffffU) / 65535.0f - 1.0f;
  const float AZIMUTH =
      static_cast<float>(hash >> 16) * (math::TWO_PI_F / 65536.0f);
  const float RADIUS = sqrtf(fmaxf(0.0f, 1.0f - Y * Y));
  return {RADIUS * cosf(AZIMUTH), Y, RADIUS * sinf(AZIMUTH)};
}

/** @brief Active ripple kernels in launch order. */
struct PreparedRipples {
  std::array<Surface::PreparedRipple, 16> waves; ///< Active kernels.
  uint8_t count = 0;                             ///< Active wave count.
};

/** @brief Applies prepared waves in launch order.
 * @param input Unit direction.
 * @param prepared Prepared wave kernels.
 * @param track_path Whether to accumulate angular travel.
 * @return Displaced direction and travel. */
HS_HOT_ITCM_MEMBER inline SurfaceResult
sequential_ripples(const math::Vector &input, const PreparedRipples &prepared,
                   bool track_path) {
  SurfaceResult result{input, 0.0f};
  for (uint8_t index = 0; index < prepared.count; ++index) {
    const math::Vector START = result.sphere;
    result.sphere = ripple_transform(START, prepared.waves[index].ripple);
    if (track_path) {
      const math::Vector CROSS = math::cross(START, result.sphere);
      result.path_length += atan2f(sqrtf(math::dot(CROSS, CROSS)),
                                   math::dot(START, result.sphere));
    }
  }
  return result;
}

/** @brief Groups consecutive waves sharing a center into one radial basis.
 * @param input Unit direction.
 * @param prepared Prepared wave kernels in launch order.
 * @param track_path Whether to accumulate angular travel.
 * @return Displaced direction and travel. */
HS_HOT_ITCM_MEMBER inline SurfaceResult
grouped_ripples(const math::Vector &input, const PreparedRipples &prepared,
                bool track_path) {
  SurfaceResult result{input, 0.0f};
  uint8_t index = 0;
  while (index < prepared.count) {
    const math::Vector START = result.sphere;
    const math::Vector CENTER = prepared.waves[index].ripple.center;
    const float DISTANCE_COS = math::dot(CENTER, START);
    float distance = math::fast_acos(hs::clamp(DISTANCE_COS, -1.0f, 1.0f));
    float total_angle = 0.0f;
    do {
      const auto &wave = prepared.waves[index++].ripple;
      const float HALF_WIDTH = wave.half_width();
      if (wave.amplitude <= 0.001f ||
          distance < wave.phase - 2.0f * HALF_WIDTH ||
          distance > wave.phase + 2.0f * HALF_WIDTH)
        continue;
      const float T = (distance - wave.phase) / HALF_WIDTH * 2.0f;
      const float THETA =
          wave.amplitude * (1.0f - T * T) *
          math::fast_expf(-0.5f * T * T - wave.decay * distance);
      if (THETA == 0.0f)
        continue;
      distance += THETA;
      if (!(distance > 0.01f && distance < math::PI_F - 0.01f))
        return sequential_ripples(input, prepared, track_path);
      total_angle += THETA;
      if (track_path)
        result.path_length += fabsf(THETA);
    } while (index < prepared.count &&
             prepared.waves[index].ripple.center.x == CENTER.x &&
             prepared.waves[index].ripple.center.y == CENTER.y &&
             prepared.waves[index].ripple.center.z == CENTER.z);
    if (total_angle != 0.0f) {
      const math::Vector CROSS = math::cross(CENTER, START);
      const float CROSS_SQ = math::dot(CROSS, CROSS);
      if (!(CROSS_SQ >= 1e-4f) ||
          !(fabsf(math::dot(START, START) - 1.0f) <= 1e-4f))
        return sequential_ripples(input, prepared, track_path);
      const math::Vector AXIS = CROSS * (1.0f / sqrtf(CROSS_SQ));
      if (fabsf(total_angle) <= RIPPLE_SMALL_ANGLE_MAX) {
        const float H = total_angle * 0.5f, H2 = H * H;
        const float S = H * (1.0f - H2 * (1.0f / 6.0f));
        const float C = 1.0f - H2 * (0.5f - H2 * (1.0f / 24.0f));
        result.sphere = math::rotate(START, math::Quaternion(C, S * AXIS));
      } else
        result.sphere =
            math::rotate(START, math::make_rotation(AXIS, total_angle));
    }
  }
  return result;
}

/** @brief Prepares the active wave kernels.
 * @tparam Params Ripple controls.
 * @param params Current controls.
 * @param state Scheduler state.
 * @return Active kernels in launch order. */
template <typename Params>
inline PreparedRipples prepare_ripples(const Params &params,
                                       const RippleState &state) {
  PreparedRipples result;
  for (uint32_t lag = 4; lag > 0; --lag) {
    const uint32_t BURST = state.burst_count - lag;
    for (size_t i = 0; i < 4; ++i) {
      const auto &wave = state.waves[(BURST % 4) * 4 + i];
      if (wave.age < 0.0f || wave.age >= params.period || wave.burst != BURST)
        continue;
      auto prepared = Surface::prepare_ripple(params, wave.age / params.period);
      if (params.random_center)
        prepared.ripple.center = ripple_center(state.seed, BURST);
      result.waves[result.count++] = prepared;
    }
  }
  return result;
}

} // namespace Pullback::SphereAnimation
