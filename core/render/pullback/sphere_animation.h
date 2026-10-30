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
