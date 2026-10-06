/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "animation/animation.h"
#include "render/pullback/material.h"
#include "render/pullback/operators/model.h"
#include "render/pullback/stage.h"
#include "render/pullback/runtime_seeds.h"

/**
 * @file common.h
 * @brief Vocabulary the per-family operator headers share: instance-owned
 *        clocks and noise fields, and the shared noise-basis value table.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Animation::RandomWalk tuning for the spin-and-wander operators. */
inline constexpr Animation::RandomWalkOptions WALK_OPTIONS{};

/**
 * @brief Instance state of the spin-and-wander operators: an accumulated
 *        orientation plus the noise field that drives its random walk.
 */
struct SpatialWalkState {
  FastNoiseLite walk_noise;
  math::Vector position;
  math::Vector direction;
  math::Quaternion wander;
  float angular_velocity = 0.0f;
  float spin_phase = 0.0f;
  uint32_t walk_time = 0;
  int32_t noise_seed = 0;
};

inline void init_walk(SpatialWalkState &state, int32_t seed) {
  state.noise_seed = seed;
  init_effect_noise(state.walk_noise, seed);
  state.walk_noise.SetFrequency(WALK_OPTIONS.noise_scale);
  state.position = math::UP;
  state.direction = math::perpendicular_axis(state.position);
}

inline void advance_walk(SpatialWalkState &state, float wander,
                         float spin_rate) {
  ++state.walk_time;
  const Animation::RandomWalkDelta delta = Animation::step_random_walk<false>(
      state.position, state.direction, state.angular_velocity, state.walk_noise,
      WALK_OPTIONS, state.walk_time);
  state.wander =
      (math::scaled_rotation_delta(delta.rotation, wander) * state.wander)
          .normalized();
  state.spin_phase = fmodf(state.spin_phase + spin_rate, math::TWO_PI_F);
}

/** @brief Instance state of the noise-driven operators: the owned field plus
    the loop phase PhaseClockModel advances. */
struct NoisePhaseState {
  FastNoiseLite noise;
  float phase = 0.0f;
  int32_t noise_seed = 0;
};

/** @brief Value-state operator with one normalized loop clock. */
template <typename StateT> struct PhaseClockModel : ValueStateModel<StateT> {
  template <typename Params>
  static void advance(StateT &state, const Params &params) {
    if constexpr (requires { params.noise_time_rate; })
      state.phase = math::wrap_t(state.phase + params.noise_time_rate);
    else
      state.phase = math::wrap_t(state.phase + params.speed);
  }
};

__attribute__((always_inline)) inline bool
noise_plane_in_domain(const math::Complex &coords) {
  return fabsf(coords.re) <= 0x1p20f && fabsf(coords.im) <= 0x1p20f;
}

inline void init_noise_phase(NoisePhaseState &state, InstanceId id) {
  state.noise_seed = static_cast<int32_t>(id.stable_hash);
  init_effect_noise(state.noise, state.noise_seed);
}

/** @brief Noise-basis topology values, in math::NoiseBasis order. */
inline constexpr const char *NOISE_BASIS_IDS[] = {"simplex", "fbm3", "ridged3"};
static_assert(std::size(NOISE_BASIS_IDS) ==
              static_cast<size_t>(math::NoiseBasis::RIDGED3) + 1);

/**
 * @brief Bounds a noise-driven operator's basis enum8.
 * @details Per-pixel basis switches rely on this check and carry no guard.
 */
inline void check_noise_basis(uint8_t basis) {
  HS_CHECK(basis <= static_cast<uint8_t>(math::NoiseBasis::RIDGED3),
           "pullback operator: invalid noise basis");
}

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
