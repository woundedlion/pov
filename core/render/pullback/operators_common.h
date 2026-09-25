/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "animation/animation.h"
#include "render/pullback/material.h"
#include "render/pullback/operator_model.h"
#include "render/pullback/stage.h"

/**
 * @file operators_common.h
 * @brief Vocabulary the per-family operator headers share: instance-owned
 *        clocks and noise fields, and the topology enum8 value tables.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Walk tuning shared with Animation::RandomWalk, which owns the
    recurrence the spin-and-wander operators step. */
inline constexpr Animation::RandomWalkOptions WALK_OPTIONS{};

inline void init_effect_noise(FastNoiseLite &noise, int32_t seed) {
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(seed);
  noise.SetFrequency(1.0f);
}

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
};

inline void init_walk(SpatialWalkState &state, int32_t seed) {
  init_effect_noise(state.walk_noise, seed);
  state.walk_noise.SetFrequency(WALK_OPTIONS.noise_scale);
  state.position = math::UP;
  state.direction = math::perpendicular_axis(state.position);
}

// Chain walks accumulate once per frame using instance seeds; composed walks
// use eased Rotation samples and effect seeds, so nonzero wander is not equivalent.
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
    the loop phase its `speed` field advances. */
struct NoisePhaseState {
  FastNoiseLite noise;
  float phase = 0.0f;
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

inline void init_noise_phase(NoisePhaseState &state, InstanceId id) {
  init_effect_noise(state.noise, static_cast<int32_t>(id.stable_hash));
}

/** @brief Noise-basis topology values, in math::NoiseBasis order. */
inline constexpr const char *NOISE_BASIS_IDS[] = {"simplex", "fbm3", "ridged3"};
static_assert(std::size(NOISE_BASIS_IDS) ==
              static_cast<size_t>(math::NoiseBasis::RIDGED3) + 1);

/**
 * @brief Bounds a noise-driven operator's basis enum8.
 * @details Called from prepare(), once per frame, so the per-pixel basis
 * switches stay total and carry no guard.
 */
inline void check_noise_basis(uint8_t basis) {
  HS_CHECK(basis <= static_cast<uint8_t>(math::NoiseBasis::RIDGED3),
           "pullback operator: invalid noise basis");
}

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
