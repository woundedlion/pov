/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/material.h"
#include "render/pullback/operator_model.h"
#include "render/pullback/runtime_seeds.h"
#include "render/pullback/stage.h"

/**
 * @file operators_common.h
 * @brief Vocabulary the per-family operator headers share: instance-owned
 *        clocks and noise fields, and the topology enum8 value tables.
 */

namespace Pullback {

namespace Interp {

namespace Op {

inline constexpr float WALK_SPEED = 0.02f;
inline constexpr float WALK_PIVOT_STRENGTH = 0.1f;
inline constexpr float WALK_NOISE_SCALE = 0.02f;
inline constexpr float WALK_SMOOTHING = 0.85f;
inline constexpr float WALK_DRIFT = 0.5f;

inline void init_effect_noise(FastNoiseLite &noise,
                              int32_t seed = EFFECT_NOISE_SEED) {
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
  Vector position;
  Vector direction;
  Quaternion wander;
  float angular_velocity = 0.0f;
  float spin_phase = 0.0f;
  uint32_t walk_time = 0;
};

inline void init_walk(SpatialWalkState &state, int32_t seed) {
  init_effect_noise(state.walk_noise, seed);
  state.walk_noise.SetFrequency(WALK_NOISE_SCALE);
  state.position = UP;
  state.direction =
      cross(state.position, least_parallel_axis(state.position)).normalized();
}

inline void advance_walk(SpatialWalkState &state, float wander,
                         float spin_rate) {
  ++state.walk_time;
  const float target_pivot =
      state.walk_noise.GetNoise(
          state.position.x * 100.0f, state.position.y * 100.0f,
          state.position.z * 100.0f +
              static_cast<float>(state.walk_time) * WALK_DRIFT) *
      WALK_PIVOT_STRENGTH;
  state.angular_velocity = state.angular_velocity * WALK_SMOOTHING +
                           target_pivot * (1.0f - WALK_SMOOTHING);
  state.direction =
      rotate(state.direction,
             make_rotation(state.position, state.angular_velocity))
          .normalized();
  const Vector axis_seed = least_parallel_axis(state.position);
  const Vector walk_axis =
      normalized_or(cross(state.position, state.direction),
                    cross(state.position, axis_seed).normalized());
  const Quaternion delta = make_rotation(walk_axis, WALK_SPEED);
  state.position = rotate(state.position, delta).normalized();
  state.direction = rotate(state.direction, delta).normalized();
  state.wander =
      (scaled_rotation_delta(delta, wander) * state.wander).normalized();
  state.spin_phase = fmodf(state.spin_phase + spin_rate, TWO_PI_F);
}

/** @brief Instance state of the noise-driven operators: the owned field plus
    the loop phase its `speed` field advances. */
struct NoisePhaseState {
  FastNoiseLite noise;
  float phase = 0.0f;
};

inline void init_noise_phase(NoisePhaseState &state, InstanceId id) {
  init_effect_noise(state.noise, static_cast<int32_t>(id.stable_hash));
}

/** @brief Per-frame source phase clocks. */
struct SourceClockState {
  float primary = 0.0f;
  float secondary = 0.0f;
  float angle = 0.0f;
};

/** @brief Noise-basis topology values, in ::NoiseBasis order. */
inline constexpr const char *NOISE_BASIS_IDS[] = {"simplex", "fbm3", "ridged3"};

enum class WeightMode : uint8_t { NONE = 0, PROJECTION = 1 };
enum class CoverageMode : uint8_t {
  NONE = 0,
  WEIGHT = 1,
  WEIGHT_SQUARED = 2,
  EDGE_FADE = 3
};

inline constexpr const char *WEIGHT_MODE_IDS[] = {"none", "projection"};
inline constexpr const char *COVERAGE_MODE_IDS[] = {
    "none", "weight", "weight-squared", "edge-fade"};

/** @brief The Sample crossing's weight switch over the shared policies. */
inline float weighted_field(uint8_t weight_mode, float raw,
                            const ProjectionProvenance &provenance,
                            const FrameContext &ctx) {
  switch (static_cast<WeightMode>(weight_mode)) {
  case WeightMode::NONE:
    return Weight::None::apply(raw, provenance, ctx);
  case WeightMode::PROJECTION:
  default:
    return Weight::Projection::apply(raw, provenance, ctx);
  }
}

/** @brief The Sample crossing's coverage switch over the shared policies. */
inline float provenance_coverage(uint8_t coverage_mode,
                                 const ProjectionProvenance &provenance,
                                 float edge_width, const FrameContext &ctx) {
  switch (static_cast<CoverageMode>(coverage_mode)) {
  case CoverageMode::NONE:
    return ProjectedCoverage::None::apply(provenance, ctx);
  case CoverageMode::WEIGHT_SQUARED:
    return ProjectedCoverage::WeightSquared::apply(provenance, ctx);
  case CoverageMode::EDGE_FADE:
    return ProjectedCoverage::edge_fade(provenance, edge_width);
  case CoverageMode::WEIGHT:
  default:
    return ProjectedCoverage::Weight::apply(provenance, ctx);
  }
}

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
