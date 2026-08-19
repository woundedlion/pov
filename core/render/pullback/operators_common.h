/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

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

/** @brief Spatial frequency of the operator-owned orientation walk. */
inline constexpr float WALK_NOISE_SCALE = 0.02f;
/** @brief Arc length of one full-amplitude walk step, in radians. */
inline constexpr float WALK_STEP_RADIANS = 0.02f;

/** @brief Deterministic per-frame walk delta from an operator-owned noise
    field; sampled per axis so the pivot direction drifts smoothly. */
inline Quaternion walk_delta(const FastNoiseLite &noise, float t) {
  const Vector axis(noise.GetNoiseSingle(t, 0.0f, 0.0f),
                    noise.GetNoiseSingle(0.0f, t, 31.7f),
                    noise.GetNoiseSingle(-17.3f, t, 0.0f));
  const float length = axis.length();
  if (length < 1e-5f)
    return Quaternion();
  return make_rotation(axis * (1.0f / length), length * WALK_STEP_RADIANS);
}

/** @brief Configures an operator-owned noise field, seeded from the instance
    hash so a fresh init of the same pair reconstructs it deterministically. */
inline void init_instance_noise(FastNoiseLite &noise, InstanceId id) {
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(static_cast<int>(id.stable_hash));
  noise.SetFrequency(1.0f);
}

/**
 * @brief Instance state of the spin-and-wander operators: an accumulated
 *        orientation plus the noise field that drives its random walk.
 */
struct SpatialWalkState {
  FastNoiseLite walk_noise;
  Quaternion wander;
  float spin_phase = 0.0f;
  float walk_time = 0.0f;
};

inline void init_walk(SpatialWalkState &state, InstanceId id) {
  init_instance_noise(state.walk_noise, id);
}

inline void advance_walk(SpatialWalkState &state, float wander,
                         float spin_rate) {
  state.walk_time += 1.0f;
  const Quaternion delta =
      walk_delta(state.walk_noise, state.walk_time * WALK_NOISE_SCALE);
  state.wander =
      (slerp(Quaternion(), delta, wander) * state.wander).normalized();
  state.spin_phase = fmodf(state.spin_phase + spin_rate, TWO_PI_F);
}

/** @brief Instance state of the noise-driven operators: the owned field plus
    the loop phase its `speed` field advances. */
struct NoisePhaseState {
  FastNoiseLite noise;
  float phase = 0.0f;
};

inline void init_noise_phase(NoisePhaseState &state, InstanceId id) {
  init_instance_noise(state.noise, id);
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
