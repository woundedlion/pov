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

/** @brief Walk tuning shared with Animation::RandomWalk, which owns the
    recurrence the spin-and-wander operators step. */
inline constexpr Animation::RandomWalkOptions WALK_OPTIONS{};

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
  state.walk_noise.SetFrequency(WALK_OPTIONS.noise_scale);
  state.position = UP;
  state.direction =
      cross(state.position, least_parallel_axis(state.position)).normalized();
}

inline void advance_walk(SpatialWalkState &state, float wander,
                         float spin_rate) {
  ++state.walk_time;
  const Animation::RandomWalkDelta delta = Animation::step_random_walk<false>(
      state.position, state.direction, state.angular_velocity, state.walk_noise,
      WALK_OPTIONS, state.walk_time);
  state.wander = (scaled_rotation_delta(delta.rotation, wander) * state.wander)
                     .normalized();
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

/**
 * @brief Bounds a noise-driven operator's basis enum8.
 * @details Called from prepare(), once per frame, so the per-pixel basis
 * switches stay total and carry no guard.
 */
inline void check_noise_basis(uint8_t basis) {
  HS_CHECK(basis <= static_cast<uint8_t>(::NoiseBasis::RIDGED3),
           "pullback operator: invalid noise basis");
}

enum class WeightMode : uint8_t { NONE = 0, PROJECTION = 1 };
using ProjectionCoverageMode = Pullback::ProjectionCoverageMode;

inline constexpr const char *WEIGHT_MODE_IDS[] = {"none", "projection"};
inline constexpr const char *COVERAGE_MODE_IDS[] = {
    "none", "weight", "weight-squared", "edge-fade"};

/**
 * @brief Bounds the Sample crossing's topology enum8s.
 * @details Called from prepare(), once per frame, so the per-pixel switches
 * below stay total and carry no guard.
 */
inline void check_sample_topology(uint8_t weight_mode, uint8_t coverage_mode) {
  HS_CHECK(weight_mode <= static_cast<uint8_t>(WeightMode::PROJECTION),
           "sample operator: invalid weight mode");
  HS_CHECK(coverage_mode <=
               static_cast<uint8_t>(ProjectionCoverageMode::EDGE_FADE),
           "sample operator: invalid projection coverage mode");
}

/** @brief The Sample crossing's weight switch over the shared policies. */
inline float weighted_field(uint8_t weight_mode, float raw,
                            const ProjectionProvenance &provenance,
                            const FrameContext &ctx) {
  switch (static_cast<WeightMode>(weight_mode)) {
  case WeightMode::NONE:
    return Weight::None::apply(raw, provenance, ctx);
  case WeightMode::PROJECTION:
    break;
  }
  return Weight::Projection::apply(raw, provenance, ctx);
}

/** @brief The Sample crossing's coverage switch over the shared policies. */
inline float projection_coverage(uint8_t coverage_mode,
                                 const ProjectionProvenance &provenance,
                                 float edge_width, const FrameContext &ctx) {
  switch (static_cast<ProjectionCoverageMode>(coverage_mode)) {
  case ProjectionCoverageMode::NONE:
    return ProjectionCoverage::None::apply(provenance, ctx);
  case ProjectionCoverageMode::WEIGHT_SQUARED:
    return ProjectionCoverage::WeightSquared::apply(provenance, ctx);
  case ProjectionCoverageMode::EDGE_FADE:
    return ProjectionCoverage::edge_fade(provenance, edge_width);
  case ProjectionCoverageMode::WEIGHT:
    break;
  }
  return ProjectionCoverage::Weight::apply(provenance, ctx);
}

/** @brief The parameters every Sample crossing family carries: the edge-fade
    width and the weight and coverage enum8s. */
struct SampleCrossingParams {
  /** Edge-fade band width; read only under edge-fade coverage. */
  float edge_width = 0.1f;
  uint8_t weight_mode = static_cast<uint8_t>(WeightMode::PROJECTION);
  uint8_t coverage_mode = static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT);
};

/** @brief The tabled field of SampleCrossingParams, retyped to the family. */
template <typename Params>
constexpr std::array<Field<Params>, 1> sample_crossing_fields() {
  return {edge_width_field<Params>(&Params::edge_width)};
}

/**
 * @brief The Sample crossing's shared weight and coverage enum8s, followed by
 *        any family-specific topology fields in @p extra.
 */
template <typename Params, typename... Extra>
constexpr std::array<TopologyField<Params>, 2 + sizeof...(Extra)>
sample_crossing_topology(const Extra &...extra) {
  return {TopologyField<Params>{"weight-mode", &Params::weight_mode,
                                WEIGHT_MODE_IDS, 2,
                                static_cast<uint8_t>(WeightMode::PROJECTION)},
          TopologyField<Params>{
              "coverage-mode", &Params::coverage_mode, COVERAGE_MODE_IDS, 4,
              static_cast<uint8_t>(ProjectionCoverageMode::WEIGHT)},
          extra...};
}

/** @brief Builds the field carrier from a raw source value under the family's
    weight and coverage enum8s. */
template <typename Params>
__attribute__((always_inline)) inline FieldSample
finish_sample(const PlaneSample &input, float raw, const Params &params,
              const FrameContext &ctx) {
  return Kernel::sample(
      input, weighted_field(params.weight_mode, raw, input.provenance, ctx),
      projection_coverage(params.coverage_mode, input.provenance,
                          params.edge_width, ctx));
}

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
