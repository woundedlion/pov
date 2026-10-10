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

/**
 * @brief Typed operator models; make_operator_descriptor() erases each into
 *        one OPERATOR_TABLE entry.
 * @details Every model declares this contract (optional members marked):
 * - `ID` — `const char *` stable versioned operator id, unique in the table.
 * - `NAME` — `const char *` display name, unique in the table.
 * - `Input`, `Output` — canonical carriers; the family rank may not decrease.
 * - `Params` — trivially copyable family block with default member values;
 *   its `FIELDS` (`Field<Params>` array) lists the float fields and the
 *   optional `TOPOLOGY` (`TopologyField<Params>` array) its structural enum8s.
 * - `State` — default-constructible instance state; usually inherited from
 *   ValueStateModel, PhaseClockModel or StatelessModel.
 * - `Prepared` — trivially destructible per-frame block built by `prepare`.
 * - `init(State &, InstanceId)` — seeds owned resources from the instance id.
 * - `migrate(State &dst, const State &src, InstanceId)` — clones into a fresh
 *   `dst`; returns Status.
 * - `advance(State &, const Params &)` — steps per-frame clocks.
 * - `prepare(const FrameContext &, const Params &, const State &)` — returns
 *   the frame's `Prepared`.
 * - `run(const Input &, const FrameContext &, const Params &,
 *   const Prepared &)` — returns the per-sample `Output`.
 * - `plane_bound(const Params &, float input_bound)` — required when `Output`
 *   is PlaneSample: output magnitude bound given the input's.
 * - `validate(const Params &)` — optional; admission warning or null. A model
 *   declaring it must SNAP every field or declare `ADMISSIBILITY_CONVEXITY`.
 * - `ADMISSIBILITY_CONVEXITY` — optional non-empty string stating the
 *   convexity contract of the admissible parameter set.
 * - `EDGE_DISTANCE_AVAILABLE` — optional, default false; a projection's plane
 *   carries edge distance downstream.
 * - `APPROXIMATE`, `ORACLE`, `METRICS`, `NON_FLOATING_FIELDS_EXACT` —
 *   optional approximation metadata; an approximate operator declares an
 *   oracle and metrics including a final framebuffer bound.
 * - `validate_frame`, `project` — projection-family hooks: frame-policy check
 *   and the per-family local-direction projection.
 */
namespace Op {

/** @brief Animation::RandomWalk tuning for the spin-and-wander operators. */
inline constexpr Animation::RandomWalkOptions WALK_OPTIONS{};

/**
 * @brief Instance state of the spin-and-wander operators: an accumulated
 *        orientation plus the noise field that drives its random walk.
 */
struct SpatialWalkState {
  FastNoiseLite walk_noise;      ///< Noise field steering the random walk.
  math::Vector position;         ///< Walker position on the unit sphere.
  math::Vector direction;        ///< Walker heading, tangent at `position`.
  math::Quaternion wander;       ///< Accumulated wander orientation.
  float angular_velocity = 0.0f; ///< Walker angular velocity.
  float spin_phase = 0.0f;       ///< Spin angle, radians, |value| < 2π.
  uint32_t walk_time = 0;        ///< Walk steps taken.
  int32_t noise_seed = 0;        ///< Seed of `walk_noise`.
};

/**
 * @brief Seeds the walk noise and places the walker at the pole.
 * @param state Walk state to initialise.
 * @param seed Noise seed.
 */
inline void init_walk(SpatialWalkState &state, int32_t seed) {
  state.noise_seed = seed;
  init_effect_noise(state.walk_noise, seed);
  state.walk_noise.SetFrequency(WALK_OPTIONS.noise_scale);
  state.position = math::UP;
  state.direction = math::perpendicular_axis(state.position);
}

/**
 * @brief Steps the random walk one frame and accumulates wander and spin.
 * @param state Walk state to advance.
 * @param wander Fraction of the walk rotation absorbed this frame, [0, 1].
 * @param spin_rate Spin angle added this frame, in radians.
 */
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
  FastNoiseLite noise;    ///< Instance-owned noise field.
  float phase = 0.0f;     ///< Normalized loop phase in [0, 1).
  int32_t noise_seed = 0; ///< Seed of `noise`.
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

/**
 * @brief Whether both plane coordinates lie within the noise input domain.
 * @param coords Plane coordinates.
 * @return True when |re| and |im| are at most 2^20.
 */
__attribute__((always_inline)) inline bool
noise_plane_in_domain(const math::Complex &coords) {
  return fabsf(coords.re) <= 0x1p20f && fabsf(coords.im) <= 0x1p20f;
}

/**
 * @brief Seeds the instance noise field from the instance hash.
 * @param state Noise state to initialise.
 * @param id Chain-entry identity; `stable_hash` becomes the seed.
 */
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
 * @param basis Raw `math::NoiseBasis` value.
 */
inline void check_noise_basis(uint8_t basis) {
  HS_CHECK(basis <= static_cast<uint8_t>(math::NoiseBasis::RIDGED3),
           "pullback operator: invalid noise basis");
}

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
