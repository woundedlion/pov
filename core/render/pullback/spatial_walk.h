/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spatial_walk.h
 * @brief Shared sphere spin and random-walk state. */

#include "core/animation/animation.h"
#include "core/render/pullback/runtime_seeds.h"
#include "core/render/pullback/fields.h"
#include "core/render/spacetime.h"

namespace Pullback::SphereAnimation {

/** @brief Controls of one keyed sphere rotation. */
struct RotationParams {
  float wander = 0.0f;     ///< Fraction of each walk delta absorbed.
  float spin_speed = 0.0f; ///< Spin about Y, radians per frame.
  /// Numeric controls and interpolation curves.
  static constexpr auto FIELDS = std::array{
      Field<RotationParams>{"wander", &RotationParams::wander, "Wander", 0.0f,
                            1.0f, FieldCurve::LERP},
      Field<RotationParams>{"spin-speed", &RotationParams::spin_speed,
                            "Spin Speed", 0.0f, 0.05f, FieldCurve::LERP}};
};
static_assert(field_ids_unique<RotationParams>());
static_assert(field_defaults_in_range<RotationParams>());

/** @brief Animation::RandomWalk tuning for the spin-and-wander operators. */
inline constexpr Animation::RandomWalkOptions WALK_OPTIONS{};

/**
 * @brief Instance state of the spin-and-wander operators: an accumulated
 *        orientation plus the noise field that drives its random walk.
 * @details Spin values are equivalent modulo 2pi; column-time samples may
 * lie outside the canonical range.
 */
struct SpatialWalkState {
  FastNoiseLite walk_noise;      ///< Noise field steering the random walk.
  math::Vector position;         ///< Walker position on the unit sphere.
  math::Vector direction;        ///< Walker heading, tangent at `position`.
  math::Quaternion wander;       ///< Accumulated wander orientation.
  float angular_velocity = 0.0f; ///< Walker angular velocity.
  float spin_phase = 0.0f;       ///< Spin angle, radians modulo 2pi.
  uint32_t walk_time = 0;        ///< Walk steps taken.
  int32_t noise_seed = 0;        ///< Seed of `walk_noise`.
};

/**
 * @brief Seeds the walk noise and places the walker at the pole.
 * @param state Walk state to initialise.
 * @param seed Noise seed.
 */
HS_FLASH_MEMBER inline void init_walk(SpatialWalkState &state, int32_t seed) {
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
 * @param speed Walk distance per frame, in radians.
 */
HS_FLASH_MEMBER inline void advance_walk(SpatialWalkState &state, float wander,
                                         float spin_rate,
                                         float speed = WALK_OPTIONS.speed) {
  ++state.walk_time;
  auto options = WALK_OPTIONS;
  options.speed = speed;
  const Animation::RandomWalkDelta delta = Animation::step_random_walk(
      state.position, state.direction, state.angular_velocity, state.walk_noise,
      options, state.walk_time);
  state.wander =
      math::absorb_rotation_delta(delta.rotation, wander, state.wander);
  state.spin_phase = fmodf(state.spin_phase + spin_rate, math::TWO_PI_F);
}

/** @brief Prepared frame clocks of SpatialWalkState. */
struct SpatialWalkInterval {
  Spacetime::TemporalRotation wander; ///< Prepared wander interval.
  Spacetime::ClockInterval spin;      ///< Prepared spin interval.
};

/** @brief Prepares the clocks between consecutive endpoints.
 * @param from Previous state.
 * @param to Current state.
 * @return Prepared interval. */
HS_FLASH_INLINE inline SpatialWalkInterval
temporal_interval(const SpatialWalkState &from, const SpatialWalkState &to) {
  return {{from.wander, to.wander},
          {from.spin_phase, to.spin_phase, math::TWO_PI_F}};
}

/** @brief Applies the prepared clocks to a column state.
 * @tparam Params Operator parameter block.
 * @param out Copy of the current state.
 * @param interval Prepared frame interval.
 * @param fraction Frame fraction in [0, 1]. */
template <typename Params>
HS_FLASH_INLINE void temporal_sample(SpatialWalkState &out,
                                     const SpatialWalkInterval &interval,
                                     const Params &, float fraction) {
  out.wander = interval.wander.sample(fraction);
  out.spin_phase = interval.spin.sample(fraction);
}

} // namespace Pullback::SphereAnimation
