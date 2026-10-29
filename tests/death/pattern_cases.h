/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file pattern_cases.h
 * @brief Pattern parameter and state admission traps. */

inline void case_pattern_sine_invalid_axis() {
  Pattern::SineParams params;
  (void)Pattern::DirectionalSine::prepare(params, {});
  params.axis = {};
  (void)Pattern::DirectionalSine::prepare(params, {});
}
inline void case_pattern_sine_invalid_state() {
  const Pattern::SineParams PARAMS;
  (void)Pattern::DirectionalSine::advance({}, PARAMS);
  (void)Pattern::DirectionalSine::advance(
      {std::numeric_limits<float>::infinity()}, PARAMS);
}
inline void case_pattern_interference_invalid_advance() {
  Pattern::InterferenceParams params;
  (void)Pattern::SphereInterference::advance({}, params);
  params.axis_b = {};
  (void)Pattern::SphereInterference::advance({}, params);
}
inline void case_pattern_interference_invalid_prepare() {
  Pattern::InterferenceParams params;
  (void)Pattern::SphereInterference::prepare(params, {});
  params.frequency_b = -1.0f;
  (void)Pattern::SphereInterference::prepare(params, {});
}
inline void case_pattern_grid_invalid_advance() {
  Pattern::GridParams params;
  (void)Pattern::SphereSmoothGrid::advance({}, params);
  params.axis_b = params.axis;
  (void)Pattern::SphereSmoothGrid::advance({}, params);
}
inline void case_pattern_grid_invalid_prepare() {
  Pattern::GridParams params;
  (void)Pattern::SphereSmoothGrid::prepare(params, {});
  params.axis_b = params.axis;
  (void)Pattern::SphereSmoothGrid::prepare(params, {});
}
inline void case_pattern_noise_invalid_advance() {
  Pattern::NoiseParams params;
  (void)Pattern::SphereNoise::advance({}, params);
  params.speed = std::numeric_limits<float>::quiet_NaN();
  (void)Pattern::SphereNoise::advance({}, params);
}
inline void case_pattern_noise_invalid_prepare() {
  Pattern::NoiseParams params;
  const FastNoiseLite GENERATOR;
  (void)Pattern::SphereNoise::prepare(params, {}, GENERATOR);
  params.scale = 0.0f;
  (void)Pattern::SphereNoise::prepare(params, {}, GENERATOR);
}
inline void case_pattern_radian_invalid_advance() {
  Pattern::PlaneRings::Params params;
  (void)Pattern::PlaneRings::advance({}, params);
  params.speed = std::numeric_limits<float>::infinity();
  (void)Pattern::PlaneRings::advance({}, params);
}
inline void case_pattern_radian_invalid_prepare() {
  Pattern::PlaneRings::Params params;
  (void)Pattern::PlaneRings::prepare(params, {});
  (void)Pattern::PlaneRings::prepare(params, {0.0f, 0.0f, 100.0f});
}
inline void case_pattern_lattice_invalid_prepare() {
  Pattern::PlanePrimitiveLattice::Params params;
  (void)Pattern::PlanePrimitiveLattice::prepare(params, {});
  params.lattice_radius = -1.0f;
  (void)Pattern::PlanePrimitiveLattice::prepare(params, {});
}
inline void case_pattern_plane_noise_invalid_advance() {
  Pattern::PlaneNoise::Params params;
  (void)Pattern::PlaneNoise::advance({}, params);
  (void)Pattern::PlaneNoise::advance({std::numeric_limits<float>::infinity()},
                                     params);
}
inline void case_pattern_plane_noise_invalid_prepare() {
  Pattern::PlaneNoise::Params params;
  const FastNoiseLite GENERATOR;
  (void)Pattern::PlaneNoise::prepare(params, {}, GENERATOR,
                                     math::NoiseBasis::SIMPLEX);
  (void)Pattern::PlaneNoise::prepare(params, {}, GENERATOR,
                                     static_cast<math::NoiseBasis>(3));
}
