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
