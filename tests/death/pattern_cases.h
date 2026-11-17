/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Death cases for pattern parameter and state guards. */

#include "tests/death/context.h"

namespace hs_test::death_tests {

/** @brief Death case: pattern sphere cells invalid prepare. */
inline void case_pattern_sphere_cells_invalid_prepare() {
  Pattern::SphereCells::Params params;
  (void)Pattern::SphereCells::prepare(params);
  death_admitted();
  params.layout = 255;
  (void)Pattern::SphereCells::prepare(params);
}

/** @brief Death case: pattern sine invalid axis. */
inline void case_pattern_sine_invalid_axis() {
  Pattern::SineParams params;
  (void)Pattern::DirectionalSine::prepare(params, {});
  death_admitted();
  params.axis = math::Vector(0.0f, 0.0f, 0.0f);
  (void)Pattern::DirectionalSine::prepare(params, {});
}
/** @brief Death case: pattern sine invalid state. */
inline void case_pattern_sine_invalid_state() {
  const Pattern::SineParams params;
  (void)Pattern::DirectionalSine::advance({}, params);
  death_admitted();
  (void)Pattern::DirectionalSine::advance(
      {std::numeric_limits<float>::infinity()}, params);
}
/** @brief Death case: pattern interference invalid advance. */
inline void case_pattern_interference_invalid_advance() {
  Pattern::InterferenceParams params;
  (void)Pattern::SphereInterference::advance({}, params);
  death_admitted();
  params.axis_b = math::Vector(0.0f, 0.0f, 0.0f);
  (void)Pattern::SphereInterference::advance({}, params);
}
/** @brief Death case: pattern interference invalid prepare. */
inline void case_pattern_interference_invalid_prepare() {
  Pattern::InterferenceParams params;
  (void)Pattern::SphereInterference::prepare(params, {});
  death_admitted();
  params.frequency_b = -1.0f;
  (void)Pattern::SphereInterference::prepare(params, {});
}
/** @brief Death case: pattern grid invalid advance. */
inline void case_pattern_grid_invalid_advance() {
  Pattern::GridParams params;
  (void)Pattern::SphereSmoothGrid::advance({}, params);
  death_admitted();
  params.axis_b = math::Vector(0.0f, 0.0f, 0.0f);
  (void)Pattern::SphereSmoothGrid::advance({}, params);
}
/** @brief Death case: pattern grid invalid prepare. */
inline void case_pattern_grid_invalid_prepare() {
  Pattern::GridParams params;
  (void)Pattern::SphereSmoothGrid::prepare(params, {});
  death_admitted();
  params.axis_b = math::Vector(0.0f, 0.0f, 0.0f);
  (void)Pattern::SphereSmoothGrid::prepare(params, {});
}
/** @brief Death case: pattern noise invalid advance. */
inline void case_pattern_noise_invalid_advance() {
  Pattern::NoiseParams params;
  (void)Pattern::SphereNoise::advance({}, params);
  death_admitted();
  params.speed = std::numeric_limits<float>::quiet_NaN();
  (void)Pattern::SphereNoise::advance({}, params);
}
/** @brief Death case: pattern noise invalid prepare. */
inline void case_pattern_noise_invalid_prepare() {
  Pattern::NoiseParams params;
  const FastNoiseLite generator;
  (void)Pattern::SphereNoise::prepare(params, {}, generator);
  death_admitted();
  params.scale = 0.0f;
  (void)Pattern::SphereNoise::prepare(params, {}, generator);
}
/** @brief Death case: pattern radian invalid advance. */
inline void case_pattern_radian_invalid_advance() {
  Pattern::PlaneInterference::Params params;
  (void)Pattern::PlaneInterference::advance({}, params);
  death_admitted();
  params.speed = std::numeric_limits<float>::infinity();
  (void)Pattern::PlaneInterference::advance({}, params);
}
/** @brief Death case: pattern radian invalid prepare. */
inline void case_pattern_radian_invalid_prepare() {
  Pattern::PlaneInterference::Params params;
  (void)Pattern::PlaneInterference::prepare(params, {});
  death_admitted();
  (void)Pattern::PlaneInterference::prepare(params, {0.0f, 0.0f, 100.0f});
}
/** @brief Death case: pattern bands invalid prepare. */
inline void case_pattern_bands_invalid_prepare() {
  Pattern::PlaneBands::Params params;
  (void)Pattern::PlaneBands::prepare(params, {});
  death_admitted();
  params.angle = std::numeric_limits<float>::infinity();
  (void)Pattern::PlaneBands::prepare(params, {});
}
/** @brief Death case: pattern spiral invalid prepare. */
inline void case_pattern_spiral_invalid_prepare() {
  Pattern::PlaneSpiral::Params params;
  (void)Pattern::PlaneSpiral::prepare(params, {});
  death_admitted();
  params.arms = -1.0f;
  (void)Pattern::PlaneSpiral::prepare(params, {});
}
/** @brief Death case: pattern tessellation invalid prepare. */
inline void case_pattern_tessellation_invalid_prepare() {
  using Cells = Pattern::PlaneTessellation<Pattern::TessellationKind::SQUARE,
                                           Pattern::TessellationFill::CELLS>;
  Cells::Params params;
  (void)Cells::prepare(params, {});
  death_admitted();
  params.cell_radius = -1.0f;
  (void)Cells::prepare(params, {});
}
/** @brief Death case: pattern plane noise invalid advance. */
inline void case_pattern_plane_noise_invalid_advance() {
  Pattern::PlaneNoise::Params params;
  (void)Pattern::PlaneNoise::advance({}, params);
  death_admitted();
  (void)Pattern::PlaneNoise::advance({std::numeric_limits<float>::infinity()},
                                     params);
}
/** @brief Death case: pattern plane noise invalid prepare. */
inline void case_pattern_plane_noise_invalid_prepare() {
  Pattern::PlaneNoise::Params params;
  const FastNoiseLite generator;
  (void)Pattern::PlaneNoise::prepare(params, {}, generator,
                                     math::NoiseBasis::SIMPLEX);
  death_admitted();
  (void)Pattern::PlaneNoise::prepare(params, {}, generator,
                                     static_cast<math::NoiseBasis>(3));
}

/** @brief Sphere direction admission rejects an infinite vector. */
inline void case_pattern_direction_infinite() {
  (void)Pattern::SphereDirection::from(math::Vector(opaque(INFINITY), 0, 0));
}

/** @brief Death case: pattern direction invalid. */
inline void case_pattern_direction_invalid() {
  (void)Pattern::SphereDirection::from(math::Vector(opaque(0.0f), 0.0f, 0.0f));
}

/** @brief Death case: pattern cube invalid direction. */
inline void case_pattern_cube_invalid_direction() {
  (void)Pattern::prepare_sphere_lut_sample(
      math::Vector(opaque(0.0f), 0.0f, 0.0f));
}

} // namespace hs_test::death_tests
