/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Pattern parameter and state admission traps.

inline void case_pattern_sphere_cells_invalid_prepare() {
  Pattern::SphereCells::Params params;
  (void)Pattern::SphereCells::prepare(params);
  params.layout = 255;
  (void)Pattern::SphereCells::prepare(params);
}

inline void case_pattern_sine_invalid_axis() {
  Pattern::SineParams params;
  (void)Pattern::DirectionalSine::prepare(params, {});
  params.axis = math::Vector(0.0f, 0.0f, 0.0f);
  (void)Pattern::DirectionalSine::prepare(params, {});
}
inline void case_pattern_sine_invalid_state() {
  const Pattern::SineParams params;
  (void)Pattern::DirectionalSine::advance({}, params);
  (void)Pattern::DirectionalSine::advance(
      {std::numeric_limits<float>::infinity()}, params);
}
inline void case_pattern_interference_invalid_advance() {
  Pattern::InterferenceParams params;
  (void)Pattern::SphereInterference::advance({}, params);
  params.axis_b = math::Vector(0.0f, 0.0f, 0.0f);
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
  params.axis_b = math::Vector(0.0f, 0.0f, 0.0f);
  (void)Pattern::SphereSmoothGrid::advance({}, params);
}
inline void case_pattern_grid_invalid_prepare() {
  Pattern::GridParams params;
  (void)Pattern::SphereSmoothGrid::prepare(params, {});
  params.axis_b = math::Vector(0.0f, 0.0f, 0.0f);
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
  const FastNoiseLite generator;
  (void)Pattern::SphereNoise::prepare(params, {}, generator);
  params.scale = 0.0f;
  (void)Pattern::SphereNoise::prepare(params, {}, generator);
}
inline void case_pattern_radian_invalid_advance() {
  Pattern::PlaneInterference::Params params;
  (void)Pattern::PlaneInterference::advance({}, params);
  params.speed = std::numeric_limits<float>::infinity();
  (void)Pattern::PlaneInterference::advance({}, params);
}
inline void case_pattern_radian_invalid_prepare() {
  Pattern::PlaneInterference::Params params;
  (void)Pattern::PlaneInterference::prepare(params, {});
  (void)Pattern::PlaneInterference::prepare(params, {0.0f, 0.0f, 100.0f});
}
inline void case_pattern_bands_invalid_prepare() {
  Pattern::PlaneBands::Params params;
  (void)Pattern::PlaneBands::prepare(params, {});
  params.angle = std::numeric_limits<float>::infinity();
  (void)Pattern::PlaneBands::prepare(params, {});
}
inline void case_pattern_spiral_invalid_prepare() {
  Pattern::PlaneSpiral::Params params;
  (void)Pattern::PlaneSpiral::prepare(params, {});
  params.arms = -1.0f;
  (void)Pattern::PlaneSpiral::prepare(params, {});
}
inline void case_pattern_tessellation_invalid_prepare() {
  using Cells = Pattern::PlaneTessellation<Pattern::TessellationKind::SQUARE,
                                           Pattern::TessellationFill::CELLS>;
  Cells::Params params;
  (void)Cells::prepare(params, {});
  params.cell_radius = -1.0f;
  (void)Cells::prepare(params, {});
}
inline void case_pattern_plane_noise_invalid_advance() {
  Pattern::PlaneNoise::Params params;
  (void)Pattern::PlaneNoise::advance({}, params);
  (void)Pattern::PlaneNoise::advance({std::numeric_limits<float>::infinity()},
                                     params);
}
inline void case_pattern_plane_noise_invalid_prepare() {
  Pattern::PlaneNoise::Params params;
  const FastNoiseLite generator;
  (void)Pattern::PlaneNoise::prepare(params, {}, generator,
                                     math::NoiseBasis::SIMPLEX);
  (void)Pattern::PlaneNoise::prepare(params, {}, generator,
                                     static_cast<math::NoiseBasis>(3));
}

/** @brief Sphere direction admission rejects an infinite vector. */
inline void case_pattern_direction_infinite() {
  (void)Pattern::SphereDirection::from(math::Vector(opaque(INFINITY), 0, 0));
}

inline void case_pattern_direction_invalid() {
  (void)Pattern::SphereDirection::from(math::Vector(opaque(0.0f), 0.0f, 0.0f));
}

inline void case_pattern_cube_invalid_direction() {
  (void)Pattern::prepare_sphere_lut_sample(
      math::Vector(opaque(0.0f), 0.0f, 0.0f));
}
