/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file fields.h
 * @brief Typed field admission, source formulas and clocks. */

#include <limits>
#include <type_traits>
#include "core/pattern/sources.h"
#include "core/pattern/noise.h"
#include "tests/test_harness.h"

namespace hs_test::pattern_tests {

/** @brief Direction tags reject invalid lengths and preserve admitted scaling. */
inline void test_pattern_direction_admission() {
  using Pattern::SphereDirection;
  static_assert(!std::is_convertible_v<math::Vector, SphereDirection>);
  static_assert(!std::is_convertible_v<Pattern::Position3D, SphereDirection>);
  static_assert(
      !std::is_convertible_v<Pattern::PlanePosition, SphereDirection>);
  const float INF = std::numeric_limits<float>::infinity();
  const float NAN_VALUE = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_TRUE(!SphereDirection::valid({0.0f, 0.0f, 0.0f}));
  HS_EXPECT_TRUE(!SphereDirection::valid({0.5e-6f, 0.0f, 0.0f}));
  HS_EXPECT_TRUE(SphereDirection::valid({1.0e-6f, 0.0f, 0.0f}));
  HS_EXPECT_TRUE(!SphereDirection::valid({INF, 1.0f, 0.0f}));
  HS_EXPECT_TRUE(!SphereDirection::valid({NAN_VALUE, 1.0f, 0.0f}));
  HS_EXPECT_TRUE(!SphereDirection::valid({1.0e30f, 1.0e30f, 0.0f}));
  const auto A = SphereDirection::from({0.3f, -0.6f, 0.9f});
  const auto B = SphereDirection::from({3.0f, -6.0f, 9.0f});
  HS_EXPECT_NEAR(A.value().x, B.value().x, 1.0e-7f);
  HS_EXPECT_NEAR(A.value().y, B.value().y, 1.0e-7f);
  HS_EXPECT_NEAR(A.value().z, B.value().z, 1.0e-7f);
}

/** @brief Sphere sources use independent admitted axes and shared turn phases. */
inline void test_pattern_sphere_formulas() {
  using namespace Pattern;
  SineParams params;
  params.frequency = 0.25f;
  const auto Y = SphereDirection::from(math::Y_AXIS);
  const auto X = SphereDirection::from(math::X_AXIS);
  const auto PREPARED = DirectionalSine::prepare(params, {});
  HS_EXPECT_EQ(DirectionalSine::sample(Y, PREPARED), 1.0f);
  HS_EXPECT_EQ(DirectionalSine::sample(X, PREPARED), 0.0f);
  HS_EXPECT_EQ(DirectSphereField<DirectionalSine>::sample(Y, PREPARED), 1.0f);
  params.offset = 0.25f;
  HS_EXPECT_EQ(DirectionalSine::sample(X, DirectionalSine::prepare(params, {})),
               1.0f);
  InterferenceParams waves;
  waves.frequency = 0.25f;
  waves.frequency_b = 0.5f;
  waves.offset_b = 0.25f;
  const auto PAIR = SphereInterference::prepare(waves, {});
  HS_EXPECT_EQ(SphereInterference::sample(Y, PAIR), 1.0f);
  HS_EXPECT_NEAR(SphereInterference::sample(X, PAIR), -0.5f, 1.0e-6f);
  GridParams grid;
  grid.frequency = 0.25f;
  grid.frequency_b = 0.5f;
  grid.offset_b = 0.25f;
  const auto GRID = SphereSmoothGrid::prepare(grid, {});
  HS_EXPECT_EQ(SphereSmoothGrid::sample(Y, GRID), 1.0f);
  HS_EXPECT_EQ(SphereSmoothGrid::sample(X, GRID), 0.0f);
  for (int i = -64; i <= 64; ++i) {
    const auto D = SphereDirection::from({0.37f, i * 0.125f, -0.91f});
    const float EXPECTED_A = std::sin(math::TWO_PI_F * 0.25f * D.value().y);
    const float EXPECTED_B =
        std::sin(math::TWO_PI_F * (0.5f * D.value().x + 0.25f));
    HS_EXPECT_NEAR(SphereInterference::sample(D, PAIR),
                   0.5f * (EXPECTED_A + EXPECTED_B), 0.0018f);
    HS_EXPECT_NEAR(SphereSmoothGrid::sample(D, GRID), EXPECTED_A * EXPECTED_B,
                   0.0035f);
    HS_EXPECT_TRUE(std::abs(SphereSmoothGrid::sample(D, GRID)) <= 1.0f);
    HS_EXPECT_NEAR(SphereInterference::sample(D, PAIR),
                   SphereInterference::sample(
                       D, SphereInterference::prepare(waves, {1.0f})),
                   1.0e-6f);
  }
  params.frequency = MAX_AXIAL_FREQUENCY + 1.0f;
  HS_EXPECT_TRUE(!DirectionalSine::valid(params));
  params.frequency = -1.0f;
  HS_EXPECT_TRUE(!DirectionalSine::valid(params));
  params.frequency = 1.0f;
  params.axis = math::Vector(0.0f, 0.0f, 0.0f);
  HS_EXPECT_TRUE(!DirectionalSine::valid(params));
  grid.axis = math::Vector(0.0f, 2.0f, 0.0f);
  grid.offset = -0.17f;
  grid.offset_b = 0.23f;
  for (const math::Vector AXIS_B :
       {math::Vector(3.0f, 3.0f, 0.0f), math::Vector(0.0f, 4.0f, 0.0f),
        math::Vector(0.0f, -4.0f, 0.0f)}) {
    grid.axis_b = AXIS_B;
    HS_EXPECT_TRUE(SphereSmoothGrid::valid(grid));
    if (!SphereSmoothGrid::valid(grid))
      return;
    const auto SKEWED = SphereSmoothGrid::prepare(grid, {0.125f});
    const float LENGTH_B = std::sqrt(AXIS_B.x * AXIS_B.x + AXIS_B.y * AXIS_B.y);
    HS_EXPECT_NEAR(SKEWED.a.axis.y, 1.0f, 1.0e-7f);
    HS_EXPECT_NEAR(SKEWED.b.axis.x, AXIS_B.x / LENGTH_B, 1.0e-7f);
    HS_EXPECT_NEAR(SKEWED.b.axis.y, AXIS_B.y / LENGTH_B, 1.0e-7f);
    for (int i = -16; i <= 16; ++i) {
      const auto D = SphereDirection::from({0.37f, i * 0.125f, -0.91f});
      const float EXPECTED_A =
          std::sin(math::TWO_PI_F * (0.25f * D.value().y + 0.125f - 0.17f));
      const float EXPECTED_B = std::sin(
          math::TWO_PI_F *
          (0.5f * (AXIS_B.x * D.value().x + AXIS_B.y * D.value().y) / LENGTH_B +
           0.125f + 0.23f));
      const float ACTUAL = SphereSmoothGrid::sample(D, SKEWED);
      HS_EXPECT_NEAR(ACTUAL, EXPECTED_A * EXPECTED_B, 0.0035f);
      HS_EXPECT_TRUE(std::abs(ACTUAL) <= 1.0f);
    }
  }
  const GridParams VALID_GRID = grid;
  for (bool second : {false, true}) {
    for (const math::Vector INVALID_AXIS :
         {math::Vector(0.0f, 0.0f, 0.0f),
          math::Vector(std::numeric_limits<float>::infinity(), 1.0f, 0.0f),
          math::Vector(std::numeric_limits<float>::quiet_NaN(), 1.0f, 0.0f)}) {
      grid = VALID_GRID;
      (second ? grid.axis_b : grid.axis) = INVALID_AXIS;
      HS_EXPECT_TRUE(!SphereSmoothGrid::valid(grid));
    }
    for (float frequency : {-1.0f, MAX_AXIAL_FREQUENCY + 1.0f}) {
      grid = VALID_GRID;
      (second ? grid.frequency_b : grid.frequency) = frequency;
      HS_EXPECT_TRUE(!SphereSmoothGrid::valid(grid));
    }
  }
}

/** @brief Turn clocks and legacy radian clocks preserve their declared rates. */
inline void test_pattern_source_clocks() {
  using namespace Pattern;
  SineParams sine;
  sine.speed = -0.125f;
  auto state = DirectionalSine::init_state(sine);
  HS_EXPECT_EQ(DirectionalSine::advance(state, sine).phase_turns, 0.875f);
  for (int i = 0; i < 8; ++i)
    state = DirectionalSine::advance(state, sine);
  HS_EXPECT_EQ(state.phase_turns, 0.0f);
  for (float speed : {-1.0417e-6f, -1.0e-7f, -1.5e-4f, 1.0417e-6f}) {
    sine.speed = speed;
    NoiseParams noise;
    noise.speed = speed;
    InterferenceParams interference;
    interference.speed = speed;
    GridParams sphere_grid;
    sphere_grid.speed = speed;
    PhaseState clock{0.375f};
    float expected = clock.phase_turns;
    for (int frame = 0; frame < 128; ++frame) {
      expected = static_cast<float>(static_cast<double>(expected) + speed);
      HS_EXPECT_EQ(SphereNoise::advance(clock, noise).phase_turns, expected);
      HS_EXPECT_EQ(SphereInterference::advance(clock, interference).phase_turns,
                   expected);
      HS_EXPECT_EQ(SphereSmoothGrid::advance(clock, sphere_grid).phase_turns,
                   expected);
      clock = DirectionalSine::advance(clock, sine);
      HS_EXPECT_EQ(clock.phase_turns, expected);
    }
    sine.offset = speed;
    HS_EXPECT_EQ(DirectionalSine::prepare(sine, {0.375f}).phase_turns,
                 static_cast<float>(0.375 + speed));
  }
  GridSourceParams grid;
  grid.speed = -0.2f;
  grid.secondary_rate = 1.25f;
  grid.angle_rate = 0.03f;
  float primary = 0.1f, secondary = 0.2f, angle = 0.4f;
  advance_clocks(grid, primary, secondary, angle);
  HS_EXPECT_NEAR(primary, -0.1f, 1.0e-8f);
  HS_EXPECT_NEAR(secondary, -0.05f, 1.0e-8f);
  HS_EXPECT_NEAR(angle, 0.43f, 1.0e-8f);
}

/** @brief Noise variants preserve configured sampling and octave contrast. */
inline void test_pattern_noise_formulas() {
  using namespace Pattern;
  FastNoiseLite generator;
  generator.SetSeed(719);
  generator.SetFrequency(0.7f);
  NoiseParams params;
  params.scale = 3.5f;
  const PhaseState STATE{0.37f};
  const auto D = SphereDirection::from({0.4f, -0.3f, 0.7f});
  const auto SINGLE = SphereNoise::prepare(params, STATE, generator);
  const auto Q = params.scale * D.value() +
                 math::noise_sphere_loop_offset(STATE.phase_turns);
  const float EXPECTED = generator.GetNoiseSingle(Q.x, Q.y, Q.z);
  HS_EXPECT_EQ(SphereNoise::sample(D, SINGLE), EXPECTED);
  generator.SetSeed(-12);
  HS_EXPECT_EQ(SphereNoise::sample(D, SINGLE), EXPECTED);
  generator.SetSeed(719);
  params.variant = NoiseVariant::OCTAVE_CONTOUR;
  params.basis = NoiseBasis::FBM3;
  params.contrast = 2.0f;
  const float OCTAVES =
      (generator.GetNoiseSingle(Q.x, Q.y, Q.z) +
       0.5f * generator.GetNoiseSingle(2.0f * Q.x, 2.0f * Q.y, 2.0f * Q.z) +
       0.25f * generator.GetNoiseSingle(4.0f * Q.x, 4.0f * Q.y, 4.0f * Q.z)) /
      1.75f;
  HS_EXPECT_NEAR(
      SphereNoise::sample(D, SphereNoise::prepare(params, STATE, generator)),
      OCTAVES * 3.0f / (1.0f + 2.0f * std::abs(OCTAVES)), 1.0e-5f);
  params.variant = NoiseVariant::CONFIGURED_SINGLE;
  HS_EXPECT_TRUE(!SphereNoise::valid(params));
}

/** @brief Extracted sources expose pure clocks and immutable prepared controls. */
inline void test_pattern_extracted_source_protocols() {
  using namespace Pattern;
  PlaneInterference::Params waves;
  waves.speed = -0.2f;
  waves.secondary_rate = 1.5f;
  waves.angle_rate = 0.03f;
  const auto ZERO = PlaneInterference::init_state(waves);
  const auto NEXT = PlaneInterference::advance(ZERO, waves);
  HS_EXPECT_EQ(ZERO.primary, 0.0f);
  HS_EXPECT_EQ(NEXT.primary, -0.2f);
  HS_EXPECT_EQ(NEXT.secondary, -0.3f);
  HS_EXPECT_EQ(NEXT.angle, 0.03f);
  const PlanePosition POINT{{0.6f, 0.25f}};
  const auto PREPARED = PlaneInterference::prepare(waves, NEXT);
  const float ROTATED =
      POINT.value.re * std::cos(0.03f) + POINT.value.im * std::sin(0.03f);
  HS_EXPECT_NEAR(PlaneInterference::sample(POINT, PREPARED),
                 0.5f * (std::sin(0.4f) + std::sin(ROTATED - 0.3f)), 0.0018f);
  const float RETAINED = PlaneInterference::sample(POINT, PREPARED);
  waves.pattern_freq = 13.0f;
  HS_EXPECT_EQ(PlaneInterference::sample(POINT, PREPARED), RETAINED);
  PlaneRings::Params rings_params;
  HS_EXPECT_NEAR(
      PlaneRings::sample(POINT,
                         PlaneRings::prepare(rings_params, {0.2f, 0.0f, 0.0f})),
      std::sin(0.45f), 0.0018f);
  PlaneSpiral::Params spiral_params;
  HS_EXPECT_NEAR(PlaneSpiral::sample(
                     {{0.5f, 0.0f}},
                     PlaneSpiral::prepare(spiral_params, {0.2f, 0.0f, 0.1f})),
                 0.0f, 0.0018f);
  PlaneCoupledGrid::Params grid_params;
  grid_params.pattern_mix = 1.0f;
  HS_EXPECT_NEAR(PlaneCoupledGrid::sample(
                     POINT, PlaneCoupledGrid::prepare(grid_params, {})),
                 std::sin(0.6f) * std::cos(0.25f), 0.0035f);
  static_assert(std::is_empty_v<PlanePrimitiveLattice::State>);
  PlanePrimitiveLattice::Params lattice;
  const auto EMPTY = PlanePrimitiveLattice::advance(
      PlanePrimitiveLattice::init_state(lattice), lattice);
  HS_EXPECT_EQ(PlanePrimitiveLattice::sample(
                   {{}}, PlanePrimitiveLattice::prepare(lattice, EMPTY)),
               1.0f);
  HS_EXPECT_EQ(
      PlaneEscapeFractal::sample({{}}, PlaneEscapeFractal::prepare({}, {})),
      1.0f);
  using Square = PlaneTessellation<TessellationKind::SQUARE>;
  HS_EXPECT_EQ(Square::sample({{}}, Square::prepare({}, {})), -1.0f);
  const auto AXIS = SphereDirection::from(math::Y_AXIS);
  const auto DIRECTION = SphereDirection::from(math::X_AXIS);
  HS_EXPECT_EQ(
      SphereRings::sample(DIRECTION, SphereRings::prepare({}, {}, AXIS)), 1.0f);
  const FastNoiseLite GENERATOR;
  PlaneNoise::Params noise;
  noise.noise_time_rate = -0.01f;
  const auto CLOCK = PlaneNoise::advance(PlaneNoise::init_state(noise), noise);
  HS_EXPECT_EQ(CLOCK.phase_turns, 0.99f);
  const auto NOISE =
      PlaneNoise::prepare(noise, CLOCK, GENERATOR, math::NoiseBasis::FBM3);
  const auto Q = math::noise_projected_coordinate(
      POINT.value, noise.noise_scale, CLOCK.phase_turns);
  HS_EXPECT_EQ(PlaneNoise::sample(POINT, NOISE),
               noise_contour(GENERATOR, math::NoiseBasis::FBM3, Q, 0.0f));
}

} // namespace hs_test::pattern_tests
