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
  const Pattern::Position3D POSITION{{3.0f, 4.0f, 0.0f}};
  HS_EXPECT_EQ(math::dot(POSITION.value, POSITION.value), 25.0f);
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
  params.axis = {};
  HS_EXPECT_TRUE(!DirectionalSine::valid(params));
  grid.axis_b = grid.axis;
  HS_EXPECT_TRUE(!SphereSmoothGrid::valid(grid));
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
  GridSourceParams grid;
  grid.speed = -0.2f;
  grid.secondary_rate = 1.25f;
  grid.angle_rate = 0.03f;
  float primary = 0.1f, secondary = 0.2f, angle = 0.4f;
  advance_clocks(grid, primary, secondary, angle);
  HS_EXPECT_EQ(primary, std::fmod(0.1f - 0.2f, math::TWO_PI_F));
  HS_EXPECT_EQ(secondary, std::fmod(0.2f - 0.2f * 1.25f, math::TWO_PI_F));
  HS_EXPECT_EQ(angle, std::fmod(0.4f + 0.03f, math::TWO_PI_F));
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

} // namespace hs_test::pattern_tests
