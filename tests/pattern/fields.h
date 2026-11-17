/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file fields.h
 * @brief Typed field admission, source formulas and clocks. */

#include "core/math/scalar_noise.h"
#include <bit>
#include <limits>
#include <numbers>
#include <type_traits>
#include "core/pattern/sources.h"
#include "core/pattern/noise.h"
#include "core/pattern/sphere_cells.h"
#include "core/pattern/spider_web.h"
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
  const auto a = SphereDirection::from({0.3f, -0.6f, 0.9f});
  const auto b = SphereDirection::from({3.0f, -6.0f, 9.0f});
  HS_EXPECT_NEAR(a.value().x, b.value().x, 1.0e-7f);
  HS_EXPECT_NEAR(a.value().y, b.value().y, 1.0e-7f);
  HS_EXPECT_NEAR(a.value().z, b.value().z, 1.0e-7f);
  for (const math::Vector &input :
       {math::Vector(0.3f, -0.6f, 0.9f), math::Vector(3.0f, -6.0f, 9.0f),
        math::Vector(1.0e-6f, 0.0f, 0.0f)}) {
    const auto expected = SphereDirection::from(input).value();
    const auto actual = SphereDirection::admit(input).value();
    HS_EXPECT_NEAR(actual.x, expected.x, 1.0e-7f);
    HS_EXPECT_NEAR(actual.y, expected.y, 1.0e-7f);
    HS_EXPECT_NEAR(actual.z, expected.z, 1.0e-7f);
  }
  const math::Vector UNIT(0.6f, 0.8f, 0.0f);
  const auto preserved = SphereDirection::from_unit(UNIT).value();
  HS_EXPECT_EQ(preserved.x, UNIT.x);
  HS_EXPECT_EQ(preserved.y, UNIT.y);
  HS_EXPECT_EQ(preserved.z, UNIT.z);
}

/** @brief Sphere sources match their independent wave formulas. */
inline void test_pattern_sphere_formulas() {
  using namespace Pattern;
  SineParams params;
  params.frequency = 0.25f;
  const auto y = SphereDirection::from(math::Y_AXIS);
  const auto x = SphereDirection::from(math::X_AXIS);
  const auto prepared = DirectionalSine::prepare(params, {});
  HS_EXPECT_EQ(DirectionalSine::sample(y, prepared), 1.0f);
  HS_EXPECT_EQ(DirectionalSine::sample(x, prepared), 0.0f);
  HS_EXPECT_EQ(DirectSphereField<DirectionalSine>::sample(y, prepared), 1.0f);
  params.offset = 0.25f;
  HS_EXPECT_EQ(DirectionalSine::sample(x, DirectionalSine::prepare(params, {})),
               1.0f);
  InterferenceParams waves;
  waves.frequency = 0.25f;
  waves.frequency_b = 0.5f;
  waves.offset_b = 0.25f;
  const auto pair = SphereInterference::prepare(waves, {});
  HS_EXPECT_EQ(SphereInterference::sample(y, pair), 1.0f);
  HS_EXPECT_NEAR(SphereInterference::sample(x, pair), -0.5f, 1.0e-6f);
  GridParams grid;
  grid.frequency = 0.25f;
  grid.frequency_b = 0.5f;
  grid.offset_b = 0.25f;
  const auto grid_value = SphereSmoothGrid::prepare(grid, {});
  HS_EXPECT_EQ(SphereSmoothGrid::sample(y, grid_value), 1.0f);
  HS_EXPECT_EQ(SphereSmoothGrid::sample(x, grid_value), 0.0f);
  for (int i = -64; i <= 64; ++i) {
    const auto d = SphereDirection::from({0.37f, i * 0.125f, -0.91f});
    const float expected_a = std::sin(math::TWO_PI_F * 0.25f * d.value().y);
    const float expected_b =
        std::sin(math::TWO_PI_F * (0.5f * d.value().x + 0.25f));
    HS_EXPECT_NEAR(SphereInterference::sample(d, pair),
                   0.5f * (expected_a + expected_b), 0.0018f);
    HS_EXPECT_NEAR(SphereSmoothGrid::sample(d, grid_value),
                   expected_a * expected_b, 0.0035f);
    HS_EXPECT_TRUE(std::abs(SphereSmoothGrid::sample(d, grid_value)) <= 1.0f);
  }
}

/** @brief Sphere interference shares phase across a complete turn. */
inline void test_pattern_sphere_shared_phase() {
  using namespace Pattern;
  InterferenceParams waves;
  waves.frequency = 0.25f;
  waves.frequency_b = 0.5f;
  waves.offset_b = 0.25f;
  const auto pair = SphereInterference::prepare(waves, {});
  for (int i = -64; i <= 64; ++i) {
    const auto d = SphereDirection::from({0.37f, i * 0.125f, -0.91f});
    HS_EXPECT_NEAR(SphereInterference::sample(d, pair),
                   SphereInterference::sample(
                       d, SphereInterference::prepare(waves, {1.0f})),
                   1.0e-6f);
  }
}

/** @brief Sphere grid axes are normalized independently, including parallel axes. */
inline void test_pattern_sphere_skewed_axes() {
  using namespace Pattern;
  GridParams grid;
  grid.frequency = 0.25f;
  grid.frequency_b = 0.5f;
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
    const auto skewed = SphereSmoothGrid::prepare(grid, {0.125f});
    const float length_b = std::sqrt(AXIS_B.x * AXIS_B.x + AXIS_B.y * AXIS_B.y);
    HS_EXPECT_NEAR(skewed.a.axis.y, 1.0f, 1.0e-7f);
    HS_EXPECT_NEAR(skewed.b.axis.x, AXIS_B.x / length_b, 1.0e-7f);
    HS_EXPECT_NEAR(skewed.b.axis.y, AXIS_B.y / length_b, 1.0e-7f);
    for (int i = -16; i <= 16; ++i) {
      const auto d = SphereDirection::from({0.37f, i * 0.125f, -0.91f});
      const float expected_a =
          std::sin(math::TWO_PI_F * (0.25f * d.value().y + 0.125f - 0.17f));
      const float expected_b = std::sin(
          math::TWO_PI_F *
          (0.5f * (AXIS_B.x * d.value().x + AXIS_B.y * d.value().y) / length_b +
           0.125f + 0.23f));
      const float actual = SphereSmoothGrid::sample(d, skewed);
      HS_EXPECT_NEAR(actual, expected_a * expected_b, 0.0035f);
      HS_EXPECT_TRUE(std::abs(actual) <= 1.0f);
    }
  }
}

/** @brief Sphere source admission rejects invalid axes and frequencies. */
inline void test_pattern_sphere_admission() {
  using namespace Pattern;
  SineParams params;
  params.frequency = MAX_AXIAL_FREQUENCY + 1.0f;
  HS_EXPECT_TRUE(!DirectionalSine::valid(params));
  params.frequency = -1.0f;
  HS_EXPECT_TRUE(!DirectionalSine::valid(params));
  params.frequency = 1.0f;
  params.axis = math::Vector(0.0f, 0.0f, 0.0f);
  HS_EXPECT_TRUE(!DirectionalSine::valid(params));
  GridParams grid;
  grid.frequency = 0.25f;
  grid.frequency_b = 0.5f;
  grid.axis = math::Vector(0.0f, 2.0f, 0.0f);
  grid.axis_b = math::Vector(0.0f, -4.0f, 0.0f);
  grid.offset = -0.17f;
  grid.offset_b = 0.23f;
  const GridParams valid_grid = grid;
  for (bool second : {false, true}) {
    for (const math::Vector INVALID_AXIS :
         {math::Vector(0.0f, 0.0f, 0.0f),
          math::Vector(std::numeric_limits<float>::infinity(), 1.0f, 0.0f),
          math::Vector(std::numeric_limits<float>::quiet_NaN(), 1.0f, 0.0f)}) {
      grid = valid_grid;
      (second ? grid.axis_b : grid.axis) = INVALID_AXIS;
      HS_EXPECT_TRUE(!SphereSmoothGrid::valid(grid));
    }
    for (float frequency : {-1.0f, MAX_AXIAL_FREQUENCY + 1.0f}) {
      grid = valid_grid;
      (second ? grid.frequency_b : grid.frequency) = frequency;
      HS_EXPECT_TRUE(!SphereSmoothGrid::valid(grid));
    }
  }
}

/** @brief Turn clocks and radian clocks preserve their declared rates. */
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
  float primary = 0.1f, secondary = 0.2f, angle = 0.4f;
  advance_clocks(grid, primary, secondary, angle);
  HS_EXPECT_NEAR(primary, -0.1f, 1.0e-8f);
  HS_EXPECT_EQ(secondary, 0.2f);
  HS_EXPECT_EQ(angle, 0.4f);
  TwinWaveSourceParams waves;
  waves.speed = -0.2f;
  waves.secondary_rate = 1.25f;
  waves.angle_rate = 0.03f;
  advance_clocks(waves, primary, secondary, angle);
  HS_EXPECT_NEAR(primary, -0.3f, 1.0e-7f);
  HS_EXPECT_NEAR(secondary, -0.05f, 1.0e-8f);
  HS_EXPECT_NEAR(angle, 0.43f, 1.0e-8f);
}

/** @brief Turn clocks reduce large phases and increments independently. */
inline void test_pattern_turn_clock_reduction() {
  const auto fraction = [](double value) {
    double reduced = std::fmod(value, 1.0);
    return reduced < 0.0 ? reduced + 1.0 : reduced;
  };
  const float LIMIT = std::numeric_limits<float>::max();
  for (const float phase : {-LIMIT, -16777216.0f, -2.75f, -0.125f, 0.0f, 0.125f,
                            0.75f, 1.0f, 2.75f, 16777216.0f, LIMIT}) {
    for (const float rate : {-LIMIT, -16777216.0f, -1.25f, -1.0f, -0.25f, 0.0f,
                             0.25f, 1.0f, 1.25f, 16777216.0f, LIMIT}) {
      const float sum =
          static_cast<float>(fraction(phase)) +
          static_cast<float>(std::fmod(static_cast<double>(rate), 1.0));
      const float expected = static_cast<float>(fraction(sum));
      HS_EXPECT_EQ(Pattern::advance_turn_phase(phase, rate), expected);
      Pattern::NoiseSourceParams params;
      params.noise_time_rate = rate;
      HS_EXPECT_EQ(Pattern::advance_noise_phase(phase, params), expected);
    }
  }
  for (uint32_t index = 1; index <= 4096; ++index) {
    const float rate = std::bit_cast<float>(index * 2654435761u);
    if (!std::isfinite(rate))
      continue;
    const float fractional_rate =
        static_cast<float>(std::fmod(static_cast<double>(rate), 1.0));
    for (const float phase : {0.0f, 0.375f, 0.99999994f}) {
      const float sum = phase + fractional_rate;
      const float reduced = static_cast<float>(fraction(sum));
      const float expected = reduced >= 1.0f ? 0.0f : reduced;
      HS_EXPECT_EQ(Pattern::advance_turn_phase(phase, rate), expected);
      Pattern::NoiseSourceParams params;
      params.noise_time_rate = rate;
      HS_EXPECT_EQ(Pattern::advance_noise_phase(phase, params), expected);
    }
  }
  for (int step = 0; step < 4096; ++step) {
    const float phase = static_cast<float>(step) / 4096.0f;
    for (const float rate :
         {-0.999f, -0.25f, -1.0e-7f, 0.0f, 1.0e-7f, 0.25f, 0.999f}) {
      const float single_wrap = math::wrap_t(phase + rate);
      HS_EXPECT_EQ(Pattern::advance_turn_phase(phase, rate), single_wrap);
      Pattern::NoiseSourceParams params;
      params.noise_time_rate = rate;
      HS_EXPECT_EQ(Pattern::advance_noise_phase(phase, params), single_wrap);
    }
  }
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
  const auto d = SphereDirection::from({0.4f, -0.3f, 0.7f});
  const auto single = SphereNoise::prepare(params, STATE, generator);
  const auto q = params.scale * d.value() +
                 math::noise_sphere_loop_offset(STATE.phase_turns);
  const float expected = generator.GetNoiseSingle(q.x, q.y, q.z);
  HS_EXPECT_EQ(SphereNoise::sample(d, single), expected);
  generator.SetSeed(-12);
  HS_EXPECT_EQ(SphereNoise::sample(d, single), expected);
  generator.SetSeed(719);
  params.variant = NoiseVariant::OCTAVE_CONTOUR;
  params.basis = NoiseBasis::FBM3;
  params.contrast = 2.0f;
  const float octaves =
      (generator.GetNoiseSingle(q.x, q.y, q.z) +
       0.5f * generator.GetNoiseSingle(2.0f * q.x, 2.0f * q.y, 2.0f * q.z) +
       0.25f * generator.GetNoiseSingle(4.0f * q.x, 4.0f * q.y, 4.0f * q.z)) /
      1.75f;
  HS_EXPECT_NEAR(
      SphereNoise::sample(d, SphereNoise::prepare(params, STATE, generator)),
      octaves * 3.0f / (1.0f + 2.0f * std::abs(octaves)), 1.0e-5f);
  params.variant = NoiseVariant::CONFIGURED_SINGLE;
  HS_EXPECT_TRUE(!SphereNoise::valid(params));
}

/** @brief Plane bands preserve frequency and rotate counterclockwise. */
inline void test_pattern_plane_bands() {
  using Source = Pattern::PlaneBands;
  static_assert(std::is_empty_v<Source::State>);
  Source::Params params;
  params.pattern_freq = 2.0f;
  const auto horizontal = Source::prepare(params, {});
  for (const float x : {-2.0f, 0.0f, 1.3f}) {
    HS_EXPECT_NEAR(Source::sample({{x, 0.3f}}, horizontal), std::sin(0.6f),
                   0.0018f);
    HS_EXPECT_NEAR(Source::sample({{x, -0.4f}}, horizontal), std::sin(-0.8f),
                   0.0018f);
  }
  params.angle = 0.5f * math::PI_F;
  const auto vertical = Source::prepare(params, {});
  for (const float y : {-1.1f, 0.0f, 0.7f})
    HS_EXPECT_NEAR(Source::sample({{0.4f, y}}, vertical), std::sin(-0.8f),
                   0.0018f);
  params.angle = 0.6f;
  const auto diagonal = Source::prepare(params, {});
  HS_EXPECT_NEAR(
      Source::sample({{0.4f, -0.7f}}, diagonal),
      std::sin(2.0f * (-0.7f * std::cos(0.6f) - 0.4f * std::sin(0.6f))),
      0.0018f);
  params.angle = 0.0f;
  params.pattern_freq = 9.0f;
  HS_EXPECT_NEAR(
      Source::sample({{0.4f, -0.7f}}, diagonal),
      std::sin(2.0f * (-0.7f * std::cos(0.6f) - 0.4f * std::sin(0.6f))),
      0.0018f);
  params.angle = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_TRUE(!Source::valid(params));
}

/** @brief Plane interference advances its clocks and samples the rotated two-wave
 * formula; prepared state ignores later param edits. */
inline void test_pattern_plane_source_protocols() {
  using namespace Pattern;
  PlaneInterference::Params waves;
  waves.speed = -0.2f;
  waves.secondary_rate = 1.5f;
  waves.angle_rate = 0.03f;
  const auto zero = PlaneInterference::init_state(waves);
  const auto next = PlaneInterference::advance(zero, waves);
  HS_EXPECT_EQ(zero.primary, 0.0f);
  HS_EXPECT_EQ(next.primary, -0.2f);
  HS_EXPECT_EQ(next.secondary, -0.3f);
  HS_EXPECT_EQ(next.angle, 0.03f);
  const PlanePosition POINT{{0.6f, 0.25f}};
  const auto prepared = PlaneInterference::prepare(waves, next);
  const float rotated =
      POINT.value.re * std::cos(0.03f) + POINT.value.im * std::sin(0.03f);
  HS_EXPECT_NEAR(PlaneInterference::sample(POINT, prepared),
                 0.5f * (std::sin(0.4f) + std::sin(rotated - 0.3f)), 0.0018f);
  const float retained = PlaneInterference::sample(POINT, prepared);
  waves.pattern_freq = 13.0f;
  HS_EXPECT_EQ(PlaneInterference::sample(POINT, prepared), retained);
}

/** @brief Plane spiral rounds arms to an integer and samples the arm-phase formula. */
inline void test_pattern_plane_source_spiral_samples() {
  using namespace Pattern;
  const PlanePosition POINT{{0.6f, 0.25f}};
  static_assert(std::is_empty_v<PlaneSpiral::State>);
  PlaneSpiral::Params spiral_params;
  spiral_params.arms = 2.6f;
  const auto spiral = PlaneSpiral::prepare(spiral_params, {});
  HS_EXPECT_EQ(spiral.arms, 3.0f);
  HS_EXPECT_NEAR(PlaneSpiral::sample({{0.0f, 0.5f}}, spiral),
                 std::sin(0.5f - 1.5f * std::numbers::pi_v<float>), 0.012f);
  spiral_params.arms = 0.0f;
  HS_EXPECT_NEAR(
      PlaneSpiral::sample(POINT, PlaneSpiral::prepare(spiral_params, {})),
      std::sin(0.65f), 0.0018f);
}

/** @brief Coupled grid samples sin(x) times cos(y) at full mix. */
inline void test_pattern_plane_source_grid_samples() {
  using namespace Pattern;
  const PlanePosition POINT{{0.6f, 0.25f}};
  PlaneCoupledGrid::Params grid_params;
  grid_params.pattern_mix = 1.0f;
  HS_EXPECT_NEAR(PlaneCoupledGrid::sample(
                     POINT, PlaneCoupledGrid::prepare(grid_params, {})),
                 std::sin(0.6f) * std::cos(0.25f), 0.0035f);
}

/** @brief A soft square cell peaks at its centre, crosses zero at the half
 *  cell and is odd about it. */
inline void test_pattern_plane_source_tessellation_samples() {
  using namespace Pattern;
  using Cells =
      PlaneTessellation<TessellationKind::SQUARE, TessellationFill::CELLS>;
  static_assert(std::is_empty_v<Cells::State>);
  Cells::Params lattice;
  const auto empty = Cells::advance(Cells::init_state(lattice), lattice);
  HS_EXPECT_EQ(Cells::sample({{}}, Cells::prepare(lattice, empty)), 1.0f);
  lattice.cell_shape = 1.0f;
  lattice.cell_softness = 0.5f;
  const auto square = Cells::prepare(lattice, empty);
  const float centre = Cells::sample({{}}, square);
  HS_EXPECT_LE(centre, 1.0f);
  for (float sign : {-1.0f, 1.0f}) {
    const float inside = Cells::sample({{sign * 0.125f, 0.0f}}, square);
    HS_EXPECT_GT(inside, 0.0f);
    HS_EXPECT_LT(inside, centre);
    HS_EXPECT_EQ(inside, Cells::sample({{0.125f, 0.0f}}, square));
    HS_EXPECT_EQ(Cells::sample({{sign * 0.25f, 0.0f}}, square), 0.0f);
    HS_EXPECT_EQ(Cells::sample({{sign * 0.375f, 0.0f}}, square), -inside);
  }
  using Edges =
      PlaneTessellation<TessellationKind::SQUARE, TessellationFill::EDGES>;
  HS_EXPECT_EQ(Edges::sample({{}}, Edges::prepare({}, {})), -1.0f);
}

/** @brief Plane noise advances its turn clock and samples projected contour noise. */
inline void test_pattern_plane_source_noise_samples() {
  using namespace Pattern;
  const PlanePosition POINT{{0.6f, 0.25f}};
  const FastNoiseLite generator;
  PlaneNoise::Params noise;
  noise.noise_time_rate = -0.01f;
  const auto clock = PlaneNoise::advance(PlaneNoise::init_state(noise), noise);
  HS_EXPECT_EQ(clock.phase_turns, 0.99f);
  const auto noise_value =
      PlaneNoise::prepare(noise, clock, generator, math::NoiseBasis::FBM3);
  const auto q = math::noise_projected_coordinate(
      POINT.value, noise.noise_scale, clock.phase_turns);
  HS_EXPECT_EQ(PlaneNoise::sample(POINT, noise_value),
               noise_contour(generator, math::NoiseBasis::FBM3, q, 0.0f));
}

/** @brief Spherical cells use the engine's canonical seed layouts. */
inline void test_pattern_sphere_cell_layouts() {
  using namespace Pattern;
  struct LayoutCase {
    SphereCellLayout layout;
    std::span<const math::Vector> vertices;
  };
  const LayoutCase cases[] = {
      {SphereCellLayout::TETRAHEDRON, Solids::Tetrahedron::vertices},
      {SphereCellLayout::CUBE, Solids::Cube::vertices},
      {SphereCellLayout::OCTAHEDRON, Solids::Octahedron::vertices},
      {SphereCellLayout::DODECAHEDRON, Solids::Dodecahedron::vertices},
      {SphereCellLayout::ICOSAHEDRON, Solids::Icosahedron::vertices}};
  SphereCellsParams params;
  for (const auto &entry : cases) {
    params.layout = static_cast<uint8_t>(entry.layout);
    const auto prepared = SphereCells::prepare(params);
    HS_EXPECT_EQ(prepared.count, entry.vertices.size());
    for (size_t index = 0; index < entry.vertices.size(); ++index) {
      HS_EXPECT_EQ(prepared.centers[index].x, entry.vertices[index].x);
      HS_EXPECT_EQ(prepared.centers[index].y, entry.vertices[index].y);
      HS_EXPECT_EQ(prepared.centers[index].z, entry.vertices[index].z);
    }
  }
  params.layout = static_cast<uint8_t>(SphereCellLayout::FIBONACCI);
  params.count = 9.6f;
  const auto ten = SphereCells::prepare(params);
  HS_EXPECT_EQ(ten.count, 10u);
  for (size_t index = 0; index < ten.count; ++index) {
    const auto expected = math::fib_spiral(static_cast<int>(ten.count), 0.5f,
                                           static_cast<int>(index));
    HS_EXPECT_EQ(ten.centers[index].x, expected.x);
    HS_EXPECT_EQ(ten.centers[index].y, expected.y);
    HS_EXPECT_EQ(ten.centers[index].z, expected.z);
  }
  params.count = 24.0f;
  const auto denser = SphereCells::prepare(params);
  HS_EXPECT_EQ(denser.count, 24u);
  HS_EXPECT_TRUE(denser.centers[0].y != ten.centers[0].y);
  HS_EXPECT_GT(denser.distance_scale, ten.distance_scale);
  params.count = static_cast<float>(SphereCells::MAX_CELLS);
  HS_EXPECT_EQ(SphereCells::prepare(params).count, SphereCells::MAX_CELLS);
}

/** @brief Scattered cells are reproducible without consuming the effect RNG. */
inline void test_pattern_sphere_cell_seed_isolation() {
  using namespace Pattern;
  SphereCellsParams params;
  params.layout = static_cast<uint8_t>(SphereCellLayout::SCATTERED);
  params.count = 31.0f;
  params.seed = 917.0f;
  auto expected_rng = hs::random();
  const auto first = SphereCells::prepare(params);
  const auto again = SphereCells::prepare(params);
  ++params.seed;
  const auto other = SphereCells::prepare(params);
  auto actual_rng = hs::random();
  HS_EXPECT_EQ(actual_rng(), expected_rng());
  bool changed = false;
  for (size_t index = 0; index < first.count; ++index) {
    HS_EXPECT_EQ(first.centers[index].x, again.centers[index].x);
    HS_EXPECT_EQ(first.centers[index].y, again.centers[index].y);
    HS_EXPECT_EQ(first.centers[index].z, again.centers[index].z);
    HS_EXPECT_NEAR(math::dot(first.centers[index], first.centers[index]), 1.0f,
                   3.0e-7f);
    changed |= first.centers[index].x != other.centers[index].x;
  }
  HS_EXPECT_TRUE(changed);
  params.layout = static_cast<uint8_t>(SphereCellLayout::OCTAHEDRON);
  params.output = static_cast<uint8_t>(SphereCellOutput::CELL_VALUE);
  const auto a = SphereCells::prepare(params);
  ++params.seed;
  const auto b = SphereCells::prepare(params);
  const auto point = SphereDirection::from_unit(math::X_AXIS);
  HS_EXPECT_TRUE(SphereCells::sample(point, a) !=
                 SphereCells::sample(point, b));
}

/** @brief Cell borders, center falloff and labels match known octahedral cells. */
inline void test_pattern_sphere_cell_outputs() {
  using namespace Pattern;
  SphereCellsParams params;
  params.layout = static_cast<uint8_t>(SphereCellLayout::OCTAHEDRON);
  params.line_thickness = 0.1f;
  params.line_softness = 0.2f;
  auto prepared = SphereCells::prepare(params);
  const auto x = SphereDirection::from_unit(math::X_AXIS);
  const auto boundary = SphereDirection::from({1.0f, 1.0f, 0.0f});
  HS_EXPECT_EQ(SphereCells::sample(x, prepared), -1.0f);
  HS_EXPECT_EQ(SphereCells::sample(boundary, prepared), 1.0f);
  const float scale = std::sqrt(6.0f / (4.0f * std::numbers::pi_v<float>));
  const float offset = 0.2f / scale;
  for (float side : {-1.0f, 1.0f}) {
    const float angle = std::numbers::pi_v<float> / 4.0f + side * offset;
    const auto point =
        SphereDirection::from({std::cos(angle), std::sin(angle), 0});
    HS_EXPECT_NEAR(SphereCells::sample(point, prepared), 0.0f, 2.0e-6f);
  }
  prepared.output = SphereCellOutput::CENTERS;
  HS_EXPECT_NEAR(SphereCells::sample(x, prepared), 1.0f, 1.0e-6f);
  for (float angle : {0.2f, 0.5f}) {
    const auto point =
        SphereDirection::from({std::cos(angle), std::sin(angle), 0});
    HS_EXPECT_NEAR(SphereCells::sample(point, prepared),
                   1.0f - 2.0f * angle * scale, 5.0e-4f);
  }
  prepared.output = SphereCellOutput::CELL_VALUE;
  prepared.seed = 391;
  HS_EXPECT_EQ(SphereCells::sample(x, prepared),
               2.0f * math::hash01(0, 391) - 1.0f);
  HS_EXPECT_EQ(SphereCells::sample(boundary, prepared),
               2.0f * math::hash01(0, 391) - 1.0f);
  HS_EXPECT_EQ(
      SphereCells::sample(SphereDirection::from_unit(math::Y_AXIS), prepared),
      2.0f * math::hash01(2, 391) - 1.0f);
}

/** @brief The closest cell border can belong to a seed beyond the nearest two. */
inline void test_pattern_sphere_cell_all_bisectors() {
  using namespace Pattern;
  SphereCells::Prepared prepared;
  prepared.count = 4;
  prepared.distance_scale = 1.0f;
  prepared.line_thickness = 0.0f;
  prepared.line_softness = 0.1f;
  const float ANGLES[] = {0.2f, 0.21f, -0.3f, 2.0f};
  for (size_t index = 0; index < prepared.count; ++index)
    prepared.centers[index] = {std::sin(ANGLES[index]), std::cos(ANGLES[index]),
                               0};
  const auto point = SphereDirection::from_unit(math::Y_AXIS);
  HS_EXPECT_NEAR(SphereCells::sample(point, prepared), 0.0f, 2.0e-6f);
}

/** @brief Cell recipes accept finite endpoints and reject invalid controls. */
inline void test_pattern_sphere_cell_admission() {
  using namespace Pattern;
  SphereCellsParams params;
  params.count = 4.0f;
  params.seed = 65535.0f;
  params.line_thickness = 0.0f;
  params.line_softness = 1.0f / 1024.0f;
  HS_EXPECT_TRUE(SphereCells::valid(params));
  for (auto member : {&SphereCellsParams::count, &SphereCellsParams::seed,
                      &SphereCellsParams::line_thickness,
                      &SphereCellsParams::line_softness}) {
    for (float invalid : {std::numeric_limits<float>::quiet_NaN(),
                          std::numeric_limits<float>::infinity(), -1.0f}) {
      auto rejected = params;
      rejected.*member = invalid;
      HS_EXPECT_TRUE(!SphereCells::valid(rejected));
    }
  }
  params.count = 129.0f;
  HS_EXPECT_TRUE(!SphereCells::valid(params));
  params.count = 42.0f;
  params.layout = 255;
  HS_EXPECT_TRUE(!SphereCells::valid(params));
  params.layout = 0;
  params.output = 255;
  HS_EXPECT_TRUE(!SphereCells::valid(params));
}

/** @brief Regular unjittered web: eight spokes, four turns, no sag or dew. */
inline Pattern::SpiderWebParams regular_spider_web() {
  Pattern::SpiderWebParams params;
  params.web_radius = 1.0f;
  params.spokes = 8.0f;
  params.turns = 4.0f;
  params.hub_radius = 0.2f;
  params.sag = 0.0f;
  params.jitter = 0.0f;
  params.dew = 0.0f;
  params.line_thickness = 0.005f;
  params.line_softness = 0.005f;
  return params;
}

/** @brief Plane point at a polar radius and angle. */
inline Pattern::PlanePosition web_point(float radius, float angle) {
  return {{radius * std::cos(angle), radius * std::sin(angle)}};
}

/** @brief Spokes and capture chords light; gaps, the hub zone and the far
 *  plane stay dark. */
inline void test_pattern_spider_web_threads() {
  using Web = Pattern::PlaneSpiderWeb;
  const auto web = Web::prepare(regular_spider_web(), {});
  const float quarter = 0.25f * math::PI_F;
  for (const float radius : {0.35f, 0.55f, 0.75f})
    HS_EXPECT_NEAR(Web::sample(web_point(radius, 2.0f * quarter), web), 1.0f,
                   1e-4f);
  const auto a = web_point(0.4f, 0.0f).value;
  const auto b = web_point(0.425f, quarter).value;
  HS_EXPECT_NEAR(
      Web::sample({{0.5f * (a.re + b.re), 0.5f * (a.im + b.im)}}, web), 1.0f,
      1e-4f);
  HS_EXPECT_EQ(Web::sample(web_point(0.473f, 0.5f * quarter), web), -1.0f);
  HS_EXPECT_EQ(Web::sample(web_point(0.14f, 0.5f * quarter), web), -1.0f);
  HS_EXPECT_EQ(Web::sample({{3.0f, 0.0f}}, web), -1.0f);

  auto hubless = regular_spider_web();
  hubless.hub_radius = 0.0f;
  HS_EXPECT_EQ(
      Web::sample(web_point(0.36f, 0.5f * quarter), Web::prepare(hubless, {})),
      -1.0f);
}

/** @brief Sag bows a capture chord's middle toward the hub. */
inline void test_pattern_spider_web_sag() {
  using Web = Pattern::PlaneSpiderWeb;
  auto params = regular_spider_web();
  params.sag = 0.3f;
  const auto web = Web::prepare(params, {});
  const auto a = web_point(0.4f, 0.0f).value;
  const auto b = web_point(0.425f, 0.25f * math::PI_F).value;
  const float dx = b.re - a.re, dy = b.im - a.im;
  const float mid_x = 0.5f * (a.re + b.re), mid_y = 0.5f * (a.im + b.im);
  const math::Complex bowed{mid_x - params.sag * dy, mid_y + params.sag * dx};
  HS_EXPECT_NEAR(Web::sample({bowed}, web), 1.0f, 1e-4f);
  HS_EXPECT_EQ(Web::sample({{mid_x, mid_y}}, web), -1.0f);
}

/** @brief Bowed capture chords stay inside their spoke sector and remain lit. */
inline void test_pattern_spider_web_sag_stays_in_sector() {
  using Web = Pattern::PlaneSpiderWeb;
  for (const float spokes : {3.0f, 4.0f, 5.0f, 8.0f}) {
    auto params = regular_spider_web();
    params.spokes = spokes;
    params.sag = 0.5f;
    const auto web = Web::prepare(params, {});
    const auto a = web_point(0.4f, 0.0f).value;
    const auto b = web_point(0.4f + web.step, math::TWO_PI_F / spokes).value;
    const float dx = b.re - a.re, dy = b.im - a.im;
    for (const float u : {0.01f, 0.25f, 0.5f, 0.75f, 0.99f}) {
      const float bow = web.sag * 4.0f * u * (1.0f - u);
      const math::Complex point{a.re + u * dx - bow * dy,
                                a.im + u * dy + bow * dx};
      HS_EXPECT_TRUE(point.im >= -1e-7f);
      HS_EXPECT_TRUE(b.im * point.re - b.re * point.im >= -1e-7f);
      HS_EXPECT_NEAR(Web::sample({point}, web), 1.0f, 1e-4f);
    }
  }
}

/** @brief A spoke's edge sits at one angular offset at every radius under
 *  gnomonic scaling. */
inline void test_pattern_spider_web_angular_width() {
  using Web = Pattern::PlaneSpiderWeb;
  const auto params = regular_spider_web();
  const auto web = Web::prepare(params, {});
  const float offset = params.line_thickness + 0.5f * params.line_softness;
  for (const float radius : {0.35f, 0.75f}) {
    const float beta =
        std::asin(offset * std::sqrt(1.0f + radius * radius) / radius);
    HS_EXPECT_NEAR(
        Web::sample(web_point(radius, 0.5f * math::PI_F - beta), web), 0.0f,
        2e-3f);
  }
}

/** @brief Thick threads fade to darkness before the web's outer cutoff. */
inline void test_pattern_spider_web_outer_fade() {
  using Web = Pattern::PlaneSpiderWeb;
  auto params = Pattern::SpiderWebParams{};
  params.line_thickness = params.line_softness = 0.1f;
  params.dew = 0.0f;
  const auto web = Web::prepare(params, {});
  const float inside = std::nextafter(web.reach, 0.0f);
  for (const float angle :
       {web.angle[0], 0.5f * (web.angle[0] + web.angle[1])}) {
    HS_EXPECT_EQ(Web::sample(web_point(inside, angle), web), -1.0f);
    HS_EXPECT_EQ(Web::sample(web_point(web.reach * 1.01f, angle), web), -1.0f);
  }
  params.web_radius = 16.0f;
  const auto horizon = Web::prepare(params, {});
  HS_EXPECT_EQ(horizon.reach, FLT_MAX);
  HS_EXPECT_TRUE(std::isfinite(Web::sample(web_point(100.0f, 0.0f), horizon)));
}

/** @brief Dew only adds coverage, and adds some. */
inline void test_pattern_spider_web_dew() {
  using Web = Pattern::PlaneSpiderWeb;
  auto params = regular_spider_web();
  const auto dry = Web::prepare(params, {});
  params.dew = 1.0f;
  const auto wet = Web::prepare(params, {});
  float gained = 0.0f;
  for (int i = 0; i < 64; ++i)
    for (int j = 0; j < 64; ++j) {
      const Pattern::PlanePosition point{
          {-1.1f + 2.2f * i / 63.0f, -1.1f + 2.2f * j / 63.0f}};
      const float before = Web::sample(point, dry);
      const float after = Web::sample(point, wet);
      HS_EXPECT_TRUE(after >= before);
      gained += after - before;
    }
  HS_EXPECT_TRUE(gained > 1.0f);
}

/** @brief Preparation admits the defaults and rejects out-of-range recipes. */
inline void test_pattern_spider_web_admission() {
  using Web = Pattern::PlaneSpiderWeb;
  HS_EXPECT_TRUE(Web::valid({}));
  auto params = Pattern::SpiderWebParams{};
  params.spokes = static_cast<float>(Web::MAX_SPOKES) + 1.0f;
  HS_EXPECT_TRUE(!Web::valid(params));
  params = {};
  params.hub_radius = 0.95f;
  HS_EXPECT_TRUE(!Web::valid(params));
  params = {};
  params.line_thickness = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_TRUE(!Web::valid(params));
}

} // namespace hs_test::pattern_tests
