/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Cached modulation test cases. */

#include "tests/pullback/color/fixture.h"

namespace hs_test::pullback_tests {

/** @brief Fused cached bases retain canonical hue, opacity and inactive fallbacks. */
inline void test_pullback_fused_cached_modulation() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  TemporalPalettes<BakedPalette::required_arena_bytes()> palettes;
  palettes.bake(CoordinatePalette{1.0f},
                ConstantPalette{Color4(Pixel(32000, 44000, 13000), 1.0f)});
  auto &a = palettes.a;
  auto &b = palettes.b;
  auto &arena = palettes.arena;
  BakedPaletteStorage translucent;
  translucent.bake(arena,
                   ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.6f)});
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  PaletteFrame frame;
  palettes.bind(frame);
  Color::ModulationState modulation{
      nullptr,
      nullptr,
      Pattern::SourceKind::DISPLACEMENT,
      {ColorModulation::Mechanic::HUE_ROTATION, 1.5958333f},
      true};
  int cache_hits = 0;
  for (bool active : {false, true})
    for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f})
      for (float value : {0.0f, 0.13f, 0.5f, 0.99f, 1.0f}) {
        frame.lab = active ? &cache : nullptr;
        frame.temporal.fraction = fraction;
        const auto prepared = Policy::prepare(frame);
        const FieldSample field{value, 0.45f, math::Y_AXIS, 0.37f};
        auto pair = Policy::apply_pair(field, frame, prepared);
        Color::apply_modulation_in_place(pair.first, modulation, true);
        Color::apply_modulation_in_place(pair.second, modulation, false);
        const auto canonical =
            blend_color_samples(pair.first, pair.second, fraction);
        const auto actual = Policy::apply_fused_modulation(
            field, frame, prepared, modulation, fraction);
        auto expected = canonical;
        if (actual.has_cached_lab) {
          ++cache_hits;
          HS_EXPECT_TRUE(active && fraction > 0.0f && fraction < 1.0f);
          HS_EXPECT_TRUE(
              Policy::apply_shared_base(field, frame, prepared, expected));
          expected.pending_hue = canonical.pending_hue;
          expected.pending_lift = canonical.pending_lift;
          expected.color.alpha = canonical.color.alpha;
          expected.brightness_gain = canonical.brightness_gain;
        }
        HS_EXPECT_EQ(actual.pending_hue, canonical.pending_hue);
        HS_EXPECT_EQ(actual.pending_lift, canonical.pending_lift);
        HS_EXPECT_EQ(actual.color.alpha, canonical.color.alpha);
        HS_EXPECT_EQ(actual.brightness_gain, canonical.brightness_gain);
        const auto finished = Color::finish(actual);
        const auto expected_value = Color::finish(expected);
        HS_EXPECT_EQ(finished.color, expected_value.color);
        HS_EXPECT_EQ(finished.alpha, expected_value.alpha);
        const auto direct = Policy::shade_fused_modulation(
            field, frame, prepared, modulation, fraction);
        HS_EXPECT_EQ(direct.color, finished.color);
        HS_EXPECT_EQ(direct.alpha, finished.alpha);
      }
  HS_EXPECT_GT(cache_hits, 0);
  cache.bake(a.view(), translucent.view());
  frame.palette = &translucent.view();
  frame.lab = &cache;
  frame.temporal.fraction = 0.5f;
  const auto prepared = Policy::prepare(frame);
  const auto fallback = Policy::apply_fused_modulation(
      FieldSample{0.5f, 0.45f, math::Y_AXIS, 0.37f}, frame, prepared,
      modulation, 0.5f);
  HS_EXPECT_FALSE(fallback.has_cached_lab);
}

/** @brief Direct Lab shading reuses endpoint taps through canonical gamut fallback. */
inline void test_pullback_fused_lab_direct_fallback() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1626, 1306, 43097), 1.0f)},
                ConstantPalette{Color4(Pixel(1451, 1516, 42130), 1.0f)});
  auto &a = palettes.a;
  auto &b = palettes.b;
  PaletteLabLut cache;
  cache.active = true;
  std::array<Pixel, HueRotationLutView::SIZE> hue;
  hue.fill(Pixel(111, 222, 333));
  PaletteFrame frame;
  frame.palette = &b.view();
  frame.lab = &cache;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  static int calls;
  const auto source = [](const Modulation::Prepared &,
                         const Pattern::SphereDirection &) {
    ++calls;
    return 0.9814453125f;
  };
  Modulation::Prepared previous{}, current{};
  previous.evaluate = current.evaluate = source;
  frame.previous_modulation = &previous;
  frame.modulation = &current;
  const auto modulation = PaletteLayerProvider::prepare(frame);
  for (bool guard : {false, true}) {
    cache.nodes.fill(guard ? OKLab{0.431782782f, 0.00217834837f, -0.24565056f}
                           : linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f));
    for (bool same_base : {false, true}) {
      frame.temporal.palette = same_base ? &b.view() : &a.view();
      frame.hue_rotation = frame.temporal.hue_rotation =
          same_base ? hue.data() : nullptr;
      for (float fraction : {0.25f, 0.5f, 0.75f}) {
        frame.temporal.fraction = fraction;
        const auto prepared = Policy::prepare(frame);
        const FieldSample FIELD{0.5f, 0.61689651f, math::Y_AXIS, 0.37f};
        const auto carrier = Policy::apply_fused_modulation(
            FIELD, frame, prepared, modulation, fraction);
        HS_EXPECT_EQ(carrier.has_cached_lab, !guard);
        const auto expected = finish_color_sample(carrier);
        calls = 0;
        const auto actual = Policy::shade_fused_modulation(
            FIELD, frame, prepared, modulation, fraction);
        HS_EXPECT_EQ(actual.color, expected.color);
        HS_EXPECT_EQ(actual.alpha, expected.alpha);
        HS_EXPECT_EQ(calls, 2);
        if (guard && same_base)
          HS_EXPECT_EQ(actual.color, hue.front());
      }
    }
  }
}

/** @brief Cached current endpoints sample one live source and retain exact alpha. */
inline void test_pullback_fused_lab_current_endpoint() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<Policy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(10000, 40000, 20000), 1.0f)},
                ConstantPalette{Color4(Pixel(1451, 1516, 42130), 1.0f)});
  auto &a = palettes.a;
  auto &b = palettes.b;
  PaletteLabLut cache;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.temporal.fraction = 1.0f;
  frame.lab = &cache;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  static int calls;
  const auto source = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  Modulation::Prepared current{};
  current.evaluate = source;
  frame.modulation = &current;
  for (bool guard : {false, true}) {
    cache.bake(a.view(), b.view());
    if (guard)
      cache.nodes.fill({0.431782782f, 0.00217834837f, -0.24565056f});
    for (float scalar : {-0.125f, 0.0f, 0.13f, 0.9814453125f}) {
      current.state.phase.phase_turns = scalar;
      const auto modulation = PaletteLayerProvider::prepare(frame);
      const auto prepared = Policy::prepare(frame);
      for (float coverage : {0.0f, 0.5f, 1.0f}) {
        const FieldSample field{0.5f, coverage, math::Y_AXIS, 0.0f};
        const auto reference = Policy::apply_fused_modulation(
            field, frame, prepared, modulation, 1.0f);
        const auto expected = finish_color_sample(reference);
        auto cached = reference;
        cached.cached_lab =
            PaletteLabLut::sample(prepared.lab, cached.palette_coordinate);
        cached.has_cached_lab = true;
        Color::restore_cached_palette_base_near_gamut_gate(cached, [] {});
        if (guard && scalar == 0.9814453125f)
          HS_EXPECT_FALSE(cached.has_cached_lab);
        const auto cached_value = finish_color_sample(cached);
        calls = 0;
        const auto actual = Policy::shade_fused_modulation(
            field, frame, prepared, modulation, 1.0f);
        HS_EXPECT_EQ(actual.color, cached_value.color);
        HS_EXPECT_EQ(actual.alpha, expected.alpha);
        HS_EXPECT_EQ(calls, 1);
        if (!guard) {
          HS_EXPECT_NEAR(actual.color.r, expected.color.r, 2);
          HS_EXPECT_NEAR(actual.color.g, expected.color.g, 2);
          HS_EXPECT_NEAR(actual.color.b, expected.color.b, 2);
        }
      }
      calls = 0;
      const auto pixel =
          Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
      HS_EXPECT_EQ(calls, 1);
      const auto direct = Policy::shade_fused_modulation(
          FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.0f}, frame, prepared,
          modulation, 1.0f);
      HS_EXPECT_EQ(pixel.color, direct.color);
      HS_EXPECT_EQ(pixel.alpha, direct.alpha);
    }
  }
  frame.temporal.fraction = 0.0f;
  frame.previous_modulation = &current;
  const auto prepared_value = Policy::prepare(frame);
  const auto modulation = PaletteLayerProvider::prepare(frame);
  const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.0f};
  const auto expected = finish_color_sample(Policy::apply_fused_modulation(
      FIELD, frame, prepared_value, modulation, 0.0f));
  calls = 0;
  const auto actual = Policy::shade_fused_modulation(
      FIELD, frame, prepared_value, modulation, 0.0f);
  HS_EXPECT_EQ(actual.color, expected.color);
  HS_EXPECT_EQ(actual.alpha, expected.alpha);
  HS_EXPECT_EQ(calls, 1);
}

/** @brief Gamut-gate crossings restore RGB-first color without normal cache fetches. */
inline void test_pullback_cached_gamut_gate_fallback() {
  using namespace Pullback;
  const float NAN_VALUE = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_TRUE(
      Color::near_palette_gamut_gate({NAN_VALUE, NAN_VALUE, NAN_VALUE}));
  ColorSample canonical{Color4(Pixel(1582, 1359, 42855), 0.61689651f), {}};
  canonical.pending_hue = 0.9814453125f;
  auto cached = canonical;
  cached.color.color = Pixel(0, 0, 0);
  cached.cached_lab = {0.431782782f, 0.00217834837f, -0.24565056f};
  cached.has_cached_lab = true;
  const auto reference = finish_color_sample(canonical);
  const auto unguarded = finish_color_sample(cached);
  auto direct = cached;
  int direct_fetches = 0;
  const auto direct_value = Color::finish_fused_color(direct, [&] {
    ++direct_fetches;
    direct.color.color = canonical.color.color;
    direct.hue_rotation = canonical.hue_rotation;
  });
  HS_EXPECT_EQ(direct_fetches, 1);
  HS_EXPECT_FALSE(direct.has_cached_lab);
  HS_EXPECT_EQ(direct_value.color, reference.color);
  HS_EXPECT_EQ(direct_value.alpha, reference.alpha);
  const Pixel expected_value = reference.color * reference.alpha;
  const Pixel wrong = unguarded.color * unguarded.alpha;
  HS_EXPECT_GT(std::max({std::abs(int(expected_value.r) - wrong.r),
                         std::abs(int(expected_value.g) - wrong.g),
                         std::abs(int(expected_value.b) - wrong.b)}),
               4000);
  int fetches = 0;
  HS_EXPECT_TRUE(
      Color::restore_cached_palette_base_near_gamut_gate(cached, [&] {
        ++fetches;
        cached.color.color = canonical.color.color;
        cached.hue_rotation = canonical.hue_rotation;
      }));
  HS_EXPECT_EQ(fetches, 1);
  HS_EXPECT_FALSE(cached.has_cached_lab);
  const auto guarded = finish_color_sample(cached);
  HS_EXPECT_EQ(guarded.color, reference.color);
  HS_EXPECT_EQ(guarded.alpha, reference.alpha);
  cached.cached_lab = linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f);
  cached.has_cached_lab = true;
  cached.pending_hue = 0.0f;
  fetches = 0;
  HS_EXPECT_FALSE(Color::restore_cached_palette_base_near_gamut_gate(
      cached, [&] { ++fetches; }));
  HS_EXPECT_EQ(fetches, 0);
  HS_EXPECT_TRUE(cached.has_cached_lab);
  for (float hue : {0.0f, 0.13f, 0.5f, 0.9814453125f})
    for (float lift : {0.0f, 0.4f, 1.0f})
      for (float gain : {0.0f, 0.3f, 1.0f})
        for (auto lab : {linear_rgb_to_oklab_fast(0.3f, 0.4f, 0.5f),
                         OKLab{0.5f, 0.3f, -0.3f}}) {
          cached.cached_lab = lab;
          cached.pending_hue = hue;
          cached.pending_lift = lift;
          cached.brightness_gain = gain;
          cached.has_cached_lab = true;
          auto expected = cached;
          int reference_fetches = 0;
          Color::restore_cached_palette_base_near_gamut_gate(expected, [&] {
            ++reference_fetches;
            expected.color.color = canonical.color.color;
          });
          const auto finished = finish_color_sample(expected);
          int actual_fetches = 0;
          const auto fast = Color::finish_fused_color(cached, [&] {
            ++actual_fetches;
            cached.color.color = canonical.color.color;
          });
          HS_EXPECT_EQ(fast.color, finished.color);
          HS_EXPECT_EQ(fast.alpha, finished.alpha);
          HS_EXPECT_EQ(actual_fetches, reference_fetches);
          HS_EXPECT_EQ(cached.has_cached_lab, expected.has_cached_lab);
        }
}

/** @brief Prepared cube-pair eligibility preserves evaluator and active-source boundaries. */
inline void test_pullback_prepared_cube_pair() {
  using namespace Pullback;
  using Layer = Color::Modulate<PaletteLayerProvider>;
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  constexpr size_t ALIGNMENT = alignof(Color::ModulationState);
  constexpr size_t SIZE_THROUGH_ENABLED =
      (offsetof(Color::ModulationState, enabled) + sizeof(bool) + ALIGNMENT -
       1) /
      ALIGNMENT * ALIGNMENT;
  static_assert(sizeof(Color::ModulationState) == SIZE_THROUGH_ENABLED);
  CubeModulationPair cubes;
  auto &a = cubes.previous;
  auto &b = cubes.current;
  std::array<float, Pattern::SphereLutView::SIZE> linear;
  std::copy(cubes.current_cube.begin(), cubes.current_cube.end(),
            linear.begin());
  TemporalPalettes<> palettes;
  palettes.bake(CoordinatePalette{0.2f},
                ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.8f)});
  PaletteFrame frame;
  palettes.bind(frame);
  frame.previous_modulation = &a;
  frame.modulation = &b;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.5958333f};
  for (bool floats : {false, true}) {
    b.cube.linear_data = floats ? linear.data() : nullptr;
    const auto prepared = Layer::prepare(frame);
    HS_EXPECT_TRUE(prepared.paired_cube_ready);
    auto dynamic = prepared;
    dynamic.paired_cube_ready = false;
    for (float fraction : {0.00001f, 0.25f, 0.75f, 0.99999f}) {
      frame.temporal.fraction = fraction;
      const auto palette = Policy::prepare(frame);
      for (int x = -3; x <= 3; ++x) {
        const FieldSample field{0.13f, 0.47f, math::Vector(x, 0.7f, -1.0f),
                                0.3f};
        const auto expected = Policy::apply_fused_modulation(
            field, frame, palette, dynamic, fraction);
        const auto actual = Policy::apply_fused_modulation(
            field, frame, palette, prepared, fraction);
        HS_EXPECT_EQ(actual.pending_hue, expected.pending_hue);
        HS_EXPECT_EQ(actual.color.color, expected.color.color);
        HS_EXPECT_EQ(actual.color.alpha, expected.color.alpha);
        const auto finished = Policy::shade_fused_modulation(
            field, frame, palette, prepared, fraction);
        const auto reference = finish_color_sample(expected);
        HS_EXPECT_EQ(finished.color, reference.color);
        HS_EXPECT_EQ(finished.alpha, reference.alpha);
      }
      const auto pipeline_frame = TemporalPalettePipeline::prepare(frame);
      HS_EXPECT_TRUE(std::get<2>(pipeline_frame.prepared).paired_cube_ready);
    }
  }
  static int calls;
  const auto custom_value = [](const Modulation::Prepared &prepared,
                               const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  a.state.phase.phase_turns = -0.3f;
  b.state.phase.phase_turns = 0.7f;
  for (int custom = 1; custom <= 3; ++custom) {
    a.evaluate = custom & 1 ? custom_value : Modulation::sample_cube;
    b.evaluate = custom & 2 ? custom_value : Modulation::sample_cube;
    const auto prepared = Layer::prepare(frame);
    HS_EXPECT_FALSE(prepared.paired_cube_ready);
    calls = 0;
    const auto palette = Policy::prepare(frame);
    (void)Policy::shade_fused_modulation(
        FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.0f}, frame, palette, prepared,
        frame.temporal.fraction);
    HS_EXPECT_EQ(calls, custom == 3 ? 2 : 1);
  }
  a.evaluate = b.evaluate = Modulation::sample_cube;
  for (int inactive = 1; inactive <= 3; ++inactive) {
    a.cube.active = !(inactive & 1);
    b.cube.active = !(inactive & 2);
    HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  }
  a.cube.active = b.cube.active = true;
  frame.modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  frame.modulation = &b;
  frame.previous_modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  frame.modulation = nullptr;
  HS_EXPECT_FALSE(Layer::prepare(frame).paired_cube_ready);
  const Color::ModulationState manual{&b, &a, Pattern::SourceKind::NOISE,
                                      frame.mechanic, true};
  HS_EXPECT_FALSE(manual.paired_cube_ready);
}

/** @brief Direct endpoint sharing respects prepared generators and opaque samplers. */
inline void test_direct_modulation_sample_identity() {
  using namespace Pullback::Modulation;
  Config config;
  config.source = Pattern::SourceKind::NOISE;
  config.evaluation = Pattern::Evaluation::DIRECT;
  Resources<false> resources;
  Prepared first = resources.prepare(config, init(config));
  auto second = first;
  second.config.mechanic.amount = first.config.mechanic.amount + 0.125f;
  HS_EXPECT_TRUE(same_direct_sample(first, second));
  auto &generator = std::get<0>(second.source).noise;
  const auto original = generator;
  generator.SetDomainWarpAmp(0.25f);
  generator.SetFractalOctaves(5);
  HS_EXPECT_TRUE(same_direct_sample(first, second));
  generator.SetFrequency(0.5f);
  HS_EXPECT_FALSE(same_direct_sample(first, second));
  generator = original;
  generator.SetRotationType3D(FastNoiseLite::RotationType3D_ImproveXYPlanes);
  HS_EXPECT_FALSE(same_direct_sample(first, second));
  generator = original;
  std::get<0>(second.source).noise.SetSeed(config.seed + 1);
  HS_EXPECT_FALSE(same_direct_sample(first, second));
  HS_EXPECT_FALSE(same_endpoint_sample(first, second));
  first.evaluate = [](const Prepared &source,
                      const Pattern::SphereDirection &) {
    return source.cube.data[0] * (1.0f / 127.0f);
  };
  const int8_t low = -96, high = 96;
  first.cube.data = &low;
  static_assert(std::is_trivially_copyable_v<Prepared>);
  std::memcpy(&second, &first, sizeof(second));
  HS_EXPECT_TRUE(same_direct_sample(first, second));
  second.cube.data = &high;
  const auto direction = Pattern::SphereDirection::from_unit(math::Y_AXIS);
  HS_EXPECT_TRUE(sample(first, direction) < sample(second, direction));
  HS_EXPECT_FALSE(same_direct_sample(first, second));
  HS_EXPECT_FALSE(same_endpoint_sample(first, second));
  CubeModulationPair cubes;
  HS_EXPECT_FALSE(same_endpoint_sample(cubes.previous, cubes.current));
  cubes.previous_cube = cubes.current_cube;
  HS_EXPECT_TRUE(same_endpoint_sample(cubes.previous, cubes.current));
}

/** @brief Paired cube samples preserve individual samplers and custom dispatch. */
inline void test_pullback_modulation_paired_samples() {
  using namespace Pullback;
  CubeModulationPair cubes;
  auto &previous = cubes.previous;
  auto &current = cubes.current;
  std::array<float, Pattern::SphereLutView::SIZE> linear_cube;
  std::copy(cubes.current_cube.begin(), cubes.current_cube.end(),
            linear_cube.begin());
  for (bool floats : {false, true}) {
    current.cube.linear_data = floats ? linear_cube.data() : nullptr;
    for (int x = -3; x <= 3; ++x)
      for (int y = -3; y <= 3; ++y)
        for (int z = -3; z <= 3; ++z) {
          if (x == 0 && y == 0 && z == 0)
            continue;
          const auto direction =
              Pattern::SphereDirection::admit(math::Vector(x, y, z));
          const auto pair =
              Modulation::sample_pair(previous, current, direction);
          HS_EXPECT_EQ(pair.first, Modulation::sample(previous, direction));
          HS_EXPECT_EQ(pair.second, Modulation::sample(current, direction));
        }
  }
  static int calls;
  calls = 0;
  const auto custom = [](const Modulation::Prepared &prepared,
                         const Pattern::SphereDirection &) {
    ++calls;
    return prepared.state.phase.phase_turns;
  };
  previous.evaluate = current.evaluate = custom;
  previous.state.phase.phase_turns = -2.0f;
  current.state.phase.phase_turns = 3.0f;
  const auto pair = Modulation::sample_pair(
      previous, current, Pattern::SphereDirection::from_unit(math::Y_AXIS));
  HS_EXPECT_EQ(calls, 2);
  HS_EXPECT_EQ(pair.first, -1.0f);
  HS_EXPECT_EQ(pair.second, 1.0f);
}

/** @brief Paired modulation keeps canonical endpoint colors and mechanic order. */
inline void test_pullback_modulation_paired_colors() {
  using namespace Pullback;
  CubeModulationPair cubes;
  auto &previous = cubes.previous;
  auto &current = cubes.current;
  Color::ModulationState state{
      &current, &previous, Pattern::SourceKind::NOISE, {}, true};
  for (auto mechanic : {ColorModulation::Mechanic::HUE_ROTATION,
                        ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                        ColorModulation::Mechanic::OPACITY})
    for (float amount : {0.0f, 0.3f, 1.5f})
      for (bool same_direction : {false, true})
        for (int x = -2; x <= 2; ++x) {
          if (amount > 1.0f &&
              mechanic != ColorModulation::Mechanic::HUE_ROTATION)
            continue;
          state.mechanic = {mechanic, amount};
          ColorSample a{Color4(Pixel(41000, 17000, 33000), 0.25f),
                        {0.3f, 0.4f, math::Vector(x, 1.0f, 0.3f), 0.7f}};
          ColorSample b{
              Color4(Pixel(27000, 53000, 9000), 0.75f),
              {0.6f, 0.8f,
               same_direction ? a.field.sphere : math::Vector(0.3f, x, 1.0f),
               0.2f}};
          a.pending_hue = 0.93f;
          b.pending_hue = 0.07f;
          a.pending_lift = 0.2f;
          b.pending_lift = 0.4f;
          const auto expected_a = Color::apply_modulation(a, state, true);
          const auto expected_b = Color::apply_modulation(b, state, false);
          Color::apply_modulation_pair_in_place(a, b, state);
          for (auto pair :
               {std::pair{a, expected_a}, std::pair{b, expected_b}}) {
            HS_EXPECT_EQ(pair.first.color.color, pair.second.color.color);
            HS_EXPECT_EQ(pair.first.color.alpha, pair.second.color.alpha);
            HS_EXPECT_EQ(pair.first.pending_hue, pair.second.pending_hue);
            HS_EXPECT_EQ(pair.first.pending_lift, pair.second.pending_lift);
            HS_EXPECT_EQ(pair.first.field.sphere, pair.second.field.sphere);
            const auto finished = Color::finish(pair.first);
            const auto expected = Color::finish(pair.second);
            HS_EXPECT_EQ(finished.color, expected.color);
            HS_EXPECT_EQ(finished.alpha, expected.alpha);
          }
        }
}

/** @brief Modulation identities retain common and selected inputs only. */
inline void test_pullback_modulation_field_identity() {
  using namespace Pullback::Modulation;
  using Pattern::SourceKind;
  Config config;
  config.version = 7;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.noise.scale = 2.0f;
  config.sine.frequency = 3.0f;
  config.interference.frequency_b = 4.0f;
  config.grid.offset_b = 0.5f;
  config.seed = 123;
  const State STATE{SourceKind::NOISE, 1, {0.25f}, true};
  for (SourceKind source :
       {SourceKind::NOISE, SourceKind::DIRECTIONAL_SINE,
        SourceKind::SPHERE_INTERFERENCE, SourceKind::SPHERE_SMOOTH_GRID,
        SourceKind::DISPLACEMENT}) {
    config.source = source;
    const FieldIdentity identity = field_identity(config, STATE);
    HS_EXPECT_EQ(identity.source, source);
    HS_EXPECT_EQ(identity.version, 7);
    HS_EXPECT_EQ(identity.evaluation, Pattern::Evaluation::DIRECT);
    HS_EXPECT_EQ(identity.phase, Pattern::PhaseState{0.25f});
    HS_EXPECT_EQ(identity.bake_layout_version,
                 Pattern::SphereLutView::LAYOUT_VERSION);
    HS_EXPECT_EQ(identity.seed, source == SourceKind::NOISE ? 123 : 0);
    HS_EXPECT_EQ(identity.noise, source == SourceKind::NOISE
                                     ? config.noise
                                     : Pattern::NoiseParams{});
    HS_EXPECT_EQ(identity.sine, source == SourceKind::DIRECTIONAL_SINE
                                    ? config.sine
                                    : Pattern::SineParams{});
    HS_EXPECT_EQ(identity.interference,
                 source == SourceKind::SPHERE_INTERFERENCE
                     ? config.interference
                     : Pattern::InterferenceParams{});
    HS_EXPECT_EQ(identity.grid, source == SourceKind::SPHERE_SMOOTH_GRID
                                    ? config.grid
                                    : Pattern::GridParams{});
  }
  config.source = static_cast<SourceKind>(255);
  const FieldIdentity unknown = field_identity(config, STATE);
  HS_EXPECT_EQ(unknown.source, static_cast<SourceKind>(255));
  HS_EXPECT_EQ(unknown.version, 7);
  HS_EXPECT_EQ(unknown.evaluation, Pattern::Evaluation::DIRECT);
  HS_EXPECT_EQ(unknown.phase, STATE.phase);
  HS_EXPECT_EQ(unknown.noise, Pattern::NoiseParams{});
  HS_EXPECT_EQ(unknown.sine, Pattern::SineParams{});
  HS_EXPECT_EQ(unknown.interference, Pattern::InterferenceParams{});
  HS_EXPECT_EQ(unknown.grid, Pattern::GridParams{});
  HS_EXPECT_EQ(unknown.seed, 0);
  HS_EXPECT_EQ(unknown.bake_layout_version,
               Pattern::SphereLutView::LAYOUT_VERSION);
}

/** @brief Modulation equality distinguishes field inputs from dormant controls. */
inline void test_pullback_modulation_field_equality() {
  using namespace Pullback::Modulation;
  using Pattern::SourceKind;
  for (SourceKind source :
       {SourceKind::NOISE, SourceKind::DIRECTIONAL_SINE,
        SourceKind::SPHERE_INTERFERENCE, SourceKind::SPHERE_SMOOTH_GRID,
        SourceKind::DISPLACEMENT}) {
    Prepared baseline{};
    baseline.config.source = source;
    baseline.config.evaluation = Pattern::Evaluation::DIRECT;
    const auto expect = [&](auto change, bool equal) {
      Prepared other = baseline;
      change(other);
      HS_EXPECT_EQ(same_field(baseline, other), equal);
      HS_EXPECT_EQ(same_field(other, baseline), equal);
      HS_EXPECT_EQ(field_identity(baseline.config, baseline.state) ==
                       field_identity(other.config, other.state),
                   equal);
    };
    HS_EXPECT_TRUE(same_field(baseline, baseline));
    expect(
        [&](auto &p) {
          p.config.source = source == SourceKind::NOISE
                                ? SourceKind::DISPLACEMENT
                                : SourceKind::NOISE;
        },
        false);
    expect([](auto &p) { ++p.config.version; }, false);
    expect(
        [](auto &p) { p.config.evaluation = Pattern::Evaluation::BAKED_CUBE; },
        false);
    expect([](auto &p) { p.state.phase.phase_turns = 0.25f; }, false);
    expect([](auto &p) { ++p.config.seed; }, source != SourceKind::NOISE);
    expect([](auto &p) { p.config.noise.scale = 2.0f; },
           source != SourceKind::NOISE);
    expect([](auto &p) { p.config.noise.speed = 0.5f; },
           source != SourceKind::NOISE);
    expect([](auto &p) { p.config.sine.frequency = 2.0f; },
           source != SourceKind::DIRECTIONAL_SINE);
    expect([](auto &p) { p.config.sine.axis.x = 0.25f; },
           source != SourceKind::DIRECTIONAL_SINE);
    expect([](auto &p) { p.config.interference.frequency = 2.0f; },
           source != SourceKind::SPHERE_INTERFERENCE);
    expect([](auto &p) { p.config.interference.axis_b.z = 0.25f; },
           source != SourceKind::SPHERE_INTERFERENCE);
    expect([](auto &p) { p.config.grid.frequency = 2.0f; },
           source != SourceKind::SPHERE_SMOOTH_GRID);
    expect([](auto &p) { p.config.grid.offset_b = 0.25f; },
           source != SourceKind::SPHERE_SMOOTH_GRID);
    expect([](auto &p) { p.config.enabled = false; }, true);
    expect(
        [](auto &p) {
          p.config.mechanic.kind = ColorModulation::Mechanic::OPACITY;
        },
        true);
    expect([](auto &p) { p.config.mechanic.amount = 0.5f; }, true);
    expect([](auto &p) { p.state.source = static_cast<SourceKind>(255); },
           true);
    expect([](auto &p) { ++p.state.version; }, true);
    expect([](auto &p) { p.state.initialized = true; }, true);
  }
  Prepared unknown{};
  unknown.config.source = static_cast<SourceKind>(255);
  HS_EXPECT_FALSE(same_field(unknown, unknown));
}

} // namespace hs_test::pullback_tests
