/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Qualification test cases. */

#include "tests/pullback/color/fixture.h"

namespace hs_test::pullback_tests {

/** @brief Prepared hue modulation meets channel limits across cycling palettes. */
inline void test_pullback_modulation_lut_qualification() {
  using namespace Pullback;
  using Layer = Color::Modulate<DisplacementLayerProvider>;
  std::array<uint8_t, GeneratedPaletteBank::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  GeneratedPaletteBank bank;
  bank.init(arena, 0.788f, math::ease_in_out_sin);
  std::array<Pixel, HueRotationLutView::SIZE> table;
  HueRotationGamutMask mask;
  PaletteFrame frame;
  frame.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 1.0f};
  const auto prepared = Layer::prepare(frame);
  uint32_t maximum = 0;
  int worst_frame = 0, worst_value = 0, worst_hue = 0;
  uint64_t total = 0;
  size_t channels = 0;
  int elapsed = 0;
  for (const int target : {0, 1, 150, 300, 600, 900, 1200, 1800}) {
    while (elapsed < target) {
      bank.step(Cycle::TRIADIC);
      ++elapsed;
    }
    const auto &palette = bank.palette(Cycle::TRIADIC);
    prepare_hue_rotation_lut(std::span(table), palette, &mask);
    constexpr int VALUES = 512, HUES = 128;
    for (int value = 0; value < VALUES; ++value) {
      const float coordinate = (value + 0.5f) / VALUES;
      const auto base = palette.get(coordinate);
      for (int hue = 0; hue < HUES; ++hue) {
        const float scalar = 2.0f * (hue + 0.5f) / HUES - 1.0f;
        ColorSample input{base, {coordinate, 0.75f, math::Y_AXIS, scalar}};
        input.hue_rotation = table.data();
        input.hue_rotation_gamut = &mask;
        input.palette_coordinate = coordinate;
        const auto actual =
            finish_color_sample(Layer::apply(input, frame, prepared));
        const auto expected =
            ColorModulation::apply(base, scalar, frame.mechanic);
        HS_EXPECT_EQ(actual.alpha, expected.alpha);
        for (const uint32_t error :
             {channel_gap(actual.color.r, expected.color.r),
              channel_gap(actual.color.g, expected.color.g),
              channel_gap(actual.color.b, expected.color.b)}) {
          if (error > maximum) {
            maximum = error;
            worst_frame = target;
            worst_value = value;
            worst_hue = hue;
          }
          total += error;
          ++channels;
        }
      }
    }
  }
  const double mean = static_cast<double>(total) / channels;
  std::printf("  modulation LUT max=%u mean=%.3f channels=%zu worst=%d/%d/%d\n",
              maximum, mean, channels, worst_frame, worst_value, worst_hue);
  HS_EXPECT_GT(maximum, 0u);
  HS_EXPECT_LE(static_cast<float>(maximum), Layer::METRICS[0].limit);
  HS_EXPECT_LE(mean, Layer::METRICS[1].limit);
}

/** @brief Shared cached bases skip both palette endpoints and preserve alpha. */
inline void test_pullback_shared_lab_dispatch() {
  using namespace Pullback;
  using Pipeline = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                      Stage::Colorize<CountingLabPalette>,
                                      Stage::Modulate<CountingPhasePolicy>>;
  using Core =
      Detail::PipelineCore<true, PaletteBinding, typename Pipeline::NodeList>;
  using Custom = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                    Stage::Colorize<ShiftedLabPalette>,
                                    Stage::Modulate<CountingPhasePolicy>>;
  using CustomCore =
      Detail::PipelineCore<true, PaletteBinding, typename Custom::NodeList>;
  TemporalPalettes<BakedPalette::required_arena_bytes()> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 1.0f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 1.0f)});
  auto &a = palettes.a;
  auto &b = palettes.b;
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  Color::PhaseModulationState layer;
  layer.source = Pattern::SourceKind::DISPLACEMENT;
  layer.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.2f};
  layer.enabled = true;
  layer.ready = true;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.phase = &layer;
  frame.lab = &cache;
  frame.temporal.fraction = 0.37f;
  const auto direct_prepared = CountingLabPalette::Base::prepare(frame);
  const FieldSample DIRECT_FIELD{0.5f, 0.5f, math::Y_AXIS, 0.7f};
  ColorSample dirty{Color4(Pixel(123, 456, 789), 0.9f),
                    {0.1f, 0.2f, math::X_AXIS, 0.3f},
                    reinterpret_cast<const Pixel *>(&cache),
                    0.4f,
                    0.5f,
                    0.6f,
                    0.7f,
                    {0.8f, 0.9f, 1.0f},
                    false};
  const auto stale_rgb = cache.view(0.0f);
  dirty.cached_rgb = &stale_rgb;
  HS_EXPECT_TRUE(CountingLabPalette::Base::apply_shared_base(
      DIRECT_FIELD, frame, direct_prepared, dirty));
  HS_EXPECT_EQ(dirty.color.color, Pixel(0, 0, 0));
  const auto direct_pair = CountingLabPalette::Base::apply_pair(
      DIRECT_FIELD, frame, direct_prepared);
  HS_EXPECT_EQ(dirty.color.alpha, direct_pair.first.color.alpha);
  HS_EXPECT_EQ(dirty.field.value, DIRECT_FIELD.value);
  HS_EXPECT_EQ(dirty.field.coverage, DIRECT_FIELD.coverage);
  HS_EXPECT_EQ(dirty.field.sphere, DIRECT_FIELD.sphere);
  HS_EXPECT_EQ(dirty.field.path_length, DIRECT_FIELD.path_length);
  HS_EXPECT_TRUE(dirty.hue_rotation == nullptr);
  HS_EXPECT_TRUE(dirty.cached_rgb == nullptr);
  HS_EXPECT_EQ(dirty.palette_coordinate, 0.5f);
  HS_EXPECT_EQ(dirty.pending_hue, 0.0f);
  HS_EXPECT_EQ(dirty.brightness_gain, 1.0f);
  HS_EXPECT_EQ(dirty.pending_lift, 0.0f);
  const auto expected_lab = PaletteLabLut::sample(cache.view(0.37f), 0.5f);
  HS_EXPECT_EQ(dirty.cached_lab.L, expected_lab.L);
  HS_EXPECT_EQ(dirty.cached_lab.a, expected_lab.a);
  HS_EXPECT_EQ(dirty.cached_lab.b, expected_lab.b);
  HS_EXPECT_TRUE(dirty.has_cached_lab);
  for (bool ready : {false, true})
    for (bool active : {false, true})
      for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
        layer.ready = ready;
        frame.lab = active ? &cache : nullptr;
        frame.temporal.fraction = fraction;
        CountingLabPalette::shared_calls = 0;
        CountingLabPalette::pair_calls = 0;
        const auto prepared = Pipeline::prepare_stages(frame);
        const auto actual = Core::template run_stage<0>(
            SphereSample{math::Y_AXIS, 0.7f}, frame, prepared);
        const bool interior = fraction > 0.0f && fraction < 1.0f;
        const bool shared = ready && fraction > 0.0f;
        HS_EXPECT_EQ(actual.has_cached_lab, active && shared);
        HS_EXPECT_EQ(CountingLabPalette::shared_calls, shared ? 1 : 0);
        HS_EXPECT_EQ(CountingLabPalette::pair_calls,
                     interior && !(active && ready) ? 1 : 0);
        const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.7f};
        const auto pair = CountingLabPalette::Base::apply_pair(
            FIELD, frame, std::get<1>(prepared));
        HS_EXPECT_EQ(actual.color.alpha,
                     fraction <= 0.0f ? pair.first.color.alpha
                     : fraction >= 1.0f
                         ? pair.second.color.alpha
                         : hs::lerp(pair.first.color.alpha,
                                    pair.second.color.alpha, fraction));
        HS_EXPECT_NEAR(actual.pending_hue, 0.14f, 1e-7f);
        if (!interior && !(active && shared))
          HS_EXPECT_EQ(actual.color.color, fraction == 0.0f
                                               ? pair.first.color.color
                                               : pair.second.color.color);
      }
  layer.ready = true;
  frame.lab = &cache;
  frame.temporal.fraction = 0.5f;
  CountingLabPalette::shared_calls = 0;
  CountingLabPalette::pair_calls = 0;
  const auto custom = CustomCore::template run_stage<0>(
      SphereSample{math::Y_AXIS, 0.7f}, frame, Custom::prepare_stages(frame));
  HS_EXPECT_FALSE(custom.has_cached_lab);
  HS_EXPECT_EQ(CountingLabPalette::shared_calls, 0);
  HS_EXPECT_EQ(CountingLabPalette::pair_calls, 1);
  HS_EXPECT_EQ(custom.color.color.r, 5123);
  BakedPaletteStorage translucent;
  translucent.bake(palettes.arena,
                   ConstantPalette{Color4(Pixel(4000, 5000, 6000), 0.6f)});
  cache.bake(a.view(), translucent.view());
  frame.palette = &translucent.view();
  const auto prepared = Pipeline::prepare_stages(frame);
  const auto fallback = Core::template run_stage<0>(
      SphereSample{math::Y_AXIS, 0.7f}, frame, prepared);
  const auto pair = CountingLabPalette::Base::apply_pair(
      FieldSample{0.5f, 0.5f, math::Y_AXIS, 0.7f}, frame,
      std::get<1>(prepared));
  HS_EXPECT_FALSE(fallback.has_cached_lab);
  HS_EXPECT_EQ(fallback.color.alpha,
               hs::lerp(pair.first.color.alpha, pair.second.color.alpha, 0.5f));
}

/** @brief Shared phase tails run once; incomplete or unready tails keep endpoints. */
inline void test_pullback_phase_tail_dispatch() {
  using namespace Pullback;
  using Phase = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<CountingPhasePolicy>>;
  using Fallback = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<CountingPhasePolicy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  using PhaseCore =
      Detail::PipelineCore<true, PaletteBinding, typename Phase::NodeList>;
  using FallbackCore =
      Detail::PipelineCore<true, PaletteBinding, typename Fallback::NodeList>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  Color::PhaseModulationState layer;
  layer.source = Pattern::SourceKind::DISPLACEMENT;
  layer.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, 0.2f};
  layer.enabled = true;
  PaletteFrame frame;
  palettes.bind(frame);
  frame.phase = &layer;
  frame.mechanic.amount = 0.0f;
  for (bool ready : {false, true})
    for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
      layer.ready = ready;
      frame.temporal.fraction = fraction;
      CountingPhasePolicy::phase_calls = 0;
      CountingPhasePolicy::endpoint_calls = 0;
      const auto prepared = Phase::prepare_stages(frame);
      const auto sample = PhaseCore::template run_stage<0>(
          SphereSample{math::Y_AXIS, 0.7f}, frame, prepared);
      const bool shared = ready && fraction > 0.0f && fraction < 1.0f;
      HS_EXPECT_EQ(CountingPhasePolicy::phase_calls, shared ? 2 : 0);
      HS_EXPECT_EQ(CountingPhasePolicy::endpoint_calls,
                   shared ? 0 : (fraction > 0.0f && fraction < 1.0f ? 4 : 2));
      HS_EXPECT_NEAR(sample.pending_hue, 0.28f, 1e-7f);
      HS_EXPECT_NEAR(sample.color.alpha, 0.04f + 0.08f * fraction, 1e-7f);
      HS_EXPECT_NEAR(sample.color.color.r, 1000 + 8000 * fraction, 1.0f);
      CountingPhasePolicy::phase_calls = 0;
      CountingPhasePolicy::endpoint_calls = 0;
      const auto other = FallbackCore::template run_stage<0>(
          SphereSample{math::Y_AXIS, 0.7f}, frame,
          Fallback::prepare_stages(frame));
      HS_EXPECT_EQ(CountingPhasePolicy::phase_calls, 0);
      HS_EXPECT_EQ(CountingPhasePolicy::endpoint_calls,
                   fraction > 0.0f && fraction < 1.0f ? 4 : 2);
      HS_EXPECT_EQ(other.color.color, sample.color.color);
      HS_EXPECT_EQ(other.color.alpha, sample.color.alpha);
      HS_EXPECT_EQ(other.pending_hue, sample.pending_hue);
    }
}

} // namespace hs_test::pullback_tests
