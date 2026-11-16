/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Temporal modulation test cases. */

#include "tests/pullback/color/fixture.h"

namespace hs_test::pullback_tests {

/** @brief Paired color dispatch skips individual endpoint calls during blends. */
inline void test_pullback_paired_color_dispatch() {
  using namespace Pullback;
  using Paired = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                    Stage::Colorize<CountingPairedPalette>>;
  using Fallback = Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                                      Stage::Colorize<CountingEndpointPalette>>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  PaletteFrame frame;
  palettes.bind(frame);
  for (float fraction : {0.0f, 0.25f, 1.0f}) {
    frame.temporal.fraction = fraction;
    CountingPairedPalette::pairs = 0;
    CountingPairedPalette::endpoints = 0;
    const auto color = Paired::shade(math::Y_AXIS, Paired::prepare(frame));
    HS_EXPECT_NEAR(color.color.r, 1000 + 8000 * fraction, 1.0f);
    HS_EXPECT_NEAR(color.alpha, 0.04f + 0.08f * fraction, 1e-6f);
    HS_EXPECT_EQ(CountingPairedPalette::pairs, fraction == 0.25f ? 1 : 0);
    HS_EXPECT_EQ(CountingPairedPalette::endpoints, fraction == 0.25f ? 0 : 1);
    CountingPairedPalette::pairs = 0;
    CountingPairedPalette::endpoints = 0;
    const auto other = Fallback::shade(math::Y_AXIS, Fallback::prepare(frame));
    HS_EXPECT_EQ(other.color, color.color);
    HS_EXPECT_EQ(other.alpha, color.alpha);
    HS_EXPECT_EQ(CountingPairedPalette::pairs, 0);
    HS_EXPECT_EQ(CountingPairedPalette::endpoints,
                 fraction > 0.0f && fraction < 1.0f ? 2 : 1);
  }
}

/** @brief Paired palettes preserve shared controls and endpoint resources. */
inline void test_pullback_paired_palette_controls() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  TemporalPalettes<> palettes;
  palettes.bake(CoordinatePalette{0.2f}, CoordinatePalette{0.6f});
  auto &a = palettes.a;
  std::array<Pixel, 2> hue{};
  PaletteFrame frame;
  palettes.bind(frame);
  auto prepared = Policy::prepare(frame);
  prepared.current.mapping_frequency = 2.0f;
  prepared.current.mapping_offset = 0.13f;
  prepared.current.brightness_envelope = Color::BrightnessEnvelope::ASCENDING;
  prepared.current.brightness_bottom = 0.2f;
  prepared.current.brightness_top = 0.8f;
  prepared.current.hue_rotation = {hue.data(), true};
  prepared.previous = prepared.current;
  prepared.previous.palette = &a.view();
  prepared.previous.hue_rotation = {hue.data() + 1, true};
  const FieldSample INPUT{0.37f, 0.45f, math::X_AXIS, 2.5f};
  const auto pair = Policy::apply_pair(INPUT, frame, prepared);
  for (const auto *sample : {&pair.first, &pair.second}) {
    HS_EXPECT_NEAR(sample->palette_coordinate, 0.87f, 1e-6f);
    HS_EXPECT_NEAR(sample->brightness_gain, 0.422f, 1e-6f);
    HS_EXPECT_EQ(sample->field.value, INPUT.value);
    HS_EXPECT_EQ(sample->field.coverage, INPUT.coverage);
    HS_EXPECT_EQ(sample->field.sphere, INPUT.sphere);
    HS_EXPECT_EQ(sample->field.path_length, INPUT.path_length);
    HS_EXPECT_EQ(sample->pending_hue, 0.0f);
    HS_EXPECT_EQ(sample->pending_lift, 0.0f);
    HS_EXPECT_EQ(sample->color.color, a.view().get(0.87f).color);
  }
  HS_EXPECT_NEAR(pair.first.color.alpha, 0.03132f, 1e-6f);
  HS_EXPECT_NEAR(pair.second.color.alpha, 0.09396f, 1e-6f);
  HS_EXPECT_TRUE(pair.first.hue_rotation == hue.data() + 1);
  HS_EXPECT_TRUE(pair.second.hue_rotation == hue.data());
}

/** @brief The complete color suffix interpolates immutable palette endpoints. */
inline void test_pullback_temporal_palette() {
  using namespace Pullback;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(1000, 2000, 3000), 0.2f)},
                ConstantPalette{Color4(Pixel(9000, 6000, 1000), 0.6f)});
  PaletteFrame frame;
  palettes.bind(frame);
  frame.mechanic.amount = 0.0f;
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto color = TemporalPalettePipeline::shade(
        math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
    HS_EXPECT_NEAR(color.color.r, 1000 + 8000 * fraction, 1.0f);
    HS_EXPECT_NEAR(color.color.g, 2000 + 4000 * fraction, 1.0f);
    HS_EXPECT_NEAR(color.color.b, 3000 - 2000 * fraction, 1.0f);
    HS_EXPECT_NEAR(color.alpha, 0.04f + 0.08f * fraction, 1e-6f);
  }
}

/** @brief Opacity modulates each endpoint before interpolation, and midpoint lightness lifts. */
inline void test_pullback_modulation_temporal_alpha() {
  using namespace Pullback;
  const Pixel RGB(1234, 4567, 8910);
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(RGB, 0.2f)},
                ConstantPalette{Color4(RGB, 0.6f)});
  Modulation::Resources<false> previous, current;
  prepare_sine_clock_pair(previous, current);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.modulation = &current.prepared;
  frame.previous_modulation = &previous.prepared;
  frame.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto color = TemporalPalettePipeline::shade(
        math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
    HS_EXPECT_EQ(color.color, RGB);
    HS_EXPECT_NEAR(color.alpha, 0.12f * fraction, 1e-6f);
    const auto check_placed = [&]<CodeEmission Emission>() {
      using Placed = TemporalPalettePlacedPipeline<Emission>;
      const auto placed = Placed::shade(math::Y_AXIS, Placed::prepare(frame));
      HS_EXPECT_EQ(placed.color, color.color);
      HS_EXPECT_NEAR(placed.alpha, color.alpha, 0.0f);
    };
    check_placed.template operator()<CodeEmission::INLINE_ONLY>();
    check_placed.template operator()<CodeEmission::OUT_OF_LINE_FLASH>();
    check_placed.template operator()<CodeEmission::OUT_OF_LINE_ITCM>();
  }
  frame.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 1.0f};
  frame.temporal.fraction = 0.5f;
  const auto color = TemporalPalettePipeline::shade(
      math::Y_AXIS, TemporalPalettePipeline::prepare(frame));
  const auto expected =
      ColorModulation::PositiveLightnessLift::apply(Color4(RGB), 0.5f, {1.0f})
          .color;
  HS_EXPECT_NEAR(color.color.r, expected.r, 1.0f);
  HS_EXPECT_NEAR(color.color.g, expected.g, 1.0f);
  HS_EXPECT_NEAR(color.color.b, expected.b, 1.0f);
}

/** @brief Perceptual finishing applies brightness before RGB16 quantization. */
inline void test_pullback_finish_brightness_quantization() {
  using namespace Pullback;
  size_t quantization_differences = 0;
  for (uint16_t r : {0, 1, 2, 1000, 32768, 65535})
    for (uint16_t g : {0, 1, 2, 1000, 32768, 65535})
      for (uint16_t b : {0, 1, 2, 1000, 32768, 65535})
        for (float hue : {0.0f, 0.17f, 0.75f})
          for (float lift : {0.0f, 0.25f, 0.9f}) {
            if (hue == 0.0f && lift == 0.0f)
              continue;
            const Color4 base(Pixel(r, g, b), 0.37f);
            const auto undimmed =
                composed_modulation_reference(base, hue, lift);
            for (float gain : {-1.0f, 0.0f, 0.001f, 0.25f, 1.0f, 2.0f}) {
              ColorSample sample{base, {0.5f, 0.7f, math::Y_AXIS, 2.0f}};
              sample.pending_hue = hue;
              sample.pending_lift = lift;
              sample.brightness_gain = gain;
              const auto actual = finish_color_sample(sample);
              const auto expected =
                  composed_modulation_reference(base, hue, lift, gain);
              HS_EXPECT_EQ(actual.color, expected.color);
              HS_EXPECT_EQ(actual.alpha, base.alpha);
              const Pixel dimmed_after_quantization = undimmed.color * gain;
              const int gap = std::max(
                  {std::abs(int(actual.color.r) - dimmed_after_quantization.r),
                   std::abs(int(actual.color.g) - dimmed_after_quantization.g),
                   std::abs(int(actual.color.b) -
                            dimmed_after_quantization.b)});
              HS_EXPECT_LE(gap, gain <= 1.0f ? 1 : 2);
              if (gain == 0.0f)
                HS_EXPECT_EQ(actual.color, Pixel(0, 0, 0));
              if (gain == 1.0f)
                HS_EXPECT_EQ(actual.color, undimmed.color);
              quantization_differences +=
                  actual.color != dimmed_after_quantization;
            }
          }
  HS_EXPECT_GT(quantization_differences, size_t{0});
}

/** @brief Hue layers compose across lightness lift without RGB round trips. */
inline void test_pullback_stacked_modulations() {
  using namespace Pullback;
  const Color4 BASE(Pixel(23000, 17000, 11000), 0.7f);
  const ColorSample INPUT{BASE, {0.4f, 0.6f, math::Y_AXIS, 2.5f}};
  Color::ModulationState first{nullptr,
                               nullptr,
                               Pattern::SourceKind::DISPLACEMENT,
                               {ColorModulation::Mechanic::HUE_ROTATION, 0.3f},
                               true};
  auto second = first;
  second.mechanic.amount = -0.3f;
  const auto cancelled = Color::finish(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), second));
  HS_EXPECT_EQ(cancelled.color, BASE.color);
  HS_EXPECT_NEAR(cancelled.alpha, BASE.alpha, 0.0f);
  std::array<Pixel, HueRotationLutView::SIZE> table;
  for (size_t index = 0; index < table.size(); ++index)
    table[index] = Pixel((index % 16) * 1000, 2000, 3000);
  auto with_table = INPUT;
  with_table.hue_rotation = table.data();
  HS_EXPECT_EQ(
      Color::finish(Color::apply_modulation(
                        Color::apply_modulation(with_table, first), second))
          .color,
      BASE.color);
  Modulation::Config lift_config;
  lift_config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  lift_config.evaluation = Pattern::Evaluation::DIRECT;
  lift_config.sine.frequency = 0.0f;
  auto clock = Modulation::init(lift_config);
  Modulation::Resources<false> lift_source;
  for (float phase : {0.0f, 0.75f}) {
    clock.phase.phase_turns = phase;
    const auto &source = lift_source.prepare(lift_config, clock);
    const Color::ModulationState identity_lift{
        &source,
        &source,
        Pattern::SourceKind::DIRECTIONAL_SINE,
        {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 1.0f},
        true};
    const auto layered = Color::apply_modulation(
        Color::apply_modulation(Color::apply_modulation(with_table, first),
                                identity_lift),
        second);
    HS_EXPECT_EQ(Color::finish(layered).color, BASE.color);
  }
  const auto rotated = Color::finish(Color::apply_modulation(INPUT, first));
  const auto expected = ColorModulation::HueRotation::apply(BASE, 2.5f, {0.3f});
  HS_EXPECT_EQ(rotated.color, expected.color);
  auto lift = first;
  lift.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.4f};
  second.mechanic.amount = 0.1f;
  auto actual = Color::apply_modulation(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), lift),
      second);
  HS_EXPECT_EQ(actual.color.color, BASE.color);
  const auto composed = composed_modulation_reference(BASE, 0.0f, 0.4f);
  HS_EXPECT_EQ(Color::finish(actual).color, composed.color);
  const auto reordered = Color::finish(Color::apply_modulation(
      Color::apply_modulation(Color::apply_modulation(INPUT, first), second),
      lift));
  HS_EXPECT_EQ(reordered.color, composed.color);
  actual.brightness_gain = 0.25f;
  const auto dimmed = Color::finish(actual).color;
  const auto dimmed_reference =
      composed_modulation_reference(BASE, 0.0f, 0.4f, 0.25f).color;
  HS_EXPECT_EQ(dimmed, dimmed_reference);
  auto opacity = first;
  opacity.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  HS_EXPECT_NEAR(Color::finish(Color::apply_modulation(INPUT, opacity)).alpha,
                 BASE.alpha, 0.0f);
  HS_EXPECT_EQ(
      Color::finish(Color::apply_modulation(INPUT, lift)).color,
      ColorModulation::PositiveLightnessLift::apply(BASE, 1.0f, {0.4f}).color);
}

/** @brief Saturated lifted colors rotate before final clipping and opacity. */
inline void test_pullback_composed_lightness_hue() {
  using namespace Pullback;
  const Color::ModulationState HUE_A{
      nullptr,
      nullptr,
      Pattern::SourceKind::DISPLACEMENT,
      {ColorModulation::Mechanic::HUE_ROTATION, 0.25f},
      true};
  auto hue_b = HUE_A;
  hue_b.mechanic.amount = 0.5f;
  auto lift_a = HUE_A;
  lift_a.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.8f};
  auto lift_b = lift_a;
  lift_b.mechanic.amount = 0.5f;
  auto opacity = HUE_A;
  opacity.mechanic = {ColorModulation::Mechanic::OPACITY, 0.5f};
  for (const Pixel BASE : {Pixel(65535, 0, 0), Pixel(0, 65535, 0),
                           Pixel(0, 0, 65535), Pixel(62000, 61000, 60000)}) {
    const ColorSample input{Color4(BASE, 0.7f),
                            {0.5f, 1.0f, math::Y_AXIS, 0.5f}};
    auto actual = Color::apply_modulation(input, HUE_A);
    actual = Color::apply_modulation(actual, lift_a);
    actual = Color::apply_modulation(actual, hue_b);
    actual = Color::apply_modulation(actual, lift_b);
    HS_EXPECT_EQ(actual.color.color, BASE);
    const auto expected =
        composed_modulation_reference(input.color, 0.375f, 0.55f);
    const auto finished = Color::finish(actual);
    HS_EXPECT_NEAR(finished.color.r, expected.color.r, 1.0f);
    HS_EXPECT_NEAR(finished.color.g, expected.color.g, 1.0f);
    HS_EXPECT_NEAR(finished.color.b, expected.color.b, 1.0f);
    HS_EXPECT_NEAR(finished.alpha, 0.7f, 0.0f);
    actual = Color::apply_modulation(actual, opacity);
    actual.brightness_gain = 0.25f;
    const auto masked = Color::finish(actual);
    HS_EXPECT_EQ(
        masked.color,
        composed_modulation_reference(input.color, 0.375f, 0.55f, 0.25f).color);
    HS_EXPECT_NEAR(masked.alpha, 0.6125f, 1e-7f);
    actual.brightness_gain = 0.0f;
    const auto dark = Color::finish(actual);
    HS_EXPECT_EQ(dark.color, Pixel(0, 0, 0));
    HS_EXPECT_EQ(dark.alpha, masked.alpha);
  }
}

/** @brief Temporal deferred color state blends before one perceptual finish. */
inline void test_pullback_temporal_composed_modulations() {
  using namespace Pullback;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage,
      Stage::Colorize<Color::GeneratedPalette<TemporalPaletteProvider>>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>,
      Stage::Modulate<Color::Modulate<PaletteHueLayerProvider>>>;
  const Color4 START(Pixel(65535, 0, 0), 0.2f);
  const Color4 END(Pixel(0, 0, 65535), 0.6f);
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{START}, ConstantPalette{END});
  Modulation::Resources<false> previous, current;
  prepare_sine_clock_pair(previous, current);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.modulation = &current.prepared;
  frame.previous_modulation = &previous.prepared;
  frame.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.4f};
  const auto previous_value =
      composed_modulation_reference(START, -0.25f, 0.0f);
  const auto current_value = composed_modulation_reference(END, 0.25f, 0.4f);
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    frame.temporal.fraction = fraction;
    const auto color = Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
    const auto base = START.color.lerp16(END.color, frac_to_q16(fraction));
    const auto expected = composed_modulation_reference(
        Color4(base), 0.75f - 0.5f * fraction, 0.4f * fraction);
    HS_EXPECT_NEAR(color.color.r, expected.color.r, 1.0f);
    HS_EXPECT_NEAR(color.color.g, expected.color.g, 1.0f);
    HS_EXPECT_NEAR(color.color.b, expected.color.b, 1.0f);
    HS_EXPECT_NEAR(color.alpha, 0.04f + 0.08f * fraction, 1e-6f);
    if (fraction > 0.0f && fraction < 1.0f) {
      const auto completed_blend = previous_value.color.lerp16(
          current_value.color, frac_to_q16(fraction));
      HS_EXPECT_GT(std::abs(int(color.color.r) - int(completed_blend.r)) +
                       std::abs(int(color.color.g) - int(completed_blend.g)) +
                       std::abs(int(color.color.b) - int(completed_blend.b)),
                   1000);
    }
  }
}

/** @brief Deferred hues take the short arc and preserve endpoint identities. */
inline void test_pullback_temporal_deferred_hue_wrap() {
  using namespace Pullback;
  ColorSample a{Color4(Pixel(65535, 0, 0), 0.2f), {}};
  ColorSample b{Color4(Pixel(0, 0, 65535), 0.8f), {}};
  a.pending_hue = 0.95f;
  b.pending_hue = 0.05f;
  a.pending_lift = 0.2f;
  b.pending_lift = 0.8f;
  a.brightness_gain = 0.25f;
  b.brightness_gain = 0.75f;
  const auto mixed = blend_color_samples(a, b, 0.5f);
  HS_EXPECT_NEAR(mixed.pending_hue, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(mixed.pending_lift, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mixed.brightness_gain, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(mixed.color.alpha, 0.5f, 1e-6f);
  HS_EXPECT_EQ(mixed.color.color,
               a.color.color.lerp16(b.color.color, frac_to_q16(0.5f)));
  HS_EXPECT_EQ(blend_color_samples(a, b, 0.0f).color.color, a.color.color);
  HS_EXPECT_EQ(blend_color_samples(a, b, 1.0f).color.color, b.color.color);
  HS_EXPECT_EQ(blend_color_samples(a, b, 0.0f).pending_hue, a.pending_hue);
  HS_EXPECT_EQ(blend_color_samples(a, b, 1.0f).pending_hue, b.pending_hue);
}

/** @brief In-place temporal blending preserves current sampling coordinates. */
inline void test_pullback_temporal_deferred_in_place() {
  using namespace Pullback;
  std::array<Pixel, 2> hue{};
  ColorSample a{Color4(Pixel(1000, 2000, 4000), 0.25f),
                {0.1f, 0.2f, math::X_AXIS, 0.3f}};
  ColorSample b{Color4(Pixel(9000, 6000, 12000), 0.75f),
                {0.6f, 0.7f, math::Y_AXIS, 0.8f}};
  a.hue_rotation = hue.data();
  b.hue_rotation = hue.data();
  a.palette_coordinate = 0.2f;
  b.palette_coordinate = 0.9f;
  a.pending_hue = 0.875f;
  b.pending_hue = 0.125f;
  a.pending_lift = 0.25f;
  b.pending_lift = 0.75f;
  a.brightness_gain = 0.5f;
  b.brightness_gain = 0.75f;
  blend_color_samples_in_place(a, b, 0.5f);
  HS_EXPECT_EQ(b.color.color, Pixel(5000, 4000, 8000));
  HS_EXPECT_EQ(b.color.alpha, 0.5f);
  HS_EXPECT_EQ(b.pending_hue, 0.0f);
  HS_EXPECT_EQ(b.pending_lift, 0.5f);
  HS_EXPECT_EQ(b.brightness_gain, 0.625f);
  HS_EXPECT_EQ(b.field.value, 0.6f);
  HS_EXPECT_EQ(b.field.coverage, 0.7f);
  HS_EXPECT_EQ(b.field.sphere, math::Y_AXIS);
  HS_EXPECT_EQ(b.field.path_length, 0.8f);
  HS_EXPECT_EQ(b.palette_coordinate, 0.9f);
  HS_EXPECT_TRUE(b.hue_rotation == nullptr);
  HS_EXPECT_EQ(a.color.color, Pixel(1000, 2000, 4000));
  HS_EXPECT_EQ(a.pending_hue, 0.875f);
  HS_EXPECT_TRUE(a.hue_rotation == hue.data());

  b = a;
  b.color.color.r = 1001;
  blend_color_samples_in_place(a, b, 0.000001f);
  HS_EXPECT_EQ(b.color.color, a.color.color);
  HS_EXPECT_TRUE(b.hue_rotation == nullptr);
  a.cached_lab = {0.4f, 0.1f, -0.2f};
  a.has_cached_lab = true;
  blend_color_samples_in_place(a, b, 0.0f);
  HS_EXPECT_TRUE(b.has_cached_lab);
  HS_EXPECT_EQ(b.cached_lab.L, 0.4f);
  HS_EXPECT_EQ(b.cached_lab.a, 0.1f);
  HS_EXPECT_EQ(b.cached_lab.b, -0.2f);
  HS_EXPECT_TRUE(b.hue_rotation == hue.data());
  blend_color_samples_in_place(a, b, 1.0f);
  HS_EXPECT_TRUE(b.has_cached_lab);
}

/** @brief A fused final layer matches the canonical temporal endpoint pipeline. */
inline void test_pullback_fused_modulation_endpoints() {
  using namespace Pullback;
  TemporalPalettes<> palettes;
  palettes.bake(CoordinatePalette{0.2f},
                ConstantPalette{Color4(Pixel(32000, 44000, 13000), 0.8f)});
  auto &a = palettes.a;
  PaletteFrame frame;
  palettes.bind(frame);
  using Policy = Color::GeneratedPalette<TemporalPaletteProvider>;
  /** @brief Non-fused palette policy used as the temporal reference. */
  struct CanonicalPalette : Policy {};
  using Actual = TemporalPalettePipeline;
  using Canonical = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<CanonicalPalette>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  using ActualCore =
      Detail::PipelineCore<true, PaletteBinding, typename Actual::NodeList>;
  using CanonicalCore =
      Detail::PipelineCore<true, PaletteBinding, typename Canonical::NodeList>;
  using Phase =
      Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                         Stage::Colorize<Policy>,
                         Stage::Modulate<Color::Modulate<PhaseLayerProvider>>>;
  using PhaseCore =
      Detail::PipelineCore<true, PaletteBinding, typename Phase::NodeList>;
  static_assert(
      ActualCore::template fused_color_tail_supported<1, FieldSample>());
  static_assert(
      !CanonicalCore::template fused_color_tail_supported<1, FieldSample>());
  static_assert(
      !PhaseCore::template fused_color_tail_supported<1, FieldSample>());

  std::array<int8_t, Pattern::SphereLutView::SIZE> ac, bc;
  for (size_t i = 0; i < ac.size(); ++i) {
    ac[i] = int(i * 37 % 255) - 127;
    bc[i] = int((i * 83 + 41) % 255) - 127;
  }
  Modulation::Prepared pa{}, pb{};
  pa.cube = {ac.data(), true};
  pb.cube = {bc.data(), true};
  pa.evaluate = pb.evaluate = Modulation::sample_cube;
  Color::ModulationState mod{&pb, &pa, Pattern::SourceKind::NOISE, {}, true};
  frame.modulation = &pb;
  frame.previous_modulation = &pa;
  for (auto mapping :
       {Color::PaletteMapping::LINEAR, Color::PaletteMapping::CUP,
        Color::PaletteMapping::BELL, Color::PaletteMapping::REVERSE})
    for (auto mechanic : {ColorModulation::Mechanic::HUE_ROTATION,
                          ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                          ColorModulation::Mechanic::OPACITY})
      for (float amount : {0.0f, 0.3f, 1.0f, 1.5958333f})
        for (float fraction :
             {0.0f, 0.000001f, 0.01f, 0.25f, 0.5f, 0.75f, 0.99999f, 1.0f})
          for (int v = 0; v <= 64; ++v) {
            if (amount > 1.0f &&
                mechanic != ColorModulation::Mechanic::HUE_ROTATION)
              continue;
            frame.temporal.fraction = fraction;
            auto prepared = Policy::prepare(frame);
            prepared.current.mapping =
                Color::PaletteMappingWeights::single(mapping);
            prepared.current.mapping_frequency = 2.0f;
            prepared.current.mapping_offset = 0.13f;
            prepared.current.brightness_envelope =
                Color::BrightnessEnvelope::ASCENDING;
            prepared.current.brightness_bottom = 0.2f;
            prepared.current.brightness_top = 0.8f;
            prepared.previous = prepared.current;
            prepared.previous.palette = &a.view();
            mod.mechanic = {mechanic, amount};
            frame.mechanic = mod.mechanic;
            const FieldSample field{v / 64.0f, 0.45f,
                                    math::Vector(v / 64.0f, 0.3f, 1.0f), 0.7f};
            auto pair = Policy::apply_pair(field, frame, prepared);
            Color::apply_modulation_in_place(pair.first, mod, true);
            Color::apply_modulation_in_place(pair.second, mod, false);
            const auto expected =
                blend_color_samples(pair.first, pair.second, fraction);
            const auto actual = Policy::apply_fused_modulation(
                field, frame, prepared, mod, fraction);
            HS_EXPECT_EQ(actual.color.color, expected.color.color);
            HS_EXPECT_EQ(actual.color.alpha, expected.color.alpha);
            HS_EXPECT_EQ(actual.pending_hue, expected.pending_hue);
            HS_EXPECT_EQ(actual.pending_lift, expected.pending_lift);
            HS_EXPECT_EQ(actual.brightness_gain, expected.brightness_gain);
            HS_EXPECT_TRUE(actual.hue_rotation == expected.hue_rotation);
            const auto x = Color::finish(actual), y = Color::finish(expected);
            HS_EXPECT_EQ(x.color, y.color);
            HS_EXPECT_EQ(x.alpha, y.alpha);
            const auto direct = Policy::shade_fused_modulation(
                field, frame, prepared, mod, fraction);
            HS_EXPECT_EQ(direct.color, y.color);
            HS_EXPECT_EQ(direct.alpha, y.alpha);
            const auto pixel =
                Actual::shade(field.sphere, Actual::prepare(frame));
            const auto reference =
                Canonical::shade(field.sphere, Canonical::prepare(frame));
            HS_EXPECT_EQ(pixel.color, reference.color);
            HS_EXPECT_EQ(pixel.alpha, reference.alpha);
          }
}

/** @brief Unit-hue wrapping matches canonical periodic interpolation bit for bit. */
inline void test_pullback_unit_hue_blending() {
  const std::array<float, 12> values{-0.0f,
                                     0.0f,
                                     std::numeric_limits<float>::denorm_min(),
                                     std::numeric_limits<float>::min(),
                                     0.125f,
                                     0.25f,
                                     std::nextafter(0.5f, 0.0f),
                                     0.5f,
                                     std::nextafter(0.5f, 1.0f),
                                     0.75f,
                                     0.875f,
                                     std::nextafter(1.0f, 0.0f)};
  const std::array<float, 12> fractions{
      -1.0f,
      -0.0f,
      0.0f,
      std::numeric_limits<float>::denorm_min(),
      std::numeric_limits<float>::min(),
      0.000001f,
      0.25f,
      0.5f,
      0.75f,
      std::nextafter(1.0f, 0.0f),
      1.0f,
      2.0f};
  const auto check = [](float a, float b, float fraction) {
    const float expected =
        a == b ? a : interp::shortest_periodic(a, b, fraction, 1.0f);
    const float actual = Pullback::blend_unit_hue(a, b, fraction);
    HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual),
                 std::bit_cast<uint32_t>(expected));
  };
  for (float a : values)
    for (float b : values)
      for (float fraction : fractions)
        check(a, b, fraction);
  uint32_t random = 0x6d2b79f5u;
  const auto next = [&] {
    random = random * 1664525u + 1013904223u;
    return float(random >> 8) * (1.0f / 16777216.0f);
  };
  for (int i = 0; i < 131072; ++i) {
    const float a = next();
    const float b = next();
    check(a, b, i % 8 == 0 ? fractions[(i / 8) % fractions.size()] : next());
  }
}

} // namespace hs_test::pullback_tests
