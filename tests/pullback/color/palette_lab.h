/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Palette lab test cases. */

#include "tests/pullback/color/fixture.h"

namespace hs_test::pullback_tests {

/** @brief Palette Lab nodes convert RGB-first blends and reject varying alpha. */
inline void test_pullback_palette_lab_nodes() {
  std::array<uint8_t, 5 * BakedPalette::required_arena_bytes()> buffer{};
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage a, b, variable;
  a.bake(arena, ConstantPalette{Color4(Pixel(0, 0, 0), 1.0f)});
  b.bake(arena, ConstantPalette{Color4(Pixel(65535, 65535, 65535), 1.0f)});
  variable.bake(arena, VariableAlphaPalette{});
  HS_EXPECT_LT(variable.view().get_alpha(1.0f / 255.0f), 1.0f);
  for (int index = 0; index < PaletteLabLut::VALUES; ++index)
    HS_EXPECT_EQ(
        variable.view().get_alpha(float(index) / (PaletteLabLut::VALUES - 1)),
        1.0f);
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  HS_EXPECT_TRUE(cache.active);
  const auto mid = PaletteLabLut::sample(cache.view(0.5f), 0.37f);
  HS_EXPECT_NEAR(mid.L, std::cbrt(32768.0 / 65535.0), 3e-5);
  HS_EXPECT_GT(mid.L, 0.75f);
  HS_EXPECT_NEAR(mid.a, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(mid.b, 0.0f, 1e-6f);
  for (int row = 0; row < PaletteLabLut::TIMES; ++row)
    for (int index = 0; index < PaletteLabLut::VALUES; ++index) {
      const float fraction = float(row) / (PaletteLabLut::TIMES - 1);
      const float value = float(index) / (PaletteLabLut::VALUES - 1);
      const Pixel rgb = a.view().get_color(value).lerp16(
          b.view().get_color(value), frac_to_q16(fraction));
      const auto linear = pixel_to_linrgb(rgb);
      const auto expected =
          linear_rgb_to_oklab_fast(linear.r, linear.g, linear.b);
      const auto actual = PaletteLabLut::sample(cache.view(fraction), value);
      HS_EXPECT_NEAR(actual.L, expected.L, 1e-7f);
      HS_EXPECT_NEAR(actual.a, expected.a, 1e-7f);
      HS_EXPECT_NEAR(actual.b, expected.b, 1e-7f);
    }
  const auto old_view = cache.view(0.25f);
  PaletteLabLut next;
  next.bake(b.view(), a.view());
  HS_EXPECT_NEAR(PaletteLabLut::sample(old_view, 0.37f).L, mid.L * 0.5f, 1e-7f);
  BakedPaletteStorage ramp, translucent;
  ramp.bake(arena, LabRampPalette{});
  translucent.bake(arena, ConstantPalette{Color4(Pixel(123, 456, 789), 0.6f)});
  next.bake(ramp.view(), b.view());
  for (float value : {0.183f, 0.677f}) {
    const float position =
        value * static_cast<float>(PaletteLabLut::VALUES - 1);
    const int index_value = static_cast<int>(position);
    const float weight = position - index_value;
    const auto node = [&](int index, float fraction) {
      const float coordinate =
          float(index) / static_cast<float>(PaletteLabLut::VALUES - 1);
      const Pixel rgb =
          ramp.view()
              .get_color(coordinate)
              .lerp16(b.view().get_color(coordinate), frac_to_q16(fraction));
      const LinRGB linear = pixel_to_linrgb(rgb);
      return linear_rgb_to_oklab_fast(linear.r, linear.g, linear.b);
    };
    const auto a_value = node(index_value, 0.0f),
               b_value = node(index_value + 1, 0.0f);
    const auto c = node(index_value, 0.5f), d = node(index_value + 1, 0.5f);
    const auto actual = PaletteLabLut::sample(next.view(0.17f), value);
    HS_EXPECT_NEAR(actual.L,
                   hs::lerp(hs::lerp(a_value.L, b_value.L, weight),
                            hs::lerp(c.L, d.L, weight),
                            (0.17f * (PaletteLabLut::TIMES - 1))),
                   1e-7f);
    HS_EXPECT_NEAR(actual.a,
                   hs::lerp(hs::lerp(a_value.a, b_value.a, weight),
                            hs::lerp(c.a, d.a, weight),
                            (0.17f * (PaletteLabLut::TIMES - 1))),
                   1e-7f);
    HS_EXPECT_NEAR(actual.b,
                   hs::lerp(hs::lerp(a_value.b, b_value.b, weight),
                            hs::lerp(c.b, d.b, weight),
                            (0.17f * (PaletteLabLut::TIMES - 1))),
                   1e-7f);
  }
  next.bake(ramp.view(), translucent.view());
  HS_EXPECT_FALSE(next.active);
  cache.bake(a.view(), variable.view());
  HS_EXPECT_FALSE(cache.active);
  HS_EXPECT_TRUE(cache.view(0.5f).low == nullptr);
}

/** @brief Owned palette copies rebind Lab nodes and retain external RGB bindings. */
inline void test_pullback_palette_lab_owned_copy() {
  const auto check = []<int Values>() {
    PaletteLabLutStorage<Values> source;
    for (size_t i = 0; i < source.nodes.size(); ++i)
      source.nodes[i] = {float(i) / source.nodes.size(), 0.1f, -0.2f};
    std::array<Pixel, HueRotationLutView::SIZE> rgb{};
    source.active = true;
    source.rgb_nodes = rgb.data();
    source.rgb_midpoint = true;
    source.rgb_midpoint_row = 0;
    const auto nodes = source.nodes;
    PaletteLabLutStorage<Values> copy(source);
    HS_EXPECT_TRUE(copy.data == copy.nodes.data() && copy.data != source.data);
    HS_EXPECT_EQ(copy.value_steps, Values);
    HS_EXPECT_TRUE(copy.active && copy.rgb_midpoint);
    HS_EXPECT_EQ(copy.rgb_midpoint_row, 0);
    HS_EXPECT_TRUE(copy.rgb_nodes == rgb.data());
    HS_EXPECT_EQ(std::memcmp(copy.nodes.data(), nodes.data(), sizeof(nodes)),
                 0);
    source.nodes[0].L = 1.0f;
    HS_EXPECT_EQ(copy.nodes[0].L, nodes[0].L);
    const auto view = copy.view(0.25f);
    HS_EXPECT_TRUE(view.low == copy.nodes.data());
    HS_EXPECT_EQ(view.value_steps, Values);
    HS_EXPECT_TRUE(view.rgb_low == rgb.data() && view.rgb_high == rgb.data());
    PaletteLabLutStorage<Values> assigned;
    HS_EXPECT_TRUE(&(assigned = copy) == &assigned);
    HS_EXPECT_TRUE(assigned.data == assigned.nodes.data());
    HS_EXPECT_TRUE(assigned.data != copy.data);
    HS_EXPECT_EQ(
        std::memcmp(assigned.nodes.data(), nodes.data(), sizeof(nodes)), 0);
    auto *same = &assigned;
    assigned = *same;
    HS_EXPECT_TRUE(assigned.data == assigned.nodes.data());
    HS_EXPECT_EQ(
        std::memcmp(assigned.nodes.data(), nodes.data(), sizeof(nodes)), 0);
    HS_EXPECT_TRUE(assigned.rgb_nodes == rgb.data() && assigned.rgb_midpoint);
    HS_EXPECT_EQ(assigned.rgb_midpoint_row, 0);
    copy.nodes[0].L = 0.75f;
    HS_EXPECT_EQ(assigned.nodes[0].L, nodes[0].L);
  };
  check.template operator()<64>();
  check.template operator()<256>();
}

/** @brief Cached neutral adjustments preserve alpha and endpoint blends keep the cached Lab. */
inline void test_pullback_cached_lab_finish() {
  using namespace Pullback;
  ColorSample sample{Color4(Pixel(0, 0, 0), 0.37f),
                     {0.5f, 0.7f, math::Y_AXIS, 2.0f}};
  sample.cached_lab = {0.7f, 0.02f, 0.03f};
  sample.has_cached_lab = true;
  const auto expected =
      linrgb_to_pixel(oklab_to_linear_rgb_lut_gamut(sample.cached_lab));
  const auto actual = finish_color_sample(sample);
  HS_EXPECT_EQ(actual.color, expected);
  HS_EXPECT_EQ(actual.alpha, 0.37f);
  HS_EXPECT_GT(actual.color.r, 0);
  sample.brightness_gain = 0.0f;
  HS_EXPECT_EQ(finish_color_sample(sample).color, Pixel(0, 0, 0));
  HS_EXPECT_EQ(finish_color_sample(sample).alpha, 0.37f);
  HS_EXPECT_TRUE(blend_color_samples(sample, sample, 0.0f).has_cached_lab);
  HS_EXPECT_TRUE(blend_color_samples(sample, sample, 1.0f).has_cached_lab);
}

/** @brief RGB-bound fused carriers retain the shared RGB and materialized-Lab finishers. */
inline void test_pullback_fused_rgb_carrier_finish() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Pipeline = Pullback::Pipeline<
      PaletteBinding, PaletteSampleStage, Stage::Colorize<Policy>,
      Stage::Modulate<Color::Modulate<PaletteLayerProvider>>>;
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(12000, 13000, 14000), 1.0f)},
                ConstantPalette{Color4(Pixel(22000, 23000, 24000), 1.0f)});
  auto &a = palettes.a;
  auto &b = palettes.b;
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  PaletteRgbHueLut rgb;
  rgb.bake(cache);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.lab = &cache;
  Modulation::Prepared previous{}, current{};
  const auto source = [](const Modulation::Prepared &p,
                         const Pattern::SphereDirection &) {
    return p.state.phase.phase_turns;
  };
  previous.evaluate = current.evaluate = source;
  previous.state.phase.phase_turns = 0.13f;
  current.state.phase.phase_turns = 0.67f;
  frame.previous_modulation = &previous;
  frame.modulation = &current;
  const FieldSample FIELD{0.5f, 0.5f, math::Y_AXIS, 0.0f};
  for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f})
    for (auto kind : {ColorModulation::Mechanic::HUE_ROTATION,
                      ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                      ColorModulation::Mechanic::OPACITY}) {
      frame.temporal.fraction = fraction;
      frame.mechanic = {kind, 0.4f};
      const auto prepared = Policy::prepare(frame);
      const auto modulation = PaletteLayerProvider::prepare(frame);
      const auto carrier = Policy::apply_fused_modulation(
          FIELD, frame, prepared, modulation, fraction);
      if (fraction > 0.0f && fraction < 1.0f) {
        HS_EXPECT_TRUE(carrier.cached_rgb == &prepared.lab);
        HS_EXPECT_TRUE(carrier.has_cached_lab);
        for (float gain : {0.0f, 0.3f, 1.0f})
          for (float lift : {0.0f, 0.4f}) {
            auto sample = carrier;
            sample.brightness_gain = gain;
            sample.pending_lift = lift;
            const auto expected = finish_color_sample(sample);
            int callbacks = 0;
            HS_EXPECT_FALSE(Color::restore_cached_palette_base_near_gamut_gate(
                sample, [&] { ++callbacks; }));
            const auto actual =
                Color::finish_fused_color(sample, [&] { ++callbacks; });
            HS_EXPECT_EQ(callbacks, 0);
            HS_EXPECT_EQ(actual.color, expected.color);
            HS_EXPECT_EQ(actual.alpha, expected.alpha);
          }
      }
      const auto expected = finish_color_sample(carrier);
      const auto direct = Policy::shade_fused_modulation(FIELD, frame, prepared,
                                                         modulation, fraction);
      HS_EXPECT_EQ(direct.color, expected.color);
      HS_EXPECT_EQ(direct.alpha, expected.alpha);
      const auto pixel =
          Pipeline::shade(math::Y_AXIS, Pipeline::prepare(frame));
      HS_EXPECT_EQ(pixel.color, expected.color);
      HS_EXPECT_EQ(pixel.alpha, expected.alpha);
    }
}

/** @brief Direct cached phase shading matches deferred finishing and its guards. */
inline void test_pullback_cached_phase_finish() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  using Layer = Color::Modulate<CachedSinePhaseProvider>;
  using Pipeline =
      Pullback::Pipeline<PaletteBinding, PaletteSampleStage,
                         Stage::Colorize<Policy>, Stage::Modulate<Layer>>;
  using Core =
      Detail::PipelineCore<true, PaletteBinding, typename Pipeline::NodeList>;
  static_assert(Core::has_direct_final_color());
  TemporalPalettes<> palettes;
  palettes.bake(ConstantPalette{Color4(Pixel(12000, 8000, 4000), 1.0f)},
                ConstantPalette{Color4(Pixel(7000, 15000, 19000), 1.0f)});
  auto &a = palettes.a;
  auto &b = palettes.b;
  PaletteLabLut cache;
  cache.bake(a.view(), b.view());
  PaletteRgbHueLut rgb;
  rgb.bake(cache);
  PaletteFrame frame;
  palettes.bind(frame);
  frame.lab = &cache;
  Color::PhaseModulationState modulation;
  modulation.source = Pattern::SourceKind::DISPLACEMENT;
  modulation.enabled = true;
  frame.phase = &modulation;
  modulation.grid.a = {math::Y_AXIS, 2.1f, 0.37f};
  Modulation::Prepared sine{};
  sine.config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  sine.source = modulation.grid.a;
  sine.evaluate = [](const Modulation::Prepared &prepared,
                     const Pattern::SphereDirection &direction) {
    return Pattern::DirectionalSine::sample(
        direction,
        std::get<Pattern::DirectionalSine::Prepared>(prepared.source));
  };
  modulation.current = modulation.previous = &sine;
  for (auto source : {Pattern::SourceKind::DISPLACEMENT,
                      Pattern::SourceKind::DIRECTIONAL_SINE}) {
    modulation.source = source;
    for (float fraction : {0.0f, 0.125f, 0.5f, 0.875f, 1.0f})
      for (float amount : {-3.7f, -0.2f, 0.0f, 0.65f})
        for (bool ready : {false, true})
          for (auto kind :
               {ColorModulation::Mechanic::HUE_ROTATION,
                ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT}) {
            if (kind == ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT &&
                amount < 0.0f)
              continue;
            modulation.ready = ready;
            modulation.mechanic = {kind, amount};
            frame.temporal.fraction = fraction;
            const auto prepared = Pipeline::prepare_stages(frame);
            const SphereSample INPUT{math::Y_AXIS, 0.7f};
            const auto expected = finish_color_sample(
                Core::template run_stage<0>(INPUT, frame, prepared));
            const auto actual =
                Core::template shade_stage<0>(INPUT, frame, prepared);
            HS_EXPECT_EQ(actual.color, expected.color);
            HS_EXPECT_EQ(actual.alpha, expected.alpha);
            const auto palette = Policy::prepare(frame);
            for (float value : {0.0f, 0.01f, 0.5f, 0.99f, 1.0f}) {
              const FieldSample field{value, 0.5f, math::Y_AXIS, 0.7f};
              ColorSample carrier;
              const bool shared =
                  Policy::apply_shared_base(field, frame, palette, carrier);
              Color4 direct;
              const bool finished = Policy::template shade_cached_phase<Layer>(
                  field, palette, modulation, direct);
              HS_EXPECT_EQ(
                  finished,
                  shared && ready &&
                      (amount == 0.0f ||
                       kind == ColorModulation::Mechanic::HUE_ROTATION));
              if (finished) {
                Layer::apply_phase_in_place(carrier, frame, modulation);
                const auto reference = finish_color_sample(carrier);
                HS_EXPECT_EQ(direct.color, reference.color);
                HS_EXPECT_EQ(direct.alpha, reference.alpha);
              }
            }
          }
  }
}

/** @brief Consecutive analogous cycle bakes meet cached channel bounds. */
inline void test_pullback_palette_lab_qualification() {
  using namespace Pullback;
  using Policy = Color::GeneratedPalette<LabPaletteProvider>;
  std::array<uint8_t, PaletteCycler::generated_arena_bytes() +
                          2 * BakedPalette::required_arena_bytes()>
      buffer{};
  Arena arena(buffer.data(), buffer.size());
  uint32_t hue = 0;
  PaletteCycler cycler;
  cycler.init_generated(
      arena,
      {.next =
           +[](void *context, uint32_t sequence, GenerativePalette &out) {
             GeneratedPaletteBank::next_palette(
                 *static_cast<uint32_t *>(context), sequence,
                 PaletteHarmony::ANALOGOUS, 0.84f, out);
           },
       .context = &hue,
       .dwell_frames = GeneratedPaletteBank::DWELL_FRAMES,
       .fade_frames = GeneratedPaletteBank::FADE_FRAMES,
       .easing = math::ease_in_out_sin,
       .display_interval_frames = 2});
  for (int frame = 0; frame < GeneratedPaletteBank::FADE_FRAMES / 2; ++frame)
    cycler.step();
  BakedPaletteStorage from, to;
  from.clone_from(cycler.palette(), arena);
  cycler.step();
  cycler.step();
  to.clone_from(cycler.palette(), arena);
  PaletteLabLut cache;
  cache.bake(from.view(), to.view());
  HS_EXPECT_TRUE(cache.active);
  PaletteFrame frame;
  frame.palette = &to.view();
  frame.temporal.palette = &from.view();
  frame.lab = &cache;
  uint32_t maximum = 0;
  uint64_t total = 0;
  size_t channels = 0;
  constexpr int VALUES = 2048;
  constexpr int FRACTIONS = 128;
  for (int time = 0; time < FRACTIONS; ++time) {
    frame.temporal.fraction = (static_cast<float>(time) + 0.5f) / FRACTIONS;
    const auto prepared = Policy::prepare(frame);
    for (int index = 0; index < VALUES; ++index) {
      const FieldSample field{(static_cast<float>(index) + 0.5f) / VALUES,
                              0.75f, math::Y_AXIS, 0.0f};
      ColorSample cached;
      HS_EXPECT_TRUE(Policy::apply_shared_base(field, frame, prepared, cached));
      const auto pair = Policy::apply_pair(field, frame, prepared);
      const auto expected = finish_color_sample(blend_color_samples(
          pair.first, pair.second, frame.temporal.fraction));
      const auto actual = finish_color_sample(cached);
      HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-7f);
      for (const uint32_t ERROR :
           {channel_gap(actual.color.r, expected.color.r),
            channel_gap(actual.color.g, expected.color.g),
            channel_gap(actual.color.b, expected.color.b)}) {
        maximum = std::max(maximum, ERROR);
        total += ERROR;
        ++channels;
      }
    }
  }
  const double mean = static_cast<double>(total) / channels;
  std::printf("  palette lab cache max=%u mean=%.3f channel codes\n", maximum,
              mean);
  HS_EXPECT_GT(maximum, 0u);
  HS_EXPECT_LE(static_cast<float>(maximum), Policy::METRICS[0].limit);
  HS_EXPECT_LE(mean, Policy::METRICS[1].limit);
}

} // namespace hs_test::pullback_tests
