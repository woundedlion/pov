/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Compiled planar source transition contracts.
 */

#include "effects/GridWarp.h"
#include "tests/composed_effect/context.h"

namespace hs_test::composed_effect_tests {

/** @brief The first and last morph frames sample the live endpoints. */
inline void test_grid_warp_source_morph() {
  using namespace Pullback;
  using FX = GridWarp<SMALL_W, SMALL_H>;
  using Selector = FX::Params::Family<"source">;
  using Planar =
      Source::Planar<FX::Binding, Source::PlanarKind::GRID,
                     Source::PlanarKind::SPIRAL, Source::PlanarKind::BANDS>;
  using Grid = Source::Grid<
      SourceProvider<FX::Binding, GridSourceParams, "source_grid">>;
  using Spiral = Source::Spiral<
      SourceProvider<FX::Binding, SpiralSourceParams, "source_spiral">>;
  const size_t GRID = Selector::index(Source::PlanarKind::GRID);
  const size_t SPIRAL = Selector::index(Source::PlanarKind::SPIRAL);
  const size_t BANDS = Selector::index(Source::PlanarKind::BANDS);

  reset_globals();
  pin_frame_clock(0);
  FX effect;
  effect.init();
  const auto start = effect.frame_for_test();
  const auto spiral = FX::test_preset_params(1);
  const auto &wanted = spiral.template get<"source_spiral">();
  HS_EXPECT_EQ(start.resources.template get<"source">().weights.exact, GRID);
  HS_EXPECT_FALSE(effect.start_document_transition(0, 2));
  HS_EXPECT_TRUE(effect.start_document_transition(0, 1));

  pin_frame_clock(1);
  effect.draw_frame();
  const auto first = ComposedFrameWhiteBox::current_frame(effect);
  const auto &first_weights = first.resources.template get<"source">().weights;
  const float first_mix = math::ease_in_out_sin(
      1.0f / Segue::Preset::frames(FX::preset_departure(0)));
  HS_EXPECT_NEAR(first_weights.values[SPIRAL], first_mix, 1e-6f);
  HS_EXPECT_EQ(first.params.template get<"source_spiral">().arms, wanted.arms);

  auto previous = start;
  previous.params.template get<"source_spiral">().arms = 0.0f;
  const auto column = FrameInterval(previous, first).at(0.5f);
  const auto &column_weights =
      column.resources.template get<"source">().weights;
  HS_EXPECT_NEAR(column_weights.values[SPIRAL], first_mix * 0.5f, 1e-6f);
  HS_EXPECT_EQ(column.params.template get<"source_spiral">().arms, wanted.arms);

  PlaneSample plane{};
  plane.coords = {0.12f, -0.38f};
  const float a = Grid::sample(plane, column, Grid::prepare(column));
  const float b = Spiral::sample(plane, column, Spiral::prepare(column));
  HS_EXPECT_NEAR(Planar::sample(plane, column, Planar::prepare(column)),
                 a * (1.0f - first_mix * 0.5f) + b * first_mix * 0.5f, 1e-6f);

  const auto frames = Segue::Preset::frames(FX::preset_departure(0));
  auto penultimate = first;
  for (int index = 2; index <= frames; ++index) {
    pin_frame_clock(index);
    effect.draw_frame();
    if (index + 1 == frames)
      penultimate = ComposedFrameWhiteBox::current_frame(effect);
  }
  auto end = ComposedFrameWhiteBox::current_frame(effect);
  const auto &end_source = end.resources.template get<"source">();
  HS_EXPECT_EQ(end_source.weights.exact, SPIRAL);
  HS_EXPECT_FALSE(end_source.snapped);
  end.params.template get<"source_grid">().pattern_freq = 2.0f;
  const auto last_column = FrameInterval(penultimate, end).at(0.5f);
  const auto &last_mix = last_column.resources.template get<"source">().weights;
  HS_EXPECT_NEAR(
      last_mix.values[GRID],
      0.5f *
          penultimate.resources.template get<"source">().weights.values[GRID],
      1e-6f);
  HS_EXPECT_EQ(last_column.params.template get<"source_grid">().pattern_freq,
               start.params.template get<"source_grid">().pattern_freq);

  HS_EXPECT_TRUE(effect.start_document_transition(0, 1));
  pin_frame_clock(frames + 1);
  effect.draw_frame();
  const auto mixed = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_EQ(mixed.resources.template get<"source">().weights.exact,
               Selector::PlanarWeights::MIXED);
  HS_EXPECT_TRUE(effect.selectPreset(2));
  const auto snapped = effect.frame_for_test();
  const auto snapped_column = FrameInterval(mixed, snapped).at(0.25f);
  HS_EXPECT_EQ(snapped_column.resources.template get<"source">().weights.exact,
               BANDS);
  HS_EXPECT_EQ(snapped_column.params.template get<"source">().source,
               Source::PlanarKind::BANDS);
}

/** @brief Every GridWarp departure is admitted and fits its show slot. */
inline void test_grid_warp_planar_admission() {
  using FX = GridWarp<SMALL_W, SMALL_H>;
  struct Probe : FX {
    using FX::transition_admitted;
  };
  constexpr size_t COUNT = FX::authored_preset_count();
  size_t frames = 0;
  for (size_t index = 0; index < COUNT; ++index) {
    HS_EXPECT_TRUE(Probe::transition_admitted(
        FX::test_preset_params(index),
        FX::test_preset_params((index + 1) % COUNT)));
    frames += FX::PRESET_DWELL_FRAMES +
              Segue::Preset::frames(FX::preset_departure(index));
  }
  HS_EXPECT_TRUE(frames <= 30 * 16);
  auto invalid = FX::test_preset_params(0);
  invalid.template get<"source">().source =
      Pullback::Source::PlanarKind::TWIN_WAVE;
  HS_EXPECT_FALSE(FX::valid_params(invalid));
}

} // namespace hs_test::composed_effect_tests
