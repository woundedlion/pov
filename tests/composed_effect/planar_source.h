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

/** @brief Rebinds the grid sample to the planar selector. */
template <typename Pipeline> struct PlanarMorphPipeline;

/// Substitutes the planar selector while preserving the original crossing policies.
template <typename B, typename Rotate, typename Displace, typename Project,
          typename Sample, typename... Tail>
struct PlanarMorphPipeline<
    Pullback::Pipeline<B, Rotate, Displace, Project, Sample, Tail...>> {
  using Type = Pullback::Pipeline<
      B, Rotate, Displace, Project,
      Pullback::Stage::Sample<
          Pullback::Source::Planar<B, Pullback::Source::PlanarKind::GRID,
                                   Pullback::Source::PlanarKind::SPIRAL,
                                   Pullback::Source::PlanarKind::BANDS>,
          typename Sample::WeightPolicy, typename Sample::CoveragePolicy>,
      Tail...>;
};

/** @brief GridWarp pipeline with the planar source selector. */
struct PlanarMorphSpec : GridWarpSpec {
  /** @brief Pipeline for the planar selector probe.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline =
      typename PlanarMorphPipeline<GridWarpSpec::Pipeline<B>>::Type;
};

/** @brief Test choreography with a grid, a spiral and bands. */
class PlanarMorphProbe
    : public Pullback::ComposedEffect<SMALL_W, SMALL_H, PlanarMorphProbe,
                                      PlanarMorphSpec> {
  using Base = Pullback::ComposedEffect<SMALL_W, SMALL_H, PlanarMorphProbe,
                                        PlanarMorphSpec>;
  using Product = GridWarp<SMALL_W, SMALL_H>;

public:
  using Params =
      Pullback::ParamsFor<PlanarMorphSpec>; ///< Probe parameter block.
  using Base::transition_admitted;

  static constexpr auto EFFECT_ID = Product::EFFECT_ID; ///< Probe identity.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 84;   ///< Hold between edges.
  static constexpr std::array<std::string_view, 3> PRESET_IDS{
      "grid", "spiral", "bands"}; ///< Probe preset IDs.

  /** @brief Supplies the grid endpoint.
   * @return Initial grid parameters. */
  static constexpr Params initial_params() {
    const auto grid = Product::initial_params();
    Params params;
    params.template get<"projection">() = grid.template get<"projection">();
    params.template get<"displace">() = grid.template get<"displace">();
    params.template get<"source_grid">() = grid.template get<"source">();
    params.template get<"color">() = grid.template get<"color">();
    params.template get<"modulate">() = grid.template get<"modulate">();
    params.template get<"modulate1">() = grid.template get<"modulate1">();
    return params;
  }

  /** @brief Returns a probe endpoint and its timed departure.
   * @param index Probe preset index.
   * @return The endpoint row. */
  static constexpr PresetEntry<Params> preset(size_t index) {
    using Pullback::Source::PlanarKind;
    Params params = initial_params();
    if (index == 1) {
      params.template get<"source">().source = PlanarKind::SPIRAL;
      params.template get<"source_spiral">().arms = 4.0f;
      params.template get<"source_spiral">().pattern_freq = 3.4f;
    } else if (index == 2) {
      params.template get<"source">().source = PlanarKind::BANDS;
      params.template get<"source_bands">().angle = 0.7f;
      params.template get<"source_bands">().pattern_freq = 4.2f;
    }
    return {params, Segue::Preset::Lerp{52, math::ease_in_out_sin}};
  }

  /** @brief Binds the probe's choreography.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Base::choreography_hooks()
        .template extend<
            PlanarMorphProbe,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &PlanarMorphProbe::initial_params>,
            Control::HookBinding<Control::ChoreographyHook::PRESET_ROW,
                                 &PlanarMorphProbe::preset>>();
  }
};

/** @brief The first and last morph frames sample the live endpoints. */
inline void test_cross_kind_source_morph() {
  using namespace Pullback;
  using FX = PlanarMorphProbe;
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
  effect.advance_display();
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
  const float column_mix = math::ease_in_out_sin(
      0.5f / Segue::Preset::frames(FX::preset_departure(0)));
  HS_EXPECT_NEAR(Planar::sample(plane, column, Planar::prepare(column)),
                 hs::lerp(a, b, column_mix), 1e-6f);

  const auto frames = Segue::Preset::frames(FX::preset_departure(0));
  auto penultimate = first;
  for (int index = 2; index <= frames; ++index) {
    pin_frame_clock(index);
    effect.draw_frame();
    effect.advance_display();
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
  effect.advance_display();
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

/** @brief Grid fields with different frequencies blend at eased column time. */
inline void test_same_kind_source_morph() {
  using namespace Pullback;
  using FX = PlanarMorphProbe;
  using Grid = Source::Grid<
      SourceProvider<FX::Binding, GridSourceParams, "source_grid">>;
  using Planar =
      Source::Planar<FX::Binding, Source::PlanarKind::GRID,
                     Source::PlanarKind::SPIRAL, Source::PlanarKind::BANDS>;
  reset_globals();
  FX effect;
  effect.init();
  auto from = effect.frame_for_test();
  auto to = from;
  to.params.template get<"source_spiral">().arms = 7.0f;
  HS_EXPECT_FALSE(Source::different_endpoint_fields(from.params, to.params));
  to.params.template get<"source_grid">().pattern_freq = 9.0f;
  HS_EXPECT_TRUE(Source::different_endpoint_fields(from.params, to.params));
  const auto a_params = from.params;
  const auto b_params = to.params;
  from.source_morph = {&a_params, &b_params, math::ease_in_out_sin, 0.25f,
                       math::ease_in_out_sin(0.25f)};
  to.source_morph = {&a_params, &b_params, math::ease_in_out_sin, 0.75f,
                     math::ease_in_out_sin(0.75f)};
  const auto column = FrameInterval(from, to).at(0.25f);
  const float mix = math::ease_in_out_sin(0.375f);
  HS_EXPECT_NEAR(column.source_morph.mix, mix, 1e-6f);
  PlaneSample plane{};
  plane.coords = {0.23f, -0.39f};
  const float a = Grid::sample(plane, from, Grid::prepare(from));
  const float b = Grid::sample(plane, to, Grid::prepare(to));
  HS_EXPECT_NEAR(Planar::sample(plane, column, Planar::prepare(column)),
                 hs::lerp(a, b, mix), 1e-6f);
}

/** @brief Enables endpoint sampling on the fixed hex-cells policy. */
struct EndpointHexSpec : RiverStonesSpec {
  /// Samples the hex field at both authored endpoints.
  static constexpr bool SOURCE_MORPH = true;
};

/** @brief Fixed hex cells blend the two preset scales and softness samples. */
inline void test_fixed_hex_source_morph() {
  using namespace Pullback;
  using Traits = ComposedTraits<SMALL_W, EndpointHexSpec>;
  constexpr auto KIND = Pattern::TessellationKind::HEXAGONAL;
  constexpr auto FILL = Pattern::TessellationFill::CELLS;
  using Policy = SourcePolicyFor<FixedTessellationParams<KIND, FILL>,
                                 Traits::Binding>::Type;
  Traits::Params from;
  Traits::Params to = from;
  auto &a = from.template get<"source">();
  auto &b = to.template get<"source">();
  a.cell_scale = 2.5352101f;
  a.cell_softness = 0.10124387f;
  b.cell_scale = 4.774648f;
  b.cell_softness = 0.32801026f;
  HS_EXPECT_TRUE(Source::different_endpoint_fields(from, to));
  HS_EXPECT_FALSE(Source::different_endpoint_fields(from, from));

  Traits::FrameState frame{};
  frame.params = Pullback::interpolate(from, to, 0.25f);
  frame.source_morph = {&from, &to, math::ease_linear, 0.25f, 0.25f};
  for (const auto point :
       {math::Complex{0.17f, 0.23f}, math::Complex{-0.39f, 0.11f},
        math::Complex{0.62f, -0.48f}}) {
    PlaneSample plane{};
    plane.coords = point;
    const Pattern::PlanePosition position{point};
    const float outgoing = Pattern::tessellation<KIND, FILL>(
        position, static_cast<const Pattern::TessellationCellParams &>(a));
    const float incoming = Pattern::tessellation<KIND, FILL>(
        position, static_cast<const Pattern::TessellationCellParams &>(b));
    HS_EXPECT_NEAR(Policy::sample(plane, frame),
                   hs::lerp(outgoing, incoming, 0.25f), 1e-6f);
  }
}

/** @brief Exercises the Teensy approximation on the native test target. */
struct TableHexSpec : RiverStonesSpec {
  /// Enables the pre-baked endpoint hex fields.
  static constexpr bool HEX_MORPH_TABLE = true;
};

/// Checks that table-backed endpoint fields agree with the analytic source.
inline void test_hex_morph_table() {
  using namespace Pullback;
  using Traits = ComposedTraits<SMALL_W, TableHexSpec>;
  using Policy =
      SourcePolicyFor<typename Traits::Params::template Family<"source">,
                      Traits::Binding>::Type;
  using FX = RiverStones<SMALL_W, SMALL_H>;
  static_assert(Policy::ORACLE == ApproximationOracleId::HEX_CELL_ENDPOINTS);
  const auto a = FX::preset(0).params, b = FX::preset(1).params;
  Source::HexMorphTables tables;
  tables.prewarm(a.template get<"source">(), b.template get<"source">());

  const auto check = [&](const auto &from, const auto &to, bool cached) {
    const auto &f = from.template get<"source">();
    const auto &t = to.template get<"source">();
    Traits::FrameState frame{};
    frame.hex_tables = tables.select(f, t);
    HS_EXPECT_EQ(frame.hex_tables.from != nullptr, cached);
    HS_EXPECT_EQ(frame.hex_tables.to != nullptr, cached);
    float worst = 0.0f;
    bool exact_endpoints = true;
    for (const float mix : {0.0f, 0.125f, 0.5f, 0.875f, 1.0f}) {
      frame.source_morph = {&from, &to, math::ease_linear, mix, mix};
      frame.params = mix == 1.0f ? to : from;
      for (int y = -24; y <= 24; ++y) {
        for (int x = -32; x <= 32; ++x) {
          const Pattern::PlanePosition coord{{x * 0.067f, y * 0.083f}};
          PlaneSample plane{};
          plane.coords = coord.value;
          const auto analytic = [&](const auto &source) {
            return Pattern::hexagonal_cells(
                coord,
                static_cast<const Pattern::TessellationCellParams &>(source));
          };
          const float ca = analytic(f), cb = analytic(t);
          const float actual = Policy::sample(plane, frame);
          if (mix == 0.0f || mix == 1.0f)
            exact_endpoints &= actual == (mix == 0.0f ? ca : cb);
          else
            worst = std::max(worst, fabsf(actual - hs::lerp(ca, cb, mix)));
        }
      }
    }
    HS_EXPECT_TRUE(exact_endpoints);
    if (cached) {
      HS_EXPECT_TRUE(worst > 0.0f);
      HS_EXPECT_TRUE(worst <= Policy::METRICS.front().limit);
    } else {
      HS_EXPECT_EQ(worst, 0.0f);
    }
  };

  check(a, b, true);
  check(b, a, true);
  auto edit = a;
  edit.template get<"source">().cell_scale += 0.75f;
  check(edit, b, true);
  for (auto field : {&Pattern::TessellationCellParams::cell_radius,
                     &Pattern::TessellationCellParams::cell_softness,
                     &Pattern::TessellationCellParams::cell_shape}) {
    edit = a;
    edit.template get<"source">().*field += 0.01f;
    check(edit, b, false);
    check(b, edit, false);
  }
}

/** @brief The planar selector admits its supported departures. */
inline void test_planar_transition_admission() {
  using FX = PlanarMorphProbe;
  constexpr size_t COUNT = FX::authored_preset_count();
  for (size_t index = 0; index < COUNT; ++index)
    HS_EXPECT_TRUE(
        FX::transition_admitted(FX::test_preset_params(index),
                                FX::test_preset_params((index + 1) % COUNT)));
  auto invalid = FX::test_preset_params(0);
  invalid.template get<"source">().source =
      Pullback::Source::PlanarKind::TWIN_WAVE;
  HS_EXPECT_FALSE(FX::valid_params(invalid));
}

} // namespace hs_test::composed_effect_tests
