/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Filigree edge-store and Fade departure tests. */

#include "tests/composed_effect/context.h"

namespace hs_test::composed_effect_tests {

/** @brief Resets globals before a probe constructs its effect. */
struct FreshGlobals {
  /** @brief Resets globals and pins the frame clock to zero. */
  FreshGlobals() {
    reset_globals();
    pin_frame_clock(0);
  }
};

/** @brief Filigree on a small canvas with frame and pixel probes. */
struct FiligreeProbe : FreshGlobals {
  using FX = Filigree<SMALL_W, SMALL_H>; ///< Probed effect.
  FX effect;                             ///< Effect under test.
  int clock = 0;                         ///< Pinned frame clock.

  /** @brief Initializes the effect on fresh globals. */
  FiligreeProbe() { effect.init(); }

  /** @brief Draws and displays one frame.
   * @return Sum of every channel of the displayed frame. */
  uint64_t draw() {
    pin_frame_clock(++clock);
    effect.draw_frame();
    effect.advance_display();
    uint64_t sum = 0;
    for (int i = 0; i < SMALL_W * SMALL_H; ++i) {
      const Pixel pixel = effect.display_buffer()[i];
      sum += pixel.r + pixel.g + pixel.b;
    }
    return sum;
  }

  /** @brief Whether the runtime's edge store holds a solid.
   * @param solid Registry solid ordinal.
   * @return True when the solid's edges are built. */
  bool stored(uint8_t solid) {
    const auto *store = effect.frame_for_test().sphere_tiling;
    return store->built_solid() == solid && !store->edges().edges.empty();
  }

  /** @brief Tiling solid of an authored preset.
   * @param index Preset number.
   * @return Registry solid ordinal. */
  static uint8_t solid(size_t index) {
    return FX::preset(index).params.template get<"source">().solid;
  }
};

/// Both Filigree departures fade through a black frame that swaps the stored solid.
inline void test_filigree_fade_swaps_solid_at_zero_opacity() {
  using FX = FiligreeProbe::FX;
  auto probe = std::make_unique<FiligreeProbe>();
  for (size_t from : {0u, 1u}) {
    HS_CONTEXT("departing preset", static_cast<int>(from));
    const auto fade = std::get<Segue::Preset::Fade>(FX::preset_departure(from));
    HS_EXPECT_EQ(fade.frames % 2, 0);
    const uint8_t outgoing = FiligreeProbe::solid(from);
    const uint8_t incoming = FiligreeProbe::solid(1 - from);
    HS_EXPECT_NE(outgoing, incoming);
    HS_EXPECT_TRUE(probe->effect.selectPreset(from));
    const uint64_t lit = probe->draw();
    HS_EXPECT_GT(lit, uint64_t{0});
    HS_EXPECT_TRUE(probe->stored(outgoing));
    HS_EXPECT_TRUE(probe->effect.start_document_transition(from, 1 - from));
    const unsigned middle = fade.frames / 2;
    for (unsigned tick = 1; tick <= fade.frames; ++tick) {
      HS_CONTEXT("tick", static_cast<int>(tick));
      const uint64_t brightness = probe->draw();
      if (tick == middle) {
        HS_EXPECT_EQ(brightness, uint64_t{0});
        HS_EXPECT_TRUE(probe->stored(incoming));
        HS_EXPECT_FALSE(probe->stored(outgoing));
      } else {
        HS_EXPECT_GT(brightness, uint64_t{0});
        HS_EXPECT_TRUE(probe->stored(tick < middle ? outgoing : incoming));
      }
      if (tick == middle - 1)
        HS_EXPECT_LT(brightness, lit / 4);
    }
  }
}

/// Filigree rebuilds solids in place, so preset changes keep the persistent offset.
inline void test_filigree_solid_changes_keep_persistent_offset() {
  auto probe = std::make_unique<FiligreeProbe>();
  probe->effect.set_clip(0, 1, 0, 1);
  probe->draw();
  const auto cycle = [&] {
    for (size_t preset : {1u, 0u}) {
      HS_EXPECT_TRUE(probe->effect.selectPreset(preset));
      probe->draw();
      HS_EXPECT_TRUE(probe->stored(FiligreeProbe::solid(preset)));
    }
  };
  cycle();
  const size_t settled = persistent_arena.get_offset();
  for (int round = 0; round < 4; ++round) {
    cycle();
    HS_EXPECT_EQ(persistent_arena.get_offset(), settled);
  }
}

/// Filigree draws one band per preset, so its stored index only changes at a
/// fade's dark frame.
inline void test_filigree_presets_keep_their_band_side() {
  using FX = FiligreeProbe::FX;
  using Store = Pullback::Source::SphereTilingEdgeStore;
  for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index) {
    HS_CONTEXT("preset", static_cast<int>(index));
    const size_t next = (index + 1) % FX::PRESET_IDS.size();
    HS_EXPECT_TRUE(
        std::holds_alternative<Segue::Preset::Fade>(
            FX::preset_departure(index)) ||
        Store::thin_band(FX::preset(index).params.template get<"source">()) ==
            Store::thin_band(FX::preset(next).params.template get<"source">()));
  }
}

/// A GUI band crossing rebuilds the stored index in place, and the drawn
/// edges match a fresh build at the new band.
inline void test_filigree_band_crossing_rebuilds_the_index() {
  using namespace Pullback::Source;
  auto probe = std::make_unique<FiligreeProbe>();
  probe->effect.set_clip(0, 1, 0, 1);
  probe->draw();
  const size_t settled = persistent_arena.get_offset();
  std::vector<uint8_t> block(1 << 20);
  std::vector<uint8_t> scratch(1 << 20);
  auto scratch_a = std::make_unique<Arena>(scratch.data(), scratch.size() / 2);
  auto scratch_b = std::make_unique<Arena>(scratch.data() + scratch.size() / 2,
                                           scratch.size() / 2);
  const auto source = [&] {
    return ComposedFrameWhiteBox::params(probe->effect)
        .template get<"source">();
  };
  const auto expect_fresh = [&] {
    const auto params = source();
    const SphereTilingEdgeBand band(params.line_thickness,
                                    params.line_softness);
    const auto stored = probe->effect.frame_for_test().sphere_tiling->edges();
    const auto fresh =
        SphereTilingGeometry::build(params.solid, block, *scratch_a, *scratch_b)
            .edge_view();
    HS_EXPECT_EQ(stored.offsets.empty(),
                 band.outer <= SphereTilingNarrowIndex::OUTER);
    for (int i = 0; i < 2048; ++i) {
      const auto direction = math::fib_spiral(2048, 0.29f, i);
      const auto a = stored.sample(direction, band);
      const auto b = fresh.sample(direction, band);
      HS_EXPECT_EQ(a.coverage, b.coverage);
      HS_EXPECT_EQ(a.value, b.value);
    }
  };
  for (int round = 0; round < 3; ++round) {
    for (float thickness : {6.0f, 0.0f}) {
      HS_CONTEXT("thickness", static_cast<int>(thickness));
      HS_EXPECT_EQ(
          probe->effect.updateParameter("Edge Thickness (deg)", thickness),
          ParamSetResult::APPLIED);
      probe->draw();
      probe->draw();
      expect_fresh();
      HS_EXPECT_EQ(persistent_arena.get_offset(), settled);
    }
  }
}

/// Filigree admits GUI tiling writes only for edge solids that fit its store.
inline void test_filigree_admits_solids_within_the_edge_store() {
  using Store = Pullback::Source::SphereTilingEdgeStore;
  auto probe = std::make_unique<FiligreeProbe>();
  for (size_t index = 0; index < FiligreeProbe::FX::PRESET_IDS.size(); ++index)
    HS_EXPECT_TRUE(Store::admits(FiligreeProbe::solid(index)));
  probe->effect.set_clip(0, 1, 0, 1);
  uint8_t last = FiligreeProbe::solid(0);
  for (unsigned solid = 0; solid < Solids::NUM_ENTRIES; ++solid) {
    HS_CONTEXT("solid", static_cast<int>(solid));
    const auto result =
        probe->effect.updateParameter("solid", static_cast<float>(solid));
    HS_EXPECT_EQ(result, Store::admits(solid) ? ParamSetResult::APPLIED
                                              : ParamSetResult::INADMISSIBLE);
    if (Store::admits(solid))
      last = static_cast<uint8_t>(solid);
  }
  HS_EXPECT_EQ(ComposedFrameWhiteBox::params(probe->effect)
                   .template get<"source">()
                   .solid,
               last);
  HS_EXPECT_EQ(probe->effect.updateParameter("output", 1.0f),
               ParamSetResult::INADMISSIBLE);
  probe->draw();
  HS_EXPECT_TRUE(probe->stored(last));
}

} // namespace hs_test::composed_effect_tests
