/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief KaleidoscopeMelt composed rendering regressions.
 */
#pragma once

#include "effects/KaleidoscopeMelt.h"
#include "tests/pixel_test_util.h"
#include "tests/composed_chain_fixture.h"

namespace hs_test {
namespace kaleidoscope_melt_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;
using KaleidoscopeMeltWhiteBox =
    ComposedWhiteBox<KaleidoscopeMelt<SMALL_W, SMALL_H>>;

inline void test_kaleidoscope_melt_identity_and_controls() {
  using WB = KaleidoscopeMeltWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "kaleidoscope-melt");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{1});
  HS_EXPECT_TRUE(FX::PRESET_IDS[0] == "kaleidoscope-melt");
  static_assert(WB::Params::template HAS<"projection"> &&
                WB::Params::template HAS<"surface"> &&
                WB::Params::template HAS<"source"> &&
                WB::Params::template HAS<"color">);
  static_assert(!WB::Params::template HAS<"lens"> &&
                !WB::Params::template HAS<"inner_warp">);
  static_assert(
      sizeof(WB::Params) ==
      sizeof(Pullback::ProjectionParams) +
          sizeof(Pullback::SurfaceNoiseParams) +
          sizeof(Pullback::LatticeSourceParams) +
          sizeof(Pullback::ColorParams) +
          sizeof(
              Pullback::Modulation::Params<Pattern::SourceKind::NOISE,
                                           Pattern::Evaluation::BAKED_CUBE>));

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Lattice Cell Scale",      "Lattice Shape",
      "Lattice Softness",        "Lattice Radius",
      "Projection Spin Speed",   "Projection Wander",
      "Camera Wander",           "Central Meridian",
      "Surface Noise Scale",     "Surface Noise Strength",
      "Surface Noise Speed",     "Singularity Fade",
      "Palette Chroma",          "Palette Mapping",
      "Mapping Frequency",       "Mapping Phase",
      "Phase Oscillation Depth", "Phase Oscillation Speed",
      "Opacity at Value 0",      "Opacity at Value 1",
      "Modulation Amount",       "Pattern Scale",
      "Pattern Speed",
  };
  HS_EXPECT_EQ(effect.getParameters().size(), std::size(CONTROL_NAMES));
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

inline void test_kaleidoscope_melt_chain_equivalence() {
  composed_chain_tests::verify_export<KaleidoscopeMeltWhiteBox::FX>(0);
}

inline void test_kaleidoscope_melt_cache_controls() {
  using FX = KaleidoscopeMeltWhiteBox::FX;
  reset_effect_globals();
  ScopedTemporalSampling temporal(true);
  pin_frame_clock(0);
  FX effect;
  effect.init();
  const auto INITIAL = KaleidoscopeMeltWhiteBox::params(effect);
  const auto draw = [&](int frame) {
    pin_frame_clock(frame);
    effect.draw_frame();
    const auto &CURRENT = ComposedFrameWhiteBox::current_frame(effect);
    const bool CACHED =
        CURRENT.palette_lab_lut && CURRENT.palette_lab_lut->active;
    effect.advance_display();
    return CACHED;
  };
  HS_EXPECT_TRUE(draw(1));
  for (int control = 0; control < 4; ++control) {
    auto changed = INITIAL;
    auto &color = changed.template get<"color">();
    if (control == 0)
      color.opacity_high = 1.0f;
    else if (control == 1)
      color.palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    else if (control == 2)
      changed.template get<"modulation">().modulation.mechanic.amount = 100.0f;
    else
      changed.template get<"modulation">().modulation.noise.scale = 1.0f;
    ComposedFrameWhiteBox::set_params(effect, changed);
    HS_EXPECT_FALSE(draw(2 + control * 3));
    HS_EXPECT_FALSE(draw(3 + control * 3));
  }
  ComposedFrameWhiteBox::set_params(effect, INITIAL);
  HS_EXPECT_FALSE(draw(14));
  HS_EXPECT_TRUE(draw(15));
}

inline int run_kaleidoscope_melt_tests() {
  ModuleFixture fixture("kaleidoscope_melt");
  HS_RUN(test_kaleidoscope_melt_identity_and_controls);
  HS_RUN(test_kaleidoscope_melt_chain_equivalence);
  HS_RUN(test_kaleidoscope_melt_cache_controls);
  return fixture.result();
}

} // namespace kaleidoscope_melt_tests
} // namespace hs_test
