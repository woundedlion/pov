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
  static_assert(sizeof(WB::Params) ==
                sizeof(Pullback::ProjectionParams) +
                    sizeof(Pullback::SurfaceNoiseParams) +
                    sizeof(Pullback::LatticeSourceParams) +
                    sizeof(Pullback::ColorParams));

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
      "Hue Shift Amount",        "Hue Noise Scale",
      "Hue Noise Speed",
  };
  HS_EXPECT_EQ(effect.getParameters().size(), std::size(CONTROL_NAMES));
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

inline void test_kaleidoscope_melt_chain_equivalence() {
  composed_chain_tests::verify_export<KaleidoscopeMeltWhiteBox::FX>(0);
}

inline int run_kaleidoscope_melt_tests() {
  ModuleFixture fixture("kaleidoscope_melt");
  test_kaleidoscope_melt_identity_and_controls();
  test_kaleidoscope_melt_chain_equivalence();
  return fixture.result();
}

} // namespace kaleidoscope_melt_tests
} // namespace hs_test
