/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief KaleidoscopeSmooth composed rendering regressions.
 */
#pragma once

#include "effects/KaleidoscopeSmooth.h"
#include "tests/pixel_test_util.h"
#include "tests/composed_chain_fixture.h"

namespace hs_test {
namespace kaleidoscope_smooth_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

using KaleidoscopeSmoothWhiteBox =
    ComposedWhiteBox<KaleidoscopeSmooth<SMALL_W, SMALL_H>>;

inline void test_kaleidoscope_smooth_identity_and_presets() {
  using WB = KaleidoscopeSmoothWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "kaleidoscope-smooth");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{4});
  static_assert(WB::Params::template HAS<"projection"> &&
                WB::Params::template HAS<"inner_warp"> &&
                WB::Params::template HAS<"source"> &&
                WB::Params::template HAS<"color">);
  static_assert(!WB::Params::template HAS<"lens"> &&
                !WB::Params::template HAS<"surface">);
  static_assert(sizeof(WB::Params) == sizeof(Pullback::ProjectionParams) +
                                          sizeof(Pullback::MirrorParams) +
                                          sizeof(Pullback::GridSourceParams) +
                                          sizeof(Pullback::ColorParams));

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Pattern Freq",
      "Speed",
      "Complexity",
      "Pattern Mix",
      "Drift",
      "Source Angle Speed",
      "Singularity Fade",
      "Projection Spin Speed",
      "Projection Wander",
      "Camera Wander",
      "Planar Warp 2 Speed",
      "Mirror Rotation",
      "Mirror Cell X",
      "Mirror Cell Y",
      "Mirror Offset X",
      "Mirror Offset Y",
      "Palette Chroma",
      "Palette Mapping",
      "Mapping Frequency",
      "Mapping Phase",
      "Phase Oscillation Depth",
      "Phase Oscillation Speed",
      "Opacity at Value 0",
      "Opacity at Value 1",
      "Hue Shift Amount",
      "Hue Noise Scale",
      "Hue Noise Speed",
  };
  HS_EXPECT_EQ(effect.getParameters().size(), std::size(CONTROL_NAMES));
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

inline void test_kaleidoscope_smooth_transition_contract() {
  using FX = KaleidoscopeSmoothWhiteBox::FX;
  check_preset_transition<FX>(
      [](const FX::Params &params) {
        return params.template get<"source">().complexity;
      },
      [](float from, float to, float amount) {
        return from + (to - from) * amount;
      });
}

inline void test_kaleidoscope_smooth_chain_equivalence() {
  for (size_t preset = 0;
       preset < KaleidoscopeSmoothWhiteBox::FX::PRESET_IDS.size(); ++preset)
    composed_chain_tests::verify_export<KaleidoscopeSmoothWhiteBox::FX>(preset);
}

inline int run_kaleidoscope_smooth_tests() {
  ModuleFixture fixture("kaleidoscope_smooth");
  test_kaleidoscope_smooth_identity_and_presets();
  test_kaleidoscope_smooth_transition_contract();
  test_kaleidoscope_smooth_chain_equivalence();
  return fixture.result();
}

} // namespace kaleidoscope_smooth_tests
} // namespace hs_test
