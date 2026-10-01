/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <limits>

#include "core/math/interpolate.h"
#include "effects/KaleidoscopeSmooth.h"
#include "tests/pixel_test_util.h"
#include "tests/test_shader_workbench.h"

namespace hs_test {
namespace kaleidoscope_smooth_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

struct KaleidoscopeSmoothWhiteBox {
  using FX = KaleidoscopeSmooth<SMALL_W, SMALL_H>;
  using Params = FX::Params;

  static constexpr size_t PARAM_CAPACITY = FX::PARAM_CAPACITY;

  static const Params &params(const FX &effect) { return effect.params; }
  static bool transition_active(const FX &effect) {
    return effect.transition.active;
  }
  static bool advance_preset(FX &effect) { return effect.advance_preset(); }
  static void drive_transition(FX &effect, float progress) {
    effect.run_transition(progress);
  }
};

inline void test_kaleidoscope_smooth_identity_and_presets() {
  using WB = KaleidoscopeSmoothWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "kaleidoscope-smooth");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{4});
  HS_EXPECT_EQ(sizeof(WB::Params), 30 * sizeof(float));

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
  using WB = KaleidoscopeSmoothWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{4});
  HS_EXPECT_TRUE(WB::advance_preset(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  uint16_t frames = 0;
  const auto drive = [&](float amount) {
    WB::drive_transition(effect, amount);
    ++frames;
  };
  drive(0.0f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset(0).params.source.complexity, 0.0f);

  drive(0.25f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 interp::linear(FX::preset(0).params.source.complexity,
                                FX::preset(1).params.source.complexity, 0.25f),
                 1e-6f);

  drive(0.5f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 interp::linear(FX::preset(0).params.source.complexity,
                                FX::preset(1).params.source.complexity, 0.5f),
                 1e-6f);

  drive(1.0f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset(1).params.source.complexity, 0.0f);
  HS_EXPECT_TRUE(WB::transition_active(effect));

  for (uint16_t frame = frames;
       frame < Segue::Preset::frames(FX::preset_departure(0)); ++frame)
    drive(0.5f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset(1).params.source.complexity, 0.0f);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

inline void test_kaleidoscope_smooth_shader_workbench_equivalence() {
  constexpr int INNER_MIRROR_GRID = 11;
  constexpr int COMPLEX_GRID = 13;
  constexpr int DOUBLE_MAPPING_GRID = 14;
  shader_workbench_tests::verify_fixed_shader_export<
      KaleidoscopeSmoothWhiteBox::FX>(INNER_MIRROR_GRID, 0);
  shader_workbench_tests::verify_fixed_shader_export<
      KaleidoscopeSmoothWhiteBox::FX>(COMPLEX_GRID, 1);
  shader_workbench_tests::verify_fixed_shader_export<
      KaleidoscopeSmoothWhiteBox::FX>(DOUBLE_MAPPING_GRID, 2);
  shader_workbench_tests::verify_fixed_shader_export<
      KaleidoscopeSmoothWhiteBox::FX>(DOUBLE_MAPPING_GRID, 3);
}

inline int run_kaleidoscope_smooth_tests() {
  ModuleFixture fixture("kaleidoscope_smooth");
  test_kaleidoscope_smooth_identity_and_presets();
  test_kaleidoscope_smooth_transition_contract();
  test_kaleidoscope_smooth_shader_workbench_equivalence();
  return fixture.result();
}

} // namespace kaleidoscope_smooth_tests
} // namespace hs_test
