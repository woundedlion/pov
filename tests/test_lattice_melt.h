/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <limits>

#include "core/math/interpolate.h"
#include "effects/LatticeMelt.h"
#include "tests/pixel_test_util.h"
#include "tests/test_shader_workbench.h"

namespace hs_test {
namespace lattice_melt_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

struct LatticeMeltWhiteBox {
  using FX = LatticeMelt<SMALL_W, SMALL_H>;
  using Params = FX::Params;

  static constexpr size_t PARAM_CAPACITY = FX::PARAM_CAPACITY;

  static const Params &params(const FX &effect) { return effect.params; }
  static bool transition_active(const FX &effect) {
    return effect.transition.active;
  }
  static bool advance_preset(FX &effect) { return effect.advance_preset(); }
  static void tick_choreography(FX &effect) { effect.step_choreography(); }
  static void saturate_timeline(FX &effect, float &sink) {
    while (Timeline::remaining() > 0)
      effect.timeline.add(
          0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  }
  static void clear_timeline(FX &effect) { effect.timeline.clear(); }
  static void drive_transition(FX &effect, float progress) {
    effect.run_blend(progress);
  }
};

inline void test_lattice_melt_identity_and_presets() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "lattice-melt");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{2});
  HS_EXPECT_EQ(sizeof(WB::Params), 28 * sizeof(float));

  // The runtime rebuilds the hue-rotation LUT on the same predicate the
  // colorizer gates its view on, so both read dead at a zero shift amount.
  Pullback::ColorParams shift;
  shift.hue_shift_amount = 0.0f;
  HS_EXPECT_FALSE(
      (Pullback::hue_rotation_active<Pullback::HueMode::NOISE>(shift)));
  shift.hue_shift_amount = 0.25f;
  HS_EXPECT_TRUE(
      (Pullback::hue_rotation_active<Pullback::HueMode::NOISE>(shift)));
  HS_EXPECT_TRUE(FX::preset_params(0).color.hue_shift_amount != 0.0f);
  HS_EXPECT_TRUE(FX::preset_params(1).color.hue_shift_amount != 0.0f);

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Lattice Cell Scale",      "Lattice Shape",
      "Lattice Softness",        "Lattice Radius",
      "Projection Spin Speed",   "Projection Wander",
      "Camera Wander",           "Central Meridian",
      "Surface Noise Scale",     "Surface Noise Strength",
      "Surface Noise Speed",     "Palette Chroma",
      "Palette Mapping",         "Mapping Frequency",
      "Mapping Phase",           "Phase Oscillation Depth",
      "Phase Oscillation Speed", "Brightness Bottom",
      "Brightness Top",          "Opacity at Value 0",
      "Opacity at Value 1",      "Hue Shift Amount",
      "Hue Noise Scale",         "Hue Noise Speed",
  };
  HS_EXPECT_EQ(effect.getParameters().size(), std::size(CONTROL_NAMES));
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
  // Folded sinusoidal returns fixed weights, so its fade slider is gated off.
  HS_EXPECT_TRUE(effect.getParameters().find("Singularity Fade") == nullptr);
}

inline void test_lattice_melt_transition_contract() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{2});
  HS_EXPECT_TRUE(WB::advance_preset(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::drive_transition(effect, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(0).surface.scale, 0.0f);

  WB::drive_transition(effect, 0.25f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(0).surface.scale *
                     powf(FX::preset_params(1).surface.scale /
                              FX::preset_params(0).surface.scale,
                          0.25f),
                 1e-6f);

  WB::drive_transition(effect, 0.5f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 sqrtf(FX::preset_params(0).surface.scale *
                       FX::preset_params(1).surface.scale),
                 1e-6f);

  WB::drive_transition(effect, 1.0f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(1).surface.scale, 0.0f);
  HS_EXPECT_TRUE(WB::transition_active(effect));

  for (uint16_t frame = 4; frame < FX::TRANSITION_DURATION; ++frame)
    WB::drive_transition(effect, 0.5f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(1).surface.scale, 0.0f);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

inline void test_lattice_melt_full_timeline_retries_transition() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  float sink = 0.0f;
  WB::saturate_timeline(effect, sink);
  const uint32_t dropped_before = Timeline::dropped_events();

  for (uint16_t f = 0; f < FX::PRESET_DWELL_FRAMES; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 1);
  HS_EXPECT_FALSE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});

  // The rejection restarts the dwell, so the retry costs one drop per dwell
  // rather than one per frame.
  for (uint16_t f = 1; f < FX::PRESET_DWELL_FRAMES; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 1);
  WB::tick_choreography(effect);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 2);

  WB::clear_timeline(effect);
  for (uint16_t f = 1; f < FX::PRESET_DWELL_FRAMES; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_FALSE(WB::transition_active(effect));
  WB::tick_choreography(effect);
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
}

inline void test_lattice_melt_overshoot_finishes_on_frame_count() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_TRUE(WB::advance_preset(effect));

  bool saw_overshoot = false;
  for (uint16_t frame = 1; frame < FX::TRANSITION_DURATION; ++frame) {
    const float progress = math::ease_out_elastic(static_cast<float>(frame) /
                                                  FX::TRANSITION_DURATION);
    WB::drive_transition(effect, progress);
    saw_overshoot |= progress > 1.0f;
    HS_EXPECT_TRUE(WB::transition_active(effect));
  }
  HS_EXPECT_TRUE(saw_overshoot);

  WB::drive_transition(effect, math::ease_out_elastic(1.0f));
  HS_EXPECT_FALSE(WB::transition_active(effect));
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(1).surface.scale, 0.0f);
}

/**
 * @brief Pins the dwell restart a manual parameter write owes the choreography.
 */
inline void test_lattice_melt_manual_write_restarts_dwell() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();

  for (uint16_t f = 0; f < FX::PRESET_DWELL_FRAMES; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::drive_transition(effect, 0.5f);
  HS_EXPECT_EQ(effect.updateParameter("Surface Noise Scale", 1.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_FALSE(WB::transition_active(effect));

  effect.setAnimationsPaused(false);
  for (int f = 0; f < FX::PRESET_DWELL_FRAMES - 1; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_FALSE(WB::transition_active(effect));
  WB::tick_choreography(effect);
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
}

inline void test_lattice_melt_shader_workbench_equivalence() {
  shader_workbench_tests::verify_fixed_shader_export<LatticeMeltWhiteBox::FX>(
      7, 0);
  shader_workbench_tests::verify_fixed_shader_export<LatticeMeltWhiteBox::FX>(
      8, 1);
}

inline int run_lattice_melt_tests() {
  ModuleFixture fixture("lattice_melt");
  test_lattice_melt_identity_and_presets();
  test_lattice_melt_transition_contract();
  test_lattice_melt_full_timeline_retries_transition();
  test_lattice_melt_overshoot_finishes_on_frame_count();
  test_lattice_melt_manual_write_restarts_dwell();
  test_lattice_melt_shader_workbench_equivalence();
  return fixture.result();
}

} // namespace lattice_melt_tests
} // namespace hs_test
