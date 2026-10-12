/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Controlled composed-effect frame fixtures.
 */
#pragma once

#include "core/render/pullback/composed_effect.h"
#include "tests/test_fixture.h"

namespace hs_test {

/**
 * @brief Composed frame preparation with controlled parameters and clocks.
 * @details Skips timeline, choreography and palette cycling; callers supply
 *          fixed parameters and drive runtime clocks explicitly.
 */
struct ComposedFrameWhiteBox {
  template <typename FX>
  static void set_params(FX &effect, const typename FX::Params &params) {
    effect.adopt_params(params);
  }

  template <typename FX> static void advance(FX &effect) {
    effect.advance_runtime();
    effect.update_spatial_frames();
    effect.update_palette_chroma();
  }

  template <typename FX> static typename FX::FrameState frame(FX &effect) {
    return effect.prepare_frame();
  }
};

template <typename Effect> struct ComposedWhiteBox {
  using FX = Effect;
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
    effect.run_transition(progress);
  }
};

template <typename FX, typename Value, typename Interpolate>
inline void check_preset_transition(Value value, Interpolate interpolate) {
  using WB = ComposedWhiteBox<FX>;
  reset_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), FX::PRESET_IDS.size());
  HS_EXPECT_TRUE(WB::advance_preset(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  const float from = value(FX::preset(0).params);
  const float to = value(FX::preset(1).params);
  uint16_t frames = 0;
  const auto drive = [&](float amount) {
    WB::drive_transition(effect, amount);
    ++frames;
  };
  for (float progress : {0.0f, 0.25f, 0.5f, 1.0f}) {
    drive(progress);
    HS_EXPECT_NEAR(value(WB::params(effect)), interpolate(from, to, progress),
                   1e-6f);
    HS_EXPECT_TRUE(WB::transition_active(effect));
  }
  for (uint16_t frame = frames;
       frame + 1 < Segue::Preset::frames(FX::preset_departure(0)); ++frame)
    drive(0.5f);
  HS_EXPECT_TRUE(WB::transition_active(effect));
  drive(0.5f);
  HS_EXPECT_NEAR(value(WB::params(effect)), to, 0.0f);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

} // namespace hs_test
