/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke_stars.h
 * @brief Rendered motion cues and stationary geometry of the Fibonacci star field.
 */

#include "effects/FourStrokeStars.h"
#include "tests/effects/fixture.h"

namespace hs_test::effects_tests {

/** @brief Direction changes shade a stationary star footprint even with the beat held. */
inline void test_four_stroke_stars_fixed_geometry() {
  constexpr int W = 96, H = 20;
  reset_globals();
  FourStrokeStars<W, H> effect;
  effect.init();
  auto render = [&](uint64_t time) {
    effect.set_presentation_interval({time, time});
    effect.draw_frame();
    effect.advance_display();
    std::array<Pixel, W * H> pixels;
    std::copy_n(effect.display_buffer(), pixels.size(), pixels.data());
    return pixels;
  };
  HS_EXPECT_EQ(effect.updateParameter("Strength", 0.0f),
               ParamSetResult::APPLIED);
  const auto unshaded = render(0);
  int occupied = 0;
  for (const auto &pixel : unshaded)
    occupied += !(pixel == FourStroke::GRAY);
  HS_EXPECT_GT(occupied, 20);
  HS_EXPECT_LT(occupied, W * H / 2);

  for (int i = 0; i < 6; ++i)
    HS_EXPECT_TRUE(render(0) == unshaded);
  HS_EXPECT_EQ(effect.updateParameter("Strength", 1.0f),
               ParamSetResult::APPLIED);
  const auto forward = render(0);
  const auto advanced = render(0);
  HS_EXPECT_TRUE(forward != advanced);
  HS_EXPECT_EQ(effect.updateParameter("Reverse Direction", true),
               ParamSetResult::APPLIED);
  const auto reversed = render(0);
  const auto beat = render(62500);
  HS_EXPECT_TRUE(forward != reversed);
  HS_EXPECT_TRUE(reversed != beat);
  for (size_t i = 0; i < unshaded.size(); ++i) {
    if (forward[i] == FourStroke::GRAY) {
      HS_EXPECT_TRUE(advanced[i] == FourStroke::GRAY);
      HS_EXPECT_TRUE(reversed[i] == FourStroke::GRAY);
      HS_EXPECT_TRUE(beat[i] == FourStroke::GRAY);
    }
  }
}

/** @brief A skipped display window advances the same global orientation. */
inline void test_four_stroke_stars_skipped_window() {
  constexpr int W = 96, H = 20;
  const auto render = [](bool skip_middle) {
    reset_globals();
    FourStrokeStars<W, H> effect;
    effect.init();
    for (uint32_t frame = 1; frame <= 3; ++frame) {
      if (skip_middle && frame == 2)
        continue;
      effect.set_scheduled_frame(frame);
      effect.set_presentation_interval({62500u * (frame - 1), 62500u * frame});
      effect.draw_frame();
      effect.advance_display();
    }
    std::array<Pixel, W * H> pixels;
    std::copy_n(effect.display_buffer(), pixels.size(), pixels.begin());
    return pixels;
  };
  HS_EXPECT_TRUE(render(false) == render(true));
}

} // namespace hs_test::effects_tests