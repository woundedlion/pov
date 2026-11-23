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

/** @brief Direction and presentation time change shading on a fixed footprint. */
inline void test_four_stroke_stars_fixed_geometry() {
  constexpr int W = 96, H = 20;
  using Frame = std::array<Pixel, W * H>;
  const auto render = [](float strength, bool reverse, uint64_t final_time,
                         int width = 1, bool black_background = false) {
    reset_globals();
    FourStrokeStars<W, H> effect;
    effect.init();
    HS_EXPECT_EQ(effect.updateParameter("Strength", strength),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(effect.updateParameter("Edge Width", width),
                 ParamSetResult::APPLIED);
    if (reverse)
      HS_EXPECT_EQ(effect.updateParameter("Reverse Direction", true),
                   ParamSetResult::APPLIED);
    if (black_background)
      HS_EXPECT_EQ(effect.updateParameter("Black Background", true),
                   ParamSetResult::APPLIED);
    std::array<Frame, 7> frames;
    for (size_t i = 0; i < frames.size(); ++i) {
      const uint64_t time = i + 1 == frames.size() ? final_time : 0;
      effect.set_presentation_interval({time, time});
      effect.draw_frame();
      effect.advance_display();
      std::copy_n(effect.display_buffer(), frames[i].size(), frames[i].begin());
    }
    return frames;
  };

  const auto stationary = render(0.0f, false, 0);
  int occupied = 0;
  for (const auto &pixel : stationary.front())
    occupied += !(pixel == FourStroke::GRAY);
  HS_EXPECT_GT(occupied, 20);
  HS_EXPECT_LT(occupied, W * H / 2);
  for (const auto &frame : stationary)
    HS_EXPECT_TRUE(frame == stationary.front());

  const auto forward = render(1.0f, false, 0);
  const auto reverse = render(1.0f, true, 0);
  const auto beat = render(1.0f, true, 62500);
  HS_EXPECT_TRUE(forward.front() != forward.back());
  HS_EXPECT_TRUE(forward.back() != reverse.back());
  HS_EXPECT_TRUE(reverse.back() != beat.back());
  const auto wide = render(1.0f, true, 62500, 2);
  HS_EXPECT_TRUE(wide.back() != beat.back());
  const auto black = render(1.0f, true, 62500, 1, true);
  HS_EXPECT_TRUE(black.front() != black.back());
  const Pixel first_fill = FourStroke::colors(0).fill;
  const Pixel second_fill = FourStroke::colors(62500).fill;
  const Pixel black_background(0, 0, 0);
  int center_pixels = 0;
  for (size_t i = 0; i < forward.back().size(); ++i) {
    if (forward.back()[i] == first_fill && reverse.back()[i] == first_fill &&
        beat.back()[i] == second_fill && wide.back()[i] == second_fill)
      ++center_pixels;
  }
  HS_EXPECT_GT(center_pixels, 10);
  for (size_t i = 0; i < forward.front().size(); ++i) {
    if (beat.back()[i] == second_fill)
      HS_EXPECT_TRUE(black.back()[i] == second_fill);
    if (beat.back()[i] == FourStroke::GRAY)
      HS_EXPECT_TRUE(black.back()[i] == black_background);
    if (forward.front()[i] == FourStroke::GRAY) {
      HS_EXPECT_TRUE(forward[1][i] == FourStroke::GRAY);
      HS_EXPECT_TRUE(reverse.back()[i] == FourStroke::GRAY);
      HS_EXPECT_TRUE(beat.back()[i] == FourStroke::GRAY);
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