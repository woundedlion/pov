/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Pixel::Feedback — Style binding and plot passthrough (no flush / no Canvas)
// ============================================================================

/**
 * @brief Verifies style() returns a live reference to the bound Style: reads,
 *        mutation, and the const overload all alias the original.
 */
inline void test_feedback_style_binding() {
  constexpr int W = 32, H = 32;
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  Filter::Pixel::Feedback<W, H> fb(style);

  HS_EXPECT_TRUE(&fb.style() == &style);

  const Filter::Pixel::Feedback<W, H> &cfb = fb;
  HS_EXPECT_TRUE(&cfb.style() == &style);

  HS_EXPECT_NEAR(fb.style().fade, 0.9f, 1e-6f);
  HS_EXPECT_EQ(fb.style().downsample, 4);

  // Mutating through the accessor is visible on the original Style.
  fb.style().fade = 0.5f;
  HS_EXPECT_NEAR(style.fade, 0.5f, 1e-6f);
}

/**
 * @brief Verifies Feedback::plot() forwards its input unchanged: one tap with
 *        the original coord, colour, age and alpha, no Canvas touched.
 * @details The feedback effect lives in flush, not plot.
 */
inline void test_feedback_plot_is_passthrough() {
  constexpr int W = 32, H = 32;
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  Filter::Pixel::Feedback<W, H> fb(style);

  int count = 0;
  float kx = -1, ky = -1, kage = -1, ka = -1;
  Pixel kc(0, 0, 0);
  Pixel src(7, 8, 9);
  fb.plot(3.0f, 4.0f, src, 1.5f, 0.75f,
          [&](float x, float y, const Pixel &c, float age, float a) {
            ++count;
            kx = x;
            ky = y;
            kc = c;
            kage = age;
            ka = a;
          });
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_NEAR(kx, 3.0f, 1e-6f);
  HS_EXPECT_NEAR(ky, 4.0f, 1e-6f);
  HS_EXPECT_PIXEL(kc, 7, 8, 9);
  HS_EXPECT_NEAR(kage, 1.5f, 1e-6f);
  HS_EXPECT_NEAR(ka, 0.75f, 1e-6f);
}
