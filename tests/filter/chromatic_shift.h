/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Pixel::ChromaticShift::plot — channel-split fan-out
// ============================================================================

/**
 * @brief Verifies ChromaticShift emits the original colour plus three
 *        single-channel copies shifted to x+1/x+2/x+3 (R, G, B respectively),
 *        and that Spread scales those offsets and the segment margin.
 */
inline void test_chromatic_shift_fanout() {
  constexpr int W = 64;
  Filter::Pixel::ChromaticShift<W> cs;

  /**
   * @brief One recorded fan-out tap: position, colour, and alpha.
   */
  struct Tap {
    float x, y;  /**< Emitted pixel coordinate. */
    Pixel c;     /**< Emitted colour. */
    float alpha; /**< Emitted alpha. */
  };
  Tap taps[8]{};
  int count = 0;
  Pixel src(100, 150, 200);
  cs.plot(10.0f, 5.0f, src, 0.0f, 1.0f,
          [&](float x, float y, const Pixel &c, float, float a) {
            if (count < 8)
              taps[count] = {x, y, c, a};
            ++count;
          });

  // Original + 3 channel-shifted copies.
  HS_EXPECT_EQ(count, 4);

  // First tap is the unmodified colour at the original coordinate.
  HS_EXPECT_NEAR(taps[0].x, 10.0f, 1e-5f);
  HS_EXPECT_NEAR(taps[0].y, 5.0f, 1e-5f);
  HS_EXPECT_EQ(taps[0].c.r, src.r);
  HS_EXPECT_EQ(taps[0].c.g, src.g);
  HS_EXPECT_EQ(taps[0].c.b, src.b);

  HS_EXPECT_EQ(taps[0].alpha, 1.0f);
  constexpr int BACKGROUND = 400;
  constexpr int MIN_RETAINED = BACKGROUND * 3 / 4 - 1;
  for (int i = 1; i < 4; ++i) {
    const Pixel lit(BACKGROUND, BACKGROUND, BACKGROUND);
    const Pixel mixed = blend_alpha(taps[i].alpha)(lit, taps[i].c);
    HS_EXPECT_GE(mixed.r, MIN_RETAINED);
    HS_EXPECT_GE(mixed.g, MIN_RETAINED);
    HS_EXPECT_GE(mixed.b, MIN_RETAINED);
  }

  // Red-only copy at x+1.
  HS_EXPECT_NEAR(taps[1].x, 11.0f, 1e-5f);
  HS_EXPECT_EQ(taps[1].c.r, src.r);
  HS_EXPECT_EQ(taps[1].c.g, 0);
  HS_EXPECT_EQ(taps[1].c.b, 0);

  // Green-only copy at x+2.
  HS_EXPECT_NEAR(taps[2].x, 12.0f, 1e-5f);
  HS_EXPECT_EQ(taps[2].c.r, 0);
  HS_EXPECT_EQ(taps[2].c.g, src.g);
  HS_EXPECT_EQ(taps[2].c.b, 0);

  // Blue-only copy at x+3.
  HS_EXPECT_NEAR(taps[3].x, 13.0f, 1e-5f);
  HS_EXPECT_EQ(taps[3].c.r, 0);
  HS_EXPECT_EQ(taps[3].c.g, 0);
  HS_EXPECT_EQ(taps[3].c.b, src.b);

  count = 0;
  cs.plot(63.0f, 5.0f, src, 0.0f, 1.0f,
          [&](float x, float y, const Pixel &c, float, float a) {
            if (count < 8)
              taps[count] = {x, y, c, a};
            ++count;
          });
  HS_EXPECT_EQ(count, 4);
  HS_EXPECT_EQ(taps[0].x, 63.0f);
  HS_EXPECT_EQ(taps[1].x, 0.0f);
  HS_EXPECT_EQ(taps[2].x, 1.0f);
  HS_EXPECT_EQ(taps[3].x, 2.0f);
  HS_EXPECT_EQ(taps[1].c.r, src.r);
  HS_EXPECT_EQ(taps[2].c.g, src.g);
  HS_EXPECT_EQ(taps[3].c.b, src.b);

  // A wider Spread scales every fringe offset and the margin that covers them.
  static_assert(Filter::Pixel::ChromaticShift<W, 3>::segment_margin == 9);
  Filter::Pixel::ChromaticShift<W, 3> wide;
  count = 0;
  wide.plot(10.0f, 5.0f, src, 0.0f, 1.0f,
            [&](float x, float y, const Pixel &c, float, float a) {
              if (count < 8)
                taps[count] = {x, y, c, a};
              ++count;
            });
  HS_EXPECT_EQ(count, 4);
  HS_EXPECT_NEAR(taps[0].x, 10.0f, 1e-5f);
  HS_EXPECT_NEAR(taps[1].x, 13.0f, 1e-5f);
  HS_EXPECT_NEAR(taps[2].x, 16.0f, 1e-5f);
  HS_EXPECT_NEAR(taps[3].x, 19.0f, 1e-5f);
}
