/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_filter.h.

// ============================================================================
// Screen::AntiAlias::plot — pure bilinear weight partition
// ============================================================================

/**
 * @brief Verifies a fractional interior sample splits into the four
 *        neighboring taps, each weighted by the quintic-eased bilinear
 *        coverage, with alphas summing back to the input alpha.
 */
inline void test_antialias_weights_partition() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  const float x = 10.3f, y = 20.6f;
  const float xs = math::quintic_kernel(x - floorf(x));
  const float ys = math::quintic_kernel(y - floorf(y));
  const float in_alpha = 0.8f;
  struct Tap {
    float x, y, a;
  };
  std::vector<Tap> taps;
  aa.plot(x, y, Pixel(1, 2, 3), 0.0f, in_alpha,
          [&](float tx, float ty, const Pixel &, float, float a) {
            taps.push_back({tx, ty, a});
          });
  HS_EXPECT_EQ(taps.size(), size_t{4});

  float sum = 0.0f;
  unsigned seen = 0;
  for (const Tap &t : taps) {
    HS_EXPECT_TRUE(t.x == 10.0f || t.x == 11.0f);
    HS_EXPECT_TRUE(t.y == 20.0f || t.y == 21.0f);
    const float wx = t.x == 10.0f ? 1.0f - xs : xs;
    const float wy = t.y == 20.0f ? 1.0f - ys : ys;
    HS_EXPECT_NEAR(t.a, in_alpha * wx * wy, 1e-6f);
    seen |= 1u << ((t.x == 11.0f ? 1u : 0u) | (t.y == 21.0f ? 2u : 0u));
    sum += t.a;
  }
  HS_EXPECT_EQ(seen, 15u);
  HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
}

/**
 * @brief Verifies integer coordinates collapse AntiAlias to a single full-weight
 *        tap at the input pixel.
 */
inline void test_antialias_integer_coord_single_tap() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  // Integer coords -> fractional parts 0 -> only the (x0,y0) tap carries weight.
  int count = 0;
  float kept_alpha = 0.0f;
  float kept_x = -1.0f, kept_y = -1.0f;
  aa.plot(12.0f, 24.0f, Pixel(0, 0, 0), 0.0f, 0.5f,
          [&](float x, float y, const Pixel &, float, float a) {
            ++count;
            kept_alpha = a;
            kept_x = x;
            kept_y = y;
          });
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_NEAR(kept_alpha, 0.5f, 1e-5f);
  HS_EXPECT_NEAR(kept_x, 12.0f, 1e-5f);
  HS_EXPECT_NEAR(kept_y, 24.0f, 1e-5f);
}

/**
 * @brief Verifies a sub-pixel sample across the theta=0 seam wraps its X taps
 *        onto columns W-1 and 0 rather than collapsing onto one unwrapped column.
 */
inline void test_antialias_seam_wraps_left_column() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  // Integer y = 32 gives zero Y fraction, isolating the two X taps.
  const float in_alpha = 1.0f;
  float sum = 0.0f;
  bool all_in_range = true, saw_w_minus_1 = false, saw_zero = false;
  aa.plot(-0.3f, 32.0f, Pixel(1, 1, 1), 0.0f, in_alpha,
          [&](float x, float y, const Pixel &, float, float a) {
            (void)y;
            sum += a;
            if (x < 0.0f || x >= static_cast<float>(W))
              all_in_range = false;
            int xi = static_cast<int>(x);
            if (xi == W - 1)
              saw_w_minus_1 = true;
            if (xi == 0)
              saw_zero = true;
          });
  HS_EXPECT_TRUE(all_in_range);
  HS_EXPECT_TRUE(saw_w_minus_1);
  HS_EXPECT_TRUE(saw_zero);
  HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
}

/**
 * @brief Verifies a sample at the upper end of the [-W, 2W) x-contract
 *        (x_floor == 2W-1) wraps both taps in range instead of deriving the
 *        right tap from an unwrapped 2W column.
 * @details At x_floor == 2W-1 the second tap's source column is 2W, one past
 *          fast_wrap's x < 2W precondition: it traps the debug assert and writes
 *          column W (out of bounds) in release. The taps must wrap onto W-1 and 0.
 */
inline void test_antialias_far_seam_wraps_both_taps() {
  constexpr int W = 64, H = 64;
  Filter::Screen::AntiAlias<W, H> aa;

  const float in_alpha = 1.0f;
  float sum = 0.0f;
  bool all_in_range = true, saw_w_minus_1 = false, saw_zero = false;
  aa.plot(2.0f * W - 0.3f, 32.0f, Pixel(1, 1, 1), 0.0f, in_alpha,
          [&](float x, float y, const Pixel &, float, float a) {
            (void)y;
            sum += a;
            if (x < 0.0f || x >= static_cast<float>(W))
              all_in_range = false;
            int xi = static_cast<int>(x);
            if (xi == W - 1)
              saw_w_minus_1 = true;
            if (xi == 0)
              saw_zero = true;
          });
  HS_EXPECT_TRUE(all_in_range);
  HS_EXPECT_TRUE(saw_w_minus_1);
  HS_EXPECT_TRUE(saw_zero);
  HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
}

/**
 * @brief Verifies AntiAlias drops fully virtual sub-pole samples (y >= H).
 * @details Boundary samples (H-1 <= y < H) fold the off-edge neighbor's weight
 *          onto the last physical row, preserving their full input alpha.
 */
inline void test_antialias_clips_virtual_subpole_row() {
  constexpr int W = 32, H = 16;
  Filter::Screen::AntiAlias<W, H> aa;

  int subpole_taps = 0;
  aa.plot(10.5f, static_cast<float>(H) + 0.3f, Pixel(40000, 0, 0), 0.0f, 1.0f,
          [&](float, float, const Pixel &, float, float) { ++subpole_taps; });
  HS_EXPECT_EQ(subpole_taps, 0);

  int row_taps = 0;
  int max_row = -1;
  float sum = 0.0f;
  aa.plot(10.5f, static_cast<float>(H - 1) + 0.3f, Pixel(40000, 0, 0), 0.0f,
          1.0f, [&](float, float y, const Pixel &, float, float alpha) {
            ++row_taps;
            sum += alpha;
            int yi = static_cast<int>(y);
            if (yi > max_row)
              max_row = yi;
          });
  HS_EXPECT_TRUE(row_taps > 0);
  HS_EXPECT_EQ(max_row, H - 1);
  HS_EXPECT_NEAR(sum, 1.0f, 1e-4f);
}
