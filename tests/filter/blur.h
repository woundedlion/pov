/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Screen::Blur::plot — kernel passthrough
// ============================================================================

/**
 * @brief Verifies factor=0 makes Blur an identity: a single full-weight tap at
 *        the input pixel.
 */
inline void test_blur_factor_zero_is_identity() {
  constexpr int W = 32, H = 32;
  // center weight c = 1.0, edge/corner weights 0 → single tap.
  Filter::Screen::Blur<W, H> blur(0.0f);

  int count = 0;
  float kept_alpha = 0.0f;
  float kept_x = -1.0f, kept_y = -1.0f;
  blur.plot(8.0f, 16.0f, Pixel(4, 5, 6), 2.0f, 1.0f,
            [&](float x, float y, const Pixel &, float, float a) {
              ++count;
              kept_alpha = a;
              kept_x = x;
              kept_y = y;
            });
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_NEAR(kept_alpha, 1.0f, 1e-5f);
  HS_EXPECT_NEAR(kept_x, 8.0f, 1e-5f);
  HS_EXPECT_NEAR(kept_y, 16.0f, 1e-5f);
}

/** @brief Blur folds sub-pole taps and drops virtual rows at all strengths. */
inline void test_blur_folds_boundary_and_drops_virtual_rows() {
  constexpr int W = 32, H = 32;
  Filter::Screen::Blur<W, H> blur(0.0f);
  for (float factor : {0.0f, 1e-6f, 1.0f}) {
    blur.update(factor);
    float total = 0.0f;
    blur.plot(8.0f, H - 0.25f, Pixel(4, 5, 6), 2.0f, 0.75f,
              [&](float, float y, const Pixel &, float, float alpha) {
                HS_EXPECT_EQ(y, static_cast<float>(H - 1));
                total += alpha;
              });
    HS_EXPECT_NEAR(total, 0.75f, 1e-5f);
  }
  int outside_taps = 0;
  for (float factor : {0.0f, 1e-6f, 1.0f}) {
    blur.update(factor);
    for (float y : {-3.0f, -1.2f, H + 0.2f})
      blur.plot(
          8.0f, y, Pixel(4, 5, 6), 2.0f, 1.0f,
          [&](float, float, const Pixel &, float, float) { ++outside_taps; });
  }
  HS_EXPECT_EQ(outside_taps, 0);
}

/**
 * @brief Verifies factor=1 gives the full 3x3 Gaussian: all 9 taps land
 *        in-bounds for an interior pixel and their alphas sum to the input alpha.
 */
inline void test_blur_full_kernel_sums_to_alpha() {
  constexpr int W = 32, H = 32;
  // Gaussian weights: center 0.25, edge 0.125, corner 0.0625 (sum to 1).
  Filter::Screen::Blur<W, H> blur(1.0f);

  float sum = 0.0f;
  int count = 0;
  const float in_alpha = 0.6f;
  blur.plot(15.0f, 16.0f, Pixel(1, 1, 1), 0.0f, in_alpha,
            [&](float, float, const Pixel &, float, float a) {
              sum += a;
              ++count;
            });
  HS_EXPECT_EQ(count, 9);
  HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
}

/**
 * @brief Verifies Blur wraps its column taps around the seam instead of
 *        emitting an out-of-range column.
 * @details The kernel reaches one column either side of the sample, so a plot
 *          at x = 0 must reach W-1 and one at x = W-1 must reach 0. Rows are
 *          interior, so all nine taps fire at their unrenormalized weights.
 */
inline void test_blur_wraps_column_taps() {
  constexpr int W = 32, H = 32;
  Filter::Screen::Blur<W, H> blur(1.0f);
  const float in_alpha = 0.5f;

  auto expect_seam = [&](int column, int left, int right) {
    HS_CONTEXT("column", column);
    bool hit[W] = {};
    bool all_in_range = true;
    int count = 0;
    float sum = 0.0f;
    blur.plot(static_cast<float>(column), 16.0f, Pixel(1, 1, 1), 0.0f, in_alpha,
              [&](float x, float, const Pixel &, float, float a) {
                ++count;
                sum += a;
                if (x < 0.0f || x >= static_cast<float>(W))
                  all_in_range = false;
                else
                  hit[static_cast<int>(x)] = true;
              });
    HS_EXPECT_EQ(count, 9);
    HS_EXPECT_TRUE(all_in_range);
    HS_EXPECT_TRUE(hit[left]);
    HS_EXPECT_TRUE(hit[column]);
    HS_EXPECT_TRUE(hit[right]);
    HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
  };

  expect_seam(0, W - 1, 1);
  expect_seam(W - 1, W - 2, 0);
}

/** @brief Constructor and update clamp strength to the nearest endpoint kernel. */
inline void test_blur_clamps_factor() {
  constexpr int W = 32, H = 32;
  for (float factor : {-3.0f, 2.0f}) {
    for (bool update : {false, true}) {
      Filter::Screen::Blur<W, H> blur(update ? 0.5f : factor);
      if (update)
        blur.update(factor);
      float weights[3][3] = {};
      blur.plot(15, 16, Pixel(1, 2, 3), 0, 1,
                [&](float x, float y, const Pixel &, float, float alpha) {
                  weights[static_cast<int>(y) - 15][static_cast<int>(x) - 14] +=
                      alpha;
                });
      for (int y = 0; y < 3; ++y)
        for (int x = 0; x < 3; ++x) {
          float expected =
              factor < 0 ? ((x == 1 && y == 1) ? 1.0f : 0.0f)
                         : ((x == 1 && y == 1)
                                ? 0.25f
                                : ((x == 1 || y == 1) ? 0.125f : 0.0625f));
          HS_EXPECT_EQ(weights[y][x], expected);
        }
    }
  }
}

/**
 * @brief Verifies update() rebuilds the kernel: update(0) collapses a full blur
 *        back to identity.
 */
inline void test_blur_update_changes_kernel() {
  constexpr int W = 32, H = 32;
  Filter::Screen::Blur<W, H> blur(1.0f);

  blur.update(0.0f);
  int count = 0;
  blur.plot(15.0f, 16.0f, Pixel(1, 1, 1), 0.0f, 1.0f,
            [&](float, float, const Pixel &, float, float) { ++count; });
  HS_EXPECT_EQ(count, 1);
}

/**
 * @brief Verifies the 3x3 kernel's center, edge and corner weights at full and
 *        half strength.
 */
inline void test_blur_kernel_weights_by_offset() {
  constexpr int W = 32, H = 32;
  constexpr float cx = 15.0f, cy = 16.0f;
  struct Case {
    float factor, center, edge, corner;
  };
  constexpr std::array<Case, 2> cases{{
      {1.0f, 0.25f, 0.125f, 0.0625f},
      {0.5f, 0.625f, 0.0625f, 0.03125f},
  }};
  for (const Case &k : cases) {
    Filter::Screen::Blur<W, H> blur(k.factor);
    int count = 0;
    float sum = 0.0f;
    blur.plot(cx, cy, Pixel(1, 1, 1), 0.0f, 1.0f,
              [&](float x, float y, const Pixel &, float, float a) {
                const int manhattan =
                    static_cast<int>(std::fabs(x - cx) + std::fabs(y - cy));
                HS_EXPECT_LE(manhattan, 2);
                const float expected = manhattan == 0   ? k.center
                                       : manhattan == 1 ? k.edge
                                                        : k.corner;
                HS_EXPECT_NEAR(a, expected, 1e-6f);
                ++count;
                sum += a;
              });
    HS_EXPECT_EQ(count, 9);
    HS_EXPECT_NEAR(sum, 1.0f, 1e-5f);
  }
}

/**
 * @brief Verifies the pole-clip renormalization: a full-blur tap on a pole row
 *        drops the off-pole neighbor row (poles don't wrap) yet still deposits
 *        the full input alpha, with no tap landing outside [0, H).
 */
inline void test_blur_pole_row_renormalizes() {
  constexpr int W = 32, H = 32;
  Filter::Screen::Blur<W, H> blur(1.0f);

  const float in_alpha = 0.6f;
  for (float py : {0.0f, static_cast<float>(H - 1)}) {
    float sum = 0.0f;
    int count = 0;
    bool all_in_bounds = true;
    blur.plot(15.0f, py, Pixel(1, 1, 1), 0.0f, in_alpha,
              [&](float, float yy, const Pixel &, float, float a) {
                sum += a;
                ++count;
                if (yy < 0.0f || yy >= static_cast<float>(H))
                  all_in_bounds = false;
              });
    HS_EXPECT_EQ(count, 6); // pole row clipped -> 2 surviving rows, 6 taps
    HS_EXPECT_TRUE(all_in_bounds);
    HS_EXPECT_NEAR(sum, in_alpha, 1e-4f);
  }
}
