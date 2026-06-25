/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the single-board POV index math (hardware/pov_single_map.h,
 * the pure arithmetic the Arduino-only pov_single.h driver derives its show_col()
 * ISR mapping from). This is the one place the single-arm index arithmetic
 * reaches the physical LEDs, so an off-by-one or hemisphere swap silently
 * mis-paints the sphere — the same failure mode the segmented tiling tests
 * (test_pov_segmented.h) were written to prevent. Covers the top/bottom strip
 * split (reversed vs straight), the bottom-half opposite-column x offset, and
 * that the S LED writes tile the two sampled canvas columns exactly once.
 */
#pragma once

#include "hardware/pov_single_map.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <vector>

namespace hs_test {
namespace pov_single_tests {

using pov::strip_bottom_led;
using pov::strip_opposite_col;
using pov::strip_top_led;

// Compile-time proof the mapping folds at compile time (the driver relies on it
// being branchless, off the ISR hot path).
static_assert(strip_top_led(0, 40) == 19);     // y=0 -> last top-half LED
static_assert(strip_top_led(19, 40) == 0);     // y=ROWS-1 -> first LED
static_assert(strip_bottom_led(0, 40) == 20);  // y=0 -> first bottom-half LED
static_assert(strip_bottom_led(19, 40) == 39); // y=ROWS-1 -> last LED
static_assert(strip_opposite_col(0, 96) == 48);
static_assert(strip_opposite_col(48, 96) == 0); // wraps back at the seam

/**
 * @brief At a fixed rotation column x, assert the strip's per-LED writes tile
 * exactly the two canvas columns the arms sample, each row covered once.
 * @param S Total physical LED count on the strip (S = 2*ROWS).
 * @param w Canvas width in columns; the bottom half samples column (x+w/2)%w.
 * @param x Rotation column in [0, w); top half samples this column directly.
 * @details Checks that S LED writes map onto 2*ROWS canvas pixels (top half at
 * column x, bottom half at (x+w/2)%w), every physical LED in [0, S) is written
 * exactly once with no gap or double-drive, and both sampled columns are fully
 * covered row for row.
 */
inline void check_strip_tiling(int S, int w, int x) {
  const int ROWS = S / 2;
  const int col_top = x;
  const int col_bot = strip_opposite_col(x, w);
  HS_EXPECT_TRUE(col_bot >= 0 && col_bot < w);
  // For even w the two sampled columns are always distinct (offset w/2 != 0).
  HS_EXPECT_TRUE(col_top != col_bot);

  std::vector<int> cover(static_cast<size_t>(w) * ROWS, 0);
  std::vector<int> led_hits(static_cast<size_t>(S), 0);
  int writes = 0;

  for (int y = 0; y < ROWS; ++y) {
    const int top_led = strip_top_led(y, S);
    const int bot_led = strip_bottom_led(y, S);
    HS_EXPECT_TRUE(top_led >= 0 && top_led < S);
    HS_EXPECT_TRUE(bot_led >= 0 && bot_led < S);
    led_hits[static_cast<size_t>(top_led)]++;
    led_hits[static_cast<size_t>(bot_led)]++;
    cover[static_cast<size_t>(col_top) * ROWS + y]++;
    cover[static_cast<size_t>(col_bot) * ROWS + y]++;
    writes += 2;
  }

  HS_EXPECT_EQ(writes, S);

  // Every physical LED written exactly once (no gap, no double-drive).
  for (int p = 0; p < S; ++p)
    HS_EXPECT_EQ(led_hits[static_cast<size_t>(p)], 1);

  int covered = 0;
  bool double_paint = false;
  bool stray_column = false;
  for (int c = 0; c < w; ++c) {
    for (int y = 0; y < ROWS; ++y) {
      const int hits = cover[static_cast<size_t>(c) * ROWS + y];
      if (hits > 1)
        double_paint = true;
      if (hits == 1)
        ++covered;
      if (hits != 0 && c != col_top && c != col_bot)
        stray_column = true;
    }
  }
  HS_EXPECT_FALSE(double_paint);
  HS_EXPECT_FALSE(stray_column);
  HS_EXPECT_EQ(covered, 2 * ROWS);

  for (int y = 0; y < ROWS; ++y) {
    HS_EXPECT_EQ(cover[static_cast<size_t>(col_top) * ROWS + y], 1);
    HS_EXPECT_EQ(cover[static_cast<size_t>(col_bot) * ROWS + y], 1);
  }
}

/**
 * @brief Pin down the top/bottom strip split for the production 96x20 config:
 * the top half [0, ROWS) maps reversed, the bottom half straight, and the two
 * halves partition the physical LED range [0, S) with no overlap.
 */
inline void test_strip_derivation() {
  // Holosphere 96x20 config: S=40 -> ROWS=20.
  const int S = 40;

  // Top half [0, 20): reversed. y=0 is the LED nearest the junction (LED 19),
  // y=ROWS-1 the pole end (LED 0).
  HS_EXPECT_EQ(strip_top_led(0, S), 19);
  HS_EXPECT_EQ(strip_top_led(19, S), 0);

  // Bottom half [20, 40): straight. y=0 -> LED 20, y=ROWS-1 -> LED 39.
  HS_EXPECT_EQ(strip_bottom_led(0, S), 20);
  HS_EXPECT_EQ(strip_bottom_led(19, S), 39);

  // The two halves partition [0, S) with no overlap.
  for (int y = 0; y < S / 2; ++y) {
    HS_EXPECT_TRUE(strip_top_led(y, S) < S / 2);
    HS_EXPECT_TRUE(strip_bottom_led(y, S) >= S / 2);
  }
}

/**
 * @brief Verify the bottom-half column offset is (x + w/2) mod w, including the
 * wrap at the seam where x >= w/2 folds back to the low columns.
 */
inline void test_opposite_col_offset() {
  const int w = 96;
  HS_EXPECT_EQ(strip_opposite_col(0, w), 48);
  HS_EXPECT_EQ(strip_opposite_col(17, w), 65);
  HS_EXPECT_EQ(strip_opposite_col(48, w), 0);  // wraps to the seam
  HS_EXPECT_EQ(strip_opposite_col(95, w), 47); // (95 + 48) % 96
}

/**
 * @brief Run the single-board POV index-math suite.
 * @return Number of failed expectations in the suite (0 on full pass).
 * @details Exercises the strip split and column offset directly, then the full
 * tiling invariant across the two production configs and a small config swept
 * over every rotation column.
 */
inline int run_pov_single_tests() {
  hs_test::ModuleFixture fixture("pov_single");

  test_strip_derivation();
  test_opposite_col_offset();

  // Holosphere 96x20: S=40, swept over representative rotation columns including
  // the x=0 and x=w/2 frame boundaries and the wrap seam.
  for (int x : {0, 1, 47, 48, 95})
    check_strip_tiling(/*S=*/40, /*w=*/96, x);

  // Holosphere 288x144: S=288.
  for (int x : {0, 1, 143, 144, 287})
    check_strip_tiling(/*S=*/288, /*w=*/288, x);

  // Small config swept over every rotation column: S=8 -> ROWS=4, w=8.
  for (int x = 0; x < 8; ++x)
    check_strip_tiling(/*S=*/8, /*w=*/8, x);

  return fixture.result();
}

} // namespace pov_single_tests
} // namespace hs_test
