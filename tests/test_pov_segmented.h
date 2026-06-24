/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the segmented-POV index math (hardware/pov_segment_map.h,
 * the pure arithmetic the Arduino-only pov_segmented.h driver derives its ISR
 * mapping from). This is the one place index arithmetic reaches the physical
 * LEDs, so an off-by-one silently mis-paints the sphere. Covers per-segment
 * derivation (arm side, top/bottom reversal), the arm-B half-image x offset,
 * and that each segment's (i -> x_col, y) writes tile the canvas exactly once.
 */
#pragma once

#include "hardware/pov_segment_map.h"
#include "tests/test_harness.h"

#include <vector>

namespace hs_test {
namespace pov_segmented_tests {

using pov::segment_map;
using pov::SegmentMap;
using pov::segment_x_col;
using pov::segment_y;

// Compile-time proof the mapping is genuinely constexpr (the driver relies on
// it folding at boot, off the ISR hot path).
static_assert(segment_map(1, 288, 4).y_base == 143);
static_assert(segment_map(1, 288, 4).y_step == -1);
static_assert(segment_x_col(true, 0, 288) == 144);

/**
 * @brief Assert the N segments tile both arm columns exactly once at column @p x.
 * @param S Total LED count across all segments.
 * @param N Number of segments.
 * @param w Canvas width in columns (rotation resolution).
 * @param x Rotation column sampled by arm A (arm B samples (x+w/2)%w).
 * @details At a fixed rotation column, every segment's per-LED write must land
 * on a distinct canvas pixel and the N segments together must tile exactly the
 * two columns the arms sample (x for arm A, (x+w/2)%w for arm B), each row
 * covered once. The S LED writes must map onto 2*ROWS canvas pixels.
 */
inline void check_tiling(int S, int N, int w, int x) {
  const int PPS = S / N;
  const int ROWS = S / 2;
  const int col_a = x;
  const int col_b = (x + w / 2) % w;

  std::vector<int> cover(static_cast<size_t>(w) * ROWS, 0);
  int writes = 0;

  for (int seg = 0; seg < N; ++seg) {
    const SegmentMap m = segment_map(seg, S, N);
    const int x_col = segment_x_col(m.arm_b, x, w);
    HS_EXPECT_TRUE(x_col >= 0 && x_col < w);

    for (int i = 0; i < PPS; ++i) {
      const int y = segment_y(m, i);
      HS_EXPECT_TRUE(y >= 0 && y < ROWS);
      cover[static_cast<size_t>(x_col) * ROWS + y]++;
      ++writes;
    }
  }

  HS_EXPECT_EQ(writes, S);

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
      if (hits != 0 && c != col_a && c != col_b)
        stray_column = true;
    }
  }
  HS_EXPECT_FALSE(double_paint);
  HS_EXPECT_FALSE(stray_column);
  HS_EXPECT_EQ(covered, 2 * ROWS);

  for (int y = 0; y < ROWS; ++y) {
    HS_EXPECT_EQ(cover[static_cast<size_t>(col_a) * ROWS + y], 1);
    HS_EXPECT_EQ(cover[static_cast<size_t>(col_b) * ROWS + y], 1);
  }
}

/**
 * @brief Verify per-segment SegmentMap fields for the 4-segment layout.
 * @details Checks arm side (A/B) and the y_base/y_step that distinguish top
 * strips (ascending from the pole) from bottom strips (descending, reversed),
 * plus that the strip endpoints map to the expected pole/junction rows.
 */
inline void test_segment_derivation() {
  // Phantasm config: N=4, S=288 -> ROWS=144, PPS=72.
  const int S = 288, N = 4;

  const SegmentMap s0 = segment_map(0, S, N); // arm A, top
  HS_EXPECT_FALSE(s0.arm_b);
  HS_EXPECT_EQ(s0.y_base, 0);
  HS_EXPECT_EQ(s0.y_step, 1);

  const SegmentMap s1 = segment_map(1, S, N); // arm A, bottom (reversed)
  HS_EXPECT_FALSE(s1.arm_b);
  HS_EXPECT_EQ(s1.y_base, 143);
  HS_EXPECT_EQ(s1.y_step, -1);

  const SegmentMap s2 = segment_map(2, S, N); // arm B, top
  HS_EXPECT_TRUE(s2.arm_b);
  HS_EXPECT_EQ(s2.y_base, 0);
  HS_EXPECT_EQ(s2.y_step, 1);

  const SegmentMap s3 = segment_map(3, S, N); // arm B, bottom (reversed)
  HS_EXPECT_TRUE(s3.arm_b);
  HS_EXPECT_EQ(s3.y_base, 143);
  HS_EXPECT_EQ(s3.y_step, -1);

  // Strip endpoints: the top strip runs N pole (y=0) -> junction (y=PPS-1);
  // the bottom strip runs S pole (y=ROWS-1) -> junction (y=ROWS-PPS), reversed.
  HS_EXPECT_EQ(segment_y(s0, 0), 0);
  HS_EXPECT_EQ(segment_y(s0, 71), 71);
  HS_EXPECT_EQ(segment_y(s1, 0), 143);
  HS_EXPECT_EQ(segment_y(s1, 71), 72);
}

/**
 * @brief Verify arm-B column offset and wrap behavior.
 * @details Arm-B columns sit half an image (w/2) ahead of the rotation column
 * and wrap modulo w; arm A samples the column unshifted.
 */
inline void test_arm_b_offset() {
  const int w = 288;
  HS_EXPECT_EQ(segment_x_col(false, 0, w), 0);
  HS_EXPECT_EQ(segment_x_col(false, 17, w), 17);
  HS_EXPECT_EQ(segment_x_col(true, 0, w), 144);
  HS_EXPECT_EQ(segment_x_col(true, 144, w), 0);   // wraps to the seam
  HS_EXPECT_EQ(segment_x_col(true, 287, w), 143); // (287 + 144) % 288
}

/**
 * @brief Module entry point: run the segmented-POV cases.
 * @return Number of failures recorded by the module.
 */
inline int run_pov_segmented_tests() {
  auto scope = begin_module("pov_segmented");

  test_segment_derivation();
  test_arm_b_offset();

  // Full-canvas tiling across representative rotation columns, including the
  // x=0 and x=w/2 frame boundaries and the wrap seam.
  for (int x : {0, 1, 143, 144, 287})
    check_tiling(/*S=*/288, /*N=*/4, /*w=*/288, x);

  // N=2: a single (non-reversed) segment owns each arm's full column.
  for (int x : {0, 1, 47, 48, 95})
    check_tiling(/*S=*/288, /*N=*/2, /*w=*/96, x);

  // Small config swept over every rotation column: S=8, N=4 -> ROWS=4, PPS=2.
  for (int x = 0; x < 8; ++x)
    check_tiling(/*S=*/8, /*N=*/4, /*w=*/8, x);

  return end_module(scope);
}

} // namespace pov_segmented_tests
} // namespace hs_test
