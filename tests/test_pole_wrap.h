/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * pole_wrap tap-reflection cases, shared by both offset builds.
 *
 * The south-pole contract differs by hs::H_OFFSET: at 0 the last rendered row
 * is the pole, at the device offset the rows below it are a virtual gap that
 * holds no data. Included by tests/test_geometry.h (offset 0) and
 * tests/test_h_offset_renorm.h (-DHS_TEST_H_OFFSET=3) so both halves run.
 */
#pragma once

#include "core/math/geometry.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace pole_wrap_tests {

/**
 * @brief Verifies in-range taps pass through pole_wrap untouched.
 */
inline void test_pole_wrap_in_range_is_identity() {
  constexpr int W = 64, H = 64;
  int rows[] = {0, 1, 30, H - 1};
  for (int r : rows) {
    int col = 17, row = r;
    HS_EXPECT_TRUE((pole_wrap<W, H>(col, row)));
    HS_EXPECT_EQ(col, 17);
    HS_EXPECT_EQ(row, r);
  }
}

/**
 * @brief Verifies a tap past the north pole reflects onto the far meridian, and
 *        one past the whole buffer reports no data.
 */
inline void test_pole_wrap_north_reflects_half_turn() {
  constexpr int W = 64, H = 64;
  int col = 5, row = -1;
  HS_EXPECT_TRUE((pole_wrap<W, H>(col, row)));
  HS_EXPECT_EQ(row, 1);
  HS_EXPECT_EQ(col, 5 + W / 2);

  col = 40;
  row = -3;
  HS_EXPECT_TRUE((pole_wrap<W, H>(col, row)));
  HS_EXPECT_EQ(row, 3);
  HS_EXPECT_EQ(col, 40 - W / 2);

  col = 5;
  row = -(H + 5);
  HS_EXPECT_TRUE(!(pole_wrap<W, H>(col, row)));
}

/**
 * @brief Verifies the south edge reflects about the true (virtual) pole row, so
 *        the sub-pole gap the device leaves unrendered reports no data.
 */
inline void test_pole_wrap_south_reflects_about_virtual_pole() {
  constexpr int W = 64, H = 64;
  int col = 9, row = H;
  const bool live = pole_wrap<W, H>(col, row);
  if constexpr (hs::H_OFFSET == 0) {
    // The south pole is the last rendered row, so row H mirrors to H - 2.
    HS_EXPECT_TRUE(live);
    HS_EXPECT_EQ(row, H - 2);
    HS_EXPECT_EQ(col, 9 + W / 2);
  } else {
    // Row H is still inside the virtual gap; nothing is ever rendered there.
    HS_EXPECT_TRUE(!live);
  }

  col = 9;
  row = 4 * H;
  HS_EXPECT_TRUE(!(pole_wrap<W, H>(col, row)));
}

/**
 * @brief Runs every pole_wrap case under the caller's module fixture.
 */
inline void run_pole_wrap_cases() {
  test_pole_wrap_in_range_is_identity();
  test_pole_wrap_north_reflects_half_turn();
  test_pole_wrap_south_reflects_about_virtual_pole();
}

} // namespace pole_wrap_tests
} // namespace hs_test
