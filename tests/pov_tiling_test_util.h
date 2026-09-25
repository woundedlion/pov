/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "tests/test_harness.h"
#include <vector>

namespace hs_test {

/** @brief Requires exactly one write per row in each sampled arm column. */
inline void check_arm_column_coverage(const std::vector<int> &cover, int width,
                                      int rows, int column_a, int column_b) {
  HS_EXPECT_SIZE_OR_RETURN(cover, static_cast<size_t>(width) * rows);
  for (int column = 0; column < width; ++column)
    for (int row = 0; row < rows; ++row) {
      HS_CONTEXT("pixel", column, row);
      HS_EXPECT_EQ(cover[static_cast<size_t>(column) * rows + row],
                   column == column_a || column == column_b ? 1 : 0);
    }
}

} // namespace hs_test
