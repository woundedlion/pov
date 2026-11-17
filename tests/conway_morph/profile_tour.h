/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Ordered-tour coverage, seed reconciliation and wrap.
 */
#pragma once

#include "tests/conway_morph/walk.h"

#include "tests/test_harness.h"
#include "core/mesh/conway_graph.h"

namespace hs_test::conway_morph_tests {

// ---------------------------------------------------------------------------
// Ordered-tour coverage, seed reconciliation and wrap.
// ---------------------------------------------------------------------------

/**
 * @brief Each ordered-tour cycle covers every node and returns to its start.
 */
inline void test_ordered_tour_full_coverage_and_wrap() {
  using namespace ConwayGraph;
  int node = TETRAHEDRON;
  int held = TETRAHEDRON;
  bool seen[NUM_NODES] = {};
  int seen_count = 1;
  int coverage_leg = -1;
  seen[node] = true;

  for (uint32_t leg = 0; leg < 2u * ORDERED_TOUR_LEN; ++leg) {
    const int e = pick_next_edge_ordered(leg);
    HS_EXPECT_TRUE(edge_touches(e, node));

    node = walk_leg(e, node, held);
    if (!seen[node]) {
      seen[node] = true;
      if (++seen_count == NUM_NODES)
        coverage_leg = static_cast<int>(leg) + 1;
    }

    // Cycle closure: every pass ends back at the registry start state.
    if ((leg + 1) % ORDERED_TOUR_LEN == 0) {
      HS_EXPECT_EQ(node, (int)TETRAHEDRON);
      HS_EXPECT_EQ(held, (int)TETRAHEDRON);
    }
  }
  HS_EXPECT_GT(coverage_leg, 0);
  HS_EXPECT_LE(coverage_leg, ORDERED_TOUR_LEN);
  std::printf("  [tour] %d legs per cycle, full coverage after %d\n",
              ORDERED_TOUR_LEN, coverage_leg);
}

} // namespace hs_test::conway_morph_tests
