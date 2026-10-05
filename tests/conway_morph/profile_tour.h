/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// test_ordered_tour_full_coverage_and_wrap checks NUM_NODES coverage, legal
// seed reconciliation, and return to the registry start state.
// ---------------------------------------------------------------------------

/**
 * @brief Simulates the walk state machine over two ordered-tour cycles.
 */
inline void test_ordered_tour_full_coverage_and_wrap() {
  using namespace ConwayGraph;
  int node = TETRAHEDRON;
  int held = TETRAHEDRON;
  int prev = -1;
  bool seen[NUM_NODES] = {};
  int seen_count = 1;
  int coverage_leg = -1;
  seen[node] = true;

  for (uint32_t leg = 0; leg < 2u * ORDERED_TOUR_LEN; ++leg) {
    const int e = pick_next_edge_ordered(node, prev, leg);
    HS_EXPECT_TRUE(edge_touches(e, node));

    reconcile_seed(e, node, held);

    const bool reverse = EDGES[e].to_node == node;
    const int arrived = reverse ? EDGES[e].from_node : EDGES[e].to_node;
    if (ConwayGraph::adopts_seed(EDGES[e], arrived, !reverse))
      held = arrived;
    node = arrived;
    prev = e;
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
