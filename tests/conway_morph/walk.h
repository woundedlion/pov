/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// Walk policy: the recency-weighted random walk visits every node within a
// bounded leg count and keeps long-run visitation balanced (no node above 3x
// the mean share, none below a quarter of it).
// ---------------------------------------------------------------------------

/** Coverage limit enforced by test_walk_policy_coverage_and_balance. */
constexpr int WALK_COVERAGE_BOUND = 250;

inline void reconcile_seed(int edge, int node, int &held) {
  ConwayGraph::SeedFix fix;
  held = ConwayGraph::reconciled_seed_identity(edge, node, held, fix);
  HS_EXPECT_TRUE(fix != ConwayGraph::SeedFix::INVALID);
}

/**
 * @brief Simulates long walks over several RNG seeds and pins coverage and
 *        per-node share bounds.
 */
inline void test_walk_policy_coverage_and_balance() {
  using namespace ConwayGraph;
  constexpr int LEGS = 10000;
  constexpr int MEAN = LEGS / NUM_NODES;

  for (uint32_t seed : {1u, 2u, 3u, 42u, 1337u}) {
    const int failed_before = hs_test::stats().failed;
    hs::random().seed(seed);
    uint8_t visits[NUM_NODES] = {};
    int counts[NUM_NODES] = {};
    int node = TETRAHEDRON;
    int held = TETRAHEDRON;
    int prev = -1;
    int in_family = 0;
    bool seen[NUM_NODES] = {};
    int seen_count = 1;
    int coverage_leg = -1;
    seen[node] = true;
    record_visit(visits, node);

    for (int leg = 0; leg < LEGS; ++leg) {
      const int e = pick_next_edge(node, prev, in_family, visits,
                                   static_cast<uint32_t>(hs::random()()));
      HS_EXPECT_TRUE(edge_touches(e, node));
      HS_EXPECT_TRUE(e != prev || node_degree(node) == 1);
      reconcile_seed(e, node, held);
      const int next = edge_other_end(e, node);
      in_family = family(next) != family(node) ? 0 : in_family + 1;
      if (ConwayGraph::adopts_seed(EDGES[e], next, EDGES[e].to_node == next))
        held = next;
      node = next;
      prev = e;
      ++counts[node];
      record_visit(visits, node);
      if (!seen[node]) {
        seen[node] = true;
        if (++seen_count == NUM_NODES)
          coverage_leg = leg + 1;
      }
    }

    HS_EXPECT_GT(coverage_leg, 0);
    HS_EXPECT_LE(coverage_leg, WALK_COVERAGE_BOUND);
    int mn = counts[0], mx = counts[0];
    for (int i = 0; i < NUM_NODES; ++i) {
      HS_EXPECT_LE(counts[i], 3 * MEAN);
      HS_EXPECT_GE(counts[i], MEAN / 4);
      mn = std::min(mn, counts[i]);
      mx = std::max(mx, counts[i]);
    }
    if (hs_test::stats().failed == failed_before)
      std::printf("  [walk] seed %u: coverage@%d legs, share max/min = "
                  "%d/%d, max/mean = %.2f\n",
                  seed, coverage_leg, mx, mn, static_cast<double>(mx) / MEAN);
    else
      std::printf("    [walk] seed %u failed (coverage@%d, max %d, min %d)\n",
                  seed, coverage_leg, mx, mn);
  }
}
