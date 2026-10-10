/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Walk policy: shortest-path routing to uniformly random destinations, with a
// recent-sweep exclusion window.
// ---------------------------------------------------------------------------

/** Sweeps within which every node must be swept at least once. */
constexpr int WALK_COVERAGE_BOUND = 120;

inline void reconcile_seed(int edge, int node, int &held) {
  ConwayGraph::SeedFix fix;
  held = ConwayGraph::reconciled_seed_identity(edge, node, held, fix);
  HS_EXPECT_TRUE(fix != ConwayGraph::SeedFix::INVALID);
}

/**
 * @brief Pins hop distances: zero on the diagonal, symmetric, one per edge.
 */
inline void test_walk_hops_metric() {
  using namespace ConwayGraph;
  for (int a = 0; a < NUM_NODES; ++a) {
    HS_EXPECT_EQ(hops(a, a), 0);
    for (int b = 0; b < NUM_NODES; ++b) {
      HS_EXPECT_GE(hops(a, b), 0);
      HS_EXPECT_EQ(hops(a, b), hops(b, a));
    }
  }
  for (const auto &edge : EDGES)
    HS_EXPECT_EQ(hops(edge.from_node, edge.to_node), 1);
  HS_EXPECT_EQ(hops(SNUB_CUBE, SNUB_DODECAHEDRON), 7);
  HS_EXPECT_EQ(hops(TRUNCATED_CUBOCTAHEDRON, TRUNCATED_ICOSIDODECAHEDRON), 5);
}

/**
 * @brief Every route from every node to every destination arrives in exactly
 *        hops() legs, whatever the tie-break draws.
 */
inline void test_walk_routes_are_shortest() {
  using namespace ConwayGraph;
  for (uint32_t seed : {1u, 7u, 42u}) {
    hs::random().seed(seed);
    for (int from = 0; from < NUM_NODES; ++from)
      for (int to = 0; to < NUM_NODES; ++to) {
        if (from == to)
          continue;
        int node = from;
        int legs = 0;
        while (node != to && legs <= NUM_NODES) {
          const int e =
              next_edge_toward(node, to, static_cast<uint32_t>(hs::random()()));
          HS_EXPECT_TRUE(edge_touches(e, node));
          node = edge_other_end(e, node);
          ++legs;
        }
        HS_EXPECT_EQ(node, to);
        HS_EXPECT_EQ(legs, hops(from, to));
      }
  }
}

/**
 * @brief Equal-length routes are both taken: truncated tetrahedron ->
 *        icosahedron runs through either the tetrahedron or the octahedron.
 */
inline void test_walk_route_tie_break_varies() {
  using namespace ConwayGraph;
  bool via[NUM_NODES] = {};
  for (uint32_t r = 0; r < 16; ++r)
    via[edge_other_end(next_edge_toward(TRUNCATED_TETRAHEDRON, ICOSAHEDRON, r),
                       TRUNCATED_TETRAHEDRON)] = true;
  HS_EXPECT_EQ(hops(TRUNCATED_TETRAHEDRON, ICOSAHEDRON), 2);
  HS_EXPECT_TRUE(via[TETRAHEDRON]);
  HS_EXPECT_TRUE(via[OCTAHEDRON]);
}

/**
 * @brief Simulates long routed walks over several RNG seeds: the walk stays
 *        seed-reconcilable, never sweeps a node inside the exclusion window,
 *        sweeps every node within a bounded count, and spreads sweeps
 *        uniformly.
 */
inline void test_walk_policy_coverage_and_balance() {
  using namespace ConwayGraph;
  constexpr int SWEEPS = 18000;
  constexpr int MEAN = SWEEPS / NUM_NODES;

  for (uint32_t seed : {1u, 2u, 3u, 42u, 1337u}) {
    const int failed_before = hs_test::stats().failed;
    hs::random().seed(seed);
    uint8_t recent[RECENT_SWEEPS];
    std::fill(std::begin(recent), std::end(recent), NO_NODE);
    int counts[NUM_NODES] = {};
    int node = TETRAHEDRON;
    int held = TETRAHEDRON;
    bool seen[NUM_NODES] = {};
    int seen_count = 1;
    int coverage_sweep = -1;
    seen[node] = true;
    record_sweep(recent, node);
    long legs = 0;

    for (int sweep = 0; sweep < SWEEPS; ++sweep) {
      const int dest =
          pick_destination(node, recent, static_cast<uint32_t>(hs::random()()));
      HS_EXPECT_NE(dest, node);
      HS_EXPECT_FALSE(recently_swept(recent, dest));
      while (node != dest) {
        const int e =
            next_edge_toward(node, dest, static_cast<uint32_t>(hs::random()()));
        reconcile_seed(e, node, held);
        const int next = edge_other_end(e, node);
        if (adopts_seed(EDGES[e], next, EDGES[e].to_node == next))
          held = next;
        node = next;
        ++legs;
      }
      ++counts[node];
      record_sweep(recent, node);
      if (!seen[node]) {
        seen[node] = true;
        if (++seen_count == NUM_NODES)
          coverage_sweep = sweep + 1;
      }
    }

    HS_EXPECT_GT(coverage_sweep, 0);
    HS_EXPECT_LE(coverage_sweep, WALK_COVERAGE_BOUND);
    int mn = counts[0], mx = counts[0];
    for (int i = 0; i < NUM_NODES; ++i) {
      HS_EXPECT_LE(counts[i], MEAN + MEAN / 10);
      HS_EXPECT_GE(counts[i], MEAN - MEAN / 10);
      mn = std::min(mn, counts[i]);
      mx = std::max(mx, counts[i]);
    }
    if (hs_test::stats().failed == failed_before)
      std::printf("  [walk] seed %u: coverage@%d sweeps, share max/min = "
                  "%d/%d, legs/sweep = %.2f\n",
                  seed, coverage_sweep, mx, mn,
                  static_cast<double>(legs) / SWEEPS);
    else
      std::printf("    [walk] seed %u failed (coverage@%d, max %d, min %d)\n",
                  seed, coverage_sweep, mx, mn);
  }
}
