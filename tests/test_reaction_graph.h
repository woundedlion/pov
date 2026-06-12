/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Integrity tests for core/reaction_graph.{h,cpp}.
 *
 * The neighbor table (core/reaction_graph.cpp) is a 92 KB machine-generated
 * K-NN adjacency array with no in-repo regeneration check. These tests guard it
 * as a trusted black box: structural invariants (range, no self-loops, no
 * duplicate rows), geometric sanity (listed neighbors are actually nearby), an
 * edge-reciprocity measurement (gross-corruption tripwire — a raw K-NN graph is
 * not required to be perfectly symmetric), and the analytic node() generator.
 * Also exercises CubemapLUT round-trip (direction -> nearest node).
 */
#pragma once

#include "core/reaction_graph.h"
#include "tests/test_3dmath.h" // for HS_EXPECT_VEC
#include "tests/test_harness.h"

namespace hs_test {
namespace reaction_graph_tests {

using ReactionGraph::RD_N;
using ReactionGraph::RD_K;
using ReactionGraph::neighbors;
using ReactionGraph::node;

/**
 * @brief Squared chord (Euclidean) distance between two points.
 * @param a First point (unit-sphere coordinates).
 * @param b Second point (unit-sphere coordinates).
 * @return Squared chord distance |a - b|^2 (dimensionless).
 * @details Monotone in arc length for points on the unit sphere, so it orders
 *          neighbors without a sqrt.
 */
static inline float chord2(const Vector &a, const Vector &b) {
  float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

// ---------------------------------------------------------------------------
// node() generator
// ---------------------------------------------------------------------------

/**
 * @brief Verifies node() places every lattice point on the unit sphere with
 *        the endpoints at the poles.
 * @details Index 0 and RD_N-1 must sit at the poles, per the Fibonacci-lattice
 *          generator's contract.
 */
inline void test_nodes_on_unit_sphere() {
  // Sample across the whole lattice; every node must lie on the unit sphere.
  for (int i = 0; i < RD_N; i += 17) {
    Vector v = node(i);
    HS_EXPECT_NEAR(v.length(), 1.0f, 1e-3f);
  }
  // Endpoints sit near the poles (y ~ +1 at i=0, y ~ -1 at i=RD_N-1).
  HS_EXPECT_GT(node(0).y, 0.999f);
  HS_EXPECT_LT(node(RD_N - 1).y, -0.999f);
}

/**
 * @brief Verifies node() is deterministic and yields distinct adjacent points.
 * @details node() is a pure function (same index -> identical point) and never
 *          collapses adjacent indices onto the same point.
 */
inline void test_node_deterministic_and_distinct() {
  HS_EXPECT_VEC(node(1234), node(1234), 0.0f);
  // Adjacent lattice indices are distinct points.
  for (int i = 0; i < 500; ++i) {
    HS_EXPECT_GT(chord2(node(i), node(i + 1)), 0.0f);
  }
}

// ---------------------------------------------------------------------------
// Table structural invariants
// ---------------------------------------------------------------------------

/**
 * @brief Verifies every table entry is the -1 sentinel or a valid node index.
 * @details An out-of-range index would index past the lattice in consumers. The
 *          shipped table contains no sentinels, but the format permits them and
 *          consumers guard `ni < 0`, so both are accepted here.
 */
inline void test_indices_in_range() {
  // Every entry is either the -1 sentinel or a valid node index. (The shipped
  // table contains no sentinels, but the format permits them and consumers
  // guard `ni < 0`, so both are accepted here.)
  int bad = 0;
  for (int i = 0; i < RD_N; ++i)
    for (int k = 0; k < RD_K; ++k) {
      int16_t ni = neighbors[i][k];
      if (!(ni == -1 || (ni >= 0 && ni < RD_N)))
        ++bad;
    }
  HS_EXPECT_EQ(bad, 0);
}

/**
 * @brief Verifies no node lists itself as a neighbor.
 * @details Self-loops would waste a slot and break consumers that assume
 *          distinct adjacency.
 */
inline void test_no_self_loops() {
  int loops = 0;
  for (int i = 0; i < RD_N; ++i)
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[i][k] == i)
        ++loops;
  HS_EXPECT_EQ(loops, 0);
}

/**
 * @brief Verifies each non-sentinel neighbor index appears at most once per row.
 * @details A duplicate would shrink the effective fan-out and hint at a corrupt
 *          table.
 */
inline void test_no_duplicate_neighbors_in_row() {
  int dupes = 0;
  for (int i = 0; i < RD_N; ++i)
    for (int k = 0; k < RD_K; ++k) {
      int16_t a = neighbors[i][k];
      if (a < 0) continue;
      for (int j = k + 1; j < RD_K; ++j)
        if (neighbors[i][j] == a)
          ++dupes;
    }
  HS_EXPECT_EQ(dupes, 0);
}

// ---------------------------------------------------------------------------
// Geometric sanity: listed neighbors must actually be nearby
// ---------------------------------------------------------------------------

/**
 * @brief Verifies every listed neighbor is geometrically nearby its node.
 * @details A shuffled or corrupted table would place neighbors far apart, well
 *          past the expected ~11 deg upper bound.
 */
inline void test_neighbors_are_local() {
  // Mean nearest-neighbor spacing on a 7680-point sphere is ~1.3 deg; even the
  // farthest of the 6 listed neighbors should be well within ~11 deg
  // (chord^2 < 0.037). A shuffled/corrupted table would blow past this.
  const float kMaxChord2 = 0.037f;
  int far = 0;
  for (int i = 0; i < RD_N; i += 7) {
    Vector p = node(i);
    for (int k = 0; k < RD_K; ++k) {
      int16_t ni = neighbors[i][k];
      if (ni < 0) continue;
      if (chord2(p, node(ni)) > kMaxChord2)
        ++far;
    }
  }
  HS_EXPECT_EQ(far, 0);
}

/**
 * @brief Verifies each listed neighbor is far closer than the lattice's far side.
 * @details Compares each neighbor's chord distance against the antipode's, which
 *          must be strictly larger.
 */
inline void test_neighbors_closer_than_antipode() {
  int violations = 0;
  for (int i = 0; i < RD_N; i += 13) {
    Vector p = node(i);
    int antipode = RD_N - 1 - i;
    float far2 = chord2(p, node(antipode));
    for (int k = 0; k < RD_K; ++k) {
      int16_t ni = neighbors[i][k];
      if (ni < 0) continue;
      if (chord2(p, node(ni)) >= far2)
        ++violations;
    }
  }
  HS_EXPECT_EQ(violations, 0);
}

// ---------------------------------------------------------------------------
// Edge reciprocity (gross-corruption tripwire, not a hard symmetry requirement)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies the fraction of reciprocated directed edges stays high.
 * @details Measures what fraction of directed edges i->ni have a return edge
 *          ni->i. A clean Fibonacci-lattice K-NN graph reciprocates the large
 *          majority of edges; the >50% threshold is a conservative tripwire that
 *          tolerates legitimate K-NN asymmetry while catching a scrambled table.
 */
inline void test_edge_reciprocity_high() {
  long total = 0, reciprocated = 0;
  for (int i = 0; i < RD_N; ++i) {
    for (int k = 0; k < RD_K; ++k) {
      int16_t ni = neighbors[i][k];
      if (ni < 0) continue;
      ++total;
      for (int j = 0; j < RD_K; ++j) {
        if (neighbors[ni][j] == i) { ++reciprocated; break; }
      }
    }
  }
  HS_EXPECT_GT(total, 0L);
  // A clean Fibonacci-lattice K-NN graph reciprocates the large majority of
  // edges. Require >50% as a conservative tripwire that tolerates legitimate
  // K-NN asymmetry while catching a scrambled table.
  float rate = total ? static_cast<float>(reciprocated) / total : 0.0f;
  std::printf("  [info] reaction_graph edge reciprocity: %.1f%%\n", rate * 100.0f);
  HS_EXPECT_GT(rate, 0.5f);
}

// ---------------------------------------------------------------------------
// CubemapLUT round-trip
// ---------------------------------------------------------------------------

/**
 * @brief Verifies CubemapLUT maps a node's own direction back to that node.
 * @details Looking up a node's own direction should return that node, or at
 *          worst a direct neighbor (cubemap texel quantization can land one cell
 *          over); the overwhelming majority must be exact or an immediate
 *          neighbor.
 */
inline void test_cubemap_lut_roundtrip() {
  static uint8_t buf[6 * ReactionGraph::CubemapLUT::RES *
                         ReactionGraph::CubemapLUT::RES * sizeof(uint16_t) +
                     64];
  Arena arena(buf, sizeof(buf));
  ReactionGraph::CubemapLUT lut;
  lut.build(arena);

  // Looking up a node's own direction should return that node, or at worst a
  // direct neighbor (cubemap texel quantization can land one cell over).
  int exact = 0, near = 0, miss = 0;
  for (int i = 0; i < RD_N; i += 23) {
    int found = lut.lookup(node(i));
    if (found == i) { ++exact; continue; }
    bool adjacent = false;
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[i][k] == found) { adjacent = true; break; }
    if (adjacent) ++near; else ++miss;
  }
  // The overwhelming majority must be exact or an immediate neighbor.
  HS_EXPECT_GT(exact + near, 0);
  HS_EXPECT_LE(miss, (exact + near) / 4);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs the full reaction_graph test suite.
 * @return The module's failure count.
 */
inline int run_reaction_graph_tests() {
  auto scope = hs_test::begin_module("reaction_graph");

  test_nodes_on_unit_sphere();
  test_node_deterministic_and_distinct();

  test_indices_in_range();
  test_no_self_loops();
  test_no_duplicate_neighbors_in_row();

  test_neighbors_are_local();
  test_neighbors_closer_than_antipode();
  test_edge_reciprocity_high();

  test_cubemap_lut_roundtrip();

  return hs_test::end_module(scope);
}

} // namespace reaction_graph_tests
} // namespace hs_test

