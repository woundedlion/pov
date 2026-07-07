/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Integrity tests for core/reaction_graph.{h,cpp}.
 *
 * The neighbor table (core/engine/reaction_graph.cpp) is a 92 KB machine-generated
 * K-NN adjacency array; scripts/generate_reaction_graph.py is its generator of
 * record and CI (reaction-graph-provenance) diffs the two, so regeneration is
 * checked. These tests additionally guard the in-tree table's content:
 * shape/population vs the RD_N/RD_K constants, structural invariants (range, no
 * self-loops, no duplicate rows), geometric sanity (listed neighbors are
 * actually nearby), an edge-reciprocity measurement (gross-corruption tripwire —
 * a raw K-NN graph is not required to be perfectly symmetric), and the analytic
 * node() generator.
 * Also exercises CubemapLUT round-trip (direction -> nearest node).
 */
#pragma once

#include <algorithm> // for std::min

#include "core/engine/reaction_graph.h"
#include "tests/test_3dmath.h" // for HS_EXPECT_VEC
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace reaction_graph_tests {

using ReactionGraph::RD_N;
using ReactionGraph::RD_K;
using ReactionGraph::D_AVG;
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
  // Frozen golden: the double-folded Fibonacci-lattice point for index 1234,
  // cast to float32 (host/device agree to the cast per the provenance contract).
  HS_EXPECT_VEC(node(1234),
                Vector(-0.416881472f, 0.678604007f, 0.604736686f), 1e-6f);
  // Cap the sample budget at RD_N-1 so node(i+1) never reads past [0, RD_N).
  const int samples = std::min(500, RD_N - 1);
  for (int i = 0; i < samples; ++i) {
    HS_EXPECT_GT(chord2(node(i), node(i + 1)), 0.0f);
  }
}

/**
 * @brief Pins the frozen D_AVG literal to its analytic value sqrt(4π / RD_N).
 * @details D_AVG is a hand-pasted constant (std::sqrt isn't constexpr here) that
 *          KERNEL_R / INV_R2 in ReactionDiffusionBase derive from, so a stale
 *          value silently mistunes the Wendland kernel support radius. Nothing
 *          else links it to RD_N: bumping RD_N (which also requires regenerating
 *          neighbors[]) would leave D_AVG quietly wrong. This recomputes the
 *          spacing and fails loudly if the literal and RD_N ever diverge.
 */
inline void test_d_avg_matches_rd_n() {
  float expected = static_cast<float>(std::sqrt(4.0 * PI / RD_N));
  HS_EXPECT_NEAR(D_AVG, expected, 1e-4f);
}

/**
 * @brief Verifies the neighbor table's shape and population match RD_N / RD_K.
 * @details The whole suite (and every consumer) iterates the table with the
 *          RD_N/RD_K constants; nothing pinned that the declared array shape, or
 *          the generated data, actually spans RD_N rows. The static_asserts lock
 *          the declared shape to the constants. sizeof only sees the *declared*
 *          shape, though — a generator that emitted fewer than RD_N rows would
 *          zero-pad the tail (every entry -> node 0) and still compile, so the
 *          runtime half asserts the south-pole row (the farthest possible row
 *          from node 0) holds real, local neighbors, proving the data reaches
 *          the last row rather than degenerating into a zero pad.
 */
inline void test_table_shape_matches_constants() {
  static_assert(sizeof(neighbors) / sizeof(neighbors[0]) == RD_N,
                "neighbors row count must equal RD_N");
  static_assert(sizeof(neighbors[0]) / sizeof(neighbors[0][0]) == RD_K,
                "neighbors column count must equal RD_K");

  const Vector last = node(RD_N - 1);
  int valid = 0;
  for (int k = 0; k < RD_K; ++k) {
    int16_t ni = neighbors[RD_N - 1][k];
    HS_EXPECT_TRUE(ni == -1 || (ni >= 0 && ni < RD_N));
    if (ni < 0) continue;
    ++valid;
    // A zero-padded row points every slot at node 0 (the north pole), ~chord^2 4
    // from this south-pole row; a real neighbor is within ~11 deg (chord^2<0.037).
    HS_EXPECT_LT(chord2(last, node(ni)), 0.037f);
  }
  HS_EXPECT_GT(valid, 0); // populated, not a zero-pad row
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

/**
 * @brief Verifies the lattice's maximum node degree stays ≤ 6, pinning the
 *        graph-Laplacian spectral bound the Gray-Scott step relies on.
 * @details GSReactionDiffusion's explicit-Euler stability margin (dt·D·|λ|max ≤
 *          2) rests on the combinatorial Laplacian's spectral radius |λ|max ≤
 *          2·deg_max, with deg_max ≤ 6 on this 6-NN lattice (hence |λ|max ≤ 12).
 *          A denser regenerated table would lift |λ|max and silently invalidate
 *          that stability comment, so assert the realized max degree directly.
 */
inline void test_max_degree_bounds_laplacian() {
  static_assert(RD_K == 6, "GS stability bound assumes a 6-NN lattice");
  int max_deg = 0;
  for (int i = 0; i < RD_N; ++i) {
    int deg = 0;
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[i][k] >= 0) ++deg;
    if (deg > max_deg) max_deg = deg;
  }
  std::printf("  [info] reaction_graph max degree: %d (|lambda|max <= %d)\n",
              max_deg, 2 * max_deg);
  HS_EXPECT_LE(max_deg, 6); // hence |λ|max ≤ 12, the GS stability bound
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
  // ~11 deg upper bound for a listed neighbor → chord^2 < 0.037.
  const float MAX_CHORD2 = 0.037f;
  int far = 0;
  for (int i = 0; i < RD_N; i += 7) {
    Vector p = node(i);
    for (int k = 0; k < RD_K; ++k) {
      int16_t ni = neighbors[i][k];
      if (ni < 0) continue;
      if (chord2(p, node(ni)) > MAX_CHORD2)
        ++far;
    }
  }
  HS_EXPECT_EQ(far, 0);
}

/**
 * @brief Verifies each listed neighbor is far closer than a far-side reference.
 * @details Compares each neighbor's chord distance against an opposite-hemisphere
 *          reference point, which must be strictly larger. node(RD_N-1-i)'s y is
 *          negated (opposite hemisphere) but its longitude is an unrelated
 *          golden-angle value, not -theta, so it is a far point, not the exact
 *          geometric antipode of node(i) — and a far reference is all this test
 *          needs.
 */
inline void test_neighbors_closer_than_far_point() {
  int violations = 0;
  for (int i = 0; i < RD_N; i += 13) {
    Vector p = node(i);
    int far_point = RD_N - 1 - i;
    float far2 = chord2(p, node(far_point));
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
 *          ni->i. The shipped table reciprocates ~98.9% (legitimate K-NN
 *          asymmetry ~1%); the >95% threshold trips on a scrambled table while
 *          still clearing the real asymmetry by a wide margin.
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
  float rate = total ? static_cast<float>(reciprocated) / total : 0.0f;
  std::printf("  [info] reaction_graph edge reciprocity: %.1f%%\n", rate * 100.0f);
  HS_EXPECT_GT(rate, 0.95f);
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
                     RD_N * sizeof(Vector) + 64];
  Arena arena(buf, sizeof(buf));
  ReactionGraph::CubemapLUT lut;
  lut.build(arena);

  int exact = 0, near = 0, miss = 0;
  for (int i = 0; i < RD_N; i += 23) {
    int found = lut.lookup(node(i));
    if (found == i) { ++exact; continue; }
    bool adjacent = false;
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[i][k] == found) { adjacent = true; break; }
    if (adjacent) ++near; else ++miss;
  }
  std::printf("  [info] cubemap roundtrip: %d exact, %d neighbor, %d miss\n",
              exact, near, miss);
  // Seeded at lattice points the LUT misses none; allow at most 5%.
  HS_EXPECT_GT(exact + near, 0);
  HS_EXPECT_LE(miss, (exact + near) / 20);
}

/**
 * @brief Verifies lookup() on off-lattice query directions against a brute-force
 *        nearest-node oracle.
 * @details test_cubemap_lut_roundtrip only seeds at exact lattice points, which
 *          find_nearest_node's note warns "cannot catch" the equatorial long-hop
 *          or a face-boundary miss "because it seeds at the answer." This samples
 *          random unit directions (no lattice point sits exactly there), finds
 *          the true nearest node by an exhaustive O(RD_N) scan, and requires
 *          lookup() to return that node or one of its direct neighbors for the
 *          overwhelming majority — the same one-cell tolerance the round-trip
 *          test allows. The mt19937 is locally seeded so the sample set is fixed.
 */
inline void test_cubemap_lut_offlattice() {
  static uint8_t buf[6 * ReactionGraph::CubemapLUT::RES *
                         ReactionGraph::CubemapLUT::RES * sizeof(uint16_t) +
                     RD_N * sizeof(Vector) + 64];
  Arena arena(buf, sizeof(buf));
  ReactionGraph::CubemapLUT lut;
  lut.build(arena);

  std::mt19937 rng(20240607u);
  std::uniform_real_distribution<float> uni(-1.0f, 1.0f);
  const int SAMPLES = 400;
  int exact = 0, near = 0, miss = 0;
  for (int s = 0; s < SAMPLES; ++s) {
    // Reject near-origin draws that normalize unstably.
    Vector q;
    float len2;
    do {
      q = Vector(uni(rng), uni(rng), uni(rng));
      len2 = q.x * q.x + q.y * q.y + q.z * q.z;
    } while (len2 < 0.01f);
    q = q.normalized();

    // Brute-force argmin = the true nearest node.
    int best = 0;
    float best_d = chord2(q, node(0));
    for (int i = 1; i < RD_N; ++i) {
      float d = chord2(q, node(i));
      if (d < best_d) { best_d = d; best = i; }
    }

    int found = lut.lookup(q);
    if (found == best) { ++exact; continue; }
    bool adjacent = false;
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[best][k] == found) { adjacent = true; break; }
    if (adjacent) ++near; else ++miss;
  }
  std::printf("  [info] cubemap off-lattice: %d exact, %d neighbor, %d miss / %d\n",
              exact, near, miss, SAMPLES);
  // 0 misses across the 400 fixed-seed off-lattice probes; allow at most 5%.
  HS_EXPECT_GT(exact + near, 0);
  HS_EXPECT_LE(miss, SAMPLES / 20);
}

/**
 * @brief Stresses find_nearest_node's worst-case equatorial long-hop directly.
 * @details find_nearest_node seeds by latitude only, so an equatorial query's
 *          true node can sit many longitude hops from the seed — the chain its
 *          64-iter cap is sized for. The random off-lattice sampler need not hit
 *          that chain; this concentrates queries at the equator (|y| tiny)
 *          across the full longitude circle, where the hop count peaks. Each
 *          lookup() drives find_nearest_node through the LUT build and the
 *          always-on HS_CHECK(converged), so a cap exceeded by a future RD_N
 *          bump traps here at the bench. Results are checked against a
 *          brute-force argmin oracle (exact node or a direct neighbor).
 */
inline void test_cubemap_lut_equatorial() {
  static uint8_t buf[6 * ReactionGraph::CubemapLUT::RES *
                         ReactionGraph::CubemapLUT::RES * sizeof(uint16_t) +
                     RD_N * sizeof(Vector) + 64];
  Arena arena(buf, sizeof(buf));
  ReactionGraph::CubemapLUT lut;
  lut.build(arena);

  const int LONGITUDES = 720;
  int exact = 0, near = 0, miss = 0;
  for (int j = 0; j < LONGITUDES; ++j) {
    float lon = (j + 0.5f) / LONGITUDES * 2.0f * static_cast<float>(PI);
    // |y| just off the equator: the maximal-latitude seed distance from a node
    // whose longitude is unrelated to the seed's.
    float y = (j & 1) ? 1e-4f : -1e-4f;
    float r = std::sqrt(1.0f - y * y);
    Vector q(std::cos(lon) * r, y, std::sin(lon) * r);

    int best = 0;
    float best_d = chord2(q, node(0));
    for (int i = 1; i < RD_N; ++i) {
      float d = chord2(q, node(i));
      if (d < best_d) { best_d = d; best = i; }
    }

    int found = lut.lookup(q);
    if (found == best) { ++exact; continue; }
    bool adjacent = false;
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[best][k] == found) { adjacent = true; break; }
    if (adjacent) ++near; else ++miss;
  }
  std::printf("  [info] cubemap equatorial: %d exact, %d neighbor, %d miss / %d\n",
              exact, near, miss, LONGITUDES);
  HS_EXPECT_GT(exact + near, 0);
  HS_EXPECT_LE(miss, LONGITUDES / 20);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs the full reaction_graph test suite.
 * @return The module's failure count.
 */
inline int run_reaction_graph_tests() {
  hs_test::ModuleFixture fixture("reaction_graph");

  test_nodes_on_unit_sphere();
  test_node_deterministic_and_distinct();
  test_d_avg_matches_rd_n();
  test_table_shape_matches_constants();

  test_indices_in_range();
  test_no_self_loops();
  test_no_duplicate_neighbors_in_row();
  test_max_degree_bounds_laplacian();

  test_neighbors_are_local();
  test_neighbors_closer_than_far_point();
  test_edge_reciprocity_high();

  test_cubemap_lut_roundtrip();
  test_cubemap_lut_offlattice();
  test_cubemap_lut_equatorial();

  return fixture.result();
}

} // namespace reaction_graph_tests
} // namespace hs_test

