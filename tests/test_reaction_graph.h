/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Integrity tests for core/spatial/reaction_graph.{h,cpp}.
 */
#pragma once

#include <cmath>
#include <cstdio>

#include "core/spatial/reaction_graph.h"
#include "tests/vec_test_util.h"
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
static inline float chord2(const math::Vector &a, const math::Vector &b) {
  float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

/**
 * @brief Upper bound on chord^2 from a node to any listed neighbor.
 * @details About 2x the shipped table's worst edge; a row shifted three rings
 *          exceeds it.
 */
constexpr float MAX_NEIGHBOR_CHORD2 = 0.008f;

// ---------------------------------------------------------------------------
// node() generator
// ---------------------------------------------------------------------------

/** @brief Compares generated flash positions with their analytic reference. */
inline void test_generated_node_positions() {
  for (int i = 0; i < RD_N; ++i) {
    const math::Vector EXPECTED = node(i);
    const math::Vector ACTUAL = ReactionGraph::node_positions[i];
    HS_EXPECT_NEAR(ACTUAL.x, EXPECTED.x, 1e-7f);
    HS_EXPECT_NEAR(ACTUAL.y, EXPECTED.y, 1e-7f);
    HS_EXPECT_NEAR(ACTUAL.z, EXPECTED.z, 1e-7f);
  }
}

/**
 * @brief Verifies node() places every lattice point on the unit sphere with
 *        the endpoints at the poles.
 */
inline void test_nodes_on_unit_sphere() {
  float worst_deviation = 0.0f;
  for (int i = 0; i < RD_N; ++i) {
    worst_deviation = hs_test::fold_worst(worst_deviation,
                                          std::fabs(node(i).length() - 1.0f));
  }
  HS_EXPECT_LT(worst_deviation, 1e-3f);
  HS_EXPECT_GT(node(0).y, 0.999f);
  HS_EXPECT_LT(node(RD_N - 1).y, -0.999f);
}

/**
 * @brief Pins node() to frozen coordinates and verifies the lattice walks
 *        strictly southward with no coincident neighbors.
 */
inline void test_node_ordered_and_distinct() {
  // Frozen double-folded goldens at indices where a float32 theta fold moves
  // x/z by more than the tolerance.
  HS_EXPECT_VEC(node(1234),
                math::Vector(-0.416881472f, 0.678604007f, 0.604736686f), 1e-6f);
  HS_EXPECT_VEC(node(5759),
                math::Vector(-0.0421082303f, -0.499934882f, -0.865038753f),
                1e-6f);
  HS_EXPECT_VEC(node(RD_N - 2),
                math::Vector(-0.00214281073f, -0.999739528f, -0.0227209534f),
                1e-6f);
  math::Vector prev = node(0);
  HS_EXPECT_TRUE(std::isfinite(prev.x) && std::isfinite(prev.y) &&
                 std::isfinite(prev.z));
  int out_of_order = 0;
  int coincident = 0;
  for (int i = 0; i < RD_N - 1; ++i) {
    math::Vector next = node(i + 1);
    HS_EXPECT_TRUE(std::isfinite(next.x) && std::isfinite(next.y) &&
                   std::isfinite(next.z));
    // y = 1 - 2i/(RD_N-1): index order is the north-to-south sweep order.
    out_of_order += next.y >= prev.y;
    coincident += chord2(prev, next) <= 0.0f;
    prev = next;
  }
  HS_EXPECT_EQ(out_of_order, 0);
  HS_EXPECT_EQ(coincident, 0);
}

/**
 * @brief Pins the frozen D_AVG literal to its analytic value sqrt(4π / RD_N).
 */
inline void test_d_avg_matches_rd_n() {
  float expected = static_cast<float>(std::sqrt(4.0 * PI / RD_N));
  HS_EXPECT_NEAR(D_AVG, expected, 1e-6f);
}

/** @brief Expands every run and compares all ordered slots to the K-NN table. */
inline void test_neighbor_runs_match_table() {
  static_assert(sizeof(ReactionGraph::NeighborRun) == 2 * (RD_K + 1));
  HS_EXPECT_GT(ReactionGraph::NEIGHBOR_RUN_COUNT, 0u);
  HS_EXPECT_LE(ReactionGraph::NEIGHBOR_RUN_COUNT, static_cast<unsigned>(RD_N));
  HS_EXPECT_LE(ReactionGraph::NEIGHBOR_RUN_COUNT, 256u);
  static_assert(sizeof(ReactionGraph::neighbor_run_index) == RD_N);
  int start = 0;
  int first_bad_slot = -1;
  for (unsigned r = 0; r < ReactionGraph::NEIGHBOR_RUN_COUNT; ++r) {
    const auto &run = ReactionGraph::neighbor_runs[r];
    HS_EXPECT_GT(run.end, start);
    HS_EXPECT_LE(run.end, RD_N);
    if (run.end <= start || run.end > RD_N)
      return;
    for (int i = start; i < run.end; ++i) {
      HS_EXPECT_EQ(ReactionGraph::neighbor_run_index[i], r);
      for (int k = 0; k < RD_K; ++k)
        if (i + run.delta[k] != neighbors[i][k] && first_bad_slot < 0)
          first_bad_slot = i * RD_K + k;
    }
    start = run.end;
  }
  HS_EXPECT_EQ(start, RD_N);
  HS_EXPECT_EQ(first_bad_slot, -1);
}

// ---------------------------------------------------------------------------
// Table structural invariants
//
// Each case walks all RD_N*RD_K slots and captures the first offending slot.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies every table entry is a valid node index.
 * @details The table has no vacant-slot sentinel.
 */
inline void test_indices_in_range() {
  int first_bad_slot = -1;
  for (int i = 0; i < RD_N; ++i)
    for (int k = 0; k < RD_K; ++k) {
      const int16_t ni = neighbors[i][k];
      if (!(ni >= 0 && ni < RD_N) && first_bad_slot < 0)
        first_bad_slot = i * RD_K + k;
    }
  HS_EXPECT_EQ(first_bad_slot, -1);
}

/**
 * @brief Verifies no node lists itself as a neighbor.
 */
inline void test_no_self_loops() {
  int first_self_slot = -1;
  for (int i = 0; i < RD_N; ++i)
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[i][k] == i && first_self_slot < 0)
        first_self_slot = i * RD_K + k;
  HS_EXPECT_EQ(first_self_slot, -1);
}

/**
 * @brief Verifies each neighbor index appears at most once per row.
 */
inline void test_no_duplicate_neighbors_in_row() {
  int first_duplicate_slot = -1;
  for (int i = 0; i < RD_N; ++i)
    for (int k = 0; k < RD_K; ++k) {
      int16_t a = neighbors[i][k];
      HS_EXPECT_TRUE(a >= 0 && a < RD_N);
      if (a < 0 || a >= RD_N)
        return;
      for (int j = k + 1; j < RD_K; ++j)
        if (neighbors[i][j] == a && first_duplicate_slot < 0)
          first_duplicate_slot = i * RD_K + j;
    }
  HS_EXPECT_EQ(first_duplicate_slot, -1);
}

// ---------------------------------------------------------------------------
// Geometric sanity: listed neighbors must actually be nearby
// ---------------------------------------------------------------------------

/**
 * @brief Verifies every listed neighbor is geometrically nearby its node.
 */
inline void test_neighbors_are_local() {
  int first_far_slot = -1;
  for (int i = 0; i < RD_N; ++i) {
    math::Vector p = node(i);
    for (int k = 0; k < RD_K; ++k) {
      int16_t ni = neighbors[i][k];
      HS_EXPECT_TRUE(ni >= 0 && ni < RD_N);
      if (ni < 0 || ni >= RD_N)
        return;
      const float distance = chord2(p, node(ni));
      if ((!std::isfinite(distance) || distance > MAX_NEIGHBOR_CHORD2) &&
          first_far_slot < 0)
        first_far_slot = i * RD_K + k;
    }
  }
  HS_EXPECT_EQ(first_far_slot, -1);
}

/**
 * @brief Verifies sampled rows are the true RD_K nearest neighbors of node().
 * @details Every STRIDE-th row is rebuilt by brute force over all other nodes
 *          under the generator's (chord^2, index) order and must match slot for
 *          slot.
 */
inline void test_neighbors_match_brute_force_knn() {
  constexpr int STRIDE = 37;
  int first_bad_row = -1;
  for (int i = 0; i < RD_N; i += STRIDE) {
    const math::Vector p = node(i);
    float best2[RD_K];
    int best[RD_K];
    int filled = 0;
    for (int j = 0; j < RD_N; ++j) {
      if (j == i)
        continue;
      const float d2 = chord2(p, node(j));
      if (filled == RD_K && d2 >= best2[RD_K - 1])
        continue;
      int slot = filled < RD_K ? filled : RD_K - 1;
      while (slot > 0 && d2 < best2[slot - 1]) {
        best2[slot] = best2[slot - 1];
        best[slot] = best[slot - 1];
        --slot;
      }
      best2[slot] = d2;
      best[slot] = j;
      if (filled < RD_K)
        ++filled;
    }
    for (int k = 0; k < RD_K && first_bad_row < 0; ++k)
      if (neighbors[i][k] != best[k]) {
        first_bad_row = i;
        std::printf("  [info] row %d slot %d: table %d, brute force %d\n", i, k,
                    neighbors[i][k], best[k]);
      }
  }
  HS_EXPECT_EQ(first_bad_row, -1);
}

// ---------------------------------------------------------------------------
// Edge reciprocity (gross-corruption tripwire, not a hard symmetry requirement)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies the fraction of reciprocated directed edges stays high.
 * @details Measures what fraction of directed edges i->ni have a return edge
 *          ni->i; a K-NN graph is not exactly symmetric.
 */
inline void test_edge_reciprocity_high() {
  long total = 0, reciprocated = 0;
  for (int i = 0; i < RD_N; ++i) {
    for (int k = 0; k < RD_K; ++k) {
      int16_t ni = neighbors[i][k];
      HS_EXPECT_TRUE(ni >= 0 && ni < RD_N);
      if (ni < 0 || ni >= RD_N)
        return;
      ++total;
      for (int j = 0; j < RD_K; ++j) {
        if (neighbors[ni][j] == i) {
          ++reciprocated;
          break;
        }
      }
    }
  }
  HS_EXPECT_GT(total, 0L);
  float rate = total ? static_cast<float>(reciprocated) / total : 0.0f;
  std::printf("  [info] reaction_graph edge reciprocity: %.1f%%\n",
              rate * 100.0f);
  HS_EXPECT_GT(rate, 0.95f);
}

// ---------------------------------------------------------------------------
// CubemapLUT round-trip
// ---------------------------------------------------------------------------

inline const ReactionGraph::CubemapLUT &built_cubemap_lut() {
  struct Fixture {
    uint8_t buffer[6 * ReactionGraph::CubemapLUT::RES *
                       ReactionGraph::CubemapLUT::RES * sizeof(uint16_t) +
                   RD_N * sizeof(math::Vector) + 64];
    Arena arena;
    ReactionGraph::CubemapLUT lut;
    Fixture() : arena(buffer, sizeof(buffer)) { lut.build(arena); }
  };
  static Fixture fixture;
  return fixture.lut;
}

enum class LookupClass { EXACT, NEIGHBOR, MISS };

inline LookupClass classify_lookup(const ReactionGraph::CubemapLUT &lut,
                                   const math::Vector &q) {
  int best = 0;
  float best_distance = chord2(q, node(0));
  for (int i = 1; i < RD_N; ++i) {
    const float DISTANCE = chord2(q, node(i));
    if (DISTANCE < best_distance) {
      best_distance = DISTANCE;
      best = i;
    }
  }
  const int FOUND = lut.lookup(q);
  if (FOUND == best)
    return LookupClass::EXACT;
  for (int k = 0; k < RD_K; ++k)
    if (neighbors[best][k] == FOUND)
      return LookupClass::NEIGHBOR;
  return LookupClass::MISS;
}

/**
 * @brief Every 23rd lattice node's direction maps back to that node or a
 *        direct neighbor.
 * @details Cubemap texel quantization can land one cell over.
 */
inline void test_cubemap_lut_roundtrip() {
  const auto &lut = built_cubemap_lut();

  int exact = 0, near = 0, miss = 0;
  for (int i = 0; i < RD_N; i += 23) {
    int found = lut.lookup(node(i));
    if (found == i) {
      ++exact;
      continue;
    }
    bool adjacent = false;
    for (int k = 0; k < RD_K; ++k)
      if (neighbors[i][k] == found) {
        adjacent = true;
        break;
      }
    if (adjacent)
      ++near;
    else
      ++miss;
  }
  std::printf("  [info] cubemap roundtrip: %d exact, %d neighbor, %d miss\n",
              exact, near, miss);
  HS_EXPECT_GT(exact + near, 0);
  HS_EXPECT_EQ(miss, 0);
}

/**
 * @brief Verifies lookup() on off-lattice query directions against a brute-force
 *        nearest-node oracle.
 * @details Allows one neighbor of error.
 */
inline void test_cubemap_lut_offlattice() {
  const auto &lut = built_cubemap_lut();

  hs::Pcg32 rng(20240607u);
  const int SAMPLES = 400;
  int exact = 0, near = 0, miss = 0;
  for (int s = 0; s < SAMPLES; ++s) {
    // Reject near-origin draws that normalize unstably.
    math::Vector q;
    float len2;
    do {
      const float x = rand_uniform(rng, -1.0f, 1.0f);
      const float y = rand_uniform(rng, -1.0f, 1.0f);
      const float z = rand_uniform(rng, -1.0f, 1.0f);
      q = math::Vector(x, y, z);
      len2 = q.x * q.x + q.y * q.y + q.z * q.z;
    } while (len2 < 0.01f);
    q = q.normalized();

    const LookupClass CLASSIFICATION = classify_lookup(lut, q);
    if (CLASSIFICATION == LookupClass::EXACT)
      ++exact;
    else if (CLASSIFICATION == LookupClass::NEIGHBOR)
      ++near;
    else
      ++miss;
  }
  std::printf(
      "  [info] cubemap off-lattice: %d exact, %d neighbor, %d miss / %d\n",
      exact, near, miss, SAMPLES);
  HS_EXPECT_GT(exact + near, 0);
  HS_EXPECT_EQ(miss, 0);
}

/**
 * @brief Checks equatorial cubemap lookups against a brute-force oracle.
 * @details Equatorial lookup quantization permits at most 12 probes outside the
 * oracle node and its direct neighbors among 720 fixed queries.
 */
inline void test_cubemap_lut_equatorial() {
  const auto &lut = built_cubemap_lut();

  const int LONGITUDES = 720;
  int exact = 0, near = 0, miss = 0;
  for (int j = 0; j < LONGITUDES; ++j) {
    float lon = (j + 0.5f) / LONGITUDES * 2.0f * static_cast<float>(PI);
    float y = (j & 1) ? 1e-4f : -1e-4f;
    float r = std::sqrt(1.0f - y * y);
    math::Vector q(std::cos(lon) * r, y, std::sin(lon) * r);

    const LookupClass CLASSIFICATION = classify_lookup(lut, q);
    if (CLASSIFICATION == LookupClass::EXACT)
      ++exact;
    else if (CLASSIFICATION == LookupClass::NEIGHBOR)
      ++near;
    else
      ++miss;
  }
  std::printf(
      "  [info] cubemap equatorial: %d exact, %d neighbor, %d miss / %d\n",
      exact, near, miss, LONGITUDES);
  HS_EXPECT_GT(exact + near, 0);
  HS_EXPECT_LE(miss, 12);
}

inline void test_cubemap_coherent_seeds() {
  const auto &LUT = built_cubemap_lut();
  constexpr int RES = ReactionGraph::CubemapLUT::RES;
  int changed = 0, misses = 0;
  for (int face = 0; face < 6; ++face)
    for (int y = 0; y < RES; ++y)
      for (int x = 0; x < RES; ++x) {
        const float U = (x + 0.5f) / RES * 2.0f - 1.0f;
        const float V = (y + 0.5f) / RES * 2.0f - 1.0f;
        math::Vector q;
        switch (face) {
        case 0:
          q = {1, V, -U};
          break;
        case 1:
          q = {-1, V, U};
          break;
        case 2:
          q = {U, 1, -V};
          break;
        case 3:
          q = {U, -1, V};
          break;
        case 4:
          q = {U, V, 1};
          break;
        default:
          q = {-U, V, -1};
          break;
        }
        q = q.normalized();
        int baseline =
            static_cast<int>(hs::clamp((1.0f - q.y) * 0.5f * (RD_N - 1) + 0.5f,
                                       0.0f, static_cast<float>(RD_N - 1)));
        float distance = chord2(q, ReactionGraph::node_positions[baseline]);
        for (int iter = 0; iter < 64; ++iter) {
          bool improved = false;
          for (int k = 0; k < RD_K; ++k) {
            const int NEIGHBOR = neighbors[baseline][k];
            const float D = chord2(q, ReactionGraph::node_positions[NEIGHBOR]);
            if (D < distance) {
              baseline = NEIGHBOR;
              distance = D;
              improved = true;
            }
          }
          if (!improved)
            break;
        }
        const int FOUND =
            LUT.lookup(ReactionGraph::CubemapLUT::Projection{face, U, V});
        changed += FOUND != baseline;
        int nearest = 0;
        float nearest_distance = chord2(q, ReactionGraph::node_positions[0]);
        for (int i = 1; i < RD_N; ++i) {
          const float D = chord2(q, ReactionGraph::node_positions[i]);
          if (D < nearest_distance) {
            nearest_distance = D;
            nearest = i;
          }
        }
        bool adjacent = FOUND == nearest;
        for (int k = 0; k < RD_K; ++k)
          adjacent |= neighbors[nearest][k] == FOUND;
        misses += !adjacent;
      }
  std::printf("  [info] coherent cubemap: %d changed, %d misses\n", changed,
              misses);
  HS_EXPECT_EQ(changed, 0);
  HS_EXPECT_EQ(misses, 0);
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
  test_generated_node_positions();
  test_node_ordered_and_distinct();
  test_d_avg_matches_rd_n();
  test_neighbor_runs_match_table();

  test_indices_in_range();
  test_no_self_loops();
  test_no_duplicate_neighbors_in_row();

  test_neighbors_are_local();
  test_neighbors_match_brute_force_knn();
  test_edge_reciprocity_high();

  test_cubemap_lut_roundtrip();
  test_cubemap_lut_offlattice();
  test_cubemap_lut_equatorial();
  test_cubemap_coherent_seeds();

  return fixture.result();
}

} // namespace reaction_graph_tests
} // namespace hs_test
