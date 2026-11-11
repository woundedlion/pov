/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for KDTree.
 */
#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include "core/spatial/kd_tree.h"
#include "tests/vec_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace spatial_tests {

// Module-scope scratch buffer, re-based by a fresh Arena in each test; nothing
// pointing into it outlives its test.
inline constexpr size_t SPATIAL_BUF_BYTES = 128 * 1024;
inline uint8_t spatial_buf[SPATIAL_BUF_BYTES];

// ============================================================================
// KDTree
// ============================================================================

/**
 * @brief Verifies building from a zero-length span leaves the tree unbuilt
 *        and queries return nothing.
 */
inline void test_kdtree_empty_input() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[1] = {math::Vector(0, 0, 0)};
  std::span<math::Vector> empty(pts, 0);
  KDTree tree(arena, empty);
  HS_EXPECT_TRUE(tree.empty());

  auto r = tree.nearest(math::Vector(0, 0, 0), 1);
  HS_EXPECT_TRUE(r.empty());
}

/**
 * @brief Verifies a one-point tree builds a root and returns that point (with
 *        its original index) for any query.
 */
inline void test_kdtree_single_point() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[1] = {math::Vector(3, 4, 5)};
  std::span<math::Vector> sp(pts, 1);
  KDTree tree(arena, sp);
  HS_EXPECT_FALSE(tree.empty());

  auto r = tree.nearest(math::Vector(0, 0, 0), 1);
  HS_EXPECT_SIZE_OR_RETURN(r, 1);
  HS_EXPECT_VEC(r[0].point, math::Vector(3, 4, 5), 1e-6f);
  HS_EXPECT_EQ(r[0].original_index, (uint16_t)0);
}

/** @brief Moving a built tree leaves the source empty and preserves queries. */
inline void test_kdtree_move() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[] = {math::Vector(3, 4, 5)};
  KDTree source(arena, std::span<math::Vector>(pts));
  KDTree destination(std::move(source));
  HS_EXPECT_TRUE(source.empty());
  HS_EXPECT_TRUE(source.nearest(math::Vector(0, 0, 0)).empty());
  auto result = destination.nearest(math::Vector(0, 0, 0));
  HS_EXPECT_SIZE_OR_RETURN(result, 1);
  HS_EXPECT_VEC(result[0].point, pts[0], 1e-6f);
}

/**
 * @brief Verifies nearest(k=1) returns the correct point and original index for
 *        queries near distinct points.
 * @details Includes a query landing exactly on a point (dist 0).
 */
inline void test_kdtree_nearest_known_set() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[6] = {
      math::Vector(0, 0, 0),  // 0
      math::Vector(10, 0, 0), // 1
      math::Vector(0, 10, 0), // 2
      math::Vector(0, 0, 10), // 3
      math::Vector(5, 5, 5),  // 4
      math::Vector(-3, 2, 1), // 5
  };
  std::span<math::Vector> sp(pts, 6);
  KDTree tree(arena, sp);

  auto r1 = tree.nearest(math::Vector(0.1f, 0.1f, 0.1f), 1);
  HS_EXPECT_SIZE_OR_RETURN(r1, 1);
  HS_EXPECT_VEC(r1[0].point, math::Vector(0, 0, 0), 1e-6f);
  HS_EXPECT_EQ(r1[0].original_index, (uint16_t)0);

  auto r2 = tree.nearest(math::Vector(10.1f, 0, 0), 1);
  HS_EXPECT_SIZE_OR_RETURN(r2, 1);
  HS_EXPECT_VEC(r2[0].point, math::Vector(10, 0, 0), 1e-6f);
  HS_EXPECT_EQ(r2[0].original_index, (uint16_t)1);

  auto r3 = tree.nearest(math::Vector(5, 5, 5), 1);
  HS_EXPECT_SIZE_OR_RETURN(r3, 1);
  HS_EXPECT_VEC(r3[0].point, math::Vector(5, 5, 5), 1e-6f);
  HS_EXPECT_EQ(r3[0].original_index, (uint16_t)4);
}

/**
 * @brief Verifies a k>1 query returns exactly k neighbors ordered
 *        nearest-first.
 */
inline void test_kdtree_k_nearest_sorted() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[5] = {
      math::Vector(0, 0, 0), // 0  d²=0
      math::Vector(1, 0, 0), // 1  d²=1
      math::Vector(2, 0, 0), // 2  d²=4
      math::Vector(3, 0, 0), // 3  d²=9
      math::Vector(4, 0, 0), // 4  d²=16
  };
  std::span<math::Vector> sp(pts, 5);
  KDTree tree(arena, sp);

  auto r = tree.nearest(math::Vector(0, 0, 0), 3);
  HS_EXPECT_SIZE_OR_RETURN(r, 3);
  HS_EXPECT_VEC(r[0].point, math::Vector(0, 0, 0), 1e-6f);
  HS_EXPECT_VEC(r[1].point, math::Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(r[2].point, math::Vector(2, 0, 0), 1e-6f);
}

/** @brief Equidistant KD-tree results prefer the lower input index. */
inline void test_kdtree_ties_prefer_source_index() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[] = {math::Vector(1, 0, 0), math::Vector(-1, 0, 0),
                        math::Vector(0, 1, 0), math::Vector(0, -1, 0)};
  KDTree tree(arena, std::span<math::Vector>(pts));

  const auto result = tree.nearest(math::Vector(0, 0, 0), 2);
  HS_EXPECT_SIZE_OR_RETURN(result, 2);
  HS_EXPECT_EQ(result[0].original_index, 0);
  HS_EXPECT_EQ(result[1].original_index, 1);
}

/**
 * @brief Oversized neighbor requests return every source point in nearest-first order.
 */
inline void test_kdtree_k_caps_at_size() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[3] = {math::Vector(0, 0, 0), math::Vector(1, 0, 0),
                         math::Vector(2, 0, 0)};
  std::span<math::Vector> sp(pts, 3);
  KDTree tree(arena, sp);

  auto r = tree.nearest(math::Vector(0, 0, 0), 5);
  HS_EXPECT_SIZE_OR_RETURN(r, 3);
  for (size_t i = 0; i < r.size(); ++i) {
    HS_EXPECT_EQ(r[i].original_index, static_cast<int>(i));
    HS_EXPECT_EQ(r[i].d_sq, static_cast<float>(i * i));
  }
}

/**
 * @brief Verifies a default-constructed (never-built) KDTree reports empty and
 *        answers queries with an empty result.
 */
inline void test_kdtree_default_unbuilt() {
  KDTree tree;
  HS_EXPECT_TRUE(tree.empty());
  auto r = tree.nearest(math::Vector(1, 2, 3), 1);
  HS_EXPECT_TRUE(r.empty());
}

/**
 * @brief Verifies the KDTree nearest neighbor matches a brute-force scan.
 */
inline void test_kdtree_matches_brute_force() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  constexpr int N = 16;
  math::Vector pts[N];
  for (int i = 0; i < N; ++i) {
    float fi = static_cast<float>(i);
    pts[i] = math::Vector(std::sin(fi) * 5.0f, std::cos(fi * 1.3f) * 4.0f,
                          std::sin(fi * 0.7f) * 3.0f);
  }

  std::span<math::Vector> sp(pts, N);
  KDTree tree(arena, sp);

  const math::Vector QUERY(0.5f, -0.25f, 1.0f);

  int best_i = -1;
  float best_d2 = FLT_MAX;
  for (int i = 0; i < N; ++i) {
    float d2 = math::distance_squared(pts[i], QUERY);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_i = i;
    }
  }

  auto r = tree.nearest(QUERY, 1);
  HS_EXPECT_SIZE_OR_RETURN(r, 1);
  HS_EXPECT_EQ((int)r[0].original_index, best_i);
  HS_EXPECT_VEC(r[0].point, pts[best_i], 1e-6f);
}

/**
 * @brief Verifies nearest() handles coincident (duplicate) points and a full
 *        k == MAX_K request.
 * @details nearest() ranks by (squared distance, source index), so ties at the
 *          k boundary are pinned index by index against a brute-force scan.
 */
inline void test_kdtree_duplicates_and_max_k() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  math::Vector pts[8] = {
      math::Vector(1, 1, 1), // 0  coincident cluster (d²=0 from query)
      math::Vector(1, 1, 1), // 1  coincident
      math::Vector(1, 1, 1), // 2  coincident
      math::Vector(2, 0, 0), // 3  d²=3
      math::Vector(0, 2,
                   0), // 4  d²=3  (boundary tie: only two of 3/4/5 fit k=5)
      math::Vector(0, 0, 2),    // 5  d²=3
      math::Vector(-1, -1, -1), // 6  d²=12
      math::Vector(5, 5, 5),    // 7  d²=48
  };
  std::span<math::Vector> sp(pts, 8);
  KDTree tree(arena, sp);

  constexpr int K = KDTree::MAX_K;
  const math::Vector QUERY(1, 1, 1); // lands on the coincident cluster

  auto r = tree.nearest(QUERY, K);
  HS_EXPECT_SIZE_OR_RETURN(r, K);

  float all_d2[8];
  for (int i = 0; i < 8; ++i)
    all_d2[i] = math::distance_squared(pts[i], QUERY);
  std::sort(all_d2, all_d2 + 8);

  // (distance, source index) is a total order, so index 5 loses the d²=3
  // boundary tie to the lower-indexed 3 and 4.
  const int EXPECTED_INDEX[K] = {0, 1, 2, 3, 4};
  float prev = -1.0f;
  for (int i = 0; i < K; ++i) {
    HS_EXPECT_EQ((int)r[i].original_index, EXPECTED_INDEX[i]);
    HS_EXPECT_GE(r[i].d_sq, prev - 1e-6f);
    prev = r[i].d_sq;
    HS_EXPECT_NEAR(r[i].d_sq, all_d2[i], 1e-5f);
    const size_t oi = r[i].original_index;
    HS_EXPECT_LT(oi, (size_t)8);
    if (oi < 8)
      HS_EXPECT_VEC(r[i].point, pts[oi], 1e-6f);
    HS_EXPECT_NEAR(math::distance_squared(r[i].point, QUERY), r[i].d_sq, 1e-5f);
  }

  HS_EXPECT_NEAR(r[0].d_sq, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(r[1].d_sq, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(r[2].d_sq, 0.0f, 1e-6f);
}

/**
 * @brief Verifies k>1 nearest matches brute force on 100 random distinct points.
 */
inline void test_kdtree_k_nearest_brute_force_random() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  constexpr int N = 100;
  math::Vector pts[N];

  hs::Pcg32 rng(20240611u);
  for (int i = 0; i < N; ++i) {
    const float x = rand_uniform(rng, -10.0f, 10.0f);
    const float y = rand_uniform(rng, -10.0f, 10.0f);
    const float z = rand_uniform(rng, -10.0f, 10.0f);
    pts[i] = math::Vector(x, y, z);
  }

  std::span<math::Vector> sp(pts, N);
  KDTree tree(arena, sp);

  const math::Vector QUERIES[] = {
      math::Vector(0, 0, 0),
      math::Vector(3.5f, -7.0f, 2.0f),
      math::Vector(-8.0f, 8.0f, -1.5f),
      math::Vector(9.9f, 9.9f, 9.9f),
      math::Vector(-2.2f, 0.3f, 5.1f),
  };
  constexpr int K = KDTree::MAX_K;

  for (const math::Vector &q : QUERIES) {
    auto r = tree.nearest(q, K);
    HS_EXPECT_SIZE_OR_RETURN(r, K);

    float all_d2[N];
    for (int i = 0; i < N; ++i)
      all_d2[i] = math::distance_squared(pts[i], q);
    std::sort(all_d2, all_d2 + N);

    float prev = -1.0f;
    for (int i = 0; i < K; ++i) {
      HS_EXPECT_GE(r[i].d_sq, prev - 1e-5f); // sorted nearest-first
      prev = r[i].d_sq;
      HS_EXPECT_NEAR(r[i].d_sq, all_d2[i], 1e-4f);
      const size_t oi = r[i].original_index;
      HS_EXPECT_LT(oi, (size_t)N);
      if (oi < (size_t)N)
        HS_EXPECT_VEC(r[i].point, pts[oi], 1e-6f);
      HS_EXPECT_NEAR(math::distance_squared(r[i].point, q), r[i].d_sq, 1e-4f);
    }
  }
}

/**
 * @brief Verifies nearest_moved() matches brute force over displaced points.
 */
inline void test_kdtree_nearest_moved_matches_brute_force() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  constexpr int N = 100;
  constexpr int K = KDTree::MAX_K;
  math::Vector built[N];
  math::Vector offsets[N];
  hs::Pcg32 rng(20261018u);
  for (int i = 0; i < N; ++i) {
    const float x = rand_uniform(rng, -10.0f, 10.0f);
    const float y = rand_uniform(rng, -10.0f, 10.0f);
    const float z = rand_uniform(rng, -10.0f, 10.0f);
    built[i] = math::Vector(x, y, z);
    const float ox = rand_uniform(rng, -1.0f, 1.0f);
    const float oy = rand_uniform(rng, -1.0f, 1.0f);
    const float oz = rand_uniform(rng, -1.0f, 1.0f);
    offsets[i] = math::Vector(ox, oy, oz).normalized();
  }
  KDTree tree(arena, std::span<math::Vector>(built, N));
  static constexpr float DRIFTS[] = {0.0f, 0.5f, 4.0f};
  for (size_t scale = 0; scale < std::size(DRIFTS); ++scale) {
    const float drift = DRIFTS[scale];
    math::Vector moved[N];
    for (int i = 0; i < N; ++i)
      moved[i] = built[i] + offsets[i] * (drift * (i % 3) / 2.0f);
    for (int query = 0; query < 200; ++query) {
      const float x = rand_uniform(rng, -11.0f, 11.0f);
      const float y = rand_uniform(rng, -11.0f, 11.0f);
      const float z = rand_uniform(rng, -11.0f, 11.0f);
      const math::Vector target(x, y, z);
      const auto result = tree.nearest_moved(target, K, drift, [&](uint16_t i) {
        return math::distance_squared(moved[i], target);
      });
      HS_EXPECT_SIZE_OR_RETURN(result, K);
      int order[N];
      for (int i = 0; i < N; ++i)
        order[i] = i;
      std::sort(order, order + N, [&](int a, int b) {
        const float da = math::distance_squared(moved[a], target);
        const float db = math::distance_squared(moved[b], target);
        return da != db ? da < db : a < b;
      });
      HS_CONTEXT("drift", static_cast<long long>(scale), query);
      for (int rank = 0; rank < K; ++rank) {
        HS_EXPECT_EQ(static_cast<int>(result[rank].original_index),
                     order[rank]);
        HS_EXPECT_EQ(result[rank].d_sq,
                     math::distance_squared(moved[order[rank]], target));
        HS_EXPECT_VEC(result[rank].point, built[order[rank]], 0.0f);
      }
    }
  }
}

/**
 * @brief Verifies bounded searches return the unbounded result when the bound
 *        is the distance of k known points.
 */
inline void test_kdtree_bounded_search_matches_unbounded() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  constexpr int N = 100;
  constexpr int K = 2;
  constexpr float DRIFT = 0.5f;
  math::Vector built[N];
  math::Vector moved[N];
  hs::Pcg32 rng(20261019u);
  for (int i = 0; i < N; ++i) {
    const float x = rand_uniform(rng, -10.0f, 10.0f);
    const float y = rand_uniform(rng, -10.0f, 10.0f);
    const float z = rand_uniform(rng, -10.0f, 10.0f);
    built[i] = math::Vector(x, y, z);
    const float ox = rand_uniform(rng, -1.0f, 1.0f);
    const float oy = rand_uniform(rng, -1.0f, 1.0f);
    const float oz = rand_uniform(rng, -1.0f, 1.0f);
    moved[i] = built[i] + math::Vector(ox, oy, oz).normalized() * DRIFT;
  }
  KDTree tree(arena, std::span<math::Vector>(built, N));
  const auto same = [](const KDTree::Neighbors &a, const KDTree::Neighbors &b) {
    HS_EXPECT_SIZE_OR_RETURN(a, b.size());
    for (size_t rank = 0; rank < a.size(); ++rank) {
      HS_EXPECT_EQ(a[rank].original_index, b[rank].original_index);
      HS_EXPECT_EQ(a[rank].d_sq, b[rank].d_sq);
    }
  };
  for (int query = 0; query < 200; ++query) {
    HS_CONTEXT("query", query);
    const float x = rand_uniform(rng, -11.0f, 11.0f);
    const float y = rand_uniform(rng, -11.0f, 11.0f);
    const float z = rand_uniform(rng, -11.0f, 11.0f);
    const math::Vector target(x, y, z);
    // Two arbitrary points bound the search; the nearest pair bounds it
    // exactly.
    const int a = query % N;
    const int b = (query * 7 + 1) % N;
    const auto unbounded = tree.nearest(target, K);
    HS_EXPECT_SIZE_OR_RETURN(unbounded, K);
    same(tree.nearest(target, K,
                      std::max(math::distance_squared(built[a], target),
                               math::distance_squared(built[b], target))),
         unbounded);
    same(tree.nearest(target, K, unbounded[K - 1].d_sq), unbounded);

    const auto moved_distance = [&](uint16_t i) {
      return math::distance_squared(moved[i], target);
    };
    const auto moved_unbounded =
        tree.nearest_moved(target, K, DRIFT, moved_distance);
    HS_EXPECT_SIZE_OR_RETURN(moved_unbounded, K);
    same(tree.nearest_moved(target, K, DRIFT, moved_distance,
                            std::max(moved_distance(a), moved_distance(b))),
         moved_unbounded);
    same(tree.nearest_moved(target, K, DRIFT, moved_distance,
                            moved_unbounded[K - 1].d_sq),
         moved_unbounded);
  }
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every spatial test case under the "spatial" module scope.
 * @return The module's failure count.
 */
inline int run_spatial_tests() {
  hs_test::ModuleFixture fixture("spatial");

  HS_RUN(test_kdtree_empty_input);
  HS_RUN(test_kdtree_single_point);
  HS_RUN(test_kdtree_move);
  HS_RUN(test_kdtree_nearest_known_set);
  HS_RUN(test_kdtree_k_nearest_sorted);
  HS_RUN(test_kdtree_ties_prefer_source_index);
  HS_RUN(test_kdtree_k_caps_at_size);
  HS_RUN(test_kdtree_default_unbuilt);
  HS_RUN(test_kdtree_matches_brute_force);
  HS_RUN(test_kdtree_duplicates_and_max_k);
  HS_RUN(test_kdtree_k_nearest_brute_force_random);
  HS_RUN(test_kdtree_nearest_moved_matches_brute_force);
  HS_RUN(test_kdtree_bounded_search_matches_unbounded);

  return fixture.result();
}

} // namespace spatial_tests
} // namespace hs_test
