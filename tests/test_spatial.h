/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/spatial.h — KDTree, MeshState.
 *
 * Tests deliberately avoid invoking the asserts in dependent types
 * (out-of-bounds, unbound access).
 */
#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <random>
#include "core/spatial.h"
#include "tests/test_3dmath.h" // re-uses approx_vec / HS_EXPECT_VEC
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace spatial {

using hs_test::math3d::approx_vec;

// Module-scope scratch buffer; each test re-bases the bump pointer by
// constructing a fresh Arena over it at entry. Do NOT retain an ArenaVector or
// pointer into this buffer past its own test scope.
inline constexpr size_t kSpatialBufBytes = 128 * 1024;
inline uint8_t spatial_buf[kSpatialBufBytes];

// Split offset for the one test needing two disjoint arenas over this buffer.
inline constexpr size_t kSpatialBufSplit = kSpatialBufBytes / 2;

// ============================================================================
// KDTree
// ============================================================================

/**
 * @brief Verifies building from a zero-length span leaves the tree unbuilt
 *        (root_index == -1) and queries return nothing.
 */
inline void test_kdtree_empty_input() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[1] = {Vector(0, 0, 0)};
  std::span<Vector> empty(pts, 0);
  KDTree tree(arena, empty);
  HS_EXPECT_EQ(tree.root_index, -1);

  auto r = tree.nearest(Vector(0, 0, 0), 1);
  HS_EXPECT_TRUE(r.is_empty());
}

/**
 * @brief Verifies a one-point tree builds a root and returns that point (with
 *        its original index) for any query.
 */
inline void test_kdtree_single_point() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[1] = {Vector(3, 4, 5)};
  std::span<Vector> sp(pts, 1);
  KDTree tree(arena, sp);
  HS_EXPECT_NE(tree.root_index, -1);

  auto r = tree.nearest(Vector(0, 0, 0), 1);
  HS_EXPECT_EQ(r.size(), (size_t)1);
  HS_EXPECT_VEC(r[0].point, Vector(3, 4, 5), 1e-6f);
  HS_EXPECT_EQ(r[0].original_index, (uint16_t)0);
}

/**
 * @brief Verifies nearest(k=1) returns the correct point and original index for
 *        queries near distinct points.
 * @details Includes a query landing exactly on a point (dist 0).
 */
inline void test_kdtree_nearest_known_set() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[6] = {
      Vector(0, 0, 0),    // 0
      Vector(10, 0, 0),   // 1
      Vector(0, 10, 0),   // 2
      Vector(0, 0, 10),   // 3
      Vector(5, 5, 5),    // 4
      Vector(-3, 2, 1),   // 5
  };
  std::span<Vector> sp(pts, 6);
  KDTree tree(arena, sp);

  auto r1 = tree.nearest(Vector(0.1f, 0.1f, 0.1f), 1);
  HS_EXPECT_EQ(r1.size(), (size_t)1);
  HS_EXPECT_VEC(r1[0].point, Vector(0, 0, 0), 1e-6f);
  HS_EXPECT_EQ(r1[0].original_index, (uint16_t)0);

  auto r2 = tree.nearest(Vector(10.1f, 0, 0), 1);
  HS_EXPECT_EQ(r2.size(), (size_t)1);
  HS_EXPECT_VEC(r2[0].point, Vector(10, 0, 0), 1e-6f);
  HS_EXPECT_EQ(r2[0].original_index, (uint16_t)1);

  auto r3 = tree.nearest(Vector(5, 5, 5), 1);
  HS_EXPECT_EQ(r3.size(), (size_t)1);
  HS_EXPECT_VEC(r3[0].point, Vector(5, 5, 5), 1e-6f);
  HS_EXPECT_EQ(r3[0].original_index, (uint16_t)4);
}

/**
 * @brief Verifies a k>1 query returns exactly k neighbors ordered
 *        nearest-first.
 */
inline void test_kdtree_k_nearest_sorted() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[5] = {
      Vector(0, 0, 0),  // 0  d²=0
      Vector(1, 0, 0),  // 1  d²=1
      Vector(2, 0, 0),  // 2  d²=4
      Vector(3, 0, 0),  // 3  d²=9
      Vector(4, 0, 0),  // 4  d²=16
  };
  std::span<Vector> sp(pts, 5);
  KDTree tree(arena, sp);

  auto r = tree.nearest(Vector(0, 0, 0), 3);
  HS_EXPECT_EQ(r.size(), (size_t)3);
  HS_EXPECT_VEC(r[0].point, Vector(0, 0, 0), 1e-6f);
  HS_EXPECT_VEC(r[1].point, Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(r[2].point, Vector(2, 0, 0), 1e-6f);
}

/**
 * @brief Verifies requesting more neighbors than exist returns all available
 *        points, not k.
 */
inline void test_kdtree_k_caps_at_size() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[3] = {Vector(0, 0, 0), Vector(1, 0, 0), Vector(2, 0, 0)};
  std::span<Vector> sp(pts, 3);
  KDTree tree(arena, sp);

  auto r = tree.nearest(Vector(0, 0, 0), 5);
  HS_EXPECT_EQ(r.size(), (size_t)3);
}

/**
 * @brief Verifies a default-constructed (never-built) KDTree reports
 *        root_index == -1 and answers queries with an empty result.
 */
inline void test_kdtree_default_unbuilt() {
  KDTree tree;
  HS_EXPECT_EQ(tree.root_index, -1);
  auto r = tree.nearest(Vector(1, 2, 3), 1);
  HS_EXPECT_TRUE(r.is_empty());
}

/**
 * @brief Verifies clear() detaches the root so the tree is logically empty
 *        again.
 * @details Node storage is reclaimed on the next build via nodes.bind.
 */
inline void test_kdtree_clear() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[2] = {Vector(0, 0, 0), Vector(1, 1, 1)};
  std::span<Vector> sp(pts, 2);
  KDTree tree(arena, sp);
  HS_EXPECT_NE(tree.root_index, -1);

  tree.clear();
  HS_EXPECT_EQ(tree.root_index, -1);
  auto r = tree.nearest(Vector(0, 0, 0), 1);
  HS_EXPECT_TRUE(r.is_empty());
}

/**
 * @brief Verifies the KDTree nearest neighbor matches a brute-force scan.
 * @details With 16 deterministic points, the nearest neighbor must match a
 *          manual scan over all distances.
 */
inline void test_kdtree_matches_brute_force() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  constexpr int N = 16;
  Vector pts[N];
  for (int i = 0; i < N; ++i) {
    float fi = static_cast<float>(i);
    pts[i] = Vector(std::sin(fi) * 5.0f, std::cos(fi * 1.3f) * 4.0f,
                    std::sin(fi * 0.7f) * 3.0f);
  }

  std::span<Vector> sp(pts, N);
  KDTree tree(arena, sp);

  const Vector query(0.5f, -0.25f, 1.0f);

  int best_i = -1;
  float best_d2 = FLT_MAX;
  for (int i = 0; i < N; ++i) {
    float d2 = distance_squared(pts[i], query);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_i = i;
    }
  }

  auto r = tree.nearest(query, 1);
  HS_EXPECT_EQ(r.size(), (size_t)1);
  HS_EXPECT_EQ((int)r[0].original_index, best_i);
  HS_EXPECT_VEC(r[0].point, pts[best_i], 1e-6f);
}

/**
 * @brief Verifies nearest() handles coincident (duplicate) points and a full
 *        k == MAX_K request, the classic degenerate cases.
 * @details Three points share a location, so a k-nearest query produces distance
 *          ties — including at the k boundary, where only some of the equidistant
 *          points fit. Tie order is unspecified, so the result is checked as a
 *          sorted multiset of squared distances against a brute-force k-smallest
 *          scan, and every returned neighbor is cross-checked against its source
 *          point and recomputed distance.
 */
inline void test_kdtree_duplicates_and_max_k() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[8] = {
      Vector(1, 1, 1),     // 0  coincident cluster (d²=0 from query)
      Vector(1, 1, 1),     // 1  coincident
      Vector(1, 1, 1),     // 2  coincident
      Vector(2, 0, 0),     // 3  d²=3
      Vector(0, 2, 0),     // 4  d²=3  (boundary tie: only one of 3/4/5 fits k=5)
      Vector(0, 0, 2),     // 5  d²=3
      Vector(-1, -1, -1),  // 6  d²=12
      Vector(5, 5, 5),     // 7  d²=48
  };
  std::span<Vector> sp(pts, 8);
  KDTree tree(arena, sp);

  constexpr int K = KDTree::MAX_K; // 5
  const Vector query(1, 1, 1);     // lands on the coincident cluster

  auto r = tree.nearest(query, K);
  HS_EXPECT_EQ(r.size(), (size_t)K);

  float all_d2[8];
  for (int i = 0; i < 8; ++i) all_d2[i] = distance_squared(pts[i], query);
  std::sort(all_d2, all_d2 + 8);

  // Compare as a sorted sequence: boundary ties make index order arbitrary.
  float prev = -1.0f;
  for (int i = 0; i < K; ++i) {
    HS_EXPECT_TRUE(r[i].d_sq >= prev - 1e-6f);
    prev = r[i].d_sq;
    HS_EXPECT_TRUE(std::fabs(r[i].d_sq - all_d2[i]) < 1e-5f);
    HS_EXPECT_VEC(r[i].point, pts[r[i].original_index], 1e-6f);
    HS_EXPECT_TRUE(std::fabs(distance_squared(r[i].point, query) - r[i].d_sq) <
                   1e-5f);
  }

  HS_EXPECT_TRUE(std::fabs(r[0].d_sq) < 1e-6f);
  HS_EXPECT_TRUE(std::fabs(r[1].d_sq) < 1e-6f);
  HS_EXPECT_TRUE(std::fabs(r[2].d_sq) < 1e-6f);
}

/**
 * @brief Verifies k>1 nearest matches brute force on 100 random distinct points.
 * @details The k>1 brute-force checks elsewhere use coincident/collinear points
 *          where bbox pruning can't err. This builds a tree over 100 randomly
 *          placed distinct points (fixed mt19937 seed for reproducibility) and,
 *          for several queries, compares the tree's k-nearest set against a
 *          brute-force k-smallest scan. Distances are well separated, so the
 *          comparison is per-rank on sorted squared distance plus a point/index
 *          cross-check.
 */
inline void test_kdtree_k_nearest_brute_force_random() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  constexpr int N = 100;
  Vector pts[N];

  std::mt19937 rng(20240611u);
  std::uniform_real_distribution<float> uni(-10.0f, 10.0f);
  for (int i = 0; i < N; ++i)
    pts[i] = Vector(uni(rng), uni(rng), uni(rng));

  std::span<Vector> sp(pts, N);
  KDTree tree(arena, sp);

  const Vector queries[] = {
      Vector(0, 0, 0),     Vector(3.5f, -7.0f, 2.0f),
      Vector(-8.0f, 8.0f, -1.5f), Vector(9.9f, 9.9f, 9.9f),
      Vector(-2.2f, 0.3f, 5.1f),
  };
  constexpr int K = KDTree::MAX_K;

  for (const Vector &q : queries) {
    auto r = tree.nearest(q, K);
    HS_EXPECT_EQ(r.size(), (size_t)K);

    float all_d2[N];
    for (int i = 0; i < N; ++i) all_d2[i] = distance_squared(pts[i], q);
    std::sort(all_d2, all_d2 + N);

    float prev = -1.0f;
    for (int i = 0; i < K; ++i) {
      HS_EXPECT_TRUE(r[i].d_sq >= prev - 1e-5f); // sorted nearest-first
      prev = r[i].d_sq;
      HS_EXPECT_TRUE(std::fabs(r[i].d_sq - all_d2[i]) < 1e-4f);
      HS_EXPECT_VEC(r[i].point, pts[r[i].original_index], 1e-6f);
      HS_EXPECT_TRUE(
          std::fabs(distance_squared(r[i].point, q) - r[i].d_sq) < 1e-4f);
    }
  }
}

// ============================================================================
// MeshState
// ============================================================================

/**
 * @brief Verifies a default-constructed MeshState owns no arena storage.
 * @details It is unbound with zero vertices and zero faces.
 */
inline void test_meshstate_default_unbound() {
  MeshState m;
  HS_EXPECT_FALSE(m.is_bound());
  HS_EXPECT_EQ(m.num_vertices(), (size_t)0);
  HS_EXPECT_EQ(m.num_faces(), (size_t)0);
}

/**
 * @brief Verifies clone() deep-copies vertices, face counts, and faces into the
 *        destination arena.
 * @details The copy is value-equal yet backed by separate storage.
 */
inline void test_meshstate_clone_deep_copies() {
  Arena src_arena(spatial_buf, kSpatialBufSplit);
  Arena dst_arena(spatial_buf + kSpatialBufSplit,
                  kSpatialBufBytes - kSpatialBufSplit);

  MeshState src;
  src.vertices.bind(src_arena, 3);
  src.vertices.push_back(Vector(1, 0, 0));
  src.vertices.push_back(Vector(0, 1, 0));
  src.vertices.push_back(Vector(0, 0, 1));

  src.face_counts.bind(src_arena, 1);
  src.face_counts.push_back(3);

  src.faces.bind(src_arena, 3);
  src.faces.push_back(0);
  src.faces.push_back(1);
  src.faces.push_back(2);

  MeshState dst;
  MeshState::clone(src, dst, dst_arena);

  HS_EXPECT_EQ(dst.num_vertices(), (size_t)3);
  HS_EXPECT_VEC(dst.vertices[0], Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(dst.vertices[2], Vector(0, 0, 1), 1e-6f);

  HS_EXPECT_EQ(dst.num_faces(), (size_t)1);
  HS_EXPECT_EQ(dst.face_counts[0], (uint8_t)3);
  HS_EXPECT_EQ(dst.get_faces_size(), (size_t)3);
  HS_EXPECT_EQ(dst.faces[0], (uint16_t)0);
  HS_EXPECT_EQ(dst.faces[1], (uint16_t)1);
  HS_EXPECT_EQ(dst.faces[2], (uint16_t)2);

  HS_EXPECT_TRUE(dst.vertices.data() != src.vertices.data());
}

/**
 * @brief Verifies clear() empties the mesh's vertex and face views, returning
 *        the counts to zero.
 */
inline void test_meshstate_clear_resets_views() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  MeshState m;
  m.vertices.bind(arena, 1);
  m.vertices.push_back(Vector(1, 1, 1));
  m.face_counts.bind(arena, 1);
  m.face_counts.push_back(1);
  HS_EXPECT_EQ(m.num_vertices(), (size_t)1);
  HS_EXPECT_EQ(m.num_faces(), (size_t)1);

  m.clear();
  HS_EXPECT_EQ(m.num_vertices(), (size_t)0);
  HS_EXPECT_EQ(m.num_faces(), (size_t)0);
}

/**
 * @brief Verifies move-constructing a MeshState transfers ownership.
 * @details The destination holds the data and the source is left unbound and
 *          empty.
 */
inline void test_meshstate_move_invalidates_source() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  MeshState src;
  src.vertices.bind(arena, 2);
  src.vertices.push_back(Vector(1, 2, 3));
  src.vertices.push_back(Vector(4, 5, 6));

  MeshState dst(std::move(src));
  HS_EXPECT_FALSE(src.vertices.is_bound());
  HS_EXPECT_EQ(src.num_vertices(), (size_t)0);
  HS_EXPECT_EQ(dst.num_vertices(), (size_t)2);
  HS_EXPECT_VEC(dst.vertices[1], Vector(4, 5, 6), 1e-6f);
}

/**
 * @brief Verifies the unified accessors fall through to the external view.
 * @details When the owned face_counts vector is empty but face_counts_view is
 *          set, accessors read from the external view.
 */
inline void test_meshstate_view_fallback() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  ArenaVector<uint8_t> owner(arena, 2);
  owner.push_back(3);
  owner.push_back(4);

  MeshState m;
  m.face_counts_view = ArenaSpan<uint8_t>(owner);

  HS_EXPECT_EQ(m.get_face_counts_size(), (size_t)2);
  HS_EXPECT_EQ(m.get_face_counts_data()[0], (uint8_t)3);
  HS_EXPECT_EQ(m.get_face_counts_data()[1], (uint8_t)4);
  HS_EXPECT_EQ(m.num_faces(), (size_t)2);
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

  test_kdtree_empty_input();
  test_kdtree_single_point();
  test_kdtree_nearest_known_set();
  test_kdtree_k_nearest_sorted();
  test_kdtree_k_caps_at_size();
  test_kdtree_default_unbuilt();
  test_kdtree_clear();
  test_kdtree_matches_brute_force();
  test_kdtree_duplicates_and_max_k();
  test_kdtree_k_nearest_brute_force_random();

  test_meshstate_default_unbound();
  test_meshstate_clone_deep_copies();
  test_meshstate_clear_resets_views();
  test_meshstate_move_invalidates_source();
  test_meshstate_view_fallback();

  return fixture.result();
}

} // namespace spatial
} // namespace hs_test

