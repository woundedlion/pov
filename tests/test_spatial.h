/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/spatial.h — AABB, KDTree, MeshState.
 *
 * Tests deliberately avoid invoking the asserts in dependent types
 * (out-of-bounds, unbound access). Zero-component (axis-aligned) ray
 * directions for AABB::intersect_ray ARE exercised here — including the
 * grazing on-face case that previously produced a 0/0 NaN — now that the
 * slab test guards parallel rays explicitly.
 */
#pragma once

#include <cfloat>
#include <cstdint>
#include "core/spatial.h"
#include "tests/test_3dmath.h" // re-uses approx_vec / HS_EXPECT_VEC
#include "tests/test_harness.h"

namespace hs_test {
namespace spatial {

using hs_test::math3d::approx_vec;

// Reuse the math3d HS_EXPECT_VEC macro (defined in test_3dmath.h)

// Dedicated arena buffer to keep KDTree/MeshState scratch isolated from the
// memory-test arenas.
inline uint8_t spatial_buf[128 * 1024];

// ============================================================================
// AABB
// ============================================================================

inline void test_aabb_default_empty() {
  AABB box;
  // Empty box: min initialized to +FLT_MAX, max to -FLT_MAX (degenerate).
  HS_EXPECT_TRUE(box.min_val.x >= FLT_MAX * 0.5f);
  HS_EXPECT_TRUE(box.max_val.x <= -FLT_MAX * 0.5f);
}

inline void test_aabb_expand_single_point() {
  AABB box;
  box.expand(Vector(1, 2, 3));
  HS_EXPECT_VEC(box.min_val, Vector(1, 2, 3), 1e-6f);
  HS_EXPECT_VEC(box.max_val, Vector(1, 2, 3), 1e-6f);
}

inline void test_aabb_expand_multiple_points() {
  AABB box;
  box.expand(Vector(1, 2, 3));
  box.expand(Vector(-1, 5, 2));
  box.expand(Vector(4, 0, -6));
  HS_EXPECT_VEC(box.min_val, Vector(-1, 0, -6), 1e-6f);
  HS_EXPECT_VEC(box.max_val, Vector(4, 5, 3), 1e-6f);
}

inline void test_aabb_union_with() {
  AABB a;
  a.expand(Vector(0, 0, 0));
  a.expand(Vector(2, 2, 2));

  AABB b;
  b.expand(Vector(-1, 1, 5));
  b.expand(Vector(3, 4, 6));

  a.union_with(b);
  HS_EXPECT_VEC(a.min_val, Vector(-1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(a.max_val, Vector(3, 4, 6), 1e-6f);
}

inline void test_aabb_union_with_empty_is_noop_when_outside() {
  AABB a;
  a.expand(Vector(0, 0, 0));
  a.expand(Vector(1, 1, 1));
  AABB before = a;

  // Union with a more-permissive box should not shrink a.
  AABB superset;
  superset.expand(Vector(-10, -10, -10));
  superset.expand(Vector(10, 10, 10));
  a.union_with(superset);
  HS_EXPECT_VEC(a.min_val, Vector(-10, -10, -10), 1e-6f);
  HS_EXPECT_VEC(a.max_val, Vector(10, 10, 10), 1e-6f);

  // Union with a contained box should not change a.
  AABB inside;
  inside.expand(Vector(-5, -5, -5));
  inside.expand(Vector(5, 5, 5));
  a.union_with(inside);
  HS_EXPECT_VEC(a.min_val, Vector(-10, -10, -10), 1e-6f);
  HS_EXPECT_VEC(a.max_val, Vector(10, 10, 10), 1e-6f);
  (void)before;
}

inline void test_aabb_ray_hit() {
  AABB box;
  box.expand(Vector(-1, -1, -1));
  box.expand(Vector(1, 1, 1));

  // Ray from +X axis pointing inward
  HS_EXPECT_TRUE(box.intersect_ray(Vector(5, 0, 0), Vector(-1, 0, 0)));
  // Ray from -Y axis pointing inward
  HS_EXPECT_TRUE(box.intersect_ray(Vector(0, -5, 0), Vector(0, 1, 0)));
  // Diagonal ray through the centre
  HS_EXPECT_TRUE(box.intersect_ray(Vector(5, 5, 5),
                                  Vector(-1, -1, -1).normalized()));
}

inline void test_aabb_ray_miss() {
  AABB box;
  box.expand(Vector(-1, -1, -1));
  box.expand(Vector(1, 1, 1));

  // Ray parallel to the box but offset in Y — passes over the top.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(-5, 5, 0), Vector(1, 0, 0)));
  // Ray parallel to box, offset in Z — passes alongside.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(0, 0, 5), Vector(1, 0, 0)));
  // Ray starts in front of the box and points AWAY from it.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(5, 0, 0), Vector(1, 0, 0)));
  // Ray starts behind the box and points AWAY from it.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(-5, 0, 0), Vector(-1, 0, 0)));
}

inline void test_aabb_ray_from_inside() {
  AABB box;
  box.expand(Vector(-1, -1, -1));
  box.expand(Vector(1, 1, 1));
  // A ray starting inside the box always hits — slab test allows negative t.
  HS_EXPECT_TRUE(box.intersect_ray(Vector(0, 0, 0), Vector(1, 0, 0)));
  HS_EXPECT_TRUE(box.intersect_ray(Vector(0, 0, 0), Vector(0, 1, 0)));
}

inline void test_aabb_ray_parallel_grazing() {
  AABB box;
  box.expand(Vector(-1, -1, -1));
  box.expand(Vector(1, 1, 1));
  // Origin lies exactly ON the +Y face and the ray runs parallel to it. The
  // naive (max_val.y - origin.y) / 0 == 0/0 == NaN; the slab guard treats a zero
  // direction component as "hit only if the origin is within the slab", so an
  // on-face grazing ray still counts as touching the box.
  HS_EXPECT_TRUE(box.intersect_ray(Vector(0, 1, 0), Vector(1, 0, 0)));
  // Same parallel direction but origin just outside the +Y face → clean miss.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(0, 1.5f, 0), Vector(1, 0, 0)));
  // Parallel to two axes, origin outside the Z slab → miss.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(0, 0, 2), Vector(1, 0, 0)));
}

// ============================================================================
// KDTree
// ============================================================================

inline void test_kdtree_empty_input() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[1] = {Vector(0, 0, 0)};
  std::span<Vector> empty(pts, 0);
  KDTree tree(arena, empty);
  HS_EXPECT_EQ(tree.root_index, -1);

  // Querying an empty tree returns an empty result.
  auto r = tree.nearest(Vector(0, 0, 0), 1);
  HS_EXPECT_TRUE(r.is_empty());
}

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

  // Query near (0.1, 0.1, 0.1) → closest is index 0 at the origin
  auto r1 = tree.nearest(Vector(0.1f, 0.1f, 0.1f), 1);
  HS_EXPECT_EQ(r1.size(), (size_t)1);
  HS_EXPECT_VEC(r1[0].point, Vector(0, 0, 0), 1e-6f);
  HS_EXPECT_EQ(r1[0].original_index, (uint16_t)0);

  // Query at (10.1, 0, 0) → closest is index 1
  auto r2 = tree.nearest(Vector(10.1f, 0, 0), 1);
  HS_EXPECT_EQ(r2.size(), (size_t)1);
  HS_EXPECT_VEC(r2[0].point, Vector(10, 0, 0), 1e-6f);
  HS_EXPECT_EQ(r2[0].original_index, (uint16_t)1);

  // Query exactly at one of the points → that point comes back, dist = 0
  auto r3 = tree.nearest(Vector(5, 5, 5), 1);
  HS_EXPECT_EQ(r3.size(), (size_t)1);
  HS_EXPECT_VEC(r3[0].point, Vector(5, 5, 5), 1e-6f);
  HS_EXPECT_EQ(r3[0].original_index, (uint16_t)4);
}

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
  // Results must be sorted by distance (closest first).
  HS_EXPECT_VEC(r[0].point, Vector(0, 0, 0), 1e-6f);
  HS_EXPECT_VEC(r[1].point, Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(r[2].point, Vector(2, 0, 0), 1e-6f);
}

inline void test_kdtree_k_caps_at_size() {
  // Request more neighbors than exist → returns all available.
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[3] = {Vector(0, 0, 0), Vector(1, 0, 0), Vector(2, 0, 0)};
  std::span<Vector> sp(pts, 3);
  KDTree tree(arena, sp);

  auto r = tree.nearest(Vector(0, 0, 0), 5);
  HS_EXPECT_EQ(r.size(), (size_t)3);
}

inline void test_kdtree_default_unbuilt() {
  // A default-constructed KDTree (no build) has root_index == -1.
  KDTree tree;
  HS_EXPECT_EQ(tree.root_index, -1);
  auto r = tree.nearest(Vector(1, 2, 3), 1);
  HS_EXPECT_TRUE(r.is_empty());
}

inline void test_kdtree_clear() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[2] = {Vector(0, 0, 0), Vector(1, 1, 1)};
  std::span<Vector> sp(pts, 2);
  KDTree tree(arena, sp);
  HS_EXPECT_NE(tree.root_index, -1);

  tree.clear();
  HS_EXPECT_EQ(tree.root_index, -1);
  HS_EXPECT_EQ(tree.node_count, (size_t)0);
  // After clear the tree returns nothing
  auto r = tree.nearest(Vector(0, 0, 0), 1);
  HS_EXPECT_TRUE(r.is_empty());
}

// Brute-force verification: with 16 deterministic points, the nearest
// neighbor must match a manual scan over all distances.
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

// ============================================================================
// MeshState
// ============================================================================

inline void test_meshstate_default_unbound() {
  MeshState m;
  HS_EXPECT_FALSE(m.is_bound());
  HS_EXPECT_EQ(m.num_vertices(), (size_t)0);
  HS_EXPECT_EQ(m.num_faces(), (size_t)0);
}

inline void test_meshstate_clone_deep_copies() {
  Arena src_arena(spatial_buf, sizeof(spatial_buf));
  Arena dst_arena(spatial_buf + 64 * 1024, sizeof(spatial_buf) - 64 * 1024);

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

  // Clone is a real copy — destination lives in dst_arena, not src_arena.
  HS_EXPECT_TRUE(dst.vertices.data() != src.vertices.data());
}

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
  HS_EXPECT_FALSE(m.cache_valid);
}

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

inline void test_meshstate_view_fallback() {
  // When face_counts is empty but face_counts_view is set, the unified
  // accessor must return the view.
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

inline int run_spatial_tests() {
  auto scope = hs_test::begin_module("spatial");

  test_aabb_default_empty();
  test_aabb_expand_single_point();
  test_aabb_expand_multiple_points();
  test_aabb_union_with();
  test_aabb_union_with_empty_is_noop_when_outside();
  test_aabb_ray_hit();
  test_aabb_ray_miss();
  test_aabb_ray_from_inside();
  test_aabb_ray_parallel_grazing();

  test_kdtree_empty_input();
  test_kdtree_single_point();
  test_kdtree_nearest_known_set();
  test_kdtree_k_nearest_sorted();
  test_kdtree_k_caps_at_size();
  test_kdtree_default_unbuilt();
  test_kdtree_clear();
  test_kdtree_matches_brute_force();

  test_meshstate_default_unbound();
  test_meshstate_clone_deep_copies();
  test_meshstate_clear_resets_views();
  test_meshstate_move_invalidates_source();
  test_meshstate_view_fallback();

  return hs_test::end_module(scope);
}

} // namespace spatial
} // namespace hs_test

