/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/spatial.h — AABB, KDTree, MeshState.
 *
 * Tests deliberately avoid invoking the asserts in dependent types
 * (out-of-bounds, unbound access). Zero-component (axis-aligned) ray
 * directions for AABB::intersect_ray are exercised here, including the
 * grazing on-face case where the slab test must guard parallel rays
 * explicitly to avoid a 0/0 NaN.
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

// Dedicated arena buffer to keep KDTree/MeshState scratch isolated from the
// memory-test arenas.
inline uint8_t spatial_buf[128 * 1024];

// ============================================================================
// AABB
// ============================================================================

// A default-constructed AABB is the degenerate "empty" box, ready to be
// widened by the first expand(): min = +FLT_MAX, max = -FLT_MAX.
inline void test_aabb_default_empty() {
  AABB box;
  HS_EXPECT_TRUE(box.min_val.x >= FLT_MAX * 0.5f);
  HS_EXPECT_TRUE(box.max_val.x <= -FLT_MAX * 0.5f);
}

// Expanding an empty box by one point collapses it to a zero-volume box at
// that point (min == max).
inline void test_aabb_expand_single_point() {
  AABB box;
  box.expand(Vector(1, 2, 3));
  HS_EXPECT_VEC(box.min_val, Vector(1, 2, 3), 1e-6f);
  HS_EXPECT_VEC(box.max_val, Vector(1, 2, 3), 1e-6f);
}

// expand() takes the component-wise min/max across all added points, so the
// final box is the tight bound of the point set.
inline void test_aabb_expand_multiple_points() {
  AABB box;
  box.expand(Vector(1, 2, 3));
  box.expand(Vector(-1, 5, 2));
  box.expand(Vector(4, 0, -6));
  HS_EXPECT_VEC(box.min_val, Vector(-1, 0, -6), 1e-6f);
  HS_EXPECT_VEC(box.max_val, Vector(4, 5, 3), 1e-6f);
}

// union_with() merges another box into this one, yielding the bound that
// encloses both.
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

// union_with() never shrinks a box: merging an empty box or a fully-contained
// box leaves it unchanged, while merging a superset grows it.
inline void test_aabb_union_with_empty_and_subset_is_noop() {
  AABB a;
  a.expand(Vector(0, 0, 0));
  a.expand(Vector(1, 1, 1));
  AABB before = a;

  // Union with a default-constructed (empty) box is a true no-op: an empty box
  // is min=+FLT_MAX, max=-FLT_MAX, so no component of `a` can be widened.
  AABB empty;
  a.union_with(empty);
  HS_EXPECT_VEC(a.min_val, before.min_val, 1e-6f);
  HS_EXPECT_VEC(a.max_val, before.max_val, 1e-6f);

  // Union with a more-permissive box should grow a to the superset.
  AABB superset;
  superset.expand(Vector(-10, -10, -10));
  superset.expand(Vector(10, 10, 10));
  a.union_with(superset);
  HS_EXPECT_VEC(a.min_val, Vector(-10, -10, -10), 1e-6f);
  HS_EXPECT_VEC(a.max_val, Vector(10, 10, 10), 1e-6f);

  // Union with a contained box should not shrink a.
  AABB inside;
  inside.expand(Vector(-5, -5, -5));
  inside.expand(Vector(5, 5, 5));
  a.union_with(inside);
  HS_EXPECT_VEC(a.min_val, Vector(-10, -10, -10), 1e-6f);
  HS_EXPECT_VEC(a.max_val, Vector(10, 10, 10), 1e-6f);
}

// Rays aimed at the box from outside report a hit: axis-aligned and diagonal
// directions all pass the slab test.
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

// Rays that pass beside the box or point away from it must report a miss,
// including the "behind and receding" case where the entry t is negative.
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

// A ray whose origin is inside the box always hits, since the slab test
// admits a negative entry t.
inline void test_aabb_ray_from_inside() {
  AABB box;
  box.expand(Vector(-1, -1, -1));
  box.expand(Vector(1, 1, 1));
  HS_EXPECT_TRUE(box.intersect_ray(Vector(0, 0, 0), Vector(1, 0, 0)));
  HS_EXPECT_TRUE(box.intersect_ray(Vector(0, 0, 0), Vector(0, 1, 0)));
}

// Rays with a zero direction component run parallel to an axis slab. The slab
// guard must treat them as "hit only if the origin is within the slab",
// avoiding the 0/0 NaN of dividing the slab distance by a zero direction.
inline void test_aabb_ray_parallel_grazing() {
  AABB box;
  box.expand(Vector(-1, -1, -1));
  box.expand(Vector(1, 1, 1));
  // Origin lies exactly ON the +Y face, ray parallel to it: the zero Y
  // direction component means (max_val.y - origin.y) / 0 would be 0/0, so the
  // guard counts this on-face grazing ray as touching the box.
  HS_EXPECT_TRUE(box.intersect_ray(Vector(0, 1, 0), Vector(1, 0, 0)));
  // Same parallel direction but origin just outside the +Y face → clean miss.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(0, 1.5f, 0), Vector(1, 0, 0)));
  // Parallel to two axes, origin outside the Z slab → miss.
  HS_EXPECT_FALSE(box.intersect_ray(Vector(0, 0, 2), Vector(1, 0, 0)));
}

// ============================================================================
// KDTree
// ============================================================================

// Building from a zero-length span leaves the tree unbuilt (root_index == -1)
// and queries return nothing.
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

// A one-point tree builds a root and returns that point (with its original
// index) for any query.
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

// nearest(k=1) returns the correct point and original index for queries near
// distinct points, including a query landing exactly on a point (dist 0).
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

// A k>1 query returns exactly k neighbors ordered nearest-first.
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

// Requesting more neighbors than exist returns all available points, not k.
inline void test_kdtree_k_caps_at_size() {
  Arena arena(spatial_buf, sizeof(spatial_buf));
  Vector pts[3] = {Vector(0, 0, 0), Vector(1, 0, 0), Vector(2, 0, 0)};
  std::span<Vector> sp(pts, 3);
  KDTree tree(arena, sp);

  auto r = tree.nearest(Vector(0, 0, 0), 5);
  HS_EXPECT_EQ(r.size(), (size_t)3);
}

// A default-constructed (never-built) KDTree reports root_index == -1 and
// answers queries with an empty result.
inline void test_kdtree_default_unbuilt() {
  KDTree tree;
  HS_EXPECT_EQ(tree.root_index, -1);
  auto r = tree.nearest(Vector(1, 2, 3), 1);
  HS_EXPECT_TRUE(r.is_empty());
}

// clear() detaches the root so the tree is logically empty again; node storage
// is reclaimed on the next build via nodes.bind.
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

// A default-constructed MeshState owns no arena storage: unbound, zero
// vertices, zero faces.
inline void test_meshstate_default_unbound() {
  MeshState m;
  HS_EXPECT_FALSE(m.is_bound());
  HS_EXPECT_EQ(m.num_vertices(), (size_t)0);
  HS_EXPECT_EQ(m.num_faces(), (size_t)0);
}

// clone() deep-copies vertices, face counts, and faces into the destination
// arena, so the copy is value-equal yet backed by separate storage.
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

// clear() empties the mesh's vertex and face views, returning the counts to
// zero.
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

// Move-constructing a MeshState transfers ownership: the destination holds the
// data and the source is left unbound and empty.
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

// When the owned face_counts vector is empty but face_counts_view is set, the
// unified accessors fall through to the external view.
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

// Runs every spatial test case under the "spatial" module scope; returns the
// module's failure count.
inline int run_spatial_tests() {
  auto scope = hs_test::begin_module("spatial");

  test_aabb_default_empty();
  test_aabb_expand_single_point();
  test_aabb_expand_multiple_points();
  test_aabb_union_with();
  test_aabb_union_with_empty_and_subset_is_noop();
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

