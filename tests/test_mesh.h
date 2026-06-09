/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/mesh.h.
 *
 * Coverage:
 *   - PolyMesh basic API (initialize, clear, accessors)
 *   - HalfEdgeMesh construction from PolyMesh and MeshState
 *     (Euler invariant V - E + F = 2, half-edge pairing, vertex/face refs)
 *   - MeshOps::compile (PolyMesh → MeshState, drops degenerate faces,
 *     populates face_offsets)
 *   - MeshOps::clone (deep copy into a target arena)
 *   - MeshOps::classify_faces_by_topology (cube → all faces same class)
 */
#pragma once

#include <cstdint>
#include "core/mesh.h"
#include "core/solids.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mesh_tests {

inline uint8_t mesh_arena_a[256 * 1024];
inline uint8_t mesh_arena_b[256 * 1024];

template <typename Solid>
inline void build_solid(PolyMesh &mesh, Arena &arena) {
  mesh.vertices.bind(arena, Solid::NUM_VERTS);
  mesh.face_counts.bind(arena, Solid::NUM_FACES);
  mesh.faces.bind(arena, Solid::faces.size());
  for (const auto &v : Solid::vertices) mesh.vertices.push_back(v);
  for (auto fc : Solid::face_counts) mesh.face_counts.push_back(fc);
  for (auto fi : Solid::faces)
    mesh.faces.push_back(static_cast<uint16_t>(fi));
}

// ---------------------------------------------------------------------------
// PolyMesh API
// ---------------------------------------------------------------------------

inline void test_polymesh_default_state() {
  PolyMesh m;
  HS_EXPECT_EQ(m.vertices.size(), (size_t)0);
  HS_EXPECT_EQ(m.face_counts.size(), (size_t)0);
  HS_EXPECT_EQ(m.faces.size(), (size_t)0);
}

inline void test_polymesh_bind_arrays() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh m;
  m.vertices.bind(arena, /*verts*/ 8);
  m.face_counts.bind(arena, /*faces*/ 6);
  m.faces.bind(arena, /*indices*/ 24);
  HS_EXPECT_TRUE(m.vertices.is_bound());
  HS_EXPECT_TRUE(m.face_counts.is_bound());
  HS_EXPECT_TRUE(m.faces.is_bound());
  HS_EXPECT_EQ(m.vertices.capacity(), (size_t)8);
  HS_EXPECT_EQ(m.face_counts.capacity(), (size_t)6);
  HS_EXPECT_EQ(m.faces.capacity(), (size_t)24);
}

inline void test_polymesh_clear_resets_data() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh m;
  build_solid<Solids::Cube>(m, arena);
  HS_EXPECT_TRUE(m.vertices.size() > 0);

  m.clear();
  HS_EXPECT_EQ(m.vertices.size(), (size_t)0);
  HS_EXPECT_EQ(m.face_counts.size(), (size_t)0);
  HS_EXPECT_EQ(m.faces.size(), (size_t)0);
}

inline void test_polymesh_unified_accessors() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh m;
  build_solid<Solids::Cube>(m, arena);

  HS_EXPECT_EQ(m.get_face_counts_size(), m.face_counts.size());
  HS_EXPECT_EQ(m.get_faces_size(), m.faces.size());
  HS_EXPECT_TRUE(m.get_face_counts_data() == m.face_counts.data());
  HS_EXPECT_TRUE(m.get_faces_data() == m.faces.data());
}

// ---------------------------------------------------------------------------
// HalfEdgeMesh
// ---------------------------------------------------------------------------

inline void test_half_edge_mesh_size_matches_input() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  HalfEdgeMesh he(arena, cube);
  HS_EXPECT_EQ(he.vertices.size(), cube.vertices.size());
  HS_EXPECT_EQ(he.faces.size(), cube.face_counts.size());
  // One half-edge per face index.
  HS_EXPECT_EQ(he.half_edges.size(), cube.faces.size());
}

inline void test_half_edge_mesh_face_loop_closes() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);
  HalfEdgeMesh he(arena, cube);

  // Walk the half-edge loop around each face and verify it closes after
  // face_counts[fi] steps.
  for (size_t fi = 0; fi < he.faces.size(); ++fi) {
    uint16_t start = he.faces[fi].half_edge;
    uint16_t curr = start;
    int steps = 0;
    do {
      HS_EXPECT_TRUE(curr != HE_NONE);
      HS_EXPECT_EQ(he.half_edges[curr].face, (uint16_t)fi);
      curr = he.half_edges[curr].next;
      steps++;
      if (steps > 100) break; // safety
    } while (curr != start);
    HS_EXPECT_EQ(steps, (int)cube.face_counts[fi]);
  }
}

inline void test_half_edge_mesh_pairs_are_symmetric() {
  // For closed manifolds (cube), every half-edge has a pair, and pairs
  // point back at each other.
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);
  HalfEdgeMesh he(arena, cube);

  for (size_t i = 0; i < he.half_edges.size(); ++i) {
    uint16_t pair = he.half_edges[i].pair;
    HS_EXPECT_TRUE(pair != HE_NONE);
    HS_EXPECT_EQ(he.half_edges[pair].pair, (uint16_t)i);
  }
}

inline void test_half_edge_mesh_euler_invariant() {
  // V - E + F = 2 for a topological sphere. E = half_edges/2 (each edge is two
  // half-edges).
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);
  HalfEdgeMesh he(arena, cube);

  int V = static_cast<int>(he.vertices.size());
  int E = static_cast<int>(he.half_edges.size()) / 2;
  int F = static_cast<int>(he.faces.size());
  HS_EXPECT_EQ(V - E + F, 2);
  HS_EXPECT_EQ(V, 8);
  HS_EXPECT_EQ(E, 12);
  HS_EXPECT_EQ(F, 6);
}

inline void test_half_edge_mesh_euler_tetrahedron() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh tet;
  build_solid<Solids::Tetrahedron>(tet, arena);
  HalfEdgeMesh he(arena, tet);

  int V = static_cast<int>(he.vertices.size());
  int E = static_cast<int>(he.half_edges.size()) / 2;
  int F = static_cast<int>(he.faces.size());
  HS_EXPECT_EQ(V, 4);
  HS_EXPECT_EQ(E, 6);
  HS_EXPECT_EQ(F, 4);
  HS_EXPECT_EQ(V - E + F, 2);
}

inline void test_half_edge_mesh_built_from_meshstate() {
  // The other HalfEdgeMesh constructor takes a MeshState.
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  MeshState ms;
  MeshOps::compile(cube, ms, arena);

  HalfEdgeMesh he(arena, ms);
  HS_EXPECT_EQ(he.vertices.size(), ms.vertices.size());
  HS_EXPECT_EQ(he.faces.size(), ms.face_counts.size());
  HS_EXPECT_EQ(he.half_edges.size(), ms.faces.size());
}

// ---------------------------------------------------------------------------
// MeshOps::compile (PolyMesh → MeshState)
// ---------------------------------------------------------------------------

inline void test_compile_polymesh_to_meshstate_basic() {
  Arena src(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, src);

  MeshState ms;
  MeshOps::compile(cube, ms, dst);

  HS_EXPECT_EQ(ms.vertices.size(), cube.vertices.size());
  HS_EXPECT_EQ(ms.face_counts.size(), cube.face_counts.size());
  HS_EXPECT_EQ(ms.faces.size(), cube.faces.size());
  // face_offsets bound and populated with running offsets
  HS_EXPECT_EQ(ms.face_offsets.size(), cube.face_counts.size());
  size_t expected_offset = 0;
  for (size_t i = 0; i < ms.face_offsets.size(); ++i) {
    HS_EXPECT_EQ(ms.face_offsets[i], (uint16_t)expected_offset);
    expected_offset += ms.face_counts[i];
  }
}

inline void test_compile_drops_degenerate_faces() {
  Arena src(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh m;
  m.vertices.bind(src, /*verts*/ 4);
  m.face_counts.bind(src, /*faces*/ 3);
  m.faces.bind(src, /*indices*/ 7);
  // 4 vertices on a tetrahedron-ish layout (positions don't matter here)
  m.vertices.push_back(Vector(1, 0, 0));
  m.vertices.push_back(Vector(0, 1, 0));
  m.vertices.push_back(Vector(0, 0, 1));
  m.vertices.push_back(Vector(-1, -1, -1));

  // Face counts: 3 (valid), 2 (degenerate), 2 (degenerate)
  m.face_counts.push_back(3);
  m.face_counts.push_back(2);
  m.face_counts.push_back(2);

  m.faces.push_back(0); m.faces.push_back(1); m.faces.push_back(2);
  m.faces.push_back(0); m.faces.push_back(3);
  m.faces.push_back(1); m.faces.push_back(3);

  MeshState ms;
  MeshOps::compile(m, ms, dst);
  // Only one face survives (the triangle).
  HS_EXPECT_EQ(ms.face_counts.size(), (size_t)1);
  HS_EXPECT_EQ(ms.faces.size(), (size_t)3);
  HS_EXPECT_EQ(ms.face_counts[0], (uint8_t)3);
}

// ---------------------------------------------------------------------------
// MeshOps::clone — deep copy
// ---------------------------------------------------------------------------

inline void test_clone_meshstate_deep_copies() {
  Arena src_arena(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst_arena(mesh_arena_b, sizeof(mesh_arena_b));

  // Build a MeshState by going PolyMesh → compile()
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, src_arena);
  MeshState src;
  MeshOps::compile(cube, src, src_arena);

  MeshState dst;
  MeshOps::clone(src, dst, dst_arena);

  HS_EXPECT_EQ(dst.vertices.size(), src.vertices.size());
  HS_EXPECT_EQ(dst.face_counts.size(), src.face_counts.size());
  HS_EXPECT_EQ(dst.faces.size(), src.faces.size());
  HS_EXPECT_EQ(dst.face_offsets.size(), src.face_offsets.size());

  // Independent storage
  HS_EXPECT_TRUE(dst.vertices.data() != src.vertices.data());

  // Values match
  for (size_t i = 0; i < src.vertices.size(); ++i) {
    HS_EXPECT_NEAR(dst.vertices[i].x, src.vertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(dst.vertices[i].y, src.vertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(dst.vertices[i].z, src.vertices[i].z, 1e-6f);
  }
  for (size_t i = 0; i < src.faces.size(); ++i)
    HS_EXPECT_EQ(dst.faces[i], src.faces[i]);
}

inline void test_clone_polymesh_deep_copies() {
  Arena src_arena(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst_arena(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh src;
  build_solid<Solids::Cube>(src, src_arena);

  PolyMesh dst;
  MeshOps::clone(src, dst, dst_arena);

  HS_EXPECT_EQ(dst.vertices.size(), src.vertices.size());
  HS_EXPECT_EQ(dst.face_counts.size(), src.face_counts.size());
  HS_EXPECT_EQ(dst.faces.size(), src.faces.size());
  HS_EXPECT_TRUE(dst.vertices.data() != src.vertices.data());
}

// ---------------------------------------------------------------------------
// MeshOps::classify_faces_by_topology
// ---------------------------------------------------------------------------

inline void test_classify_faces_cube_uniform_topology() {
  // All 6 cube faces are topologically equivalent (square, 4 right angles,
  // same neighbor signature) — must map to the same class.
  Arena geom(mesh_arena_a, sizeof(mesh_arena_a));
  Arena scratch_a(mesh_arena_b, sizeof(mesh_arena_b) / 2);
  Arena scratch_b(mesh_arena_b + sizeof(mesh_arena_b) / 2,
                  sizeof(mesh_arena_b) / 2);

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);

  MeshOps::classify_faces_by_topology(cube, scratch_a, scratch_b, geom);
  HS_EXPECT_EQ(cube.topology.size(), cube.face_counts.size());
  // All 6 faces should share class 0.
  for (size_t i = 0; i < cube.topology.size(); ++i)
    HS_EXPECT_EQ(cube.topology[i], 0);
}

inline void test_classify_faces_tetrahedron_uniform_topology() {
  Arena geom(mesh_arena_a, sizeof(mesh_arena_a));
  Arena scratch_a(mesh_arena_b, sizeof(mesh_arena_b) / 2);
  Arena scratch_b(mesh_arena_b + sizeof(mesh_arena_b) / 2,
                  sizeof(mesh_arena_b) / 2);

  PolyMesh tet;
  build_solid<Solids::Tetrahedron>(tet, geom);

  MeshOps::classify_faces_by_topology(tet, scratch_a, scratch_b, geom);
  HS_EXPECT_EQ(tet.topology.size(), tet.face_counts.size());
  for (size_t i = 0; i < tet.topology.size(); ++i)
    HS_EXPECT_EQ(tet.topology[i], 0);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

inline int run_mesh_tests() {
  auto scope = hs_test::begin_module("mesh");

  test_polymesh_default_state();
  test_polymesh_bind_arrays();
  test_polymesh_clear_resets_data();
  test_polymesh_unified_accessors();

  test_half_edge_mesh_size_matches_input();
  test_half_edge_mesh_face_loop_closes();
  test_half_edge_mesh_pairs_are_symmetric();
  test_half_edge_mesh_euler_invariant();
  test_half_edge_mesh_euler_tetrahedron();
  test_half_edge_mesh_built_from_meshstate();

  test_compile_polymesh_to_meshstate_basic();
  test_compile_drops_degenerate_faces();

  test_clone_meshstate_deep_copies();
  test_clone_polymesh_deep_copies();

  test_classify_faces_cube_uniform_topology();
  test_classify_faces_tetrahedron_uniform_topology();

  return hs_test::end_module(scope);
}

} // namespace mesh_tests
} // namespace hs_test

