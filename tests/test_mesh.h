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
#include "tests/mesh_test_util.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mesh_tests {

inline uint8_t mesh_arena_a[256 * 1024];
inline uint8_t mesh_arena_b[256 * 1024];

// build_solid() lives in tests/mesh_test_util.h.

// ---------------------------------------------------------------------------
// PolyMesh API
// ---------------------------------------------------------------------------

// A freshly constructed PolyMesh owns no storage: all arrays are empty.
inline void test_polymesh_default_state() {
  PolyMesh m;
  HS_EXPECT_EQ(m.vertices.size(), (size_t)0);
  HS_EXPECT_EQ(m.face_counts.size(), (size_t)0);
  HS_EXPECT_EQ(m.faces.size(), (size_t)0);
}

// Binding each array to an arena reserves the requested capacity and marks
// the array bound.
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

// clear() empties all arrays (sizes return to zero) without unbinding them.
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

// The get_*_size()/get_*_data() accessors are thin views onto the underlying
// arrays and must agree with direct member access.
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

// HalfEdgeMesh derived counts track the source PolyMesh: one vertex per
// vertex, one face per face, and one half-edge per face index.
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

// Each face's next-pointer ring is closed and consistent: following next from
// the face's start half-edge returns to it in exactly face_counts[fi] steps,
// and every half-edge on the ring reports that same face.
inline void test_half_edge_mesh_face_loop_closes() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);
  HalfEdgeMesh he(arena, cube);

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

// For a closed manifold (cube), every half-edge has a pair and the pairing is
// reciprocal: pair(pair(i)) == i.
inline void test_half_edge_mesh_pairs_are_symmetric() {
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

// Euler characteristic V - E + F = 2 holds for the cube (a topological sphere),
// with the expected V=8, E=12, F=6. E = half_edges/2 since each edge is two
// half-edges.
inline void test_half_edge_mesh_euler_invariant() {
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

// Same Euler check on a tetrahedron: V=4, E=6, F=4 satisfy V - E + F = 2.
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

// HalfEdgeMesh can also be built from a compiled MeshState; the derived counts
// must track that MeshState just as they do for the PolyMesh constructor.
inline void test_half_edge_mesh_built_from_meshstate() {
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

// compile() copies a clean PolyMesh into a MeshState verbatim and additionally
// fills face_offsets with the running prefix sum of face_counts.
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
  HS_EXPECT_EQ(ms.face_offsets.size(), cube.face_counts.size());
  size_t expected_offset = 0;
  for (size_t i = 0; i < ms.face_offsets.size(); ++i) {
    HS_EXPECT_EQ(ms.face_offsets[i], (uint16_t)expected_offset);
    expected_offset += ms.face_counts[i];
  }
}

// compile() discards faces with fewer than three vertices: of one triangle and
// two 2-vertex faces, only the triangle survives into the MeshState.
inline void test_compile_drops_degenerate_faces() {
  Arena src(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh m;
  m.vertices.bind(src, /*verts*/ 4);
  m.face_counts.bind(src, /*faces*/ 3);
  m.faces.bind(src, /*indices*/ 7);
  // Positions are arbitrary; only the topology (face vertex counts) matters.
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

// clone() of a MeshState yields equal sizes and values but independent storage
// in the destination arena (no aliasing of the source buffers).
inline void test_clone_meshstate_deep_copies() {
  Arena src_arena(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst_arena(mesh_arena_b, sizeof(mesh_arena_b));

  // Build a MeshState via PolyMesh -> compile().
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

  HS_EXPECT_TRUE(dst.vertices.data() != src.vertices.data());

  for (size_t i = 0; i < src.vertices.size(); ++i) {
    HS_EXPECT_NEAR(dst.vertices[i].x, src.vertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(dst.vertices[i].y, src.vertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(dst.vertices[i].z, src.vertices[i].z, 1e-6f);
  }
  for (size_t i = 0; i < src.faces.size(); ++i)
    HS_EXPECT_EQ(dst.faces[i], src.faces[i]);
}

// clone() of a PolyMesh likewise copies all arrays into independent storage.
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

// All 6 cube faces are topologically equivalent (square, 4 right angles, same
// neighbor signature), so classify_faces_by_topology must assign them all the
// same class (0). scratch_a/scratch_b are working arenas for the classifier.
inline void test_classify_faces_cube_uniform_topology() {
  Arena geom(mesh_arena_a, sizeof(mesh_arena_a));
  Arena scratch_a(mesh_arena_b, sizeof(mesh_arena_b) / 2);
  Arena scratch_b(mesh_arena_b + sizeof(mesh_arena_b) / 2,
                  sizeof(mesh_arena_b) / 2);

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, geom);

  MeshOps::classify_faces_by_topology(cube, scratch_a, scratch_b, geom);
  HS_EXPECT_EQ(cube.topology.size(), cube.face_counts.size());
  for (size_t i = 0; i < cube.topology.size(); ++i)
    HS_EXPECT_EQ(cube.topology[i], 0);
}

// A tetrahedron's 4 faces are also mutually equivalent triangles, so they too
// must all classify into the same class (0).
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

// Run every mesh test case in order; returns the module's failure count.
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

