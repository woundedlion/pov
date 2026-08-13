/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/mesh/mesh.h and core/mesh/triangular_bitset.h.
 *
 * Coverage:
 *   - PolyMesh basic API (initialize, clear, accessors)
 *   - HalfEdgeMesh construction from PolyMesh and MeshState
 *     (Euler invariant V - E + F = 2, half-edge pairing, vertex/face refs)
 *   - MeshOps::compile (PolyMesh → MeshState, drops degenerate faces,
 *     populates face_offsets)
 *   - MeshOps::clone (deep copy into a target arena)
 *   - MeshOps::classify_faces_by_topology (cube → all faces same class)
 *   - TriangularBitset pair dedup (sizing, index bijection, bit independence)
 */
#pragma once

#include <cstdint>
#include <span>
#include "core/mesh/mesh.h"
#include "core/mesh/solids.h"
#include "core/mesh/triangular_bitset.h"
#include "tests/mesh_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mesh_tests {

inline uint8_t mesh_arena_a[1024 * 1024];
inline uint8_t mesh_arena_b[1024 * 1024];
inline uint8_t mesh_arena_c[1024 * 1024];

// ---------------------------------------------------------------------------
// PolyMesh API
// ---------------------------------------------------------------------------

/**
 * @brief Verifies a freshly constructed PolyMesh owns no storage (all arrays empty).
 */
inline void test_polymesh_default_state() {
  PolyMesh m;
  HS_EXPECT_EQ(m.vertices.size(), (size_t)0);
  HS_EXPECT_EQ(m.face_counts.size(), (size_t)0);
  HS_EXPECT_EQ(m.faces.size(), (size_t)0);
}

/**
 * @brief Verifies binding each array to an arena reserves the requested
 *        capacity and marks the array bound.
 */
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

/**
 * @brief Verifies clear() empties all arrays (sizes return to zero) without
 *        unbinding them.
 */
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

/**
 * @brief Verifies the get_*_size()/get_*_data() accessors are thin views onto
 *        the underlying arrays and agree with direct member access.
 */
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

/**
 * @brief Verifies HalfEdgeMesh derived counts track the source PolyMesh: one
 *        face per face, and one half-edge per face index.
 */
inline void test_half_edge_mesh_size_matches_input() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  HalfEdgeMesh he(arena, cube);
  HS_EXPECT_EQ(he.faces.size(), cube.face_counts.size());
  HS_EXPECT_EQ(he.half_edges.size(), cube.faces.size());
}

/**
 * @brief Verifies each face's next-pointer ring is closed and consistent.
 * @details Following next from the face's start half-edge returns to it in
 *          exactly face_counts[fi] steps, and every half-edge on the ring
 *          reports that same face.
 */
inline void test_half_edge_mesh_face_loop_closes() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);
  HalfEdgeMesh he(arena, cube);

  // A valid ring visits each half-edge at most once, so anything longer than the
  // total half-edge count cannot close — a mesh-derived bound, not a constant.
  const int max_steps = static_cast<int>(he.half_edges.size()) + 1;
  for (size_t fi = 0; fi < he.faces.size(); ++fi) {
    uint16_t start = he.faces[fi].half_edge;
    uint16_t curr = start;
    int steps = 0;
    do {
      HS_EXPECT_TRUE(curr != HE_NONE);
      if (curr == HE_NONE)
        break; // don't index half_edges with the sentinel
      HS_EXPECT_EQ(he.half_edges[curr].face, (uint16_t)fi);
      curr = he.half_edges[curr].next;
      steps++;
      if (steps > max_steps)
        break; // ring failed to close
    } while (curr != start);
    HS_EXPECT_EQ(steps, (int)cube.face_counts[fi]);
  }
}

/**
 * @brief Verifies that for a closed manifold (cube), every half-edge has a pair,
 *        the pairing is reciprocal (pair(pair(i)) == i), AND the twin is the
 *        geometric opposite — it traverses the same undirected edge in reverse.
 * @details Reciprocity alone is too weak: a pairing that consistently joined the
 *          wrong edges (e.g. swapped two twins) could still satisfy
 *          pair(pair(i)) == i. The opposite-direction check pins the geometry:
 *          for a half-edge u→v, its twin must run v→u, so head(pair) == tail(i)
 *          and tail(pair) == head(i). `he.vertex` is the head (destination); the
 *          tail is the head of the previous half-edge in the face loop.
 */
inline void test_half_edge_mesh_pairs_are_symmetric() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);
  HalfEdgeMesh he(arena, cube);

  for (size_t i = 0; i < he.half_edges.size(); ++i) {
    uint16_t pair = he.half_edges[i].pair;
    HS_EXPECT_TRUE(pair != HE_NONE);
    HS_EXPECT_EQ(he.half_edges[pair].pair, (uint16_t)i);

    // Same undirected edge, reversed endpoints: twin of u->v must be v->u.
    const uint16_t head_i = he.half_edges[i].vertex;
    const uint16_t tail_i = he.half_edges[he.half_edges[i].prev].vertex;
    const uint16_t head_p = he.half_edges[pair].vertex;
    const uint16_t tail_p = he.half_edges[he.half_edges[pair].prev].vertex;
    HS_EXPECT_EQ(head_p, tail_i);
    HS_EXPECT_EQ(tail_p, head_i);
  }
}

/**
 * @brief Verifies the half-edge builder's boundary path: an open mesh leaves its
 *        outer edges unpaired (pair == HE_NONE).
 * @details Every closed-solid fixture above pairs *every* half-edge, so the
 *          boundary branch of pair_half_edges (a run of one half-edge on an
 *          undirected edge) is otherwise unexercised. Two triangles (0,1,2) and
 *          (0,2,3) share edge 0-2: that single interior edge pairs reciprocally;
 *          the other four half-edges border the open boundary and stay HE_NONE.
 */
inline void test_half_edge_mesh_open_boundary_edges() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh open;
  open.vertices.bind(arena, 4);
  open.face_counts.bind(arena, 2);
  open.faces.bind(arena, 6);
  // Positions are irrelevant to half-edge topology; any distinct points work.
  open.vertices.push_back(Vector(0, 0, 1));
  open.vertices.push_back(Vector(1, 0, 0));
  open.vertices.push_back(Vector(0, 1, 0));
  open.vertices.push_back(Vector(-1, 0, 0));
  open.face_counts.push_back(3);
  open.face_counts.push_back(3);
  const uint16_t idx[] = {0, 1, 2, 0, 2, 3}; // tri A: 0->1->2, tri B: 0->2->3
  for (uint16_t i : idx)
    open.faces.push_back(i);

  HalfEdgeMesh he(arena, open);

  HS_EXPECT_EQ(he.half_edges.size(), (size_t)6);
  int paired = 0, boundary = 0;
  for (size_t i = 0; i < he.half_edges.size(); ++i) {
    uint16_t pair = he.half_edges[i].pair;
    if (pair == HE_NONE) {
      boundary++;
    } else {
      paired++;
      HS_EXPECT_EQ(he.half_edges[pair].pair, (uint16_t)i);
    }
  }
  HS_EXPECT_EQ(paired, 2);   // the shared edge 0-2, both directions
  HS_EXPECT_EQ(boundary, 4); // the four outer edges have no twin
}

/**
 * @brief Verifies the Euler characteristic V - E + F = 2 holds for the cube (a
 *        topological sphere), with the expected V=8, E=12, F=6.
 * @details E = half_edges/2 since each edge is two half-edges.
 */
inline void test_half_edge_mesh_euler_invariant() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);
  HalfEdgeMesh he(arena, cube);

  int V = static_cast<int>(cube.vertices.size());
  int E = static_cast<int>(he.half_edges.size()) / 2;
  int F = static_cast<int>(he.faces.size());
  HS_EXPECT_EQ(V - E + F, 2);
  HS_EXPECT_EQ(V, 8);
  HS_EXPECT_EQ(E, 12);
  HS_EXPECT_EQ(F, 6);
}

/**
 * @brief Verifies the Euler check on a tetrahedron: V=4, E=6, F=4 satisfy
 *        V - E + F = 2.
 */
inline void test_half_edge_mesh_euler_tetrahedron() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh tet;
  build_solid<Solids::Tetrahedron>(tet, arena);
  HalfEdgeMesh he(arena, tet);

  int V = static_cast<int>(tet.vertices.size());
  int E = static_cast<int>(he.half_edges.size()) / 2;
  int F = static_cast<int>(he.faces.size());
  HS_EXPECT_EQ(V, 4);
  HS_EXPECT_EQ(E, 6);
  HS_EXPECT_EQ(F, 4);
  HS_EXPECT_EQ(V - E + F, 2);
}

/**
 * @brief Verifies HalfEdgeMesh can be built from a compiled MeshState, with the
 *        derived counts tracking that MeshState as for the PolyMesh constructor.
 */
inline void test_half_edge_mesh_built_from_meshstate() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  MeshState ms;
  MeshOps::compile(cube, ms, arena, scratch_arena_a);

  HalfEdgeMesh he(arena, ms);
  HS_EXPECT_EQ(he.faces.size(), ms.face_counts.size());
  HS_EXPECT_EQ(he.half_edges.size(), ms.faces.size());
}

// ---------------------------------------------------------------------------
// MeshOps::compile (PolyMesh → MeshState)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies compile() copies a clean PolyMesh into a MeshState verbatim
 *        and fills face_offsets with the running prefix sum of face_counts.
 */
inline void test_compile_polymesh_to_meshstate_basic() {
  Arena src(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, src);

  MeshState ms;
  MeshOps::compile(cube, ms, dst, scratch_arena_a);

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

/**
 * @brief Verifies compile() discards faces with fewer than three vertices and
 *        compacts away the vertices they orphaned.
 * @details Of one triangle and two 2-vertex faces, only the triangle survives
 *          into the MeshState. Vertex 3 is referenced only by the dropped faces,
 *          so it is compacted out while the survivor keeps its remapped indices.
 */
inline void test_compile_drops_degenerate_faces() {
  Arena src(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh m;
  m.vertices.bind(src, /*verts*/ 4);
  m.face_counts.bind(src, /*faces*/ 3);
  m.faces.bind(src, /*indices*/ 7);
  // The triangle references vertices 0-2; vertex 3 only the degenerate faces.
  m.vertices.push_back(Vector(1, 0, 0));
  m.vertices.push_back(Vector(0, 1, 0));
  m.vertices.push_back(Vector(0, 0, 1));
  m.vertices.push_back(Vector(-1, -1, -1));

  // Face counts: 3 (valid), then two degenerate 2-vertex faces.
  m.face_counts.push_back(3);
  m.face_counts.push_back(2);
  m.face_counts.push_back(2);

  m.faces.push_back(0);
  m.faces.push_back(1);
  m.faces.push_back(2);
  m.faces.push_back(0);
  m.faces.push_back(3);
  m.faces.push_back(1);
  m.faces.push_back(3);

  MeshState ms;
  MeshOps::compile(m, ms, dst, scratch_arena_a);
  HS_EXPECT_EQ(ms.face_counts.size(), (size_t)1);
  HS_EXPECT_EQ(ms.faces.size(), (size_t)3);
  HS_EXPECT_EQ(ms.face_counts[0], (uint8_t)3);

  // Orphan vertex 3 is compacted out; the kept vertices retain source order so
  // the surviving triangle's indices are unchanged.
  HS_EXPECT_EQ(ms.vertices.size(), (size_t)3);
  HS_EXPECT_EQ(ms.faces[0], (uint16_t)0);
  HS_EXPECT_EQ(ms.faces[1], (uint16_t)1);
  HS_EXPECT_EQ(ms.faces[2], (uint16_t)2);
}

// ---------------------------------------------------------------------------
// MeshOps::clone — deep copy
// ---------------------------------------------------------------------------

/**
 * @brief Verifies clone() of a MeshState yields equal sizes and values but
 *        independent storage in the destination arena (no aliasing of source).
 * @details compile() leaves topology empty, so it is filled with distinct
 *          per-face values first: clone() is the only path that carries it.
 */
inline void test_clone_meshstate_deep_copies() {
  Arena src_arena(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst_arena(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, src_arena);
  MeshState src;
  MeshOps::compile(cube, src, src_arena, scratch_arena_a);

  src.topology.bind(src_arena, src.face_counts.size());
  for (size_t i = 0; i < src.face_counts.size(); ++i)
    src.topology.push_back(static_cast<uint16_t>(100 + i));

  MeshState dst;
  MeshOps::clone(src, dst, dst_arena);

  HS_EXPECT_EQ(dst.vertices.size(), src.vertices.size());
  HS_EXPECT_EQ(dst.face_counts.size(), src.face_counts.size());
  HS_EXPECT_EQ(dst.faces.size(), src.faces.size());
  HS_EXPECT_EQ(dst.face_offsets.size(), src.face_offsets.size());
  HS_EXPECT_EQ(dst.topology.size(), src.topology.size());

  HS_EXPECT_TRUE(dst.vertices.data() != src.vertices.data());
  HS_EXPECT_TRUE(dst.topology.data() != src.topology.data());

  for (size_t i = 0; i < src.vertices.size(); ++i) {
    HS_EXPECT_NEAR(dst.vertices[i].x, src.vertices[i].x, 1e-6f);
    HS_EXPECT_NEAR(dst.vertices[i].y, src.vertices[i].y, 1e-6f);
    HS_EXPECT_NEAR(dst.vertices[i].z, src.vertices[i].z, 1e-6f);
  }
  for (size_t i = 0; i < src.face_counts.size(); ++i)
    HS_EXPECT_EQ(dst.face_counts[i], src.face_counts[i]);
  for (size_t i = 0; i < src.faces.size(); ++i)
    HS_EXPECT_EQ(dst.faces[i], src.faces[i]);
  for (size_t i = 0; i < src.face_offsets.size(); ++i)
    HS_EXPECT_EQ(dst.face_offsets[i], src.face_offsets[i]);
  for (size_t i = 0; i < src.topology.size(); ++i)
    HS_EXPECT_EQ(dst.topology[i], src.topology[i]);
}

/**
 * @brief Verifies clone() of a PolyMesh likewise copies all arrays into
 *        independent storage.
 * @details No Conway operator propagates topology, so clone() is the only path
 *          that carries it; distinct per-face values pin the copy element-wise.
 */
inline void test_clone_polymesh_deep_copies() {
  Arena src_arena(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst_arena(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh src;
  build_solid<Solids::Cube>(src, src_arena);
  src.topology.bind(src_arena, src.face_counts.size());
  for (size_t i = 0; i < src.face_counts.size(); ++i)
    src.topology.push_back(static_cast<int>(100 + i));

  PolyMesh dst;
  MeshOps::clone(src, dst, dst_arena);

  HS_EXPECT_EQ(dst.vertices.size(), src.vertices.size());
  HS_EXPECT_EQ(dst.face_counts.size(), src.face_counts.size());
  HS_EXPECT_EQ(dst.faces.size(), src.faces.size());
  HS_EXPECT_EQ(dst.topology.size(), src.topology.size());
  HS_EXPECT_TRUE(dst.vertices.data() != src.vertices.data());
  HS_EXPECT_TRUE(dst.topology.data() != src.topology.data());
  for (size_t i = 0; i < src.topology.size(); ++i)
    HS_EXPECT_EQ(dst.topology[i], src.topology[i]);
}

// ---------------------------------------------------------------------------
// MeshOps::classify_faces_by_topology
// ---------------------------------------------------------------------------

/**
 * @brief Verifies all 6 cube faces classify into the same class (0).
 * @details The faces are topologically equivalent (square, 4 right angles, same
 *          neighbor signature). scratch_a/scratch_b are working arenas for the
 *          classifier.
 */
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

/**
 * @brief Verifies classify_faces_by_topology on an UNCOMPILED PolyMesh with a
 *        degenerate 2-gon neither self-pairs it nor trips the non-manifold trap.
 * @details compile() drops sub-triangle faces, but the public PolyMesh overload
 *          runs on raw input. A 2-gon reusing the triangle's edge would put a
 *          third half-edge on that edge and trap without the record-loop
 *          side-count guard; the guard leaves the degenerate edges unpaired so
 *          the triangle still classifies.
 */
inline void test_classify_faces_uncompiled_degenerate() {
  Arena geom(mesh_arena_a, sizeof(mesh_arena_a));
  Arena scratch_a(mesh_arena_b, sizeof(mesh_arena_b) / 2);
  Arena scratch_b(mesh_arena_b + sizeof(mesh_arena_b) / 2,
                  sizeof(mesh_arena_b) / 2);

  PolyMesh m;
  m.vertices.bind(geom, 3);
  m.face_counts.bind(geom, 2);
  m.faces.bind(geom, 5);
  m.vertices.push_back(Vector(1, 0, 0));
  m.vertices.push_back(Vector(0, 1, 0));
  m.vertices.push_back(Vector(0, 0, 1));
  m.face_counts.push_back(3);
  m.face_counts.push_back(2);
  m.faces.push_back(0);
  m.faces.push_back(1);
  m.faces.push_back(2);
  m.faces.push_back(0);
  m.faces.push_back(1);

  MeshOps::classify_faces_by_topology(m, scratch_a, scratch_b, geom);
  HS_EXPECT_SIZE_OR_RETURN(m.topology, m.face_counts.size());
  // The triangle still classifies: it lands in its own class rather than being
  // merged with the degenerate 2-gon, and the two ids stay dense.
  HS_EXPECT_NE(m.topology[0], m.topology[1]);
  HS_EXPECT_EQ(std::min(m.topology[0], m.topology[1]),
               static_cast<uint16_t>(0));
  HS_EXPECT_EQ(std::max(m.topology[0], m.topology[1]),
               static_cast<uint16_t>(1));
}

/** @brief Verifies the degenerate sentinel cannot collide with a real edge. */
inline void test_degenerate_edge_records_never_pair() {
  HalfEdgePairRecord records[2];
  fill_edge_record(records[0], HE_NONE, HE_NONE, 0);
  fill_edge_record(records[1], 0, 0, 1);
  int pairs = 0;
  pair_half_edges(records, 2, [&](uint16_t, uint16_t) { ++pairs; });
  HS_EXPECT_EQ(pairs, 0);
}

/**
 * @brief Verifies a tetrahedron's 4 faces, being mutually equivalent triangles,
 *        all classify into the same class (0).
 */
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

/**
 * @brief Verifies classify_faces_by_topology accepts one arena in all three
 *        roles, matching the three-arena classification exactly.
 * @details HankinSolids classifies its bookend mesh with scratch_arena_a passed
 *          for both scratch parameters and for persistent, so the aliasing is a
 *          supported contract, not an accident of the current allocation order.
 */
inline void test_classify_faces_aliased_arenas() {
  Arena target(mesh_arena_a, sizeof(mesh_arena_a));
  Arena one(mesh_arena_c, sizeof(mesh_arena_c));

  PolyMesh tr, tr_aliased;
  {
    Arena temp(mesh_arena_b, sizeof(mesh_arena_b));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, temp);
    tr = MeshOps::truncate(cube, target, temp, 0.25f);
    tr_aliased = MeshOps::truncate(cube, one, temp, 0.25f);
  }

  Arena scratch_a(mesh_arena_b, sizeof(mesh_arena_b) / 2);
  Arena scratch_b(mesh_arena_b + sizeof(mesh_arena_b) / 2,
                  sizeof(mesh_arena_b) / 2);
  MeshOps::classify_faces_by_topology(tr, scratch_a, scratch_b, target);
  MeshOps::classify_faces_by_topology(tr_aliased, one, one, one);

  HS_EXPECT_EQ(tr.topology.size(), tr.face_counts.size());
  HS_EXPECT_EQ(tr_aliased.topology.size(), tr.topology.size());
  int max_id = 0;
  for (size_t i = 0; i < tr.topology.size(); ++i) {
    HS_EXPECT_EQ(tr_aliased.topology[i], tr.topology[i]);
    max_id = std::max(max_id, static_cast<int>(tr.topology[i]));
  }
  HS_EXPECT_EQ(max_id, 1); // 6 octagons + 8 triangles
}

/**
 * @brief Verifies the classifier separates genuinely distinct face classes.
 * @details truncate(cube) is a mixed-face solid: 6 octagons + 8 triangles. The
 *          cube/tetrahedron tests above are all-equivalent faces, so the
 *          expected output is trivially all-zeros and would pass even if the
 *          neighbor-hash folding, HashNode sort, and dense-id assignment were
 *          broken. This exercises that discriminating logic: it must yield
 *          exactly two classes, group all octagons under one id and all
 *          triangles under another, and assign dense ids {0, 1}.
 */
inline void test_classify_faces_truncated_cube_distinct_topology() {
  Arena target(mesh_arena_a, sizeof(mesh_arena_a));
  Arena temp(mesh_arena_b, sizeof(mesh_arena_b));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);
  PolyMesh tr = MeshOps::truncate(cube, target, temp, 0.25f);
  HS_EXPECT_EQ(tr.face_counts.size(), (size_t)14); // 6 octagons + 8 triangles

  // truncate's temp working set is no longer referenced; reuse mesh_arena_b for
  // the classifier scratch. topology grows into `target`, where tr lives.
  Arena scratch_a(mesh_arena_b, sizeof(mesh_arena_b) / 2);
  Arena scratch_b(mesh_arena_b + sizeof(mesh_arena_b) / 2,
                  sizeof(mesh_arena_b) / 2);
  MeshOps::classify_faces_by_topology(tr, scratch_a, scratch_b, target);

  HS_EXPECT_EQ(tr.topology.size(), tr.face_counts.size());

  // Faces of the same side count must share an id; the two shapes must differ.
  int octagon_id = -1, triangle_id = -1, max_id = 0;
  for (size_t i = 0; i < tr.face_counts.size(); ++i) {
    int sides = tr.face_counts[i];
    int id = tr.topology[i];
    if (id > max_id)
      max_id = id;
    if (sides == 8) {
      if (octagon_id < 0)
        octagon_id = id;
      HS_EXPECT_EQ(id, octagon_id);
    } else if (sides == 3) {
      if (triangle_id < 0)
        triangle_id = id;
      HS_EXPECT_EQ(id, triangle_id);
    }
  }
  HS_EXPECT_TRUE(octagon_id >= 0 && triangle_id >= 0);
  HS_EXPECT_TRUE(octagon_id != triangle_id);
  // Dense ids: two classes -> {0, 1}.
  HS_EXPECT_EQ(max_id, 1);
}

/** Largest face side count the topology sweep's per-face buffers hold. */
inline constexpr int MAX_TOPO_SIDES = 64;

/**
 * @brief A face topology key plus the base hash the classifier derives from it.
 */
struct FaceTopoRecord {
  int count;                  /**< Side count. */
  int angles[MAX_TOPO_SIDES]; /**< Sorted whole-degree interior angles. */
  uint32_t hash;              /**< MeshOps base topology hash for this key. */
};

/**
 * @brief Recomputes a face's canonical key and asks the classifier for its base
 *        topology hash.
 * @details The interior angle uses the classifier's own formulation (unnormalized
 *   edge dot over the length product, fast_acos, rounded to whole degrees) on
 *   purpose: a normalize-first variant lands on the other side of a degree
 *   boundary on relaxed meshes, which would split classes the classifier merges
 *   without either result being wrong. The key's identity, the neighbour
 *   canonicalization and the class partition stay independent.
 */
inline FaceTopoRecord face_topo_record(const PolyMesh &mesh,
                                       const uint16_t *idx, int count) {
  FaceTopoRecord rec;
  rec.count = count;
  for (int k = 0; k < MAX_TOPO_SIDES; ++k)
    rec.angles[k] = 0;
  HS_CHECK(count <= MAX_TOPO_SIDES,
           "face_topo_record: face sides overrun angles[]");
  if (count >= 3) {
    for (int k = 0; k < count; ++k) {
      const Vector &prev = mesh.vertices[idx[(k - 1 + count) % count]];
      const Vector &curr = mesh.vertices[idx[k]];
      const Vector &next = mesh.vertices[idx[(k + 1) % count]];
      const Vector e1 = prev - curr;
      const Vector e2 = next - curr;
      const float m1 = dot(e1, e1), m2 = dot(e2, e2);
      float ang = 0.0f;
      if (m1 > math::EPS_LEN_SQ && m2 > math::EPS_LEN_SQ)
        ang = fast_acos(hs::clamp(dot(e1, e2) / sqrtf(m1 * m2), -1.0f, 1.0f));
      rec.angles[k] = static_cast<int>(std::round(ang * 180.0f / PI_F));
    }
    std::sort(rec.angles, rec.angles + count);
  }
  rec.hash = MeshOps::face_topology_base_hash(count, rec.angles);
  return rec;
}

/**
 * @brief 64-bit FNV-1a over a byte span.
 * @param data Bytes to hash.
 * @param n Byte count.
 * @param h Seed, letting a key be built from several spans.
 * @return The accumulated hash.
 */
inline uint64_t fnv1a64(const void *data, size_t n,
                        uint64_t h = 0xcbf29ce484222325ull) {
  const uint8_t *p = static_cast<const uint8_t *>(data);
  for (size_t i = 0; i < n; ++i) {
    h ^= p[i];
    h *= 0x100000001b3ull;
  }
  return h;
}

/**
 * @brief Reference identity of a face's canonical key (count + sorted angles).
 */
inline uint64_t topo_key_id(const FaceTopoRecord &rec) {
  const uint64_t h = fnv1a64(&rec.count, sizeof(rec.count));
  return fnv1a64(rec.angles, sizeof(rec.angles[0]) * rec.count, h);
}

/** Upper bound on the dense topology ids one roster solid can carry. */
inline constexpr int MAX_TOPO_CLASSES = 256;
/** Slots in the collision sweep's hash table; a power of two. */
inline constexpr int TOPO_HASH_SLOTS = 8192;
/** Faces the collision sweep holds per mesh. */
inline constexpr size_t MAX_SWEEP_FACES = 8192;

/**
 * @brief Open-addressed classifier-hash -> reference-key table.
 * @details Every insert expects an earlier entry under the same 32-bit
 *          classifier hash to carry the same reference key; a mismatch is a
 *          collision that merges two distinct topologies into one class.
 */
struct TopoHashTable {
  uint32_t hashes[TOPO_HASH_SLOTS] = {};
  uint64_t keys[TOPO_HASH_SLOTS] = {};
  bool used[TOPO_HASH_SLOTS] = {};
  int n = 0; /**< Occupied slots. */

  /** @brief Empties the table so a repeated module run starts clean. */
  void clear() {
    for (int i = 0; i < TOPO_HASH_SLOTS; ++i)
      used[i] = false;
    n = 0;
  }

  /** @brief Records one (classifier hash, reference key) pair. */
  void insert(uint32_t hash, uint64_t key) {
    HS_EXPECT_TRUE(n < TOPO_HASH_SLOTS);
    if (n >= TOPO_HASH_SLOTS)
      return;
    for (int i = static_cast<int>(hash % TOPO_HASH_SLOTS);;
         i = (i + 1) % TOPO_HASH_SLOTS) {
      if (!used[i]) {
        used[i] = true;
        hashes[i] = hash;
        keys[i] = key;
        ++n;
        return;
      }
      if (hashes[i] == hash) {
        HS_EXPECT_EQ(keys[i], key);
        return;
      }
    }
  }
};

/**
 * @brief Verifies both classifier hash stages are collision-free across all
 *        three solid registries.
 * @details A fmix32 collision merges two distinct face topologies into one
 *          class. Builds every roster solid and sweeps the pre-fold hash
 *          (count + sorted angles) and the neighbour-folded hash the topology
 *          ids are actually derived from, asserting each maps to one key.
 *
 *          The reference fold above is not the classifier, so it is also run
 *          against the real classify_faces_by_topology output: within each mesh
 *          the production ids and the reference keys must induce the SAME
 *          partition of faces, in both directions. Agreement between the hash
 *          primitives alone would survive a classifier that wired its stages
 *          up differently or skipped the fold entirely.
 */
inline void test_classify_faces_roster_hash_collision_free() {
  static TopoHashTable pre_fold;
  static TopoHashTable folded;
  static uint32_t face_hashes[MAX_SWEEP_FACES];
  static uint64_t face_keys[MAX_SWEEP_FACES];
  static uint64_t folded_keys[MAX_SWEEP_FACES];
  pre_fold.clear();
  folded.clear();
  int swept_meshes = 0;
  int classified_meshes = 0;

  for (std::span<const Solids::Entry> reg : Solids::all_registries()) {
    for (const Solids::Entry &entry : reg) {
      Arena a(mesh_arena_a, sizeof(mesh_arena_a));
      Arena b(mesh_arena_b, sizeof(mesh_arena_b));
      Arena c(mesh_arena_c, sizeof(mesh_arena_c));
      PolyMesh mesh = entry.generate(a, b);

      const uint8_t *fc = mesh.get_face_counts_data();
      const uint16_t *faces = mesh.get_faces_data();
      const size_t F = mesh.get_face_counts_size();
      HS_EXPECT_TRUE(F <= MAX_SWEEP_FACES);
      if (F > MAX_SWEEP_FACES)
        continue;

      size_t off = 0;
      for (size_t f = 0; f < F; ++f) {
        const int count = fc[f];
        const FaceTopoRecord rec = face_topo_record(mesh, faces + off, count);
        off += count;
        face_hashes[f] = rec.hash;
        face_keys[f] = topo_key_id(rec);
        pre_fold.insert(rec.hash, face_keys[f]);
      }

      // Neighbour fold: own hash mixed with the sum of the paired faces' mixed
      // hashes, matching classify_faces_impl's second stage.
      HalfEdgeMesh he(c, mesh);
      size_t he_off = 0;
      for (size_t f = 0; f < F; ++f) {
        const int count = fc[f];
        uint64_t neighbor_keys[MAX_TOPO_SIDES];
        int n_neighbors = 0;
        uint32_t neighbor_acc = 0;
        for (int k = 0; k < count; ++k) {
          const uint16_t pair = he.half_edges[he_off + k].pair;
          if (pair == HE_NONE)
            continue;
          const uint16_t neighbor = he.half_edges[pair].face;
          uint32_t neigh_h = face_hashes[neighbor];
          MeshOps::hash_combine(neigh_h, 0);
          neighbor_acc += MeshOps::fmix32(neigh_h);
          neighbor_keys[n_neighbors++] = face_keys[neighbor];
        }
        std::sort(neighbor_keys, neighbor_keys + n_neighbors);
        folded_keys[f] =
            fnv1a64(neighbor_keys, sizeof(neighbor_keys[0]) * n_neighbors,
                    face_keys[f]);
        folded.insert(
            MeshOps::fold_face_topology_hash(face_hashes[f], neighbor_acc),
            folded_keys[f]);
        he_off += count;
      }
      ++swept_meshes;

      // Run the real classifier and require it to induce the same partition.
      // `c` held only the half-edge mesh, which is dead now, so it backs the
      // classifier scratch; topology grows into `a`, where the mesh lives.
      Arena scratch_a(mesh_arena_c, sizeof(mesh_arena_c) / 2);
      Arena scratch_b(mesh_arena_c + sizeof(mesh_arena_c) / 2,
                      sizeof(mesh_arena_c) / 2);
      MeshOps::classify_faces_by_topology(mesh, scratch_a, scratch_b, a);
      HS_EXPECT_SIZE_OR_RETURN(mesh.topology, F);

      // Dense ids, so a per-id slot table is enough for one direction; the
      // class count is small, so the reverse direction is a linear scan.
      uint64_t id_key[MAX_TOPO_CLASSES];
      bool id_used[MAX_TOPO_CLASSES] = {};
      int classes = 0;
      for (size_t f = 0; f < F; ++f) {
        const int id = mesh.topology[f];
        HS_CONTEXT(entry.name, static_cast<long long>(f));
        HS_EXPECT_TRUE(id >= 0 && id < MAX_TOPO_CLASSES);
        if (id < 0 || id >= MAX_TOPO_CLASSES)
          continue;
        if (!id_used[id]) {
          id_used[id] = true;
          id_key[id] = folded_keys[f];
          ++classes;
          continue;
        }
        // Same production class must mean the same reference topology.
        HS_EXPECT_EQ(id_key[id], folded_keys[f]);
      }
      for (size_t f = 0; f < F; ++f) {
        HS_CONTEXT(entry.name, static_cast<long long>(f));
        int matching = -1;
        for (int id = 0; id < MAX_TOPO_CLASSES; ++id)
          if (id_used[id] && id_key[id] == folded_keys[f]) {
            matching = id;
            break;
          }
        // Same reference topology must mean the same production class.
        HS_EXPECT_EQ(matching, static_cast<int>(mesh.topology[f]));
      }
      HS_EXPECT_GT(classes, 0);
      ++classified_meshes;
    }
  }

  HS_EXPECT_EQ(swept_meshes, Solids::NUM_ENTRIES);
  HS_EXPECT_EQ(classified_meshes, swept_meshes);
  HS_EXPECT_TRUE(pre_fold.n > 0 && folded.n > 0);
}

// ---------------------------------------------------------------------------
// TriangularBitset
// ---------------------------------------------------------------------------

/**
 * @brief Verifies BITS/BYTES are the compile-time size of the upper-triangular
 *        pair matrix.
 * @details MAX_V*(MAX_V-1)/2 bits, ceil(bits/8) bytes.
 */
inline void test_tribitset_sizes() {
  // MAX_V=10 → 10*9/2 = 45 bits → ceil(45/8) = 6 bytes
  HS_EXPECT_EQ(TriangularBitset<10>::BITS, 45);
  HS_EXPECT_EQ(TriangularBitset<10>::BYTES, 6);
  // MAX_V=16 → 16*15/2 = 120 bits = 15 bytes
  HS_EXPECT_EQ(TriangularBitset<16>::BITS, 120);
  HS_EXPECT_EQ(TriangularBitset<16>::BYTES, 15);
}

/**
 * @brief Verifies that after clear(), every pair reads unset.
 */
inline void test_tribitset_clear_and_test() {
  TriangularBitset<8> bs;
  bs.clear();
  for (int a = 0; a < 8; ++a)
    for (int b = a + 1; b < 8; ++b)
      HS_EXPECT_FALSE(bs.test(a, b));
}

/**
 * @brief Verifies test_and_set() returns the prior bit (false on first insert,
 *        true thereafter) and leaves neighboring pairs untouched.
 */
inline void test_tribitset_test_and_set() {
  TriangularBitset<10> bs;
  bs.clear();
  HS_EXPECT_FALSE(bs.test_and_set(2, 5));
  HS_EXPECT_TRUE(bs.test_and_set(2, 5));
  HS_EXPECT_TRUE(bs.test(2, 5));
  HS_EXPECT_FALSE(bs.test(2, 6));
  HS_EXPECT_FALSE(bs.test(3, 5));
}

/**
 * @brief Verifies index(a,b) is a bijection from the MAX_V*(MAX_V-1)/2 unordered
 *        pairs onto the bit range [0, BITS).
 * @details No two pairs collide, and all indices are covered.
 */
inline void test_tribitset_index_uniqueness() {
  constexpr int N = 8;
  int seen[64] = {};
  int unique = 0;
  for (int a = 0; a < N; ++a) {
    for (int b = a + 1; b < N; ++b) {
      int idx = TriangularBitset<N>::index(a, b);
      HS_EXPECT_TRUE(idx >= 0 && idx < TriangularBitset<N>::BITS);
      HS_EXPECT_FALSE(seen[idx]);
      seen[idx] = 1;
      ++unique;
    }
  }
  HS_EXPECT_EQ(unique, TriangularBitset<N>::BITS);
}

/**
 * @brief Verifies setting one set of pairs leaves disjoint pairs unset (bits are
 *        independent).
 */
inline void test_tribitset_all_pairs_independent() {
  TriangularBitset<8> bs;
  bs.clear();
  for (int i = 0; i + 2 < 8; ++i)
    bs.test_and_set(i, i + 2);
  for (int i = 0; i + 2 < 8; ++i)
    HS_EXPECT_TRUE(bs.test(i, i + 2));
  for (int i = 0; i + 1 < 8; ++i)
    HS_EXPECT_FALSE(bs.test(i, i + 1));
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs every mesh test case in order.
 * @return The module's failure count.
 */
inline int run_mesh_tests() {
  hs_test::ModuleFixture fixture("mesh");

  test_polymesh_default_state();
  test_polymesh_bind_arrays();
  test_polymesh_clear_resets_data();
  test_polymesh_unified_accessors();

  test_half_edge_mesh_size_matches_input();
  test_half_edge_mesh_face_loop_closes();
  test_half_edge_mesh_pairs_are_symmetric();
  test_half_edge_mesh_open_boundary_edges();
  test_half_edge_mesh_euler_invariant();
  test_half_edge_mesh_euler_tetrahedron();
  test_half_edge_mesh_built_from_meshstate();

  test_compile_polymesh_to_meshstate_basic();
  test_compile_drops_degenerate_faces();

  test_clone_meshstate_deep_copies();
  test_clone_polymesh_deep_copies();

  test_classify_faces_cube_uniform_topology();
  test_classify_faces_uncompiled_degenerate();
  test_degenerate_edge_records_never_pair();
  test_classify_faces_tetrahedron_uniform_topology();
  test_classify_faces_aliased_arenas();
  test_classify_faces_truncated_cube_distinct_topology();
  test_classify_faces_roster_hash_collision_free();

  test_tribitset_sizes();
  test_tribitset_clear_and_test();
  test_tribitset_test_and_set();
  test_tribitset_index_uniqueness();
  test_tribitset_all_pairs_independent();

  return fixture.result();
}

} // namespace mesh_tests
} // namespace hs_test
