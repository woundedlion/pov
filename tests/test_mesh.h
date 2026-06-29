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
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace mesh_tests {

inline uint8_t mesh_arena_a[256 * 1024];
inline uint8_t mesh_arena_b[256 * 1024];

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
 *        vertex per vertex, one face per face, and one half-edge per face index.
 */
inline void test_half_edge_mesh_size_matches_input() {
  Arena arena(mesh_arena_a, sizeof(mesh_arena_a));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  HalfEdgeMesh he(arena, cube);
  HS_EXPECT_EQ(he.vertices.size(), cube.vertices.size());
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
      if (curr == HE_NONE) break; // don't index half_edges with the sentinel
      HS_EXPECT_EQ(he.half_edges[curr].face, (uint16_t)fi);
      curr = he.half_edges[curr].next;
      steps++;
      if (steps > max_steps) break; // ring failed to close
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
  for (uint16_t i : idx) open.faces.push_back(i);

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

  int V = static_cast<int>(he.vertices.size());
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

  int V = static_cast<int>(he.vertices.size());
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
  MeshOps::compile(cube, ms, arena);

  HalfEdgeMesh he(arena, ms);
  HS_EXPECT_EQ(he.vertices.size(), ms.vertices.size());
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

  m.faces.push_back(0); m.faces.push_back(1); m.faces.push_back(2);
  m.faces.push_back(0); m.faces.push_back(3);
  m.faces.push_back(1); m.faces.push_back(3);

  MeshState ms;
  MeshOps::compile(m, ms, dst);
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
 */
inline void test_clone_meshstate_deep_copies() {
  Arena src_arena(mesh_arena_a, sizeof(mesh_arena_a));
  Arena dst_arena(mesh_arena_b, sizeof(mesh_arena_b));

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
  for (size_t i = 0; i < src.face_counts.size(); ++i)
    HS_EXPECT_EQ(dst.face_counts[i], src.face_counts[i]);
  for (size_t i = 0; i < src.faces.size(); ++i)
    HS_EXPECT_EQ(dst.faces[i], src.faces[i]);
  for (size_t i = 0; i < src.face_offsets.size(); ++i)
    HS_EXPECT_EQ(dst.face_offsets[i], src.face_offsets[i]);
}

/**
 * @brief Verifies clone() of a PolyMesh likewise copies all arrays into
 *        independent storage.
 */
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

/**
 * @brief A face topology key plus the base hash the classifier derives from it.
 */
struct FaceTopoRecord {
  int count;          /**< Side count. */
  int angles[24];     /**< Sorted whole-degree interior angles. */
  uint32_t hash;      /**< MeshOps base topology hash for this key. */
};

/**
 * @brief Recomputes a face's canonical key (count, sorted angles) and the base
 *        topology hash, mirroring classify_faces_impl.
 */
inline FaceTopoRecord face_topo_record(const PolyMesh &mesh, const uint16_t *idx,
                                       int count) {
  FaceTopoRecord rec;
  rec.count = count;
  for (int k = 0; k < 24; ++k)
    rec.angles[k] = 0;
  HS_EXPECT_TRUE(count <= 24);
  uint32_t h = 0x12345678;
  MeshOps::hash_combine(h, static_cast<uint32_t>(count));
  if (count >= 3) {
    for (int k = 0; k < count; ++k) {
      const Vector &prev = mesh.vertices[idx[(k - 1 + count) % count]];
      const Vector &curr = mesh.vertices[idx[k]];
      const Vector &next = mesh.vertices[idx[(k + 1) % count]];
      Vector v1 = (prev - curr).normalized();
      Vector v2 = (next - curr).normalized();
      rec.angles[k] =
          (int)std::round(angle_between(v1, v2) * 180.0f / PI_F);
    }
    std::sort(rec.angles, rec.angles + count);
    for (int k = 0; k < count; ++k)
      MeshOps::hash_combine(h, static_cast<uint32_t>(rec.angles[k]));
  }
  rec.hash = MeshOps::fmix32(h);
  return rec;
}

/**
 * @brief True when two records share the same canonical key (count + angles).
 */
inline bool same_topo_key(const FaceTopoRecord &a, const FaceTopoRecord &b) {
  if (a.count != b.count)
    return false;
  for (int k = 0; k < a.count; ++k)
    if (a.angles[k] != b.angles[k])
      return false;
  return true;
}

/**
 * @brief Verifies the per-face topology hash is collision-free across the
 *        Platonic/Archimedean/Catalan roster.
 * @details A fmix32 collision between two distinct face topologies would merge
 *          them into one palette class. Builds every roster solid, collects the
 *          distinct (count, sorted-angle) keys, and asserts each maps to a
 *          unique hash.
 */
inline void test_classify_faces_roster_hash_collision_free() {
  static FaceTopoRecord keys[256];
  int n = 0;

  auto collect = [&](const Solids::Entry *reg, size_t reg_n) {
    for (size_t r = 0; r < reg_n; ++r) {
      Arena a(mesh_arena_a, sizeof(mesh_arena_a));
      Arena b(mesh_arena_b, sizeof(mesh_arena_b));
      PolyMesh mesh = reg[r].generate(a, b);

      const uint8_t *fc = mesh.get_face_counts_data();
      const uint16_t *faces = mesh.get_faces_data();
      size_t F = mesh.get_face_counts_size();
      size_t off = 0;
      for (size_t f = 0; f < F; ++f) {
        int count = fc[f];
        FaceTopoRecord rec = face_topo_record(mesh, faces + off, count);
        off += count;

        bool seen = false;
        for (int i = 0; i < n; ++i) {
          if (same_topo_key(keys[i], rec)) {
            seen = true;
            break;
          }
          // Distinct keys must not share a hash.
          HS_EXPECT_TRUE(keys[i].hash != rec.hash);
        }
        if (!seen) {
          HS_EXPECT_TRUE(n < (int)(sizeof(keys) / sizeof(keys[0])));
          keys[n++] = rec;
        }
      }
    }
  };

  collect(Solids::simple_registry,
          sizeof(Solids::simple_registry) / sizeof(Solids::simple_registry[0]));
  collect(Solids::catalan_registry,
          sizeof(Solids::catalan_registry) /
              sizeof(Solids::catalan_registry[0]));

  HS_EXPECT_TRUE(n > 0);
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
  test_classify_faces_tetrahedron_uniform_topology();
  test_classify_faces_truncated_cube_distinct_topology();
  test_classify_faces_roster_hash_collision_free();

  return fixture.result();
}

} // namespace mesh_tests
} // namespace hs_test

