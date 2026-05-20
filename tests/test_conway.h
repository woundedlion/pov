/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/conway.h.
 *
 * Builds a small input PolyMesh (a unit-sphere cube from Solids::Cube),
 * runs each Conway operator, and verifies structural invariants:
 *   - non-empty output
 *   - all vertices on (approximately) the unit sphere
 *   - face_counts/faces arrays internally consistent (Σ face_counts == |faces|)
 *   - all face indices reference valid vertices
 *   - operator-specific counts (e.g. |dual.vertices| == |input.faces|)
 */
#pragma once
#ifndef HOLOSPHERE_TESTS_TEST_CONWAY_H_
#define HOLOSPHERE_TESTS_TEST_CONWAY_H_

#include <algorithm>
#include <cstdint>
#include <vector>
#include "core/conway.h"
#include "core/solids.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_tests {

// Large scratch buffers — Conway ops on a cube need a comfortable budget.
inline uint8_t conway_target_buf[256 * 1024];
inline uint8_t conway_temp_buf[256 * 1024];

// Build a PolyMesh from a Solids::* descriptor into the given arena.
template <typename Solid>
inline void build_solid(PolyMesh &mesh, Arena &arena) {
  mesh.initialize(arena, Solid::NUM_VERTS, Solid::NUM_FACES,
                  Solid::faces.size());
  for (const auto &v : Solid::vertices) mesh.vertices.push_back(v);
  for (auto fc : Solid::face_counts) mesh.face_counts.push_back(fc);
  for (auto fi : Solid::faces)
    mesh.faces.push_back(static_cast<uint16_t>(fi));
}

// ---------------------------------------------------------------------------
// Structural invariants — used by every Conway-op test
// ---------------------------------------------------------------------------

inline void check_face_counts_consistent(const PolyMesh &m) {
  size_t total = 0;
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    total += m.face_counts[i];
  HS_EXPECT_EQ(total, m.faces.size());
}

inline void check_indices_in_range(const PolyMesh &m) {
  size_t V = m.vertices.size();
  for (size_t i = 0; i < m.faces.size(); ++i) {
    HS_EXPECT_TRUE(m.faces[i] < V);
  }
}

inline void check_all_unit_vertices(const PolyMesh &m) {
  for (size_t i = 0; i < m.vertices.size(); ++i) {
    float len = m.vertices[i].length();
    HS_EXPECT_NEAR(len, 1.0f, 1e-3f);
  }
}

// Newell's method face normal — robust for non-planar faces (curved faces
// on the unit sphere). Returns the unnormalised normal vector.
inline Vector face_newell_normal(const PolyMesh &m, size_t face_idx_offset,
                                 int count) {
  Vector n(0, 0, 0);
  for (int k = 0; k < count; ++k) {
    const Vector &curr = m.vertices[m.faces[face_idx_offset + k]];
    const Vector &next =
        m.vertices[m.faces[face_idx_offset + (k + 1) % count]];
    n.x += (curr.y - next.y) * (curr.z + next.z);
    n.y += (curr.z - next.z) * (curr.x + next.x);
    n.z += (curr.x - next.x) * (curr.y + next.y);
  }
  return n;
}

inline Vector face_centroid_pos(const PolyMesh &m, size_t face_idx_offset,
                                int count) {
  Vector c(0, 0, 0);
  for (int k = 0; k < count; ++k)
    c = c + m.vertices[m.faces[face_idx_offset + k]];
  return c * (1.0f / static_cast<float>(count));
}

// Returns the number of faces whose winding produces an INWARD-pointing
// normal (relative to the mesh origin). On a unit-sphere mesh, correct
// CCW-from-outside winding gives a normal that agrees in direction with
// the face centroid.
inline int count_inward_winding(const PolyMesh &m) {
  int bad = 0;
  size_t offset = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    int count = m.face_counts[fi];
    if (count >= 3) {
      Vector n = face_newell_normal(m, offset, count);
      Vector c = face_centroid_pos(m, offset, count);
      // Skip degenerate normals (colinear vertices) — caller can decide.
      if (n.length() > 1e-6f && dot(n, c) <= 0.0f) ++bad;
    }
    offset += count;
  }
  return bad;
}

// Verifies that no two faces share an oriented edge in the same direction.
// In a manifold mesh with consistent winding, each undirected edge appears
// in two faces with OPPOSITE orientations. Same-direction reuse means at
// least one of the two faces is wound backwards.
inline int count_same_direction_edge_violations(const PolyMesh &m) {
  // Collect every directed edge (u -> v) from every face.
  struct DirEdge { uint16_t u, v; };
  size_t total_edges = 0;
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    total_edges += m.face_counts[i];

  // Use stack arrays via std::vector replacement: ArenaVector requires arena.
  // We're test code — std::vector is fine here.
  std::vector<DirEdge> edges;
  edges.reserve(total_edges);

  size_t offset = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    int count = m.face_counts[fi];
    for (int k = 0; k < count; ++k) {
      uint16_t u = m.faces[offset + k];
      uint16_t v = m.faces[offset + (k + 1) % count];
      edges.push_back({u, v});
    }
    offset += count;
  }

  // Sort lexicographically by (u, v).
  std::sort(edges.begin(), edges.end(),
            [](const DirEdge &a, const DirEdge &b) {
              return a.u != b.u ? a.u < b.u : a.v < b.v;
            });

  // Count adjacent duplicates (same direction).
  int dup = 0;
  for (size_t i = 1; i < edges.size(); ++i) {
    if (edges[i].u == edges[i - 1].u && edges[i].v == edges[i - 1].v) ++dup;
  }
  return dup;
}

inline void check_consistent_winding(const PolyMesh &m) {
  int inward = count_inward_winding(m);
  if (inward != 0)
    std::printf("    [winding] %d / %zu faces have inward-facing normals\n",
                inward, m.face_counts.size());
  HS_EXPECT_EQ(inward, 0);

  int dups = count_same_direction_edge_violations(m);
  if (dups != 0)
    std::printf("    [winding] %d directed edges reused in same direction\n",
                dups);
  HS_EXPECT_EQ(dups, 0);
}

inline void check_basic_invariants(const PolyMesh &m) {
  HS_EXPECT_TRUE(m.vertices.size() > 0);
  HS_EXPECT_TRUE(m.face_counts.size() > 0);
  HS_EXPECT_TRUE(m.faces.size() > 0);
  check_face_counts_consistent(m);
  check_indices_in_range(m);
  check_all_unit_vertices(m);
  check_consistent_winding(m);
}

// ---------------------------------------------------------------------------
// Input fixture: the unit-sphere cube. Vertex magnitudes are 1/√3 each, so
// |v|=1 — confirms our cube data is correctly normalised.
// ---------------------------------------------------------------------------

inline void test_input_cube_is_well_formed() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  HS_EXPECT_EQ(cube.vertices.size(), (size_t)8);
  HS_EXPECT_EQ(cube.face_counts.size(), (size_t)6);
  HS_EXPECT_EQ(cube.faces.size(), (size_t)24);
  check_face_counts_consistent(cube);
  check_indices_in_range(cube);
  check_all_unit_vertices(cube);
  check_consistent_winding(cube);
}

inline void test_input_tetrahedron_winding() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh tet;
  build_solid<Solids::Tetrahedron>(tet, arena);
  check_consistent_winding(tet);
}

// ---------------------------------------------------------------------------
// normalize — exposed helper inside MeshOps namespace.
// ---------------------------------------------------------------------------

inline void test_normalize_pushes_to_unit_sphere() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh m;
  m.initialize(arena, 3, 1, 3);
  m.vertices.push_back(Vector(3, 4, 0));     // length 5
  m.vertices.push_back(Vector(0, 0, 7));     // length 7
  m.vertices.push_back(Vector(2, 2, 2));     // length √12
  m.face_counts.push_back(3);
  m.faces.push_back(0);
  m.faces.push_back(1);
  m.faces.push_back(2);

  MeshOps::normalize(m);
  check_all_unit_vertices(m);
}

// ---------------------------------------------------------------------------
// dual(cube) → octahedron-topology (6 vertices, 8 faces, 24 face indices)
// ---------------------------------------------------------------------------

inline void test_dual_cube_has_octahedral_topology() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh d = MeshOps::dual(cube, target, temp);

  check_basic_invariants(d);
  // Each original face → one dual vertex. Cube has 6 faces.
  HS_EXPECT_EQ(d.vertices.size(), (size_t)6);
  // Each original vertex → one dual face. Cube has 8 vertices → 8 faces,
  // each triangular (3 indices) for an octahedron.
  HS_EXPECT_EQ(d.face_counts.size(), (size_t)8);
  HS_EXPECT_EQ(d.faces.size(), (size_t)24);
}

// ---------------------------------------------------------------------------
// kis(cube) — pyramidalize each face. Adds 1 center vertex per face,
// each n-gon face becomes n triangles.
// ---------------------------------------------------------------------------

inline void test_kis_cube_pyramidalizes() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh k = MeshOps::kis(cube, target, temp);

  check_basic_invariants(k);
  HS_EXPECT_EQ(k.vertices.size(), (size_t)(8 + 6)); // original V + one per F
  // 6 quad faces × 4 triangles each = 24 triangles
  HS_EXPECT_EQ(k.face_counts.size(), (size_t)24);
  HS_EXPECT_EQ(k.faces.size(), (size_t)72);
  // Every face must be a triangle
  for (size_t i = 0; i < k.face_counts.size(); ++i)
    HS_EXPECT_EQ(k.face_counts[i], (uint8_t)3);
}

// ---------------------------------------------------------------------------
// ambo(cube) → cuboctahedron-topology. Each original edge → one new vertex.
// Cube has 12 edges → 12 vertices, 6+8 = 14 faces.
// ---------------------------------------------------------------------------

inline void test_ambo_cube_has_cuboctahedral_topology() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh a = MeshOps::ambo(cube, target, temp);

  check_basic_invariants(a);
  HS_EXPECT_EQ(a.vertices.size(), (size_t)12);
  HS_EXPECT_EQ(a.face_counts.size(), (size_t)14); // 6 squares + 8 triangles
}

// ---------------------------------------------------------------------------
// truncate(cube) — cuts each corner. Result: 24 vertices (2 per edge),
// 6 octagons + 8 triangles = 14 faces.
// ---------------------------------------------------------------------------

inline void test_truncate_cube_has_truncated_topology() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh tr = MeshOps::truncate(cube, target, temp, 0.25f);

  check_basic_invariants(tr);
  HS_EXPECT_EQ(tr.vertices.size(), (size_t)24);
  HS_EXPECT_EQ(tr.face_counts.size(), (size_t)14);

  // 6 octagonal faces (8 sides) + 8 triangles = 6*8 + 8*3 = 72 indices
  HS_EXPECT_EQ(tr.faces.size(), (size_t)(6 * 8 + 8 * 3));
}

// truncate(t = 0.5) is equivalent to ambo().
inline void test_truncate_t_half_is_ambo() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube1;
  build_solid<Solids::Cube>(cube1, temp);
  PolyMesh tr = MeshOps::truncate(cube1, target, temp, 0.5f);

  // Match topology of ambo()
  HS_EXPECT_EQ(tr.vertices.size(), (size_t)12);
  HS_EXPECT_EQ(tr.face_counts.size(), (size_t)14);
}

// ---------------------------------------------------------------------------
// expand(cube) — rhombicuboctahedron-topology. Each face shrinks, each edge
// becomes a quad, each vertex becomes a triangle.
// ---------------------------------------------------------------------------

inline void test_expand_cube() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh e = MeshOps::expand(cube, target, temp);

  check_basic_invariants(e);
  // Vertices: one per (face, vertex) pair → 24 for the cube.
  HS_EXPECT_EQ(e.vertices.size(), (size_t)24);
  // Faces: 6 original (shrunken squares) + 12 edge quads + 8 vertex triangles
  HS_EXPECT_EQ(e.face_counts.size(), (size_t)(6 + 12 + 8));
}

// ---------------------------------------------------------------------------
// chamfer(cube) — replaces each edge with a hexagon. Result: V + 2E vertices
// = 8 + 24 = 32, F + E faces = 6 + 12 = 18.
// ---------------------------------------------------------------------------

inline void test_chamfer_cube() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh c = MeshOps::chamfer(cube, target, temp, 0.5f);

  check_basic_invariants(c);
  HS_EXPECT_EQ(c.vertices.size(), (size_t)(8 + 2 * 12));
  HS_EXPECT_EQ(c.face_counts.size(), (size_t)(6 + 12));
}

// ---------------------------------------------------------------------------
// canonicalize — relaxation, output must keep the same topology and all
// vertices remain on the unit sphere.
// ---------------------------------------------------------------------------

inline void test_canonicalize_preserves_topology() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh c = MeshOps::canonicalize(cube, target, temp, /*iterations*/ 3);

  HS_EXPECT_EQ(c.vertices.size(), cube.vertices.size());
  HS_EXPECT_EQ(c.face_counts.size(), cube.face_counts.size());
  HS_EXPECT_EQ(c.faces.size(), cube.faces.size());
  check_all_unit_vertices(c);
  check_consistent_winding(c);
}

// ---------------------------------------------------------------------------
// snub(cube) — chiral operator; structural sanity only.
// ---------------------------------------------------------------------------

inline void test_snub_cube_is_well_formed() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh s = MeshOps::snub(cube, target, temp);

  check_basic_invariants(s);
}

// ---------------------------------------------------------------------------
// transform — variadic vertex transformer pipeline. Operates on MeshState
// (which has the face_counts_view / faces_view borrowed members).
// ---------------------------------------------------------------------------

inline void test_transform_applies_translation_chain() {
  Arena src_arena(conway_target_buf, sizeof(conway_target_buf) / 2);
  Arena dst_arena(conway_target_buf + sizeof(conway_target_buf) / 2,
                  sizeof(conway_target_buf) / 2);

  MeshState src;
  src.vertices.bind(src_arena, 4);
  src.vertices.push_back(Vector(1, 0, 0));
  src.vertices.push_back(Vector(0, 1, 0));
  src.vertices.push_back(Vector(0, 0, 1));
  src.vertices.push_back(Vector(0.577f, 0.577f, 0.577f));
  src.face_counts.bind(src_arena, 1);
  src.face_counts.push_back(3);
  src.faces.bind(src_arena, 3);
  src.faces.push_back(0);
  src.faces.push_back(1);
  src.faces.push_back(2);

  MeshState dst;
  auto scale = [](const Vector &v) { return v * 2.0f; };
  auto shift = [](const Vector &v) { return v + Vector(1, 0, 0); };
  MeshOps::transform(src, dst, dst_arena, scale, shift);

  HS_EXPECT_EQ(dst.vertices.size(), src.vertices.size());
  for (size_t i = 0; i < src.vertices.size(); ++i) {
    Vector expected = src.vertices[i] * 2.0f + Vector(1, 0, 0);
    HS_EXPECT_NEAR(dst.vertices[i].x, expected.x, 1e-5f);
    HS_EXPECT_NEAR(dst.vertices[i].y, expected.y, 1e-5f);
    HS_EXPECT_NEAR(dst.vertices[i].z, expected.z, 1e-5f);
  }

  // Topology is borrowed (view) from the source — accessible via unified
  // accessors even though the destination didn't bind its own buffer.
  HS_EXPECT_EQ(dst.get_face_counts_size(), (size_t)1);
  HS_EXPECT_EQ(dst.get_faces_size(), (size_t)3);
}

// ---------------------------------------------------------------------------
// face_centroid + vertex_orbit (low-level half-edge helpers)
// ---------------------------------------------------------------------------

inline void test_face_centroid_for_cube_top_face() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  HalfEdgeMesh he(arena, cube);
  int count = 0;
  Vector c = MeshOps::face_centroid(he, cube, 0, count);
  HS_EXPECT_EQ(count, 4); // cube face is a quad
  // Centroid magnitude is bounded by 1 (mean of unit vectors)
  HS_EXPECT_TRUE(c.length() <= 1.0f + 1e-4f);
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

inline int run_conway_tests() {
  auto scope = hs_test::begin_module("conway");

  test_input_cube_is_well_formed();
  test_input_tetrahedron_winding();
  test_normalize_pushes_to_unit_sphere();
  test_dual_cube_has_octahedral_topology();
  test_kis_cube_pyramidalizes();
  test_ambo_cube_has_cuboctahedral_topology();
  test_truncate_cube_has_truncated_topology();
  test_truncate_t_half_is_ambo();
  test_expand_cube();
  test_chamfer_cube();
  test_canonicalize_preserves_topology();
  test_snub_cube_is_well_formed();
  test_transform_applies_translation_chain();
  test_face_centroid_for_cube_top_face();

  return hs_test::end_module(scope);
}

} // namespace conway_tests
} // namespace hs_test

#endif // HOLOSPHERE_TESTS_TEST_CONWAY_H_
