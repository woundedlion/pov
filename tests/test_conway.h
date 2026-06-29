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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <vector>
#include "core/conway.h"
#include "core/solids.h"
#include "tests/mesh_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace conway_tests {

/** @brief Scratch arena for the Conway operator's primary output mesh. */
inline uint8_t conway_target_buf[256 * 1024]; /**< Comfortable budget for cube-scale ops. */
/** @brief Scratch arena for the Conway operator's temporary/intermediate mesh. */
inline uint8_t conway_temp_buf[256 * 1024]; /**< Comfortable budget for cube-scale ops. */

// build_solid(), check_all_unit_vertices(), check_face_counts_consistent(),
// and check_indices_in_range() live in tests/mesh_test_util.h.

// ---------------------------------------------------------------------------
// Structural invariants
// ---------------------------------------------------------------------------

/**
 * @brief Computes a face normal via Newell's method.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Unnormalised normal vector for the face.
 * @details Newell's method is robust for non-planar faces (e.g. curved faces
 *          on the unit sphere) where a simple cross product would be ambiguous.
 */
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

/**
 * @brief Computes the unweighted centroid of a face's vertex positions.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Arithmetic mean of the face's vertex positions.
 */
inline Vector face_centroid_pos(const PolyMesh &m, size_t face_idx_offset,
                                int count) {
  Vector c(0, 0, 0);
  for (int k = 0; k < count; ++k)
    c = c + m.vertices[m.faces[face_idx_offset + k]];
  return c * (1.0f / static_cast<float>(count));
}

/**
 * @brief Counts faces whose winding produces an inward-pointing normal.
 * @param m Mesh to inspect; normals are taken relative to the mesh origin.
 * @return Number of faces whose normal disagrees with the centroid direction.
 * @details On a unit-sphere mesh, correct CCW-from-outside winding gives a
 *          normal that agrees in direction with the face centroid; a face
 *          failing this is wound backwards. Degenerate normals are ignored.
 */
inline int count_inward_winding(const PolyMesh &m) {
  int bad = 0;
  size_t offset = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    int count = m.face_counts[fi];
    if (count >= 3) {
      Vector n = face_newell_normal(m, offset, count);
      Vector c = face_centroid_pos(m, offset, count);
      // Skip degenerate normals (colinear vertices).
      if (n.length() > 1e-6f && dot(n, c) <= 0.0f) ++bad;
    }
    offset += count;
  }
  return bad;
}

/**
 * @brief Counts directed edges reused in the same direction by two faces.
 * @param m Mesh to inspect.
 * @return Number of duplicate same-direction directed edges.
 * @details In a manifold mesh with consistent winding, each undirected edge
 *          appears in two faces with OPPOSITE orientations. Same-direction
 *          reuse means at least one of the two faces is wound backwards.
 */
inline int count_same_direction_edge_violations(const PolyMesh &m) {
  /** @brief A single directed edge from source vertex u to destination v. */
  struct DirEdge { uint16_t u, v; };
  size_t total_edges = 0;
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    total_edges += m.face_counts[i];

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

  std::sort(edges.begin(), edges.end(),
            [](const DirEdge &a, const DirEdge &b) {
              return a.u != b.u ? a.u < b.u : a.v < b.v;
            });

  int dup = 0;
  for (size_t i = 1; i < edges.size(); ++i) {
    if (edges[i].u == edges[i - 1].u && edges[i].v == edges[i - 1].v) ++dup;
  }
  return dup;
}

/**
 * @brief Asserts the mesh is wound consistently.
 * @param m Mesh to validate.
 * @details Requires every face to point outward and no directed edge to be
 *          reused in the same direction; reports counts before each assertion.
 */
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

/**
 * @brief Asserts the structural invariants shared by every Conway-op test.
 * @param m Mesh to validate.
 * @details Checks non-empty arrays, internally consistent face_counts/faces,
 *          in-range indices, vertices on the unit sphere, and consistent
 *          winding.
 */
inline void check_basic_invariants(const PolyMesh &m) {
  HS_EXPECT_TRUE(m.vertices.size() > 0);
  HS_EXPECT_TRUE(m.face_counts.size() > 0);
  HS_EXPECT_TRUE(m.faces.size() > 0);
  check_face_counts_consistent(m);
  check_indices_in_range(m);
  check_all_unit_vertices(m, 1e-3f);
  check_consistent_winding(m);
}

/**
 * @brief Asserts the mesh is a closed genus-0 2-manifold (V - E + F == 2).
 * @param m Mesh to validate.
 * @details Conway operators map a genus-0 seed to another genus-0 closed
 *          polyhedron, so every undirected edge must be shared by exactly two
 *          faces. A boundary edge (shared once) or a non-manifold edge
 *          (shared >2) breaks the topology the renderer assumes; this catches
 *          both, then verifies Euler's formula.
 */
inline void check_euler_genus0(const PolyMesh &m) {
  std::vector<std::pair<uint16_t, uint16_t>> edges;
  edges.reserve(m.faces.size());
  size_t offset = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    int count = m.face_counts[fi];
    for (int k = 0; k < count; ++k) {
      uint16_t u = m.faces[offset + k];
      uint16_t v = m.faces[offset + (k + 1) % count];
      if (u > v) std::swap(u, v);
      edges.push_back({u, v});
    }
    offset += count;
  }
  std::sort(edges.begin(), edges.end());

  int E = 0;
  for (size_t i = 0; i < edges.size();) {
    size_t j = i;
    while (j < edges.size() && edges[j] == edges[i]) ++j;
    HS_EXPECT_EQ((int)(j - i), 2); // each edge bounds exactly two faces
    ++E;
    i = j;
  }

  int V = (int)m.vertices.size();
  int F = (int)m.face_counts.size();
  HS_EXPECT_EQ(V - E + F, 2);
}

/**
 * @brief Histogram of vertex degree (incident faces) → number of such vertices.
 * @param m Mesh to inspect.
 * @details In a closed manifold a vertex's incident-face count equals its edge
 *          degree, so this pins the local connectivity a bare vertex/face count
 *          cannot. Every vertex is tallied, so an unreferenced (degree-0) vertex
 *          shows up as a {0, n} bucket.
 */
inline std::map<int, int> vertex_degree_histogram(const PolyMesh &m) {
  std::vector<int> degree(m.vertices.size(), 0);
  for (uint16_t idx : m.faces)
    ++degree[idx];
  std::map<int, int> hist;
  for (int d : degree)
    ++hist[d];
  return hist;
}

/**
 * @brief Histogram of face side-count (triangle, quad, ...) → number of faces.
 * @param m Mesh to inspect.
 */
inline std::map<int, int> face_type_histogram(const PolyMesh &m) {
  std::map<int, int> hist;
  for (uint8_t c : m.face_counts)
    ++hist[(int)c];
  return hist;
}

/**
 * @brief Runs every primitive Conway operator on one seed and checks Euler.
 * @tparam Solid Platonic seed solid (e.g. Solids::Cube) to build and operate on.
 * @details Asserts each operator's result is a closed genus-0 manifold.
 *          Templated on the seed so one body covers triangle-, quad-, and
 *          pentagon-faced Platonic solids. Each op gets a fresh target/temp
 *          pair with the seed rebuilt into temp, mirroring the per-op tests.
 */
template <typename Solid>
inline void check_euler_for_seed() {
#define HS_EULER_OP(CALL)                                                      \
  do {                                                                         \
    Arena target(conway_target_buf, sizeof(conway_target_buf));               \
    Arena temp(conway_temp_buf, sizeof(conway_temp_buf));                      \
    PolyMesh seed;                                                             \
    build_solid<Solid>(seed, temp);                                           \
    check_euler_genus0(MeshOps::CALL);                                        \
  } while (0)

  {
    Arena arena(conway_target_buf, sizeof(conway_target_buf));
    PolyMesh seed;
    build_solid<Solid>(seed, arena);
    check_euler_genus0(seed);
  }
  HS_EULER_OP(dual(seed, target, temp));
  HS_EULER_OP(kis(seed, target, temp));
  HS_EULER_OP(ambo(seed, target, temp));
  HS_EULER_OP(truncate(seed, target, temp));
  HS_EULER_OP(expand(seed, target, temp));
  HS_EULER_OP(chamfer(seed, target, temp));
  HS_EULER_OP(snub(seed, target, temp));
  HS_EULER_OP(gyro(seed, target, temp));
  HS_EULER_OP(meta(seed, target, temp));
  HS_EULER_OP(needle(seed, target, temp));
  HS_EULER_OP(zip(seed, target, temp));
  HS_EULER_OP(bevel(seed, target, temp));
#undef HS_EULER_OP
}

/**
 * @brief Verifies every Conway operator preserves the Euler characteristic.
 * @details Exercises check_euler_for_seed across triangle-, quad-, and
 *          pentagon-faced Platonic seeds.
 */
inline void test_conway_ops_preserve_euler_characteristic() {
  check_euler_for_seed<Solids::Tetrahedron>();
  check_euler_for_seed<Solids::Cube>();
  check_euler_for_seed<Solids::Icosahedron>();
  check_euler_for_seed<Solids::Dodecahedron>();
}

// ---------------------------------------------------------------------------
// Input fixture: the unit-sphere cube.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies the cube seed fixture is well formed and on the unit sphere.
 * @details Checks 8 vertices, 6 quad faces, 24 indices, consistent face_counts,
 *          in-range indices, unit-sphere magnitudes, and consistent winding.
 */
inline void test_input_cube_is_well_formed() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  HS_EXPECT_EQ(cube.vertices.size(), (size_t)8);
  HS_EXPECT_EQ(cube.face_counts.size(), (size_t)6);
  HS_EXPECT_EQ(cube.faces.size(), (size_t)24);
  check_face_counts_consistent(cube);
  check_indices_in_range(cube);
  check_all_unit_vertices(cube, 1e-3f);
  check_consistent_winding(cube);
}

/** @brief Verifies the tetrahedron seed fixture is wound consistently. */
inline void test_input_tetrahedron_winding() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh tet;
  build_solid<Solids::Tetrahedron>(tet, arena);
  check_consistent_winding(tet);
}

// ---------------------------------------------------------------------------
// normalize — exposed helper inside MeshOps namespace.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies MeshOps::normalize projects all vertices onto the unit sphere.
 * @details Builds a mesh with off-sphere vertices (lengths 5, 7, √12) and
 *          asserts they land at |v|=1 after normalization.
 */
inline void test_normalize_pushes_to_unit_sphere() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh m;
  m.vertices.bind(arena, 3);
  m.face_counts.bind(arena, 1);
  m.faces.bind(arena, 3);
  m.vertices.push_back(Vector(3, 4, 0));     // length 5
  m.vertices.push_back(Vector(0, 0, 7));     // length 7
  m.vertices.push_back(Vector(2, 2, 2));     // length √12
  m.face_counts.push_back(3);
  m.faces.push_back(0);
  m.faces.push_back(1);
  m.faces.push_back(2);

  MeshOps::normalize(m);
  check_all_unit_vertices(m, 1e-3f);
}

// ---------------------------------------------------------------------------
// dual(cube)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies dual(cube) yields octahedral topology.
 * @details Each cube face becomes a dual vertex (6) and each cube vertex
 *          becomes a triangular dual face (8 faces, 24 indices).
 */
inline void test_dual_cube_has_octahedral_topology() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh d = MeshOps::dual(cube, target, temp);

  check_basic_invariants(d);
  HS_EXPECT_EQ(d.vertices.size(), (size_t)6);
  HS_EXPECT_EQ(d.face_counts.size(), (size_t)8);
  HS_EXPECT_EQ(d.faces.size(), (size_t)24);

  // Octahedral connectivity: all 6 vertices have degree 4, all 8 faces triangles.
  auto deg = vertex_degree_histogram(d);
  HS_EXPECT_EQ(deg.size(), (size_t)1);
  HS_EXPECT_EQ(deg[4], 6);
  auto faces = face_type_histogram(d);
  HS_EXPECT_EQ(faces.size(), (size_t)1);
  HS_EXPECT_EQ(faces[3], 8);
}

// ---------------------------------------------------------------------------
// kis(cube)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies kis(cube) pyramidalizes every face into triangles.
 * @details Adds one center vertex per face (8 + 6 vertices) and turns each
 *          quad into 4 triangles (24 triangular faces, 72 indices).
 */
inline void test_kis_cube_pyramidalizes() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh k = MeshOps::kis(cube, target, temp);

  check_basic_invariants(k);
  HS_EXPECT_EQ(k.vertices.size(), (size_t)(8 + 6)); // original V + one per F
  HS_EXPECT_EQ(k.face_counts.size(), (size_t)24);
  HS_EXPECT_EQ(k.faces.size(), (size_t)72);

  // kis lifts each original corner to degree 6 (its 3 faces each split in two)
  // and adds a degree-4 apex per face; every face is a triangle.
  auto deg = vertex_degree_histogram(k);
  HS_EXPECT_EQ(deg.size(), (size_t)2);
  HS_EXPECT_EQ(deg[6], 8);
  HS_EXPECT_EQ(deg[4], 6);
  auto faces = face_type_histogram(k);
  HS_EXPECT_EQ(faces.size(), (size_t)1);
  HS_EXPECT_EQ(faces[3], 24);
}

// ---------------------------------------------------------------------------
// ambo(cube)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies ambo(cube) yields cuboctahedral topology.
 * @details Each cube edge becomes a vertex (12) and the result has 6 squares
 *          plus 8 triangles (14 faces).
 */
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
// truncate(cube)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies truncate(cube) cuts each corner into the expected topology.
 * @details Yields 24 vertices (2 per edge) and 14 faces (6 octagons + 8
 *          triangles, 72 indices).
 */
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

/**
 * @brief Verifies truncate(t = 0.5) reproduces ambo() topology.
 * @details At t = 0.5 the truncation midpoints coincide, collapsing to 12
 *          vertices and 14 faces — the cuboctahedron ambo() also produces.
 */
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
// expand(cube)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies expand(cube) yields rhombicuboctahedral topology.
 * @details One vertex per (face, vertex) pair (24) and 26 faces: 6 shrunken
 *          squares, 12 edge quads, and 8 vertex triangles.
 */
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
// chamfer(cube)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies chamfer(cube) replaces each edge with a hexagon.
 * @details Yields V + 2E = 8 + 24 = 32 vertices and F + E = 6 + 12 = 18 faces.
 */
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
// relax
// ---------------------------------------------------------------------------

/**
 * @brief Verifies relax preserves topology and keeps vertices on the sphere.
 * @details Spring-based edge-length relaxation must not change vertex/face/
 *          index counts and must leave every vertex at |v|=1 with consistent
 *          winding.
 */
inline void test_relax_preserves_topology() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh c = MeshOps::relax(cube, target, temp, /*iterations*/ 3);

  HS_EXPECT_EQ(c.vertices.size(), cube.vertices.size());
  HS_EXPECT_EQ(c.face_counts.size(), cube.face_counts.size());
  HS_EXPECT_EQ(c.faces.size(), cube.faces.size());
  check_all_unit_vertices(c, 1e-3f);
  check_consistent_winding(c);
}

/**
 * @brief Population variance of a mesh's edge lengths.
 * @details Sums every per-face edge (each shared edge of a closed mesh is
 *          counted exactly twice, so the uniform double-count leaves the variance
 *          unchanged). Doubles accumulate to keep the two-pass mean/variance
 *          stable. This is the quantity relax() drives down by pulling every edge
 *          toward the mesh-wide mean length.
 */
template <typename MeshT>
inline float edge_length_variance(const MeshT &m) {
  const auto *fc = m.get_face_counts_data();
  const auto *faces = m.get_faces_data();
  const size_t F = m.get_face_counts_size();

  double sum = 0.0;
  int n = 0;
  size_t idx = 0;
  for (size_t f = 0; f < F; ++f) {
    const int c = fc[f];
    for (int k = 0; k < c; ++k) {
      const int a = faces[idx + k];
      const int b = faces[idx + (k + 1) % c];
      sum += distance_between(m.vertices[a], m.vertices[b]);
      ++n;
    }
    idx += c;
  }
  const float mean = static_cast<float>(sum / n);

  double var = 0.0;
  idx = 0;
  for (size_t f = 0; f < F; ++f) {
    const int c = fc[f];
    for (int k = 0; k < c; ++k) {
      const int a = faces[idx + k];
      const int b = faces[idx + (k + 1) % c];
      const float d = distance_between(m.vertices[a], m.vertices[b]) - mean;
      var += static_cast<double>(d) * d;
    }
    idx += c;
  }
  return static_cast<float>(var / n);
}

/**
 * @brief Verifies relax actually reduces edge-length variance — its purpose.
 * @details A cube has uniform edges, so it cannot show relax working. Truncating
 *          a cube yields corner triangles and face octagons whose edges differ in
 *          length; relax's spring system pulls every edge toward the mesh mean, so
 *          the post-relax edge-length variance must be strictly below the input's.
 */
inline void test_relax_reduces_edge_variance() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  // Build the uneven-edge source in `target`; the truncation scratch in `temp`
  // is freed once the truncated output (in `target`) is built.
  PolyMesh uneven;
  {
    ScratchScope tscope(temp);
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, temp);
    uneven = MeshOps::truncate(cube, target, temp);
  }

  const float var_before = edge_length_variance(uneven);
  HS_EXPECT_GT(var_before, 1e-5f);

  // `target` still holds `uneven`, so relax writes its output and scratch into
  // the two halves of the now-free `temp` buffer.
  Arena relax_out(conway_temp_buf, sizeof(conway_temp_buf) / 2);
  Arena relax_scratch(conway_temp_buf + sizeof(conway_temp_buf) / 2,
                      sizeof(conway_temp_buf) / 2);
  PolyMesh relaxed = MeshOps::relax(uneven, relax_out, relax_scratch,
                                    /*iterations*/ 12);

  const float var_after = edge_length_variance(relaxed);
  HS_EXPECT_LT(var_after, var_before);
}

/**
 * @brief Verifies relax runs on an open (boundary) mesh without tripping a trap
 *        and leaves every output vertex finite.
 * @details Two triangles sharing one edge: the outer edges are boundary edges
 *          (shared once), so the per-vertex orbit must fall back to the
 *          boundary-tolerant path rather than HS_CHECK-trapping.
 */
inline void test_relax_open_mesh_finite() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh strip;
  strip.vertices.bind(target, 4);
  strip.vertices.push_back(Vector(1, 0, 0).normalized());
  strip.vertices.push_back(Vector(0, 1, 0).normalized());
  strip.vertices.push_back(Vector(0, 0, 1).normalized());
  strip.vertices.push_back(Vector(1, 1, 0).normalized());
  strip.face_counts.bind(target, 2);
  strip.face_counts.push_back(3);
  strip.face_counts.push_back(3);
  strip.faces.bind(target, 6);
  strip.faces.push_back(0);
  strip.faces.push_back(1);
  strip.faces.push_back(2);
  strip.faces.push_back(0);
  strip.faces.push_back(3);
  strip.faces.push_back(1);

  PolyMesh r = MeshOps::relax(strip, target, temp, /*iterations*/ 5);

  HS_EXPECT_EQ(r.vertices.size(), strip.vertices.size());
  for (size_t i = 0; i < r.vertices.size(); ++i) {
    HS_EXPECT_TRUE(std::isfinite(r.vertices[i].x));
    HS_EXPECT_TRUE(std::isfinite(r.vertices[i].y));
    HS_EXPECT_TRUE(std::isfinite(r.vertices[i].z));
  }
}

// ---------------------------------------------------------------------------
// Compositional + standalone operators — structural invariants only.
// ---------------------------------------------------------------------------

/** @brief Verifies meta(cube) satisfies the basic structural invariants. */
inline void test_meta_cube() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);
  PolyMesh m = MeshOps::meta(cube, target, temp);
  check_basic_invariants(m);
}

/** @brief Verifies needle(cube) satisfies the basic structural invariants. */
inline void test_needle_cube() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);
  PolyMesh n = MeshOps::needle(cube, target, temp);
  check_basic_invariants(n);
}

/** @brief Verifies zip(cube) satisfies the basic structural invariants. */
inline void test_zip_cube() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);
  PolyMesh z = MeshOps::zip(cube, target, temp);
  check_basic_invariants(z);
}

/** @brief Verifies bevel(cube) satisfies the basic structural invariants. */
inline void test_bevel_cube() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);
  PolyMesh b = MeshOps::bevel(cube, target, temp);
  check_basic_invariants(b);
}

// ---------------------------------------------------------------------------
// Composition polarity: a primitive operator returns its output in `target`; a
// composed operator (gyro/meta/needle/zip/bevel) returns it in `temp`.
// ---------------------------------------------------------------------------

/**
 * @brief Reports whether a pointer lands inside a byte buffer's address range.
 * @param p Pointer to test (e.g. the address of an output mesh's first vertex).
 * @param buf Base of the arena's backing byte buffer.
 * @param n Size of the buffer in bytes.
 * @return True iff p lies within [buf, buf + n).
 * @details Compares via uintptr_t so the cross-object pointer test is fully
 *          defined (no unspecified relational pointer comparison under UBSan).
 */
inline bool ptr_in_buffer(const void *p, const uint8_t *buf, size_t n) {
  auto addr = reinterpret_cast<uintptr_t>(p);
  auto base = reinterpret_cast<uintptr_t>(buf);
  return addr >= base && addr < base + n;
}

/**
 * @brief Verifies the primitive-vs-composed arena polarity in conway.h.
 * @details A primitive (dual) must return its output in `target`; each composed
 *          operator (gyro/meta/needle/zip/bevel) must return its output in
 *          `temp`. The seed is built in `temp` for every case, matching the
 *          per-operator tests and SolidBuilder's calling convention. Each output
 *          also satisfies the basic structural invariants, so the polarity check
 *          is asserted on a genuinely valid mesh.
 */
inline void test_conway_composition_polarity() {
  const uint8_t *tgt = conway_target_buf;
  const uint8_t *tmp = conway_temp_buf;
  const size_t tgt_n = sizeof(conway_target_buf);
  const size_t tmp_n = sizeof(conway_temp_buf);

  // A primitive returns in `target`.
  {
    Arena target(conway_target_buf, sizeof(conway_target_buf));
    Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, temp);
    PolyMesh d = MeshOps::dual(cube, target, temp);
    check_basic_invariants(d);
    HS_EXPECT_TRUE(ptr_in_buffer(&d.vertices[0], tgt, tgt_n));
    HS_EXPECT_FALSE(ptr_in_buffer(&d.vertices[0], tmp, tmp_n));
  }

  // Each composed operator returns in `temp`.
#define HS_POLARITY_COMPOSED(CALL)                                             \
  do {                                                                         \
    Arena target(conway_target_buf, sizeof(conway_target_buf));               \
    Arena temp(conway_temp_buf, sizeof(conway_temp_buf));                      \
    PolyMesh seed;                                                             \
    build_solid<Solids::Cube>(seed, temp);                                     \
    PolyMesh out = MeshOps::CALL;                                              \
    check_basic_invariants(out);                                              \
    HS_EXPECT_TRUE(ptr_in_buffer(&out.vertices[0], tmp, tmp_n));              \
    HS_EXPECT_FALSE(ptr_in_buffer(&out.vertices[0], tgt, tgt_n));             \
  } while (0)

  HS_POLARITY_COMPOSED(gyro(seed, target, temp));
  HS_POLARITY_COMPOSED(meta(seed, target, temp));
  HS_POLARITY_COMPOSED(needle(seed, target, temp));
  HS_POLARITY_COMPOSED(zip(seed, target, temp));
  HS_POLARITY_COMPOSED(bevel(seed, target, temp));
#undef HS_POLARITY_COMPOSED
}

// ---------------------------------------------------------------------------
// snub(cube)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies snub(cube) yields snub-cube topology with the expected counts.
 * @details One inset vertex per half-edge (I = 24); faces are the 6 twisted
 *          quad primaries, 8 vertex-orbit triangles, and 2 triangles per edge
 *          (2 * E = 24), for 38 faces and 6*4 + 32*3 = 120 indices. The chiral
 *          snub leaves the primaries as quads (no triangulation at twist = 0).
 */
inline void test_snub_cube_is_well_formed() {
  Arena target(conway_target_buf, sizeof(conway_target_buf));
  Arena temp(conway_temp_buf, sizeof(conway_temp_buf));

  PolyMesh cube;
  build_solid<Solids::Cube>(cube, temp);

  PolyMesh s = MeshOps::snub(cube, target, temp);

  check_basic_invariants(s);
  HS_EXPECT_EQ(s.vertices.size(), (size_t)24);             // one per half-edge
  HS_EXPECT_EQ(s.face_counts.size(), (size_t)38);          // 6 + 8 + 2*12
  HS_EXPECT_EQ(s.faces.size(), (size_t)(6 * 4 + 32 * 3));  // 6 quads + 32 triangles

  auto faces = face_type_histogram(s);
  HS_EXPECT_EQ(faces.size(), (size_t)2);
  HS_EXPECT_EQ(faces[3], 32); // 8 vertex-orbit + 24 edge triangles
  HS_EXPECT_EQ(faces[4], 6);  // twisted square primaries
}

// ---------------------------------------------------------------------------
// transform — variadic vertex transformer pipeline.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies transform applies a variadic vertex-transformer chain.
 * @details Runs a scale-then-shift pipeline over a MeshState and checks each
 *          vertex equals v*2 + (1,0,0); also confirms the destination borrows
 *          the source topology via the unified accessors.
 */
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

  // Topology is borrowed (view) from the source; the dst bound no buffer.
  HS_EXPECT_EQ(dst.get_face_counts_size(), (size_t)1);
  HS_EXPECT_EQ(dst.get_faces_size(), (size_t)3);
}

/**
 * @brief Verifies transform resets stale owned topology when reusing a dst.
 * @details transform() yields a BORROWED-mode mesh (owned vertices, topology
 *          via *_view). Accessors discriminate on is_bound() and clear() does
 *          not unbind, so a destination reused from a prior owned-mode life
 *          carries still-bound owned topology that would shadow the freshly set
 *          views. This pins that transform() resets the owned topology, so the
 *          source views (sizes 1 / 3) win over the stale owned sizes (2 / 9).
 */
inline void test_transform_unbinds_stale_owned_topology_on_reuse() {
  Arena src_arena(conway_target_buf, sizeof(conway_target_buf) / 2);
  Arena dst_arena(conway_target_buf + sizeof(conway_target_buf) / 2,
                  sizeof(conway_target_buf) / 2);

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

  // Destination arrives in OWNED mode with stale topology (2 faces / 9 indices).
  MeshState dst;
  dst.face_counts.bind(dst_arena, 2);
  dst.face_counts.push_back(4);
  dst.face_counts.push_back(4);
  dst.faces.bind(dst_arena, 9);
  for (int i = 0; i < 9; ++i)
    dst.faces.push_back(9);

  MeshOps::transform(src, dst, dst_arena);

  // Accessors must reflect the source view, not the stale owned buffers.
  HS_EXPECT_EQ(dst.get_face_counts_size(), (size_t)1);
  HS_EXPECT_EQ(dst.get_faces_size(), (size_t)3);
  HS_EXPECT_EQ((int)dst.get_face_counts_data()[0], 3);
  HS_EXPECT_EQ((int)dst.get_faces_data()[0], 0);
  HS_EXPECT_EQ((int)dst.get_faces_data()[1], 1);
  HS_EXPECT_EQ((int)dst.get_faces_data()[2], 2);
}

// ---------------------------------------------------------------------------
// face_centroid + vertex_orbit (low-level half-edge helpers)
// ---------------------------------------------------------------------------

/**
 * @brief Verifies face_centroid reports a quad and a bounded centroid.
 * @details Builds a HalfEdgeMesh over the cube and checks face 0 has 4 sides
 *          and a centroid magnitude bounded by 1 (mean of unit vectors).
 */
inline void test_face_centroid_for_cube_top_face() {
  Arena arena(conway_target_buf, sizeof(conway_target_buf));
  PolyMesh cube;
  build_solid<Solids::Cube>(cube, arena);

  HalfEdgeMesh he(arena, cube);
  int count = 0;
  Vector c = MeshOps::face_centroid(he, cube, 0, count);
  HS_EXPECT_EQ(count, 4);
  // Centroid magnitude is bounded by 1 (mean of unit vectors).
  HS_EXPECT_TRUE(c.length() <= 1.0f + 1e-4f);
}

// ---------------------------------------------------------------------------
// Degenerate-face checks
// ---------------------------------------------------------------------------

/**
 * @brief Asserts every face in the mesh has at least 3 sides.
 * @param m Mesh to validate.
 */
inline void check_no_degenerate_faces(const PolyMesh &m) {
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    HS_EXPECT_TRUE(m.face_counts[i] >= 3);
}

/**
 * @brief Builds a minimal malformed digon: 2 vertices forming one 2-gon face.
 * @param m Mesh to populate (its arrays are bound from arena).
 * @param arena Arena providing backing storage for the mesh arrays.
 * @details HalfEdgeMesh accepts this (it has no <3-sided guard), so it reaches
 *          the primary-face loop of expand/chamfer/snub with count == 2.
 */
inline void build_degenerate_digon(PolyMesh &m, Arena &arena) {
  m.vertices.bind(arena, 2);
  m.face_counts.bind(arena, 1);
  m.faces.bind(arena, 2);
  m.vertices.push_back(Vector(1.0f, 0.0f, 0.0f));
  m.vertices.push_back(Vector(-1.0f, 0.0f, 0.0f));
  m.face_counts.push_back(2);
  m.faces.push_back(0);
  m.faces.push_back(1);
}

/**
 * @brief Verifies expand/chamfer/snub drop degenerate primary faces.
 * @details Feeding each operator a digon must not leak a sub-triangular primary
 *          face: each drops the degenerate primary face (matching ambo/truncate
 *          and the vertex-orbit emitter) so no output face has fewer than 3
 *          sides.
 */
inline void test_conway_ops_drop_degenerate_primary_faces() {
  {
    Arena target(conway_target_buf, sizeof(conway_target_buf));
    Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
    PolyMesh m;
    build_degenerate_digon(m, temp);
    check_no_degenerate_faces(MeshOps::expand(m, target, temp));
  }
  {
    Arena target(conway_target_buf, sizeof(conway_target_buf));
    Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
    PolyMesh m;
    build_degenerate_digon(m, temp);
    check_no_degenerate_faces(MeshOps::chamfer(m, target, temp, 0.5f));
  }
  {
    Arena target(conway_target_buf, sizeof(conway_target_buf));
    Arena temp(conway_temp_buf, sizeof(conway_temp_buf));
    PolyMesh m;
    build_degenerate_digon(m, temp);
    check_no_degenerate_faces(MeshOps::snub(m, target, temp));
  }
}

// ---------------------------------------------------------------------------

/**
 * @brief Runs every Conway test under a "conway" module scope.
 * @return Failure count reported by end_module for the "conway" module.
 */
inline int run_conway_tests() {
  hs_test::ModuleFixture fixture("conway");

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
  test_relax_preserves_topology();
  test_relax_reduces_edge_variance();
  test_relax_open_mesh_finite();
  test_meta_cube();
  test_needle_cube();
  test_zip_cube();
  test_bevel_cube();
  test_conway_composition_polarity();
  test_snub_cube_is_well_formed();
  test_conway_ops_preserve_euler_characteristic();
  test_conway_ops_drop_degenerate_primary_faces();
  test_transform_applies_translation_chain();
  test_transform_unbinds_stale_owned_topology_on_reuse();
  test_face_centroid_for_cube_top_face();

  return fixture.result();
}

} // namespace conway_tests
} // namespace hs_test

