/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared mesh fixtures, vertex checks, morph recipes, probes, and sliver gates.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <vector>
#include "core/mesh/mesh.h"
#include "core/mesh/solids.h"
#include "tests/test_harness.h"

namespace hs_test {

/** The identity truncate depth the smooth path sweeps to (the "uniform" Conway
 * depth). */
inline constexpr float RECONCILE_TRUNCATE_T = 1.0f / 3.0f;

/**
 * @brief Interior angle at corner `k` of a face, in degrees, before rounding.
 * @param m Mesh owning the vertices.
 * @param idx Start of the face's index run.
 * @param count Sides in the face.
 * @param k Corner to measure.
 * @return The value classify_faces_impl rounds to a whole degree.
 * @details Mirrors classify_faces_impl (core/mesh/mesh.h) term for term:
 *          unnormalized edges, the same degenerate-edge cutoff, acosf.
 */
inline float classifier_angle_deg(const PolyMesh &m, const uint16_t *idx,
                                  int count, int k) {
  const int prev_k = k == 0 ? count - 1 : k - 1;
  const int next_k = k + 1 == count ? 0 : k + 1;
  const math::Vector e1 = m.vertices[idx[prev_k]] - m.vertices[idx[k]];
  const math::Vector e2 = m.vertices[idx[next_k]] - m.vertices[idx[k]];
  const float m1 = math::dot(e1, e1);
  const float m2 = math::dot(e2, e2);
  if (!(m1 > math::EPS_LEN_SQ && m2 > math::EPS_LEN_SQ))
    return 0.0f;
  const float d = hs::clamp(math::dot(e1, e2) / sqrtf(m1 * m2), -1.0f, 1.0f);
  return acosf(d) * 180.0f / math::PI_F;
}

/**
 * @brief Builds a PolyMesh from a Solids::* descriptor into the given arena.
 * @tparam Solid Solids::* descriptor type providing NUM_VERTS, NUM_FACES,
 *               vertices, face_counts, and faces.
 * @param mesh Destination mesh; its vertex/face arrays are bound and filled.
 * @param arena Arena from which the mesh arrays are allocated.
 */
template <typename Solid>
inline void build_solid(PolyMesh &mesh, Arena &arena) {
  mesh.vertices.bind(arena, Solid::NUM_VERTS);
  mesh.face_counts.bind(arena, Solid::NUM_FACES);
  mesh.faces.bind(arena, Solid::faces.size());
  for (const auto &v : Solid::vertices)
    mesh.vertices.push_back(v);
  for (auto fc : Solid::face_counts)
    mesh.face_counts.push_back(fc);
  for (auto fi : Solid::faces)
    mesh.faces.push_back(static_cast<uint16_t>(fi));
}

/** @brief Builds an icosahedron into caller-owned arenas. */
inline void build_icosahedron_meshstate(Arena &seed_a, Arena &seed_b,
                                        Arena &geometry, MeshState &mesh) {
  PolyMesh base = Solids::Platonic::icosahedron(seed_a, seed_b);
  mesh.vertices.bind(geometry, base.vertices.size());
  for (const math::Vector &vertex : base.vertices)
    mesh.vertices.push_back(vertex);
  mesh.faces.bind(geometry, base.faces.size());
  mesh.face_counts.bind(geometry, base.face_counts.size());
  for (size_t i = 0; i < base.face_counts.size(); ++i)
    mesh.face_counts.push_back(static_cast<uint8_t>(base.face_counts[i]));
  for (size_t i = 0; i < base.faces.size(); ++i)
    mesh.faces.push_back(base.faces[i]);
}

/** Step table for the ambo/relax/hk54/needle test-fixture recipe. */
inline constexpr Solids::OpStep
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_STEPS[] = {
        {Solids::Op::AMBO},
        {.op = Solids::Op::RELAX,
         .bake = &Solids::RelaxBakes::truncated_icosahedron_ambo_converged},
        {Solids::Op::HANKIN, 54.0f * Solids::IslamicStarPatterns::D2R},
        {Solids::Op::NEEDLE}};
inline constexpr Solids::Recipe
    TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_RECIPE = {
        Solids::SEED_TRUNCATED_ICOSAHEDRON,
        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_STEPS,
        static_cast<uint8_t>(
            std::size(TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_STEPS))};

/**
 * @brief Adapts the plain dodecahedron to the seed function-pointer signature
 *        the mesh probe site tables hold.
 * @param a Output arena for the built mesh.
 * @param b Scratch arena for the intermediate meshes.
 * @return The dodecahedron, allocated in @p a.
 */
inline PolyMesh probe_dodecahedron(Arena &a, Arena &b) {
  return Solids::Platonic::dodecahedron(a, b);
}

/**
 * @brief Builds the needle partition-probe seed: the ambo /
 *        relax(100) / hankin(54 deg) prefix of
 *        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_STEPS.
 * @param a Output arena for the built mesh.
 * @param b Scratch arena for the intermediate meshes.
 * @return The hankin(54 deg) arrival, allocated in @p a.
 * @details The relax runs live at 100 iterations, not from the fixture step's
 *          baked payload.
 */
inline PolyMesh build_ticosa_ambo_relax100_hk54(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::truncated_icosahedron(a, b),
                              a, b)
      .ambo()
      .relax(100)
      .hankin(54.0f * D2R)
      .build();
}

/**
 * @brief Verifies the mesh has vertices and that every one lies on the unit
 *        sphere to within tol.
 * @param m Mesh whose vertices are checked.
 * @param tol Absolute tolerance on the deviation of each vertex length from 1.0.
 * @details Reports the worst deviation across all vertices.
 */
inline void check_all_unit_vertices(const PolyMesh &m, float tol) {
  HS_EXPECT_TRUE(m.vertices.size() > 0);
  float worst = 0.0f;
  for (size_t i = 0; i < m.vertices.size(); ++i)
    worst = fold_worst(worst, std::fabs(m.vertices[i].length() - 1.0f));
  HS_EXPECT_LE(worst, tol);
}

/**
 * @brief Verifies the sum of face_counts equals the flat face-index array length.
 * @param m Mesh whose face_counts and faces arrays are checked.
 * @details The mesh must carry at least one face.
 */
inline void check_face_counts_consistent(const PolyMesh &m) {
  HS_EXPECT_TRUE(m.face_counts.size() > 0);
  size_t total = 0;
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    total += m.face_counts[i];
  HS_EXPECT_EQ(total, m.faces.size());
}

/**
 * @brief Verifies every face index references a valid vertex.
 * @param m Mesh whose face indices are checked against the vertex count.
 * @details The index array must be non-empty. Reports the largest index found
 *          rather than asserting per entry.
 */
inline void check_indices_in_range(const PolyMesh &m) {
  HS_EXPECT_TRUE(m.faces.size() > 0);
  size_t V = m.vertices.size();
  size_t max_index = 0;
  for (size_t i = 0; i < m.faces.size(); ++i)
    max_index = std::max<size_t>(max_index, m.faces[i]);
  HS_EXPECT_LT(max_index, V);
}

/** Longest geodesic edge a healthy solid reaches, as a multiple of its median
 * edge. */
inline constexpr float MAX_SLIVER_EDGE_RATIO = 6.0f;

/**
 * @brief Verifies the longest geodesic edge is within MAX_SLIVER_EDGE_RATIO
 *        times the median edge.
 * @param m Mesh whose face edges are measured as arcs on the unit sphere.
 */
inline void check_no_sliver_edges(const PolyMesh &m) {
  std::vector<float> edges;
  size_t off = 0;
  for (size_t f = 0; f < m.face_counts.size(); ++f) {
    const int n = m.face_counts[f];
    for (int i = 0; i < n; ++i) {
      const math::Vector u = m.vertices[m.faces[off + i]].normalized();
      const math::Vector v =
          m.vertices[m.faces[off + (i + 1) % n]].normalized();
      edges.push_back(
          std::acos(std::max(-1.0f, std::min(1.0f, math::dot(u, v)))));
    }
    off += n;
  }
  HS_EXPECT_TRUE(!edges.empty());
  if (edges.empty())
    return;
  std::sort(edges.begin(), edges.end());
  const float median = edges[edges.size() / 2];
  HS_EXPECT_LE(edges.back(), MAX_SLIVER_EDGE_RATIO * median);
}

/** @brief Newell vector of a polygon supplied in cyclic vertex order. */
template <typename VertexAt>
inline math::Vector newell_normal(int count, VertexAt vertex_at) {
  math::Vector n(0, 0, 0);
  for (int k = 0; k < count; ++k) {
    const math::Vector curr = vertex_at(k);
    const math::Vector next = vertex_at((k + 1) % count);
    n.x += (curr.y - next.y) * (curr.z + next.z);
    n.y += (curr.z - next.z) * (curr.x + next.x);
    n.z += (curr.x - next.x) * (curr.y + next.y);
  }
  return n;
}

/** @brief Index-run offset of a PolyMesh face. */
inline size_t poly_face_offset(const PolyMesh &mesh, size_t face) {
  size_t offset = 0;
  for (size_t i = 0; i < face; ++i)
    offset += mesh.face_counts[i];
  return offset;
}

/** @brief Finds a directed graph edge, or -1 if it is absent. */
template <typename Edges, typename Node>
inline int find_directed_edge(const Edges &edges, Node from, Node to) {
  for (size_t e = 0; e < std::size(edges); ++e)
    if (edges[e].from_node == from && edges[e].to_node == to)
      return static_cast<int>(e);
  return -1;
}

/**
 * @brief Computes a face normal via Newell's method.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Unnormalised normal vector for the face; its magnitude is twice the
 *         planar face area.
 */
inline math::Vector face_newell_normal(const PolyMesh &m,
                                       size_t face_idx_offset, int count) {
  return newell_normal(
      count, [&](int k) { return m.vertices[m.faces[face_idx_offset + k]]; });
}

/**
 * @brief Computes a face's Newell area vector.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Face normal scaled to the planar face area — half
 *         face_newell_normal(), which carries twice the area.
 */
inline math::Vector face_area_vector(const PolyMesh &m, size_t face_idx_offset,
                                     int count) {
  return face_newell_normal(m, face_idx_offset, count) * 0.5f;
}

/**
 * @brief Computes the unweighted centroid of a face's vertex positions.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Arithmetic mean of the face's vertex positions.
 */
inline math::Vector face_centroid_pos(const PolyMesh &m, size_t face_idx_offset,
                                      int count) {
  math::Vector c(0, 0, 0);
  for (int k = 0; k < count; ++k)
    c = c + m.vertices[m.faces[face_idx_offset + k]];
  return c * (1.0f / static_cast<float>(count));
}

/**
 * @brief Projects a face's vertex-average centroid onto the unit sphere.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Normalised centroid direction.
 */
inline math::Vector face_centroid_unit(const PolyMesh &m,
                                       size_t face_idx_offset, int count) {
  math::Vector c(0.0f, 0.0f, 0.0f);
  for (int k = 0; k < count; ++k)
    c = c + m.vertices[m.faces[face_idx_offset + k]];
  return c.normalized();
}

/** @brief Copies a fixed solid into mutable renderer storage, optionally omitting its last face. */
template <typename Solid>
inline void build_meshstate_solid(MeshState &mesh, Arena &arena,
                                  bool open = false) {
  mesh.vertices.bind(arena, Solid::vertices.size());
  for (const auto &vertex : Solid::vertices)
    mesh.vertices.push_back(vertex);
  const size_t face_count = Solid::face_counts.size() - (open ? 1 : 0);
  mesh.face_counts.bind(arena, face_count);
  size_t index_count = 0;
  for (size_t face = 0; face < face_count; ++face) {
    mesh.face_counts.push_back(Solid::face_counts[face]);
    index_count += Solid::face_counts[face];
  }
  mesh.faces.bind(arena, index_count);
  for (size_t i = 0; i < index_count; ++i)
    mesh.faces.push_back(static_cast<uint16_t>(Solid::faces[i]));
}

/**
 * @brief Asserts two meshes have equal vertex, face_counts and faces sizes.
 * @param a First mesh.
 * @param b Second mesh.
 * @return True when all three sizes match.
 */
inline bool expect_same_mesh_sizes(const PolyMesh &a, const PolyMesh &b) {
  HS_EXPECT_EQ(a.vertices.size(), b.vertices.size());
  HS_EXPECT_EQ(a.face_counts.size(), b.face_counts.size());
  HS_EXPECT_EQ(a.faces.size(), b.faces.size());
  return a.vertices.size() == b.vertices.size() &&
         a.face_counts.size() == b.face_counts.size() &&
         a.faces.size() == b.faces.size();
}

/**
 * @brief Asserts two meshes are bitwise identical: same counts and same bytes
 *        for vertices, face_counts, and faces.
 * @param m1,m2 Meshes to compare.
 */
inline void check_bitwise_equal_meshes(const PolyMesh &m1, const PolyMesh &m2) {
  if (!expect_same_mesh_sizes(m1, m2))
    return;
  HS_EXPECT_EQ(std::memcmp(m1.vertices.data(), m2.vertices.data(),
                           m1.vertices.size() * sizeof(math::Vector)),
               0);
  HS_EXPECT_EQ(std::memcmp(m1.face_counts.data(), m2.face_counts.data(),
                           m1.face_counts.size() * sizeof(uint8_t)),
               0);
  HS_EXPECT_EQ(std::memcmp(m1.faces.data(), m2.faces.data(),
                           m1.faces.size() * sizeof(uint16_t)),
               0);
}

/**
 * @brief Asserts exact element equality (IEEE ==) between two meshes.
 * @param a First mesh.
 * @param b Second mesh.
 * @details Unlike check_bitwise_equal_meshes, float comparison accepts signed
 *          zero differences and rejects NaNs.
 */
inline void check_meshes_identical(const PolyMesh &a, const PolyMesh &b) {
  if (!expect_same_mesh_sizes(a, b))
    return;
  for (size_t i = 0; i < a.vertices.size(); ++i) {
    HS_EXPECT_EQ(a.vertices[i].x, b.vertices[i].x);
    HS_EXPECT_EQ(a.vertices[i].y, b.vertices[i].y);
    HS_EXPECT_EQ(a.vertices[i].z, b.vertices[i].z);
  }
  for (size_t i = 0; i < a.face_counts.size(); ++i)
    HS_EXPECT_EQ(a.face_counts[i], b.face_counts[i]);
  for (size_t i = 0; i < a.faces.size(); ++i)
    HS_EXPECT_EQ(a.faces[i], b.faces[i]);
}

/**
 * @brief Newell-sum area of face fi.
 */
inline float poly_face_area(const PolyMesh &m, size_t fi) {
  const size_t off = poly_face_offset(m, fi);
  const int n = m.face_counts[fi];
  math::Vector s(0.0f, 0.0f, 0.0f);
  for (int k = 1; k + 1 < n; ++k) {
    const math::Vector e1 =
        m.vertices[m.faces[off + k]] - m.vertices[m.faces[off]];
    const math::Vector e2 =
        m.vertices[m.faces[off + k + 1]] - m.vertices[m.faces[off]];
    s = s + math::cross(e1, e2);
  }
  return 0.5f * s.length();
}

namespace conway_tests {
/**
 * @brief Counts faces whose winding produces an inward-pointing normal.
 * @param m Mesh to inspect; normals are taken relative to the mesh origin.
 * @return Number of faces whose normal disagrees with the centroid direction.
 * @details On a unit-sphere mesh, correct CCW-from-outside winding gives a
 *          normal that agrees in direction with the face centroid; a face
 *          failing this is wound backwards. Finite degenerate normals are
 *          ignored; non-finite normal lengths count as bad.
 */
inline int count_inward_winding(const PolyMesh &m) {
  int bad = 0;
  size_t offset = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    int count = m.face_counts[fi];
    if (count >= 3) {
      math::Vector n = face_newell_normal(m, offset, count);
      math::Vector c = face_centroid_pos(m, offset, count);
      const float LENGTH = n.length();
      if (!std::isfinite(LENGTH) || (LENGTH > 1e-6f && math::dot(n, c) <= 0.0f))
        ++bad;
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
  struct DirEdge {
    uint16_t u, v;
  };
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

  std::sort(edges.begin(), edges.end(), [](const DirEdge &a, const DirEdge &b) {
    return a.u != b.u ? a.u < b.u : a.v < b.v;
  });

  int dup = 0;
  for (size_t i = 1; i < edges.size(); ++i) {
    if (edges[i].u == edges[i - 1].u && edges[i].v == edges[i - 1].v)
      ++dup;
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
 * @brief Checks two-face edge incidence and Euler characteristic V - E + F == 2.
 * @param m Mesh to validate.
 * @details Every undirected edge must be shared by exactly two faces; the scan
 *          reports the extreme fan-out across all edges. Vertex-fan
 *          connectivity and connectedness are not checked.
 */
inline void check_euler_characteristic_two(const PolyMesh &m) {
  std::vector<std::pair<uint16_t, uint16_t>> edges;
  edges.reserve(m.faces.size());
  size_t offset = 0;
  for (size_t fi = 0; fi < m.face_counts.size(); ++fi) {
    int count = m.face_counts[fi];
    for (int k = 0; k < count; ++k) {
      uint16_t u = m.faces[offset + k];
      uint16_t v = m.faces[offset + (k + 1) % count];
      if (u > v)
        std::swap(u, v);
      edges.push_back({u, v});
    }
    offset += count;
  }
  std::sort(edges.begin(), edges.end());

  int E = 0;
  int min_fanout = edges.empty() ? 2 : (int)edges.size();
  int max_fanout = edges.empty() ? 2 : 0;
  for (size_t i = 0; i < edges.size();) {
    size_t j = i;
    while (j < edges.size() && edges[j] == edges[i])
      ++j;
    min_fanout = std::min(min_fanout, (int)(j - i));
    max_fanout = std::max(max_fanout, (int)(j - i));
    ++E;
    i = j;
  }
  // Each undirected edge must bound exactly two faces.
  HS_EXPECT_EQ(min_fanout, 2);
  HS_EXPECT_EQ(max_fanout, 2);

  int V = (int)m.vertices.size();
  int F = (int)m.face_counts.size();
  HS_EXPECT_EQ(V - E + F, 2);
}

/**
 * @brief Histogram of vertex degree (incident faces) → number of such vertices.
 * @param m Mesh to inspect.
 * @details In a closed manifold a vertex's incident-face count equals its edge
 *          degree. An unreferenced (degree-0) vertex shows up as a {0, n}
 *          bucket. Out-of-range face indices occupy the {-1, n} bucket.
 */
inline std::map<int, int> vertex_degree_histogram(const PolyMesh &m) {
  std::vector<int> degree(m.vertices.size(), 0);
  std::map<int, int> hist;
  for (uint16_t idx : m.faces) {
    if (idx < degree.size())
      ++degree[idx];
    else
      ++hist[-1];
  }
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

} // namespace conway_tests

} // namespace hs_test
