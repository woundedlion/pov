/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared mesh test fixtures. The conway/mesh/hankin/solids suites all build a
 * PolyMesh from a Solids::* descriptor and check that its vertices land on the
 * unit sphere.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>
#include "core/mesh/mesh.h"
#include "core/mesh/solids.h"
#include "tests/test_harness.h"

namespace hs_test {

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

/** Step table for the ambo/relax/hk54/needle recipe: a test fixture, not in
 * islamic_registry; the reconcile tests' canonical needle-ending recipe. */
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
 * @details Shared by the morph and opchain probe suites, whose site tables both
 *          name it as an un-transformed seed.
 */
inline PolyMesh probe_dodecahedron(Arena &a, Arena &b) {
  return Solids::Platonic::dodecahedron(a, b);
}

/**
 * @brief Builds the seed of the needle recipe's gated swaps: the ambo /
 *        relax(100) / hankin(54 deg) prefix of
 *        TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_STEPS.
 * @param a Output arena for the built mesh.
 * @param b Scratch arena for the intermediate meshes.
 * @return The hankin(54 deg) arrival, allocated in @p a.
 * @details Shared by the morph and opchain probe suites, which both pin shapes
 *          against this chain; a spec change must move exactly one build here.
 *          The relax is live at 100 iterations rather than the shipping step's
 *          baked payload, so the probes measure the operator, not the bake.
 */
inline PolyMesh build_ticosa_ambo_relax100_hk54(Arena &a, Arena &b) {
  using Solids::IslamicStarPatterns::D2R;
  return Solids::SolidBuilder(Solids::Archimedean::truncatedIcosahedron(a, b),
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
    worst = std::max(worst, std::fabs(m.vertices[i].length() - 1.0f));
  HS_EXPECT_LE(worst, tol);
}

/**
 * @brief Verifies the sum of face_counts equals the flat face-index array length.
 * @param m Mesh whose face_counts and faces arrays are checked.
 * @details Σ face_counts must equal m.faces.size() for the flat index layout to
 *          be self-consistent. The mesh must carry at least one face.
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
 * @details Each entry of m.faces must be strictly less than m.vertices.size().
 *          The index array must be non-empty. Reports the largest index found
 *          rather than asserting per entry, for the same floor reason as
 *          check_all_unit_vertices.
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
 * edge: registry recipes measure at most ~3.4x, a hankin resonance sling ~24x. */
inline constexpr float MAX_SLIVER_EDGE_RATIO = 6.0f;

/**
 * @brief Verifies the mesh has no sliver faces: its longest geodesic edge stays
 *        within MAX_SLIVER_EDGE_RATIO times the median edge.
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

/**
 * @brief Computes a face normal via Newell's method.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Unnormalised normal vector for the face; its magnitude is twice the
 *         planar face area.
 * @details Newell's method is robust for non-planar faces (e.g. curved faces
 *          on the unit sphere) where a simple cross product would be ambiguous.
 */
inline math::Vector face_newell_normal(const PolyMesh &m,
                                       size_t face_idx_offset, int count) {
  math::Vector n(0, 0, 0);
  for (int k = 0; k < count; ++k) {
    const math::Vector &curr = m.vertices[m.faces[face_idx_offset + k]];
    const math::Vector &next =
        m.vertices[m.faces[face_idx_offset + (k + 1) % count]];
    n.x += (curr.y - next.y) * (curr.z + next.z);
    n.y += (curr.z - next.z) * (curr.x + next.x);
    n.z += (curr.x - next.x) * (curr.y + next.y);
  }
  return n;
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
 * @return Normalised centroid direction; the vertex sum is normalised directly,
 *         so the division by count of face_centroid_pos() is skipped.
 */
inline math::Vector face_centroid_unit(const PolyMesh &m,
                                       size_t face_idx_offset, int count) {
  math::Vector c(0.0f, 0.0f, 0.0f);
  for (int k = 0; k < count; ++k)
    c = c + m.vertices[m.faces[face_idx_offset + k]];
  return c.normalized();
}

} // namespace hs_test
