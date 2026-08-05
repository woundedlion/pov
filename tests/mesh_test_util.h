/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared mesh test fixtures. The conway/mesh/hankin/solids suites all build a
 * PolyMesh from a Solids::* descriptor and check that its vertices land on the
 * unit sphere.
 */
#pragma once

#include "core/mesh/mesh.h"
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

/**
 * @brief Verifies the mesh has vertices and that every one lies on the unit
 *        sphere to within tol.
 * @param m Mesh whose vertices are checked.
 * @param tol Absolute tolerance on the deviation of each vertex length from 1.0.
 */
inline void check_all_unit_vertices(const PolyMesh &m, float tol) {
  HS_EXPECT_TRUE(m.vertices.size() > 0);
  for (size_t i = 0; i < m.vertices.size(); ++i) {
    float len = m.vertices[i].length();
    HS_EXPECT_NEAR(len, 1.0f, tol);
  }
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
 *          The index array must be non-empty.
 */
inline void check_indices_in_range(const PolyMesh &m) {
  HS_EXPECT_TRUE(m.faces.size() > 0);
  size_t V = m.vertices.size();
  for (size_t i = 0; i < m.faces.size(); ++i)
    HS_EXPECT_TRUE(m.faces[i] < V);
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
inline Vector face_newell_normal(const PolyMesh &m, size_t face_idx_offset,
                                 int count) {
  Vector n(0, 0, 0);
  for (int k = 0; k < count; ++k) {
    const Vector &curr = m.vertices[m.faces[face_idx_offset + k]];
    const Vector &next = m.vertices[m.faces[face_idx_offset + (k + 1) % count]];
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
 * @brief Projects a face's vertex-average centroid onto the unit sphere.
 * @param m Mesh owning the vertices and face-index array.
 * @param face_idx_offset Offset into m.faces where this face's indices begin.
 * @param count Number of vertices (sides) in the face.
 * @return Normalised centroid direction; the vertex sum is normalised directly,
 *         so the division by count of face_centroid_pos() is skipped.
 */
inline Vector face_centroid_unit(const PolyMesh &m, size_t face_idx_offset,
                                 int count) {
  Vector c(0.0f, 0.0f, 0.0f);
  for (int k = 0; k < count; ++k)
    c = c + m.vertices[m.faces[face_idx_offset + k]];
  return c.normalized();
}

} // namespace hs_test
