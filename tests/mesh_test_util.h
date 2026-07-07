/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
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
  for (const auto &v : Solid::vertices) mesh.vertices.push_back(v);
  for (auto fc : Solid::face_counts) mesh.face_counts.push_back(fc);
  for (auto fi : Solid::faces)
    mesh.faces.push_back(static_cast<uint16_t>(fi));
}

/**
 * @brief Verifies every vertex lies on the unit sphere to within tol.
 * @param m Mesh whose vertices are checked.
 * @param tol Absolute tolerance on the deviation of each vertex length from 1.0.
 */
inline void check_all_unit_vertices(const PolyMesh &m, float tol) {
  for (size_t i = 0; i < m.vertices.size(); ++i) {
    float len = m.vertices[i].length();
    HS_EXPECT_NEAR(len, 1.0f, tol);
  }
}

/**
 * @brief Verifies the sum of face_counts equals the flat face-index array length.
 * @param m Mesh whose face_counts and faces arrays are checked.
 * @details Σ face_counts must equal m.faces.size() for the flat index layout to
 *          be self-consistent.
 */
inline void check_face_counts_consistent(const PolyMesh &m) {
  size_t total = 0;
  for (size_t i = 0; i < m.face_counts.size(); ++i)
    total += m.face_counts[i];
  HS_EXPECT_EQ(total, m.faces.size());
}

/**
 * @brief Verifies every face index references a valid vertex.
 * @param m Mesh whose face indices are checked against the vertex count.
 * @details Each entry of m.faces must be strictly less than m.vertices.size().
 */
inline void check_indices_in_range(const PolyMesh &m) {
  size_t V = m.vertices.size();
  for (size_t i = 0; i < m.faces.size(); ++i)
    HS_EXPECT_TRUE(m.faces[i] < V);
}

} // namespace hs_test
