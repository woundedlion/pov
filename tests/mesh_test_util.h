/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Shared mesh test fixtures. The conway/mesh/hankin/solids suites all build a
 * PolyMesh from a Solids::* descriptor and check that its vertices land on the
 * unit sphere; these helpers live here so the fixture cannot drift between
 * suites (the unit-sphere tolerance was previously 1e-3 in one file and 1e-2
 * in another). The tolerance is an explicit parameter so each caller states
 * the precision it actually needs.
 */
#pragma once

#include "core/mesh.h"
#include "tests/test_harness.h"

namespace hs_test {

// Build a PolyMesh from a Solids::* descriptor into the given arena.
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

// Verify every vertex lies on the unit sphere to within `tol`.
inline void check_all_unit_vertices(const PolyMesh &m, float tol) {
  for (size_t i = 0; i < m.vertices.size(); ++i) {
    float len = m.vertices[i].length();
    HS_EXPECT_NEAR(len, 1.0f, tol);
  }
}

} // namespace hs_test
