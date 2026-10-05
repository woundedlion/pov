/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// §7.4 Bridge convergence: the tetrahedral edges that cross symmetry families.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies snub(tetrahedron, 0.5, SNUB_BRIDGE_TWIST).relax(50) is the
 *        regular icosahedron: 12 vertices, 20 triangles, equal edges on the
 *        unit sphere (relax supplies the canonical form, as the registry snub
 *        chains rely on), at the bridge's tabled arrival twist.
 */
inline void test_snub_tetrahedron_relax_converges_to_icosahedron() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh snubbed =
      MeshOps::snub(tetra, target, temp, 0.5f, ConwayGraph::SNUB_BRIDGE_TWIST);
  PolyMesh relaxed = MeshOps::relax(snubbed, aux, temp, 50);

  HS_EXPECT_EQ(relaxed.vertices.size(), (size_t)12);
  HS_EXPECT_EQ(relaxed.face_counts.size(), (size_t)20);
  for (size_t fi = 0; fi < relaxed.face_counts.size(); ++fi)
    HS_EXPECT_EQ((int)relaxed.face_counts[fi], 3);
  check_face_counts_consistent(relaxed);
  check_indices_in_range(relaxed);
  check_all_unit_vertices(relaxed, 1e-3f);

  // Regular icosahedron: every edge equals the mean (chord ~1.0515).
  HS_EXPECT_LE(max_edge_length_deviation(relaxed), 0.02f);
}

/**
 * @brief Verifies ambo(tetrahedron) is the regular octahedron: 6 vertices,
 *        8 triangles, equal edges, and a vertex-set bijection onto the
 *        Octahedron seed (normalized tetra edge midpoints are the ±axes).
 */
inline void test_ambo_tetrahedron_is_regular_octahedron() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh a = MeshOps::ambo(tetra, target, temp);

  HS_EXPECT_EQ(a.vertices.size(), (size_t)6);
  HS_EXPECT_EQ(a.face_counts.size(), (size_t)8);
  for (size_t fi = 0; fi < a.face_counts.size(); ++fi)
    HS_EXPECT_EQ((int)a.face_counts[fi], 3);
  check_face_counts_consistent(a);
  check_indices_in_range(a);
  check_all_unit_vertices(a, 1e-3f);
  HS_EXPECT_LE(max_edge_length_deviation(a), 1e-4f);

  bool used[Solids::Octahedron::NUM_VERTS] = {};
  for (size_t i = 0; i < a.vertices.size(); ++i) {
    int match = -1;
    for (int j = 0; j < Solids::Octahedron::NUM_VERTS; ++j) {
      if (!used[j] &&
          (a.vertices[i] - Solids::Octahedron::vertices[j]).length() <= 1e-4f) {
        match = j;
        break;
      }
    }
    HS_EXPECT_TRUE(match >= 0);
    if (match >= 0)
      used[match] = true;
  }
}
