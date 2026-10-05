/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// §7.3 Settle correspondence: relax output vertex order is the identity over
// its input, so a relaxed endpoint is per-vertex slerpable.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies relax(50) on a registry form (expand(dodecahedron), the
 *        rhombicosidodecahedron chain) preserves vertex order and topology.
 * @details Counts equal, face_counts/faces byte-identical, and each relaxed
 *          vertex stays strictly nearest to its own input vertex — a relax
 *          rewrite that reorders vertices fails here loudly.
 */
inline void test_relax_is_vertex_order_identity() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh dodeca;
  build_solid<Solids::Dodecahedron>(dodeca, temp);
  PolyMesh unrelaxed = MeshOps::expand(dodeca, target, temp);
  PolyMesh relaxed = MeshOps::relax(unrelaxed, aux, temp, 50);

  HS_EXPECT_EQ(relaxed.vertices.size(), unrelaxed.vertices.size());
  HS_EXPECT_EQ(relaxed.face_counts.size(), unrelaxed.face_counts.size());
  HS_EXPECT_EQ(relaxed.faces.size(), unrelaxed.faces.size());
  HS_EXPECT_EQ(std::memcmp(relaxed.face_counts.data(),
                           unrelaxed.face_counts.data(),
                           relaxed.face_counts.size() * sizeof(uint8_t)),
               0);
  HS_EXPECT_EQ(std::memcmp(relaxed.faces.data(), unrelaxed.faces.data(),
                           relaxed.faces.size() * sizeof(uint16_t)),
               0);

  for (size_t i = 0; i < relaxed.vertices.size(); ++i) {
    size_t nearest = 0;
    float best = 1e9f;
    for (size_t j = 0; j < unrelaxed.vertices.size(); ++j) {
      const float d = (relaxed.vertices[i] - unrelaxed.vertices[j]).length();
      if (d < best) {
        best = d;
        nearest = j;
      }
    }
    HS_EXPECT_EQ(nearest, i);
  }
}
