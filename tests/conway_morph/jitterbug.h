/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// Jitterbug bridge (icosahedron <-> octahedron on the tetra snub family):
// both endpoint parameter pins plus the clamped-leg topology sweep.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies snub(tetrahedron, T_JITTERBUG_ICOSA, TWIST_JITTERBUG_ICOSA)
 *        is the regular icosahedron directly — 12 vertices, 20 triangles, all
 *        30 edges equal on the unit sphere with no relax — pinning the
 *        jitterbug bridge's icosa endpoint parameters.
 */
inline void test_jitterbug_icosa_point_is_regular() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh s =
      MeshOps::snub(tetra, target, temp, ConwayGraph::T_JITTERBUG_ICOSA,
                    ConwayGraph::TWIST_JITTERBUG_ICOSA);

  HS_EXPECT_EQ(s.vertices.size(), (size_t)12);
  HS_EXPECT_EQ(s.face_counts.size(), (size_t)20);
  for (size_t fi = 0; fi < s.face_counts.size(); ++fi)
    HS_EXPECT_EQ((int)s.face_counts[fi], 3);
  check_face_counts_consistent(s);
  check_indices_in_range(s);
  check_all_unit_vertices(s, 1e-3f);
  HS_EXPECT_LE(max_edge_length_deviation(s), 1e-5f);
}

/**
 * @brief Verifies the jitterbug octa endpoint snub(tetrahedron, 0.5, -pi/3):
 *        the 12 vertices merge pairwise onto the registry octahedron's 6 and
 *        exactly the 12 edge-orbit faces are zero-area (the SDF zero-area cull
 *        hides them, so the clean swap to the held octahedron changes no
 *        pixels).
 */
inline void test_jitterbug_octa_end_covers_octahedron() {
  Arena target(morph_target_buf, sizeof(morph_target_buf));
  Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, temp);
  PolyMesh s = MeshOps::snub(tetra, target, temp, 0.5f,
                             ConwayGraph::TWIST_JITTERBUG_OCTA);
  PolyMesh octa;
  build_solid<Solids::Octahedron>(octa, aux);

  check_pairwise_vertex_cover(s, octa, 1e-4f);

  // Emission order: 4 primary + 4 vertex-orbit faces (the octahedron's 8),
  // then the 12 collapsed edge-orbit faces.
  HS_EXPECT_EQ(s.face_counts.size(), (size_t)20);
  int zero_area = 0;
  for (size_t fi = 0; fi < s.face_counts.size(); ++fi) {
    const float a = poly_face_area(s, fi);
    if (a < 1e-6f)
      ++zero_area;
    else
      HS_EXPECT_GT(a, 0.5f); // equilateral sqrt(2)-side triangle: ~0.866
    if (fi < 8)
      HS_EXPECT_GT(a, 0.5f);
  }
  HS_EXPECT_EQ(zero_area, 12);
}

/**
 * @brief Verifies the jitterbug leg exactly as OpLeg runs it — t from
 *        the icosa point to the T_JITTERBUG_OCTA_MIN clamp with the tabled twist
 *        endpoints — holds V12/F20/E30 with two-face edge incidence,
 *        >= 3-side faces, and unit vertices, with the collapsing edge never
 *        shorter than the clamp chord (spec section 7.2 for the edge).
 */
inline void test_jitterbug_sweep_holds_topology() {
  constexpr int SAMPLES = 17;
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
  PolyMesh tetra;
  build_solid<Solids::Tetrahedron>(tetra, aux);

  for (int s = 0; s < SAMPLES; ++s) {
    const float k = static_cast<float>(s) / (SAMPLES - 1);
    const float t =
        ConwayGraph::T_JITTERBUG_ICOSA +
        (ConwayGraph::T_JITTERBUG_OCTA_MIN - ConwayGraph::T_JITTERBUG_ICOSA) *
            k;
    const float twist = ConwayGraph::TWIST_JITTERBUG_ICOSA +
                        (ConwayGraph::TWIST_JITTERBUG_OCTA -
                         ConwayGraph::TWIST_JITTERBUG_ICOSA) *
                            k;

    Arena target(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
    PolyMesh out = MeshOps::snub(tetra, target, temp, t, twist);

    HS_EXPECT_EQ(out.vertices.size(), (size_t)12);
    HS_EXPECT_EQ(out.face_counts.size(), (size_t)20);
    HS_EXPECT_EQ(out.faces.size(), (size_t)60); // E = I / 2 = 30
    for (size_t fi = 0; fi < out.face_counts.size(); ++fi)
      HS_EXPECT_TRUE(out.face_counts[fi] >= 3);
    check_face_counts_consistent(out);
    check_indices_in_range(out);
    check_all_unit_vertices(out, 1e-3f);
    conway_tests::check_euler_characteristic_two(out);

    // The clamp keeps the shortest edge above the sliver threshold.
    float min_edge = 1e9f;
    size_t off = 0;
    for (size_t fi = 0; fi < out.face_counts.size(); ++fi) {
      const int c = out.face_counts[fi];
      for (int j = 0; j < c; ++j)
        min_edge = std::min(
            min_edge,
            math::distance_between(out.vertices[out.faces[off + j]],
                                   out.vertices[out.faces[off + (j + 1) % c]]));
      off += c;
    }
    HS_EXPECT_GE(min_edge, 0.019f);
  }
}
