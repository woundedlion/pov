/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// §7.1 Endpoint exactness: sweeping to an edge endpoint arrives at the
// registry generator's output.
// ---------------------------------------------------------------------------

/** How an edge endpoint is compared against its node's registry output. */
enum class EndRegime {
  EXACT,        /**< Same code path: bitwise vertices, identical topology. */
  VERTEX_MATCH, /**< Same geometry, different vertex order (dual-family ambo,
                     ambo(tetra) bridge). */
  REGULAR,      /**< Relax-canonical arrival in a walk-dependent orientation
                     (tetra -> icosa bridge). */
  PAIR_COVER,   /**< Jitterbug octa end: vertices merge pairwise onto the node
                     mesh's. */
  BAKED_RELAX,  /**< Registry node ends in relax_baked: identical topology, and
                     vertices within the relax convergence gate. */
};

/**
 * @brief Whether a simple-registry node's generator ends in relax_baked.
 * @param node Simple-registry index.
 * @return True for the nodes carrying a flash bake.
 * @details The baked mesh holds host-IEEE relax() bits; a live relax agrees
 *   bitwise only under matching float semantics.
 */
inline bool is_relax_baked_node(uint8_t node) {
  return node == ConwayGraph::TRUNCATED_CUBOCTAHEDRON ||
         node == ConwayGraph::SNUB_CUBE ||
         node == ConwayGraph::RHOMBICOSIDODECAHEDRON ||
         node == ConwayGraph::SNUB_DODECAHEDRON ||
         node == ConwayGraph::TRUNCATED_ICOSIDODECAHEDRON;
}

/**
 * @brief Comparison regime of an edge's to_node end.
 * @param e Edge to classify.
 * @return EXACT when op(seed, t_to) [+ relax] is the to_node registry chain;
 *         BAKED_RELAX when that chain ends in a flash bake rather than the live
 *         relax the leg runs; VERTEX_MATCH for arrivals off the registry seed
 *         (dual-family ambo, non-settle bridges); REGULAR for the settling
 *         bridge, whose relax orientation tracks the seed frame, not the
 *         registry icosahedron; PAIR_COVER for the jitterbug bridge's collapsed
 *         octa end.
 */
inline EndRegime to_end_regime(const ConwayGraph::EdgeSpec &e) {
  using namespace ConwayGraph;
  if (is_jitterbug_edge(e))
    return EndRegime::PAIR_COVER;
  if (e.to_node == CUBOCTAHEDRON && e.seed_solid == OCTAHEDRON)
    return EndRegime::VERTEX_MATCH;
  if (e.to_node == ICOSIDODECAHEDRON && e.seed_solid == ICOSAHEDRON)
    return EndRegime::VERTEX_MATCH;
  if (e.bridge)
    return e.settle ? EndRegime::REGULAR : EndRegime::VERTEX_MATCH;
  if (e.settle && is_relax_baked_node(e.to_node))
    return EndRegime::BAKED_RELAX;
  return EndRegime::EXACT;
}

/**
 * @brief Asserts two meshes share a bitwise-identical face list and agree
 *        vertex-for-vertex within the relax convergence gate.
 * @param got Mesh whose settle end ran a live relax.
 * @param want Registry mesh whose chain ends in a flash bake.
 * @details Vertex tolerance is sqrt(RELAX_CONVERGE_EPS_SQ); vertex order,
 *   face_counts and faces stay exact.
 */
inline void check_equal_within_relax_gate(const PolyMesh &got,
                                          const PolyMesh &want) {
  const float tol = std::sqrt(MeshOps::RELAX_CONVERGE_EPS_SQ);
  HS_EXPECT_EQ(got.vertices.size(), want.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), want.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), want.faces.size());
  if (got.vertices.size() != want.vertices.size() ||
      got.face_counts.size() != want.face_counts.size() ||
      got.faces.size() != want.faces.size())
    return;
  for (size_t i = 0; i < got.vertices.size(); ++i)
    HS_EXPECT_LE(math::distance_between(got.vertices[i], want.vertices[i]),
                 tol);
  for (size_t i = 0; i < got.face_counts.size(); ++i)
    HS_EXPECT_EQ((int)got.face_counts[i], (int)want.face_counts[i]);
  for (size_t i = 0; i < got.faces.size(); ++i)
    HS_EXPECT_EQ(got.faces[i], want.faces[i]);
}

/**
 * @brief Asserts two meshes carry the same geometry up to vertex order:
 *        equal counts, equal face-type histograms, and a vertex-set bijection
 *        within tol.
 */
inline void check_equal_up_to_vertex_order(const PolyMesh &got,
                                           const PolyMesh &want, float tol) {
  HS_EXPECT_EQ(got.vertices.size(), want.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), want.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), want.faces.size());
  HS_EXPECT_TRUE(conway_tests::face_type_histogram(got) ==
                 conway_tests::face_type_histogram(want));
  if (got.vertices.size() != want.vertices.size())
    return;
  std::vector<bool> used(want.vertices.size(), false);
  for (size_t i = 0; i < got.vertices.size(); ++i) {
    bool matched = false;
    for (size_t j = 0; j < want.vertices.size(); ++j) {
      if (!used[j] && (got.vertices[i] - want.vertices[j]).length() <= tol) {
        used[j] = true;
        matched = true;
        break;
      }
    }
    HS_EXPECT_TRUE(matched);
  }
}

/**
 * @brief Asserts a mesh is the registry solid's regular form in an arbitrary
 *        orientation: equal counts, equal face-type histograms, unit vertices,
 *        and near-equal edge lengths.
 */
inline void check_regular_form(const PolyMesh &got, const PolyMesh &want,
                               float edge_dev_tol) {
  HS_EXPECT_EQ(got.vertices.size(), want.vertices.size());
  HS_EXPECT_EQ(got.face_counts.size(), want.face_counts.size());
  HS_EXPECT_EQ(got.faces.size(), want.faces.size());
  HS_EXPECT_TRUE(conway_tests::face_type_histogram(got) ==
                 conway_tests::face_type_histogram(want));
  check_all_unit_vertices(got, 1e-3f);
  HS_EXPECT_LE(max_edge_length_deviation(got), edge_dev_tol);
}

/**
 * @brief Verifies every ConwayGraph edge endpoint against its node's registry
 *        generator: exact on the registry code path, geometric tolerance for
 *        the t = 0 ends, the off-registry-seed arrivals and the flash-baked
 *        relax arrivals.
 * @details Seeds are built via the registry generators, so the DERIVE_AMBO
 *          rows run the exact bevel decomposition of their to_node chains.
 */
inline void test_edge_endpoints_match_registry() {
  constexpr size_t TARGET_HALF = sizeof(morph_target_buf) / 2;
  constexpr size_t TEMP_HALF = sizeof(morph_temp_buf) / 2;
  constexpr size_t AUX_HALF = sizeof(morph_aux_buf) / 2;
  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];
    const int failed_before = hs_test::stats().failed;

    Arena sa(morph_aux_buf, AUX_HALF);
    Arena sb(morph_aux_buf + AUX_HALF, AUX_HALF);
    PolyMesh seed = Solids::simple_registry[e.seed_solid].generate(sa, sb);

    // from end: t = 0 emits expanded topology, so compare op(seed, T_EPS)
    // primaries against the seed (= the from_node registry mesh); a non-zero
    // t_from is the from_node registry chain itself, except the jitterbug
    // icosa point, which is regular in the tetra frame, not the registry
    // orientation.
    {
      Arena oa(morph_temp_buf, TEMP_HALF);
      Arena ob(morph_temp_buf + TEMP_HALF, TEMP_HALF);
      if (e.t_from == 0.0f) {
        PolyMesh got = run_edge_op(e, seed, oa, ob, T_EPS, e.twist_from);
        const int per_corner = e.op == ConwayGraph::MorphOp::TRUNCATE ? 2 : 1;
        const float tol = e.op == ConwayGraph::MorphOp::TRUNCATE
                              ? PRIMARY_CORNER_TOL_TRUNCATE
                              : PRIMARY_CORNER_TOL_SINGLE;
        check_primary_faces_match_seed(seed, got, per_corner, tol);
      } else {
        Arena ra(morph_target_buf, TARGET_HALF);
        Arena rb(morph_target_buf + TARGET_HALF, TARGET_HALF);
        PolyMesh want = Solids::simple_registry[e.from_node].generate(ra, rb);
        PolyMesh got = run_edge_op(e, seed, oa, ob, e.t_from, e.twist_from);
        if (ConwayGraph::is_jitterbug_edge(e))
          check_regular_form(got, want, 1e-4f);
        else
          hs_test::check_meshes_identical(got, want);
      }
    }

    // to end.
    {
      Arena ra(morph_target_buf, TARGET_HALF);
      Arena rb(morph_target_buf + TARGET_HALF, TARGET_HALF);
      PolyMesh want = Solids::simple_registry[e.to_node].generate(ra, rb);

      Arena oa(morph_temp_buf, TEMP_HALF);
      Arena ob(morph_temp_buf + TEMP_HALF, TEMP_HALF);
      PolyMesh got = run_edge_op(e, seed, oa, ob, e.t_to, e.twist_to);
      if (e.settle)
        got = MeshOps::relax(got, ob, oa, ConwayGraph::SETTLE_RELAX_ITERATIONS);

      switch (to_end_regime(e)) {
      case EndRegime::EXACT:
        hs_test::check_meshes_identical(got, want);
        break;
      case EndRegime::BAKED_RELAX:
        check_equal_within_relax_gate(got, want);
        break;
      case EndRegime::VERTEX_MATCH:
        check_equal_up_to_vertex_order(got, want, 1e-4f);
        break;
      case EndRegime::REGULAR:
        check_regular_form(got, want, 0.02f);
        break;
      case EndRegime::PAIR_COVER:
        check_pairwise_vertex_cover(got, want, 1e-4f);
        break;
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [endpoint] edge %d: %s -> %s\n", ei,
                  Solids::simple_registry[e.from_node].name,
                  Solids::simple_registry[e.to_node].name);
  }
}
