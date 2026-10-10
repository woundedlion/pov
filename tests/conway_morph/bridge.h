/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Bridge convergence: the tetrahedral edges that cross symmetry families.
// ---------------------------------------------------------------------------

/**
 * @brief Distinct hankin star-face classes of a mesh at interlace angle 0.
 * @param mesh Node base mesh.
 * @param a Scratch arena.
 * @param b Scratch arena.
 * @return Number of distinct classes over the star faces (the first
 *         face-count faces of the hankin mesh).
 */
inline int hankin_star_class_count(const PolyMesh &mesh, Arena &a, Arena &b) {
  CompiledHankin ch;
  MeshOps::compile_hankin(mesh, ch, a, b);
  MeshState hk;
  MeshOps::update_hankin(ch, hk, a, 0.0f);
  MeshOps::classify_faces_by_topology(hk, b, b, b);
  int classes = 0;
  for (size_t f = 0; f < mesh.face_counts.size(); ++f)
    classes = std::max(classes, static_cast<int>(hk.topology[f]) + 1);
  return classes;
}

/**
 * @brief Verifies every edge swept on the tetra -> icosa bridge's arrival
 *        classifies its hankin star faces like the registry node, so the
 *        icosahedral family adopted from the bridge colors symmetrically.
 */
inline void test_snub_bridge_icosahedron_classifies_like_registry() {
  constexpr size_t TARGET_HALF = sizeof(morph_target_buf) / 2;
  constexpr size_t TEMP_HALF = sizeof(morph_temp_buf) / 2;
  Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

  const ConwayGraph::EdgeSpec *bridge = nullptr;
  for (const auto &e : ConwayGraph::EDGES)
    if (e.from_node == ConwayGraph::TETRAHEDRON &&
        e.to_node == ConwayGraph::ICOSAHEDRON)
      bridge = &e;
  HS_EXPECT_TRUE(bridge != nullptr);
  if (!bridge)
    return;

  PolyMesh icosa;
  {
    Arena ta(morph_target_buf, sizeof(morph_target_buf));
    Arena tb(morph_temp_buf, sizeof(morph_temp_buf));
    PolyMesh tetra =
        Solids::simple_registry[ConwayGraph::TETRAHEDRON].generate(ta, tb);
    PolyMesh arrival =
        run_edge_op(*bridge, tetra, ta, tb, bridge->t_to, bridge->twist_to);
    MeshOps::clone(arrival, icosa, aux);
  }

  int checked = 0;
  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];
    if (e.seed_solid != ConwayGraph::ICOSAHEDRON)
      continue;
    for (const bool to_end : {false, true}) {
      Arena ra(morph_target_buf, TARGET_HALF);
      Arena rb(morph_target_buf + TARGET_HALF, TARGET_HALF);
      Arena oa(morph_temp_buf, TEMP_HALF);
      Arena ob(morph_temp_buf + TEMP_HALF, TEMP_HALF);
      const int node = to_end ? e.to_node : e.from_node;
      const float t = to_end ? e.t_to : e.t_from;
      PolyMesh got;
      if (t > 0.0f)
        got = run_edge_op(e, icosa, oa, ob, t,
                          to_end ? e.twist_to : e.twist_from);
      else
        MeshOps::clone(icosa, got, oa);
      PolyMesh want = Solids::simple_registry[node].generate(ra, rb);
      const int got_classes = hankin_star_class_count(got, ob, rb);
      const int want_classes = hankin_star_class_count(want, ob, rb);
      if (got_classes != want_classes)
        std::printf("    [bridge] edge %d %s end: %d hankin classes, registry "
                    "'%s' has %d\n",
                    ei, to_end ? "to" : "from", got_classes,
                    Solids::simple_registry[node].name, want_classes);
      HS_EXPECT_EQ(got_classes, want_classes);
      ++checked;
    }
  }
  HS_EXPECT_GT(checked, 0);
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
