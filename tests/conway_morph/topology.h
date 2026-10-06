/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// §7.2 Topology-constancy sweep: connectivity is fixed on the open interval,
// so classification and palette assignment can hoist to once per leg.
// ---------------------------------------------------------------------------

/** Samples per edge sweep. */
constexpr int SWEEP_SAMPLES = 16;

/**
 * @brief Per-leg accumulator for the invariants every OpLeg draw callback
 *        shares: a compiled face count matching the shading and constant
 *        across the leg, in-range ramp indices, and a vertex list that keeps
 *        its count and order frame to frame.
 */
struct LegDrawProbe {
  int ramp_count = 0;      /**< LUT array length for this leg. */
  size_t drawn = 0;        /**< Frames handed to the callback. */
  size_t faces = 0;        /**< Face count latched on the first frame. */
  float worst_step = 0.0f; /**< Largest per-vertex motion between frames. */
  std::vector<math::Vector> prev_v; /**< Previous frame's vertices. */

  /**
   * @brief Folds one drawn frame in.
   * @param m Compiled mesh handed to the callback.
   * @param sh Shading handed alongside it.
   */
  void observe(const MeshState &m, const Animation::OpLeg::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    if (drawn == 0)
      faces = sh.faces;
    else
      HS_EXPECT_EQ(sh.faces, faces);
    for (size_t f = 0; f < sh.faces; ++f)
      HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]), ramp_count);
    if (!prev_v.empty()) {
      HS_EXPECT_EQ(m.vertices.size(), prev_v.size());
      for (size_t i = 0; i < m.vertices.size(); ++i)
        worst_step = fold_worst(
            worst_step, math::distance_between(m.vertices[i], prev_v[i]));
    }
    prev_v.assign(m.vertices.begin(), m.vertices.end());
    ++drawn;
  }
};

/** Raw mesh counts a Conway operator is required to produce. */
struct OpCounts {
  size_t v; /**< Vertices. */
  size_t f; /**< Faces. */
  size_t i; /**< Flat face indices; on a closed mesh this is 2E. */
};

/**
 * @brief Reads a closed mesh's operand counts.
 * @param m Closed manifold mesh.
 * @return {V, F, I} of @p m.
 */
inline OpCounts mesh_op_counts(const PolyMesh &m) {
  return {m.vertices.size(), m.face_counts.size(), m.faces.size()};
}

/**
 * @brief Counts a closed mesh's degree-2 vertices.
 * @param m Closed manifold mesh; each index is one face incidence.
 * @return Vertices with exactly two incident faces.
 * @details Truncating one produces a 2-gon the op drops, so this is the only
 *          term the Conway forms need beyond (V, E, F). Hankin seeds carry
 *          degree-2 star points; the registry solids carry none.
 */
inline size_t degree_two_vertex_count(const PolyMesh &m) {
  std::vector<uint32_t> degree(m.vertices.size(), 0u);
  for (size_t k = 0; k < m.faces.size(); ++k)
    ++degree[m.faces[k]];
  size_t n = 0;
  for (uint32_t d : degree)
    n += d == 2u ? 1u : 0u;
  return n;
}

/**
 * @brief Closed-form output counts of each Conway operator the graph sweeps.
 * @param op Operator applied.
 * @param seed Seed mesh.
 * @return The counts @p op must produce from @p seed, at every sweep parameter
 *         on the open interval.
 * @details Standard Conway arithmetic in the seed's (V, E, F): truncate
 *          (2E, F+V, 3E), expand (2E, F+V+E, 4E), snub (2E, F+V+2E, 5E),
 *          chamfer (V+2E, F+E, 4E) as (V', F', E'); the returned index count
 *          is 2E'. Truncate additionally drops the 2-gon each degree-2 vertex
 *          would raise, costing one face and one edge apiece.
 */
inline OpCounts morph_op_counts(ConwayGraph::MorphOp op, const PolyMesh &seed) {
  const size_t v = seed.vertices.size();
  const size_t f = seed.face_counts.size();
  const size_t e = seed.faces.size() / 2;
  switch (op) {
  case ConwayGraph::MorphOp::TRUNCATE: {
    const size_t d2 = degree_two_vertex_count(seed);
    return {2 * e, f + v - d2, 6 * e - 2 * d2};
  }
  case ConwayGraph::MorphOp::EXPAND:
    return {2 * e, f + v + e, 8 * e};
  case ConwayGraph::MorphOp::SNUB:
    return {2 * e, f + v + 2 * e, 10 * e};
  case ConwayGraph::MorphOp::CHAMFER:
    return {v + 2 * e, f + e, 8 * e};
  }
  return {v, f, seed.faces.size()};
}

/**
 * @brief Closed-form output counts of MeshOps::hankin, at any contact angle.
 * @param seed Seed mesh.
 * @return {3E, F+V, 8E}.
 * @details compile_hankin lays down one midpoint per edge and one star point
 *          per half-edge (3E vertices), one star face per seed face and one
 *          rosette per seed vertex (F+V), and each face contributes twice its
 *          degree in indices on both sides (4I = 8E). Degree-2 vertices raise
 *          a quad rosette like any other, so no correction applies.
 */
inline OpCounts hankin_op_counts(const PolyMesh &seed) {
  const size_t e = seed.faces.size() / 2;
  return {3 * e, seed.face_counts.size() + seed.vertices.size(), 8 * e};
}

/**
 * @brief Asserts a sweep sample's raw counts equal a derived expectation.
 * @param v Sample vertex count.
 * @param f Sample face count.
 * @param i Sample flat index count.
 * @param d Expected counts.
 */
inline void expect_op_counts(size_t v, size_t f, size_t i, const OpCounts &d) {
  HS_EXPECT_EQ(v, d.v);
  HS_EXPECT_EQ(f, d.f);
  HS_EXPECT_EQ(i, d.i);
}

/**
 * @brief Returns the production clamp of both edge sweep endpoints.
 * @param e Edge to clamp.
 * @param t_lo Out: clamped t_from.
 * @param t_hi Out: clamped t_to.
 */
inline void edge_sweep_interval(const ConwayGraph::EdgeSpec &e, float &t_lo,
                                float &t_hi) {
  t_lo = ConwayGraph::clamp_edge_endpoint(e, e.t_from);
  t_hi = ConwayGraph::clamp_edge_endpoint(e, e.t_to);
}

/**
 * @brief Verifies every edge holds constant topology across its sweep:
 *        fixed V/F/I, two-face edge incidence, Euler characteristic 2,
 *        all faces >= 3 sides, near-unit vertices, no traps.
 * @details Snub twist and clamped t interpolate with the same leg progress.
 */
inline void test_edge_sweeps_hold_topology() {
  constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];
    const int failed_before = hs_test::stats().failed;

    Arena sa(morph_aux_buf, HALF);
    Arena sb(morph_aux_buf + HALF, HALF);
    PolyMesh seed = Solids::simple_registry[e.seed_solid].generate(sa, sb);

    float t_lo, t_hi;
    edge_sweep_interval(e, t_lo, t_hi);

    size_t v0 = 0, f0 = 0, i0 = 0;
    for (int s = 0; s < SWEEP_SAMPLES; ++s) {
      const float u = static_cast<float>(s) / (SWEEP_SAMPLES - 1);
      const float t = t_lo + (t_hi - t_lo) * u;
      const float twist = e.twist_from + (e.twist_to - e.twist_from) * u;

      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      PolyMesh out = run_edge_op(e, seed, target, temp, t, twist);

      if (s == 0) {
        v0 = out.vertices.size();
        f0 = out.face_counts.size();
        i0 = out.faces.size();
        expect_op_counts(v0, f0, i0, morph_op_counts(e.op, seed));
      } else {
        HS_EXPECT_EQ(out.vertices.size(), v0);
        HS_EXPECT_EQ(out.face_counts.size(), f0);
        HS_EXPECT_EQ(out.faces.size(), i0);
      }
      for (size_t fi = 0; fi < out.face_counts.size(); ++fi)
        HS_EXPECT_TRUE(out.face_counts[fi] >= 3);
      check_face_counts_consistent(out);
      check_indices_in_range(out);
      check_all_unit_vertices(out, 1e-3f);
      conway_tests::check_euler_characteristic_two(out);
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [sweep] edge %d: %s -> %s\n", ei,
                  Solids::simple_registry[e.from_node].name,
                  Solids::simple_registry[e.to_node].name);
  }
}
