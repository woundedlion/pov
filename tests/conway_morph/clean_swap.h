/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Clean-swap invisibility: the boundary swaps exchange geometrically
// matching meshes.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies truncate(seed, 0.5 - eps) vertices pairwise merge onto
 *        ambo(seed) vertices for every sweep seed: each ambo vertex has
 *        exactly 2 truncate vertices within tolerance.
 */
inline void test_truncate_near_half_merges_onto_ambo() {
  constexpr float MERGE_TOL = 0.02f;
  for (MorphSeed s : MORPH_SEEDS) {
    HS_CONTEXT(seed_name(s));
    Arena target(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
    Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

    PolyMesh seed = build_morph_seed(s, aux, temp);
    PolyMesh tr =
        MeshOps::truncate(seed, target, temp, 0.5f - ConwayGraph::T_EPS_AMBO);
    PolyMesh am = MeshOps::ambo(seed, temp, target);

    check_pairwise_vertex_cover(tr, am, MERGE_TOL);
  }
}

/**
 * @brief Verifies each parameterized op at t = T_EPS emits primary faces that
 *        geometrically match the seed's faces, for every sweep seed.
 * @details truncate contributes two cut corners per seed corner; expand,
 *          snub (zero twist) and chamfer contribute one inset corner each.
 *          Tolerances bound the T_EPS displacement plus the unit-sphere
 *          renormalization.
 */
inline void test_ops_at_t_eps_primary_faces_match_seed() {
  for (MorphSeed s : MORPH_SEEDS) {
    HS_CONTEXT(seed_name(s));
    auto check = [&](auto op, int corners, float tolerance) {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      const PolyMesh seed = build_morph_seed(s, aux, temp);
      const PolyMesh out = op(seed, target, temp);
      check_primary_faces_match_seed(seed, out, corners, tolerance);
    };
    check(
        [](const PolyMesh &seed, Arena &target, Arena &temp) {
          return MeshOps::truncate(seed, target, temp, T_EPS);
        },
        2, PRIMARY_CORNER_TOL_TRUNCATE);
    check(
        [](const PolyMesh &seed, Arena &target, Arena &temp) {
          return MeshOps::expand(seed, target, temp, T_EPS);
        },
        1, PRIMARY_CORNER_TOL_SINGLE);
    check(
        [](const PolyMesh &seed, Arena &target, Arena &temp) {
          return MeshOps::snub(seed, target, temp, T_EPS, 0.0f);
        },
        1, PRIMARY_CORNER_TOL_SINGLE);
    check(
        [](const PolyMesh &seed, Arena &target, Arena &temp) {
          return MeshOps::chamfer(seed, target, temp, T_EPS);
        },
        1, PRIMARY_CORNER_TOL_SINGLE);
  }
}
