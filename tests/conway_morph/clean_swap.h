/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// §7.5 Clean-swap invisibility: the boundary swaps exchange geometrically
// matching meshes.
// ---------------------------------------------------------------------------

/**
 * @brief Verifies truncate(seed, 0.5 - eps) vertices pairwise merge onto
 *        ambo(seed) vertices for every sweep seed: each ambo vertex has
 *        exactly 2 truncate vertices within tolerance.
 */
inline void test_truncate_near_half_merges_onto_ambo() {
  constexpr float MERGE_EPS = 1e-3f;
  constexpr float MERGE_TOL = 1e-2f;
  for (MorphSeed s : MORPH_SEEDS) {
    Arena target(morph_target_buf, sizeof(morph_target_buf));
    Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
    Arena aux(morph_aux_buf, sizeof(morph_aux_buf));

    PolyMesh seed = build_morph_seed(s, aux, temp);
    PolyMesh tr = MeshOps::truncate(seed, target, temp, 0.5f - MERGE_EPS);
    PolyMesh am = MeshOps::ambo(seed, temp, target);

    HS_CONTEXT(seed_name(s));
    check_pairwise_vertex_cover(tr, am, MERGE_TOL);
  }
}

/**
 * @brief Verifies each parameterized op at t = T_EPS emits primary faces that
 *        geometrically match the seed's faces, for every sweep seed.
 * @details truncate contributes two cut corners per seed corner; expand and
 *          snub (zero twist) contribute one inset corner each. Tolerances
 *          bound the T_EPS displacement plus the unit-sphere renormalization.
 */
inline void test_ops_at_t_eps_primary_faces_match_seed() {
  for (MorphSeed s : MORPH_SEEDS) {
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::truncate(seed, target, temp, T_EPS);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 2,
                                     PRIMARY_CORNER_TOL_TRUNCATE);
    }
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::expand(seed, target, temp, T_EPS);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 1,
                                     PRIMARY_CORNER_TOL_SINGLE);
    }
    {
      Arena target(morph_target_buf, sizeof(morph_target_buf));
      Arena temp(morph_temp_buf, sizeof(morph_temp_buf));
      Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
      PolyMesh seed = build_morph_seed(s, aux, temp);
      PolyMesh out = MeshOps::snub(seed, target, temp, T_EPS, 0.0f);
      check_primary_faces_match_seed(seed, out, /*corners_per_source*/ 1,
                                     PRIMARY_CORNER_TOL_SINGLE);
    }
  }
}
