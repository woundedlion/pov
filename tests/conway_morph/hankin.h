/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// MeshOps::hankin topology probe over sampled recipe prefixes.
// ---------------------------------------------------------------------------

/** @brief One hankin-sweep leg seed and its arrival angle. */
struct HankinSweepSite {
  const char *name;                     /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b); /**< Chain prefix up to the hankin. */
  float theta_star;                     /**< Arrival contact angle, radians. */
};

template <const Solids::Recipe &RECIPE, size_t OCCURRENCE>
inline constexpr HankinSweepSite hankin_sweep_site(const char *name) {
  size_t occurrence = 0;
  for (size_t i = 0; i < RECIPE.count; ++i)
    if (RECIPE.steps[i].op == Solids::Op::HANKIN && occurrence++ == OCCURRENCE)
      return {name, recipe_step_seed<RECIPE, Solids::Op::HANKIN, OCCURRENCE>,
              RECIPE.steps[i].param};
  throw "hankin_sweep_site: recipe has no matching HANKIN step";
}

inline constexpr HankinSweepSite HANKIN_SWEEP_SITES[] = {
    hankin_sweep_site<Solids::DODECAHEDRON_HK62_AMBO_HK62_RECIPE, 0>(
        "dodecahedron"),
    hankin_sweep_site<Solids::DODECAHEDRON_HK62_AMBO_HK62_RECIPE, 1>(
        "dodecahedron_hk62_ambo"),
    hankin_sweep_site<Solids::OCTAHEDRON_HK17_AMBO_HK73_RECIPE, 0>(
        "octahedron"),
    hankin_sweep_site<Solids::OCTAHEDRON_HK17_AMBO_HK73_RECIPE, 1>(
        "octahedron_hk17_ambo"),
};

/**
 * @brief Checks topology over the sampled dodecahedron, dodecahedron_hk62_ambo,
 *        octahedron and octahedron_hk17_ambo prefixes.
 * @details Asserts constant raw and compiled face counts, two-face edge incidence
 *        and Euler characteristic 2 from THETA_EPS to each arrival angle.
 */
inline void test_hankin_sweep_on_islamic_seeds_holds_topology() {
  constexpr int SAMPLES = 17;
  constexpr float THETA_EPS = Animation::OpLeg::THETA_EPS;

  for (const HankinSweepSite &site : HANKIN_SWEEP_SITES) {
    const int failed_before = hs_test::stats().failed;

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(site.seed(ga, gb), persist);
    }

    size_t v0 = 0, f0 = 0, i0 = 0, compiled0 = 0;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float theta =
          THETA_EPS + (site.theta_star - THETA_EPS) *
                          (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept = MeshOps::hankin(seed, a, b, theta);
      MeshState compiled;
      MeshOps::compile(swept, compiled, a, b);
      if (s == 0) {
        v0 = swept.vertices.size();
        f0 = swept.face_counts.size();
        i0 = swept.faces.size();
        compiled0 = compiled.face_counts.size();
        expect_op_counts(v0, f0, i0, hankin_op_counts(seed));
      } else {
        HS_EXPECT_EQ(swept.vertices.size(), v0);
        HS_EXPECT_EQ(swept.face_counts.size(), f0);
        HS_EXPECT_EQ(swept.faces.size(), i0);
        HS_EXPECT_EQ(compiled.face_counts.size(), compiled0);
      }
      check_face_counts_consistent(swept);
      check_indices_in_range(swept);
      check_all_unit_vertices(swept, 1e-3f);
      conway_tests::check_euler_characteristic_two(swept);
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [hankin-sweep] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, f0, compiled0);
    else
      std::printf("  [hankin-sweep] %s: F=%zu compiled=%zu across %d samples\n",
                  site.name, f0, compiled0, SAMPLES);
  }
}
