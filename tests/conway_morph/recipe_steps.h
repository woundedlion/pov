/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Recipe-step leg kinds: truncate, snub and relax legs, sampled on the chain
// prefixes the shipping recipes reach.
// ---------------------------------------------------------------------------

template <const Solids::Recipe &RECIPE, Solids::Op OP>
inline PolyMesh recipe_step_seed(Arena &a, Arena &b) {
  constexpr size_t CAPACITY = Solids::lowered_step_count(RECIPE);
  Solids::OpStep lowered[CAPACITY];
  const size_t count = Solids::expand_to_primitives(RECIPE, lowered, CAPACITY);
  size_t prefix = 0;
  while (prefix < count && lowered[prefix].op != OP)
    ++prefix;
  HS_EXPECT_LT(prefix, count);
  return Solids::build_steps(RECIPE.seed, lowered, prefix, a, b);
}

/** @brief One recipe-step leg site: the chain prefix the step sweeps on. */
struct StepLegSite {
  const char *name;                         /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b);     /**< Chain prefix up to the step. */
  float param;                              /**< Arrival t. */
  const MeshOps::RelaxBake *bake = nullptr; /**< Relax arrival bake. */
  bool via_dt_macro = false; /**< Bridge follows the dt truncate leg. */
  double min_face_area = 1e-3;
  float seam_quiet_ratio = .70f;
};

inline PolyMesh probe_icosahedron(Arena &a, Arena &b) {
  return Solids::Platonic::icosahedron(a, b);
}
inline PolyMesh probe_icosa_ambo(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::icosahedron(a, b), a, b)
      .ambo()
      .build();
}
inline PolyMesh probe_icosa_snub(Arena &a, Arena &b) {
  return Solids::SolidBuilder(Solids::Platonic::icosahedron(a, b), a, b)
      .snub()
      .build();
}
inline PolyMesh probe_icosa_snub_relax(Arena &a, Arena &b) {
  return recipe_step_seed<
      Solids::ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_RECIPE,
      Solids::Op::TRUNCATE>(a, b);
}
inline PolyMesh probe_ticosa_ambo_relax_converged(Arena &a, Arena &b) {
  return recipe_step_seed<
      Solids::TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_RECIPE,
      Solids::Op::TRUNCATE>(a, b);
}

/** @brief Representative chain prefixes for truncate-leg sweeps. */
inline constexpr StepLegSite TRUNCATE_LEG_SITES[] = {
    {"icosahedron_ambo", probe_icosa_ambo, 0.33f},
    {"truncatedIcosahedron_ambo_relax_converged",
     probe_ticosa_ambo_relax_converged, 0.33f},
    {"icosahedron_snub_relax", probe_icosa_snub_relax, 0.33f},
};

/** @brief Representative snub-leg seed at the snub() defaults. */
inline constexpr StepLegSite SNUB_LEG_SITES[] = {
    {"icosahedron", probe_icosahedron, 0.5f},
};

template <const Solids::Recipe &RECIPE>
inline constexpr StepLegSite relax_leg_site(const char *name) {
  for (size_t i = 0; i < RECIPE.count; ++i)
    if (RECIPE.steps[i].op == Solids::Op::RELAX)
      return {name, recipe_step_seed<RECIPE, Solids::Op::RELAX>, 0.0f,
              RECIPE.steps[i].bake};
  throw "relax_leg_site: recipe has no RELAX step";
}

/** Unique baked-relax seed sites in the shipping recipes. */
inline constexpr StepLegSite RELAX_LEG_SITES[] = {
    relax_leg_site<Solids::ICOSAHEDRON_SNUB_RELAX_TRUNCATE033_HANKIN62_RECIPE>(
        "icosahedron_snub"),
    relax_leg_site<Solids::DODECAHEDRON_AMBO_BEVEL33_RELAX_HK66_RECIPE>(
        "dodecahedron_ambo_bevel33"),
    relax_leg_site<
        Solids::TRUNCATED_ICOSAHEDRON_AMBO_RELAX_TRUNCATE33_HK64_RECIPE>(
        "truncatedIcosahedron_ambo"),
    relax_leg_site<Solids::DODECAHEDRON_BEVEL2_RELAX_GYRO_RECIPE>(
        "dodecahedron_bevel20"),
    relax_leg_site<
        Solids::TRUNCATED_ICOSIDODECAHEDRON_BEVEL5_RELAX_HK77_RECIPE>(
        "truncatedIcosidodecahedron_bevel50"),
    relax_leg_site<Solids::DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE>(
        "dodecahedron_hk35_ambo_hk62_ambo"),
};

/** @brief Topology fingerprint of one sweep sample. */
struct SweepFingerprint {
  size_t v = 0;        /**< Raw vertex count. */
  size_t f = 0;        /**< Raw face count. */
  size_t i = 0;        /**< Raw face-index count. */
  size_t compiled = 0; /**< Compiled face count. */
};

/**
 * @brief Fingerprints a swept mesh and runs the structural checks.
 * @param swept Sweep sample.
 * @param a Arena receiving the compiled mesh.
 * @param b Compile scratch arena.
 * @return The sample's fingerprint.
 */
inline SweepFingerprint check_sweep_sample(const PolyMesh &swept, Arena &a,
                                           Arena &b) {
  MeshState compiled;
  MeshOps::compile(swept, compiled, a, b);
  check_face_counts_consistent(swept);
  check_indices_in_range(swept);
  check_all_unit_vertices(swept, 1e-3f);
  conway_tests::check_euler_characteristic_two(swept);
  return {swept.vertices.size(), swept.face_counts.size(), swept.faces.size(),
          compiled.face_counts.size()};
}

/**
 * @brief Asserts a sweep sample fingerprint matches the opening one.
 * @param s Sample fingerprint.
 * @param first Opening-sample fingerprint.
 */
inline void expect_same_fingerprint(const SweepFingerprint &s,
                                    const SweepFingerprint &first) {
  HS_EXPECT_EQ(s.v, first.v);
  HS_EXPECT_EQ(s.f, first.f);
  HS_EXPECT_EQ(s.i, first.i);
  HS_EXPECT_EQ(s.compiled, first.compiled);
}

/**
 * @brief Builds a leg site's seed into a persistent arena.
 * @param site Site whose chain prefix is built.
 * @param persist Arena receiving the finalized seed.
 * @return The seed mesh in @p persist.
 */
inline PolyMesh build_step_leg_seed(const StepLegSite &site, Arena &persist) {
  constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
  Arena ga(morph_aux_buf, HALF);
  Arena gb(morph_aux_buf + HALF, HALF);
  PolyMesh seed = site.seed(ga, gb);
  if (site.via_dt_macro)
    seed = MeshOps::truncate(seed, gb, ga, RECONCILE_TRUNCATE_T);
  return Solids::finalize_solid(seed, persist);
}

/**
 * @brief Steps truncate sweeps on TRUNCATE_LEG_SITES, asserting
 *        constant raw and compiled face counts, two-face edge incidence, and Euler
 *        characteristic 2 at every sampled parameter.
 */
inline void test_truncate_leg_on_recipe_seeds_holds_topology() {
  constexpr int SAMPLES = 33;

  for (const StepLegSite &site : TRUNCATE_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed = build_step_leg_seed(site, persist);

    SweepFingerprint first;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float t = T_EPS + (site.param - T_EPS) *
                                  (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      const SweepFingerprint fp =
          check_sweep_sample(MeshOps::truncate(seed, a, b, t), a, b);
      if (s == 0) {
        first = fp;
        expect_op_counts(fp.v, fp.f, fp.i,
                         morph_op_counts(ConwayGraph::MorphOp::TRUNCATE, seed));
      } else {
        expect_same_fingerprint(fp, first);
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [truncate-leg] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, first.f, first.compiled);
    else
      std::printf("  [truncate-leg] %s: t*=%.2f F=%zu compiled=%zu across %d "
                  "samples\n",
                  site.name, (double)site.param, first.f, first.compiled,
                  SAMPLES);
  }
}

/**
 * @brief Steps snub sweeps on SNUB_LEG_SITES, asserting constant
 *        raw and compiled face counts, two-face edge incidence, and Euler
 *        characteristic 2 at every sampled parameter.
 */
inline void test_snub_leg_on_recipe_seeds_holds_topology() {
  constexpr int SAMPLES = 33;

  for (const StepLegSite &site : SNUB_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed = build_step_leg_seed(site, persist);

    SweepFingerprint first;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    for (int s = 0; s < SAMPLES; ++s) {
      const float t = T_EPS + (site.param - T_EPS) *
                                  (static_cast<float>(s) / (SAMPLES - 1));
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      const SweepFingerprint fp =
          check_sweep_sample(MeshOps::snub(seed, a, b, t, 0.0f), a, b);
      if (s == 0) {
        first = fp;
        expect_op_counts(fp.v, fp.f, fp.i,
                         morph_op_counts(ConwayGraph::MorphOp::SNUB, seed));
      } else {
        expect_same_fingerprint(fp, first);
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [snub-leg] %s failed (raw F=%zu, compiled F=%zu)\n",
                  site.name, first.f, first.compiled);
    else
      std::printf("  [snub-leg] %s: t*=%.2f F=%zu compiled=%zu across %d "
                  "samples\n",
                  site.name, (double)site.param, first.f, first.compiled,
                  SAMPLES);
  }
}

/**
 * @brief Steps a baked-relax slerp on every unique shipping recipe seed,
 *        asserting the vertex count/order identity the kind rests on plus
 *        constant compiled face counts and unit vertices across the slerp.
 */
inline void test_relax_leg_on_recipe_seeds_holds_topology() {
  constexpr int SAMPLES = 33;

  for (const StepLegSite &site : RELAX_LEG_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed = build_step_leg_seed(site, persist);

    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    PolyMesh relaxed = MeshOps::relax_baked(seed, a, *site.bake);

    // The precondition of a standalone relax leg: same vertex count, same
    // topology bytes, and vertex i still nearest its own seed vertex, so the
    // leg slerps per-vertex with no correspondence pass.
    HS_EXPECT_EQ(relaxed.vertices.size(), seed.vertices.size());
    HS_EXPECT_EQ(relaxed.face_counts.size(), seed.face_counts.size());
    HS_EXPECT_EQ(relaxed.faces.size(), seed.faces.size());
    if (relaxed.vertices.size() != seed.vertices.size() ||
        relaxed.face_counts.size() != seed.face_counts.size() ||
        relaxed.faces.size() != seed.faces.size())
      continue;
    HS_EXPECT_EQ(std::memcmp(relaxed.face_counts.data(),
                             seed.face_counts.data(),
                             seed.face_counts.size() * sizeof(uint8_t)),
                 0);
    HS_EXPECT_EQ(std::memcmp(relaxed.faces.data(), seed.faces.data(),
                             seed.faces.size() * sizeof(uint16_t)),
                 0);
    for (size_t i = 0; i < relaxed.vertices.size(); ++i) {
      size_t nearest = 0;
      float best = 1e9f;
      for (size_t j = 0; j < seed.vertices.size(); ++j) {
        const float d =
            math::distance_between(relaxed.vertices[i], seed.vertices[j]);
        if (d < best) {
          best = d;
          nearest = j;
        }
      }
      HS_EXPECT_EQ(nearest, i);
    }

    size_t first_compiled = 0;
    for (int s = 0; s < SAMPLES; ++s) {
      const float K = static_cast<float>(s) / (SAMPLES - 1);
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept;
      MeshOps::clone(seed, swept, a);
      for (size_t i = 0; i < swept.vertices.size(); ++i) {
        swept.vertices[i] =
            math::slerp(seed.vertices[i], relaxed.vertices[i], K);
        HS_EXPECT_NEAR(math::dot(swept.vertices[i], swept.vertices[i]), 1.0f,
                       1e-3f);
      }
      MeshState compiled;
      MeshOps::compile(swept, compiled, a, b);
      if (s == 0)
        first_compiled = compiled.face_counts.size();
      else
        HS_EXPECT_EQ(compiled.face_counts.size(), first_compiled);
    }

    std::printf("  [relax-leg] %s: baked compiled=%zu across %d samples%s\n",
                site.name, first_compiled, SAMPLES,
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}
