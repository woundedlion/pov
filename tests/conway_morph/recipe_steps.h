/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Recipe-step leg kinds: truncate, snub and relax legs, sampled on the chain
// prefixes the shipping recipes reach.
// ---------------------------------------------------------------------------

template <const Solids::Recipe &RECIPE, Solids::Op OP, size_t OCCURRENCE = 0>
inline PolyMesh recipe_step_seed(Arena &a, Arena &b) {
  constexpr size_t CAPACITY = Solids::lowered_step_count(RECIPE);
  Solids::OpStep lowered[CAPACITY];
  const size_t count = Solids::expand_to_primitives(RECIPE, lowered, CAPACITY);
  size_t prefix = 0;
  size_t occurrence = 0;
  while (prefix < count) {
    if (lowered[prefix].op == OP && occurrence++ == OCCURRENCE)
      break;
    ++prefix;
  }
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

/** @brief Recipe-step leg kind driven by the smoke test. */
enum class StepLegKind { TRUNCATE, SNUB, RELAX };

/** Paused redraws check_step_leg_smoke issues at its pause frame. */
inline constexpr int PAUSED_REDRAWS = 2;

/**
 * @brief Drives one recipe-step leg to completion through OpLeg, gating
 *        compiled-face-count constancy, ramp indices and per-frame motion.
 * @param kind Leg kind under test.
 * @param site Seed site and arrival parameter.
 * @param frames Leg length in frames.
 * @param max_step_chord Per-frame max-vertex-motion bound.
 * @param easing Easing driving the sweep parameter.
 * @param frames_out Optional sink receiving each drawn frame's vertices.
 * @param pause_after Leg frame after which two step_paused() redraws run; 0
 *        drives the leg unpaused.
 * @details The departed palettes come from the seed's own classification and
 * the bookend from the step's clean endpoint.
 */
inline void check_step_leg_smoke(
    StepLegKind kind, const StepLegSite &site, int frames, float max_step_chord,
    EasingFn easing = math::ease_in_out_sin,
    std::vector<std::vector<math::Vector>> *frames_out = nullptr,
    int pause_after = 0) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;

  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  hs::random().seed(2026u);

  Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh seed = build_step_leg_seed(site, leg_arena);
  PolyMesh endpoint;
  {
    constexpr size_t HALF = sizeof(morph_temp_buf) / 2;
    Arena ea(morph_temp_buf, HALF);
    Arena eb(morph_temp_buf + HALF, HALF);
    PolyMesh raw;
    switch (kind) {
    case StepLegKind::TRUNCATE:
      raw = MeshOps::truncate(seed, ea, eb, site.param);
      break;
    case StepLegKind::SNUB:
      raw = MeshOps::snub(seed, ea, eb, site.param, 0.0f);
      break;
    case StepLegKind::RELAX:
      raw = MeshOps::relax_baked(seed, ea, *site.bake);
      break;
    }
    endpoint = Solids::finalize_solid(raw, leg_arena);
  }

  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  MeshOps::classify_faces_by_topology(seed, scratch_arena_a, scratch_arena_b,
                                      leg_arena);
  MeshOps::classify_faces_by_topology(endpoint, scratch_arena_a,
                                      scratch_arena_b, leg_arena);

  std::array<int, OpLeg::PALETTES> slots;
  MeshPaletteBank::shuffle_indices(slots);

  const size_t prev_faces = seed.face_counts.size();
  uint8_t *prev_pal = leg_arena.allocate_n<uint8_t>(prev_faces);
  math::Vector *prev_centroid = leg_arena.allocate_n<math::Vector>(prev_faces);
  size_t off = 0;
  for (size_t f = 0; f < prev_faces; ++f) {
    prev_pal[f] = static_cast<uint8_t>(
        slots[math::wrap(static_cast<int>(seed.topology[f]), OpLeg::PALETTES)]);
    const int n = seed.face_counts[f];
    prev_centroid[f] = face_centroid_unit(seed, off, n);
    off += n;
  }

  OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                .prev_face_palette = prev_pal,
                                .prev_faces = prev_faces,
                                .prev_face_centroid = prev_centroid};
  OpLeg::BookendClasses bookend{.topology = endpoint.topology.data(),
                                .faces = endpoint.face_counts.size()};

  LegDrawProbe probe;
  auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
    probe.observe(m, sh);
    if (frames_out)
      frames_out->push_back(probe.prev_v);
  };

  hs_test::StubEffect fx(288, 144);

  auto run = [&](OpLeg &&leg) {
    const OpLeg::Landing &landing = leg.landing();
    probe.ramp_count = landing.blend_pairs;
    HS_EXPECT_EQ(landing.primary_faces, seed.face_counts.size());
    HS_EXPECT_TRUE(landing.topology != nullptr);
    for (int f = 0; f < frames; ++f) {
      {
        Canvas c(fx);
        leg.step(c);
      }
      fx.advance_display();
      if (f + 1 != pause_after)
        continue;
      for (int p = 0; p < PAUSED_REDRAWS; ++p) {
        {
          Canvas c(fx);
          leg.step_paused(c);
        }
        fx.advance_display();
      }
    }
    HS_EXPECT_EQ(probe.drawn,
                 (size_t)(frames + (pause_after > 0 ? PAUSED_REDRAWS : 0)));
    HS_EXPECT_EQ(probe.faces, landing.faces);
  };

  switch (kind) {
  case StepLegKind::TRUNCATE:
    run(OpLeg(seed,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                    .t_start = 0.0f,
                                    .t_end = site.param,
                                    .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  case StepLegKind::SNUB:
    run(OpLeg(seed,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::SNUB,
                                    .t_start = 0.0f,
                                    .t_end = site.param,
                                    .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  case StepLegKind::RELAX:
    run(OpLeg(seed, OpLeg::RelaxSpec{.bake = site.bake, .sweep_frames = frames},
              leg_arena, cb, handoff, bookend, OpLeg::classic_blend, easing));
    break;
  }

  HS_EXPECT_LT(probe.worst_step, max_step_chord);
  const char *label = kind == StepLegKind::TRUNCATE ? "truncate"
                      : kind == StepLegKind::SNUB   ? "snub"
                                                    : "relax";
  std::printf("  [opleg %s] %s: %zu faces, worst per-frame vertex step %.4f "
              "chord (bound %.2f)%s\n",
              label, site.name, probe.faces, (double)probe.worst_step,
              (double)max_step_chord,
              hs_test::stats().failed != failed_before ? " FAILED" : "");
}

/** Seed fixtures for kis and dual gated-swap legs. */
inline constexpr StepLegSite GATE_LEG_SITES[] = {
    {"icosahedron", probe_icosahedron, 0.0f},
    {"icosahedron_snub", probe_icosa_snub, 0.0f},
};

/**
 * @brief Drives one gated-swap leg to completion through OpLeg, gating the
 *        per-side compiled-face-count constancy, the single swap, the ramp
 *        indices.
 * @param op Partition operator the leg swaps to.
 * @param site Seed site.
 * @param gate Half-gate length in frames.
 * @details Derives this test-local handoff and bookend: the
 * departed palettes from the seed's own classification, the bookend from the
 * op's clean endpoint.
 */
inline void check_gated_leg_smoke(Animation::OpLeg::SwapOp op,
                                  const StepLegSite &site, int gate) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;
  const bool is_kis = op == OpLeg::SwapOp::KIS;

  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  hs::random().seed(2026u);

  Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  PolyMesh seed = build_step_leg_seed(site, leg_arena);
  PolyMesh endpoint;
  {
    constexpr size_t HALF = sizeof(morph_temp_buf) / 2;
    Arena ea(morph_temp_buf, HALF);
    Arena eb(morph_temp_buf + HALF, HALF);
    endpoint = Solids::finalize_solid(is_kis ? MeshOps::kis(seed, ea, eb)
                                             : MeshOps::dual(seed, ea, eb),
                                      leg_arena);
  }

  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  MeshOps::classify_faces_by_topology(seed, scratch_arena_a, scratch_arena_b,
                                      leg_arena);
  MeshOps::classify_faces_by_topology(endpoint, scratch_arena_a,
                                      scratch_arena_b, leg_arena);

  std::array<int, OpLeg::PALETTES> slots;
  MeshPaletteBank::shuffle_indices(slots);

  const size_t prev_faces = seed.face_counts.size();
  uint8_t *prev_pal = leg_arena.allocate_n<uint8_t>(prev_faces);
  for (size_t f = 0; f < prev_faces; ++f)
    prev_pal[f] = static_cast<uint8_t>(
        slots[math::wrap(static_cast<int>(seed.topology[f]), OpLeg::PALETTES)]);

  OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                .prev_face_palette = prev_pal,
                                .prev_faces = prev_faces};
  OpLeg::BookendClasses bookend{.topology = endpoint.topology.data(),
                                .faces = endpoint.face_counts.size()};

  const int frames = 2 * gate + 1;
  int drawn = 0, swaps = 0, swap_frame = -1;
  size_t side_faces = 0;
  int target_ramps = 0;
  std::array<bool, OpLeg::PALETTES> seed_palettes{};
  for (size_t f = 0; f < prev_faces; ++f)
    seed_palettes[prev_pal[f]] = true;
  const int SEED_RAMPS = static_cast<int>(
      std::count(seed_palettes.begin(), seed_palettes.end(), true));
  auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    for (size_t f = 0; f < sh.faces; ++f)
      HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]),
                   drawn < gate ? SEED_RAMPS : target_ramps);

    if (drawn == 0) {
      side_faces = sh.faces;
    } else if (sh.faces != side_faces) {
      ++swaps;
      swap_frame = drawn;
      side_faces = sh.faces;
    }
    ++drawn;
  };

  hs_test::StubEffect fx(288, 144);

  OpLeg leg(seed, OpLeg::GatedSwapSpec{.op = op, .gate_frames = gate},
            leg_arena, cb, handoff, bookend);
  const OpLeg::Landing &landing = leg.landing();
  target_ramps = landing.blend_pairs;
  for (int f = 0; f < frames; ++f) {
    {
      Canvas c(fx);
      leg.step(c);
    }
    fx.advance_display();
  }

  HS_EXPECT_EQ(drawn, frames);
  // Topology changes exactly once, at the swap frame, and the leg lands on the
  // arrival face count.
  HS_EXPECT_EQ(swaps, 1);
  HS_EXPECT_EQ(swap_frame, gate);
  HS_EXPECT_EQ(side_faces, landing.faces);
  HS_EXPECT_EQ(landing.faces, endpoint.face_counts.size());
  HS_EXPECT_EQ(landing.primary_faces, prev_faces);
  HS_EXPECT_TRUE(landing.from_palette != nullptr);

  std::printf("  [opleg %s] %s: F %zu -> %zu, swap at frame %d of %d%s\n",
              is_kis ? "kis" : "dual", site.name, prev_faces, landing.faces,
              swap_frame, frames,
              hs_test::stats().failed != failed_before ? " FAILED" : "");
}

/**
 * @brief Smoke-tests a gated-swap leg for both partition ops on every gate
 *        site, with a test-local 6-frame half-gate.
 */
inline void test_opleg_gated_swap_smoke() {
  using Animation::OpLeg;
  for (const StepLegSite &site : GATE_LEG_SITES) {
    check_gated_leg_smoke(OpLeg::SwapOp::KIS, site, 6);
    check_gated_leg_smoke(OpLeg::SwapOp::DUAL, site, 6);
  }
}

/**
 * @brief Smoke-tests the first truncate and snub sites and every baked relax site end to end.
 */
inline void test_opleg_step_leg_smoke() {
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0],
                       RecipeLegLengths::SWEEP_LEG_FRAMES, 0.15f);
  check_step_leg_smoke(StepLegKind::SNUB, SNUB_LEG_SITES[0],
                       RecipeLegLengths::SWEEP_LEG_FRAMES, 0.15f);
  for (const StepLegSite &site : RELAX_LEG_SITES)
    check_step_leg_smoke(StepLegKind::RELAX, site,
                         RecipeLegLengths::RELAX_LEG_FRAMES, 0.15f);
}

/**
 * @brief Pins a paused leg's redraw: the held frame, drawn again, with the
 *        sweep clock untouched.
 * @details Both paused draws must reproduce the held frame's vertices bitwise,
 * and the resumed frames must match the unpaused run frame for frame.
 */
inline void test_opleg_step_paused_holds_frame() {
  constexpr int FRAMES = 8;
  constexpr int PAUSE_AFTER = 4;
  // Loose chord bound: motion is not under test.
  constexpr float CHORD_MAX = 1.0f;
  std::vector<std::vector<math::Vector>> unpaused, held;
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, math::ease_in_out_sin, &unpaused);
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, math::ease_in_out_sin, &held, PAUSE_AFTER);
  HS_EXPECT_SIZE_OR_RETURN(unpaused, (size_t)FRAMES);
  HS_EXPECT_EQ(held.size(), (size_t)(FRAMES + PAUSED_REDRAWS));
  if (held.size() != (size_t)(FRAMES + PAUSED_REDRAWS))
    return;

  auto identical = [](const std::vector<math::Vector> &a,
                      const std::vector<math::Vector> &b) {
    if (a.size() != b.size() || a.empty())
      return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i].x != b[i].x || a[i].y != b[i].y || a[i].z != b[i].z)
        return false;
    return true;
  };

  // Every paused draw repeats the frame the leg was holding.
  for (int p = 0; p < PAUSED_REDRAWS; ++p)
    HS_EXPECT_TRUE(identical(held[PAUSE_AFTER - 1], held[PAUSE_AFTER + p]));

  // The clock never moved: the resumed frames are the unpaused run's.
  size_t drifted = 0;
  for (int f = 0; f < FRAMES; ++f) {
    const size_t k = f < PAUSE_AFTER ? (size_t)f : (size_t)(f + PAUSED_REDRAWS);
    if (!identical(unpaused[f], held[k]))
      ++drifted;
  }
  HS_EXPECT_EQ(drifted, (size_t)0);
  std::printf("  [opleg paused] truncate: %d frames, %d paused redraws at "
              "frame %d, %zu resumed frames drifted\n",
              FRAMES, PAUSED_REDRAWS, PAUSE_AFTER, drifted);
}

/**
 * @brief Drives swept legs under an easing whose range leaves [0, 1], gating
 *        the sweep parameter against extrapolation past the arrival.
 * @details ease_out_elastic's first overshoot lobe spans x in [0.075, 0.225],
 * peaking near 1.37. Frames 2-5 of 24 fall in it and clamp to the arrival (the
 * frame-4-vs-last comparison); frame 1 is still mid-sweep. The chord bound is
 * loose because elastic covers half the sweep in one frame.
 */
inline void test_opleg_step_leg_overshooting_easing() {
  constexpr StepLegSite NEAR_AMBO{"icosahedron_ambo_truncate049",
                                  probe_icosa_ambo, 0.49f};
  constexpr int FRAMES = 24;
  for (int k = 0; k < 2; ++k) {
    std::vector<std::vector<math::Vector>> drawn;
    const bool truncate = k == 0;
    check_step_leg_smoke(truncate ? StepLegKind::TRUNCATE : StepLegKind::SNUB,
                         truncate ? NEAR_AMBO : SNUB_LEG_SITES[0], FRAMES, 2.0f,
                         math::ease_out_elastic, &drawn);
    HS_EXPECT_SIZE_OR_RETURN(drawn, FRAMES);
    const std::vector<math::Vector> &arrival = drawn.back();
    const std::vector<math::Vector> &peak = drawn[3];
    const std::vector<math::Vector> &opening = drawn[0];
    HS_EXPECT_SIZE_OR_RETURN(peak, arrival.size());
    HS_EXPECT_SIZE_OR_RETURN(opening, arrival.size());
    float worst_peak = 0.0f, worst_opening = 0.0f;
    for (size_t i = 0; i < arrival.size(); ++i) {
      worst_peak =
          fold_worst(worst_peak, math::distance_between(peak[i], arrival[i]));
      worst_opening = fold_worst(
          worst_opening, math::distance_between(opening[i], arrival[i]));
    }
    HS_EXPECT_LT(worst_peak, 1e-5f);
    HS_EXPECT_GT(worst_opening, 1e-3f);
    std::printf("  [opleg elastic] %s: overshoot frame %.2e from arrival, "
                "opening frame %.4f\n",
                truncate ? "truncate" : "snub", (double)worst_peak,
                (double)worst_opening);
  }
}

/**
 * @brief Pins the crossfade colour model on a ConwayGraph edge leg under the
 *        classic_blend default: exact `from` at frame 1, moved off `from` by
 *        mid-leg, exact `to` at the arrival frame.
 */
inline void test_opleg_edge_leg_crossfade() {
  using Animation::OpLeg;
  reset_globals();
  const ScopedArenaSplit split(
      IslamicStars<288, 144>::GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_a,
      IslamicStars<288, 144>::GENERATED_BUDGET.scratch_b);
  hs::random().seed(2026u);

  Arena bank_arena(morph_bank_buf, sizeof(morph_bank_buf));
  MeshPaletteBank bank;
  bank.bake_all(bank_arena);

  hs_test::StubEffect fx(288, 144);

  // LUT grid-aligned sample coordinates for exact ramp-color comparisons.
  constexpr float SAMPLES[] = {0.0f, 0.5f, 1.0f};
  const OpLeg::Landing *lp = nullptr;
  std::vector<char> all_from, all_to; // per drawn frame: exact from/to palettes
  auto matches = [&](const OpLeg::Shading &sh, size_t f, uint8_t pal) {
    for (float t : SAMPLES) {
      const Color4 got = sh.ramps[sh.face_ramp[f]].get(t);
      const Color4 exp = bank.bank.entries[pal].get(t);
      if (got.color.r != exp.color.r || got.color.g != exp.color.g ||
          got.color.b != exp.color.b)
        return false;
    }
    return true;
  };
  auto cb = [&](Canvas &, const MeshState &, const OpLeg::Shading &sh) {
    bool from_ok = true, to_ok = true;
    for (size_t f = 0; f < sh.faces; ++f) {
      const uint8_t from = lp->from_palette[f];
      const uint8_t to = lp->to_palette[math::wrap(
          static_cast<int>(lp->topology[f]), OpLeg::PALETTES)];
      from_ok = from_ok && matches(sh, f, from);
      to_ok = to_ok && matches(sh, f, to);
    }
    all_from.push_back(from_ok);
    all_to.push_back(to_ok);
  };
  auto divergent_faces = [&]() {
    int n = 0;
    for (size_t f = 0; f < lp->faces; ++f)
      if (lp->from_palette[f] !=
          lp->to_palette[math::wrap(static_cast<int>(lp->topology[f]),
                                    OpLeg::PALETTES)])
        ++n;
    return n;
  };
  auto run_frames = [&](OpLeg &leg, int frames) {
    all_from.clear();
    all_to.clear();
    lp = &leg.landing();
    for (int f = 0; f < frames; ++f) {
      {
        Canvas c(fx);
        leg.step(c);
      }
      fx.advance_display();
    }
    HS_EXPECT_EQ(all_from.size(), (size_t)frames);
  };
  {
    const int failed_before = hs_test::stats().failed;
    Arena leg_arena(morph_target_buf, sizeof(morph_target_buf));
    PolyMesh cube;
    build_solid<Solids::Cube>(cube, leg_arena);
    uint8_t pal[16];
    for (size_t f = 0; f < cube.face_counts.size(); ++f)
      pal[f] = static_cast<uint8_t>(f % OpLeg::PALETTES);
    const int edge = find_directed_edge(ConwayGraph::EDGES, ConwayGraph::CUBE,
                                        ConwayGraph::TRUNCATED_CUBE);
    HS_EXPECT_GE(edge, 0);
    if (edge < 0)
      return;
    OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                  .prev_face_palette = pal,
                                  .prev_faces = cube.face_counts.size()};
    constexpr int EDGE_FRAMES = 24;
    OpLeg leg(cube,
              OpLeg::EdgeSweepSpec{.edge = &ConwayGraph::EDGES[edge],
                                   .sweep_frames = EDGE_FRAMES},
              leg_arena, cb, handoff);
    run_frames(leg, EDGE_FRAMES);
    HS_EXPECT_SIZE_OR_RETURN(all_from, EDGE_FRAMES);
    HS_EXPECT_SIZE_OR_RETURN(all_to, EDGE_FRAMES);
    HS_EXPECT_GT(divergent_faces(), 0);
    HS_EXPECT_TRUE(all_from[0]);
    HS_EXPECT_TRUE(!all_from[EDGE_FRAMES / 2 - 1]);
    HS_EXPECT_TRUE(all_to[EDGE_FRAMES - 1]);
    std::printf("  [crossfade] edge leg: F=%zu mid-leg blend intact%s\n",
                lp->faces,
                hs_test::stats().failed != failed_before ? " FAILED" : "");
  }
}

/**
 * @brief Pins recipe-step acceptance at the supported sweep boundaries.
 */
inline void test_unsweepable_recipe_steps_are_gated() {
  using Solids::Op;
  using Solids::IslamicStarPatterns::D2R;
  using ConwayGraph::T_TRUNCATE_ARRIVAL_MIN;
  using ConwayGraph::T_TRUNCATE_FAR_MAX;

  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.33f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.49f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::TRUNCATE, 0.5f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.51f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::SNUB, 0.5f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::RELAX, 8.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::HANKIN, 62.0f * D2R}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::AMBO}));
  HS_EXPECT_TRUE(
      Solids::is_morphable_step({Op::TRUNCATE, T_TRUNCATE_ARRIVAL_MIN}));
  HS_EXPECT_FALSE(Solids::is_morphable_step(
      {Op::TRUNCATE, std::nextafter(T_TRUNCATE_ARRIVAL_MIN, 0.0f)}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, T_TRUNCATE_FAR_MAX}));
  HS_EXPECT_FALSE(Solids::is_morphable_step(
      {Op::TRUNCATE, std::nextafter(T_TRUNCATE_FAR_MAX, 1.0f)}));
  HS_EXPECT_TRUE(
      Solids::is_morphable_step({Op::CHAMFER, Solids::CHAMFER_T_MAX}));
  HS_EXPECT_FALSE(Solids::is_morphable_step(
      {Op::CHAMFER, std::nextafter(Solids::CHAMFER_T_MAX, 1.0f)}));
  for (Op op : {Op::EXPAND, Op::BEVEL, Op::GYRO, Op::META, Op::NEEDLE, Op::ZIP})
    HS_EXPECT_FALSE(Solids::is_morphable_step({op}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::KIS}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::DUAL}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::CHAMFER, 0.001f}));
  // Zero snub inset, zero hankin angle and bake-less relax below one iteration
  // are not morphable.
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::SNUB, 0.0f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::HANKIN, 0.0f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::RELAX, 0.0f}));
  HS_EXPECT_TRUE(!Solids::is_morphable_step({Op::RELAX, 0.5f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::RELAX, 1.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step(
      {.op = Op::RELAX,
       .bake = &Solids::RelaxBakes::dodecahedron_bevel20_converged}));
}
