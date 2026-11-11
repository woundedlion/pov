/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Recipe-step leg kinds: truncate, snub and relax legs, sampled on the chain
// prefixes the shipping recipes reach.
// ---------------------------------------------------------------------------

/** @brief One recipe-step leg site: the chain prefix the step sweeps on. */
struct StepLegSite {
  const char *name;                         /**< Diagnostic label. */
  PolyMesh (*seed)(Arena &a, Arena &b);     /**< Chain prefix up to the step. */
  float param;                              /**< Arrival t. */
  const MeshOps::RelaxBake *bake = nullptr; /**< Relax arrival bake. */
  bool via_dt_macro = false;     /**< Bridge follows the dt truncate leg. */
  double min_face_area = 1e-3;   /**< Floor on the bridge sweep's face area. */
  float seam_quiet_ratio = .70f; /**< Floor on seam-to-control quiet pixels. */
  size_t opening_retraced_faces =
      0; /**< Doubled arcs at the opening endpoint. */
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
    relax_leg_site<
        Solids::ICOSAHEDRON_TRUNCATE25_HK60_AMBO_RELAX100_HK60_RECIPE>(
        "icosahedron_truncate25_hk60_ambo"),
    relax_leg_site<Solids::DODECAHEDRON_HK35_AMBO_HK62_AMBO_RELAX_HK42_RECIPE>(
        "dodecahedron_hk35_ambo_hk62_ambo"),
};

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

/** @brief Recipe-step leg kind driven by the sweep and smoke tests. */
enum class StepLegKind { TRUNCATE, SNUB, RELAX };

/**
 * @brief Steps a truncate or snub sweep on each site, asserting constant raw
 *        and compiled face counts, two-face edge incidence, and Euler
 *        characteristic 2 at every sampled parameter.
 * @param kind `StepLegKind::TRUNCATE` or `StepLegKind::SNUB`.
 * @param sites Seed sites and arrival parameters.
 */
inline void check_step_leg_sweep(StepLegKind kind,
                                 std::span<const StepLegSite> sites) {
  constexpr int SAMPLES = 33;
  const bool truncate = kind == StepLegKind::TRUNCATE;
  HS_EXPECT(truncate || kind == StepLegKind::SNUB,
            "sweep takes truncate or snub");
  const char *tag = truncate ? "truncate-leg" : "snub-leg";
  const ConwayGraph::MorphOp op =
      truncate ? ConwayGraph::MorphOp::TRUNCATE : ConwayGraph::MorphOp::SNUB;

  for (const StepLegSite &site : sites) {
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
          check_sweep_sample(truncate ? MeshOps::truncate(seed, a, b, t)
                                      : MeshOps::snub(seed, a, b, t, 0.0f),
                             a, b);
      if (s == 0) {
        first = fp;
        expect_op_counts(fp.v, fp.f, fp.i, morph_op_counts(op, seed));
      } else {
        expect_same_fingerprint(fp, first);
      }
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [%s] %s failed (raw F=%zu, compiled F=%zu)\n", tag,
                  site.name, first.f, first.compiled);
    else
      std::printf("  [%s] %s: t*=%.2f F=%zu compiled=%zu across %d "
                  "samples\n",
                  tag, site.name, (double)site.param, first.f, first.compiled,
                  SAMPLES);
  }
}

/** @brief Steps truncate sweeps on TRUNCATE_LEG_SITES. */
inline void test_truncate_leg_on_recipe_seeds_holds_topology() {
  check_step_leg_sweep(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES);
}

/** @brief Steps snub sweeps on SNUB_LEG_SITES. */
inline void test_snub_leg_on_recipe_seeds_holds_topology() {
  check_step_leg_sweep(StepLegKind::SNUB, SNUB_LEG_SITES);
}

/** @brief Baked recipe relaxes preserve vertex count and order. */
inline void test_relax_recipe_seeds_preserve_vertex_order() {
  for (const StepLegSite &site : RELAX_LEG_SITES) {
    HS_CONTEXT(site.name);
    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed = build_step_leg_seed(site, persist);
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    PolyMesh relaxed = MeshOps::relax_baked(seed, a, *site.bake);
    check_vertex_order_identity(relaxed, seed);
  }
}

/** Paused redraws check_step_leg_smoke issues at its pause frame. */
inline constexpr int PAUSED_REDRAWS = 2;

/** @brief Arena-owned seed, endpoint and palette handoff for leg smoke checks. */
struct LegFixture {
  ScopedArenaSplit split =
      islamic_split(IslamicStars<288, 144>::GENERATED_BUDGET);
  Arena leg_arena{morph_target_buf, sizeof(morph_target_buf)};
  Arena bank_arena{morph_bank_buf, sizeof(morph_bank_buf)};
  ScratchScope a_guard{scratch_arena_a};
  ScratchScope b_guard{scratch_arena_b};
  MeshPaletteBank bank;
  PolyMesh seed, endpoint;
  std::array<int, Animation::OpLeg::PALETTES> slots;
  size_t prev_faces = 0;
  uint8_t *prev_pal = nullptr;

  template <typename BuildEndpoint>
  LegFixture(const StepLegSite &site, BuildEndpoint build_endpoint) {
    hs::random().seed(2026u);
    bank.bake_all(bank_arena);
    seed = build_step_leg_seed(site, leg_arena);
    constexpr size_t HALF = sizeof(morph_temp_buf) / 2;
    Arena ea(morph_temp_buf, HALF);
    Arena eb(morph_temp_buf + HALF, HALF);
    endpoint = Solids::finalize_solid(build_endpoint(seed, ea, eb), leg_arena);
    MeshOps::classify_faces_by_topology(seed, scratch_arena_a, scratch_arena_b,
                                        leg_arena);
    MeshOps::classify_faces_by_topology(endpoint, scratch_arena_a,
                                        scratch_arena_b, leg_arena);
    MeshPaletteBank::shuffle_indices(slots);
    prev_faces = seed.face_counts.size();
    prev_pal = leg_arena.allocate_n<uint8_t>(prev_faces);
    for (size_t face = 0; face < prev_faces; ++face)
      prev_pal[face] = static_cast<uint8_t>(
          slots[MeshPaletteBank::slot_of(seed.topology[face])]);
  }
};

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
 * @param weights_out Optional sink receiving each drawn frame's crossfade weight.
 * @details The departed palettes come from the seed's own classification and
 * the bookend from the step's clean endpoint.
 */
inline void check_step_leg_smoke(
    StepLegKind kind, const StepLegSite &site, int frames, float max_step_chord,
    EasingFn easing = math::ease_in_out_sin,
    std::vector<std::vector<math::Vector>> *frames_out = nullptr,
    int pause_after = 0, std::vector<float> *weights_out = nullptr) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;

  LegFixture fixture(site, [&](const PolyMesh &seed, Arena &ea, Arena &eb) {
    switch (kind) {
    case StepLegKind::TRUNCATE:
      return MeshOps::truncate(seed, ea, eb, site.param);
    case StepLegKind::SNUB:
      return MeshOps::snub(seed, ea, eb, site.param, 0.0f);
    case StepLegKind::RELAX:
      return MeshOps::relax_baked(seed, ea, *site.bake);
    }
    return PolyMesh{};
  });
  auto &leg_arena = fixture.leg_arena;
  auto &bank = fixture.bank;
  const auto &seed = fixture.seed;
  const auto &endpoint = fixture.endpoint;
  const size_t prev_faces = fixture.prev_faces;
  const auto *prev_pal = fixture.prev_pal;
  math::Vector *prev_centroid = leg_arena.allocate_n<math::Vector>(prev_faces);
  size_t off = 0;
  for (size_t f = 0; f < prev_faces; ++f) {
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
    if (weights_out)
      weights_out->push_back(sh.weight);
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
 * @brief A gated leg preserves each side's face count, swaps once and emits
 * the expected ramp indices.
 * @param op Partition operator the leg swaps to.
 * @param site Seed site.
 * @param gate Half-gate length in frames.
 */
inline void check_gated_leg_smoke(Animation::OpLeg::SwapOp op,
                                  const StepLegSite &site, int gate) {
  using Animation::OpLeg;
  const int failed_before = hs_test::stats().failed;
  const bool is_kis = op == OpLeg::SwapOp::KIS;

  LegFixture fixture(site, [&](const PolyMesh &seed, Arena &ea, Arena &eb) {
    return is_kis ? MeshOps::kis(seed, ea, eb) : MeshOps::dual(seed, ea, eb);
  });
  auto &leg_arena = fixture.leg_arena;
  auto &bank = fixture.bank;
  const auto &seed = fixture.seed;
  const auto &endpoint = fixture.endpoint;
  const size_t prev_faces = fixture.prev_faces;
  const auto *prev_pal = fixture.prev_pal;

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
  for (size_t f = 0; f < prev_faces; ++f) {
    HS_EXPECT_LT(prev_pal[f], OpLeg::PALETTES);
    if (prev_pal[f] < OpLeg::PALETTES)
      seed_palettes[prev_pal[f]] = true;
  }
  const int seed_ramps = static_cast<int>(
      std::count(seed_palettes.begin(), seed_palettes.end(), true));
  float last_weight = 0.0f;
  auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
    HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
    for (size_t f = 0; f < sh.faces; ++f)
      HS_EXPECT_LT(static_cast<int>(sh.face_ramp[f]),
                   drawn < gate ? seed_ramps : target_ramps);
    HS_EXPECT_EQ(sh.bank == nullptr, drawn < gate);
    const float previous = drawn == 0 ? sh.weight : last_weight;
    HS_EXPECT_EQ(sh.previous_weight, previous);
    last_weight = sh.weight;

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
  hs_test::drive_leg(leg, fx, frames);

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
  constexpr float SWEEP_MAX_STEP_CHORD = 0.04f;
  constexpr float RELAX_MAX_STEP_CHORD = 0.04f;
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0],
                       RecipeLegLengths::SWEEP_LEG_FRAMES,
                       SWEEP_MAX_STEP_CHORD);
  check_step_leg_smoke(StepLegKind::SNUB, SNUB_LEG_SITES[0],
                       RecipeLegLengths::SWEEP_LEG_FRAMES,
                       SWEEP_MAX_STEP_CHORD);
  for (const StepLegSite &site : RELAX_LEG_SITES)
    check_step_leg_smoke(StepLegKind::RELAX, site,
                         RecipeLegLengths::RELAX_LEG_FRAMES,
                         RELAX_MAX_STEP_CHORD);
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
  std::vector<float> unpaused_weights, held_weights;
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, math::ease_in_out_sin, &unpaused, 0,
                       &unpaused_weights);
  check_step_leg_smoke(StepLegKind::TRUNCATE, TRUNCATE_LEG_SITES[0], FRAMES,
                       CHORD_MAX, math::ease_in_out_sin, &held, PAUSE_AFTER,
                       &held_weights);
  HS_EXPECT_SIZE_OR_RETURN(unpaused, (size_t)FRAMES);
  HS_EXPECT_EQ(held.size(), (size_t)(FRAMES + PAUSED_REDRAWS));
  if (held.size() != (size_t)(FRAMES + PAUSED_REDRAWS))
    return;
  HS_EXPECT_SIZE_OR_RETURN(unpaused_weights, unpaused.size());
  HS_EXPECT_SIZE_OR_RETURN(held_weights, held.size());

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
  for (int p = 0; p < PAUSED_REDRAWS; ++p) {
    HS_EXPECT_TRUE(identical(held[PAUSE_AFTER - 1], held[PAUSE_AFTER + p]));
    HS_EXPECT_EQ(held_weights[PAUSE_AFTER + p], held_weights[PAUSE_AFTER - 1]);
  }

  // The clock never moved: the resumed frames are the unpaused run's.
  size_t drifted = 0;
  for (int f = 0; f < FRAMES; ++f) {
    const size_t k = f < PAUSE_AFTER ? (size_t)f : (size_t)(f + PAUSED_REDRAWS);
    if (!identical(unpaused[f], held[k]))
      ++drifted;
    HS_EXPECT_EQ(unpaused_weights[f], held_weights[k]);
  }
  HS_EXPECT_EQ(drifted, (size_t)0);
  std::printf("  [opleg paused] truncate: %d frames, %d paused redraws at "
              "frame %d, %zu resumed frames drifted\n",
              FRAMES, PAUSED_REDRAWS, PAUSE_AFTER, drifted);
}

/**
 * @brief Drives swept legs under an easing whose range leaves [0, 1], gating
 *        the sweep parameter against extrapolation past the arrival.
 * @details Frames inside ease_out_elastic's first overshoot lobe clamp to the
 * arrival. The chord bound is loose because elastic covers half the sweep in
 * one frame.
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
 * @brief Pins recipe-step acceptance at the supported sweep boundaries.
 */
inline void test_recipe_step_morphable_boundaries() {
  using Solids::Op;
  using Solids::IslamicStarPatterns::D2R;
  using ConwayGraph::T_TRUNCATE_ARRIVAL_MIN;
  using ConwayGraph::T_TRUNCATE_FAR_MAX;

  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.33f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, 0.49f}));
  constexpr float NEAR_PINCH = 0.5f - ConwayGraph::T_EPS_AMBO;
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::TRUNCATE, NEAR_PINCH}));
  HS_EXPECT_FALSE(Solids::is_morphable_step(
      {Op::TRUNCATE, std::nextafter(NEAR_PINCH, 0.5f)}));
  HS_EXPECT_FALSE(Solids::is_morphable_step({Op::TRUNCATE, 0.5f}));
  HS_EXPECT_TRUE(
      Solids::is_morphable_step({Op::TRUNCATE, std::nextafter(0.5f, 1.0f)}));
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
  HS_EXPECT_FALSE(Solids::is_morphable_step({Op::CHAMFER, 0.001f}));
  // Zero snub inset, zero hankin angle and bake-less relax below one iteration
  // are not morphable.
  HS_EXPECT_FALSE(Solids::is_morphable_step({Op::SNUB, 0.0f}));
  HS_EXPECT_FALSE(Solids::is_morphable_step({Op::HANKIN, 0.0f}));
  HS_EXPECT_FALSE(Solids::is_morphable_step({Op::RELAX, 0.0f}));
  HS_EXPECT_FALSE(Solids::is_morphable_step({Op::RELAX, 0.5f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step({Op::RELAX, 1.0f}));
  HS_EXPECT_TRUE(Solids::is_morphable_step(
      {.op = Op::RELAX,
       .bake = &Solids::RelaxBakes::dodecahedron_bevel20_converged}));
}
