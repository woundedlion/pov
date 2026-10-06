/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// Recipe chain build replay (docs/specs/opchain_morph_spec.md, "Validation
// contract"): test-local partition chains are lowered and replayed leg by leg
// using individual OpLegs. DUAL/KIS use gated swaps here; RecipeBuild instead
// schedules bridges/macros, so those chains do not establish production budgets.
// test_opchain_arena_survey.h covers the registry.
// Each leg is stepped frame by frame with a recording draw callback.
// Gates per-leg compiled-face-count constancy, the crossfade colour model
// (every frame of every leg draws every face's (from, to) ramp bit-exact at
// the leg's blend weight; each leg departs from the palette the previous leg
// landed on), the final per-face sprite handoff, the per-final-class palette
// symmetry of the finished shape, and the
// persistent/scratch high-water against IslamicStars' configured split.
// ---------------------------------------------------------------------------

/** The arena split is canvas-independent; this instantiation names it. */
using IslamicFx = IslamicStars<288, 144>;

constexpr size_t ISLAMIC_SCRATCH_A_BUDGET =
    IslamicFx::RECIPE_BUDGET.scratch_a; /**< IslamicStars recipe scratch_a. */
constexpr size_t ISLAMIC_SCRATCH_B_BUDGET =
    IslamicFx::RECIPE_BUDGET.scratch_b; /**< IslamicStars build scratch_b. */
/** Device persistent budget of IslamicStars' arena split. */
constexpr size_t ISLAMIC_PERSISTENT_BUDGET =
    IslamicFx::RECIPE_BUDGET.device_persistent();
/** Scratch capacity the replay runs with, above every budget it gates. */
constexpr size_t REPLAY_SCRATCH_CAPACITY = 512 * 1024;

/**
 * @brief Arena high-waters and shape of one replayed chain.
 */
struct ChainPeaks {
  size_t persistent = 0;  /**< persistent_arena high-water, bytes. */
  size_t scratch_a = 0;   /**< scratch_arena_a high-water, bytes. */
  size_t scratch_b = 0;   /**< scratch_arena_b high-water, bytes. */
  size_t legs = 0;        /**< Lowered primitive step count. */
  size_t faces = 0;       /**< Face count of the finished solid. */
  int palettes = 0;       /**< Distinct palettes on the finished shape. */
  int final_classes = 0;  /**< Distinct newborn classes on the final leg. */
  int blend_pairs = 0;    /**< Max distinct (from, to) pairs on any leg. */
  bool supported = false; /**< Every lowered step has a replay leg kind. */
  bool production_schedule =
      true; /**< No approximated DUAL/KIS bridge scheduling. */
};

/**
 * @brief Replays recipe legs and gates continuity and classification.
 * @details DUAL/KIS are gated-swap approximations of production bridges/macros;
 *          arena budgets are checked only for chains without those steps.
 * @param name Diagnostic label.
 * @param recipe Recipe replayed.
 * @return The chain's arena high-waters.
 */
inline ChainPeaks replay_build_chain(const char *name,
                                     const Solids::Recipe &recipe) {
  ChainPeaks peaks;
  using Animation::OpLeg;
  constexpr int GATE_HALF_FRAMES = 6;
  constexpr size_t MAX_FACES = 1152;
  constexpr size_t MAX_STEPS = 8;

  {
    const int failed_before = hs_test::stats().failed;

    reset_globals();
    // Capacities are the host's, not the device split: a chain that overruns a
    // budget must report its high-water, not OOM-trap the replay before the
    // measurement. The budgets below are what the peaks are gated against.
    const ScopedArenaSplit split(
        GLOBAL_ARENA_SIZE - 2 * REPLAY_SCRATCH_CAPACITY,
        REPLAY_SCRATCH_CAPACITY, REPLAY_SCRATCH_CAPACITY);
    hs::random().seed(2026u);

    MeshPaletteBank bank;
    bank.bake_all(persistent_arena);

    // Seed solid: persistent PolyMesh plus the compiled/classified carousel
    // slot, exactly as spawn_shape prepares them.
    PolyMesh cur;
    MeshState seed_slot;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      cur = Solids::finalize_solid(
          Solids::simple_registry[recipe.seed].generate(a, b), target);
      MeshOps::compile(cur, seed_slot, target, a);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(seed_slot, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // The spawned seed's shuffled palette order, consumed by class ordinal.
    std::array<int, OpLeg::PALETTES> slots;
    MeshPaletteBank::shuffle_indices(slots);
    std::array<uint8_t, OpLeg::PALETTES> order;
    for (int i = 0; i < OpLeg::PALETTES; ++i)
      order[i] = static_cast<uint8_t>(slots[i]);

    Solids::OpStep steps[MAX_STEPS];
    const size_t count = Solids::expand_to_primitives(recipe, steps, MAX_STEPS);
    HS_EXPECT_GT(count, (size_t)0);
    if (count == 0)
      return peaks;
    bool supported = true;
    int leg_frames[MAX_STEPS] = {};
    bool gated[MAX_STEPS] = {};
    for (size_t k = 0; k < count; ++k) {
      if (!Solids::is_morphable_step(steps[k]))
        supported = false;
      gated[k] =
          steps[k].op == Solids::Op::KIS || steps[k].op == Solids::Op::DUAL;
      peaks.production_schedule &= !gated[k];
      leg_frames[k] = gated[k] ? 2 * GATE_HALF_FRAMES + 1
                      : steps[k].op == Solids::Op::HANKIN
                          ? RecipeLegLengths::HANKIN_LEG_FRAMES
                      : steps[k].op == Solids::Op::RELAX
                          ? RecipeLegLengths::RELAX_LEG_FRAMES
                          : RecipeLegLengths::SWEEP_LEG_FRAMES;
    }
    peaks.legs = count;
    peaks.supported = supported;
    HS_EXPECT_TRUE(supported);
    if (!supported) {
      std::printf("    [chain] %s: unsupported lowered op\n", name);
      return peaks;
    }

    hs_test::StubEffect fx(288, 144);
    uint8_t prev_pal_buf[MAX_FACES];
    math::Vector prev_centroid[MAX_FACES];
    uint8_t carried_to[MAX_FACES] = {};
    std::vector<int> full_topo;
    const OpLeg::Landing *prev_landing = nullptr;
    PolyMesh next;
    // A hankin leg's endpoint is rebuilt into scratch_a and stays there until
    // the boundary evacuates it into the fresh persistent arena, as
    // finish_build_leg does. One is live at a time, so each rebuild rewinds to
    // this mark first.
    const size_t endpoint_mark = scratch_arena_a.get_offset();

    for (size_t k = 0; k < count; ++k) {
      const size_t prev_faces = cur.face_counts.size();
      HS_EXPECT_LE(prev_faces, MAX_FACES);
      // The replay threads state leg to leg, so an over-cap leg cannot be
      // skipped: abandon the chain instead of writing past the handoff arrays.
      if (prev_faces > MAX_FACES)
        return peaks;
      {
        size_t off = 0;
        for (size_t f = 0; f < prev_faces; ++f) {
          const int n = cur.face_counts[f];
          prev_centroid[f] = face_centroid_unit(cur, off, n);
          off += n;
        }
      }
      const uint8_t *prev_pal;
      if (k == 0) {
        // Seed keying, as spawn_shape does it: class ids are dense, so class
        // c is the c-th cohort.
        for (size_t f = 0; f < prev_faces; ++f) {
          const int cls = seed_slot.topology[f];
          prev_pal_buf[f] =
              static_cast<uint8_t>(order[math::wrap(cls, OpLeg::PALETTES)]);
        }
        prev_pal = prev_pal_buf;
      } else {
        for (size_t f = 0; f < prev_faces; ++f)
          prev_pal_buf[f] = carried_to[f];
        prev_pal = prev_pal_buf;
      }

      // Eager clean endpoint, exactly as start_build_leg derives it (for the
      // ambo leg: the AMBO mesh, never the swept truncate form). A hankin step
      // builds none: its leg rebuilds its own arrival once it has run.
      const bool hankin_step = steps[k].op == Solids::Op::HANKIN;
      OpLeg::BookendClasses bookend;
      if (!hankin_step) {
        hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
          PolyMesh nx;
          switch (steps[k].op) {
          case Solids::Op::AMBO:
            nx = MeshOps::ambo(cur, a, b);
            break;
          case Solids::Op::TRUNCATE:
            nx = MeshOps::truncate(cur, a, b, steps[k].param);
            break;
          case Solids::Op::SNUB:
            nx = MeshOps::snub(cur, a, b, steps[k].param, steps[k].twist);
            break;
          case Solids::Op::KIS:
            nx = MeshOps::kis(cur, a, b);
            break;
          case Solids::Op::DUAL:
            nx = MeshOps::dual(cur, a, b);
            break;
          case Solids::Op::CHAMFER:
            nx = MeshOps::chamfer(cur, a, b, steps[k].param);
            break;
          default:
            nx = steps[k].bake
                     ? MeshOps::relax_baked(cur, a, *steps[k].bake)
                     : MeshOps::relax(cur, a, b,
                                      static_cast<int>(steps[k].param));
            break;
          }
          next = Solids::finalize_solid(nx, target);
        });
        {
          ScratchScope a_guard(scratch_arena_a);
          ScratchScope b_guard(scratch_arena_b);
          MeshOps::classify_faces_by_topology(
              next, scratch_arena_a, scratch_arena_b, persistent_arena);
        }
        const size_t bookend_faces = next.face_counts.size();
        HS_EXPECT_LE(bookend_faces, MAX_FACES);
        bookend = {.topology = next.topology.data(), .faces = bookend_faces};
      }

      OpLeg::PaletteHandoff handoff{.bank = &bank.bank,
                                    .prev_face_palette = prev_pal,
                                    .prev_faces = prev_faces,
                                    .prev_face_centroid = prev_centroid};

      size_t drawn = 0;
      size_t leg_faces = 0;
      int off_palette_frames = 0;
      const OpLeg::Landing *lp = nullptr;
      // LUT grid-aligned sample coordinates for exact ramp-color comparisons.
      constexpr float PROBE_T[] = {0.0f, 0.5f, 1.0f};
      Arena blend_arena(morph_temp_buf, sizeof(morph_temp_buf));
      // A gated leg's face count is constant per side and changes once, at the
      // swap; every other kind holds one count for the whole leg. Every frame
      // must draw every face's (from, to) ramp bit-exact at the frame's blend
      // weight: endpoint weights alias the bank LUTs, interior weights are
      // rebuilt here through the same bake_palette_blend the leg uses.
      auto cb = [&](Canvas &, const MeshState &m, const OpLeg::Shading &sh) {
        HS_EXPECT_EQ(m.face_counts.size(), sh.faces);
        const bool side_start =
            drawn == 0 ||
            (gated[k] && drawn == static_cast<size_t>(GATE_HALF_FRAMES));
        if (side_start)
          leg_faces = sh.faces;
        else
          HS_EXPECT_EQ(sh.faces, leg_faces);
        const bool seed_side =
            gated[k] && drawn < static_cast<size_t>(GATE_HALF_FRAMES);
        const int frame = static_cast<int>(drawn) + 1;
        const float w = gated[k] ? OpLeg::trailing_blend(frame, leg_frames[k])
                                 : OpLeg::classic_blend(frame, leg_frames[k]);
        blend_arena.reset();
        BakedPalette expected[OpLeg::PALETTES][OpLeg::PALETTES];
        bool baked[OpLeg::PALETTES][OpLeg::PALETTES] = {};
        bool frame_ok = true;
        for (size_t f = 0; f < sh.faces && lp; ++f) {
          const uint8_t from =
              seed_side ? prev_pal_buf[f] : lp->from_palette[f];
          const uint8_t to = seed_side ? from : lp->landed_palette(f);
          const BakedPalette *want = &bank.bank.entries[to].view();
          if (from != to) {
            if (!baked[from][to]) {
              expected[from][to] =
                  bake_palette_blend(blend_arena, bank.bank.entries[from],
                                     bank.bank.entries[to], w);
              baked[from][to] = true;
            }
            want = &expected[from][to];
          }
          for (float t : PROBE_T) {
            const Color4 got = sh.ramps[sh.face_ramp[f]].get(t);
            const Color4 exp = want->get(t);
            if (got.color.r != exp.color.r || got.color.g != exp.color.g ||
                got.color.b != exp.color.b)
              frame_ok = false;
          }
        }
        if (!frame_ok)
          ++off_palette_frames;
        ++drawn;
      };

      auto make_leg = [&]() {
        switch (steps[k].op) {
        case Solids::Op::KIS:
          return OpLeg(cur,
                       OpLeg::GatedSwapSpec{.op = OpLeg::SwapOp::KIS,
                                            .gate_frames = GATE_HALF_FRAMES},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::DUAL:
          return OpLeg(cur,
                       OpLeg::GatedSwapSpec{.op = OpLeg::SwapOp::DUAL,
                                            .gate_frames = GATE_HALF_FRAMES},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::HANKIN:
          return OpLeg(cur,
                       OpLeg::HankinSweepSpec{.theta_start = 0.0f,
                                              .theta_end = steps[k].param,
                                              .sweep_frames = leg_frames[k]},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::AMBO:
          return OpLeg(
              cur,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                    .t_start = 0.0f,
                                    .t_end = 0.5f,
                                    .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        case Solids::Op::TRUNCATE:
          return OpLeg(
              cur,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                    .t_start = 0.0f,
                                    .t_end = steps[k].param,
                                    .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        case Solids::Op::SNUB:
          return OpLeg(cur,
                       OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::SNUB,
                                             .t_start = 0.0f,
                                             .t_end = steps[k].param,
                                             .twist_end = steps[k].twist,
                                             .sweep_frames = leg_frames[k]},
                       persistent_arena, cb, handoff, bookend);
        case Solids::Op::CHAMFER:
          return OpLeg(
              cur,
              OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::CHAMFER,
                                    .t_start = 0.0f,
                                    .t_end = steps[k].param,
                                    .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        default:
          return OpLeg(
              cur,
              OpLeg::RelaxSpec{.iterations = static_cast<int>(steps[k].param),
                               .bake = steps[k].bake,
                               .sweep_frames = leg_frames[k]},
              persistent_arena, cb, handoff, bookend);
        }
      };
      OpLeg leg = make_leg();
      const OpLeg::Landing &landing = leg.landing();
      lp = &landing;

      for (int f = 0; f < leg_frames[k]; ++f) {
        {
          Canvas c(fx);
          leg.step(c);
        }
        fx.advance_display();
      }
      HS_EXPECT_EQ(drawn, (size_t)leg_frames[k]);
      HS_EXPECT_EQ(landing.faces, leg_faces);
      HS_EXPECT_TRUE(landing.from_palette != nullptr);
      // Every frame drew every face's (from, to) ramp bit-exact.
      HS_EXPECT_EQ(off_palette_frames, 0);

      // The landed carry: a surviving face departs the next leg from the
      // palette this leg landed it on. A partition op keeps no emission-order
      // correspondence to its seed, so its from-palettes are the leg's own
      // geometric provenance.
      HS_EXPECT_LE(landing.faces, MAX_FACES);
      if (landing.faces > MAX_FACES)
        return peaks;
      if (k > 0 && !gated[k]) {
        // Report the first offending face across the replay's full face sweep.
        int first_broken_carry = -1;
        for (size_t f = 0; f < prev_faces; ++f)
          if ((int)landing.from_palette[f] != (int)carried_to[f] &&
              first_broken_carry < 0)
            first_broken_carry = (int)f;
        HS_EXPECT_EQ(first_broken_carry, -1);
      }
      {
        // The leg's fresh target set is a permutation of the bank.
        bool seen_to[OpLeg::PALETTES] = {};
        for (int i = 0; i < OpLeg::PALETTES; ++i) {
          HS_EXPECT_LE((int)landing.to_palette[i], OpLeg::PALETTES - 1);
          seen_to[landing.to_palette[i]] = true;
        }
        for (bool s : seen_to)
          HS_EXPECT_TRUE(s);
        // Final-leg birth grain: distinct newborn classes (the perceptual
        // grouping report).
        if (k + 1 == count) {
          const size_t births_from = gated[k] ? 0 : prev_faces;
          int classes = 0;
          int prev_cls = -1;
          for (;;) {
            int cls = -1;
            for (size_t f = births_from; f < landing.faces; ++f) {
              const int c = landing.topology[f];
              if (c > prev_cls && (cls < 0 || c < cls))
                cls = c;
            }
            if (cls < 0)
              break;
            ++classes;
            prev_cls = cls;
          }
          peaks.final_classes = classes;
        }
      }
      peaks.blend_pairs = std::max(peaks.blend_pairs, landing.blend_pairs);
      for (size_t f = 0; f < landing.faces; ++f)
        carried_to[f] = landing.landed_palette(f);

      // Landing lives in the leg's arena-backed Transients; outlives `leg`.
      prev_landing = &landing;

      if (hankin_step) {
        next = PolyMesh();
        scratch_arena_a.set_offset(endpoint_mark);
        OpLeg::arrival_mesh(landing, next, scratch_arena_a);
      }

      // Full-precision final classification for the symmetry pin, rebuilt as
      // the leg's constructor builds its arrival (before the snorm16 pack).
      // Test-local arenas keep the replay's gated peaks untouched; a
      // non-hankin final leg lands on the clean endpoint, whose full-precision
      // classification the closing compile below provides.
      if (k + 1 == count && hankin_step) {
        Arena pa(morph_target_buf, sizeof(morph_target_buf));
        Arena pb(morph_temp_buf, sizeof(morph_temp_buf));
        Arena pc(morph_aux_buf, sizeof(morph_aux_buf));
        CompiledHankin hk;
        MeshOps::compile_hankin(cur, hk, pa, pb);
        PolyMesh full;
        MeshOps::update_hankin(hk, full, pa,
                               std::max(steps[k].param, OpLeg::THETA_EPS));
        MeshOps::classify_faces_by_topology(full, pb, pc, pa);
        full_topo.assign(full.topology.begin(), full.topology.end());
      }

      // Mirror finish_build_leg's boundary compaction: the finished leg's
      // transients are reclaimed and only the endpoint the next leg sweeps
      // from crosses the reset. Without it the replay accumulates every leg
      // and reports peaks the effect never reaches. The last leg's landing is
      // consumed by the final swap below, so it is left standing exactly as
      // the effect leaves it for the next shape's compaction.
      if (k + 1 < count) {
        Persist<PolyMesh> pn(next, scratch_arena_b, persistent_arena);
        cur = PolyMesh();
        seed_slot = MeshState();
        prev_landing = nullptr;
        persistent_arena.reset();
        bank.bake_all(persistent_arena);
      }
      cur = std::move(next);
    }

    // Final sprite handoff, mirroring finish_build: the finished solid's
    // per-face palettes are the last landing's landed palettes, snapshotted
    // before the closing compaction. The compiled slot must match them face
    // for face, so the sprite's first frame draws exactly what the last leg
    // frame drew (the closing w = 1 plateau).
    const size_t landed_faces = cur.face_counts.size();
    HS_EXPECT_LE(landed_faces, prev_landing->faces);
    HS_EXPECT_LE(landed_faces, MAX_FACES);
    if (landed_faces > MAX_FACES)
      return peaks;
    uint8_t sprite_pal[MAX_FACES];
    for (size_t f = 0; f < landed_faces; ++f)
      sprite_pal[f] = prev_landing->landed_palette(f);
    prev_landing = nullptr;
    MeshState final_slot;
    {
      ScratchScope a_guard(scratch_arena_a);
      {
        ScratchScope seed_guard(scratch_arena_b);
        PolyMesh built;
        MeshOps::clone(cur, built, scratch_arena_b);
        cur = PolyMesh();
        seed_slot = MeshState();
        persistent_arena.reset();
        bank.bake_all(persistent_arena);
        MeshOps::compile(built, final_slot, persistent_arena, scratch_arena_a);
      }
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(final_slot, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }
    {
      HS_EXPECT_EQ(landed_faces, final_slot.topology.size());
      int first_broken_carry = -1;
      for (size_t f = 0; f < landed_faces; ++f)
        if ((int)sprite_pal[f] != (int)carried_to[f] && first_broken_carry < 0)
          first_broken_carry = (int)f;
      HS_EXPECT_EQ(first_broken_carry, -1);
    }
    if (full_topo.empty())
      full_topo.assign(final_slot.topology.begin(), final_slot.topology.end());

    // Variety: distinct palettes on the finished shape.
    peaks.palettes = 0;
    {
      bool seen[OpLeg::PALETTES] = {};
      for (size_t f = 0; f < landed_faces; ++f)
        if (sprite_pal[f] < OpLeg::PALETTES)
          seen[sprite_pal[f]] = true;
      for (bool s : seen)
        peaks.palettes += s;
    }

    // Symmetry pin: on the finished shape any two faces with the same
    // full-precision final class wear the same palette — the landed palette
    // is a function of the final classification alone. A
    // quantized-classification regression shatters classes and breaks this.
    {
      HS_EXPECT_EQ(full_topo.size(), landed_faces);
      int asymmetric = 0;
      if (full_topo.size() == landed_faces) {
        for (size_t f = 0; f < landed_faces; ++f)
          for (size_t g = f + 1; g < landed_faces; ++g)
            if (full_topo[f] == full_topo[g] && sprite_pal[f] != sprite_pal[g])
              ++asymmetric;
      }
      HS_EXPECT_EQ(asymmetric, 0);
    }

    peaks.persistent = persistent_arena.get_high_water_mark();
    peaks.scratch_a = scratch_arena_a.get_high_water_mark();
    peaks.scratch_b = scratch_arena_b.get_high_water_mark();
    peaks.faces = final_slot.face_counts.size();
    std::printf("  [chain] %s: %zu legs, %d/%d palettes, final leg %d "
                "classes, %d/%d blend pairs, persistent=%zu B / %zu B, "
                "scratch a=%zu B / %zu B, b=%zu B / %zu B\n",
                name, count, peaks.palettes, OpLeg::PALETTES,
                peaks.final_classes, peaks.blend_pairs, OpLeg::MAX_BLEND_PAIRS,
                peaks.persistent, (size_t)ISLAMIC_PERSISTENT_BUDGET,
                peaks.scratch_a, (size_t)ISLAMIC_SCRATCH_A_BUDGET,
                peaks.scratch_b, (size_t)ISLAMIC_SCRATCH_B_BUDGET);
    if (peaks.production_schedule) {
      HS_EXPECT_LE(peaks.persistent, ISLAMIC_PERSISTENT_BUDGET);
      HS_EXPECT_LE(peaks.scratch_a, ISLAMIC_SCRATCH_A_BUDGET);
      HS_EXPECT_LE(peaks.scratch_b, ISLAMIC_SCRATCH_B_BUDGET);
    } else {
      std::printf(
          "    [chain] gated-swap approximation; production bridge budget excluded\n");
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [chain] %s FAILED\n", name);
  }
  return peaks;
}

/** Partition chains replayed as gated swaps, without production bridges/macros.
 * Registry partition recipes are surveyed in test_opchain_arena_survey.h. */
inline constexpr Solids::OpStep CHAIN_KIS[] = {{Solids::Op::KIS}};
inline constexpr Solids::OpStep CHAIN_KIS_DUAL[] = {{Solids::Op::KIS},
                                                    {Solids::Op::DUAL}};
inline constexpr Solids::OpStep CHAIN_AMBO_DUAL[] = {{Solids::Op::AMBO},
                                                     {Solids::Op::DUAL}};
inline constexpr Solids::OpStep CHAIN_HK62_DUAL[] = {
    {Solids::Op::HANKIN, 62.0f * Solids::IslamicStarPatterns::D2R},
    {Solids::Op::DUAL}};
inline constexpr Solids::Recipe DODECAHEDRON_KIS_RECIPE = {
    Solids::SEED_DODECAHEDRON, CHAIN_KIS, std::size(CHAIN_KIS)};
inline constexpr Solids::Recipe CUBE_KIS_DUAL_RECIPE = {
    ConwayGraph::CUBE, CHAIN_KIS_DUAL, std::size(CHAIN_KIS_DUAL)};
inline constexpr Solids::Recipe ICOSAHEDRON_AMBO_DUAL_RECIPE = {
    Solids::SEED_ICOSAHEDRON, CHAIN_AMBO_DUAL, std::size(CHAIN_AMBO_DUAL)};
inline constexpr Solids::Recipe DODECAHEDRON_HK62_DUAL_RECIPE = {
    Solids::SEED_DODECAHEDRON, CHAIN_HK62_DUAL, std::size(CHAIN_HK62_DUAL)};

/** @brief Identity measurements for one generated relax source. */
struct RelaxSourceIdentity {
  uint32_t hash;
  float quantization_margin;
};

/** @brief Measures the pre-relax source of a dodecahedron bevel recipe. */
inline RelaxSourceIdentity dodecahedron_bevel_source_identity(float depth) {
  Arena a(morph_target_buf, sizeof(morph_target_buf));
  Arena b(morph_temp_buf, sizeof(morph_temp_buf));
  PolyMesh source =
      Solids::SolidBuilder(Solids::Platonic::dodecahedron(a, b), a, b)
          .bevel(depth)
          .build();
  return {MeshOps::relax_source_hash(source),
          MeshOps::relax_source_quantization_margin(source)};
}

/** @brief Source identity separates the bakes whose topology hashes collide. */
inline void test_relax_source_hash_separates_bevel_inputs() {
  const RelaxSourceIdentity truncated =
      dodecahedron_bevel_source_identity(Solids::T_TRUNC_ICOS);
  const RelaxSourceIdentity bevel20 = dodecahedron_bevel_source_identity(0.2f);

  HS_EXPECT_EQ(
      truncated.hash,
      Solids::RelaxBakes::truncated_icosidodecahedron_converged.source_hash);
  HS_EXPECT_EQ(bevel20.hash,
               Solids::RelaxBakes::dodecahedron_bevel20_converged.source_hash);
  HS_EXPECT_TRUE(truncated.hash != bevel20.hash);
  HS_EXPECT_GE(truncated.quantization_margin, MeshOps::RELAX_SOURCE_MIN_MARGIN);
  HS_EXPECT_GE(bevel20.quantization_margin, MeshOps::RELAX_SOURCE_MIN_MARGIN);
}

inline void test_recipe_chain_build_replay() {
  replay_build_chain("dodecahedron_kis", DODECAHEDRON_KIS_RECIPE);
  replay_build_chain("cube_kis_dual", CUBE_KIS_DUAL_RECIPE);
  replay_build_chain("icosahedron_ambo_dual", ICOSAHEDRON_AMBO_DUAL_RECIPE);
  replay_build_chain("dodecahedron_hk62_dual", DODECAHEDRON_HK62_DUAL_RECIPE);
}
