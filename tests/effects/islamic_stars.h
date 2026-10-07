/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// IslamicStars build chain, recipes and arena budgets.
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for IslamicStars' private build-chain state.
 * @details Build bookkeeping is resolution-independent.
 */
struct IslamicBuildProbe {
  using IS = IslamicStars<SMALL_W, SMALL_H>;
  static_assert(IS::MACRO_TRUNCATE_T == RECONCILE_TRUNCATE_T);
  static void check_build_budget(IS &e, size_t budget) {
    e.device_persistent_budget = budget;
    e.check_build_budget();
  }
  static void invalid_bridge_continuation(IS &e) {
    e.schedule_dual_bridge(IS::BuildContinuation::DUAL_DONE);
  }

  template <int W, int H>
  static void set_trans_speed(IslamicStars<W, H> &e, float v) {
    e.params.trans_speed = v;
  }
  template <int W, int H>
  static bool build_active(const IslamicStars<W, H> &e) {
    return e.build_active;
  }
  template <int W, int H> static int solid_idx(const IslamicStars<W, H> &e) {
    return e.solid_idx;
  }
  template <int W, int H> static int dual_bridges(const IslamicStars<W, H> &e) {
    return e.dual_bridges_built;
  }
  template <int W, int H>
  static size_t persistent_budget(const IslamicStars<W, H> &e) {
    return e.device_persistent_budget;
  }
  template <int W, int H> static int front_slot(IslamicStars<W, H> &e) {
    return e.carousel.front_index();
  }
  template <int W, int H>
  static const uint8_t *slot_palette(const IslamicStars<W, H> &e, int slot) {
    return e.slot_face_palette[slot];
  }
  template <int W, int H>
  static size_t slot_faces(IslamicStars<W, H> &e, int slot) {
    return e.carousel.slot(slot).topology.size();
  }
  static constexpr size_t bridge_scratch_a() {
    return IS::BRIDGE_BUDGET.scratch_a;
  }
  static constexpr size_t bridge_scratch_b() {
    return IS::BRIDGE_BUDGET.scratch_b;
  }
  static constexpr int sprite_fade_frames() { return IS::SPRITE_FADE_FRAMES; }
  /**
   * @brief Spawns @p entry out of band, whatever the registry cycle is on.
   * @details Clears the timeline first: the pending scheduled spawn and the
   * outgoing shape's sprite would otherwise run against the injected shape.
   */
  template <int W, int H>
  static void spawn_entry(IslamicStars<W, H> &e, const Solids::Entry &entry) {
    e.timeline.clear();
    e.spawn_entry(entry);
  }
  template <int W, int H>
  static void set_burst_size(IslamicStars<W, H> &e, int size) {
    e.params.burst_size = size;
  }
  template <int W, int H>
  static int cached_burst_size(const IslamicStars<W, H> &e) {
    return e.burst_size_eff;
  }
  template <int W, int H>
  static int cached_burst_window(const IslamicStars<W, H> &e) {
    return (e.burst_size_eff - 1) * e.ripple_stagger_eff + e.ripple_dur_eff;
  }
  template <int W, int H> static int fire_ripple(IslamicStars<W, H> &e) {
    int count = 0;
    {
      Canvas canvas(e);
      e.ripple(canvas);
      count = e.ripple_gen.active_count();
    }
    e.advance_display();
    return count;
  }
  template <int W, int H>
  static PolyMesh clean_endpoint(IslamicStars<W, H> &e,
                                 const Solids::OpStep &step, Arena &a,
                                 Arena &b) {
    return e.clean_endpoint(step, a, b);
  }
};

/** Whole-solid generator of the needle-ending recipe. */
inline PolyMesh generate_needle_recipe_solid(Arena &a, Arena &b) {
  return Solids::build_recipe(
      TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_RECIPE, a, b);
}

/** The needle-ending recipe as a spawnable entry. */
inline constexpr Solids::Entry NEEDLE_ENTRY = {
    "truncatedIcosahedron_ambo_relax_hk54_needle", generate_needle_recipe_solid,
    Solids::Category::Complex,
    &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_RECIPE};

/**
 * @brief Verifies the first recipe seed uses the Sprite envelope's 16-frame
 *        fade-in and starts its build at the full-opacity boundary.
 */
inline void test_islamicstars_seed_sprite_fade_in() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  effect.init();

  constexpr int FADE_FRAMES = IslamicBuildProbe::sprite_fade_frames();
  HS_EXPECT_EQ(FADE_FRAMES, 16);
  for (int frame = 1; frame < FADE_FRAMES; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_FALSE(IslamicBuildProbe::build_active(effect));
  }

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(IslamicBuildProbe::build_active(effect));
}

inline void test_islamicstars_burst_size_is_snapshotted_per_spawn() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  effect.init();
  const auto entry = Solids::Collections::get_islamic_solids().front();

  IslamicBuildProbe::set_burst_size(effect, 1);
  IslamicBuildProbe::spawn_entry(effect, entry);
  const int single_window = IslamicBuildProbe::cached_burst_window(effect);

  IslamicBuildProbe::set_burst_size(effect, 4);
  HS_EXPECT_EQ(IslamicBuildProbe::cached_burst_size(effect), 1);
  HS_EXPECT_EQ(IslamicBuildProbe::cached_burst_window(effect), single_window);
  HS_EXPECT_EQ(IslamicBuildProbe::fire_ripple(effect), 1);

  IslamicBuildProbe::spawn_entry(effect, entry);
  HS_EXPECT_EQ(IslamicBuildProbe::cached_burst_size(effect), 4);
  HS_EXPECT_TRUE(IslamicBuildProbe::cached_burst_window(effect) >
                 single_window);
  HS_EXPECT_EQ(IslamicBuildProbe::fire_ripple(effect), 4);
}

/** @brief Small recipe builds finish every bridge and retain landed palettes. */
inline void test_islamicstars_smooth_recipe_completion() {
  struct Case {
    Solids::Op op;
    size_t faces;
    int bridges;
  };
  for (const Case &c :
       {Case{Solids::Op::DUAL, 8, 1}, Case{Solids::Op::NEEDLE, 24, 1},
        Case{Solids::Op::KIS, 24, 2}}) {
    reset_effect_globals();
    const Solids::OpStep steps[] = {{c.op}};
    const Solids::Recipe recipe = Solids::make_recipe(1, steps);
    HS_EXPECT_TRUE(
        std::string_view(Solids::simple_registry[recipe.seed].name) == "cube");
    const Solids::Entry entry = {"cube_bridge", Solids::Platonic::cube,
                                 Solids::Category::Complex, &recipe};
    IslamicBuildProbe::IS effect;
    IslamicBuildProbe::set_trans_speed(effect, 2.0f);
    effect.init();
    IslamicBuildProbe::spawn_entry(effect, entry);
    bool started = false;
    bool finished = false;
    for (int frame = 0; frame < 256; ++frame) {
      effect.draw_frame();
      effect.advance_display();
      HS_EXPECT_LE(persistent_arena.get_offset(),
                   IslamicBuildProbe::persistent_budget(effect));
      if (IslamicBuildProbe::build_active(effect))
        started = true;
      else if (started) {
        finished = true;
        break;
      }
    }
    HS_EXPECT_TRUE(started);
    HS_EXPECT_TRUE(finished);
    HS_EXPECT_EQ(IslamicBuildProbe::dual_bridges(effect), c.bridges);
    const int slot = IslamicBuildProbe::front_slot(effect);
    HS_EXPECT_EQ(IslamicBuildProbe::slot_faces(effect, slot), c.faces);
    const uint8_t *palette = IslamicBuildProbe::slot_palette(effect, slot);
    std::vector<uint8_t> landed(palette, palette + c.faces);
    for (uint8_t value : landed)
      HS_EXPECT_LT(value, MeshPaletteBank::N);
    for (int frame = 0; frame < 4; ++frame) {
      effect.draw_frame();
      effect.advance_display();
      HS_EXPECT_FALSE(IslamicBuildProbe::build_active(effect));
      HS_EXPECT_EQ(IslamicBuildProbe::front_slot(effect), slot);
      HS_EXPECT_TRUE(std::equal(landed.begin(), landed.end(),
                                IslamicBuildProbe::slot_palette(effect, slot)));
    }
    HS_EXPECT_GT((frame_energy<SMALL_W, SMALL_H>(effect)), uint64_t(0));
  }
}

/**
 * @brief Drives IslamicStars through the first registry entry's full build at
 *        max trans speed: no trap, stable per-face colours through its
 *        display, and lit frames once entry 2 starts.
 */
inline void test_islamicstars_recipe_build_smoke() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  IslamicBuildProbe::set_trans_speed(effect, 8.0f);
  effect.init();

  // The snapshot captures entry 0's first completed recipe build.
  constexpr int MAX_FRAMES = 400;
  int frames = 0;
  int build_frames = 0;
  bool was_building = false;
  int built_slot = -1;
  int snap_solid = -1;
  int changed_after_build = 0;
  int constant_frames = 0;
  std::vector<uint8_t> built_pal;
  while (frames < MAX_FRAMES && IslamicBuildProbe::solid_idx(effect) < 2) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
    const bool building = IslamicBuildProbe::build_active(effect);
    if (building)
      ++build_frames;
    // Snapshot the built shape's per-face colours the moment its build
    // completes; they must stay byte-identical through its still, ripple, and
    // fade phases (the next spawn retires the shape and may reuse its array).
    if (was_building && !building && built_slot < 0) {
      built_slot = IslamicBuildProbe::front_slot(effect);
      snap_solid = IslamicBuildProbe::solid_idx(effect);
      const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, built_slot);
      built_pal.assign(pal,
                       pal + IslamicBuildProbe::slot_faces(effect, built_slot));
    } else if (built_slot >= 0 &&
               IslamicBuildProbe::solid_idx(effect) == snap_solid) {
      const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, built_slot);
      for (size_t f = 0; f < built_pal.size(); ++f)
        if (pal[f] != built_pal[f])
          ++changed_after_build;
      ++constant_frames;
    }
    was_building = building;
  }
  HS_EXPECT_LT(frames, MAX_FRAMES);
  HS_EXPECT_GT(build_frames, 0);
  HS_EXPECT_TRUE(!IslamicBuildProbe::build_active(effect));
  HS_EXPECT_GE(built_slot, 0);
  HS_EXPECT_GT(built_pal.size(), (size_t)0);
  HS_EXPECT_GT(constant_frames, 0);
  HS_EXPECT_EQ(changed_after_build, 0);

  // Entry 2 has started; its frames may overlap entry 1's fade.
  for (int f = 0; f < 12; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  const uint64_t acc = frame_energy<SMALL_W, SMALL_H>(effect);
  HS_EXPECT_GT(acc, (uint64_t)0);
}

/**
 * @brief Drives IslamicStars through every registry entry and then through the
 *        needle recipe, pinning the persistent arena against the effect's own
 *        budget.
 * @details An arena overrun traps. The needle sets the scratch_a-heavy split,
 *          so it is measured separately.
 */
inline void test_islamicstars_roster_cycle_fits_budget() {
  reset_effect_globals();
  // Host scratch uses the device caps; persistent usage is checked against the
  // live per-shape device budget. Resplitting rebases scratch high-water marks.
  {
    IslamicStars<288, 144> effect;
    IslamicBuildProbe::set_trans_speed(effect, 8.0f);
    effect.init();

    auto solids = Solids::Collections::get_islamic_solids();
    const int entries = static_cast<int>(solids.size());
    constexpr int MAX_FRAMES = 20000;
    size_t a_peak = 0, b_peak = 0, persist_peak = 0;
    size_t worst_p = 0, worst_p_budget = 1;
    int worst_p_idx = -1; // shape at the worst persistent/budget ratio
    int frames = 0, shapes = 0, builds = 0;
    bool was_building = false;
    int last = IslamicBuildProbe::solid_idx(effect);
    while (frames < MAX_FRAMES && shapes <= entries) {
      effect.draw_frame();
      effect.advance_display();
      ++frames;
      const int cur = IslamicBuildProbe::solid_idx(effect);
      // Palette variety of each finished build: distinct landed palettes on the
      // shape the real leg chain just landed.
      const bool building = IslamicBuildProbe::build_active(effect);
      if (was_building && !building) {
        const int front = IslamicBuildProbe::front_slot(effect);
        const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, front);
        const size_t nf = IslamicBuildProbe::slot_faces(effect, front);
        bool seen[MeshPaletteBank::N] = {};
        int distinct = 0;
        for (size_t f = 0; f < nf; ++f)
          if (pal[f] < MeshPaletteBank::N && !seen[pal[f]]) {
            seen[pal[f]] = true;
            ++distinct;
          }
        std::printf("  [built] %s: %d/%d palettes on %zu faces\n",
                    (cur >= 0 && cur < entries) ? solids[cur].name : "?",
                    distinct, MeshPaletteBank::N, nf);
      }
      was_building = building;
      const size_t p = persistent_arena.get_offset();
      const size_t p_budget = IslamicBuildProbe::persistent_budget(effect);
      a_peak = std::max(a_peak, scratch_arena_a.get_high_water_mark());
      b_peak = std::max(b_peak, scratch_arena_b.get_high_water_mark());
      persist_peak = std::max(persist_peak, p);
      if (p_budget &&
          uint64_t(p) * worst_p_budget > uint64_t(worst_p) * p_budget) {
        worst_p = p;
        worst_p_budget = p_budget;
        worst_p_idx = cur;
      }
      // Per-shape persistent budget (device figure); scratch is trap-enforced.
      HS_EXPECT_LE(p, p_budget);
      if (building)
        ++builds;
      if (cur != last) {
        last = cur;
        ++shapes;
      }
    }

    const char *worst_name = (worst_p_idx >= 0 && worst_p_idx < entries)
                                 ? solids[worst_p_idx].name
                                 : "?";
    std::printf(
        "  [roster] %d shapes over %d frames, %d build frames: scratch_a "
        "peak=%zu B, scratch_b peak=%zu B, persistent peak=%zu B; tightest "
        "persistent %zu/%zu at %s\n",
        shapes, frames, builds, a_peak, b_peak, persist_peak, worst_p,
        worst_p_budget, worst_name);
    HS_EXPECT_GT(shapes, entries - 1);
    HS_EXPECT_GT(builds, 0);
  }

  // The needle build. Trans Speed 2, not the roster's 8: a compressed stage can
  // drop the closing bridge leg before it runs, which is the peak.
  reset_effect_globals();
  size_t na_peak = 0, nb_peak = 0, np_peak = 0;
  int needle_frames = 0;
  bool needle_built = false;
  {
    constexpr int NEEDLE_MAX_FRAMES = 4000;
    IslamicStars<288, 144> effect;
    IslamicBuildProbe::set_trans_speed(effect, 2.0f);
    effect.init();
    IslamicBuildProbe::spawn_entry(effect, NEEDLE_ENTRY);
    bool was_building = false;
    while (needle_frames < NEEDLE_MAX_FRAMES) {
      effect.draw_frame();
      effect.advance_display();
      ++needle_frames;
      const bool building = IslamicBuildProbe::build_active(effect);
      const size_t p = persistent_arena.get_offset();
      na_peak = std::max(na_peak, scratch_arena_a.get_high_water_mark());
      nb_peak = std::max(nb_peak, scratch_arena_b.get_high_water_mark());
      np_peak = std::max(np_peak, p);
      HS_EXPECT_LE(p, IslamicBuildProbe::persistent_budget(effect));
      if (building)
        was_building = true;
      else if (was_building) {
        needle_built = true;
        break;
      }
    }
  }
  std::printf("  [needle] smooth-path peaks over %d frames: scratch_a=%zu/%zu "
              "scratch_b=%zu/%zu persistent=%zu/%zu\n",
              needle_frames, na_peak, IslamicBuildProbe::bridge_scratch_a(),
              nb_peak, IslamicBuildProbe::bridge_scratch_b(), np_peak,
              DEVICE_GLOBAL_ARENA_SIZE - IslamicBuildProbe::bridge_scratch_a() -
                  IslamicBuildProbe::bridge_scratch_b());
  HS_EXPECT_TRUE(needle_built);
  // The needle reached its scratch_a-heavy split.
  HS_EXPECT_GT(na_peak, 120u * 1024u);
}

/**
 * @brief Drives IslamicStars until TARGET_BRIDGES dual bridges complete,
 *        pinning the scratch peaks against the effect's budget.
 * @details The bridge's leg 3 rebuilds the medial for its handoff centroids,
 *          whose scratch must not co-reside with the leg's own arrival mesh.
 */
inline void test_islamicstars_dual_bridge_fits_budget() {
  reset_effect_globals();
  IslamicStars<288, 144> effect;
  IslamicBuildProbe::set_trans_speed(effect, 2.0f);
  effect.init();

  constexpr int TARGET_BRIDGES = 5;
  constexpr int MAX_FRAMES = 40000;
  size_t a_peak = 0, b_peak = 0, persist_peak = 0;
  int frames = 0;
  // Scratch is hard-capped per-shape (spawn_shape's resplit), so a leg over its
  // split traps here; completing without a trap proves every full bridge fit.
  // Persistent is host-inflated, so it is checked per frame against the effect's
  // live per-shape device budget.
  while (frames < MAX_FRAMES &&
         IslamicBuildProbe::dual_bridges(effect) < TARGET_BRIDGES) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
    a_peak = std::max(a_peak, scratch_arena_a.get_high_water_mark());
    b_peak = std::max(b_peak, scratch_arena_b.get_high_water_mark());
    const size_t p = persistent_arena.get_offset();
    persist_peak = std::max(persist_peak, p);
    HS_EXPECT_LE(p, IslamicBuildProbe::persistent_budget(effect));
  }
  std::printf(
      "  [dual-bridge] %d bridges over %d frames: scratch_a peak=%zu B, "
      "scratch_b peak=%zu B, persistent peak=%zu B\n",
      IslamicBuildProbe::dual_bridges(effect), frames, a_peak, b_peak,
      persist_peak);
  HS_EXPECT_GE(IslamicBuildProbe::dual_bridges(effect), TARGET_BRIDGES);
}
