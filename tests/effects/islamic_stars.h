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
  /** @brief Enables or bypasses early immutable bridge preparation. */
  static void preparation(IS &e, bool enabled) {
    e.test_prepare_build = enabled;
  }
  /** @return Opening, endpoint and bookend preparation counts. */
  static std::array<int, 3> preparations(const IS &e) {
    return {e.prepared_openings, e.prepared_duals, e.prepared_bookends};
  }
  /** @brief Selects the reference ripple capture loop and clears its counters. */
  static void reference_ripples(bool enabled) {
    IS::test_reference_ripples = enabled;
    IS::test_vertex_passes = 0;
    IS::test_ripple_candidates = 0;
    IS::test_ripple_interpolations = 0;
  }
  /**
   * @brief Vertex path evaluations of an effect at any resolution.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @return Evaluations so far, at the frame's start and at captures.
   */
  template <int W, int H> static size_t vertex_passes() {
    return IslamicStars<W, H>::test_vertex_passes;
  }
  /** @return Vertex passes, visited bands and interpolated ripples. */
  static std::array<size_t, 3> ripple_counts() {
    return {IS::test_vertex_passes, IS::test_ripple_candidates,
            IS::test_ripple_interpolations};
  }
  /** @brief Selects the ungrouped endpoint chain and resets its counters. */
  static void reference_endpoints(bool enabled) {
    IS::test_reference_endpoints = enabled;
    IS::test_endpoint_groups = 0;
    IS::test_endpoint_bands = 0;
  }
  /** @return Group and individual endpoint band checks. */
  static std::array<size_t, 2> endpoint_counts() {
    return {IS::test_endpoint_groups, IS::test_endpoint_bands};
  }
  static constexpr const auto &presets() { return IS::PRESETS; }
  /**
   * @brief Returns the active build entry.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @param e Effect to inspect.
   * @return Current shape entry.
   */
  template <int W, int H>
  static const Solids::Entry *selected_entry(const IslamicStars<W, H> &e) {
    return e.build_entry;
  }
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
    return static_cast<int>(e.getPresetIndex());
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
  /** @return One default registry cycle using the runtime build planner. */
  static size_t runtime_cycle_frames() {
    IS effect;
    size_t frames = 0;
    for (const auto &preset : presets()) {
      const auto &entry = *preset.params;
      effect.build_step_count =
          entry.recipe ? Solids::expand_to_primitives(*entry.recipe,
                                                      effect.build_step_chain,
                                                      IS::MAX_BUILD_STEPS)
                       : 0;
      bool morphable = true;
      for (size_t k = 0; k < effect.build_step_count; ++k)
        morphable &= Solids::is_morphable_step(effect.build_step_chain[k]);
      frames += IS::PRESET_DWELL_FRAMES +
                (morphable ? effect.plan_build_legs(1.0f) : 0);
    }
    return frames;
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
  template <int W, int H>
  static int active_ripples(const IslamicStars<W, H> &e) {
    return e.ripple_gen.active_count();
  }
  /** @brief One build-leg draw as the draw callback received it. */
  struct BuildDraw {
    std::vector<math::Vector> vertices; ///< The drawn mesh's vertices.
    std::vector<math::Vector> start;    ///< Its leg path at u = 0.
    std::vector<math::Vector> end;      ///< Its leg path at u = 1.
    bool on_leg_path = false;           ///< Whether the leg handed its path.
  };
  /** @brief Build-leg draws recorded by record_build_draws(). */
  struct BuildRecorder {
    IS *effect = nullptr;         ///< Effect whose draws are forwarded.
    std::vector<BuildDraw> draws; ///< Recorded draws, oldest first.
  };
  /**
   * @brief Routes an effect's build-leg draws through a recorder.
   * @param e Effect whose draw callback is replaced.
   * @param recorder Receives every draw; must outlive the effect's legs.
   */
  static void record_build_draws(IS &e, BuildRecorder &recorder) {
    recorder.effect = &e;
    e.draw_build_fn = decltype(e.draw_build_fn)(
        [r = &recorder](Canvas &c, MeshState &m,
                        const Animation::OpLeg::Shading &sh) {
          BuildDraw draw;
          draw.vertices.assign(m.vertices.begin(), m.vertices.end());
          draw.on_leg_path = sh.leg != nullptr;
          if (draw.on_leg_path)
            for (size_t i = 0; i < m.vertices.size(); ++i) {
              draw.start.push_back(sh.vertex_at(i, 0.0f));
              draw.end.push_back(sh.vertex_at(i, 1.0f));
            }
          r->draws.push_back(std::move(draw));
          r->effect->draw_build_mesh(c, m, sh);
        });
  }

  /**
   * @brief An active ripple's params, writable.
   * @param e Effect whose ripple pool is read.
   * @param k Active index.
   * @return The ripple's live params.
   */
  template <int W, int H>
  static Animation::RippleParams &ripple(IslamicStars<W, H> &e, int k) {
    return const_cast<Animation::RippleParams &>(e.ripple_gen.active_params(k));
  }
  /**
   * @brief Vertices of the front slot's mesh.
   * @param e Effect whose front mesh is read.
   * @return The mesh's model-space vertices.
   */
  template <int W, int H>
  static const ArenaVector<math::Vector> &
  front_vertices(IslamicStars<W, H> &e) {
    return e.carousel.slot(e.carousel.front_index()).vertices;
  }
  template <int W, int H>
  static float first_ripple_phase(const IslamicStars<W, H> &e) {
    return e.ripple_gen.active_params(0).phase;
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

  template <int W, int H>
  static float record_phase(IslamicStars<W, H> &e, int slot, float phase) {
    return e.record_phase(slot, phase);
  }
  template <int W, int H>
  static bool phase_drawn(const IslamicStars<W, H> &e, int slot) {
    return e.phase_history[slot].initialized;
  }
  /**
   * @brief Draws the front slot's shape at a camera held still.
   * @param e Effect whose front shape is drawn.
   * @param phase Segue phase now.
   * @param previous_phase Segue phase a frame ago.
   * @return Displayed frame.
   */
  template <int W, int H>
  static std::vector<Pixel> draw_front(IslamicStars<W, H> &e, float phase,
                                       float previous_phase) {
    const int front = e.carousel.front_index();
    {
      Canvas canvas(e);
      e.draw_shape(canvas, phase, previous_phase, e.carousel.slot(front),
                   e.slot_face_palette[front]);
    }
    e.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(e, out);
    return out;
  }

  /**
   * @brief Camera of an effect at its current endpoint.
   * @param e Effect whose orientation is read.
   * @return The orientation's latest quaternion.
   */
  template <int W, int H>
  static math::Quaternion camera(const IslamicStars<W, H> &e) {
    return e.orientation.get();
  }
  /**
   * @brief The active ripples' params.
   * @param e Effect whose ripple pool is read.
   * @return One params per active ripple, in spawn order.
   */
  template <int W, int H>
  static std::vector<Animation::RippleParams>
  ripples(const IslamicStars<W, H> &e) {
    std::vector<Animation::RippleParams> out;
    for (int k = 0; k < e.ripple_gen.active_count(); ++k)
      out.push_back(e.ripple_gen.active_params(k));
    return out;
  }
  /**
   * @brief A front-shape vertex on a frame's ripple path, model space.
   * @param e Effect whose front shape is read.
   * @param frame The frame's ripples, as ripples() gives them.
   * @param i Vertex index.
   * @param u Frame fraction.
   * @return The rippled vertex.
   */
  template <int W, int H>
  static math::Vector
  ripple_path(const IslamicStars<W, H> &e,
              const std::vector<Animation::RippleParams> &frame, size_t i,
              float u) {
    std::vector<Animation::RippleParams> during;
    for (const Animation::RippleParams &p : frame)
      during.push_back(p.during_step(u));
    return e.ripple_gen.transform(
        e.carousel.slot(e.carousel.front_index()).vertices[i], during.data());
  }
  /**
   * @brief Shades a front-shape face from its palette as draw_shape() does on
   *        the segue plateau.
   * @param e Effect whose palettes shade the face.
   * @param face Face index.
   * @param size The face's inradius.
   * @param shader Receives the face's palette and gradient scale.
   */
  template <int W, int H>
  static void shade_face(IslamicStars<W, H> &e, size_t face, float size,
                         FacePaletteShader &shader) {
    const int slot = e.carousel.front_index();
    shader.set_palette(&e.palette_bank[e.slot_face_palette[slot][face]].view());
    shader.scale = size > math::TOLERANCE ? 1.0f / size : 0.0f;
  }
  /**
   * @brief Draws the front shape's rippling frame through a moving camera,
   *        recording the faces scanned against each arm boundary.
   * @param e Effect whose front shape is drawn; its camera is replaced.
   * @param from Camera at the previous endpoint.
   * @param to Camera at the current endpoint.
   * @param copies Receives the arm sides' copies.
   */
  template <int W, int H>
  static void draw_copies(IslamicStars<W, H> &e, const math::Quaternion &from,
                          const math::Quaternion &to, FaceCopies &copies) {
    e.orientation.set(from).push(to);
    {
      Canvas canvas(e);
      ScratchScope scope(scratch_arena_a);
      const MeshState &source = e.carousel.slot(e.carousel.front_index());
      const typename IslamicStars<W, H>::ModelPath path{source.vertices.data()};
      const float bound = e.motion_bound(path, source.vertices.size());
      const auto time = Scan::prepare_column_basis<W, H>(scratch_arena_a,
                                                         e.orientation, true);
      Scan::ColumnBasisScope columns(canvas, time);
      MeshState mesh = e.transform_shape(source);
      const auto plan = e.plan_motion(canvas, mesh, path, bound,
                                      scratch_arena_a, scratch_arena_b);
      HS_EXPECT_TRUE(plan.active);
      FacePaletteShader fragment;
      int window = -1;
      Scan::deform_mesh<W, H>(
          canvas, mesh, plan.active, plan.columns, plan.moving, plan.start,
          scratch_arena_a, bound,
          [&](size_t i, float u) {
            return e.capture_position(plan, path, i, u);
          },
          [&]() {
            window = record_face_copies(canvas, mesh, e.orientation, copies);
            auto setup = [&](size_t f, float size) {
              shade_face(e, f, size, fragment);
              if (window >= 0)
                copies[window][f].size = size;
            };
            Scan::Mesh::draw_specialized<W, H>(e.filters, canvas, mesh,
                                               fragment, scratch_arena_a,
                                               nullptr, setup);
          });
    }
    e.advance_display();
  }
};

/** @brief Channel bound of an IslamicStars column against its lerped
 *  oracle, where two faces share an edge pixel. */
inline constexpr int ISLAMIC_OPACITY_BOUND = 32;

/** Whole-solid generator of the needle-ending recipe. */
inline PolyMesh generate_needle_recipe_solid(Arena &a, Arena &b) {
  return Solids::build_recipe(
      TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_RECIPE, a, b);
}

/** The needle-ending recipe as a spawnable entry. */
inline constexpr Solids::Entry NEEDLE_ENTRY = {
    "truncatedIcosahedron_ambo_relax_hk54_needle", generate_needle_recipe_solid,
    Solids::Category::COMPLEX,
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
                                 Solids::Category::COMPLEX, &recipe};
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

/** @brief Index of the first playback preset with a recipe. */
inline size_t first_recipe_entry() {
  const auto &presets = IslamicBuildProbe::presets();
  for (size_t i = 0; i < presets.size(); ++i)
    if (presets[i].params->recipe)
      return i;
  HS_CHECK(false, "no Islamic solid has a recipe");
  return 0;
}

/** @brief The shipping cycle budget includes every shape's planned build. */
inline void test_islamicstars_default_cycle_frames() {
  reset_effect_globals();
  using IS = IslamicBuildProbe::IS;
  HS_EXPECT_EQ(IS::cycle_frames(), 4896u);
  HS_EXPECT_EQ(IslamicBuildProbe::runtime_cycle_frames(), IS::cycle_frames());
}

/**
 * @brief Pins the playback roster order, the boot preset and selection
 *        mapping, and which controls are preset-scoped versus global.
 */
inline void test_islamicstars_presets_cover_roster() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  effect.init();
  constexpr std::string_view EXPECTED_NAMES[] = {
      "dodecahedron_hk62_ambo_hk62",
      "truncatedIcosahedron_hk58_chamfer63",
      "dodecahedron_bevel2_relax_gyro",
      "truncatedCube_hk90_meta",
      "truncatedIcosahedron_ambo_relax_truncate001_hankin59",
      "icosahedron_ambo_truncate033_hankin59",
      "dodecahedron_hk35_ambo_hk62_ambo_relax_hk42",
      "octahedron_hk17_ambo_hk73",
      "truncatedIcosidodecahedron_truncate50d_ambo_dual",
      "icosidodecahedron_truncate5d_ambo_dual",
      "snubDodecahedron_truncate5d_ambo_dual",
      "rhombicuboctahedron_hk63_ambo_hk63",
      "dodecahedron_hk54_ambo_hk72",
      "dodecahedron_hk72_ambo_dual_hk20",
      "icosahedron_snub_relax_truncate033_hankin62",
      "icosahedron_truncate25_hk60_ambo_relax100_hk60"};
  HS_EXPECT_EQ(effect.getPresetCount(), std::size(EXPECTED_NAMES));
  HS_EXPECT_EQ(effect.getPresetIndex(), first_recipe_entry());
  for (size_t i = 0; i < std::size(EXPECTED_NAMES); ++i) {
    HS_EXPECT_TRUE(effect.selectPreset(i));
    HS_EXPECT_EQ(
        std::string_view(IslamicBuildProbe::selected_entry(effect)->name),
        EXPECTED_NAMES[i]);
  }
  HS_EXPECT_FALSE(effect.selectPreset(std::size(EXPECTED_NAMES)));
  for (const char *name :
       {"Face Fade Lo", "Face Fade Hi", "Ripp Amp", "Ripp Decay"}) {
    const auto *def = effect.getParameters().find(name);
    HS_EXPECT_TRUE(def && !def->preset);
  }
  for (const char *name : {"Burst", "Ripp Dur", "Trans Speed"}) {
    const auto *def = effect.getParameters().find(name);
    HS_EXPECT_TRUE(def && def->preset);
  }
}

/**
 * @brief A manual select mid-build cuts the build; the selected recipe shape
 *        segues in and builds while paused.
 */
inline void test_islamicstars_manual_select_cuts_build() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  IslamicBuildProbe::set_trans_speed(effect, 8.0f);
  effect.init();
  int frames = 0;
  while (!IslamicBuildProbe::build_active(effect) && frames++ < 64) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(IslamicBuildProbe::build_active(effect));

  const size_t selected =
      (effect.getPresetIndex() + 1) % effect.getPresetCount();
  HS_EXPECT_TRUE(IslamicBuildProbe::presets()[selected].params->recipe !=
                 nullptr);
  HS_EXPECT_TRUE(effect.selectPreset(selected));
  HS_EXPECT_FALSE(IslamicBuildProbe::build_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), selected);
  HS_EXPECT_TRUE(effect.animations_paused());
  bool built = false;
  for (int frame = 0; frame < 256 && !built; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_LE(persistent_arena.get_offset(),
                 IslamicBuildProbe::persistent_budget(effect));
    built = IslamicBuildProbe::build_active(effect);
  }
  HS_EXPECT_TRUE(built);
  HS_EXPECT_EQ(effect.getPresetIndex(), selected);
  HS_EXPECT_GT((frame_energy<SMALL_W, SMALL_H>(effect)), uint64_t(0));
}

/**
 * @brief Pause holds the built preset with moving, repeating ripples;
 *        unpausing resumes advances.
 */
inline void test_islamicstars_pause_holds_shape() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  IslamicBuildProbe::set_trans_speed(effect, 8.0f);
  IslamicBuildProbe::set_burst_size(effect, 1);
  effect.init();
  const size_t held = effect.getPresetIndex();
  effect.draw_frame();
  effect.advance_display();
  effect.setAnimationsPaused(true);
  bool built = false;
  for (int frame = 0; frame < 400; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    built |= IslamicBuildProbe::build_active(effect);
  }
  HS_EXPECT_TRUE(built);
  HS_EXPECT_FALSE(IslamicBuildProbe::build_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), held);
  HS_EXPECT_GT((frame_energy<SMALL_W, SMALL_H>(effect)), uint64_t(0));

  const int window = IslamicBuildProbe::cached_burst_window(effect);
  int moving_frames = 0;
  int repeated_bursts = 0;
  float previous_phase = 0.0f;
  for (int frame = 0; frame < 3 * window + 3; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_EQ(effect.getPresetIndex(), held);
    if (IslamicBuildProbe::active_ripples(effect) > 0) {
      const float phase = IslamicBuildProbe::first_ripple_phase(effect);
      moving_frames += phase > previous_phase;
      repeated_bursts += phase < previous_phase;
      previous_phase = phase;
    }
  }
  HS_EXPECT_GT(moving_frames, window);
  HS_EXPECT_GE(repeated_bursts, 2);

  effect.setAnimationsPaused(false);
  int frames = 0;
  while (effect.getPresetIndex() == held && frames++ < 400) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(effect.getPresetIndex(), (held + 1) % effect.getPresetCount());
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
  hs_test::run_frames(effect, 12);
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

    const auto &solids = IslamicBuildProbe::presets();
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
                    (cur >= 0 && cur < entries) ? solids[cur].params->name
                                                : "?",
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
                                 ? solids[worst_p_idx].params->name
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

/** @brief Early shape-8 preparation preserves every build frame and fits its arena. */
inline void test_islamicstars_prepared_build() {
  using P = IslamicBuildProbe;
  std::vector<uint64_t> reference;
  std::array<int, 3> counts{};
  for (int mode = 0; mode < 2; ++mode) {
    reset_effect_globals();
    pin_frame_clock(0);
    P::IS effect;
    P::set_trans_speed(effect, 4.0f);
    effect.init();
    P::preparation(effect, mode != 0);
    P::spawn_entry(effect, *P::presets()[8].params);
    bool saw_build = false, completed = false;
    for (int frame = 0; frame < 120; ++frame) {
      pin_frame_clock(frame);
      effect.draw_frame();
      saw_build = saw_build || P::build_active(effect);
      completed = completed || (saw_build && !P::build_active(effect));
      effect.advance_display();
      HS_EXPECT_LE(persistent_arena.get_offset(), P::persistent_budget(effect));
      const uint64_t hash = frame_hash<SMALL_W, SMALL_H>(effect);
      if (mode == 0)
        reference.push_back(hash);
      else
        HS_EXPECT_EQ(hash, reference[frame]);
    }
    if (mode == 0) {
      HS_EXPECT_TRUE(saw_build);
      HS_EXPECT_TRUE(completed);
    }
    counts = P::preparations(effect);
    for (const int count : counts)
      HS_EXPECT_EQ(count, mode);
  }

  hs::clear_mock_time();
}

/** @brief Prepared shape-6 ripples preserve frames while rejecting interpolation work. */
inline void test_islamicstars_prepared_ripples() {
  using P = IslamicBuildProbe;
  std::vector<uint64_t> reference;
  std::array<size_t, 3> counts[2];
  for (int mode = 0; mode < 2; ++mode) {
    reset_effect_globals();
    pin_frame_clock(0);
    P::IS effect;
    P::set_trans_speed(effect, 4.0f);
    effect.init();
    P::set_burst_size(effect, 8);
    P::reference_ripples(mode == 0);
    P::spawn_entry(effect, *P::presets()[6].params);
    for (int frame = 0; frame < 120; ++frame) {
      pin_frame_clock(frame);
      effect.draw_frame();
      effect.advance_display();
      const uint64_t hash = frame_hash<SMALL_W, SMALL_H>(effect);
      if (mode == 0)
        reference.push_back(hash);
      else
        HS_EXPECT_EQ(hash, reference[frame]);
    }
    counts[mode] = P::ripple_counts();
  }
  HS_EXPECT_GT(counts[0][0], size_t{0});
  HS_EXPECT_EQ(counts[0][0], counts[1][0]);
  HS_EXPECT_EQ(counts[0][1], counts[1][1]);
  HS_EXPECT_LT(counts[1][2], counts[0][2]);

  P::reference_ripples(false);
  hs::clear_mock_time();
}

/** @brief Group rejection preserves shape-6 and shape-8 frames. */
inline void test_islamicstars_endpoint_groups() {
  using P = IslamicBuildProbe;
  for (const int shape : {6, 8}) {
    std::vector<uint64_t> reference;
    std::array<size_t, 2> counts[2];
    for (int mode = 0; mode < 2; ++mode) {
      reset_effect_globals();
      pin_frame_clock(0);
      P::IS effect;
      P::set_trans_speed(effect, 4.0f);
      effect.init();
      P::set_burst_size(effect, 8);
      P::reference_endpoints(mode == 0);
      P::spawn_entry(effect, *P::presets()[shape].params);
      for (int frame = 0; frame < 120; ++frame) {
        pin_frame_clock(frame);
        effect.draw_frame();
        effect.advance_display();
        const uint64_t hash = frame_hash<SMALL_W, SMALL_H>(effect);
        if (mode == 0)
          reference.push_back(hash);
        else
          HS_EXPECT_EQ(hash, reference[frame]);
      }
      counts[mode] = P::endpoint_counts();
    }
    HS_EXPECT_GT(counts[1][0], size_t{0});
    HS_EXPECT_LT(counts[1][1], counts[0][1]);
    std::printf(
        "  [shape %d endpoint groups] individual bands %zu -> %zu, group checks %zu\n",
        shape, counts[0][1], counts[1][1], counts[1][0]);
  }
  P::reference_endpoints(false);
  hs::clear_mock_time();
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

/**
 * @brief Each column of a segue draw shades every face's opacity at its sweep
 *        time: the frame lerped between the two endpoint frames, on both
 *        arms.
 * @details Faces keep their colours, so a lone face's pixel is linear in its
 * opacity; a pixel two faces share blends their products, within the bound.
 */
inline void test_islamicstars_opacity_shades_at_column_time() {
  using P = IslamicBuildProbe;
  reset_effect_globals();
  P::IS effect;
  effect.init();
  constexpr float PREVIOUS = 0.35f, PHASE = 0.55f;
  const std::vector<Pixel> drawn = P::draw_front(effect, PHASE, PREVIOUS);
  const std::vector<Pixel> from = P::draw_front(effect, PREVIOUS, PREVIOUS);
  const std::vector<Pixel> to = P::draw_front(effect, PHASE, PHASE);
  const std::vector<Pixel> expected =
      column_time_oracle<SMALL_W, SMALL_H>([&](float u) {
        std::vector<Pixel> frame(from.size());
        for (size_t i = 0; i < frame.size(); ++i) {
          const auto lerp = [u](uint16_t a, uint16_t b) {
            return static_cast<uint16_t>(lroundf(a + (b - a) * u));
          };
          frame[i] = Pixel(lerp(from[i].r, to[i].r), lerp(from[i].g, to[i].g),
                           lerp(from[i].b, to[i].b));
        }
        return frame;
      });
  HS_EXPECT_TRUE(from != to);
  HS_EXPECT_TRUE(expect_frames_within("islamic opacity", drawn, expected,
                                      ISLAMIC_OPACITY_BOUND) > 200);
}

/**
 * @brief A slot's previous segue phase is the one it drew last, held across
 *        a repeat draw, reset by a spawn into the slot and replaced by the
 *        current phase with temporal sampling off.
 */
inline void test_islamicstars_opacity_history() {
  using P = IslamicBuildProbe;
  reset_effect_globals();
  P::IS effect;
  effect.init();
  const int front = P::front_slot(effect);
  HS_EXPECT_FALSE(P::phase_drawn(effect, front));
  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(P::phase_drawn(effect, front));

  HS_EXPECT_EQ(P::record_phase(effect, 1 - front, 0.2f), 0.2f);
  HS_EXPECT_EQ(P::record_phase(effect, 1 - front, 0.3f), 0.2f);
  HS_EXPECT_EQ(P::record_phase(effect, 1 - front, 0.3f), 0.3f);
  Spacetime::set_flat_surface(true);
  HS_EXPECT_EQ(P::record_phase(effect, 1 - front, 0.5f), 0.5f);
  Spacetime::set_flat_surface(false);
  HS_EXPECT_EQ(P::record_phase(effect, 1 - front, 0.6f), 0.5f);

  const size_t selected =
      (effect.getPresetIndex() + 1) % effect.getPresetCount();
  HS_EXPECT_TRUE(effect.selectPreset(selected));
  HS_EXPECT_FALSE(P::phase_drawn(effect, P::front_slot(effect)));
}

/**
 * @brief The columns of a held shape show its ripple at their sweep time,
 *        across both arm boundaries.
 * @details One ripple steps once at the shipping speed and amplitude, its
 * wavefront at the steepest slope of the wavelet on the shape's vertices, so
 * they move about two pixels of arc over the frame.
 */
inline void test_islamicstars_ripple_motion_is_seam_free() {
  using P = IslamicBuildProbe;
  constexpr int W = 288, H = 144;
  reset_effect_globals();
  IslamicStars<W, H> effect;
  effect.init();
  HS_EXPECT_GT(P::fire_ripple(effect), 0);
  Animation::RippleParams &ripple = P::ripple(effect, 0);
  const auto &vertices = P::front_vertices(effect);
  // The wavelet is steepest 0.37 of its half-width ahead of a vertex.
  const float reach = 0.37f * ripple.half_width();
  const float distance = math::fast_acos(hs::clamp(
      math::dot(vertices[0].normalized(), ripple.center), -1.0f, 1.0f));
  ripple.amplitude = RIPPLE_SMALL_ANGLE_MAX;
  ripple.previous_amplitude = RIPPLE_SMALL_ANGLE_MAX;
  ripple.phase = distance - reach;
  ripple.previous_phase = ripple.phase - math::PI_F / 80.0f;
  ripple.sync();
  const std::vector<Pixel> drawn = P::draw_front(effect, 1.0f, 1.0f);
  const Animation::RippleParams moving = ripple;
  const auto at = [&](float u) {
    ripple = moving.during_step(u);
    ripple.previous_phase = ripple.phase;
    ripple.previous_amplitude = ripple.amplitude;
    return P::draw_front(effect, 1.0f, 1.0f);
  };
  const std::vector<Pixel> held = at(1.0f);
  expect_seam_resolved<W, H>("islamic ripple", drawn, held, at);
}

/**
 * @brief An IslamicStars face across the arm boundary scans as one polygon,
 *        shaded alike, in one ripple frame's late copy and the next frame's
 *        early copy.
 * @details The roster's truncatedCube_hk90_meta, generated whole, under a
 * ripple centred on the arm boundary's pole so that its wavefront runs along
 * the boundary; the ripple steps as a decaying ripple does while the camera
 * turns.
 */
inline void test_islamicstars_seam_copies_are_one_polygon() {
  using P = IslamicBuildProbe;
  constexpr int W = 288, H = 144;
  reset_effect_globals();
  IslamicStars<W, H> effect;
  effect.init();
  Solids::Entry entry{};
  for (const auto &preset : P::presets())
    if (std::string_view(preset.params->name) == "truncatedCube_hk90_meta")
      entry = *preset.params;
  HS_EXPECT_TRUE(entry.name != nullptr);
  entry.recipe = nullptr;
  P::spawn_entry(effect, entry);
  HS_EXPECT_GT(P::fire_ripple(effect), 0);
  Animation::RippleParams &ripple = P::ripple(effect, 0);
  const math::Quaternion step =
      math::make_rotation(math::Vector(0.3f, 1.0f, 0.2f).normalized(), 0.01f);
  std::array<math::Quaternion, 3> cameras{P::camera(effect)};
  for (int j = 1; j < 3; ++j)
    cameras[j] = (step * cameras[j - 1]).normalized();
  ripple.center =
      math::rotate(math::Vector(0.0f, 0.0f, 1.0f), cameras[1].conjugate());
  constexpr float STEP = math::PI_F / 40.0f;
  const float steepest = math::PI_F / 2.0f - 0.37f * ripple.half_width();
  const float phases[3] = {steepest - STEP, steepest, steepest + STEP};
  const float amplitudes[3] = {RIPPLE_SMALL_ANGLE_MAX,
                               0.97f * RIPPLE_SMALL_ANGLE_MAX,
                               0.94f * RIPPLE_SMALL_ANGLE_MAX};
  std::array<FaceCopies, 2> frames;
  std::array<std::vector<Animation::RippleParams>, 2> ripples;
  for (int j = 0; j < 2; ++j) {
    ripple.previous_phase = phases[j];
    ripple.phase = phases[j + 1];
    ripple.previous_amplitude = amplitudes[j];
    ripple.amplitude = amplitudes[j + 1];
    ripple.sync();
    ripples[j] = P::ripples(effect);
    P::draw_copies(effect, cameras[j], cameras[j + 1], frames[j]);
  }
  const size_t count = P::front_vertices(effect).size();
  float arc = 0.0f;
  for (int j = 0; j < 2; ++j)
    for (size_t i = 0; i < count; ++i)
      arc = fmaxf(
          arc,
          acosf(fminf(1.0f,
                      math::dot(P::ripple_path(effect, ripples[j], i, 0.0f),
                                P::ripple_path(effect, ripples[j], i, 1.0f)))));
  HS_EXPECT_TRUE(Scan::mesh_motion_active<W>(arc));
  const float turn = quaternion_angle(cameras[0], cameras[1]);
  const size_t faces = expect_seams_match<W, H>(
      effect, "islamic seam copies", frames, cameras[1],
      [&](uint16_t i) {
        return seam_vertex_tolerance<W>(
            [&](float u) { return P::ripple_path(effect, ripples[0], i, u); },
            [&](float u) { return P::ripple_path(effect, ripples[1], i, u); },
            cameras[1], turn);
      },
      [&](size_t f, float size, FacePaletteShader &shader) {
        P::shade_face(effect, f, size, shader);
      });
  HS_EXPECT_GE(faces, size_t{4});
}

/**
 * @brief The playback's moving mesh shows no excess colour step at the arm
 *        boundaries or mid-arm.
 * @details Its faces are small beside a half turn, so a seam step stays
 * within a few percent of its neighbours'. Faces reaching a pole draw once
 * per half-arm window, as on master, whose steps bound these: 1.025 and
 * 1.062.
 */
inline void test_islamicstars_playback_motion_is_seam_free() {
  constexpr int W = 288, H = 144;
  reset_effect_globals();
  pin_frame_clock(0);
  IslamicStars<W, H> effect;
  effect.init();
  std::vector<Pixel> previous, current;
  MeshSeamSteps steps;
  for (int f = 0; f < 240; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    capture_frame<W, H>(effect, current);
    if (f > 0)
      steps.add(previous, current, W, H);
    previous.swap(current);
  }
  expect_mesh_seam_continuous(steps, "islamic playback", 1.025, 1.062);
  hs::clear_mock_time();
}

/**
 * @brief A hankin, relax or medial build leg hands each draw its own vertex
 *        path: the frame drawn last at u = 0, the drawn mesh at u = 1.
 */
inline void test_islamicstars_build_legs_hand_their_path() {
  using P = IslamicBuildProbe;
  reset_effect_globals();
  P::IS effect;
  effect.init();
  P::BuildRecorder recorder;
  P::record_build_draws(effect, recorder);
  size_t on_path = 0;
  for (int f = 0; f < 600 && on_path < 12; ++f) {
    effect.draw_frame();
    effect.advance_display();
    on_path = 0;
    for (const P::BuildDraw &draw : recorder.draws)
      on_path += draw.on_leg_path;
  }
  HS_EXPECT_GE(on_path, size_t{12});
  bool moved = false;
  for (size_t d = 1; d < recorder.draws.size(); ++d) {
    const P::BuildDraw &draw = recorder.draws[d];
    if (!draw.on_leg_path)
      continue;
    HS_CONTEXT("draw", static_cast<long long>(d));
    HS_EXPECT_TRUE(draw.end == draw.vertices);
    HS_EXPECT_TRUE(draw.start == recorder.draws[d - 1].vertices);
    moved |= draw.start != draw.end;
  }
  HS_EXPECT_TRUE(moved);
}
