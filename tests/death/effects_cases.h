/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Effects death fixtures and guard cases.

/** @brief Death case: a GS stencil must fit the delayed-write history. */
inline void case_gs_neighbor_exceeds_history() {
  ReactionGraph::NeighborRun run{};
  run.delta[0] = opaque<int16_t>(-145);
  hs_test::effects_tests::GSWhiteBox::validate_physics_neighbors(&run, 1);
}

/** @brief Death case: GS rejects a zero color-noise scale. */
inline void case_gs_color_noise_zero_scale() {
  using WB = hs_test::effects_tests::GSWhiteBox;
  WB::GS gs;
  gs.init();
  WB::set_color_params(gs, 0.0f, opaque(0.0f), 0.2f, 0.4f);
  WB::advance_color_noise(gs);
}

/** @brief Death case: GS rejects a non-finite color-noise scale. */
inline void case_gs_color_noise_nan_scale() {
  using WB = hs_test::effects_tests::GSWhiteBox;
  WB::GS gs;
  gs.init();
  WB::set_color_params(
      gs, 0.0f, opaque(std::numeric_limits<float>::quiet_NaN()), 0.2f, 0.4f);
  WB::advance_color_noise(gs);
}

inline void case_lattice_shells_oob() {
  SDF::Lattice::Settings settings;
  settings.shells = static_cast<SDF::Lattice::ShellCount>(3);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: zero lattice softness would divide by zero in shading. */
inline void case_lattice_zero_softness() {
  SDF::Lattice::Settings settings;
  settings.softness = opaque(0.0f);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: zero lattice cell size has no inverse transform. */
inline void case_lattice_zero_cell_size() {
  SDF::Lattice::Settings settings;
  settings.cell_size = opaque(0.0f);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: a negative AA strength gives invalid crossing widths. */
inline void case_lattice_negative_aa() {
  SDF::Lattice::Settings settings;
  settings.aa_strength = opaque(-1.0f);
  (void)SDF::Lattice::prepare(settings, {}, math::Mat4::identity(), 4.0f,
                              0.01f);
}

/** @brief Death case: a HyperLattice frame without crossing scratch traps. */
inline void case_hyperlattice_frame_without_crossings() {
  const HyperLatticeDetail::FrameState frame{};
  (void)HyperLatticeDetail::prepare_trace(frame);
}

inline void case_mindsplatter_profile_preset_oob() {
  MindSplatter<96, 20> effect;
  effect.profile_select_preset(opaque<size_t>(SIZE_MAX));
}

/** @brief Death case: contour preparation past its table capacity traps. */
inline void case_shapeshifter_count_over_capacity() {
  using namespace shapeshifter_oracle_tests;
  OracleEffect effect;
  ShapeShifterWhiteBox::prepare_count(effect,
                                      opaque(OracleEffect::MAX_SHAPES + 1));
}

/** @brief Death case: a woven edge whose start vertex is absent must trap. */
inline void case_dreamballs_woven_owner_vertex_oob() {
  using WB = effects_tests::DreamBallsWhiteBox;
  static uint8_t arena_buf[64];
  Arena arena(arena_buf, sizeof(arena_buf));
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(arena, 1);
  edges.push_back({opaque<uint16_t>(1), 0});
  uint16_t owners[1];
  WB::assign_woven_start_owners(edges, owners, 1);
}

/** @brief Death case: a woven-edge ownership query past the list must trap. */
inline void case_dreamballs_woven_owner_edge_oob() {
  using WB = effects_tests::DreamBallsWhiteBox;
  static uint8_t arena_buf[64];
  Arena arena(arena_buf, sizeof(arena_buf));
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(arena, 1);
  edges.push_back({0, 0});
  const std::vector<uint16_t> owners{0};
  if (WB::owns_woven_start(edges, owners, opaque<size_t>(1)))
    std::printf("x");
}

/** @brief Death case: a Raymarch placement-solid index past its table traps. */
inline void case_raymarch_placement_solid_oob() {
  using WB = effects_tests::RaymarchWhiteBox;
  effects_tests::reset_effect_globals();
  Raymarch<effects_tests::SMALL_W, effects_tests::SMALL_H> effect;
  WB::set_base_solid(effect, RaymarchPlacementSolid::COUNT);
  WB::build_points(effect);
}

/** @brief Death case: a harmonic morph cannot synchronize an invalid mode. */
inline void case_spherical_harmonics_invalid_morph_mode() {
  using WB = effects_tests::SphericalHarmonicsWhiteBox;
  effects_tests::reset_effect_globals();
  WB::SH effect;
  effect.init();
  WB::set_next_idx(effect, WB::max_mode_idx() + 1);
  Canvas canvas(effect);
  for (int frame = 0; frame < 64; ++frame)
    WB::step_timeline(effect, canvas);
}

inline void case_hankinsolids_missing_topology() {
  using WB = effects_tests::HankinPauseWhiteBox;
  effects_tests::reset_effect_globals();
  WB::EffectT effect;
  effect.init();
  WB::draw_without_topology(effect);
}

inline void case_islamicstars_build_budget() {
  configure_arenas_default();
  persistent_arena.allocate_n<uint8_t>(1);
  effects_tests::IslamicBuildProbe::IS effect;
  effects_tests::IslamicBuildProbe::check_build_budget(effect, 0);
}

inline void case_islamicstars_bridge_continuation() {
  effects_tests::IslamicBuildProbe::IS effect;
  effects_tests::IslamicBuildProbe::invalid_bridge_continuation(effect);
}

/** @brief Death case: a Hankin step has no eagerly generated endpoint. */
inline void case_islamicstars_hankin_eager_endpoint() {
  using WB = effects_tests::IslamicBuildProbe;
  effects_tests::reset_effect_globals();
  WB::IS effect;
  static uint8_t a_buf[64];
  static uint8_t b_buf[64];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep step{Solids::Op::HANKIN};
  PolyMesh out = WB::clean_endpoint(effect, step, a, b);
  if (out.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/** @brief Rejects a unsupported HyperLattice pattern ID. */
inline void case_hyperlattice_pattern_defaults_invalid() {
  using Effect = HyperLattice<32, 16>;
  for (const auto &config : Effect::CONFIGURATIONS) {
    const auto valid = Effect::pattern_defaults(config.pattern, config.domain);
    HS_EXPECT_EQ(valid.pattern, config.pattern);
    HS_EXPECT_EQ(valid.mode, config.domain);
  }
  Effect::pattern_defaults(static_cast<Effect::Pattern>(opaque(uint8_t{2})),
                           Effect::LatticeMode::THREE_D);
}
