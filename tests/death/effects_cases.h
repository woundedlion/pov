/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Effects death fixtures and guard cases.

inline void check_shader_chain_colorize_companion() {
  using WB = hs_test::shader_chain_tests::ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(WB::colorize_index(effect), 3);
  reset_globals();
}

inline void case_shader_chain_missing_colorize() {
  using WB = hs_test::shader_chain_tests::ShaderChainWhiteBox;
  WB::FX effect;
  (void)WB::colorize_index(effect);
}

/** @brief PetalFlow traps when every ring slot is occupied. */
inline void case_petalflow_ring_pool_full() {
  using WB = hs_test::effects_tests::PetalFlowWhiteBox;
  WB::PF effect;
  effect.init();
  WB::fill_and_spawn(effect, true);
}

/** @brief Death case: a GS stencil must fit the delayed-write history. */
inline void case_gs_neighbor_exceeds_history() {
  using WB = hs_test::effects_tests::GSWhiteBox;
  ReactionGraph::NeighborRun run{};
  run.delta[0] = opaque<int16_t>(-(WB::PHYSICS_NEIGHBOR_REACH + 1));
  WB::validate_physics_neighbors(&run, 1);
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

/** @brief Traps on lattice shells oob. */
inline void case_lattice_shells_oob() {
  SDF::Lattice::Settings settings;
  settings.shells =
      static_cast<SDF::Lattice::ShellCount>(opaque(SDF::Lattice::MAX_SHELLS));
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

/** @brief Death case: cubic-wire patterns have no lattice trace geometry. */
inline void case_hyperlattice_cubic_trace_geometry() {
  HyperLatticeDetail::FrameState frame{};
  frame.params.pattern = opaque(HyperLatticeDetail::Pattern::CUBIC_WIRE);
  (void)HyperLatticeDetail::trace_settings(frame, {});
}

inline HyperLatticeDetail::FrameState hyperlattice_octet_frame(
    const hs_test::hyper_lattice_tests::HyperLatticeWhiteBox::Effect &effect) {
  using WB = hs_test::hyper_lattice_tests::HyperLatticeWhiteBox;
  auto frame = WB::frame(effect);
  frame.params =
      WB::Effect::pattern_defaults(HyperLatticeDetail::Pattern::OCTET,
                                   HyperLatticeDetail::LatticeMode::SPATIAL_3D);
  return frame;
}

inline void check_hyperlattice_trace_companion() {
  using WB = hs_test::hyper_lattice_tests::HyperLatticeWhiteBox;
  reset_globals();
  {
    WB::Effect effect;
    effect.init();
    HS_EXPECT_TRUE(
        WB::prepare_traced(effect, hyperlattice_octet_frame(effect)).valid);
  }
  reset_globals();
}

inline void case_hyperlattice_invalid_trace_settings() {
  using WB = hs_test::hyper_lattice_tests::HyperLatticeWhiteBox;
  WB::Effect effect;
  effect.init();
  auto frame = hyperlattice_octet_frame(effect);
  frame.pixel_half_angle = opaque(-1.0f);
  (void)WB::prepare_traced(effect, frame);
}

/** @brief Traps on mindsplatter profile preset oob. */
inline void case_mindsplatter_profile_preset_oob() {
  MindSplatter<96, 20> effect;
  effect.profile_select_preset(opaque<size_t>(SIZE_MAX));
}

/** @brief Traps on mindsplatter profile preset uninitialized. */
inline void case_mindsplatter_profile_preset_uninitialized() {
  MindSplatter<96, 20> effect;
  effect.profile_select_preset(opaque<size_t>(0));
}

/** @brief Death case: contour preparation past its table capacity traps. */
inline void case_shapeshifter_count_over_capacity() {
  using namespace shapeshifter_oracle_tests;
  OracleEffect effect;
  ShapeShifterWhiteBox::prepare_count(effect,
                                      opaque(OracleEffect::MAX_SHAPES + 1));
}

/** @brief Death case: profiling rejects an out-of-range displacement preset. */
inline void case_displacement_field_profile_preset_oob() {
  DisplacementField<96, 20> effect;
  effect.profile_select_preset(opaque<size_t>(SIZE_MAX));
}

/** @brief Death case: profiling requires initialized displacement presets. */
inline void case_displacement_field_profile_preset_uninitialized() {
  DisplacementField<96, 20> effect;
  effect.profile_select_preset(opaque<size_t>(0));
}

/** @brief Death case: Balls departure must commit its Noise preset. */
inline void case_displacement_field_noise_advance_refused() {
  DisplacementField<96, 20> effect;
  effect.init();
  effects_tests::DisplacementFieldWhiteBox::deny_phase_advance(effect, false);
  effect.draw_frame();
}

/** @brief Death case: Noise departure must commit its Balls preset. */
inline void case_displacement_field_balls_advance_refused() {
  DisplacementField<96, 20> effect;
  effect.init();
  effects_tests::DisplacementFieldWhiteBox::deny_phase_advance(effect, true);
  effect.draw_frame();
}

/** @brief Death case: preview warming requires the active frame canvas. */
inline void case_displacement_field_preview_missing_canvas() {
  DisplacementField<96, 20> effect;
  effect.init();
  effects_tests::DisplacementFieldWhiteBox::spawn_preview_without_canvas(
      effect);
}

/** @brief Commits @p n flat rings into a DisplacementField ring pool. */
template <typename Pool> inline void commit_flat_rings(Pool &pool, int n) {
  const math::Basis basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  for (int i = 0; i < n; ++i)
    pool.commit(i, 1.0f, 16, basis, 0.5f, 0.03f,
                ScalarFn([](float) { return 0.0f; }), 0.0f, 0.0f);
}

/** @brief Death case: a ring pool frame cannot start over live rings. */
inline void case_displacement_field_ring_pool_live_begin() {
  using WB = effects_tests::DisplacementFieldWhiteBox;
  DisplacementField<96, 20> effect;
  effect.init();
  auto &pool = WB::ring_pool(effect);
  pool.begin_frame();
  commit_flat_rings(pool, opaque(1));
  pool.begin_frame();
}

/** @brief Death case: a full ring pool has no next slot. */
inline void case_displacement_field_ring_pool_full() {
  using WB = effects_tests::DisplacementFieldWhiteBox;
  DisplacementField<96, 20> effect;
  effect.init();
  auto &pool = WB::ring_pool(effect);
  pool.begin_frame();
  commit_flat_rings(pool, opaque(WB::ring_slots(effect)));
  static_cast<void>(pool.next());
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
  (void)WB::owns_woven_start(edges, owners, opaque<size_t>(1));
}

/** @brief Death case: a Raymarch placement-solid index past its table traps. */
inline void case_raymarch_placement_solid_oob() {
  using WB = effects_tests::RaymarchWhiteBox;
  hs_test::reset_globals();
  Raymarch<effects_tests::SMALL_W, effects_tests::SMALL_H> effect;
  WB::set_base_solid(effect, RaymarchPlacementSolid::COUNT);
  WB::build_points(effect);
}

/** @brief Death case: a harmonic morph cannot synchronize an invalid mode. */
inline void case_spherical_harmonics_invalid_morph_mode() {
  using WB = effects_tests::SphericalHarmonicsWhiteBox;
  hs_test::reset_globals();
  WB::SH effect;
  effect.init();
  WB::set_next_idx(effect, WB::max_mode_idx() + 1);
  Canvas canvas(effect);
  for (int frame = 0; frame < 64; ++frame)
    WB::step_timeline(effect, canvas);
}

/** @brief Traps on hankinsolids missing topology. */
inline void case_hankinsolids_missing_topology() {
  using WB = effects_tests::HankinSolidsWhiteBox;
  hs_test::reset_globals();
  WB::EffectT effect;
  effect.init();
  WB::draw_without_topology(effect);
}

/** @brief Traps on islamicstars build budget. */
inline void case_islamicstars_build_budget() {
  configure_arenas_default();
  persistent_arena.allocate_n<uint8_t>(1);
  effects_tests::IslamicBuildProbe::IS effect;
  effects_tests::IslamicBuildProbe::check_build_budget(effect, 0);
}

/** @brief Traps on islamicstars bridge continuation. */
inline void case_islamicstars_bridge_continuation() {
  effects_tests::IslamicBuildProbe::IS effect;
  effects_tests::IslamicBuildProbe::invalid_bridge_continuation(effect);
}

/** @brief Death case: a Hankin step has no eagerly generated endpoint. */
inline void case_islamicstars_hankin_eager_endpoint() {
  using WB = effects_tests::IslamicBuildProbe;
  hs_test::reset_globals();
  WB::IS effect;
  static uint8_t a_buf[64];
  static uint8_t b_buf[64];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep step{Solids::Op::HANKIN};
  PolyMesh out = WB::clean_endpoint(effect, step, a, b);
  (void)out.vertices.size();
}

/** @brief Death case: an unsupported HyperLattice pattern ID must trap. */
inline void case_hyperlattice_pattern_defaults_invalid() {
  using Effect = HyperLattice<32, 16>;
  Effect::pattern_defaults(static_cast<Effect::Pattern>(opaque(uint8_t{2})),
                           Effect::LatticeMode::SPATIAL_3D);
}

struct UninitializedReactionLattice
    : ReactionDiffusionBase<UninitializedReactionLattice, 32, 16> {
  void init() override {}
  void render(Canvas &) {}
  using ReactionDiffusionBase::orient_lattice;
};

/** @brief Traps on reaction lattice uninitialized. */
inline void case_reaction_lattice_uninitialized() {
  UninitializedReactionLattice effect;
  auto lattice = effect.orient_lattice();
  (void)lattice;
}
