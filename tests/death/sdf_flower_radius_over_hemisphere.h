/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

/**
 * @brief Death case: a flower wider than a hemisphere must trap.
 * @details SDF surface — the petal cap bound is taken about the antipode, so a
 *          radius past the hemisphere inverts the band it derives.
 */
inline void case_sdf_flower_radius_over_hemisphere() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Flower flower(b, opaque(1.5f), opaque(5), opaque(0.0f));
  if (flower.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a zero-radius flower must trap.
 * @details SDF surface — the petal parameter divides by the circumradius, so a
 *          zero radius hands every probe a non-finite distance.
 */
inline void case_sdf_flower_zero_radius() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Flower flower(b, opaque(0.0f), opaque(5), opaque(0.0f));
  if (flower.circumradius == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: baking a class LUT for a degenerate polygon must trap.
 * @details SDF class-LUT surface — fewer than three vertices leaves no closed
 *          boundary for the crossing test, so every sample would read as
 *          outside.
 */
inline void case_sdf_class_lut_too_few_vertices() {
  static const float poly_xy[4] = {-0.5f, -0.5f, 0.5f, -0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(poly_xy, opaque(2), opaque(4), grid, lut);
  if (lut.n == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: baking a class LUT on a single-cell grid must trap.
 * @details SDF class-LUT surface — the cell step divides by (n - 1), so a
 *          resolution below 2 makes the whole domain non-finite.
 */
inline void case_sdf_class_lut_grid_too_small() {
  static const float poly_xy[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  static int16_t grid[16];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(poly_xy, opaque(3), opaque(1), grid, lut);
  if (lut.n == opaque(42))
    std::printf("x");
}

/**
 * @brief Death case: binding a class LUT at a vertex offset outside the face
 *        must trap.
 * @details SDF class-LUT surface — the offset indexes the canonical polygon
 *          cyclically, so an out-of-range one correlates the face against
 *          storage past the shape.
 */
inline void case_sdf_bind_class_lut_offset_out_of_range() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[3];
  uint16_t idx[3];
  for (int i = 0; i < 3; ++i) {
    float a = (2.0f * math::PI_F * i) / 3.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }
  static SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const math::Vector>(verts, 3),
                 std::span<const uint16_t>(idx, 3), scratch, HV, H);
  static const float canon_xy[6] = {-0.5f, -0.5f, 0.5f, -0.5f, 0.0f, 0.5f};
  SDF::ClassLut lut;
  if (face.bind_class_lut(&lut, canon_xy, opaque(7), false))
    std::printf("x");
}

/**
 * @brief Death case: a ring wider than the antipode must trap.
 * @details SDF ring surface — target_angle is radius * PI/2, so past 2 the
 *          band's cosine limits wrap and the stroke lands at the wrong
 *          latitude.
 */
inline void case_sdf_ring_radius_past_antipode() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Ring ring(b, opaque(2.5f), opaque(0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a ring with a negative stroke half-width must trap.
 * @details SDF ring surface — a negative thickness inverts the angular band,
 *          so every probe returns the far sentinel and the ring renders
 *          nothing.
 */
inline void case_sdf_ring_negative_thickness() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Ring ring(b, opaque(1.0f), opaque(-0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a distorted ring wider than the antipode must trap.
 * @details SDF ring surface — the shared ring geometry derives its band from
 *          radius * PI/2, which past 2 wraps its cosine limits.
 */
inline void case_sdf_distorted_ring_radius_past_antipode() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::FlatDistortedRing ring(b, opaque(2.5f), opaque(0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a distorted ring with a negative half-width must trap.
 * @details SDF ring surface — a negative thickness inverts the angular band,
 *          so every probe returns the far sentinel and the ring renders
 *          nothing.
 */
inline void case_sdf_distorted_ring_negative_thickness() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::FlatDistortedRing ring(b, opaque(0.5f), opaque(-0.05f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a line with a negative stroke half-width must trap.
 * @details SDF line surface — a negative thickness inverts the angular band
 *          and shrinks the bounding cap below the arc's own half-length, so
 *          the cull drops rows the arc covers.
 */
inline void case_sdf_line_negative_thickness() {
  SDF::Line line(math::Vector(1, 0, 0), math::Vector(0, 0, 1), opaque(-0.05f));
  if (line.thickness == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a distorted ring built with a null shift callback must
 *        trap.
 * @details SDF ring surface — the callback is invoked per azimuth on every
 *          probe, so a null one faults deep inside the rasterizer instead.
 */
inline void case_sdf_distorted_ring_null_shift() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  ScalarFn shift; // default-constructed -> empty
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), shift, opaque(0.1f),
                          opaque(0.0f));
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

/** @brief Death case: contour preparation past its table capacity traps. */
inline void case_shapeshifter_count_over_capacity() {
  using namespace shapeshifter_oracle_tests;
  OracleEffect effect;
  ShapeShifterWhiteBox::prepare_count(effect,
                                      opaque(OracleEffect::MAX_SHAPES + 1));
}

/** @brief Death case: planar chord storage for no vertices traps. */
inline void case_planar_chords_empty_storage() {
  static uint8_t arena_buf[256];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, opaque(0));
}

/** @brief Death case: a polyline past the bound chord storage traps. */
inline void case_planar_chords_over_capacity() {
  static uint8_t arena_buf[256];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.prepare(canvas);
  chords.prepare(canvas.clip());
  Fragments points;
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, opaque(3), math::Basis{}, Color4{},
                     shader);
}

/** @brief Death case: a chord draw before prepare traps. */
inline void case_planar_chords_unprepared() {
  static uint8_t arena_buf[256];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.prepare(canvas);
  Fragments points;
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, 2, math::Basis{}, Color4{}, shader);
}

inline void case_planar_chords_stale_clip() {
  static uint8_t arena_buf[256];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarChords<16, 8> chords;
  chords.init_storage(arena, 2);
  static hs_test::StubEffect fx(16, 8);
  Canvas canvas(fx);
  Filter::Screen::DirectAntiAliasSink<16, 8> sink;
  sink.prepare(canvas);
  ClipRegion clip = canvas.clip();
  clip.y_start = opaque(1);
  chords.prepare(clip);
  Fragments points;
  auto shader = [](const math::Vector &, Fragment &) {};
  chords.draw_closed(sink, canvas, points, 2, math::Basis{}, Color4{}, shader);
}

/** @brief Death case: band-split storage for under two points traps. */
inline void case_planar_band_split_empty_storage() {
  static uint8_t arena_buf[64];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, opaque(1));
}

/** @brief Death case: a band split cannot append to a populated destination. */
inline void case_planar_band_split_nonempty_output() {
  static uint8_t arena_buf[256];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out;
  out.bind(arena, 5);
  out.push_back(Fragment{});
  Fragments ring;
  band_split.split(out, ring, opaque(1), 1, math::Basis{},
                   Plot::ClipBand<16, 8>{});
}

/** @brief Death case: a band split past its bound storage traps. */
inline void case_planar_band_split_over_capacity() {
  static uint8_t arena_buf[64];
  Arena arena(arena_buf, sizeof(arena_buf));
  Plot::PlanarBandSplit<16, 8> band_split;
  band_split.init_storage(arena, 5);
  Fragments out;
  Fragments ring;
  band_split.split(out, ring, opaque(2), 4, math::Basis{},
                   Plot::ClipBand<16, 8>{});
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

/** @brief Death case: reconcile endpoints with different sizes must trap. */
inline void case_reconcile_vertices_size_mismatch() {
  static uint8_t identity_buf[64];
  static uint8_t target_buf[64];
  static uint8_t scratch_buf[64];
  Arena identity_arena(identity_buf, sizeof(identity_buf));
  Arena target(target_buf, sizeof(target_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh identity;
  identity.vertices.bind(identity_arena, 1);
  identity.vertices.push_back(math::Vector{});
  PolyMesh authored;
  PolyMesh out;
  MeshOps::reconcile_vertices(identity, authored, out, target, scratch);
}

/** @brief Rejects output aliasing a reconciliation input. */
inline void case_reconcile_vertices_aliased_output() {
  static uint8_t target_bytes[128], scratch_bytes[128];
  Arena target(target_bytes, sizeof(target_bytes));
  Arena scratch(scratch_bytes, sizeof(scratch_bytes));
  PolyMesh identity, authored;
  MeshOps::reconcile_vertices(identity, authored, identity, target, scratch);
}

/** @brief Death case: reconciling a vertex-less endpoint pair must trap. */
inline void case_reconcile_vertices_empty() {
  static uint8_t target_buf[64];
  static uint8_t scratch_buf[64];
  Arena target(target_buf, sizeof(target_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));
  PolyMesh identity;
  PolyMesh authored;
  PolyMesh out;
  MeshOps::reconcile_vertices(identity, authored, out, target, scratch);
}

/**
 * @brief Death case: an operator table whose entry decreases carrier family
 *        rank must trap.
 * @details Interpreter surface — compile() only matches adjacent carriers, so
 *          a rank-decreasing entry would run a chain the type system forbids.
 */
inline void case_chain_table_rank_decreases() {
  static Pullback::Interp::ChainProgram program;
  static Pullback::Interp::OperatorDescriptor descriptor{};
  descriptor.input = Pullback::Interp::CarrierId::COLOR;
  descriptor.output = Pullback::Interp::CarrierId::SPHERE;
  alignas(std::max_align_t) static uint8_t block_a[64];
  alignas(std::max_align_t) static uint8_t block_b[64];
  program.bind_storage(
      block_a, block_b, opaque<size_t>(sizeof(block_a)),
      std::span<const Pullback::Interp::OperatorDescriptor>(&descriptor, 1));
}

inline void case_star_mismatched_chart() {
  using Star = Plot::Star<Plot::PlanarProjection>;
  Fragments points;
  float x[10], y[10];
  const auto basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  const auto wrong = math::make_basis(math::Quaternion(), math::Y_AXIS);
  Star::sample_chart_positions(points, x, y, basis, 0.5f, 5, 0.0f,
                               Star::radius_trig(0.5f), Star::step_trig(5),
                               wrong);
}

inline void case_sdf_distorted_ring_negative_distortion() {
  const math::Basis basis{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  ScalarFn shift = [](float) { return 0.0f; };
  SDF::DistortedRing ring(basis, opaque(0.5f), opaque(0.05f), shift,
                          opaque(-0.1f), 0.0f);
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

inline void case_chain_overlapping_storage() {
  Pullback::Interp::ChainProgram program;
  alignas(std::max_align_t) uint8_t block[128];
  program.bind_storage(block, block + alignof(std::max_align_t), 64);
}

inline void case_chain_identical_storage() {
  Pullback::Interp::ChainProgram program;
  alignas(std::max_align_t) uint8_t block[64];
  program.bind_storage(block, block, sizeof(block));
}

inline void chain_invalid_layout(unsigned variant) {
  using namespace Pullback::Interp;
  ChainProgram program;
  OperatorDescriptor descriptor = OPERATOR_TABLE[0];
  alignas(std::max_align_t) uint8_t block_a[512], block_b[512];
  size_t capacity = sizeof(block_a);
  switch (variant) {
  case 0:
    descriptor.runtime.param.align = 0;
    break;
  case 1:
    descriptor.runtime.prepared.align = 3;
    break;
  case 2:
    descriptor.runtime.state.align = 2 * alignof(std::max_align_t);
    break;
  case 3:
    descriptor.runtime.param.size = 0;
    break;
  case 4:
    descriptor.runtime.state.size = 3;
    break;
  case 5:
    capacity = std::numeric_limits<uint32_t>::max();
    break;
  }
  program.bind_storage(block_a, block_b, capacity, {&descriptor, 1});
}
inline void case_chain_zero_alignment() { chain_invalid_layout(0); }
inline void case_chain_non_power_alignment() { chain_invalid_layout(1); }
inline void case_chain_overaligned_block() { chain_invalid_layout(2); }
inline void case_chain_zero_size() { chain_invalid_layout(3); }
inline void case_chain_misaligned_size() { chain_invalid_layout(4); }
inline void case_chain_capacity_overflow() { chain_invalid_layout(5); }

inline void case_pullback_project_nonunit_direction() {
  using namespace hs_test::pullback_tests;
  using Project =
      Pullback::Stage::Project<CountingProjectionPolicy>::Bind<TestBinding>;
  const TestFrame frame;
  const Pullback::SphereSample input{math::Vector(opaque(2.0f), 0.0f, 0.0f),
                                     0.0f};
  const auto result = Project::run(input, frame, Project::prepare(frame));
  if (result.coords.re != 0.0f)
    std::printf("x");
}

/** @brief Death case: a Sample operator rejects an unknown coverage mode. */
inline void case_pullback_operator_invalid_coverage_mode() {
  Pullback::Interp::Op::GridSampleParams params;
  params.coverage_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::SourceClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::SampleGridV3::prepare(context, params, state)
          .primary != 0.0f)
    std::printf("x");
}

/** @brief Death case: a Sample operator rejects an unknown weight mode. */
inline void case_pullback_operator_invalid_weight_mode() {
  Pullback::Interp::Op::GridSampleParams params;
  params.weight_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::SourceClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::SampleGridV3::prepare(context, params, state)
          .primary != 0.0f)
    std::printf("x");
}

/** @brief Death case: a warp operator rejects an unknown envelope. */
inline void case_pullback_operator_invalid_warp_envelope() {
  Pullback::Interp::Op::WaveShearWarpParams params;
  params.envelope = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::WarpPhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpWaveShear::prepare(context, params, state)
          .phase != 0.0f)
    std::printf("x");
}

/** @brief Death case: the polar chart rejects an unknown polar mode. */
inline void case_pullback_operator_invalid_polar_mode() {
  Pullback::Interp::Op::PolarChartParams params;
  params.mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::WarpPhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpPolarChart::prepare(context, params, state)
          .phase != 0.0f)
    std::printf("x");
}

/** @brief Death case: the polar chart rejects an out-of-range harmonic. */
inline void case_pullback_operator_invalid_polar_harmonic() {
  Pullback::Interp::Op::PolarChartParams params;
  params.harmonic = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::WarpPhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpPolarChart::prepare(context, params, state)
          .phase != 0.0f)
    std::printf("x");
}

/** @brief Death case: the curl-flow operator rejects an unknown integrator. */
inline void case_pullback_operator_invalid_curl_integrator() {
  Pullback::Interp::Op::CurlFlowWarpParams params;
  params.integrator = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::NoisePhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpCurlFlow::prepare(context, params, state)
          .intervals != 0)
    std::printf("x");
}

/** @brief Death case: the curl displacement rejects an unknown integrator. */
inline void case_pullback_operator_invalid_surface_integrator() {
  Pullback::Interp::Op::CurlDisplaceParams params;
  params.integrator = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::NoisePhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::DisplaceCurl::prepare(context, params, state)
          .noise != nullptr)
    std::printf("x");
}

/** @brief Rejects an unknown Bonne hemisphere. */
inline void case_pullback_operator_invalid_bonne_hemisphere() {
  Pullback::Interp::Op::ProjectBonneV3::Params params;
  params.hemisphere = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ProjectBonneV3::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectBonneV3::prepare(context, params, state);
}

/** @brief Death case: the gnomonic projection rejects an unknown hemisphere. */
inline void case_pullback_operator_invalid_gnomonic_hemisphere() {
  Pullback::Interp::Op::GnomonicChainParams params;
  params.hemisphere = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ProjectGnomonic::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectGnomonic::prepare(context, params, state);
}

inline void case_pullback_operator_invalid_airocean_layout() {
  Pullback::Interp::Op::ProjectAiroceanV3::Params params;
  params.layout = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ProjectAiroceanV3::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectAiroceanV3::prepare(context, params, state);
}

inline void case_pullback_operator_invalid_peirce_layout() {
  Pullback::Interp::Op::ProjectPeirceV3::Params params;
  params.layout =
      opaque<uint8_t>(std::size(Pullback::Interp::Op::PEIRCE_LAYOUT_IDS));
  Pullback::Interp::Op::ProjectPeirceV3::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectPeirceV3::prepare(context, params, state);
}

/** @brief Death case: a noise-driven operator rejects an unknown basis. */
inline void case_pullback_operator_invalid_noise_basis() {
  Pullback::Interp::Op::CurlDisplaceParams params;
  params.basis = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::NoisePhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::DisplaceCurl::prepare(context, params, state)
          .noise != nullptr)
    std::printf("x");
}

/** @brief Death case: the tessellation source rejects an unknown kind. */
inline void case_pullback_operator_invalid_tessellation_kind() {
  Pullback::Interp::Op::TessellationSampleParams params;
  params.kind = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::SourceClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::SampleTessellation::prepare(context, params, state)
          .primary != 0.0f)
    std::printf("x");
}

/** @brief Death case: the kaleidoscope lens rejects an unknown symmetry. */
inline void case_pullback_operator_invalid_kaleidoscope_symmetry() {
  Pullback::Interp::Op::KaleidoscopeChainParams params;
  params.symmetry = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::LensKaleidoscope::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::LensKaleidoscope::prepare(context, params, state);
}

/** @brief Death case: the generated-palette operator rejects an unknown hue mode. */
inline void case_pullback_operator_invalid_hue_mode() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.hue_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
}

/** @brief Death case: the generated-palette operator rejects an unknown palette
    mode. */
inline void case_pullback_operator_invalid_palette_mode() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.palette_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
}

/** @brief Death case: the generated-palette operator rejects an unknown palette
    mapping. */
inline void case_pullback_operator_invalid_palette_mapping() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.mapping_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
}
