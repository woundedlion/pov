/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

/**
 * @brief Death case: a class name equal to another effect's stable ID must
 *        trap.
 */
inline void case_effect_registry_name_matches_stable_id() {
  EffectRegistration first{};
  first.name = "DeathFirst";
  first.stable_id = "DeathPersistedAlias";

  EffectRegistration second{};
  second.name = "DeathPersistedAlias";
  second.stable_id = "death-second";
  validate_effect_registrations(std::array{first, second});
}

/**
 * @brief Death case: a Flywheel period of zero must trap at construction.
 * @details POV-sync surface — position() divides the int32 elapsed window by the
 *          period, so a zero divides by zero and an over-large one voids the
 *          signed-safe coast window; the constructor rejects both before the
 *          driver ever schedules a column.
 */
inline void case_flywheel_period_zero() {
  pov::sync::Config cfg;
  cfg.cycles_per_half_rev = opaque<uint32_t>(0);
  pov::sync::Flywheel fw(cfg); // period 0 -> HS_CHECK
  (void)fw;
}

/**
 * @brief Death case: a virtual height of one row must trap in the phi mapping.
 * @details Geometry surface — the row-to-angle scale divides by (h_virt - 1),
 *          so a single-row canvas would map every row to a non-finite phi.
 */
inline void case_y_to_phi_degenerate_height() {
  float phi =
      math::y_to_phi_virtual(opaque(0.0f), opaque(1)); // divisor 0 -> trap
  if (phi == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: reading an orientation frame past the history must trap.
 * @details Geometry surface — the motion-blur history is a fixed array whose
 *          live prefix is num_frames long, so an index past it would read a
 *          stale or never-written quaternion instead of failing.
 */
inline void case_orientation_frame_index_oob() {
  math::Orientation<> orientation; // constructed with one frame
  const math::Quaternion &q = orientation.get(opaque(3));
  if (q.r == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: make_basis with a non-unit quaternion must trap.
 * @details Geometry surface — the rotation assumes a unit quaternion, so a
 *          finite but over-long one would scale and shear the frame rather than
 *          rotate it; the guard fires before the axes are built.
 */
inline void case_make_basis_nonunit_quaternion() {
  math::Quaternion q(opaque(2.0f), opaque(0.0f), opaque(0.0f), opaque(0.0f));
  math::Vector normal{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Basis b = math::make_basis(q, normal); // |q| = 2 -> HS_CHECK
  if (b.u.x == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: parallel transport between antipodal endpoints must trap.
 * @details Geometry surface — the great circle through antipodes is
 *          ill-determined and the transport divides by 1 + dot, so the guard
 *          fires before the tangent is amplified.
 */
inline void case_parallel_transport_antipodal() {
  math::Vector from{opaque(1.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector to{opaque(-1.0f), opaque(0.0f), opaque(0.0f)};
  math::Vector tangent{opaque(0.0f), opaque(1.0f), opaque(0.0f)};
  math::Vector t =
      math::parallel_transport(from, to, tangent); // dot = -1 -> HS_CHECK
  if (t.x == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a polyhedral fold that never converges must trap.
 * @details Lens surface — two opposed mirrors are not a chamber: each pass
 *          reflects the direction back across the other, so the bounded
 *          reflection loop exhausts its passes and fires the guard.
 */
inline void case_polyhedral_kaleidoscope_no_converge() {
  const std::array<math::Vector, 3> mirrors = {
      math::Vector(opaque(1.0f), 0.0f, 0.0f),
      math::Vector(opaque(-1.0f), 0.0f, 0.0f),
      math::Vector(0.0f, opaque(1.0f), 0.0f)};
  math::Vector v = lenses::polyhedral_kaleidoscope_lens(
      math::Vector(opaque(0.5f), opaque(0.5f), 0.0f), mirrors);
  if (v.x == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a polygon with fewer than three sides must trap.
 * @details SDF surface — the sector fold divides a full turn by the side count,
 *          so a 2-gon has no interior for the distance to be measured against.
 */
inline void case_sdf_polygon_side_count() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::PlanarPolygon poly(b, opaque(0.5f), opaque(2), opaque(0.0f));
  if (poly.apothem == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: an angular repeat around a non-unit axis must trap.
 * @details SDF surface — the sector fold rotates the query point about the
 *          axis, so a non-unit one scales every folded copy off the sphere.
 */
inline void case_sdf_angular_repeat_nonunit_axis() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::Ring ring(b, opaque(1.0f), opaque(0.1f));
  math::Vector axis{opaque(0.0f), opaque(2.0f), opaque(0.0f)};
  SDF::AngularRepeat<SDF::Ring> rep(ring, opaque(4), axis); // non-unit -> trap
  if (rep.sector == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a knot ring with no cells must trap.
 * @details SDF surface — the per-pixel cell index divides the azimuth by
 *          2π/n, so n == 0 wraps to knots[-1] on every probe.
 */
inline void case_sdf_distorted_ring_zero_knots() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[1] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), knots, opaque(0),
                          opaque(0.0f), pf); // no knot cells -> trap
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

inline void case_sdf_distorted_ring_one_knots() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[1] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), knots, opaque(1),
                          opaque(0.0f), pf);
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

inline void case_sdf_distorted_ring_two_knots() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[2] = {0.0f};
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, opaque(0.5f), opaque(0.05f), knots, opaque(2),
                          opaque(0.0f), pf);
  if (ring.thickness == opaque(42.0f))
    std::printf("x");
}

inline void case_gamut_lut_scratch_a() {
  init_gamut_lut(scratch_arena_a, GAMUT_LUT_MIN_ANGLE_STEPS,
                 GAMUT_LUT_MIN_L_STEPS);
}

inline void case_gamut_lut_scratch_b() {
  init_gamut_lut(scratch_arena_b, GAMUT_LUT_MIN_ANGLE_STEPS,
                 GAMUT_LUT_MIN_L_STEPS);
}

inline void case_scan_ring_stack_too_many_slots() {
  constexpr int W = 32, H = 16;
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[4] = {0.0f, 0.01f, 0.0f, -0.01f};
  SDF::DistortedRing ring(b, 1.0f, 0.05f, knots, 4, 0.0f, nullptr);
  const int8_t slot_by_ring[1] = {0};
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, 1, &ring, slot_by_ring, opaque(INT8_MAX + 1), table,
      [](int, const math::Vector &, Fragment &) {});
}

inline void case_scan_ring_stack_too_many_rings() {
  constexpr int W = 32, H = 16;
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  const float knots[4] = {0.0f, 0.01f, 0.0f, -0.01f};
  SDF::DistortedRing ring(b, 1.0f, 0.05f, knots, 4, 0.0f, nullptr);
  int8_t slot_by_ring[256];
  for (int8_t &s : slot_by_ring)
    s = -1;
  slot_by_ring[0] = 0;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, opaque(256), &ring, slot_by_ring, 1, table,
      [](int, const math::Vector &, Fragment &) {}); // 256 > 255 -> trap
}

inline void case_scan_ring_stack_callback_ring() {
  constexpr int W = 32, H = 16;
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::DistortedRing ring(
      b, 1.0f, 0.05f, [](float) { return 0.0f; }, opaque(0.01f), 0.0f);
  const int8_t slot_by_ring[1] = {0};
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  Canvas canvas(fx);
  static Scan::DistortedRingStack::CandidateTable<W, H> table;
  Scan::DistortedRingStack::draw<W, H>(
      pipeline, canvas, 1, &ring, slot_by_ring, 1, table,
      [](int, const math::Vector &, Fragment &) {}); // no knots -> trap
}

/**
 * @brief Death case: a twist warp around a zero-radius torus must trap.
 * @details SDF warp surface — the Lipschitz bound scales by 2/R, so a zero
 *          major radius hands the rasterizer a non-finite step bound.
 */
inline void case_sdf_twist_zero_major_radius() {
  SDF::Warp::Twist tw(opaque(2), opaque(0.1f), opaque(0.0f)); // R = 0 -> trap
  if (tw.two_over_r == opaque(42.0f))
    std::printf("x");
}

inline void case_transformed_torus_invalid_minor_radius() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
      {opaque(1.0f), opaque(0.6f)}, {2, 0.1f, 1.0f}};
  Scan::TransformedVolume volume(torus, math::Vector(), math::Quaternion());
  volume.check_trace_preconditions();
}

/** @brief Draw callback for the OpLeg construction death cases; never runs. */
inline void death_opleg_draw(Canvas &, MeshState &,
                             const Animation::OpLeg::Shading &) {}

/**
 * @brief A palette handoff complete enough to clear the OpLeg handoff guard.
 * @return A handoff naming a default bank and a one-face departed palette.
 * @details Lets a case reach a guard that sits behind the handoff check; the
 *          LUTs are never sampled, since every such case traps in the
 *          constructor before the first frame.
 */
inline Animation::OpLeg::PaletteHandoff death_opleg_handoff() {
  static const BakedPaletteBank bank;
  static const uint8_t face_palette[1] = {0};
  return {.bank = &bank,
          .prev_face_palette = face_palette,
          .prev_faces = 1,
          .prev_face_centroid = nullptr,
          .correspondence = Animation::OpLeg::FaceCorrespondence::GEOMETRIC};
}

inline void case_opleg_rewind_refill() {
  static uint8_t seed_storage[64 * 1024], leg_storage[128 * 1024];
  Arena seed_arena(seed_storage, sizeof(seed_storage));
  Arena leg_arena(leg_storage, sizeof(leg_storage));
  PolyMesh seed;
  build_solid<Solids::Cube>(seed, seed_arena);
  static const BakedPaletteBank bank;
  static const uint8_t palettes[6]{};
  const Animation::OpLeg::PaletteHandoff HANDOFF{
      .bank = &bank, .prev_face_palette = palettes, .prev_faces = 6};
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::GatedSwapSpec{.op = Animation::OpLeg::SwapOp::KIS,
                                      .gate_frames = 1},
      leg_arena, death_opleg_draw, HANDOFF);
  const size_t BYTES = leg_arena.get_offset();
  leg_arena.set_offset(0);
  leg_arena.allocate(BYTES);
  (void)leg.landing();
}

/**
 * @brief Death case: choosing an edge from a node outside the graph must trap.
 * @details ConwayGraph surface — no EDGES row touches such a node, so the
 *          weighted pick would have nothing to divide by.
 */
inline void case_pick_next_edge_unknown_node() {
  uint8_t visits[ConwayGraph::NUM_NODES] = {};
  const int e = ConwayGraph::pick_next_edge(opaque<int>(ConwayGraph::NUM_NODES),
                                            -1, 0, visits, 0u);
  if (e == opaque<int>(-42))
    std::printf("x");
}

/** @brief A well-formed non-settling graph edge for the OpLeg death cases. */
inline constexpr ConwayGraph::EdgeSpec death_opleg_edge{
    .from_node = 0,
    .to_node = 0,
    .seed_solid = 0,
    .op = ConwayGraph::MorphOp::TRUNCATE,
    .t_from = 0.0f,
    .t_to = 0.4f,
    .twist_from = 0.0f,
    .twist_to = 0.0f,
    .settle = false,
    .reseed = ConwayGraph::Reseed::NONE,
    .bridge = false};

/**
 * @brief Death case: an edge-sweep leg without a graph edge must trap.
 * @details OpLeg surface — the constructor reads the edge's operator and settle
 *          flag on its first line, so a null edge is a null dereference.
 */
inline void case_opleg_edge_sweep_no_edge() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed, Animation::OpLeg::EdgeSweepSpec{}, arena,
                       death_opleg_draw, handoff); // null edge -> HS_CHECK
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a leg with a non-positive sweep length must trap.
 * @details OpLeg surface — the per-frame sweep parameter divides by the frame
 *          count, and a zero-frame leg would also complete before drawing.
 */
inline void case_opleg_zero_sweep_frames() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                       .sweep_frames = opaque(0)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a leg built without a palette handoff must trap.
 * @details OpLeg surface — the departed node's per-face palette keys every
 *          blend ramp the leg bakes, so an absent bank leaves the whole
 *          crossfade unresolvable rather than merely uncolored.
 */
inline void case_opleg_incomplete_palette_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff; // no bank, no per-face palette
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                       .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: settle frames that contradict the edge must trap.
 * @details OpLeg surface — the edge's settle flag decides whether the leg
 *          computes a relaxed endpoint at all, so a settle window on a
 *          non-settling edge would slerp toward vertices nothing produced.
 */
inline void case_opleg_edge_settle_mismatch() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::EdgeSweepSpec{
                           .edge = &death_opleg_edge,
                           .reverse = false,
                           .sweep_frames = opaque(1),
                           .settle_frames = opaque(1)}, // edge does not settle
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: an edge-sweep leg without a palette handoff must trap.
 * @details Pins the edge-sweep constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_edge_sweep_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::EdgeSweepSpec{.edge = &death_opleg_edge,
                                      .reverse = false,
                                      .sweep_frames = opaque(1),
                                      .settle_frames = opaque(0)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a hankin leg without a palette handoff must trap.
 * @details Pins the hankin constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_hankin_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::HankinSweepSpec{.theta_start = opaque(0.1f),
                                        .theta_end = opaque(0.5f),
                                        .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a hankin leg sweeping to a smaller angle must trap.
 * @details OpLeg surface — the leg sweeps the slerp fraction outward from the
 *          collapsed corner, which is monotone only while the arrival angle is
 *          the larger of the two.
 */
inline void case_opleg_hankin_backward_theta() {
  static uint8_t buf[8192];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  const Animation::OpLeg::PaletteHandoff handoff = death_opleg_handoff();
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::HankinSweepSpec{.theta_start = opaque(0.9f),
                                        .theta_end = opaque(0.3f),
                                        .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a relax leg with neither a bake nor iterations must trap.
 * @details OpLeg surface — the leg needs a relaxed endpoint to slerp to, which
 *          is either the shipped bake or the result of live iterations; with
 *          neither it would slerp the seed onto itself for its whole duration.
 */
inline void case_opleg_relax_no_iterations() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::RelaxSpec{.iterations = opaque(0),
                                                   .bake = nullptr,
                                                   .sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a relax leg without a palette handoff must trap.
 * @details Pins the relax constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_relax_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::RelaxSpec{.iterations = opaque(1),
                                                   .bake = nullptr,
                                                   .sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a medial leg without a palette handoff must trap.
 * @details Pins the medial constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_medial_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::MedialSpec{.sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a reconcile leg without endpoints must trap.
 * @details OpLeg surface — the leg slerps every seed vertex to an authored
 *          position, so an absent endpoint array is the whole leg's target.
 */
inline void case_opleg_reconcile_no_endpoints() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(seed,
                       Animation::OpLeg::ReconcileSpec{
                           .to_positions = nullptr, .sweep_frames = opaque(1)},
                       arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/** @brief Death case: reconcile endpoint and seed counts must agree. */
inline void case_opleg_reconcile_endpoint_count() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  const math::Vector endpoint{};
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::ReconcileSpec{.to_positions = &endpoint,
                                      .to_count = opaque<size_t>(1),
                                      .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a reconcile leg without a palette handoff must trap.
 * @details Pins the reconcile constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_reconcile_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  static const math::Vector endpoints[1] = {math::Vector(0, 0, 1)};
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::ReconcileSpec{.to_positions = endpoints,
                                      .sweep_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a gated-swap leg with no gate window must trap.
 * @details OpLeg surface — the leg runs 2*gate_frames + 1 frames around the
 *          swap, so a zero gate leaves the swap frame with no approach or
 *          departure to blend across.
 */
inline void case_opleg_gated_swap_zero_gate_frames() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::GatedSwapSpec{.op = Animation::OpLeg::SwapOp::KIS,
                                      .gate_frames = opaque(0)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a gated-swap leg without a palette handoff must trap.
 * @details Pins the gated-swap constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_gated_swap_incomplete_handoff() {
  static uint8_t buf[1024];
  Arena arena(buf, sizeof(buf));
  PolyMesh seed;
  Animation::OpLeg::PaletteHandoff handoff;
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::GatedSwapSpec{.op = Animation::OpLeg::SwapOp::KIS,
                                      .gate_frames = opaque(1)},
      arena, death_opleg_draw, handoff);
  if (leg.landing().faces == opaque<size_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: a shading lookup past the leg's face table must trap.
 * @details OpLeg surface — the ramp index is read straight from the per-face
 *          table, so an out-of-range face would shade through whatever follows
 *          it instead of failing.
 */
inline void case_opleg_shading_face_out_of_range() {
  static BakedPalette ramps[1];
  static const uint8_t face_ramp[1] = {0};
  const Animation::OpLeg::Shading shading{
      .ramps = ramps, .face_ramp = face_ramp, .faces = opaque<size_t>(1)};
  if (&shading.ramp_for(opaque<size_t>(3)) == ramps)
    std::printf("x");
}

/**
 * @brief Death case: a Motion over an empty path must trap on the first step.
 * @details An unfilled Path samples the origin. Motion rejects its zero vectors
 *          before computing an angle or updating Orientation.
 */
inline void case_motion_empty_path_origin_sample() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  math::Orientation<4> orientation;
  static Path<32> path; // never appended -> get_point returns the origin
  Animation::Motion<W, 4> motion(orientation, path, opaque(10));
  motion.step(c);
}

/**
 * @brief Death case: a negative equator sample count must trap.
 * @details Spherical-field surface — the count sizes every ring's longitude
 *          walk, so a negative one underflows the per-ring sample allocation.
 */
inline void case_spherical_field_negative_equator_samples() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4, 0, 0, opaque(-1));
  if (layout.sample_count() == 42)
    std::printf("x");
}

/**
 * @brief Death case: a spherical polygon wider than a hemisphere must trap.
 * @details SDF surface — beyond the hemisphere the cap fold changes sign, so
 *          the shape must be built inverted about its antipode instead.
 */
inline void case_sdf_spherical_polygon_radius_over_hemisphere() {
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                      math::Vector(0, 0, 1)};
  SDF::SphericalPolygon poly(b, opaque(1.5f), opaque(5), opaque(0.0f));
  if (poly.circumradius == opaque(42.0f))
    std::printf("x");
}
