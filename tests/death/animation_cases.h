/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Animation death cases.

/** @brief Draw callback for the OpLeg construction death cases; never runs. */
inline void death_opleg_draw(Canvas &, MeshState &,
                             const Animation::OpLeg::Shading &) {}

/**
 * @brief Constructs an OpLeg over an empty seed for death probes.
 * @tparam CAPACITY Scratch arena capacity in bytes.
 * @param spec Constructor specification.
 * @param handoff Palette handoff supplied to the constructor.
 */
template <size_t CAPACITY = 1024, typename Spec>
inline void
death_construct_opleg(const Spec &spec,
                      Animation::OpLeg::PaletteHandoff handoff = {}) {
  static uint8_t buffer[CAPACITY];
  Arena arena(buffer, sizeof(buffer));
  PolyMesh seed;
  Animation::OpLeg leg(seed, spec, arena, death_opleg_draw, handoff);
  (void)leg.landing();
}

/**
 * @brief A palette handoff complete enough to clear the OpLeg handoff guard.
 * @return A handoff naming a default bank and a one-face departed palette.
 * @details Lets a case reach a guard that sits behind the handoff check.
 */
inline Animation::OpLeg::PaletteHandoff death_opleg_handoff() {
  static const BakedPaletteBank bank;
  static const uint8_t FACE_PALETTE[1] = {0};
  return {.bank = &bank,
          .prev_face_palette = FACE_PALETTE,
          .prev_faces = 1,
          .prev_face_centroid = nullptr,
          .correspondence = Animation::OpLeg::FaceCorrespondence::GEOMETRIC};
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

/** @brief Adds an event from a clear hook, violating the hook contract. */
inline void add_event_from_clear_hook(void *ctx) {
  static float value = 0.0f;
  static_cast<Timeline *>(ctx)->add(
      0, Animation::Transition(value, 1.0f, 1, math::ease_linear));
}

/** @brief Lens clocks reject nonfinite multipliers. */
inline void case_lens_phase_nonfinite_scale() {
  float angle = 0.0f, magnitude = 1.0f;
  Animation::LensPhase phase(angle, magnitude, opaque(NAN), 16);
}

/**
 * @brief Death case: relocating a retained (pinned) add_get() handle must trap.
 * @details Animation surface — TimelineEvent::move_into traps when the event
 *          was handed out via add_get(Pin::PINNED).
 */
inline void case_timeline_pinned_relocation() {
  TimelineEvent src;
  src.pinned = opaque(true); // as if handed out by add_get(Pin::PINNED)
  TimelineEvent dst;
  src.move_into(dst); // HS_CHECK(!pinned) -> trap
}

/**
 * @brief Death case: relocating into a slot that still owns an animation must
 *        trap.
 * @details Animation surface — move_into overwrites dst.manager/dst.iface, so a
 *          live destination would lose its animation's destructor.
 */
inline void case_timeline_move_into_live_destination() {
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 10, math::ease_linear));
  tl.add(0, Animation::Transition(v, 1.0f, 10, math::ease_linear));
  global_timeline_events[opaque(0)].move_into(global_timeline_events[1]);
}

/**
 * @brief Death case: a negative timeline delay must trap.
 */
inline void case_timeline_negative_delay() {
  Timeline tl;
  float value = 0.0f;
  tl.add(opaque(-1), Animation::Transition(value, 1.0f, 1, math::ease_linear));
}

/**
 * @brief Death case: a pinned animation that COMPLETES must trap.
 * @details Animation surface — the pin contract is pinned => infinite, so
 *          step()'s completion branch traps a pinned animation that completes.
 *          cancel() is exempt (is_canceled()), so this case completes naturally.
 */
inline void case_timeline_pinned_completion() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  // add_get(Pin::PINNED) rejects a finite non-repeating animation, so the event
  // is marked pinned directly. A 1-frame Transition as the sole event goes
  // through completion/destroy.
  tl.add(0, Animation::Transition(v, 1.0f, 1, math::ease_linear));
  global_timeline_events[0].pinned = opaque(true);
  tl.step(canvas); // t=1: done() && !repeats() && !canceled, keep=false -> trap
}

/**
 * @brief Death case: pinning a finite, non-repeating animation must trap.
 * @details Animation surface — add_get(Pin::PINNED) promises the caller a pointer
 *          valid across frames, which only holds for an animation that never
 *          completes on its own.
 */
inline void case_timeline_pinned_finite_animation() {
  Timeline tl;
  float v = 0.0f;
  tl.add_get(0, Animation::Transition(v, 1.0f, opaque(1), math::ease_linear),
             Timeline::Pin::PINNED);
}

/** @brief Pinning from a callback rejects a destroyed predecessor. */
inline void case_timeline_pinned_after_destroyed_predecessor() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  float a = 0.0f, b = 0.0f;
  tl.add(0, Animation::Transition(a, 1.0f, 1, math::ease_linear));
  tl.add(0, Animation::Transition(b, 1.0f, 1, math::ease_linear).then([&]() {
    tl.add_get(0,
               Animation::PeriodicTimer(1, [](Canvas &) {}, {.repeat = true}),
               Timeline::Pin::PINNED);
  }));
  tl.step(canvas);
}

inline void case_timeline_pinned_after_cancelled_perpetual() {
  Timeline tl;
  auto *predecessor = tl.add_get(
      0, Animation::PeriodicTimer(1, [](Canvas &) {}, {.repeat = true}),
      Timeline::Pin::UNPINNED);
  tl.add_get(0, Animation::PeriodicTimer(1, [](Canvas &) {}, {.repeat = true}),
             Timeline::Pin::PINNED);
  death_admitted();
  predecessor->cancel();
  tl.add_get(0, Animation::PeriodicTimer(1, [](Canvas &) {}, {.repeat = true}),
             Timeline::Pin::PINNED);
}

/**
 * @brief Death case: dropping a pinned add on a full timeline must trap.
 * @details Animation surface — the capacity guard returns nullptr, which an
 *          add_get(Pin::PINNED) caller would retain across frames.
 */
inline void case_timeline_pinned_add_on_full_timeline() {
  Timeline tl;
  float sink = 0.0f;
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  tl.add_get(0, Animation::PeriodicTimer(1, [](Canvas &) {}, {.repeat = true}),
             opaque(Timeline::Pin::PINNED));
}

/**
 * @brief Death case: a pinned one-shot timer must trap when it fires.
 * @details Animation surface — a one-shot timer ends itself via finish(), not
 *          cancel(), so a pinned one hits step()'s completion guard.
 */
inline void case_timeline_pinned_one_shot_timer() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  tl.add(0, Animation::PeriodicTimer(1, [](Canvas &) {}, {.repeat = false}));
  global_timeline_events[0].pinned = opaque(true);
  tl.step(canvas); // t=1: fires, finish() -> done() && !canceled -> trap
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
  (void)q;
}

/** @brief Invalid orientation upsample count below one traps. */
inline void case_orientation_upsample_count_below_one() {
  math::Orientation<> orientation;
  orientation.upsample(opaque(0));
}

/** @brief Invalid orientation orient index oob traps. */
inline void case_orientation_orient_index_oob() {
  math::Orientation<> orientation;
  (void)orientation.orient(math::X_AXIS, opaque(3));
}

/** @brief Invalid orientation unorient index oob traps. */
inline void case_orientation_unorient_index_oob() {
  math::Orientation<> orientation;
  (void)orientation.unorient(math::X_AXIS, opaque(3));
}

/** @brief Invalid orientation set nonunit traps. */
inline void case_orientation_set_nonunit() {
  math::Orientation<> orientation;
  orientation.set(math::Quaternion(opaque(2.0f), 0, 0, 0));
}

/** @brief Invalid orientation push nonunit traps. */
inline void case_orientation_push_nonunit() {
  math::Orientation<> orientation;
  orientation.push(math::Quaternion(opaque(2.0f), 0, 0, 0));
}

/** @brief Invalid orientation history full traps. */
inline void case_orientation_history_full() {
  math::Orientation<> orientation;
  orientation.upsample(4);
  orientation.push(math::Quaternion());
}

/** @brief Invalid orientation set at index oob traps. */
inline void case_orientation_set_at_index_oob() {
  math::Orientation<> orientation;
  orientation.set_at(opaque(3), math::Quaternion());
}

/** @brief Invalid orientation set at nonunit traps. */
inline void case_orientation_set_at_nonunit() {
  math::Orientation<> orientation;
  orientation.set_at(0, math::Quaternion(opaque(2.0f), 0, 0, 0));
}

/** @brief Invalid orientation set at normalized index oob traps. */
inline void case_orientation_set_at_normalized_index_oob() {
  math::Orientation<> orientation;
  orientation.set_at_normalized(opaque(3), math::Quaternion());
}

/** @brief Reclaiming a live phase block traps before its next write. */
inline void case_mobiuswarp_evolving_rewind() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  math::MobiusParams params;
  const size_t start = persistent_arena.get_offset();
  Animation::MobiusWarpEvolving warp(params);
  warp.step(canvas);
  persistent_arena.set_offset(start);
  warp.step(canvas);
}

/** @brief A moved phase block retains its reclaim stamp after refill. */
inline void case_mobiuswarp_evolving_move_rewind_refill() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  math::MobiusParams params;
  const size_t start = persistent_arena.get_offset();
  Animation::MobiusWarpEvolving source(params);
  Animation::MobiusWarpEvolving moved(std::move(source));
  Animation::MobiusWarpEvolving assigned(params);
  assigned = std::move(moved);
  assigned.step(canvas);
  const size_t end = persistent_arena.get_offset();
  persistent_arena.set_offset(start);
  persistent_arena.allocate(end - start);
  assigned.step(canvas);
}

/** @brief A copied phase block retains its reclaim stamp after refill. */
inline void case_mobiuswarp_evolving_copy_rewind_refill() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  math::MobiusParams params;
  Animation::MobiusWarpEvolving source(params);
  const size_t start = persistent_arena.get_offset();
  Animation::MobiusWarpEvolving copied(source);
  Animation::MobiusWarpEvolving assigned(params);
  assigned = copied;
  assigned.step(canvas);
  const size_t end = persistent_arena.get_offset();
  persistent_arena.set_offset(start);
  persistent_arena.allocate(end - start);
  assigned.step(canvas);
}

inline void case_opleg_rewind_refill() {
  static uint8_t seed_storage[64 * 1024], leg_storage[128 * 1024];
  Arena seed_arena(seed_storage, sizeof(seed_storage));
  Arena leg_arena(leg_storage, sizeof(leg_storage));
  PolyMesh seed;
  build_solid<Solids::Cube>(seed, seed_arena);
  static const BakedPaletteBank bank;
  static const uint8_t PALETTES[6]{};
  const Animation::OpLeg::PaletteHandoff handoff{
      .bank = &bank, .prev_face_palette = PALETTES, .prev_faces = 6};
  Animation::OpLeg leg(
      seed,
      Animation::OpLeg::GatedSwapSpec{.op = Animation::OpLeg::SwapOp::KIS,
                                      .gate_frames = 1},
      leg_arena, death_opleg_draw, handoff);
  const size_t bytes = leg_arena.get_offset();
  leg_arena.set_offset(0);
  leg_arena.allocate(bytes);
  (void)leg.landing();
}

/**
 * @brief Death case: routing from a node outside the graph must trap.
 */
inline void case_next_edge_toward_unknown_node() {
  const int e = ConwayGraph::next_edge_toward(
      opaque<int>(ConwayGraph::NUM_NODES), ConwayGraph::TETRAHEDRON, 0u);
  (void)e;
}

/**
 * @brief Death case: an edge-sweep leg without a graph edge must trap.
 */
inline void case_opleg_edge_sweep_no_edge() {
  death_construct_opleg(Animation::OpLeg::EdgeSweepSpec{});
}

/**
 * @brief Death case: a leg with a non-positive sweep length must trap.
 * @details OpLeg surface — the per-frame sweep parameter divides by the frame
 *          count, and a zero-frame leg would also complete before drawing.
 */
inline void case_opleg_zero_sweep_frames() {
  death_construct_opleg(Animation::OpLeg::ParamSweepSpec{
      .op = ConwayGraph::MorphOp::TRUNCATE, .sweep_frames = opaque(0)});
}

/**
 * @brief Death case: a leg built without a palette handoff must trap.
 * @details OpLeg surface — the departed node's per-face palette keys every
 *          blend ramp the leg bakes, so an absent bank leaves the whole
 *          crossfade unresolvable rather than merely uncolored.
 */
inline void case_opleg_incomplete_palette_handoff() {
  death_construct_opleg(Animation::OpLeg::ParamSweepSpec{
      .op = ConwayGraph::MorphOp::TRUNCATE, .sweep_frames = opaque(1)});
}

/** @brief Prepared bridge data requires the opening geometric handoff. */
inline void case_opleg_incompatible_prepared_bridge() {
  using Animation::OpLeg;
  static uint8_t seed_storage[64 * 1024], leg_storage[128 * 1024];
  Arena seed_arena(seed_storage, sizeof(seed_storage));
  Arena leg_arena(leg_storage, sizeof(leg_storage));
  PolyMesh seed;
  build_solid<Solids::Cube>(seed, seed_arena);
  OpLeg::BridgePreparation prepared;
  const uint8_t PALETTES[6]{};
  auto handoff = death_opleg_handoff();
  handoff.prev_face_palette = PALETTES;
  handoff.prev_faces = std::size(PALETTES);
  handoff.prev_face_centroid = nullptr;
  OpLeg leg(OpLeg::SweepSeed::borrow(seed),
            OpLeg::ParamSweepSpec{.t_end = 0.5f,
                                  .sweep_frames = 1,
                                  .bridge_provenance = true,
                                  .prepared = &prepared},
            leg_arena, death_opleg_draw, handoff);
  (void)leg.landing();
}

/** @brief Prepared bridge hosts must name faces in the departed mesh. */
inline void case_opleg_invalid_prepared_bridge_host() {
  using Animation::OpLeg;
  static uint8_t seed_storage[64 * 1024], leg_storage[128 * 1024];
  Arena seed_arena(seed_storage, sizeof(seed_storage));
  Arena leg_arena(leg_storage, sizeof(leg_storage));
  PolyMesh seed;
  build_solid<Solids::Cube>(seed, seed_arena);
  OpLeg::BridgePreparation prepared;
  const uint16_t CLASSES[7]{};
  prepared.topology.bind(seed_arena, std::size(CLASSES));
  prepared.topology.append_bulk(CLASSES, std::size(CLASSES));
  prepared.hosts.fill(UINT16_MAX);

  const uint8_t PALETTES[6]{};
  const math::Vector CENTROIDS[6]{};
  auto handoff = death_opleg_handoff();
  handoff.prev_face_palette = PALETTES;
  handoff.prev_faces = std::size(PALETTES);
  handoff.prev_face_centroid = CENTROIDS;
  OpLeg leg(OpLeg::SweepSeed::borrow(seed),
            OpLeg::ParamSweepSpec{.t_end = 0.5f,
                                  .sweep_frames = 1,
                                  .bridge_provenance = true,
                                  .prepared = &prepared},
            leg_arena, death_opleg_draw, handoff);
  (void)leg.landing();
}

/**
 * @brief Death case: settle frames that contradict the edge must trap.
 * @details OpLeg surface — the edge's settle flag decides whether the leg
 *          computes a relaxed endpoint at all, so a settle window on a
 *          non-settling edge would slerp toward vertices nothing produced.
 */
inline void case_opleg_edge_settle_mismatch() {
  death_construct_opleg(
      Animation::OpLeg::EdgeSweepSpec{.edge = &death_opleg_edge,
                                      .reverse = false,
                                      .sweep_frames = opaque(1),
                                      .settle_frames = opaque(1)});
}

/**
 * @brief Death case: an edge-sweep leg without a palette handoff must trap.
 * @details Pins the edge-sweep constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_edge_sweep_incomplete_handoff() {
  death_construct_opleg(
      Animation::OpLeg::EdgeSweepSpec{.edge = &death_opleg_edge,
                                      .reverse = false,
                                      .sweep_frames = opaque(1),
                                      .settle_frames = opaque(0)});
}

/**
 * @brief Death case: a hankin leg without a palette handoff must trap.
 * @details Pins the hankin constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_hankin_incomplete_handoff() {
  death_construct_opleg(
      Animation::OpLeg::HankinSweepSpec{.theta_start = opaque(0.1f),
                                        .theta_end = opaque(0.5f),
                                        .sweep_frames = opaque(1)});
}

/**
 * @brief Death case: a hankin leg sweeping to a smaller angle must trap.
 * @details OpLeg surface — the leg sweeps the slerp fraction outward from the
 *          collapsed corner, which is monotone only while the arrival angle is
 *          the larger of the two.
 */
inline void case_opleg_hankin_backward_theta() {
  death_construct_opleg<8192>(
      Animation::OpLeg::HankinSweepSpec{.theta_start = opaque(0.9f),
                                        .theta_end = opaque(0.3f),
                                        .sweep_frames = opaque(1)},
      death_opleg_handoff());
}

/**
 * @brief Death case: a relax leg with neither a bake nor iterations must trap.
 * @details OpLeg surface — the leg needs a relaxed endpoint to slerp to, which
 *          is either the shipped bake or the result of live iterations; with
 *          neither it would slerp the seed onto itself for its whole duration.
 */
inline void case_opleg_relax_no_iterations() {
  death_construct_opleg(Animation::OpLeg::RelaxSpec{
      .iterations = opaque(0), .bake = nullptr, .sweep_frames = opaque(1)});
}

/**
 * @brief Death case: a relax leg without a palette handoff must trap.
 * @details Pins the relax constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_relax_incomplete_handoff() {
  death_construct_opleg(Animation::OpLeg::RelaxSpec{
      .iterations = opaque(1), .bake = nullptr, .sweep_frames = opaque(1)});
}

/**
 * @brief Death case: a medial leg without a palette handoff must trap.
 * @details Pins the medial constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_medial_incomplete_handoff() {
  death_construct_opleg(
      Animation::OpLeg::MedialSpec{.sweep_frames = opaque(1)});
}

/**
 * @brief Death case: a reconcile leg without endpoints must trap.
 * @details OpLeg surface — the leg slerps every seed vertex to an authored
 *          position, so an absent endpoint array is the whole leg's target.
 */
inline void case_opleg_reconcile_no_endpoints() {
  death_construct_opleg(Animation::OpLeg::ReconcileSpec{
      .to_positions = nullptr, .sweep_frames = opaque(1)});
}

/** @brief Death case: reconcile endpoint and seed counts must agree. */
inline void case_opleg_reconcile_endpoint_count() {
  const math::Vector ENDPOINT{};
  death_construct_opleg(
      Animation::OpLeg::ReconcileSpec{.to_positions = &ENDPOINT,
                                      .to_count = opaque<size_t>(1),
                                      .sweep_frames = opaque(1)});
}

/**
 * @brief Death case: a reconcile leg without a palette handoff must trap.
 * @details Pins the reconcile constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_reconcile_incomplete_handoff() {
  static const math::Vector ENDPOINTS[1] = {math::Vector(0, 0, 1)};
  death_construct_opleg(Animation::OpLeg::ReconcileSpec{
      .to_positions = ENDPOINTS, .sweep_frames = opaque(1)});
}

/**
 * @brief Death case: a gated-swap leg with no gate window must trap.
 * @details OpLeg surface — the leg runs 2*gate_frames + 1 frames around the
 *          swap, so a zero gate leaves the swap frame with no approach or
 *          departure to blend across.
 */
inline void case_opleg_gated_swap_zero_gate_frames() {
  death_construct_opleg(Animation::OpLeg::GatedSwapSpec{
      .op = Animation::OpLeg::SwapOp::KIS, .gate_frames = opaque(0)});
}

/**
 * @brief Death case: a gated-swap leg without a palette handoff must trap.
 * @details Pins the gated-swap constructor routing through the shared
 * init_transients palette-handoff guard.
 */
inline void case_opleg_gated_swap_incomplete_handoff() {
  death_construct_opleg(Animation::OpLeg::GatedSwapSpec{
      .op = Animation::OpLeg::SwapOp::KIS, .gate_frames = opaque(1)});
}

/**
 * @brief Death case: a crossfade lookup past the leg's face table must trap.
 */
inline void case_opleg_shading_blend_face_out_of_range() {
  static BakedPalette ramps[1];
  static const uint8_t FACE_RAMP[1] = {0};
  const Animation::OpLeg::Shading shading{
      .ramps = ramps, .face_ramp = FACE_RAMP, .faces = opaque<size_t>(1)};
  FacePaletteShader shader;
  shading.blend_face(shader, opaque<size_t>(3));
}

/**
 * @brief Death case: a shading lookup past the leg's face table must trap.
 */
inline void case_opleg_shading_face_out_of_range() {
  static BakedPalette ramps[1];
  static const uint8_t FACE_RAMP[1] = {0};
  const Animation::OpLeg::Shading shading{
      .ramps = ramps, .face_ramp = FACE_RAMP, .faces = opaque<size_t>(1)};
  (void)&shading.ramp_for(opaque<size_t>(3));
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
 * @brief Death case: a live-source Driver built with a null speed pointer must trap.
 * @details The guard traps before the constructor body reads *speed_src.
 */
inline void case_driver_null_speed_src() {
  static float mutant = 0.0f;
  Animation::Driver d(mutant, opaque<const float *>(nullptr),
                      1.0f); // -> HS_CHECK
  (void)d;
  (void)mutant;
}

/**
 * @brief Death case: Path::append_segment with zero samples must trap.
 * @details Animation surface — a zero sample count divides by zero in the
 *          t / samples term (easing(0/0) = NaN).
 */
inline void case_path_append_zero_samples() {
  Path<32> path;
  path.append_segment([](float s) { return math::Vector(s, 0.0f, 0.0f); }, 1.0f,
                      opaque(0),
                      [](float t) { return t; }); // samples < 1 -> HS_CHECK
}

/**
 * @brief Death case: a RandomTimer with min > max must trap.
 * @details The constructor traps an inverted or negative range.
 */
inline void case_random_timer_inverted_range() {
  Animation::RandomTimer timer({.min = opaque(5), .max = opaque(2)},
                               [](Canvas &) {}); // min > max -> HS_CHECK
  (void)timer;
}

/**
 * @brief Death case: clear()ing a pinned event must trap.
 * @details Animation surface — clear() would otherwise free an event whose
 *          animation pointer the caller still holds.
 */
inline void case_timeline_clear_pinned() {
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 1, math::ease_linear));
  global_timeline_events[0].pinned = opaque(true);
  tl.clear(); // HS_CHECK(!pinned) -> trap
}

/**
 * @brief Death case: clear()ing from a completion callback must trap.
 * @details A clear() inside the callback would free the executing callable.
 */
inline void case_timeline_clear_during_step() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline tl;
  float v = 0.0f;
  tl.add(0, Animation::Transition(v, 1.0f, 1, math::ease_linear).then([&tl]() {
    tl.clear();
  }));
  tl.step(canvas); // t=1: completes -> callback -> clear() while stepping
}

/** @brief Rejects owner cancellation from a completion callback. */
inline void case_timeline_cancel_owner_during_step() {
  static hs_test::StubEffect fx(8, 8);
  static Canvas canvas(fx);
  Timeline timeline;
  float value = 0.0f;
  timeline.add(
      0, Animation::Transition(value, 1.0f, 1, math::ease_linear)
             .then([&timeline, &value]() { timeline.cancel_owner(&value); }));
  timeline.step(canvas);
}

/** @brief Death case: finite parameter animations reject the -1 sentinel. */
inline void case_finite_param_perpetual_duration() {
  float value = 0.0f;
  Animation::Transition transition(value, 1.0f, opaque(-1), math::ease_linear);
  (void)transition.done();
}

/** @brief Death case: a Transition target must be finite. */
inline void case_transition_nonfinite_target() {
  float value = 0.0f;
  Animation::Transition transition(
      value, opaque(std::numeric_limits<float>::quiet_NaN()), 1,
      math::ease_linear);
  (void)transition.done();
}

/** @brief Death case: clear hooks must not mutate timeline event storage. */
inline void case_timeline_clear_hook_adds_event() {
  Timeline tl;
  tl.add_clear_hook(&tl, add_event_from_clear_hook);
  tl.clear();
}

/**
 * @brief Death case: scheduling a segue sprite with no free timeline slot must
 *        trap.
 * @details Animation surface — a dropped sprite add would leave the sphere
 *          dark for a whole transition.
 */
inline void case_segue_sprite_no_slot() {
  Timeline tl;
  float sink = 0.0f;
  while (Timeline::remaining() > 0)
    tl.add(0, Animation::Transition(sink, 1.0f, 1000, math::ease_linear));
  Segue::schedule_faded_sprite(tl, [](Canvas &, float) {}, 4, 1);
}

/** @brief Death case: a segue must target the already-flipped front slot. */
inline void case_mesh_carousel_unflipped_slot() {
  Timeline tl;
  MeshCarousel<> carousel;
  carousel.schedule_segue(tl, 1, [](Canvas &, float) {}, 4, 1);
}

/** @brief Death case: fill_face_phases on an unclassified mesh must trap. */
inline void case_mesh_carousel_unclassified_face_phases() {
  static uint8_t buf[1 << 14];
  Arena arena(buf, sizeof(buf));
  PolyMesh poly;
  build_solid<Solids::Octahedron>(poly, arena);
  MeshState mesh;
  MeshOps::compile(poly, mesh, persistent_arena, scratch_arena_a);
  ArenaVector<float> phases;
  phases.bind(persistent_arena, mesh.num_faces());
  MeshCarousel<Segue::Breakdown> carousel;
  carousel.fill_face_phases(mesh, mesh, opaque(0.5f), phases);
}

/**
 * @brief Death case: a second simultaneously-live Timeline must trap.
 * @details Animation surface — every Timeline shares the single global event
 *          array, so a second live instance would stomp the first's events.
 */
inline void case_timeline_double_construct() {
  Timeline a;
  Timeline b; // second live ctor -> HS_CHECK(!global_timeline_live) -> trap
  (void)global_timeline_num_events;
}

inline void case_random_walk_nonfinite_options() {
  math::Orientation<> orientation;
  FastNoiseLite noise;
  Animation::RandomWalkOptions options;
  options.drift = opaque(std::numeric_limits<float>::quiet_NaN());
  Animation::RandomWalk<32> walk(orientation, math::Vector(0, 0, 1), noise,
                                 options);
}

/** @brief A required event traps on a full timeline. */
inline void case_timeline_required_no_event_slot() {
  Timeline tl;
  float sink = 0.0f;
  for (int i = 0; i < Timeline::MAX_EVENTS; ++i)
    tl.add(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
  tl.add_required(0, Animation::Transition(sink, 1.0f, 10, math::ease_linear));
}

/** @brief Traps on a nonfinite captured ripple peak. */
inline void case_ripple_nonfinite_peak() {
  Animation::RippleParams params;
  params.amplitude = opaque(std::numeric_limits<float>::quiet_NaN());
  Animation::Ripple ripple(params, math::Y_AXIS, 0.2f, 4);
}

/** @brief Traps on an infinite ball-drop radius. */
inline void case_ball_drop_nonfinite_radius() {
  Animation::BumpParams params;
  params.radius = opaque(std::numeric_limits<float>::infinity());
  math::Orientation<4> orientation;
  Animation::BallDrop drop(params, orientation, math::Y_AXIS, 0.0f, 4);
}

/** @brief Rejects non-finite ball-drop azimuth. */
inline void case_ball_drop_nonfinite_azimuth() {
  Animation::BumpParams params;
  params.radius = 0.5f;
  math::Orientation<> orientation;
  Animation::BallDrop<> drop(params, orientation, math::Y_AXIS,
                             opaque(std::numeric_limits<float>::quiet_NaN()),
                             10);
}

/** @brief A cached bump offset must agree with the sample's cap distance. */
inline void case_bump_offset_outside_cap_distance() {
  Animation::BumpParams params;
  params.center = math::Y_AXIS;
  params.axis = math::Y_AXIS;
  params.radius = 0.5f;
  params.amplitude = 1.0f;
  params.envelope = 1.0f;
  params.sync();
  (void)bump_field_with_y(math::Y_AXIS, params, opaque(0.25f));
}

inline void case_timeline_add_into_live_slot() {
  Timeline timeline;
  float value = 0;
  timeline.add(0, Animation::Transition(value, 1.0f, 10, math::ease_linear));
  global_timeline_num_events = 0;
  timeline.add(0, Animation::Transition(value, 1.0f, 10, math::ease_linear));
}

/** @brief Rejects a random-timer maximum whose inclusive bound overflows. */
inline void case_random_timer_max_int() {
  Animation::RandomTimer timer(
      {.min = 0, .max = opaque(std::numeric_limits<int>::max())},
      [](Canvas &) {});
  (void)timer;
}

/** @brief Plateau boundaries reject a full timeline. */
inline void case_plateau_boundary_full_timeline() {
  Timeline timeline;
  Segue::PlateauHold hold;
  hold.arm_plateau_end(timeline, 1);
  death_admitted();
  while (Timeline::remaining() > 0)
    (void)timeline.add(0, Animation::PeriodicTimer(1, [](Canvas &) {}));
  hold.arm_plateau_end(timeline, 1);
}
