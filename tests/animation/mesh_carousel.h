/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// MeshCarousel segues
// ============================================================================

/**
 * @brief Verifies Segue::Crossfade schedules one fading sprite through the
 * carousel seam and returns a next-transition delay that overlaps consecutive
 * sprites by exactly the fade window.
 */
inline void test_crossfade_segue_schedules_overlapping_sprite() {
  Timeline tl;
  std::vector<float> ops;
  const int dur = 10, window = 3;
  MeshCarousel<Segue::Crossfade> carousel;
  int next_delay = carousel.schedule_segue(
      tl, carousel.front_index(), [&](Canvas &, float o) { ops.push_back(o); },
      dur, window);
  HS_EXPECT_EQ(next_delay, dur - window);
  HS_EXPECT_EQ(tl.event_count(), 1);

  for (int i = 0; i < dur; ++i)
    tl.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_SIZE_OR_RETURN(ops, dur);
  // Fade-in ramp, full-opacity plateau, transparent on the final frame.
  HS_EXPECT_LT(ops[0], 1.0f);
  HS_EXPECT_NEAR(ops[4], 1.0f, 1e-3f);
  HS_EXPECT_NEAR(ops[dur - 1], 0.0f, 1e-3f);
  HS_EXPECT_EQ(tl.event_count(), 0);
}

/**
 * @brief Verifies Segue::Crossfade clamps the fade window to half the duration.
 */
inline void test_crossfade_segue_clamps_fade_to_half_duration() {
  Timeline tl;
  const int dur = 10, window = 9; // > dur/2, clamps to 5
  MeshCarousel<Segue::Crossfade> carousel;
  int next_delay = carousel.schedule_segue(
      tl, carousel.front_index(), [](Canvas &, float) {}, dur, window);
  HS_EXPECT_EQ(next_delay, dur - dur / 2);
}

/**
 * @brief Verifies Segue::Crossfade's overlap parameter: zero overlap returns
 * the full duration, an oversized overlap clamps to the fade window, and an
 * in-range overlap is honored exactly.
 */
inline void test_crossfade_segue_overlap_is_configurable() {
  Timeline tl;
  const int dur = 10, window = 3;
  MeshCarousel<Segue::Crossfade> carousel;
  auto noop = [](Canvas &, float) {};
  carousel.segue().overlap = 0;
  HS_EXPECT_EQ(
      carousel.schedule_segue(tl, carousel.front_index(), noop, dur, window),
      dur);
  carousel.segue().overlap = 100;
  HS_EXPECT_EQ(
      carousel.schedule_segue(tl, carousel.front_index(), noop, dur, window),
      dur - window);
  carousel.segue().overlap = 2;
  HS_EXPECT_EQ(
      carousel.schedule_segue(tl, carousel.front_index(), noop, dur, window),
      dur - 2);
}

/**
 * @brief Verifies the default (Base) scheduling is sequential: the returned
 * delay equals the full duration.
 */
inline void test_sequential_segue_never_overlaps_sprites() {
  Timeline tl;
  std::vector<float> ops;
  const int dur = 10, window = 3;
  MeshCarousel<Segue::SpinFlip> carousel;
  int next_delay = carousel.schedule_segue(
      tl, carousel.front_index(), [&](Canvas &, float p) { ops.push_back(p); },
      dur, window);
  HS_EXPECT_EQ(next_delay, dur);
  HS_EXPECT_EQ(tl.event_count(), 1);

  for (int i = 0; i < dur; ++i)
    tl.step(fake_canvas()); // observed at t = 1..10

  HS_EXPECT_SIZE_OR_RETURN(ops, dur);
  // Phase ramps up, plateaus at 1, and returns to 0 on the final frame.
  HS_EXPECT_LT(ops[0], 1.0f);
  HS_EXPECT_NEAR(ops[4], 1.0f, 1e-3f);
  HS_EXPECT_NEAR(ops[dur - 1], 0.0f, 1e-3f);
}

/**
 * @brief Verifies Segue::Dissolve's complementary masks partition the key
 *        domain: every key is owned by exactly one of the two meshes, and the
 *        incoming share tracks the phase.
 */
inline void test_dissolve_segue_masks_partition_keys() {
  constexpr int KA = 64, KB = 32;
  Segue::Dissolve dissolve;
  constexpr float PHASES[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  for (int pi = 0; pi < 5; ++pi) {
    HS_CONTEXT("phase", pi);
    const float phase = PHASES[pi];
    auto masks = dissolve.mask_pair(phase, 7u);
    int owned_in = 0;
    for (int kb = 0; kb < KB; ++kb)
      for (int ka = 0; ka < KA; ++ka) {
        HS_CONTEXT("key", ka, kb);
        bool a = masks.incoming.owns(ka, kb), b = masks.outgoing.owns(ka, kb);
        HS_EXPECT_NE(a, b);
        owned_in += a ? 1 : 0;
      }
    // Hash spread is not perfectly uniform, so allow a few percent of slack.
    float share = static_cast<float>(owned_in) / (KA * KB);
    HS_EXPECT_NEAR(share, phase, 0.05f);
  }
}

/**
 * @brief Verifies Segue::Dissolve re-rolls its pattern every frame and
 *        re-seeds per transition, so successive frames dither rather than
 *        freezing one stochastic partition on screen.
 */
inline void test_dissolve_segue_reseeds_per_frame_and_transition() {
  Segue::Dissolve dissolve;
  auto f0 = dissolve.mask_pair(0.5f, 0u);
  auto f1 = dissolve.mask_pair(0.5f, 1u);
  HS_EXPECT_NE(f0.incoming.salt, f1.incoming.salt);
  HS_EXPECT_EQ(f0.incoming.threshold, f1.incoming.threshold);
  // Both halves of a pair share the salt and threshold that make them partition.
  HS_EXPECT_EQ(f0.incoming.salt, f0.outgoing.salt);
  HS_EXPECT_EQ(f0.incoming.threshold, f0.outgoing.threshold);
  HS_EXPECT_TRUE(f0.incoming.invert != f0.outgoing.invert);

  uint32_t before = dissolve.mask_pair(0.5f, 0u).incoming.salt;
  dissolve.retarget(math::Vector(0, 1, 0));
  HS_EXPECT_NE(dissolve.mask_pair(0.5f, 0u).incoming.salt, before);
}

/**
 * @brief Verifies Segue::Dissolve schedules the two meshes co-resident for the
 *        whole fade window: its negative overlap selects the full window, so
 *        every frame of the transition has both halves of the mask pair
 *        drawing.
 */
inline void test_dissolve_segue_overlaps_the_full_fade_window() {
  Timeline tl;
  const int dur = 10, window = 3; // fade = min(window, dur/2) = window
  MeshCarousel<Segue::Dissolve> carousel;
  bool drew_out = false, drew_in = false;
  int co_resident = 0;
  int outgoing_frames = 0;
  int next_delay = carousel.schedule_segue(
      tl, carousel.front_index(), [&](Canvas &, float) { drew_out = true; },
      dur, window);
  HS_EXPECT_EQ(next_delay, dur - window);

  for (int i = 0; i < next_delay; ++i) {
    drew_out = drew_in = false;
    tl.step(fake_canvas());
    outgoing_frames += drew_out;
  }
  HS_EXPECT_EQ(outgoing_frames, next_delay);

  carousel.schedule_segue(
      tl, carousel.front_index(), [&](Canvas &, float) { drew_in = true; }, dur,
      window);
  for (int i = 0; i < dur; ++i) {
    drew_out = drew_in = false;
    tl.step(fake_canvas());
    co_resident += drew_out && drew_in;
  }
  HS_EXPECT_EQ(co_resident, window);
}

/**
 * @brief Verifies Breakdown's degenerate-input guards: a single-class mesh
 * fades as one unit instead of dividing by a zero rank span, an out-of-range
 * class id resolves through rank[0] rather than reading past the array, and
 * face_phase saturates outside a class's window.
 */
inline void test_breakdown_guards_degenerate_class_inputs() {
  const auto saved_rng = hs::random();
  Segue::Breakdown bd; // default: one class, no reorder() yet
  const math::Vector any(0.0f, 1.0f, 0.0f);
  HS_EXPECT_EQ(bd.num_classes, 1);
  HS_EXPECT_NEAR(bd.face_offset(any, 0, 0), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(bd.face_offset(any, 0, 5), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(bd.face_phase(0.0f, 0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(bd.face_phase(1.0f, 0.0f), 1.0f, 1e-6f);

  hs::random().seed(11u);
  const std::vector<uint8_t> face_classes{0, 1, 2};
  bd.reorder(face_classes);
  HS_EXPECT_EQ(bd.num_classes, 3);
  const float class0 = bd.face_offset(any, 0, 0);
  HS_EXPECT_NEAR(bd.face_offset(any, 0, -1), class0, 1e-6f);
  HS_EXPECT_NEAR(bd.face_offset(any, 0, bd.num_classes), class0, 1e-6f);
  for (int c = 0; c < bd.num_classes; ++c) {
    const float o = bd.face_offset(any, 0, c);
    HS_EXPECT_NEAR(bd.face_phase(0.0f, o), 0.0f, 1e-6f);
    HS_EXPECT_NEAR(bd.face_phase(1.0f, o), 1.0f, 1e-6f);
  }
  hs::random() = saved_rng;
}

/**
 * @brief Verifies the Base shading hooks are identities: full opacity and
 * coverage, unmodified edge distance and color, never culled.
 */
inline void test_segue_base_hooks_are_identity() {
  Segue::Base base;
  HS_EXPECT_NEAR(base.opacity(0.3f), 1.0f, 1e-6f);
  float t = 0.42f;
  HS_EXPECT_NEAR(base.fill(t, 0.3f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(t, 0.42f, 1e-6f);
  Color4 c(Pixel(1000, 2000, 3000), 0.5f);
  Color4 g = base.grade(c, 0.3f);
  HS_EXPECT_EQ(g.color.r, c.color.r);
  HS_EXPECT_EQ(g.color.g, c.color.g);
  HS_EXPECT_EQ(g.color.b, c.color.b);
  HS_EXPECT_NEAR(g.alpha, c.alpha, 1e-6f);
  HS_EXPECT_TRUE(base.visible(0.5f));
  HS_EXPECT_TRUE(base.visible(0.0f));
}

/**
 * @brief Peak shading weight a policy can put on screen at one global phase.
 * @details Per-face policies resolve opacity from their own face-local phase,
 * so they are sampled across the offset range; fragment policies are sampled
 * across the edge-distance range through fill.
 */
template <typename S>
inline float segue_peak_weight(const S &seg, float phase) {
  float peak = 0.0f;
  for (int i = 0; i <= 20; ++i) {
    if constexpr (Segue::PerFace<S>) {
      const float offset = static_cast<float>(i) / 20.0f;
      const float local = seg.face_phase(phase, offset, seg.face_fade_frac(i));
      peak = std::max(peak, seg.opacity(local));
    } else {
      float t = static_cast<float>(i) / 20.0f;
      peak = std::max(peak, seg.fill(t, phase) * seg.opacity(phase));
    }
  }
  return peak;
}

/**
 * @brief Verifies every segue policy's visible() gate agrees with what that
 * policy actually shades, in both directions.
 * @details visible() is a whole-draw cull: it must pass every phase the
 * policy shades and cull phases where it shades nothing.
 */
inline void test_segue_visible_gate_culls_only_dark_phases() {
  auto check = [](const auto &seg) {
    for (int i = 0; i <= 1000; ++i) {
      const float phase = static_cast<float>(i) / 1000.0f;
      if (seg.visible(phase))
        HS_EXPECT_GT(segue_peak_weight(seg, phase), 0.0f);
      else
        HS_EXPECT_LT(segue_peak_weight(seg, phase), 0.02f);
    }
  };
  Segue::AllPolicies::for_each(check);

  // The gate still culls for the policies that do fade to black.
  HS_EXPECT_FALSE(Segue::Crossfade().visible(0.0f));
  HS_EXPECT_TRUE(Segue::Crossfade().visible(0.5f));
  HS_EXPECT_FALSE(Segue::TerminatorSweep().visible(0.0f));
  HS_EXPECT_FALSE(Segue::Shockwave().visible(0.0f));
  HS_EXPECT_TRUE(Segue::Shockwave().visible(0.5f));
  // Breakdown holds its whole BLACK_DWELL slice at face phase 0.
  HS_EXPECT_FALSE(Segue::Breakdown().visible(Segue::Breakdown::BLACK_DWELL));
  HS_EXPECT_TRUE(Segue::Breakdown().visible(0.5f));
  // ...and never culls the policies that keep shading at phase 0.
  HS_EXPECT_TRUE(Segue::GoldConvergence().visible(0.0f));
  HS_EXPECT_TRUE(Segue::SpinFlip().visible(0.0f));
  HS_EXPECT_TRUE(Segue::IrisBloom().visible(0.0f));
  HS_EXPECT_TRUE(Segue::Lace().visible(0.0f));
  HS_EXPECT_TRUE(Segue::Dissolve().visible(0.0f));
}

/**
 * @brief Verifies IrisBloom's fill contracts faces toward their centers: at
 * full phase everything survives, at mid phase only fragments deeper than the
 * inset do, and the surviving core renormalizes to the full gradient.
 */
inline void test_iris_bloom_fill_contracts_to_face_centers() {
  Segue::IrisBloom iris;
  // Full phase: everything covered, t unchanged.
  float t = 0.3f;
  HS_EXPECT_NEAR(iris.fill(t, 1.0f), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(t, 0.3f, 1e-3f);
  // Mid phase: shallow fragments culled...
  t = 0.3f;
  HS_EXPECT_NEAR(iris.fill(t, 0.5f), 0.0f, 1e-6f);
  // ...deep fragments survive with t renormalized over the shrunken core.
  t = 0.9f;
  HS_EXPECT_NEAR(iris.fill(t, 0.5f), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(t, 0.8f, 1e-3f);
}

/**
 * @brief Verifies Lace's fill is the inverse mask: only fragments within the
 * phase-driven band of an edge survive, renormalized across the band.
 */
inline void test_lace_fill_keeps_edge_band() {
  Segue::Lace lace;
  // Full phase: everything covered.
  float t = 0.9f;
  HS_EXPECT_NEAR(lace.fill(t, 1.0f), 1.0f, 1e-3f);
  // Narrow band: deep fragments culled, near-edge fragments survive.
  t = 0.5f;
  HS_EXPECT_NEAR(lace.fill(t, 0.3f), 0.0f, 1e-6f);
  t = 0.15f;
  HS_EXPECT_NEAR(lace.fill(t, 0.3f), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(t, 0.5f, 1e-3f);
}

/**
 * @brief Verifies the shared sweep front: full everywhere at phase 1, dark
 * everywhere at phase 0, monotone in phase, and higher offsets extinguish
 * earlier.
 */
inline void test_sweep_phase_front_ordering() {
  const float band = 0.25f;
  for (float o : {0.0f, 0.5f, 1.0f}) {
    HS_EXPECT_NEAR(Segue::sweep_phase(1.0f, o, band), 1.0f, 1e-3f);
    HS_EXPECT_NEAR(Segue::sweep_phase(0.0f, o, band), 0.0f, 1e-3f);
  }
  // Monotone in phase at fixed offset, sampled inside the front's ramp
  // (phase 0.14..0.39 for this offset/band under the sqrt ease).
  HS_EXPECT_GT(Segue::sweep_phase(0.3f, 0.5f, band),
               Segue::sweep_phase(0.2f, 0.5f, band));
  // At a fixed phase, a higher offset is further extinguished.
  HS_EXPECT_GT(Segue::sweep_phase(0.5f, 0.2f, band),
               Segue::sweep_phase(0.5f, 0.8f, band));
}

/** @brief Pins sweep coordinate selection and transformed topology slots. */
inline void test_meshcarousel_face_phases_use_sweep_frame_and_slots() {
  hs_test::reset_globals();
  static uint8_t polybuf[1 << 14];
  Arena polyarena(polybuf, sizeof(polybuf));
  PolyMesh poly;
  build_solid<Solids::Octahedron>(poly, polyarena);
  MeshState base, transformed;
  MeshOps::compile(poly, base, persistent_arena, scratch_arena_a);
  MeshOps::compile(poly, transformed, persistent_arena, scratch_arena_a);
  for (auto &vertex : transformed.vertices) {
    vertex.y = -vertex.y;
    vertex.z = -vertex.z;
  }
  base.topology.bind(persistent_arena, base.num_faces());
  transformed.topology.bind(persistent_arena, transformed.num_faces());
  base.topology.clear();
  transformed.topology.clear();
  for (size_t face = 0; face < base.num_faces(); ++face) {
    base.topology.push_back(0);
    transformed.topology.push_back(
        static_cast<uint16_t>(MeshPaletteBank::N + 1));
  }
  ArenaVector<float> phases;
  phases.bind(persistent_arena, base.num_faces());
  auto check = [&]<typename Policy>(const MeshState &expected,
                                    const MeshState &other) {
    MeshCarousel<Policy> carousel;
    if constexpr (requires { carousel.segue().retarget(math::Y_AXIS); })
      carousel.segue().retarget(math::Y_AXIS);
    else {
      carousel.segue().num_classes = 2;
      for (int i = 0; i < 16; ++i)
        carousel.segue().rank[i] = static_cast<uint8_t>(i);
    }
    carousel.fill_face_phases(base, transformed, 0.5f, phases);
    HS_EXPECT_EQ(phases.size(), expected.num_faces());
    bool distinguishes_frame = false;
    for (size_t face = 0; face < expected.num_faces(); ++face) {
      auto center = [&](const MeshState &mesh) {
        return math::normalized_or(
            MeshOps::face_vertex_sum(mesh.vertices.data(),
                                     mesh.get_faces_data(),
                                     mesh.get_face_offsets_data()[face],
                                     mesh.get_face_counts_data()[face]),
            math::UP);
      };
      const int CLS = MeshPaletteBank::slot_of(transformed.topology[face]);
      const auto &policy = carousel.segue();
      const float WANT = policy.face_phase(
          0.5f, policy.face_offset(center(expected), face, CLS),
          policy.face_fade_frac(face));
      const float OTHER =
          policy.face_phase(0.5f, policy.face_offset(center(other), face, CLS),
                            policy.face_fade_frac(face));
      HS_EXPECT_NEAR(phases[face], WANT, 1e-6f);
      distinguishes_frame |= fabsf(WANT - OTHER) > 0.1f;
    }
    if constexpr (!std::is_same_v<Policy, Segue::Breakdown>)
      HS_EXPECT_TRUE(distinguishes_frame);
    else {
      const auto &policy = carousel.segue();
      HS_EXPECT_NEAR(phases[0],
                     policy.face_phase(0.5f, policy.face_offset(
                                                 math::UP, 0,
                                                 MeshPaletteBank::slot_of(
                                                     transformed.topology[0]))),
                     1e-6f);
      HS_EXPECT_TRUE(phases[0] != policy.face_phase(0.5f, policy.face_offset(
                                                              math::UP, 0, 0)));
    }
  };
  check.template operator()<Segue::TerminatorSweep>(base, transformed);
  check.template operator()<Segue::Shockwave>(transformed, base);
  check.template operator()<Segue::Breakdown>(transformed, base);
}

/**
 * @brief Verifies TerminatorSweep orders faces along its axis: the axis pole
 * extinguishes first (offset 1), the antipode last (offset 0), and the
 * per-face fade alpha is the perceptual square of the face-local phase with
 * exact endpoints.
 */
inline void test_terminator_sweep_orders_by_axis() {
  Segue::TerminatorSweep term;
  math::Vector raw_axis(1.0f, 2.0f, -0.5f);
  math::Vector axis = raw_axis.normalized();
  term.retarget(raw_axis);
  HS_EXPECT_NEAR(term.axis.length(), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(term.face_offset(axis, 0, 0), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(term.face_offset(-axis, 0, 0), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(
      term.face_offset(math::cross(axis, math::X_AXIS).normalized(), 0, 0),
      0.5f, 1e-2f);
  HS_EXPECT_NEAR(term.opacity(0.4f), 0.16f, 1e-6f);
  HS_EXPECT_NEAR(term.opacity(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(term.opacity(1.0f), 1.0f, 1e-6f);
}

/**
 * @brief Verifies TerminatorSweep's per-face fade is time-based: the fade
 * length divides the window schedule() recorded, so a face ramps over its fade
 * length once the front reaches it, with exact 0/1 window endpoints, and a
 * window shorter than the fade length degrades to one whole-sphere fade.
 */
inline void test_terminator_sweep_fades_faces_over_fixed_frames() {
  Timeline tl;
  const int dur = 200, window = 32;
  MeshCarousel<Segue::TerminatorSweep> carousel;
  carousel.segue().fade_frames_min = 8.0f;
  carousel.segue().fade_frames_max = 8.0f;
  int next_delay = carousel.schedule_segue(
      tl, carousel.front_index(), [](Canvas &, float) {}, dur, window);
  HS_EXPECT_EQ(next_delay, dur); // sequential: one mesh per frame
  const Segue::TerminatorSweep &term = carousel.segue();
  const float f = 8.0f / window;
  HS_EXPECT_NEAR(term.face_fade_frac(0), f, 1e-6f);
  HS_EXPECT_NEAR(term.face_fade_frac(97), f, 1e-6f);
  for (float o : {0.0f, 0.5f, 1.0f}) {
    float touch = o * (1.0f - f); // phase at which the front reaches the face
    HS_EXPECT_NEAR(term.face_phase(touch, o, f), 0.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(touch + 0.5f * f, o, f), 0.5f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(touch + f, o, f), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(1.0f, o, f), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(0.0f, o, f), 0.0f, 1e-5f);
    // Independent of the phase algebra: a clamped monotone ramp over the whole
    // transition, so no sampled phase can back up or leave [0, 1].
    float prev = 0.0f;
    for (int i = 0; i <= 40; ++i) {
      const float ph = term.face_phase(i / 40.0f, o, f);
      HS_EXPECT_GE(ph, prev - 1e-5f);
      HS_EXPECT_GE(ph, 0.0f);
      HS_EXPECT_LE(ph, 1.0f);
      prev = ph;
    }
  }
  carousel.schedule_segue(
      tl, carousel.front_index(), [](Canvas &, float) {}, 8, 2);
  HS_EXPECT_NEAR(carousel.segue().face_fade_frac(0), 1.0f, 1e-6f);
}

/**
 * @brief Verifies a mid-transition Face Fade slider move lands on the next
 * frame: face_fade_frac reads the frame bounds live.
 */
inline void test_terminator_sweep_fade_sliders_apply_without_reschedule() {
  Timeline tl;
  MeshCarousel<Segue::TerminatorSweep> carousel;
  carousel.segue().retarget(math::Y_AXIS);
  carousel.segue().fade_frames_min = 4.0f;
  carousel.segue().fade_frames_max = 4.0f;
  carousel.schedule_segue(
      tl, carousel.front_index(), [](Canvas &, float) {}, 400, 64);
  HS_EXPECT_NEAR(carousel.segue().face_fade_frac(11), 4.0f / 64.0f, 1e-6f);
  // Slider drag with no reschedule.
  carousel.segue().fade_frames_min = 16.0f;
  carousel.segue().fade_frames_max = 16.0f;
  HS_EXPECT_NEAR(carousel.segue().face_fade_frac(11), 16.0f / 64.0f, 1e-6f);
}

/**
 * @brief Verifies TerminatorSweep draws each face's fade length randomly from
 * [fade_frames_min, fade_frames_max]: the fractions stay in range, differ
 * across faces, and preserve exact 0/1 window endpoints for every per-face fade
 * length.
 */
inline void test_terminator_sweep_per_face_fade_random_in_range() {
  const auto saved_rng = hs::random();
  hs::random().seed(1337);
  Timeline tl;
  const int dur = 400, window = 64;
  MeshCarousel<Segue::TerminatorSweep> carousel;
  carousel.segue().retarget(math::Y_AXIS); // rolls the per-face fade seed
  carousel.segue().fade_frames_min = 4.0f;
  carousel.segue().fade_frames_max = 16.0f;
  carousel.schedule_segue(
      tl, carousel.front_index(), [](Canvas &, float) {}, dur, window);
  const Segue::TerminatorSweep &term = carousel.segue();
  const float lo = 4.0f / window, hi = 16.0f / window;
  const float first = term.face_fade_frac(0);
  bool varied = false;
  for (int i = 0; i < 256; ++i) {
    float ff = term.face_fade_frac(i);
    HS_EXPECT_TRUE(ff >= lo - 1e-6f && ff <= hi + 1e-6f);
    HS_EXPECT_NEAR(term.face_phase(1.0f, 0.7f, ff), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(term.face_phase(0.0f, 0.7f, ff), 0.0f, 1e-5f);
    if (std::fabs(ff - first) > 1e-4f)
      varied = true;
  }
  HS_EXPECT_TRUE(varied);
  hs::random() = saved_rng;
}

/**
 * @brief Verifies Shockwave orders faces by angular distance from its origin:
 * nearest faces extinguish first, the antipode last.
 */
inline void test_shockwave_orders_by_distance_from_origin() {
  Segue::Shockwave wave;
  math::Vector origin = math::Vector(0.3f, -1.0f, 0.7f).normalized();
  wave.retarget(origin);
  HS_EXPECT_NEAR(wave.face_offset(origin, 0, 0), 1.0f, 1e-2f);
  HS_EXPECT_NEAR(wave.face_offset(-origin, 0, 0), 0.0f, 1e-2f);
  // Equidistant ring sits mid-order.
  HS_EXPECT_NEAR(
      wave.face_offset(math::cross(origin, math::X_AXIS).normalized(), 0, 0),
      0.5f, 2e-2f);
}

/** @brief The per-face hook set the mesh draw calls once face_offset resolves.
 */
template <typename SegueT>
concept PerFaceSegueDrawable =
    requires(const SegueT &s, const math::Vector &c) {
      s.face_offset(c, 0, 0);
      s.face_fade_frac(0);
      s.face_phase(0.5f, 0.5f, 0.1f);
    };

/** @brief A policy whose face_phase takes two arguments: the authoring slip
 * MeshCarousel's per-face contract assert rejects. */
struct TwoArgFacePhaseSegue : Segue::Base {
  float face_offset(const math::Vector &, int, int) const { return 0.0f; }
  float face_phase(float phase, float) const { return phase; }
};

/** @brief A policy whose face_offset drops the palette-class argument: the
 * per-face draw path's call site no longer resolves. */
struct DriftedFaceOffsetSegue : Segue::Base {
  float face_offset(const math::Vector &, int) const { return 0.0f; }
  float face_phase(float phase, float, float) const { return phase; }
};

/** @brief Policies whose optional hooks carry drifted signatures: each is named
 * but uncallable at the contract's argument list. */
struct DriftedWarpSegue : Segue::Base {
  math::Vector warp(const math::Vector &v, float, int) const { return v; }
};
struct DriftedRetargetSegue : Segue::Base {
  int retarget(const math::Vector &) { return 0; }
};
struct DriftedReorderSegue : Segue::Base {
  template <typename Classes> void reorder(const Classes &, int) {}
};
struct DriftedMaskPairSegue : Segue::Base {
  int mask_pair(float, uint32_t) const { return 0; }
};
struct DriftedLocalSweepSegue : Segue::Base {
  static constexpr int LOCAL_SWEEP = 1;
};

/** @brief A final policy: no name carrier can be merged into it, so every
 * Declares* probe answers false, its drifted warp included. */
struct FinalSegue final : Segue::Base {
  math::Vector warp(const math::Vector &v, float, int) const { return v; }
};

/** @brief A policy shadowing Base's visible() with a float: nonzero phases
 * convert to true. */
struct DriftedVisibleSegue : Segue::Base {
  float visible(float phase) const { return phase; }
};

/** @brief A policy taking fill()'s edge distance by value instead of a mutable
 * reference. */
struct DriftedFillSegue : Segue::Base {
  float fill(float t, float) const { return t; }
};

/** @brief A policy taking grade()'s Color4 by reference instead of by value. */
struct DriftedGradeSegue : Segue::Base {
  Color4 grade(Color4 &c, float) const { return c; }
};

/**
 * @brief Pins every per-face segue against PerFaceSegueDrawable, plus the
 * segue trait probes against conforming and drifted policies.
 * @details A per-face draw path never calls fill/grade, so shadowing either
 * alongside face_offset would drop it silently.
 */
inline void test_per_face_segues_satisfy_draw_contract() {
  static_assert(PerFaceSegueDrawable<Segue::TerminatorSweep>);
  static_assert(PerFaceSegueDrawable<Segue::Shockwave>);
  static_assert(PerFaceSegueDrawable<Segue::Breakdown>);
  static_assert(!Segue::SHADOWS_FRAGMENT_HOOKS<Segue::TerminatorSweep>);
  static_assert(!Segue::SHADOWS_FRAGMENT_HOOKS<Segue::Shockwave>);
  static_assert(!Segue::SHADOWS_FRAGMENT_HOOKS<Segue::Breakdown>);
  static_assert(Segue::SHADOWS_FRAGMENT_HOOKS<Segue::IrisBloom>);
  static_assert(Segue::SHADOWS_FRAGMENT_HOOKS<Segue::GoldConvergence>);
  static_assert(!Segue::PerFace<Segue::IrisBloom>);
  static_assert(!Segue::PerFace<Segue::GoldConvergence>);
  static_assert(!Segue::HasFaceOffset<Segue::IrisBloom>);
  static_assert(Segue::HasFaceOffset<TwoArgFacePhaseSegue>);
  static_assert(!Segue::PerFace<TwoArgFacePhaseSegue>);
  static_assert(Segue::NeedsClasses<Segue::Breakdown>);
  static_assert(!Segue::NeedsClasses<Segue::TerminatorSweep>);
  static_assert(Segue::Masked<Segue::Dissolve>);
  static_assert(!Segue::Masked<Segue::Crossfade>);
  static_assert(Segue::HasWarp<Segue::SpinFlip>);
  static_assert(!Segue::HasWarp<Segue::Crossfade>);
  static_assert(Segue::HasRetarget<Segue::SpinFlip>);
  static_assert(Segue::HasRetarget<Segue::Dissolve>);
  static_assert(!Segue::HasRetarget<Segue::Crossfade>);
  static_assert(Segue::DeclaresWarp<Segue::SpinFlip>);
  static_assert(!Segue::DeclaresWarp<Segue::Crossfade>);
  static_assert(Segue::DeclaresRetarget<Segue::TerminatorSweep>);
  static_assert(!Segue::DeclaresRetarget<Segue::Crossfade>);
  static_assert(Segue::DeclaresReorder<Segue::Breakdown>);
  static_assert(!Segue::DeclaresReorder<Segue::TerminatorSweep>);
  static_assert(Segue::DeclaresMaskPair<Segue::Dissolve>);
  static_assert(!Segue::DeclaresMaskPair<Segue::Crossfade>);
  static_assert(Segue::DeclaresFaceOffset<Segue::TerminatorSweep>);
  static_assert(!Segue::DeclaresFaceOffset<Segue::Crossfade>);
  static_assert(Segue::DeclaresFacePhase<Segue::Shockwave>);
  static_assert(!Segue::DeclaresFacePhase<Segue::Crossfade>);
  static_assert(Segue::DeclaresLocalSweep<Segue::TerminatorSweep> &&
                Segue::LocalSweeps<Segue::TerminatorSweep>);
  static_assert(!Segue::DeclaresLocalSweep<Segue::Crossfade>);
  // A drifted hook is seen by name and rejected by signature.
  static_assert(Segue::DeclaresWarp<DriftedWarpSegue> &&
                !Segue::HasWarp<DriftedWarpSegue>);
  static_assert(Segue::DeclaresRetarget<DriftedRetargetSegue> &&
                !Segue::HasRetarget<DriftedRetargetSegue>);
  static_assert(Segue::DeclaresReorder<DriftedReorderSegue> &&
                !Segue::NeedsClasses<DriftedReorderSegue>);
  static_assert(Segue::DeclaresMaskPair<DriftedMaskPairSegue> &&
                !Segue::Masked<DriftedMaskPairSegue>);
  static_assert(Segue::DeclaresFaceOffset<DriftedFaceOffsetSegue> &&
                !Segue::HasFaceOffset<DriftedFaceOffsetSegue>);
  static_assert(Segue::DeclaresFacePhase<TwoArgFacePhaseSegue> &&
                !Segue::HasFacePhase<TwoArgFacePhaseSegue>);
  static_assert(Segue::DeclaresLocalSweep<DriftedLocalSweepSegue> &&
                !Segue::LocalSweeps<DriftedLocalSweepSegue>);
  static_assert(!Segue::PolicyList<DriftedLocalSweepSegue>::LOCAL_SWEEPS_TYPED);
  static_assert(!Segue::DeclaresWarp<FinalSegue> &&
                !Segue::PolicyList<FinalSegue>::MERGEABLE);
  static_assert(Segue::AllPolicies::MERGEABLE);
  static_assert(Segue::AllPolicies::CONFORMING);
  static_assert(Segue::AllPolicies::LOCAL_SWEEPS_TYPED);
  static_assert(!Segue::HasPhaseHooks<DriftedVisibleSegue>);
  static_assert(!Segue::HasPhaseHooks<DriftedFillSegue>);
  static_assert(!Segue::HasPhaseHooks<DriftedGradeSegue>);

  HS_EXPECT_NEAR(Segue::Base().face_fade_frac(3), 1.0f, 1e-6f);
}

/**
 * @brief Verifies every segue policy hands schedule()'s pause gate to the
 * sprite it schedules.
 * @details A gated sprite holds its envelope: the phase reaching the draw
 * callback must not move while the flag is set and must climb again once it
 * clears.
 */
inline void test_segue_policies_forward_pause_gate() {
  struct Probe {
    float phase = -1.0f;
    int draws = 0;
  };
  auto holds_under_pause = [](auto policy) {
    Timeline tl;
    Probe probe;
    bool paused = false;
    policy.schedule(
        tl,
        [&probe](Canvas &, float phase) {
          probe.phase = phase;
          probe.draws++;
        },
        /*duration=*/60, /*window=*/20, &paused);
    tl.step(fake_canvas());
    tl.step(fake_canvas());
    const int drawn = probe.draws;
    const float held = probe.phase;
    HS_EXPECT_GT(drawn, 0);
    HS_EXPECT_LT(held, 1.0f); // mid fade-in phase

    paused = true;
    for (int i = 0; i < 5; ++i)
      tl.step(fake_canvas());
    HS_EXPECT_EQ(probe.draws, drawn + 5); // held frames still draw
    HS_EXPECT_NEAR(probe.phase, held, 1e-6f);

    paused = false;
    tl.step(fake_canvas());
    HS_EXPECT_GT(probe.phase, held);
  };
  Segue::AllPolicies::for_each(holds_under_pause);
}

/**
 * @brief Verifies Breakdown fades classes sequentially: reorder() yields a
 * permutation of the class ranks, offsets follow the ranks, and each class's
 * fade window is an abutting 1/n slice of the phase range — fully faded
 * before the next class starts.
 */
inline void test_breakdown_fades_classes_sequentially() {
  const auto saved_rng = hs::random();
  hs::random().seed(7u);
  Segue::Breakdown bd;
  constexpr int n = 5;
  // Per-face classes (dense [0, n), out of order with repeats); reorder derives
  // num_classes = max + 1 = n from them rather than taking a declared count.
  const std::vector<uint8_t> face_classes{2, 0, 4, 1, 3, 4, 0};
  bd.reorder(face_classes);
  HS_EXPECT_EQ(bd.num_classes, n);
  bool seen[n] = {};
  for (int c = 0; c < n; ++c) {
    HS_EXPECT_LT(static_cast<int>(bd.rank[c]), n);
    seen[bd.rank[c]] = true;
  }
  for (int r = 0; r < n; ++r)
    HS_EXPECT_TRUE(seen[r]); // a permutation: every rank assigned once

  math::Vector any(0.0f, 1.0f, 0.0f);
  for (int c = 0; c < n; ++c) {
    float o = bd.face_offset(any, 0, c);
    int r = bd.rank[c];
    HS_EXPECT_NEAR(o, static_cast<float>(n - 1 - r) / (n - 1), 1e-6f);
    // Class rank r fades linearly over one band of [BLACK_DWELL, 1].
    float band = (1.0f - Segue::Breakdown::BLACK_DWELL) / n;
    float floor_p = Segue::Breakdown::BLACK_DWELL + (n - 1 - r) * band;
    HS_EXPECT_NEAR(bd.face_phase(floor_p, o), 0.0f, 1e-5f);
    HS_EXPECT_NEAR(bd.face_phase(floor_p + band, o), 1.0f, 1e-5f);
    HS_EXPECT_NEAR(bd.face_phase(Segue::Breakdown::BLACK_DWELL, o), 0.0f,
                   1e-5f);
    HS_EXPECT_NEAR(bd.face_phase(0.0f, o), 0.0f, 1e-6f);
  }

  // Independent of the band algebra: every class is a clamped monotone ramp,
  // and at any phase a later-ranked class is never behind an earlier one.
  for (int step = 0; step <= 40; ++step) {
    const float t = step / 40.0f;
    // Rank n-1 owns the earliest window, so descending rank never fades later.
    float prev_phase = 1.0f;
    for (int r = n - 1; r >= 0; --r) {
      int cls = -1;
      for (int c = 0; c < n; ++c)
        if (bd.rank[c] == r)
          cls = c;
      HS_EXPECT_GE(cls, 0);
      if (cls < 0)
        continue;
      const float ph = bd.face_phase(t, bd.face_offset(any, 0, cls));
      HS_EXPECT_GE(ph, 0.0f);
      HS_EXPECT_LE(ph, 1.0f);
      HS_EXPECT_LE(ph, prev_phase + 1e-5f);
      prev_phase = ph;
    }
  }
  hs::random() = saved_rng;
}

/**
 * @brief Verifies SpinFlip's warp is rigid: pairwise angles are preserved at
 * every phase, the mid-phase warp really spins about the axis, and the plateau
 * is the identity.
 */
inline void test_spin_flip_warp_is_rigid() {
  Segue::SpinFlip spin;
  spin.retarget(math::Vector(0.5f, 0.5f, -0.7f));
  HS_EXPECT_NEAR(spin.axis.length(), 1.0f, 1e-6f);
  math::Vector a = math::Vector(1.0f, 0.2f, 0.1f).normalized();
  math::Vector b = math::Vector(-0.3f, 0.9f, 0.4f).normalized();
  for (int i = 0; i <= 10; ++i) {
    const float phase = static_cast<float>(i) / 10.0f;
    math::Vector wa = spin.warp(a, phase), wb = spin.warp(b, phase);
    HS_EXPECT_NEAR(math::dot(wa, wb), math::dot(a, b), 1e-3f);
    HS_EXPECT_NEAR(wa.length(), 1.0f, 1e-3f);
    HS_EXPECT_NEAR(wb.length(), 1.0f, 1e-3f);
  }

  // Winding is (1 - phase)^2 * REVS revolutions; this phase makes it a quarter
  // turn, so an axis-perpendicular input lands perpendicular to where it began.
  float quarter = 1.0f - std::sqrt(0.25f / Segue::SpinFlip::REVS);
  math::Vector perp = math::cross(spin.axis, a).normalized();
  math::Vector wperp = spin.warp(perp, quarter);
  HS_EXPECT_NEAR(math::dot(wperp, perp), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(math::dot(wperp, spin.axis), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(math::dot(wperp, math::cross(spin.axis, perp)), 1.0f, 1e-3f);

  HS_EXPECT_NEAR((spin.warp(a, 1.0f) - a).length(), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(spin.opacity(0.0f), 1.0f,
                 1e-6f); // never fades: blur hides the swap
}

/**
 * @brief Verifies GoldConvergence grades toward its gold at the swap and is
 * the identity on the plateau, with the mild opacity dip.
 */
inline void test_gold_convergence_grades_to_gold() {
  Segue::GoldConvergence gc;
  Color4 c(Pixel(1000, 2000, 3000), 0.8f);
  Color4 plateau = gc.grade(c, 1.0f);
  HS_EXPECT_EQ(plateau.color.r, c.color.r);
  HS_EXPECT_EQ(plateau.color.g, c.color.g);
  HS_EXPECT_EQ(plateau.color.b, c.color.b);
  HS_EXPECT_NEAR(plateau.alpha, c.alpha, 1e-6f);
  Color4 swap = gc.grade(c, 0.0f);
  HS_EXPECT_EQ(swap.color.r, gc.gold.r);
  HS_EXPECT_EQ(swap.color.g, gc.gold.g);
  HS_EXPECT_EQ(swap.color.b, gc.gold.b);
  HS_EXPECT_NEAR(gc.opacity(0.0f), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(gc.opacity(1.0f), 1.0f, 1e-6f);
}
