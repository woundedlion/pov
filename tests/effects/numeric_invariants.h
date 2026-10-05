/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_effects.h.

// ---------------------------------------------------------------------------
// In-code-flagged numeric invariants with no oracle in the smoke harness
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for Comets' Lissajous-loop closing snap.
 * @details Befriended in effects/Comets.h. Reaches the private closing_domain()
 *          snap and the authored function table to verify every entry closes —
 *          path_fn(domain) == path_fn(0) — so the per-cycle drift reset never
 *          teleports the head. A wrong snap still renders a (discontinuous)
 *          curve, invisible to the smoke harness.
 */
struct CometsWhiteBox {
  /**
   * @brief Verifies every authored function table entry closes its loop.
   * @details For each entry the snapped endpoint must coincide with the t=0
   *          start (0,1,0), and the snap must stay positive (the floor-at-1
   *          guard keeps the head moving rather than freezing at path_fn(0)).
   */
  static void check_paths_close() {
    using C = Comets<DEFAULT_W, DEFAULT_H>;
    int idx = 0;
    for (const PresetEntry<math::LissajousParams> &entry : C::PRESETS) {
      const math::LissajousParams &cfg = entry.params;
      const float cd = C::closing_domain(cfg);
      HS_EXPECT_GT(cd, 0.0f); // floor-at-1 keeps the head moving
      // Every authored entry must clear the floor (m2*domain >= PI rounds to >= 1
      // closing cycle) so the floor never silently rewrites an authored domain.
      HS_EXPECT_GE(cfg.m2 * cfg.domain, math::PI_F);
      const math::Vector start = math::lissajous(cfg.m1, cfg.m2, cfg.a, 0.0f);
      const math::Vector end = math::lissajous(cfg.m1, cfg.m2, cfg.a, cd);
      const float gap = (end - start).magnitude();
      if (gap > 1e-3f)
        std::printf("  COMETS entry %d does not close: gap=%.6f\n", idx,
                    static_cast<double>(gap));
      HS_EXPECT_LT(gap, 1e-3f);
      ++idx;
    }
  }

  using C = Comets<SMALL_W, SMALL_H>;
  static constexpr int wipe_frames() { return C::WIPE_FRAMES; }
  static void roll_palette_over(C &c) { c.update_palette(); }
  static int wipe_frames_remaining(const C &c) {
    return c.wipe.frames_remaining;
  }
  static const GenerativePalette::Snapshot &palette_start(const C &c) {
    return c.wipe.start;
  }
  static const GenerativePalette::Snapshot &palette_target(const C &c) {
    return c.wipe.target;
  }
  static const math::Quaternion &node_orientation(const C &c) {
    return c.node->orientation.get();
  }
  static size_t trail_length(const C &c) { return c.node->trail.length(); }
  static int dwell_frames() { return C::PRESET_DWELL_FRAMES; }
};

/**
 * @brief Field-wise equality for two generative-palette snapshots.
 * @param a First snapshot.
 * @param b Second snapshot.
 * @return true when every authored key and axis range matches.
 */
inline bool palette_snapshots_equal(const GenerativePalette::Snapshot &a,
                                    const GenerativePalette::Snapshot &b) {
  if (a.key_count != b.key_count || a.lightness_low != b.lightness_low ||
      a.lightness_high != b.lightness_high || a.chroma_low != b.chroma_low ||
      a.chroma_high != b.chroma_high ||
      a.lightness_curve != b.lightness_curve ||
      a.chroma_curve != b.chroma_curve)
    return false;
  for (size_t i = 0; i < a.keys.size(); ++i)
    if (a.keys[i].bytes != b.keys[i].bytes)
      return false;
  return true;
}

/**
 * @brief Pins Comets' mid-wipe rollover skip.
 * @details A future dwell shorter than WIPE_FRAMES could trigger a rollover
 *          while a ColorWipe is still animating.
 *          The guard drops that rollover; without it the second wipe would
 *          overwrite palette_start/palette_target, which the live ColorWipe
 *          still holds references to. Both wipes still render, so the smoke
 *          pass cannot see the difference — assert the snapshots survive the
 *          mid-wipe call and that a rollover arms again once the wipe drains.
 */
inline void test_comets_rollover_skipped_mid_wipe() {
  using WB = CometsWhiteBox;
  reset_effect_globals();

  WB::C effect;
  effect.init();

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
  const GenerativePalette::Snapshot start = WB::palette_start(effect);
  const GenerativePalette::Snapshot target = WB::palette_target(effect);

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
  HS_EXPECT_TRUE(palette_snapshots_equal(start, WB::palette_start(effect)));
  HS_EXPECT_TRUE(palette_snapshots_equal(target, WB::palette_target(effect)));

  // step_wipe_rebake burns one armed frame plus WIPE_FRAMES stepped ones.
  for (int f = 0; f <= WB::wipe_frames(); ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), 0);

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
}

/** @brief Manual Comets selection restarts the authored path deterministically. */
inline void test_comets_manual_preset_restarts_path() {
  using WB = CometsWhiteBox;
  reset_effect_globals();

  WB::C effect;
  effect.init();
  for (int f = 0; f < 7; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }

  constexpr size_t PRESET = 4;
  HS_EXPECT_TRUE(effect.selectPreset(PRESET));
  HS_EXPECT_EQ(WB::node_orientation(effect), math::Quaternion());
  HS_EXPECT_EQ(WB::trail_length(effect), 0u);

  effect.draw_frame();
  effect.advance_display();
  const math::Quaternion first_step = WB::node_orientation(effect);
  HS_EXPECT_EQ(WB::trail_length(effect), 1u);

  for (int f = 0; f < 7; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(effect.selectPreset(PRESET));
  HS_EXPECT_EQ(WB::node_orientation(effect), math::Quaternion());
  HS_EXPECT_EQ(WB::trail_length(effect), 0u);

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(WB::node_orientation(effect), first_step);

  // Manual selection pauses the automatic dwell countdown as well as
  // selecting the requested path.
  for (int f = 0; f < WB::dwell_frames(); ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(effect.getPresetIndex(), PRESET);
}

/**
 * @brief Pins that AshCloud's value cutout reaches the rendered frame.
 * @details AshCloud is the roster's only FieldCoverageKind::VALUE_CUTOUT
 *          tenant. The kernel has a parity oracle against the chain
 *          interpreter (tests/test_shader_chain.h), but that exercises the
 *          stage standalone, not the composed wiring that feeds it: a pipeline
 *          that dropped the coverage stage would still render, and every other
 *          check would stay green. Sweeping cutout-threshold across its whole
 *          authored range must open the frame at one end and close it at the
 *          other.
 */
inline void test_ash_cloud_value_cutout_gates_the_frame() {
  using FX = AshCloud<SMALL_W, SMALL_H>;

  auto lit_pixels = [](float threshold) {
    reset_effect_globals();
    FX effect;
    effect.init();

    auto snapshot = effect.serialize_parameters();
    snapshot.params.template get<"value">().cutout_threshold = threshold;
    snapshot.params.template get<"value">().cutout_softness = 1.0f / 1024.0f;
    HS_EXPECT_TRUE(effect.restore_parameters(snapshot));

    effect.draw_frame();
    effect.advance_display();

    size_t lit = 0;
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x) {
        const Pixel &p = effect.get_pixel(x, y);
        lit += (p.r != 0 || p.g != 0 || p.b != 0) ? 1u : 0u;
      }
    return lit;
  };

  HS_EXPECT_GT(lit_pixels(0.0f), size_t{0});
  HS_EXPECT_EQ(lit_pixels(1.0f), size_t{0});
}

/**
 * @brief White-box accessor for the Thrusters warp curve and fire path.
 * @details Befriended in effects/Thrusters.h. Reaches the private warp_decay()
 *          curve to pin its shift-and-renormalized endpoints: a bare
 *          0.7*exp(-2t) would bottom out at ~0.095 and leave a residual wobble,
 *          which still renders and so passes the smoke harness. Also drives
 *          on_fire_thruster() directly: the roster smoke sweep never reaches a
 *          fire at the local 8-frame default, so the opposed pair, the
 *          spin-axis fallback and the FIFO pairing are otherwise unobserved.
 */
struct ThrustersWhiteBox {
  using FX = Thrusters<DEFAULT_W, DEFAULT_H>;
  using Slot = typename FX::ThrusterContext;

  /**
   * @brief Verifies warp_decay peaks at 0.7 at t=0 and relaxes to exactly 0 at t=1.
   */
  static void check_warp_endpoints() {
    HS_EXPECT_NEAR(FX::warp_decay(0.0f), 0.7f, 1e-6f); // peak at fire
    HS_EXPECT_NEAR(FX::warp_decay(1.0f), 0.0f, 1e-6f); // full relaxation by end
  }

  /**
   * @brief Verifies one fire spawns two thrusters a half turn apart on the ring.
   * @details The pair's polar angles differ (the warp shifts each radially), so
   *          the opposition is in azimuth about the ring axis, not antipodality.
   */
  static void check_fire_spawns_opposed_pair() {
    reset_effect_globals();
    FX effect;
    effect.init();

    effect.on_fire_thruster();
    HS_EXPECT_EQ(effect.thrusters.size(), size_t{2});

    const math::Vector a =
        azimuth_dir(effect.thrusters[0].point, effect.ring_vec);
    const math::Vector b =
        azimuth_dir(effect.thrusters[1].point, effect.ring_vec);
    HS_EXPECT_NEAR(math::dot(a, b), -1.0f, 1e-4f);
  }

  /**
   * @brief Verifies a collapsed ring still fires.
   * @details theta_eq is radius * pi/2, so a sub-microradian ring puts both
   *          thrust points on ring_vec and drops the spin-axis cross product
   *          below EPS_NORMALIZE_SQ. Without the normalized_or() fallback the
   *          normalize() inside the fire traps, so reaching the assertions is
   *          the check.
   */
  static void check_collapsed_ring_falls_back_to_a_spin_axis() {
    reset_effect_globals();
    FX effect;
    effect.init();

    effect.params.radius = 1e-7f;
    effect.on_fire_thruster();
    HS_EXPECT_EQ(effect.thrusters.size(), size_t{2});
    HS_EXPECT_NEAR(math::dot(effect.thrusters[0].point, effect.ring_vec), 1.0f,
                   1e-6f);
    HS_EXPECT_NEAR(math::dot(effect.thrusters[1].point, effect.ring_vec), 1.0f,
                   1e-6f);
  }

  /**
   * @brief Verifies the slot pool saturates and evicts whole fires in order.
   */
  static void check_fifo_evicts_the_oldest_pair() {
    reset_effect_globals();
    FX effect;
    effect.init();

    constexpr size_t CAPACITY = effect.thrusters.capacity();
    constexpr int FIRES = CAPACITY / 2 + 1;
    math::Vector lead[FIRES];
    for (int fire = 0; fire < FIRES; ++fire) {
      effect.on_fire_thruster();
      const size_t spawned = 2u * static_cast<size_t>(fire + 1);
      const size_t expected = spawned < CAPACITY ? spawned : CAPACITY;
      HS_EXPECT_EQ(effect.thrusters.size(), expected);
      // The fire's first spawn: back() is its opposed twin.
      lead[fire] = effect.thrusters[effect.thrusters.size() - 2].point;
    }
    // The last fire pushed two slots into a full pool, retiring fire 0's pair
    // and leaving fire 1's first slot at the front.
    HS_EXPECT_NEAR(math::dot(effect.thrusters.front().point, lead[1]), 1.0f,
                   1e-6f);
  }

  /**
   * @brief Verifies draw_frame retires the expired prefix and ages the rest.
   */
  static void check_expired_slots_retire_by_pair() {
    reset_effect_globals();
    FX effect;
    effect.init();

    effect.on_fire_thruster();
    effect.on_fire_thruster();
    HS_EXPECT_EQ(effect.thrusters.size(), size_t{4});

    // The first fire's pair is spent; the second's has one frame left.
    effect.thrusters[0].age = Slot::LIFE;
    effect.thrusters[1].age = Slot::LIFE;
    effect.thrusters[2].age = Slot::LIFE - 1;
    effect.thrusters[3].age = Slot::LIFE - 1;

    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_EQ(effect.thrusters.size(), size_t{2});
    HS_EXPECT_EQ(effect.thrusters.front().age, Slot::LIFE);

    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_TRUE(effect.thrusters.is_empty());
  }

private:
  /** @brief Unit azimuth of @p p about @p axis, with the radial warp removed. */
  static math::Vector azimuth_dir(const math::Vector &p,
                                  const math::Vector &axis) {
    return (p - axis * math::dot(p, axis)).normalized();
  }
};

/**
 * @brief White-box accessor for RingShower's radius easing endpoints.
 * @details Befriended in effects/RingShower.h. Reaches the private Ring type to
 *          pin its age+1 convention: the ring must reach RADIUS_MAX on its final
 *          visible frame (age+1 == life) and render a non-zero first step rather
 *          than radius 0. An off-by-one in the convention still renders a ring.
 */
struct RingShowerWhiteBox {
  /**
   * @brief Verifies radius_at hits RADIUS_MAX on the final frame and is non-zero
   *        on the first.
   */
  static void check_radius_endpoints() {
    using RS = RingShower<DEFAULT_W, DEFAULT_H>;
    typename RS::Ring ring;
    ring.life = 50;
    ring.age = ring.life - 1; // final visible frame: age+1 == life -> t == 1
    HS_EXPECT_NEAR(ring.radius_at(), RS::Ring::RADIUS_MAX, 1e-5f);
    ring.age = 0; // first frame: one eased step in, not radius 0
    HS_EXPECT_GT(ring.radius_at(), 0.0f);
    HS_EXPECT_LT(ring.radius_at(), RS::Ring::RADIUS_MAX);
  }
};

/**
 * @brief White-box accessor for Dynamo's overlapping-wipe band ordering.
 * @details Stages inverted overlapping-wipe boundaries and pins finite alpha
 *          in [0,1] across the angular span. Palette indexing is bounded by the
 *          boundary container's capacity.
 */
struct DynamoWhiteBox {
  /**
   * @brief Verifies color() is memory-safe and in-range under inverted bands.
   */
  static void check_overlapping_wipes_stay_in_range() {
    // Dynamo::init() bakes from persistent_arena and schedules on the shared
    // global timeline, so reset the shared globals as smoke_one does.
    reset_effect_globals();

    Dynamo<DEFAULT_W, DEFAULT_H> effect;
    effect.init();

    const Color4 opening = effect.color(math::Z_AXIS, 0.5f);
    effect.color_wipe();
    HS_EXPECT_EQ(effect.color(math::Z_AXIS, 0.5f).color, opening.color);
    effect.palette_boundaries.front() = math::PI_F + effect.WIPE_BLEND_WIDTH;
    HS_EXPECT_EQ(effect.color(math::Z_AXIS * -1.0f, 0.5f).color,
                 effect.baked_palettes[0].get(0.5f).color);
    effect.color_wipe();
    HS_EXPECT_EQ(effect.palette_boundaries.size(), static_cast<size_t>(2));

    // Force the documented worst case: the newer band (index 0) has overtaken
    // the older (index 1) -> non-monotonic order. The chosen magnitudes also push
    // the scan past the first iteration into the second boundary and the
    // baked_palettes[i+1] access for part of the sweep, exercising the bounds
    // path the safety claim rests on.
    effect.palette_boundaries[0] = 1.5f; // newer, overtaken ahead of the older
    effect.palette_boundaries[1] = 0.5f; // older, left behind

    // PALETTE_NORMAL is Z_AXIS, so v = (sin theta, 0, cos theta) sweeps
    // angle_between(v, normal) across the full [0, PI] band span.
    constexpr int STEPS = 256;
    for (int i = 0; i <= STEPS; ++i) {
      float theta =
          math::PI_F * static_cast<float>(i) / static_cast<float>(STEPS);
      math::Vector v(std::sin(theta), 0.0f, std::cos(theta));
      Color4 c = effect.color(v, 0.5f);
      HS_EXPECT_TRUE(std::isfinite(c.alpha));
      HS_EXPECT_GE(c.alpha, 0.0f);
      HS_EXPECT_LE(c.alpha, 1.0f);
    }
  }

  using D = Dynamo<DEFAULT_W, DEFAULT_H>;
  using Ring = Filter::World::Trails<D::TRAIL_CAPACITY>;
  static constexpr int trail_capacity() { return D::TRAIL_CAPACITY; }
  static constexpr int trail_len_max() { return D::TRAIL_LEN_MAX; }
  static void set_speed(D &d, float v) { d.params.speed = v; }
  static void set_trail_length(D &d, float v) { d.params.trail_length = v; }
  static float trail_ceiling(const D &d) { return d.params.trail_ceiling; }
  static uint32_t emitted_points(const D &d) { return d.emitted_points; }
  static uint32_t points_per_emission(const D &d) {
    return d.points_per_emission;
  }
  static size_t trail_points(D &d) { return d.filters.get<Ring>().size(); }
};

/**
 * @brief Pins Dynamo's trail-ring ceiling as the thing that keeps the ring from
 *        evicting.
 * @details trail_length_ceiling() caps the live trail at what the ring can hold
 *          for the current emission rate; past it the ring overruns and evicts
 *          points of arbitrary age (flush()'s compaction leaves it unordered),
 *          punching holes in the tail rather than shortening it — corruption
 *          that still renders, so the smoke pass never sees it. Drive the
 *          worst case both sliders allow and assert the requested trail really
 *          would have overrun while the ceiling keeps steady-state occupancy
 *          inside the ring.
 */
inline void test_dynamo_trail_ceiling_bounds_the_ring() {
  using WB = DynamoWhiteBox;
  reset_effect_globals();

  WB::D effect;
  effect.init();

  constexpr float MAX_SPEED = 10.0f; // "Speed" slider bound
  constexpr float MAX_TRAIL = WB::trail_len_max();
  WB::set_speed(effect, MAX_SPEED);
  WB::set_trail_length(effect, MAX_TRAIL);

  const auto capacity = static_cast<uint32_t>(WB::trail_capacity());
  for (int f = 0; f < static_cast<int>(MAX_TRAIL) + 2; ++f) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_LT(WB::trail_points(effect), static_cast<size_t>(capacity));
  }

  // Points buffered per frame at this speed: points_per_emission per whole step,
  // and speed_accumulator carries at most one extra step into a frame.
  const uint32_t per_frame =
      WB::points_per_emission(effect) * (static_cast<uint32_t>(MAX_SPEED) + 1);
  HS_EXPECT_GT(per_frame, 0u);
  HS_EXPECT_GT(per_frame * static_cast<uint32_t>(MAX_TRAIL), capacity);
  HS_EXPECT_GT(WB::trail_points(effect), size_t{0});
  HS_EXPECT_LT(WB::trail_points(effect), static_cast<size_t>(capacity));
}

/**
 * @brief Pins what Dynamo's emitted_points counter actually measures.
 * @details draw_nodes() bumps it from inside the strand's fragment shader, so
 *          the trail-ring ceiling derived from it is only a valid bound while
 *          the rasterizer runs that shader exactly once per point it hands to
 *          the pipeline. Assert that directly: with a trail long enough that
 *          nothing ages out over the window, the ring's growth across a frame
 *          must equal the frame's emitted_points. A rasterizer that shaded a
 *          fragment it then dropped, or plotted one it never shaded, would
 *          silently rescale the bound.
 */
inline void test_dynamo_emitted_points_counts_ring_seeds() {
  using WB = DynamoWhiteBox;
  reset_effect_globals();

  WB::D effect;
  effect.init();
  WB::set_trail_length(effect, 100.0f);
  WB::set_speed(effect, 1.0f);
  constexpr int FRAMES = 3;

  for (int f = 0; f < FRAMES; ++f) {
    const size_t before = WB::trail_points(effect);
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_GT(WB::trail_ceiling(effect), static_cast<float>(FRAMES));
    const size_t seeded = WB::trail_points(effect) - before;
    HS_EXPECT_GT(WB::emitted_points(effect), 0u);
    HS_EXPECT_EQ(seeded, static_cast<size_t>(WB::emitted_points(effect)));
  }
}

/**
 * @brief White-box accessor for HopfFibration's S3-lift + stereographic
 *        projection (befriended in effects/HopfFibration.h).
 * @details The smoke/determinism harness only proves the effect renders and
 *          reproduces; it never pins hopf_project()'s numeric output. This seam
 *          sets the private per-frame cache (tumble sines/cosines, fold/flow/
 *          tumble-y phases) and a fiber's base coordinates, then calls the real
 *          projection so a test can compare it to the closed form.
 */
struct HopfWhiteBox {
  using HF = HopfFibration<DEFAULT_W, DEFAULT_H>;
  static size_t fiber_count() { return HF::ACTUAL_FIBERS; }
  static math::Vector project(HF &fx, size_t i) { return fx.hopf_project(i); }
  static void set_fiber(HF &fx, size_t i, float azimuth, float polar) {
    fx.fibers[i] = math::Spherical(azimuth, polar);
  }
  static void set_cache(HF &fx, float cx, float sx, float cy, float sy,
                        float fold_base, float flow_rad, float ty_rad) {
    fx.cx = cx;
    fx.sx = sx;
    fx.cy = cy;
    fx.sy = sy;
    fx.fold_base = fold_base;
    fx.flow_rad = flow_rad;
    fx.ty_rad = ty_rad;
  }
  static void set_folding(HF &fx, float v) { fx.params.folding = v; }
  static void set_twist(HF &fx, float v) { fx.params.twist = v; }
  static size_t trim_start(size_t len, float alpha) {
    return HF::trail_trim_start(len, alpha);
  }
};

/**
 * @brief Pins HopfFibration::hopf_project against a closed form and asserts its
 *        unit-direction/finite invariant under a nontrivial 4D tumble.
 * @details Identity tumble with zero folding and twist reduces fiber 0
 *          (beta == 0) to a plain S3 lift, whose stereographic image is the
 *          normalized (q0, q1, q2); the check uses the same fast trig the effect
 *          does so it pins the projection, not the trig approximation. The second
 *          pass exercises every fiber under active tumble/fold/twist and requires
 *          each result be finite, unit-length (normalized_or), and deterministic.
 */
inline void test_hopf_projection_math() {
  // init() bakes from persistent_arena and schedules on the shared timeline, so
  // reset the shared globals as smoke_one does.
  reset_effect_globals();
  using WB = HopfWhiteBox;

  HopfFibration<DEFAULT_W, DEFAULT_H> fx;
  fx.init();

  // Identity tumble + no folding/twist: fiber 0 has beta == 0, so q3 == 0 and the
  // projection collapses to the normalized (q0, q1, q2).
  WB::set_folding(fx, 0.0f);
  WB::set_twist(fx, 0.0f);
  WB::set_cache(fx, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);

  const float polar = 1.2f, azimuth = 0.7f;
  WB::set_fiber(fx, 0, azimuth, polar);
  const math::Vector v = WB::project(fx, 0);

  const float eta = polar * 0.5f;
  const float ce = math::fast_cosf(eta), se = math::fast_sinf(eta);
  const math::Vector expected = math::Vector(ce * math::fast_cosf(azimuth),
                                             ce * math::fast_sinf(azimuth), se)
                                    .normalized();
  HS_EXPECT_NEAR(v.x, expected.x, 1e-4f);
  HS_EXPECT_NEAR(v.y, expected.y, 1e-4f);
  HS_EXPECT_NEAR(v.z, expected.z, 1e-4f);
  HS_EXPECT_NEAR(v.magnitude(), 1.0f, 1e-4f);

  // Every fiber projects to a finite unit direction under active tumble/fold/twist.
  const float ax = 0.9f, ay = 0.4f;
  WB::set_cache(fx, math::fast_cosf(ax), math::fast_sinf(ax),
                math::fast_cosf(ay), math::fast_sinf(ay),
                math::fast_sinf(ax * 0.5f) * 0.5f, 0.3f, ay);
  WB::set_folding(fx, 0.5f);
  WB::set_twist(fx, 2.0f);
  for (size_t i = 0; i < WB::fiber_count(); ++i) {
    const math::Vector p = WB::project(fx, i);
    HS_EXPECT_VEC(WB::project(fx, i), p, 0.0f);
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.magnitude(), 1.0f, 1e-3f);
  }
}

/**
 * @brief Checks HopfFibration's supported trail lengths at the alpha samples
 *        listed below.
 * @details render_trails() stages points [first, len) and rasterizes
 *          len - first fragments, so a trim that reached len would bind an
 *          empty polyline; the visibility gate ahead of it guarantees
 *          alpha >= MIN_VISIBLE_ALPHA, which clears the trim's own
 *          MIN_ENCODABLE_ALPHA floor and so bounds first at len - 2.
 *          The trim must match its per-sample cutoff and be monotone in alpha:
 *          a brighter trail can never show less of its tail.
 */
inline void test_hopf_trail_trim_keeps_a_segment() {
  using WB = HopfWhiteBox;
  using HF = HopfFibration<DEFAULT_W, DEFAULT_H>;
  // Ascending, and all at or above the gate render_trails() applies first.
  const float alphas[] = {MIN_VISIBLE_ALPHA, 0.005f, 0.01f, 0.05f, 0.4f, 1.0f};
  for (size_t len = 2; len <= static_cast<size_t>(HF::TRAIL_LEN); ++len) {
    size_t brighter = len;
    for (float alpha : alphas) {
      const size_t first = WB::trim_start(len, alpha);
      HS_EXPECT_LT(first + 1, len);
      HS_EXPECT_LE(first, brighter);
      // The kept segment meets the cutoff; its predecessor falls below it.
      HS_EXPECT_GE(static_cast<float>(first + 1) / (len - 1) * alpha,
                   MIN_ENCODABLE_ALPHA);
      if (first > 0)
        HS_EXPECT_LT(static_cast<float>(first) / (len - 1) * alpha,
                     MIN_ENCODABLE_ALPHA);
      brighter = first;
    }
  }
}
