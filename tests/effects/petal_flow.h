/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief PetalFlow white-box cases. */

/**
 * @brief White-box accessor for PetalFlow's spawn-gap accumulator and hue
 *        cursor.
 */
struct PetalFlowWhiteBox {
  using PF = PetalFlow<DEFAULT_W, DEFAULT_H>;
  /** @brief Exercises the final free slot or the fully occupied pool. */
  static void fill_and_spawn(PF &pf, bool full) {
    for (auto &ring : pf.rings)
      ring.active = true;
    pf.rings[PF::MAX_RINGS - 1].active = full;
    pf.spawn_ring_at_pos(PF::START_RHO);
  }
  /** @brief Renders one production ring or its independently frozen geometry.
   * @param timed Uses the production temporal draw.
   * @param fraction Frozen frame fraction.
   * @param displacement Vertical frozen-reference offset in pixels.
   * @param rho Current ring travel.
   * @return Row-major linear luminance. */
  static std::vector<float> render_temporal_ring(bool timed, float fraction,
                                                 float displacement,
                                                 float rho) {
    reset_effect_globals();
    PF fx;
    fx.init();
    const float previous_rho = rho - fx.move_dist();
    const math::Vector axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
    const auto start_pose = math::make_rotation(math::X_AXIS, 0.27f);
    fx.orientation.set(start_pose);
    fx.orientation.push(math::make_rotation(axis, 0.012f) * start_pose);
    typename PF::Ring ring{rho, 0.25f, true, false};
    {
      Canvas canvas(fx);
      if (timed) {
        Spacetime::geometry_counters = {};
        fx.draw_ring(canvas, ring, previous_rho);
        HS_EXPECT_EQ(Spacetime::geometry_counters.column_draw_setups, 1u);
      } else {
        ScratchScope guard(scratch_arena_a);
        const float travel = previous_rho + (rho - previous_rho) * fraction;
        const float first_twist =
            previous_rho / fx.spacing() * fx.params.twist_factor;
        const float last_twist = rho / fx.spacing() * fx.params.twist_factor;
        const float twist = first_twist + (last_twist - first_twist) * fraction;
        const auto pose = fx.orientation.get(fraction);
        fx.orientation.set(pose);
        auto *vertices = scratch_arena_a.make_n<math::Vector>(PF::NUM_SAMPLES);
        for (int i = 0; i < PF::NUM_SAMPLES; ++i) {
          vertices[i] = fx.ring_vertex(i, expf(travel), twist);
          if (displacement != 0.0f) {
            const auto world = math::rotate(vertices[i], pose);
            vertices[i] =
                math::rotate(hs_test::displace_reference_rows<DEFAULT_H>(
                                 world, displacement),
                             pose.conjugate());
          }
        }
        Color4 shade = fx.palette.get(ring.hue).color;
        shade.alpha = PF::pole_opacity(travel) * fx.params.alpha;
        auto shader = [&](const math::Vector &, Fragment &f) {
          f.color = shade;
        };
        Plot::Mesh::draw_columns<DEFAULT_W, DEFAULT_H>(
            fx.filters, canvas,
            Plot::Mesh::VertexSpan{{vertices, size_t(PF::NUM_SAMPLES)}},
            fx.ring_edges, shader, PF::RING_GEOMETRY, 0.0f);
      }
    }
    fx.advance_display();
    return capture_luminance<DEFAULT_W, DEFAULT_H>(fx);
  }
  static float gap(const PF &pf) { return pf.gap_accumulator; }
  static float start_rho() { return PF::START_RHO; }
  static float speed_max() { return PF::SPEED_MAX; }
  static float density_max() { return PF::DENSITY_MAX; }
  static float youngest_rho(const PF &pf) {
    float youngest = PF::END_RHO;
    for (const auto &ring : pf.rings)
      if (ring.active && ring.rho < youngest)
        youngest = ring.rho;
    return youngest;
  }
  static float next_hue(const PF &pf) { return pf.next_hue; }
  static float live_spacing(const PF &pf) { return pf.spacing(); }
  static float move_dist(const PF &pf) { return pf.move_dist(); }
  static int spawns(const PF &pf) { return pf.spawns; }
  static int max_rings() { return PF::MAX_RINGS; }
  static int rings_on_path() { return PF::RINGS_ON_PATH; }
  static int active_rings(const PF &pf) {
    int active = 0;
    for (int i = 0; i < PF::MAX_RINGS; ++i)
      if (pf.rings[i].active)
        ++active;
    return active;
  }
};

/** @brief PetalFlow can spawn into the final free ring slot. */
inline void test_petalflow_last_free_slot() {
  reset_effect_globals();
  using WB = PetalFlowWhiteBox;
  WB::PF effect;
  effect.init();
  const int before = WB::spawns(effect);
  WB::fill_and_spawn(effect, false);
  HS_EXPECT_EQ(WB::spawns(effect), before + 1);
  HS_EXPECT_EQ(WB::active_rings(effect), WB::max_rings());
}

/**
 * @brief Verifies the spawn-gap accumulator drains every frame, the hue cursor
 *        stays wrapped, and the ring pool absorbs the worst slider corner.
 * @details After any frame the residue must satisfy 0 <= gap < spacing() and
 *          next_hue must stay in [0, 1).
 */
inline void test_petalflow_spawn_gap_bounded() {
  using WB = PetalFlowWhiteBox;
  reset_effect_globals();
  WB::PF pf;
  pf.init();
  HS_EXPECT_NEAR(WB::youngest_rho(pf) - WB::start_rho(), WB::gap(pf), 1e-5f);
  HS_EXPECT_EQ(pf.updateParameter("Speed", WB::speed_max()),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(pf.updateParameter("Density", WB::density_max()),
               ParamSetResult::APPLIED);

  const float spacing = WB::live_spacing(pf);
  // Below this the while-loop emits at most one ring per frame, whatever the
  // frame count.
  HS_EXPECT_GT(WB::move_dist(pf), spacing);

  // Long enough for the path prefilled at the default density to be fully
  // replaced at the tightest one, so the pool bound is read at steady state
  // rather than mid-fill.
  const int frames = smoke_frames() < 128 ? 128 : smoke_frames();
  int worst_active = 0;
  int worst_burst = 0;
  int previous_spawns = WB::spawns(pf);
  for (int f = 0; f < frames; ++f) {
    pf.draw_frame();
    pf.advance_display();
    const float gap = WB::gap(pf);
    HS_EXPECT_GE(gap, 0.0f);
    HS_EXPECT_LT(gap, spacing); // accumulator never runs away past one spacing
    const float hue = WB::next_hue(pf);
    HS_EXPECT_GE(hue, 0.0f);
    HS_EXPECT_LT(hue, 1.0f); // hue cursor stays wrapped
    const int spawns = WB::spawns(pf);
    worst_burst = std::max(worst_burst, spawns - previous_spawns);
    previous_spawns = spawns;
    worst_active = std::max(worst_active, WB::active_rings(pf));
  }
  HS_EXPECT_GE(worst_burst, 2); // the multi-spawn branch actually ran
  HS_EXPECT_LE(worst_active, WB::rings_on_path());
  HS_EXPECT_LE(WB::rings_on_path(), WB::max_rings());
}
