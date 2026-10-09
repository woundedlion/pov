/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Scan::Ring::draw — SDF rasterize() path through a Pipeline sink
// ============================================================================

/**
 * @brief Verifies the SDF rasterize() path plots a complete hollow ring.
 * @details Every pixel inside the stroke core is lit and every lit pixel lies
 * within the stroke plus its AA fringe.
 */
inline void test_ring_rasterize_produces_bounded_output() {
  constexpr int W = 64, H = 48;
  constexpr float radius = 0.5f, thickness = 0.4f;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe; // bare 2D sink (no filters)

  {
    Canvas c(fx);
    math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
    Scan::Ring::draw<W, H, false>(pipe, c, basis, radius, thickness,
                                  [](const math::Vector &, Fragment &f) {
                                    f.color = Color4(Pixel(60000, 60000, 60000),
                                                     1.0f);
                                  });
  }
  fx.advance_display();

  const float target = radius * (math::PI_F / 2.0f); // ring-centre polar angle
  const float row = math::PI_F / (H - 1); // angular height of one row
  // Two rows of quantization either side of the nominal stroke: inside the
  // inner bound the ring is solid, outside the outer one it is clear.
  const float core = thickness - 2.0f * row;
  const float band = thickness + 2.0f * row;
  size_t plotted = 0, solid = 0;
  float widest_lit = 0.0f, nearest_dark = math::PI_F;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      HS_CONTEXT("px", x, y);
      const math::Vector v = math::pixel_to_vector<W, H>(x, y);
      const float offset = fabsf(acosf(hs::clamp(v.y, -1.0f, 1.0f)) - target);
      const bool lit = !is_black(fx.get_pixel(x, y));
      if (lit) {
        ++plotted;
        widest_lit = std::max(widest_lit, offset);
        HS_EXPECT_LE(offset, band);
      } else {
        nearest_dark = std::min(nearest_dark, offset);
      }
      if (offset <= core) {
        ++solid;
        HS_EXPECT_TRUE(lit);
      }
    }
  std::printf("  [ring mask] lit %zu/%d, core %zu, widest lit %.4f, nearest "
              "dark %.4f of half-width %.4f\n",
              plotted, W * H, solid, (double)widest_lit, (double)nearest_dark,
              (double)thickness);
  HS_EXPECT_GT(solid, (size_t)0);
  HS_EXPECT_LT(plotted, (size_t)(W * H));
}

/**
 * @brief Verifies a radius past 1 shades the azimuth of the caller's own frame.
 * @details SDF::Ring covers the whole [0, 2] radius range directly.
 */
inline void test_ring_long_radius_azimuth_unflipped() {
  constexpr int W = 96, H = 64;
  constexpr float radius = 1.4f, thickness = 0.15f;
  const math::Basis basis = math::make_basis(
      math::Quaternion(), math::Vector(0.3f, 0.8f, -0.5f).normalized());

  size_t lit = 0;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Scan::Ring::draw<W, H>(
        pipe, c, basis, radius, thickness,
        [&](const math::Vector &p, Fragment &f) {
          const double u = static_cast<double>(p.x) * basis.u.x +
                           static_cast<double>(p.y) * basis.u.y +
                           static_cast<double>(p.z) * basis.u.z;
          const double w = static_cast<double>(p.x) * basis.w.x +
                           static_cast<double>(p.y) * basis.w.y +
                           static_cast<double>(p.z) * basis.w.z;
          const double turn = std::atan2(w, u) / (2.0 * 3.14159265358979323846);
          double delta = std::fabs(f.v0 - (turn < 0 ? turn + 1 : turn));
          delta = std::min(delta, 1.0 - delta);
          HS_EXPECT_LT(delta, 0.00061);
          ++lit;
          f.color = Color4(Pixel(60000, 60000, 60000), f.v2);
        });
  }
  HS_EXPECT_GT(lit, (size_t)0);
}

/**
 * @brief Verifies a degenerate clip band makes rasterize plot nothing.
 */
inline void test_ring_rasterize_empty_clip_draws_nothing() {
  constexpr int W = 64, H = 48;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;

  fx.set_clip(30, 30, 0, W);
  fx.set_margin(0);

  {
    Canvas c(fx);
    math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
    Scan::Ring::draw<W, H, false>(
        pipe, c, basis, 0.5f, 0.4f, [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();

  const size_t plotted = count_lit_region<W, H>(fx);
  HS_EXPECT_EQ(plotted, (size_t)0);
}

/**
 * @brief Verifies flat-ring rasterization matches zero knots within one code.
 */
inline void test_distorted_ring_flat_matches_zero_knot_raster() {
  constexpr int W = 96, H = 64, LUT_N = 16;
  float knots[LUT_N + 1] = {};
  math::Basis basis = math::make_basis(
      math::Quaternion(), math::Vector(0.3f, 0.8f, -0.5f).normalized());

  auto check = [&](bool partial_clip) {
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(51000, 27000, 9000), 0.8f * f.v2);
    };
    std::vector<Pixel> expected(W * H);
    {
      hs_test::StubEffect reference(W, H);
      if (partial_clip) {
        reference.set_clip(9, 53, 17, 81);
        reference.set_margin(0);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(reference);
        Scan::DistortedRing::draw<W, H>(pipeline, canvas, basis, 1.7f, 0.12f,
                                        knots, LUT_N, shader);
      }
      reference.advance_display();
      capture_frame<W, H>(reference, expected);
    }

    hs_test::StubEffect flat(W, H);
    if (partial_clip) {
      flat.set_clip(9, 53, 17, 81);
      flat.set_margin(0);
    }
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(flat);
      Scan::DistortedRing::draw_flat<W, H, false>(pipeline, canvas, basis, 1.7f,
                                                  0.12f, shader);
    }
    flat.advance_display();

    size_t lit = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = flat.get_pixel(x, y);
        if (!is_black(a))
          ++lit;
        HS_EXPECT_NEAR(static_cast<int>(a.r), static_cast<int>(b.r), 1);
        HS_EXPECT_NEAR(static_cast<int>(a.g), static_cast<int>(b.g), 1);
        HS_EXPECT_NEAR(static_cast<int>(a.b), static_cast<int>(b.b), 1);
      }
    }
    // Both paths drawing nothing would satisfy every comparison above.
    HS_EXPECT_GT(lit, (size_t)0);
  };

  check(false);
  check(true);
}

// AA-dust ceiling of the fused ring-group scan under -ffast-math: a handful of
// stroke-edge pixels at 1/128 of full scale.
constexpr int GROUP_MAX_CHANNEL_DELTA = 512;
constexpr int GROUP_MAX_DIFF_PIXELS = 16;

/**
 * @brief Compares fused and sequential ring rendering, bit-identical under IEEE.
 * @details Under -ffast-math, per-loop reassociation differences are bounded by
 *          GROUP_MAX_* tolerances.
 */
inline void test_ring_group_matches_sequential() {
  constexpr int W = 96, H = 64;
  constexpr int N = 4;
  const ScopedPoleLod lod(0.0f);

  auto run_case = [&](const math::Vector &normal, bool partial_clip) {
    const float ths[N] = {0.08f, 0.04f, 0.04f, 0.08f};
    const Color4 colors[N] = {Color4(Pixel(60000, 10000, 5000), 0.9f),
                              Color4(Pixel(5000, 60000, 10000), 0.6f),
                              Color4(Pixel(10000, 5000, 60000), 0.4f),
                              Color4(Pixel(30000, 30000, 30000), 0.7f)};
    SDF::Ring shapes[N];
    for (int s = 0; s < N; ++s) {
      math::Quaternion q = math::make_rotation(
          math::Vector(0.2f, 0.5f, 0.8f).normalized(), 0.02f * s);
      shapes[s] = SDF::Ring(math::make_basis(q, normal), 1.0f, ths[s]);
    }

    std::vector<Pixel> expected(W * H);
    {
      hs_test::StubEffect seq(W, H);
      if (partial_clip) {
        seq.set_clip(9, 53, 17, 81);
        seq.set_margin(0);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(seq);
        for (int s = 0; s < N; ++s) {
          auto shader = [&](const math::Vector &, Fragment &f) {
            f.color = colors[s];
          };
          Scan::rasterize<W, H, false>(pipeline, canvas, shapes[s], shader);
        }
      }
      seq.advance_display();
      capture_frame<W, H>(seq, expected);
    }

    hs_test::StubEffect fused(W, H);
    if (partial_clip) {
      fused.set_clip(9, 53, 17, 81);
      fused.set_margin(0);
    }
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fused);
      Scan::RingGroup::draw<W, H>(pipeline, canvas, shapes, N,
                                  [&](int s, const math::Vector &,
                                      Fragment &f) { f.color = colors[s]; });
    }
    fused.advance_display();

    int diff_px = 0;
    int worst_delta = 0;
    size_t lit = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = fused.get_pixel(x, y);
        if (a.r != 0 || a.g != 0 || a.b != 0)
          ++lit;
        if (a.r == b.r && a.g == b.g && a.b == b.b)
          continue;
        ++diff_px;
        const int deltas[3] = {std::abs(static_cast<int>(a.r) - b.r),
                               std::abs(static_cast<int>(a.g) - b.g),
                               std::abs(static_cast<int>(a.b) - b.b)};
        for (int d : deltas) {
          worst_delta = std::max(worst_delta, d);
          HS_EXPECT_LE(d, GROUP_MAX_CHANNEL_DELTA);
        }
      }
    }
    std::printf("  [ring-group] diff_px=%d worst_delta=%d\n", diff_px,
                worst_delta);
    HS_EXPECT_LE(diff_px, GROUP_MAX_DIFF_PIXELS);
    HS_EXPECT_GT(lit, size_t(200));
#if !defined(HS_TEST_FAST_MATH)
    HS_EXPECT_EQ(diff_px, 0);
#endif
  };

  run_case(math::Vector(0.3f, 0.8f, -0.5f).normalized(), false);
  run_case(math::Vector(0.3f, 0.8f, -0.5f).normalized(), true);
  // Near-pole axis: r_val under the horizontal-projection floor forces the
  // group's full-row-scan fallback.
  run_case(math::Vector(0.005f, 1.0f, 0.0f).normalized(), false);
}

/** @brief Distorted rings wholly past a pole add no candidate cells. */
inline void test_distorted_ring_candidates_outside_poles() {
  constexpr int W = 32, H = 16, KNOTS = 8;
  const math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  const int8_t slots[] = {0};
  for (const float offset : {-4.0f, 4.0f}) {
    float knots[KNOTS + 1];
    std::fill_n(knots, KNOTS + 1, offset);
    const SDF::DistortedRing shapes[] = {
        SDF::DistortedRing(basis, 1.0f, 0.01f, knots, KNOTS, 0.0f, nullptr)};
    Scan::DistortedRingStack::CandidateTable<W, H> table;
    Scan::DistortedRingStack::build_candidate_table<W, H>(1, shapes, slots,
                                                          table);
    for (const auto &cell : table.cells)
      HS_EXPECT_GT(cell.lo, cell.hi);
  }
}

/**
 * @brief Verifies DistortedRingStack::draw matches rasterizing the stack's
 *        rings one by one.
 * @details The per-ring path uses suppress_pole_fill. Scoped to
 * Render::pole_lod_aggressiveness 0. Lit pixels match exactly; under
 * -ffast-math each composited blend may move a channel by one 16-bit step.
 */
inline void test_distorted_ring_stack_matches_sequential() {
  constexpr int W = 96, H = 64;
  constexpr int N_RINGS = 5, LUT_N = 32;
  const ScopedPoleLod lod(0.0f);

  const float ths[N_RINGS] = {0.07f, 0.05f, 0.09f, 0.05f, 0.07f};
  const Color4 colors[N_RINGS] = {Color4(Pixel(60000, 0, 5000), 0.9f),
                                  Color4(Pixel(5000, 0, 10000), 0.6f),
                                  Color4(Pixel(10000, 0, 60000), 0.4f),
                                  Color4(Pixel(30000, 0, 30000), 0.7f),
                                  Color4(Pixel(45000, 0, 20000), 1.0f)};

  auto run_case = [&](const math::Vector &normal, bool partial_clip,
                      int culled) {
    math::Basis basis = math::make_basis(
        math::make_rotation(math::Vector(0.2f, 0.5f, 0.8f).normalized(), 0.3f),
        normal);

    float knots[N_RINGS][LUT_N + 1];
    build_ring_knots<N_RINGS, LUT_N>(knots);
    // Ring i's centerline colatitude must be PI * (i + 1) / (N_RINGS + 1);
    // target_angle is radius * PI/2.
    auto ring_radius = [](int i) { return 2.0f * (i + 1) / (N_RINGS + 1); };

    alignas(SDF::DistortedRing) unsigned char
        mem[N_RINGS * sizeof(SDF::DistortedRing)];
    auto *shapes = reinterpret_cast<SDF::DistortedRing *>(mem);
    int8_t slot_by_ring[N_RINGS];
    Color4 slot_color[N_RINGS];
    int n_slots = 0;
    for (int i = 0; i < N_RINGS; ++i) {
      if (i == culled) {
        slot_by_ring[i] = -1;
        continue;
      }
      new (&shapes[n_slots]) SDF::DistortedRing(basis, ring_radius(i), ths[i],
                                                knots[i], LUT_N, 0.0f, nullptr);
      slot_color[n_slots] = colors[i];
      slot_by_ring[i] = static_cast<int8_t>(n_slots);
      ++n_slots;
    }

    auto shade = [&](int s, Fragment &f) {
      const uint16_t g = static_cast<uint16_t>(math::wrap_t(f.v0) * 60000.0f);
      f.color = Color4(Pixel(slot_color[s].color.r, g, slot_color[s].color.b),
                       slot_color[s].alpha * f.v2);
    };

    std::vector<Pixel> expected(W * H);
    {
      hs_test::StubEffect seq(W, H);
      if (partial_clip) {
        seq.set_clip(9, 53, 17, 81);
        seq.set_margin(0);
      }
      Pipeline<W, H> pipeline;
      {
        Canvas canvas(seq);
        for (int i = 0; i < N_RINGS; ++i) {
          const int s = slot_by_ring[i];
          if (s < 0)
            continue;
          auto shader = [&](const math::Vector &, Fragment &f) { shade(s, f); };
          Scan::DistortedRing::draw<W, H>(pipeline, canvas, basis,
                                          ring_radius(i), ths[i], knots[i],
                                          LUT_N, shader, 0.0f, false,
                                          /*suppress_pole_fill=*/true);
        }
      }
      seq.advance_display();
      capture_frame<W, H>(seq, expected);
    }

    hs_test::StubEffect fused(W, H);
    if (partial_clip) {
      fused.set_clip(9, 53, 17, 81);
      fused.set_margin(0);
    }
    Pipeline<W, H> pipeline;
    Scan::DistortedRingStack::CandidateTable<W, H> table;
    {
      Canvas canvas(fused);
      Scan::DistortedRingStack::draw<W, H>(
          pipeline, canvas, N_RINGS, shapes, slot_by_ring, n_slots, table,
          [&](int s, const math::Vector &, Fragment &f) { shade(s, f); });
    }
    fused.advance_display();

    // ScalarFn's inplace_function member is not trivially destructible.
    for (int s = 0; s < n_slots; ++s)
      shapes[s].~DistortedRing();

#if defined(HS_TEST_FAST_MATH)
    constexpr int CHANNEL_TOL = N_RINGS; // one 16-bit step per composited blend
#else
    constexpr int CHANNEL_TOL = 0;
#endif
    size_t lit = 0;
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const Pixel &a = expected[y * W + x];
        const Pixel &b = fused.get_pixel(x, y);
        if (!is_black(a))
          ++lit;
        HS_EXPECT_EQ(is_black(a), is_black(b));
        HS_EXPECT_NEAR(static_cast<int>(a.r), static_cast<int>(b.r),
                       CHANNEL_TOL);
        HS_EXPECT_NEAR(static_cast<int>(a.g), static_cast<int>(b.g),
                       CHANNEL_TOL);
        HS_EXPECT_NEAR(static_cast<int>(a.b), static_cast<int>(b.b),
                       CHANNEL_TOL);
      }
    }
    // Guard against both paths drawing nothing.
    HS_EXPECT_GT(lit, (size_t)200);
  };

  const math::Vector axis = math::Vector(0.3f, 0.8f, -0.5f).normalized();
  run_case(axis, false, -1);
  run_case(axis, true, -1);
  run_case(axis, false, 2);
  // Near-pole axis: r_val under the horizontal-projection floor forces the
  // full-row-scan fallback on both paths.
  run_case(math::Vector(0.005f, 1.0f, 0.0f).normalized(), false, -1);
}

/** @brief Empty clips skip distorted ring candidate tables. */
inline void test_distorted_ring_stack_empty_clip_skips_table() {
  constexpr int W = 64, H = 48;
  hs_test::StubEffect effect(W, H);
  effect.set_clip(30, 30, 0, W);
  effect.set_margin(0);
  Pipeline<W, H> pipeline;
  const auto basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  float knots[9]{};
  SDF::DistortedRing ring(basis, 1.0f, 0.18f, knots, 8, 0.0f, nullptr);
  const int8_t slots[] = {0};
  Scan::DistortedRingStack::CandidateTable<W, H> table;
  for (auto &cell : table.cells)
    cell = {17, 23};
  int shaded = 0;
  {
    Canvas canvas(effect);
    Scan::DistortedRingStack::draw<W, H>(
        pipeline, canvas, 1, &ring, slots, 1, table,
        [&](int, const math::Vector &, Fragment &) { ++shaded; });
  }
  HS_EXPECT_EQ(shaded, 0);
  for (const auto &cell : table.cells) {
    HS_EXPECT_EQ(cell.lo, 17);
    HS_EXPECT_EQ(cell.hi, 23);
  }
}

/**
 * @brief Verifies the fused RingGroup and DistortedRingStack walks ignore
 *        Render::pole_lod_aggressiveness.
 * @details Frames at aggressiveness 0 and 4 are bit-identical.
 */
inline void test_fused_walks_ignore_pole_lod() {
  constexpr int W = 96, H = 64;
  constexpr int N = 4, LUT_N = 32;
  const ScopedPoleLod scoped_lod(0.0f);
  const math::Vector normal = math::Vector(0.3f, 0.8f, -0.5f).normalized();
  const math::Basis basis = math::make_basis(
      math::make_rotation(math::Vector(0.2f, 0.5f, 0.8f).normalized(), 0.3f),
      normal);
  const float ths[N] = {0.08f, 0.04f, 0.06f, 0.05f};
  const Color4 colors[N] = {Color4(Pixel(60000, 10000, 5000), 0.9f),
                            Color4(Pixel(5000, 60000, 10000), 0.6f),
                            Color4(Pixel(10000, 5000, 60000), 0.4f),
                            Color4(Pixel(30000, 30000, 30000), 0.7f)};
  auto shader = [&](int s, const math::Vector &, Fragment &f) {
    f.color = colors[s];
  };

  auto draw_group = [&](std::vector<Pixel> &out) {
    SDF::Ring shapes[N];
    for (int s = 0; s < N; ++s) {
      shapes[s] = SDF::Ring(
          math::make_basis(
              math::make_rotation(math::Vector(0.2f, 0.5f, 0.8f).normalized(),
                                  0.02f * s),
              normal),
          1.0f, ths[s]);
    }
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      Scan::RingGroup::draw<W, H>(pipeline, canvas, shapes, N, shader);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, out);
  };

  auto draw_stack = [&](std::vector<Pixel> &out) {
    float knots[N][LUT_N + 1];
    build_ring_knots<N, LUT_N>(knots);
    alignas(
        SDF::DistortedRing) unsigned char mem[N * sizeof(SDF::DistortedRing)];
    auto *shapes = reinterpret_cast<SDF::DistortedRing *>(mem);
    int8_t slot_by_ring[N];
    for (int i = 0; i < N; ++i) {
      new (&shapes[i])
          SDF::DistortedRing(basis, 2.0f * (i + 1) / (N + 1), ths[i], knots[i],
                             LUT_N, 0.0f, nullptr);
      slot_by_ring[i] = static_cast<int8_t>(i);
    }
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    Scan::DistortedRingStack::CandidateTable<W, H> table;
    {
      Canvas canvas(fx);
      Scan::DistortedRingStack::draw<W, H>(pipeline, canvas, N, shapes,
                                           slot_by_ring, N, table, shader);
    }
    fx.advance_display();
    for (int i = 0; i < N; ++i)
      shapes[i].~DistortedRing();
    capture_frame<W, H>(fx, out);
  };

  auto expect_knob_independent = [&](const char *label, auto &&draw) {
    HS_CONTEXT(label);
    std::vector<Pixel> undecimated, decimated;
    Render::pole_lod_aggressiveness = 0.0f;
    draw(undecimated);
    Render::pole_lod_aggressiveness = 4.0f;
    HS_EXPECT_GT(Scan::pole_lod_run(1.0f), 1);
    draw(decimated);
    size_t lit = 0;
    for (size_t i = 0; i < undecimated.size(); ++i) {
      if (!is_black(undecimated[i]))
        ++lit;
      HS_EXPECT_EQ(undecimated[i].r, decimated[i].r);
      HS_EXPECT_EQ(undecimated[i].g, decimated[i].g);
      HS_EXPECT_EQ(undecimated[i].b, decimated[i].b);
    }
    HS_EXPECT_GT(lit, (size_t)200);
  };
  expect_knob_independent("ring group", draw_group);
  expect_knob_independent("distorted ring stack", draw_stack);
}
