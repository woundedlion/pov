/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Accepts points in the row band whose owner or either neighbour is unclipped. */
template <int W>
inline bool reaches_band(const math::Vector &p,
                         const ClipRegion::XClip &columns, float top,
                         float bottom) {
  const float column = math::vector_to_theta<W>(p);
  const int owner = math::fast_wrap(int(floorf(column + 0.5f)), W);
  return p.y <= top && p.y >= bottom &&
         (!columns.active || !columns.clipped(owner) ||
          !columns.clipped(math::fast_wrap(owner + 1, W)) ||
          !columns.clipped(math::fast_wrap(owner - 1, W)));
}

/** @brief The swept edge cull never rejects an edge whose arcs over the frame
 * reach the draw's band. */
inline void test_column_edge_cull_is_conservative() {
  constexpr int W = 288, H = 144;
  ScratchScope scope(scratch_arena_a);
  Spacetime::Worldlines<> motion;
  motion.bind(scratch_arena_a, 2);
  uint32_t rejected = 0, checked = 0;
  hs::random().seed(20261021u);
  for (const auto band :
       {std::array<int, 4>{0, 72, 0, 144}, std::array<int, 4>{0, 72, 144, 288},
        std::array<int, 4>{36, 72, 100, 136}, std::array<int, 4>{0, H, 30, 230},
        std::array<int, 4>{0, 30, 0, W}, std::array<int, 4>{110, H, 0, 144}}) {
    ClipRegion clip;
    clip.w = W;
    clip.h = H;
    clip.y_start = band[0];
    clip.y_end = band[1];
    clip.x_start = band[2];
    clip.x_end = band[3];
    const float top = cosf(
        fmaxf(0.0f, math::DisplayGeometry<H>::row_to_phi(float(band[0] - 2))));
    const float bottom = cosf(
        math::DisplayGeometry<H>::row_to_phi(float(std::min(H, band[1] + 2))));
    const auto columns = clip.x_clip();
    const Plot::SweptEdgeCull<W> cull(columns, top, bottom);
    for (float radians : {0.0f, 0.01f, 0.06f}) {
      for (int trial = 0; trial < 4000; ++trial) {
        const auto a = hs_test::rand_unit();
        const auto edge_axis =
            math::cross(a, hs_test::rand_unit()).normalized();
        const float edge_angle = hs::rand_f(0.0f, 0.7f);
        const auto b =
            math::rotate(a, math::make_rotation(edge_axis, edge_angle));
        const auto a_axis = hs_test::rand_unit();
        const float a_angle = hs::rand_f(0.0f, radians);
        const auto a1 = math::rotate(a, math::make_rotation(a_axis, a_angle));
        const auto b_axis = hs_test::rand_unit();
        const float b_angle = hs::rand_f(0.0f, radians);
        const auto b1 = math::rotate(b, math::make_rotation(b_axis, b_angle));
        if (trial % 100 == 0) {
          motion.begin_frame(2, false);
          motion.capture_end(0, a);
          motion.capture_end(1, b);
          motion.begin_frame(2, true);
          motion.capture_end(0, a1);
          motion.capture_end(1, b1);
          Plot::EdgeWorldsheet reference, gated;
          reference.prepare<W, H>(motion, 0, 1, nullptr, 0.0f, false);
          const bool HIDDEN =
              cull.rejects({reference.bases[0].v1, reference.bases[1].v1,
                            reference.ends[0], reference.ends[1]},
                           reference.deviation + 1e-4f);
          Plot::g_geodesic_edge_span_builds = 0;
          const bool PREPARED = gated.prepare<W, H>(motion, 0, 1, nullptr, 0.0f,
                                                    false, true, &cull);
          HS_EXPECT_EQ(PREPARED, !HIDDEN);
          HS_EXPECT_EQ(Plot::g_geodesic_edge_span_builds, HIDDEN ? 0u : 2u);
        }
        const bool rejects = cull.rejects({a, a1, b, b1}, 1e-4f);
        rejected += rejects;
        bool reaches = false;
        for (int move = 0; move <= 8 && !reaches; ++move) {
          const float u = float(move) / 8.0f;
          const auto ma = math::nlerp_unit(a, a1, u);
          const auto mb = math::nlerp_unit(b, b1, u);
          for (int k = 0; k <= 32 && !reaches; ++k) {
            const auto p = math::slerp(ma, mb, float(k) / 32.0f);
            reaches = reaches_band<W>(p, columns, top, bottom);
          }
        }
        checked += reaches;
        if (reaches)
          HS_EXPECT(!rejects, "a reachable edge is kept");
      }
    }
  }
  std::printf("  [column edge cull] rejected %u, reachable %u\n", rejected,
              checked);
  HS_EXPECT_GT(rejected, 1000u);
  HS_EXPECT_GT(checked, 20000u);
}

/** @brief A canvas's column cull keeps every padded polyline, or the padded
 *         cap around it, that can reach the walk's rows and band. */
inline void test_column_edge_cull_keeps_padded_point_sets() {
  const Plot::SweptEdgeCull<288> NORTHERN_ROWS({}, 1.0f, 0.25f);
  const math::Vector ANTIPODAL[] = {math::X_AXIS, -math::X_AXIS};
  HS_EXPECT_FALSE(NORTHERN_ROWS.rejects(ANTIPODAL, 2, 0.0f));
  const math::Vector SOUTHERN[] = {
      math::Vector(1.0f, -1.0f, 0.0f).normalized(),
      math::Vector(1.0f, -1.0f, 0.1f).normalized()};
  HS_EXPECT_TRUE(NORTHERN_ROWS.rejects(SOUTHERN, 2, 0.0f));
  const auto poles = []<int HEIGHT>() {
    constexpr int WIDTH = 288;
    hs_test::StubEffect effect(WIDTH, HEIGHT);
    Canvas canvas(effect);
    const auto cull = Plot::SweptEdgeCull<WIDTH>::of<HEIGHT>(canvas);
    const math::Vector north{0.0f, 1.0f, 0.0f};
    const math::Vector south{0.0f, -1.0f, 0.0f};
    if constexpr (HEIGHT == 20 || HS_DISPLAY_PROFILE == 0) {
      HS_EXPECT_EQ(Plot::SweptEdgeCull<WIDTH>::top_of<HEIGHT>(canvas), 1.0f);
      HS_EXPECT_EQ(Plot::SweptEdgeCull<WIDTH>::bottom_of<HEIGHT>(canvas),
                   -1.0f);
      HS_EXPECT_FALSE(cull.rejects(&north, 1, 0.0f));
      HS_EXPECT_FALSE(cull.rejects(&south, 1, 0.0f));
    }
  };
  poles.template operator()<20>();
  poles.template operator()<144>();
  constexpr int W = 288, H = 144;
  uint32_t rejected = 0, cap_rejected = 0, reachable = 0;
  hs::random().seed(20261022u);
  for (const auto band :
       {std::array<int, 4>{36, 72, 0, 144}, std::array<int, 4>{0, 72, 144, 288},
        std::array<int, 4>{110, H, 30, 230}}) {
    hs_test::StubEffect effect(W, H);
    effect.set_clip(band[0], band[1], band[2], band[3]);
    Canvas canvas(effect);
    const auto cull = Plot::SweptEdgeCull<W>::of<H>(canvas);
    const float top = Plot::SweptEdgeCull<W>::top_of<H>(canvas);
    const float bottom = Plot::SweptEdgeCull<W>::bottom_of<H>(canvas);
    HS_EXPECT_NEAR(top,
                   cosf(fmaxf(0.0f, math::DisplayGeometry<H>::row_to_phi(
                                        float(band[0] - 2)))),
                   1e-6f);
    const auto &columns = canvas.render_x_clip();
    for (float pad : {0.0f, 0.02f, 0.06f}) {
      for (int trial = 0; trial < 1500; ++trial) {
        const int count = 2 + trial % 8;
        math::Vector points[9];
        points[0] = hs_test::rand_unit();
        for (int i = 1; i < count; ++i) {
          const auto axis =
              math::cross(points[i - 1], hs_test::rand_unit()).normalized();
          const float angle = hs::rand_f(0.0f, 0.08f);
          points[i] =
              math::rotate(points[i - 1], math::make_rotation(axis, angle));
        }
        float reach = 1.0f;
        for (int i = 0; i < count; ++i)
          reach = fminf(reach, math::dot(points[count / 2], points[i]));
        const bool rejects = cull.rejects(points, count, pad);
        const bool cap_rejects =
            cull.rejects(&points[count / 2], 1,
                         pad + sqrtf(fmaxf(0.0f, 2.0f - 2.0f * reach)));
        rejected += rejects;
        cap_rejected += cap_rejects;
        bool reaches = false;
        for (int i = 0; i + 1 < count && !reaches; ++i)
          for (int k = 0; k <= 16 && !reaches; ++k) {
            const math::Vector at =
                math::slerp(points[i], points[i + 1], float(k) / 16.0f);
            const math::Vector west =
                math::cross(math::Y_AXIS, at).magnitude() > 1e-4f
                    ? math::cross(math::Y_AXIS, at).normalized()
                    : math::X_AXIS;
            const math::Vector north = math::cross(at, west);
            for (int d = 0; d <= 8 && !reaches; ++d) {
              const float angle = math::TWO_PI_F * float(d) / 8.0f;
              const math::Vector p =
                  d == 8
                      ? at
                      : (at + (west * cosf(angle) + north * sinf(angle)) * pad)
                            .normalized();
              reaches = reaches_band<W>(p, columns, top, bottom);
            }
          }
        reachable += reaches;
        if (reaches) {
          HS_EXPECT(!rejects, "a reachable polyline is kept");
          HS_EXPECT(!cap_rejects, "a reachable polyline's cap is kept");
        }
      }
    }
  }
  std::printf("  [column point-set cull] rejected %u, cap %u, reachable %u\n",
              rejected, cap_rejected, reachable);
  HS_EXPECT_GT(rejected, 1000u);
  HS_EXPECT_GT(cap_rejected, 1000u);
  HS_EXPECT_GT(reachable, 2800u);
}

/** @brief Column culls retain points inside widened filter margins. */
inline void test_column_edge_cull_honors_render_margin() {
  constexpr int W = 288, H = 144;
  for (const int margin : {1, 3, 5}) {
    hs_test::StubEffect effect(W, H);
    effect.set_clip(40, 60, 0, W);
    effect.set_margin(margin);
    Canvas canvas(effect);
    const auto cull = Plot::SweptEdgeCull<W>::of<H>(canvas);
    for (const float row : {37.0f, 63.0f}) {
      const float phi = math::DisplayGeometry<H>::row_to_phi(row);
      const math::Vector point(sinf(phi), cosf(phi), 0.0f);
      HS_EXPECT_EQ(cull.rejects(&point, 1, 0.0f), margin == 1);
    }
    for (const float row : {32.0f, 68.0f}) {
      const float phi = math::DisplayGeometry<H>::row_to_phi(row);
      const math::Vector point(sinf(phi), cosf(phi), 0.0f);
      HS_EXPECT_TRUE(cull.rejects(&point, 1, 0.0f));
    }
  }
}

/** @brief Physical northern AA support keeps visible points and swept arcs. */
inline void test_column_physical_north_cap() {
  constexpr int W = 288, H = 144;
  constexpr math::LatitudeGeometry GEOMETRY(H, 0.02f * math::PI_F,
                                            0.98f * math::PI_F);
  const float TOP =
      Plot::SweptEdgeCull<W>::physical_top(GEOMETRY.row_to_phi(-1.0f));
  const auto point = [&](float row, float longitude) {
    const float phi = GEOMETRY.row_to_phi(row);
    return math::Vector(sinf(phi) * cosf(longitude), cosf(phi),
                        sinf(phi) * sinf(longitude));
  };
  const auto reaches = [&](const math::Vector &a, const math::Vector &b,
                           const math::Vector &a1, const math::Vector &b1) {
    Plot::EdgeWorldsheet frame{};
    frame.bases[0].v1 = a;
    frame.bases[1].v1 = a1;
    frame.ends[0] = b;
    frame.ends[1] = b1;
    return frame.reaches_rows<W>(TOP, -1.0f);
  };
  for (const float row : {-0.5f, -2.5f}) {
    const auto a = point(row, 0.0f);
    const auto b = point(row, 0.02f);
    const auto taps = Filter::Screen::splat_taps<W, H>(
        40.25f, GEOMETRY.phi_to_row(math::fast_acos(a.y)));
    const bool visible = taps.y0_physical || taps.y1_physical;
    HS_CONTEXT("physical row tenths", int(row * 10.0f));
    HS_EXPECT_EQ(visible, row == -0.5f);
    HS_EXPECT_EQ(reaches(a, a, a, a), visible);
    HS_EXPECT_EQ(reaches(a, b, a, b), visible);
  }
  uint32_t visible = 0, rejected = 0;
  for (int history = 0; history < 160; ++history) {
    const float row = -2.7f + 0.03f * float(history);
    const auto a = point(row, 0.0f);
    const auto b = point(row + 0.1f, 0.04f);
    const auto a1 = point(row + 0.2f, 0.03f);
    const auto b1 = point(row - 0.1f, 0.07f);
    const bool kept = reaches(a, b, a1, b1);
    rejected += !kept;
    for (int time = 0; time <= 8; ++time) {
      const float t = float(time) / 8.0f;
      const auto start = math::nlerp_unit(a, a1, t);
      const auto end = math::nlerp_unit(b, b1, t);
      for (int sample = 0; sample <= 16; ++sample) {
        const auto p = math::slerp(start, end, float(sample) / 16.0f);
        const float y = GEOMETRY.phi_to_row(math::fast_acos(p.y));
        const auto taps = Filter::Screen::splat_taps<W, H>(40.25f, y);
        if (taps.y0_physical || taps.y1_physical) {
          ++visible;
          HS_CONTEXT("physical cap history", history);
          HS_EXPECT(kept, "a swept arc with physical AA support is kept");
        }
      }
    }
  }
  HS_EXPECT_GT(visible, 10000u);
  HS_EXPECT_GT(rejected, 20u);
}
