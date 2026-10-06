/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// Plot::Multiline::draw — sampling plus rasterization, held against an
// analytic point-to-geodesic-arc oracle.
// ============================================================================

/** @brief Four non-coplanar control directions used by the Multiline cases. */
inline void multiline_control_points(std::vector<math::Vector> &out) {
  out = {math::Vector(0.95f, -0.1f, 0.3f).normalized(),
         math::Vector(0.3f, 0.9f, 0.2f).normalized(),
         math::Vector(-0.5f, 0.2f, 0.84f).normalized(),
         math::Vector(-0.2f, -0.85f, 0.49f).normalized()};
}

/** @brief Rejects a great-circle foot on the complementary arc. */
inline void test_arc_angular_distance_clamps_to_minor_arc() {
  const auto equator = [](float degrees) {
    const float angle = degrees * math::PI_F / 180.0f;
    return math::Vector(cosf(angle), 0.0f, sinf(angle));
  };
  HS_EXPECT_NEAR(
      arc_angular_distance(equator(270.0f), equator(0.0f), equator(150.0f)),
      math::PI_F * 0.5f, 1e-5f);
}

/**
 * @brief Verifies Multiline::draw paints its geodesic edges and nothing else.
 * @details Every plotted position must lie on one of the polyline's geodesic
 * arcs, every control point must be reached, and the walk must be gap-free.
 * The tolerance is one screen row.
 */
inline void test_multiline_draw_covers_only_its_geodesic_edges() {
  constexpr int W = 128, H = 64;
  const float row = math::PI_F / (H - 1);

  std::vector<math::Vector> control;
  multiline_control_points(control);

  ScratchScope sc(plot_arena());
  Fragments verts;
  verts.bind(plot_arena(), control.size());
  for (const math::Vector &v : control) {
    Fragment f;
    f.pos = v;
    verts.push_back(f);
  }

  hs_test::StubEffect fx(W, H);
  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::Multiline::draw<W, H>(pipe, c, verts, noop_shader);
  }
  fx.advance_display();

  HS_EXPECT_GT(pipe.plotted.size(), (size_t)0);
  float worst_off_path = 0.0f;
  for (const math::Vector &p : pipe.plotted) {
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
    float nearest = math::PI_F;
    for (size_t i = 1; i < control.size(); ++i)
      nearest = std::min(nearest,
                         arc_angular_distance(p, control[i - 1], control[i]));
    worst_off_path = std::max(worst_off_path, nearest);
  }
  HS_EXPECT_LT(worst_off_path, row);

  // Every control point is reached, so no edge was skipped.
  for (const math::Vector &c : control) {
    float nearest = math::PI_F;
    for (const math::Vector &p : pipe.plotted)
      nearest = std::min(nearest, math::angle_between(p, c));
    HS_EXPECT_LT(nearest, row);
  }

  // Consecutive samples in the open walk stay within a row.
  HS_EXPECT_LT(max_consecutive_gap(pipe.plotted, /*wrap=*/false), row);
}

/**
 * @brief Verifies the closed flag draws the last->first seam edge.
 * @details Closing routes a loop_seam fragment through draw_fragments. The
 * seam's own midpoint separates the two renders: it is off every open edge and
 * on the closed one.
 */
inline void test_multiline_draw_closed_adds_the_seam_edge() {
  constexpr int W = 128, H = 64;
  const float row = math::PI_F / (H - 1);

  std::vector<math::Vector> control;
  multiline_control_points(control);
  const math::Vector seam_midpoint =
      (control.back() + control.front()).normalized();

  // The midpoint is only a seam witness if it is off every drawn open edge.
  float midpoint_to_open_path = math::PI_F;
  for (size_t i = 1; i < control.size(); ++i)
    midpoint_to_open_path = std::min(
        midpoint_to_open_path,
        arc_angular_distance(seam_midpoint, control[i - 1], control[i]));
  HS_EXPECT_GT(midpoint_to_open_path, 4.0f * row);

  auto draw_path = [&](bool closed) {
    ScratchScope sc(plot_arena());
    Fragments verts;
    verts.bind(plot_arena(), control.size());
    for (const math::Vector &v : control) {
      Fragment f;
      f.pos = v;
      verts.push_back(f);
    }
    hs_test::StubEffect fx(W, H);
    CapturePipeline pipe;
    {
      Canvas c(fx);
      Plot::Multiline::draw<W, H>(pipe, c, verts, noop_shader, closed);
    }
    fx.advance_display();
    float nearest = math::PI_F;
    for (const math::Vector &p : pipe.plotted)
      nearest = std::min(nearest, math::angle_between(p, seam_midpoint));
    return std::pair<size_t, float>{pipe.plotted.size(), nearest};
  };

  const auto open_path = draw_path(false);
  const auto closed_path = draw_path(true);

  HS_EXPECT_GT(closed_path.first, open_path.first);
  HS_EXPECT_LT(closed_path.second, row);
  HS_EXPECT_GT(open_path.second, 4.0f * row);
}

/** @brief Replay parameters stay within the interpolation domain at antipodes. */
inline void test_plot_line_antipodal_replay_parameter() {
  constexpr int W = 288, H = 144;
  for (const auto &start : {math::Y_AXIS, math::X_AXIS,
                            math::Vector(1.0f, 1.0f, 1.0f).normalized()}) {
    hs_test::StubEffect fx(W, H);
    CapturePipeline pipe;
    Fragment from, to;
    from.pos = start;
    to.pos = start * -1.0f;
    from.v0 = 0.0f;
    to.v0 = 1.0f;
    int samples = 0;
    {
      Canvas canvas(fx);
      Plot::Line::draw<W, H>(pipe, canvas, from, to,
                             [&](const math::Vector &, Fragment &f) {
                               HS_EXPECT_GE(f.v0, 0.0f);
                               HS_EXPECT_LE(f.v0, 1.0f);
                               ++samples;
                               f.color = Color4(Pixel(65535, 0, 0), 1.0f);
                             });
    }
    fx.advance_display();
    HS_EXPECT_GT(samples, 100);
  }
}

/**
 * @brief Verifies a geodesic line through the north pole plots the pole row.
 * @details Interpolated points are up to 1.7e-3 non-unit, and acos's infinite
 * slope at y=1 turns that into a multi-row shift unless the drawing phase
 * re-normalizes with newton_unit().
 */
inline void test_plot_line_over_pole_reaches_row0() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe; // bare sink (no AA) so we see raw sample placement

  // Geodesic from 0.4 rad down the +Z side of the N pole to 0.4 rad down the
  // -Z side; its midpoint is the pole (row 0).
  Fragment f1, f2;
  f1.pos = math::Vector(0.0f, cosf(0.4f), sinf(0.4f));
  f2.pos = math::Vector(0.0f, cosf(0.4f), -sinf(0.4f));
  {
    Canvas c(fx);
    Plot::Line::draw<W, H>(pipe, c, f1, f2,
                           [](const math::Vector &, Fragment &f) {
                             f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
                           });
  }
  fx.advance_display();

  size_t row0 = 0;
  for (int x = 0; x < W; ++x)
    if (!is_black(fx.get_pixel(x, 0)))
      ++row0;
  HS_EXPECT_GT(row0, (size_t)0);
}
