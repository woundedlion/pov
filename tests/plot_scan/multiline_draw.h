/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

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
 * @details The seam's midpoint separates the two renders: it is off every open
 * edge and on the closed one.
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
