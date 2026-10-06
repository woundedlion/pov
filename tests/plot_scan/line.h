/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Plot::Line::sample
// ============================================================================

/**
 * @brief Verifies Line::sample emits density+1 great-circle fragments with exact
 *        endpoints, unit-length positions, v0 progress spanning 0..1, and v1 arc
 *        length ending at the segment's subtended angle.
 */
inline void test_line_sample_endpoints_and_unit_length() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = math::Vector(0.6f, 0.0f, 0.8f);
  b.pos = math::Vector(-0.64f, 0.6f, 0.48f);
  const int density = 8;
  Plot::Line::sample(points, a, b, density);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(density + 1));

  HS_EXPECT_EQ(points[0].pos.x, a.pos.x);
  HS_EXPECT_EQ(points[0].pos.y, a.pos.y);
  HS_EXPECT_EQ(points[0].pos.z, a.pos.z);
  HS_EXPECT_EQ(points[density].pos.x, b.pos.x);
  HS_EXPECT_EQ(points[density].pos.y, b.pos.y);
  HS_EXPECT_EQ(points[density].pos.z, b.pos.z);

  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }

  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points[density].v0, 1.0f, 1e-6f);
  float total_angle = math::angle_between(a.pos, b.pos);
  HS_EXPECT_NEAR(points[density].v1, total_angle, 1e-4f);
  HS_EXPECT_NEAR(total_angle, math::PI_F * 0.5f, 1e-4f);
}

/** Angular slack on a Line::sample position: fast_sincosf_0_pi error leaves
 * about 2.4e-3 rad of directional error after renormalization. */
constexpr float LINE_SAMPLE_ANGLE_TOL = 4e-3f;

/**
 * @brief Verifies interior Line::sample fragments lie on the minor arc itself,
 *        at even angular spacing, not merely inside a cone bounding it.
 * @details A point is on the minor arc iff its angles to the two endpoints sum
 * to the whole span. Each sample's angle from the start must match its share of
 * the span, and the arc-length register must agree with the geometry.
 */
inline void test_line_sample_interior_between_endpoints() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = math::Vector(1, 0, 0);
  b.pos = math::Vector(0, 1, 0);
  const int density = 4;
  Plot::Line::sample(points, a, b, density);
  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(density + 1));

  const float total = math::angle_between(a.pos, b.pos);
  for (size_t i = 1; i + 1 < points.size(); ++i) {
    HS_CONTEXT("sample", (long long)i);
    const math::Vector &p = points[i].pos;
    const float from_a = math::angle_between(a.pos, p);
    const float from_b = math::angle_between(b.pos, p);
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
    HS_EXPECT_NEAR(from_a + from_b, total, LINE_SAMPLE_ANGLE_TOL);
    HS_EXPECT_NEAR(from_a, total * static_cast<float>(i) / density,
                   LINE_SAMPLE_ANGLE_TOL);
    HS_EXPECT_NEAR(points[i].v0, static_cast<float>(i) / density, 1e-6f);
    HS_EXPECT_NEAR(points[i].v1, from_a, LINE_SAMPLE_ANGLE_TOL);
  }
}

/**
 * @brief Verifies a zero-length segment (coincident endpoints) emits a dot: two
 *        fragments sitting exactly on the shared endpoint, whatever density was
 *        asked for, with the registers pinned to the start of the arc.
 */
inline void test_line_sample_degenerate_segment() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = math::Vector(0, 0, 1);
  b.pos = math::Vector(0, 0, 1);

  Plot::Line::sample(points, a, b, 8);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)2);
  for (size_t i = 0; i < points.size(); ++i) {
    HS_CONTEXT("sample", (long long)i);
    HS_EXPECT_EQ(points[i].pos.x, a.pos.x);
    HS_EXPECT_EQ(points[i].pos.y, a.pos.y);
    HS_EXPECT_EQ(points[i].pos.z, a.pos.z);
    HS_EXPECT_EQ(points[i].v0, 0.0f);
    HS_EXPECT_EQ(points[i].v1, 0.0f);
    HS_EXPECT_EQ(points[i].v2, 0.0f);
  }
}

/**
 * @brief Verifies Line::sample picks a stable perpendicular axis for antipodal
 *        endpoints so the arc stays finite, unit-length, and passes through a
 *        real ~90deg midpoint.
 * @details Antipodal endpoints (angle == pi) make cross(a, b) == 0, so the
 *          rotation axis needs a perpendicular fallback.
 */
inline void test_line_sample_antipodal_stable_axis() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = math::Vector(1, 0, 0);
  b.pos = math::Vector(-1, 0, 0);
  const int density = 8;
  Plot::Line::sample(points, a, b, density);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(density + 1));
  HS_EXPECT_NEAR(points[0].pos.x, a.pos.x, 1e-6f);
  HS_EXPECT_NEAR(points[density].pos.x, b.pos.x, 1e-6f);

  for (size_t i = 0; i < points.size(); ++i) {
    const math::Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }

  // The midpoint is ~90deg from each antipodal endpoint.
  const math::Vector &mid = points[density / 2].pos;
  HS_EXPECT_NEAR(math::angle_between(a.pos, mid), math::PI_F * 0.5f, 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(b.pos, mid), math::PI_F * 0.5f, 1e-3f);
}

/**
 * @brief Antipodal endpoints one ULP apart in length still pick a stable axis
 *        through make_geodesic_edge_span.
 * @details acos' derivative diverges at ±1, so a single-ULP perturbation of the
 *          normalized dot moves angle_between off π by ~5e-4 rad while
 *          cross(a, b) stays exactly zero. The miss comes from correctly-rounded
 *          mul/div/sqrt alone, so it holds on every target.
 */
inline void test_line_sample_near_antipodal_ulp_stable_axis() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = math::Vector(-0.28f, 0.96f, 0.0f);
  b.pos = a.pos * -(1.0f + 0x1p-23f); // one ULP longer than -a

  // The arc pole cannot come from the cross product: it is unnormalizable.
  const math::Vector pole = math::cross(a.pos, b.pos);
  HS_EXPECT_LT(math::dot(pole, pole), math::EPS_NORMALIZE_SQ);

  const int density = 8;
  Plot::Line::sample(points, a, b, density);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)(density + 1));
  for (size_t i = 0; i < points.size(); ++i) {
    const math::Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }
  const math::Vector &mid = points[density / 2].pos;
  HS_EXPECT_NEAR(math::angle_between(a.pos, mid), math::PI_F * 0.5f, 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(b.pos, mid), math::PI_F * 0.5f, 1e-3f);

  // The span setup resolves the same edge to a unit arc pole perpendicular to a.
  const Plot::GeodesicEdgeSpan es = Plot::make_geodesic_edge_span(a.pos, b.pos);
  HS_EXPECT_TRUE(es.have_axis);
  HS_EXPECT_TRUE(es.antipodal);
  HS_EXPECT_NEAR(es.axis.length(), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(math::dot(es.axis, a.pos.normalized()), 0.0f, 1e-5f);

  // And so does rasterize_geodesic_strategy, reached through the rasterizer.
  constexpr int W = 128, H = 64;
  hs_test::StubEffect fx(W, H);
  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader);
  }
  fx.advance_display();
  HS_EXPECT_GT(pipe.plotted.size(), (size_t)0);
  for (const math::Vector &p : pipe.plotted) {
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }
}
