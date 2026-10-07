/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Interval-cull conservativeness: the cull may over-visit but must never drop
// a pixel belonging to the shape.
// ============================================================================

/**
 * @brief Records every pixel scan_region visits for a shape, as rasterize() drives it.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Shape SDF shape type providing the cull and interval interface.
 * @param shape Shape whose culled coverage is being captured.
 * @param visited Output flag grid (W*H), set to 1 for each visited pixel.
 * @details Drives the full canvas with no clip.
 */
template <int W, int H, typename Shape>
inline void cull_visited(const Shape &shape, std::vector<uint8_t> &visited) {
  visited.assign(static_cast<size_t>(W) * H, 0);
  auto bounds = shape.template get_vertical_bounds<H>();
  int y_lo = std::max(0, bounds.y_min);
  int y_hi = std::min(H - 1, bounds.y_max);
  if (y_lo > y_hi)
    return;
  Scan::scan_region<W, H>(
      y_lo, y_hi,
      [&](int y, auto &&out) {
        return shape.template get_horizontal_intervals<W, H>(y, out);
      },
      [&](int wx, int y, const math::Vector &, int run) {
        for (int i = 0; i < run; ++i)
          if (wx + i >= 0 && wx + i < W && y >= 0 && y < H)
            visited[static_cast<size_t>(y) * W + wx + i] = 1;
        return run;
      });
}

/**
 * @brief Asserts no pixel clearly inside the shape is dropped by the cull.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Shape SDF shape type providing the cull and distance interface.
 * @param shape Shape under test.
 * @param label Caller-identifying label reported with each failing pixel.
 * @return Count of interior pixels (dist < -pixel_width) found.
 * @details Asserts at least one interior pixel.
 */
template <int W, int H, typename Shape>
inline int expect_cull_covers_interior(const Shape &shape, const char *label) {
  HS_CONTEXT(label);
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  std::vector<uint8_t> visited;
  cull_visited<W, H>(shape, visited);

  const float pixel_width = 2.0f * math::PI_F / W;
  int interior = 0;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      if (SDF::distance_of(shape, p).dist < -pixel_width) {
        ++interior;
        HS_CONTEXT("interior px", x, y);
        HS_EXPECT_TRUE(visited[static_cast<size_t>(y) * W + x]);
      }
    }
  }
  HS_EXPECT_GT(interior, 0);
  return interior;
}

/**
 * @brief Asserts no pixel the AA band would shade is dropped by the cull.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Shape SDF shape type providing the cull and distance interface.
 * @param shape Shape under test.
 * @param label Caller-identifying label reported with each failing pixel.
 * @return Count of paintable pixels (dist < pixel_width) found.
 * @details Every pixel under pixel_width carries non-zero coverage.
 */
template <int W, int H, typename Shape>
inline int expect_cull_covers_fringe(const Shape &shape, const char *label) {
  HS_CONTEXT(label);
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  std::vector<uint8_t> visited;
  cull_visited<W, H>(shape, visited);

  const float pixel_width = 2.0f * math::PI_F / W;
  int paintable = 0;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      if (SDF::distance_of(shape, p).dist < pixel_width) {
        ++paintable;
        HS_CONTEXT("paintable px", x, y);
        HS_EXPECT_TRUE(visited[static_cast<size_t>(y) * W + x]);
      }
    }
  }
  HS_EXPECT_GT(paintable, 0);
  return paintable;
}

/**
 * @brief Verifies the Star / PlanarPolygon / SphericalPolygon cull covers the
 *   whole AA fringe.
 * @details All three read a shallow radial gradient at the tips, which the
 *   circumscribed-disc clamp in distance() bounds.
 */
inline void test_star_polygon_cull_covers_aa_fringe() {
  constexpr int W = 96, H = 48;
  const math::Vector axes[] = {math::Vector(0, 1, 0), math::Vector(1, 0, 0),
                               math::Vector(0.3f, -0.8f, 0.5f)};
  int total = 0;
  for (const math::Vector &axis : axes) {
    math::Basis basis = math::make_basis(math::Quaternion(), axis);
    for (float radius : {0.03f, 0.3f, 0.9f}) {
      for (int sides : {5, 8}) {
        SDF::Star star(basis, radius, sides, 0.0f);
        total += expect_cull_covers_fringe<W, H>(star, "star");

        SDF::PlanarPolygon poly(basis, radius, sides, 0.0f);
        total += expect_cull_covers_fringe<W, H>(poly, "planar polygon");

        SDF::SphericalPolygon sphpoly(basis, radius, sides, 0.0f);
        total += expect_cull_covers_fringe<W, H>(sphpoly, "spherical polygon");
      }
    }
  }
  HS_EXPECT_GT(total, 1000);
}

/** @brief Annular angle bounds contain the reference angle before pixel rounding. */
inline void test_annular_angles_bound_reference() {
  static_assert(SDF::ANNULAR_ANGLE_PAD >= 5.1e-5f);
  for (int i = 0; i <= 1000; ++i) {
    const float cosine = 2.0f * i / 1000 - 1.0f;
    float angle_min, angle_max;
    HS_EXPECT_TRUE(SDF::annular_band_angles(cosine, cosine, 0, 0, 1, angle_min,
                                            angle_max));
    const float reference = std::acos(cosine);
    HS_EXPECT_LE(angle_min, reference);
    HS_EXPECT_GE(angle_max, reference);
  }
  float angle_min, angle_max;
  HS_EXPECT_FALSE(
      SDF::annular_band_angles(2, 3, 0, 0, 1, angle_min, angle_max));
}

/** @brief Verifies interval cull coverage over the sampled orientation/radius grid. */
inline void test_cull_covers_interior_over_orientation_grid() {
  constexpr int W = 96, H = 48;

  const math::Vector axes[] = {
      math::Vector(0, 1, 0),          math::Vector(0, -1, 0),
      math::Vector(1, 0, 0),          math::Vector(0, 0, 1),
      math::Vector(1, 1, 0.4f),       math::Vector(-0.5f, 0.7f, -0.6f),
      math::Vector(0.3f, -0.8f, 0.5f)};

  for (const math::Vector &axis : axes) {
    math::Basis basis = math::make_basis(math::Quaternion(), axis);
    for (float radius : {0.3f, 0.6f, 0.9f}) {
      SDF::Ring ring(basis, radius, /*thickness=*/0.25f);
      expect_cull_covers_interior<W, H>(ring, "ring");
      SDF::FlatDistortedRing flat_ring(basis, radius, 0.25f);
      expect_cull_covers_interior<W, H>(flat_ring, "flat distorted ring");
      expect_cull_covers_fringe<W, H>(flat_ring, "flat distorted ring");

      SDF::SphericalPolygon spoly(basis, radius, /*sides=*/5, 0.0f);
      expect_cull_covers_interior<W, H>(spoly, "spherical polygon");

      SDF::Star star(basis, radius, /*sides=*/5, 0.0f);
      expect_cull_covers_interior<W, H>(star, "star");

      SDF::PlanarPolygon ppoly(basis, /*radius=*/radius / (math::PI_F / 2.0f),
                               /*sides=*/6, 0.0f);
      expect_cull_covers_interior<W, H>(ppoly, "planar polygon");
      const SDF::Union combined(ppoly, star);
      expect_cull_covers_interior<W, H>(combined, "polygon-star union");
      expect_cull_covers_fringe<W, H>(combined, "polygon-star union");

      SDF::Flower flower(basis, radius, /*sides=*/5, 0.0f);
      expect_cull_covers_interior<W, H>(flower, "flower");
    }
  }
}

/**
 * @brief Verifies a pole-axis ring's row band skips the pole rows its stroke
 *        never reaches, while still covering every interior pixel.
 * @details A ring about +Y is the worst case for a loose latitude band: its
 *   axis has no horizontal projection, so every row inside the band is
 *   full-width scanned.
 */
inline void test_pole_axis_ring_bounds_skip_pole_rows() {
  constexpr int W = 96, H = 48;
  math::Basis basis = equator_basis();
  SDF::Ring ring(basis, /*radius=*/0.6f, /*thickness=*/0.25f);

  auto bounds = ring.get_vertical_bounds<H>();
  HS_EXPECT_GT(bounds.y_min, 0);
  HS_EXPECT_LT(bounds.y_max, H - 1);
  expect_cull_covers_interior<W, H>(ring, "pole-axis ring");
}

/** @brief Linearized ring bounds retain every row with visible stroke alpha. */
inline void test_linearized_ring_bounds_cover_visible_rows() {
  constexpr int W = 288, H = 144;
  const math::Basis basis = equator_basis();
  int linearized = 0;
  for (int radius_index = 1; radius_index < 100; ++radius_index) {
    for (float thickness : {0.445f, 0.7f, 1.1f, 1.5f}) {
      const SDF::Ring ring(basis, radius_index * 0.02f, thickness);
      if (ring.inv_sin_target == 0.0f)
        continue;
      ++linearized;
      const auto bounds = ring.get_vertical_bounds<H>();
      for (int y = 0; y < H; ++y) {
        const auto point = math::pixel_to_vector<W, H>(0, y);
        if (ring.stroke_alpha(math::dot(point, basis.v)) > Scan::MIN_ALPHA)
          HS_EXPECT_TRUE(y >= bounds.y_min && y <= bounds.y_max);
      }
    }
  }
  HS_EXPECT_GT(linearized, 0);
}

/**
 * @brief Verifies the Intersection interval cull covers every interior pixel of
 *        a real leaf pair.
 * @details The last pose centers both polygons on +X, so their overlap
 *   straddles theta = 0 and each child can emit its span in a different wrap
 *   frame.
 */
inline void test_intersection_cull_covers_interior_over_polygon_pairs() {
  constexpr int W = 96, H = 48;
  using Poly = SDF::PlanarPolygon;

  struct Pose {
    math::Vector axis_a, axis_b;
    float radius_a, radius_b;
  };
  const Pose poses[] = {
      {math::Vector(0, 0, 1), math::Vector(0.3f, 0.2f, 1.0f), 0.8f, 0.5f},
      {math::Vector(0, 1, 0), math::Vector(0.25f, 1.0f, -0.15f), 0.7f, 0.45f},
      {math::Vector(-0.4f, 0.6f, 0.7f), math::Vector(-0.2f, 0.75f, 0.6f), 0.9f,
       0.6f},
      {math::Vector(1, 0, 0), math::Vector(1.0f, 0.15f, 0.2f), 0.8f, 0.5f},
  };

  for (const Pose &pose : poses) {
    math::Basis basis_a = math::make_basis(math::Quaternion(), pose.axis_a);
    math::Basis basis_b = math::make_basis(math::Quaternion(), pose.axis_b);
    Poly poly_a(basis_a, pose.radius_a / (math::PI_F / 2.0f), /*sides=*/6,
                0.0f);
    Poly poly_b(basis_b, pose.radius_b / (math::PI_F / 2.0f), /*sides=*/5,
                0.4f);

    SDF::Intersection<Poly, Poly> both(poly_a, poly_b);
    expect_cull_covers_interior<W, H>(both, "intersection leaf pair");
  }
}

/**
 * @brief Verifies the Subtract interval cull covers every interior pixel of a
 *        real leaf pair.
 * @details The last pose centers both on +X, so the minuend's spans straddle
 *   theta = 0.
 */
inline void test_subtract_cull_covers_interior_over_leaf_pairs() {
  constexpr int W = 96, H = 48;
  using Poly = SDF::PlanarPolygon;
  using Star = SDF::Star;

  struct Pose {
    math::Vector axis_a, axis_b;
    float radius_a; /**< Polygon circumradius in radians. */
    float radius_b; /**< Star tip radius, normalized (x PI/2 for radians). */
  };
  const Pose poses[] = {
      {math::Vector(0, 0, 1), math::Vector(0.2f, 0.1f, 1.0f), 0.9f, 0.35f},
      {math::Vector(0, 1, 0), math::Vector(0.2f, 1.0f, -0.1f), 0.8f, 0.3f},
      {math::Vector(-0.4f, 0.6f, 0.7f), math::Vector(-0.3f, 0.65f, 0.7f), 1.0f,
       0.4f},
      {math::Vector(1, 0, 0), math::Vector(1.0f, 0.1f, 0.15f), 0.9f, 0.35f},
  };

  for (const Pose &pose : poses) {
    math::Basis basis_a = math::make_basis(math::Quaternion(), pose.axis_a);
    math::Basis basis_b = math::make_basis(math::Quaternion(), pose.axis_b);
    Poly poly(basis_a, pose.radius_a / (math::PI_F / 2.0f), /*sides=*/6, 0.0f);
    Star star(basis_b, pose.radius_b, /*sides=*/5, 0.4f);

    SDF::Subtract<Poly, Star> carved(poly, star);
    expect_cull_covers_interior<W, H>(carved, "subtract leaf pair");
  }
}

/**
 * @brief Verifies the SmoothUnion interval cull covers every AA-fringe pixel of a
 *        real leaf pair.
 * @details The weld bulges the surface outside both children, so the cull rests
 *   on the k pad. The last pose
 *   centers both on +X so the padded spans straddle theta = 0.
 */
inline void test_smooth_union_cull_covers_fringe_over_leaf_pairs() {
  constexpr int W = 288, H = 144;
  using Poly = SDF::PlanarPolygon;

  struct Pose {
    math::Vector axis_a, axis_b;
    float radius_a, radius_b;
  };
  const Pose poses[] = {
      {math::Vector(0, 0, 1), math::Vector(0.6f, 0.2f, 1.0f), 0.9f, 0.8f},
      {math::Vector(0, 0, 1), math::Vector(0.6f, 0.2f, 1.0f), 0.5f, 0.4f},
      {math::Vector(0, 1, 0), math::Vector(0.5f, 1.0f, -0.2f), 0.45f, 0.35f},
      {math::Vector(-0.4f, 0.6f, 0.7f), math::Vector(0.1f, 0.9f, 0.4f), 0.6f,
       0.5f},
      {math::Vector(1, 0, 0), math::Vector(1.0f, 0.35f, 0.2f), 0.5f, 0.4f},
  };

  for (const Pose &pose : poses) {
    math::Basis basis_a = math::make_basis(math::Quaternion(), pose.axis_a);
    math::Basis basis_b = math::make_basis(math::Quaternion(), pose.axis_b);
    Poly poly_a(basis_a, pose.radius_a / (math::PI_F / 2.0f), /*sides=*/6,
                0.0f);
    Poly poly_b(basis_b, pose.radius_b / (math::PI_F / 2.0f), /*sides=*/5,
                0.4f);

    SDF::SmoothUnion<Poly, Poly> welded(poly_a, poly_b, /*k=*/0.25f);
    expect_cull_covers_fringe<W, H>(welded, "smooth union leaf pair");
  }
}

/**
 * @brief Verifies the SmoothUnion cull scans the rows past both children's
 *        spans.
 * @details Welding a polygon to itself dilates the surface by k/6, past the
 *   last row either child emits a span for; those rows must request a full
 *   scan. A row past the blend reach stays culled.
 */
inline void test_smooth_union_scans_rows_past_both_children() {
  constexpr int W = 96, H = 48;
  math::init_geometry_luts<W, H>();
  using Poly = SDF::PlanarPolygon;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0.2f, 1.0f, 0.1f));
  Poly poly(basis, /*radius=*/0.5f / (math::PI_F / 2.0f), /*sides=*/6, 0.0f);

  // k/6 = 0.2 rad of dilation, several pixel widths past the polygon's own cap.
  SDF::SmoothUnion<Poly, Poly> welded(poly, poly, /*k=*/1.2f);
  expect_cull_covers_interior<W, H>(welded, "smooth union self weld");

  SDF::SmoothUnion<Poly, Poly> tight(poly, poly, /*k=*/0.05f);
  const int far_row = poly.get_vertical_bounds<H>().y_max + 5;
  HS_EXPECT_LT(far_row, H);
  int spans = 0;
  bool handled = tight.get_horizontal_intervals<W, H>(
      far_row, [&](float, float) { ++spans; });
  HS_EXPECT_TRUE(handled);
  HS_EXPECT_EQ(spans, 0);
}

/**
 * @brief Verifies AngularRepeat around a non-Y axis culls in the full canvas, covering all copies.
 * @details A non-Y axis sweeps the folded copies through latitudes the
 *   un-repeated child never occupies, so get_vertical_bounds must fall back to
 *   the full canvas.
 */
inline void test_angular_repeat_non_y_axis_cull_covers_copies() {
  constexpr int W = 96, H = 48;
  SDF::Line ln(math::Vector(0.25f, 1, 0).normalized(),
               math::Vector(-0.25f, 1, 0).normalized(), /*thickness=*/0.12f);
  SDF::AngularRepeat<SDF::Line> rep(ln, /*reps=*/4, math::Vector(1, 0, 0));
  int interior = expect_cull_covers_interior<W, H>(rep, "angular repeat");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies a Y-axis AngularRepeat culls to its copies' columns without
 *        dropping any of them.
 * @details A Y-axis fold shifts azimuth by a whole sector and holds latitude,
 *   so the child's spans replayed once per copy bound every copy.
 */
inline void test_angular_repeat_y_axis_cull_narrows_rows() {
  constexpr int W = 288, H = 144;
  constexpr int REPS = 5;
  math::init_geometry_luts<W, H>();
  // Small star on the equator at azimuth 0, the sector the fold emits.
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(1, 0, 0));
  SDF::Star star(basis, /*radius=*/0.15f, /*sides=*/5, 0.0f);
  SDF::AngularRepeat<SDF::Star> rep(star, REPS, math::Vector(0, 1, 0));
  int fringe = expect_cull_covers_fringe<W, H>(rep, "angular repeat y axis");
  HS_EXPECT_GT(fringe, 0);

  int spans = 0;
  bool handled =
      rep.get_horizontal_intervals<W, H>(H / 2, [&](float, float) { ++spans; });
  HS_EXPECT_TRUE(handled);
  HS_EXPECT_EQ(spans, REPS);

  std::vector<uint8_t> visited;
  cull_visited<W, H>(rep, visited);
  int visited_px = 0;
  for (uint8_t v : visited)
    visited_px += v;
  auto bounds = rep.get_vertical_bounds<H>();
  const int rows =
      std::min(H - 1, bounds.y_max) - std::max(0, bounds.y_min) + 1;
  HS_EXPECT_GT(rows, 0);
  // Five narrow copies cannot reach half the band a full-width scan visits.
  HS_EXPECT_LT(visited_px, rows * W / 2);
}

/**
 * @brief Verifies a slightly tilted fold axis forfeits the Y-fold cull.
 * @details The copies of a tilted fold drift in latitude and azimuth by the
 *   tilt, past what the fold slop pads, so the child's band and spans no longer
 *   bound them.
 */
inline void test_angular_repeat_tilted_axis_forfeits_cull() {
  constexpr int W = 288, H = 144;
  constexpr int REPS = 5;
  math::init_geometry_luts<W, H>();
  // Exceeds the Y-fold axis tolerance (`ANGULAR_REPEAT_Y_AXIS_TOL_SQ`).
  const float tilt = 5e-3f;
  SDF::Line ln(math::Vector(0.25f, 1, 0).normalized(),
               math::Vector(-0.25f, 1, 0).normalized(), /*thickness=*/0.12f);
  SDF::AngularRepeat<SDF::Line> rep(ln, REPS,
                                    math::Vector(tilt, 1, 0).normalized());
  auto bounds = rep.get_vertical_bounds<H>();
  HS_EXPECT_EQ(bounds.y_min, 0);
  HS_EXPECT_EQ(bounds.y_max, H - 1);
  int spans = 0;
  bool handled =
      rep.get_horizontal_intervals<W, H>(H / 2, [&](float, float) { ++spans; });
  HS_EXPECT_FALSE(handled);
  HS_EXPECT_EQ(spans, 0);
}

/**
 * @brief Verifies the arc-extrema cull widens phi to a Line's great-circle bulge.
 * @details The Line's two endpoints share a latitude but its great-circle arc
 *   bulges to a pole between them.
 */
inline void test_line_arc_bulge_cull_covers_interior() {
  constexpr int W = 96, H = 48;
  // Endpoints at phi≈0.4 either side of +Y; the arc bulges through the north
  // pole (phi=0), above either endpoint's latitude.
  SDF::Line ln(math::Vector(0, cosf(0.4f), sinf(0.4f)),
               math::Vector(0, cosf(0.4f), -sinf(0.4f)), /*thickness=*/0.15f);
  int interior = expect_cull_covers_interior<W, H>(ln, "line arc bulge");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies the cull covers a Line whose antipodal endpoints select no arc.
 * @details With antipodal endpoints distance() measures the whole great
 *   circle, which the bounds must cover.
 */
inline void test_line_antipodal_cull_covers_interior() {
  const math::Vector a(-0.21973225f, -0.52185529f, -0.82424802f);
  const math::Vector b(0.21973109f, 0.52185678f, 0.82424736f);
  SDF::Line jittered(a, b, 0.01f);
  const auto n = math::cross(a, math::Vector(0, 1, 0)).normalized();
  const auto p = math::cross(n, a).normalized();
  SDF::DistanceResult result;
  jittered.distance(p, result);
  HS_EXPECT_LT(result.dist, 0.0f);
  constexpr int W = 96, H = 48;
  const math::Vector ENDPOINT = math::Vector(0.4f, 0.6f, 0.69f).normalized();
  SDF::Line ln(ENDPOINT, -ENDPOINT, /*thickness=*/0.15f);
  int interior = expect_cull_covers_interior<W, H>(ln, "line antipodal");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies the cull covers a Line whose bounding cap radius exceeds pi.
 * @details A quarter arc with a stroke this wide gives half-length + thickness
 *   ≈ 3.39 rad, past the pi where cos turns back up.
 */
inline void test_line_thick_cap_past_pi_cull_covers_interior() {
  constexpr int W = 96, H = 48;
  SDF::Line ln(math::Vector(1, 0, 0), math::Vector(0, 0, 1),
               /*thickness=*/2.6f);
  int interior = expect_cull_covers_interior<W, H>(ln, "line thick cap");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies the Ring interval cull covers interior pixels for thin rings
 *        whose band wraps a pole while the centerline still takes the fast path.
 */
inline void test_ring_pole_wrap_cull_covers_interior() {
  constexpr int W = 256, H = 128;
  struct Cfg {
    float tilt, radius, thickness;
  };
  const Cfg cfgs[] = {
      {0.15f, 0.22f, 0.13f},
      {0.16f, 0.25f, 0.14f},
      {0.16f, 0.33f, 0.15f},
      {0.16f, 0.39f, 0.15f},
  };
  for (const Cfg &c : cfgs) {
    math::Basis basis_n =
        math::make_basis(math::Quaternion(), math::Vector(c.tilt, 1.0f, 0.0f));
    SDF::Ring ring_n(basis_n, c.radius, c.thickness);
    expect_cull_covers_interior<W, H>(ring_n, "ring north pole");

    math::Basis basis_s =
        math::make_basis(math::Quaternion(), math::Vector(c.tilt, -1.0f, 0.0f));
    SDF::Ring ring_s(basis_s, c.radius, c.thickness);
    expect_cull_covers_interior<W, H>(ring_s, "ring south pole");
  }
}

/**
 * @brief Verifies the DistortedRing cull drops no interior arc column under
 *        high-frequency centerline shifts with exact max_distortion bounds.
 * @details max_distortion widens the row/column/per-pixel reject bands; each
 *   shift_fn passes its exact analytic peak.
 */
inline void test_distorted_ring_cull_covers_interior_high_freq() {
  constexpr int W = 256, H = 128;
  struct Cfg {
    float amp;
    int harmonic;
    float phase_frac;
  };
  const Cfg cfgs[] = {
      {0.18f, 127, 0.5f / 256.0f},
      {0.20f, 255, 0.5f / 256.0f},
      {0.15f, 384, 0.25f / 256.0f},
      {0.22f, 200, 0.5f / 256.0f},
  };
  const math::Vector axes[] = {math::Vector(0, 1, 0),
                               math::Vector(0.3f, 1.0f, 0.2f),
                               math::Vector(1, 0, 0)};
  for (const math::Vector &axis : axes) {
    math::Basis basis = math::make_basis(math::Quaternion(), axis);
    for (const Cfg &c : cfgs) {
      float amp = c.amp;
      int harmonic = c.harmonic;
      float ph = c.phase_frac;
      auto shift = [amp, harmonic, ph](float t) {
        return amp * std::sin(2.0f * math::PI_F * (harmonic * t + ph));
      };
      SDF::DistortedRing ring(basis, /*radius=*/0.6f, /*thickness=*/0.12f,
                              shift,
                              /*max_distortion=*/amp, /*phase=*/0.0f);
      expect_cull_covers_interior<W, H>(ring, "distorted ring");

      // The knot overload's exact distance widens the interior toward the
      // band edges.
      constexpr int LUT_N = 256;
      float knots[LUT_N + 1];
      for (int k = 0; k <= LUT_N; ++k)
        knots[k] = shift(static_cast<float>(k % LUT_N) / LUT_N);
      SDF::KnotPrefilter poly_pf;
      SDF::DistortedRing poly(basis, /*radius=*/0.6f, /*thickness=*/0.12f,
                              knots, LUT_N, /*phase=*/0.0f, poly_pf);
      expect_cull_covers_interior<W, H>(poly, "distorted polygon");
    }
  }

  constexpr int ASYM_LUT_N = 64;
  float asymmetric[ASYM_LUT_N + 1];
  for (int k = 0; k <= ASYM_LUT_N; ++k)
    asymmetric[k] =
        0.075f + 0.1f * sinf(6.0f * math::PI_F * (k % ASYM_LUT_N) / ASYM_LUT_N);
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0.3f, 1.0f, 0.2f));
  SDF::KnotPrefilter asymmetric_pf;
  SDF::DistortedRing asymmetric_ring(basis, 0.6f, 0.08f, asymmetric, ASYM_LUT_N,
                                     0.0f, asymmetric_pf);
  expect_cull_covers_interior<W, H>(asymmetric_ring, "asymmetric ring");
}

/**
 * @brief Verifies the Face azimuth-interval cull covers every paintable pixel
 *        (including the outer AA fringe column at a silhouette edge).
 * @return The paintable-pixel count, so the caller confirms the case is non-trivial.
 * @details A pixel is paintable when its exact distance < pixel_width.
 */
template <int W, int H>
inline int expect_face_cull_covers_fringe(int sides, float rho,
                                          const math::Vector &axis) {
  HS_CONTEXT("face fringe", sides, rho);
  HS_EXPECT_TRUE(sides >= 3 && sides <= 8);
  if (sides < 3 || sides > 8)
    return 0;
  constexpr int HV = H + hs::H_OFFSET;
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  math::Basis basis = math::make_basis(math::Quaternion(), axis);
  math::Vector verts3d[8];
  uint16_t idx[8];
  for (int i = 0; i < sides; ++i) {
    float a = (2.0f * math::PI_F * i) / sides + 0.37f;
    verts3d[i] = (basis.v * cosf(rho) +
                  (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(rho))
                     .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }
  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const math::Vector>(verts3d, sides),
                 std::span<const uint16_t>(idx, sides), scratch, HV, H);

  const auto bounds = face.get_vertical_bounds<H>();
  for (int y : {bounds.y_min - 1, bounds.y_max + 1}) {
    int emitted = 0;
    HS_EXPECT_TRUE((face.get_horizontal_intervals<W, H>(
        y, [&](float, float) { ++emitted; })));
    HS_EXPECT_EQ(emitted, 0);
  }

  std::vector<uint8_t> visited;
  cull_visited<W, H>(face, visited);

  const float pixel_width = 2.0f * math::PI_F / W;
  int paintable = 0;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      if (SDF::distance_of(face, p).dist < pixel_width) {
        ++paintable;
        HS_EXPECT_TRUE(visited[static_cast<size_t>(y) * W + x]);
      }
    }
  }
  return paintable;
}

/** @brief Pins azimuth culling to the columns emitted at a rounded boundary. */
inline void test_face_azimuth_cull_matches_boundary_column() {
  constexpr int W = 288, H = 144, Y = 72;
  const math::Vector VERTS[] = {math::Vector(1.0f, 0.1f, 0.0f).normalized(),
                                math::Vector(1.0f, -0.1f, 0.1f).normalized(),
                                math::Vector(1.0f, -0.1f, -0.1f).normalized()};
  const uint16_t INDICES[] = {0, 1, 2};
  SDF::FaceScratchBuffer scratch;
  SDF::Face face(VERTS, INDICES, scratch, math::LatitudeGeometry(H), H);
  std::array<float, H> pads;
  pads.fill(0.05f);
  SDF::Interval interval{0x1.e3930ap-1f, 1.1f};
  face.intervals = std::span<SDF::Interval>(&interval, 1);
  face.full_width = false;
  face.y_min = Y;
  face.y_max = Y;
  face.build_azimuth_pads = pads.data();
  int start = -1, end = -1;
  HS_EXPECT_TRUE((face.get_horizontal_intervals<W, H>(Y, [&](float a, float b) {
    start = static_cast<int>(a);
    end = static_cast<int>(b);
  })));
  HS_EXPECT_EQ(start, 40);
  HS_EXPECT_GT(end, start);
  ClipRegion clip{.y_start = Y,
                  .y_end = Y + 1,
                  .x_start = start,
                  .x_end = start + 1,
                  .margin = 0,
                  .w = W,
                  .h = H};
  HS_EXPECT_TRUE(!face.clip_rejects_azimuth(clip, Y, Y));
  clip.x_start = start - 1;
  clip.x_end = start;
  HS_EXPECT_TRUE(face.clip_rejects_azimuth(clip, Y, Y));
}

/** @brief Verifies Face culling includes the AA fringe near column boundaries. */
inline void test_face_cull_covers_aa_fringe() {
  constexpr int W = 256, H = 128;
  struct Cfg {
    int sides;
    float rho;
    math::Vector axis;
  };
  const Cfg cfgs[] = {
      {3, 0.15f, math::Vector(0, 0, 1)},
      {3, 0.30f, math::Vector(0, 0, 1)},
      {3, 0.48f, math::Vector(0, 0, 1)},
      {3, 0.18f, math::Vector(1, 0, 0)},
  };
  int total_paintable = 0;
  for (const Cfg &c : cfgs) {
    const int paintable =
        expect_face_cull_covers_fringe<W, H>(c.sides, c.rho, c.axis);
    HS_EXPECT_GT(paintable, 0);
    total_paintable += paintable;
  }
  HS_EXPECT_GT(total_paintable, 1000);
}

/** @brief Vertical bounds cover wide AA bands on tall, low-width grids. */
inline void test_face_vertical_margin_tracks_pixel_width() {
  constexpr int W = 64, H = 144, HV = H + hs::H_OFFSET;
  const float tilt = math::PI_F / 60.0f;
  const math::Vector axis(sinf(tilt) * cosf(0.37f), cosf(tilt),
                          sinf(tilt) * sinf(0.37f));
  const math::Basis basis = math::make_basis(math::Quaternion(), axis);
  math::Vector vertices[3];
  const uint16_t indices[] = {0, 1, 2};
  for (int i = 0; i < 3; ++i) {
    const float angle = math::TWO_PI_F * i / 3.0f + 0.37f;
    vertices[i] =
        (basis.v * cosf(0.025f) +
         (basis.u * cosf(angle) + basis.w * sinf(angle)) * sinf(0.025f))
            .normalized();
  }
  SDF::FaceScratchBuffer original_scratch, widened_scratch;
  SDF::Face original(vertices, indices, original_scratch, HV, H);
  SDF::Face widened(vertices, indices, widened_scratch, HV, H, nullptr, nullptr,
                    std::max(SDF::BOUNDS_MARGIN, math::TWO_PI_F / W));
  const auto original_bounds = original.get_vertical_bounds<H>();
  const auto widened_bounds = widened.get_vertical_bounds<H>();
  int original_misses = 0;
  int covered = 0;
  for (int y = 0; y < H; ++y) {
    const float phi = math::PI_F * y / (HV - 1);
    for (int x = 0; x < W; ++x) {
      const float theta = math::TWO_PI_F * x / W;
      const math::Vector point(sinf(phi) * cosf(theta), cosf(phi),
                               sinf(phi) * sinf(theta));
      if (SDF::distance_of(widened, point).dist >= math::TWO_PI_F / W)
        continue;
      ++covered;
      original_misses += y < original_bounds.y_min || y > original_bounds.y_max;
      HS_EXPECT_TRUE(y >= widened_bounds.y_min && y <= widened_bounds.y_max);
    }
  }
  HS_EXPECT_GT(covered, 0);
  HS_EXPECT_GT(original_misses, 0);

  SDF::FaceScratchBuffer shipping_scratch;
  SDF::Face shipping(vertices, indices, shipping_scratch, HV, H, nullptr,
                     nullptr,
                     std::max(SDF::BOUNDS_MARGIN, math::TWO_PI_F / 288));
  const auto shipping_bounds = shipping.get_vertical_bounds<H>();
  HS_EXPECT_EQ(shipping_bounds.y_min, original_bounds.y_min);
  HS_EXPECT_EQ(shipping_bounds.y_max, original_bounds.y_max);
}

/** @brief Records the columns emitted by Face's fixed equatorial pad. */
template <int W, int H>
inline void face_fixed_pad_visited(const SDF::Face &face,
                                   std::vector<uint8_t> &visited) {
  visited.assign(static_cast<size_t>(W) * H, 0);
  auto bounds = face.template get_vertical_bounds<H>();
  const int y_lo = std::max(0, bounds.y_min);
  const int y_hi = std::min(H - 1, bounds.y_max);
  if (y_lo > y_hi)
    return;
  const float pad = SDF::face_azimuth_pad(W);
  Scan::scan_region<W, H>(
      y_lo, y_hi,
      [&](int, auto &&out) {
        if (face.full_width)
          return false;
        for (const auto &iv : face.intervals) {
          out(floorf((iv.start - pad) * W / math::TWO_PI_F),
              ceilf((iv.end + pad) * W / math::TWO_PI_F));
        }
        return true;
      },
      [&](int wx, int y, const math::Vector &, int run) {
        for (int i = 0; i < run; ++i)
          if (wx + i >= 0 && wx + i < W)
            visited[static_cast<size_t>(y) * W + wx + i] = 1;
        return run;
      });
}

/** @brief Counts paintable pixels missed before and after latitude widening. */
template <int W, int H>
inline std::pair<int, int> face_fringe_misses(const SDF::Face &face) {
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  std::vector<uint8_t> fixed_visited;
  std::vector<uint8_t> widened_visited;
  face_fixed_pad_visited<W, H>(face, fixed_visited);
  cull_visited<W, H>(face, widened_visited);

  const float pixel_width = math::TWO_PI_F / W;
  int fixed_misses = 0;
  int widened_misses = 0;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      if (SDF::distance_of(face, p).dist >= pixel_width)
        continue;
      const size_t px = static_cast<size_t>(y) * W + x;
      fixed_misses += !fixed_visited[px];
      widened_misses += !widened_visited[px];
    }
  }
  return {fixed_misses, widened_misses};
}

/**
 * @brief Pins the reduction from converting Face's AA reach to azimuth by row.
 * @details The widened mask closes at least 92% of the fixed-pad misses in this
 * deterministic sample and no configuration in it regresses.
 */
inline void test_face_latitude_pad_reduces_fringe_drops() {
  constexpr int W = 288, H = 144, HV = H + hs::H_OFFSET;
  struct Cfg {
    int sides;
    float rho;
    float phi;
  };
  const Cfg cfgs[] = {
      {3, 0.14f, 2.58f}, {3, 0.14f, 2.76f}, {3, 0.20f, 2.38f},
      {3, 0.20f, 2.58f}, {3, 0.26f, 2.76f}, {3, 0.26f, 0.38f},
      {3, 0.26f, 0.56f}, {5, 0.14f, 0.38f}, {4, 0.20f, 0.38f},
      {3, 0.08f, 0.38f},
  };
  int fixed_misses = 0;
  int widened_misses = 0;
  for (const Cfg &cfg : cfgs) {
    const float theta = 0.31f * cfg.sides + 1.7f * cfg.rho + 0.23f * cfg.phi;
    const math::Vector axis(sinf(cfg.phi) * cosf(theta), cosf(cfg.phi),
                            sinf(cfg.phi) * sinf(theta));
    math::Basis basis = math::make_basis(math::Quaternion(), axis);
    math::Vector verts[6];
    uint16_t idx[6];
    for (int i = 0; i < cfg.sides; ++i) {
      const float a = math::TWO_PI_F * i / cfg.sides + 0.37f;
      verts[i] = (basis.v * cosf(cfg.rho) +
                  (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(cfg.rho))
                     .normalized();
      idx[i] = static_cast<uint16_t>(i);
    }
    SDF::FaceScratchBuffer scratch;
    SDF::Face face(std::span<const math::Vector>(verts, cfg.sides),
                   std::span<const uint16_t>(idx, cfg.sides), scratch, HV, H);
    const auto misses = face_fringe_misses<W, H>(face);
    HS_EXPECT_LE(misses.second, misses.first);
    fixed_misses += misses.first;
    widened_misses += misses.second;
  }
  HS_EXPECT_GE(fixed_misses, 100);
  HS_EXPECT_LT(widened_misses, fixed_misses);
  HS_EXPECT_LE(widened_misses * 100, fixed_misses * 8);
}

/**
 * @brief Rasterizes a wedge face apexed on a pole and compares it against a
 *        full-canvas distance scan.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param pole_y +1 for the north pole, -1 for the south pole.
 * @return Count of pixels the scan expects painted.
 * @details The apex projects onto a polygon vertex, so the face classifies as
 *   PoleHit::BOUNDARY. The expectation covers the AA fringe, whose azimuth
 *   reach grows as the rows approach the apex.
 */
template <int W, int H>
inline int expect_pole_vertex_face_matches_full_scan(float pole_y) {
  constexpr int HV = H + hs::H_OFFSET;
  constexpr int N_VERTS = 4;
  const float rho = 0.7f;

  math::Vector verts[N_VERTS];
  uint16_t idx[N_VERTS];
  verts[0] = math::Vector(0, pole_y, 0);
  for (int i = 1; i < N_VERTS; ++i) {
    const float a = pole_y * 0.55f * static_cast<float>(i - 2);
    verts[i] = math::Vector(sinf(rho) * cosf(a), pole_y * cosf(rho),
                            sinf(rho) * sinf(a));
  }
  for (int i = 0; i < N_VERTS; ++i)
    idx[i] = static_cast<uint16_t>(i);

  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const math::Vector>(verts, N_VERTS),
                 std::span<const uint16_t>(idx, N_VERTS), scratch, HV, H);
  // Pins the setup on the BOUNDARY branch: !full_width alone also holds for a
  // face that misses the pole entirely.
  const float inv_c = 1.0f / std::abs(face.center.y);
  const float sgn = face.center.y > 0 ? 1.0f : -1.0f;
  HS_EXPECT_TRUE(face.pole_hit(sgn * face.basis_u.y * inv_c,
                               sgn * face.basis_w.y * inv_c) ==
                 SDF::Face::PoleHit::BOUNDARY);
  HS_EXPECT_TRUE(!face.full_width);

  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(fx);
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
    };
    Scan::rasterize_face<W, H>(pipeline, canvas, face, shader);
  }
  fx.advance_display();

  const float pixel_width = 2.0f * math::PI_F / W;
  int painted = 0;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      const float d = SDF::distance_of(face, p).dist;
      const Pixel px = fx.get_pixel(x, y);
      const bool lit = px.r != 0 || px.g != 0 || px.b != 0;
      if (d >= pixel_width) {
        HS_EXPECT_TRUE(!lit);
        continue;
      }
      // Dead band around Scan::MIN_ALPHA: an alpha that rounds the shade to
      // black is neither required nor forbidden.
      const float alpha = Scan::solid_coverage(d, pixel_width);
      if (alpha > 0.05f) {
        ++painted;
        HS_EXPECT_TRUE(lit);
      }
    }
  }
  return painted;
}

/** @brief Covers the pole-boundary classification and its per-row raster at both poles. */
inline void test_face_pole_vertex_matches_full_scan() {
  constexpr int W = 96, H = 64;
  const int north = expect_pole_vertex_face_matches_full_scan<W, H>(1.0f);
  const int south = expect_pole_vertex_face_matches_full_scan<W, H>(-1.0f);
  HS_EXPECT_GT(north, 100);
  HS_EXPECT_GT(south, 100);
}
