/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// Plot::screen_step — adaptive-density sub-step
// ============================================================================

/**
 * @brief Pins screen_step against an independent screen-velocity oracle in the
 *        unclamped regime.
 * @details Reconstructs the pixel speed by finite-differencing the canvas map
 *          (x = longitude·W/2π, y = colatitude·(H_VIRT-1)/π) along a small
 *          geodesic step in the tangent direction, then asserts screen_step
 *          returns SCREEN_STEP_PX/|v_screen|. Inputs are chosen so the result
 *          lands strictly between the pole and equator clamps. Because the speed
 *          squares the velocity components, flipping the tangent's sign leaves
 *          the step unchanged — checked here.
 */
inline void test_screen_step_matches_analytic_unclamped() {
  constexpr int W = 288, H = 144;
  constexpr int H_VIRT = H + hs::H_OFFSET;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  const float KX = W / (2.0f * math::PI_F);
  const float KY = (H_VIRT - 1) / math::PI_F;

  auto screen_speed = [&](const math::Vector &pos, const math::Vector &tan) {
    const float ds = 1e-4f;
    auto xy = [&](const math::Vector &p) {
      float phi = std::acos(hs::clamp(p.y, -1.0f, 1.0f));
      float lam = std::atan2(p.z, p.x);
      return std::pair<float, float>(lam * KX, phi * KY);
    };
    math::Vector pp = (pos * std::cos(ds) + tan * std::sin(ds)).normalized();
    math::Vector pm = (pos * std::cos(ds) - tan * std::sin(ds)).normalized();
    auto a = xy(pp);
    auto b = xy(pm);
    float dx = (a.first - b.first) / (2.0f * ds);
    float dy = (a.second - b.second) / (2.0f * ds);
    return std::sqrt(dx * dx + dy * dy);
  };

  struct Case {
    math::Vector pos, tan;
  };
  // The tangent must be perpendicular to pos (a genuine arc-length tangent), or
  // the geodesic oracle and screen_step's formula parametrize differently.
  const math::Vector pos_off(std::sqrt(0.75f), 0.5f, 0.0f);
  const math::Vector raw(0.0f, 1.0f,
                         1.0f); // arbitrary; projected onto the tangent plane
  const math::Vector tan_mixed =
      (raw - pos_off * math::dot(raw, pos_off)).normalized();
  // Equatorial longitudinal, off-equator longitudinal, and a mixed (colatitude +
  // longitude) tangent — all verified below to land inside the clamp window.
  const Case cases[] = {
      {math::Vector(1.0f, 0.0f, 0.0f), math::Vector(0.0f, 0.0f, 1.0f)},
      {pos_off, math::Vector(0.0f, 0.0f, 1.0f)},
      {pos_off, tan_mixed},
  };

  const float lo = base_step * Plot::MIN_POLE_SCALE;
  for (const Case &c : cases) {
    float speed = screen_speed(c.pos, c.tan);
    float expected = Plot::SCREEN_STEP_PX / speed;
    // Unclamped regime guard: the analytic step sits strictly inside the window.
    HS_EXPECT_GT(expected, lo);
    HS_EXPECT_LT(expected, base_step);

    float got = Plot::screen_step<W, H>(c.pos, c.tan, base_step);
    HS_EXPECT_NEAR(got, expected, expected * 2e-3f);

    // Speed squares the tangent, so the sign drops out.
    float flipped = Plot::screen_step<W, H>(c.pos, c.tan * -1.0f, base_step);
    HS_EXPECT_NEAR(flipped, got, 1e-6f);
  }
}

/**
 * @brief Verifies edge_fits_one_dot is a strict tightening of the rasterizer
 *        fast-path test, so a routed edge always renders as the same single
 *        dot the full path would emit.
 * @details For every accepted edge: theta >= EPS_GEOMETRIC (the routed edge is
 *          not one process_segment treats as degenerate) and theta <=
 *          screen_step at the edge start with the geodesic tangent built
 *          exactly as rasterize_geodesic_strategy builds it. Sweeps random
 *          headings across latitudes (poles included) and log-spaced arc
 *          lengths around the one-pixel scale, and checks the predicate
 *          accepts a healthy share of genuinely sub-pixel edges.
 */
inline void test_edge_fits_one_dot_is_conservative() {
  constexpr int W = 288, H = 144;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  hs::random().seed(20260720);

  int accepted = 0;
  for (int n = 0; n < 20000; ++n) {
    math::Vector a;
    if (n % 7 == 0) {
      // Pole-proximal start: predicate must stay safe as sin(phi) -> 0.
      float e = hs::rand_f(0.0f, 0.05f);
      float az = hs::rand_f(0.0f, 2.0f * math::PI_F);
      float s = std::sin(e);
      a = math::Vector(s * std::cos(az),
                       (n % 14 == 0) ? std::cos(e) : -std::cos(e),
                       s * std::sin(az));
    } else {
      for (;;) {
        const float rx = hs::rand_f(-1, 1);
        const float ry = hs::rand_f(-1, 1);
        const float rz = hs::rand_f(-1, 1);
        math::Vector r(rx, ry, rz);
        if (r.length() > 0.1f) {
          a = r.normalized();
          break;
        }
      }
    }
    const float raw_x = hs::rand_f(-1, 1);
    const float raw_y = hs::rand_f(-1, 1);
    const float raw_z = hs::rand_f(-1, 1);
    math::Vector raw(raw_x, raw_y, raw_z);
    math::Vector t = raw - a * math::dot(raw, a);
    if (t.length() < 1e-3f)
      continue;
    t = t.normalized();
    float theta = base_step * std::pow(10.0f, hs::rand_f(-4.0f, 0.5f));
    math::Vector b = (a * std::cos(theta) + t * std::sin(theta)).normalized();

    if (!Plot::edge_fits_one_dot<W, H>(a, b))
      continue;
    ++accepted;

    float total = math::angle_between(a, b);
    HS_EXPECT_GE(total, math::EPS_GEOMETRIC);
    math::Vector axis = math::cross(a, b).normalized();
    math::Vector v_perp = math::cross(axis, a);
    float first_step = Plot::screen_step<W, H>(a, v_perp, base_step);
    HS_EXPECT_LE(total, first_step);
  }
  // The predicate must actually fire on the sub-pixel population it targets.
  HS_EXPECT_GT(accepted, 500);
}

/**
 * @brief Pins the one-dot gate to AntiAlias's seam, row, and tap-cutoff rules.
 */
inline void test_antialiased_dot_clip_footprint() {
  constexpr int W = 96, H = 48;
  ClipRegion cr;
  cr.w = W;
  cr.h = H;
  cr.margin = 0;
  cr.y_start = 10;
  cr.y_end = 11;
  cr.x_start = 0;
  cr.x_end = 1;
  auto xc = cr.x_clip();

  // A fractional dot just left of the cylindrical seam splats into column 0.
  HS_EXPECT_TRUE(
      (Plot::antialiased_dot_visible_in_clip<W, H>(cr, xc, 10.25f, W - 0.25f)));
  // An integral seam-neighbor emits only at W-1; its zero-weight x1 tap at 0
  // must not make the dot visible.
  HS_EXPECT_FALSE(
      (Plot::antialiased_dot_visible_in_clip<W, H>(cr, xc, 10.0f, W - 1.0f)));

  cr.y_start = H - 1;
  cr.y_end = H;
  cr.x_start = 20;
  cr.x_end = 21;
  xc = cr.x_clip();
  // AntiAlias renormalizes a splat straddling the physical bottom row onto the
  // in-range tap.
  HS_EXPECT_TRUE(
      (Plot::antialiased_dot_visible_in_clip<W, H>(cr, xc, H - 0.25f, 20.0f)));
  HS_EXPECT_FALSE(
      (Plot::antialiased_dot_visible_in_clip<W, H>(cr, xc, H - 2.0f, 20.0f)));
}

/**
 * @brief Keeps a geodesic whose upper AntiAlias tap reaches the render margin.
 */
inline void test_geodesic_edge_gate_keeps_upper_antialias_tap() {
  constexpr int W = 96, H = 48;
  constexpr float EDGE_ROW = 10.25f;
  const float phi = math::y_to_phi_virtual(EDGE_ROW, H + hs::H_OFFSET);
  const float sp = sinf(phi);
  const float cp = cosf(phi);
  const float theta = 0.2f;
  const math::Vector a(sp * cosf(theta), cp, sp * sinf(theta));
  const math::Vector b(sp * cosf(theta), cp, -sp * sinf(theta));
  const Plot::GeodesicEdgeSpan es = Plot::make_geodesic_edge_span(a, b);

  ClipRegion cr;
  cr.w = W;
  cr.h = H;
  cr.margin = 1;
  cr.y_start = 12;
  cr.y_end = 20;
  cr.x_start = 0;
  cr.x_end = W;
  const auto xc = cr.x_clip();

  const float ra = Plot::y_to_screen_row<H>(a.y);
  const float rb = Plot::y_to_screen_row<H>(b.y);
  float row_lo, row_hi;
  Plot::geodesic_row_span_rows<H>(ra, rb, a, b, es, row_lo, row_hi);
  HS_EXPECT_FALSE(cr.could_intersect_y(row_lo, row_hi));
  HS_EXPECT_GT(row_hi, cr.render_y_start() - Plot::GEODESIC_ROW_AA_PAD);
  HS_EXPECT_TRUE((Plot::exact_geodesic_edge_visible<W, H>(
      cr, xc, ra, rb, a, b, es, [](int &, int &) { return false; })));
  HS_EXPECT_EQ((Plot::raw_geodesic_edge_gate<W, H>(
                   cr, xc, ra, rb, math::vector_to_theta<W>(a),
                   math::vector_to_theta<W>(b), a, b)),
               Plot::RawGeodesicGateResult::VISIBLE);

  Filter::Screen::AntiAlias<W, H> aa;
  bool emitted_in_margin = false;
  bool emitted_in_display = false;
  aa.plot(20.25f, row_hi, Pixel(1, 2, 3), 0.0f, 1.0f,
          [&](float, float y, const Pixel &, float, float) {
            const int row = static_cast<int>(y);
            emitted_in_margin |= cr.contains_y(row);
            emitted_in_display |= row >= cr.y_start && row < cr.y_end;
          });
  HS_EXPECT_TRUE(emitted_in_margin);
  HS_EXPECT_FALSE(emitted_in_display);
}

/**
 * @brief Pins the one-dot gate to the taps Screen::AntiAlias actually emits.
 * @details The gate hand-mirrors the filter's splat geometry, so a sweep of
 * sub-pixel positions against several clip bands compares the gate's verdict
 * with the filter's own tap set: visible iff some emitted tap lands inside the
 * band.
 */
inline void test_antialiased_dot_gate_matches_antialias_taps() {
  constexpr int W = 32, H = 24;
  Filter::Screen::AntiAlias<W, H> aa;

  struct Band {
    int y_start, y_end, x_start, x_end, margin;
  };
  const Band bands[] = {
      {0, H, 0, W, 1},   {10, 14, 0, W, 1},    {10, 14, 4, 9, 0},
      {10, 14, 4, 9, 1}, {H - 2, H, 6, 7, 0},  {0, 2, 0, W, 0},
      {8, 12, 0, 3, 2},  {8, 12, W - 3, W, 2},
  };

  int mismatches = 0, visible = 0, hidden = 0;
  for (const Band &b : bands) {
    ClipRegion cr;
    cr.w = W;
    cr.h = H;
    cr.margin = b.margin;
    cr.y_start = b.y_start;
    cr.y_end = b.y_end;
    cr.x_start = b.x_start;
    cr.x_end = b.x_end;
    const ClipRegion::XClip xc = cr.x_clip();

    for (int ry = -4; ry <= 4 * H + 4; ++ry) {
      const float row = 0.25f * static_cast<float>(ry);
      for (int rx = 0; rx < 8 * W; ++rx) {
        const float col = 0.125f * static_cast<float>(rx);
        bool any_tap_in_band = false;
        aa.plot(col, row, Pixel(1, 2, 3), 0.0f, 1.0f,
                [&](float tx, float ty, const Pixel &, float, float) {
                  if (cr.contains_y(static_cast<int>(ty)) &&
                      !xc.clipped(static_cast<int>(tx)))
                    any_tap_in_band = true;
                });
        const bool gate =
            Plot::antialiased_dot_visible_in_clip<W, H>(cr, xc, row, col);
        if (gate != any_tap_in_band)
          ++mismatches;
        if (any_tap_in_band)
          ++visible;
        else
          ++hidden;
      }
    }
  }
  HS_EXPECT_EQ(mismatches, 0);
  // Both verdicts must occur, or the sweep proves nothing.
  HS_EXPECT_GT(visible, 0);
  HS_EXPECT_GT(hidden, 0);
}
