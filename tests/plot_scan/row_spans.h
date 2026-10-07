/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Plot::geodesic_row_span / Plot::planar_row_span — arc-aware clip cull
// ============================================================================

/**
 * @brief Verifies the row-span helpers conservatively cover the rendered arc's
 *        screen-row extent, including the interior latitude bulge where the arc
 *        reaches rows beyond both endpoints.
 */
inline void test_row_span_covers_arc_bulge() {
  constexpr int TW = 288, TH = 144;
  auto row_of = [](const math::Vector &v) {
    return math::vector_to_pixel<TW, TH>(v.normalized()).y;
  };

  hs::random().seed(20260609);
  int bulge_cases = 0; // edges whose arc bulges past both endpoints

  for (int trial = 0; trial < 6000; ++trial) {
    const bool planar = (trial & 1);
    math::Vector a, b;
    math::Basis basis;
    const math::Basis *pb = nullptr;

    if (planar) {
      const float cx = hs::rand_f(-1, 1);
      const float cy = hs::rand_f(-1, 1);
      const float cz = hs::rand_f(-1, 1);
      math::Vector center(cx, cy, cz);
      if (center.length() < 0.1f)
        continue;
      basis = basis_from_normal(center.normalized());
      pb = &basis;
      random_disk_edge(basis, a, b);
      // Antipodal-seam edges fall back to the geodesic strategy; skip them.
      if (math::dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
          math::dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
        continue;
    } else {
      const float rax = hs::rand_f(-1, 1);
      const float ray = hs::rand_f(-1, 1);
      const float raz = hs::rand_f(-1, 1);
      const float rbx = hs::rand_f(-1, 1);
      const float rby = hs::rand_f(-1, 1);
      const float rbz = hs::rand_f(-1, 1);
      math::Vector ra(rax, ray, raz);
      math::Vector rb(rbx, rby, rbz);
      if (ra.length() < 0.1f || rb.length() < 0.1f)
        continue;
      a = ra.normalized();
      b = rb.normalized();
      if (math::angle_between(a, b) < 0.05f)
        continue;
    }

    // Dense ground truth for the rendered arc's row extent.
    float t_lo = 1e9f, t_hi = -1e9f;
    constexpr int N = 2000;
    std::pair<float, float> p1{}, p2{};
    float ang = 0.0f;
    math::Vector vperp = a;
    if (planar) {
      p1 = Plot::azimuthal_project(a, basis);
      p2 = Plot::azimuthal_project(b, basis);
    } else {
      ang = math::angle_between(a, b);
      vperp = math::cross(math::cross(a, b).normalized(), a);
    }
    for (int i = 0; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      math::Vector p = planar
                           ? Plot::azimuthal_unproject(
                                 p1.first + (p2.first - p1.first) * t,
                                 p1.second + (p2.second - p1.second) * t, basis)
                           : (a * cosf(ang * t) + vperp * sinf(ang * t));
      float r = row_of(p);
      t_lo = std::min(t_lo, r);
      t_hi = std::max(t_hi, r);
    }

    float n_lo, n_hi;
    if (pb != nullptr)
      Plot::planar_row_span<TH>(a, b, Plot::make_planar_edge_span(a, b, *pb),
                                n_lo, n_hi);
    else
      Plot::geodesic_row_span<TH>(a, b, Plot::make_geodesic_edge_span(a, b),
                                  n_lo, n_hi);

    // Span must conservatively contain the arc (sub-pixel tolerance for fast-math).
    HS_EXPECT_LE(n_lo, t_lo + 0.25f);
    HS_EXPECT_GE(n_hi, t_hi - 0.25f);

    // Count edges whose interior bulges past the endpoints.
    float e_lo = std::min(row_of(a), row_of(b));
    float e_hi = std::max(row_of(a), row_of(b));
    if (t_lo < e_lo - 1.0f || t_hi > e_hi + 1.0f)
      bulge_cases++;
  }

  HS_EXPECT_GT(bulge_cases, 500);

  // Exact-antipodal geodesic edges: ground truth is the semicircle about
  // stable_perpendicular_axis.
  for (int trial = 0; trial < 500; ++trial) {
    const float rax = hs::rand_f(-1, 1);
    const float ray = hs::rand_f(-1, 1);
    const float raz = hs::rand_f(-1, 1);
    math::Vector ra(rax, ray, raz);
    if (ra.length() < 0.1f)
      continue;
    math::Vector a = ra.normalized();
    math::Vector b = a * -1.0f;
    math::Vector axis = Plot::stable_perpendicular_axis(a);
    math::Vector vperp = math::cross(axis, a);

    float t_lo = 1e9f, t_hi = -1e9f;
    constexpr int N = 2000;
    for (int i = 0; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      math::Vector p = a * cosf(math::PI_F * t) + vperp * sinf(math::PI_F * t);
      float r = row_of(p);
      t_lo = std::min(t_lo, r);
      t_hi = std::max(t_hi, r);
    }

    float n_lo, n_hi;
    Plot::geodesic_row_span<TH>(a, b, Plot::make_geodesic_edge_span(a, b), n_lo,
                                n_hi);
    HS_EXPECT_LE(n_lo, t_lo + 0.25f);
    HS_EXPECT_GE(n_hi, t_hi - 0.25f);
  }
}

/**
 * @brief Verifies cap_may_touch_clip never rejects a cap that reaches the
 *        clip's render region.
 * @details Only false positives are admissible.
 */
inline void test_cap_may_touch_clip_is_conservative() {
  constexpr int W = 288, H = 144;
  hs::random().seed(0xCA9C);

  const int clips[][5] = {
      {0, H / 2, 0, W / 2, 0},  {0, H / 2, W / 2, W, 2},
      {H / 2, H, 0, W / 2, 0},  {H / 2, H, W / 2, W, 1},
      {20, H - 20, 40, 200, 0}, {0, H, 0, W, 0},
  };

  int rejects = 0, reaches = 0;
  for (const auto &bounds : clips) {
    ClipRegion cr;
    cr.w = W;
    cr.h = H;
    cr.y_start = bounds[0];
    cr.y_end = bounds[1];
    cr.x_start = bounds[2];
    cr.x_end = bounds[3];
    cr.margin = bounds[4];

    for (int trial = 0; trial < 300; ++trial) {
      const float dx = hs::rand_f(-1, 1);
      const float dy = hs::rand_f(-1, 1);
      const float dz = hs::rand_f(-1, 1);
      math::Vector d(dx, dy, dz);
      if (d.length() < 0.1f)
        continue;
      const math::Vector dir = d.normalized();
      const float half_angle = hs::rand_f(0.01f, 1.2f);
      const bool passed = Plot::cap_may_touch_clip<H>(cr, dir, half_angle);
      if (!passed)
        rejects++;

      const math::Basis basis = basis_from_normal(dir);
      bool reached = false;
      for (int i = 0; i <= 24 && !reached; ++i) {
        const float t = half_angle * static_cast<float>(i) / 24.0f;
        for (int j = 0; j < 64; ++j) {
          const float az = 2.0f * math::PI_F * static_cast<float>(j) / 64.0f;
          const math::Vector rim = basis.u * cosf(az) + basis.w * sinf(az);
          const math::Vector p = (dir * cosf(t) + rim * sinf(t)).normalized();
          const float row =
              math::phi_to_y<H>(acosf(hs::clamp(p.y, -1.0f, 1.0f)));
          float lam = atan2f(p.z, p.x);
          if (lam < 0.0f)
            lam += 2.0f * math::PI_F;
          const int col =
              std::min(W - 1, static_cast<int>(lam * W / (2.0f * math::PI_F)));
          if (row >= cr.render_y_start() && row < cr.render_y_end() &&
              cr.contains_x(col)) {
            reached = true;
            break;
          }
        }
      }
      if (reached) {
        reaches++;
        HS_EXPECT_TRUE(passed);
      }
    }
  }

  HS_EXPECT_GT(rejects, 100);
  HS_EXPECT_GT(reaches, 100);
}
