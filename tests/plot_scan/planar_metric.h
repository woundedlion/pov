/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Azimuthal-equidistant projection + dual-metric planar arc length, pinned
// against independent oracles (libm great-circle reconstruction, fine on-sphere
// quadrature).
// ============================================================================

/**
 * @brief Full-precision azimuthal unprojection (libm), an oracle independent of
 *        plot.h's LUT-based fast-trig path.
 */
inline math::Vector az_unproject_exact(float Px, float Py,
                                       const math::Basis &b) {
  float R = std::sqrt(Px * Px + Py * Py);
  if (R < math::EPS_GEOMETRIC)
    return b.v;
  float th = std::atan2(Py, Px);
  math::Vector axis = b.u * std::cos(th) + b.w * std::sin(th);
  return b.v * std::cos(R) + axis * std::sin(R);
}

/**
 * @brief Great-circle angle between two unit vectors, accurate for tiny angles.
 * @details atan2(|p x q|, p.q) stays well-conditioned near 0 where acos is
 *          flat; fast_acos (via angle_between) collapses sub-milliradian steps
 *          to zero.
 */
inline float az_arc_exact(const math::Vector &p, const math::Vector &q) {
  return std::atan2(math::cross(p, q).length(), math::dot(p, q));
}

/**
 * @brief azimuthal_project's radius equals the great-circle angle from center.
 * @details Radius is checked against an independent libm great-circle angle.
 */
inline void test_azimuthal_project_radius_is_geodesic_angle() {
  hs::random().seed(0xA21E);
  int mid = 0;
  for (int trial = 0; trial < 4000; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());
    math::Vector p = rand_unit();
    float geo = az_arc_exact(p, basis.v);
    auto proj = Plot::azimuthal_project(p, basis);
    float r = std::hypot(proj.first, proj.second);
    HS_EXPECT_NEAR(r, geo, 5e-3f * (geo + 1.0f));
    if (geo > 0.3f && geo < math::PI_F - 0.3f)
      ++mid;
  }
  HS_EXPECT_GT(mid, 1000);
}

/** Relative allowance on a plane->sphere->plane roundtrip, scaled by R + 1. */
constexpr float AZ_ROUNDTRIP_REL_TOL = 2e-2f;

/**
 * @brief azimuthal_project and azimuthal_unproject invert each other.
 * @details plane->sphere->plane and sphere->plane->sphere both return the
 *          input, away from the antipodal band where the azimuth is unstable.
 */
inline void test_azimuthal_roundtrip_identity() {
  hs::random().seed(0xB33F);
  int fwd = 0;
  for (int trial = 0; trial < 4000; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());

    float R = hs::rand_f(0.05f, math::PI_F - 0.05f);
    float th = hs::rand_f(-math::PI_F, math::PI_F);
    float Px = R * std::cos(th), Py = R * std::sin(th);
    math::Vector s = Plot::azimuthal_unproject(Px, Py, basis);
    auto rp = Plot::azimuthal_project(s, basis);
    HS_EXPECT_NEAR(rp.first, Px, AZ_ROUNDTRIP_REL_TOL * (R + 1.0f));
    HS_EXPECT_NEAR(rp.second, Py, AZ_ROUNDTRIP_REL_TOL * (R + 1.0f));

    math::Vector p = rand_unit();
    if (math::dot(p, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;
    auto proj = Plot::azimuthal_project(p, basis);
    math::Vector back =
        Plot::azimuthal_unproject(proj.first, proj.second, basis);
    HS_EXPECT_NEAR(math::angle_between(p, back), 0.0f, 1.5e-2f);
    ++fwd;
  }
  HS_EXPECT_GT(fwd, 3000);
}

/**
 * @brief azimuthal_unproject lands on the great-circle point at (R, theta).
 * @details Oracle is an independent libm reconstruction
 *          v*cos(R) + (u*cos(th)+w*sin(th))*sin(R).
 */
inline void test_azimuthal_unproject_hits_great_circle_point() {
  hs::random().seed(0xC0DE);
  int n = 0;
  for (int trial = 0; trial < 4000; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());
    float R = hs::rand_f(0.02f, math::PI_F - 0.02f);
    float th = hs::rand_f(-math::PI_F, math::PI_F);
    math::Vector got =
        Plot::azimuthal_unproject(R * std::cos(th), R * std::sin(th), basis);
    math::Vector axis = basis.u * std::cos(th) + basis.w * std::sin(th);
    math::Vector want = basis.v * std::cos(R) + axis * std::sin(R);
    HS_EXPECT_NEAR(math::angle_between(got, want), 0.0f, 1e-2f);
    // The unprojection landed off both poles of the chart, so the oracle
    // compared a point with a defined azimuth.
    float got_R = math::angle_between(got, basis.v);
    if (got_R > 0.05f && got_R < math::PI_F - 0.05f)
      ++n;
  }
  HS_EXPECT_GT(n, 3900);
}

/**
 * @brief planar_arc_length matches a fine libm quadrature of the edge.
 * @details Compares the 4-panel table against a 2000-panel libm reference on
 *          short polygon edges; most edges must bow past their great-circle
 *          chord. Long edges sweeping near the chart center straddle the azimuth
 *          singularity and are excluded.
 */
inline void test_planar_arc_length_matches_fine_quadrature() {
  hs::random().seed(0xD41A);
  int bows = 0;
  float max_rel_err = 0.0f;
  for (int trial = 0; trial < 3000; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());
    float R1 = hs::rand_f(0.2f, 1.2f), R2 = hs::rand_f(0.2f, 1.2f);
    float t1 = hs::rand_f(-math::PI_F, math::PI_F),
          t2 = t1 + hs::rand_f(0.15f, 0.8f);
    math::Vector a =
        Plot::azimuthal_unproject(R1 * std::cos(t1), R1 * std::sin(t1), basis);
    math::Vector b =
        Plot::azimuthal_unproject(R2 * std::cos(t2), R2 * std::sin(t2), basis);
    if (math::dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
        math::dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;

    auto p1 = Plot::azimuthal_project(a, basis);
    auto p2 = Plot::azimuthal_project(b, basis);
    constexpr int N = 2000;
    float fine = 0.0f;
    math::Vector prev = az_unproject_exact(p1.first, p1.second, basis);
    for (int i = 1; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      math::Vector cur =
          az_unproject_exact(p1.first + (p2.first - p1.first) * t,
                             p1.second + (p2.second - p1.second) * t, basis);
      fine += az_arc_exact(prev, cur);
      prev = cur;
    }
    float got = Plot::planar_arc_length(a, b, basis);
    float geo = az_arc_exact(a, b);

    HS_EXPECT_NEAR(got, fine, 0.05f * fine + 5e-3f);
    if (fine > geo * 1.005f)
      ++bows;
    max_rel_err = std::max(max_rel_err, std::abs(got - fine) / (fine + 1e-4f));
  }
  HS_EXPECT_GT(bows, 500);
  HS_EXPECT_LT(max_rel_err, 0.1f);
}

/**
 * @brief The dual metric: radial edges are isometric, azimuthal edges bow.
 * @details A constant-azimuth edge is a meridian great circle, so its planar
 *          arc length equals the geodesic angle exactly; a constant-radius edge
 *          bows strictly past its chord, and the bow grows with radius as the
 *          azimuthal stretch R/sin R rises.
 */
inline void test_dual_metric_radial_vs_azimuthal() {
  hs::random().seed(0xE1A5);
  int radial = 0, azi = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());

    float th = hs::rand_f(-math::PI_F, math::PI_F);
    float Ra = hs::rand_f(0.1f, 0.6f), Rb = hs::rand_f(0.8f, 1.4f);
    math::Vector a =
        Plot::azimuthal_unproject(Ra * std::cos(th), Ra * std::sin(th), basis);
    math::Vector b =
        Plot::azimuthal_unproject(Rb * std::cos(th), Rb * std::sin(th), basis);
    HS_EXPECT_NEAR(Plot::planar_arc_length(a, b, basis),
                   math::angle_between(a, b), 1.2e-2f);
    // The two radii unprojected to a genuine radial edge, not a collapsed one.
    if (math::angle_between(a, b) > 0.1f)
      ++radial;

    float a1 = hs::rand_f(-math::PI_F, math::PI_F), a2 = a1 + 2.4f;
    float chord = 0.0f;
    auto bow = [&](float rad) {
      math::Vector p = Plot::azimuthal_unproject(rad * std::cos(a1),
                                                 rad * std::sin(a1), basis);
      math::Vector q = Plot::azimuthal_unproject(rad * std::cos(a2),
                                                 rad * std::sin(a2), basis);
      chord = math::angle_between(p, q);
      return Plot::planar_arc_length(p, q, basis) - chord;
    };
    float lo = bow(1.0f), hi = bow(1.3f);
    HS_EXPECT_GT(lo, 1.5e-2f);
    HS_EXPECT_GT(hi, lo);
    // The azimuthal separation survived the unprojection.
    if (chord > 0.1f)
      ++azi;
  }
  HS_EXPECT_GT(radial, 1900);
  HS_EXPECT_GT(azi, 1900);
}

/**
 * @brief planar_arc_cumul is monotone and totals what the rasterizer walks.
 * @details Starts at 0, rises strictly, totals what the span-based edge
 *          sampler accumulates from its cached interior points, and is bounded
 *          below by the geodesic angle.
 */
inline void test_planar_arc_cumul_monotone_and_endpoints() {
  hs::random().seed(0xF00D);
  int checked = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    math::Basis basis = basis_from_normal(rand_unit());
    float R1 = hs::rand_f(0.1f, 1.3f), R2 = hs::rand_f(0.1f, 1.3f);
    float t1 = hs::rand_f(-math::PI_F, math::PI_F),
          t2 = t1 + hs::rand_f(0.3f, 1.5f);
    math::Vector a =
        Plot::azimuthal_unproject(R1 * std::cos(t1), R1 * std::sin(t1), basis);
    math::Vector b =
        Plot::azimuthal_unproject(R2 * std::cos(t2), R2 * std::sin(t2), basis);
    if (math::dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
        math::dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;

    auto p1 = Plot::azimuthal_project(a, basis);
    auto p2 = Plot::azimuthal_project(b, basis);
    std::array<float, Plot::PLANAR_LEN_SAMPLES + 1> cumul;
    Plot::planar_arc_cumul(p1, p2.first - p1.first, p2.second - p1.second,
                           basis, cumul);

    HS_EXPECT_NEAR(cumul[0], 0.0f, 1e-6f);
    for (int k = 1; k <= Plot::PLANAR_LEN_SAMPLES; ++k)
      HS_EXPECT_GT(cumul[k], cumul[k - 1]);
    // The per-segment sampler rebuilds the table from the cull span's cached
    // interior points, and the pre-pass takes planar_arc_length. All three must
    // total the same length.
    const Plot::PlanarEdgeSpan span = Plot::make_planar_edge_span(a, b, basis);
    const math::Vector span_end = Plot::azimuthal_unproject(
        span.p1.first + span.dX, span.p1.second + span.dY, basis);
    const Plot::PlanarEdgeSampler sampler =
        Plot::make_planar_edge_sampler(span, span_end, basis);
    constexpr float total_tol = 2e-5f;
    HS_EXPECT_NEAR(cumul[Plot::PLANAR_LEN_SAMPLES], sampler.dist, total_tol);
    HS_EXPECT_NEAR(Plot::planar_arc_length(a, b, basis), sampler.dist,
                   total_tol);
    // Spherical triangle inequality against the endpoints the table actually
    // joins: a total below their separation violates this lower bound.
    const math::Vector chart_start =
        Plot::azimuthal_unproject(p1.first, p1.second, basis);
    const float endpoint_cos = math::dot(chart_start, span_end) /
                               sqrtf(math::dot(chart_start, chart_start) *
                                     math::dot(span_end, span_end));
    constexpr float CHORD_ERROR_BUDGET = Plot::PLANAR_LEN_SAMPLES * 5.1e-5f;
    HS_EXPECT_GE(cumul[Plot::PLANAR_LEN_SAMPLES],
                 acosf(hs::clamp(endpoint_cos, -1.0f, 1.0f)) -
                     CHORD_ERROR_BUDGET);
    ++checked;
  }
  HS_EXPECT_GT(checked, 1500);
}
