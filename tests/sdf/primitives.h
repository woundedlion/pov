/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// clamp_phi
// ============================================================================

/** @brief Verifies phi already in [0, π] passes through clamp_phi unchanged. */
inline void test_clamp_phi_in_range() {
  HS_EXPECT_NEAR(SDF::clamp_phi(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(0.5f), 0.5f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(math::PI_F), math::PI_F, 1e-6f);
}

/** @brief Verifies negative phi reflects across the north pole (|phi|). */
inline void test_clamp_phi_negative_reflects() {
  HS_EXPECT_NEAR(SDF::clamp_phi(-0.3f), 0.3f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(-1.2f), 1.2f, 1e-6f);
}

/** @brief Verifies phi above π reflects across the south pole (2π - phi). */
inline void test_clamp_phi_above_pi_reflects() {
  HS_EXPECT_NEAR(SDF::clamp_phi(math::PI_F + 0.2f), math::PI_F - 0.2f, 1e-5f);
  HS_EXPECT_NEAR(SDF::clamp_phi(2.0f * math::PI_F), 0.0f, 1e-5f);
}

/**
 * @brief Verifies inputs outside [-π, 2π] still fold into [0, π] (full-range
 *        acosf(cosf(x)) equivalence).
 */
inline void test_clamp_phi_full_range() {
  // 2π + 0.2 folds to 0.2.
  HS_EXPECT_NEAR(SDF::clamp_phi(2.0f * math::PI_F + 0.2f), 0.2f, 1e-5f);
  // acosf(cosf(3π)) = π.
  HS_EXPECT_NEAR(SDF::clamp_phi(3.0f * math::PI_F), math::PI_F, 1e-5f);
  // -1.5π folds to 0.5π.
  HS_EXPECT_NEAR(SDF::clamp_phi(-1.5f * math::PI_F), 0.5f * math::PI_F, 1e-5f);
}

/**
 * @brief Verifies clamp_phi_band reports the exact colatitude extent of the
 *        circle it bounds, against a brute-force sweep of that circle.
 * @details The extremes sit at psi = 0 and psi = π, both sampled exactly.
 */
inline void test_clamp_phi_band_matches_circle_extent() {
  constexpr int SAMPLES = 512;
  auto circle_extent = [](float c, float t) {
    float lo = math::PI_F, hi = 0.0f;
    for (int i = 0; i < SAMPLES; ++i) {
      float psi = math::TWO_PI_F * i / SAMPLES;
      float y =
          std::cos(t) * std::cos(c) - std::sin(t) * std::sin(c) * std::cos(psi);
      float phi = std::acos(std::max(-1.0f, std::min(1.0f, y)));
      lo = std::min(lo, phi);
      hi = std::max(hi, phi);
    }
    return std::pair<float, float>(lo, hi);
  };

  const float centers[] = {0.0f, 0.2f, 0.9f, 1.5f, 2.4f, 3.0f, math::PI_F};
  // Past π and negative: a complement-radius ring and a mirrored one still
  // trace a real circle, so the band must order its folded endpoints.
  const float radii[] = {0.0f, 0.1f,       0.5f, 1.0f, 2.0f,
                         2.5f, math::PI_F, 4.0f, -0.7f};
  for (int ci = 0; ci < static_cast<int>(std::size(centers)); ++ci) {
    for (int ti = 0; ti < static_cast<int>(std::size(radii)); ++ti) {
      HS_CONTEXT("center/radius index", ci, ti);
      const float c = centers[ci];
      const float t = radii[ti];
      SDF::PhiBand band = SDF::clamp_phi_band(c, t);
      auto expected = circle_extent(c, t);
      HS_EXPECT_NEAR(band.phi_min, expected.first, 1e-4f);
      HS_EXPECT_NEAR(band.phi_max, expected.second, 1e-4f);
      HS_EXPECT_LE(band.phi_min, band.phi_max);
      HS_EXPECT_GE(band.phi_min, 0.0f);
      HS_EXPECT_LE(band.phi_max, math::PI_F);
    }
  }
}

/**
 * @brief Verifies the two degenerate poses report their real band rather than
 *        opening to a pole: a pole-axis circle collapses to one latitude, and a
 *        circle whose far edge runs past the south pole reflects back.
 */
inline void test_clamp_phi_band_pole_crossing_poses() {
  SDF::PhiBand pole_axis = SDF::clamp_phi_band(0.0f, 0.5f);
  HS_EXPECT_NEAR(pole_axis.phi_min, 0.5f, 1e-5f);
  HS_EXPECT_NEAR(pole_axis.phi_max, 0.5f, 1e-5f);

  SDF::PhiBand wrapped = SDF::clamp_phi_band(3.0f, 2.5f);
  HS_EXPECT_NEAR(wrapped.phi_min, 0.5f, 1e-5f);
  HS_EXPECT_NEAR(wrapped.phi_max, math::TWO_PI_F - 5.5f, 1e-5f);
}

/** @brief Verifies the reciprocal sector fold matches the general wrap path. */
inline void test_centered_sector_angle_matches_wrap() {
  const float angles[] = {-13.7f, -2.3f, -0.7f, 0.21f, 1.8f, 9.4f};
  for (int sides : {3, 5, 8, 17}) {
    float sector = 2.0f * math::PI_F / sides;
    float reciprocal_sector = static_cast<float>(sides) / (2.0f * math::PI_F);
    for (float angle : angles) {
      const float folded =
          SDF::centered_sector_angle(angle, sector, reciprocal_sector);
      float expected =
          math::wrap(angle + sector * 0.5f, sector) - sector * 0.5f;
      HS_EXPECT_NEAR(folded, expected, 1e-5f);
      HS_EXPECT_LE(folded, sector * 0.5f + 1e-5f);
      HS_EXPECT_GE(folded, -sector * 0.5f - 1e-5f);
      const float turns = (angle - folded) / sector;
      HS_EXPECT_NEAR(turns, std::round(turns), 1e-4f);
    }
  }
}

// ============================================================================
// Ring
// ============================================================================

/** @brief Ring distances remain finite at the exact center. */
inline void test_ring_roundoff_at_exact_center() {
  auto basis = equator_basis();
  basis.v = math::Vector(0.0f, std::nextafter(1.0f, 2.0f), 0.0f);
  for (float radius : {0.0f, 2.0f}) {
    const float direction = radius == 0.0f ? 1.0f : -1.0f;
    SDF::Ring dot(basis, radius, 0.1f);
    const auto p = math::Vector(0.0f, direction, 0.0f);
    const auto result = SDF::distance_of(dot, p);
    HS_EXPECT_NEAR(result.dist, -0.1f, 1e-6f);
    HS_EXPECT_NEAR(result.raw_dist, 0.0f, 1e-6f);
    HS_EXPECT_NEAR(dot.stroke_alpha(math::dot(p, basis.v)), 1.0f, 1e-6f);
  }
}

/** @brief Verifies a point on the ring centerline reads raw_dist 0 and dist = -thickness. */
inline void test_ring_on_centerline() {
  math::Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  auto r = SDF::distance_of(ring, math::Vector(1, 0, 0));
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
}

/** @brief Verifies a point within the ring band reads negative dist. */
inline void test_ring_inside_band() {
  math::Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  float off = 0.05f;
  math::Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = SDF::distance_of(ring, p);
  HS_EXPECT_TRUE(r.dist < 0.0f);
  HS_EXPECT_NEAR(r.raw_dist, 0.05f, 1e-3f);

  SDF::Ring thin_ring(b, 0.5f, 0.05f);
  for (float offset : {-0.03f, 0.03f}) {
    const float ANGLE = math::PI_F / 4.0f + offset;
    const auto SAMPLE = SDF::distance_of(
        thin_ring, math::Vector(sinf(ANGLE), cosf(ANGLE), 0.0f));
    HS_EXPECT_NEAR(SAMPLE.raw_dist, 0.03f, 1e-3f);
  }
}

/** @brief Verifies a point far outside the band reads the cull sentinel rather than a real dist. */
inline void test_ring_outside_band_returns_sentinel() {
  math::Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  auto r = SDF::distance_of(ring, math::Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist > 50.0f);
}

/** @brief Verifies a point just past the band edge trips the sentinel. */
inline void test_ring_just_outside_band() {
  math::Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.05f);

  // off (0.07) > thickness (0.05): point is 0.02 past the band edge.
  float off = 0.07f;
  math::Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = SDF::distance_of(ring, p);
  HS_EXPECT_TRUE(r.dist > 50.0f);
}

/**
 * @brief Verifies a small-radius ring thick enough to break the linearized
 *        centerline distance still reports the geodesic distance,
 *        symmetric about the centerline.
 */
inline void test_ring_small_radius_distance_symmetric() {
  const float RADIUS = 0.04f;
  const float THICKNESS = 0.05f;
  math::Basis b = equator_basis();
  SDF::Ring ring(b, RADIUS, THICKNESS);
  float target = RADIUS * (math::PI_F / 2.0f);

  // Ring axis is +Y, so a point at polar angle a off the axis is
  // (sin a, cos a, 0).
  auto at_polar = [](float a) {
    return math::Vector(std::sin(a), std::cos(a), 0.0f);
  };

  for (float off : {-0.9f, -0.4f, 0.4f, 0.9f}) {
    float delta = off * THICKNESS;
    auto r = SDF::distance_of(ring, at_polar(target + delta));
    HS_EXPECT_NEAR(r.raw_dist, std::abs(delta), 1e-4f);
    HS_EXPECT_NEAR(r.dist, std::abs(delta) - THICKNESS, 1e-4f);
  }

  // The AA ramp reads the same either side of the centerline.
  float d_in = math::dot(at_polar(target - 0.9f * THICKNESS), b.v);
  float d_out = math::dot(at_polar(target + 0.9f * THICKNESS), b.v);
  HS_EXPECT_NEAR(ring.stroke_alpha(d_in), ring.stroke_alpha(d_out), 2e-4f);
}

// ============================================================================
// DistortedRing  (per-azimuth centerline shift)
// ============================================================================

/**
 * @brief Verifies a constant shift_fn moves the centerline by exactly that
 *        offset, and that max_distortion widens the early-reject band so the
 *        shifted on-centerline point is not falsely culled.
 */
inline void test_distorted_ring_constant_shift_moves_centerline() {
  math::Basis b =
      equator_basis(); // v=+Y, u=+X, w=+Z; radius=1 → target_angle = π/2
  const float shift = 0.2f;
  const float thickness = 0.05f;

  // Azimuth 0 (along +X) on the shifted centerline: polar angle from +Y is
  // π/2 + shift.
  math::Vector p(std::sin(math::PI_F / 2 + shift),
                 std::cos(math::PI_F / 2 + shift), 0.0f);

  SDF::DistortedRing shifted(
      b, 1.0f, thickness, [shift](float) { return shift; },
      /*max_distortion=*/shift, /*phase=*/0.0f);
  auto rs = SDF::distance_of(shifted, p);
  HS_EXPECT_TRUE(rs.dist < 50.0f);
  HS_EXPECT_NEAR(rs.raw_dist, 0.0f, 5e-4f);
  HS_EXPECT_NEAR(rs.dist, -thickness, 5e-4f);
  HS_EXPECT_NEAR(rs.t, 0.0f, 5e-4f);

  // Same point, no shift: the centerline stays at π/2, so it now sits `shift`
  // radians off (raw_dist ≈ shift).
  SDF::DistortedRing plain(
      b, 1.0f, thickness, [](float) { return 0.0f; },
      /*max_distortion=*/shift, /*phase=*/0.0f);
  auto rp = SDF::distance_of(plain, p);
  HS_EXPECT_NEAR(rp.raw_dist, shift, 5e-4f);
}

/**
 * @brief Verifies a sinusoidal shift_fn moves the centerline by a per-azimuth
 *        amount, so the t parameter feeding shift_fn is wired correctly.
 */
inline void test_distorted_ring_sin_shift_varies_by_azimuth() {
  math::Basis b = equator_basis();
  const float amp = 0.2f;
  SDF::DistortedRing ring(
      b, 1.0f, 0.05f,
      [amp](float t) { return amp * std::sin(2 * math::PI_F * t); }, amp, 0.0f);

  // Azimuth π/2 (along +Z) → t = 0.25 → shift = amp; centerline polar angle is
  // π/2 + amp.
  math::Vector on(0.0f, std::cos(math::PI_F / 2 + amp),
                  std::sin(math::PI_F / 2 + amp));
  auto r_on = SDF::distance_of(ring, on);
  HS_EXPECT_NEAR(r_on.t, 0.25f, 5e-4f);
  HS_EXPECT_NEAR(r_on.raw_dist, 0.0f, 5e-4f);

  // Same azimuth on the unshifted equator (+Z): centerline moved by amp here.
  auto r_off = SDF::distance_of(ring, math::Vector(0, 0, 1));
  HS_EXPECT_NEAR(r_off.t, 0.25f, 5e-4f);
  HS_EXPECT_NEAR(r_off.raw_dist, amp, 5e-4f);
}

/**
 * @brief Verifies flat mode preserves the exact zero-knot polar distance.
 */
inline void test_distorted_ring_flat_matches_zero_knots() {
  constexpr int LUT_N = 16;
  float knots[LUT_N + 1] = {};

  auto check = [&](const math::Basis &basis, float radius) {
    constexpr float thickness = 0.08f;
    SDF::FlatDistortedRing flat(basis, radius, thickness);
    SDF::KnotPrefilter pf;
    SDF::DistortedRing polyline(basis, radius, thickness, knots, LUT_N, 0.0f,
                                pf);
    const float target = radius * (math::PI_F / 2.0f);
    const float azimuths[] = {0.0f, 1e-5f, math::PI_F / 2.0f, math::PI_F,
                              2.0f * math::PI_F - 1e-5f};
    const float offsets[] = {-0.04f, 0.0f, 0.04f};
    for (float azimuth : azimuths) {
      for (float offset : offsets) {
        float polar = hs::clamp(target + offset, 0.0f, math::PI_F);
        math::Vector p =
            basis.v * cosf(polar) +
            (basis.u * cosf(azimuth) + basis.w * sinf(azimuth)) * sinf(polar);
        auto actual = SDF::distance_of(flat, p);
        auto expected = SDF::distance_of(polyline, p);
        HS_EXPECT_NEAR(actual.dist, expected.dist, 1e-5f);
        HS_EXPECT_NEAR(actual.raw_dist, expected.raw_dist, 1e-5f);
        HS_EXPECT_NEAR(actual.t, expected.t, 1e-5f);
      }
    }

    SDF::DistanceResult no_uv;
    flat.distance<false>(basis.v, no_uv);
    HS_EXPECT_EQ(no_uv.t, 0.0f);
  };

  check(math::make_basis(math::Quaternion(), math::Y_AXIS), 0.01f);
  check(math::make_basis(math::Quaternion(), math::X_AXIS), 1.0f);
  check(math::make_basis(math::Quaternion(),
                         math::Vector(0.3f, 0.8f, -0.5f).normalized()),
        1.99f);
}

/**
 * @brief Verifies the knot overload's raw_dist matches a brute-force geodesic
 *        minimum over the densely sampled polyline, within the stroke reach.
 * @details Accuracy is contracted only for true distances below thickness,
 *   the outward search's reach; every probe is placed inside it.
 */
template <int LUT_N> inline void expect_polyline_distance_matches_bruteforce() {
  math::Basis b =
      math::make_basis(math::Quaternion(), math::Vector(0.3f, 1.0f, 0.2f));
  const float amp = 0.2f;
  const int harmonic = 5;
  const float radius = 0.5f; // target_angle = π/4: curved chart, tilted axis
  const float thickness = 0.06f;
  float knots[LUT_N + 1];
  for (int k = 0; k <= LUT_N; ++k)
    knots[k] =
        amp * std::sin(2.0f * math::PI_F * harmonic * (k % LUT_N) / LUT_N);
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, radius, thickness, knots, LUT_N, 0.0f, pf);

  const float target = radius * (math::PI_F / 2.0f);
  auto on_sphere = [&](float t, float dv) {
    float theta =
        target + amp * std::sin(2.0f * math::PI_F * harmonic * t) + dv;
    float a = 2.0f * math::PI_F * t;
    return (b.v * std::cos(theta)) +
           ((b.u * std::cos(a)) + (b.w * std::sin(a))) * std::sin(theta);
  };
  auto brute = [&](const math::Vector &p) {
    constexpr int SAMPLES = LUT_N * 64;
    float best = 100.0f;
    for (int s = 0; s < SAMPLES; ++s) {
      int k = s / 64;
      float f = (s % 64) / 64.0f;
      float theta = target + knots[k] + f * (knots[k + 1] - knots[k]);
      float a = 2.0f * math::PI_F * (k + f) / LUT_N;
      math::Vector q =
          (b.v * std::cos(theta)) +
          ((b.u * std::cos(a)) + (b.w * std::sin(a))) * std::sin(theta);
      best = std::min(best, std::acos(hs::clamp(math::dot(p, q), -1.0f, 1.0f)));
    }
    return best;
  };

  // t = 0.05: crest (slope zero, curvature max); t = 0.1: steep flank;
  // t = 0.998: wrap seam. Offsets stay within the thickness reach.
  const float probes[][2] = {{0.05f, 0.04f},   {0.05f, -0.05f}, {0.1f, 0.05f},
                             {0.1f, -0.04f},   {0.998f, 0.04f}, {0.25f, 0.05f},
                             {0.375f, -0.05f}, {0.6f, 0.03f}};
  for (const auto &pr : probes) {
    math::Vector p = on_sphere(pr[0], pr[1]);
    auto r = SDF::distance_of(ring, p);
    float expected = brute(p);
    HS_EXPECT_TRUE(expected < thickness);
    HS_EXPECT_NEAR(r.raw_dist, expected, 3e-3f);
  }
}

/** @brief Closing-segment bounds depend on the first knot, not a sentinel. */
inline void test_distorted_ring_closes_without_a_sentinel() {
  const math::Basis basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  const float knots[] = {0.1f, -0.2f, 0.3f, 0.05f};
  const float stale_sentinel[] = {0.1f, -0.2f, 0.3f, 0.05f, 10.0f};
  SDF::KnotPrefilter exact, extra;
  SDF::DistortedRing a(basis, 0.5f, 0.02f, knots, 4, 0.0f, exact);
  SDF::DistortedRing b(basis, 0.5f, 0.02f, stale_sentinel, 4, 0.0f, extra);
  for (int c = 0; c < SDF::KnotPrefilter::CHUNKS; ++c) {
    HS_EXPECT_EQ(exact.lo[c], extra.lo[c]);
    HS_EXPECT_EQ(exact.hi[c], extra.hi[c]);
  }
  HS_EXPECT_EQ(a.max_distortion, b.max_distortion);
}

/** @brief Compares distance at aligned and chunk-straddling knot counts. */
inline void test_distorted_ring_polyline_distance_matches_bruteforce() {
  expect_polyline_distance_matches_bruteforce<96>();
  expect_polyline_distance_matches_bruteforce<97>();
}

/** @brief Verifies the knot reject band follows asymmetric extrema, not ±max_distortion. */
inline void test_distorted_ring_knot_extrema_tighten_band() {
  constexpr int LUT_N = 8;
  constexpr float RADIUS = 0.8f;
  constexpr float THICKNESS = 0.05f;
  constexpr float TARGET = RADIUS * (math::PI_F / 2.0f);
  float knots[LUT_N + 1] = {0.18f, 0.12f, 0.04f, -0.01f, -0.03f,
                            0.02f, 0.09f, 0.16f, 0.18f};
  math::Basis basis = equator_basis();
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(basis, RADIUS, THICKNESS, knots, LUT_N, 0.0f, pf);

  HS_EXPECT_NEAR(ring.max_distortion, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(ring.max_thickness, 0.23f, 1e-6f);
  for (int k : {0, 4, 8}) {
    float t = static_cast<float>(k % LUT_N) / LUT_N;
    float azimuth = 2.0f * math::PI_F * t;
    float polar = TARGET + knots[k];
    math::Vector p(cosf(azimuth) * sinf(polar), cosf(polar),
                   sinf(azimuth) * sinf(polar));
    HS_EXPECT_NEAR(SDF::distance_of(ring, p).raw_dist, 0.0f, 2e-4f);
  }

  float below = TARGET - 0.03f - THICKNESS - 1e-3f;
  float above = TARGET + 0.18f + THICKNESS + 1e-3f;
  HS_EXPECT_GT(
      SDF::distance_of(ring, math::Vector(sinf(below), cosf(below), 0.0f)).dist,
      50.0f);
  HS_EXPECT_GT(
      SDF::distance_of(ring, math::Vector(sinf(above), cosf(above), 0.0f)).dist,
      50.0f);
}

/**
 * @brief Verifies a knot ring reports the far sentinel past its stroke reach.
 * @details A pixel no segment comes within thickness of has no distance to
 *   report; reporting the reach itself would put dist on the surface.
 */
inline void test_distorted_ring_past_reach_reports_far_sentinel() {
  constexpr int LUT_N = 32;
  constexpr float RADIUS = 0.8f;
  constexpr float THICKNESS = 0.05f;
  constexpr float TARGET = RADIUS * (math::PI_F / 2.0f);
  constexpr float SPIKE = 0.5f;
  // One tall spike at azimuth 0, every other cell on the centerline. A pixel in
  // the next chunk shares the spike's prefilter window, so the segment search
  // runs and terminates without a point inside the reach.
  float knots[LUT_N + 1] = {};
  knots[0] = SPIKE;
  knots[LUT_N] = SPIKE;
  math::Basis basis = equator_basis();
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(basis, RADIUS, THICKNESS, knots, LUT_N, 0.0f, pf);

  const float azimuth = 2.0f * math::PI_F * 1.5f / LUT_N;
  const float polar = TARGET + SPIKE - 0.05f;
  math::Vector p(cosf(azimuth) * sinf(polar), cosf(polar),
                 sinf(azimuth) * sinf(polar));
  HS_EXPECT_GT(SDF::distance_of(ring, p).dist, 50.0f);

  math::Basis poly_basis = math::make_basis(math::Quaternion(), p);
  SDF::PlanarPolygon poly(poly_basis, /*radius=*/0.3f / (math::PI_F / 2.0f),
                          /*sides=*/6, 0.0f);
  SDF::Subtract<SDF::PlanarPolygon, SDF::DistortedRing> carved(poly, ring);
  const float solid = SDF::distance_of(poly, p).dist;
  HS_EXPECT_LT(solid, -0.1f);
  HS_EXPECT_NEAR(SDF::distance_of(carved, p).dist, solid, 1e-6f);
}

/**
 * @brief Verifies distance_from_frame() lights exactly where distance() does,
 *        at the same distance.
 * @details Where either lights (dist < 0) both must, with raw distances equal
 *          to float rounding.
 */
inline void test_distorted_ring_frame_distance_matches_distance() {
  const math::Basis basis = math::make_basis(
      math::Quaternion(), math::Vector(0.3f, 1.0f, 0.2f).normalized());
  size_t lit = 0;
  for (int lut_n : {16, 97, 288}) {
    std::vector<float> knots(lut_n + 1);
    for (int k = 0; k <= lut_n; ++k) {
      const float t = static_cast<float>(k % lut_n) / lut_n;
      knots[k] = 0.08f * sinf(2.0f * math::PI_F * 3.0f * t) +
                 0.03f * cosf(2.0f * math::PI_F * 7.0f * t);
    }
    for (float thickness : {0.01f, 0.03f, 0.12f}) {
      for (float radius : {0.02f, 0.5f, 1.0f, 1.9f}) {
        SDF::KnotPrefilter prefilter;
        const SDF::DistortedRing ring(basis, radius, thickness, knots.data(),
                                      lut_n, 0.0f, prefilter);
        const float target = radius * (math::PI_F / 2.0f);
        for (int i = 0; i <= 64; ++i) {
          const float polar =
              hs::clamp(target - 0.25f + 0.5f * i / 64.0f, 0.0f, math::PI_F);
          for (int k = 0; k < 256; ++k) {
            const float a = 2.0f * math::PI_F * (k + 0.37f) / 256.0f;
            const math::Vector p =
                (basis.v * cosf(polar)) +
                ((basis.u * cosf(a)) + (basis.w * sinf(a))) * sinf(polar);
            SDF::DistanceResult full;
            ring.distance<true>(p, full);
            const float d = math::dot(p, ring.normal);
            const float frame_polar =
                math::fast_acos(hs::clamp(d, -1.0f, 1.0f));
            const float t_norm = math::wrap_t(
                SDF::basis_azimuth(p, ring.u, ring.w, 0.0f) / math::TWO_PI_F);
            const float sin_polar = sqrtf(
                std::max(1.0f - d * d, SDF::DistortedRing::POLE_SIN2_FLOOR));
            SDF::DistanceResult frame;
            ring.distance_from_frame(d, frame_polar, sin_polar, t_norm, frame);
            if (full.dist < 0.0f || frame.dist < 0.0f) {
              HS_CONTEXT("lut_n / sample", lut_n, (i << 8) | k);
              HS_EXPECT_LT(full.dist, 0.0f);
              HS_EXPECT_LT(frame.dist, 0.0f);
              HS_EXPECT_NEAR(frame.raw_dist, full.raw_dist, 1e-6f);
              ++lit;
            }
          }
        }
      }
    }
  }
  HS_EXPECT_GT(lit, static_cast<size_t>(10000));
}

// ============================================================================
// PlanarPolygon  (Basis at top of sphere; distance to nearest edge)
// ============================================================================

/** @brief Verifies the polygon center is inside, with dist equal to the negated apothem. */
inline void test_polygon_at_center_inside() {
  math::Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, /*radius*/ 0.5f / (math::PI_F / 2.0f), /*sides*/ 6,
                          /*phase*/ 0.0f);

  auto r = SDF::distance_of(poly, math::Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
  float apothem = 0.5f * std::cos(math::PI_F / 6.0f);
  HS_EXPECT_NEAR(r.dist, -apothem, 1e-3f);
}

/** @brief Verifies the antipode of the polygon center is outside (positive dist). */
inline void test_polygon_far_point_outside() {
  math::Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, 0.3f / (math::PI_F / 2.0f), 6, 0.0f);

  auto r = SDF::distance_of(poly, math::Vector(0, -1, 0));
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

// ============================================================================
// SphericalPolygon — great-circle edges
// ============================================================================

/** @brief Verifies the spherical-polygon center is strictly inside. */
inline void test_spherical_polygon_center_inside() {
  math::Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, /*radius*/ 0.5f, /*sides*/ 5, /*phase*/ 0.0f);
  auto r = SDF::distance_of(sp, math::Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies the antipode of the spherical-polygon center is outside. */
inline void test_spherical_polygon_far_outside() {
  math::Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, 0.3f, 6, 0.0f);
  auto r = SDF::distance_of(sp, math::Vector(0, -1, 0));
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

/**
 * @brief Verifies the spherical-polygon center reads the negated geodesic
 *        inradius, and a point on an edge midpoint reads dist 0 there.
 * @details Regular-spherical-polygon inradius r and circumradius R relate by
 *   tan(r) = tan(R)·cos(π/n) (Napier, right triangle center-midpoint-vertex).
 *   The edge bisector lies along +u, so the point at polar angle r on +u sits on
 *   an edge great circle: dist 0, raw_dist = r.
 */
inline void test_spherical_polygon_center_and_edge_magnitude() {
  math::Basis b = equator_basis();
  const int sides = 5;
  const float radius = 0.5f;
  SDF::SphericalPolygon sp(b, radius, sides, 0.0f);

  const float R = radius * (math::PI_F / 2.0f);
  const float inradius = std::atan(std::tan(R) * std::cos(math::PI_F / sides));

  auto center = SDF::distance_of(sp, math::Vector(0, 1, 0));
  HS_EXPECT_NEAR(center.dist, -inradius, 5e-4f);
  HS_EXPECT_NEAR(center.raw_dist, 0.0f, 1e-3f);

  // Edge midpoint: polar angle = inradius along +u (the sector bisector).
  math::Vector edge_mid(std::sin(inradius), std::cos(inradius), 0.0f);
  auto em = SDF::distance_of(sp, edge_mid);
  HS_EXPECT_NEAR(em.dist, 0.0f, 5e-4f);
  HS_EXPECT_NEAR(em.raw_dist, inradius, 5e-4f);
}

/**
 * @brief Bounds sine-domain distance error across the device-width AA band.
 * @details Where the edge dot wins, the paths differ only by sin(x) - x. Where
 *   the circumscribed-disc clamp wins, distance() uses fast_acos, so the gap
 *   widens to fast_acos's error.
 */
inline void test_spherical_polygon_sine_distance_aa_error() {
  constexpr int W = 288;
  constexpr int H = 144;
  constexpr float PIXEL_WIDTH = 2.0f * math::PI_F / W;
  math::Basis basis = math::make_basis(
      math::make_rotation(math::Vector(0.3f, -0.8f, 0.5f).normalized(), 0.71f),
      math::Y_AXIS);
  struct Case {
    float radius;
    int sides;
    float phase;
  };
  const Case cases[] = {{0.22f, 3, -2.7f},
                        {0.72f, 5, 0.31f},
                        {0.98f, 12, 5.4f},
                        {1.42f, 7, -0.9f},
                        {1.999f, 5, 0.0f}};

  float max_error = 0.0f;
  int edge_samples = 0;
  for (const Case &c : cases) {
    auto folded = math::get_antipode(basis, c.radius);
    SDF::SphericalPolygon shape(folded.first, folded.second, c.sides, c.phase,
                                c.radius > 1.0f);
    for (int y = 0; y < H; ++y) {
      float polar = math::PI_F * (static_cast<float>(y) + 0.5f) / H;
      float sin_p = sinf(polar);
      float cos_p = cosf(polar);
      for (int x = 0; x < W; ++x) {
        float azimuth = 2.0f * math::PI_F * (static_cast<float>(x) + 0.5f) / W;
        math::Vector p(sin_p * cosf(azimuth), cos_p, sin_p * sinf(azimuth));
        SDF::DistanceResult exact;
        shape.distance<false>(p, exact);
        float sine = shape.sine_distance(p);
        HS_EXPECT_TRUE(exact.dist == 0.0f || sine == 0.0f ||
                       std::signbit(exact.dist) == std::signbit(sine));
        if (std::abs(exact.dist) <= PIXEL_WIDTH) {
          max_error = fold_worst(max_error, std::abs(exact.dist - sine));
          ++edge_samples;
        }
      }
    }
  }

  std::printf("spherical sine AA samples=%d max_error=%.9g rad\n", edge_samples,
              max_error);
  HS_EXPECT_GT(edge_samples, 1000);
  HS_EXPECT_LE(max_error, 1.5e-4f);
}

/** @brief Sine-domain polygon distance retains the full interior. */
inline void test_spherical_polygon_sine_full_interior() {
  constexpr int W = 288;
  constexpr int H = 144;
  constexpr float PIXEL_WIDTH = 2.0f * math::PI_F / W;
  const math::Basis basis = equator_basis();
  for (float radius : {0.0001f, 0.001f, 0.01f, 0.5f, 1.0f}) {
    for (bool invert : {false, true}) {
      SDF::SphericalPolygon shape(basis, radius, 5, 0.0f, invert);
      for (float angle : {0.0f, 0.1f, 0.3f}) {
        const math::Vector point(std::sin(angle), -std::cos(angle), 0.0f);
        const float distance = shape.sine_distance(point);
        HS_EXPECT_EQ(Scan::solid_coverage(distance, PIXEL_WIDTH),
                     invert ? 1.0f : 0.0f);
      }
    }
  }
  StubEffect effect(W, H);
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(effect);
    Scan::SphericalPolygon::draw_solid<W, H, true>(
        pipeline, canvas, basis, 1.999f, 5,
        Color4(Pixel(60000, 60000, 60000), 1.0f));
  }
  effect.advance_display();
  for (int y = 0; y < H / 2; ++y)
    for (int x = 0; x < W; ++x)
      HS_EXPECT_EQ(effect.get_pixel(x, y).r, 60000);
}

/**
 * @brief Verifies SphericalPolygon satisfies SDFShape and composes under CSG.
 * @details Two disjoint polygons must emit their two arcs through the
 *   span-bounded scanline path.
 */
inline void test_spherical_polygon_composes_under_csg() {
  static_assert(SDF::SDFShape<SDF::SphericalPolygon>,
                "SphericalPolygon must satisfy the CSG child contract");
  using U = SDF::Union<SDF::SphericalPolygon, SDF::SphericalPolygon>;
  static_assert(SDF::sdf_max_spans<SDF::SphericalPolygon>::value == 1);
  static_assert(SDF::sdf_max_spans<U>::value == 2);
  static_assert(SDF::sdf_max_spans<SDF::Intersection<U, U>>::value == 8);

  math::Basis b = equator_basis();
  SDF::SphericalPolygon inner(b, 0.3f, 5, 0.0f);
  SDF::SphericalPolygon outer(b, 0.7f, 5, 0.0f);

  U coaxial(inner, outer);
  // Sector bisector (+u) at polar 0.7: past inner's inradius, short of outer's.
  math::Vector p(std::sin(0.7f), std::cos(0.7f), 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(inner, p).dist > 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(outer, p).dist < 0.0f);
  HS_EXPECT_NEAR(SDF::distance_of(coaxial, p).dist,
                 SDF::distance_of(outer, p).dist, 1e-6f);

  // Poles on +X and +Z: both cross the equatorial row, a quarter turn apart.
  constexpr int W = 256, H = 128;
  const math::Basis bx{math::Vector(0, 1, 0), math::Vector(1, 0, 0),
                       math::Vector(0, 0, 1)};
  const math::Basis bz{math::Vector(0, 1, 0), math::Vector(0, 0, 1),
                       math::Vector(-1, 0, 0)};
  SDF::SphericalPolygon px(bx, 0.3f, 5, 0.0f), pz(bz, 0.3f, 5, 0.0f);
  U disjoint(px, pz);

  std::vector<std::pair<float, float>> out;
  bool ok = disjoint.get_horizontal_intervals<W, H>(
      H / 2, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
}

// ============================================================================
// Star
// ============================================================================

/** @brief Verifies the star center is interior (negative dist). */
inline void test_star_center_inside() {
  math::Basis b = equator_basis();
  SDF::Star star(b, /*radius*/ 0.6f, /*sides*/ 5, /*phase*/ 0.0f);
  auto r = SDF::distance_of(star, math::Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies the antipode of the star center is outside. */
inline void test_star_far_outside() {
  math::Basis b = equator_basis();
  SDF::Star star(b, 0.4f, 5, 0.0f);
  auto r = SDF::distance_of(star, math::Vector(0, -1, 0));
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

/**
 * @brief Verifies a star point tip lies on the boundary (dist 0) at the outer
 *        radius along a sector bisector.
 * @details A tip sits at polar angle outer_radius = radius·π/2 along +u (with
 *   phase 0 the sector folds azimuth 0 to the tip bisector), so its distance to
 *   the point edge is exactly 0 and raw_dist equals the outer radius.
 */
inline void test_star_tip_on_boundary() {
  math::Basis b = equator_basis();
  const float radius = 0.6f;
  SDF::Star star(b, radius, /*sides=*/5, 0.0f);

  const float outer = radius * (math::PI_F / 2.0f);
  math::Vector tip(std::sin(outer), std::cos(outer), 0.0f);
  auto r = SDF::distance_of(star, tip);
  HS_EXPECT_NEAR(r.dist, 0.0f, 5e-4f);
  HS_EXPECT_NEAR(r.raw_dist, outer, 5e-4f);
}

// ============================================================================
// Flower  (petals radiating from the axis antipode)
// ============================================================================

/**
 * @brief Verifies an interior point on a petal bisector reads dist = scan_dist −
 *        outer radius.
 * @details The flower is scanned from the antipode of basis.v. Along a petal
 *   bisector (+u from the antipode) at scan distance s the petal-edge distance is
 *   s − outer (negative for s < outer). The exact antipode is avoided: azimuth is
 *   undefined at the pole, so the sector fold there is degenerate.
 */
inline void test_flower_interior_along_petal() {
  math::Basis b = equator_basis();
  const float radius = 0.6f;
  SDF::Flower flower(b, radius, /*sides=*/6, 0.0f);

  const float outer = radius * (math::PI_F / 2.0f);
  const float s = 0.1f;
  // From the antipode (-Y), step s toward +u (+X): interior of a petal.
  math::Vector p(std::sin(s), -std::cos(s), 0.0f);
  auto r = SDF::distance_of(flower, p);
  HS_EXPECT_NEAR(r.dist, s - outer, 5e-4f);
  HS_EXPECT_NEAR(r.raw_dist, s, 1e-3f);
}

/**
 * @brief Verifies a petal tip sits on the boundary (dist 0) at the outer radius
 *        along a petal bisector.
 * @details A petal bisector runs along +u from the antipode; at scan distance
 *   outer = radius·π/2 the petal-edge distance is exactly 0 and raw_dist equals
 *   that scan distance.
 */
inline void test_flower_petal_tip_on_boundary() {
  math::Basis b = equator_basis();
  const float radius = 0.6f;
  SDF::Flower flower(b, radius, /*sides=*/6, 0.0f);

  const float outer = radius * (math::PI_F / 2.0f);
  // From the antipode (-Y), step `outer` toward +u (+X).
  math::Vector tip(std::sin(outer), -std::cos(outer), 0.0f);
  auto r = SDF::distance_of(flower, tip);
  HS_EXPECT_NEAR(r.dist, 0.0f, 5e-4f);
  HS_EXPECT_NEAR(r.raw_dist, outer, 5e-4f);
}

/** @brief Verifies the solid-shape unit-vector and no-UV distance paths. */
inline void test_solid_shape_unit_angle_and_no_uv_paths() {
  math::Basis b = math::make_basis(
      math::Quaternion(), math::Vector(0.3f, 0.8f, -0.5f).normalized());
  float polar = 0.73f;
  float azimuth = -1.17f;
  math::Vector p = b.v * cosf(polar) +
                   (b.u * cosf(azimuth) + b.w * sinf(azimuth)) * sinf(polar);

  auto check = [&](const auto &shape, float expected_raw) {
    SDF::DistanceResult with_uv;
    SDF::DistanceResult no_uv;
    shape.template distance<true>(p, with_uv);
    shape.template distance<false>(p, no_uv);
    HS_EXPECT_NEAR(no_uv.dist, with_uv.dist, 1e-6f);
    HS_EXPECT_NEAR(no_uv.raw_dist, expected_raw, 2e-4f);
    HS_EXPECT_EQ(no_uv.t, 0.0f);
  };

  check(SDF::PlanarPolygon(b, 0.8f, 7, -2.4f), math::angle_between(p, b.v));
  check(SDF::SphericalPolygon(b, 0.8f, 7, -2.4f), math::angle_between(p, b.v));
  check(SDF::Star(b, 0.8f, 7, -2.4f), math::angle_between(p, b.v));
  check(SDF::Flower(b, 0.8f, 7, -2.4f), math::angle_between(p, -b.v));
}

// ============================================================================
// Inverted (complement) fill — radius past the hemisphere
// ============================================================================

/**
 * @brief Verifies the inverted fill keeps a radius > 1 shape centered on its
 *        original axis instead of jumping to the antipode.
 * @details The fill must cover the shape's own center side and exclude the
 *   folded (small) shape.
 */
inline void test_inverted_fill_stays_centered() {
  math::Basis b = equator_basis();
  const float radius = 1.5f;
  auto res = math::get_antipode(b, radius);
  const math::Basis &fb = res.first;
  const float fr = res.second;
  HS_EXPECT_NEAR(fr, 0.5f, 1e-6f);

  const math::Vector center(0, 1, 0);
  const math::Vector far_side(0, -1, 0);

  SDF::SphericalPolygon sp(fb, fr, 5, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(SDF::distance_of(sp, center).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(sp, far_side).dist > 0.0f);

  SDF::Star star(fb, fr, 5, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(SDF::distance_of(star, center).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(star, far_side).dist > 0.0f);

  SDF::PlanarPolygon pp(fb, fr, 6, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(SDF::distance_of(pp, center).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(pp, far_side).dist > 0.0f);

  // Flower fills around the antipode of its axis, so the sides swap; sample
  // off the exact poles (azimuth is degenerate there).
  SDF::Flower fl(fb, fr, 6, 0.0f, /*invert=*/true);
  math::Vector near_far(std::sin(0.1f), -std::cos(0.1f), 0.0f);
  math::Vector near_center(std::sin(0.1f), std::cos(0.1f), 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(fl, near_far).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(fl, near_center).dist > 0.0f);
}

/**
 * @brief Verifies an inverted shape scans the whole sphere: full row bounds and
 *        a full-row (unbounded) interval request per scanline.
 */
inline void test_inverted_fill_scans_full_sphere() {
  math::Basis b = equator_basis();
  auto res = math::get_antipode(b, 1.5f);

  SDF::SphericalPolygon sp(res.first, res.second, 5, 0.0f, /*invert=*/true);
  auto bounds = sp.get_vertical_bounds<144>();
  HS_EXPECT_EQ(bounds.y_min, 0);
  HS_EXPECT_EQ(bounds.y_max, 143);
  int spans = 0;
  bool handled =
      sp.get_horizontal_intervals<288, 144>(72, [&](float, float) { ++spans; });
  HS_EXPECT_TRUE(!handled);
  HS_EXPECT_EQ(spans, 0);
}

// ============================================================================
// Line
// ============================================================================

/** @brief Verifies a point on the line's arc reads raw_dist 0 and dist = -thickness. */
inline void test_line_on_arc_is_inside() {
  math::Vector a(1, 0, 0);
  math::Vector bv(0, 0, 1);
  SDF::Line ln(a, bv, /*thickness*/ 0.1f);

  math::Vector mid = ((a + bv) * 0.5f).normalized();
  auto r = SDF::distance_of(ln, mid);
  HS_EXPECT_NEAR(r.dist, -0.1f, 5e-4f);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 5e-4f);
}

/** @brief Verifies an endpoint counts as on the line (raw_dist 0, dist = -thickness). */
inline void test_line_endpoint_is_on_line() {
  math::Vector a(1, 0, 0);
  math::Vector b(0, 0, 1);
  SDF::Line ln(a, b, 0.1f);
  auto r = SDF::distance_of(ln, a);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);
}

/** @brief Verifies a point off the line's great-circle plane reads positive dist. */
inline void test_line_perpendicular_off() {
  math::Vector a(1, 0, 0);
  math::Vector b(0, 0, 1);
  SDF::Line ln(a, b, 0.05f);

  // Off the arc in +Y (perpendicular to the great-circle plane of a and b).
  math::Vector p = math::Vector(0.5f, 0.7f, 0.5f).normalized();
  auto r = SDF::distance_of(ln, p);
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

/** @brief Line distances remain stable just above the cross-product threshold. */
inline void test_line_just_above_cross_threshold() {
  const auto a = math::Vector(0.3f, -0.5f, 0.8f).normalized();
  const auto tangent = math::cross(a, math::Y_AXIS).normalized();
  for (float direction : {-1.0f, 1.0f}) {
    const auto b = (a * direction + tangent * 0.000102f).normalized();
    const auto cross = math::cross(a, b);
    HS_EXPECT_GT(math::dot(cross, cross), math::EPS_CROSS_SQ);
    const SDF::Line line(a, b, 0.01f);
    const auto midpoint = (a + b).normalized();
    const auto result = SDF::distance_of(line, midpoint);
    HS_EXPECT_LT(result.raw_dist, 0.003f);
    HS_EXPECT_LT(result.dist, 0.0f);
  }
}

/**
 * @brief Verifies a zero-length line degenerates to a point.
 * @details On the point reads dist = -thickness; a quarter-turn away reads
 *   raw_dist π/2 (positive dist).
 */
inline void test_line_degenerate_zero_length() {
  math::Vector a(1, 0, 0);
  SDF::Line ln(a, a, 0.1f);
  auto r = SDF::distance_of(ln, a);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);

  auto r2 = SDF::distance_of(ln, math::Vector(0, 1, 0));
  HS_EXPECT_NEAR(r2.raw_dist, math::PI_F * 0.5f, 5e-4f);
  HS_EXPECT_TRUE(r2.dist > 0.0f);
}

/**
 * @brief Verifies near-coincident endpoints stay point-like, not antipodal.
 * @details These two unit vectors are 1.2e-5 rad apart, so |cross|² lands under
 *   EPS_CROSS_SQ, the same band a π separation occupies.
 */
inline void test_line_near_coincident_endpoints_stay_point_like() {
  math::Vector a(0.30058673f, -0.500977874f, 0.811584115f);
  math::Vector b(0.300576329f, -0.500984192f, 0.811584175f);
  SDF::Line ln(a, b, 0.1f);

  auto bounds = ln.get_vertical_bounds<144>();
  HS_EXPECT_TRUE(bounds.y_min > 0);
  HS_EXPECT_TRUE(bounds.y_max < 143);

  int row = (bounds.y_min + bounds.y_max) / 2;
  bool handled =
      ln.get_horizontal_intervals<288, 144>(row, [](float, float) {});
  HS_EXPECT_TRUE(handled);
}

// ============================================================================
// Torus (3D volumetric)
// ============================================================================

/** @brief Verifies on the ring centerline the torus is maximally inside: dist = -minor radius. */
inline void test_torus_on_centerline_is_inside() {
  SDF::Torus t{2.0f, 0.5f};
  HS_EXPECT_NEAR(t.distance(math::Vector(2, 0, 0)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(math::Vector(0, 0, 2)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(math::Vector(-2, 0, 0)), -0.5f, 1e-5f);
}

/** @brief Verifies inner/outer rim and top-of-tube points all read dist 0 (on the surface). */
inline void test_torus_on_surface() {
  SDF::Torus t{2.0f, 0.5f};
  HS_EXPECT_NEAR(t.distance(math::Vector(2.5f, 0, 0)), 0.0f,
                 1e-5f); // outer rim, R+r
  HS_EXPECT_NEAR(t.distance(math::Vector(1.5f, 0, 0)), 0.0f,
                 1e-5f); // inner rim, R-r
  HS_EXPECT_NEAR(t.distance(math::Vector(2.0f, 0.5f, 0)), 0.0f,
                 1e-5f); // top of tube
}

/** @brief Verifies the donut-hole center is outside, at distance R - r from the tube. */
inline void test_torus_origin_is_outside_hole() {
  SDF::Torus t{2.0f, 0.5f};
  // Donut-hole center: distance = R - r = 1.5.
  HS_EXPECT_NEAR(t.distance(math::Vector(0, 0, 0)), 1.5f, 1e-5f);
}

/** @brief Verifies the surface normal on the outer rim points radially outward (+X here). */
inline void test_torus_normal_points_outward_on_outer_rim() {
  SDF::Torus t{2.0f, 0.5f};
  math::Vector n = t.normal(math::Vector(2.5f, 0, 0));
  HS_EXPECT_VEC(n, math::Vector(1, 0, 0), 1e-4f);
}

/** @brief Verifies the surface normal at the top of the tube points +Y. */
inline void test_torus_normal_points_outward_on_top() {
  SDF::Torus t{2.0f, 0.5f};
  math::Vector n = t.normal(math::Vector(2.0f, 0.5f, 0));
  HS_EXPECT_VEC(n, math::Vector(0, 1, 0), 1e-4f);
}
