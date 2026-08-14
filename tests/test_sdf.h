/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/sdf.h.
 *
 * Coverage:
 *   - clamp_phi utility
 *   - Spherical SDF primitives (Ring, PlanarPolygon, SphericalPolygon, Star, Line)
 *   - 3D Torus
 *   - WarpedVolume with Warp::Twist
 *   - CSG operators (Union, SmoothUnion, Subtract, Intersection)
 *   - AngularRepeat
 *
 * Tests focus on the distance() interface — the rendering pipeline
 * (get_vertical_bounds<H> / get_horizontal_intervals<W,H>) is exercised
 * indirectly by feeding known points and verifying signed distance.
 */
#pragma once

#include "core/render/sdf.h"
#include "core/render/sdf_volume.h"
#include "core/mesh/mesh_classes.h"
#include "core/render/scan.h"
#include "core/math/geometry.h"
#include "tests/vec_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cfloat>
#include <cmath>
#include <cstring>
#include <iterator>
#include <type_traits>
#include <utility>
#include <vector>

namespace hs_test {
namespace sdf_tests {

using hs_test::approx_vec;

/**
 * @brief Builds the canonical equator-facing basis: v = +Y, u = +X, w = +Z.
 * @return A Basis oriented so its pole points along +Y.
 */
inline Basis equator_basis() {
  return Basis{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
}

// ============================================================================
// clamp_phi
// ============================================================================

/** @brief Verifies phi already in [0, π] passes through clamp_phi unchanged. */
inline void test_clamp_phi_in_range() {
  HS_EXPECT_NEAR(SDF::clamp_phi(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(0.5f), 0.5f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(PI_F), PI_F, 1e-6f);
}

/** @brief Verifies negative phi reflects across the north pole (|phi|). */
inline void test_clamp_phi_negative_reflects() {
  HS_EXPECT_NEAR(SDF::clamp_phi(-0.3f), 0.3f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(-1.2f), 1.2f, 1e-6f);
}

/** @brief Verifies phi above π reflects across the south pole (2π - phi). */
inline void test_clamp_phi_above_pi_reflects() {
  HS_EXPECT_NEAR(SDF::clamp_phi(PI_F + 0.2f), PI_F - 0.2f, 1e-5f);
  HS_EXPECT_NEAR(SDF::clamp_phi(2.0f * PI_F), 0.0f, 1e-5f);
}

/**
 * @brief Verifies inputs outside [-π, 2π] still fold into [0, π] (full-range
 *        acosf(cosf(x)) equivalence), not the out-of-range values the old
 *        single-reflection body returned.
 */
inline void test_clamp_phi_full_range() {
  // 2π + 0.2 folds to 0.2.
  HS_EXPECT_NEAR(SDF::clamp_phi(2.0f * PI_F + 0.2f), 0.2f, 1e-5f);
  // acosf(cosf(3π)) = π.
  HS_EXPECT_NEAR(SDF::clamp_phi(3.0f * PI_F), PI_F, 1e-5f);
  // -1.5π folds to 0.5π.
  HS_EXPECT_NEAR(SDF::clamp_phi(-1.5f * PI_F), 0.5f * PI_F, 1e-5f);
}

/**
 * @brief Verifies clamp_phi_band reports the exact colatitude extent of the
 *        circle it bounds, against a brute-force sweep of that circle.
 * @details The extremes sit at psi = 0 and psi = π, both sampled exactly, so
 *   the sweep is an exact oracle rather than an approximation.
 */
inline void test_clamp_phi_band_matches_circle_extent() {
  constexpr int SAMPLES = 512;
  auto circle_extent = [](float c, float t) {
    float lo = PI_F, hi = 0.0f;
    for (int i = 0; i < SAMPLES; ++i) {
      float psi = TWO_PI_F * i / SAMPLES;
      float y =
          std::cos(t) * std::cos(c) - std::sin(t) * std::sin(c) * std::cos(psi);
      float phi = std::acos(std::max(-1.0f, std::min(1.0f, y)));
      lo = std::min(lo, phi);
      hi = std::max(hi, phi);
    }
    return std::pair<float, float>(lo, hi);
  };

  const float centers[] = {0.0f, 0.2f, 0.9f, 1.5f, 2.4f, 3.0f, PI_F};
  // Past π and negative: a complement-radius ring and a mirrored one still
  // trace a real circle, so the band must order its folded endpoints.
  const float radii[] = {0.0f, 0.1f, 0.5f, 1.0f, 2.0f, 2.5f, PI_F, 4.0f, -0.7f};
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
      HS_EXPECT_LE(band.phi_max, PI_F);
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
  HS_EXPECT_NEAR(wrapped.phi_max, TWO_PI_F - 5.5f, 1e-5f);
}

/** @brief Verifies the reciprocal sector fold matches the general wrap path. */
inline void test_centered_sector_angle_matches_wrap() {
  const float angles[] = {-13.7f, -2.3f, -0.7f, 0.21f, 1.8f, 9.4f};
  for (int sides : {3, 5, 8, 17}) {
    float sector = 2.0f * PI_F / sides;
    float reciprocal_sector = static_cast<float>(sides) / (2.0f * PI_F);
    for (float angle : angles) {
      const float folded =
          SDF::centered_sector_angle(angle, sector, reciprocal_sector);
      float expected = wrap(angle + sector * 0.5f, sector) - sector * 0.5f;
      HS_EXPECT_NEAR(folded, expected, 1e-5f);
      // Independent of both folds: the result lands in the centered sector and
      // differs from the input by a whole number of sectors.
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

/** @brief Verifies a point on the ring centerline reads raw_dist 0 and dist = -thickness. */
inline void test_ring_on_centerline() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  auto r = SDF::distance_of(ring, Vector(1, 0, 0));
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
}

/** @brief Verifies a point within the ring band reads negative dist. */
inline void test_ring_inside_band() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  float off = 0.05f;
  Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = SDF::distance_of(ring, p);
  HS_EXPECT_TRUE(r.dist < 0.0f);
  HS_EXPECT_TRUE(r.raw_dist <= 0.1f + 1e-3f);
}

/** @brief Verifies a point far outside the band reads the cull sentinel rather than a real dist. */
inline void test_ring_outside_band_returns_sentinel() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  auto r = SDF::distance_of(ring, Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist > 50.0f);
}

/** @brief Verifies just past the band edge still trips the sentinel (band edge is exclusive). */
inline void test_ring_just_outside_band() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.05f);

  // off (0.07) > thickness (0.05): point is 0.02 past the band edge.
  float off = 0.07f;
  Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = SDF::distance_of(ring, p);
  HS_EXPECT_TRUE(r.dist > 50.0f);
}

/**
 * @brief Verifies a small-radius ring thick enough to break the linearized
 *        centerline distance still reports the exact geodesic distance,
 *        symmetric about the centerline.
 */
inline void test_ring_small_radius_distance_symmetric() {
  const float RADIUS = 0.04f;
  const float THICKNESS = 0.05f;
  Basis b = equator_basis();
  SDF::Ring ring(b, RADIUS, THICKNESS);
  float target = RADIUS * (PI_F / 2.0f);

  // Ring axis is +Y, so a point at polar angle a off the axis is
  // (sin a, cos a, 0).
  auto at_polar = [](float a) {
    return Vector(std::sin(a), std::cos(a), 0.0f);
  };

  for (float off : {-0.9f, -0.4f, 0.4f, 0.9f}) {
    float delta = off * THICKNESS;
    auto r = SDF::distance_of(ring, at_polar(target + delta));
    HS_EXPECT_NEAR(r.raw_dist, std::abs(delta), 1e-4f);
    HS_EXPECT_NEAR(r.dist, std::abs(delta) - THICKNESS, 1e-4f);
  }

  // The AA ramp reads the same either side of the centerline.
  float d_in = dot(at_polar(target - 0.9f * THICKNESS), b.v);
  float d_out = dot(at_polar(target + 0.9f * THICKNESS), b.v);
  HS_EXPECT_NEAR(ring.stroke_alpha(d_in), ring.stroke_alpha(d_out), 1e-4f);
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
  Basis b = equator_basis(); // v=+Y, u=+X, w=+Z; radius=1 → target_angle = π/2
  const float shift = 0.2f;
  const float thickness = 0.05f;

  // Azimuth 0 (along +X) on the shifted centerline: polar angle from +Y is
  // π/2 + shift.
  Vector p(std::sin(PI_F / 2 + shift), std::cos(PI_F / 2 + shift), 0.0f);

  SDF::DistortedRing shifted(
      b, 1.0f, thickness, [shift](float) { return shift; },
      /*max_distortion=*/shift, /*phase=*/0.0f);
  auto rs = SDF::distance_of(shifted, p);
  HS_EXPECT_TRUE(rs.dist < 50.0f);
  HS_EXPECT_NEAR(rs.raw_dist, 0.0f, 1e-2f);
  HS_EXPECT_NEAR(rs.dist, -thickness, 1e-2f);
  HS_EXPECT_NEAR(rs.t, 0.0f, 1e-2f);

  // Same point, no shift: the centerline stays at π/2, so it now sits `shift`
  // radians off (raw_dist ≈ shift) — the shift moved the centerline.
  SDF::DistortedRing plain(
      b, 1.0f, thickness, [](float) { return 0.0f; },
      /*max_distortion=*/shift, /*phase=*/0.0f);
  auto rp = SDF::distance_of(plain, p);
  HS_EXPECT_NEAR(rp.raw_dist, shift, 1e-2f);
}

/**
 * @brief Verifies a sinusoidal shift_fn moves the centerline by a per-azimuth
 *        amount, so the t parameter feeding shift_fn is wired correctly.
 */
inline void test_distorted_ring_sin_shift_varies_by_azimuth() {
  Basis b = equator_basis();
  const float amp = 0.2f;
  SDF::DistortedRing ring(
      b, 1.0f, 0.05f, [amp](float t) { return amp * std::sin(2 * PI_F * t); },
      amp, 0.0f);

  // Azimuth π/2 (along +Z) → t = 0.25 → shift = amp; centerline polar angle is
  // π/2 + amp.
  Vector on(0.0f, std::cos(PI_F / 2 + amp), std::sin(PI_F / 2 + amp));
  auto r_on = SDF::distance_of(ring, on);
  HS_EXPECT_NEAR(r_on.t, 0.25f, 1e-2f);
  HS_EXPECT_NEAR(r_on.raw_dist, 0.0f, 1e-2f);

  // Same azimuth on the unshifted equator (+Z): centerline moved by amp here.
  auto r_off = SDF::distance_of(ring, Vector(0, 0, 1));
  HS_EXPECT_NEAR(r_off.t, 0.25f, 1e-2f);
  HS_EXPECT_NEAR(r_off.raw_dist, amp, 1e-2f);
}

/**
 * @brief Verifies flat mode preserves the exact zero-knot polar distance.
 */
inline void test_distorted_ring_flat_matches_zero_knots() {
  constexpr int LUT_N = 16;
  float knots[LUT_N + 1] = {};

  auto check = [&](const Basis &basis, float radius) {
    constexpr float thickness = 0.08f;
    SDF::FlatDistortedRing flat(basis, radius, thickness);
    SDF::KnotPrefilter pf;
    SDF::DistortedRing polyline(basis, radius, thickness, knots, LUT_N, 0.0f,
                                pf);
    const float target = radius * (PI_F / 2.0f);
    const float azimuths[] = {0.0f, 1e-5f, PI_F / 2.0f, PI_F,
                              2.0f * PI_F - 1e-5f};
    const float offsets[] = {-0.04f, 0.0f, 0.04f};
    for (float azimuth : azimuths) {
      for (float offset : offsets) {
        float polar = hs::clamp(target + offset, 0.0f, PI_F);
        Vector p =
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

  check(make_basis(Quaternion(), Y_AXIS), 0.01f);
  check(make_basis(Quaternion(), X_AXIS), 1.0f);
  check(make_basis(Quaternion(), Vector(0.3f, 0.8f, -0.5f).normalized()),
        1.99f);
}

/**
 * @brief Verifies the knot overload's raw_dist matches a brute-force geodesic
 *        minimum over the densely sampled polyline, within the stroke reach.
 * @details Probes pixels around crests, steep flanks, and the wrap seam of a
 *   high-harmonic knot polyline on a tilted ring. Accuracy is only contracted
 *   for true distances below thickness (the outward search stops at that
 *   reach); every probe is placed inside it.
 */
template <int LUT_N> inline void expect_polyline_distance_matches_bruteforce() {
  Basis b = make_basis(Quaternion(), Vector(0.3f, 1.0f, 0.2f));
  const float amp = 0.2f;
  const int harmonic = 5;
  const float radius = 0.5f; // target_angle = π/4: curved chart, tilted axis
  const float thickness = 0.06f;
  float knots[LUT_N + 1];
  for (int k = 0; k <= LUT_N; ++k)
    knots[k] = amp * std::sin(2.0f * PI_F * harmonic * (k % LUT_N) / LUT_N);
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(b, radius, thickness, knots, LUT_N, 0.0f, pf);

  const float target = radius * (PI_F / 2.0f);
  auto on_sphere = [&](float t, float dv) {
    float theta = target + amp * std::sin(2.0f * PI_F * harmonic * t) + dv;
    float a = 2.0f * PI_F * t;
    return (b.v * std::cos(theta)) +
           ((b.u * std::cos(a)) + (b.w * std::sin(a))) * std::sin(theta);
  };
  auto brute = [&](const Vector &p) {
    constexpr int SAMPLES = LUT_N * 64;
    float best = 100.0f;
    for (int s = 0; s < SAMPLES; ++s) {
      int k = s / 64;
      float f = (s % 64) / 64.0f;
      float theta = target + knots[k] + f * (knots[k + 1] - knots[k]);
      float a = 2.0f * PI_F * (k + f) / LUT_N;
      Vector q = (b.v * std::cos(theta)) +
                 ((b.u * std::cos(a)) + (b.w * std::sin(a))) * std::sin(theta);
      best = std::min(best, std::acos(hs::clamp(dot(p, q), -1.0f, 1.0f)));
    }
    return best;
  };

  // t = 0.05: crest (slope zero, curvature max); t = 0.1: steep flank;
  // t = 0.998: wrap seam. Offsets stay within the thickness reach.
  const float probes[][2] = {{0.05f, 0.04f},   {0.05f, -0.05f}, {0.1f, 0.05f},
                             {0.1f, -0.04f},   {0.998f, 0.04f}, {0.25f, 0.05f},
                             {0.375f, -0.05f}, {0.6f, 0.03f}};
  for (const auto &pr : probes) {
    Vector p = on_sphere(pr[0], pr[1]);
    auto r = SDF::distance_of(ring, p);
    float expected = brute(p);
    HS_EXPECT_TRUE(expected < thickness);
    HS_EXPECT_NEAR(r.raw_dist, expected, 3e-3f);
  }
}

/**
 * @brief Runs the brute-force comparison at a chunk-aligned and a
 *        chunk-straddling knot count (97 is coprime to the prefilter chunks,
 *        so segments cross chunk boundaries).
 */
inline void test_distorted_ring_polyline_distance_matches_bruteforce() {
  expect_polyline_distance_matches_bruteforce<96>();
  expect_polyline_distance_matches_bruteforce<97>();
}

/** @brief Verifies knot extrema tighten an asymmetric loose caller band. */
inline void test_distorted_ring_knot_extrema_tighten_band() {
  constexpr int LUT_N = 8;
  constexpr float RADIUS = 0.8f;
  constexpr float THICKNESS = 0.05f;
  constexpr float TARGET = RADIUS * (PI_F / 2.0f);
  float knots[LUT_N + 1] = {0.18f, 0.12f, 0.04f, -0.01f, -0.03f,
                            0.02f, 0.09f, 0.16f, 0.18f};
  Basis basis = equator_basis();
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(basis, RADIUS, THICKNESS, knots, LUT_N, 0.0f, pf);

  HS_EXPECT_NEAR(ring.max_distortion, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(ring.max_thickness, 0.23f, 1e-6f);
  for (int k : {0, 4, 8}) {
    float t = static_cast<float>(k % LUT_N) / LUT_N;
    float azimuth = 2.0f * PI_F * t;
    float polar = TARGET + knots[k];
    Vector p(cosf(azimuth) * sinf(polar), cosf(polar),
             sinf(azimuth) * sinf(polar));
    HS_EXPECT_NEAR(SDF::distance_of(ring, p).raw_dist, 0.0f, 2e-4f);
  }

  float below = TARGET - 0.03f - THICKNESS - 1e-3f;
  float above = TARGET + 0.18f + THICKNESS + 1e-3f;
  HS_EXPECT_GT(
      SDF::distance_of(ring, Vector(sinf(below), cosf(below), 0.0f)).dist,
      50.0f);
  HS_EXPECT_GT(
      SDF::distance_of(ring, Vector(sinf(above), cosf(above), 0.0f)).dist,
      50.0f);
}

/**
 * @brief Verifies a knot ring reports the far sentinel past its stroke reach.
 * @details The outward search stops at the reach, so a pixel the prefilter lets
 *   through but no segment comes within thickness of has no distance to report.
 *   Reporting the reach itself puts dist on exactly 0 — the surface — and a
 *   Subtract then carves its minuend across the whole bounding annulus.
 */
inline void test_distorted_ring_past_reach_reports_far_sentinel() {
  constexpr int LUT_N = 32;
  constexpr float RADIUS = 0.8f;
  constexpr float THICKNESS = 0.05f;
  constexpr float TARGET = RADIUS * (PI_F / 2.0f);
  constexpr float SPIKE = 0.5f;
  // One tall spike at azimuth 0, every other cell on the centerline. A pixel in
  // the next chunk shares the spike's prefilter window, so the segment search
  // runs and terminates without a point inside the reach.
  float knots[LUT_N + 1] = {};
  knots[0] = SPIKE;
  knots[LUT_N] = SPIKE;
  Basis basis = equator_basis();
  SDF::KnotPrefilter pf;
  SDF::DistortedRing ring(basis, RADIUS, THICKNESS, knots, LUT_N, 0.0f, pf);

  const float azimuth = 2.0f * PI_F * 1.5f / LUT_N;
  const float polar = TARGET + SPIKE - 0.05f;
  Vector p(cosf(azimuth) * sinf(polar), cosf(polar),
           sinf(azimuth) * sinf(polar));
  HS_EXPECT_GT(SDF::distance_of(ring, p).dist, 50.0f);

  Basis poly_basis = make_basis(Quaternion(), p);
  SDF::PlanarPolygon poly(poly_basis, /*radius=*/0.3f, /*sides=*/6, 0.0f);
  SDF::Subtract<SDF::PlanarPolygon, SDF::DistortedRing> carved(poly, ring);
  const float solid = SDF::distance_of(poly, p).dist;
  HS_EXPECT_LT(solid, -0.1f);
  HS_EXPECT_NEAR(SDF::distance_of(carved, p).dist, solid, 1e-6f);
}

// ============================================================================
// PlanarPolygon  (Basis at top of sphere; distance to nearest edge)
// ============================================================================

/** @brief Verifies the polygon center is inside, with dist equal to the negated apothem. */
inline void test_polygon_at_center_inside() {
  Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, /*circumradius*/ 0.5f, /*sides*/ 6,
                          /*phase*/ 0.0f);

  auto r = SDF::distance_of(poly, Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
  float apothem = 0.5f * std::cos(PI_F / 6.0f);
  HS_EXPECT_NEAR(r.dist, -apothem, 1e-3f);
}

/** @brief Verifies the antipode of the polygon center is outside (positive dist). */
inline void test_polygon_far_point_outside() {
  Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, 0.3f, 6, 0.0f);

  auto r = SDF::distance_of(poly, Vector(0, -1, 0));
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

// ============================================================================
// SphericalPolygon — great-circle edges
// ============================================================================

/** @brief Verifies the spherical-polygon center is strictly inside. */
inline void test_spherical_polygon_center_inside() {
  Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, /*radius*/ 0.5f, /*sides*/ 5, /*phase*/ 0.0f);
  auto r = SDF::distance_of(sp, Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies the antipode of the spherical-polygon center is outside. */
inline void test_spherical_polygon_far_outside() {
  Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, 0.3f, 6, 0.0f);
  auto r = SDF::distance_of(sp, Vector(0, -1, 0));
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
  Basis b = equator_basis();
  const int sides = 5;
  const float radius = 0.5f;
  SDF::SphericalPolygon sp(b, radius, sides, 0.0f);

  const float R = radius * (PI_F / 2.0f);
  const float inradius = std::atan(std::tan(R) * std::cos(PI_F / sides));

  auto center = SDF::distance_of(sp, Vector(0, 1, 0));
  HS_EXPECT_NEAR(center.dist, -inradius, 1e-2f);
  HS_EXPECT_NEAR(center.raw_dist, 0.0f, 1e-3f);

  // Edge midpoint: polar angle = inradius along +u (the sector bisector).
  Vector edge_mid(std::sin(inradius), std::cos(inradius), 0.0f);
  auto em = SDF::distance_of(sp, edge_mid);
  HS_EXPECT_NEAR(em.dist, 0.0f, 1e-2f);
  HS_EXPECT_NEAR(em.raw_dist, inradius, 1e-2f);
}

/**
 * @brief Bounds sine-domain distance error across the device-width AA band.
 * @details Where the edge dot wins, the two paths share it verbatim and differ
 *   only by sin(x) - x (~1.7e-6 at one pixel_width). Where the
 *   circumscribed-disc clamp wins, distance() spends fast_acos on polar while
 *   sine_distance forms sin(polar - circumradius) from exact trig, so the gap
 *   widens to fast_acos' ~1.3e-4 rad peak -- under 1% of a pixel_width, and the
 *   sine path is the more accurate of the two there.
 */
inline void test_spherical_polygon_sine_distance_aa_error() {
  constexpr int W = 288;
  constexpr int H = 144;
  constexpr float PIXEL_WIDTH = 2.0f * PI_F / W;
  Basis basis = make_basis(
      make_rotation(Vector(0.3f, -0.8f, 0.5f).normalized(), 0.71f), Y_AXIS);
  struct Case {
    float radius;
    int sides;
    float phase;
  };
  const Case cases[] = {{0.22f, 3, -2.7f},
                        {0.72f, 5, 0.31f},
                        {0.98f, 12, 5.4f},
                        {1.42f, 7, -0.9f}};

  float max_error = 0.0f;
  int edge_samples = 0;
  for (const Case &c : cases) {
    auto folded = get_antipode(basis, c.radius);
    SDF::SphericalPolygon shape(folded.first, folded.second, c.sides, c.phase,
                                c.radius > 1.0f);
    for (int y = 0; y < H; ++y) {
      float polar = PI_F * (static_cast<float>(y) + 0.5f) / H;
      float sin_p = sinf(polar);
      float cos_p = cosf(polar);
      for (int x = 0; x < W; ++x) {
        float azimuth = 2.0f * PI_F * (static_cast<float>(x) + 0.5f) / W;
        Vector p(sin_p * cosf(azimuth), cos_p, sin_p * sinf(azimuth));
        SDF::DistanceResult exact;
        shape.distance<false>(p, exact);
        float sine = shape.sine_distance(p);
        HS_EXPECT_TRUE(exact.dist == 0.0f || sine == 0.0f ||
                       std::signbit(exact.dist) == std::signbit(sine));
        if (std::abs(exact.dist) <= PIXEL_WIDTH) {
          max_error = std::max(max_error, std::abs(exact.dist - sine));
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

/**
 * @brief Verifies SphericalPolygon satisfies SDFShape and composes under CSG.
 * @details Every combinator static_asserts SDFShape on its children, so a leaf
 *   that fails it is barred from all of them and its sdf_max_spans
 *   specialization never applies. Two disjoint polygons must emit their two arcs
 *   through the span-bounded scanline path.
 */
inline void test_spherical_polygon_composes_under_csg() {
  static_assert(SDF::SDFShape<SDF::SphericalPolygon>,
                "SphericalPolygon must satisfy the CSG child contract");
  using U = SDF::Union<SDF::SphericalPolygon, SDF::SphericalPolygon>;
  static_assert(SDF::sdf_max_spans<SDF::SphericalPolygon>::value == 1);
  static_assert(SDF::sdf_max_spans<U>::value == 2);
  static_assert(SDF::sdf_max_spans<SDF::Intersection<U, U>>::value == 8);

  Basis b = equator_basis();
  SDF::SphericalPolygon inner(b, 0.3f, 5, 0.0f);
  SDF::SphericalPolygon outer(b, 0.7f, 5, 0.0f);

  U coaxial(inner, outer);
  // Sector bisector (+u) at polar 0.7: past inner's inradius, short of outer's.
  Vector p(std::sin(0.7f), std::cos(0.7f), 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(inner, p).dist > 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(outer, p).dist < 0.0f);
  HS_EXPECT_NEAR(SDF::distance_of(coaxial, p).dist,
                 SDF::distance_of(outer, p).dist, 1e-6f);

  // Poles on +X and +Z: both cross the equatorial row, a quarter turn apart.
  constexpr int W = 256, H = 128;
  const Basis bx{Vector(0, 1, 0), Vector(1, 0, 0), Vector(0, 0, 1)};
  const Basis bz{Vector(0, 1, 0), Vector(0, 0, 1), Vector(-1, 0, 0)};
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
  Basis b = equator_basis();
  SDF::Star star(b, /*radius*/ 0.6f, /*sides*/ 5, /*phase*/ 0.0f);
  auto r = SDF::distance_of(star, Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies the antipode of the star center is outside. */
inline void test_star_far_outside() {
  Basis b = equator_basis();
  SDF::Star star(b, 0.4f, 5, 0.0f);
  auto r = SDF::distance_of(star, Vector(0, -1, 0));
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
  Basis b = equator_basis();
  const float radius = 0.6f;
  SDF::Star star(b, radius, /*sides=*/5, 0.0f);

  const float outer = radius * (PI_F / 2.0f);
  Vector tip(std::sin(outer), std::cos(outer), 0.0f);
  auto r = SDF::distance_of(star, tip);
  HS_EXPECT_NEAR(r.dist, 0.0f, 1e-2f);
  HS_EXPECT_NEAR(r.raw_dist, outer, 1e-2f);
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
  Basis b = equator_basis();
  const float radius = 0.6f;
  SDF::Flower flower(b, radius, /*sides=*/6, 0.0f);

  const float outer = radius * (PI_F / 2.0f);
  const float s = 0.1f;
  // From the antipode (-Y), step s toward +u (+X): interior of a petal.
  Vector p(std::sin(s), -std::cos(s), 0.0f);
  auto r = SDF::distance_of(flower, p);
  HS_EXPECT_NEAR(r.dist, s - outer, 1e-2f);
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
  Basis b = equator_basis();
  const float radius = 0.6f;
  SDF::Flower flower(b, radius, /*sides=*/6, 0.0f);

  const float outer = radius * (PI_F / 2.0f);
  // From the antipode (-Y), step `outer` toward +u (+X).
  Vector tip(std::sin(outer), -std::cos(outer), 0.0f);
  auto r = SDF::distance_of(flower, tip);
  HS_EXPECT_NEAR(r.dist, 0.0f, 1e-2f);
  HS_EXPECT_NEAR(r.raw_dist, outer, 1e-2f);
}

/** @brief Verifies the solid-shape unit-vector and no-UV distance paths. */
inline void test_solid_shape_unit_angle_and_no_uv_paths() {
  Basis b = make_basis(Quaternion(), Vector(0.3f, 0.8f, -0.5f).normalized());
  float polar = 0.73f;
  float azimuth = -1.17f;
  Vector p = b.v * cosf(polar) +
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

  check(SDF::PlanarPolygon(b, 0.8f, 7, -2.4f), angle_between(p, b.v));
  check(SDF::SphericalPolygon(b, 0.8f, 7, -2.4f), angle_between(p, b.v));
  check(SDF::Star(b, 0.8f, 7, -2.4f), angle_between(p, b.v));
  check(SDF::Flower(b, 0.8f, 7, -2.4f), angle_between(p, -b.v));
}

// ============================================================================
// Inverted (complement) fill — radius past the hemisphere
// ============================================================================

/**
 * @brief Verifies the inverted fill keeps a radius > 1 shape centered on its
 *        original axis instead of jumping to the antipode.
 * @details Builds each solid shape the way the Scan wrappers do for radius 1.5:
 *   antipode-folded basis, folded radius 0.5, invert = true. The fill must
 *   cover the shape's own center side and exclude the folded (small) shape;
 *   Flower's fill is centered on the antipode of its axis, so its sides swap.
 */
inline void test_inverted_fill_stays_centered() {
  Basis b = equator_basis();
  const float radius = 1.5f;
  auto res = get_antipode(b, radius);
  const Basis &fb = res.first;
  const float fr = res.second;
  HS_EXPECT_NEAR(fr, 0.5f, 1e-6f);

  const Vector center(0, 1, 0);
  const Vector far_side(0, -1, 0);

  SDF::SphericalPolygon sp(fb, fr, 5, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(SDF::distance_of(sp, center).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(sp, far_side).dist > 0.0f);

  SDF::Star star(fb, fr, 5, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(SDF::distance_of(star, center).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(star, far_side).dist > 0.0f);

  SDF::PlanarPolygon pp(fb, fr * (PI_F / 2.0f), 6, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(SDF::distance_of(pp, center).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(pp, far_side).dist > 0.0f);

  // Flower fills around the antipode of its axis, so the sides swap; sample
  // off the exact poles (azimuth is degenerate there).
  SDF::Flower fl(fb, fr, 6, 0.0f, /*invert=*/true);
  Vector near_far(std::sin(0.1f), -std::cos(0.1f), 0.0f);
  Vector near_center(std::sin(0.1f), std::cos(0.1f), 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(fl, near_far).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(fl, near_center).dist > 0.0f);
}

/**
 * @brief Verifies an inverted shape scans the whole sphere: full row bounds and
 *        a full-row (unbounded) interval request per scanline.
 */
inline void test_inverted_fill_scans_full_sphere() {
  Basis b = equator_basis();
  auto res = get_antipode(b, 1.5f);

  SDF::SphericalPolygon sp(res.first, res.second, 5, 0.0f, /*invert=*/true);
  auto bounds = sp.get_vertical_bounds<144>();
  HS_EXPECT_EQ(bounds.y_min, 0);
  HS_EXPECT_EQ(bounds.y_max, 143);
  bool handled = sp.get_horizontal_intervals<288, 144>(72, [](float, float) {});
  HS_EXPECT_TRUE(!handled);
}

// ============================================================================
// Line
// ============================================================================

/** @brief Verifies a point on the line's arc reads raw_dist 0 and dist = -thickness. */
inline void test_line_on_arc_is_inside() {
  Vector a(1, 0, 0);
  Vector bv(0, 0, 1);
  SDF::Line ln(a, bv, /*thickness*/ 0.1f);

  Vector mid = ((a + bv) * 0.5f).normalized();
  auto r = SDF::distance_of(ln, mid);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-2f);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-2f);
}

/** @brief Verifies an endpoint counts as on the line (raw_dist 0, dist = -thickness). */
inline void test_line_endpoint_is_on_line() {
  Vector a(1, 0, 0);
  Vector b(0, 0, 1);
  SDF::Line ln(a, b, 0.1f);
  auto r = SDF::distance_of(ln, a);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);
}

/** @brief Verifies a point off the line's great-circle plane reads positive dist. */
inline void test_line_perpendicular_off() {
  Vector a(1, 0, 0);
  Vector b(0, 0, 1);
  SDF::Line ln(a, b, 0.05f);

  // Off the arc in +Y (perpendicular to the great-circle plane of a and b).
  Vector p = Vector(0.5f, 0.7f, 0.5f).normalized();
  auto r = SDF::distance_of(ln, p);
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

/**
 * @brief Verifies a zero-length line degenerates to a point.
 * @details On the point reads dist = -thickness; a quarter-turn away reads
 *   raw_dist π/2 (positive dist).
 */
inline void test_line_degenerate_zero_length() {
  Vector a(1, 0, 0);
  SDF::Line ln(a, a, 0.1f);
  auto r = SDF::distance_of(ln, a);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);

  auto r2 = SDF::distance_of(ln, Vector(0, 1, 0));
  HS_EXPECT_NEAR(r2.raw_dist, PI_F * 0.5f, 1e-2f);
  HS_EXPECT_TRUE(r2.dist > 0.0f);
}

/**
 * @brief Verifies near-coincident endpoints stay point-like, not antipodal.
 * @details These two unit vectors are 1.2e-5 rad apart, so |cross|² lands under
 *   EPS_CROSS_SQ — the same band a π separation occupies — while
 *   angle_between's acos reports 4.9e-4 (quantisation, not geometry). Read as
 *   antipodal, the arc becomes a substituted great circle: every canvas row in
 *   bounds and no horizontal cull, i.e. a full-width scan per row for one dot.
 */
inline void test_line_near_coincident_endpoints_stay_point_like() {
  Vector a(0.30058673f, -0.500977874f, 0.811584115f);
  Vector b(0.300576329f, -0.500984192f, 0.811584175f);
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
  HS_EXPECT_NEAR(t.distance(Vector(2, 0, 0)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(Vector(0, 0, 2)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(Vector(-2, 0, 0)), -0.5f, 1e-5f);
}

/** @brief Verifies inner/outer rim and top-of-tube points all read dist 0 (on the surface). */
inline void test_torus_on_surface() {
  SDF::Torus t{2.0f, 0.5f};
  HS_EXPECT_NEAR(t.distance(Vector(2.5f, 0, 0)), 0.0f, 1e-5f); // outer rim, R+r
  HS_EXPECT_NEAR(t.distance(Vector(1.5f, 0, 0)), 0.0f, 1e-5f); // inner rim, R-r
  HS_EXPECT_NEAR(t.distance(Vector(2.0f, 0.5f, 0)), 0.0f, 1e-5f); // top of tube
}

/** @brief Verifies the donut-hole center is outside, at distance R - r from the tube. */
inline void test_torus_origin_is_outside_hole() {
  SDF::Torus t{2.0f, 0.5f};
  // Donut-hole center: distance = R - r = 1.5.
  HS_EXPECT_NEAR(t.distance(Vector(0, 0, 0)), 1.5f, 1e-5f);
}

/** @brief Verifies the surface normal on the outer rim points radially outward (+X here). */
inline void test_torus_normal_points_outward_on_outer_rim() {
  SDF::Torus t{2.0f, 0.5f};
  Vector n = t.normal(Vector(2.5f, 0, 0));
  HS_EXPECT_VEC(n, Vector(1, 0, 0), 1e-4f);
}

/** @brief Verifies the surface normal at the top of the tube points +Y. */
inline void test_torus_normal_points_outward_on_top() {
  SDF::Torus t{2.0f, 0.5f};
  Vector n = t.normal(Vector(2.0f, 0.5f, 0));
  HS_EXPECT_VEC(n, Vector(0, 1, 0), 1e-4f);
}

// ============================================================================
// WarpedVolume + Warp::Twist (domain warp, Lipschitz bound, normal correction)
// ============================================================================

/** @brief Verifies Twist::apply displaces Y by amplitude·sin(twist·θ), θ=atan2(z,x). */
inline void test_twist_apply_displaces_y() {
  SDF::Warp::Twist tw{/*twist=*/1, /*amplitude=*/0.3f, /*R=*/1.0f};

  // θ = atan2(0, 1) = 0 → no displacement.
  Vector a(1.0f, 0.5f, 0.0f);
  Vector ra = tw.apply(a, tw.make_ctx(a));
  HS_EXPECT_VEC(ra, Vector(1.0f, 0.5f, 0.0f), 1e-3f);

  // θ = π/2 → sin(twist·π/2) = 1 → Y drops by amplitude.
  Vector b(0.0f, 0.5f, 1.0f);
  Vector rb = tw.apply(b, tw.make_ctx(b));
  HS_EXPECT_VEC(rb, Vector(0.0f, 0.5f - 0.3f, 1.0f), 1e-2f);
}

/** @brief Verifies Twist::lipschitz is 1 for twist 0 and matches the closed form otherwise. */
inline void test_twist_lipschitz_identity_and_closed_form() {
  SDF::Warp::Twist flat{0, 0.5f, 1.0f};
  HS_EXPECT_NEAR(
      flat.lipschitz(Vector(2, 0, 0), flat.make_ctx(Vector(2, 0, 0))), 1.0f,
      1e-6f);

  // twist=2, amplitude=0.5 at s=2: γ = 0.5, bound = γ/2 + √(1 + γ²/4).
  SDF::Warp::Twist tw{2, 0.5f, 1.0f};
  Vector p(2, 0, 0);
  float s = tw.make_ctx(p);
  HS_EXPECT_NEAR(s, 2.0f, 1e-6f);
  float gamma = 0.5f;
  float expected = 0.5f * gamma + std::sqrt(1.0f + 0.25f * gamma * gamma);
  HS_EXPECT_NEAR(tw.lipschitz(p, s), expected, 1e-5f);
}

/** @brief Verifies Twist::bounding_inflation returns the displacement amplitude. */
inline void test_twist_bounding_inflation() {
  SDF::Warp::Twist tw{3, 0.42f, 1.0f};
  HS_EXPECT_NEAR(tw.bounding_inflation(), 0.42f, 1e-6f);
}

/**
 * @brief Verifies WarpedVolume::distance never over-estimates the warped
 *        distance (sphere-trace safety): on every sampled point the returned
 *        march distance is <= the raw warped distance, on both the bounding
 *        fast-path and the Lipschitz-corrected path.
 */
inline void test_warped_volume_distance_is_sphere_trace_safe() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{
      SDF::Torus{1.0f, 0.3f}, SDF::Warp::Twist{3, 0.2f, 1.0f}};

  for (float x = -2.0f; x <= 2.0f; x += 0.5f)
    for (float y = -1.0f; y <= 1.0f; y += 0.5f)
      for (float z = -2.0f; z <= 2.0f; z += 0.5f) {
        Vector p(x, y, z);
        float d = wv.distance(p);
        float raw = wv.raw_distance(p);
        HS_EXPECT_TRUE(d <= raw + 1e-4f);
      }
}

/**
 * @brief Exact distance from a point to a twisted torus surface, by brute force.
 * @param p Query point.
 * @param R Major radius.
 * @param r Minor radius.
 * @param n Twist count.
 * @param A Twist amplitude.
 * @return The unsigned distance to the warped surface.
 * @details The surface is the union over theta of tube circles of radius r
 * centered at (R cos t, A sin(n t), R sin t) in the plane spanned by the radial
 * direction at t and the y axis, so the tube angle is solved in closed form and
 * only theta is sampled.
 */
constexpr double PI_DBL = 3.14159265358979323846;

inline double twisted_torus_distance(const Vector &p, double R, double r, int n,
                                     double A, int steps) {
  double best = 1e30;
  for (int i = 0; i < steps; ++i) {
    double t = 2.0 * PI_DBL * i / steps;
    double ct = std::cos(t), st = std::sin(t);
    double dx = p.x - R * ct, dy = p.y - A * std::sin(n * t), dz = p.z - R * st;
    double u = dx * ct + dz * st;
    double w = -dx * st + dz * ct;
    double rad = std::sqrt(u * u + dy * dy) - r;
    best = std::min(best, std::sqrt(rad * rad + w * w));
  }
  return best;
}

/**
 * @brief Verifies WarpedVolume::bounding_distance never over-estimates the
 *        distance to the warped surface, over randomized points and torus/twist
 *        parameter sets.
 * @details The fast path returns this bound directly, so an over-estimate would
 * let a sphere trace step through the surface. Covers twist 0, amplitude 0,
 * amplitude larger than the major radius, points on the y axis (s == 0), points
 * inside the tube and points hugging the surface.
 */
inline void test_warped_volume_bounding_distance_never_over_estimates() {
  struct Case {
    double R, r, A;
    int n;
  };
  const Case cases[] = {
      {1.0, 0.3, 0.2, 3},    // nominal
      {1.0, 0.3, 0.0, 3},    // zero amplitude
      {1.0, 0.3, 0.2, 0},    // zero twist
      {1.0, 0.31, 2.5, 5},   // amplitude well past the major radius
      {0.45, 0.14, 0.35, 2}, // Raymarch proportions, default twist
      {0.45, 0.14, 0.35, 8}, // Raymarch proportions, max twist
      {2.0, 0.05, 0.9, 7},   // thin tube
      {0.5, 0.45, 0.1, 1},   // fat tube
  };

  int violations = 0;
  hs::Pcg32 rng(0x5eed1337u);
  const auto next = [&rng]() {
    return static_cast<double>(rand_uniform(rng, 0.0f, 1.0f));
  };

  for (const Case &c : cases) {
    SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{
        SDF::Torus{static_cast<float>(c.R), static_cast<float>(c.r)},
        SDF::Warp::Twist{c.n, static_cast<float>(c.A),
                         static_cast<float>(c.R)}};
    const double reach = c.R + c.r + c.A + 1.0;

    for (int k = 0; k < 240; ++k) {
      Vector p;
      const int kind = k % 4;
      if (kind == 0) {
        p = Vector(0.0f, static_cast<float>((next() * 2 - 1) * reach), 0.0f);
      } else if (kind == 1 || kind == 2) {
        // On or inside the tube: radius scaled to at most the minor radius.
        const double t = next() * 2 * PI_DBL, ph = next() * 2 * PI_DBL;
        const double rr = c.r * (kind == 1 ? next() * 0.9 : 1.0);
        const double X = c.R + rr * std::cos(ph);
        p = Vector(
            static_cast<float>(X * std::cos(t)),
            static_cast<float>(rr * std::sin(ph) + c.A * std::sin(c.n * t)),
            static_cast<float>(X * std::sin(t)));
      } else {
        const double px = (next() * 2 - 1) * reach;
        const double py = (next() * 2 - 1) * reach;
        const double pz = (next() * 2 - 1) * reach;
        p = Vector(static_cast<float>(px), static_cast<float>(py),
                   static_cast<float>(pz));
      }

      const double bd = wv.bounding_distance(p);
      const double truth = twisted_torus_distance(p, c.R, c.r, c.n, c.A, 20000);
      if (bd - truth > 1e-4)
        ++violations;
      // The bound must also never exceed the slow path's raw warped distance.
      if (bd - static_cast<double>(wv.raw_distance(p)) > 1e-4)
        ++violations;
    }
  }
  HS_EXPECT_EQ(violations, 0);
}

/**
 * @brief Verifies the Lipschitz-corrected path returns raw/lipschitz on a
 *        near-surface outside point (not on the bounding fast-path).
 */
inline void test_warped_volume_distance_matches_lipschitz_correction() {
  SDF::Torus torus{1.0f, 0.3f};
  SDF::Warp::Twist tw{3, 0.2f, 1.0f};
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{torus, tw};

  // Just outside the outer rim at θ=0: small positive base distance lands off
  // the bounding fast-path and triggers the Lipschitz divide.
  Vector p(1.4f, 0.1f, 0.0f);
  float raw = wv.raw_distance(p);
  HS_EXPECT_TRUE(raw > 0.0f);
  auto ctx = tw.make_ctx(p);
  float lip = tw.lipschitz(p, ctx);
  HS_EXPECT_TRUE(lip > 1.0f);
  HS_EXPECT_NEAR(wv.distance(p), raw / lip, 1e-4f);
}

/** @brief Verifies Twist::correct_normal returns a unit vector and is identity at twist 0. */
inline void test_twist_correct_normal_unit_length() {
  Vector base_n = Vector(0.6f, 0.8f, 0.0f); // already unit

  SDF::Warp::Twist flat{0, 0.3f, 1.0f};
  Vector cf = flat.correct_normal(Vector(1, 0.2f, 0.5f), base_n,
                                  flat.make_ctx(Vector(1, 0.2f, 0.5f)));
  HS_EXPECT_VEC(cf, base_n, 1e-6f);

  SDF::Warp::Twist tw{4, 0.25f, 1.0f};
  for (float x = -1.0f; x <= 1.0f; x += 0.5f)
    for (float z = -1.0f; z <= 1.0f; z += 0.5f) {
      Vector p(x, 0.3f, z);
      Vector c = tw.correct_normal(p, base_n, tw.make_ctx(p));
      HS_EXPECT_NEAR(c.length(), 1.0f, 1e-4f);
    }
}

// ============================================================================
// Union — min of distances
// ============================================================================

/** @brief Verifies Union returns the min of member distances, picking the closest shape. */
inline void test_union_picks_closest_shape() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.1f);

  SDF::Union<SDF::Line, SDF::Line> u(la, lb);

  Vector mid_a = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = SDF::distance_of(u, mid_a);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-2f);

  Vector mid_b = ((Vector(-1, 0, 0) + Vector(0, 0, -1)) * 0.5f).normalized();
  auto r2 = SDF::distance_of(u, mid_b);
  HS_EXPECT_NEAR(r2.dist, -0.1f, 1e-2f);
}

// ============================================================================
// Subtract — max(A, -B)
// ============================================================================

/** @brief Verifies a point inside A but outside B stays inside the difference A - B. */
inline void test_subtract_inside_a_outside_b_remains_inside() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.2f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, lb);

  Vector mid_a = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = SDF::distance_of(s, mid_a);
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies a point inside both A and B becomes outside the difference A - B. */
inline void test_subtract_inside_both_becomes_outside() {
  // Same line for A and B → A - A is empty everywhere.
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, la);
  Vector mid = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = SDF::distance_of(s, mid);
  // max(dist(A), -dist(B)) = max(-0.1, 0.1) = 0.1 → outside.
  HS_EXPECT_NEAR(r.dist, 0.1f, 1e-3f);
}

/** @brief Verifies the AA size metric stays the minuend's when the subtrahend wins. */
inline void test_subtract_keeps_minuend_size_when_b_wins() {
  Vector p(1, 0, 0), q(0, 0, 1);
  SDF::Line la(p, q, 0.1f);
  SDF::Line lb(p, q, 0.4f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, lb);
  auto r = SDF::distance_of(s, ((p + q) * 0.5f).normalized());
  HS_EXPECT_NEAR(r.dist, 0.4f, 1e-3f);
  HS_EXPECT_NEAR(r.size, 0.1f, 1e-6f);
}

namespace sdf_subtract_detail {
/**
 * @brief Mock SDF shape that emits a fixed (possibly unsorted, multi-) interval list.
 * @details Exercises the CSG scanline interval paths independently of any real
 *   shape. Minimal surface: only is_solid and get_horizontal_intervals are
 *   touched by a combinator's ctor + interval path.
 */
struct MockIntervalShape {
  const std::vector<std::pair<float, float>>
      *ivs; /**< Interval list this mock replays. */
  static constexpr bool is_solid =
      true; /**< Marks the mock as a solid fill shape. */
  /**
   * @brief Emits the stored intervals to the scanline sink.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Out Interval-sink callable type taking (start, end).
   * @param out Sink invoked once per stored interval.
   * @return Always true: this mock definitively produces intervals.
   */
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out out) const {
    for (const auto &p : *ivs)
      out(p.first, p.second);
    return true;
  }
  /**
   * @brief Claims every row, matching the unconditional interval emission.
   * @tparam H Canvas height in rows.
   * @return Row bounds spanning the canvas.
   */
  template <int H> SDF::Bounds get_vertical_bounds() const {
    return {0, H - 1};
  }
};

/**
 * @brief Mock that falls back to a full-row scan.
 * @details Returning false means "I cannot produce intervals — caller must scan
 *   the whole row with distance()", NOT that the shape covers the row.
 */
struct MockFullWidthShape {
  static constexpr bool is_solid =
      true; /**< Marks the mock as a solid fill shape. */
  /**
   * @brief Declines to emit intervals, forcing a full-row distance scan.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Out Interval-sink callable type (unused).
   * @return Always false: caller must scan the whole row with distance().
   */
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out) const {
    return false;
  }
};

/**
 * @brief Mock that violates the interval protocol by emitting before falling back.
 * @details The protocol requires a false return to emit nothing. Nothing in the
 *   shape library does this; the mock exists so a CSG parent's fallback path can
 *   be checked for propagating the violation upward.
 */
struct MockEmitThenFullWidthShape {
  static constexpr bool is_solid =
      true; /**< Marks the mock as a solid fill shape. */
  /**
   * @brief Emits one span and then requests a full-row distance scan anyway.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Out Interval-sink callable type taking (start, end).
   * @param out Sink invoked once before the fallback.
   * @return Always false.
   */
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out out) const {
    out(10.0f, 30.0f);
    return false;
  }
};
} // namespace sdf_subtract_detail

/**
 * @brief Verifies a solid B's spans never carve the minuend's scanline emission.
 * @details A child's spans BOUND its coverage (a polygon or star emits its
 *   circumscribed cap, padded a pixel for AA), so differencing them would cull
 *   columns B does not cover. Every subtrahend is carved per pixel by
 *   max(A, -B) instead, leaving A's spans intact whatever B emits.
 */
inline void test_subtract_solid_b_leaves_the_minuend_uncarved() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{0.0f, 40.0f}, {60.0f, 100.0f}};
  std::vector<P> b_ivs = {{60.0f, 70.0f}, {20.0f, 30.0f}}; // unsorted, multi
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 40.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 60.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 100.0f, 1e-4f);
}

/**
 * @brief Verifies a star notch inside the subtrahend's bounding cap still gets scanned.
 * @details Star emits its circumscribed disc as one span, so carving that span
 *   out of the polygon would drop every column between the star's points -- all
 *   of them inside the difference -- plus the carve-edge AA fringe. The notch
 *   point below is outside the star, inside the polygon, and inside the
 *   difference, so its column must lie in an emitted span.
 */
inline void test_subtract_star_notch_columns_survive_the_carve() {
  using P = std::pair<float, float>;
  constexpr int W = 256, H = 128;
  // Axis +Z puts the shared center at column W/4, clear of the theta=0 seam.
  const Basis b{Vector(1, 0, 0), Vector(0, 0, 1), Vector(0, 1, 0)};
  constexpr float outer = 0.5f;              // star tip radius (radians)
  SDF::PlanarPolygon poly(b, 0.7f, 5, 0.0f); // apothem 0.566 > outer
  SDF::Star star(b, outer / (PI_F / 2.0f), 5, 0.0f);
  SDF::Subtract<SDF::PlanarPolygon, SDF::Star> s(poly, star);

  // Half a sector off a tip, just short of the tip radius: a notch.
  const float polar = 0.9f * outer, az = PI_F / 5.0f;
  Vector p = (b.v * std::cos(polar) +
              (b.u * std::cos(az) + b.w * std::sin(az)) * std::sin(polar))
                 .normalized();
  HS_EXPECT_TRUE(SDF::distance_of(star, p).dist > 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(poly, p).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(s, p).dist < 0.0f);

  float phi = std::acos(hs::clamp(p.y, -1.0f, 1.0f));
  int y = static_cast<int>(phi * (H + hs::H_OFFSET - 1) / PI_F + 0.5f);
  float theta = std::atan2(p.z, p.x);
  if (theta < 0.0f)
    theta += TWO_PI_F;
  const float col = theta * W / TWO_PI_F;

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<W, H>(
      y, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  bool covered = false;
  for (const auto &iv : out)
    covered = covered || (col >= iv.first && col <= iv.second);
  HS_EXPECT_TRUE(covered);
}

/**
 * @brief Verifies that when B removes nothing, every A interval passes through untouched.
 * @details Passthrough replays A's emission verbatim — no span dropped, merged
 *   or reordered. Consumers (scan_region's coalescer, a CSG parent's merge)
 *   order the list themselves.
 */
inline void test_subtract_empty_b_passes_a_through_verbatim() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{50.0f, 60.0f}, {0.0f, 10.0f}}; // unsorted
  std::vector<P> b_ivs = {};                              // empty → passthrough
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
  HS_EXPECT_NEAR(out[0].first, 50.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 60.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 10.0f, 1e-4f);
}

/**
 * @brief Verifies a subtrahend that cannot produce intervals costs the minuend nothing.
 * @details Subtract never consults B's intervals, so B's fallback neither widens
 *   the row to a full scan nor erases A: the row is handled with A's own spans
 *   and scan_region evaluates distance()=max(A,-B) inside them.
 */
inline void test_subtract_full_width_b_still_emits_the_minuend() {
  using P = std::pair<float, float>;
  using MockA = sdf_subtract_detail::MockIntervalShape;
  using MockB = sdf_subtract_detail::MockFullWidthShape;
  std::vector<P> a_ivs = {{0.0f, 100.0f}};
  MockA A{&a_ivs};
  MockB B;
  SDF::Subtract<MockA, MockB> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 100.0f, 1e-4f);
}

/**
 * @brief Verifies a minuend band straddling θ=0 is emitted seam-split into [0, W).
 * @details A emits the seam band in a negative wrap frame as [-10, 10]. The
 *   emission is normalized into [0, W), so the band reaches the sink as its two
 *   in-frame pieces rather than as a span the consumer must wrap itself.
 */
inline void test_subtract_seam_straddle_splits_the_minuend_into_frame() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-10.0f, 10.0f}};  // seam band, negative frame
  std::vector<P> b_ivs = {{246.0f, 266.0f}}; // same band, [W,2W) frame (W=256)
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
  HS_EXPECT_NEAR(out[0].first, 246.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 256.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 10.0f, 1e-4f);
}

/**
 * @brief Verifies a many-arc seam-straddling Subtract stays within the span-count bound.
 * @details normalize_intervals_to_range seam-splits the minuend into [0, W); a span
 *   crossing θ=0 becomes two, so the norm buffer is sized 2x. push_interval traps
 *   (fail-fast) if the post-split count exceeds that cap. Drive A with many disjoint
 *   arcs — two straddling the seam in a wrapped frame so they split — and B with arcs
 *   overlapping two of them, then assert every arc survives in [0, W) and that no
 *   trap fired (reaching this line at all).
 */
inline void test_subtract_many_arc_seam_split_within_bound() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  constexpr float W = 256.0f;

  // Twelve disjoint A arcs; the first two straddle the seam (negative / over-W
  // frame) so each splits into two on normalization (14 norm spans, well under 64).
  std::vector<P> a_ivs = {{-6.0f, 4.0f},    {250.0f, 262.0f}, // seam straddlers
                          {20.0f, 30.0f},   {40.0f, 50.0f},   {60.0f, 70.0f},
                          {80.0f, 90.0f},   {100.0f, 110.0f}, {120.0f, 130.0f},
                          {140.0f, 150.0f}, {160.0f, 170.0f}, {180.0f, 190.0f},
                          {200.0f, 210.0f}};
  std::vector<P> b_ivs = {{45.0f, 55.0f}, {125.0f, 135.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  // Twelve arcs, two of them split at the seam. Reaching here means
  // push_interval never trapped.
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(14));
  for (size_t i = 0; i < out.size(); ++i) {
    HS_EXPECT_TRUE(out[i].first >= 0.0f && out[i].second <= W);
    HS_EXPECT_TRUE(out[i].first < out[i].second);
  }
  // The arcs B overlaps come through whole.
  HS_EXPECT_NEAR(out[5].first, 40.0f, 1e-4f);
  HS_EXPECT_NEAR(out[5].second, 50.0f, 1e-4f);
  HS_EXPECT_NEAR(out[9].first, 120.0f, 1e-4f);
  HS_EXPECT_NEAR(out[9].second, 130.0f, 1e-4f);
}

// ============================================================================
// Intersection — max(A, B)
// ============================================================================

/** @brief Verifies Intersection is inside only where both children are inside. */
inline void test_intersection_requires_both_inside() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.3f);
  SDF::Line lb(Vector(1, 0, 0), Vector(0, 1, 0), 0.3f);
  SDF::Intersection<SDF::Line, SDF::Line> inter(la, lb);

  // Endpoint a is on both arcs.
  Vector a(1, 0, 0);
  auto r = SDF::distance_of(inter, a);
  HS_EXPECT_TRUE(r.dist < 0.0f);

  Vector far_pt(-1, 0, 0);
  auto r2 = SDF::distance_of(inter, far_pt);
  HS_EXPECT_TRUE(r2.dist > 0.0f);
}

/**
 * @brief Verifies an UNSORTED multi-interval child still yields a start-sorted intersection.
 * @details Intersection's merge-sweep advances the child lists in start order,
 *   so it sorts the seam-split lists itself before sweeping them.
 */
inline void test_intersection_unsorted_child_yields_sorted_result() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{0.0f, 100.0f}};
  std::vector<P> b_ivs = {{60.0f, 80.0f}, {20.0f, 40.0f}}; // unsorted, multi
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Intersection<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  // [0,100] ∩ {[20,40],[60,80]} = [20,40],[60,80], start-sorted.
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
  HS_EXPECT_NEAR(out[0].first, 20.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 40.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 60.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 80.0f, 1e-4f);
}

/**
 * @brief Verifies a full-width child intersected with the other replays the other's intervals.
 * @details When one child falls back to a full-width scan, the intersection is
 *   just the other child's intervals (replayed from the buffer already collected).
 *   Pins the equivalence in both orientations and the both-fall-back full-scan case.
 */
inline void test_intersection_full_width_child_replays_other() {
  using P = std::pair<float, float>;
  using MockI = sdf_subtract_detail::MockIntervalShape;
  using MockF = sdf_subtract_detail::MockFullWidthShape;
  std::vector<P> ivs = {{60.0f, 80.0f},
                        {20.0f, 40.0f}}; // multi, emission order
  MockI shape{&ivs};
  MockF full;

  {
    SDF::Intersection<MockF, MockI> s(full, shape);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_TRUE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
    HS_EXPECT_NEAR(out[0].first, 60.0f, 1e-4f);
    HS_EXPECT_NEAR(out[1].first, 20.0f, 1e-4f);
  }

  // Symmetric: A has intervals, B full-width.
  {
    SDF::Intersection<MockI, MockF> s(shape, full);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_TRUE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
    HS_EXPECT_NEAR(out[0].first, 60.0f, 1e-4f);
    HS_EXPECT_NEAR(out[1].first, 20.0f, 1e-4f);
  }

  // Both fall back -> full-scan fallback (return false).
  {
    SDF::Intersection<MockF, MockF> s(full, full);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_FALSE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
  }
}

/**
 * @brief Verifies Intersection emits nothing whenever it requests a full-row scan.
 * @details The interval protocol pairs a false return with an empty emission: the
 *   caller walks every column and never reads the buffer, so a span emitted first
 *   is dropped or shaded twice. Intersection replays a buffered child, which makes
 *   it the one combinator that could pair both. Driven by a child that itself
 *   violates the protocol, so the check does not rest on leaf behaviour.
 */
inline void test_intersection_full_scan_emits_no_spans() {
  using P = std::pair<float, float>;
  using MockE = sdf_subtract_detail::MockEmitThenFullWidthShape;
  using MockF = sdf_subtract_detail::MockFullWidthShape;
  MockE emitting;
  MockF full;

  {
    SDF::Intersection<MockF, MockE> s(full, emitting);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_FALSE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
  }

  {
    SDF::Intersection<MockE, MockF> s(emitting, full);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_FALSE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
  }
}

/**
 * @brief Verifies a seam-straddling band shared by both children is intersected across wrap frames.
 * @details A emits the seam band as [-10, 10]; B emits the SAME physical band as
 *   [W-10, W+10]. A raw-coordinate overlap test sees no shared columns and
 *   under-reports the seam. After normalizing both into [0, W) the bands coincide
 *   and the intersection is the full shared band, split at the seam.
 */
inline void test_intersection_seam_straddle_overlaps_across_wrap_frames() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-10.0f, 10.0f}};  // seam band, negative frame
  std::vector<P> b_ivs = {{246.0f, 266.0f}}; // same band, [W,2W) frame (W=256)
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Intersection<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  // Both normalize to {[246,256],[0,10]}; the intersection is the same band.
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(2));
}

// ============================================================================
// SmoothUnion — blends at the boundary
// ============================================================================

/** @brief Verifies that away from the blend zone, SmoothUnion's distance equals the hard Union's. */
inline void test_smooth_union_matches_union_far_from_boundary() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.1f);
  SDF::Union<SDF::Line, SDF::Line> u(la, lb);
  SDF::SmoothUnion<SDF::Line, SDF::Line> su(la, lb, /*k*/ 0.05f);

  Vector mid_a = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r_hard = SDF::distance_of(u, mid_a);
  auto r_soft = SDF::distance_of(su, mid_a);
  HS_EXPECT_NEAR(r_hard.dist, r_soft.dist, 1e-3f);
}

/**
 * @brief Verifies the cubic smin blend term inside the blend band.
 * @details The far-from-boundary test only re-checks the hard-union min (the
 *   blend term m is zero there). This exercises the smin core: at a point
 *   roughly equidistant from both children (|dA - dB| < k) the smooth distance
 *   must dip strictly below min(dA, dB) by the expected m, and outside the band
 *   it must collapse back to the hard min (m == 0).
 */
inline void test_smooth_union_blends_inside_band() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.1f);
  const float k = 0.5f;
  SDF::SmoothUnion<SDF::Line, SDF::Line> su(la, lb, k);

  // Both arcs lie in the y=0 plane, so the north pole is equidistant from both
  // (|dA - dB| ≈ 0), maximizing the cubic blend (h == 1, m == k/6).
  Vector p(0, 1, 0);
  float dA = SDF::distance_of(la, p).dist;
  float dB = SDF::distance_of(lb, p).dist;
  HS_EXPECT_TRUE(std::abs(dA - dB) < k);

  float h = std::max(k - std::abs(dA - dB), 0.0f) / k;
  float m = h * h * h * k * (1.0f / 6.0f);
  float soft = SDF::distance_of(su, p).dist;
  HS_EXPECT_NEAR(soft, std::min(dA, dB) - m, 1e-4f);
  HS_EXPECT_LT(soft, std::min(dA, dB) - 1e-4f);
  HS_EXPECT_NEAR(std::min(dA, dB) - soft, m, 1e-4f);
  // Independent of the cubic: the smooth min never rises above the hard min,
  // and k/6 is the deepest it can dip below it anywhere.
  HS_EXPECT_LE(soft, std::min(dA, dB) + 1e-5f);
  HS_EXPECT_GE(soft, std::min(dA, dB) - k * (1.0f / 6.0f) - 1e-5f);

  // Outside the band the blend vanishes and collapses to the hard min.
  Vector q = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  float qA = SDF::distance_of(la, q).dist;
  float qB = SDF::distance_of(lb, q).dist;
  HS_EXPECT_TRUE(std::abs(qA - qB) >= k);
  HS_EXPECT_NEAR(SDF::distance_of(su, q).dist, std::min(qA, qB), 1e-5f);
}

/**
 * @brief Verifies SmoothUnion solidity follows its children.
 * @details Children must share solidity (enforced by static_assert): two strokes
 *   take the soft falloff path, two solids the hard 1-px silhouette path.
 */
inline void test_smooth_union_solidity_follows_children() {
  static_assert(!SDF::SmoothUnion<SDF::Line, SDF::Line>::is_solid,
                "two strokes -> falloff path");
  static_assert(
      SDF::SmoothUnion<SDF::PlanarPolygon, SDF::PlanarPolygon>::is_solid,
      "two solids -> silhouette path");
}

/**
 * @brief Verifies the blendability trait tracks which shapes clamp to the far
 *   sentinel, and that a combinator blends only when every child does.
 * @details Ring, DistortedRing, FlatDistortedRing and Face report dist = 100
 *   outside their reject band, so both children read the sentinel across the
 *   weld and SmoothUnion collapses to Union; SmoothUnion static_asserts the
 *   trait to reject those instantiations at compile time.
 */
inline void test_sentinel_clampers_are_not_blendable() {
  static_assert(!SDF::blends_smoothly<SDF::Ring>);
  static_assert(!SDF::blends_smoothly<SDF::DistortedRing>);
  static_assert(!SDF::blends_smoothly<SDF::FlatDistortedRing>);
  static_assert(!SDF::blends_smoothly<SDF::Face>);
  static_assert(SDF::blends_smoothly<SDF::Line>);
  static_assert(SDF::blends_smoothly<SDF::PlanarPolygon>);
  static_assert(SDF::blends_smoothly<SDF::SphericalPolygon>);
  static_assert(SDF::blends_smoothly<SDF::Star>);
  static_assert(SDF::blends_smoothly<SDF::Flower>);
  static_assert(!SDF::blends_smoothly<SDF::Union<SDF::Line, SDF::Ring>>);
  static_assert(SDF::blends_smoothly<SDF::Union<SDF::Line, SDF::Line>>);
  static_assert(!SDF::blends_smoothly<SDF::AngularRepeat<SDF::Ring>>);
  static_assert(SDF::blends_smoothly<SDF::AngularRepeat<SDF::Line>>);
}

/**
 * @brief Verifies the CSG combinators reject a temporary child.
 * @details Every combinator holds its children by reference and reads them on
 *   each pixel, so a temporary bound at construction dangles for the whole
 *   render. The leaves delete the rvalue-Basis overloads for the same reason;
 *   an accepted rvalue here is the same defect one level up.
 */
inline void test_csg_combinators_reject_temporary_children() {
  using L = SDF::Line;
  using Poly = SDF::PlanarPolygon;
  static_assert(std::is_constructible_v<SDF::Union<L, L>, L &, L &>);
  static_assert(!std::is_constructible_v<SDF::Union<L, L>, L &&, L &>);
  static_assert(!std::is_constructible_v<SDF::Union<L, L>, L &, L &&>);
  static_assert(!std::is_constructible_v<SDF::Union<L, L>, L &&, L &&>);

  static_assert(
      std::is_constructible_v<SDF::SmoothUnion<L, L>, L &, L &, float>);
  static_assert(
      !std::is_constructible_v<SDF::SmoothUnion<L, L>, L &&, L &, float>);
  static_assert(
      !std::is_constructible_v<SDF::SmoothUnion<L, L>, L &, L &&, float>);

  static_assert(std::is_constructible_v<SDF::Subtract<Poly, L>, Poly &, L &>);
  static_assert(!std::is_constructible_v<SDF::Subtract<Poly, L>, Poly &&, L &>);
  static_assert(!std::is_constructible_v<SDF::Subtract<Poly, L>, Poly &, L &&>);

  static_assert(std::is_constructible_v<SDF::Intersection<L, L>, L &, L &>);
  static_assert(!std::is_constructible_v<SDF::Intersection<L, L>, L &&, L &>);
  static_assert(!std::is_constructible_v<SDF::Intersection<L, L>, L &, L &&>);

  // AngularRepeat copies its axis, so only the child rejects a temporary.
  static_assert(std::is_constructible_v<SDF::AngularRepeat<L>, L &, int>);
  static_assert(!std::is_constructible_v<SDF::AngularRepeat<L>, L &&, int>);
  static_assert(
      std::is_constructible_v<SDF::AngularRepeat<L>, L &, int, Vector &&>);
  static_assert(
      !std::is_constructible_v<SDF::AngularRepeat<L>, L &&, int, Vector &&>);

  HS_EXPECT_FALSE((std::is_constructible_v<SDF::Union<L, L>, L &&, L &>));
  HS_EXPECT_FALSE((std::is_constructible_v<SDF::AngularRepeat<L>, L &&, int>));
}

/**
 * @brief Verifies Union coalesces two overlapping child intervals into one span.
 * @details merge_intervals sorts by start and welds spans that touch or overlap;
 *   the children emit overlapping bands (A [0,40], B [30,70]) that must collapse
 *   to a single [0,70] — emitting both would double-paint the [30,40] overlap.
 */
inline void test_union_merges_overlapping_intervals() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{0.0f, 40.0f}};
  std::vector<P> b_ivs = {{30.0f, 70.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Union<Mock, Mock> u(A, B);

  std::vector<P> out;
  bool ok = u.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 70.0f, 1e-4f);
}

/**
 * @brief Verifies Union welds two overlapping seam-straddling spans into one.
 * @details Both children emit bands straddling θ=0 in the same negative frame
 *   (A [-10,6], B [2,12]); their overlap must coalesce to a single [-10,12] span,
 *   not two duplicate spans the downstream seam-split would then double-paint.
 */
inline void test_union_seam_straddle_merges_overlapping_intervals() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-10.0f, 6.0f}};
  std::vector<P> b_ivs = {{2.0f, 12.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Union<Mock, Mock> u(A, B);

  std::vector<P> out;
  bool ok = u.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(out[0].first, -10.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 12.0f, 1e-4f);
}

/**
 * @brief Verifies three- and four-way nested Unions of real leaves compile and
 *        emit every child arc.
 * @details The nesting depth a combinator admits is gated by sdf_max_spans; a
 *   leaf pinned to its true emission (2 for the annular-band rings) keeps a
 *   four-deep union at 8 spans, far under the 2*INTERVAL_SPAN_CAP accumulator.
 *   Sizing the types alone would prove the static_asserts pass, so the test also
 *   runs the scanline path and counts the merged arcs: four coaxial rings at
 *   disjoint radii cross an equatorial row in 8 disjoint spans, which is also
 *   the bound the trait reports.
 */
inline void test_nested_union_emits_every_child_arc() {
  using P = std::pair<float, float>;
  constexpr int W = 256, H = 128;
  // Ring axis along +X so the row math has a non-degenerate horizontal
  // projection (a +Y axis full-row scans instead).
  const Basis b{Vector(0, 1, 0), Vector(1, 0, 0), Vector(0, 0, 1)};
  SDF::Ring r1(b, 0.4f, 0.05f), r2(b, 0.6f, 0.05f);
  SDF::Ring r3(b, 0.8f, 0.05f), r4(b, 1.0f, 0.05f);

  using U2 = SDF::Union<SDF::Ring, SDF::Ring>;
  using U3 = SDF::Union<U2, SDF::Ring>;
  using U4 = SDF::Union<U2, U2>;
  static_assert(SDF::sdf_max_spans<SDF::Ring>::value == 2);
  static_assert(SDF::sdf_max_spans<U3>::value == 6);
  static_assert(SDF::sdf_max_spans<U4>::value == 8);
  // Nesting under the frame-sensitive combinators too: each child must fit one
  // IntervalBuffer, which a 32-span leaf bound could never satisfy.
  static_assert(sizeof(SDF::Subtract<U2, SDF::Ring>) > 0);
  static_assert(sizeof(SDF::Intersection<U2, U2>) > 0);

  U2 u_lo(r1, r2), u_hi(r3, r4);
  U3 u3(u_lo, r3);
  U4 u4(u_lo, u_hi);

  std::vector<P> out3, out4;
  bool ok3 = u3.get_horizontal_intervals<W, H>(
      H / 2, [&](float st, float en) { out3.push_back({st, en}); });
  bool ok4 = u4.get_horizontal_intervals<W, H>(
      H / 2, [&](float st, float en) { out4.push_back({st, en}); });
  HS_EXPECT_TRUE(ok3);
  HS_EXPECT_TRUE(ok4);

  HS_EXPECT_EQ(out3.size(), static_cast<size_t>(6));
  HS_EXPECT_EQ(out4.size(), static_cast<size_t>(8));
  for (size_t i = 1; i < out4.size(); ++i)
    HS_EXPECT_TRUE(out4[i - 1].second < out4[i].first);
}

/**
 * @brief Verifies SmoothUnion's k-padded union welds overlapping seam-straddling spans.
 * @details Each child interval is inflated by pad_px = k·W/(2π·sinφ) before the
 *   merge, so a tiny k keeps the bounds near the raw spans while still routing
 *   through the pad path; two overlapping seam-straddling bands must coalesce to
 *   one padded span rather than two overlapping spans. The pad mirrors the
 *   latitude-scaled production formula at the sampled row.
 */
inline void test_smooth_union_seam_straddle_merges_padded_intervals() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  constexpr int W = 256, H = 128;
  init_geometry_luts<W,
                     H>(); // fill sin_phi; scan_region does this in production
  const int row = H / 2;   // equatorial row: sinφ ≈ 1, pad ≈ k·W/(2π)
  const float k = 0.02f;
  const float sin_phi = TrigLUT<W, H>::sin_phi[row];
  const float pad =
      std::min(k * W / (2.0f * PI_F) / sin_phi, static_cast<float>(W));
  std::vector<P> a_ivs = {{-10.0f, 6.0f}};
  std::vector<P> b_ivs = {{2.0f, 12.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::SmoothUnion<Mock, Mock> su(A, B, k);

  std::vector<P> out;
  bool ok = su.get_horizontal_intervals<W, H>(
      row, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(1));
  HS_EXPECT_NEAR(out[0].first, -10.0f - pad, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 12.0f + pad, 1e-4f);
  // Independent of the pad formula: the weld must cover both raw spans and
  // widen past their union in both directions.
  HS_EXPECT_LT(out[0].first, -10.0f);
  HS_EXPECT_GT(out[0].second, 12.0f);
}

/**
 * @brief Verifies the weld pad widens toward the poles by the 1/sinφ factor.
 * @details A point interval's emitted span is exactly twice the row pad. Sampling
 *   a near-pole row (small sinφ) against the equator (sinφ ≈ 1) must yield a
 *   strictly wider span, so near-pole welds are not clipped to the equatorial pad.
 */
inline void test_smooth_union_pad_widens_toward_pole() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  constexpr int W = 256, H = 128;
  init_geometry_luts<W,
                     H>(); // fill sin_phi; scan_region does this in production
  const float k = 0.05f;
  std::vector<P> ivs = {{100.0f, 100.0f}}; // a point; only the pad sets width
  Mock A{&ivs}, B{&ivs};
  SDF::SmoothUnion<Mock, Mock> su(A, B, k);

  auto span_at = [&](int row) {
    std::vector<P> out;
    su.get_horizontal_intervals<W, H>(
        row, [&](float st, float en) { out.push_back({st, en}); });
    return out.empty() ? 0.0f : out[0].second - out[0].first;
  };
  HS_EXPECT_GT(span_at(1), span_at(H / 2));
}

// ============================================================================
// AngularRepeat — folds azimuth into N sectors
// ============================================================================

/** @brief Verifies AngularRepeat agrees with the base shape in the canonical (zero-angle) sector. */
inline void test_angular_repeat_matches_base_at_zero_angle() {
  SDF::Line ln(Vector(1, 0, 0), Vector(0.7071f, 0, 0.7071f), 0.05f);
  SDF::AngularRepeat<SDF::Line> rep(ln, /*reps*/ 4, Vector(0, 1, 0));

  Vector mid =
      ((Vector(1, 0, 0) + Vector(0.7071f, 0, 0.7071f)) * 0.5f).normalized();
  auto r_base = SDF::distance_of(ln, mid);
  auto r_rep = SDF::distance_of(rep, mid);
  HS_EXPECT_TRUE(r_rep.dist < 0.0f);
  HS_EXPECT_NEAR(r_rep.dist, r_base.dist, 1e-3f);
}

/** @brief Verifies AngularRepeat folds a line in the canonical sector into a folded copy. */
inline void test_angular_repeat_creates_copies() {
  SDF::Line ln(Vector(1, 0, 0), Vector(0.7071f, 0, 0.7071f), 0.05f);
  SDF::AngularRepeat<SDF::Line> rep(ln, 4, Vector(0, 1, 0));

  // Rotate the midpoint by one sector (90° around Y) onto a folded copy.
  Vector mid =
      ((Vector(1, 0, 0) + Vector(0.7071f, 0, 0.7071f)) * 0.5f).normalized();
  Quaternion q90 = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  Vector mid_rot = rotate(mid, q90);

  auto r = SDF::distance_of(rep, mid_rot);
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/**
 * @brief Verifies AngularRepeat's child UV (t) is sector-local, not global.
 * @details distance() folds p into one sector before evaluating the child, so
 *   the child's azimuthal t resets every sector. A point and its copy one full
 *   sector away fold to the same point and thus share a t, even though the
 *   un-repeated ring reports two distinct global azimuths for them. This pins
 *   the documented sector-local convention so a future change is caught.
 */
inline void test_angular_repeat_t_is_sector_local() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);
  const int reps = 4;
  SDF::AngularRepeat<SDF::Ring> rep(ring, reps, Vector(0, 1, 0));

  // A point at 30° azimuth (inside sector 0, off a boundary) and its copy one
  // sector away.
  float az = PI_F / 6.0f;
  Vector p(cosf(az), 0.0f, sinf(az));
  Quaternion q_sector = make_rotation(Vector(0, 1, 0), 2 * PI_F / reps);
  Vector p2 = rotate(p, q_sector);

  // The un-repeated ring sees two global azimuths one sector apart (sign is
  // handedness-dependent, so compare the wrapped gap).
  float t_global_1 = SDF::distance_of(ring, p).t;
  float t_global_2 = SDF::distance_of(ring, p2).t;
  float dg = t_global_2 - t_global_1;
  dg -= floorf(dg);
  HS_EXPECT_NEAR(std::min(dg, 1.0f - dg), 1.0f / reps, 1e-3f);

  // The repeated shape folds both into the same sector → identical sector-local t.
  float t_rep_1 = SDF::distance_of(rep, p).t;
  float t_rep_2 = SDF::distance_of(rep, p2).t;
  HS_EXPECT_NEAR(t_rep_1, t_rep_2, 1e-3f);
  HS_EXPECT_NEAR(t_rep_1, t_global_1, 1e-3f);
}

// ============================================================================
// Interval-cull conservativeness  ("interval cull == full-row scan")
//
// The cull must be conservative: it may over-visit but must never drop a pixel
// belonging to the shape. Invariant: a pixel clearly inside (distance().dist <
// -pixel_width, one pixel deep into body or stroke band) is among those
// scan_region visits. The one-pixel margin keeps the test off two soft edges the
// cull does NOT promise — a stroke's outer ~5% rim and a solid's sub-pixel AA
// halo — and clears the ~0.004 rad fast-trig noise.
// ============================================================================

/**
 * @brief Records every pixel scan_region visits for a shape, as rasterize() drives it.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Shape SDF shape type providing the cull and interval interface.
 * @param shape Shape whose culled coverage is being captured.
 * @param visited Output flag grid (W*H), set to 1 for each visited pixel.
 * @details Drives the full canvas with no clip, mirroring the rasterizer.
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
      [&](int wx, int y, const Vector &, int run) {
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
 * @param label Caller-identifying label; every failure here reports it plus the
 *   pixel, since __func__ names this shared helper for all of them.
 * @return Count of interior pixels (dist < -pixel_width) found.
 * @details Interior pixels are found by a brute-force full-canvas exact distance
 *   scan. Asserts the case is non-trivial (at least one interior pixel) so no
 *   shape silently contributes zero coverage.
 */
template <int W, int H, typename Shape>
inline int expect_cull_covers_interior(const Shape &shape, const char *label) {
  HS_CONTEXT(label);
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();
  std::vector<uint8_t> visited;
  cull_visited<W, H>(shape, visited);

  const float *cos_theta =
      TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  const float pixel_width = 2.0f * PI_F / W;
  int interior = 0;
  for (int y = 0; y < H; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    for (int x = 0; x < W; ++x) {
      Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
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
 * @param label Caller-identifying label; every failure here reports it plus the
 *   pixel, since __func__ names this shared helper for all of them.
 * @return Count of paintable pixels (dist < pixel_width) found.
 * @details The stricter sibling of expect_cull_covers_interior: the shade path
 *   returns early only at dist >= pixel_width, so every pixel under that
 *   threshold carries non-zero coverage and dropping one steps the fringe to
 *   zero mid-ramp.
 */
template <int W, int H, typename Shape>
inline int expect_cull_covers_fringe(const Shape &shape, const char *label) {
  HS_CONTEXT(label);
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();
  std::vector<uint8_t> visited;
  cull_visited<W, H>(shape, visited);

  const float *cos_theta =
      TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  const float pixel_width = 2.0f * PI_F / W;
  int paintable = 0;
  for (int y = 0; y < H; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    for (int x = 0; x < W; ++x) {
      Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
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
 * @details All three fold a sector onto one edge and read a shallow radial
 *   gradient at the tips (|nx| = 0.309 at 5 star points, cos(PI/sides) at a
 *   polygon vertex), which without the circumscribed-disc clamp in distance()
 *   ramps coverage well past the one-pixel cap pad and steps it to zero there.
 *   The grid spans pole, equator and oblique axes at sub-pixel through
 *   near-hemisphere radii.
 */
inline void test_star_polygon_cull_covers_aa_fringe() {
  constexpr int W = 96, H = 48;
  const Vector axes[] = {Vector(0, 1, 0), Vector(1, 0, 0),
                         Vector(0.3f, -0.8f, 0.5f)};
  int total = 0;
  for (const Vector &axis : axes) {
    Basis basis = make_basis(Quaternion(), axis);
    for (float radius : {0.03f, 0.3f, 0.9f}) {
      for (int sides : {5, 8}) {
        SDF::Star star(basis, radius, sides, 0.0f);
        total += expect_cull_covers_fringe<W, H>(star, "star");

        SDF::PlanarPolygon poly(basis, radius * (PI_F / 2.0f), sides, 0.0f);
        total += expect_cull_covers_fringe<W, H>(poly, "planar polygon");

        SDF::SphericalPolygon sphpoly(basis, radius, sides, 0.0f);
        total += expect_cull_covers_fringe<W, H>(sphpoly, "spherical polygon");
      }
    }
  }
  HS_EXPECT_GT(total, 1000);
}

/** @brief Verifies the interval cull covers every interior pixel across an orientation/radius grid. */
inline void test_cull_covers_interior_over_orientation_grid() {
  constexpr int W = 96, H = 48;

  // Poles, equator, and oblique tilts.
  const Vector axes[] = {Vector(0, 1, 0),          Vector(0, -1, 0),
                         Vector(1, 0, 0),          Vector(0, 0, 1),
                         Vector(1, 1, 0.4f),       Vector(-0.5f, 0.7f, -0.6f),
                         Vector(0.3f, -0.8f, 0.5f)};

  for (const Vector &axis : axes) {
    Basis basis = make_basis(Quaternion(), axis);
    for (float radius : {0.3f, 0.6f, 0.9f}) {
      SDF::Ring ring(basis, radius, /*thickness=*/0.25f);
      expect_cull_covers_interior<W, H>(ring, "ring");

      SDF::SphericalPolygon spoly(basis, radius, /*sides=*/5, 0.0f);
      expect_cull_covers_interior<W, H>(spoly, "spherical polygon");

      SDF::Star star(basis, radius, /*sides=*/5, 0.0f);
      expect_cull_covers_interior<W, H>(star, "star");

      SDF::PlanarPolygon ppoly(basis, /*circumradius=*/radius, /*sides=*/6,
                               0.0f);
      expect_cull_covers_interior<W, H>(ppoly, "planar polygon");

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
  Basis basis = equator_basis();
  SDF::Ring ring(basis, /*radius=*/0.6f, /*thickness=*/0.25f);

  auto bounds = ring.get_vertical_bounds<H>();
  HS_EXPECT_GT(bounds.y_min, 0);
  HS_EXPECT_LT(bounds.y_max, H - 1);
  expect_cull_covers_interior<W, H>(ring, "pole-axis ring");
}

/**
 * @brief Verifies the Intersection interval cull covers every interior pixel of
 *        a real leaf pair.
 * @details The seam-split normalization and the merge sweep are otherwise driven
 *   only by mock span lists. The last pose centers both polygons on +X, so their
 *   overlap straddles theta = 0 and each child can emit its span in a different
 *   wrap frame.
 */
inline void test_intersection_cull_covers_interior_over_polygon_pairs() {
  constexpr int W = 96, H = 48;
  using Poly = SDF::PlanarPolygon;

  struct Pose {
    Vector axis_a, axis_b;
    float radius_a, radius_b;
  };
  const Pose poses[] = {
      {Vector(0, 0, 1), Vector(0.3f, 0.2f, 1.0f), 0.8f, 0.5f},
      {Vector(0, 1, 0), Vector(0.25f, 1.0f, -0.15f), 0.7f, 0.45f},
      {Vector(-0.4f, 0.6f, 0.7f), Vector(-0.2f, 0.75f, 0.6f), 0.9f, 0.6f},
      {Vector(1, 0, 0), Vector(1.0f, 0.15f, 0.2f), 0.8f, 0.5f},
  };

  for (const Pose &pose : poses) {
    Basis basis_a = make_basis(Quaternion(), pose.axis_a);
    Basis basis_b = make_basis(Quaternion(), pose.axis_b);
    Poly poly_a(basis_a, pose.radius_a, /*sides=*/6, 0.0f);
    Poly poly_b(basis_b, pose.radius_b, /*sides=*/5, 0.4f);

    SDF::Intersection<Poly, Poly> both(poly_a, poly_b);
    expect_cull_covers_interior<W, H>(both, "intersection leaf pair");
  }
}

/**
 * @brief Verifies the Subtract interval cull covers every interior pixel of a
 *        real leaf pair.
 * @details The rest of the Subtract suite drives mock span lists, which stand in
 *   for coverage exactly; a real leaf emits its circumscribed cap padded a pixel
 *   for AA, so its spans only bound it. Subtract answers that by emitting the
 *   minuend's spans and carving per pixel, and this sweeps the whole canvas to
 *   pin that no pixel the carve leaves solid falls outside them. A star
 *   subtrahend is the hostile case: its cap covers the notches the difference
 *   keeps. The last pose centers both on +X, so the minuend's spans straddle
 *   theta = 0.
 */
inline void test_subtract_cull_covers_interior_over_leaf_pairs() {
  constexpr int W = 96, H = 48;
  using Poly = SDF::PlanarPolygon;
  using Star = SDF::Star;

  struct Pose {
    Vector axis_a, axis_b;
    float radius_a; /**< Polygon circumradius in radians. */
    float radius_b; /**< Star tip radius, normalized (x PI/2 for radians). */
  };
  const Pose poses[] = {
      {Vector(0, 0, 1), Vector(0.2f, 0.1f, 1.0f), 0.9f, 0.35f},
      {Vector(0, 1, 0), Vector(0.2f, 1.0f, -0.1f), 0.8f, 0.3f},
      {Vector(-0.4f, 0.6f, 0.7f), Vector(-0.3f, 0.65f, 0.7f), 1.0f, 0.4f},
      {Vector(1, 0, 0), Vector(1.0f, 0.1f, 0.15f), 0.9f, 0.35f},
  };

  for (const Pose &pose : poses) {
    Basis basis_a = make_basis(Quaternion(), pose.axis_a);
    Basis basis_b = make_basis(Quaternion(), pose.axis_b);
    Poly poly(basis_a, pose.radius_a, /*sides=*/6, 0.0f);
    Star star(basis_b, pose.radius_b, /*sides=*/5, 0.4f);

    SDF::Subtract<Poly, Star> carved(poly, star);
    expect_cull_covers_interior<W, H>(carved, "subtract leaf pair");
  }
}

/**
 * @brief Verifies the SmoothUnion interval cull covers every interior pixel of a
 *        real leaf pair.
 * @details The weld bulges the surface outside both children, so the cull rests
 *   entirely on the k pad — and the pad is an equatorial column count divided by
 *   sin(phi), which the mock span lists exercise only at fixed rows. Real leaves
 *   put the weld at every latitude the poses reach. The last pose centers both
 *   on +X so the padded spans straddle theta = 0.
 */
inline void test_smooth_union_cull_covers_interior_over_leaf_pairs() {
  constexpr int W = 96, H = 48;
  using Poly = SDF::PlanarPolygon;

  struct Pose {
    Vector axis_a, axis_b;
    float radius_a, radius_b;
  };
  const Pose poses[] = {
      {Vector(0, 0, 1), Vector(0.6f, 0.2f, 1.0f), 0.5f, 0.4f},
      {Vector(0, 1, 0), Vector(0.5f, 1.0f, -0.2f), 0.45f, 0.35f},
      {Vector(-0.4f, 0.6f, 0.7f), Vector(0.1f, 0.9f, 0.4f), 0.6f, 0.5f},
      {Vector(1, 0, 0), Vector(1.0f, 0.35f, 0.2f), 0.5f, 0.4f},
  };

  for (const Pose &pose : poses) {
    Basis basis_a = make_basis(Quaternion(), pose.axis_a);
    Basis basis_b = make_basis(Quaternion(), pose.axis_b);
    Poly poly_a(basis_a, pose.radius_a, /*sides=*/6, 0.0f);
    Poly poly_b(basis_b, pose.radius_b, /*sides=*/5, 0.4f);

    SDF::SmoothUnion<Poly, Poly> welded(poly_a, poly_b, /*k=*/0.15f);
    expect_cull_covers_interior<W, H>(welded, "smooth union leaf pair");
  }
}

/**
 * @brief Verifies the SmoothUnion cull scans the rows past both children's
 *        spans.
 * @details Welding a polygon to itself puts the smin at its peak everywhere, so
 *   the surface dilates by k/6 in every direction — past the last row either
 *   child emits a span for. Those rows carry no interval to pad, and skipping
 *   them clips the weld's polar fringe; they must request a full scan instead.
 *   A row past the blend reach still holds no surface and stays culled.
 */
inline void test_smooth_union_scans_rows_past_both_children() {
  constexpr int W = 96, H = 48;
  init_geometry_luts<W, H>();
  using Poly = SDF::PlanarPolygon;
  Basis basis = make_basis(Quaternion(), Vector(0.2f, 1.0f, 0.1f));
  Poly poly(basis, /*radius=*/0.5f, /*sides=*/6, 0.0f);

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
 * @details A non-Y axis sweeps the folded copies through latitudes the un-repeated
 *   child never occupies, so the child's vertical band no longer bounds them.
 *   get_vertical_bounds must fall back to the full canvas; forwarding the child's
 *   narrow band would row-clip every off-band copy out of the scan and drop it.
 */
inline void test_angular_repeat_non_y_axis_cull_covers_copies() {
  constexpr int W = 96, H = 48;
  // Near-pole (+Y) arc folded into 4 sectors around X: copies rotate to
  // +Z / -Y / -Z, far below the child's near-pole band.
  SDF::Line ln(Vector(0.25f, 1, 0).normalized(),
               Vector(-0.25f, 1, 0).normalized(), /*thickness=*/0.12f);
  SDF::AngularRepeat<SDF::Line> rep(ln, /*reps=*/4, Vector(1, 0, 0));
  int interior = expect_cull_covers_interior<W, H>(rep, "angular repeat");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies the arc-extrema cull widens phi to a Line's great-circle bulge.
 * @details The Line's two endpoints share a latitude but its great-circle arc
 *   bulges to a pole between them. Endpoint-only vertical bounds clip the polar
 *   portion of the stroke; the arc-extrema test must widen phi to the bulge.
 */
inline void test_line_arc_bulge_cull_covers_interior() {
  constexpr int W = 96, H = 48;
  // Endpoints at phi≈0.4 either side of +Y; the arc bulges through the north
  // pole (phi=0), above either endpoint's latitude.
  SDF::Line ln(Vector(0, cosf(0.4f), sinf(0.4f)),
               Vector(0, cosf(0.4f), -sinf(0.4f)), /*thickness=*/0.15f);
  int interior = expect_cull_covers_interior<W, H>(ln, "line arc bulge");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies the cull covers a Line whose antipodal endpoints select no arc.
 * @details With antipodal endpoints the arc plane is undefined and distance()
 *   measures the whole great circle, not a segment. Endpoint-derived phi bounds
 *   and a half-circle bounding cap both describe a segment, so they clip that
 *   circle to a sliver; the bounds must cover the circle the shape actually
 *   renders.
 */
inline void test_line_antipodal_cull_covers_interior() {
  constexpr int W = 96, H = 48;
  Vector a = Vector(0.4f, 0.6f, 0.69f).normalized();
  SDF::Line ln(a, -a, /*thickness=*/0.15f);
  int interior = expect_cull_covers_interior<W, H>(ln, "line antipodal");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies the cull covers a Line whose bounding cap radius exceeds pi.
 * @details A quarter arc with a stroke this wide gives half-length + thickness
 *   ≈ 3.39 rad, past the pi where cos turns back up: an unclamped cos reports a
 *   cap tighter than the whole sphere the stroke actually covers, and the
 *   horizontal cull then drops interior columns.
 */
inline void test_line_thick_cap_past_pi_cull_covers_interior() {
  constexpr int W = 96, H = 48;
  SDF::Line ln(Vector(1, 0, 0), Vector(0, 0, 1), /*thickness=*/2.6f);
  int interior = expect_cull_covers_interior<W, H>(ln, "line thick cap");
  HS_EXPECT_GT(interior, 0);
}

/**
 * @brief Verifies the Ring interval cull covers interior pixels for thin rings
 *        whose band wraps a pole while the centerline still takes the fast path.
 * @details Regression for the centerline fast path emitting two *unmerged* arcs
 *   on a pole-wrap row, leaving a one-band seam gap. The axis is tilted just off
 *   the canvas pole so r_val clears MIN_HORIZONTAL_PROJ (the fast path is
 *   eligible) yet the thin band encircles the pole, where emit_annular_band
 *   merges its two arcs into one span. Sweeps a small near-pole tilt/radius/
 *   thickness grid (both poles); expect_cull_covers_interior asserts no interior
 *   pixel is dropped.
 */
inline void test_ring_pole_wrap_cull_covers_interior() {
  // 256x128 — coarser resolutions can hide the sub-pixel gap. The triples below
  // were search-found to drop 1-4 interior pixels at the pole seam pre-fix.
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
    Basis basis_n = make_basis(Quaternion(), Vector(c.tilt, 1.0f, 0.0f));
    SDF::Ring ring_n(basis_n, c.radius, c.thickness);
    expect_cull_covers_interior<W, H>(ring_n, "ring north pole");

    Basis basis_s = make_basis(Quaternion(), Vector(c.tilt, -1.0f, 0.0f));
    SDF::Ring ring_s(basis_s, c.radius, c.thickness);
    expect_cull_covers_interior<W, H>(ring_s, "ring south pole");
  }
}

/**
 * @brief Verifies the DistortedRing cull drops no interior arc column under
 *        high-frequency centerline shifts with exact max_distortion bounds.
 * @details max_distortion widens the row/column/per-pixel reject bands; an
 *   underestimate silently culls genuine arcs. This is the pin that the caller's
 *   analytic amplitude is a true bound: sweeps spiky high-harmonic, off-grid-phase
 *   shift_fns at their exact analytic peak and asserts expect_cull_covers_interior
 *   — a full-canvas distance() reference — drops no interior column.
 */
inline void test_distorted_ring_cull_covers_interior_high_freq() {
  constexpr int W = 256, H = 128;
  // Spiky harmonic / off-grid-phase shift_fns; the bound passed is the exact
  // analytic peak.
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
  const Vector axes[] = {Vector(0, 1, 0), Vector(0.3f, 1.0f, 0.2f),
                         Vector(1, 0, 0)};
  for (const Vector &axis : axes) {
    Basis basis = make_basis(Quaternion(), axis);
    for (const Cfg &c : cfgs) {
      float amp = c.amp;
      int harmonic = c.harmonic;
      float ph = c.phase_frac;
      auto shift = [amp, harmonic, ph](float t) {
        return amp * std::sin(2.0f * PI_F * (harmonic * t + ph));
      };
      SDF::DistortedRing ring(basis, /*radius=*/0.6f, /*thickness=*/0.12f,
                              shift,
                              /*max_distortion=*/amp, /*phase=*/0.0f);
      expect_cull_covers_interior<W, H>(ring, "distorted ring");

      // The knot overload's exact distance widens the interior toward the
      // band edges; every extra pixel must still fall inside the emitted
      // intervals.
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
        0.075f + 0.1f * sinf(6.0f * PI_F * (k % ASYM_LUT_N) / ASYM_LUT_N);
  Basis basis = make_basis(Quaternion(), Vector(0.3f, 1.0f, 0.2f));
  SDF::KnotPrefilter asymmetric_pf;
  SDF::DistortedRing asymmetric_ring(basis, 0.6f, 0.08f, asymmetric, ASYM_LUT_N,
                                     0.0f, asymmetric_pf);
  expect_cull_covers_interior<W, H>(asymmetric_ring, "asymmetric ring");
}

/**
 * @brief Verifies the Face azimuth-interval cull covers every paintable pixel
 *        (including the outer AA fringe column at a silhouette edge).
 * @return The paintable-pixel count, so the caller confirms the case is non-trivial.
 * @details A pixel is paintable when its exact distance < pixel_width — the AA
 *   reach the polygon family's one-pixel cap pad is sized for. The pre-fix Face
 *   cull floor/ceil'd its azimuth intervals with no pad, so an edge falling near
 *   an integer column lost its outer AA column. Brute-forces the full canvas and
 *   asserts each paintable pixel is among those scan_region visits.
 */
template <int W, int H>
inline int expect_face_cull_covers_fringe(int sides, float rho,
                                          const Vector &axis) {
  constexpr int HV = H + hs::H_OFFSET;
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();
  Basis basis = make_basis(Quaternion(), axis);
  Vector verts3d[8];
  uint16_t idx[8];
  for (int i = 0; i < sides; ++i) {
    float a = (2.0f * PI_F * i) / sides + 0.37f;
    verts3d[i] = (basis.v * cosf(rho) +
                  (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(rho))
                     .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }
  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const Vector>(verts3d, sides),
                 std::span<const uint16_t>(idx, sides), scratch, HV, H);

  std::vector<uint8_t> visited;
  cull_visited<W, H>(face, visited);

  const float *cos_theta =
      TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  const float pixel_width = 2.0f * PI_F / W;
  int paintable = 0;
  for (int y = 0; y < H; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    for (int x = 0; x < W; ++x) {
      Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
      if (SDF::distance_of(face, p).dist < pixel_width) {
        ++paintable;
        HS_EXPECT_TRUE(visited[static_cast<size_t>(y) * W + x]);
      }
    }
  }
  return paintable;
}

/** @brief Regresses the Face AA-fringe pad over configs whose edges fall near a column. */
inline void test_face_cull_covers_aa_fringe() {
  // 256x128 — coarser resolutions can hide the dropped fringe. The triples below
  // were search-found to drop 3-6 AA-fringe pixels with the un-padded floor/ceil.
  constexpr int W = 256, H = 128;
  struct Cfg {
    int sides;
    float rho;
    Vector axis;
  };
  const Cfg cfgs[] = {
      {3, 0.15f, Vector(0, 0, 1)},
      {3, 0.30f, Vector(0, 0, 1)},
      {3, 0.48f, Vector(0, 0, 1)},
      {3, 0.18f, Vector(1, 0, 0)},
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

/**
 * @brief Rasterizes a wedge face apexed on a pole and compares it against a
 *        full-canvas distance scan.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param pole_y +1 for the north pole, -1 for the south pole.
 * @return Count of pixels the scan expects painted.
 * @details The apex projects onto a polygon vertex, so the face classifies as
 *   PoleHit::BOUNDARY: it keeps its azimuth wedge, widens the pad per row, and
 *   rasterize_face rebuilds its runs per row rather than once for the band. The
 *   expectation covers the AA fringe, whose azimuth reach is what grows as the
 *   rows approach the apex.
 */
template <int W, int H>
inline int expect_pole_vertex_face_matches_full_scan(float pole_y) {
  constexpr int HV = H + hs::H_OFFSET;
  constexpr int N_VERTS = 4;
  const float rho = 0.7f;

  Vector verts[N_VERTS];
  uint16_t idx[N_VERTS];
  verts[0] = Vector(0, pole_y, 0);
  for (int i = 1; i < N_VERTS; ++i) {
    const float a = pole_y * 0.55f * static_cast<float>(i - 2);
    verts[i] =
        Vector(sinf(rho) * cosf(a), pole_y * cosf(rho), sinf(rho) * sinf(a));
  }
  for (int i = 0; i < N_VERTS; ++i)
    idx[i] = static_cast<uint16_t>(i);

  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const Vector>(verts, N_VERTS),
                 std::span<const uint16_t>(idx, N_VERTS), scratch, HV, H);
  HS_EXPECT_TRUE(face.pole_touch);
  HS_EXPECT_TRUE(!face.full_width);

  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(fx);
    auto shader = [](const Vector &, Fragment &f) {
      f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
    };
    Scan::rasterize_face<W, H>(pipeline, canvas, face, shader);
  }
  fx.advance_display();

  const float *cos_theta =
      TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  const float pixel_width = 2.0f * PI_F / W;
  int painted = 0;
  for (int y = 0; y < H; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    for (int x = 0; x < W; ++x) {
      Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
      const float d = SDF::distance_of(face, p).dist;
      const Pixel px = fx.get_pixel(x, y);
      const bool lit = px.r != 0 || px.g != 0 || px.b != 0;
      if (d >= pixel_width) {
        HS_EXPECT_TRUE(!lit);
        continue;
      }
      // Dead band around Scan::MIN_ALPHA: an alpha that rounds the shade to
      // black is neither required nor forbidden.
      const float alpha = d <= -pixel_width
                              ? 1.0f
                              : quintic_kernel(0.5f - d / (2.0f * pixel_width));
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

// ============================================================================
// Face distance vs an independent exact oracle
//
// Face::distance computes the per-edge point-to-polygon distance in the
// gnomonic tangent plane. This pins it against an independently coded oracle
// over the face's whole gnomonic box, to float precision.
// ============================================================================

/**
 * @brief Builds one tilted N-gon face and scans its gnomonic box against an exact oracle.
 * @param sample_total Accumulates the number of non-culled samples evaluated.
 * @param sides Number of polygon points (must be <= 8).
 * @param rho Angular circumradius of the face, in radians.
 * @param axis Pole direction the face's basis is built around.
 * @param rho_inner Inner-vertex radius; > 0 interleaves star points (concave).
 * @details A convex face must match the oracle exactly inside and stay within
 * [0, oracle] outside (the half-plane path underestimates in vertex cones); a
 * concave face must match the oracle everywhere via the exact walk.
 */
inline void check_face_distance_oracle(int &sample_total, int sides, float rho,
                                       const Vector &axis,
                                       float rho_inner = 0.0f) {
  constexpr int H = 144;
  constexpr int HV = H + hs::H_OFFSET;
  HS_EXPECT_TRUE(sides <= 8);

  Basis basis = make_basis(Quaternion(), axis);
  Vector verts3d[16];
  uint16_t idx[16];
  const int n_verts = rho_inner > 0.0f ? 2 * sides : sides;
  for (int i = 0; i < n_verts; ++i) {
    float a = (2.0f * PI_F * i) / n_verts + 0.37f;
    float r = (rho_inner > 0.0f && (i & 1)) ? rho_inner : rho;
    verts3d[i] =
        (basis.v * cosf(r) + (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(r))
            .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }

  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const Vector>(verts3d, n_verts),
                 std::span<const uint16_t>(idx, n_verts), scratch, HV, H);
  HS_EXPECT_EQ(face.convex, rho_inner <= 0.0f);
  const uint32_t probe_flags = face.probe_flags();

  int samples = 0;
  // Gnomonic point normalize(center + u*px + w*py): distance() recovers (px,py)
  // exactly after dividing by dot(p,center).
  const float reach = face.max_dist * 0.98f;
  constexpr int G = 64;
  for (int gi = 0; gi <= G; ++gi) {
    for (int gj = 0; gj <= G; ++gj) {
      float px = -reach + (2.0f * reach) * gi / G;
      float py = -reach + (2.0f * reach) * gj / G;
      Vector p =
          (face.basis_v + face.basis_u * px + face.basis_w * py).normalized();

      hs::g_scan_metrics.exact_hits = 0;
      SDF::DistanceResult res = SDF::distance_of(face, p);
      SDF::DistanceResult cached_res;
      face.distance_with_flags(p, cached_res, FLT_MAX, probe_flags);
      HS_EXPECT_EQ(std::memcmp(&res, &cached_res, sizeof(res)), 0);
      if (hs::g_scan_metrics.exact_hits == 0)
        continue; // culled (outside max_dist / behind the center)

      // Independent exact oracle: point-to-polygon in the same tangent plane.
      float dmin = FLT_MAX;
      bool inside = false;
      for (int i = 0; i < face.count; ++i) {
        const auto &ep = face.packed_edges[i];
        float wx = px - ep.vx, wy = py - ep.vy;
        float t = (wx * ep.ex + wy * ep.ey) * ep.inv_len_sq;
        float cv = hs::clamp(t, 0.0f, 1.0f);
        float bx = wx - ep.ex * cv, by = wy - ep.ey * cv;
        float dsq = bx * bx + by * by;
        if (dsq < dmin)
          dmin = dsq;
        if ((ep.vy > py) != (face.poly_2d[i + 1].y > py)) {
          float isx = ep.vx + (py - ep.vy) * ep.ex * ep.inv_ej;
          if (px < isx)
            inside = !inside;
        }
      }
      float plane_exact = (inside ? -1.0f : 1.0f) * sqrtf(dmin);
      float expected =
          face.linear_dist ? plane_exact : fast_atan2(plane_exact, 1.0f);

      if (face.convex && plane_exact > 0.0f) {
        // Outside a convex face the half-plane max is a lower bound (line
        // distance, not vertex distance) that never crosses zero.
        HS_EXPECT_TRUE(res.raw_dist >= -1e-4f);
        HS_EXPECT_TRUE(res.raw_dist <= expected + 1e-4f);
      } else {
        HS_EXPECT_NEAR(res.raw_dist, expected, 1e-4f);
      }
      ++samples;
    }
  }
  sample_total += samples;
}

/**
 * @brief Verifies Face::distance reproduces the exact point-to-polygon oracle.
 * @details Drives check_face_distance_oracle across a spread of polygons
 *   (triangle, pentagon, hexagon) and tilts.
 */
inline void test_face_distance_matches_exact_oracle() {
  int samples = 0;
  check_face_distance_oracle(samples, /*sides=*/3, 0.45f,
                             Vector(0.4f, 0.3f, 1.0f));
  check_face_distance_oracle(samples, /*sides=*/5, 0.50f,
                             Vector(0.4f, 0.3f, 1.0f));
  check_face_distance_oracle(samples, /*sides=*/6, 0.40f,
                             Vector(-0.6f, 0.5f, 0.7f));
  check_face_distance_oracle(samples, /*sides=*/3, 0.12f,
                             Vector(0.4f, 0.3f, 1.0f));
  check_face_distance_oracle(samples, /*sides=*/4, 0.50f,
                             Vector(0.4f, 0.3f, 1.0f), /*rho_inner=*/0.25f);
  // Concave star: convexity detection must reject it and the exact walk must
  // reproduce the oracle everywhere.
  check_face_distance_oracle(samples, /*sides=*/6, 0.50f,
                             Vector(0.4f, 0.3f, 1.0f), /*rho_inner=*/0.25f);
  // The grid actually exercised the distance path.
  HS_EXPECT_GT(samples, 1000);
}

// ============================================================================
// Face congruence-class LUT vs the exact oracle
//
// Face::distance with a bound ClassLut serves sign-pure probes >= one cell
// diagonal from the boundary via a bilinear lookup in the canonical class
// frame; everything else falls back to the exact walk on the true edges. This
// pins the invariants of the deleted per-face check_face_lut, now across the
// canonical alignment (cyclic vertex offset, 3D rotation, mirror family):
//   1. SIGN is always correct on the LUT path (sign-purity guard).
//   2. The LUT never serves a near-boundary magnitude (>= safe_dist floor).
//   3. LUT-served values stay within the interpolation bound of the oracle.
//   4. Fallback samples still match the oracle to float precision.
// ============================================================================

/**
 * @brief Rotates a vector about a unit axis by an angle (Rodrigues).
 * @param v Vector to rotate.
 * @param k Unit rotation axis.
 * @param theta Rotation angle (radians).
 * @return The rotated vector.
 */
inline Vector rotate_about(const Vector &v, const Vector &k, float theta) {
  float c = cosf(theta), s = sinf(theta);
  return v * c + cross(k, v) * s + k * (dot(k, v) * (1.0f - c));
}

/**
 * @brief Bakes a canonical LUT from one concave star face and sweeps a
 *        transformed congruent copy against the exact oracle.
 * @param lut_total Accumulates the number of samples served by the LUT path.
 * @param cyc Cyclic shift applied to the copy's vertex order.
 * @param reflected Mirror the copy (and reverse its order, preserving winding).
 * @param rot_angle 3D rotation applied to the copy's vertices.
 */
inline void check_face_class_lut(int &lut_total, int cyc, bool reflected,
                                 float rot_angle) {
  constexpr int H = 144;
  constexpr int HV = H + hs::H_OFFSET;
  constexpr int sides = 6, n_verts = 2 * sides;
  constexpr float rho = 0.40f, rho_inner = 0.20f;
  const Vector axis(0.4f, 0.3f, 1.0f);

  Basis basis = make_basis(Quaternion(), axis);
  Vector orig[n_verts];
  for (int i = 0; i < n_verts; ++i) {
    float a = (2.0f * PI_F * i) / n_verts + 0.37f;
    float r = (i & 1) ? rho_inner : rho;
    orig[i] =
        (basis.v * cosf(r) + (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(r))
            .normalized();
  }

  // Canonical polygon: the untransformed face's own centered 2D projection.
  SDF::FaceScratchBuffer canon_scratch;
  uint16_t canon_idx[n_verts];
  for (int i = 0; i < n_verts; ++i)
    canon_idx[i] = static_cast<uint16_t>(i);
  SDF::Face canon_face(std::span<const Vector>(orig, n_verts),
                       std::span<const uint16_t>(canon_idx, n_verts),
                       canon_scratch, HV, H);
  HS_EXPECT_TRUE(!canon_face.convex);
  float canon[2 * n_verts];
  float mx = 0.0f, my = 0.0f;
  for (int i = 0; i < n_verts; ++i) {
    mx += canon_face.poly_2d[i].x;
    my += canon_face.poly_2d[i].y;
  }
  mx /= n_verts;
  my /= n_verts;
  for (int i = 0; i < n_verts; ++i) {
    canon[2 * i] = canon_face.poly_2d[i].x - mx;
    canon[2 * i + 1] = canon_face.poly_2d[i].y - my;
  }

  static int16_t lut_data[64 * 64];
  SDF::ClassLut lut;
  SDF::build_canonical_distance_lut(canon, n_verts, 64, lut_data, lut);
  // Cell diagonal of the 64x64 grid over this star's box, ~0.023. Pinned from
  // above because the sweep tolerance below is a multiple of it.
  HS_EXPECT_GT(lut.safe_dist, 0.0f);
  HS_EXPECT_LT(lut.safe_dist, 0.03f);

  // Congruent copy: rotate the sphere, then reindex. A cyclic shift by cyc
  // maps canonical vertex k to copy index (n - cyc + k) % n; the mirror family
  // reverses the order (keeping winding consistent) and binds with offset 0.
  Vector verts[n_verts];
  const Vector rot_axis = Vector(0.3f, -0.8f, 0.52f).normalized();
  for (int j = 0; j < n_verts; ++j) {
    int src = reflected ? (n_verts - j) % n_verts : (j + cyc) % n_verts;
    Vector v = rotate_about(orig[src], rot_axis, rot_angle);
    if (reflected)
      v.x = -v.x;
    verts[j] = v;
  }
  int off = reflected ? 0 : (n_verts - cyc) % n_verts;

  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const Vector>(verts, n_verts),
                 std::span<const uint16_t>(canon_idx, n_verts), scratch, HV, H);
  HS_EXPECT_TRUE(face.bind_class_lut(&lut, canon, off, reflected));
  const uint32_t probe_flags = face.probe_flags();

  // A degenerate canonical shape must be rejected by the correlation guard.
  {
    SDF::FaceScratchBuffer reject_scratch;
    static const float zeros[2 * n_verts] = {};
    SDF::Face reject_face(std::span<const Vector>(verts, n_verts),
                          std::span<const uint16_t>(canon_idx, n_verts),
                          reject_scratch, HV, H);
    HS_EXPECT_TRUE(!reject_face.bind_class_lut(&lut, zeros, 0, false));
  }

  int lut_samples = 0, sign_mismatches = 0;
  float min_lut_mag = FLT_MAX;
  const float reach = face.max_dist * 0.98f;
  constexpr int G = 96;
  for (int gi = 0; gi <= G; ++gi) {
    for (int gj = 0; gj <= G; ++gj) {
      float px = -reach + (2.0f * reach) * gi / G;
      float py = -reach + (2.0f * reach) * gj / G;
      Vector p =
          (face.basis_v + face.basis_u * px + face.basis_w * py).normalized();

      hs::g_scan_metrics.lut_hits = 0;
      hs::g_scan_metrics.exact_hits = 0;
      SDF::DistanceResult res = SDF::distance_of(face, p);
      SDF::DistanceResult cached_res;
      face.distance_with_flags(p, cached_res, FLT_MAX, probe_flags);
      HS_EXPECT_EQ(std::memcmp(&res, &cached_res, sizeof(res)), 0);
      bool took_lut = hs::g_scan_metrics.lut_hits > 0;
      if (!took_lut && hs::g_scan_metrics.exact_hits == 0)
        continue; // culled

      // Exact oracle on the copy's true edges.
      float dmin = FLT_MAX;
      bool inside = false;
      for (int i = 0; i < face.count; ++i) {
        const auto &ep = face.packed_edges[i];
        float wx = px - ep.vx, wy = py - ep.vy;
        float t = (wx * ep.ex + wy * ep.ey) * ep.inv_len_sq;
        float cv = hs::clamp(t, 0.0f, 1.0f);
        float bx = wx - ep.ex * cv, by = wy - ep.ey * cv;
        float dsq = bx * bx + by * by;
        if (dsq < dmin)
          dmin = dsq;
        if ((ep.vy > py) != (face.poly_2d[i + 1].y > py)) {
          float isx = ep.vx + (py - ep.vy) * ep.ex * ep.inv_ej;
          if (px < isx)
            inside = !inside;
        }
      }
      float plane_exact = (inside ? -1.0f : 1.0f) * sqrtf(dmin);
      float expected =
          face.linear_dist ? plane_exact : fast_atan2(plane_exact, 1.0f);

      if (took_lut) {
        ++lut_samples;
        if ((res.raw_dist < 0.0f) != (expected < 0.0f))
          ++sign_mismatches;
        float mag = std::abs(res.raw_dist);
        if (mag < min_lut_mag)
          min_lut_mag = mag;
        // Bilinear over a sign-pure cell of exact corner samples stays within
        // ~2 cell diagonals of the 1-Lipschitz field; slack covers the int16
        // quantization and the fast_atan2 conversion.
        HS_EXPECT_NEAR(res.raw_dist, expected, 3.0f * lut.safe_dist);
      } else {
        HS_EXPECT_NEAR(res.raw_dist, expected, 1e-4f);
      }
    }
  }
  HS_EXPECT_EQ(sign_mismatches, 0);
  // The sign-purity guard keeps every served magnitude a cell diagonal from
  // zero — outside the AA ramp.
  if (lut_samples > 0) {
    float floor_mag =
        face.linear_dist ? lut.safe_dist : fast_atan2(lut.safe_dist, 1.0f);
    HS_EXPECT_GT(min_lut_mag, floor_mag - 0.01f);
  }
  lut_total += lut_samples;
}

/**
 * @brief Verifies the class-LUT path across identity, cyclic-offset + 3D
 *        rotation, and mirror-family bindings.
 */
inline void test_face_class_lut_matches_oracle() {
  int lut_samples = 0;
  check_face_class_lut(lut_samples, /*cyc=*/0, /*reflected=*/false, 0.0f);
  check_face_class_lut(lut_samples, /*cyc=*/5, /*reflected=*/false, 1.1f);
  check_face_class_lut(lut_samples, /*cyc=*/0, /*reflected=*/true, 0.7f);
  // The sweep actually fired the LUT path.
  HS_EXPECT_GT(lut_samples, 300);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every sdf test case.
 * @return The module's failure count.
 */
inline int run_sdf_tests() {
  hs_test::ModuleFixture fixture("sdf");

  test_clamp_phi_in_range();
  test_clamp_phi_negative_reflects();
  test_clamp_phi_above_pi_reflects();
  test_clamp_phi_full_range();
  test_clamp_phi_band_matches_circle_extent();
  test_clamp_phi_band_pole_crossing_poses();
  test_centered_sector_angle_matches_wrap();

  test_ring_on_centerline();
  test_ring_inside_band();
  test_ring_outside_band_returns_sentinel();
  test_ring_just_outside_band();
  test_ring_small_radius_distance_symmetric();

  test_distorted_ring_constant_shift_moves_centerline();
  test_distorted_ring_sin_shift_varies_by_azimuth();
  test_distorted_ring_flat_matches_zero_knots();
  test_distorted_ring_polyline_distance_matches_bruteforce();
  test_distorted_ring_knot_extrema_tighten_band();
  test_distorted_ring_past_reach_reports_far_sentinel();

  test_polygon_at_center_inside();
  test_polygon_far_point_outside();

  test_spherical_polygon_center_inside();
  test_spherical_polygon_far_outside();
  test_spherical_polygon_center_and_edge_magnitude();
  test_spherical_polygon_sine_distance_aa_error();
  test_spherical_polygon_composes_under_csg();

  test_star_center_inside();
  test_star_far_outside();
  test_star_tip_on_boundary();

  test_flower_interior_along_petal();
  test_flower_petal_tip_on_boundary();
  test_solid_shape_unit_angle_and_no_uv_paths();

  test_inverted_fill_stays_centered();
  test_inverted_fill_scans_full_sphere();

  test_line_on_arc_is_inside();
  test_line_endpoint_is_on_line();
  test_line_perpendicular_off();
  test_line_degenerate_zero_length();
  test_line_near_coincident_endpoints_stay_point_like();

  test_torus_on_centerline_is_inside();
  test_torus_on_surface();
  test_torus_origin_is_outside_hole();
  test_torus_normal_points_outward_on_outer_rim();
  test_torus_normal_points_outward_on_top();

  test_twist_apply_displaces_y();
  test_twist_lipschitz_identity_and_closed_form();
  test_twist_bounding_inflation();
  test_warped_volume_distance_is_sphere_trace_safe();
  test_warped_volume_distance_matches_lipschitz_correction();
  test_twist_correct_normal_unit_length();

  test_union_picks_closest_shape();

  test_subtract_inside_a_outside_b_remains_inside();
  test_subtract_inside_both_becomes_outside();
  test_subtract_keeps_minuend_size_when_b_wins();
  test_subtract_solid_b_leaves_the_minuend_uncarved();
  test_subtract_star_notch_columns_survive_the_carve();
  test_subtract_empty_b_passes_a_through_verbatim();
  test_subtract_full_width_b_still_emits_the_minuend();
  test_subtract_seam_straddle_splits_the_minuend_into_frame();
  test_subtract_many_arc_seam_split_within_bound();

  test_intersection_requires_both_inside();
  test_intersection_unsorted_child_yields_sorted_result();
  test_intersection_full_width_child_replays_other();
  test_intersection_full_scan_emits_no_spans();
  test_intersection_seam_straddle_overlaps_across_wrap_frames();

  test_smooth_union_matches_union_far_from_boundary();
  test_smooth_union_blends_inside_band();
  test_smooth_union_solidity_follows_children();
  test_sentinel_clampers_are_not_blendable();
  test_csg_combinators_reject_temporary_children();

  test_union_merges_overlapping_intervals();
  test_union_seam_straddle_merges_overlapping_intervals();
  test_nested_union_emits_every_child_arc();
  test_smooth_union_seam_straddle_merges_padded_intervals();
  test_smooth_union_pad_widens_toward_pole();

  test_angular_repeat_matches_base_at_zero_angle();
  test_angular_repeat_creates_copies();
  test_angular_repeat_t_is_sector_local();

  test_warped_volume_bounding_distance_never_over_estimates();

  test_cull_covers_interior_over_orientation_grid();
  test_pole_axis_ring_bounds_skip_pole_rows();
  test_intersection_cull_covers_interior_over_polygon_pairs();
  test_subtract_cull_covers_interior_over_leaf_pairs();
  test_smooth_union_cull_covers_interior_over_leaf_pairs();
  test_smooth_union_scans_rows_past_both_children();
  test_angular_repeat_non_y_axis_cull_covers_copies();
  test_line_arc_bulge_cull_covers_interior();
  test_line_antipodal_cull_covers_interior();
  test_line_thick_cap_past_pi_cull_covers_interior();
  test_ring_pole_wrap_cull_covers_interior();
  test_distorted_ring_cull_covers_interior_high_freq();
  test_face_cull_covers_aa_fringe();
  test_star_polygon_cull_covers_aa_fringe();
  test_face_pole_vertex_matches_full_scan();
  test_face_distance_matches_exact_oracle();
  test_face_class_lut_matches_oracle();

  return fixture.result();
}

} // namespace sdf_tests
} // namespace hs_test
