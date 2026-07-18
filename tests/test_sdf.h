/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
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
#include "core/render/scan.h"
#include "core/math/geometry.h"
#include "tests/test_3dmath.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cmath>
#include <utility>
#include <vector>

namespace hs_test {
namespace sdf {

using hs_test::math3d::approx_vec;

/**
 * @brief Builds the canonical equator-facing basis: v = +Y, u = +X, w = +Z.
 * @return A Basis oriented so its pole points along +Y.
 */
inline Basis equator_basis() { return Basis{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)}; }

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

// ============================================================================
// Ring
// ============================================================================

/** @brief Verifies a point on the ring centerline reads raw_dist 0 and dist = -thickness. */
inline void test_ring_on_centerline() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  auto r = ring.distance(Vector(1, 0, 0));
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
}

/** @brief Verifies a point within the ring band reads negative dist. */
inline void test_ring_inside_band() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  float off = 0.05f;
  Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = ring.distance(p);
  HS_EXPECT_TRUE(r.dist < 0.0f);
  HS_EXPECT_TRUE(r.raw_dist <= 0.1f + 1e-3f);
}

/** @brief Verifies a point far outside the band reads the cull sentinel rather than a real dist. */
inline void test_ring_outside_band_returns_sentinel() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  auto r = ring.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist > 50.0f);
}

/** @brief Verifies just past the band edge still trips the sentinel (band edge is exclusive). */
inline void test_ring_just_outside_band() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.05f);

  // off (0.07) > thickness (0.05): point is 0.02 past the band edge.
  float off = 0.07f;
  Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = ring.distance(p);
  HS_EXPECT_TRUE(r.dist > 50.0f);
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

  SDF::DistortedRing shifted(b, 1.0f, thickness,
                             [shift](float) { return shift; },
                             /*max_distortion=*/shift, /*phase=*/0.0f);
  auto rs = shifted.distance(p);
  HS_EXPECT_TRUE(rs.dist < 50.0f);
  HS_EXPECT_NEAR(rs.raw_dist, 0.0f, 1e-2f);
  HS_EXPECT_NEAR(rs.dist, -thickness, 1e-2f);
  HS_EXPECT_NEAR(rs.t, 0.0f, 1e-2f);

  // Same point, no shift: the centerline stays at π/2, so it now sits `shift`
  // radians off (raw_dist ≈ shift) — the shift moved the centerline.
  SDF::DistortedRing plain(b, 1.0f, thickness, [](float) { return 0.0f; },
                           /*max_distortion=*/shift, /*phase=*/0.0f);
  auto rp = plain.distance(p);
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
  auto r_on = ring.distance(on);
  HS_EXPECT_NEAR(r_on.t, 0.25f, 1e-2f);
  HS_EXPECT_NEAR(r_on.raw_dist, 0.0f, 1e-2f);

  // Same azimuth on the unshifted equator (+Z): centerline moved by amp here.
  auto r_off = ring.distance(Vector(0, 0, 1));
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
    SDF::DistortedRing polyline(basis, radius, thickness, knots, LUT_N, 0.0f,
                                0.0f);
    const float target = radius * (PI_F / 2.0f);
    const float azimuths[] = {0.0f, 1e-5f, PI_F / 2.0f, PI_F,
                              2.0f * PI_F - 1e-5f};
    const float offsets[] = {-0.04f, 0.0f, 0.04f};
    for (float azimuth : azimuths) {
      for (float offset : offsets) {
        float polar = hs::clamp(target + offset, 0.0f, PI_F);
        Vector p = basis.v * cosf(polar) +
                   (basis.u * cosf(azimuth) + basis.w * sinf(azimuth)) *
                       sinf(polar);
        auto actual = flat.distance(p);
        auto expected = polyline.distance(p);
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
template <int LUT_N>
inline void expect_polyline_distance_matches_bruteforce() {
  Basis b = make_basis(Quaternion(), Vector(0.3f, 1.0f, 0.2f));
  const float amp = 0.2f;
  const int harmonic = 5;
  const float radius = 0.5f; // target_angle = π/4: curved chart, tilted axis
  const float thickness = 0.06f;
  float knots[LUT_N + 1];
  for (int k = 0; k <= LUT_N; ++k)
    knots[k] = amp * std::sin(2.0f * PI_F * harmonic * (k % LUT_N) / LUT_N);
  SDF::DistortedRing ring(b, radius, thickness, knots, LUT_N, amp, 0.0f);

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
  const float probes[][2] = {{0.05f, 0.04f},  {0.05f, -0.05f}, {0.1f, 0.05f},
                             {0.1f, -0.04f},  {0.998f, 0.04f}, {0.25f, 0.05f},
                             {0.375f, -0.05f}, {0.6f, 0.03f}};
  for (const auto &pr : probes) {
    Vector p = on_sphere(pr[0], pr[1]);
    auto r = ring.distance(p);
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
  SDF::DistortedRing ring(basis, RADIUS, THICKNESS, knots, LUT_N, 0.8f, 0.0f);

  HS_EXPECT_NEAR(ring.max_distortion, 0.18f, 1e-6f);
  HS_EXPECT_NEAR(ring.max_thickness, 0.23f, 1e-6f);
  for (int k : {0, 4, 8}) {
    float t = static_cast<float>(k % LUT_N) / LUT_N;
    float azimuth = 2.0f * PI_F * t;
    float polar = TARGET + knots[k];
    Vector p(cosf(azimuth) * sinf(polar), cosf(polar),
             sinf(azimuth) * sinf(polar));
    HS_EXPECT_NEAR(ring.distance(p).raw_dist, 0.0f, 2e-4f);
  }

  float below = TARGET - 0.03f - THICKNESS - 1e-3f;
  float above = TARGET + 0.18f + THICKNESS + 1e-3f;
  HS_EXPECT_GT(ring.distance(Vector(sinf(below), cosf(below), 0.0f)).dist,
               50.0f);
  HS_EXPECT_GT(ring.distance(Vector(sinf(above), cosf(above), 0.0f)).dist,
               50.0f);
}

// ============================================================================
// PlanarPolygon  (Basis at top of sphere; distance to nearest edge)
// ============================================================================

/** @brief Verifies the polygon center is inside, with dist equal to the negated apothem. */
inline void test_polygon_at_center_inside() {
  Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, /*thickness*/ 0.5f, /*sides*/ 6, /*phase*/ 0.0f);

  auto r = poly.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
  float apothem = 0.5f * std::cos(PI_F / 6.0f);
  HS_EXPECT_NEAR(r.dist, -apothem, 1e-3f);
}

/** @brief Verifies the antipode of the polygon center is outside (positive dist). */
inline void test_polygon_far_point_outside() {
  Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, 0.3f, 6, 0.0f);

  auto r = poly.distance(Vector(0, -1, 0));
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

// ============================================================================
// SphericalPolygon — great-circle edges
// ============================================================================

/** @brief Verifies the spherical-polygon center is strictly inside. */
inline void test_spherical_polygon_center_inside() {
  Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, /*radius*/ 0.5f, /*sides*/ 5, /*phase*/ 0.0f);
  auto r = sp.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies the antipode of the spherical-polygon center is outside. */
inline void test_spherical_polygon_far_outside() {
  Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, 0.3f, 6, 0.0f);
  auto r = sp.distance(Vector(0, -1, 0));
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

  auto center = sp.distance(Vector(0, 1, 0));
  HS_EXPECT_NEAR(center.dist, -inradius, 1e-2f);
  HS_EXPECT_NEAR(center.raw_dist, 0.0f, 1e-3f);

  // Edge midpoint: polar angle = inradius along +u (the sector bisector).
  Vector edge_mid(std::sin(inradius), std::cos(inradius), 0.0f);
  auto em = sp.distance(edge_mid);
  HS_EXPECT_NEAR(em.dist, 0.0f, 1e-2f);
  HS_EXPECT_NEAR(em.raw_dist, inradius, 1e-2f);
}

// ============================================================================
// Star
// ============================================================================

/** @brief Verifies the star center is interior (negative dist). */
inline void test_star_center_inside() {
  Basis b = equator_basis();
  SDF::Star star(b, /*radius*/ 0.6f, /*sides*/ 5, /*phase*/ 0.0f);
  auto r = star.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies the antipode of the star center is outside. */
inline void test_star_far_outside() {
  Basis b = equator_basis();
  SDF::Star star(b, 0.4f, 5, 0.0f);
  auto r = star.distance(Vector(0, -1, 0));
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
  auto r = star.distance(tip);
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
  auto r = flower.distance(p);
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
  auto r = flower.distance(tip);
  HS_EXPECT_NEAR(r.dist, 0.0f, 1e-2f);
  HS_EXPECT_NEAR(r.raw_dist, outer, 1e-2f);
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
  HS_EXPECT_TRUE(sp.distance(center).dist < 0.0f);
  HS_EXPECT_TRUE(sp.distance(far_side).dist > 0.0f);

  SDF::Star star(fb, fr, 5, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(star.distance(center).dist < 0.0f);
  HS_EXPECT_TRUE(star.distance(far_side).dist > 0.0f);

  SDF::PlanarPolygon pp(fb, fr * (PI_F / 2.0f), 6, 0.0f, /*invert=*/true);
  HS_EXPECT_TRUE(pp.distance(center).dist < 0.0f);
  HS_EXPECT_TRUE(pp.distance(far_side).dist > 0.0f);

  // Flower fills around the antipode of its axis, so the sides swap; sample
  // off the exact poles (azimuth is degenerate there).
  SDF::Flower fl(fb, fr, 6, 0.0f, /*invert=*/true);
  Vector near_far(std::sin(0.1f), -std::cos(0.1f), 0.0f);
  Vector near_center(std::sin(0.1f), std::cos(0.1f), 0.0f);
  HS_EXPECT_TRUE(fl.distance(near_far).dist < 0.0f);
  HS_EXPECT_TRUE(fl.distance(near_center).dist > 0.0f);
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
  bool handled =
      sp.get_horizontal_intervals<288, 144>(72, [](float, float) {});
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
  auto r = ln.distance(mid);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-2f);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-2f);
}

/** @brief Verifies an endpoint counts as on the line (raw_dist 0, dist = -thickness). */
inline void test_line_endpoint_is_on_line() {
  Vector a(1, 0, 0);
  Vector b(0, 0, 1);
  SDF::Line ln(a, b, 0.1f);
  auto r = ln.distance(a);
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
  auto r = ln.distance(p);
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
  auto r = ln.distance(a);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);

  auto r2 = ln.distance(Vector(0, 1, 0));
  HS_EXPECT_NEAR(r2.raw_dist, PI_F * 0.5f, 1e-2f);
  HS_EXPECT_TRUE(r2.dist > 0.0f);
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
  HS_EXPECT_NEAR(flat.lipschitz(Vector(2, 0, 0), flat.make_ctx(Vector(2, 0, 0))),
                 1.0f, 1e-6f);

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
  uint32_t rs = 0x5eed1337u;
  const auto next = [&rs]() {
    rs = rs * 1664525u + 1013904223u;
    return static_cast<double>(rs >> 8) / 16777216.0;
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
        p = Vector(static_cast<float>((next() * 2 - 1) * reach),
                   static_cast<float>((next() * 2 - 1) * reach),
                   static_cast<float>((next() * 2 - 1) * reach));
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
  auto r = u.distance(mid_a);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-2f);

  Vector mid_b = ((Vector(-1, 0, 0) + Vector(0, 0, -1)) * 0.5f).normalized();
  auto r2 = u.distance(mid_b);
  HS_EXPECT_NEAR(r2.dist, -0.1f, 1e-2f);
}

/** @brief Verifies Union exposes the max of its children's thicknesses. */
inline void test_union_thickness_is_max() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.3f);
  SDF::Union<SDF::Line, SDF::Line> u(la, lb);
  HS_EXPECT_NEAR(u.thickness, 0.3f, 1e-6f);
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
  auto r = s.distance(mid_a);
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies a point inside both A and B becomes outside the difference A - B. */
inline void test_subtract_inside_both_becomes_outside() {
  // Same line for A and B → A - A is empty everywhere.
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, la);
  Vector mid = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = s.distance(mid);
  // max(dist(A), -dist(B)) = max(-0.1, 0.1) = 0.1 → outside.
  HS_EXPECT_NEAR(r.dist, 0.1f, 1e-3f);
}

namespace sdf_subtract_detail {
/**
 * @brief Mock SDF shape that emits a fixed (possibly unsorted, multi-) interval list.
 * @details Exercises Subtract's scanline set-difference independently of any real
 *   shape. Minimal surface: only thickness, is_solid, and get_horizontal_intervals
 *   are touched by Subtract's ctor + interval path.
 */
struct MockIntervalShape {
  const std::vector<std::pair<float, float>> *ivs; /**< Interval list this mock replays. */
  float thickness = 0.1f;                          /**< Stroke half-width the parent reads. */
  static constexpr bool is_solid = true;           /**< Marks the mock as a solid fill shape. */
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
};

/**
 * @brief Mock that falls back to a full-row scan.
 * @details Returning false means "I cannot produce intervals — caller must scan
 *   the whole row with distance()", NOT that the shape covers the row.
 */
struct MockFullWidthShape {
  float thickness = 0.1f;                /**< Stroke half-width the parent reads. */
  static constexpr bool is_solid = true; /**< Marks the mock as a solid fill shape. */
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
} // namespace sdf_subtract_detail

/**
 * @brief Verifies UNSORTED, multi-piece B intervals still yield the correct set difference.
 * @details The result must be emitted in start-sorted order; scan_region's
 *   coalescer silently drops any out-of-order interval.
 */
inline void test_subtract_unsorted_b_yields_sorted_set_difference() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{0.0f, 100.0f}};
  std::vector<P> b_ivs = {{60.0f, 70.0f}, {20.0f, 30.0f}}; // unsorted, multi
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  // [0,100] - {[20,30],[60,70]} = [0,20],[30,60],[70,100].
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(3));
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 20.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 30.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 60.0f, 1e-4f);
  HS_EXPECT_NEAR(out[2].first, 70.0f, 1e-4f);
  HS_EXPECT_NEAR(out[2].second, 100.0f, 1e-4f);
  for (size_t i = 1; i < out.size(); ++i)
    HS_EXPECT_TRUE(out[i].first >= out[i - 1].first);
}

/**
 * @brief Verifies that when B removes nothing, unsorted A intervals pass through start-sorted.
 * @details Sorted passthrough keeps the coalescer from dropping the earlier interval.
 */
inline void test_subtract_unsorted_a_passthrough_is_sorted() {
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
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 50.0f, 1e-4f);
}

/**
 * @brief Verifies that when B cannot produce intervals, Subtract requests a full-row scan.
 * @details B returning false means Subtract cannot compute the set difference, so
 *   it must return false — like Union/SmoothUnion — letting scan_region evaluate
 *   distance()=max(A,-B) per pixel. Returning true while emitting nothing would
 *   make scan_region SKIP the row and silently erase all of A.
 */
inline void test_subtract_full_width_b_requests_full_row_scan() {
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
  HS_EXPECT_TRUE(!ok);
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
}

/**
 * @brief Verifies a B band straddling θ=0 in a different wrap frame still carves A.
 * @details A emits the seam band as [-10, 10]; B emits the SAME physical band as
 *   [W-10, W+10]. A raw-coordinate disjointness test sees no overlap and leaves A
 *   uncarved (under-carve at the seam). After normalizing both into [0, W) the
 *   bands coincide and the difference is empty.
 */
inline void test_subtract_seam_straddle_carves_across_wrap_frames() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-10.0f, 10.0f}};   // seam band, negative frame
  std::vector<P> b_ivs = {{246.0f, 266.0f}};  // same band, [W,2W) frame (W=256)
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  // Both normalize to {[246,256],[0,10]}; A - B is empty.
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
}

/**
 * @brief Verifies a many-arc seam-straddling Subtract stays within the span-count bound.
 * @details normalize_intervals_to_range seam-splits each child into [0, W); a span
 *   crossing θ=0 becomes two, so the norm buffers are sized 2x. push_interval traps
 *   (fail-fast) if the post-split count exceeds that cap. Drive A with many disjoint
 *   arcs — several straddling the seam in a wrapped frame so they split — and a few
 *   carving B arcs, then assert the result is valid (in [0, W), start-sorted) and
 *   that no trap fired (reaching this line at all).
 */
inline void test_subtract_many_arc_seam_split_within_bound() {
  using P = std::pair<float, float>;
  using Mock = sdf_subtract_detail::MockIntervalShape;
  constexpr float W = 256.0f;

  // Twelve disjoint A arcs; the first two straddle the seam (negative / over-W
  // frame) so each splits into two on normalization (14 norm spans, well under 64).
  std::vector<P> a_ivs = {
      {-6.0f, 4.0f},   {250.0f, 262.0f},                          // seam straddlers
      {20.0f, 30.0f},  {40.0f, 50.0f},   {60.0f, 70.0f},
      {80.0f, 90.0f},  {100.0f, 110.0f}, {120.0f, 130.0f},
      {140.0f, 150.0f},{160.0f, 170.0f}, {180.0f, 190.0f},
      {200.0f, 210.0f}};
  std::vector<P> b_ivs = {{45.0f, 55.0f}, {125.0f, 135.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  // Every emitted span is a valid in-frame arc, start-sorted (the coalescer
  // requirement). Reaching here means push_interval never trapped.
  HS_EXPECT_TRUE(!out.empty());
  for (size_t i = 0; i < out.size(); ++i) {
    HS_EXPECT_TRUE(out[i].first >= 0.0f && out[i].second <= W);
    HS_EXPECT_TRUE(out[i].first < out[i].second);
    if (i > 0)
      HS_EXPECT_TRUE(out[i].first >= out[i - 1].first);
  }
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
  auto r = inter.distance(a);
  HS_EXPECT_TRUE(r.dist < 0.0f);

  Vector far_pt(-1, 0, 0);
  auto r2 = inter.distance(far_pt);
  HS_EXPECT_TRUE(r2.dist > 0.0f);
}

/** @brief Verifies Intersection exposes the min of its children's thicknesses. */
inline void test_intersection_thickness_is_min() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.3f);
  SDF::Intersection<SDF::Line, SDF::Line> inter(la, lb);
  HS_EXPECT_NEAR(inter.thickness, 0.1f, 1e-6f);
}

/**
 * @brief Verifies an UNSORTED multi-interval child still yields a start-sorted intersection.
 * @details Intersection's merge-sweep assumes both child interval lists are
 *   start-sorted; out-of-order output would be dropped by scan_region's coalescer.
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
  std::vector<P> ivs = {{60.0f, 80.0f}, {20.0f, 40.0f}}; // multi, emission order
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
  auto r_hard = u.distance(mid_a);
  auto r_soft = su.distance(mid_a);
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
  float dA = la.distance(p).dist;
  float dB = lb.distance(p).dist;
  HS_EXPECT_TRUE(std::abs(dA - dB) < k);

  float h = std::max(k - std::abs(dA - dB), 0.0f) / k;
  float m = h * h * h * k * (1.0f / 6.0f);
  float soft = su.distance(p).dist;
  HS_EXPECT_NEAR(soft, std::min(dA, dB) - m, 1e-4f);
  HS_EXPECT_LT(soft, std::min(dA, dB) - 1e-4f);
  HS_EXPECT_NEAR(std::min(dA, dB) - soft, m, 1e-4f);

  // Outside the band the blend vanishes and collapses to the hard min.
  Vector q = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  float qA = la.distance(q).dist;
  float qB = lb.distance(q).dist;
  HS_EXPECT_TRUE(std::abs(qA - qB) >= k);
  HS_EXPECT_NEAR(su.distance(q).dist, std::min(qA, qB), 1e-5f);
}

/**
 * @brief Verifies SmoothUnion solidity follows its children.
 * @details Children must share solidity (enforced by static_assert): two strokes
 *   take the soft falloff path, two solids the hard 1-px silhouette path.
 */
inline void test_smooth_union_solidity_follows_children() {
  static_assert(!SDF::SmoothUnion<SDF::Ring, SDF::Ring>::is_solid,
                "two strokes -> falloff path");
  static_assert(
      SDF::SmoothUnion<SDF::PlanarPolygon, SDF::PlanarPolygon>::is_solid,
      "two solids -> silhouette path");
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
  init_geometry_luts<W, H>(); // fill sin_phi; scan_region does this in production
  const int row = H / 2; // equatorial row: sinφ ≈ 1, pad ≈ k·W/(2π)
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
  init_geometry_luts<W, H>(); // fill sin_phi; scan_region does this in production
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
  auto r_base = ln.distance(mid);
  auto r_rep = rep.distance(mid);
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

  auto r = rep.distance(mid_rot);
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
  float t_global_1 = ring.distance(p).t;
  float t_global_2 = ring.distance(p2).t;
  float dg = t_global_2 - t_global_1;
  dg -= floorf(dg);
  HS_EXPECT_NEAR(std::min(dg, 1.0f - dg), 1.0f / reps, 1e-3f);

  // The repeated shape folds both into the same sector → identical sector-local t.
  float t_rep_1 = rep.distance(p).t;
  float t_rep_2 = rep.distance(p2).t;
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
      [&](int wx, int y, const Vector &) {
        if (wx >= 0 && wx < W && y >= 0 && y < H)
          visited[static_cast<size_t>(y) * W + wx] = 1;
      });
}

/**
 * @brief Asserts no pixel clearly inside the shape is dropped by the cull.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Shape SDF shape type providing the cull and distance interface.
 * @param shape Shape under test.
 * @return Count of interior pixels (dist < -pixel_width) found.
 * @details Interior pixels are found by a brute-force full-canvas exact distance
 *   scan. Asserts the case is non-trivial (at least one interior pixel) so no
 *   shape silently contributes zero coverage.
 */
template <int W, int H, typename Shape>
inline int expect_cull_covers_interior(const Shape &shape) {
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();
  std::vector<uint8_t> visited;
  cull_visited<W, H>(shape, visited);

  const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  const float pixel_width = 2.0f * PI_F / W;
  int interior = 0;
  for (int y = 0; y < H; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    for (int x = 0; x < W; ++x) {
      Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
      if (shape.distance(p).dist < -pixel_width) {
        ++interior;
        HS_EXPECT_TRUE(visited[static_cast<size_t>(y) * W + x]);
      }
    }
  }
  HS_EXPECT_GT(interior, 0);
  return interior;
}

/** @brief Verifies the interval cull covers every interior pixel across an orientation/radius grid. */
inline void test_cull_covers_interior_over_orientation_grid() {
  constexpr int W = 96, H = 48;

  // Poles, equator, and oblique tilts.
  const Vector axes[] = {
      Vector(0, 1, 0),         Vector(0, -1, 0),        Vector(1, 0, 0),
      Vector(0, 0, 1),         Vector(1, 1, 0.4f),      Vector(-0.5f, 0.7f, -0.6f),
      Vector(0.3f, -0.8f, 0.5f)};

  for (const Vector &axis : axes) {
    Basis basis = make_basis(Quaternion(), axis);
    for (float radius : {0.3f, 0.6f, 0.9f}) {
      SDF::Ring ring(basis, radius, /*thickness=*/0.25f);
      expect_cull_covers_interior<W, H>(ring);

      SDF::SphericalPolygon spoly(basis, radius, /*sides=*/5, 0.0f);
      expect_cull_covers_interior<W, H>(spoly);

      SDF::Star star(basis, radius, /*sides=*/5, 0.0f);
      expect_cull_covers_interior<W, H>(star);

      // PlanarPolygon's `thickness` is its angular circumradius.
      SDF::PlanarPolygon ppoly(basis, /*thickness=*/radius, /*sides=*/6, 0.0f);
      expect_cull_covers_interior<W, H>(ppoly);
    }
  }
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
  int interior = expect_cull_covers_interior<W, H>(rep);
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
  int interior = expect_cull_covers_interior<W, H>(ln);
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
  struct Cfg { float tilt, radius, thickness; };
  const Cfg cfgs[] = {
      {0.15f, 0.22f, 0.13f}, {0.16f, 0.25f, 0.14f},
      {0.16f, 0.33f, 0.15f}, {0.16f, 0.39f, 0.15f},
  };
  for (const Cfg &c : cfgs) {
    Basis basis_n = make_basis(Quaternion(), Vector(c.tilt, 1.0f, 0.0f));
    SDF::Ring ring_n(basis_n, c.radius, c.thickness);
    expect_cull_covers_interior<W, H>(ring_n);

    Basis basis_s = make_basis(Quaternion(), Vector(c.tilt, -1.0f, 0.0f));
    SDF::Ring ring_s(basis_s, c.radius, c.thickness);
    expect_cull_covers_interior<W, H>(ring_s);
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
  struct Cfg { float amp; int harmonic; float phase_frac; };
  const Cfg cfgs[] = {
      {0.18f, 127, 0.5f / 256.0f}, {0.20f, 255, 0.5f / 256.0f},
      {0.15f, 384, 0.25f / 256.0f}, {0.22f, 200, 0.5f / 256.0f},
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
      SDF::DistortedRing ring(basis, /*radius=*/0.6f, /*thickness=*/0.12f, shift,
                              /*max_distortion=*/amp, /*phase=*/0.0f);
      expect_cull_covers_interior<W, H>(ring);

      // The knot overload's exact distance widens the interior toward the
      // band edges; every extra pixel must still fall inside the emitted
      // intervals.
      constexpr int LUT_N = 256;
      float knots[LUT_N + 1];
      for (int k = 0; k <= LUT_N; ++k)
        knots[k] = shift(static_cast<float>(k % LUT_N) / LUT_N);
      SDF::DistortedRing poly(basis, /*radius=*/0.6f, /*thickness=*/0.12f,
                              knots, LUT_N, /*max_distortion=*/amp,
                              /*phase=*/0.0f);
      expect_cull_covers_interior<W, H>(poly);
    }
  }

  constexpr int ASYM_LUT_N = 64;
  float asymmetric[ASYM_LUT_N + 1];
  for (int k = 0; k <= ASYM_LUT_N; ++k)
    asymmetric[k] = 0.075f + 0.1f * sinf(6.0f * PI_F * (k % ASYM_LUT_N) /
                                        ASYM_LUT_N);
  Basis basis = make_basis(Quaternion(), Vector(0.3f, 1.0f, 0.2f));
  SDF::DistortedRing asymmetric_ring(basis, 0.6f, 0.08f, asymmetric,
                                     ASYM_LUT_N, 0.8f, 0.0f);
  expect_cull_covers_interior<W, H>(asymmetric_ring);
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
                 std::span<const uint16_t>(idx, sides), /*thickness=*/0.0f,
                 scratch, HV, H);

  std::vector<uint8_t> visited;
  cull_visited<W, H>(face, visited);

  const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  const float pixel_width = 2.0f * PI_F / W;
  int paintable = 0;
  for (int y = 0; y < H; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    for (int x = 0; x < W; ++x) {
      Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
      if (face.distance(p).dist < pixel_width) {
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
  struct Cfg { int sides; float rho; Vector axis; };
  const Cfg cfgs[] = {
      {3, 0.15f, Vector(0, 0, 1)}, {3, 0.30f, Vector(0, 0, 1)},
      {3, 0.48f, Vector(0, 0, 1)}, {3, 0.18f, Vector(1, 0, 0)},
  };
  int total_paintable = 0;
  for (const Cfg &c : cfgs)
    total_paintable += expect_face_cull_covers_fringe<W, H>(c.sides, c.rho, c.axis);
  HS_EXPECT_GT(total_paintable, 1000);
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
 * @param sides Number of polygon points (must be <= 8).
 * @param rho Angular circumradius of the face, in radians.
 * @param axis Pole direction the face's basis is built around.
 * @param rho_inner Inner-vertex radius; > 0 interleaves star points (concave).
 * @return The number of non-culled samples evaluated.
 * @details A convex face must match the oracle exactly inside and stay within
 * [0, oracle] outside (the half-plane path underestimates in vertex cones); a
 * concave face must match the oracle everywhere via the exact walk.
 */
inline int check_face_distance_oracle(int sides, float rho, const Vector &axis,
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
    verts3d[i] = (basis.v * cosf(r) +
                  (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(r))
                     .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }

  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const Vector>(verts3d, n_verts),
                 std::span<const uint16_t>(idx, n_verts), /*thickness=*/0.0f,
                 scratch, HV, H);
  HS_EXPECT_EQ(face.convex, rho_inner <= 0.0f);

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
      SDF::DistanceResult res = face.distance(p);
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
        if ((ep.vy > py) != (ep.next_vy > py)) {
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
  return samples;
}

/**
 * @brief Verifies Face::distance reproduces the exact point-to-polygon oracle.
 * @details Drives check_face_distance_oracle across a spread of polygons
 *   (triangle, pentagon, hexagon) and tilts.
 */
inline void test_face_distance_matches_exact_oracle() {
  int samples = 0;
  samples += check_face_distance_oracle(/*sides=*/3, 0.45f, Vector(0.4f, 0.3f, 1.0f));
  samples += check_face_distance_oracle(/*sides=*/5, 0.50f, Vector(0.4f, 0.3f, 1.0f));
  samples += check_face_distance_oracle(/*sides=*/6, 0.40f, Vector(-0.6f, 0.5f, 0.7f));
  // Concave star: convexity detection must reject it and the exact walk must
  // reproduce the oracle everywhere.
  samples += check_face_distance_oracle(/*sides=*/6, 0.50f, Vector(0.4f, 0.3f, 1.0f),
                                        /*rho_inner=*/0.25f);
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
 * @param cyc Cyclic shift applied to the copy's vertex order.
 * @param reflected Mirror the copy (and reverse its order, preserving winding).
 * @param rot_angle 3D rotation applied to the copy's vertices.
 * @return The number of samples served by the LUT path.
 */
inline int check_face_class_lut(int cyc, bool reflected, float rot_angle) {
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
    orig[i] = (basis.v * cosf(r) +
               (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(r))
                  .normalized();
  }

  // Canonical polygon: the untransformed face's own centered 2D projection.
  SDF::FaceScratchBuffer canon_scratch;
  uint16_t canon_idx[n_verts];
  for (int i = 0; i < n_verts; ++i)
    canon_idx[i] = static_cast<uint16_t>(i);
  SDF::Face canon_face(std::span<const Vector>(orig, n_verts),
                       std::span<const uint16_t>(canon_idx, n_verts), 0.0f,
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
  HS_EXPECT_GT(lut.safe_dist, 0.0f);

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
                 std::span<const uint16_t>(canon_idx, n_verts), 0.0f, scratch,
                 HV, H);
  HS_EXPECT_TRUE(face.bind_class_lut(&lut, canon, off, reflected));

  // A degenerate canonical shape must be rejected by the correlation guard.
  {
    SDF::FaceScratchBuffer reject_scratch;
    static const float zeros[2 * n_verts] = {};
    SDF::Face reject_face(std::span<const Vector>(verts, n_verts),
                          std::span<const uint16_t>(canon_idx, n_verts), 0.0f,
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
      SDF::DistanceResult res = face.distance(p);
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
        if ((ep.vy > py) != (ep.next_vy > py)) {
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
    float floor_mag = face.linear_dist ? lut.safe_dist
                                       : fast_atan2(lut.safe_dist, 1.0f);
    HS_EXPECT_GT(min_lut_mag, floor_mag - 0.01f);
  }
  return lut_samples;
}

/**
 * @brief Verifies the class-LUT path across identity, cyclic-offset + 3D
 *        rotation, and mirror-family bindings.
 */
inline void test_face_class_lut_matches_oracle() {
  int lut_samples = 0;
  lut_samples += check_face_class_lut(/*cyc=*/0, /*reflected=*/false, 0.0f);
  lut_samples += check_face_class_lut(/*cyc=*/5, /*reflected=*/false, 1.1f);
  lut_samples += check_face_class_lut(/*cyc=*/0, /*reflected=*/true, 0.7f);
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

  test_ring_on_centerline();
  test_ring_inside_band();
  test_ring_outside_band_returns_sentinel();
  test_ring_just_outside_band();

  test_distorted_ring_constant_shift_moves_centerline();
  test_distorted_ring_sin_shift_varies_by_azimuth();
  test_distorted_ring_flat_matches_zero_knots();
  test_distorted_ring_polyline_distance_matches_bruteforce();
  test_distorted_ring_knot_extrema_tighten_band();

  test_polygon_at_center_inside();
  test_polygon_far_point_outside();

  test_spherical_polygon_center_inside();
  test_spherical_polygon_far_outside();
  test_spherical_polygon_center_and_edge_magnitude();

  test_star_center_inside();
  test_star_far_outside();
  test_star_tip_on_boundary();

  test_flower_interior_along_petal();
  test_flower_petal_tip_on_boundary();

  test_inverted_fill_stays_centered();
  test_inverted_fill_scans_full_sphere();

  test_line_on_arc_is_inside();
  test_line_endpoint_is_on_line();
  test_line_perpendicular_off();
  test_line_degenerate_zero_length();

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
  test_union_thickness_is_max();

  test_subtract_inside_a_outside_b_remains_inside();
  test_subtract_inside_both_becomes_outside();
  test_subtract_unsorted_b_yields_sorted_set_difference();
  test_subtract_unsorted_a_passthrough_is_sorted();
  test_subtract_full_width_b_requests_full_row_scan();
  test_subtract_seam_straddle_carves_across_wrap_frames();
  test_subtract_many_arc_seam_split_within_bound();

  test_intersection_requires_both_inside();
  test_intersection_thickness_is_min();
  test_intersection_unsorted_child_yields_sorted_result();
  test_intersection_full_width_child_replays_other();
  test_intersection_seam_straddle_overlaps_across_wrap_frames();

  test_smooth_union_matches_union_far_from_boundary();
  test_smooth_union_blends_inside_band();
  test_smooth_union_solidity_follows_children();

  test_union_merges_overlapping_intervals();
  test_union_seam_straddle_merges_overlapping_intervals();
  test_smooth_union_seam_straddle_merges_padded_intervals();
  test_smooth_union_pad_widens_toward_pole();

  test_angular_repeat_matches_base_at_zero_angle();
  test_angular_repeat_creates_copies();
  test_angular_repeat_t_is_sector_local();

  test_warped_volume_bounding_distance_never_over_estimates();

  test_cull_covers_interior_over_orientation_grid();
  test_angular_repeat_non_y_axis_cull_covers_copies();
  test_line_arc_bulge_cull_covers_interior();
  test_ring_pole_wrap_cull_covers_interior();
  test_distorted_ring_cull_covers_interior_high_freq();
  test_face_cull_covers_aa_fringe();
  test_face_distance_matches_exact_oracle();
  test_face_class_lut_matches_oracle();

  return fixture.result();
}

} // namespace sdf
} // namespace hs_test

