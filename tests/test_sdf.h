/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/sdf.h.
 *
 * Coverage:
 *   - clamp_phi utility
 *   - Spherical SDF primitives (Ring, PlanarPolygon, SphericalPolygon, Star, Line)
 *   - 3D Torus
 *   - CSG operators (Union, SmoothUnion, Subtract, Intersection)
 *   - AngularRepeat
 *
 * Tests focus on the distance() interface — the rendering pipeline
 * (get_vertical_bounds<H> / get_horizontal_intervals<W,H>) is exercised
 * indirectly by feeding known points and verifying signed distance.
 */
#pragma once

#include "core/sdf.h"
#include "core/scan.h"
#include "core/geometry.h"
#include "tests/test_3dmath.h"
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
 * @brief Verifies SmoothUnion derives its solidity from its children.
 * @details Like Union/Subtract/Intersection, a stroke child keeps its soft falloff
 *   instead of collapsing to a hard 1-px silhouette in process_pixel.
 */
inline void test_smooth_union_solidity_follows_children() {
  static_assert(!SDF::SmoothUnion<SDF::Ring, SDF::Ring>::is_solid,
                "two strokes -> falloff path");
  static_assert(
      SDF::SmoothUnion<SDF::PlanarPolygon, SDF::PlanarPolygon>::is_solid,
      "two solids -> silhouette path");
  static_assert(!SDF::SmoothUnion<SDF::PlanarPolygon, SDF::Ring>::is_solid,
                "mixed -> falloff path");
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
 * @return Count of interior pixels (dist < -pixel_width) found, so the caller can
 *   confirm the case was non-trivial.
 * @details Interior pixels are found by a brute-force full-canvas exact distance scan.
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

  int total_interior = 0;
  for (const Vector &axis : axes) {
    Basis basis = make_basis(Quaternion(), axis);
    for (float radius : {0.3f, 0.6f, 0.9f}) {
      SDF::Ring ring(basis, radius, /*thickness=*/0.25f);
      total_interior += expect_cull_covers_interior<W, H>(ring);

      SDF::SphericalPolygon spoly(basis, radius, /*sides=*/5, 0.0f);
      total_interior += expect_cull_covers_interior<W, H>(spoly);

      SDF::Star star(basis, radius, /*sides=*/5, 0.0f);
      total_interior += expect_cull_covers_interior<W, H>(star);

      // PlanarPolygon's `thickness` is its angular circumradius.
      SDF::PlanarPolygon ppoly(basis, /*thickness=*/radius, /*sides=*/6, 0.0f);
      total_interior += expect_cull_covers_interior<W, H>(ppoly);
    }
  }
  HS_EXPECT_GT(total_interior, 1000);
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
  int total_interior = 0;
  for (const Cfg &c : cfgs) {
    Basis basis_n = make_basis(Quaternion(), Vector(c.tilt, 1.0f, 0.0f));
    SDF::Ring ring_n(basis_n, c.radius, c.thickness);
    total_interior += expect_cull_covers_interior<W, H>(ring_n);

    Basis basis_s = make_basis(Quaternion(), Vector(c.tilt, -1.0f, 0.0f));
    SDF::Ring ring_s(basis_s, c.radius, c.thickness);
    total_interior += expect_cull_covers_interior<W, H>(ring_s);
  }
  HS_EXPECT_GT(total_interior, 1000);
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
// Face distance LUT vs exact  (LUT bilinear vs an independent exact oracle)
//
// Face::distance interpolates a 32x32 LUT in sign-pure cells >= one cell-diagonal
// from any edge, else falls back to the exact per-edge scan. The bilinear error
// is NOT bounded by the cell-diagonal (in medial-axis creases it can reach ~4x),
// so this pins two stronger invariants against an exact point-to-polygon oracle:
//   1. SIGN is always correct (sign-purity guard).
//   2. The LUT never serves a near-boundary magnitude — every LUT result is
//      >= cell-diagonal from zero, outside the AA ramp.
// The exact fallback path must reproduce the oracle to float precision.
// ============================================================================

/**
 * @brief Builds one tilted N-gon face and scans its gnomonic box against an exact oracle.
 * @param sides Number of polygon sides (must be <= 8).
 * @param rho Angular circumradius of the face, in radians.
 * @param axis Pole direction the face's basis is built around.
 * @return The number of samples that took the LUT path.
 * @details Asserts sign agreement and bounded magnitude error on the LUT path.
 */
inline int check_face_lut(int sides, float rho, const Vector &axis) {
  constexpr int H = 144;
  constexpr int HV = H + hs::H_OFFSET;
  HS_EXPECT_TRUE(sides <= 8);

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

  const float cell_diag = face.lut_safe_dist;
  HS_EXPECT_GT(cell_diag, 0.0f);

  int lut_samples = 0, sign_mismatches = 0;
  float min_lut_mag = FLT_MAX;
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

      hs::g_scan_metrics.lut_hits = 0;
      hs::g_scan_metrics.exact_hits = 0;
      SDF::DistanceResult res = face.distance(p);
      bool took_lut = hs::g_scan_metrics.lut_hits > 0;
      bool took_exact = hs::g_scan_metrics.exact_hits > 0;
      if (!took_lut && !took_exact)
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
      float exact_angular = fast_atan2(plane_exact, 1.0f);

      if (took_lut) {
        ++lut_samples;
        if ((res.raw_dist < 0.0f) != (exact_angular < 0.0f))
          ++sign_mismatches;
        float mag = std::abs(res.raw_dist);
        if (mag < min_lut_mag)
          min_lut_mag = mag;
      } else {
        HS_EXPECT_NEAR(res.raw_dist, exact_angular, 1e-4f);
      }
    }
  }
  // (1) sign always correct on the LUT path.
  HS_EXPECT_EQ(sign_mismatches, 0);
  // (2) LUT magnitude floor ~atan(cell_diagonal); slack for fast_atan2.
  if (lut_samples > 0) {
    float floor_mag = fast_atan2(cell_diag, 1.0f) - 0.01f;
    HS_EXPECT_GT(min_lut_mag, floor_mag);
  }
  return lut_samples;
}

/**
 * @brief Verifies Face's distance LUT never mis-signs and never serves a near-boundary magnitude.
 * @details Drives check_face_lut across a spread of polygons (triangle, pentagon,
 *   hexagon) and tilts.
 */
inline void test_face_lut_matches_exact_within_cell_diagonal() {
  int lut_samples = 0;
  lut_samples += check_face_lut(/*sides=*/3, 0.45f, Vector(0.4f, 0.3f, 1.0f));
  lut_samples += check_face_lut(/*sides=*/5, 0.50f, Vector(0.4f, 0.3f, 1.0f));
  lut_samples += check_face_lut(/*sides=*/6, 0.40f, Vector(-0.6f, 0.5f, 0.7f));
  // The grid actually fired the LUT path.
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
  auto scope = hs_test::begin_module("sdf");

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

  test_polygon_at_center_inside();
  test_polygon_far_point_outside();

  test_spherical_polygon_center_inside();
  test_spherical_polygon_far_outside();

  test_star_center_inside();
  test_star_far_outside();

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

  test_intersection_requires_both_inside();
  test_intersection_thickness_is_min();
  test_intersection_unsorted_child_yields_sorted_result();
  test_intersection_full_width_child_replays_other();

  test_smooth_union_matches_union_far_from_boundary();
  test_smooth_union_blends_inside_band();
  test_smooth_union_solidity_follows_children();

  test_angular_repeat_matches_base_at_zero_angle();
  test_angular_repeat_creates_copies();
  test_angular_repeat_t_is_sector_local();

  test_cull_covers_interior_over_orientation_grid();
  test_angular_repeat_non_y_axis_cull_covers_copies();
  test_line_arc_bulge_cull_covers_interior();
  test_ring_pole_wrap_cull_covers_interior();
  test_face_cull_covers_aa_fringe();
  test_face_lut_matches_exact_within_cell_diagonal();

  return hs_test::end_module(scope);
}

} // namespace sdf
} // namespace hs_test

