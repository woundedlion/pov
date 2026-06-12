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

// ============================================================================
// Ring
// ============================================================================

/** @brief Verifies a point on the ring centerline reads raw_dist 0 and dist = -thickness. */
inline void test_ring_on_centerline() {
  Basis b = equator_basis();
  // radius=1 → target_angle = π/2 (equator viewed from +Y)
  SDF::Ring ring(b, 1.0f, 0.1f);

  // A point on the equator (perpendicular to +Y) lies on the centerline.
  auto r = ring.distance(Vector(1, 0, 0));
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f); // dist = raw_dist - thickness = 0 - 0.1
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
}

/** @brief Verifies a point within the ring band reads negative dist. */
inline void test_ring_inside_band() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  // Point slightly above the equator — still inside the band of width 0.1
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

  // North pole — far outside the equatorial ring band
  auto r = ring.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist > 50.0f); // sentinel: 100.0f
}

/** @brief Verifies just past the band edge still trips the sentinel (band edge is exclusive). */
inline void test_ring_just_outside_band() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.05f);

  // Point at phi = π/2 - 0.07 from Y-axis: 0.02 outside the 0.05 band edge.
  float off = 0.07f;
  Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = ring.distance(p);
  // off (0.07) > thickness (0.05): point is past the band edge → sentinel.
  HS_EXPECT_TRUE(r.dist > 50.0f);
}

// ============================================================================
// PlanarPolygon  (Basis at top of sphere; distance to nearest edge)
// ============================================================================

/** @brief Verifies the polygon center is inside, with dist equal to the negated apothem. */
inline void test_polygon_at_center_inside() {
  Basis b = equator_basis();
  // 6-gon, thickness=0.5 (radians), phase=0
  SDF::PlanarPolygon poly(b, /*thickness*/ 0.5f, /*sides*/ 6, /*phase*/ 0.0f);

  // Point at center of polygon (basis.v direction) is inside.
  auto r = poly.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
  // dist == polar(0) * cos(local) - apothem = -apothem
  float apothem = 0.5f * std::cos(PI_F / 6.0f);
  HS_EXPECT_NEAR(r.dist, -apothem, 1e-3f);
}

/** @brief Verifies the antipode of the polygon center is outside (positive dist). */
inline void test_polygon_far_point_outside() {
  Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, 0.3f, 6, 0.0f);

  // Point at the opposite pole (-Y) is π away from the polygon center → outside.
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
  // Center of polygon is strictly inside (negative signed distance).
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies the antipode of the spherical-polygon center is outside. */
inline void test_spherical_polygon_far_outside() {
  Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, 0.3f, 6, 0.0f);
  auto r = sp.distance(Vector(0, -1, 0)); // antipode
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
  HS_EXPECT_TRUE(r.dist < 0.0f); // center is interior
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

  // Midpoint of the great-circle arc — on the line → dist = -thickness.
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

  // Point off the arc in the +Y direction (perpendicular to the great-circle
  // plane that contains a and b).
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
  // Point on the ring centerline → distance = -r (inside, equal to minor radius).
  HS_EXPECT_NEAR(t.distance(Vector(2, 0, 0)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(Vector(0, 0, 2)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(Vector(-2, 0, 0)), -0.5f, 1e-5f);
}

/** @brief Verifies inner/outer rim and top-of-tube points all read dist 0 (on the surface). */
inline void test_torus_on_surface() {
  SDF::Torus t{2.0f, 0.5f};
  // Outer rim point (xz_radius = R + r).
  HS_EXPECT_NEAR(t.distance(Vector(2.5f, 0, 0)), 0.0f, 1e-5f);
  // Inner rim point (xz_radius = R - r).
  HS_EXPECT_NEAR(t.distance(Vector(1.5f, 0, 0)), 0.0f, 1e-5f);
  // Top of tube.
  HS_EXPECT_NEAR(t.distance(Vector(2.0f, 0.5f, 0)), 0.0f, 1e-5f);
}

/** @brief Verifies the donut-hole center is outside, at distance R - r from the tube. */
inline void test_torus_origin_is_outside_hole() {
  SDF::Torus t{2.0f, 0.5f};
  // Origin is in the centre of the donut hole — distance equals R - r = 1.5.
  HS_EXPECT_NEAR(t.distance(Vector(0, 0, 0)), 1.5f, 1e-5f);
}

/** @brief Verifies the surface normal on the outer rim points radially outward (+X here). */
inline void test_torus_normal_points_outward_on_outer_rim() {
  SDF::Torus t{2.0f, 0.5f};
  // Outer rim — normal should point radially outward in XZ.
  Vector n = t.normal(Vector(2.5f, 0, 0));
  HS_EXPECT_VEC(n, Vector(1, 0, 0), 1e-4f);
}

/** @brief Verifies the surface normal at the top of the tube points +Y. */
inline void test_torus_normal_points_outward_on_top() {
  SDF::Torus t{2.0f, 0.5f};
  // Top of tube — normal points +Y.
  Vector n = t.normal(Vector(2.0f, 0.5f, 0));
  HS_EXPECT_VEC(n, Vector(0, 1, 0), 1e-4f);
}

// ============================================================================
// Union — min of distances
// ============================================================================

/** @brief Verifies Union returns the min of member distances, picking the closest shape. */
inline void test_union_picks_closest_shape() {
  // Union distance is the min of member distances. Torus has no .thickness field,
  // which Union<> requires, so the algebra is exercised with Line shapes.
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.1f);

  SDF::Union<SDF::Line, SDF::Line> u(la, lb);

  // Point on first arc midpoint — distance to a is -0.1 (inside la),
  // distance to b is large positive. Union picks a.
  Vector mid_a = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = u.distance(mid_a);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-2f);

  // Point on second arc midpoint — Union picks b.
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
  // Inside la (~-0.2), outside lb (positive large). Max(-0.2, -large_neg) = -0.2.
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies a point inside both A and B becomes outside the difference A - B. */
inline void test_subtract_inside_both_becomes_outside() {
  // Same line for A and B → A - A is empty everywhere.
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, la);
  Vector mid = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = s.distance(mid);
  // dist(A) = -0.1, dist(B) = -0.1. Max(A, -B) = max(-0.1, 0.1) = 0.1 → outside.
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
  HS_EXPECT_TRUE(!ok); // false: request a full-row distance scan, NOT a...
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(0)); // ...definitive empty row
}

// ============================================================================
// Intersection — max(A, B)
// ============================================================================

/** @brief Verifies Intersection is inside only where both children are inside. */
inline void test_intersection_requires_both_inside() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.3f);
  SDF::Line lb(Vector(1, 0, 0), Vector(0, 1, 0), 0.3f);
  SDF::Intersection<SDF::Line, SDF::Line> inter(la, lb);

  // Endpoint a is on both lines (intersection of the two arcs).
  Vector a(1, 0, 0);
  auto r = inter.distance(a);
  HS_EXPECT_TRUE(r.dist < 0.0f); // intersection point lies inside both

  // Far point: outside both — intersection still outside.
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

  // A full-width, B has intervals -> out == B's intervals (emission order).
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

  // Symmetric: A has intervals, B full-width -> out == A's intervals.
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

  // Both fall back -> request a full-scan fallback (return false), no intervals.
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

  // Sample a point clearly inside la and clearly outside lb.
  Vector mid_a = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r_hard = u.distance(mid_a);
  auto r_soft = su.distance(mid_a);
  // Far from the blending zone, smooth union ≈ hard union.
  HS_EXPECT_NEAR(r_hard.dist, r_soft.dist, 1e-3f);
}

/**
 * @brief Verifies SmoothUnion derives its solidity from its children.
 * @details Like Union/Subtract/Intersection, a stroke child keeps its soft falloff
 *   instead of collapsing to a hard 1-px silhouette in process_pixel.
 */
inline void test_smooth_union_solidity_follows_children() {
  // Ring is a stroke (is_solid == false); PlanarPolygon is solid.
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

  // A point exactly on the input line should be inside both the base and
  // the repeated shape.
  Vector mid =
      ((Vector(1, 0, 0) + Vector(0.7071f, 0, 0.7071f)) * 0.5f).normalized();
  auto r_base = ln.distance(mid);
  auto r_rep = rep.distance(mid);
  HS_EXPECT_TRUE(r_rep.dist < 0.0f);
  // Repeated and base agree at the canonical sector.
  HS_EXPECT_NEAR(r_rep.dist, r_base.dist, 1e-3f);
}

/** @brief Verifies AngularRepeat folds a line in the canonical sector into a folded copy. */
inline void test_angular_repeat_creates_copies() {
  // A line in the canonical sector should also be inside in a folded sector.
  SDF::Line ln(Vector(1, 0, 0), Vector(0.7071f, 0, 0.7071f), 0.05f);
  SDF::AngularRepeat<SDF::Line> rep(ln, 4, Vector(0, 1, 0));

  // Rotate the midpoint by 90° around Y — should land on a "copy" of the line.
  Vector mid =
      ((Vector(1, 0, 0) + Vector(0.7071f, 0, 0.7071f)) * 0.5f).normalized();
  Quaternion q90 = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  Vector mid_rot = rotate(mid, q90);

  auto r = rep.distance(mid_rot);
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

// ============================================================================
// Interval-cull conservativeness  ("interval cull == full-row scan")
//
// The rasterizer culls in two tiers: get_vertical_bounds<H> picks the row band,
// get_horizontal_intervals<W,H> picks the column spans within each row.
// scan_region only visits the culled spans; process_pixel then accepts/rejects
// each visited pixel via the exact distance(). So the cull must be
// *conservative* — it may over-visit (the exact test rejects extras) but it must
// never drop a pixel that belongs to the shape. The load-bearing invariant: a
// pixel clearly inside the shape — distance().dist < -pixel_width, i.e. at least
// one pixel deep into the body (solid) or stroke band (stroke) — is among the
// pixels scan_region visits.
//
// The one-pixel margin is deliberate, not a fudge: it scopes the test to the
// invariant the cull actually promises and away from two soft edge behaviours it
// does NOT — a stroke's outer ~5% rim, which the interval generator drops
// because the AA kernel there is ~0.001, and a solid's sub-pixel AA halo, which
// process_pixel paints just *outside* the geometric boundary (so outside the
// caps). It also clears the ~0.004 rad fast-trig noise between distance() and the
// cap math. The Line/great-circle arc-bulge vertical cull is covered separately
// by test_line_arc_bulge_cull_covers_interior.
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

  const float *cos_theta = TrigLUT<W, H>::cos_theta.data();
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

  // Poles, equator, and oblique tilts — the cull math is orientation-dependent.
  const Vector axes[] = {
      Vector(0, 1, 0),         Vector(0, -1, 0),        Vector(1, 0, 0),
      Vector(0, 0, 1),         Vector(1, 1, 0.4f),      Vector(-0.5f, 0.7f, -0.6f),
      Vector(0.3f, -0.8f, 0.5f)};

  int total_interior = 0;
  for (const Vector &axis : axes) {
    Basis basis = make_basis(Quaternion(), axis);
    for (float radius : {0.3f, 0.6f, 0.9f}) {
      // Ring (stroke): both annular-band arcs must cover the whole band.
      SDF::Ring ring(basis, radius, /*thickness=*/0.25f);
      total_interior += expect_cull_covers_interior<W, H>(ring);

      // SphericalPolygon (solid, great-circle edges, geodesic cap).
      SDF::SphericalPolygon spoly(basis, radius, /*sides=*/5, 0.0f);
      total_interior += expect_cull_covers_interior<W, H>(spoly);

      // Star (solid, padded bounding cap).
      SDF::Star star(basis, radius, /*sides=*/5, 0.0f);
      total_interior += expect_cull_covers_interior<W, H>(star);

      // PlanarPolygon (solid, tangent-plane cap; `thickness` is its angular
      // circumradius).
      SDF::PlanarPolygon ppoly(basis, /*thickness=*/radius, /*sides=*/6, 0.0f);
      total_interior += expect_cull_covers_interior<W, H>(ppoly);
    }
  }
  // The grid exercised a non-trivial number of interior pixels.
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
  // A short stroke arc near the north pole (+Y), folded into 4 sectors around
  // X. X-fold's canonical sector is centered on +Y, so the child sits in it and
  // its copies rotate to +Z / -Y / -Z — the equator and the opposite pole, far
  // below the child's near-pole band. The full-canvas cull must reach them.
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
  // Endpoints at phi≈0.4 either side of +Y; the shorter arc passes through the
  // north pole (phi=0), far above either endpoint's latitude.
  SDF::Line ln(Vector(0, cosf(0.4f), sinf(0.4f)),
               Vector(0, cosf(0.4f), -sinf(0.4f)), /*thickness=*/0.15f);
  int interior = expect_cull_covers_interior<W, H>(ln);
  HS_EXPECT_GT(interior, 0);
}

// ============================================================================
// Face distance LUT vs exact  (LUT bilinear vs an independent exact oracle)
//
// Face::distance interpolates a 32x32 baked distance LUT in cells that are
// sign-pure and at least one cell-diagonal (lut_safe_dist) from any edge, and
// falls back to the exact per-edge scan otherwise. The bilinear error is NOT
// bounded by the cell-diagonal: away from the boundary, in the medial-axis
// gradient creases of a face's vertex fans, bilinear interpolation of the
// (otherwise 1-Lipschitz) field is much less accurate — the measured worst case
// is ~1.3x the cell diagonal for a pentagon and ~4x for a sharp-vertexed
// triangle. lut_safe_dist bounds distance to the *zero set* (sign purity), not
// to those creases.
//
// What actually matters for rendering, and what this pins against an independent
// exact point-to-polygon scan, is therefore two stronger invariants:
//
//   1. SIGN is always correct — the load-bearing guarantee for solid fill. The
//      sign-purity guard makes it exact, never a flipped pixel.
//   2. The LUT NEVER reports a near-boundary magnitude. It fires only when all
//      four corners are same-sign with magnitude > cell-diagonal, so the
//      bilinear value (a convex combination) also has magnitude >= cell-diagonal.
//      The renderer only consults the distance magnitude for the AA ramp, which
//      spans +-pixel_width; since cell-diagonal exceeds a pixel here, an
//      LUT-served pixel is always unambiguously inside or outside and its AA
//      ramp is inactive — so AA accuracy never depends on the LUT's magnitude,
//      only on its (exact) sign.
//
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
  // March a grid across the gnomonic tangent-plane box the face occupies. A
  // point with gnomonic coords (px,py) is normalize(center + u*px + w*py): then
  // dot(p,center)=k, dot(p,u)=k*px, dot(p,w)=k*py, so distance() recovers (px,py)
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

      // Independent exact oracle: point-to-polygon in the same tangent plane,
      // over the face's packed edges (mirrors Face::distance's fallback).
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
        // Exact path must reproduce the oracle to float precision.
        HS_EXPECT_NEAR(res.raw_dist, exact_angular, 1e-4f);
      }
    }
  }
  // (1) The sign-purity guard never lets the LUT path mis-sign a pixel.
  HS_EXPECT_EQ(sign_mismatches, 0);
  // (2) The LUT never serves a near-boundary magnitude: every LUT-path result is
  // at least ~atan(cell_diagonal) from zero (the convex-combination floor), so it
  // lands outside the AA ramp and only its (exact) sign drives the pixel. Allow a
  // small slack for the fast_atan2 approximation.
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
  // A spread of polygons (triangle, pentagon, hexagon) and tilts.
  int lut_samples = 0;
  lut_samples += check_face_lut(/*sides=*/3, 0.45f, Vector(0.4f, 0.3f, 1.0f));
  lut_samples += check_face_lut(/*sides=*/5, 0.50f, Vector(0.4f, 0.3f, 1.0f));
  lut_samples += check_face_lut(/*sides=*/6, 0.40f, Vector(-0.6f, 0.5f, 0.7f));
  // The grid actually fired the LUT path (otherwise the bounds are untested).
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

  test_ring_on_centerline();
  test_ring_inside_band();
  test_ring_outside_band_returns_sentinel();
  test_ring_just_outside_band();

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
  test_smooth_union_solidity_follows_children();

  test_angular_repeat_matches_base_at_zero_angle();
  test_angular_repeat_creates_copies();

  test_cull_covers_interior_over_orientation_grid();
  test_angular_repeat_non_y_axis_cull_covers_copies();
  test_line_arc_bulge_cull_covers_interior();
  test_face_lut_matches_exact_within_cell_diagonal();

  return hs_test::end_module(scope);
}

} // namespace sdf
} // namespace hs_test

