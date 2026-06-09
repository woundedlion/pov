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
#include "core/geometry.h"
#include "tests/test_3dmath.h"
#include "tests/test_harness.h"

#include <utility>
#include <vector>

namespace hs_test {
namespace sdf {

using hs_test::math3d::approx_vec;

// Canonical equator-facing basis: v = +Y, u = +X, w = +Z.
inline Basis equator_basis() { return Basis{Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)}; }

// ============================================================================
// clamp_phi
// ============================================================================

inline void test_clamp_phi_in_range() {
  HS_EXPECT_NEAR(SDF::clamp_phi(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(0.5f), 0.5f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(PI_F), PI_F, 1e-6f);
}

inline void test_clamp_phi_negative_reflects() {
  HS_EXPECT_NEAR(SDF::clamp_phi(-0.3f), 0.3f, 1e-6f);
  HS_EXPECT_NEAR(SDF::clamp_phi(-1.2f), 1.2f, 1e-6f);
}

inline void test_clamp_phi_above_pi_reflects() {
  HS_EXPECT_NEAR(SDF::clamp_phi(PI_F + 0.2f), PI_F - 0.2f, 1e-5f);
  HS_EXPECT_NEAR(SDF::clamp_phi(2.0f * PI_F), 0.0f, 1e-5f);
}

// ============================================================================
// Ring
// ============================================================================

inline void test_ring_on_centerline() {
  Basis b = equator_basis();
  // radius=1 → target_angle = π/2 (equator viewed from +Y)
  SDF::Ring ring(b, 1.0f, 0.1f);

  // A point on the equator (perpendicular to +Y) lies on the centerline.
  auto r = ring.distance(Vector(1, 0, 0));
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f); // dist = raw_dist - thickness = 0 - 0.1
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
}

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

inline void test_ring_outside_band_returns_sentinel() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);

  // North pole — far outside the equatorial ring band
  auto r = ring.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist > 50.0f); // sentinel: 100.0f
}

inline void test_ring_just_outside_band() {
  Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.05f);

  // Point at phi = π/2 - 0.07 from Y-axis: 0.02 outside the 0.05 band edge.
  float off = 0.07f;
  Vector p(std::cos(off), std::sin(off), 0.0f);
  auto r = ring.distance(p);
  // Still within cos_min/cos_max because thickness margin... actually if
  // off > thickness, the point is outside the band. Check returned sentinel.
  HS_EXPECT_TRUE(r.dist > 50.0f);
}

// ============================================================================
// PlanarPolygon  (Basis at top of sphere; distance to nearest edge)
// ============================================================================

inline void test_polygon_at_center_inside() {
  Basis b = equator_basis();
  // 6-gon, thickness=0.5 (radians), phase=0, h_virt=144, height=144
  SDF::PlanarPolygon poly(b, /*r*/ 0.0f, /*thickness*/ 0.5f, /*sides*/ 6,
                          /*phase*/ 0.0f, 144, 144);

  // Point at center of polygon (basis.v direction) is inside.
  auto r = poly.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f);
  // dist == polar(0) * cos(local) - apothem = -apothem
  float apothem = 0.5f * std::cos(PI_F / 6.0f);
  HS_EXPECT_NEAR(r.dist, -apothem, 1e-3f);
}

inline void test_polygon_far_point_outside() {
  Basis b = equator_basis();
  SDF::PlanarPolygon poly(b, 0.0f, 0.3f, 6, 0.0f, 144, 144);

  // Point at the opposite pole (-Y) is π away from the polygon center → outside.
  auto r = poly.distance(Vector(0, -1, 0));
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

// ============================================================================
// SphericalPolygon — great-circle edges
// ============================================================================

inline void test_spherical_polygon_center_inside() {
  Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, /*radius*/ 0.5f, /*sides*/ 5, /*phase*/ 0.0f,
                           144, 144);
  auto r = sp.distance(Vector(0, 1, 0));
  // Center of polygon is strictly inside (negative signed distance).
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

inline void test_spherical_polygon_far_outside() {
  Basis b = equator_basis();
  SDF::SphericalPolygon sp(b, 0.3f, 6, 0.0f, 144, 144);
  auto r = sp.distance(Vector(0, -1, 0)); // antipode
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

// ============================================================================
// Star
// ============================================================================

inline void test_star_center_inside() {
  Basis b = equator_basis();
  SDF::Star star(b, /*radius*/ 0.6f, /*sides*/ 5, /*phase*/ 0.0f, 144, 144);
  auto r = star.distance(Vector(0, 1, 0));
  HS_EXPECT_TRUE(r.dist < 0.0f); // center is interior
}

inline void test_star_far_outside() {
  Basis b = equator_basis();
  SDF::Star star(b, 0.4f, 5, 0.0f, 144, 144);
  auto r = star.distance(Vector(0, -1, 0));
  HS_EXPECT_TRUE(r.dist > 0.0f);
}

// ============================================================================
// Line
// ============================================================================

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

inline void test_line_endpoint_is_on_line() {
  Vector a(1, 0, 0);
  Vector b(0, 0, 1);
  SDF::Line ln(a, b, 0.1f);
  auto r = ln.distance(a);
  HS_EXPECT_NEAR(r.raw_dist, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.dist, -0.1f, 1e-3f);
}

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

inline void test_torus_on_centerline_is_inside() {
  SDF::Torus t{2.0f, 0.5f};
  // Point on the ring centerline → distance = -r (inside, equal to minor radius).
  HS_EXPECT_NEAR(t.distance(Vector(2, 0, 0)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(Vector(0, 0, 2)), -0.5f, 1e-5f);
  HS_EXPECT_NEAR(t.distance(Vector(-2, 0, 0)), -0.5f, 1e-5f);
}

inline void test_torus_on_surface() {
  SDF::Torus t{2.0f, 0.5f};
  // Outer rim point (xz_radius = R + r).
  HS_EXPECT_NEAR(t.distance(Vector(2.5f, 0, 0)), 0.0f, 1e-5f);
  // Inner rim point (xz_radius = R - r).
  HS_EXPECT_NEAR(t.distance(Vector(1.5f, 0, 0)), 0.0f, 1e-5f);
  // Top of tube.
  HS_EXPECT_NEAR(t.distance(Vector(2.0f, 0.5f, 0)), 0.0f, 1e-5f);
}

inline void test_torus_origin_is_outside_hole() {
  SDF::Torus t{2.0f, 0.5f};
  // Origin is in the centre of the donut hole — distance equals R - r = 1.5.
  HS_EXPECT_NEAR(t.distance(Vector(0, 0, 0)), 1.5f, 1e-5f);
}

inline void test_torus_normal_points_outward_on_outer_rim() {
  SDF::Torus t{2.0f, 0.5f};
  // Outer rim — normal should point radially outward in XZ.
  Vector n = t.normal(Vector(2.5f, 0, 0));
  HS_EXPECT_VEC(n, Vector(1, 0, 0), 1e-4f);
}

inline void test_torus_normal_points_outward_on_top() {
  SDF::Torus t{2.0f, 0.5f};
  // Top of tube — normal points +Y.
  Vector n = t.normal(Vector(2.0f, 0.5f, 0));
  HS_EXPECT_VEC(n, Vector(0, 1, 0), 1e-4f);
}

// ============================================================================
// Union — min of distances
// ============================================================================

inline void test_union_picks_closest_shape() {
  // Two disjoint torii. Distance to Union = min of distances.
  SDF::Torus a{2.0f, 0.5f};
  SDF::Torus b{5.0f, 0.5f};
  // Note: Torus does not satisfy the .thickness field required by Union<>.
  // Use Line shapes (which do have .thickness) for the algebra tests instead.

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

inline void test_union_thickness_is_max() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.3f);
  SDF::Union<SDF::Line, SDF::Line> u(la, lb);
  HS_EXPECT_NEAR(u.thickness, 0.3f, 1e-6f);
}

// ============================================================================
// Subtract — max(A, -B)
// ============================================================================

inline void test_subtract_inside_a_outside_b_remains_inside() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.2f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, lb);

  Vector mid_a = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = s.distance(mid_a);
  // Inside la (~-0.2), outside lb (positive large). Max(-0.2, -large_neg) = -0.2.
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

inline void test_subtract_inside_both_becomes_outside() {
  // Same line for A and B → A - A is empty everywhere.
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, la);
  Vector mid = ((Vector(1, 0, 0) + Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = s.distance(mid);
  // dist(A) = -0.1, dist(B) = -0.1. Max(A, -B) = max(-0.1, 0.1) = 0.1 → outside.
  HS_EXPECT_NEAR(r.dist, 0.1f, 1e-3f);
}

// Mock SDF shape that emits a fixed (possibly unsorted, multi-) interval list,
// to exercise Subtract's scanline set-difference independently of any real
// shape. Minimal surface: only thickness, is_solid, and get_horizontal_intervals
// are touched by Subtract's ctor + interval path.
namespace sdf_subtract_detail {
struct MockIntervalShape {
  const std::vector<std::pair<float, float>> *ivs;
  float thickness = 0.1f;
  static constexpr bool is_solid = true;
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out out) const {
    for (const auto &p : *ivs)
      out(p.first, p.second);
    return true;
  }
};

// Mock that falls back to full-width coverage: returning false means "no
// interval restriction", i.e. the shape covers the entire row.
struct MockFullWidthShape {
  float thickness = 0.1f;
  static constexpr bool is_solid = true;
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out) const {
    return false;
  }
};
} // namespace sdf_subtract_detail

// Regression: a child emitting UNSORTED, multi-piece B intervals must still
// yield the correct set difference, emitted in start-sorted order (scan_region's
// coalescer drops any out-of-order interval). Before the fix, Subtract processed
// B in emission order and produced a wrong, unsorted result.
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

// Regression: when B removes nothing (empty), unsorted A intervals must pass
// through in start-sorted order so the coalescer doesn't drop the earlier one.
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

// Regression: when B falls back to full-width (covers the whole row), A - B is
// empty. Before the fix a full-width B was conflated with an empty B, so A was
// passed through, drawing geometry the subtraction should have removed.
inline void test_subtract_full_width_b_yields_empty() {
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
  HS_EXPECT_TRUE(ok); // definitive result (not a full-width fallback)...
  HS_EXPECT_EQ(out.size(), static_cast<size_t>(0)); // ...and it is empty
}

// ============================================================================
// Intersection — max(A, B)
// ============================================================================

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

inline void test_intersection_thickness_is_min() {
  SDF::Line la(Vector(1, 0, 0), Vector(0, 0, 1), 0.1f);
  SDF::Line lb(Vector(-1, 0, 0), Vector(0, 0, -1), 0.3f);
  SDF::Intersection<SDF::Line, SDF::Line> inter(la, lb);
  HS_EXPECT_NEAR(inter.thickness, 0.1f, 1e-6f);
}

// Regression: Intersection's merge-sweep assumes both child interval lists are
// start-sorted. With an UNSORTED multi-interval child it must still produce the
// correct, start-sorted intersection (before the fix it emitted out-of-order
// intervals that scan_region's coalescer then dropped).
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

// ============================================================================
// SmoothUnion — blends at the boundary
// ============================================================================

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

// ============================================================================
// AngularRepeat — folds azimuth into N sectors
// ============================================================================

inline void test_angular_repeat_matches_base_at_zero_angle() {
  Basis b = equator_basis();
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
// Runner
// ============================================================================

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
  test_subtract_full_width_b_yields_empty();

  test_intersection_requires_both_inside();
  test_intersection_thickness_is_min();
  test_intersection_unsorted_child_yields_sorted_result();

  test_smooth_union_matches_union_far_from_boundary();

  test_angular_repeat_matches_base_at_zero_angle();
  test_angular_repeat_creates_copies();

  return hs_test::end_module(scope);
}

} // namespace sdf
} // namespace hs_test

