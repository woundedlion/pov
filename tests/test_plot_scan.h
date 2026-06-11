/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/plot.h (and the pure helpers reached via scan.h's
 * dependency on constants.h ClipRegion).
 *
 * Focus: PURE sampling / geometry paths that produce Fragments from geometry
 * WITHOUT a live Canvas. The draw()/rasterize() paths (which require an
 * Effect + Canvas) and the entire Scan:: rasterizer are intentionally NOT
 * exercised here — see the SKIPPED note at the bottom of this file.
 *
 * Coverage:
 *   - Plot::Line::sample        : geodesic endpoints unit-length, span, v0/v1;
 *                                 degenerate (zero-length) segment handled.
 *   - ClipRegion::could_intersect_y : y-range cull culling (constants.h).
 *   - Plot::Ring::calcPoint / Ring::sample : unit-length, angular progress.
 *   - Plot::DistortedRing::sample : angle-addition identity (LUT) matches
 *                                   direct cos/sin within tolerance.
 *   - Plot::Spiral::sample      : unit-length, monotone arc length.
 *   - Plot::Multiline::sample   : arc-length parameterization, v0 in [0,1].
 *   - Plot::Bezier::sample      : unit-length-ish, monotone arc length,
 *                                 endpoints land on control points.
 *   - Plot::Star::sample / Flower::sample : unit-length, closed loop.
 */
#pragma once

#include "core/plot.h"
#include "core/scan.h"
#include "core/geometry.h"
#include "core/canvas.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace plot_scan_tests {

// ---------------------------------------------------------------------------
// Local arena for sampling. Plot::*::sample takes a caller-bound Fragments,
// so we provide our own arena rather than relying on the global scratch arena.
// ---------------------------------------------------------------------------

inline uint8_t plot_scan_arena_buf[256 * 1024];

inline Arena &plot_arena() {
  static Arena a(plot_scan_arena_buf, sizeof(plot_scan_arena_buf));
  return a;
}

// ============================================================================
// Plot::Line::sample
// ============================================================================

inline void test_line_sample_endpoints_and_unit_length() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  b.pos = Vector(0, 1, 0);
  const int density = 8;
  Plot::Line::sample(points, a, b, density);

  HS_EXPECT_EQ(points.size(), (size_t)(density + 1));

  // Endpoints land exactly on the inputs.
  HS_EXPECT_NEAR(points[0].pos.x, a.pos.x, 1e-6f);
  HS_EXPECT_NEAR(points[0].pos.y, a.pos.y, 1e-6f);
  HS_EXPECT_NEAR(points[density].pos.x, b.pos.x, 1e-6f);
  HS_EXPECT_NEAR(points[density].pos.y, b.pos.y, 1e-6f);

  // All sampled positions are unit length (great-circle interpolation).
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }

  // v0 (progress) spans 0..1; v1 (arc length) ends at the segment angle.
  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points[density].v0, 1.0f, 1e-6f);
  float total_angle = angle_between(a.pos, b.pos);
  HS_EXPECT_NEAR(points[density].v1, total_angle, 1e-4f);
  HS_EXPECT_NEAR(total_angle, PI_F * 0.5f, 1e-4f); // x-axis to y-axis = 90deg
}

inline void test_line_sample_interior_between_endpoints() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  b.pos = Vector(0, 1, 0);
  Plot::Line::sample(points, a, b, 4);

  // Interior samples must be finite and lie on the minor arc: their angle to
  // each endpoint never exceeds the total span.
  float total = angle_between(a.pos, b.pos);
  for (size_t i = 1; i + 1 < points.size(); ++i) {
    const Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_LE(angle_between(a.pos, p), total + 1e-3f);
    HS_EXPECT_LE(angle_between(b.pos, p), total + 1e-3f);
  }
}

inline void test_line_sample_degenerate_segment() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = Vector(0, 0, 1);
  b.pos = Vector(0, 0, 1); // identical -> zero-length

  Plot::Line::sample(points, a, b, 8);

  // Degenerate path emits a dot (two coincident fragments), no NaN, no crash.
  HS_EXPECT_EQ(points.size(), (size_t)2);
  for (size_t i = 0; i < points.size(); ++i) {
    const Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }
}

inline void test_line_sample_antipodal_stable_axis() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  // Antipodal endpoints: angle == pi, so cross(a, b) == 0. Without a stable-axis
  // fallback this normalizes a ~zero vector into a garbage axis, collapsing the
  // interior to the start point (or NaN). The fix picks a perpendicular axis.
  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  b.pos = Vector(-1, 0, 0);
  const int density = 8;
  Plot::Line::sample(points, a, b, density);

  HS_EXPECT_EQ(points.size(), (size_t)(density + 1));
  HS_EXPECT_NEAR(points[0].pos.x, a.pos.x, 1e-6f);
  HS_EXPECT_NEAR(points[density].pos.x, b.pos.x, 1e-6f);

  // Every sample is finite and unit length.
  for (size_t i = 0; i < points.size(); ++i) {
    const Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }

  // Discriminating check: the midpoint is a real point ~90deg from each
  // antipodal endpoint (without the fix it collapses to the start, angle 0).
  const Vector &mid = points[density / 2].pos;
  HS_EXPECT_NEAR(angle_between(a.pos, mid), PI_F * 0.5f, 1e-3f);
  HS_EXPECT_NEAR(angle_between(b.pos, mid), PI_F * 0.5f, 1e-3f);
}

// ============================================================================
// ClipRegion::could_intersect_y  (constants.h — pure clip culling)
// ============================================================================

inline void test_clip_could_intersect_y() {
  // A clip window covering the middle band of a 144-row canvas.
  ClipRegion cr;
  cr.y_start = 40;
  cr.y_end = 80;
  cr.x_start = 0;
  cr.x_end = MAX_W;
  cr.margin = 0;
  cr.h = MAX_H;
  cr.w = MAX_W;

  HS_EXPECT_FALSE(cr.is_full());

  // Render band (margin 0) is [40, 80).
  HS_EXPECT_EQ(cr.render_y_start(), 40);
  HS_EXPECT_EQ(cr.render_y_end(), 80);

  // Segment fully inside the band -> NOT culled.
  HS_EXPECT_TRUE(cr.could_intersect_y(50.0f, 60.0f));

  // Segment straddling the lower edge -> NOT culled.
  HS_EXPECT_TRUE(cr.could_intersect_y(30.0f, 45.0f));

  // Segment fully ABOVE the band (smaller y) -> culled.
  HS_EXPECT_FALSE(cr.could_intersect_y(0.0f, 39.0f));

  // Segment fully BELOW the band (larger y) -> culled.
  HS_EXPECT_FALSE(cr.could_intersect_y(80.0f, 120.0f));

  // Order independence: arguments swapped give the same answer.
  HS_EXPECT_TRUE(cr.could_intersect_y(60.0f, 50.0f));
  HS_EXPECT_FALSE(cr.could_intersect_y(120.0f, 80.0f));
}

// ============================================================================
// Plot::edge_row_span — arc-aware clip cull (CODE_REVIEW #2)
// ============================================================================

// An orthonormal basis from a unit normal (mirrors make_basis's construction
// without needing a quaternion). v is the normal; u, w span the tangent plane.
inline Basis basis_from_normal(const Vector &n) {
  Vector v = n.normalized();
  Vector ref = std::abs(dot(v, X_AXIS)) > math::COS_AXIS_PARALLEL ? Y_AXIS : X_AXIS;
  Vector u = cross(v, ref).normalized();
  Vector w = cross(v, u).normalized();
  return {u, v, w};
}

// edge_row_span must conservatively cover the actual rendered arc's screen-row
// extent — including the interior latitude bulge that an endpoint-only test
// silently drops (the curve-clip-cull gap). We densely sample the true arc for
// both the geodesic and planar strategies and assert the span contains it, and
// that the bulge case (where endpoints alone would under-cover) is exercised.
inline void test_edge_row_span_covers_arc_bulge() {
  constexpr int TW = 288, TH = 144;
  auto row_of = [](const Vector &v) {
    return vector_to_pixel<TW, TH>(v.normalized()).y;
  };

  hs::random().seed(20260609);
  int bulge_cases = 0; // edges where endpoints alone under-cover the true arc

  for (int trial = 0; trial < 6000; ++trial) {
    const bool planar = (trial & 1);
    Vector a, b;
    Basis basis;
    const Basis *pb = nullptr;

    if (planar) {
      // A planar-polygon edge: two points on a disk of angular radius `radius`
      // about a random center, joined by an azimuthal-equidistant straight line.
      Vector center(hs::rand_f(-1, 1), hs::rand_f(-1, 1), hs::rand_f(-1, 1));
      if (center.length() < 0.1f) continue;
      basis = basis_from_normal(center.normalized());
      pb = &basis;
      float radius = hs::rand_f(0.2f, 1.4f);
      auto on_disk = [&](float ang) {
        Vector dir = basis.u * cosf(ang) + basis.w * sinf(ang);
        return (basis.v * cosf(radius) + dir * sinf(radius)).normalized();
      };
      float a0 = hs::rand_f(0, 2 * PI_F);
      a = on_disk(a0);
      b = on_disk(a0 + hs::rand_f(0.3f, 2.3f));
      // Antipodal-seam edges fall back to the geodesic strategy; skip them.
      if (dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
          dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
        continue;
    } else {
      Vector ra(hs::rand_f(-1, 1), hs::rand_f(-1, 1), hs::rand_f(-1, 1));
      Vector rb(hs::rand_f(-1, 1), hs::rand_f(-1, 1), hs::rand_f(-1, 1));
      if (ra.length() < 0.1f || rb.length() < 0.1f) continue;
      a = ra.normalized();
      b = rb.normalized();
      if (angle_between(a, b) < 0.05f) continue;
    }

    // Dense ground truth: the true rendered arc's row extent.
    float t_lo = 1e9f, t_hi = -1e9f;
    constexpr int N = 2000;
    std::pair<float, float> p1{}, p2{};
    float ang = 0.0f;
    Vector vperp = a;
    if (planar) {
      p1 = Plot::azimuthal_project(a, basis);
      p2 = Plot::azimuthal_project(b, basis);
    } else {
      ang = angle_between(a, b);
      vperp = cross(cross(a, b).normalized(), a);
    }
    for (int i = 0; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      Vector p = planar ? Plot::azimuthal_unproject(
                              p1.first + (p2.first - p1.first) * t,
                              p1.second + (p2.second - p1.second) * t, basis)
                        : (a * cosf(ang * t) + vperp * sinf(ang * t));
      float r = row_of(p);
      t_lo = std::min(t_lo, r);
      t_hi = std::max(t_hi, r);
    }

    float n_lo, n_hi;
    Plot::edge_row_span<TW, TH>(a, b, pb, n_lo, n_hi);

    // The span must conservatively contain the true arc (sub-pixel tolerance for
    // fast-math; the real render band additionally expands by >=1px margin).
    HS_EXPECT_LE(n_lo, t_lo + 0.25f);
    HS_EXPECT_GE(n_hi, t_hi - 0.25f);

    // Count edges whose interior bulges past the endpoints — the case the old
    // endpoint-only cull dropped. Confirms the test actually exercises it.
    float e_lo = std::min(row_of(a), row_of(b));
    float e_hi = std::max(row_of(a), row_of(b));
    if (t_lo < e_lo - 1.0f || t_hi > e_hi + 1.0f) bulge_cases++;
  }

  // The randomized sweep must have produced many genuine bulge cases, else the
  // test would pass vacuously without ever exercising the fix.
  HS_EXPECT_GT(bulge_cases, 500);
}

// ============================================================================
// Plot::Ring::calcPoint / Ring::sample
// ============================================================================

inline void test_ring_calc_point_unit_length() {
  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  for (float a = 0.0f; a < 6.28f; a += 0.4f) {
    Vector p = Plot::Ring::calcPoint(a, 0.5f, b.u, b.v, b.w);
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-4f);
  }
}

inline void test_ring_sample_unit_length_and_progress() {
  ScratchScope sc(plot_arena());
  Fragments points;
  const int N = 32;
  points.bind(plot_arena(), N + 2);

  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  Plot::Ring::sample(points, b, 0.5f, N, 0.0f);

  // N samples + 1 manual-close overlap fragment.
  HS_EXPECT_EQ(points.size(), (size_t)(N + 1));

  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
  // Angular progress v0 runs 0 .. 1.
  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points[points.size() - 1].v0, 1.0f, 1e-6f);

  // Closing fragment coincides with the first sample (closed ring).
  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);
}

// Ring::sample<W,H> (Ring::draw's full-resolution path) builds the ring from
// the TrigLUT angle-addition identity instead of a per-sample libm cos/sin.
// Verify each fragment — position and the analytic arc-length register — matches
// a direct cos/sin(theta+phase) construction of the same ring.
inline void test_ring_sample_lut_matches_direct() {
  constexpr int W = 64;
  constexpr int H = 64;

  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), W + 2);

  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  const float radius = 0.5f;
  const float phase = 0.7f;
  Plot::Ring::sample<W, H>(points, b, radius, phase);

  HS_EXPECT_EQ(points.size(), (size_t)(W + 1));

  auto res = get_antipode(b, radius);
  const Basis &wb = res.first;
  float wr = res.second;
  const float theta_eq = wr * (PI_F / 2.0f);
  const float r_val = sinf(theta_eq);
  const float d_val = cosf(theta_eq);
  const float arc_scale = sinf(wr);
  const float step = 2.0f * PI_F / W;

  for (int i = 0; i < W; ++i) {
    float theta = i * step;
    float t = theta + phase;
    Vector u_temp = (wb.u * cosf(t)) + (wb.w * sinf(t));
    Vector expected = ((wb.v * d_val) + (u_temp * r_val)).normalized();
    HS_EXPECT_NEAR(points[i].pos.x, expected.x, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.y, expected.y, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.z, expected.z, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
    HS_EXPECT_NEAR(points[i].v1, theta * arc_scale, 2e-3f);
  }

  // Close vertex coincides with sample 0 and carries v0 == 1.
  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-6f);
}

// ============================================================================
// Plot::DistortedRing::sample  — angle-addition identity (LUT) vs direct
// ============================================================================

inline void test_distorted_ring_sample_angle_addition_identity() {
  // With a zero shift function, DistortedRing::sample reduces to a plain ring
  // built via the TrigLUT angle-addition identity cos/sin(theta+phase). We
  // verify each fragment matches a direct cos/sin construction of the same
  // ring, confirming the precomputed-LUT identity is correct.
  constexpr int W = 64;
  constexpr int H = 64;

  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), W + 2);

  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  const float radius = 0.5f;
  const float phase = 0.7f;
  ScalarFn zero_shift = [](float) { return 0.0f; };

  Plot::DistortedRing::sample<W, H>(points, b, radius, zero_shift, phase);

  HS_EXPECT_EQ(points.size(), (size_t)(W + 1));

  // Reconstruct the ring directly (no LUT) and compare.
  auto res = get_antipode(b, radius);
  const Basis &wb = res.first;
  float wr = res.second;
  const float theta_eq = wr * (PI_F / 2.0f);
  const float r_val = sinf(theta_eq);
  const float d_val = cosf(theta_eq);
  const float step = 2.0f * PI_F / W;

  for (int i = 0; i < W; ++i) {
    float theta = i * step;
    float t = theta + phase;
    Vector u_temp = (wb.u * cosf(t)) + (wb.w * sinf(t));
    Vector expected = ((wb.v * d_val) + (u_temp * r_val)).normalized();

    HS_EXPECT_NEAR(points[i].pos.x, expected.x, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.y, expected.y, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.z, expected.z, 2e-3f);
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
}

// ============================================================================
// Plot::Spiral::sample
// ============================================================================

inline void test_spiral_sample_unit_length_and_monotone_arc() {
  ScratchScope sc(plot_arena());
  Fragments frags;
  const int N = 48;
  frags.bind(plot_arena(), N);

  Plot::Spiral::sample(frags, N, 0.5f);
  HS_EXPECT_EQ(frags.size(), (size_t)N);

  float last_v1 = -1.0f;
  for (size_t i = 0; i < frags.size(); ++i) {
    HS_EXPECT_NEAR(frags[i].pos.length(), 1.0f, 1e-3f);
    // Cumulative arc length is non-decreasing.
    HS_EXPECT_GE(frags[i].v1, last_v1);
    last_v1 = frags[i].v1;
  }
  HS_EXPECT_NEAR(frags[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(frags[N - 1].v0, 1.0f, 1e-6f);
}

// ============================================================================
// Plot::Multiline::sample
// ============================================================================

inline void test_multiline_sample_arclength_param() {
  ScratchScope sc(plot_arena());

  // Source vertices as Fragments (sample reads .pos).
  Fragments verts;
  verts.bind(plot_arena(), 4);
  Fragment v;
  v.pos = Vector(1, 0, 0);  verts.push_back(v);
  v.pos = Vector(0, 1, 0);  verts.push_back(v);
  v.pos = Vector(-1, 0, 0); verts.push_back(v);

  Fragments points;
  points.bind(plot_arena(), 8);
  Plot::Multiline::sample(points, verts, /*closed=*/false);

  HS_EXPECT_EQ(points.size(), (size_t)3);

  // v0 normalized progress: first 0, last 1, monotone.
  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-4f);
  float last = -1.0f;
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_GE(points[i].v0, last - 1e-6f);
    last = points[i].v0;
    // v2 carries the vertex index.
    HS_EXPECT_NEAR(points[i].v2, static_cast<float>(i), 1e-6f);
  }
  // v1 cumulative arc length equals the sum of the two equal 90deg hops.
  HS_EXPECT_NEAR(points.back().v1, PI_F, 1e-3f);
}

// ============================================================================
// Plot::Bezier::sample
// ============================================================================

inline void test_bezier_sample_endpoints_and_monotone_arc() {
  ScratchScope sc(plot_arena());
  Fragments points;
  const int N = 24;
  points.bind(plot_arena(), N + 1);

  Vector p0(1, 0, 0);
  Vector p1 = Vector(0.7f, 0.7f, 0.0f).normalized();
  Vector p2 = Vector(0.0f, 0.7f, 0.7f).normalized();
  Vector p3(0, 0, 1);

  Plot::Bezier::sample(points, p0, p1, p2, p3, N);
  HS_EXPECT_EQ(points.size(), (size_t)(N + 1));

  // Endpoints land on the outer control points.
  HS_EXPECT_NEAR(points[0].pos.x, p0.x, 1e-3f);
  HS_EXPECT_NEAR(points[0].pos.z, p0.z, 1e-3f);
  HS_EXPECT_NEAR(points[N].pos.x, p3.x, 1e-3f);
  HS_EXPECT_NEAR(points[N].pos.z, p3.z, 1e-3f);

  // Geodesic-mode samples stay on the unit sphere; arc length is monotone.
  float last_v1 = -1.0f;
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 5e-3f);
    HS_EXPECT_GE(points[i].v1, last_v1);
    last_v1 = points[i].v1;
  }
  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points[N].v0, 1.0f, 1e-6f);
}

// ============================================================================
// Plot::Star::sample / Plot::Flower::sample
// ============================================================================

inline void test_star_sample_unit_length_closed() {
  ScratchScope sc(plot_arena());
  Fragments points;
  const int sides = 5;
  points.bind(plot_arena(), sides * 2 + 2);

  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  Plot::Star::sample(points, b, 0.5f, sides, 0.0f);

  // 2*sides vertices + 1 close fragment.
  HS_EXPECT_EQ(points.size(), (size_t)(sides * 2 + 1));
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
  // Closing fragment matches the first vertex.
  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-6f);
}

inline void test_flower_sample_unit_length_closed() {
  ScratchScope sc(plot_arena());
  Fragments points;
  const int sides = 6;
  points.bind(plot_arena(), sides * 2 + 2);

  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
  Plot::Flower::sample(points, b, 0.5f, sides, 0.0f);

  HS_EXPECT_EQ(points.size(), (size_t)(sides * 2 + 1));
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }
  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
}

// ============================================================================
// Runner
// ============================================================================

inline int run_plot_scan_tests() {
  auto scope = hs_test::begin_module("plot_scan");

  hs::random().seed(1337);

  test_line_sample_endpoints_and_unit_length();
  test_line_sample_interior_between_endpoints();
  test_line_sample_degenerate_segment();
  test_line_sample_antipodal_stable_axis();

  test_clip_could_intersect_y();
  test_edge_row_span_covers_arc_bulge();

  test_ring_calc_point_unit_length();
  test_ring_sample_unit_length_and_progress();
  test_ring_sample_lut_matches_direct();

  test_distorted_ring_sample_angle_addition_identity();

  test_spiral_sample_unit_length_and_monotone_arc();
  test_multiline_sample_arclength_param();
  test_bezier_sample_endpoints_and_monotone_arc();

  test_star_sample_unit_length_closed();
  test_flower_sample_unit_length_closed();

  return hs_test::end_module(scope);
}

} // namespace plot_scan_tests
} // namespace hs_test

