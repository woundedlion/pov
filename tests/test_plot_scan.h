/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/render/plot.h (and the pure helpers reached via scan.h's
 * dependency on constants.h ClipRegion).
 *
 * Focus: PURE sampling / geometry paths that produce Fragments from geometry
 * WITHOUT a live Canvas. The Scan:: rasterizer is covered in test_scan.h.
 *
 * Coverage:
 *   - Plot::Line::sample        : geodesic endpoints unit-length, span, v0/v1;
 *                                 degenerate (zero-length) segment handled.
 *   - ClipRegion::could_intersect_y : y-range cull (constants.h).
 *   - ClipRegion x-clip (render_x_*, contains_x, x_clip/XClip) : cylindrical
 *                                 band topologies + contains_x/XClip parity.
 *   - Plot::Ring::sample : unit-length, angular progress.
 *   - Plot::DistortedRing::sample : angle-addition identity (LUT) matches
 *                                   direct cos/sin within tolerance.
 *   - Plot::Spiral::sample      : unit-length, monotone arc length.
 *   - Plot::Multiline::sample   : arc-length parameterization, v0 in [0,1].
 *   - Plot::Star::sample / Flower::sample : unit-length, closed loop.
 */
#pragma once

#include "core/render/plot.h"
#include "core/render/scan.h"
#include "core/render/filter.h"
#include "core/math/geometry.h"
#include "core/render/canvas.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <vector>

namespace hs_test {
namespace plot_scan_tests {

// ---------------------------------------------------------------------------
// Local arena for sampling.
// ---------------------------------------------------------------------------

/** @brief Backing storage for the module-local sampling arena. */
inline uint8_t plot_scan_arena_buf[256 * 1024];

/**
 * @brief Returns the module-local arena used to bind Fragments for sampling.
 * @return Reference to a function-static Arena over plot_scan_arena_buf.
 * @details Plot::*::sample takes a caller-bound Fragments, so this provides a
 *          dedicated arena rather than relying on the global scratch arena.
 *          The arena is not reset between tests: every caller must wrap its
 *          allocations in a ScratchScope, or allocations leak across tests and
 *          couple correctness to run order.
 */
inline Arena &plot_arena() {
  static Arena a(plot_scan_arena_buf, sizeof(plot_scan_arena_buf));
  return a;
}

// ---------------------------------------------------------------------------
// Headless rasterize() harness. A stub pipeline records plotted world positions;
// a no-op Effect supplies a Canvas with a full (unclipped) clip band.
// ---------------------------------------------------------------------------

/**
 * @brief Pipeline stub: records each plotted position; ignores color/age.
 * @details Carries both plot() overloads so it can also back a type-erased
 * PipelineRef (the Plot::ParticleSystem draw() entry points take one); only
 * the 3D world-space overload records — the 2D screen-space form is unused by
 * these paths.
 */
struct CapturePipeline {
  std::vector<Vector> plotted; /**< World positions handed to plot(), in order. */
  void plot(Canvas &, const Vector &v, const Pixel &, float, float) {
    plotted.push_back(v);
  }
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
};

/** @brief Minimal Effect backing a Canvas (no per-frame draw, no background). */
struct RasterFx : public Effect {
  RasterFx(int W, int H) : Effect(W, H) {}
  void draw_frame() override {}
};

/** @brief Identity fragment shader (leaves the fragment untouched). */
inline void noop_shader(const Vector &, Fragment &) {}

/**
 * @brief Largest angular gap (radians) between consecutive recorded positions.
 * @param pts Plotted positions in plot order.
 * @param wrap When true, also measures the gap from the last point back to the
 *             first (closed-loop seam continuity).
 */
inline float max_consecutive_gap(const std::vector<Vector> &pts, bool wrap) {
  float worst = 0.0f;
  for (size_t i = 1; i < pts.size(); ++i)
    worst = std::max(worst, angle_between(pts[i - 1], pts[i]));
  if (wrap && pts.size() >= 2)
    worst = std::max(worst, angle_between(pts.back(), pts.front()));
  return worst;
}

// ============================================================================
// Plot::Line::sample
// ============================================================================

/**
 * @brief Verifies Line::sample emits density+1 great-circle fragments with exact
 *        endpoints, unit-length positions, v0 progress spanning 0..1, and v1 arc
 *        length ending at the segment's subtended angle.
 */
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

  HS_EXPECT_NEAR(points[0].pos.x, a.pos.x, 1e-6f);
  HS_EXPECT_NEAR(points[0].pos.y, a.pos.y, 1e-6f);
  HS_EXPECT_NEAR(points[density].pos.x, b.pos.x, 1e-6f);
  HS_EXPECT_NEAR(points[density].pos.y, b.pos.y, 1e-6f);

  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_NEAR(points[i].pos.length(), 1.0f, 1e-3f);
  }

  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points[density].v0, 1.0f, 1e-6f);
  float total_angle = angle_between(a.pos, b.pos);
  HS_EXPECT_NEAR(points[density].v1, total_angle, 1e-4f);
  HS_EXPECT_NEAR(total_angle, PI_F * 0.5f, 1e-4f);
}

/**
 * @brief Verifies interior Line::sample fragments stay finite and lie on the
 *        minor arc, never farther from either endpoint than the total span.
 */
inline void test_line_sample_interior_between_endpoints() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  b.pos = Vector(0, 1, 0);
  Plot::Line::sample(points, a, b, 4);

  float total = angle_between(a.pos, b.pos);
  for (size_t i = 1; i + 1 < points.size(); ++i) {
    const Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_LE(angle_between(a.pos, p), total + 1e-3f);
    HS_EXPECT_LE(angle_between(b.pos, p), total + 1e-3f);
  }
}

/**
 * @brief Verifies a zero-length segment (coincident endpoints) emits a dot: two
 *        coincident, finite, unit-length fragments rather than NaN or a crash.
 */
inline void test_line_sample_degenerate_segment() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = Vector(0, 0, 1);
  b.pos = Vector(0, 0, 1);

  Plot::Line::sample(points, a, b, 8);

  HS_EXPECT_EQ(points.size(), (size_t)2);
  for (size_t i = 0; i < points.size(); ++i) {
    const Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }
}

/**
 * @brief Verifies Line::sample picks a stable perpendicular axis for antipodal
 *        endpoints so the arc stays finite, unit-length, and passes through a
 *        real ~90deg midpoint.
 * @details Antipodal endpoints (angle == pi) make cross(a, b) == 0, so the
 *          rotation axis is degenerate; a perpendicular fallback is required to
 *          avoid collapsing to the start point or NaN.
 */
inline void test_line_sample_antipodal_stable_axis() {
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 16);

  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  b.pos = Vector(-1, 0, 0);
  const int density = 8;
  Plot::Line::sample(points, a, b, density);

  HS_EXPECT_EQ(points.size(), (size_t)(density + 1));
  HS_EXPECT_NEAR(points[0].pos.x, a.pos.x, 1e-6f);
  HS_EXPECT_NEAR(points[density].pos.x, b.pos.x, 1e-6f);

  for (size_t i = 0; i < points.size(); ++i) {
    const Vector &p = points[i].pos;
    HS_EXPECT_TRUE(std::isfinite(p.x) && std::isfinite(p.y) &&
                   std::isfinite(p.z));
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
  }

  // The midpoint is ~90deg from each antipodal endpoint.
  const Vector &mid = points[density / 2].pos;
  HS_EXPECT_NEAR(angle_between(a.pos, mid), PI_F * 0.5f, 1e-3f);
  HS_EXPECT_NEAR(angle_between(b.pos, mid), PI_F * 0.5f, 1e-3f);
}

// ============================================================================
// ClipRegion::could_intersect_y  (constants.h — pure clip culling)
// ============================================================================

/**
 * @brief Verifies ClipRegion::could_intersect_y culls by screen-row range:
 *        segments overlapping the clip band pass, fully-above/below segments are
 *        rejected, and the test is order-independent in its two y arguments.
 */
inline void test_clip_could_intersect_y() {
  ClipRegion cr;
  cr.y_start = 40;
  cr.y_end = 80;
  cr.x_start = 0;
  cr.x_end = MAX_W;
  cr.margin = 0;
  cr.h = MAX_H;
  cr.w = MAX_W;

  HS_EXPECT_FALSE(cr.is_full());

  HS_EXPECT_EQ(cr.render_y_start(), 40);
  HS_EXPECT_EQ(cr.render_y_end(), 80);

  HS_EXPECT_TRUE(cr.could_intersect_y(50.0f, 60.0f));
  HS_EXPECT_TRUE(cr.could_intersect_y(30.0f, 45.0f));

  HS_EXPECT_FALSE(cr.could_intersect_y(0.0f, 39.0f));
  HS_EXPECT_FALSE(cr.could_intersect_y(80.0f, 120.0f));

  // Order independence: swapped arguments give the same answer.
  HS_EXPECT_TRUE(cr.could_intersect_y(60.0f, 50.0f));
  HS_EXPECT_FALSE(cr.could_intersect_y(120.0f, 80.0f));
}

/**
 * @brief Asserts contains_x() and !XClip::clipped() agree on every column.
 * @details The runtime hot loops pick one predicate or the other (filters query
 *          contains_x() per fragment; the scanline rasterizer precomputes XClip
 *          once per draw), so any divergence would blank or leak pixels on one
 *          path but not the other. Parity over [0, w) pins them together.
 */
inline void expect_xclip_parity(const ClipRegion &cr) {
  const ClipRegion::XClip xc = cr.x_clip();
  for (int x = 0; x < cr.w; ++x) {
    HS_EXPECT_EQ(cr.contains_x(x), !xc.clipped(x));
  }
}

/**
 * @brief Exercises the cylindrical x-clip predicates (render_x_*, contains_x,
 *        x_clip/XClip) across every documented band topology.
 * @details Covers a non-wrapping sub-band, a partial seam-crossing band
 *          (rs > re), the wrap-to-full fold (rs == re forcing active == false),
 *          an over-wrap band whose display width plus both margins exceeds w
 *          (must read as full coverage, not a thin wrap sliver), and the
 *          explicit full-width band (x_end - x_start >= w). Each case checks the
 *          XClip flags and membership at the band edges, then asserts
 *          contains_x()/XClip parity over all columns. Uses the hardware-like
 *          w = 96 so a small margin exercises a real seam cross.
 */
inline void test_clip_x_band_topologies() {
  constexpr int W = 96;

  auto make = [](int x0, int x1, int margin) {
    ClipRegion cr;
    cr.y_start = 0;
    cr.y_end = MAX_H;
    cr.x_start = x0;
    cr.x_end = x1;
    cr.margin = margin;
    cr.w = W;
    cr.h = MAX_H;
    return cr;
  };

  // 1) Non-wrapping sub-band: [20,60) expanded by 2 -> render band [18,62).
  {
    ClipRegion cr = make(20, 60, 2);
    HS_EXPECT_EQ(cr.render_x_start(), 18);
    HS_EXPECT_EQ(cr.render_x_end(), 62);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_TRUE(xc.active);
    HS_EXPECT_FALSE(xc.wrap);
    HS_EXPECT_EQ(xc.rs, 18);
    HS_EXPECT_EQ(xc.re, 62);
    HS_EXPECT_FALSE(cr.contains_x(17));
    HS_EXPECT_TRUE(cr.contains_x(18));
    HS_EXPECT_TRUE(cr.contains_x(61));
    HS_EXPECT_FALSE(cr.contains_x(62)); // re exclusive
    expect_xclip_parity(cr);
  }

  // 2) Partial seam-crossing band: [2,90) expanded by 3 -> render band wraps to
  //    [95, w) U [0, 93), i.e. rs (95) > re (93); the 2-column gap {93,94} stays
  //    clipped (display width + both margins = 94 < w, so it is a true sub-arc).
  {
    ClipRegion cr = make(2, 90, 3);
    HS_EXPECT_EQ(cr.render_x_start(), 95);
    HS_EXPECT_EQ(cr.render_x_end(), 93);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_TRUE(xc.active);
    HS_EXPECT_TRUE(xc.wrap);
    HS_EXPECT_TRUE(cr.contains_x(0));
    HS_EXPECT_TRUE(cr.contains_x(92));  // last column before the gap
    HS_EXPECT_FALSE(cr.contains_x(93)); // re exclusive
    HS_EXPECT_FALSE(cr.contains_x(94)); // interior of the clipped gap
    HS_EXPECT_TRUE(cr.contains_x(95));  // rs inclusive
    expect_xclip_parity(cr);
  }

  // 3) Wrap-to-full fold: [10,50) with margin 28 makes both edges land on the
  //    same column (rs == re == 78), so the margin expansion has wrapped to
  //    cover the full width and XClip must deactivate rather than read as empty.
  {
    ClipRegion cr = make(10, 50, 28);
    HS_EXPECT_EQ(cr.render_x_start(), 78);
    HS_EXPECT_EQ(cr.render_x_end(), 78);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_FALSE(xc.active); // rs == re folds to "no clipping", not empty band
    HS_EXPECT_TRUE(cr.contains_x(0));
    HS_EXPECT_TRUE(cr.contains_x(78));
    HS_EXPECT_TRUE(cr.contains_x(95));
    expect_xclip_parity(cr);
  }

  // 4) Explicit full-width band: x_end - x_start >= w covers everything
  //    regardless of margin.
  {
    ClipRegion cr = make(0, W, 0);
    HS_EXPECT_TRUE(cr.is_full());
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_FALSE(xc.active);
    HS_EXPECT_TRUE(cr.contains_x(0));
    HS_EXPECT_TRUE(cr.contains_x(W - 1));
    expect_xclip_parity(cr);
  }

  // 5) Over-wrap to full coverage: [2,90) expanded by 8 gives display width 88 +
  //    both margins = 104 >= w, so every column renders. render_x_start (90) !=
  //    render_x_end (2) would otherwise fold to an 8-column sliver and wrongly
  //    clip the display band's own interior; the full-coverage test forbids it.
  {
    ClipRegion cr = make(2, 90, 8);
    HS_EXPECT_EQ(cr.render_x_start(), 90);
    HS_EXPECT_EQ(cr.render_x_end(), 2);
    const ClipRegion::XClip xc = cr.x_clip();
    HS_EXPECT_FALSE(xc.active);
    HS_EXPECT_TRUE(cr.contains_x(2));  // display left edge
    HS_EXPECT_TRUE(cr.contains_x(50)); // display interior
    HS_EXPECT_TRUE(cr.contains_x(89)); // display right edge
    HS_EXPECT_TRUE(cr.contains_x(90)); // gap, but reached from both margins
    expect_xclip_parity(cr);
  }
}

// ============================================================================
// Plot::edge_row_span — arc-aware clip cull
// ============================================================================

/**
 * @brief Builds an orthonormal basis from a unit normal without a quaternion.
 * @param n Direction used as the basis normal; need not be pre-normalized.
 * @return Basis whose v is the normalized normal and whose u, w span the
 *         tangent plane.
 * @details Mirrors make_basis's construction; picks a reference axis that is not
 *          near-parallel to n to keep the cross products well-conditioned.
 */
inline Basis basis_from_normal(const Vector &n) {
  Vector v = n.normalized();
  Vector ref = std::abs(dot(v, X_AXIS)) > math::COS_AXIS_PARALLEL ? Y_AXIS : X_AXIS;
  Vector u = cross(v, ref).normalized();
  Vector w = cross(v, u).normalized();
  return {u, v, w};
}

/**
 * @brief Verifies Plot::edge_row_span conservatively covers the rendered arc's
 *        screen-row extent, including the interior latitude bulge where the arc
 *        reaches rows beyond both endpoints.
 * @details Densely samples the true arc for both the geodesic and planar
 *          strategies, asserts the span contains it, and confirms the randomized
 *          sweep produces many genuine bulge cases so the check is not vacuous.
 */
inline void test_edge_row_span_covers_arc_bulge() {
  constexpr int TW = 288, TH = 144;
  auto row_of = [](const Vector &v) {
    return vector_to_pixel<TW, TH>(v.normalized()).y;
  };

  hs::random().seed(20260609);
  int bulge_cases = 0; // edges whose arc bulges past both endpoints

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

    // Dense ground truth for the rendered arc's row extent.
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

    // Span must conservatively contain the arc (sub-pixel tolerance for fast-math).
    HS_EXPECT_LE(n_lo, t_lo + 0.25f);
    HS_EXPECT_GE(n_hi, t_hi - 0.25f);

    // Count edges whose interior bulges past the endpoints.
    float e_lo = std::min(row_of(a), row_of(b));
    float e_hi = std::max(row_of(a), row_of(b));
    if (t_lo < e_lo - 1.0f || t_hi > e_hi + 1.0f) bulge_cases++;
  }

  // Non-vacuity guard: the sweep must produce many genuine bulge cases.
  HS_EXPECT_GT(bulge_cases, 500);

  // Exact-antipodal geodesic edges: cross(a, b) collapses, so the renderer
  // slerps the semicircle about stable_perpendicular_axis. Ground truth is
  // built about that same axis; an endpoint-only span would cull the arc.
  for (int trial = 0; trial < 500; ++trial) {
    Vector ra(hs::rand_f(-1, 1), hs::rand_f(-1, 1), hs::rand_f(-1, 1));
    if (ra.length() < 0.1f) continue;
    Vector a = ra.normalized();
    Vector b = a * -1.0f;
    Vector axis = Plot::stable_perpendicular_axis(a);
    Vector vperp = cross(axis, a);

    float t_lo = 1e9f, t_hi = -1e9f;
    constexpr int N = 2000;
    for (int i = 0; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      Vector p = a * cosf(PI_F * t) + vperp * sinf(PI_F * t);
      float r = row_of(p);
      t_lo = std::min(t_lo, r);
      t_hi = std::max(t_hi, r);
    }

    float n_lo, n_hi;
    Plot::edge_row_span<TW, TH>(a, b, nullptr, n_lo, n_hi);
    HS_EXPECT_LE(n_lo, t_lo + 0.25f);
    HS_EXPECT_GE(n_hi, t_hi - 0.25f);
  }
}

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
  constexpr float base_step = (2.0f * PI_F) / W;
  const float KX = W / (2.0f * PI_F);
  const float KY = (H_VIRT - 1) / PI_F;

  auto screen_speed = [&](const Vector &pos, const Vector &tan) {
    const float ds = 1e-4f;
    auto xy = [&](const Vector &p) {
      float phi = std::acos(hs::clamp(p.y, -1.0f, 1.0f));
      float lam = std::atan2(p.z, p.x);
      return std::pair<float, float>(lam * KX, phi * KY);
    };
    Vector pp = (pos * std::cos(ds) + tan * std::sin(ds)).normalized();
    Vector pm = (pos * std::cos(ds) - tan * std::sin(ds)).normalized();
    auto a = xy(pp);
    auto b = xy(pm);
    float dx = (a.first - b.first) / (2.0f * ds);
    float dy = (a.second - b.second) / (2.0f * ds);
    return std::sqrt(dx * dx + dy * dy);
  };

  struct Case {
    Vector pos, tan;
  };
  // The tangent must be perpendicular to pos (a genuine arc-length tangent), or
  // the geodesic oracle and screen_step's formula parametrize differently.
  const Vector pos_off(std::sqrt(0.75f), 0.5f, 0.0f);
  const Vector raw(0.0f, 1.0f, 1.0f); // arbitrary; projected onto the tangent plane
  const Vector tan_mixed = (raw - pos_off * dot(raw, pos_off)).normalized();
  // Equatorial longitudinal, off-equator longitudinal, and a mixed (colatitude +
  // longitude) tangent — all verified below to land inside the clamp window.
  const Case cases[] = {
      {Vector(1.0f, 0.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f)},
      {pos_off, Vector(0.0f, 0.0f, 1.0f)},
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

// ============================================================================
// Plot::Ring::sample
// ============================================================================

/**
 * @brief Verifies Ring::sample emits N unit-length fragments plus one closing
 *        overlap fragment, with v0 progress running 0..1 and the closing
 *        fragment coinciding with sample 0.
 */
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
  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points[points.size() - 1].v0, 1.0f, 1e-6f);

  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);
}

/**
 * @brief Verifies Ring::sample<W,H> built from the TrigLUT angle-addition
 *        identity matches a direct cos/sin(theta+phase) construction of the same
 *        ring, in both position and the analytic arc-length register.
 * @details Covers Ring::draw's full-resolution path, which uses the LUT identity
 *          instead of a per-sample libm cos/sin.
 */
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
    HS_EXPECT_NEAR(points[i].v1, theta * r_val, 2e-3f);
  }

  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.z, points[0].pos.z, 1e-3f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-6f);

  // Close vertex carries the full-perimeter arc length: the ring's true
  // great-circle circumference 2*pi*sin(theta_eq).
  HS_EXPECT_NEAR(points.back().v1, 2.0f * PI_F * r_val, 2e-3f);
}

// ============================================================================
// Plot::DistortedRing::sample  — angle-addition identity (LUT) vs direct
// ============================================================================

/**
 * @brief Verifies DistortedRing::sample with a zero shift function reduces to a
 *        plain ring built via the TrigLUT angle-addition identity, matching a
 *        direct cos/sin construction of the same ring per fragment.
 */
inline void test_distorted_ring_sample_angle_addition_identity() {
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

/**
 * @brief Verifies DistortedRing::sample_arc emits vertices bit-identical
 *        (position, v0, v2) to the closed sample() at the same indices, up to
 *        and including the i == W overlap vertex, with v1 accumulating from
 *        the arc start.
 * @details This identity is what lets a clip-culling caller decompose a ring
 *          into open arcs (splitting at the seam) without changing any pixel
 *          the closed draw produces inside the clip.
 */
inline void test_distorted_ring_arc_matches_closed() {
  constexpr int W = 64;
  constexpr int H = 64;

  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0.3f, 0.8f, 0.5f));
  const float phase = 0.4f;
  ScalarFn shift = [](float t) { return 0.2f * sinf(2.0f * PI_F * 3.0f * t); };

  // radius < 1 and > 1 cover both sides of the antipode remap.
  for (float radius : {0.7f, 1.4f}) {
    ScratchScope sc(plot_arena());
    Fragments closed;
    closed.bind(plot_arena(), W + 2);
    Plot::DistortedRing::sample<W, H>(closed, b, radius, shift, phase);

    // Interior arc and an arc ending on the seam's overlap vertex.
    const int ranges[2][2] = {{5, 20}, {W - 8, W}};
    for (auto &r : ranges) {
      Fragments arc;
      arc.bind(plot_arena(), r[1] - r[0] + 2);
      Plot::DistortedRing::sample_arc<W, H>(arc, b, radius, shift, r[0], r[1],
                                            phase);
      HS_EXPECT_EQ(arc.size(), (size_t)(r[1] - r[0] + 1));
      HS_EXPECT_EQ(arc[0].v1, 0.0f);
      for (int k = 0; k < (int)arc.size(); ++k) {
        int i = r[0] + k;
        HS_EXPECT_EQ(arc[k].pos.x, closed[i].pos.x);
        HS_EXPECT_EQ(arc[k].pos.y, closed[i].pos.y);
        HS_EXPECT_EQ(arc[k].pos.z, closed[i].pos.z);
        HS_EXPECT_EQ(arc[k].v0, closed[i].v0);
        HS_EXPECT_EQ(arc[k].v2, closed[i].v2);
      }
      // v1 re-accumulates the same great-circle steps from the arc's start.
      float expect_len = 0.0f;
      for (int k = 1; k < (int)arc.size(); ++k) {
        expect_len += angle_between(arc[k - 1].pos, arc[k].pos);
        HS_EXPECT_NEAR(arc[k].v1, expect_len, 1e-5f);
      }
    }
  }
}

/**
 * @brief Verifies DistortedRing::draw_culled renders pixels bit-identical to
 *        the closed draw inside a partial clip, requests LUT bakes only
 *        through bake_run, and skips whole a ring that cannot reach the clip.
 * @details The shift LUT starts at an out-of-reach sentinel, so any draw that
 *          samples a range bake_run was never given displaces the ring far off
 *          its band and breaks the pixel identity.
 */
inline void test_distorted_ring_draw_culled_matches_closed() {
  constexpr int W = 96, H = 48;
  const float pad = 3.0f * PI_F / H;
  const float MAX_SHIFT = 0.15f;
  const float SENTINEL = 0.9f; // unbaked reads displace far outside reach
  const float reach = MAX_SHIFT + pad;

  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0.3f, 0.8f, 0.5f));

  float lut[W + 1];
  auto reset_lut = [&] {
    for (int i = 0; i <= W; ++i)
      lut[i] = SENTINEL;
  };
  int bake_calls = 0;
  auto bake_run = [&](int i0, int i1) {
    ++bake_calls;
    // -1/+1 pads cover a parameter rounding just below i0 and the lerp
    // neighbor above i1; modulo folds a seam-ending run onto the wrap
    // vertices, mirroring a LUT-baking caller.
    for (int m = i0 - 1; m <= i1 + 1; ++m) {
      int idx = ((m % W) + W) % W;
      lut[idx] = MAX_SHIFT * sinf(2.0f * PI_F * 3.0f * idx / W);
    }
    lut[W] = lut[0];
  };
  ScalarFn shift = [&](float t) {
    float x = wrap_t(t) * W;
    int j = static_cast<int>(x);
    return lut[j] + (x - j) * (lut[j + 1] - lut[j]);
  };
  auto shade = [](const Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.8f);
  };

  for (float radius : {0.7f, 1.4f}) {
    // Reference: full clip routes through the closed-draw path.
    std::vector<Pixel> ref(static_cast<size_t>(W) * H);
    {
      RasterFx fx(W, H);
      reset_lut();
      Pipeline<W, H> filters;
      {
        Canvas c(fx);
        Plot::DistortedRing::draw_culled<W, H>(filters, c, fx.clip(), b,
                                               radius, shift, shade, reach,
                                               bake_run);
      }
      fx.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          ref[static_cast<size_t>(y) * W + x] = fx.get_pixel(x, y);
    }

    // A row band and a column wedge: interior pixels must match exactly.
    const int clips[2][4] = {{H / 4, H / 2, 0, W}, {0, H, 10, 34}};
    for (auto &cl : clips) {
      RasterFx fx(W, H);
      fx.set_clip(cl[0], cl[1], cl[2], cl[3]);
      reset_lut();
      bake_calls = 0;
      Pipeline<W, H> filters;
      {
        Canvas c(fx);
        Plot::DistortedRing::draw_culled<W, H>(filters, c, fx.clip(), b,
                                               radius, shift, shade, reach,
                                               bake_run);
      }
      fx.advance_display();
      int lit = 0, diff = 0;
      for (int y = cl[0]; y < cl[1]; ++y)
        for (int x = cl[2]; x < cl[3]; ++x) {
          Pixel p = fx.get_pixel(x, y);
          const Pixel &r = ref[static_cast<size_t>(y) * W + x];
          if (r.r | r.g | r.b)
            ++lit;
          if (p.r != r.r || p.g != r.g || p.b != r.b)
            ++diff;
        }
      HS_EXPECT_GT(lit, 0);
      HS_EXPECT_EQ(diff, 0);
      HS_EXPECT_GT(bake_calls, 0);
    }
  }

  // A ring whose widened cap cannot reach the clip is skipped whole: no bake.
  {
    RasterFx fx(W, H);
    fx.set_clip(H - H / 8, H, 0, W);
    reset_lut();
    bake_calls = 0;
    Pipeline<W, H> filters;
    {
      Canvas c(fx);
      Basis axis_y = make_basis(Quaternion(1, 0, 0, 0), Vector(0, 1, 0));
      Plot::DistortedRing::draw_culled<W, H>(filters, c, fx.clip(), axis_y,
                                             0.3f, shift, shade, reach,
                                             bake_run);
    }
    fx.advance_display();
    HS_EXPECT_EQ(bake_calls, 0);
  }
}

/**
 * @brief Verifies draw_culled's run decomposition and its remaining knobs: a
 *        wedge the ring crosses twice yields multiple bake runs that are
 *        ascending, non-overlapping, and chunk-aligned; a nonzero phase and an
 *        alternate CHUNKS instantiation both preserve the in-clip pixel
 *        identity with the closed draw.
 */
inline void test_distorted_ring_draw_culled_runs_phase_chunks() {
  constexpr int W = 96, H = 48;
  const float pad = 3.0f * PI_F / H;
  const float MAX_SHIFT = 0.1f;
  const float reach = MAX_SHIFT + pad;
  const float phase = 0.6f;
  const float radius = 1.0f;

  // Near-equatorial axis: the radius-1 ring is a near-polar great circle that
  // enters and leaves a mid-canvas column wedge in two separate azimuth spans.
  Basis b = make_basis(Quaternion(1, 0, 0, 0), Vector(0.9f, 0.1f, 0.42f));

  float lut[W + 1];
  auto reset_lut = [&] {
    for (int i = 0; i <= W; ++i)
      lut[i] = 0.9f;
  };
  std::vector<std::pair<int, int>> runs;
  auto bake_run = [&](int i0, int i1) {
    runs.push_back({i0, i1});
    for (int m = i0 - 1; m <= i1 + 1; ++m) {
      int idx = ((m % W) + W) % W;
      lut[idx] = MAX_SHIFT * sinf(2.0f * PI_F * 2.0f * idx / W);
    }
    lut[W] = lut[0];
  };
  ScalarFn shift = [&](float t) {
    float x = wrap_t(t) * W;
    int j = static_cast<int>(x);
    return lut[j] + (x - j) * (lut[j + 1] - lut[j]);
  };
  auto shade = [](const Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.8f);
  };

  // Reference: full clip, same phase.
  std::vector<Pixel> ref(static_cast<size_t>(W) * H);
  {
    RasterFx fx(W, H);
    reset_lut();
    Pipeline<W, H> filters;
    {
      Canvas c(fx);
      Plot::DistortedRing::draw_culled<W, H>(filters, c, fx.clip(), b, radius,
                                             shift, shade, reach, bake_run,
                                             phase);
    }
    fx.advance_display();
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        ref[static_cast<size_t>(y) * W + x] = fx.get_pixel(x, y);
  }

  const int cx0 = 40, cx1 = 56;
  for (int chunks_case = 0; chunks_case < 2; ++chunks_case) {
    const int V = chunks_case == 0 ? W / 24 : W / 8;
    RasterFx fx(W, H);
    fx.set_clip(0, H, cx0, cx1);
    reset_lut();
    runs.clear();
    Pipeline<W, H> filters;
    {
      Canvas c(fx);
      if (chunks_case == 0)
        Plot::DistortedRing::draw_culled<W, H, 24>(filters, c, fx.clip(), b,
                                                   radius, shift, shade, reach,
                                                   bake_run, phase);
      else
        Plot::DistortedRing::draw_culled<W, H, 8>(filters, c, fx.clip(), b,
                                                  radius, shift, shade, reach,
                                                  bake_run, phase);
    }
    fx.advance_display();

    // Partial cull engaged (not the full-mask path), and at 24 chunks the two
    // wedge crossings decompose into separate runs.
    HS_EXPECT_GT(runs.size(), (size_t)0);
    if (chunks_case == 0)
      HS_EXPECT_GT(runs.size(), (size_t)1);
    HS_EXPECT_TRUE(!(runs.size() == 1 && runs[0].first == 0 &&
                     runs[0].second == W));
    for (size_t k = 0; k < runs.size(); ++k) {
      HS_EXPECT_TRUE(runs[k].first >= 0 && runs[k].second <= W &&
                     runs[k].first < runs[k].second);
      HS_EXPECT_EQ(runs[k].first % V, 0);
      HS_EXPECT_EQ(runs[k].second % V, 0);
      if (k > 0)
        HS_EXPECT_TRUE(runs[k - 1].second < runs[k].first);
    }

    int lit = 0, diff = 0;
    for (int y = 0; y < H; ++y)
      for (int x = cx0; x < cx1; ++x) {
        Pixel p = fx.get_pixel(x, y);
        const Pixel &r = ref[static_cast<size_t>(y) * W + x];
        if (r.r | r.g | r.b)
          ++lit;
        if (p.r != r.r || p.g != r.g || p.b != r.b)
          ++diff;
      }
    HS_EXPECT_GT(lit, 0);
    HS_EXPECT_EQ(diff, 0);
  }
}

// ============================================================================
// Plot::Spiral::sample
// ============================================================================

/**
 * @brief Verifies Spiral::sample emits N unit-length fragments with
 *        non-decreasing cumulative arc length (v1) and v0 progress spanning 0..1,
 *        and that the result is genuinely a pole-to-pole Fibonacci spiral.
 * @details Shape-discriminating check: the latitude (y) sweeps strictly
 *          monotonically from pole to pole while the azimuth winds many full
 *          turns — neither holds for a ring (constant latitude, one turn) or a
 *          meridian arc (no azimuthal winding), so this distinguishes the spiral
 *          from the other closed/curve shapes rather than merely confirming the
 *          points lie on the sphere.
 */
inline void test_spiral_sample_unit_length_and_monotone_arc() {
  ScratchScope sc(plot_arena());
  Fragments frags;
  const int N = 48;
  frags.bind(plot_arena(), N);

  Plot::Spiral::sample(frags, N, 0.5f);
  HS_EXPECT_EQ(frags.size(), (size_t)N);

  float last_v1 = -1.0f;
  float last_y = 2.0f;        // above any unit-sphere y, so the first compare passes
  float last_theta = 0.0f;
  float winding = 0.0f;       // total absolute azimuthal travel (radians)
  for (size_t i = 0; i < frags.size(); ++i) {
    const Vector &p = frags[i].pos;
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
    HS_EXPECT_GE(frags[i].v1, last_v1);
    last_v1 = frags[i].v1;

    HS_EXPECT_LT(p.y, last_y);
    last_y = p.y;

    float theta = std::atan2(p.z, p.x);
    if (i > 0) {
      float d = theta - last_theta;
      while (d > PI_F) d -= 2.0f * PI_F;
      while (d < -PI_F) d += 2.0f * PI_F;
      winding += std::fabs(d);
    }
    last_theta = theta;
  }
  HS_EXPECT_NEAR(frags[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(frags[N - 1].v0, 1.0f, 1e-6f);

  // Many azimuthal turns — a single ring would wind only ~2π, an arc ~0.
  HS_EXPECT_GT(winding, 4.0f * PI_F);
}

// ============================================================================
// Plot::Multiline::sample
// ============================================================================

/**
 * @brief Verifies Multiline::sample arc-length-parameterizes the polyline: v0
 *        progress is monotone 0..1, v2 carries the vertex index, and v1
 *        cumulative arc length sums the per-edge angles.
 */
inline void test_multiline_sample_arclength_param() {
  ScratchScope sc(plot_arena());

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

  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-4f);
  float last = -1.0f;
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_GE(points[i].v0, last - 1e-6f);
    last = points[i].v0;
    HS_EXPECT_NEAR(points[i].v2, static_cast<float>(i), 1e-6f);
  }
  // v1 cumulative arc length sums the two equal 90deg hops.
  HS_EXPECT_NEAR(points.back().v1, PI_F, 1e-3f);
}

// ============================================================================
// Plot::Star::sample / Plot::Flower::sample
// ============================================================================

/**
 * @brief Verifies Star::sample emits 2*sides unit-length vertices plus one
 *        closing fragment matching the first vertex (closed loop, v0 == 1 at
 *        close), and that the vertices form a genuine star — not a flower.
 * @details Shape-discriminating check: the colatitude (angle from the shape's
 *          center axis) alternates between an outer and an inner radius, with
 *          inner/outer == STAR_INNER_RATIO. The Flower test below pins the
 *          opposite invariant (constant colatitude), so the pair would catch a
 *          regression that swapped the two shapes or dropped the star's notches.
 */
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
  HS_EXPECT_NEAR(points.back().pos.x, points[0].pos.x, 1e-3f);
  HS_EXPECT_NEAR(points.back().pos.y, points[0].pos.y, 1e-3f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-6f);

  // Alternating outer/inner colatitude about the center axis.
  const Vector axis = get_antipode(b, 0.5f).first.v;
  const float outer = angle_between(points[0].pos, axis); // even index -> outer
  const float inner = angle_between(points[1].pos, axis); // odd index  -> inner
  HS_EXPECT_GT(outer, inner + 1e-3f);
  HS_EXPECT_NEAR(inner / outer, Plot::STAR_INNER_RATIO, 1e-2f);
  for (int i = 0; i < sides * 2; ++i) {
    const float colat = angle_between(points[i].pos, axis);
    HS_EXPECT_NEAR(colat, (i % 2 == 0) ? outer : inner, 1e-3f);
  }
}

/**
 * @brief Verifies Flower::sample emits 2*sides unit-length vertices plus one
 *        closing fragment matching the first vertex (closed loop), and that the
 *        vertices form a genuine flower — constant radius, not a star.
 * @details Shape-discriminating check: every vertex sits at the SAME colatitude
 *          about the center axis (a constant polar radius), the opposite of the
 *          star's alternating outer/inner notches. The pair pins the two shapes
 *          apart at the sample level.
 */
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

  // Constant colatitude about the center axis.
  const Vector axis = get_antipode(b, 0.5f).first.v;
  const float colat0 = angle_between(points[0].pos, axis);
  for (int i = 0; i < sides * 2; ++i) {
    const float colat = angle_between(points[i].pos, axis);
    HS_EXPECT_NEAR(colat, colat0, 1e-3f);
  }
}

// ============================================================================
// Plot::rasterize — control-flow coverage through a capturing pipeline
// ============================================================================

/**
 * @brief A sub-base_step open segment takes the fast path and plots BOTH
 *        endpoints (start, then the final vertex since the loop is open).
 */
inline void test_rasterize_subpixel_open_segment_plots_both_endpoints() {
  constexpr int W = 128, H = 64;
  RasterFx fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  // ~0.02 rad apart, well under base_step (2*pi/W = ~0.049 rad).
  b.pos = Vector(1.0f, 0.02f, 0.0f).normalized();
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader, /*close_loop=*/false);
  }
  fx.advance_display();

  // Fast path on an open last segment plots curr and next.
  HS_EXPECT_EQ(pipe.plotted.size(), (size_t)2);
  HS_EXPECT_NEAR(angle_between(pipe.plotted.front(), a.pos), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(angle_between(pipe.plotted.back(), b.pos), 0.0f, 1e-4f);
  for (const Vector &p : pipe.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
}

/**
 * @brief A normal-length open segment is sampled densely enough to be gap-free
 *        (no inter-sample gap exceeds one pixel column) and lands on both
 *        endpoints.
 */
inline void test_rasterize_open_segment_gap_free() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * PI_F) / W;
  RasterFx fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  b.pos = Vector(0, 1, 0);
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader, /*close_loop=*/false);
  }
  fx.advance_display();

  HS_EXPECT_GT(pipe.plotted.size(), (size_t)10);
  // Adaptive sub-stepping caps each advance at ~base_step (slack for quantization).
  HS_EXPECT_LE(max_consecutive_gap(pipe.plotted, /*wrap=*/false),
               1.5f * base_step);
  HS_EXPECT_NEAR(angle_between(pipe.plotted.front(), a.pos), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(angle_between(pipe.plotted.back(), b.pos), 0.0f, 1e-3f);
  for (const Vector &p : pipe.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
}

/**
 * @brief A closed loop omits each segment's terminal vertex, so the plotted
 *        sequence wraps continuously (gap-free across the seam) with no shared
 *        vertex plotted twice.
 */
inline void test_rasterize_closed_loop_gap_free_no_dup() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * PI_F) / W;
  RasterFx fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // A spherical triangle with well-separated vertices.
  Fragment v0, v1, v2;
  v0.pos = Vector(1, 0, 0);
  v1.pos = Vector(0, 1, 0);
  v2.pos = Vector(0, 0, 1);
  points.push_back(v0);
  points.push_back(v1);
  points.push_back(v2);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader, /*close_loop=*/true);
  }
  fx.advance_display();

  HS_EXPECT_GT(pipe.plotted.size(), (size_t)20);
  // Continuity all the way around, including the last->first seam.
  HS_EXPECT_LE(max_consecutive_gap(pipe.plotted, /*wrap=*/true),
               1.5f * base_step);
  // No vertex plotted twice: consecutive samples stay distinct.
  for (size_t i = 1; i < pipe.plotted.size(); ++i)
    HS_EXPECT_GT(angle_between(pipe.plotted[i - 1], pipe.plotted[i]), 1e-5f);
}

/**
 * @brief A planar segment whose endpoint sits at the basis antipode falls back
 *        to the geodesic strategy: rasterizing with that planar_basis must
 *        produce exactly the same samples as rasterizing geodesically.
 */
inline void test_rasterize_antipodal_seam_planar_falls_back_geodesic() {
  constexpr int W = 128, H = 64;
  RasterFx fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // planar basis centered on +Z; the second endpoint sits within the antipodal
  // margin of -Z (dot < -COS_PLANAR_ANTIPODE), tripping the seam guard.
  Basis basis = basis_from_normal(Vector(0, 0, 1));
  Fragment a, b;
  a.pos = Vector(1, 0, 0);
  b.pos = Vector(0.02f, 0.0f, -1.0f).normalized();
  HS_EXPECT_LT(dot(b.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  points.push_back(a);
  points.push_back(b);

  CapturePipeline planar_pipe, geo_pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(planar_pipe, c, points, noop_shader,
                          /*close_loop=*/false, &basis);
  }
  fx.advance_display();
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(geo_pipe, c, points, noop_shader,
                          /*close_loop=*/false, /*planar_basis=*/nullptr);
  }
  fx.advance_display();

  HS_EXPECT_EQ(planar_pipe.plotted.size(), geo_pipe.plotted.size());
  size_t n = std::min(planar_pipe.plotted.size(), geo_pipe.plotted.size());
  for (size_t i = 0; i < n; ++i)
    HS_EXPECT_NEAR(angle_between(planar_pipe.plotted[i], geo_pipe.plotted[i]),
                   0.0f, 1e-5f);
}

/**
 * @brief A non-seam planar segment renders gap-free in ARC length: the
 *        arc-length parameterization (map_planar's cumulative-arc inversion)
 *        keeps every plotted step near one pixel column and lands on both
 *        endpoints, with no clustering or gaps the projection-linear chord would
 *        otherwise leave.
 * @details Exercises the planar strategy path (rasterize_planar_strategy +
 *          map_planar), which the other rasterize tests do not reach — the
 *          antipodal-seam case falls back to geodesic. This locks in the
 *          end-to-end arc-uniform sampling the LEN_SAMPLES table provides; it
 *          does not isolate the table's contribution from the rasterizer's
 *          adaptive (sin-phi) sub-stepping, which also shapes local density.
 */
inline void test_rasterize_planar_segment_gap_free_arclength() {
  constexpr int W = 128, H = 64;
  constexpr float base_step = (2.0f * PI_F) / W;
  RasterFx fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // Planar disk about +Y; endpoints sweep colatitude 0.3 -> 1.3 across azimuths
  // so the chord crosses regions of differing azimuthal stretch (r / sin r).
  Basis basis = basis_from_normal(Vector(0, 1, 0));
  auto on_disk = [&](float colat, float az) {
    Vector dir = basis.u * cosf(az) + basis.w * sinf(az);
    return (basis.v * cosf(colat) + dir * sinf(colat)).normalized();
  };
  Fragment a, b;
  a.pos = on_disk(0.3f, 0.0f);
  b.pos = on_disk(1.3f, 1.0f);
  // Stays out of the antipodal-seam fallback so the planar strategy is used.
  HS_EXPECT_GT(dot(a.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  HS_EXPECT_GT(dot(b.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, noop_shader, /*close_loop=*/false,
                          &basis);
  }
  fx.advance_display();

  HS_EXPECT_GT(pipe.plotted.size(), (size_t)10);
  HS_EXPECT_LE(max_consecutive_gap(pipe.plotted, /*wrap=*/false),
               1.5f * base_step);
  // Endpoints land within map_planar's project/unproject round-trip error.
  HS_EXPECT_NEAR(angle_between(pipe.plotted.front(), a.pos), 0.0f, 1e-2f);
  HS_EXPECT_NEAR(angle_between(pipe.plotted.back(), b.pos), 0.0f, 1e-2f);
  for (const Vector &p : pipe.plotted)
    HS_EXPECT_NEAR(p.length(), 1.0f, 1e-3f);
}

/**
 * @brief A planar segment's v0/v1 registers measure the RENDERED azimuthal arc,
 *        not the geodesic chord: v1 rises monotonically to the planar arc length
 *        (strictly longer than the endpoints' great-circle distance) and v0 spans
 *        0..1 over that same arc. Locks in the rasterizer's arc-register
 *        re-derivation for the always-planar primitives (PlanarPolygon / Star /
 *        Flower), which the sample()-level tests cannot observe.
 */
inline void test_rasterize_planar_arc_registers_track_drawn_arc() {
  constexpr int W = 128, H = 64;
  RasterFx fx(W, H);
  ScratchScope sc(plot_arena());
  Fragments points;
  points.bind(plot_arena(), 4);

  // Same non-seam planar disk edge as the gap-free test: colatitude 0.3 -> 1.3
  // across azimuths, so the rendered edge bows well clear of its chord.
  Basis basis = basis_from_normal(Vector(0, 1, 0));
  auto on_disk = [&](float colat, float az) {
    Vector dir = basis.u * cosf(az) + basis.w * sinf(az);
    return (basis.v * cosf(colat) + dir * sinf(colat)).normalized();
  };
  Fragment a, b;
  a.pos = on_disk(0.3f, 0.0f);
  b.pos = on_disk(1.3f, 1.0f);
  HS_EXPECT_GT(dot(a.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  HS_EXPECT_GT(dot(b.pos, basis.v), -Plot::COS_PLANAR_ANTIPODE);
  // Bare control points default v0/v1 to 0, so any nonzero arc below comes
  // solely from the rasterizer's rendered-arc override.
  points.push_back(a);
  points.push_back(b);

  CapturePipeline pipe;
  std::vector<float> v0s, v1s;
  auto capture = [&](const Vector &, Fragment &f) {
    v0s.push_back(f.v0);
    v1s.push_back(f.v1);
  };
  {
    Canvas c(fx);
    Plot::rasterize<W, H>(pipe, c, points, capture, /*close_loop=*/false, &basis);
  }
  fx.advance_display();

  HS_EXPECT_GT(v1s.size(), (size_t)10);

  for (size_t i = 1; i < v1s.size(); ++i) {
    HS_EXPECT_GE(v1s[i], v1s[i - 1] - 1e-6f);
    HS_EXPECT_GE(v0s[i], v0s[i - 1] - 1e-6f);
  }

  const float geo = angle_between(a.pos, b.pos);
  const float planar = Plot::planar_arc_length(a.pos, b.pos, basis);
  // The rendered planar edge bows strictly longer than the geodesic chord...
  HS_EXPECT_GT(planar, geo);
  // ...and v1 tracks that rendered arc (start ~0, end ~planar length, > chord).
  HS_EXPECT_NEAR(v1s.front(), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(v1s.back(), planar, 2e-2f);
  HS_EXPECT_GT(v1s.back(), geo);

  // v0 is v1 normalized by the single-segment total arc: 0 at the start, ~1 end.
  HS_EXPECT_NEAR(v0s.front(), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(v0s.back(), 1.0f, 2e-2f);
}

// ============================================================================
// Plot::ParticleSystem — trail rasterization
// ============================================================================

/**
 * @brief Angular distance from a unit point to a geodesic ARC a→b (not the full
 *        great circle).
 * @param p Query point (unit length).
 * @param a Arc start (unit length).
 * @param b Arc end (unit length).
 * @return Perpendicular distance to the arc when the foot of the perpendicular
 *         falls within it, else the nearer endpoint distance.
 * @details Arc-clamped (not plane-clamped) so a point bulging off its own segment
 *          is not scored as "near" a different segment whose infinite great circle
 *          it happens to pass close to.
 */
inline float dist_to_arc(const Vector &p, const Vector &a, const Vector &b) {
  Vector n = cross(a, b);
  float nlen = n.length();
  if (nlen < 1e-6f) // degenerate (coincident/antipodal) arc
    return angle_between(p, a);
  n = n * (1.0f / nlen);
  Vector foot = p - n * dot(p, n); // project onto the great-circle plane
  float flen = foot.length();
  if (flen > 1e-6f) {
    foot = foot * (1.0f / flen);
    float ab = angle_between(a, b);
    if (angle_between(a, foot) <= ab + 1e-4f &&
        angle_between(b, foot) <= ab + 1e-4f)
      return std::fabs(std::asin(hs::clamp(dot(p, n), -1.0f, 1.0f)));
  }
  return std::min(angle_between(p, a), angle_between(p, b));
}

/** @brief Minimal particle for the ParticleSystem draw concept: trail + life. */
struct StubParticle {
  Animation::VectorTrail<16> history; /**< World-space trail positions. */
  uint16_t life = 0;                  /**< Remaining life (drives v3). */
};

/** @brief Minimal pool/active_count/max_life triple the draw concept reads. */
struct StubSystem {
  std::vector<StubParticle> pool; /**< Pool; pool[i] read for i < active_count. */
  int active_count = 0;           /**< Live prefix length. */
  uint16_t max_life = 0;          /**< Life normaliser for v3. */
  int active() const { return active_count; }
};

/**
 * @brief Verifies ParticleSystem::draw rasterizes only the active prefix's trails
 *        and stamps the per-particle registers (v2 source index, v3 life ratio).
 * @details Pins three things the smoke loop never checks: (1) only particles in
 * [0, active_count) are drawn — an inactive particle parked on the ±Y poles never
 * contributes a point; (2) the drawn trail follows the particle's recorded
 * history (an equatorial +X→+Z arc); (3) every emitted fragment carries the
 * source index in v2 and life/max_life in v3, constant across the trail (the
 * register convention MindSplatter's invariant assert depends on).
 */
inline void test_particle_system_draws_active_trails_with_registers() {
  constexpr int W = 288, H = 144;
  RasterFx fx(W, H);

  StubSystem sys;
  sys.max_life = 100;
  sys.active_count = 1; // pool holds 2; only particle 0 is active

  StubParticle p0;
  p0.life = 60;
  const Vector t0[3] = {Vector(1, 0, 0), Vector(0.7071f, 0.0f, 0.7071f),
                        Vector(0, 0, 1)}; // equatorial +X -> +Z arc
  for (const Vector &v : t0)
    p0.history.record(v);

  StubParticle p1; // inactive: parked on the poles, must NOT be drawn
  p1.life = 99;
  p1.history.record(Vector(0, 1, 0));
  p1.history.record(Vector(0, -1, 0));

  sys.pool.push_back(p0);
  sys.pool.push_back(p1);

  CapturePipeline pipe;
  float v2_lo = 1e9f, v2_hi = -1e9f, v3_lo = 1e9f, v3_hi = -1e9f;
  {
    Canvas c(fx);
    Plot::ParticleSystem::draw<W, H>(
        pipe, c, sys, [&](const Vector &, Fragment &f) {
          v2_lo = std::min(v2_lo, f.v2);
          v2_hi = std::max(v2_hi, f.v2);
          v3_lo = std::min(v3_lo, f.v3);
          v3_hi = std::max(v3_hi, f.v3);
        });
  }
  fx.advance_display();

  // (1)+(2) Active trail follows its recorded arc; the inactive ±Y particle is absent.
  HS_EXPECT_GT(pipe.plotted.size(), (size_t)2);
  for (const Vector &v : pipe.plotted) {
    HS_EXPECT_LE(dist_to_arc(v, t0[0], t0[2]), 0.05f);
    HS_EXPECT_LT(std::fabs(v.y), 0.1f);
  }
  // (3) Registers: v2 == source index 0; v3 == life/max_life == 0.6, constant.
  HS_EXPECT_NEAR(v2_lo, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(v2_hi, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(v3_lo, 0.6f, 1e-3f);
  HS_EXPECT_NEAR(v3_hi, 0.6f, 1e-3f);
}

/**
 * @brief Verifies an empty system needs no lifetime normalizer.
 */
inline void test_particle_system_empty_zero_lifetime_is_noop() {
  constexpr int W = 288, H = 144;
  RasterFx fx(W, H);
  StubSystem sys;
  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::ParticleSystem::draw<W, H>(pipe, c, sys, noop_shader);
  }
  fx.advance_display();
  HS_EXPECT_TRUE(pipe.plotted.empty());
}

/**
 * @brief The segment cull follows a filter-chain orientation: an edge the
 *        World::Orient stage rotates into a clip band is drawn, not culled.
 * @details When orientation lives in the filter chain (FlowField) the rasterizer
 *          must bound the edge by its RENDERED latitude, so the cull is routed
 *          through the pipeline. 90° about X maps the equatorial +X->+Z arc onto
 *          a polar meridian, so the rendered arc reaches the bottom band the
 *          source arc never touches; a band-clipped worker there must match the
 *          full render. Without the pipeline-routed cull the edge is bounded by
 *          its un-rotated (equatorial) latitude and dropped — the segmented-mode
 *          FlowField clipping bug (docs/segmented_stateful_effects_spec.md).
 */
inline void test_rasterize_cull_follows_filter_orientation() {
  constexpr int W = 128, H = 64;
  constexpr int BAND = H / 4; // bottom band [H-BAND, H) the rotated arc enters

  Orientation<> orientation(make_rotation(X_AXIS, PI_F / 2.0f));
  auto shade = [](const Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 1.0f);
  };

  auto band_lit = [&](int cy0, int cy1) -> int {
    RasterFx fx(W, H);
    fx.set_clip(cy0, cy1, 0, W);
    Pipeline<W, H, Filter::World::Orient<W>> filters{
        Filter::World::Orient<W>(orientation)};
    {
      ScratchScope sc(plot_arena());
      Fragments pts;
      pts.bind(plot_arena(), 4);
      Fragment a, b;
      a.pos = Vector(1, 0, 0); // equator (+X)
      b.pos = Vector(0, 0, 1); // equator (+Z); 90° about X sends it to -Y (pole)
      pts.push_back(a);
      pts.push_back(b);
      Canvas c(fx);
      Plot::rasterize<W, H>(filters, c, pts, shade, /*close_loop=*/false);
    }
    fx.advance_display();
    int lit = 0;
    for (int y = H - BAND; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        Pixel p = fx.get_pixel(x, y);
        if (p.r | p.g | p.b)
          ++lit;
      }
    return lit;
  };

  int full = band_lit(0, H);          // full canvas: rotated arc reaches the band
  int banded = band_lit(H - BAND, H); // worker clipped to that band
  HS_EXPECT_GT(full, 0);
  HS_EXPECT_EQ(banded, full); // routed cull keeps the edge; identical in the band
}

// ============================================================================
// Azimuthal-equidistant projection + dual-metric planar arc length
//
// These pin plot.h's dual-metric core directly against independent oracles (a
// libm great-circle reconstruction, a fine on-sphere quadrature) instead of
// against the primitives themselves, which the sample()/rasterize() tests reuse
// as their own ground truth. Radial displacement is isometric on the sphere;
// azimuthal displacement stretches by R/sin R, and the arc integrator must
// track that anisotropy.
// ============================================================================

/** @brief Random unit vector, rejecting near-zero draws. */
inline Vector az_rand_unit() {
  for (;;) {
    Vector r(hs::rand_f(-1, 1), hs::rand_f(-1, 1), hs::rand_f(-1, 1));
    if (r.length() > 0.1f) return r.normalized();
  }
}

/**
 * @brief Full-precision azimuthal unprojection (libm), an oracle independent of
 *        plot.h's LUT-based fast-trig path.
 */
inline Vector az_unproject_exact(float Px, float Py, const Basis &b) {
  float R = std::sqrt(Px * Px + Py * Py);
  if (R < math::EPS_GEOMETRIC)
    return b.v;
  float th = std::atan2(Py, Px);
  Vector axis = b.u * std::cos(th) + b.w * std::sin(th);
  return b.v * std::cos(R) + axis * std::sin(R);
}

/**
 * @brief Great-circle angle between two unit vectors, accurate for tiny angles.
 * @details atan2(|p x q|, p.q) stays well-conditioned near 0 where acos is flat;
 *          fast_acos (via angle_between) collapses sub-milliradian steps to zero,
 *          which would corrupt a fine-quadrature reference.
 */
inline float az_arc_exact(const Vector &p, const Vector &q) {
  return std::atan2(cross(p, q).length(), dot(p, q));
}

/**
 * @brief azimuthal_project's radius equals the great-circle angle from center.
 * @details The projection's defining property: |proj| == angle_between(p,
 *          center). Checked against an independent angle_between, so an axis or
 *          scale error in the projection would break it.
 */
inline void test_azimuthal_project_radius_is_geodesic_angle() {
  hs::random().seed(0xA21E);
  int mid = 0;
  for (int trial = 0; trial < 4000; ++trial) {
    Basis basis = basis_from_normal(az_rand_unit());
    Vector p = az_rand_unit();
    float geo = angle_between(p, basis.v);
    auto proj = Plot::azimuthal_project(p, basis);
    float r = std::hypot(proj.first, proj.second);
    HS_EXPECT_NEAR(r, geo, 5e-3f * (geo + 1.0f));
    if (geo > 0.3f && geo < PI_F - 0.3f) ++mid;
  }
  HS_EXPECT_GT(mid, 1000);
}

/**
 * @brief azimuthal_project and azimuthal_unproject invert each other.
 * @details plane->sphere->plane and sphere->plane->sphere both return the
 *          input, away from the antipodal band where the azimuth is unstable.
 */
inline void test_azimuthal_roundtrip_identity() {
  hs::random().seed(0xB33F);
  int fwd = 0, inv = 0;
  for (int trial = 0; trial < 4000; ++trial) {
    Basis basis = basis_from_normal(az_rand_unit());

    float R = hs::rand_f(0.05f, PI_F - 0.05f);
    float th = hs::rand_f(-PI_F, PI_F);
    float Px = R * std::cos(th), Py = R * std::sin(th);
    Vector s = Plot::azimuthal_unproject(Px, Py, basis);
    auto rp = Plot::azimuthal_project(s, basis);
    HS_EXPECT_NEAR(rp.first, Px, 2e-2f * (R + 1.0f));
    HS_EXPECT_NEAR(rp.second, Py, 2e-2f * (R + 1.0f));
    ++inv;

    Vector p = az_rand_unit();
    if (dot(p, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;
    auto proj = Plot::azimuthal_project(p, basis);
    Vector back = Plot::azimuthal_unproject(proj.first, proj.second, basis);
    HS_EXPECT_NEAR(angle_between(p, back), 0.0f, 1.5e-2f);
    ++fwd;
  }
  HS_EXPECT_GT(inv, 3000);
  HS_EXPECT_GT(fwd, 3000);
}

/**
 * @brief azimuthal_unproject lands on the great-circle point at (R, theta).
 * @details Oracle is an independent libm reconstruction
 *          v*cos(R) + (u*cos(th)+w*sin(th))*sin(R); a sign or axis swap in the
 *          fast-trig unprojection would diverge from it.
 */
inline void test_azimuthal_unproject_hits_great_circle_point() {
  hs::random().seed(0xC0DE);
  int n = 0;
  for (int trial = 0; trial < 4000; ++trial) {
    Basis basis = basis_from_normal(az_rand_unit());
    float R = hs::rand_f(0.02f, PI_F - 0.02f);
    float th = hs::rand_f(-PI_F, PI_F);
    Vector got =
        Plot::azimuthal_unproject(R * std::cos(th), R * std::sin(th), basis);
    Vector axis = basis.u * std::cos(th) + basis.w * std::sin(th);
    Vector want = basis.v * std::cos(R) + axis * std::sin(R);
    HS_EXPECT_NEAR(angle_between(got, want), 0.0f, 1e-2f);
    ++n;
  }
  HS_EXPECT_GT(n, 3900);
}

/**
 * @brief planar_arc_length matches a fine libm quadrature of the edge.
 * @details Compares the 4-panel table against a 2000-panel libm reference (a
 *          full-precision unprojection summed with a small-angle-robust arc) on
 *          the short polygon-edge regime the primitive is built for. Non-vacuity:
 *          most edges bow past their great-circle chord. Long edges sweeping near
 *          the chart center are out of the primitive's domain — 4 samples
 *          straddle the azimuth singularity — and are excluded, as real polygon
 *          edges are.
 */
inline void test_planar_arc_length_matches_fine_quadrature() {
  hs::random().seed(0xD41A);
  int bows = 0;
  float max_rel_err = 0.0f;
  for (int trial = 0; trial < 3000; ++trial) {
    Basis basis = basis_from_normal(az_rand_unit());
    float R1 = hs::rand_f(0.2f, 1.2f), R2 = hs::rand_f(0.2f, 1.2f);
    float t1 = hs::rand_f(-PI_F, PI_F), t2 = t1 + hs::rand_f(0.15f, 0.8f);
    Vector a =
        Plot::azimuthal_unproject(R1 * std::cos(t1), R1 * std::sin(t1), basis);
    Vector b =
        Plot::azimuthal_unproject(R2 * std::cos(t2), R2 * std::sin(t2), basis);
    if (dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
        dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;

    auto p1 = Plot::azimuthal_project(a, basis);
    auto p2 = Plot::azimuthal_project(b, basis);
    constexpr int N = 2000;
    float fine = 0.0f;
    Vector prev = az_unproject_exact(p1.first, p1.second, basis);
    for (int i = 1; i <= N; ++i) {
      float t = static_cast<float>(i) / N;
      Vector cur = az_unproject_exact(p1.first + (p2.first - p1.first) * t,
                                      p1.second + (p2.second - p1.second) * t,
                                      basis);
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
 *          azimuthal stretch R/sin R rises. This separates the two metrics
 *          directly, which the end-to-end tests cannot.
 */
inline void test_dual_metric_radial_vs_azimuthal() {
  hs::random().seed(0xE1A5);
  int radial = 0, azi = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    Basis basis = basis_from_normal(az_rand_unit());

    float th = hs::rand_f(-PI_F, PI_F);
    float Ra = hs::rand_f(0.1f, 0.6f), Rb = hs::rand_f(0.8f, 1.4f);
    Vector a =
        Plot::azimuthal_unproject(Ra * std::cos(th), Ra * std::sin(th), basis);
    Vector b =
        Plot::azimuthal_unproject(Rb * std::cos(th), Rb * std::sin(th), basis);
    HS_EXPECT_NEAR(Plot::planar_arc_length(a, b, basis), angle_between(a, b),
                   1.2e-2f);
    ++radial;

    float a1 = hs::rand_f(-PI_F, PI_F), a2 = a1 + 2.4f;
    auto bow = [&](float rad) {
      Vector p = Plot::azimuthal_unproject(rad * std::cos(a1),
                                           rad * std::sin(a1), basis);
      Vector q = Plot::azimuthal_unproject(rad * std::cos(a2),
                                           rad * std::sin(a2), basis);
      return Plot::planar_arc_length(p, q, basis) - angle_between(p, q);
    };
    float lo = bow(1.0f), hi = bow(1.3f);
    HS_EXPECT_GT(lo, 1.5e-2f);
    HS_EXPECT_GT(hi, lo);
    ++azi;
  }
  HS_EXPECT_GT(radial, 1900);
  HS_EXPECT_GT(azi, 1900);
}

/**
 * @brief planar_arc_cumul is monotone and totals planar_arc_length.
 * @details Locks the table shared by the rasterizer's pre-pass and per-segment
 *          accumulator: it starts at 0, rises strictly, and its last entry is
 *          exactly planar_arc_length, so both consumers sum identical lengths.
 */
inline void test_planar_arc_cumul_monotone_and_endpoints() {
  hs::random().seed(0xF00D);
  int checked = 0;
  for (int trial = 0; trial < 2000; ++trial) {
    Basis basis = basis_from_normal(az_rand_unit());
    float R1 = hs::rand_f(0.1f, 1.3f), R2 = hs::rand_f(0.1f, 1.3f);
    float t1 = hs::rand_f(-PI_F, PI_F), t2 = t1 + hs::rand_f(0.3f, 1.5f);
    Vector a =
        Plot::azimuthal_unproject(R1 * std::cos(t1), R1 * std::sin(t1), basis);
    Vector b =
        Plot::azimuthal_unproject(R2 * std::cos(t2), R2 * std::sin(t2), basis);
    if (dot(a, basis.v) < -Plot::COS_PLANAR_ANTIPODE ||
        dot(b, basis.v) < -Plot::COS_PLANAR_ANTIPODE)
      continue;

    auto p1 = Plot::azimuthal_project(a, basis);
    auto p2 = Plot::azimuthal_project(b, basis);
    std::array<float, Plot::PLANAR_LEN_SAMPLES + 1> cumul;
    Plot::planar_arc_cumul(p1, p2.first - p1.first, p2.second - p1.second, basis,
                           cumul);

    HS_EXPECT_NEAR(cumul[0], 0.0f, 1e-6f);
    for (int k = 1; k <= Plot::PLANAR_LEN_SAMPLES; ++k)
      HS_EXPECT_GT(cumul[k], cumul[k - 1]);
    HS_EXPECT_NEAR(cumul[Plot::PLANAR_LEN_SAMPLES],
                   Plot::planar_arc_length(a, b, basis), 1e-6f);
    ++checked;
  }
  HS_EXPECT_GT(checked, 1500);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every plot/scan sampling test in this module.
 * @return Number of failed assertions reported across the module's tests.
 */
inline int run_plot_scan_tests() {
  hs_test::ModuleFixture fixture("plot_scan");

  hs::random().seed(1337);

  test_line_sample_endpoints_and_unit_length();
  test_line_sample_interior_between_endpoints();
  test_line_sample_degenerate_segment();
  test_line_sample_antipodal_stable_axis();

  test_clip_could_intersect_y();
  test_clip_x_band_topologies();
  test_edge_row_span_covers_arc_bulge();
  test_screen_step_matches_analytic_unclamped();

  test_ring_sample_unit_length_and_progress();
  test_ring_sample_lut_matches_direct();

  test_distorted_ring_sample_angle_addition_identity();
  test_distorted_ring_arc_matches_closed();
  test_distorted_ring_draw_culled_matches_closed();
  test_distorted_ring_draw_culled_runs_phase_chunks();

  test_spiral_sample_unit_length_and_monotone_arc();
  test_multiline_sample_arclength_param();

  test_star_sample_unit_length_closed();
  test_flower_sample_unit_length_closed();

  // rasterize() control-flow coverage. The 2*W steps_cache backstop is a
  // defensive path unreachable through any realistic single segment, so it is
  // not asserted here.
  test_rasterize_subpixel_open_segment_plots_both_endpoints();
  test_rasterize_open_segment_gap_free();
  test_rasterize_closed_loop_gap_free_no_dup();
  test_rasterize_antipodal_seam_planar_falls_back_geodesic();
  test_rasterize_planar_segment_gap_free_arclength();
  test_rasterize_planar_arc_registers_track_drawn_arc();
  test_rasterize_cull_follows_filter_orientation();

  test_particle_system_draws_active_trails_with_registers();
  test_particle_system_empty_zero_lifetime_is_noop();

  test_azimuthal_project_radius_is_geodesic_angle();
  test_azimuthal_roundtrip_identity();
  test_azimuthal_unproject_hits_great_circle_point();
  test_planar_arc_length_matches_fine_quadrature();
  test_dual_metric_radial_vs_azimuthal();
  test_planar_arc_cumul_monotone_and_endpoints();

  return fixture.result();
}

} // namespace plot_scan_tests
} // namespace hs_test

