/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file test_radial_surface.h
 * @brief First-arrival queries against independent swept-edge intersections. */

#include <algorithm>
#include <array>
#include <bit>
#include <cfloat>
#include <cmath>
#include <vector>
#include "core/render/scan/radial_surface.h"
#include "core/render/scan/radial_surface_edge_aa.h"
#include "core/render/scan/radial_surface_math.h"
#include "core/render/scan/radial_surface_renderer.h"
#include "core/render/scan/radial_surface_boundary.h"
#include "tests/test_harness.h"
#include "tests/test_fixture.h"

namespace hs_test::radial_surface_tests {

using Surface = Scan::RadialSurface<64, 9>;

/** @brief First-arrival grayscale ramp for a configured interval count. */
inline auto arrival_ramp(int intervals, bool inverted = false) {
  return [intervals, inverted](const auto &hit) {
    const float q = (hit.interval + hit.fraction) / intervals;
    const float value = inverted ? 1 - q : q;
    return Scan::SurfaceColor{value, value, value};
  };
}

/** @brief Double-precision chart point for the independent intersection oracle. */
struct ReferencePoint {
  double x, y;
};

/** @param a,b Chart vectors.
 * @return Signed determinant. */
inline double reference_cross(ReferencePoint a, ReferencePoint b) {
  return a.x * b.y - a.y * b.x;
}

/** @param c Similarity multiplier.
 * @param v Unit-star vertex.
 * @return Vertex after the similarity transform. */
inline ReferencePoint reference_multiply(Surface::Knot c, ReferencePoint v) {
  return {c.real * v.x - c.imag * v.y, c.imag * v.x + c.real * v.y};
}

/** @brief Computes a midpoint area-quadrature color reference.
 * @tparam Shade First-arrival shader.
 * @param surface Prepared radial surface.
 * @param p Pixel center in chart coordinates.
 * @param dx,dy Full pixel footprint axes.
 * @param n Samples along each footprint axis.
 * @param shade Maps a first arrival to its color.
 * @return Mean red channel over the footprint. */
template <typename Shade>
inline double area_reference(const Surface &surface, Surface::Point p,
                             Surface::Point dx, Surface::Point dy, int n,
                             const Shade &shade) {
  double total = 0;
  for (int y = 0; y < n; ++y)
    for (int x = 0; x < n; ++x) {
      const float U = (x + .5f) / n - .5f, V = (y + .5f) / n - .5f;
      const Surface::Point Q{p.x + U * dx.x + V * dy.x,
                             p.y + U * dx.y + V * dy.y};
      const auto hit = surface.template sample<false, true, false>(
          Q, Scan::radial_atan2(Q.y, Q.x));
      total += hit.found() ? shade(hit).r : 0;
    }
  return total / (n * n);
}

/** @brief Finds the first point on any swept edge without angular sectors.
 * @param p Chart point.
 * @param knots Similarity knots from the first swept frame.
 * @param sides Outer tip count.
 * @param ratio Inner-to-outer radius ratio.
 * @return Earliest knot coordinate, or -1 for an uncovered point. */
inline double reference_hit(ReferencePoint p,
                            std::span<const Surface::Knot> knots, int sides,
                            double ratio) {
  constexpr double TAU = 6.28318530717958647692;
  for (size_t i = 0; i + 1 < knots.size(); ++i) {
    double earliest = 2.0;
    for (int edge = 0; edge < 2 * sides; ++edge) {
      const double A = TAU * edge / (2 * sides);
      const double B = TAU * (edge + 1) / (2 * sides);
      const double RA = (edge & 1) ? ratio : 1.0;
      const double RB = (edge & 1) ? 1.0 : ratio;
      const ReferencePoint UNIT_A{RA * std::cos(A), RA * std::sin(A)};
      const ReferencePoint UNIT_B{RB * std::cos(B), RB * std::sin(B)};
      const auto a = reference_multiply(knots[i], UNIT_A);
      const auto b = reference_multiply(knots[i], UNIT_B);
      const auto end_a = reference_multiply(knots[i + 1], UNIT_A);
      const auto end_b = reference_multiply(knots[i + 1], UNIT_B);
      const ReferencePoint da{end_a.x - a.x, end_a.y - a.y};
      const ReferencePoint e{b.x - a.x, b.y - a.y};
      const ReferencePoint de{end_b.x - end_a.x - e.x, end_b.y - end_a.y - e.y};
      const ReferencePoint pa{p.x - a.x, p.y - a.y};
      const double C = reference_cross(pa, e);
      const double L = reference_cross(pa, de) - reference_cross(da, e);
      const double Q = -reference_cross(da, de);
      std::array<double, 2> roots{-1, -1};
      if (std::abs(Q) < 1e-25) {
        if (L != 0)
          roots[0] = -C / L;
      } else {
        const double DISCRIMINANT = L * L - 4 * Q * C;
        if (DISCRIMINANT >= 0) {
          roots[0] = (-L - std::sqrt(DISCRIMINANT)) / (2 * Q);
          roots[1] = (-L + std::sqrt(DISCRIMINANT)) / (2 * Q);
        }
      }
      for (double t : roots) {
        if (t < -1e-8 || t > 1 + 1e-8)
          continue;
        const ReferencePoint EDGE{e.x + t * de.x, e.y + t * de.y};
        const double LENGTH2 = EDGE.x * EDGE.x + EDGE.y * EDGE.y;
        if (LENGTH2 < 1e-24)
          continue;
        const double U =
            ((pa.x - t * da.x) * EDGE.x + (pa.y - t * da.y) * EDGE.y) / LENGTH2;
        if (U >= -1e-8 && U <= 1 + 1e-8)
          earliest = std::min(earliest, std::clamp(t, 0.0, 1.0));
      }
    }
    if (earliest <= 1)
      return double(i) + earliest;
  }
  return -1;
}

/** @brief Largest disagreement between two queries about one root's fraction.
 * @param hit Query result carrying the edge normal and predicate slope.
 * @param p Chart point.
 * @return Zero in IEEE builds, where every query evaluates the edge predicate
 * with the same operations. Under fast-math, four ulps of the predicate's
 * largest term |p|·|normal| over its interval slope: reassociation perturbs
 * the predicate by ulps of that term per call site, and the root moves by the
 * perturbation over the slope. */
inline float root_rounding_tolerance(const Surface::Hit &hit,
                                     Surface::Point p) {
#if defined(HS_TEST_FAST_MATH)
  const float TERM =
      sqrtf(p.x * p.x + p.y * p.y) *
      sqrtf(hit.normal.x * hit.normal.x + hit.normal.y * hit.normal.y);
  return TERM > 0.0f ? 4.0f * FLT_EPSILON * TERM / fabsf(hit.derivative) : 0.0f;
#else
  (void)hit;
  (void)p;
  return 0.0f;
#endif
}

/** @brief A tip sweep covers points outside both endpoint stars. */
inline void test_radial_surface_tip_sweep() {
  ScratchScope scope(scratch_arena_a);
  Surface surface;
  const std::array<Surface::Knot, 2> KNOTS{
      Surface::Knot{cosf(-.1f), sinf(-.1f), -.1f},
      Surface::Knot{cosf(.1f), sinf(.1f), .1f}};
  surface.prepare(scratch_arena_a, 7, .382f, KNOTS);
  HS_EXPECT_FALSE(surface.strictly_increasing_radius());
  Surface::Stats stats;
  const auto hit = surface.sample<true>({.99f, 0}, 0, 0, &stats);
  const auto expected = reference_hit({.99f, 0}, KNOTS, 7, .382f);
  HS_EXPECT_TRUE(hit.found());
  HS_EXPECT_NEAR(hit.fraction, expected, 2e-5);
  HS_EXPECT_EQ(stats.intervals, 1u);
  HS_EXPECT_EQ(stats.events, 1u);
  HS_EXPECT_EQ(stats.roots, 1u);
  HS_EXPECT_GT(hit.fraction, 0.0f);
  HS_EXPECT_LT(hit.fraction, .5f);
}

/** @brief A descending root survives cancellation at the start edge. */
inline void test_radial_surface_descending_root_cancellation() {
  ScratchScope scope(scratch_arena_a);
  Surface surface;
  const float PHASE = -0.3f;
  const std::array KNOTS{
      Surface::Knot{1, 0, 0},
      Surface::Knot{30 * cosf(PHASE), 30 * sinf(PHASE), PHASE}};
  surface.prepare(scratch_arena_a, 4, 0.01f, KNOTS);
  const auto VERTEX = surface.vertex(0);
  const auto EDGE = surface.edge(0);
  Surface::Point p{VERTEX.x + 0.5f * EDGE.x, VERTEX.y + 0.5f * EDGE.y};
  const ReferencePoint DELTA{double(KNOTS[1].real) - 1, KNOTS[1].imag};
  const ReferencePoint EDGE_DELTA{DELTA.x * EDGE.x - DELTA.y * EDGE.y,
                                  DELTA.y * EDGE.x + DELTA.x * EDGE.y};
  const double DETERMINANT =
      reference_cross({VERTEX.x, VERTEX.y}, {EDGE.x, EDGE.y});
  const double QUADRATIC =
      DETERMINANT * (DELTA.x * DELTA.x + DELTA.y * DELTA.y);
  for (int i = 0; i < 64; ++i)
    p.x = nextafterf(p.x, 0.0f);
  int descending = 0;
  for (int i = 0; i < 256; ++i) {
    const auto hit = surface.sample<false, false, false>(p, atan2f(p.y, p.x));
    HS_EXPECT_TRUE(hit.found());
    if (hit.fraction > 0.001f) {
      const double VALUE =
          reference_cross({p.x, p.y}, {EDGE.x, EDGE.y}) - DETERMINANT;
      const double LINEAR =
          reference_cross({p.x, p.y}, EDGE_DELTA) - 2 * DETERMINANT * DELTA.x;
      const double EXPECTED =
          (LINEAR + std::sqrt(LINEAR * LINEAR + 4 * QUADRATIC * VALUE)) /
          (2 * QUADRATIC);
      HS_EXPECT_NEAR(hit.fraction, EXPECTED, 1e-6);
      ++descending;
    }
    p.x = nextafterf(p.x, 2.0f);
  }
  HS_EXPECT_GT(descending, 0);
}

/** @brief Rotating contour gradients match central differences of arrival fraction. */
inline void test_radial_surface_fraction_gradient() {
  ScratchScope scope(scratch_arena_a);
  Surface surface;
  std::array<Surface::Knot, 9> knots;
  for (int i = 0; i <= 8; ++i) {
    const float radius = i / 8.0f;
    const float phase = sinf(math::PI_F * radius + .4f);
    knots[i] = {radius * cosf(phase), radius * sinf(phase), phase};
  }
  surface.prepare(scratch_arena_a, 7, .382f, knots);
  constexpr float STEP = .002f;
  std::vector<double> errors;
  auto sample = [&](Surface::Point point) {
    return surface.sample(point, atan2f(point.y, point.x));
  };
  for (int y = -40; y <= 40; ++y)
    for (int x = -40; x <= 40; ++x) {
      const Surface::Point point{x * .025f, y * .025f};
      const auto hit = sample(point);
      if (!hit.found() || hit.fraction < .1f || hit.fraction > .9f ||
          fabsf(hit.derivative) < 1e-5f)
        continue;
      const std::array adjacent{
          sample({point.x - STEP, point.y}), sample({point.x + STEP, point.y}),
          sample({point.x, point.y - STEP}), sample({point.x, point.y + STEP})};
      if (!std::all_of(adjacent.begin(), adjacent.end(),
                       [&](const auto &neighbor) {
                         return neighbor.interval == hit.interval &&
                                neighbor.edge == hit.edge;
                       }))
        continue;
      const double gx =
          (adjacent[1].fraction - adjacent[0].fraction) / (2 * STEP);
      const double gy =
          (adjacent[3].fraction - adjacent[2].fraction) / (2 * STEP);
      const double magnitude = std::max(1e-6, std::hypot(gx, gy));
      errors.push_back(std::hypot(hit.gradient.x - gx, hit.gradient.y - gy) /
                       magnitude);
      HS_EXPECT_NEAR(hit.gradient.x, -hit.normal.x / hit.derivative,
                     1e-5f * std::max(1.0f, fabsf(hit.gradient.x)));
      HS_EXPECT_NEAR(hit.gradient.y, -hit.normal.y / hit.derivative,
                     1e-5f * std::max(1.0f, fabsf(hit.gradient.y)));
    }
  HS_EXPECT_GT(errors.size(), size_t{100});
  if (errors.empty())
    return;
  std::sort(errors.begin(), errors.end());
  std::printf("  fraction gradient median %.8g p99 %.8g samples %zu\n",
              errors[errors.size() / 2], errors[errors.size() * 99 / 100],
              errors.size());
  HS_EXPECT_LT(errors[errors.size() / 2], .001);
}

/** @brief Degenerate intervals retain membership and finite first arrivals. */
inline void test_radial_surface_degenerate_intervals() {
  ScratchScope scope(scratch_arena_a);
  Surface surface;
  std::array<Surface::Knot, 2> knots{Surface::Knot{0, 0, 0},
                                     Surface::Knot{1, 0, 0}};
  surface.prepare(scratch_arena_a, 7, .382f, knots);
  const auto half = surface.sample({.5f, 0}, 0);
  HS_EXPECT_TRUE(half.found());
  HS_EXPECT_NEAR(half.fraction, .5f, 2e-6);
  HS_EXPECT_NEAR(half.gradient.x, 1.0f, 2e-5);
  HS_EXPECT_TRUE(surface.sample({0, 0}, 0).found());
  HS_EXPECT_FALSE(surface.sample({1.01f, 0}, 0).found());
  knots[0] = knots[1];
  surface.set_knots(knots);
  HS_EXPECT_TRUE(surface.sample({.1f, .1f}, math::PI_F / 4).found());
  HS_EXPECT_FALSE(surface.sample({2, 0}, 0).found());
  HS_EXPECT_FALSE(surface.sample({.1f, 0}, 0, 1).found());
  Surface thin;
  const std::array TINY_KNOTS{Surface::Knot{0, 0, 0},
                              Surface::Knot{1e-10f, 0, 0}};
  thin.prepare(scratch_arena_a, 7, 1e-30f, TINY_KNOTS);
  const std::array TINY_HITS{thin.sample<false, true, false>({5e-11f, 0}, 0),
                             thin.sample<false, false, false>({5e-11f, 0}, 0)};
  for (const auto &hit : TINY_HITS) {
    HS_EXPECT_TRUE(hit.found());
    HS_EXPECT_TRUE(std::isfinite(hit.fraction));
    HS_EXPECT_GE(hit.fraction, 0.0f);
    HS_EXPECT_LE(hit.fraction, 1.0f);
  }
  Surface near_corner;
  constexpr float RADIUS = 1e-13f, PHASE = 1e-20f;
  const std::array CORNER_KNOTS{
      Surface::Knot{RADIUS, 0, 0},
      Surface::Knot{2 * RADIUS, 2 * RADIUS * PHASE, PHASE}};
  near_corner.prepare(scratch_arena_a, 7, .382f, CORNER_KNOTS);
  const auto corner_hit =
      near_corner.sample<false, true, false>({1.5f * RADIUS, 0}, 0);
  HS_EXPECT_TRUE(corner_hit.found());
  HS_EXPECT_TRUE(std::isfinite(corner_hit.fraction));
  HS_EXPECT_GE(corner_hit.fraction, 0.0f);
  HS_EXPECT_LE(corner_hit.fraction, 1.0f);
}

/** @brief Coarse rejection retains fine radial bulges and intersecting disks. */
inline void test_radial_surface_block_bounds() {
  ScratchScope scope(scratch_arena_a);
  std::array<Surface::Knot, 9> knots;
  for (int i = 0; i <= 8; ++i)
    knots[size_t(i)] = {1.0f + .01f * (1.0f - std::abs(i - 4) / 4.0f), 0, 0};
  Surface surface;
  surface.prepare(scratch_arena_a, 7, .382f, knots);
  const auto a = surface.vertex(0), b = surface.vertex(1);
  const Surface::Point EDGE{b.x - a.x, b.y - a.y};
  const float LENGTH = sqrtf(EDGE.x * EDGE.x + EDGE.y * EDGE.y);
  const Surface::Point CENTER{.5f * (a.x + b.x), .5f * (a.y + b.y)};
  const Surface::Point BULGE{CENTER.x * 1.005f, CENTER.y * 1.005f};
  const float THETA = atan2f(BULGE.y, BULGE.x);
  HS_EXPECT_EQ(surface.skip_outside(BULGE, 0, 0), 0);
  const auto hit = surface.sample(BULGE, THETA);
  const double EXPECTED = reference_hit({BULGE.x, BULGE.y}, knots, 7, .382f);
  HS_EXPECT_TRUE(hit.found());
  HS_EXPECT_NEAR(float(hit.interval) + hit.fraction, EXPECTED, 1e-4);
  for (auto &k : knots)
    k = {1, 0, 0};
  surface.set_knots(knots);
  const Surface::Point NEAR{CENTER.x + .01f * EDGE.y / LENGTH,
                            CENTER.y - .01f * EDGE.x / LENGTH};
  HS_EXPECT_EQ(surface.skip_outside(NEAR, .02f, 0), 0);
  HS_EXPECT_EQ(surface.skip_outside(NEAR, .001f, 0), 2);
  HS_EXPECT_EQ(surface.skip_outside(NEAR, .001f, 3), 4);
  HS_EXPECT_EQ(surface.skip_outside(NEAR, .001f, 8), 8);
  HS_EXPECT_EQ(surface.skip_outside({0, 0}, .02f, 0), 0);
  HS_EXPECT_EQ(surface.skip_outside({2, 0}, .02f, 0), 2);
}

/** @brief Sector-rejected scans match all swept-edge roots in either turn. */
inline void test_radial_surface_matches_all_edges() {
  ScratchScope scope(scratch_arena_a);
  constexpr int N = 64;
  auto *knots = scratch_arena_a.allocate_n<Surface::Knot>(N + 1);
  Surface surface;
  surface.bind(scratch_arena_a, 7, .382f, N);
  float worst = 0;
  uint64_t coarse_checks = 0, skipped = 0, fine_visits = 0;
  uint64_t bounded_checks = 0, bounded_skipped = 0, bounded_fine = 0;
  for (float direction : {-1.0f, 1.0f}) {
    for (int i = 0; i <= N; ++i) {
      const float RADIUS = float(i) / N;
      const float PHASE = direction * sinf(3.0f * RADIUS + .4f);
      knots[i] = {RADIUS * cosf(PHASE), RADIUS * sinf(PHASE), PHASE};
    }
    surface.set_knots({knots, N + 1});
    for (int angle = 0; angle < 97; ++angle) {
      const float THETA = (angle + .371f) * math::TWO_PI_F / 97;
      for (int radius = 0; radius < 25; ++radius) {
        const float RHO = (radius + .217f) / 25;
        const Surface::Point P{RHO * cosf(THETA), RHO * sinf(THETA)};
        const double EXPECTED =
            reference_hit({P.x, P.y}, {knots, N + 1}, 7, .382f);
        Surface::Stats stats;
        const auto hit = surface.sample<true>(P, THETA, 0, &stats);
        const float ROOT_TOLERANCE = root_rounding_tolerance(hit, P);
        const auto unskipped = surface.sample<false, false>(P, THETA);
        const auto parameter_only =
            surface.sample<false, true, false>(P, THETA);
        HS_EXPECT_EQ(parameter_only.interval, hit.interval);
        HS_EXPECT_NEAR(parameter_only.fraction, hit.fraction, ROOT_TOLERANCE);
        HS_EXPECT_EQ(parameter_only.edge, hit.edge);
        Surface::Stats bounded_stats;
        const auto bounded =
            surface.sample<true>(P, THETA, int(RHO * N), &bounded_stats);
        coarse_checks += stats.coarse_checks;
        skipped += stats.skipped;
        fine_visits += stats.intervals;
        bounded_checks += bounded_stats.coarse_checks;
        bounded_skipped += bounded_stats.skipped;
        bounded_fine += bounded_stats.intervals;
        HS_EXPECT_EQ(hit.found(), EXPECTED >= 0);
        HS_EXPECT_EQ(hit.found(), bounded.found());
        HS_EXPECT_EQ(hit.interval, unskipped.interval);
        HS_EXPECT_NEAR(hit.fraction, unskipped.fraction, ROOT_TOLERANCE);
        HS_EXPECT_LE(stats.intervals, uint32_t(N));
        HS_EXPECT_LE(stats.roots, 1u);
        if (hit.found() && EXPECTED >= 0) {
          const float ERROR =
              fabsf(float(hit.interval) + hit.fraction - float(EXPECTED));
          worst = std::max(worst, ERROR);
          HS_EXPECT_LE(ERROR, 2e-3f);
          HS_EXPECT_EQ(hit.interval, bounded.interval);
          HS_EXPECT_NEAR(hit.fraction, bounded.fraction,
                         std::max(1e-6f, ROOT_TOLERANCE));
        }
      }
    }
  }
  std::printf("  radial surface worst knot-coordinate error %.9f\n", worst);
  std::printf("  radial blocks: checks=%llu skipped=%llu fine=%llu\n",
              static_cast<unsigned long long>(coarse_checks),
              static_cast<unsigned long long>(skipped),
              static_cast<unsigned long long>(fine_visits));
  std::printf("  radius-bounded blocks: checks=%llu skipped=%llu fine=%llu\n",
              static_cast<unsigned long long>(bounded_checks),
              static_cast<unsigned long long>(bounded_skipped),
              static_cast<unsigned long long>(bounded_fine));
  HS_EXPECT_GE(skipped, uint64_t{93000});
  HS_EXPECT_LT(bounded_fine, fine_visits);
}

/** @brief Chart-angle error across the float range and octant boundaries. */
inline void test_radial_surface_atan2() {
  double maximum = 0;
  float worst_x = 0, worst_y = 0;
  uint32_t failures = 0, cases = 0;
  auto check = [&](float y, float x) {
    const float ACTUAL = Scan::radial_atan2(y, x);
    const double EXPECTED = std::atan2(double(y), double(x));
    const double ERROR = std::abs(double(ACTUAL) - EXPECTED);
    if (!std::isfinite(ACTUAL) || ERROR > 5e-7)
      ++failures;
    if (ERROR > maximum) {
      maximum = ERROR;
      worst_x = x;
      worst_y = y;
    }
    ++cases;
  };
  for (int x = -150; x <= 150; ++x)
    for (int y = -150; y <= 150; ++y)
      check(y / 37.0f, x / 53.0f);
  uint32_t random = 0x7319acbd;
  auto next_float = [&] {
    random = random * 1664525u + 1013904223u;
    uint32_t bits = random;
    if ((bits & 0x7f800000u) == 0x7f800000u)
      bits ^= 0x00800000u;
    return std::bit_cast<float>(bits);
  };
  for (int i = 0; i < 500000; ++i) {
    const float X = next_float(), Y = next_float();
    check(Y, X);
  }
  for (float ratio : {0.0f, .4142135623730950488f, 1.0f})
    for (int exponent = -120; exponent <= 120; ++exponent)
      for (int dx : {-1, 1})
        for (int dy : {-1, 1}) {
          const float SCALE = std::ldexp(1.0f, exponent);
          const float BASE = ratio * SCALE;
          for (float small : {BASE, std::nextafter(BASE, 0.0f),
                              std::nextafter(BASE, INFINITY)}) {
            check(dy * small, dx * SCALE);
            check(dy * SCALE, dx * small);
          }
        }
  for (float x : {-0.0f, 0.0f, -1.0f, 1.0f})
    for (float y : {-0.0f, 0.0f}) {
      const float ACTUAL = Scan::radial_atan2(y, x);
      const float EXPECTED = std::atan2(y, x);
      HS_EXPECT_EQ(std::signbit(ACTUAL), std::signbit(EXPECTED));
      HS_EXPECT_NEAR(ACTUAL, EXPECTED, 5e-7);
    }
  HS_EXPECT_EQ(failures, 0u);
  std::printf("  radial atan2: cases=%u max_error=%.10g x=%.8g y=%.8g\n", cases,
              maximum, worst_x, worst_y);
}

/** @brief Footprint witnesses imply independent corner winding containment. */
inline void test_radial_surface_footprint_witness() {
  ScratchScope scope(scratch_arena_a);
  Surface surface;
  const float PHASE = .231f;
  const std::array KNOTS{Surface::Knot{0, 0, PHASE},
                         Surface::Knot{cosf(PHASE), sinf(PHASE), PHASE}};
  surface.prepare(scratch_arena_a, 7, .382f, KNOTS);
  auto inside = [&](ReferencePoint p) {
    bool result = false;
    for (int i = 0; i < 14; ++i) {
      const auto va = surface.vertex(i), vb = surface.vertex((i + 1) % 14);
      const auto A = reference_multiply(KNOTS[1], {va.x, va.y});
      const auto B = reference_multiply(KNOTS[1], {vb.x, vb.y});
      if ((A.y > p.y) != (B.y > p.y) &&
          p.x < (B.x - A.x) * (p.y - A.y) / (B.y - A.y) + A.x)
        result = !result;
    }
    return result;
  };
  uint32_t certificates = 0, failures = 0;
  auto check = [&](Surface::Point p, Surface::Point dx, Surface::Point dy) {
    const float THETA = Scan::radial_atan2(p.y, p.x);
    if (!surface.contains_footprint_at_knot(p, dx, dy, THETA, 1))
      return;
    ++certificates;
    for (int x : {-1, 1})
      for (int y : {-1, 1})
        if (!inside({p.x + .5 * (x * dx.x + y * dy.x),
                     p.y + .5 * (x * dx.y + y * dy.y)}))
          ++failures;
  };
  for (int x = -100; x <= 100; ++x)
    for (int y = -100; y <= 100; ++y) {
      const float ANGLE = .31f * x + .17f * y;
      check({x / 100.0f, y / 100.0f},
            {.017f * cosf(ANGLE), .017f * sinf(ANGLE)},
            {-.011f * sinf(ANGLE), .011f * cosf(ANGLE)});
    }
  for (int edge = 0; edge < 14; ++edge) {
    const auto va = surface.vertex(edge), ve = surface.edge(edge);
    const auto A = reference_multiply(KNOTS[1], {va.x, va.y});
    const auto E = reference_multiply(KNOTS[1], {ve.x, ve.y});
    for (int i = 0; i <= 100; ++i)
      for (float epsilon : {-1e-5f, -1e-7f, 0.0f, 1e-7f, 1e-5f}) {
        const float T = i / 100.0f;
        check({float(A.x + T * E.x + epsilon * E.y),
               float(A.y + T * E.y - epsilon * E.x)},
              {1e-7f, 0}, {0, 1e-7f});
      }
  }
  hs::Pcg32 random(0x42da5819);
  auto coordinate = [&] { return rand_uniform(random, -1.0f, 1.0f); };
  for (int i = 0; i < 20000; ++i) {
    const Surface::Point P{coordinate(), coordinate()};
    const Surface::Point DX{.05f * coordinate(), .05f * coordinate()};
    const Surface::Point DY{.05f * coordinate(), .05f * coordinate()};
    check(P, DX, DY);
  }
  HS_EXPECT_GT(certificates, 1000u);
  HS_EXPECT_EQ(failures, 0u);
  HS_EXPECT_FALSE(
      surface.contains_footprint_at_knot({.2f, 0}, {.01f, 0}, {0, .01f}, 0, 0));
  std::printf("  radial footprint: certificates=%u false=%u\n", certificates,
              failures);
}

/** @brief Two-color edge AA matches the area reference within its error bound. */
inline void test_radial_surface_edge_aa() {
  ScratchScope scope(scratch_arena_a);
  using AA = Scan::RadialSurfaceEdgeAA<Surface>;
  HS_EXPECT_NEAR(AA::coverage(1, 0, 0), .5f, 1e-7f);
  HS_EXPECT_NEAR(AA::coverage(1, 1, .5f), .875f, 1e-7f);
  HS_EXPECT_NEAR(AA::coverage(1, 1, -.5f), .125f, 1e-7f);
  HS_EXPECT_NEAR(AA::coverage(1, 0, .25f), .75f, 1e-7f);
  constexpr int N = 64;
  Surface surface;
  auto *knots = scratch_arena_a.allocate_n<Surface::Knot>(N + 1);
  for (int i = 0; i <= N; ++i) {
    const float R = float(i) / N, P = sinf(3 * R + .4f);
    knots[i] = {R * cosf(P), R * sinf(P), P};
  }
  surface.prepare(scratch_arena_a, 7, .382f, {knots, N + 1});
  AA aa;
  Scan::SurfaceAreaStats stats;
  double aa_error = 0, point_error = 0;
  float worst_error = 0, worst_x = 0, worst_y = 0;
  uint32_t worse = 0;
  const auto shade = arrival_ramp(N);
  for (int x = -20; x <= 20; ++x)
    for (int y = -20; y <= 20; ++y) {
      const Surface::Point P{x * .05f, y * .05f};
      const auto color =
          aa.sample(surface, P, {.05f, 0}, {0, .05f},
                    Scan::radial_atan2(P.y, P.x), 0, shade, &stats);
      HS_EXPECT_TRUE(std::isfinite(color.r));
      HS_EXPECT_GE(color.r, 0);
      HS_EXPECT_LE(color.r, 1);
      const double reference =
          area_reference(surface, P, {.05f, 0}, {0, .05f}, 12, shade);
      const auto center = surface.template sample<false, true, false>(
          P, Scan::radial_atan2(P.y, P.x));
      const double POINT = center.found() ? shade(center).r : 0;
      const float ERROR = float(std::abs(color.r - reference));

      aa_error += ERROR;
      point_error += std::abs(POINT - reference);
      worse += ERROR > std::abs(POINT - reference) + .02;
      if (ERROR > worst_error) {
        worst_error = ERROR;
        worst_x = P.x;
        worst_y = P.y;
      }
    }
  HS_EXPECT_GT(stats.edge_pixels, 0u);
  HS_EXPECT_GT(stats.clips, 0u);
  HS_EXPECT_LT(worst_error, .1f);
  HS_EXPECT_LT(aa_error, point_error);
  HS_EXPECT_LT(aa_error / 1681, 0.005);
  HS_EXPECT_LE(worse, 10u);
  std::printf(
      "  radial edge quality: aa_mae=%.8f point_mae=%.8f worst=%.6f at=(%.4f,%.4f) worse=%u\n",
      aa_error / 1681, point_error / 1681, worst_error, worst_x, worst_y,
      worse);
  std::printf(
      "  radial edge AA: pixels=%llu edge=%llu point=%llu outside=%llu total=%llu\n",
      static_cast<unsigned long long>(stats.pixels),
      static_cast<unsigned long long>(stats.edge_pixels),
      static_cast<unsigned long long>(stats.point_intervals),
      static_cast<unsigned long long>(stats.outside_intervals),
      static_cast<unsigned long long>(stats.intervals));
  std::printf("  radial paired pixels=%llu\n",
              static_cast<unsigned long long>(stats.clips));
}

/** @brief First-arrival knot envelopes retain their visible extent over a full phase. */
inline void test_radial_surface_fold_visibility() {
  ScratchScope scope(scratch_arena_a);
  constexpr int N = 64;
  Surface surface;
  surface.bind(scratch_arena_a, 7, .382f, N);
  auto *knots = scratch_arena_a.allocate_n<Surface::Knot>(N + 1);
  constexpr float DENSITY_FLOOR = .5f;
  constexpr float BREAK_ANGLE = math::PI_F / 6;
  constexpr float DENSITY_INTEGRAL = 1.1278247916f;
  constexpr float BREAK_QUANTILE =
      DENSITY_FLOOR * BREAK_ANGLE / DENSITY_INTEGRAL;
  uint32_t segments = 0, visible = 0;
  float longest = 0, max_gap = 0;
  for (int phase = 0; phase < 64; ++phase) {
    for (int i = 0; i <= N; ++i) {
      const float U = float(i) / N;
      const float R = U < BREAK_QUANTILE ? U * DENSITY_INTEGRAL / DENSITY_FLOOR
                                         : acosf(DENSITY_INTEGRAL * (1 - U));
      const float P = sinf(math::PI_F * U + phase * math::TWO_PI_F / 64);
      knots[i] = {R * cosf(P), R * sinf(P), P};
    }
    surface.set_knots({knots, N + 1});
    HS_EXPECT_TRUE(surface.strictly_increasing_radius());
    Scan::RadialSurfaceBoundaryCandidates::visit_candidates(
        surface, [&](const auto &candidate) {
          if (!candidate.knot)
            return;
          ++segments;
          bool any = false;
          for (int j = 1; j <= 5; ++j) {
            const float T = j / 6.0f;
            const Surface::Point P{
                candidate.a.x + T * (candidate.b.x - candidate.a.x),
                candidate.a.y + T * (candidate.b.y - candidate.a.y)};
            const auto hit = surface.template sample<false, true, false>(
                P, Scan::radial_atan2(P.y, P.x));
            const bool VISIBLE =
                !hit.found() ||
                hit.interval + hit.fraction >= candidate.origin_interval - .01f;
            any |= VISIBLE;
            if (VISIBLE) {
              const float R = hypotf(P.x, P.y);
              max_gap =
                  fmaxf(max_gap, hypotf(knots[candidate.origin_interval].real,
                                        knots[candidate.origin_interval].imag) -
                                     R);
            }
          }
          if (any) {
            ++visible;
            longest = fmaxf(longest, hypotf(candidate.a.x - candidate.b.x,
                                            candidate.a.y - candidate.b.y));
          }
        });
  }
  HS_EXPECT_GE(segments, 900u);
  HS_EXPECT_LE(segments, 960u);
  HS_EXPECT_GE(visible, 580u);
  HS_EXPECT_LE(visible, 630u);
  HS_EXPECT_GE(longest, 0.50f);
  HS_EXPECT_LE(longest, 0.55f);
  HS_EXPECT_GE(max_gap, 0.68f);
  HS_EXPECT_LE(max_gap, 0.70f);
  std::printf(
      "  radial fold visibility: segments=%u visible=%u longest=%g radial_gap=%g\n",
      segments, visible, longest, max_gap);
}

/** @brief Edge shading matches the area reference at a partially covered pixel. */
inline void test_radial_surface_fixture_edge_trace() {
  ScratchScope scope(scratch_arena_a);
  const std::array KNOTS{
      Surface::Knot{float(0), float(0), 0.72896862f},
      Surface::Knot{0.0270809401f, 0.0225565769f, 0.6945014f},
      Surface::Knot{0.0557565615f, 0.0431266911f, 0.658361077f},
      Surface::Knot{0.0860152766f, 0.0614895262f, 0.620634735f},
      Surface::Knot{0.117813617f, 0.0774259195f, 0.58141315f},
      Surface::Knot{0.151075974f, 0.0907218978f, 0.540790915f},
      Surface::Knot{0.185694739f, 0.101172216f, 0.498865962f},
      Surface::Knot{0.221531272f, 0.108584233f, 0.455739081f},
      Surface::Knot{0.258417428f, 0.1127818f, 0.411514282f},
      Surface::Knot{0.296157509f, 0.113609128f, 0.366298229f},
      Surface::Knot{0.334531397f, 0.110934265f, 0.320199609f},
      Surface::Knot{0.373297751f, 0.104652569f, 0.273329675f},
      Surface::Knot{0.412198186f, 0.0946896374f, 0.2258012f},
      Surface::Knot{0.45096156f, 0.0810035169f, 0.177728802f},
      Surface::Knot{0.489308983f, 0.0635868534f, 0.129228204f},
      Surface::Knot{0.526937485f, 0.0424659103f, 0.0804163814f},
      Surface::Knot{0.562327325f, 0.0176688246f, 0.0314107239f},
      Surface::Knot{0.594732523f, -0.0105103729f, -0.0176706016f},
      Surface::Knot{0.624187768f, -0.0417007692f, -0.0667092502f},
      Surface::Knot{0.650704622f, -0.0755498782f, -0.115587324f},
      Surface::Knot{0.67428875f, -0.111715004f, -0.164186954f},
      Surface::Knot{0.694952309f, -0.149861753f, -0.212390959f},
      Surface::Knot{0.712721169f, -0.189662695f, -0.260083288f},
      Surface::Knot{0.72764045f, -0.230797902f, -0.307149023f},
      Surface::Knot{0.739777207f, -0.272956789f, -0.353474855f},
      Surface::Knot{0.749223769f, -0.315839022f, -0.398949146f},
      Surface::Knot{0.756096363f, -0.359158784f, -0.443462372f},
      Surface::Knot{0.760536909f, -0.402644694f, -0.486907214f},
      Surface::Knot{0.762710333f, -0.446044207f, -0.529179096f},
      Surface::Knot{0.762804151f, -0.489122957f, -0.570176005f},
      Surface::Knot{0.761025906f, -0.531668007f, -0.609799385f},
      Surface::Knot{0.757600129f, -0.573488235f, -0.647953749f},
      Surface::Knot{0.75276643f, -0.61441505f, -0.684547186f},
      Surface::Knot{0.746776521f, -0.654302001f, -0.719491422f},
      Surface::Knot{0.739891231f, -0.69302398f, -0.752702296f},
      Surface::Knot{0.732377052f, -0.730477929f, -0.784099877f},
      Surface::Knot{0.724504232f, -0.766579449f, -0.813608468f},
      Surface::Knot{0.716543794f, -0.801262558f, -0.841157019f},
      Surface::Knot{0.708764195f, -0.834477246f, -0.866679192f},
      Surface::Knot{0.7014305f, -0.86618638f, -0.890113413f},
      Surface::Knot{0.69480139f, -0.896364093f, -0.911403358f},
      Surface::Knot{0.689127088f, -0.924992561f, -0.930497527f},
      Surface::Knot{0.684648991f, -0.952058792f, -0.947350144f},
      Surface::Knot{0.681596339f, -0.977553189f, -0.961920381f},
      Surface::Knot{0.680185676f, -1.00146544f, -0.974173427f},
      Surface::Knot{0.68062067f, -1.0237813f, -0.984079599f},
      Surface::Knot{0.68309021f, -1.04448164f, -0.991614938f},
      Surface::Knot{0.687765598f, -1.06353951f, -0.996761382f},
      Surface::Knot{0.694802642f, -1.08091712f, -0.999506593f},
      Surface::Knot{0.704337895f, -1.09656513f, -0.999843895f},
      Surface::Knot{0.716490448f, -1.11041999f, -0.997772455f},
      Surface::Knot{0.731357992f, -1.12240326f, -0.993297338f},
      Surface::Knot{0.749015927f, -1.13242221f, -0.986429214f},
      Surface::Knot{0.769518554f, -1.14036703f, -0.977184832f},
      Surface::Knot{0.792895377f, -1.14611232f, -0.965586185f},
      Surface::Knot{0.819150507f, -1.14951777f, -0.951661408f},
      Surface::Knot{0.848262131f, -1.15042794f, -0.935443997f},
      Surface::Knot{0.880178809f, -1.1486764f, -0.916973054f},
      Surface::Knot{0.914821804f, -1.14408469f, -0.896293044f},
      Surface::Knot{0.95208174f, -1.1364665f, -0.873453617f},
      Surface::Knot{0.991816998f, -1.12563121f, -0.848510206f},
      Surface::Knot{1.03385639f, -1.11138487f, -0.821522474f},
      Surface::Knot{1.07799459f, -1.09353769f, -0.792555749f},
      Surface::Knot{1.12399709f, -1.07190359f, -0.76167953f},
      Surface::Knot{1.171597f, -1.0463084f, -0.728968561f}};
  Surface surface;
  surface.prepare(scratch_arena_a, 7, .382f, KNOTS);
  const Surface::Point P{-.435921282f, -.661391735f},
      DX{.0171549805f, -.00341231585f}, DY{.00249151839f, .0213687979f};
  Scan::RadialSurfaceEdgeAA<Surface> aa;
  const auto shade = arrival_ramp(64);
  const auto result = aa.sample(surface, P, DX, DY, -2.15354013f, 23, shade);
  const auto zero = aa.sample(surface, P, DX, DY, -2.15354013f, 0, shade);
  const float RADIUS = sqrtf(P.x * P.x + P.y * P.y);
  const float EXTENT = .5f * (sqrtf(DX.x * DX.x + DX.y * DX.y) +
                              sqrtf(DY.x * DY.x + DY.y * DY.y));
  const auto prepared = aa.sample(surface, P, DX, DY, -2.15354013f, 23, shade,
                                  nullptr, RADIUS, EXTENT);
  HS_EXPECT_NEAR(prepared.r, result.r, 1e-6);
  Scan::SurfaceAreaStats disabled_stats;
  const auto uncounted =
      aa.template sample<false>(surface, P, DX, DY, -2.15354013f, 23, shade,
                                &disabled_stats, RADIUS, EXTENT);
  HS_EXPECT_NEAR(uncounted.r, result.r, 1e-6);
  HS_EXPECT_EQ(disabled_stats.pixels, 0u);
  const double reference = area_reference(surface, P, DX, DY, 32, shade);
  HS_EXPECT_GT(result.r, 0.0f);
  HS_EXPECT_NEAR(result.r, reference, 0.005);
  std::printf("  fixture edge trace: q=%.8f first0=%.8f reference=%.8f\n",
              result.r, zero.r, reference);
  HS_EXPECT_NEAR(result.r, zero.r, 1e-6);
}

/** @brief A pixel spanning two terminal tip edges uses both half-planes. */
inline void test_radial_surface_edge_pair_corner() {
  ScratchScope scope(scratch_arena_a);
  Surface surface;
  const std::array KNOTS{Surface::Knot{0, 0, 0}, Surface::Knot{1, 0, 0}};
  surface.prepare(scratch_arena_a, 4, .382f, KNOTS);
  Scan::RadialSurfaceEdgeAA<Surface> aa;
  Scan::SurfaceAreaStats stats;
  const auto color = aa.sample(
      surface, {.98f, 0}, {.1f, 0}, {0, .1f}, 0, 0,
      [](const Surface::Hit &) { return Scan::SurfaceColor{1, 1, 1}; }, &stats);
  const float INNER = .382f * sqrtf(.5f);
  const float HALF_PLANE = .5f + .02f * INNER / (.1f * (1 - INNER));
  HS_EXPECT_EQ(stats.clips, 1u);
  const auto shade = [](const auto &) { return Scan::SurfaceColor{1, 1, 1}; };
  const double reference =
      area_reference(surface, {.98f, 0}, {.1f, 0}, {0, .1f}, 128, shade);
  std::printf("  terminal corner coverage %.6f area %.6f\n", color.r,
              reference);
  HS_EXPECT_NEAR(color.r, reference, .16);
  HS_EXPECT_LT(color.r, HALF_PLANE - .1f);
}

/** @brief Checks a terminal-boundary overlap against subpixel area.
 * @param knots Surface contour.
 * @param p Pixel center.
 * @param dx,dy Pixel frame.
 * @param theta Pixel azimuth.
 * @param lower_bound First candidate interval.
 * @param tolerance Maximum color error.
 * @param label Diagnostic label. */
inline void check_terminal_overlap(std::span<const Surface::Knot> knots,
                                   Surface::Point p, Surface::Point dx,
                                   Surface::Point dy, float theta,
                                   uint16_t lower_bound, double tolerance,
                                   const char *label) {
  HS_CONTEXT(label);
  ScratchScope scope(scratch_arena_a);
  Surface surface;
  surface.prepare(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, knots);
  Scan::RadialSurfaceEdgeAA<Surface> aa;
  const auto shade = arrival_ramp(static_cast<int>(knots.size() - 1));
  Scan::SurfaceAreaStats stats;
  const auto result =
      aa.sample(surface, p, dx, dy, theta, lower_bound, shade, &stats);
  const double reference = area_reference(surface, p, dx, dy, 128, shade);
  HS_EXPECT_EQ(stats.clips, 1u);
  HS_EXPECT_NEAR(result.r, reference, tolerance);
  std::printf("  terminal overlap %s: aa=%g reference=%g\n", label, result.r,
              reference);
}

/** @brief Checks overlapping early and terminal boundaries against subpixel area. */
inline void test_radial_surface_terminal_background() {
  {
    const std::array KNOTS{
        Surface::Knot{float(0), float(0), 0.72896862f},
        Surface::Knot{0.075755164f, 0.0556274317f, 0.63338083f},
        Surface::Knot{0.162471071f, 0.0945313275f, 0.526955783f},
        Surface::Knot{0.258417428f, 0.1127818f, 0.411514282f},
        Surface::Knot{0.360347688f, 0.107152462f, 0.289031774f},
        Surface::Knot{0.463804126f, 0.075611867f, 0.161603808f},
        Surface::Knot{0.562327325f, 0.0176688246f, 0.0314107239f},
        Surface::Knot{0.642191768f, -0.0639929622f, -0.099319756f},
        Surface::Knot{0.70119518f, -0.162961051f, -0.228350967f},
        Surface::Knot{0.739777207f, -0.272956789f, -0.353474855f},
        Surface::Knot{0.75931716f, -0.388147056f, -0.47255075f},
        Surface::Knot{0.762408674f, -0.503373921f, -0.583541274f},
        Surface::Knot{0.75276643f, -0.61441505f, -0.684547186f},
        Surface::Knot{0.734934986f, -0.718139827f, -0.773840308f},
        Surface::Knot{0.713917315f, -0.812499404f, -0.849892735f},
        Surface::Knot{0.69480139f, -0.896364093f, -0.911403358f},
        Surface::Knot{0.682441533f, -0.969230711f, -0.957319498f},
        Surface::Knot{0.681209445f, -1.03086174f, -0.986855686f},
        Surface::Knot{0.694802642f, -1.08091712f, -0.999506593f},
        Surface::Knot{0.7260952f, -1.11862218f, -0.995055556f},
        Surface::Knot{0.776990771f, -1.14253318f, -0.97357893f},
        Surface::Knot{0.848262131f, -1.15042794f, -0.935443997f},
        Surface::Knot{0.939379215f, -1.13935375f, -0.88130337f},
        Surface::Knot{1.04834688f, -1.10584521f, -0.812083423f},
        Surface::Knot{1.171597f, -1.0463084f, -0.728968561f}};
    const Surface::Point P{0.126707062f, -0.879523873f};
    const Surface::Point DX{0.0154463407f, -0.00674017007f};
    const Surface::Point DY{0.0113697099f, 0.0190199912f};
    check_terminal_overlap(KNOTS, P, DX, DY, -1.42771757f, 10, 0.08,
                           "folded terminal overlap");
  }
  {
    const std::array KNOTS{
        Surface::Knot{float(0), float(0), 0.368096858f},
        Surface::Knot{0.0912108868f, 0.0226678196f, 0.243586123f},
        Surface::Knot{0.186731219f, 0.0215517394f, 0.114907548f},
        Surface::Knot{0.281921297f, -0.00443695905f, -0.0157371461f},
        Surface::Knot{0.371935785f, -0.0547344834f, -0.146112502f},
        Surface::Knot{0.452398509f, -0.127149373f, -0.273987889f},
        Surface::Knot{0.518810213f, -0.217624038f, -0.397175252f},
        Surface::Knot{0.562117755f, -0.317063063f, -0.513566792f},
        Surface::Knot{0.58540678f, -0.418962836f, -0.621171176f},
        Surface::Knot{0.593782067f, -0.518843412f, -0.718147039f},
        Surface::Knot{0.592395544f, -0.613423228f, -0.802835226f},
        Surface::Knot{0.586460471f, -0.700511515f, -0.873786807f},
        Surface::Knot{0.581069171f, -0.778795123f, -0.929787457f},
        Surface::Knot{0.580974758f, -0.847539246f, -0.969879329f},
        Surface::Knot{0.590400279f, -0.906234384f, -0.993376136f},
        Surface::Knot{0.612883449f, -0.954249144f, -0.999876201f},
        Surface::Knot{0.651133478f, -0.99053508f, -0.989267945f},
        Surface::Knot{0.706889272f, -1.01342487f, -0.961733162f},
        Surface::Knot{0.780762732f, -1.02055955f, -0.917742789f},
        Surface::Knot{0.872070372f, -1.00897121f, -0.858049512f},
        Surface::Knot{0.978691936f, -0.975324929f, -0.783675015f},
        Surface::Knot{1.09699905f, -0.916311085f, -0.695891261f},
        Surface::Knot{1.22190893f, -0.829155505f, -0.596200883f},
        Surface::Knot{1.34712458f, -0.712165892f, -0.48630923f},
        Surface::Knot{1.46557486f, -0.565235972f, -0.368096769f}};
    const Surface::Point P{-0.349607706f, -0.830267727f};
    const Surface::Point DX{0.0242360439f, 0.00149365724f};
    const Surface::Point DY{-0.00394339021f, 0.0215454288f};
    check_terminal_overlap(KNOTS, P, DX, DY, -1.96934056f, 10, 0.13,
                           "rotated terminal overlap");
  }
  {
    const std::array KNOTS{
        Surface::Knot{float(0), float(0), 0.99361074f},
        Surface::Knot{0.0507930592f, 0.0790779442f, 0.99984163f},
        Surface::Knot{0.103300296f, 0.157041624f, 0.988964915f},
        Surface::Knot{0.161437929f, 0.231164649f, 0.961166739f},
        Surface::Knot{0.228672355f, 0.298397481f, 0.916922748f},
        Surface::Knot{0.307668567f, 0.355206162f, 0.85698998f},
        Surface::Knot{0.399015129f, 0.396624714f, 0.782393813f},
        Surface::Knot{0.495923966f, 0.412994742f, 0.694410741f},
        Surface::Knot{0.596352816f, 0.403229833f, 0.59454608f},
        Surface::Knot{0.697771132f, 0.367275268f, 0.484508574f},
        Surface::Knot{0.796234846f, 0.305336982f, 0.366181046f},
        Surface::Knot{0.887061715f, 0.218572244f, 0.241587907f},
        Surface::Knot{0.965498745f, 0.109432288f, 0.112861246f},
        Surface::Knot{1.027385f, -0.0182857905f, -0.0177966207f},
        Surface::Knot{1.06974053f, -0.159651712f, -0.148149729f},
        Surface::Knot{1.09120262f, -0.309022009f, -0.275968194f},
        Surface::Knot{1.09224212f, -0.460588574f, -0.39906469f},
        Surface::Knot{1.07513607f, -0.608937204f, -0.515332997f},
        Surface::Knot{1.04372191f, -0.749517918f, -0.622783899f},
        Surface::Knot{1.00298965f, -0.878943443f, -0.719578981f},
        Surface::Knot{0.95860666f, -0.995072722f, -0.804061472f},
        Surface::Knot{0.916440845f, -1.09689069f, -0.874786615f},
        Surface::Knot{0.882160842f, -1.1842097f, -0.930543721f},
        Surface::Knot{0.860919476f, -1.25727582f, -0.970378935f},
        Surface::Knot{0.857132673f, -1.31633019f, -0.99361074f}};
    const Surface::Point P{0.851393998f, 0.705435812f};
    const Surface::Point DX{0.0215828046f, -0.0114967246f};
    const Surface::Point DY{0.0069451537f, 0.0210930854f};
    check_terminal_overlap(KNOTS, P, DX, DY, 0.691918015f, 13, 0.04,
                           "offset terminal overlap");
  }
}

/** @brief Whole-footprint occlusion takes the point path while wider pixels retain AA. */
inline void test_radial_surface_hidden_terminal() {
  ScratchScope scope(scratch_arena_a);
  const std::array KNOTS{
      Surface::Knot{0.0f, 0.0f, 0.702642024f},
      Surface::Knot{0.0773696452f, 0.0533590987f, 0.603755772f},
      Surface::Knot{0.16544956f, 0.08921583f, 0.494539082f},
      Surface::Knot{0.262169778f, 0.103760779f, 0.376860648f},
      Surface::Knot{0.36399886f, 0.0940049887f, 0.252734095f},
      Surface::Knot{0.466302335f, 0.0582537539f, 0.124283135f},
      Surface::Knot{0.562593699f, -0.00354114221f, -0.00629431009f},
      Surface::Knot{0.639346004f, -0.0879887938f, -0.13676405f},
      Surface::Knot{0.694773316f, -0.188470274f, -0.26489377f},
      Surface::Knot{0.729767799f, -0.298688233f, -0.388491005f},
      Surface::Knot{0.746142507f, -0.412906766f, -0.50544107f},
      Surface::Knot{0.746860683f, -0.526166737f, -0.613743007f},
      Surface::Knot{0.735907078f, -0.634510756f, -0.71154356f},
      Surface::Knot{0.717983186f, -0.735087991f, -0.797169447f},
      Surface::Knot{0.69813484f, -0.826099873f, -0.869155407f},
      Surface::Knot{0.681399107f, -0.906594098f, -0.926270008f},
      Surface::Knot{0.672504127f, -0.976152003f, -0.967535794f},
      Surface::Knot{0.675642073f, -1.0345192f, -0.992246747f},
      Surface::Knot{0.694290698f, -1.08124602f, -0.999980152f},
      Surface::Knot{0.731068313f, -1.11537838f, -0.990603626f},
      Surface::Knot{0.787583947f, -1.13525689f, -0.964277625f},
      Surface::Knot{0.864274681f, -1.1384474f, -0.921452403f},
      Surface::Knot{0.96023047f, -1.12183678f, -0.862861037f},
      Surface::Knot{1.07304466f, -1.08189654f, -0.789505839f},
      Surface::Knot{1.19873369f, -1.01510525f, -0.702641964f}};
  Surface surface;
  surface.prepare(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, KNOTS);
  Scan::RadialSurfaceEdgeAA<Surface> aa;
  for (int sector = 0; sector < 7; ++sector) {
    const float ANGLE = sector * math::TWO_PI_F / 7;
    auto rotate = [&](Surface::Point p) {
      return Surface::Point{p.x * cosf(ANGLE) - p.y * sinf(ANGLE),
                            p.x * sinf(ANGLE) + p.y * cosf(ANGLE)};
    };
    const auto P = rotate({.220684409f, -.662453413f});
    const auto DX = rotate({.00910583977f, .0197949838f});
    const auto DY = rotate({-.021678308f, .00790116936f});
    double reference = 0;
    for (int y = 0; y < 16; ++y)
      for (int x = 0; x < 16; ++x) {
        const double U = (x + .5) / 16 - .5, V = (y + .5) / 16 - .5;
        const double HIT = reference_hit(
            {P.x + U * DX.x + V * DY.x, P.y + U * DX.y + V * DY.y}, KNOTS, 7,
            Plot::STAR_INNER_RATIO);
        HS_EXPECT_GT(HIT, 8.0);
        HS_EXPECT_LT(HIT, 10.0);
        reference += 1 - HIT / 24;
      }
    reference /= 256;
    const auto shade = arrival_ramp(24, true);
    Scan::SurfaceAreaStats stats;
    int calls = 0;
    const float THETA = Scan::radial_atan2(P.y, P.x);
    const auto POINT = shade(surface.sample<false, true, false>(P, THETA, 7));
    const auto COLOR = aa.sample(
        surface, P, DX, DY, THETA, 7,
        [&](const Surface::Hit &hit) {
          ++calls;
          return shade(hit);
        },
        &stats);
    HS_EXPECT_EQ(stats.edge_pixels, uint64_t{0});
    HS_EXPECT_EQ(stats.outside_intervals, uint64_t{0});
    HS_EXPECT_EQ(calls, 1);
    HS_EXPECT_EQ(COLOR.r, POINT.r);
    Scan::SurfaceAreaStats wide_stats;
    aa.sample(surface, P, {2 * DX.x, 2 * DX.y}, {2 * DY.x, 2 * DY.y}, THETA, 7,
              shade, &wide_stats);
    HS_EXPECT_EQ(wide_stats.edge_pixels, uint64_t{1});
    std::printf("  hidden terminal sector=%d aa=%g reference=%g\n", sector,
                COLOR.r, reference);
    HS_EXPECT_NEAR(COLOR.r, reference, .005);
  }
}

/** @brief AA colors do not fall below the footprint's first-arrival minimum. */
inline void test_radial_surface_footprint_colours() {
  ScratchScope scope(scratch_arena_a);
  const std::array KNOTS{
      Surface::Knot{0.0f, 0.0f, 0.770516276f},
      Surface::Knot{0.0730376616f, 0.0591502786f, 0.680724382f},
      Surface::Knot{0.157304168f, 0.102900043f, 0.579285026f},
      Surface::Knot{0.251646489f, 0.127174452f, 0.467933953f},
      Surface::Knot{0.35333246f, 0.128406689f, 0.348576486f},
      Surface::Knot{0.458264321f, 0.104044087f, 0.22325477f},
      Surface::Knot{0.560115099f, 0.0528703257f, 0.0941130817f},
      Surface::Knot{0.644939125f, -0.0236403141f, -0.0366389155f},
      Surface::Knot{0.70989567f, -0.119494982f, -0.16676411f},
      Surface::Knot{0.754685581f, -0.228528902f, -0.294035733f},
      Surface::Knot{0.779946387f, -0.344825029f, -0.416276395f},
      Surface::Knot{0.787609756f, -0.462950408f, -0.531394541f},
      Surface::Knot{0.780875385f, -0.578270555f, -0.637420356f},
      Surface::Knot{0.763959527f, -0.687182665f, -0.732539713f},
      Surface::Knot{0.741728902f, -0.787192047f, -0.815124989f},
      Surface::Knot{0.719308317f, -0.876819968f, -0.883763433f},
      Surface::Knot{0.701726139f, -0.955361187f, -0.937280357f},
      Surface::Knot{0.693628013f, -1.02254701f, -0.974760115f},
      Surface::Knot{0.699061275f, -1.0781678f, -0.995561481f},
      Surface::Knot{0.721308768f, -1.12171459f, -0.999328554f},
      Surface::Knot{0.762743473f, -1.15209341f, -0.985996842f},
      Surface::Knot{0.824676335f, -1.16745102f, -0.955794394f},
      Surface::Knot{0.907189548f, -1.16514707f, -0.90923804f},
      Surface::Knot{1.00896096f, -1.14189422f, -0.847124398f},
      Surface::Knot{1.12712681f, -1.09406877f, -0.770516217f},
  };
  Surface surface;
  surface.prepare(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, KNOTS);
  Scan::RadialSurfaceEdgeAA<Surface> aa;
  const Surface::Point P{-.700527012f, -.636006296f};
  const Surface::Point DX{.00336223631f, .00599108543f};
  const Surface::Point DY{-.0201944206f, .0136090377f};
  const auto shade = arrival_ramp(24, true);
  double minimum = 1;
  for (int y = 0; y <= 64; ++y)
    for (int x = 0; x <= 64; ++x) {
      const double U = x / 64.0 - .5, V = y / 64.0 - .5;
      const double hit =
          reference_hit({P.x + U * DX.x + V * DY.x, P.y + U * DX.y + V * DY.y},
                        KNOTS, 7, Plot::STAR_INNER_RATIO);
      HS_EXPECT_GT(hit, 0.0);
      minimum = std::min(minimum, 1 - hit / 24);
    }
  Scan::SurfaceAreaStats stats;
  const auto result = aa.sample(surface, P, DX, DY,
                                Scan::radial_atan2(P.y, P.x), 0, shade, &stats);
  HS_EXPECT_EQ(stats.edge_pixels, uint64_t{1});
  HS_EXPECT_GE(result.r, minimum - .00001);
  std::printf("  footprint first-arrival range: aa=%.9f minimum=%.9f\n",
              result.r, minimum);
}

/** @brief Runs the generic point and edge-coverage oracles. */
inline void run_radial_surface_checks() {
  HS_RUN(test_radial_surface_tip_sweep);
  HS_RUN(test_radial_surface_descending_root_cancellation);
  HS_RUN(test_radial_surface_fraction_gradient);
  HS_RUN(test_radial_surface_degenerate_intervals);
  HS_RUN(test_radial_surface_block_bounds);
  HS_RUN(test_radial_surface_matches_all_edges);
  HS_RUN(test_radial_surface_atan2);
  HS_RUN(test_radial_surface_footprint_witness);
  HS_RUN(test_radial_surface_edge_aa);
  HS_RUN(test_radial_surface_fold_visibility);
  HS_RUN(test_radial_surface_fixture_edge_trace);
  HS_RUN(test_radial_surface_edge_pair_corner);
  HS_RUN(test_radial_surface_terminal_background);
  HS_RUN(test_radial_surface_hidden_terminal);
  HS_RUN(test_radial_surface_footprint_colours);
}
} // namespace hs_test::radial_surface_tests
