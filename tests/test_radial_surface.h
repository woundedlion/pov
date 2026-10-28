/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file test_radial_surface.h
 * @brief First-arrival queries against independent swept-edge intersections. */

#include <array>
#include <bit>
#include <cfloat>
#include <vector>
#include <cmath>
#include "core/render/scan/radial_surface.h"
#include "core/render/scan/radial_surface_edge_aa.h"
#include "core/render/scan/radial_surface_math.h"
#include "core/render/scan/radial_surface_renderer.h"
#include "core/render/scan/radial_surface_boundary.h"
#include "tests/test_harness.h"

namespace hs_test::radial_surface_tests {

using Surface = Scan::RadialSurface<64, 9>;

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

/** @brief Finds the first point on any swept edge without angular sectors.
 * @param p Chart point.
 * @param knots Similarity knots, starting at zero scale.
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

/** @brief Static and zero-scale intervals retain their solid membership. */
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
  HS_EXPECT_GT(skipped, 0u);
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
  uint32_t random = 0x42da5819;
  auto coordinate = [&] {
    random = random * 1664525u + 1013904223u;
    return float(int32_t(random)) / 2147483648.0f;
  };
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

/** @brief Two-color jump coverage is finite and bounded. */
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
  auto shade = [](const Surface::Hit &hit) {
    const float Q = (hit.interval + hit.fraction) / N;
    return Scan::SurfaceColor{Q, Q, Q};
  };
  for (int x = -20; x <= 20; ++x)
    for (int y = -20; y <= 20; ++y) {
      const Surface::Point P{x * .05f, y * .05f};
      const auto color =
          aa.sample(surface, P, {.05f, 0}, {0, .05f},
                    Scan::radial_atan2(P.y, P.x), 0, shade, &stats);
      HS_EXPECT_TRUE(std::isfinite(color.r));
      HS_EXPECT_GE(color.r, 0);
      HS_EXPECT_LE(color.r, 1);
      double reference = 0;
      for (int sx = 0; sx < 12; ++sx)
        for (int sy = 0; sy < 12; ++sy) {
          const Surface::Point Q{P.x + .05f * ((sx + .5f) / 12 - .5f),
                                 P.y + .05f * ((sy + .5f) / 12 - .5f)};
          const auto hit = surface.template sample<false, true, false>(
              Q, Scan::radial_atan2(Q.y, Q.x));
          reference += hit.found() ? shade(hit).r : 0;
        }
      reference /= 144;
      const auto center = surface.template sample<false, true, false>(
          P, Scan::radial_atan2(P.y, P.x));
      const double POINT = center.found() ? shade(center).r : 0;
      const float ERROR = float(std::abs(color.r - reference));
      if (x == 8 && y == -18)
        std::printf("  corner replay: aa=%g center=%g reference=%g\n", color.r,
                    POINT, reference);
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
  uint32_t segments = 0, visible = 0;
  float longest = 0, max_gap = 0;
  for (int phase = 0; phase < 64; ++phase) {
    for (int i = 0; i <= N; ++i) {
      const float U = float(i) / N;
      const float R = U < (.5f * math::PI_F / 6 / 1.1278247916f)
                          ? 2 * U * 1.1278247916f
                          : acosf(1.1278247916f * (1 - U));
      const float P = sinf(math::PI_F * U + phase * math::TWO_PI_F / 64);
      knots[i] = {R * cosf(P), R * sinf(P), P};
    }
    surface.set_knots({knots, N + 1});
    HS_EXPECT_TRUE(surface.strictly_increasing_radius());
    Scan::RadialSurfaceBoundaryCandidates::visit_candidates_detailed(
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
            if (!hit.found() ||
                hit.interval + hit.fraction >= candidate.origin_interval - .01f)
              any = true;
            if (any) {
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
  HS_EXPECT_GE(max_gap, 0.69f);
  HS_EXPECT_LE(max_gap, 0.74f);
  std::printf(
      "  radial fold visibility: segments=%u visible=%u longest=%g radial_gap=%g\n",
      segments, visible, longest, max_gap);
}

/** @brief Edge shading matches the area reference at a partially covered pixel. */
inline void test_radial_surface_fixture_edge_trace() {
  ScratchScope scope(scratch_arena_a);
  const std::array KNOTS{
      Surface::Knot{float(0), float(0), float(0.72896862)},
      Surface::Knot{float(0.0270809401), float(0.0225565769), float(0.6945014)},
      Surface::Knot{float(0.0557565615), float(0.0431266911),
                    float(0.658361077)},
      Surface::Knot{float(0.0860152766), float(0.0614895262),
                    float(0.620634735)},
      Surface::Knot{float(0.117813617), float(0.0774259195), float(0.58141315)},
      Surface::Knot{float(0.151075974), float(0.0907218978),
                    float(0.540790915)},
      Surface::Knot{float(0.185694739), float(0.101172216), float(0.498865962)},
      Surface::Knot{float(0.221531272), float(0.108584233), float(0.455739081)},
      Surface::Knot{float(0.258417428), float(0.1127818), float(0.411514282)},
      Surface::Knot{float(0.296157509), float(0.113609128), float(0.366298229)},
      Surface::Knot{float(0.334531397), float(0.110934265), float(0.320199609)},
      Surface::Knot{float(0.373297751), float(0.104652569), float(0.273329675)},
      Surface::Knot{float(0.412198186), float(0.0946896374), float(0.2258012)},
      Surface::Knot{float(0.45096156), float(0.0810035169), float(0.177728802)},
      Surface::Knot{float(0.489308983), float(0.0635868534),
                    float(0.129228204)},
      Surface::Knot{float(0.526937485), float(0.0424659103),
                    float(0.0804163814)},
      Surface::Knot{float(0.562327325), float(0.0176688246),
                    float(0.0314107239)},
      Surface::Knot{float(0.594732523), float(-0.0105103729),
                    float(-0.0176706016)},
      Surface::Knot{float(0.624187768), float(-0.0417007692),
                    float(-0.0667092502)},
      Surface::Knot{float(0.650704622), float(-0.0755498782),
                    float(-0.115587324)},
      Surface::Knot{float(0.67428875), float(-0.111715004),
                    float(-0.164186954)},
      Surface::Knot{float(0.694952309), float(-0.149861753),
                    float(-0.212390959)},
      Surface::Knot{float(0.712721169), float(-0.189662695),
                    float(-0.260083288)},
      Surface::Knot{float(0.72764045), float(-0.230797902),
                    float(-0.307149023)},
      Surface::Knot{float(0.739777207), float(-0.272956789),
                    float(-0.353474855)},
      Surface::Knot{float(0.749223769), float(-0.315839022),
                    float(-0.398949146)},
      Surface::Knot{float(0.756096363), float(-0.359158784),
                    float(-0.443462372)},
      Surface::Knot{float(0.760536909), float(-0.402644694),
                    float(-0.486907214)},
      Surface::Knot{float(0.762710333), float(-0.446044207),
                    float(-0.529179096)},
      Surface::Knot{float(0.762804151), float(-0.489122957),
                    float(-0.570176005)},
      Surface::Knot{float(0.761025906), float(-0.531668007),
                    float(-0.609799385)},
      Surface::Knot{float(0.757600129), float(-0.573488235),
                    float(-0.647953749)},
      Surface::Knot{float(0.75276643), float(-0.61441505), float(-0.684547186)},
      Surface::Knot{float(0.746776521), float(-0.654302001),
                    float(-0.719491422)},
      Surface::Knot{float(0.739891231), float(-0.69302398),
                    float(-0.752702296)},
      Surface::Knot{float(0.732377052), float(-0.730477929),
                    float(-0.784099877)},
      Surface::Knot{float(0.724504232), float(-0.766579449),
                    float(-0.813608468)},
      Surface::Knot{float(0.716543794), float(-0.801262558),
                    float(-0.841157019)},
      Surface::Knot{float(0.708764195), float(-0.834477246),
                    float(-0.866679192)},
      Surface::Knot{float(0.7014305), float(-0.86618638), float(-0.890113413)},
      Surface::Knot{float(0.69480139), float(-0.896364093),
                    float(-0.911403358)},
      Surface::Knot{float(0.689127088), float(-0.924992561),
                    float(-0.930497527)},
      Surface::Knot{float(0.684648991), float(-0.952058792),
                    float(-0.947350144)},
      Surface::Knot{float(0.681596339), float(-0.977553189),
                    float(-0.961920381)},
      Surface::Knot{float(0.680185676), float(-1.00146544),
                    float(-0.974173427)},
      Surface::Knot{float(0.68062067), float(-1.0237813), float(-0.984079599)},
      Surface::Knot{float(0.68309021), float(-1.04448164), float(-0.991614938)},
      Surface::Knot{float(0.687765598), float(-1.06353951),
                    float(-0.996761382)},
      Surface::Knot{float(0.694802642), float(-1.08091712),
                    float(-0.999506593)},
      Surface::Knot{float(0.704337895), float(-1.09656513),
                    float(-0.999843895)},
      Surface::Knot{float(0.716490448), float(-1.11041999),
                    float(-0.997772455)},
      Surface::Knot{float(0.731357992), float(-1.12240326),
                    float(-0.993297338)},
      Surface::Knot{float(0.749015927), float(-1.13242221),
                    float(-0.986429214)},
      Surface::Knot{float(0.769518554), float(-1.14036703),
                    float(-0.977184832)},
      Surface::Knot{float(0.792895377), float(-1.14611232),
                    float(-0.965586185)},
      Surface::Knot{float(0.819150507), float(-1.14951777),
                    float(-0.951661408)},
      Surface::Knot{float(0.848262131), float(-1.15042794),
                    float(-0.935443997)},
      Surface::Knot{float(0.880178809), float(-1.1486764), float(-0.916973054)},
      Surface::Knot{float(0.914821804), float(-1.14408469),
                    float(-0.896293044)},
      Surface::Knot{float(0.95208174), float(-1.1364665), float(-0.873453617)},
      Surface::Knot{float(0.991816998), float(-1.12563121),
                    float(-0.848510206)},
      Surface::Knot{float(1.03385639), float(-1.11138487), float(-0.821522474)},
      Surface::Knot{float(1.07799459), float(-1.09353769), float(-0.792555749)},
      Surface::Knot{float(1.12399709), float(-1.07190359), float(-0.76167953)},
      Surface::Knot{float(1.171597), float(-1.0463084), float(-0.728968561)}};
  Surface surface;
  surface.prepare(scratch_arena_a, 7, .382f, KNOTS);
  const Surface::Point P{-.435921282f, -.661391735f},
      DX{.0171549805f, -.00341231585f}, DY{.00249151839f, .0213687979f};
  Scan::RadialSurfaceEdgeAA<Surface> aa;
  auto shade = [](const Surface::Hit &h) {
    const float Q = (h.interval + h.fraction) / 64;
    return Scan::SurfaceColor{Q, Q, Q};
  };
  std::array<Surface::Hit, 2> hits{};
  int count = 0;
  auto capture = [&](const Surface::Hit &h) {
    if (count < 2)
      hits[count++] = h;
    return shade(h);
  };
  const auto result = aa.sample(surface, P, DX, DY, -2.15354013f, 23, capture);
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
  double reference = 0;
  for (int x = 0; x < 32; ++x)
    for (int y = 0; y < 32; ++y) {
      const float U = (x + .5f) / 32 - .5f, V = (y + .5f) / 32 - .5f;
      const Surface::Point Q{P.x + U * DX.x + V * DY.x,
                             P.y + U * DX.y + V * DY.y};
      const auto hit = surface.template sample<false, true, false>(
          Q, Scan::radial_atan2(Q.y, Q.x));
      reference += hit.found() ? shade(hit).r : 0;
    }
  reference /= 1024;
  HS_EXPECT_GT(result.r, 0.0f);
  HS_EXPECT_NEAR(result.r, reference, 0.005);
  std::printf(
      "  fixture edge trace: q=%.8f first0=%.8f reference=%.8f calls=%d h0=%.8f h1=%.8f\n",
      result.r, zero.r, reference, count, hits[0].interval + hits[0].fraction,
      hits[1].interval + hits[1].fraction);
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
  HS_EXPECT_NEAR(color.r, HALF_PLANE * HALF_PLANE, 2e-6);
  HS_EXPECT_LT(color.r, HALF_PLANE - .1f);
}

/** @brief Checks overlapping early and terminal boundaries against subpixel area. */
inline void test_radial_surface_terminal_background() {
  {
    ScratchScope scope(scratch_arena_a);
    const std::array KNOTS{
        Surface::Knot{float(0), float(0), float(0.72896862)},
        Surface::Knot{float(0.075755164), float(0.0556274317),
                      float(0.63338083)},
        Surface::Knot{float(0.162471071), float(0.0945313275),
                      float(0.526955783)},
        Surface::Knot{float(0.258417428), float(0.1127818), float(0.411514282)},
        Surface::Knot{float(0.360347688), float(0.107152462),
                      float(0.289031774)},
        Surface::Knot{float(0.463804126), float(0.075611867),
                      float(0.161603808)},
        Surface::Knot{float(0.562327325), float(0.0176688246),
                      float(0.0314107239)},
        Surface::Knot{float(0.642191768), float(-0.0639929622),
                      float(-0.099319756)},
        Surface::Knot{float(0.70119518), float(-0.162961051),
                      float(-0.228350967)},
        Surface::Knot{float(0.739777207), float(-0.272956789),
                      float(-0.353474855)},
        Surface::Knot{float(0.75931716), float(-0.388147056),
                      float(-0.47255075)},
        Surface::Knot{float(0.762408674), float(-0.503373921),
                      float(-0.583541274)},
        Surface::Knot{float(0.75276643), float(-0.61441505),
                      float(-0.684547186)},
        Surface::Knot{float(0.734934986), float(-0.718139827),
                      float(-0.773840308)},
        Surface::Knot{float(0.713917315), float(-0.812499404),
                      float(-0.849892735)},
        Surface::Knot{float(0.69480139), float(-0.896364093),
                      float(-0.911403358)},
        Surface::Knot{float(0.682441533), float(-0.969230711),
                      float(-0.957319498)},
        Surface::Knot{float(0.681209445), float(-1.03086174),
                      float(-0.986855686)},
        Surface::Knot{float(0.694802642), float(-1.08091712),
                      float(-0.999506593)},
        Surface::Knot{float(0.7260952), float(-1.11862218),
                      float(-0.995055556)},
        Surface::Knot{float(0.776990771), float(-1.14253318),
                      float(-0.97357893)},
        Surface::Knot{float(0.848262131), float(-1.15042794),
                      float(-0.935443997)},
        Surface::Knot{float(0.939379215), float(-1.13935375),
                      float(-0.88130337)},
        Surface::Knot{float(1.04834688), float(-1.10584521),
                      float(-0.812083423)},
        Surface::Knot{float(1.171597), float(-1.0463084), float(-0.728968561)}};
    Surface surface;
    surface.prepare(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, KNOTS);
    const Surface::Point P{float(0.126707062), float(-0.879523873)};
    const Surface::Point DX{float(0.0154463407), float(-0.00674017007)};
    const Surface::Point DY{float(0.0113697099), float(0.0190199912)};
    Scan::RadialSurfaceEdgeAA<Surface> aa;
    auto shade = [](const Surface::Hit &h) {
      const float Q = (h.interval + h.fraction) / 24;
      return Scan::SurfaceColor{Q, Q, Q};
    };
    Scan::SurfaceAreaStats stats;
    const auto RESULT =
        aa.sample(surface, P, DX, DY, float(-1.42771757), 10, shade, &stats);
    double reference = 0;
    for (int y = 0; y < 128; ++y)
      for (int x = 0; x < 128; ++x) {
        const float U = (x + .5f) / 128 - .5f, V = (y + .5f) / 128 - .5f;
        const Surface::Point P_SAMPLE{P.x + U * DX.x + V * DY.x,
                                      P.y + U * DX.y + V * DY.y};
        const auto HIT = surface.template sample<false, true, false>(
            P_SAMPLE, Scan::radial_atan2(P_SAMPLE.y, P_SAMPLE.x));
        if (HIT.found())
          reference += shade(HIT).r;
      }
    reference /= 16384;
    HS_EXPECT_EQ(stats.clips, 1u);
    HS_EXPECT_LT(fabs(RESULT.r - reference), 0.08);
    std::printf("  terminal overlap 20,34: aa=%g reference=%g\n", RESULT.r,
                reference);
  }
  {
    ScratchScope scope(scratch_arena_a);
    const std::array KNOTS{
        Surface::Knot{float(0), float(0), float(0.368096858)},
        Surface::Knot{float(0.0912108868), float(0.0226678196),
                      float(0.243586123)},
        Surface::Knot{float(0.186731219), float(0.0215517394),
                      float(0.114907548)},
        Surface::Knot{float(0.281921297), float(-0.00443695905),
                      float(-0.0157371461)},
        Surface::Knot{float(0.371935785), float(-0.0547344834),
                      float(-0.146112502)},
        Surface::Knot{float(0.452398509), float(-0.127149373),
                      float(-0.273987889)},
        Surface::Knot{float(0.518810213), float(-0.217624038),
                      float(-0.397175252)},
        Surface::Knot{float(0.562117755), float(-0.317063063),
                      float(-0.513566792)},
        Surface::Knot{float(0.58540678), float(-0.418962836),
                      float(-0.621171176)},
        Surface::Knot{float(0.593782067), float(-0.518843412),
                      float(-0.718147039)},
        Surface::Knot{float(0.592395544), float(-0.613423228),
                      float(-0.802835226)},
        Surface::Knot{float(0.586460471), float(-0.700511515),
                      float(-0.873786807)},
        Surface::Knot{float(0.581069171), float(-0.778795123),
                      float(-0.929787457)},
        Surface::Knot{float(0.580974758), float(-0.847539246),
                      float(-0.969879329)},
        Surface::Knot{float(0.590400279), float(-0.906234384),
                      float(-0.993376136)},
        Surface::Knot{float(0.612883449), float(-0.954249144),
                      float(-0.999876201)},
        Surface::Knot{float(0.651133478), float(-0.99053508),
                      float(-0.989267945)},
        Surface::Knot{float(0.706889272), float(-1.01342487),
                      float(-0.961733162)},
        Surface::Knot{float(0.780762732), float(-1.02055955),
                      float(-0.917742789)},
        Surface::Knot{float(0.872070372), float(-1.00897121),
                      float(-0.858049512)},
        Surface::Knot{float(0.978691936), float(-0.975324929),
                      float(-0.783675015)},
        Surface::Knot{float(1.09699905), float(-0.916311085),
                      float(-0.695891261)},
        Surface::Knot{float(1.22190893), float(-0.829155505),
                      float(-0.596200883)},
        Surface::Knot{float(1.34712458), float(-0.712165892),
                      float(-0.48630923)},
        Surface::Knot{float(1.46557486), float(-0.565235972),
                      float(-0.368096769)}};
    Surface surface;
    surface.prepare(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, KNOTS);
    const Surface::Point P{float(-0.349607706), float(-0.830267727)};
    const Surface::Point DX{float(0.0242360439), float(0.00149365724)};
    const Surface::Point DY{float(-0.00394339021), float(0.0215454288)};
    Scan::RadialSurfaceEdgeAA<Surface> aa;
    auto shade = [](const Surface::Hit &h) {
      const float Q = (h.interval + h.fraction) / 24;
      return Scan::SurfaceColor{Q, Q, Q};
    };
    Scan::SurfaceAreaStats stats;
    const auto RESULT =
        aa.sample(surface, P, DX, DY, float(-1.96934056), 10, shade, &stats);
    double reference = 0;
    for (int y = 0; y < 128; ++y)
      for (int x = 0; x < 128; ++x) {
        const float U = (x + .5f) / 128 - .5f, V = (y + .5f) / 128 - .5f;
        const Surface::Point P_SAMPLE{P.x + U * DX.x + V * DY.x,
                                      P.y + U * DX.y + V * DY.y};
        const auto HIT = surface.template sample<false, true, false>(
            P_SAMPLE, Scan::radial_atan2(P_SAMPLE.y, P_SAMPLE.x));
        if (HIT.found())
          reference += shade(HIT).r;
      }
    reference /= 16384;
    HS_EXPECT_EQ(stats.clips, 1u);
    HS_EXPECT_LT(fabs(RESULT.r - reference), 0.13);
    std::printf("  terminal overlap 53,74: aa=%g reference=%g\n", RESULT.r,
                reference);
  }
  {
    ScratchScope scope(scratch_arena_a);
    const std::array KNOTS{Surface::Knot{float(0), float(0), float(0.99361074)},
                           Surface::Knot{float(0.0507930592),
                                         float(0.0790779442),
                                         float(0.99984163)},
                           Surface::Knot{float(0.103300296), float(0.157041624),
                                         float(0.988964915)},
                           Surface::Knot{float(0.161437929), float(0.231164649),
                                         float(0.961166739)},
                           Surface::Knot{float(0.228672355), float(0.298397481),
                                         float(0.916922748)},
                           Surface::Knot{float(0.307668567), float(0.355206162),
                                         float(0.85698998)},
                           Surface::Knot{float(0.399015129), float(0.396624714),
                                         float(0.782393813)},
                           Surface::Knot{float(0.495923966), float(0.412994742),
                                         float(0.694410741)},
                           Surface::Knot{float(0.596352816), float(0.403229833),
                                         float(0.59454608)},
                           Surface::Knot{float(0.697771132), float(0.367275268),
                                         float(0.484508574)},
                           Surface::Knot{float(0.796234846), float(0.305336982),
                                         float(0.366181046)},
                           Surface::Knot{float(0.887061715), float(0.218572244),
                                         float(0.241587907)},
                           Surface::Knot{float(0.965498745), float(0.109432288),
                                         float(0.112861246)},
                           Surface::Knot{float(1.027385), float(-0.0182857905),
                                         float(-0.0177966207)},
                           Surface::Knot{float(1.06974053), float(-0.159651712),
                                         float(-0.148149729)},
                           Surface::Knot{float(1.09120262), float(-0.309022009),
                                         float(-0.275968194)},
                           Surface::Knot{float(1.09224212), float(-0.460588574),
                                         float(-0.39906469)},
                           Surface::Knot{float(1.07513607), float(-0.608937204),
                                         float(-0.515332997)},
                           Surface::Knot{float(1.04372191), float(-0.749517918),
                                         float(-0.622783899)},
                           Surface::Knot{float(1.00298965), float(-0.878943443),
                                         float(-0.719578981)},
                           Surface::Knot{float(0.95860666), float(-0.995072722),
                                         float(-0.804061472)},
                           Surface::Knot{float(0.916440845), float(-1.09689069),
                                         float(-0.874786615)},
                           Surface::Knot{float(0.882160842), float(-1.1842097),
                                         float(-0.930543721)},
                           Surface::Knot{float(0.860919476), float(-1.25727582),
                                         float(-0.970378935)},
                           Surface::Knot{float(0.857132673), float(-1.31633019),
                                         float(-0.99361074)}};
    Surface surface;
    surface.prepare(scratch_arena_a, 7, Plot::STAR_INNER_RATIO, KNOTS);
    const Surface::Point P{float(0.851393998), float(0.705435812)};
    const Surface::Point DX{float(0.0215828046), float(-0.0114967246)};
    const Surface::Point DY{float(0.0069451537), float(0.0210930854)};
    Scan::RadialSurfaceEdgeAA<Surface> aa;
    auto shade = [](const Surface::Hit &h) {
      const float Q = (h.interval + h.fraction) / 24;
      return Scan::SurfaceColor{Q, Q, Q};
    };
    Scan::SurfaceAreaStats stats;
    const auto RESULT =
        aa.sample(surface, P, DX, DY, float(0.691918015), 13, shade, &stats);
    double reference = 0;
    for (int y = 0; y < 128; ++y)
      for (int x = 0; x < 128; ++x) {
        const float U = (x + .5f) / 128 - .5f, V = (y + .5f) / 128 - .5f;
        const Surface::Point P_SAMPLE{P.x + U * DX.x + V * DY.x,
                                      P.y + U * DX.y + V * DY.y};
        const auto HIT = surface.template sample<false, true, false>(
            P_SAMPLE, Scan::radial_atan2(P_SAMPLE.y, P_SAMPLE.x));
        if (HIT.found())
          reference += shade(HIT).r;
      }
    reference /= 16384;
    HS_EXPECT_EQ(stats.clips, 1u);
    HS_EXPECT_LT(fabs(RESULT.r - reference), 0.04);
    std::printf("  terminal overlap 89,88: aa=%g reference=%g\n", RESULT.r,
                reference);
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
    auto shade = [](const Surface::Hit &hit) {
      const float C = 1 - (hit.interval + hit.fraction) / 24;
      return Scan::SurfaceColor{C, C, C};
    };
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

/** @brief AA colors stay within the footprint's first-arrival range. */
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
  auto shade = [](const Surface::Hit &hit) {
    const float C = 1 - (hit.interval + hit.fraction) / 24;
    return Scan::SurfaceColor{C, C, C};
  };
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
