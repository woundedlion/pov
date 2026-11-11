/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file test_radial_surface_mask.h
 * @brief Conservative boundary-mask coverage and work census. */

#include <array>
#include <bit>
#include "core/render/scan/radial_surface.h"
#include "core/render/scan/radial_surface_candidate_mask.h"
#include "core/render/scan/radial_surface_edge_aa.h"
#include "core/render/scan/radial_surface_renderer.h"
#include "tests/test_harness.h"
#include "tests/test_radial_surface.h"

namespace hs_test::radial_surface_mask_tests {

using AreaTestSurface = Scan::RadialSurface<4, 7>;
inline constexpr std::array NEAR_AREA_KNOTS{
    AreaTestSurface::Knot{0, 0, 0}, AreaTestSurface::Knot{.2f, 0, 0},
    AreaTestSurface::Knot{.5f, .1f, .2f}, AreaTestSurface::Knot{.6f, .3f, .5f},
    AreaTestSurface::Knot{.8f, .1f, .1f}};
inline constexpr std::array FAR_AREA_KNOTS{
    AreaTestSurface::Knot{0, 0, 0}, AreaTestSurface::Knot{.1f, 0, 0},
    AreaTestSurface::Knot{.3f, -.1f, -.3f},
    AreaTestSurface::Knot{.5f, -.2f, -.4f}, AreaTestSurface::Knot{.6f, 0, 0}};

/** @brief Knot radius profile: monotone expansion or folded radius. */
enum class MaskProfile { EXPANDING, FOLDED };

/** @brief Checks padded segment coverage independently of mask rectangles. */
template <int RADIAL, int ANGULAR, int INTERVALS = 24,
          MaskProfile PROFILE = MaskProfile::EXPANDING>
inline void check_radial_surface_mask(float padding) {
  using Surface = Scan::RadialSurface<INTERVALS, 9>;
  constexpr int SIDES = PROFILE == MaskProfile::EXPANDING ? 7 : 9;
  using Mask = Scan::RadialSurfaceCandidateMask<RADIAL, ANGULAR>;
  using Boundary = Scan::RadialSurfaceBoundaryCandidates;
  using Point = Boundary::Point;
  const float pad = padding;
  uint32_t probes = 0, misses = 0, marked = 0, total = 0;
  uint64_t interval_bits = 0;
  uint32_t candidates = 0, pieces = 0, writes = 0;
  for (int phase = 0; phase < 6; ++phase) {
    ScratchScope scope(scratch_arena_a);
    std::array<typename Surface::Knot, INTERVALS + 1> knots;
    for (int i = 0; i <= INTERVALS; ++i) {
      const float q = float(i) / INTERVALS;
      const float r = PROFILE == MaskProfile::EXPANDING
                          ? math::PI_F * .5f * (1.44f * q - .44f * q * q)
                          : .7f + .25f * sinf(math::TWO_PI_F * q);
      const float p =
          (PROFILE == MaskProfile::EXPANDING ? 1 : 3) *
          sinf((PROFILE == MaskProfile::EXPANDING ? math::PI_F
                                                  : math::TWO_PI_F) *
                   q +
               phase * math::PI_F / 3);
      knots[i] = {r * cosf(p), r * sinf(p), p};
    }
    float maximum_radius = 0;
    for (const auto &knot : knots)
      maximum_radius = std::max(maximum_radius, hypotf(knot.real, knot.imag));
    Surface surface;
    surface.prepare(scratch_arena_a, SIDES, .382f, knots);
    Mask mask;
    mask.bind(scratch_arena_a);
    const auto stats = mask.prepare(surface, pad);
    candidates += stats.candidates;
    pieces += stats.pieces;
    writes += stats.row_writes;
    Boundary::visit_candidates(surface, [&](const auto &candidate) {
      for (int along = 0; along <= 32; ++along) {
        const float t = along / 32.0f;
        const Point p{candidate.a.x + t * (candidate.b.x - candidate.a.x),
                      candidate.a.y + t * (candidate.b.y - candidate.a.y)};
        for (int angle = 0; angle < 16; ++angle) {
          const float a = angle * math::TWO_PI_F / 16;
          const Point q{p.x + pad * cosf(a), p.y + pad * sinf(a)};
          const float r = sqrtf(q.x * q.x + q.y * q.y);
          const float theta = atan2f(q.y, q.x);
          ++probes;
          misses += (mask.interval_bits(r, theta) &
                     (uint32_t(1) << (candidate.origin_interval >>
                                      mask.interval_shift()))) == 0;
          ++probes;
          misses += (mask.interval_bits(r, theta + math::TWO_PI_F / SIDES) &
                     (uint32_t(1) << (candidate.origin_interval >>
                                      mask.interval_shift()))) == 0;
        }
      }
    });
    for (int r = 0; r < 100; ++r)
      for (int a = 0; a < 100; ++a) {
        const float r_value = (r + .5f) * (maximum_radius + pad) / 100;
        const float theta = (a + .5f) * math::TWO_PI_F / (SIDES * 100);
        const uint32_t bits = mask.interval_bits(r_value, theta);
        marked += bits != 0;
        interval_bits += std::popcount(bits);
        ++total;
      }
  }
  HS_EXPECT_EQ(misses, 0u);
  if constexpr (PROFILE == MaskProfile::EXPANDING)
    HS_EXPECT_LT(marked, total * 3 / 5);
  else
    HS_EXPECT_LT(marked, total * 7 / 8);
  HS_EXPECT_GT(probes, 0u);
  HS_EXPECT_GT(candidates, 0u);
  HS_EXPECT_LT(interval_bits, static_cast<uint64_t>(marked) * 4u);
  std::printf("  mean interval groups per marked cell: %.6f\n",
              static_cast<double>(interval_bits) / marked);
  HS_EXPECT_LE(candidates, 6u * (3u * INTERVALS + 2u));
  std::printf(
      "  mask %dx%d N=%d pad=%g: probes=%u misses=%u marked=%u/%u candidates=%u pieces=%u row_writes=%u bytes=%zu\n",
      RADIAL, ANGULAR, INTERVALS, pad, probes, misses, marked, total,
      candidates, pieces, writes, Mask::STORAGE_BYTES);
}

/** @brief Covers candidate births throughout a certified temporal cache. */
inline void test_radial_surface_mask_motion() {
  using Surface = Scan::RadialSurface<24, 7>;
  using Mask = Scan::RadialSurfaceCandidateMask<64, 16>;
  using Boundary = Scan::RadialSurfaceBoundaryCandidates;
  constexpr float PIXEL = .02425f, DRIFT = .006f, RADIUS = math::PI_F * .5f;
  uint32_t probes = 0, misses = 0;
  for (int fixed = 0; fixed < 2; ++fixed)
    for (int phase = 0; phase < 8; ++phase) {
      ScratchScope scope(scratch_arena_a);
      std::array<Surface::Knot, 25> reference, moved;
      for (int i = 0; i <= 24; ++i) {
        const float q = float(i) / 24;
        const float r = RADIUS * (1.44f * q - .44f * q * q);
        const float p = sinf(math::PI_F * q + phase * math::PI_F / 4);
        reference[i] = {r * cosf(p), r * sinf(p), p};
      }
      Surface surface;
      surface.prepare(scratch_arena_a, 7, .382f, reference);
      Mask mask;
      mask.bind(scratch_arena_a);
      mask.prepare(surface, PIXEL, (fixed ? 1 : 3) * RADIUS * DRIFT + 1e-6f,
                   fixed ? DRIFT * math::PI_F / 24 : -1);
      for (int time = -4; time <= 4; ++time) {
        const float t = time / 4.0f;
        for (int i = 0; i <= 24; ++i) {
          const float q = float(i) / 24;
          float r = RADIUS * (1.44f * q - .44f * q * q);
          float p;
          if (fixed)
            p = sinf(math::PI_F * q + phase * math::PI_F / 4 + DRIFT * t);
          else {
            r *= 1 + DRIFT * t * sinf(i * .7f);
            p = reference[i].phase + DRIFT * t * cosf(i * .8f);
          }
          moved[i] = {r * cosf(p), r * sinf(p), p};
        }
        surface.set_knots(moved);
        Boundary::visit_candidates(surface, [&](const auto &candidate) {
          for (int along = 0; along <= 16; ++along)
            for (int angle = 0; angle < 8; ++angle) {
              const float u = along / 16.0f, a = angle * math::TWO_PI_F / 8;
              const float x = candidate.a.x +
                              u * (candidate.b.x - candidate.a.x) +
                              PIXEL * cosf(a);
              const float y = candidate.a.y +
                              u * (candidate.b.y - candidate.a.y) +
                              PIXEL * sinf(a);
              ++probes;
              const auto bits =
                  mask.interval_bits(sqrtf(x * x + y * y), atan2f(y, x));
              misses += ((bits >> (candidate.origin_interval >>
                                   mask.interval_shift())) &
                         1u) == 0;
            }
        });
      }
    }
  HS_EXPECT_GT(probes, 100000u);
  HS_EXPECT_EQ(misses, 0u);
  std::printf("  temporal mask: probes=%u misses=%u\n", probes, misses);
}

/** @brief Candidate filtering preserves two-color area samples. */
inline void test_radial_surface_candidate_parity() {
  using Surface = Scan::RadialSurface<24, 7>;
  using Mask = Scan::RadialSurfaceCandidateMask<64, 16>;
  uint32_t tested = 0, different = 0;
  for (int phase = 0; phase < 8; ++phase) {
    ScratchScope scope(scratch_arena_a);
    std::array<Surface::Knot, 25> knots;
    for (int i = 0; i <= 24; ++i) {
      const float q = float(i) / 24, r = math::PI_F * .5f * q;
      const float p = sinf(math::PI_F * q + phase * math::PI_F / 4);
      knots[i] = {r * cosf(p), r * sinf(p), p};
    }
    Surface surface;
    surface.prepare(scratch_arena_a, 7, .382f, knots);
    Mask mask;
    mask.bind(scratch_arena_a);
    mask.prepare(surface, .025f, .006f * math::PI_F * .5f + 1e-6f,
                 .006f * math::PI_F / 24);
    if (phase & 1) {
      for (int i = 0; i <= 24; ++i) {
        const float q = float(i) / 24, r = math::PI_F * .5f * q;
        const float p = sinf(math::PI_F * q + phase * math::PI_F / 4 + .006f);
        knots[i] = {r * cosf(p), r * sinf(p), p};
      }
      surface.set_knots(knots);
    }
    Scan::RadialSurfaceEdgeAA<Surface> aa;
    const auto shade = radial_surface_tests::arrival_ramp(24);
    for (int y = -30; y <= 30; ++y)
      for (int x = -30; x <= 30; ++x) {
        const Surface::Point p{x * .05f, y * .05f}, DX{.025f, .001f},
            DY{-.002f, .025f};
        const float r = sqrtf(p.x * p.x + p.y * p.y),
                    a = Scan::radial_atan2(p.y, p.x);
        const auto reference =
            aa.template sample<false>(surface, p, DX, DY, a, 0, shade);
        const auto masked = aa.template sample<false>(
            surface, p, DX, DY, a, 0, shade, nullptr,
            {.radius = -1,
             .extent = -1,
             .candidate_bits = mask.interval_bits(r, a),
             .candidate_shift = mask.interval_shift()});
        different += std::bit_cast<uint32_t>(reference.r) !=
                     std::bit_cast<uint32_t>(masked.r);
        ++tested;
      }
  }
  HS_EXPECT_EQ(different, 0u);
  std::printf("  candidate mask parity: pixels=%u differences=%u\n", tested,
              different);
}

/** @brief Retains newborn tip intervals when reference knots coincide. */
inline void test_radial_surface_mask_degenerate_motion() {
  using Surface = Scan::RadialSurface<4, 7>;
  using Mask = Scan::RadialSurfaceCandidateMask<64, 16>;
  using Boundary = Scan::RadialSurfaceBoundaryCandidates;
  ScratchScope scope(scratch_arena_a);
  std::array knots{Surface::Knot{0, 0, 0}, Surface::Knot{.6f, 0, 0},
                   Surface::Knot{.5f, 0, 0}, Surface::Knot{.5f, 0, 0},
                   Surface::Knot{.8f, 0, 0}};
  Surface surface;
  surface.prepare(scratch_arena_a, 7, .382f, knots);
  Mask mask;
  mask.bind(scratch_arena_a);
  mask.prepare(surface, .0005f, .0011f, .004f);
  knots[2] = {.5f * cosf(.002f), -.5f * sinf(.002f), -.002f};
  knots[3] = {.5f * cosf(.002f), .5f * sinf(.002f), .002f};
  surface.set_knots(knots);
  uint32_t probes = 0, misses = 0;
  Boundary::visit_candidates(surface, [&](const auto &candidate) {
    if (candidate.origin_interval != 2 || candidate.knot)
      return;
    for (int i = 0; i <= 16; ++i) {
      const float t = i / 16.0f;
      const float x = candidate.a.x + t * (candidate.b.x - candidate.a.x);
      const float y = candidate.a.y + t * (candidate.b.y - candidate.a.y);
      ++probes;
      misses +=
          (mask.interval_bits(sqrtf(x * x + y * y), atan2f(y, x)) & 4u) == 0;
    }
  });
  HS_EXPECT_GT(probes, 0u);
  HS_EXPECT_EQ(misses, 0u);
}

/** @brief Mask routing defers footprints and preserves point and edge-AA results. */
inline void test_radial_surface_area() {
  using Surface = Scan::RadialSurface<4, 7>;
  using Mask = Scan::RadialSurfaceCandidateMask<64, 16>;
  ScratchScope scope(scratch_arena_a);
  Surface near, far;
  near.prepare(scratch_arena_a, 7, .382f, NEAR_AREA_KNOTS);
  far.prepare(scratch_arena_a, 7, .382f, FAR_AREA_KNOTS);
  Mask near_mask, far_mask;
  near_mask.bind(scratch_arena_a);
  far_mask.bind(scratch_arena_a);
  near_mask.prepare(near, .02f);
  far_mask.prepare(far, .02f);
  Scan::RadialSurfaceArea<Surface, true, true, Mask> area(near, near_mask,
                                                          far_mask);
  Scan::RadialSurfaceArea<Surface, false, false, Mask> unmasked(near, near_mask,
                                                                far_mask);
  Scan::RadialSurfaceEdgeAA<Surface> aa;
  auto lower_bound = [](float) { return 0; };
  auto shade = [](const Surface::Hit &hit) {
    const float q = hit.interval + hit.fraction;
    return Scan::SurfaceColor{q, q * 2, q * 3};
  };
  uint32_t points = 0, edges = 0, hemispheres = 0;
  for (const auto *surface : {&near, &far}) {
    const auto &mask = surface == &near ? near_mask : far_mask;
    for (int r = 1; r <= 100; ++r)
      for (int a = 0; a < 16; ++a) {
        const float r_value = r * .01f, theta = a * math::TWO_PI_F / 112;
        const Surface::Point p{r_value * cosf(theta), r_value * sinf(theta)};
        auto footprint = [&] {
          return Scan::radial_chart_footprint<Surface::Point>(
              {sinf(r_value) * cosf(theta), cosf(r_value),
               sinf(r_value) * sinf(theta)},
              {.01f, 0, 0}, {0, 0, .01f});
        };
        const uint32_t bits = mask.interval_bits(r_value, theta);
        const auto before = Scan::RadialChartWork::footprints;
        Scan::SurfaceAreaStats stats;
        const auto actual = area(*surface, p, theta, r_value, footprint,
                                 lower_bound, shade, &stats);
        HS_EXPECT_EQ(Scan::RadialChartWork::footprints, before + (bits != 0));
        HS_EXPECT_EQ(stats.mask_queries, 1u);
        HS_EXPECT_EQ(stats.mask_hits, uint64_t(bits != 0));
        Scan::SurfaceColor expected;
        if (!bits) {
          ++points;
          Surface::Stats work;
          const auto hit =
              surface->sample<true, true, false>(p, theta, 0, &work);
          expected = hit.found() ? shade(hit) : Scan::SurfaceColor{};
          HS_EXPECT_EQ(stats.pixels, 1u);
          HS_EXPECT_EQ(stats.intervals, work.intervals);
          HS_EXPECT_EQ(stats.point_intervals, work.intervals);
          HS_EXPECT_EQ(stats.max_intervals, int(work.intervals));
          HS_EXPECT_EQ(stats.max_point_intervals, int(work.intervals));
          HS_EXPECT_EQ(stats.fragments, 0u);
        } else {
          ++edges;
          const auto query = footprint();
          Scan::SurfaceAreaStats reference;
          expected = aa.sample<true>(*surface, p, query.dx, query.dy, theta, 0,
                                     shade, &reference,
                                     {.radius = r_value,
                                      .extent = query.extent,
                                      .candidate_bits = bits,
                                      .candidate_shift = mask.interval_shift(),
                                      .point_first = 0});
          HS_EXPECT_EQ(stats.intervals, reference.intervals);
          HS_EXPECT_EQ(stats.point_intervals, reference.point_intervals);
          HS_EXPECT_EQ(stats.outside_intervals, reference.outside_intervals);
          HS_EXPECT_EQ(stats.fragments, reference.fragments);
          HS_EXPECT_EQ(stats.edge_pixels, reference.edge_pixels);
          HS_EXPECT_EQ(stats.hidden_pixels, reference.hidden_pixels);
        }
        HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual.r),
                     std::bit_cast<uint32_t>(expected.r));
        HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual.g),
                     std::bit_cast<uint32_t>(expected.g));
        HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual.b),
                     std::bit_cast<uint32_t>(expected.b));
        const auto unmasked_before = Scan::RadialChartWork::footprints;
        const auto all = unmasked(*surface, p, theta, r_value, footprint,
                                  lower_bound, shade, nullptr);
        HS_EXPECT_EQ(Scan::RadialChartWork::footprints, unmasked_before + 1);
        const auto query = footprint();
        const auto reference = aa.sample<false>(*surface, p, query.dx, query.dy,
                                                theta, 0, shade, nullptr,
                                                {.radius = r_value,
                                                 .extent = query.extent,
                                                 .candidate_bits = UINT32_MAX,
                                                 .candidate_shift = 0,
                                                 .point_first = 0});
        HS_EXPECT_EQ(std::bit_cast<uint32_t>(all.r),
                     std::bit_cast<uint32_t>(reference.r));
      }
    ++hemispheres;
    HS_EXPECT_GT(points, 0u);
    HS_EXPECT_GT(edges, 0u);
    points = edges = 0;
  }
  HS_EXPECT_EQ(hemispheres, 2u);
}

/** @brief Radius-based first intervals preserve all boundaries touching a wide footprint. */
inline void test_radial_surface_area_footprint_first_interval() {
  using Surface = AreaTestSurface;
  using Mask = Scan::RadialSurfaceCandidateMask<64, 16>;
  ScratchScope scope(scratch_arena_a);
  size_t positive = 0;
  for (const auto &knots : {NEAR_AREA_KNOTS, FAR_AREA_KNOTS}) {
    Surface surface;
    surface.prepare(scratch_arena_a, 7, .382f, knots);
    Mask mask;
    Scan::RadialSurfaceArea<Surface, false, false, Mask> area(surface, mask,
                                                              mask);
    auto first_interval = [&](float radius) {
      int first = 0;
      for (int index = 1; index < static_cast<int>(knots.size()); ++index)
        if (hypotf(knots[index].real, knots[index].imag) <= radius)
          first = index;
      return first;
    };
    auto zero_bound = [](float) { return 0; };
    const auto shade = [](const auto &hit) {
      const float value = hit.interval + hit.fraction;
      return Scan::SurfaceColor{value, value * 2, value * 3};
    };
    for (int radial = 1; radial <= 100; ++radial)
      for (int angular = 0; angular < 16; ++angular) {
        const float radius = radial * .01f;
        const float theta = angular * math::TWO_PI_F / 112;
        const Surface::Point point{radius * cosf(theta), radius * sinf(theta)};
        const auto footprint = [&] {
          return Scan::radial_chart_footprint<Surface::Point>(
              {sinf(radius) * cosf(theta), cosf(radius),
               sinf(radius) * sinf(theta)},
              {.06f, 0, 0}, {0, 0, .06f});
        };
        const auto bounded = area(surface, point, theta, radius, footprint,
                                  first_interval, shade, nullptr);
        const auto reference = area(surface, point, theta, radius, footprint,
                                    zero_bound, shade, nullptr);
        HS_EXPECT_EQ(bounded.r > 0, reference.r > 0);
        HS_EXPECT_EQ(bounded.r, reference.r);
        HS_EXPECT_EQ(bounded.g, reference.g);
        HS_EXPECT_EQ(bounded.b, reference.b);
        positive += reference.r > 0;
      }
  }
  HS_EXPECT_GT(positive, size_t{0});
}

/** @brief Checks conservative interval masks and temporal candidate parity. */
inline void run_radial_surface_mask_checks() {
  HS_RUN((check_radial_surface_mask<64, 32>), .035f);
  HS_RUN((check_radial_surface_mask<64, 32>), .02425f);
  HS_RUN((check_radial_surface_mask<64, 32, 96, MaskProfile::FOLDED>), .02425f);
  HS_RUN((check_radial_surface_mask<64, 16, 24>), .02425f);
  HS_RUN((check_radial_surface_mask<64, 16, 96, MaskProfile::FOLDED>), .02425f);
  HS_RUN(test_radial_surface_mask_motion);
  HS_RUN(test_radial_surface_candidate_parity);
  HS_RUN(test_radial_surface_mask_degenerate_motion);
  HS_RUN(test_radial_surface_area);
  HS_RUN(test_radial_surface_area_footprint_first_interval);
}

} // namespace hs_test::radial_surface_mask_tests
