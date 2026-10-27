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
#include "tests/test_harness.h"

namespace hs_test::radial_surface_mask_tests {

/** @brief Checks padded segment coverage independently of mask rectangles. */
template <int RADIAL, int ANGULAR, int INTERVALS = 24>
inline void check_radial_surface_mask(float padding) {
  using Surface = Scan::RadialSurface<INTERVALS, 9>;
  constexpr int SIDES = INTERVALS == 24 ? 7 : 9;
  using Mask = Scan::RadialSurfaceCandidateMask<RADIAL, ANGULAR>;
  using Boundary = Scan::RadialSurfaceBoundaryCandidates;
  using Point = Boundary::Point;
  const float PAD = padding;
  uint32_t probes = 0, misses = 0, marked = 0, total = 0;
  uint32_t candidates = 0, pieces = 0, writes = 0;
  for (int phase = 0; phase < 6; ++phase) {
    ScratchScope scope(scratch_arena_a);
    std::array<typename Surface::Knot, INTERVALS + 1> knots;
    for (int i = 0; i <= INTERVALS; ++i) {
      const float Q = float(i) / INTERVALS;
      const float R = INTERVALS == 24
                          ? math::PI_F * .5f * (1.44f * Q - .44f * Q * Q)
                          : .7f + .25f * sinf(math::TWO_PI_F * Q);
      const float P = (INTERVALS == 24 ? 1 : 3) *
                      sinf((INTERVALS == 24 ? math::PI_F : math::TWO_PI_F) * Q +
                           phase * math::PI_F / 3);
      knots[i] = {R * cosf(P), R * sinf(P), P};
    }
    Surface surface;
    surface.prepare(scratch_arena_a, SIDES, .382f, knots);
    Mask mask;
    mask.bind(scratch_arena_a);
    const auto stats = mask.prepare(surface, PAD);
    candidates += stats.candidates;
    pieces += stats.pieces;
    writes += stats.row_writes;
    Boundary::visit_candidates_detailed(surface, [&](const auto &candidate) {
      for (int along = 0; along <= 32; ++along) {
        const float T = along / 32.0f;
        const Point P{candidate.a.x + T * (candidate.b.x - candidate.a.x),
                      candidate.a.y + T * (candidate.b.y - candidate.a.y)};
        for (int angle = 0; angle < 16; ++angle) {
          const float A = angle * math::TWO_PI_F / 16;
          const Point Q{P.x + PAD * cosf(A), P.y + PAD * sinf(A)};
          const float R = sqrtf(Q.x * Q.x + Q.y * Q.y);
          const float THETA = atan2f(Q.y, Q.x);
          ++probes;
          misses += (mask.interval_bits(R, THETA) &
                     (uint32_t(1) << (candidate.origin_interval >>
                                      mask.interval_shift()))) == 0;
          ++probes;
          misses += (mask.interval_bits(R, THETA + math::TWO_PI_F / SIDES) &
                     (uint32_t(1) << (candidate.origin_interval >>
                                      mask.interval_shift()))) == 0;
        }
      }
    });
    for (int r = 0; r < 100; ++r)
      for (int a = 0; a < 100; ++a) {
        const float R = (r + .5f) * (math::PI_F * .5f + PAD) / 100;
        const float THETA = (a + .5f) * math::TWO_PI_F / (SIDES * 100);
        marked += mask.near_boundary(R, THETA);
        ++total;
      }
  }
  HS_EXPECT_EQ(misses, 0u);
  HS_EXPECT_LT(marked, total);
  HS_EXPECT_LE(candidates, 6u * (3u * INTERVALS + 2u));
  std::printf(
      "  mask %dx%d N=%d pad=%g: probes=%u misses=%u marked=%u/%u candidates=%u pieces=%u row_writes=%u bytes=%zu\n",
      RADIAL, ANGULAR, INTERVALS, PAD, probes, misses, marked, total,
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
        const float Q = float(i) / 24;
        const float R = RADIUS * (1.44f * Q - .44f * Q * Q);
        const float P = sinf(math::PI_F * Q + phase * math::PI_F / 4);
        reference[i] = {R * cosf(P), R * sinf(P), P};
      }
      Surface surface;
      surface.prepare(scratch_arena_a, 7, .382f, reference);
      Mask mask;
      mask.bind(scratch_arena_a);
      mask.prepare(surface, PIXEL, (fixed ? 1 : 3) * RADIUS * DRIFT + 1e-6f,
                   fixed ? DRIFT * math::PI_F / 24 : -1);
      for (int time = -4; time <= 4; ++time) {
        const float T = time / 4.0f;
        for (int i = 0; i <= 24; ++i) {
          const float Q = float(i) / 24;
          float R = RADIUS * (1.44f * Q - .44f * Q * Q);
          float P;
          if (fixed)
            P = sinf(math::PI_F * Q + phase * math::PI_F / 4 + DRIFT * T);
          else {
            R *= 1 + DRIFT * T * sinf(i * .7f);
            P = reference[i].phase + DRIFT * T * cosf(i * .8f);
          }
          moved[i] = {R * cosf(P), R * sinf(P), P};
        }
        surface.set_knots(moved);
        Boundary::visit_candidates_detailed(
            surface, [&](const auto &candidate) {
              for (int along = 0; along <= 16; ++along)
                for (int angle = 0; angle < 8; ++angle) {
                  const float U = along / 16.0f, A = angle * math::TWO_PI_F / 8;
                  const float X = candidate.a.x +
                                  U * (candidate.b.x - candidate.a.x) +
                                  PIXEL * cosf(A);
                  const float Y = candidate.a.y +
                                  U * (candidate.b.y - candidate.a.y) +
                                  PIXEL * sinf(A);
                  ++probes;
                  misses +=
                      !mask.near_boundary(sqrtf(X * X + Y * Y), atan2f(Y, X));
                }
            });
      }
    }
  HS_EXPECT_GT(probes, 100000u);
  HS_EXPECT_EQ(misses, 0u);
  std::printf("  temporal mask: probes=%u misses=%u\n", probes, misses);
}

/** @brief Candidate filtering preserves the reviewed two-color sampler. */
inline void test_radial_surface_candidate_parity() {
  using Surface = Scan::RadialSurface<24, 7>;
  using Mask = Scan::RadialSurfaceCandidateMask<64, 16>;
  uint32_t tested = 0, different = 0;
  for (int phase = 0; phase < 8; ++phase) {
    ScratchScope scope(scratch_arena_a);
    std::array<Surface::Knot, 25> knots;
    for (int i = 0; i <= 24; ++i) {
      const float Q = float(i) / 24, R = math::PI_F * .5f * Q;
      const float P = sinf(math::PI_F * Q + phase * math::PI_F / 4);
      knots[i] = {R * cosf(P), R * sinf(P), P};
    }
    Surface surface;
    surface.prepare(scratch_arena_a, 7, .382f, knots);
    Mask mask;
    mask.bind(scratch_arena_a);
    mask.prepare(surface, .025f, .006f * math::PI_F * .5f + 1e-6f,
                 .006f * math::PI_F / 24);
    if (phase & 1) {
      for (int i = 0; i <= 24; ++i) {
        const float Q = float(i) / 24, R = math::PI_F * .5f * Q;
        const float P = sinf(math::PI_F * Q + phase * math::PI_F / 4 + .006f);
        knots[i] = {R * cosf(P), R * sinf(P), P};
      }
      surface.set_knots(knots);
    }
    Scan::RadialSurfaceEdgeAA<Surface> aa;
    auto shade = [](const auto &hit) {
      const float Q = (hit.interval + hit.fraction) / 24;
      return Scan::SurfaceColor{Q, Q, Q};
    };
    float worst = 0;
    for (int y = -30; y <= 30; ++y)
      for (int x = -30; x <= 30; ++x) {
        const Surface::Point P{x * .05f, y * .05f}, DX{.025f, .001f},
            DY{-.002f, .025f};
        const float R = sqrtf(P.x * P.x + P.y * P.y),
                    A = Scan::radial_atan2(P.y, P.x);
        const auto reference =
            aa.template sample<false>(surface, P, DX, DY, A, 0, shade);
        const auto masked = aa.template sample<false>(
            surface, P, DX, DY, A, 0, shade, nullptr, -1, -1,
            mask.interval_bits(R, A), mask.interval_shift());
        const float ERROR = fabsf(reference.r - masked.r);
        worst = fmaxf(worst, ERROR);
        different += std::bit_cast<uint32_t>(reference.r) !=
                     std::bit_cast<uint32_t>(masked.r);
        ++tested;
      }
    HS_EXPECT_LT(worst, 1e-6f);
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
  Boundary::visit_candidates_detailed(surface, [&](const auto &candidate) {
    if (candidate.origin_interval != 2 || candidate.knot)
      return;
    for (int i = 0; i <= 16; ++i) {
      const float T = i / 16.0f;
      const float X = candidate.a.x + T * (candidate.b.x - candidate.a.x);
      const float Y = candidate.a.y + T * (candidate.b.y - candidate.a.y);
      ++probes;
      misses +=
          (mask.interval_bits(sqrtf(X * X + Y * Y), atan2f(Y, X)) & 4u) == 0;
    }
  });
  HS_EXPECT_GT(probes, 0u);
  HS_EXPECT_EQ(misses, 0u);
}

/** @brief Compares conservative masks at two arena sizes. */
inline void test_radial_surface_mask() {
  check_radial_surface_mask<64, 32>(.035f);
  check_radial_surface_mask<64, 32>(.02425f);
  check_radial_surface_mask<64, 32, 96>(.02425f);
  check_radial_surface_mask<64, 16, 24>(.02425f);
  check_radial_surface_mask<64, 16, 96>(.02425f);
  test_radial_surface_mask_motion();
  test_radial_surface_candidate_parity();
  test_radial_surface_mask_degenerate_motion();
}

} // namespace hs_test::radial_surface_mask_tests
