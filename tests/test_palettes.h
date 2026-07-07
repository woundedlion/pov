/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/color/palettes.h — the named generative/OKLCH palette layer.
 *
 * Coverage:
 *   - Named ProceduralPalette endpoints: pinned 16-bit linear colors at t=0/1,
 *     including the cos(0)=1 channels derivable from the cosine coefficients.
 *   - Shortest-arc hue direction across the +/-PI seam, driven by a named
 *     palette whose endpoint hues straddle it.
 *   - MeshPaletteBank: slot 0 reproduces its embers source, every slot bakes a
 *     distinct LUT, and shuffle_indices is a permutation of [0, N).
 */
#pragma once

#include <array>
#include <cmath>

#include "core/color/color.h"
#include "core/color/palettes.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace palettes_tests {

/**
 * @brief Converts a realized Pixel color back to OKLCH for hue inspection.
 * @param c The color to convert.
 * @return The OKLCH of c (via the exact forward transform).
 */
inline OKLCH color_to_oklch(const Color4 &c) {
  float r = c.color.r / 65535.0f, g = c.color.g / 65535.0f,
        b = c.color.b / 65535.0f;
  return oklab_to_oklch(linear_rgb_to_oklab(r, g, b));
}

/**
 * @brief Wraps a hue difference into (-PI, PI] for circular comparison.
 * @param dh Raw hue difference in radians.
 * @return The equivalent difference in (-PI, PI].
 */
inline float wrap_hue(float dh) {
  while (dh > PI_F) dh -= 2.0f * PI_F;
  while (dh < -PI_F) dh += 2.0f * PI_F;
  return dh;
}

/**
 * @brief Pins named ProceduralPalette endpoints against golden 16-bit colors.
 * @details The cosine formula C = a + b*cos(2*PI*(c*t+d)) is exact at the
 *          channels where the argument is an integer multiple of 2*PI: for
 *          darkRainbow (c=1, d_r=0) the red channel reads a+b at both t=0 and
 *          t=1, and for mauveFade (d_r=d_b=0) red and blue saturate at 1.0. The
 *          remaining channels carry fast_cosf, so they are pinned as goldens of
 *          the current deterministic output — any drift in the named constants
 *          or the sRGB->linear interpolation fails here.
 */
inline void test_named_procedural_palette_endpoints() {
  // darkRainbow: c={1,1,1}, d={0,0.33,0.67}. Red arg = 2*PI*(t) hits cos=1 at
  // both ends, so t=0 and t=1 must produce the identical pinned color.
  Color4 dr0 = Palettes::darkRainbow.get(0.0f);
  Color4 dr1 = Palettes::darkRainbow.get(1.0f);
  HS_EXPECT_EQ(dr0.color.r, 47426); // 0.367 + 0.5*cos(0) = 0.867 sRGB
  HS_EXPECT_EQ(dr0.color.g, 954);
  HS_EXPECT_EQ(dr0.color.b, 954);
  HS_EXPECT_EQ(dr1.color.r, dr0.color.r);
  HS_EXPECT_EQ(dr1.color.g, dr0.color.g);
  HS_EXPECT_EQ(dr1.color.b, dr0.color.b);
  HS_EXPECT_NEAR(dr0.alpha, 1.0f, 1e-6f);

  // mauveFade: a_r=0.583 b_r=1, a_b=0.583 b_b=1 -> 1.583 clamped to 1.0 at t=0
  // (cos=1). Green has b_g=0, so it is the constant a_g=0 -> black.
  Color4 mf0 = Palettes::mauveFade.get(0.0f);
  HS_EXPECT_EQ(mf0.color.r, 65535);
  HS_EXPECT_EQ(mf0.color.g, 0);
  HS_EXPECT_EQ(mf0.color.b, 65535);

  // fireGlow / peachPop: fully golden-pinned (fast_cosf channels).
  Color4 fg0 = Palettes::fireGlow.get(0.0f);
  Color4 fg1 = Palettes::fireGlow.get(1.0f);
  HS_EXPECT_EQ(fg0.color.r, 108);
  HS_EXPECT_EQ(fg0.color.g, 0);
  HS_EXPECT_EQ(fg0.color.b, 0);
  HS_EXPECT_EQ(fg1.color.r, 17340);
  HS_EXPECT_EQ(fg1.color.g, 9961);
  HS_EXPECT_EQ(fg1.color.b, 0);

  Color4 pp0 = Palettes::peachPop.get(0.0f);
  HS_EXPECT_EQ(pp0.color.r, 65535);
  HS_EXPECT_EQ(pp0.color.g, 28156);
  HS_EXPECT_EQ(pp0.color.b, 0);
}

/**
 * @brief Verifies hue interpolates along the short arc for a seam-straddling
 *        named-palette pair.
 * @details undersea's endpoints sit at h ~= -2.00 and +2.49 rad: ~1.80 rad
 *          apart the short way (through the +/-PI seam), ~4.48 rad apart the
 *          long way (through 0). lerp_oklch must drive the midpoint across the
 *          seam (a hue on the negative/|h|>PI/2 side), NOT through ~+0.25 where
 *          a naive average of the two angles would land. This is the hue-arc
 *          contract the named generative palettes rely on.
 */
inline void test_named_palette_hue_short_arc() {
  OKLCH a = color_to_oklch(Palettes::undersea.get(0.0f));
  OKLCH b = color_to_oklch(Palettes::undersea.get(1.0f));
  // Precondition: endpoints really do straddle the seam (numeric gap > PI).
  HS_EXPECT_GT(std::fabs(b.h - a.h), PI_F);

  OKLCH mid = lerp_oklch(a, b, 0.5f);
  // The short-arc midpoint is the circular average of the two hues.
  float short_mid = a.h + 0.5f * wrap_hue(b.h - a.h);
  HS_EXPECT_NEAR(wrap_hue(mid.h - short_mid), 0.0f, 1e-3f);
  // ...and it is NOT the naive through-zero average.
  float naive_mid = 0.5f * (a.h + b.h);
  HS_EXPECT_GT(std::fabs(wrap_hue(mid.h - naive_mid)), 1.0f);
}

/**
 * @brief Verifies MeshPaletteBank bakes its sources and indexes them distinctly.
 * @details Slot 0's baked endpoints must match the embers source (the first
 *          entry in sources()), proving the bank baked the right palette into
 *          the right slot; the five slots must produce distinct LUTs so a
 *          slot-mapping bug can't collapse them.
 */
inline void test_mesh_palette_bank_lookup() {
  alignas(std::max_align_t) static uint8_t buf[MeshPaletteBank::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  MeshPaletteBank bank;
  bank.bake_all(arena);

  // Slot 0 == embers (sources()[0]) baked: pinned golden endpoints.
  Color4 s0 = bank[0].get(0.0f);
  Color4 s1 = bank[0].get(1.0f);
  HS_EXPECT_EQ(s0.color.r, 307);
  HS_EXPECT_EQ(s0.color.g, 180);
  HS_EXPECT_EQ(s0.color.b, 1906);
  HS_EXPECT_EQ(s1.color.r, 36642);
  HS_EXPECT_EQ(s1.color.g, 9703);
  HS_EXPECT_EQ(s1.color.b, 157);
  // The baked slot reproduces the embers source it was baked from.
  Color4 e0 = Palettes::embers.get(0.0f);
  HS_EXPECT_EQ(s0.color.r, e0.color.r);
  HS_EXPECT_EQ(s0.color.g, e0.color.g);
  HS_EXPECT_EQ(s0.color.b, e0.color.b);

  // Every slot bakes a distinct LUT (no two share a t=0 color).
  for (int i = 0; i < MeshPaletteBank::N; ++i)
    for (int j = i + 1; j < MeshPaletteBank::N; ++j) {
      Color4 ci = bank[i].get(0.0f), cj = bank[j].get(0.0f);
      HS_EXPECT_TRUE(ci.color.r != cj.color.r || ci.color.g != cj.color.g ||
                     ci.color.b != cj.color.b);
    }
}

/**
 * @brief Verifies shuffle_indices yields a permutation of [0, N).
 * @details Each of 0..N-1 must appear exactly once: a shuffle that dropped or
 *          duplicated a slot would leave a palette unassigned or doubled.
 */
inline void test_mesh_palette_bank_shuffle_is_permutation() {
  std::array<int, MeshPaletteBank::N> idx{};
  MeshPaletteBank::shuffle_indices(idx);
  std::array<int, MeshPaletteBank::N> seen{};
  for (int v : idx) {
    HS_EXPECT_GE(v, 0);
    HS_EXPECT_LT(v, MeshPaletteBank::N);
    seen[v]++;
  }
  for (int count : seen)
    HS_EXPECT_EQ(count, 1);
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every palettes-module test and reports the aggregate result.
 * @return 0 on success, non-zero on any failure.
 */
inline int run_palettes_tests() {
  hs_test::ModuleFixture fixture("palettes");

  test_named_procedural_palette_endpoints();
  test_named_palette_hue_short_arc();
  test_mesh_palette_bank_lookup();
  test_mesh_palette_bank_shuffle_is_permutation();

  return fixture.result();
}

} // namespace palettes_tests
} // namespace hs_test
