/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/color.h.
 *
 * Coverage:
 *   - Pixel16::lerp16 (endpoints + midpoint)
 *   - Blend modes (over/under/max/mean/add) — identity & boundedness invariants
 *   - OKLab / OKLCH round-trips (sRGB -> OKLab -> sRGB, sRGB -> OKLCH -> sRGB)
 *     for grays and saturated primaries; achromatic hue handling in lerp_oklch
 *   - sRGB <-> linear LUTs vs. float reference functions; known endpoints
 *   - Gradient::get and BakedPalette::get: in-range validity + endpoint match
 */
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

#include "core/color.h"
#include "core/util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace color_tests {

// ============================================================================
// lerp16
// ============================================================================

/**
 * @brief Verifies frac=0 and frac=65535 recover the two endpoints exactly.
 */
inline void test_lerp16_endpoints() {
  Pixel16 a(1000, 2000, 3000);
  Pixel16 b(40000, 50000, 60000);

  Pixel16 at0 = a.lerp16(b, 0);
  HS_EXPECT_EQ(at0.r, a.r);
  HS_EXPECT_EQ(at0.g, a.g);
  HS_EXPECT_EQ(at0.b, a.b);

  Pixel16 at1 = a.lerp16(b, 65535);
  HS_EXPECT_EQ(at1.r, b.r);
  HS_EXPECT_EQ(at1.g, b.g);
  HS_EXPECT_EQ(at1.b, b.b);
}

/**
 * @brief Verifies frac~0.5 lands each channel on (a+b)/2 within rounding error.
 * @details Each channel lands on (a+b)/2 within the /65535 rounding error;
 *          equal endpoints stay put.
 */
inline void test_lerp16_midpoint() {
  Pixel16 a(0, 100, 65535);
  Pixel16 b(65535, 300, 65535);
  Pixel16 mid = a.lerp16(b, 32768); // ~0.5

  HS_EXPECT_NEAR(static_cast<float>(mid.r), 32767.0f, 2.0f);
  HS_EXPECT_NEAR(static_cast<float>(mid.g), 200.0f, 2.0f);
  HS_EXPECT_EQ(mid.b, 65535); // both endpoints equal -> stays put
}

/**
 * @brief Verifies lerp16 rounds to nearest, not floor.
 * @details The reconstruction tail (x + (x>>16) + 32768) >> 16 adds half a
 *          quantum so the divide rounds; at frac = 49152 (~0.75) round-to-nearest
 *          and floor disagree on every channel, pinning the rounding behavior for
 *          both the portable and smlad paths.
 */
inline void test_lerp16_rounds_to_nearest() {
  Pixel16 a(0, 0, 0);
  Pixel16 b(1, 2, 4);
  Pixel16 m = a.lerp16(b, 49152); // 0.75
  // True values 0.75 / 1.5 / 3.0 -> round-to-nearest 1 / 2 / 3 (floor: 0 / 1 / 2).
  HS_EXPECT_EQ(m.r, 1);
  HS_EXPECT_EQ(m.g, 2);
  HS_EXPECT_EQ(m.b, 3);
}

/**
 * @brief Verifies every interpolated channel lies within the endpoint envelope.
 * @details Each channel must lie within the [min,max] envelope of the two
 *          endpoints (allowing +/- 1 LSB of rounding slack).
 */
inline void test_lerp16_bounded() {
  Pixel16 a(123, 45678, 65535);
  Pixel16 b(65535, 0, 12345);
  for (uint32_t f = 0; f <= 65535; f += 4095) {
    Pixel16 m = a.lerp16(b, static_cast<uint16_t>(f));
    HS_EXPECT_LE(m.r, static_cast<uint16_t>(std::max(a.r, b.r)));
    HS_EXPECT_GE(m.r + 1, static_cast<uint16_t>(std::min(a.r, b.r)));
    HS_EXPECT_LE(m.g, static_cast<uint16_t>(std::max(a.g, b.g)));
    HS_EXPECT_GE(m.g + 1, static_cast<uint16_t>(std::min(a.g, b.g)));
    HS_EXPECT_LE(m.b, static_cast<uint16_t>(std::max(a.b, b.b)));
    HS_EXPECT_GE(m.b + 1, static_cast<uint16_t>(std::min(a.b, b.b)));
  }
}

/**
 * @brief Round-to-nearest div-by-65535 lerp reference computed in double.
 * @param a First endpoint channel value in [0, 65535].
 * @param b Second endpoint channel value in [0, 65535].
 * @param frac Interpolation fraction in [0, 65535] mapping to [0, 1].
 * @return Interpolated channel value in [0, 65535], rounded to nearest.
 */
inline uint16_t lerp16_reference(uint16_t a, uint16_t b, uint16_t frac) {
  double t = static_cast<double>(frac) / 65535.0;
  return static_cast<uint16_t>(a * (1.0 - t) + b * t + 0.5);
}

/**
 * @brief Verifies lerp16 is correct across the full 0..65535 operand range.
 * @details A signed 16x16 multiply would misread any operand >= 32768 (a frac,
 *          an inverse-frac, or a bright channel) as negative and corrupt the
 *          whole upper half by up to a full 65535; this pins unsigned-correct
 *          results so the device's MAC path stays bit-exact with the double
 *          reference. The native build can't run the ARM asm, so this guards the
 *          behavior, not the instruction.
 */
inline void test_lerp16_full_range_correct() {
  // Midpoint of two maximal channels stays maximal — the canonical case a signed
  // multiply collapses (65535 read as -1 -> product ~0).
  Pixel16 hi(65535, 65535, 65535), lo(0, 0, 0);
  Pixel16 mid = hi.lerp16(lo, 32768);
  HS_EXPECT_NEAR(static_cast<float>(mid.r), 32768.0f, 2.0f);
  HS_EXPECT_NEAR(static_cast<float>(mid.g), 32768.0f, 2.0f);

  // Bright endpoints recovered exactly.
  Pixel16 a(65535, 49152, 40000), b(32768, 60000, 33000);
  HS_EXPECT_TRUE(a.lerp16(b, 0) == a);
  HS_EXPECT_TRUE(a.lerp16(b, 65535) == b);

  // Sweep high-operand pairs (all >= 32768) against the double reference — the
  // regime a signed multiply would corrupt.
  const uint16_t vals[] = {32768, 40000, 49152, 60000, 65535};
  for (uint16_t av : vals)
    for (uint16_t bv : vals)
      for (uint32_t f = 0; f <= 65535; f += 8191) {
        Pixel16 pa(av, 0, 0), pb(bv, 0, 0);
        uint16_t got = pa.lerp16(pb, static_cast<uint16_t>(f)).r;
        uint16_t ref = lerp16_reference(av, bv, static_cast<uint16_t>(f));
        HS_EXPECT_TRUE(std::abs(static_cast<int>(got) - static_cast<int>(ref)) <= 1);
      }
}

// ============================================================================
// Blend modes
// ============================================================================

/**
 * @brief Verifies blend_over picks the source and blend_under the destination.
 */
inline void test_blend_over_under() {
  Pixel16 dst(100, 200, 300);
  Pixel16 src(4000, 5000, 6000);
  HS_EXPECT_TRUE(blend_over(dst, src) == src);
  HS_EXPECT_TRUE(blend_under(dst, src) == dst);
}

/**
 * @brief Verifies blend_max takes the per-channel maximum.
 */
inline void test_blend_max() {
  Pixel16 a(100, 5000, 300);
  Pixel16 b(4000, 200, 6000);
  Pixel16 m = blend_max(a, b);
  HS_EXPECT_EQ(m.r, 4000);
  HS_EXPECT_EQ(m.g, 5000);
  HS_EXPECT_EQ(m.b, 6000);
}

/**
 * @brief Verifies blend_mean averages the two pixels per channel.
 * @details blend_mean rounds to nearest ((a+b+1)/2), matching this file's
 *          +0.5f / +32768 round-to-nearest parity, not floor. The even-sum
 *          cases below can't tell the two modes apart, so an explicit odd-sum
 *          case pins the rounding contract (mirroring test_lerp16_rounds_to_nearest).
 */
inline void test_blend_mean() {
  Pixel16 a(0, 100, 65534);
  Pixel16 b(65534, 300, 0);
  Pixel16 m = blend_mean(a, b);
  // Even sums: floor and round-to-nearest agree.
  HS_EXPECT_EQ(m.r, 32767);
  HS_EXPECT_EQ(m.g, 200);
  HS_EXPECT_EQ(m.b, 32767);

  // Odd sums: blend_mean rounds up to (a+b+1)/2, not (a+b)/2.
  Pixel16 odd = blend_mean(Pixel16(1, 3, 5), Pixel16(2, 4, 6));
  HS_EXPECT_EQ(odd.r, 2); // (1+2+1)/2 = 2  (floor would give 1)
  HS_EXPECT_EQ(odd.g, 4); // (3+4+1)/2 = 4  (floor would give 3)
  HS_EXPECT_EQ(odd.b, 6); // (5+6+1)/2 = 6  (floor would give 5)
}

/**
 * @brief Verifies blend_add with black is the identity in either order.
 */
inline void test_blend_add_identity_with_black() {
  Pixel16 c(1234, 5678, 9012);
  Pixel16 black(0, 0, 0);
  Pixel16 r = blend_add(c, black);
  HS_EXPECT_TRUE(r == c);
  Pixel16 r2 = blend_add(black, c);
  HS_EXPECT_TRUE(r2 == c);
}

/**
 * @brief Verifies blend_add sums per channel and clamps at 65535, not wrapping.
 */
inline void test_blend_add_saturates() {
  Pixel16 a(60000, 60000, 100);
  Pixel16 b(60000, 1000, 50);
  Pixel16 r = blend_add(a, b);
  HS_EXPECT_EQ(r.r, 65535);       // saturated
  HS_EXPECT_EQ(r.g, 61000);       // no saturation
  HS_EXPECT_EQ(r.b, 150);
  HS_EXPECT_LE(r.r, 65535);
}

/**
 * @brief Pins the device's packed (uqadd16) saturating-add lane layout.
 * @details The device path of operator+= / blend_add packs g|b into one 32-bit
 *          uqadd16 lane and r into another, then unpacks. That asm path never
 *          runs on the host, so a transposed g/b lane or a wrong unpack shift
 *          would ship silently (wrong colors, not a crash). pixel16_blend_add_packed
 *          shares the exact lane layout with the device and compiles natively via
 *          the software uqadd16, so this checks it against an independent
 *          per-channel saturating reference across cases that stress each lane.
 */
inline void test_blend_add_packed_lane_layout() {
  auto ref = [](uint32_t x, uint32_t y) -> uint16_t {
    uint32_t s = x + y;
    return (uint16_t)(s > 65535 ? 65535 : s);
  };
  const Pixel16 cases[][2] = {
      {Pixel16(60000, 1000, 40000), Pixel16(10000, 200, 40000)}, // r,b sat; g not
      {Pixel16(0, 65535, 0), Pixel16(65535, 0, 65535)},          // each lane to max
      {Pixel16(123, 45678, 9000), Pixel16(40000, 30000, 50)},    // g sat only
      {Pixel16(0, 0, 0), Pixel16(0, 0, 0)},                      // zero
  };
  for (const auto &c : cases) {
    Pixel16 got = pixel16_blend_add_packed(c[0], c[1]);
    HS_EXPECT_EQ(got.r, ref(c[0].r, c[1].r));
    HS_EXPECT_EQ(got.g, ref(c[0].g, c[1].g));
    HS_EXPECT_EQ(got.b, ref(c[0].b, c[1].b));
    // Host add operators must agree with the packed device layout.
    Pixel16 acc = c[0];
    acc += c[1];
    HS_EXPECT_TRUE(acc == got);
    HS_EXPECT_TRUE(blend_add(c[0], c[1]) == got);
  }
}

/**
 * @brief Verifies Color4::operator*= scales both pixel and alpha.
 * @details The SSAA averaging algebra (scale each sample by 1/N, then sum)
 *          relies on *= touching alpha alongside color; scaling color but not
 *          alpha would silently corrupt the running average's weight.
 */
inline void test_color4_scale_affects_color_and_alpha() {
  Color4 c(Pixel16(1000, 2000, 4000), 0.8f);
  c *= 0.5f;
  HS_EXPECT_EQ(c.color.r, 500);
  HS_EXPECT_EQ(c.color.g, 1000);
  HS_EXPECT_EQ(c.color.b, 2000);
  HS_EXPECT_NEAR(c.alpha, 0.4f, 1e-6f);
}

/**
 * @brief Verifies Color4::operator+= sums color and clamps alpha at 1.0.
 * @details Alpha must saturate at 1.0 so the sum stays a valid blend weight; a
 *          missing clamp would let summed coverage ride past 1. Also exercises
 *          the SSAA pattern: N samples each pre-scaled by 1/N sum back to the
 *          average with alpha <= 1.
 */
inline void test_color4_add_clamps_alpha_and_sums_color() {
  // Two near-opaque samples sum past 1.0 -> clamped.
  Color4 a(Pixel16(10000, 0, 0), 0.7f);
  a += Color4(Pixel16(20000, 100, 0), 0.6f);
  HS_EXPECT_EQ(a.color.r, 30000);
  HS_EXPECT_EQ(a.color.g, 100);
  HS_EXPECT_NEAR(a.alpha, 1.0f, 1e-6f); // 1.3 clamped

  // SSAA averaging: sum of (sample * 1/N) reproduces the average, alpha <= 1.
  const int N = 4;
  Color4 samples[N] = {
      Color4(Pixel16(40000, 0, 0), 1.0f),
      Color4(Pixel16(0, 40000, 8000), 1.0f),
      Color4(Pixel16(0, 8000, 40000), 0.5f),
      Color4(Pixel16(0, 8000, 24000), 0.5f),
  };
  Color4 acc(Pixel16(0, 0, 0), 0.0f);
  for (int i = 0; i < N; ++i) {
    Color4 s = samples[i];
    s *= 1.0f / N;
    acc += s;
  }
  HS_EXPECT_EQ(acc.color.r, 10000); // 40000/4
  HS_EXPECT_EQ(acc.color.g, 14000); // 56000/4
  HS_EXPECT_EQ(acc.color.b, 18000); // 72000/4
  HS_EXPECT_NEAR(acc.alpha, 0.75f, 1e-6f); // (1+1+0.5+0.5)/4
  HS_EXPECT_LE(acc.alpha, 1.0f);
}

/**
 * @brief Verifies blend_max against black is the identity.
 */
inline void test_blend_max_with_black_identity() {
  Pixel16 c(111, 222, 333);
  Pixel16 black(0, 0, 0);
  HS_EXPECT_TRUE(blend_max(c, black) == c);
}

/**
 * @brief Verifies blend_alpha clamps its alpha to [0,1] before the float->int cast.
 * @details In-range values interpolate with round-to-nearest (+0.5f) weight
 *          quantization and out-of-range/overflowing/NaN alphas saturate to an
 *          endpoint instead of invoking cast UB.
 */
inline void test_blend_alpha_clamps_before_cast() {
  Pixel16 a(0, 0, 0);
  Pixel16 b(60000, 40000, 20000);

  HS_EXPECT_TRUE(blend_alpha(0.0f)(a, b) == a); // fully a
  HS_EXPECT_TRUE(blend_alpha(1.0f)(a, b) == b); // fully b

  // Alpha rounds to nearest (+0.5f): 0.5 -> weight 32768, not 32767.
  HS_EXPECT_TRUE(blend_alpha(0.5f)(a, b) == a.lerp16(b, 32768));

  // Out-of-range alpha saturates: a >= 1 -> full b; a <= 0 -> full a.
  HS_EXPECT_TRUE(blend_alpha(1000.0f)(a, b) == b);
  HS_EXPECT_TRUE(blend_alpha(-5.0f)(a, b) == a);
  // Large enough to overflow int in an unclamped (int)(a*65535).
  HS_EXPECT_TRUE(blend_alpha(1e9f)(a, b) == b);
  // NaN folds to the hi bound via hs::clamp.
  Pixel16 nan_res = blend_alpha(NAN)(a, b);
  HS_EXPECT_TRUE(nan_res == b);
}

/**
 * @brief Verifies Pixel16 * float clamps each scaled channel into [0,65535]
 *        before the cast.
 * @details Overflowing scales saturate, negatives clamp to 0, and NaN folds to
 *          the hi bound (matching blend_alpha) rather than invoking cast UB.
 */
inline void test_pixel16_scale_clamps_before_cast() {
  Pixel16 c(100, 2000, 30000);

  HS_EXPECT_TRUE(c * 0.0f == Pixel16(0, 0, 0));
  HS_EXPECT_TRUE(c * 1.0f == c);
  HS_EXPECT_TRUE(c * 2.0f == Pixel16(200, 4000, 60000));

  // Half-LSB results round up: odd channels * 0.5 land on .5 and carry up.
  HS_EXPECT_TRUE(Pixel16(3, 5, 7) * 0.5f == Pixel16(2, 3, 4));

  // Overflowing scale saturates at 65535.
  HS_EXPECT_TRUE(c * 1e9f == Pixel16(65535, 65535, 65535));
  // Negative scale clamps to zero.
  HS_EXPECT_TRUE(c * -3.0f == Pixel16(0, 0, 0));
  // NaN folds to the hi bound.
  HS_EXPECT_TRUE(c * NAN == Pixel16(65535, 65535, 65535));
}

// ============================================================================
// OKLab / OKLCH round-trips
// ============================================================================

/**
 * @brief Round-trips one sRGB color through OKLab and asserts recovery.
 * @param r Source red channel in [0, 255].
 * @param g Source green channel in [0, 255].
 * @param b Source blue channel in [0, 255].
 * @param tol255 Allowed absolute error in 8-bit sRGB levels.
 * @details Path: sRGB[0-255] -> linear float -> OKLab -> linear float ->
 *          sRGB[0-255].
 */
inline void roundtrip_oklab(uint8_t r, uint8_t g, uint8_t b, float tol255) {
  float rf = srgb_to_linear_float(r / 255.0f);
  float gf = srgb_to_linear_float(g / 255.0f);
  float bf = srgb_to_linear_float(b / 255.0f);

  OKLab lab = linear_rgb_to_oklab(rf, gf, bf);
  float r2, g2, b2;
  oklab_to_linear_rgb(lab, r2, g2, b2);

  float rs = linear_to_srgb_float(hs::clamp(r2, 0.0f, 1.0f)) * 255.0f;
  float gs = linear_to_srgb_float(hs::clamp(g2, 0.0f, 1.0f)) * 255.0f;
  float bs = linear_to_srgb_float(hs::clamp(b2, 0.0f, 1.0f)) * 255.0f;

  HS_EXPECT_NEAR(rs, static_cast<float>(r), tol255);
  HS_EXPECT_NEAR(gs, static_cast<float>(g), tol255);
  HS_EXPECT_NEAR(bs, static_cast<float>(b), tol255);
}

/**
 * @brief Verifies sRGB -> OKLab -> sRGB recovers a spread of colors.
 * @details Recovers grays, primaries/secondaries, and an arbitrary color to
 *          within tol255 of the original 8-bit value.
 */
inline void test_oklab_roundtrip() {
  // Grays
  roundtrip_oklab(0, 0, 0, 1.0f);
  roundtrip_oklab(128, 128, 128, 1.0f);
  roundtrip_oklab(255, 255, 255, 1.0f);
  // Saturated primaries / secondaries
  roundtrip_oklab(255, 0, 0, 1.0f);
  roundtrip_oklab(0, 255, 0, 1.0f);
  roundtrip_oklab(0, 0, 255, 1.0f);
  roundtrip_oklab(255, 255, 0, 1.0f);
  roundtrip_oklab(0, 255, 255, 1.0f);
  roundtrip_oklab(255, 0, 255, 1.0f);
  // Arbitrary
  roundtrip_oklab(37, 142, 211, 1.0f);
}

/**
 * @brief Round-trips one sRGB color through OKLCH and asserts recovery.
 * @param r Source red channel in [0, 255].
 * @param g Source green channel in [0, 255].
 * @param b Source blue channel in [0, 255].
 * @param tol255 Allowed absolute error in 8-bit sRGB levels.
 * @details Path: sRGB[0-255] -> OKLCH -> OKLab -> linear -> sRGB[0-255].
 */
inline void roundtrip_oklch(uint8_t r, uint8_t g, uint8_t b, float tol255) {
  OKLCH lch = srgb_to_oklch(r, g, b);
  float r2, g2, b2;
  oklab_to_linear_rgb(oklch_to_oklab(lch), r2, g2, b2);
  float rs = linear_to_srgb_float(hs::clamp(r2, 0.0f, 1.0f)) * 255.0f;
  float gs = linear_to_srgb_float(hs::clamp(g2, 0.0f, 1.0f)) * 255.0f;
  float bs = linear_to_srgb_float(hs::clamp(b2, 0.0f, 1.0f)) * 255.0f;
  HS_EXPECT_NEAR(rs, static_cast<float>(r), tol255);
  HS_EXPECT_NEAR(gs, static_cast<float>(g), tol255);
  HS_EXPECT_NEAR(bs, static_cast<float>(b), tol255);
}

/**
 * @brief Verifies sRGB -> OKLCH -> sRGB recovers grays, primaries, and a sample.
 */
inline void test_oklch_roundtrip() {
  roundtrip_oklch(0, 0, 0, 1.0f);
  roundtrip_oklch(128, 128, 128, 1.0f);
  roundtrip_oklch(255, 255, 255, 1.0f);
  roundtrip_oklch(255, 0, 0, 1.0f);
  roundtrip_oklch(0, 255, 0, 1.0f);
  roundtrip_oklch(0, 0, 255, 1.0f);
  roundtrip_oklch(64, 180, 75, 1.0f);
}

/**
 * @brief Pins sRGB -> OKLab/OKLCH against published reference coordinates.
 * @details Round-trip and palette-golden tests only check the forward and
 *          inverse transforms agree; a transposed or otherwise wrong-but-
 *          invertible M1/M2 pair would round-trip cleanly. These triples are the
 *          canonical Ottosson sRGB references (white, pure red/green/blue) in the
 *          engine's units: L in [0,1], a/b Cartesian, h in radians. The tolerance
 *          (4e-3 on L/a/b, ~0.3deg on hue) catches a swapped matrix row/column
 *          while absorbing cbrtf rounding.
 */
inline void test_oklab_reference_triples() {
  struct Ref {
    uint8_t r, g, b;
    float L, a, bb; // expected OKLab
  };
  const Ref refs[] = {
      {255, 255, 255, 1.0000f, 0.0000f, 0.0000f},
      {255, 0, 0, 0.6279f, 0.2249f, 0.1258f},
      {0, 255, 0, 0.8664f, -0.2339f, 0.1795f},
      {0, 0, 255, 0.4520f, -0.0324f, -0.3115f},
  };
  const float tol = 4e-3f;
  for (const Ref &x : refs) {
    float rf = srgb_to_linear_float(x.r / 255.0f);
    float gf = srgb_to_linear_float(x.g / 255.0f);
    float bf = srgb_to_linear_float(x.b / 255.0f);
    OKLab lab = linear_rgb_to_oklab(rf, gf, bf);
    HS_EXPECT_NEAR(lab.L, x.L, tol);
    HS_EXPECT_NEAR(lab.a, x.a, tol);
    HS_EXPECT_NEAR(lab.b, x.bb, tol);
  }

  // Pure red in OKLCH: L=0.6279, C=0.2577, h=29.23deg.
  OKLCH red = srgb_to_oklch(255, 0, 0);
  HS_EXPECT_NEAR(red.L, 0.6279f, tol);
  HS_EXPECT_NEAR(red.C, 0.2577f, tol);
  HS_EXPECT_NEAR(red.h, 29.23f * PI_F / 180.0f, 5e-3f);
}

/**
 * @brief Verifies pure grays have ~zero chroma in OKLCH.
 */
inline void test_oklch_gray_is_achromatic() {
  OKLCH g0 = srgb_to_oklch(0, 0, 0);
  OKLCH g1 = srgb_to_oklch(128, 128, 128);
  OKLCH g2 = srgb_to_oklch(255, 255, 255);
  HS_EXPECT_NEAR(g0.C, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(g1.C, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(g2.C, 0.0f, 1e-3f);
}

/**
 * @brief Verifies lerp_oklch hue handling for achromatic endpoints.
 * @details Two grays force hue to 0, and a gray/chromatic pair adopts the
 *          chromatic side's hue.
 */
inline void test_lerp_oklch_achromatic_hue() {
  // Both achromatic -> hue forced to 0, endpoints recovered in L.
  OKLCH a = srgb_to_oklch(0, 0, 0);
  OKLCH b = srgb_to_oklch(255, 255, 255);
  OKLCH mid = lerp_oklch(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.h, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(mid.C, 0.0f, 1e-3f);
  HS_EXPECT_GT(mid.L, a.L);
  HS_EXPECT_LT(mid.L, b.L);

  // One achromatic endpoint -> hue taken from the chromatic side.
  OKLCH red = srgb_to_oklch(255, 0, 0);
  OKLCH from_gray = lerp_oklch(a, red, 0.5f);
  HS_EXPECT_NEAR(from_gray.h, red.h, 1e-6f);
}

/**
 * @brief Verifies lerp_oklch interpolates hue along the short arc across the seam.
 * @details atan2f puts h in [-PI, PI], so two hues straddling that seam are only
 *          a small arc apart even though their numeric difference is ~2*PI. A
 *          naive linear lerp would pass through h=0 (the opposite side of the
 *          wheel); the shortest-arc lerp must pass through the seam at +/-PI
 *          instead. This underpins the marquee color feature.
 */
inline void test_lerp_oklch_shortest_arc_midpoint() {
  const float L = 0.6f, C = 0.15f;

  // Straddle the +/-PI seam: 2.8 and -2.8 rad are ~0.68 rad apart the short way
  // (through +/-PI), ~5.6 rad apart the long way (through 0).
  OKLCH a{L, C, 2.8f};
  OKLCH b{L, C, -2.8f};
  OKLCH mid = lerp_oklch(a, b, 0.5f);
  // Short arc midpoint sits at the seam (+/-PI), not at 0.
  HS_EXPECT_NEAR(std::fabs(mid.h), PI_F, 1e-4f);
  HS_EXPECT_NEAR(mid.L, L, 1e-5f);
  HS_EXPECT_NEAR(mid.C, C, 1e-5f);

  // Same seam, opposite winding.
  OKLCH c{L, C, 3.0f};
  OKLCH d{L, C, -3.0f};
  OKLCH mid2 = lerp_oklch(c, d, 0.5f);
  HS_EXPECT_NEAR(std::fabs(mid2.h), PI_F, 1e-4f);

  // A non-seam-crossing pair interpolates directly.
  OKLCH e{L, C, 0.5f};
  OKLCH f{L, C, 1.5f};
  HS_EXPECT_NEAR(lerp_oklch(e, f, 0.5f).h, 1.0f, 1e-4f);
}

/**
 * @brief Verifies lerp_oklch at amount 0/1 reproduces the endpoints' L and C.
 */
inline void test_lerp_oklch_endpoints() {
  OKLCH a = srgb_to_oklch(200, 30, 30);
  OKLCH b = srgb_to_oklch(30, 30, 200);
  OKLCH at0 = lerp_oklch(a, b, 0.0f);
  OKLCH at1 = lerp_oklch(a, b, 1.0f);
  HS_EXPECT_NEAR(at0.L, a.L, 1e-5f);
  HS_EXPECT_NEAR(at0.C, a.C, 1e-5f);
  HS_EXPECT_NEAR(at1.L, b.L, 1e-5f);
  HS_EXPECT_NEAR(at1.C, b.C, 1e-5f);
}

/**
 * @brief Verifies extrapolating amounts still yield a valid OKLCH.
 * @details Extrapolating amounts (reachable via unbounded
 *          GenerativePalette::lerp / ColorWipe paths) must clamp L to [0,1] and
 *          keep C non-negative, so an overshoot can't flip the hue 180deg or
 *          render near-black.
 */
inline void test_lerp_oklch_extrapolation_clamped() {
  OKLCH dark{0.1f, 0.05f, 0.0f};
  OKLCH bright{0.9f, 0.20f, 1.0f};

  OKLCH under = lerp_oklch(dark, bright, -2.0f); // overshoots below dark
  HS_EXPECT_GE(under.L, 0.0f);
  HS_EXPECT_GE(under.C, 0.0f);

  OKLCH over = lerp_oklch(bright, dark, 3.0f); // overshoots past dark toward 0
  HS_EXPECT_GE(over.L, 0.0f);
  HS_EXPECT_GE(over.C, 0.0f);

  OKLCH high = lerp_oklch(dark, bright, 5.0f); // overshoots above bright
  HS_EXPECT_LE(high.L, 1.0f);
}

/**
 * @brief Verifies an out-of-gamut OKLCH clamps into [0,65535] per channel.
 */
inline void test_oklch_to_pixel_bounded() {
  OKLCH vivid{0.7f, 0.4f, 1.0f};
  Pixel p = oklch_to_pixel(vivid);
  HS_EXPECT_LE(p.r, 65535);
  HS_EXPECT_LE(p.g, 65535);
  HS_EXPECT_LE(p.b, 65535);
}

// ============================================================================
// Chroma-reduction gamut mapping
// ============================================================================

/**
 * @brief Wraps a hue difference into (-PI, PI] for circular comparison.
 * @param dh Raw hue difference in radians.
 * @return The equivalent difference in (-PI, PI].
 */
inline float wrap_hue_delta(float dh) {
  while (dh > PI_F) dh -= 2.0f * PI_F;
  while (dh < -PI_F) dh += 2.0f * PI_F;
  return dh;
}

/**
 * @brief Verifies the chroma-reduction map holds hue and lightness in-gamut.
 * @details A deeply out-of-gamut OKLCH (chroma far past the sRGB cusp at this L)
 *          must map back inside the cube by SHRINKING chroma — not by the
 *          per-channel RGB clip that twists hue on saturated colors. For a spread
 *          of hues: assert the direct conversion really is out of gamut (so the
 *          map is exercised), then assert the mapped color is in gamut, its L and
 *          hue match the input within tolerance, and its chroma is strictly
 *          smaller but still positive.
 */
inline void test_gamut_clip_preserves_hue() {
  const float L = 0.65f, C = 0.40f;
  for (float h : {0.3f, 1.2f, 2.5f, -1.0f, -2.6f}) {
    OKLab lab = oklch_to_oklab({L, C, h});

    // Precondition: this chroma is unreachable in sRGB at this L.
    float r0, g0, b0;
    oklab_to_linear_rgb(lab, r0, g0, b0);
    HS_EXPECT_FALSE(linear_rgb_in_gamut(r0, g0, b0));

    OKLab mapped = gamut_clip_preserve_chroma(lab);
    float r, g, b;
    oklab_to_linear_rgb(mapped, r, g, b);
    HS_EXPECT_TRUE(linear_rgb_in_gamut(r, g, b));

    OKLCH out = oklab_to_oklch(mapped);
    HS_EXPECT_NEAR(out.L, L, 1e-4f);
    HS_EXPECT_NEAR(wrap_hue_delta(out.h - h), 0.0f, 1e-3f);
    HS_EXPECT_LT(out.C, C);   // chroma reduced...
    HS_EXPECT_GT(out.C, 0.0f); // ...but not crushed
  }
}

/**
 * @brief Verifies oklch_to_pixel routes out-of-gamut colors through the
 *        chroma-reduction map, holding hue where a per-channel clip would not.
 * @details Realizes a past-cusp OKLCH as a Pixel, reads the realized color back
 *          through the exact forward transform, and checks the hue survived the
 *          16-bit quantization. A per-channel clip on this color would swing the
 *          hue well past the tolerance toward the nearest primary.
 */
inline void test_oklch_to_pixel_holds_hue_out_of_gamut() {
  const float L = 0.62f, C = 0.42f, h = 0.9f;
  Pixel p = oklch_to_pixel({L, C, h});

  float r = p.r / 65535.0f, g = p.g / 65535.0f, b = p.b / 65535.0f;
  OKLCH got = oklab_to_oklch(linear_rgb_to_oklab(r, g, b));
  HS_EXPECT_NEAR(wrap_hue_delta(got.h - h), 0.0f, 2e-2f);
}

// ============================================================================
// fast_cbrt + perceptual hue_rotate (OKLab)
// ============================================================================

/**
 * @brief Verifies fast_cbrt accuracy against cbrtf and its domain guard.
 * @details Matches cbrtf to ~1e-4 relative error over the linear-RGB range, is
 *          exact at 0, and the negative/zero-domain guard returns 0.
 */
inline void test_fast_cbrt_accuracy() {
  HS_EXPECT_EQ(fast_cbrt(0.0f), 0.0f);
  HS_EXPECT_NEAR(fast_cbrt(1.0f), 1.0f, 1e-4f);
  HS_EXPECT_NEAR(fast_cbrt(8.0f), 2.0f, 1e-3f);
  // Negative / zero domain guard returns 0.
  HS_EXPECT_EQ(fast_cbrt(-1.0f), 0.0f);
  for (int k = 1; k <= 800; ++k) {
    float x = 8.0f * static_cast<float>(k) / 800.0f;
    float approx = fast_cbrt(x);
    float exact = cbrtf(x);
    HS_EXPECT_TRUE(std::fabs(approx - exact) / exact < 1e-4f);
  }
}

/**
 * @brief Verifies a perceptual hue rotation leaves a gray unchanged.
 * @details A gray (achromatic) color has zero chroma in OKLab, so rotating the
 *          (a,b) plane is a no-op and the color must come back unchanged for any
 *          amount, with alpha preserved. The OKLab matrices round-trip a gray
 *          channel exactly here (measured delta 0 of 65535), so the tolerance is
 *          1 LSB of the 16-bit linear channel — precision-grade, not the former
 *          48-LSB smoke band — leaving room only for a stray rounding ULP.
 */
inline void test_hue_rotate_preserves_gray() {
  Color4 gray(128, 128, 128, 0.5f);
  for (float amt : {0.1f, 0.25f, 0.5f, 0.8f}) {
    Color4 out = hue_rotate(gray, amt);
    HS_EXPECT_NEAR(static_cast<float>(out.color.r),
                   static_cast<float>(gray.color.r), 1.0f);
    HS_EXPECT_NEAR(static_cast<float>(out.color.g),
                   static_cast<float>(gray.color.g), 1.0f);
    HS_EXPECT_NEAR(static_cast<float>(out.color.b),
                   static_cast<float>(gray.color.b), 1.0f);
    HS_EXPECT_NEAR(out.alpha, 0.5f, 1e-5f);
  }
}

/**
 * @brief Verifies a full-turn rotation returns to the original color.
 * @details A full-turn rotation (amount = 1.0) should be the identity. The
 *          residual is the combined error of fast_cosf/fast_sinf at 2*PI (the
 *          turn never lands on an exact cos=1, sin=0) and the fast_cbrt OKLab
 *          round-trip: measured at most 4 LSB of the 16-bit linear channel on
 *          this saturated sample. The tolerance is 12 LSB — a ~3x margin over
 *          that measured error, so it stays precision-grade (vs the former
 *          64-LSB smoke band) while tolerating minor float-rounding drift.
 */
inline void test_hue_rotate_full_turn_identity() {
  Color4 c(200, 60, 30, 1.0f);
  Color4 out = hue_rotate(c, 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(out.color.r),
                 static_cast<float>(c.color.r), 12.0f);
  HS_EXPECT_NEAR(static_cast<float>(out.color.g),
                 static_cast<float>(c.color.g), 12.0f);
  HS_EXPECT_NEAR(static_cast<float>(out.color.b),
                 static_cast<float>(c.color.b), 12.0f);
}

// ============================================================================
// sRGB <-> linear LUTs vs. float reference
// ============================================================================

/**
 * @brief Verifies sRGB 0 maps to linear 0 and sRGB 255 to (near) max linear.
 */
inline void test_srgb_to_linear_endpoints() {
  HS_EXPECT_EQ(srgb_to_linear(0), 0);
  HS_EXPECT_EQ(srgb_to_linear(255), 65535);
}

/**
 * @brief Verifies the inverse LUT maps linear 0 to sRGB 0 and max to sRGB 255.
 */
inline void test_linear_to_srgb_endpoints() {
  HS_EXPECT_EQ(linear_to_srgb_lut[0], 0);
  HS_EXPECT_EQ(linear_to_srgb_lut[65535], 255);
}

/**
 * @brief Verifies the 8-bit -> 16-bit linear LUT matches the float reference.
 * @details Matches the float reference (scaled to 16-bit) across all 256 entries
 *          within rounding tolerance.
 */
inline void test_srgb_linear_lut_vs_float_reference() {
  for (int s = 0; s <= 255; ++s) {
    float ref = srgb_to_linear_float(s / 255.0f) * 65535.0f;
    float lut = static_cast<float>(srgb_to_linear(static_cast<uint8_t>(s)));
    HS_EXPECT_NEAR(lut, ref, 2.0f);
  }
}

/**
 * @brief Verifies sRGB -> linear (LUT) -> sRGB (LUT) recovers the 8-bit value.
 */
inline void test_srgb_linear_roundtrip_lut() {
  for (int s = 0; s <= 255; ++s) {
    uint16_t lin = srgb_to_linear(static_cast<uint8_t>(s));
    uint8_t back = linear_to_srgb_lut[lin];
    HS_EXPECT_EQ(static_cast<int>(back), s);
  }
}

/**
 * @brief Verifies srgb_to_linear_interp recovers sub-pixel precision.
 * @details Interpolates between LUT entries so sub-8-bit fractions resolve to
 *          distinct, ordered, monotonic 16-bit values rather than collapsing to
 *          one bucket. Endpoints and 1/255 steps still match the integer LUT.
 */
inline void test_srgb_to_linear_interp_recovers_subpixel_precision() {
  // Endpoints match the integer LUT.
  HS_EXPECT_EQ(srgb_to_linear_interp(0.0f), 0);
  HS_EXPECT_EQ(srgb_to_linear_interp(1.0f), 65535);

  // At exact 1/255 steps the interpolated value equals the integer LUT entry
  // (within float rounding).
  for (int s = 0; s <= 255; ++s) {
    int interp = srgb_to_linear_interp(s / 255.0f);
    int lut = srgb_to_linear(static_cast<uint8_t>(s));
    HS_EXPECT_TRUE(std::abs(interp - lut) <= 2);
  }

  // Two sRGB values in the SAME 8-bit bucket (200) but at different fractions
  // must map to DIFFERENT, ordered 16-bit linear values.
  uint16_t a = srgb_to_linear_interp(200.2f / 255.0f);
  uint16_t b = srgb_to_linear_interp(200.8f / 255.0f);
  HS_EXPECT_TRUE(b > a);
  // ...and both lie between the two LUT entries being interpolated.
  uint16_t e200 = srgb_to_linear(200), e201 = srgb_to_linear(201);
  HS_EXPECT_TRUE(a >= e200 && a <= e201);
  HS_EXPECT_TRUE(b >= e200 && b <= e201);

  // Monotonic non-decreasing across the full range.
  uint16_t prev = 0;
  for (int k = 0; k <= 1000; ++k) {
    uint16_t v = srgb_to_linear_interp(k / 1000.0f);
    HS_EXPECT_TRUE(v >= prev);
    prev = v;
  }
}

/**
 * @brief Verifies the float reference sRGB<->linear round-trip recovers values.
 * @details sRGB -> linear -> sRGB through the float reference functions recovers
 *          the original 8-bit value within half a level.
 */
inline void test_srgb_linear_roundtrip_float() {
  for (int s = 0; s <= 255; ++s) {
    float f = s / 255.0f;
    float lin = srgb_to_linear_float(f);
    float back = linear_to_srgb_float(lin);
    HS_EXPECT_NEAR(back * 255.0f, static_cast<float>(s), 0.5f);
  }
}

// ============================================================================
// Gradient::get
// ============================================================================

/**
 * @brief Verifies a black->white gradient yields black at t=0 and white at t~1.
 */
inline void test_gradient_endpoints() {
  // Black at t=0, white at t=1, linear sRGB ramp between.
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  Color4 c0 = grad.get(0.0f);
  HS_EXPECT_EQ(c0.color.r, 0);
  HS_EXPECT_EQ(c0.color.g, 0);
  HS_EXPECT_EQ(c0.color.b, 0);

  // t near 1.0 lands on the last (white) entry.
  Color4 c1 = grad.get(0.999f);
  HS_EXPECT_GT(c1.color.r, 60000);
  HS_EXPECT_EQ(c1.alpha, 1.0f);
}

/**
 * @brief Verifies the in-range ramp yields a non-decreasing red channel.
 * @details Walking the in-range ramp of a black->white gradient yields a
 *          non-decreasing red channel. Out-of-range clamping is covered by the
 *          *_clamps_* test.
 */
inline void test_gradient_in_range_valid_and_monotone() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};
  uint16_t prev = 0;
  for (int i = 0; i <= 100; ++i) {
    float t = i / 100.0f;
    if (t > 0.999f) t = 0.999f; // index = uint8_t(t*255); keep within [0,255]
    Color4 c = grad.get(t);
    HS_EXPECT_GE(c.color.r, prev);
    prev = c.color.r;
  }
}

/**
 * @brief Verifies a gradient with identical stops returns one color for any t.
 */
inline void test_gradient_solid_color() {
  Gradient grad{{0.0f, CPixel(10u, 20u, 30u)}, {1.0f, CPixel(10u, 20u, 30u)}};
  Color4 a = grad.get(0.0f);
  Color4 b = grad.get(0.5f);
  HS_EXPECT_EQ(a.color.r, b.color.r);
  HS_EXPECT_EQ(a.color.g, b.color.g);
  HS_EXPECT_EQ(a.color.b, b.color.b);
}

/**
 * @brief Verifies Gradient::get interpolates between LUT entries.
 * @details Two t values inside the same cell (both truncate to index 200) yield
 *          distinct colors, each bracketed by its neighbouring entries. A
 *          nearest-index lookup would band them together.
 */
inline void test_gradient_interpolates_between_entries() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};
  Color4 lo = grad.get(200.25f / 255.0f);
  Color4 hi = grad.get(200.75f / 255.0f);
  HS_EXPECT_GT(hi.color.r, lo.color.r);

  Color4 e_lo = grad.get(200.0f / 255.0f);
  Color4 e_hi = grad.get(201.0f / 255.0f);
  HS_EXPECT_GE(lo.color.r, e_lo.color.r);
  HS_EXPECT_LE(hi.color.r, e_hi.color.r);
}

/**
 * @brief Verifies Gradient::get clamps t to [0,1] for out-of-range input.
 * @details Out-of-range input saturates to an endpoint and never indexes past
 *          the 256-entry table (GenerativePalette::lerp can pass t slightly
 *          outside the unit interval). NaN folds to the hi bound.
 */
inline void test_gradient_get_clamps_out_of_range() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  // t < 0 saturates to the first entry (black).
  Color4 lo_end = grad.get(0.0f);
  Color4 below = grad.get(-0.5f);
  HS_EXPECT_EQ(below.color.r, lo_end.color.r);
  HS_EXPECT_EQ(below.color.g, lo_end.color.g);
  HS_EXPECT_EQ(below.color.b, lo_end.color.b);

  // t > 1 saturates to the last entry (white), same as the t==1 endpoint.
  Color4 hi_end = grad.get(1.0f);
  Color4 above = grad.get(1.5f);
  HS_EXPECT_EQ(above.color.r, hi_end.color.r);
  HS_EXPECT_EQ(above.color.g, hi_end.color.g);
  HS_EXPECT_EQ(above.color.b, hi_end.color.b);

  // NaN folds to the hi bound via hs::clamp -> last entry.
  Color4 nan_res = grad.get(NAN);
  HS_EXPECT_EQ(nan_res.color.r, hi_end.color.r);
}

/**
 * @brief Verifies a first stop at pos>0 flat-fills the LUT prefix with its color.
 * @details The constructor fills entries[0..firstStop] with the first stop's
 *          color, so a gradient whose first stop sits at 0.25 returns that color
 *          for all t in [0, 0.25]; only past the stop do the flanks interpolate.
 */
inline void test_gradient_first_stop_offset_flat_fills_prefix() {
  // First stop (pure red) at 0.25; second (pure blue) at 1.0.
  Gradient grad{{0.25f, CPixel(255u, 0u, 0u)}, {1.0f, CPixel(0u, 0u, 255u)}};

  Color4 at0 = grad.get(0.0f);
  Color4 flat = grad.get(0.1f); // still inside the [0,0.25] flat prefix
  // Pure red.
  HS_EXPECT_GT(at0.color.r, 60000);
  HS_EXPECT_EQ(at0.color.g, 0);
  HS_EXPECT_EQ(at0.color.b, 0);
  HS_EXPECT_EQ(flat.color.r, at0.color.r);
  HS_EXPECT_EQ(flat.color.b, at0.color.b);
  // Past the first stop the flank interpolates toward blue.
  Color4 ramp = grad.get(0.7f);
  HS_EXPECT_LT(ramp.color.r, at0.color.r);
  HS_EXPECT_GT(ramp.color.b, at0.color.b);
}

/**
 * @brief Verifies a >=3-stop gradient places the interior stop and interpolates flanks.
 * @details An interior stop not at 0/1 is the segment-join the 2-stop tests never
 *          exercise: the interior color must appear near its position and each
 *          flanking segment must blend between its bracketing stops.
 */
inline void test_gradient_three_stops_interior_and_flanks() {
  // red -> green (interior, 0.5) -> blue.
  Gradient grad{{0.0f, CPixel(255u, 0u, 0u)},
                {0.5f, CPixel(0u, 255u, 0u)},
                {1.0f, CPixel(0u, 0u, 255u)}};
  // Endpoints are the pure stops.
  Color4 a = grad.get(0.0f);
  Color4 b = grad.get(1.0f);
  HS_EXPECT_GT(a.color.r, 60000);
  HS_EXPECT_EQ(a.color.g, 0);
  HS_EXPECT_GT(b.color.b, 60000);
  HS_EXPECT_EQ(b.color.r, 0);
  // Interior green dominates at its stop.
  Color4 mid = grad.get(0.5f);
  HS_EXPECT_GT(mid.color.g, mid.color.r);
  HS_EXPECT_GT(mid.color.g, mid.color.b);
  // First flank (red->green): both red and green present, blue absent.
  Color4 f1 = grad.get(0.25f);
  HS_EXPECT_GT(f1.color.r, 0);
  HS_EXPECT_GT(f1.color.g, 0);
  HS_EXPECT_EQ(f1.color.b, 0);
  // Second flank (green->blue): green and blue present, red absent.
  Color4 f2 = grad.get(0.75f);
  HS_EXPECT_GT(f2.color.g, 0);
  HS_EXPECT_GT(f2.color.b, 0);
  HS_EXPECT_EQ(f2.color.r, 0);
}

/**
 * @brief Verifies two stops at the same quantized index produce a hard stop.
 * @details Coincident stop positions leave end==start, so the segment is skipped
 *          and the LUT jumps abruptly rather than interpolating. A smooth two-stop
 *          red->blue ramp would read as a red/blue mix on both sides of 0.5; the
 *          hard stop instead stays near-pure red below the boundary and near-pure
 *          blue above it.
 */
inline void test_gradient_hard_stop_is_abrupt() {
  Gradient grad{{0.0f, CPixel(255u, 0u, 0u)},
                {0.5f, CPixel(255u, 0u, 0u)},
                {0.5f, CPixel(0u, 0u, 255u)},
                {1.0f, CPixel(0u, 0u, 255u)}};
  // Below the boundary: essentially pure red.
  Color4 lo = grad.get(0.4f);
  HS_EXPECT_GT(lo.color.r, 60000);
  HS_EXPECT_LT(lo.color.b, 100);
  // Above the boundary: essentially pure blue. No red/blue blend in between.
  Color4 hi = grad.get(0.6f);
  HS_EXPECT_GT(hi.color.b, 60000);
  HS_EXPECT_LT(hi.color.r, 100);
}

// ============================================================================
// BakedPalette::get  (requires an Arena)
// ============================================================================

/**
 * @brief Verifies baking a solid-color source reproduces it at every sample.
 * @details Baking a solid-color source reproduces that color (and alpha) at
 *          every sample.
 */
inline void test_baked_palette_matches_source_endpoints() {
  // Source: solid color palette so every entry is identical.
  Color4 target(Pixel(1000, 2000, 3000), 1.0f);
  SolidColorPalette src(target);

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPalette baked;
  baked.bake(arena, src);

  Color4 c0 = baked.get(0.0f);
  Color4 c1 = baked.get(1.0f);
  HS_EXPECT_EQ(c0.color.r, target.color.r);
  HS_EXPECT_EQ(c0.color.g, target.color.g);
  HS_EXPECT_EQ(c0.color.b, target.color.b);
  HS_EXPECT_EQ(c1.color.r, target.color.r);
  HS_EXPECT_NEAR(c0.alpha, 1.0f, 1e-6f);
}

/**
 * @brief Verifies baking a black->white gradient preserves the envelope.
 * @details Preserves the endpoints and keeps interior samples within the
 *          [black,white] envelope at full alpha.
 */
inline void test_baked_palette_in_range() {
  // Ramp source via Gradient (black->white), bake, then sample.
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  alignas(std::max_align_t) static uint8_t
      buf[BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPalette baked;
  baked.bake(arena, grad);

  // Endpoints: t=0 -> first entry (black), t=1 -> last entry (white).
  Color4 c0 = baked.get(0.0f);
  Color4 c1 = baked.get(1.0f);
  HS_EXPECT_EQ(c0.color.r, 0);
  HS_EXPECT_GT(c1.color.r, 60000);

  // Interior samples stay within the [c0, c1] envelope and alpha == 1.
  for (int i = 0; i <= 50; ++i) {
    float t = i / 50.0f;
    Color4 c = baked.get(t);
    HS_EXPECT_GE(c.color.r, c0.color.r);
    HS_EXPECT_LE(c.color.r, c1.color.r);
    HS_EXPECT_NEAR(c.alpha, 1.0f, 1e-6f);
  }
}

// ============================================================================
// Palette layer: source palettes, modifiers, compositions, and wrappers
// ============================================================================

/**
 * @brief Verifies ProceduralPalette evaluates its cosine color formula.
 * @details C(t) = a + b*cos(2*PI*(c*t + d)) in sRGB, then to linear. A
 *          {a=.5,b=.5,c=1,d=0} channel is cos-driven: t=0 -> 1.0 (full),
 *          t=0.5 -> 0.0.
 */
inline void test_procedural_palette_cosine() {
  ProceduralPalette pp({0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f},
                       {1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f});
  Color4 c0 = pp.get(0.0f);
  Color4 chalf = pp.get(0.5f);
  HS_EXPECT_EQ(c0.color.r, 65535); // 0.5 + 0.5*cos(0) = 1.0
  HS_EXPECT_EQ(chalf.color.r, 0);  // 0.5 + 0.5*cos(PI) = 0.0
  HS_EXPECT_NEAR(c0.alpha, 1.0f, 1e-6f);
}

/**
 * @brief Verifies MutatingPalette blends between its two endpoint palettes.
 * @details mutate(0) reproduces palette #1, mutate(1) reproduces #2, and an
 *          interior amount lands between them.
 */
inline void test_mutating_palette_blends_endpoints() {
  // #1 -> white at t=0 ; #2 -> black everywhere.
  MutatingPalette m({0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, {1.0f, 1.0f, 1.0f},
                    {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f},
                    {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f});
  m.mutate(0.0f);
  HS_EXPECT_EQ(m.get(0.0f).color.r, 65535);
  m.mutate(1.0f);
  HS_EXPECT_EQ(m.get(0.0f).color.r, 0);
  m.mutate(0.5f); // a=.25, b=.25 -> 0.25 + 0.25*cos(0) = 0.5 sRGB
  uint16_t mid = m.get(0.0f).color.r;
  HS_EXPECT_GT(mid, 0);
  HS_EXPECT_LT(mid, 65535);
}

/**
 * @brief Verifies GenerativePalette is deterministic for a manual seed.
 * @details With a manual seed and rand-free profiles (FLAT/VIBRANT/TRIADIC use
 *          no RNG) the palette is fully deterministic.
 */
inline void test_generative_palette_deterministic() {
  auto make = [](int seed) {
    return GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                             BrightnessProfile::FLAT, SaturationProfile::VIBRANT,
                             seed);
  };
  GenerativePalette g1 = make(0), g2 = make(0), g3 = make(128);
  bool g3_differs = false;
  for (int i = 0; i <= 16; ++i) {
    float t = i / 16.0f;
    Color4 a = g1.get(t), b = g2.get(t), c = g3.get(t);
    // Same seed/params -> bit-identical.
    HS_EXPECT_EQ(a.color.r, b.color.r);
    HS_EXPECT_EQ(a.color.g, b.color.g);
    HS_EXPECT_EQ(a.color.b, b.color.b);
    if (c.color.r != a.color.r || c.color.g != a.color.g) g3_differs = true;
  }
  HS_EXPECT_TRUE(g3_differs); // a different seed must change the palette
}

/**
 * @brief Verifies GenerativePalette::get clamps t to [0,1].
 * @details t < 0 saturates to the first stop and t > 1 to the last, with no
 *          discontinuity at the ends. A NaN t folds to 1.0 (house clamp
 *          contract) -> last stop.
 */
inline void test_generative_palette_get_clamps_out_of_range() {
  GenerativePalette g(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                      BrightnessProfile::FLAT, SaturationProfile::VIBRANT, 0);
  Color4 first = g.get(0.0f);
  Color4 below = g.get(-0.5f);
  HS_EXPECT_EQ(below.color.r, first.color.r); // clamps to the first stop
  HS_EXPECT_EQ(below.color.g, first.color.g);
  HS_EXPECT_EQ(below.color.b, first.color.b);

  Color4 last = g.get(1.0f);
  Color4 nan = g.get(std::numeric_limits<float>::quiet_NaN());
  HS_EXPECT_EQ(nan.color.r, last.color.r); // NaN -> hi -> last stop
  HS_EXPECT_EQ(nan.color.g, last.color.g);
  HS_EXPECT_EQ(nan.color.b, last.color.b);
}

/**
 * @brief Verifies the Mobius longitude singularity saturates to a palette endpoint.
 * @details Pins the load-bearing NaN propagation in MobiusGrid::draw_longitudes
 *          (effects/MobiusGrid.h): at the stereographic singularity z = +/-1 the
 *          conformal radius R = sqrtf((1+z)/(1-z)) blows up to +inf, logf carries
 *          the inf through, and wrap() folds it to NaN, which palette.get clamps
 *          to a palette endpoint so the longitude saturates to its terminal color
 *          rather than misbehaving. The behavior is deliberate (the WASM release
 *          build keeps -fno-finite-math-only specifically to preserve it; see
 *          CMakeLists.txt). Locks the whole inf -> NaN -> endpoint chain so a
 *          future math nudge or a stray -ffinite-math-only that quietly breaks it
 *          fails here instead of on the sphere.
 */
inline void test_mobius_longitude_singularity_saturates_to_endpoint() {
  const float z = 1.0f; // sinf(0.25 * 2*PI) at the t_line = 0.25 singularity
  float R = std::sqrt((1.0f + z) / (1.0f - z));
  HS_EXPECT_TRUE(std::isinf(R)); // division by (1 - z) == 0 blows R up

  float log_r = std::log(R);
  const float log_min = -2.5f, log_max = 2.5f;
  float t = (log_r - log_min) / (log_max - log_min);
  float wrapped = wrap(t, 1.0f); // phase = 0; inf -> fmod -> NaN
  HS_EXPECT_TRUE(std::isnan(wrapped));

  // Same palette type/shape MobiusGrid uses (default GenerativePalette).
  GenerativePalette palette(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                            BrightnessProfile::FLAT, SaturationProfile::VIBRANT,
                            0);
  Color4 endpoint = palette.get(1.0f);
  Color4 singular = palette.get(wrapped);
  HS_EXPECT_EQ(singular.color.r, endpoint.color.r);
  HS_EXPECT_EQ(singular.color.g, endpoint.color.g);
  HS_EXPECT_EQ(singular.color.b, endpoint.color.b);
}

/**
 * @brief Verifies the auto-seed cursor advances and reset_hue_seed restores it.
 * @details The auto-seed cursor advances per construction, so two consecutive
 *          auto-seeded palettes differ; reset_hue_seed restores the cursor for
 *          reproducibility.
 */
inline void test_generative_palette_auto_seed_advances() {
  GenerativePalette::reset_hue_seed(0);
  GenerativePalette a(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                      BrightnessProfile::FLAT, SaturationProfile::VIBRANT);
  GenerativePalette b(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                      BrightnessProfile::FLAT, SaturationProfile::VIBRANT);
  HS_EXPECT_TRUE(a.get(0.0f).color.r != b.get(0.0f).color.r ||
                 a.get(0.0f).color.g != b.get(0.0f).color.g);

  // Resetting the cursor reproduces the first palette exactly.
  GenerativePalette::reset_hue_seed(0);
  GenerativePalette a2(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                       BrightnessProfile::FLAT, SaturationProfile::VIBRANT);
  HS_EXPECT_EQ(a.get(0.3f).color.r, a2.get(0.3f).color.r);
  HS_EXPECT_EQ(a.get(0.3f).color.g, a2.get(0.3f).color.g);
}

/**
 * @brief Verifies GenerativePalette::lerp cross-fades between two snapshots.
 * @details lerp(snapshot, snapshot, amount) drives a cross-fade: at amount 0/1
 *          it reproduces the endpoint palettes (within the OKLCH round-trip).
 */
inline void test_generative_palette_snapshot_lerp() {
  GenerativePalette from(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                         BrightnessProfile::FLAT, SaturationProfile::VIBRANT, 0);
  GenerativePalette to(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                       BrightnessProfile::FLAT, SaturationProfile::VIBRANT, 128);
  // Endpoints must be clearly distinct for the test to mean anything.
  HS_EXPECT_TRUE(std::abs(static_cast<int>(from.get(0.0f).color.r) -
                          static_cast<int>(to.get(0.0f).color.r)) > 4000 ||
                 std::abs(static_cast<int>(from.get(0.0f).color.g) -
                          static_cast<int>(to.get(0.0f).color.g)) > 4000);

  GenerativePalette mixed(GradientShape::STRAIGHT, HarmonyType::TRIADIC,
                          BrightnessProfile::FLAT, SaturationProfile::VIBRANT, 0);
  // amount 0 -> ~from, amount 1 -> ~to (OKLCH round-trip tolerance in linear).
  const float tol = 2500.0f;
  mixed.lerp(from.snapshot(), to.snapshot(), 0.0f);
  HS_EXPECT_NEAR(static_cast<float>(mixed.get(0.0f).color.r),
                 static_cast<float>(from.get(0.0f).color.r), tol);
  mixed.lerp(from.snapshot(), to.snapshot(), 1.0f);
  HS_EXPECT_NEAR(static_cast<float>(mixed.get(0.0f).color.r),
                 static_cast<float>(to.get(0.0f).color.r), tol);
}

/**
 * @brief Verifies each coordinate modifier's modify() in isolation.
 * @details Modifiers transform the palette coordinate; this exercises Scale,
 *          Cycle, Quantize, Fold, Breathe, Ripple, and Pinch one at a time.
 */
inline void test_palette_modifiers() {
  // Scale multiplies the coordinate.
  HS_EXPECT_NEAR(ScaleModifier(2.0f).modify(0.3f), 0.6f, 1e-5f);
  float dyn_scale = 3.0f;
  HS_EXPECT_NEAR(ScaleModifier(1.0f, &dyn_scale).modify(0.2f), 0.6f, 1e-5f);

  // Cycle adds the driver offset; a null driver is a deliberate pass-through.
  float off = 0.25f;
  HS_EXPECT_NEAR(CycleModifier(&off).modify(0.5f), 0.75f, 1e-5f);
  HS_EXPECT_NEAR(CycleModifier(nullptr).modify(0.5f), 0.5f, 1e-5f);

  // Quantize snaps to the nearest step: round(0.3*4)/4 = 1/4 = 0.25.
  HS_EXPECT_NEAR(QuantizeModifier(4.0f).modify(0.3f), 0.25f, 1e-5f);

  // Fold is a triangle wave: 0->1, 0.25->0.5, 0.5->0.
  FoldModifier fold(2.0f);
  HS_EXPECT_NEAR(fold.modify(0.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(fold.modify(0.25f), 0.5f, 1e-5f);
  HS_EXPECT_NEAR(fold.modify(0.5f), 0.0f, 1e-5f);

  // Breathe with a zero-phase driver is identity; a quarter-turn shifts by amp.
  float phase0 = 0.0f;
  HS_EXPECT_NEAR(BreatheModifier(&phase0, 0.1f).modify(0.5f), 0.5f, 1e-3f);

  // Ripple at t=0, phase=0 leaves the coordinate fixed (sin(0) = 0).
  float rphase = 0.0f;
  HS_EXPECT_NEAR(RippleModifier(&rphase, 3.0f, 0.1f).modify(0.0f), 0.0f, 1e-5f);

  // Pinch (positive tension) pulls an off-center coordinate toward 0.5.
  float tension = 0.5f;
  float pinched = PinchModifier(&tension).modify(0.25f);
  HS_EXPECT_GT(pinched, 0.25f);
  HS_EXPECT_LT(pinched, 0.5f);
  // A null tension driver passes through.
  HS_EXPECT_NEAR(PinchModifier(nullptr).modify(0.3f), 0.3f, 1e-5f);
}

/**
 * @brief Verifies StaticPalette folds its modifier chain then queries the source.
 * @details Applies the modifier chain in order before sampling the source.
 */
inline void test_static_palette_composition() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  // Single modifier: scale=2 turns get(0.25) into source.get(0.5).
  ScaleModifier scale(2.0f);
  StaticPalette<Gradient, Coords<ScaleModifier>> sp;
  sp.bind(&grad, &scale);
  HS_EXPECT_EQ(sp.get(0.25f).color.r, grad.get(0.5f).color.r);

  // Two modifiers apply in tuple order (scale THEN cycle): 0.2 -> 0.4 -> 0.5.
  float off = 0.1f;
  CycleModifier cycle(&off);
  StaticPalette<Gradient, Coords<ScaleModifier, CycleModifier>> sp2;
  sp2.bind(&grad, &scale, &cycle);
  HS_EXPECT_EQ(sp2.get(0.2f).color.r, grad.get(0.5f).color.r);
}

/**
 * @brief Constant 0.5 falloff function for AlphaFalloffShade.
 * @param Unused normalized coordinate (the falloff is constant).
 * @return Constant falloff factor 0.5.
 * @details Must be a plain function pointer to bind to AlphaFalloffShade.
 */
inline float test_half_falloff(float) { return 0.5f; }

/**
 * @brief Verifies each modifier composed via StaticPalette produces its remap.
 * @details Exercises Reverse, Mirror, Inset+EdgeFade, Inset+EdgeAlpha,
 *          AlphaFalloff, PaletteFacade, and SolidColor. Bounded remaps use
 *          Wrap=false so a coordinate landing exactly on 1.0 reaches the source's
 *          last stop rather than wrapping to 0.
 */
inline void test_palette_wrappers() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  // ReverseModifier: t -> 1-t, so t=0 samples the white end and t=1 the black.
  ReverseModifier rev;
  StaticPalette<Gradient, Coords<ReverseModifier>, Colors<>, /*Wrap=*/false> revp;
  revp.bind(&grad, &rev);
  HS_EXPECT_GT(revp.get(0.0f).color.r, 60000); // source's t=1 (white)
  HS_EXPECT_EQ(revp.get(1.0f).color.r, 0);     // source's t=0 (black)

  // MirrorModifier: [0,1] -> [0,1,0]; the midpoint reaches the far end.
  MirrorModifier mir;
  StaticPalette<Gradient, Coords<MirrorModifier>, Colors<>, /*Wrap=*/false> mirp;
  mirp.bind(&grad, &mir);
  HS_EXPECT_EQ(mirp.get(0.0f).color.r, 0);
  HS_EXPECT_GT(mirp.get(0.5f).color.r, 60000);
  HS_EXPECT_EQ(mirp.get(1.0f).color.r, 0);

  // Opaque vignette = InsetModifier + EdgeFadeShade: fades to black at the
  // edges, source color in the middle band.
  InsetModifier inset;
  EdgeFadeShade fade;
  StaticPalette<Gradient, Coords<InsetModifier>, Colors<EdgeFadeShade>,
                /*Wrap=*/false>
      vig;
  vig.bind(&grad, &inset, &fade);
  HS_EXPECT_LT(vig.get(0.0f).color.r, 1000); // edge -> ~black
  uint16_t vmid = vig.get(0.5f).color.r;     // middle -> source.get(0.5)
  HS_EXPECT_GT(vmid, 1000);
  HS_EXPECT_LT(vmid, 64000);

  // Transparent vignette = InsetModifier + EdgeAlphaShade: alpha (not color)
  // fades at the edges.
  EdgeAlphaShade alpha_fade;
  StaticPalette<Gradient, Coords<InsetModifier>, Colors<EdgeAlphaShade>,
                /*Wrap=*/false>
      tv;
  tv.bind(&grad, &inset, &alpha_fade);
  HS_EXPECT_NEAR(tv.get(0.0f).alpha, 0.0f, 1e-3f); // edge -> transparent
  HS_EXPECT_NEAR(tv.get(0.1f).alpha, 0.5f, 1e-2f); // quintic(0.5) = 0.5
  HS_EXPECT_NEAR(tv.get(0.5f).alpha, 1.0f, 1e-3f); // middle -> opaque

  // AlphaFalloffShade scales alpha by the falloff function.
  AlphaFalloffShade afs(test_half_falloff);
  StaticPalette<Gradient, Coords<>, Colors<AlphaFalloffShade>, /*Wrap=*/false>
      afp;
  afp.bind(&grad, &afs);
  HS_EXPECT_NEAR(afp.get(0.5f).alpha, 0.5f, 1e-5f);
  HS_EXPECT_EQ(afp.get(0.5f).color.r, grad.get(0.5f).color.r);

  // PaletteFacade exposes a composition through the polymorphic Palette API.
  PaletteFacade<decltype(afp)> facade(&afp);
  const Palette &as_palette = facade;
  HS_EXPECT_NEAR(as_palette.get(0.5f).alpha, 0.5f, 1e-5f);

  // Solid color ignores t.
  SolidColorPalette solid(Color4(Pixel(111, 222, 333), 0.7f));
  HS_EXPECT_EQ(solid.get(0.9f).color.g, 222);
  HS_EXPECT_NEAR(solid.get(0.1f).alpha, 0.7f, 1e-6f);
}

// Clamp-before-cast / NaN-saturation checks whose contract must also hold under
// the shipping WASM fast-math codegen. fastmath_clamp_check.cpp iterates this
// same list, so adding a case here automatically extends both the default-IEEE
// run and the -ffast-math -fno-finite-math-only pass.
#define HS_FASTMATH_CLAMP_TESTS(X)                                             \
  X(test_blend_alpha_clamps_before_cast)                                       \
  X(test_pixel16_scale_clamps_before_cast)                                     \
  X(test_gradient_get_clamps_out_of_range)                                     \
  X(test_generative_palette_get_clamps_out_of_range)                           \
  X(test_mobius_longitude_singularity_saturates_to_endpoint)

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every color-module test and reports the aggregate result.
 * @return Process exit code from hs_test::end_module: 0 on success, non-zero on
 *         any failure.
 */
inline int run_color_tests() {
  hs_test::ModuleFixture fixture("color");

  test_lerp16_endpoints();
  test_lerp16_midpoint();
  test_lerp16_rounds_to_nearest();
  test_lerp16_full_range_correct();
  test_lerp16_bounded();

  test_blend_over_under();
  test_blend_max();
  test_blend_mean();
  test_blend_add_identity_with_black();
  test_blend_add_saturates();
  test_blend_add_packed_lane_layout();
  test_color4_scale_affects_color_and_alpha();
  test_color4_add_clamps_alpha_and_sums_color();
  test_blend_max_with_black_identity();

  test_oklab_roundtrip();
  test_oklch_roundtrip();
  test_oklab_reference_triples();
  test_oklch_gray_is_achromatic();
  test_lerp_oklch_achromatic_hue();
  test_lerp_oklch_shortest_arc_midpoint();
  test_lerp_oklch_endpoints();
  test_lerp_oklch_extrapolation_clamped();
  test_oklch_to_pixel_bounded();
  test_gamut_clip_preserves_hue();
  test_oklch_to_pixel_holds_hue_out_of_gamut();

  test_fast_cbrt_accuracy();
  test_hue_rotate_preserves_gray();
  test_hue_rotate_full_turn_identity();

  test_srgb_to_linear_endpoints();
  test_linear_to_srgb_endpoints();
  test_srgb_linear_lut_vs_float_reference();
  test_srgb_linear_roundtrip_lut();
  test_srgb_to_linear_interp_recovers_subpixel_precision();
  test_srgb_linear_roundtrip_float();

  test_gradient_endpoints();
  test_gradient_in_range_valid_and_monotone();
  test_gradient_solid_color();
  test_gradient_interpolates_between_entries();
  test_gradient_first_stop_offset_flat_fills_prefix();
  test_gradient_three_stops_interior_and_flanks();
  test_gradient_hard_stop_is_abrupt();

  test_baked_palette_matches_source_endpoints();
  test_baked_palette_in_range();

  test_procedural_palette_cosine();
  test_mutating_palette_blends_endpoints();
  test_generative_palette_deterministic();
  test_generative_palette_auto_seed_advances();
  test_generative_palette_snapshot_lerp();
  test_palette_modifiers();
  test_static_palette_composition();
  test_palette_wrappers();

  // Clamp-before-cast / NaN-saturation checks, shared with the fast-math pass.
#define HS_RUN_CLAMP_TEST(fn) fn();
  HS_FASTMATH_CLAMP_TESTS(HS_RUN_CLAMP_TEST)
#undef HS_RUN_CLAMP_TEST

  return fixture.result();
}

} // namespace color_tests
} // namespace hs_test

