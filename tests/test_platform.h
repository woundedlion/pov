/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for the host-side FastLED / Arduino mocks in core/engine/platform.h.
 *
 * These mocks stand in for <FastLED.h> and the Arduino runtime in the native
 * and WASM builds, so any divergence from the device's integer semantics makes
 * the simulator lie about the hardware — and the determinism contract requires
 * the sim to match the device bit-for-bit. We pin golden vectors for the FastLED
 * sine/scale primitives (sin8, sin16, scale8/scale16, beatsin16): hand-traced
 * cardinal anchors plus a full-period accuracy bound against the true sine, so a
 * regression in sin16's byte-truncated secoffset8 (or any LUT/slope drift) is
 * caught even where it still hits the anchors. We also pin the FastLED map8
 * mapping, the degenerate-range guard on map(), and Serial.printf varargs
 * expansion — the seams most prone to diverging from FastLED.
 */
#pragma once

#include "core/engine/platform.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>

namespace hs_test {
namespace platform_tests {

static constexpr double TWO_PI = 6.283185307179586;

/**
 * @brief Verifies sin8 is bit-exact with FastLED's sin8_C LUT and tracks the
 *        true sine to within a few LSBs.
 * @details Pins the cardinal anchors (centred on 128, peak at theta=64, trough
 *          near theta=192), then sweeps the full period and bounds the LUT
 *          error against the true sine, so a regression in the interleaved-slope
 *          math is caught even where it still hits the anchors.
 */
inline void test_sin8_golden() {
  // Exact anchors (FastLED sin8_C).
  HS_EXPECT_EQ(sin8(0), 128);
  HS_EXPECT_EQ(sin8(64), 255); // quarter cycle -> peak
  HS_EXPECT_EQ(sin8(128), 128);
  HS_EXPECT_EQ(sin8(192), 1); // three-quarter -> trough

  int max_err = 0;
  for (int t = 0; t < 256; ++t) {
    double truth = 128.0 + 127.0 * std::sin(TWO_PI * t / 256.0);
    int err = std::abs(static_cast<int>(sin8(static_cast<uint8_t>(t))) -
                       static_cast<int>(std::lround(truth)));
    if (err > max_err) max_err = err;
  }
  // sin8_C tracks the true sine to within a few LSBs (measured worst case 3).
  HS_EXPECT_LT(max_err, 4);
}

/**
 * @brief Verifies sin16 matches FastLED's sin16_C and bounds its full-period
 *        LUT error against the true sine.
 * @details sin16_C is signed -32767..32767 with a byte-truncated secoffset8 that
 *          is not self-evidently correct. Pins hand-traced cardinal anchors AND
 *          bounds the full-period LUT error, so a truncation/slope regression
 *          that still happens to hit the anchors is caught.
 */
inline void test_sin16_golden() {
  // Exact anchors traced through sin16_C.
  HS_EXPECT_EQ(sin16(0), 0);
  HS_EXPECT_EQ(sin16(8192), 23170);   // 45 deg
  HS_EXPECT_EQ(sin16(16384), 32645);  // 90 deg (LUT peak, just under 32767)
  HS_EXPECT_EQ(sin16(32768), 0);      // 180 deg
  HS_EXPECT_EQ(sin16(49152), -32645); // 270 deg

  int max_err = 0;
  for (int t = 0; t < 65536; t += 3) {
    double truth = 32767.0 * std::sin(TWO_PI * t / 65536.0);
    int err = std::abs(static_cast<int>(sin16(static_cast<uint16_t>(t))) -
                       static_cast<int>(std::lround(truth)));
    if (err > max_err) max_err = err;
  }
  // sin16_C's 8-section linear LUT tracks the true sine to within 226 of 32767
  // (~0.7%); a broken secoffset8 truncation or slope table would blow far past.
  HS_EXPECT_LT(max_err, 227);
}

/**
 * @brief Verifies scale8/scale16 implement FastLED's SCALE8_FIXED fades at the
 *        full-scale, half-scale, and zero corners.
 * @details The fades are (i*(1+sc))>>8 and (i*(1+sc))>>16, so a full scale
 *          (sc==max) is the identity. Pins the corners that effects rely on.
 */
inline void test_scale_golden() {
  HS_EXPECT_EQ(scale8(255, 255), 255); // 255*256>>8 — full scale is identity
  HS_EXPECT_EQ(scale8(255, 128), 128); // 255*129>>8
  HS_EXPECT_EQ(scale8(128, 255), 128); // 128*256>>8
  HS_EXPECT_EQ(scale8(200, 100), 78);  // 200*101>>8
  HS_EXPECT_EQ(scale8(1, 255), 1);     // 1*256>>8
  HS_EXPECT_EQ(scale8(0, 255), 0);

  HS_EXPECT_EQ(scale16(65535, 65535), 65535); // full scale is identity
  HS_EXPECT_EQ(scale16(65535, 32768), 32768); // 65535*32769>>16
  HS_EXPECT_EQ(scale16(32768, 65535), 32768); // 32768*65536>>16
  HS_EXPECT_EQ(scale16(1000, 5000), 76); // 1000*5001>>16
  HS_EXPECT_EQ(scale16(0, 65535), 0);
}

/**
 * @brief Verifies scale8 at fractional scales.
 * @details Confirms zero in -> zero, full scale is identity, and the
 *          SCALE8_FIXED rounding (i*(1+sc))>>8 holds at the half point.
 */
inline void test_scale8_fractional() {
  HS_EXPECT_EQ(scale8(0, 200), 0);
  HS_EXPECT_EQ(scale8(255, 255), 255); // 255*256>>8 — full scale is identity
  HS_EXPECT_EQ(scale8(128, 128), 64);  // 128*129>>8
}

/**
 * @brief Verifies FastLED map8(in, start, end) semantics.
 * @details map8 maps the FULL 0..255 input onto [start, end] via scale8, so
 *          every input stays in range. Confirms the floor/ceil endpoints and
 *          that the mapping is monotone non-decreasing.
 */
inline void test_map8_fastled_semantics() {
  HS_EXPECT_EQ(map8(0, 100, 140), 100);   // input floor -> rangeStart
  HS_EXPECT_EQ(map8(255, 100, 140), 140); // input ceil  -> rangeEnd
  // Every input stays inside [start, end] and is monotone non-decreasing.
  int prev = -1;
  for (int in = 0; in <= 255; ++in) {
    uint8_t out = map8(static_cast<uint8_t>(in), 100, 140);
    HS_EXPECT_GE(out, 100);
    HS_EXPECT_LE(out, 140);
    HS_EXPECT_GE(out, prev);
    prev = out;
  }
}

/**
 * @brief Verifies map() does not divide by zero on a degenerate input range.
 * @details The device's Cortex-M7 yields out_min when in_min==in_max, so the
 *          host must match that rather than SIGFPE; normal mapping still holds.
 */
inline void test_map_degenerate_range() {
  HS_EXPECT_EQ(map(5, 10, 10, 0, 100), 0);    // in_min == in_max -> out_min
  HS_EXPECT_EQ(map(42, 7, 7, 3, 9), 3);
  HS_EXPECT_EQ(map(5, 0, 10, 0, 100), 50);    // normal mapping still correct
}

/**
 * @brief Verifies random()/random8() do not divide by zero on a degenerate
 *        range and that a normal range stays in bounds.
 * @details Two cases with different provenance:
 *          - random(0) / random(<=0) / random(min, min) are device-faithful: the
 *            empty range reduces to Arduino's random(0) -> 0 (i.e. min), which is
 *            well-defined on the device.
 *          - random(min, max) with min > max is NOT something the Arduino runtime
 *            guarantees (its long-range modulo of a negative span is effectively
 *            unspecified). The host mock defines "inverted -> min" purely to dodge
 *            the modulo-by-zero/UB; this test pins that mock contract, not a
 *            cross-device behavior.
 */
inline void test_random_degenerate_range() {
  HS_EXPECT_EQ(random(0), 0);      // Arduino random(0) -> 0 (device-faithful)
  HS_EXPECT_EQ(random(-3), 0);     // non-positive bound -> 0 (device-faithful)
  HS_EXPECT_EQ(random(5, 5), 5);   // empty range -> min (device-faithful)
  HS_EXPECT_EQ(random(9, 4), 9);   // inverted range -> min (host-mock contract)
  HS_EXPECT_EQ(random8(0), 0);     // FastLED random8(0) -> 0
  // A normal range still falls inside [min, max).
  for (int i = 0; i < 64; ++i) {
    int v = random(3, 7);
    HS_EXPECT_GE(v, 3);
    HS_EXPECT_LE(v, 6);
  }
}

/**
 * @brief Verifies beatsin8 oscillates within [lowest, highest], is
 *        deterministic, applies phase_offset, and matches FastLED's LUT phase.
 * @details Confirms the scale8 range fit, determinism under the injected clock,
 *          the phase_offset (5th arg) shift, and the LUT phase at known times.
 */
inline void test_beatsin8_faithful() {
  // 60 BPM == one cycle per 1000 ms. Quarter cycle (250 ms) -> peak, three-
  // quarter (750 ms) -> trough.
  hs::set_mock_time(250, 250000);
  HS_EXPECT_EQ(beatsin8(60, 0, 255), 255); // sin8(64)=255 -> scale8(255,255)=255
  uint8_t again = beatsin8(60, 0, 255);
  HS_EXPECT_EQ(again, 255); // deterministic under the same mock time

  hs::set_mock_time(750, 750000);
  HS_EXPECT_EQ(beatsin8(60, 0, 255), 1); // sin8(192)=1 -> scale8(1,255)=1

  // phase_offset (5th arg) shifts the wave: at t=0 the bare beat is the LUT
  // midpoint, a +64 offset advances it to the peak.
  hs::set_mock_time(0, 0);
  HS_EXPECT_EQ(beatsin8(60, 0, 255, 0, 0), 128);   // sin8(0)=128 -> scale8(128,255)=128
  HS_EXPECT_EQ(beatsin8(60, 0, 255, 0, 64), 255);  // sin8(64)=255 -> scale8(255,255)=255

  // Range fit holds across the whole cycle for a non-trivial [lowest, highest].
  for (unsigned long ms = 0; ms < 1000; ms += 37) {
    hs::set_mock_time(ms, ms * 1000);
    uint8_t v = beatsin8(120, 40, 90);
    HS_EXPECT_GE(v, 40);
    HS_EXPECT_LE(v, 90);
  }
  hs::clear_mock_time();
}

/**
 * @brief Verifies beatsin16 is value-exact at known phases and stays in range
 *        across a full cycle.
 * @details beatsin16 = lowest + scale16(sin16(beat16(...)) + 32768,
 *          highest-lowest). Pins it value-exact at t=0 (beat phase 0) while
 *          sweeping phase_offset onto the sin16 anchors, catching an arg swap, a
 *          dropped offset, or a missing +32768, then confirms the range fit.
 */
inline void test_beatsin16_golden() {
  hs::set_mock_time(0, 0);
  // phase 0:      sin16(0)=0      -> +32768 -> scale16(32768,4000)=2000 -> 3000
  HS_EXPECT_EQ(beatsin16(90, 1000, 5000, 0, 0), 3000);
  // phase 16384:  sin16=32645     -> 65413  -> scale16(65413,4000)=3993 -> 4993
  HS_EXPECT_EQ(beatsin16(90, 1000, 5000, 0, 16384), 4993);
  // phase 49152:  sin16=-32645    -> 123    -> scale16(123,4000)=7      -> 1007
  HS_EXPECT_EQ(beatsin16(90, 1000, 5000, 0, 49152), 1007);

  for (unsigned long ms = 0; ms < 1000; ms += 53) {
    hs::set_mock_time(ms, ms * 1000);
    uint16_t v = beatsin16(90, 1000, 5000);
    HS_EXPECT_GE(v, 1000);
    HS_EXPECT_LE(v, 5000);
  }
  hs::clear_mock_time();
}

/**
 * @brief Verifies Serial.printf expands its varargs (Arduino behaviour) rather
 *        than emitting the raw format string.
 */
inline void test_serial_printf_formats_varargs() {
  std::ostringstream cap;
  std::streambuf *saved = std::cout.rdbuf(cap.rdbuf());
  Serial.printf("req %u / cap %u", 12u, 48u);
  std::cout.rdbuf(saved);
  HS_EXPECT(cap.str() == std::string("req 12 / cap 48"), "printf expands args");
}

/**
 * @brief Verifies rand_f() stays in the half-open interval [0, 1) and that
 *        random_to_unit() clamps only the top band.
 * @details The top band of RNG draws rounds value/max to exactly 1.0f, which
 *          would make (int)(rand_f()*N) index N out of bounds; random_to_unit()
 *          clamps those to the float just below 1.0f and leaves the rest
 *          unchanged.
 */
inline void test_rand_f_half_open() {
  constexpr uint32_t MAX = 0xFFFFFFFFu; // hs::Pcg32::max()
  // The boundary case: the naive ratio is 1.0f, the clamp pulls it to the float
  // just below.
  HS_EXPECT_LT(hs::random_to_unit(MAX, MAX), 1.0f);
  HS_EXPECT_EQ(hs::random_to_unit(MAX, MAX), 0x1.fffffep-1f);
  HS_EXPECT_LT(hs::random_to_unit(MAX - 64, MAX), 1.0f); // rest of the top band
  // Ordinary draws pass through unchanged.
  HS_EXPECT_EQ(hs::random_to_unit(0, MAX), 0.0f);
  HS_EXPECT_NEAR(hs::random_to_unit(MAX / 2, MAX), 0.5f, 1e-6f);
  // And the live generator never escapes [0, 1) across a long sweep.
  for (int i = 0; i < 100000; ++i) {
    float r = hs::rand_f();
    HS_EXPECT_GE(r, 0.0f);
    HS_EXPECT_LT(r, 1.0f);
  }
}

/**
 * @brief Verifies CRGB's single-argument constructor decodes a 0xRRGGBB
 *        colorcode like FastLED rather than as a grayscale fill.
 * @details CRGB(0xFF8000) must be orange, not black.
 */
inline void test_crgb_colorcode_constructor() {
  CRGB orange(0xFF8000u);
  HS_EXPECT_EQ(orange.r, 0xFF);
  HS_EXPECT_EQ(orange.g, 0x80);
  HS_EXPECT_EQ(orange.b, 0x00);
  CRGB teal(0x008080u);
  HS_EXPECT_EQ(teal.r, 0x00);
  HS_EXPECT_EQ(teal.g, 0x80);
  HS_EXPECT_EQ(teal.b, 0x80);
}

/**
 * @brief Runs every platform-mock test under the "platform" module scope.
 * @return The module's failure count.
 */
inline int run_platform_tests() {
  hs_test::ModuleFixture fixture("platform");

  test_sin8_golden();
  test_sin16_golden();
  test_scale8_fractional();
  test_scale_golden();
  test_map8_fastled_semantics();
  test_map_degenerate_range();
  test_random_degenerate_range();
  test_rand_f_half_open();
  test_crgb_colorcode_constructor();
  test_beatsin8_faithful();
  test_beatsin16_golden();
  test_serial_printf_formats_varargs();

  return fixture.result();
}

} // namespace platform_tests
} // namespace hs_test
