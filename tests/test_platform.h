/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for the host-side FastLED / Arduino mocks in core/platform.h.
 *
 * These mocks stand in for <FastLED.h> and the Arduino runtime in the native
 * and WASM builds, so any divergence from the device's integer semantics makes
 * the simulator lie about the hardware — and the determinism contract requires
 * the sim to match the device bit-for-bit. We pin golden vectors for the FastLED
 * sine/scale primitives (sin8, sin16, scale8/scale16, beatsin16): hand-traced
 * cardinal anchors plus a full-period accuracy bound against the true sine, so a
 * regression in sin16's byte-truncated secoffset8 (or any LUT/slope drift) is
 * caught even where it still hits the anchors. We also pin the FastLED map8
 * mapping, the degenerate-range guard on map(), and the now-faithful
 * Serial.printf — the seams the old "rough approximation" stubs got wrong.
 */
#pragma once

#include "core/platform.h"
#include "tests/test_harness.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>

namespace hs_test {
namespace platform_tests {

static constexpr double kTwoPi = 6.283185307179586;

// sin8 must be bit-exact with FastLED's sin8_C LUT: centred on 128, peak at
// theta=64, trough near theta=192. Beyond the cardinal anchors, sweep the full
// period and bound the LUT error against the true sine, so a regression in the
// interleaved-slope math is caught even where it still hits the anchors.
inline void test_sin8_golden() {
  // Exact anchors (FastLED sin8_C).
  HS_EXPECT_EQ(sin8(0), 128);
  HS_EXPECT_EQ(sin8(64), 255); // quarter cycle -> peak
  HS_EXPECT_EQ(sin8(128), 128);
  HS_EXPECT_EQ(sin8(192), 1); // three-quarter -> trough

  int max_err = 0;
  for (int t = 0; t < 256; ++t) {
    double truth = 128.0 + 127.0 * std::sin(kTwoPi * t / 256.0);
    int err = std::abs(static_cast<int>(sin8(static_cast<uint8_t>(t))) -
                       static_cast<int>(std::lround(truth)));
    if (err > max_err) max_err = err;
  }
  // sin8_C tracks the true sine to within a few LSBs (measured worst case 3).
  HS_EXPECT_LT(max_err, 6);
}

// sin16 is FastLED's sin16_C: signed -32767..32767, with the byte-truncated
// secoffset8 the review flagged as not self-evidently correct. Pin hand-traced
// cardinal anchors AND bound the full-period LUT error against the true sine, so
// a truncation/slope regression that still happens to hit the anchors is caught.
inline void test_sin16_golden() {
  // Exact anchors traced through sin16_C.
  HS_EXPECT_EQ(sin16(0), 0);
  HS_EXPECT_EQ(sin16(8192), 23170);   // 45 deg
  HS_EXPECT_EQ(sin16(16384), 32645);  // 90 deg (LUT peak, just under 32767)
  HS_EXPECT_EQ(sin16(32768), 0);      // 180 deg
  HS_EXPECT_EQ(sin16(49152), -32645); // 270 deg

  int max_err = 0;
  for (int t = 0; t < 65536; t += 3) {
    double truth = 32767.0 * std::sin(kTwoPi * t / 65536.0);
    int err = std::abs(static_cast<int>(sin16(static_cast<uint16_t>(t))) -
                       static_cast<int>(std::lround(truth)));
    if (err > max_err) max_err = err;
  }
  // sin16_C's 8-section linear LUT tracks the true sine to within ~226 of 32767
  // (~0.7%); a broken secoffset8 truncation or slope table would blow far past.
  HS_EXPECT_LT(max_err, 300);
}

// scale8/scale16 are exact fixed-point fades: (i*sc)>>8 and (i*sc)>>16. Pin the
// full-scale, half-scale, and zero corners that effects rely on.
inline void test_scale_golden() {
  HS_EXPECT_EQ(scale8(255, 255), 254); // 255*255>>8
  HS_EXPECT_EQ(scale8(255, 128), 127);
  HS_EXPECT_EQ(scale8(128, 255), 127);
  HS_EXPECT_EQ(scale8(200, 100), 78); // 20000>>8
  HS_EXPECT_EQ(scale8(1, 255), 0);
  HS_EXPECT_EQ(scale8(0, 255), 0);

  HS_EXPECT_EQ(scale16(65535, 65535), 65534);
  HS_EXPECT_EQ(scale16(65535, 32768), 32767);
  HS_EXPECT_EQ(scale16(32768, 65535), 32767);
  HS_EXPECT_EQ(scale16(1000, 5000), 76); // 5,000,000>>16
  HS_EXPECT_EQ(scale16(0, 65535), 0);
}

inline void test_scale8_fractional() {
  HS_EXPECT_EQ(scale8(0, 200), 0);
  HS_EXPECT_EQ(scale8(255, 255), 254); // 255*255>>8
  HS_EXPECT_EQ(scale8(128, 128), 64);
}

// FastLED map8(in, start, end) maps the FULL 0..255 input onto [start, end] via
// scale8 — the inverse of the old mock, which remapped [start,end]->[0,255] and
// produced an out-of-range (wrapped) palette index for in < start.
inline void test_map8_fastled_semantics() {
  HS_EXPECT_EQ(map8(0, 100, 140), 100);   // input floor -> rangeStart
  HS_EXPECT_EQ(map8(255, 100, 140), 139); // input ceil  -> ~rangeEnd
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

// map() must not divide by zero on a degenerate input range. The device's
// Cortex-M7 yields out_min there; the host now matches instead of SIGFPE-ing.
inline void test_map_degenerate_range() {
  HS_EXPECT_EQ(map(5, 10, 10, 0, 100), 0);    // in_min == in_max -> out_min
  HS_EXPECT_EQ(map(42, 7, 7, 3, 9), 3);
  HS_EXPECT_EQ(map(5, 0, 10, 0, 100), 50);    // normal mapping still correct
}

// beatsin8 oscillates within [lowest, highest] (scale8 fit), is deterministic
// under the injected clock, applies phase_offset (the old mock dropped it), and
// reproduces FastLED's LUT phase at known times.
inline void test_beatsin8_faithful() {
  // 60 BPM == one cycle per 1000 ms. Quarter cycle (250 ms) -> peak, three-
  // quarter (750 ms) -> trough.
  hs::set_mock_time(250, 250000);
  HS_EXPECT_EQ(beatsin8(60, 0, 255), 254); // sin8(64)=255 -> scale8(255,255)
  uint8_t again = beatsin8(60, 0, 255);
  HS_EXPECT_EQ(again, 254); // deterministic under the same mock time

  hs::set_mock_time(750, 750000);
  HS_EXPECT_EQ(beatsin8(60, 0, 255), 0); // sin8(192)=1 -> scale8(1,255)=0

  // phase_offset shifts the wave: at t=0 the bare beat is the LUT midpoint, a
  // +64 offset advances it to the peak. The old mock ignored this 5th argument.
  hs::set_mock_time(0, 0);
  HS_EXPECT_EQ(beatsin8(60, 0, 255, 0, 0), 127);   // sin8(0)=128 -> scale8(128,255)
  HS_EXPECT_EQ(beatsin8(60, 0, 255, 0, 64), 254);  // sin8(64)=255

  // Range fit holds across the whole cycle for a non-trivial [lowest, highest].
  for (unsigned long ms = 0; ms < 1000; ms += 37) {
    hs::set_mock_time(ms, ms * 1000);
    uint8_t v = beatsin8(120, 40, 90);
    HS_EXPECT_GE(v, 40);
    HS_EXPECT_LE(v, 90);
  }
  hs::clear_mock_time();
}

// beatsin16 = lowest + scale16(sin16(beat16(...)) + 32768, highest-lowest). Pin
// it value-exact at t=0 (beat phase 0) while sweeping phase_offset onto the
// sin16 anchors above — this catches the arg-swap / dropped-offset / missing
// +32768 regressions the comment in platform.h warns about — then confirm the
// range fit holds across a full cycle.
inline void test_beatsin16_golden() {
  hs::set_mock_time(0, 0);
  // phase 0:      sin16(0)=0      -> +32768 -> scale16(32768,4000)=2000 -> 3000
  HS_EXPECT_EQ(beatsin16(90, 1000, 5000, 0, 0), 3000);
  // phase 16384:  sin16=32645     -> 65413  -> scale16(65413,4000)=3992 -> 4992
  HS_EXPECT_EQ(beatsin16(90, 1000, 5000, 0, 16384), 4992);
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

// Serial.printf must expand its varargs (Arduino behaviour). The old stub wrote
// the raw format string and dropped every argument.
inline void test_serial_printf_formats_varargs() {
  std::ostringstream cap;
  std::streambuf *saved = std::cout.rdbuf(cap.rdbuf());
  Serial.printf("req %u / cap %u", 12u, 48u);
  std::cout.rdbuf(saved);
  HS_EXPECT(cap.str() == std::string("req 12 / cap 48"), "printf expands args");
}

inline int run_platform_tests() {
  auto scope = hs_test::begin_module("platform");

  test_sin8_golden();
  test_sin16_golden();
  test_scale8_fractional();
  test_scale_golden();
  test_map8_fastled_semantics();
  test_map_degenerate_range();
  test_beatsin8_faithful();
  test_beatsin16_golden();
  test_serial_printf_formats_varargs();

  return hs_test::end_module(scope);
}

} // namespace platform_tests
} // namespace hs_test
