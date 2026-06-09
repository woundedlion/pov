/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for the host-side FastLED / Arduino mocks in core/platform.h.
 *
 * These mocks stand in for <FastLED.h> and the Arduino runtime in the native
 * and WASM builds, so any divergence from the device's integer semantics makes
 * the simulator lie about the hardware. We pin the LUT-exact sine, the scale8
 * range fit, the FastLED map8 mapping, the degenerate-range guard on map(), and
 * the now-faithful Serial.printf — the seams the old "rough approximation"
 * stubs got wrong.
 */
#pragma once

#include "core/platform.h"
#include "tests/test_harness.h"

#include <sstream>
#include <string>

namespace hs_test {
namespace platform_tests {

// sin8 must be bit-exact with FastLED's sin8_C LUT: centred on 128, peak at
// theta=64, trough near theta=192.
inline void test_sin8_lut_exact() {
  HS_EXPECT_EQ(sin8(0), 128);
  HS_EXPECT_EQ(sin8(64), 255); // quarter cycle -> peak
  HS_EXPECT_EQ(sin8(128), 128);
  HS_EXPECT_EQ(sin8(192), 1); // three-quarter -> trough
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

inline void test_beatsin16_in_range() {
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

  test_sin8_lut_exact();
  test_scale8_fractional();
  test_map8_fastled_semantics();
  test_map_degenerate_range();
  test_beatsin8_faithful();
  test_beatsin16_in_range();
  test_serial_printf_formats_varargs();

  return hs_test::end_module(scope);
}

} // namespace platform_tests
} // namespace hs_test
