/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the HD107S protocol buffer + color correction
 * (hardware/hd107s_frame.h), split out of dma_led.h precisely so this
 * wire-format and arithmetic is testable without a Teensy. Covers the bytes
 * that actually go on the SPI wire: frame layout, [0xFF][B][G][R] channel
 * order, the linear-space correction pipeline, and load()/packPixel() parity.
 * The DMA/SPI driver itself (register access, eDMA, ISR) is hardware-only and
 * out of host-test scope.
 */
#pragma once

#include "hardware/hd107s_frame.h"
#include "tests/test_harness.h"

#include <cstdint>

namespace hs_test {
namespace hd107s_tests {

constexpr int N = 40; // small strip: END_FRAME_BYTES=3, COMPOSITE=334 bytes

using Frame = HD107SFrame<N>;

// Pointer to pixel i's 4-byte record [0xFF][B][G][R] in the image frame.
inline const uint8_t *pixel(const Frame &f, int i) {
  return f.data() + 4 + i * 4;
}

// Restore the shared static correction state to unity (255 = scale ×255/256).
inline void reset_correction() {
  Frame::setCorrection(255, 255, 255);
  Frame::setTemperature(255, 255, 255);
  Frame::setBrightness(255);
}

inline void test_layout_constants() {
  // The end-frame latch and buffer sizing are spec-derived; pin the formulas.
  HS_EXPECT_EQ(Frame::END_FRAME_BYTES, (N + 15) / 16);
  HS_EXPECT_EQ(Frame::BUFFER_SIZE, 4 + N * 4 + (N + 15) / 16);
  HS_EXPECT_EQ(Frame::COMPOSITE_SIZE, Frame::BUFFER_SIZE * 2);

  Frame f;
  HS_EXPECT_EQ(static_cast<int>(f.size()), Frame::BUFFER_SIZE);
  HS_EXPECT_EQ(static_cast<int>(f.sizeWithBg()), Frame::COMPOSITE_SIZE);
}

inline void test_fresh_frame_skeleton() {
  Frame f;
  const uint8_t *d = f.data();

  // Start frame: 4 bytes of 0x00.
  for (int i = 0; i < 4; ++i)
    HS_EXPECT_EQ(d[i], 0);

  // Every pixel slot carries the fixed 0xFF brightness byte and zero color.
  for (int i = 0; i < N; ++i) {
    const uint8_t *p = pixel(f, i);
    HS_EXPECT_EQ(p[0], 0xFF);
    HS_EXPECT_EQ(p[1], 0);
    HS_EXPECT_EQ(p[2], 0);
    HS_EXPECT_EQ(p[3], 0);
  }

  // End frame: zeros (SK9822/HD107S latch is 0x00, not 0xFF).
  for (int i = 0; i < Frame::END_FRAME_BYTES; ++i)
    HS_EXPECT_EQ(d[4 + N * 4 + i], 0);

  // Trailing black frame mirrors the image frame's brightness bytes, zero color.
  const uint8_t *bg = d + Frame::BUFFER_SIZE;
  for (int i = 0; i < N; ++i) {
    HS_EXPECT_EQ(bg[4 + i * 4], 0xFF);
    HS_EXPECT_EQ(bg[4 + i * 4 + 1], 0);
  }
}

inline void test_correct_pipeline() {
  reset_correction();
  Frame f;

  // Zero in → zero out, regardless of factors.
  uint32_t r = 0, g = 0, b = 0;
  f.correct(r, g, b);
  HS_EXPECT_EQ(r, 0u);
  HS_EXPECT_EQ(g, 0u);
  HS_EXPECT_EQ(b, 0u);

  // Unity factors are ×255/256 applied three times (corr·temp·brightness), so
  // the pipeline is scale-down-only — output never exceeds input, never clamps.
  auto unity = [](uint32_t v) {
    v = (v * 255) >> 8;
    v = (v * 255) >> 8;
    v = (v * 255) >> 8;
    return v;
  };
  r = g = b = 65535;
  f.correct(r, g, b);
  HS_EXPECT_EQ(r, unity(65535));
  HS_EXPECT_LE(r, 65535u); // clamp invariant holds (and is never reached here)
  HS_EXPECT_EQ(g, unity(65535));
  HS_EXPECT_EQ(b, unity(65535));

  // Brightness scales down monotonically; 0 zeroes the channel.
  Frame::setBrightness(128);
  uint32_t half = 65535;
  uint32_t dummy = 0;
  f.correct(half, dummy, dummy);
  HS_EXPECT_LT(half, unity(65535)); // dimmer than full brightness

  Frame::setBrightness(0);
  uint32_t off = 65535;
  uint32_t z2 = 65535, z3 = 65535;
  f.correct(off, z2, z3);
  HS_EXPECT_EQ(off, 0u);
  reset_correction();
}

inline void test_packpixel_wire_order() {
  reset_correction();
  Frame f;

  // Linear primaries via the sRGB→linear conversion CRGB carries.
  const Pixel16 red(CRGB(255, 0, 0));
  const Pixel16 green(CRGB(0, 255, 0));
  const Pixel16 blue(CRGB(0, 0, 255));

  // Wire record is [0xFF][B][G][R]: a primary must light exactly its own slot.
  f.packPixel(0, red);
  HS_EXPECT_EQ(pixel(f, 0)[0], 0xFF);
  HS_EXPECT_EQ(pixel(f, 0)[1], 0);    // B
  HS_EXPECT_EQ(pixel(f, 0)[2], 0);    // G
  HS_EXPECT_GT(pixel(f, 0)[3], 0);    // R

  f.packPixel(1, green);
  HS_EXPECT_EQ(pixel(f, 1)[1], 0);    // B
  HS_EXPECT_GT(pixel(f, 1)[2], 0);    // G
  HS_EXPECT_EQ(pixel(f, 1)[3], 0);    // R

  f.packPixel(2, blue);
  HS_EXPECT_GT(pixel(f, 2)[1], 0);    // B
  HS_EXPECT_EQ(pixel(f, 2)[2], 0);    // G
  HS_EXPECT_EQ(pixel(f, 2)[3], 0);    // R

  // Packing pixel 2 must not have disturbed pixel 0.
  HS_EXPECT_GT(pixel(f, 0)[3], 0);
  HS_EXPECT_EQ(pixel(f, 0)[1], 0);
}

inline void test_load_matches_packpixel() {
  reset_correction();

  // load() (sRGB CRGB → linear → correct → sRGB) and packPixel() (linear
  // Pixel16 → correct → sRGB) share correct() and must produce identical wire
  // bytes for equivalent inputs — the "cannot drift" contract.
  const CRGB src[3] = {CRGB(200, 100, 50), CRGB(0, 0, 0), CRGB(255, 255, 255)};

  Frame fl;
  fl.load(src, 3);

  Frame fp;
  for (int i = 0; i < 3; ++i)
    fp.packPixel(i, Pixel16(src[i]));

  for (int i = 0; i < 3; ++i)
    for (int k = 0; k < 4; ++k)
      HS_EXPECT_EQ(pixel(fl, i)[k], pixel(fp, i)[k]);

  // And load() honors CRGB channel identity: pure red → only the R wire slot.
  Frame fr;
  const CRGB red_only[1] = {CRGB(255, 0, 0)};
  fr.load(red_only, 1);
  HS_EXPECT_GT(pixel(fr, 0)[3], 0); // R
  HS_EXPECT_EQ(pixel(fr, 0)[2], 0); // G
  HS_EXPECT_EQ(pixel(fr, 0)[1], 0); // B
}

inline void test_load_clamps_count() {
  reset_correction();
  Frame f;
  // count > N must clamp, not overrun the buffer (the end frame stays zero).
  CRGB big[N + 8];
  for (int i = 0; i < N + 8; ++i)
    big[i] = CRGB(255, 255, 255);
  f.load(big, N + 8);
  for (int i = 0; i < Frame::END_FRAME_BYTES; ++i)
    HS_EXPECT_EQ(f.data()[4 + N * 4 + i], 0);
}

inline int run_hd107s_tests() {
  auto scope = hs_test::begin_module("hd107s_frame");
  test_layout_constants();
  test_fresh_frame_skeleton();
  test_correct_pipeline();
  test_packpixel_wire_order();
  test_load_matches_packpixel();
  test_load_clamps_count();
  reset_correction(); // leave shared static state clean for any later module
  return hs_test::end_module(scope);
}

} // namespace hd107s_tests
} // namespace hs_test
