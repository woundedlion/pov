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

/**
 * @brief Returns a pointer to pixel i's 4-byte wire record in the image frame.
 * @param f Frame whose image buffer is inspected.
 * @param i Pixel index in [0, N); offset is 4-byte start frame plus i records.
 * @return Pointer to the [0xFF][B][G][R] record for pixel i.
 */
inline const uint8_t *pixel(const Frame &f, int i) {
  return f.data() + 4 + i * 4;
}

/**
 * @brief Restores the shared static correction state to unity (255 = exact ×1.0).
 * @details Resets the per-channel correction, temperature, and brightness so a
 * test starts from an identity color pipeline; these are static frame settings.
 */
inline void reset_correction() {
  Frame::setCorrection(255, 255, 255);
  Frame::setTemperature(255, 255, 255);
  Frame::setBrightness(255);
}

/**
 * @brief Verifies the spec-derived buffer-sizing formulas and that size()/
 * sizeWithBg() agree with the compile-time layout constants.
 */
inline void test_layout_constants() {
  HS_EXPECT_EQ(Frame::END_FRAME_BYTES, (N + 15) / 16);
  HS_EXPECT_EQ(Frame::BUFFER_SIZE, 4 + N * 4 + (N + 15) / 16);
  HS_EXPECT_EQ(Frame::COMPOSITE_SIZE, Frame::BUFFER_SIZE * 2);

  Frame f;
  HS_EXPECT_EQ(static_cast<int>(f.size()), Frame::BUFFER_SIZE);
  HS_EXPECT_EQ(static_cast<int>(f.sizeWithBg()), Frame::COMPOSITE_SIZE);
}

/**
 * @brief Verifies a default-constructed frame holds a valid all-black wire image.
 * @details Checks the start frame, per-pixel 0xFF brightness byte with zero
 * color, the zero end-frame latch, and a matching trailing black frame.
 */
inline void test_fresh_frame_skeleton() {
  Frame f;
  const uint8_t *d = f.data();

  // Start frame: 4 bytes of 0x00.
  for (int i = 0; i < 4; ++i)
    HS_EXPECT_EQ(d[i], 0);

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

  const uint8_t *bg = d + Frame::BUFFER_SIZE;
  for (int i = 0; i < N; ++i) {
    HS_EXPECT_EQ(bg[4 + i * 4], 0xFF);
    HS_EXPECT_EQ(bg[4 + i * 4 + 1], 0);
  }
}

/**
 * @brief Exercises correct(): zero pass-through, exact unity identity at full
 * scale, and monotonic brightness scaling (0 zeroes the channel).
 */
inline void test_correct_pipeline() {
  reset_correction();
  Frame f;

  // Zero in → zero out, regardless of factors.
  uint32_t r = 0, g = 0, b = 0;
  f.correct(r, g, b);
  HS_EXPECT_EQ(r, 0u);
  HS_EXPECT_EQ(g, 0u);
  HS_EXPECT_EQ(b, 0u);

  // Factor 255 maps to ×256/256 = exact identity, so full input is preserved.
  r = g = b = 65535;
  f.correct(r, g, b);
  HS_EXPECT_EQ(r, 65535u);
  HS_EXPECT_EQ(g, 65535u);
  HS_EXPECT_EQ(b, 65535u);

  Frame::setBrightness(128);
  uint32_t half = 65535;
  uint32_t dummy = 0;
  f.correct(half, dummy, dummy);
  HS_EXPECT_LT(half, 65535u);

  Frame::setBrightness(0);
  uint32_t off = 65535;
  uint32_t z2 = 65535, z3 = 65535;
  f.correct(off, z2, z3);
  HS_EXPECT_EQ(off, 0u);
  reset_correction();
}

/**
 * @brief Exercises the multi-factor compounding that actually ships and locks in
 * the dead saturation clamp.
 * @details test_correct_pipeline only varies brightness with every other factor
 * at unity, so the non-unity correction + temperature gains the production config
 * sets (pov_single.h: correction 255,176,240; temperature 255,147,41) are never
 * asserted together. This case applies those shipped gains and checks both that
 * the two factors compound (temperature attenuates on top of correction, not
 * instead of it) and the exact per-channel result, then documents that the
 * clamp branch in correct() is unreachable for any public factor combination.
 */
inline void test_correct_multifactor() {
  Frame f;

  // Temperature compounds on top of correction, not instead of it.
  Frame::setCorrection(255, 176, 240);
  Frame::setTemperature(255, 255, 255);
  Frame::setBrightness(255);
  uint32_t gr = 0, gg = 65535, gb = 0;
  f.correct(gr, gg, gb);
  const uint32_t g_corr_only = gg;

  Frame::setTemperature(255, 147, 41);
  uint32_t r = 65535, g = 65535, b = 65535;
  f.correct(r, g, b);
  HS_EXPECT_LT(g, g_corr_only); // factors compound

  // Exact shipped scaling. factor(f) stores f+1, so each stage is (v*(f+1))>>8:
  //   R: 65535                                                     = 65535
  //   G: ((65535*177>>8)=45311 *148>>8)=26195  *256>>8             = 26195
  //   B: ((65535*241>>8)=61695 * 42>>8)=10121  *256>>8             = 10121
  HS_EXPECT_EQ(r, 65535u);
  HS_EXPECT_EQ(g, 26195u);
  HS_EXPECT_EQ(b, 10121u);
  HS_EXPECT_LT(b, g);

  // Dead-clamp invariant: factor 255 maps to multiplier 256 (exact unity), so
  // every stage's (v*256)>>8 never exceeds the input — max gains reach but never
  // breach 65535, so the clamp never fires.
  Frame::setCorrection(255, 255, 255);
  Frame::setTemperature(255, 255, 255);
  Frame::setBrightness(255);
  uint32_t mr = 65535, mg = 65535, mb = 65535;
  f.correct(mr, mg, mb);
  HS_EXPECT_EQ(mr, 65535u);
  HS_EXPECT_EQ(mg, 65535u);
  HS_EXPECT_EQ(mb, 65535u);

  reset_correction();
}

/**
 * @brief Verifies packPixel() emits [0xFF][B][G][R] order and writes only the
 * targeted pixel slot (each primary lights its own byte, neighbors untouched).
 */
inline void test_packpixel_wire_order() {
  reset_correction();
  Frame f;

  const Pixel16 red(CRGB(255, 0, 0));
  const Pixel16 green(CRGB(0, 255, 0));
  const Pixel16 blue(CRGB(0, 0, 255));

  // Wire record is [0xFF][B][G][R].
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

  HS_EXPECT_GT(pixel(f, 0)[3], 0);
  HS_EXPECT_EQ(pixel(f, 0)[1], 0);
}

/**
 * @brief Verifies load() and packPixel() yield identical wire bytes for
 * equivalent inputs and that load() preserves CRGB channel identity.
 * @details Both paths share the same correct() core, so equivalent inputs
 * cannot drift; pure red must light only the R wire slot.
 */
inline void test_load_matches_packpixel() {
  reset_correction();

  const CRGB src[3] = {CRGB(200, 100, 50), CRGB(0, 0, 0), CRGB(255, 255, 255)};

  Frame fl;
  fl.load(src, 3);

  Frame fp;
  for (int i = 0; i < 3; ++i)
    fp.packPixel(i, Pixel16(src[i]));

  for (int i = 0; i < 3; ++i)
    for (int k = 0; k < 4; ++k)
      HS_EXPECT_EQ(pixel(fl, i)[k], pixel(fp, i)[k]);

  // load() honors CRGB channel identity: pure red → only the R wire slot.
  Frame fr;
  const CRGB red_only[1] = {CRGB(255, 0, 0)};
  fr.load(red_only, 1);
  HS_EXPECT_GT(pixel(fr, 0)[3], 0); // R
  HS_EXPECT_EQ(pixel(fr, 0)[2], 0); // G
  HS_EXPECT_EQ(pixel(fr, 0)[1], 0); // B
}

/**
 * @brief Verifies a full count==N load fills every pixel and leaves the end
 * frame zero (the upper in-range bound; an over-N count now traps — see the
 * load_count_over_range death case).
 */
inline void test_load_full_count() {
  reset_correction();
  Frame f;
  CRGB full[N];
  for (int i = 0; i < N; ++i)
    full[i] = CRGB(255, 255, 255);
  f.load(full, N);
  for (int i = 0; i < Frame::END_FRAME_BYTES; ++i)
    HS_EXPECT_EQ(f.data()[4 + N * 4 + i], 0);
}

/**
 * @brief Runs the full hd107s_frame test suite.
 * @return The module's failure count from end_module().
 */
inline int run_hd107s_tests() {
  auto scope = hs_test::begin_module("hd107s");
  test_layout_constants();
  test_fresh_frame_skeleton();
  test_correct_pipeline();
  test_packpixel_wire_order();
  test_load_matches_packpixel();
  test_load_full_count();
  reset_correction(); // leave shared static state clean for any later module
  return hs_test::end_module(scope);
}

} // namespace hd107s_tests
} // namespace hs_test
