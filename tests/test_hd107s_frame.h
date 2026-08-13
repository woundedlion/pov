/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the HD107S protocol buffer + color correction
 * (hardware/hd107s_frame.h), split out of dma_led.h precisely so this
 * wire-format and arithmetic is testable without a Teensy. Covers the bytes
 * that actually go on the SPI wire: frame layout, [0xFF][B][G][R] channel
 * order, and the linear-space correction pipeline.
 * The DMA/SPI driver itself (register access, eDMA, ISR) is hardware-only and
 * out of host-test scope.
 */
#pragma once

#include "hardware/hd107s_frame.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cstdint>

namespace hs_test {
namespace hd107s_tests {

constexpr int N = 40; // small strip: END_FRAME_BYTES=4, COMPOSITE=336 bytes

// Shipping Phantasm segment: TOTAL_PIXELS 288 over NUM_SEGMENTS 4. First size
// whose ceil(N/16) end frame clears the 4-byte pad, so END_FRAME_BYTES=8.
constexpr int PHANTASM_PPS = 72;

using Frame = HD107SFrame<N>;
using PhantasmFrame = HD107SFrame<PHANTASM_PPS>;

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
 * @brief Restores HD107SFrame<S>'s shared static correction state to unity
 *        (255 = exact ×1.0).
 * @tparam S Pixel count of the frame instantiation to reset; the correction
 *           state is a static of HD107SFrame<S>, so each S carries its own.
 * @details Resets the per-channel correction, temperature, and brightness so a
 * test starts from an identity color pipeline; these are static frame settings.
 * Shared with the DMA controller suite, which drives a different S.
 */
template <int S> inline void reset_correction() {
  HD107SFrame<S>::setCorrection(255, 255, 255);
  HD107SFrame<S>::setTemperature(255, 255, 255);
  HD107SFrame<S>::setBrightness(255);
}

/**
 * @brief Verifies the spec-derived buffer-sizing formulas and that size()/
 * sizeWithBg() agree with the compile-time layout constants.
 */
inline void test_layout_constants() {
  HS_EXPECT_EQ(Frame::END_FRAME_BYTES, 4);
  HS_EXPECT_EQ(Frame::BUFFER_SIZE, 4 + N * 4 + Frame::END_FRAME_BYTES);
  HS_EXPECT_EQ(Frame::COMPOSITE_SIZE, Frame::BUFFER_SIZE * 2);

  Frame f;
  HS_EXPECT_EQ(static_cast<int>(f.size()), Frame::BUFFER_SIZE);
  HS_EXPECT_EQ(static_cast<int>(f.sizeWithBg()), Frame::COMPOSITE_SIZE);
  HS_EXPECT_EQ(reinterpret_cast<uintptr_t>(f.data() + Frame::BUFFER_SIZE) %
                   alignof(uint32_t),
               0u);
}

/**
 * @brief Verifies the layout constants at the shipping Phantasm segment size.
 * @details N=40 and the dma_controller module's N=8 both land on the 4-byte
 * end-frame floor, so a ceil-div slip that only shows past the first
 * multiple-of-4 boundary stays invisible. 72 pixels is the size the device
 * actually instantiates and the first one whose end frame is 8 bytes.
 */
inline void test_layout_constants_phantasm() {
  HS_EXPECT_EQ(PhantasmFrame::END_FRAME_BYTES, 8);
  HS_EXPECT_EQ(PhantasmFrame::BUFFER_SIZE,
               4 + PHANTASM_PPS * 4 + PhantasmFrame::END_FRAME_BYTES);
  HS_EXPECT_EQ(PhantasmFrame::COMPOSITE_SIZE, PhantasmFrame::BUFFER_SIZE * 2);

  PhantasmFrame f;
  HS_EXPECT_EQ(static_cast<int>(f.size()), PhantasmFrame::BUFFER_SIZE);
  HS_EXPECT_EQ(static_cast<int>(f.sizeWithBg()), PhantasmFrame::COMPOSITE_SIZE);
  HS_EXPECT_EQ(
      reinterpret_cast<uintptr_t>(f.data() + PhantasmFrame::BUFFER_SIZE) %
          alignof(uint32_t),
      0u);
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
 * scale, and exact brightness scaling at half and zero.
 */
inline void test_correct_pipeline() {
  reset_correction<N>();
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

  // factor(128) is 129, so each channel is (65535*129)>>8.
  Frame::setBrightness(128);
  uint32_t hr = 65535, hg = 65535, hb = 65535;
  f.correct(hr, hg, hb);
  HS_EXPECT_EQ(hr, 33023u);
  HS_EXPECT_EQ(hg, 33023u);
  HS_EXPECT_EQ(hb, 33023u);

  // factor(0) is 0, so every channel zeroes regardless of input.
  Frame::setBrightness(0);
  uint32_t zr = 65535, zg = 65535, zb = 65535;
  f.correct(zr, zg, zb);
  HS_EXPECT_EQ(zr, 0u);
  HS_EXPECT_EQ(zg, 0u);
  HS_EXPECT_EQ(zb, 0u);
  reset_correction<N>();
}

/**
 * @brief Exercises the multi-factor compounding that actually ships, then the
 * no-overflow invariant at maximum gain.
 * @details test_correct_pipeline only varies brightness with every other factor
 * at unity, so the non-unity correction + temperature gains the production config
 * sets (pov_single.h: correction 255,176,240; temperature 255,147,41) are never
 * asserted together. This case applies those shipped gains and checks both that
 * the two factors compound (temperature attenuates on top of correction, not
 * instead of it) and the exact per-channel result, then that the largest public
 * factor combination leaves a full-scale input unchanged.
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

  // No-overflow invariant: factor 255 maps to multiplier 256 (exact unity), so
  // every stage's (v*256)>>8 returns the input untouched — max gains reach but
  // never breach 65535, keeping every stage a valid linear_to_srgb_lut index.
  Frame::setCorrection(255, 255, 255);
  Frame::setTemperature(255, 255, 255);
  Frame::setBrightness(255);
  uint32_t mr = 65535, mg = 65535, mb = 65535;
  f.correct(mr, mg, mb);
  HS_EXPECT_EQ(mr, 65535u);
  HS_EXPECT_EQ(mg, 65535u);
  HS_EXPECT_EQ(mb, 65535u);

  reset_correction<N>();
}

/**
 * @brief Verifies packPixel() emits [0xFF][B][G][R] order and writes only the
 * targeted pixel slot (each primary lights its own byte, neighbors untouched).
 */
inline void test_packpixel_wire_order() {
  reset_correction<N>();
  Frame f;

  const Pixel red(CRGB(255, 0, 0));
  const Pixel green(CRGB(0, 255, 0));
  const Pixel blue(CRGB(0, 0, 255));

  // Wire record is [0xFF][B][G][R]. Under unity correction the lit channel
  // round-trips the sRGB<->linear LUTs back to 255; the rest stay 0.
  f.packPixel(0, red);
  HS_EXPECT_EQ(pixel(f, 0)[0], 0xFF);
  HS_EXPECT_EQ(pixel(f, 0)[1], 0);   // B
  HS_EXPECT_EQ(pixel(f, 0)[2], 0);   // G
  HS_EXPECT_EQ(pixel(f, 0)[3], 255); // R

  f.packPixel(1, green);
  HS_EXPECT_EQ(pixel(f, 1)[1], 0);   // B
  HS_EXPECT_EQ(pixel(f, 1)[2], 255); // G
  HS_EXPECT_EQ(pixel(f, 1)[3], 0);   // R

  f.packPixel(2, blue);
  HS_EXPECT_EQ(pixel(f, 2)[1], 255); // B
  HS_EXPECT_EQ(pixel(f, 2)[2], 0);   // G
  HS_EXPECT_EQ(pixel(f, 2)[3], 0);   // R

  HS_EXPECT_EQ(pixel(f, 0)[3], 255);
  HS_EXPECT_EQ(pixel(f, 0)[1], 0);
}

/**
 * @brief Verifies packPixel()'s wire bytes under the shipped factor set at
 * non-unity brightness.
 * @details Full-scale white through correction 255,176,240, temperature
 * 255,147,41 and brightness 128 leaves linear (33023, 13199, 5100), which
 * linear_to_srgb8 encodes as R=188, G=124, B=79.
 */
inline void test_packpixel_shipped_brightness() {
  Frame::setCorrection(255, 176, 240);
  Frame::setTemperature(255, 147, 41);
  Frame::setBrightness(128);
  Frame f;

  const Pixel white(CRGB(255, 255, 255));
  f.packPixel(3, white);
  HS_EXPECT_EQ(pixel(f, 3)[0], 0xFF);
  HS_EXPECT_EQ(pixel(f, 3)[1], 79);  // B
  HS_EXPECT_EQ(pixel(f, 3)[2], 124); // G
  HS_EXPECT_EQ(pixel(f, 3)[3], 188); // R

  reset_correction<N>();
}

/**
 * @brief Runs the full hd107s_frame test suite.
 * @return The module's failure count from end_module().
 */
inline int run_hd107s_tests() {
  hs_test::ModuleFixture fixture("hd107s");
  test_layout_constants();
  test_layout_constants_phantasm();
  test_fresh_frame_skeleton();
  test_correct_pipeline();
  test_correct_multifactor();
  test_packpixel_wire_order();
  test_packpixel_shipped_brightness();
  reset_correction<N>(); // leave shared static state clean for any later module
  return fixture.result();
}

} // namespace hd107s_tests
} // namespace hs_test
