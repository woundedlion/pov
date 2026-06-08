/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file hd107s_frame.h
 * @brief Pre-formatted HD107S (APA102-compatible) protocol buffer + color
 *        correction. Pure CPU logic, no Teensy peripherals — split out of
 *        dma_led.h so the wire-format and correction math are host-unit-testable
 *        (see tests/test_hd107s_frame.h). dma_led.h #includes this and adds the
 *        DMA/SPI hardware driver on top (Teensy-only).
 *
 * Color correction uses the PROGMEM sRGB<->Linear LUTs from color.h; all
 * corrections are applied in linear 16-bit space.
 */

#include <cstddef>
#include <cstdint>
#include <cstring>

#include "color.h" // Pixel16, CRGB, srgb_to_linear_lut / linear_to_srgb_lut

// arm_dcache_flush_delete is a Teensy cache intrinsic (Arduino.h). It only
// matters when a real DMA engine reads the buffer; on host builds (unit tests,
// WASM sim) there is no DMA, so the flush is a no-op.
#ifndef ARDUINO
static inline void arm_dcache_flush_delete(void *, size_t) {}
#endif

// ============================================================================
// HD107SFrame — Pre-formatted DMA buffer for the HD107S protocol
// ============================================================================

/**
 * @brief Pre-formatted DMA buffer for HD107S (APA102-compatible) LEDs.
 * @tparam N Maximum number of pixels.
 *
 * HD107S frame layout:
 *   Start frame : 4 bytes of 0x00
 *   Per pixel   : [0xFF] [B] [G] [R]   (brightness byte fixed at max)
 *   End frame   : ceil(N/16) bytes of 0x00. Each LED re-clocks the data one
 *                 half-cycle later, so the last pixel needs ceil(N/2) extra
 *                 clocks to latch; at 8 clocks/byte that is ceil(N/16) =
 *                 (N+15)/16 bytes. 0x00 (not 0xFF) per the SK9822/HD107S latch.
 *
 * Color correction pipeline (all in linear 16-bit space):
 *   1. sRGB 8-bit → linear 16-bit   (srgb_to_linear_lut, PROGMEM)
 *   2. Color correction multiply     (TypicalLEDStrip equivalent)
 *   3. Temperature correction multiply (Candle equivalent)
 *   4. Brightness scaling
 *   5. Linear 16-bit → sRGB 8-bit   (linear_to_srgb_lut, PROGMEM)
 */
template <int N>
class HD107SFrame {
public:
  /// End-frame latch per the SK9822/HD107S spec: ceil(N/2) extra clocks to
  /// push data through the chain → ceil(N/16) = (N+15)/16 bytes of 0x00.
  static constexpr int END_FRAME_BYTES = (N + 15) / 16;
  /// Single-frame buffer size in bytes.
  static constexpr int BUFFER_SIZE = 4 + (N * 4) + END_FRAME_BYTES;
  /// Composite size: image frame + trailing black frame (for show_bg).
  static constexpr int COMPOSITE_SIZE = BUFFER_SIZE * 2;

  // The whole composite buffer can be handed to a single DMA transfer
  // (submitFrame(withBg=true) → transmitAsync(data(), COMPOSITE_SIZE)). Teensy
  // 4's eDMA encodes a transfer's major-loop count in the 15-bit CITER/BITER
  // field (minor-loop linking disabled), so one transfer tops out at 32767
  // bytes — past that the count silently truncates and the strip tail goes
  // dark. Trap at compile time if a future pixel count would overflow it.
  static_assert(COMPOSITE_SIZE <= 32767,
      "HD107SFrame composite buffer exceeds the 15-bit eDMA single-transfer "
      "limit (CITER/BITER); split the transfer or reduce N");

  HD107SFrame() {
    memset(buffer_, 0, COMPOSITE_SIZE);
    // Image frame: set brightness bytes
    for (int i = 0; i < N; ++i) {
      buffer_[4 + i * 4] = 0xFF;
    }
    // Trailing black frame: set brightness bytes (colors stay zero)
    for (int i = 0; i < N; ++i) {
      buffer_[BUFFER_SIZE + 4 + i * 4] = 0xFF;
    }
  }

  /**
   * @brief Apply color → temperature → brightness correction in linear 16-bit
   * space, in place, then clamp to [0, 65535].
   *
   * Shared by load() (CRGB path) and packPixel() (Pixel16 path) so the two
   * cannot drift. r/g/b are already-linear channels. FASTRUN/inline because
   * packPixel() calls it on the per-column ISR hot path; with scale-down-only
   * factors the clamp never triggers but is kept for safety.
   */
  FASTRUN inline void correct(uint32_t& r, uint32_t& g, uint32_t& b) const {
    // Color correction
    r = (r * corrR_) >> 8;
    g = (g * corrG_) >> 8;
    b = (b * corrB_) >> 8;

    // Temperature correction
    r = (r * tempR_) >> 8;
    g = (g * tempG_) >> 8;
    b = (b * tempB_) >> 8;

    // Brightness scaling
    r = (r * brightness_) >> 8;
    g = (g * brightness_) >> 8;
    b = (b * brightness_) >> 8;

    if (r > 65535) r = 65535;
    if (g > 65535) g = 65535;
    if (b > 65535) b = 65535;
  }

  /**
   * @brief Loads pixel data into the buffer with full color correction.
   * @param pixels Source CRGB array (sRGB 8-bit).
   * @param count  Number of pixels to load (clamped to N).
   *
   * All corrections are applied in linear 16-bit space using the PROGMEM
   * sRGB↔Linear LUTs. Result is written in HD107S byte order: [0xFF][B][G][R].
   */
  void load(const CRGB* pixels, int count) {
    if (count > N) count = N;

    uint8_t* dest = buffer_ + 4; // skip start frame
    for (int i = 0; i < count; ++i) {
      const CRGB& c = pixels[i];

      // 1. sRGB 8-bit → Linear 16-bit (PROGMEM LUT)
      uint32_t r = srgb_to_linear_lut[c.r];
      uint32_t g = srgb_to_linear_lut[c.g];
      uint32_t b = srgb_to_linear_lut[c.b];

      // 2-4. Color + temperature + brightness correction (linear), clamped.
      correct(r, g, b);

      // 5. Linear 16-bit → sRGB 8-bit (PROGMEM LUT)
      uint8_t r8 = linear_to_srgb_lut[r];
      uint8_t g8 = linear_to_srgb_lut[g];
      uint8_t b8 = linear_to_srgb_lut[b];

      // HD107S byte order: [0xFF][B][G][R]
      dest[0] = 0xFF;
      dest[1] = b8;
      dest[2] = g8;
      dest[3] = r8;
      dest += 4;
    }

    // Flush data cache so DMA sees the updated buffer.
    arm_dcache_flush_delete(buffer_, COMPOSITE_SIZE);
  }

  /**
   * @brief Packs a single Pixel16 directly into the buffer with corrections.
   * @param index LED index (0-based).
   * @param p     Linear 16-bit pixel.
   *
   * Applies color/temperature/brightness corrections in linear 16-bit space
   * then converts to sRGB 8-bit in a single pass (no intermediate CRGB).
   */
  FASTRUN inline void packPixel(int index, const Pixel16& p) {
    uint8_t* dest = buffer_ + 4 + index * 4;

    uint32_t r = p.r;
    uint32_t g = p.g;
    uint32_t b = p.b;

    // Color + temperature + brightness correction (all linear 16-bit), clamped.
    correct(r, g, b);

    // Single linear 16-bit → sRGB 8-bit conversion, packed in wire order
    dest[0] = 0xFF;
    dest[1] = linear_to_srgb_lut[b];
    dest[2] = linear_to_srgb_lut[g];
    dest[3] = linear_to_srgb_lut[r];
  }

  /**
   * @brief Flushes data cache so DMA sees the buffer. Call after packPixel().
   */
  void flush() {
    arm_dcache_flush_delete(buffer_, COMPOSITE_SIZE);
  }

  const uint8_t* data() const { return buffer_; }
  constexpr size_t size() const { return BUFFER_SIZE; }
  /// Size including trailing black frame (for show_bg composite DMA).
  constexpr size_t sizeWithBg() const { return COMPOSITE_SIZE; }

  // --- Static correction configuration (shared across all frames) -----------

  static void setTemperature(uint8_t r, uint8_t g, uint8_t b) {
    tempR_ = r; tempG_ = g; tempB_ = b;
  }

  static void setCorrection(uint8_t r, uint8_t g, uint8_t b) {
    corrR_ = r; corrG_ = g; corrB_ = b;
  }

  static void setBrightness(uint8_t brightness) {
    brightness_ = brightness;
  }

private:
  uint8_t buffer_[COMPOSITE_SIZE] __attribute__((aligned(32)));

  // Shared correction state — 8-bit scale factors (255 = 1.0)
  static uint8_t tempR_, tempG_, tempB_;
  static uint8_t corrR_, corrG_, corrB_;
  static uint8_t brightness_;
};

// Static member definitions
template <int N> uint8_t HD107SFrame<N>::tempR_ = 255;
template <int N> uint8_t HD107SFrame<N>::tempG_ = 255;
template <int N> uint8_t HD107SFrame<N>::tempB_ = 255;
template <int N> uint8_t HD107SFrame<N>::corrR_ = 255;
template <int N> uint8_t HD107SFrame<N>::corrG_ = 255;
template <int N> uint8_t HD107SFrame<N>::corrB_ = 255;
template <int N> uint8_t HD107SFrame<N>::brightness_ = 255;
