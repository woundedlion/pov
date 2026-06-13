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
 *        (see tests/test_hd107s_frame.h). dma_led.h `#includes` this and adds the
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
  /** End-frame latch per the SK9822/HD107S spec: ceil(N/2) extra clocks to
       push data through the chain → ceil(N/16) = (N+15)/16 bytes of 0x00. */
  static constexpr int END_FRAME_BYTES = (N + 15) / 16;
  /** Single-frame buffer size in bytes. */
  static constexpr int BUFFER_SIZE = 4 + (N * 4) + END_FRAME_BYTES;
  /** Composite size: image frame + trailing black frame (for show_bg). */
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

  /**
   * @brief Constructs the frame: zeroes the composite buffer, then primes every
   *        pixel's brightness byte to 0xFF (max).
   * @details Primes brightness bytes in both the image frame and the trailing
   *          black frame. Only the B/G/R color bytes change at runtime;
   *          start/end frames stay zero.
   */
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
   * @brief Applies color → temperature → brightness correction in linear 16-bit
   *        space, in place, then clamps to [0, 65535].
   * @param r In/out red channel, already-linear 16-bit (0..65535).
   * @param g In/out green channel, already-linear 16-bit (0..65535).
   * @param b In/out blue channel, already-linear 16-bit (0..65535).
   * @details Shared by load() (CRGB path) and packPixel() (Pixel16 path) so the
   *          two cannot drift. FASTRUN/inline because packPixel() calls it on
   *          the per-column ISR hot path. The final clamp to 65535 is a hard
   *          memory-safety bound, not a fidelity tweak: the result indexes
   *          linear_to_srgb_lut[r] (a 65536-entry table), so an out-of-range r
   *          would read past the LUT. With the shipped scale-down-only factors
   *          it never triggers, but a >1.0 correction factor would, and the
   *          clamp keeps the index in bounds rather than corrupting memory.
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
   * @details All corrections are applied in linear 16-bit space using the
   *          PROGMEM sRGB↔Linear LUTs. Result is written in HD107S byte order:
   *          [0xFF][B][G][R].
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

    // Flush data cache so DMA sees the updated buffer. Full COMPOSITE_SIZE on
    // purpose even for an image-only transfer — see flush()'s note.
    arm_dcache_flush_delete(buffer_, COMPOSITE_SIZE);
  }

  /**
   * @brief Packs a single Pixel16 directly into the buffer with corrections.
   * @param index LED index (0-based).
   * @param p     Linear 16-bit pixel.
   * @details Applies color/temperature/brightness corrections in linear 16-bit
   *          space then converts to sRGB 8-bit in a single pass (no intermediate
   *          CRGB).
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
   * @details Flushes the full COMPOSITE_SIZE (image + trailing black frame)
   *          even though a non-withBg submitFrame DMAs only the image half:
   *          flush()/load() run before submitFrame and do not know which range
   *          it will pick, so flushing the superset keeps the cache coherent for
   *          either transfer. The extra cost is one dcache pass over the bg half;
   *          under-flushing would instead risk DMAing a stale composite buffer.
   */
  void flush() {
    arm_dcache_flush_delete(buffer_, COMPOSITE_SIZE);
  }

  /**
   * @brief Returns a pointer to the start of the composite DMA buffer.
   * @return Read-only pointer to the first byte of the composite buffer.
   */
  const uint8_t* data() const { return buffer_; }
  /**
   * @brief Returns the size of a single image frame in bytes.
   * @return Image-frame size in bytes (excludes the trailing black frame).
   */
  constexpr size_t size() const { return BUFFER_SIZE; }
  /**
   * @brief Returns the composite size including the trailing black frame.
   * @return Composite size in bytes (image frame + black frame, for show_bg DMA).
   */
  constexpr size_t sizeWithBg() const { return COMPOSITE_SIZE; }

  // --- Static correction configuration (shared across all frames) -----------

  /**
   * @brief Sets the white-balance temperature factors (shared across frames).
   * @param r Red temperature factor, 8-bit scale (255 = ×1.0, 0 = off).
   * @param g Green temperature factor, 8-bit scale (255 = ×1.0, 0 = off).
   * @param b Blue temperature factor, 8-bit scale (255 = ×1.0, 0 = off).
   */
  static void setTemperature(uint8_t r, uint8_t g, uint8_t b) {
    tempR_ = factor(r); tempG_ = factor(g); tempB_ = factor(b);
  }

  /**
   * @brief Sets the per-channel color-correction factors (shared across frames).
   * @param r Red correction factor, 8-bit scale (255 = ×1.0, 0 = off).
   * @param g Green correction factor, 8-bit scale (255 = ×1.0, 0 = off).
   * @param b Blue correction factor, 8-bit scale (255 = ×1.0, 0 = off).
   */
  static void setCorrection(uint8_t r, uint8_t g, uint8_t b) {
    corrR_ = factor(r); corrG_ = factor(g); corrB_ = factor(b);
  }

  /**
   * @brief Sets the global brightness factor (shared across frames).
   * @param brightness Brightness factor, 8-bit scale (255 = full, 0 = off).
   */
  static void setBrightness(uint8_t brightness) {
    brightness_ = factor(brightness);
  }

private:
  uint8_t buffer_[COMPOSITE_SIZE] __attribute__((aligned(32)));

  /**
   * @brief Converts a public 8-bit scale factor into the internal multiplier
   *        used by correct()'s `(v * f) >> 8`.
   * @param f Public 8-bit scale factor (255 = ×1.0, 0 = off).
   * @return Internal multiplier (256 = exact unity, 0 = off).
   * @details Storing factor+1 makes 255 resolve to ×256/256 = exact unity (a
   *          bare `(v*255)>>8` loses ~0.4% per stage, ~1.2% compounded over the
   *          three stages — full brightness would be unreachable in a pipeline
   *          whose selling point is exact 16-bit fidelity). 0 is kept at 0 so it
   *          still zeroes the channel exactly; the cost lives here, off the ISR
   *          hot path.
   */
  static constexpr uint16_t factor(uint8_t f) {
    return f == 0 ? 0u : static_cast<uint16_t>(f) + 1u;
  }

  // Shared correction state — internal multipliers (256 = ×1.0, 0 = off).
  static uint16_t tempR_, tempG_, tempB_;
  static uint16_t corrR_, corrG_, corrB_;
  static uint16_t brightness_;
};

// Static member definitions (256 = unity; see factor()).
template <int N> uint16_t HD107SFrame<N>::tempR_ = 256;
template <int N> uint16_t HD107SFrame<N>::tempG_ = 256;
template <int N> uint16_t HD107SFrame<N>::tempB_ = 256;
template <int N> uint16_t HD107SFrame<N>::corrR_ = 256;
template <int N> uint16_t HD107SFrame<N>::corrG_ = 256;
template <int N> uint16_t HD107SFrame<N>::corrB_ = 256;
template <int N> uint16_t HD107SFrame<N>::brightness_ = 256;
