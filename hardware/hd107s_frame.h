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

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>


#include "core/platform.h" // HS_CHECK, hs::log, FASTRUN — used directly below;
                           // included explicitly rather than relying on color.h
                           // pulling it (this header is independently host-tested)
#include "core/color.h" // Pixel16, CRGB, srgb_to_linear_lut / linear_to_srgb_lut

// arm_dcache_flush (Arduino.h) cleans dirty D-cache lines without invalidating:
// a TX-only buffer must reach RAM but stay resident for the next frame's write.
// No DMA on host builds, so the flush is a no-op there.
#ifndef ARDUINO
static inline void arm_dcache_flush(void *, size_t) {}
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
 *                 clocks to latch. At 8 clocks/byte that is ceil(ceil(N/2)/8)
 *                 bytes, which equals ceil(N/16) = (N+15)/16 *exactly* by the
 *                 nested-ceiling identity ceil(ceil(x/a)/b) = ceil(x/(a·b)) for
 *                 integer a,b>0 — not a coincidence of the two ceilings. The only
 *                 slack is byte granularity: a whole-byte end frame can carry up
 *                 to 7 more zero clocks than the ceil(N/2) minimum at small N,
 *                 which is harmless (extra 0x00 clocks just keep latching zeros).
 *                 0x00 (not 0xFF) per the SK9822/HD107S latch.
 *
 * Color correction pipeline (all in linear 16-bit space):
 *   1. sRGB 8-bit → linear 16-bit   (srgb_to_linear_lut, PROGMEM)
 *   2. Color correction multiply     (TypicalLEDStrip equivalent)
 *   3. Temperature correction multiply (Candle equivalent)
 *   4. Brightness scaling
 *   5. Linear 16-bit → sRGB 8-bit   (linear_to_srgb_lut, PROGMEM)
 *
 * PROGMEM is a placement no-op on Teensy 4.x (Cortex-M7, unified flash/RAM
 * address space) and empty on the PC/wasm sim; the LUTs are read with plain
 * operator[], not pgm_read_*. Direct subscripting is intentional.
 */
template <int N>
class HD107SFrame {
public:
  /** End-frame latch per the SK9822/HD107S spec: ceil(N/2) extra clocks to push
       data through the chain, packed 8 clocks/byte → ceil(ceil(N/2)/8) =
       ceil(N/16) = (N+15)/16 bytes of 0x00 (exact via the nested-ceiling
       identity; see the layout note above for the byte-granularity slack). */
  static constexpr int END_FRAME_BYTES = (N + 15) / 16;
  /** Single-frame buffer size in bytes. */
  static constexpr int BUFFER_SIZE = 4 + (N * 4) + END_FRAME_BYTES;
  /** Composite size: image frame + trailing black frame (for strobe_columns). */
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
    // Prime the 0xFF brightness byte of every pixel slot in the image frame
    // (base 4) and the trailing black frame (base BUFFER_SIZE+4). load() and
    // packPixel() never rewrite it, so this is its sole writer.
    const int bases[2] = {4, BUFFER_SIZE + 4};
    for (int base : bases) {
      for (int i = 0; i < N; ++i) {
        buffer_[base + i * 4] = 0xFF;
      }
    }
  }

  /**
   * @brief Applies color → temperature → brightness correction in linear 16-bit
   *        space, in place.
   * @param r In/out red channel, already-linear 16-bit (0..65535).
   * @param g In/out green channel, already-linear 16-bit (0..65535).
   * @param b In/out blue channel, already-linear 16-bit (0..65535).
   * @details Shared by load() (CRGB path) and packPixel() (Pixel16 path) so the
   *          two cannot drift. FASTRUN/inline because packPixel() calls it on
   *          the per-column ISR hot path. No output clamp is needed: factor()
   *          caps every multiplier at 256 (×1.0), so each (v * f) >> 8 stage with
   *          v ≤ 65535 yields ≤ v. The result therefore never grows past the
   *          16-bit input and stays a valid index into linear_to_srgb_lut (a
   *          65536-entry table). The public API cannot express a >1.0 factor; if
   *          one were ever added (by widening factor()'s domain), a clamp to
   *          65535 would become live here and must be reinstated.
   */
  FASTRUN inline void correct(uint32_t& r, uint32_t& g, uint32_t& b) const {
    r = (r * corrR_) >> 8;
    g = (g * corrG_) >> 8;
    b = (b * corrB_) >> 8;

    r = (r * tempR_) >> 8;
    g = (g * tempG_) >> 8;
    b = (b * tempB_) >> 8;

    r = (r * brightness_) >> 8;
    g = (g * brightness_) >> 8;
    b = (b * brightness_) >> 8;
  }

  /**
   * @brief Loads pixel data into the buffer with full color correction.
   * @param pixels Source CRGB array (sRGB 8-bit).
   * @param count  Number of pixels to load; must be in [0, N].
   * @details All corrections are applied in linear 16-bit space using the
   *          PROGMEM sRGB↔Linear LUTs. Result is written in HD107S byte order:
   *          [0xFF][B][G][R]. This cold bind seam traps an out-of-range count
   *          rather than silently clamping (count > N) or no-op'ing (count < 0),
   *          mirroring packPixel()'s index guard.
   *
   *          Slots [count, N) are blanked to [0xFF][0][0][0], so a partial load()
   *          leaves a clean frame with no stale tail.
   */
  void load(const CRGB* pixels, int count) {
    HS_CHECK(count >= 0 && count <= N, "load count out of range");

    uint8_t* dest = buffer_ + 4; // skip start frame
    for (int i = 0; i < count; ++i) {
      const CRGB& c = pixels[i];

      uint32_t r = srgb_to_linear_lut[c.r];
      uint32_t g = srgb_to_linear_lut[c.g];
      uint32_t b = srgb_to_linear_lut[c.b];

      correct(r, g, b);

      uint8_t r8 = linear_to_srgb_lut[r];
      uint8_t g8 = linear_to_srgb_lut[g];
      uint8_t b8 = linear_to_srgb_lut[b];

      // HD107S byte order: [0xFF][B][G][R]; brightness primed at construction
      dest[1] = b8;
      dest[2] = g8;
      dest[3] = r8;
      dest += 4;
    }
    for (int i = count; i < N; ++i) {
      dest[1] = dest[2] = dest[3] = 0;
      dest += 4;
    }

    arm_dcache_flush(buffer_, COMPOSITE_SIZE);
  }

  /**
   * @brief Packs a single Pixel16 directly into the buffer with corrections.
   * @param index LED index, must be in [0, N). Unchecked on the device hot path:
   *              the write is `buffer_ + 4 + index*4` with no clamp, so an
   *              out-of-range index corrupts the trailing black frame or runs off
   *              the composite buffer (UB). Callers own the bound; the cold
   *              effect-height HS_CHECK at the bind site backstops the shipped
   *              caller. The assert below is a host-build trip-wire for future
   *              callers — stripped on the device by NDEBUG (zero ISR cost).
   *              Deliberately a stripped `assert`, NOT the always-on `HS_CHECK`
   *              that `load()` uses for its cold one-shot bind: this is the
   *              per-pixel ISR hot path, where an always-on trap would cost
   *              cycles on every LED. Do not "harmonize" the two onto HS_CHECK.
   * @param p     Linear 16-bit pixel.
   * @details Applies color/temperature/brightness corrections in linear 16-bit
   *          space then converts to sRGB 8-bit in a single pass (no intermediate
   *          CRGB).
   */
  FASTRUN inline void packPixel(int index, const Pixel16& p) {
    assert(index >= 0 && index < N);
    uint8_t* dest = buffer_ + 4 + index * 4;

    uint32_t r = p.r;
    uint32_t g = p.g;
    uint32_t b = p.b;

    correct(r, g, b);

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
   *          Uses arm_dcache_flush (clean, no invalidate): the buffer is TX-only
   *          (DMA reads, never writes it), so the lines need only be written back
   *          to RAM, not evicted — keeping them resident saves the next frame's
   *          re-read/re-write a cache miss.
   */
  void flush() {
    arm_dcache_flush(buffer_, COMPOSITE_SIZE);
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
   * @return Composite size in bytes (image frame + black frame, for strobe_columns DMA).
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
   * @details Maps f to (f+1)/256, not the nominal f/255: exact at the endpoints
   *          (0→off, 255→×1.0) and ~1/256 high at intermediate values. Storing
   *          f+1 keeps full brightness reachable as exact ×256/256 unity; 0 stays
   *          0 to zero the channel exactly.
   */
  static constexpr uint16_t factor(uint8_t f) {
    return f == 0 ? 0u : static_cast<uint16_t>(f) + 1u;
  }

  static_assert(factor(255) <= 256u,
                "factor() must cap at unity (256) so correct()'s (v*f)>>8 stages "
                "cannot grow v past 16 bits and over-read linear_to_srgb_lut");

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
