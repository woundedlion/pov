/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file dma_led.h
 * @brief Non-blocking DMA LED controller for HD107S (APA102-compatible) LEDs.
 *
 * Replaces FastLED's blocking SPI with an async DMA pipeline:
 *   HD107SFrame  — pre-formatted protocol buffer with inline color correction
 *   TeensySPIDMA — low-level DMA+SPI hardware driver
 *   DMALEDController — double-buffered high-level controller
 *
 * Color correction uses the existing PROGMEM sRGB↔Linear LUTs from
 * color_luts.h. All corrections are applied in linear 16-bit space.
 *
 * Only compiles on Teensy 4.x (ARDUINO defined). On WASM/sim builds this
 * header is a no-op.
 */

#ifdef ARDUINO

#include <Arduino.h>
#include <SPI.h>
#include <DMAChannel.h>
#include <atomic>
#include "color_luts.h"

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
 *   End frame   : ceil((N+1)/2) bytes of 0x00
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
  /// End-frame size per SK9822/HD107S convention.
  static constexpr int END_FRAME_BYTES = (N / 2) + 1;
  /// Total buffer size in bytes.
  static constexpr int BUFFER_SIZE = 4 + (N * 4) + END_FRAME_BYTES;

  HD107SFrame() {
    memset(buffer_, 0, BUFFER_SIZE);
    // Pre-fill the brightness byte for every pixel slot (0xFF = max).
    for (int i = 0; i < N; ++i) {
      buffer_[4 + i * 4] = 0xFF;
    }
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

      // 2. Color correction in linear space (0-255 = 0.0-1.0 scale)
      r = (r * corrR_) >> 8;
      g = (g * corrG_) >> 8;
      b = (b * corrB_) >> 8;

      // 3. Temperature correction in linear space
      r = (r * tempR_) >> 8;
      g = (g * tempG_) >> 8;
      b = (b * tempB_) >> 8;

      // 4. Brightness scaling in linear space
      r = (r * brightness_) >> 8;
      g = (g * brightness_) >> 8;
      b = (b * brightness_) >> 8;

      // Clamp to 16-bit range (should not exceed with scale-down only)
      if (r > 65535) r = 65535;
      if (g > 65535) g = 65535;
      if (b > 65535) b = 65535;

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
    arm_dcache_flush_delete(buffer_, BUFFER_SIZE);
  }

  const uint8_t* data() const { return buffer_; }
  constexpr size_t size() const { return BUFFER_SIZE; }

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
  uint8_t buffer_[BUFFER_SIZE] __attribute__((aligned(32)));

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

// ============================================================================
// TeensySPIDMA — Low-level async DMA+SPI driver for Teensy 4.x
// ============================================================================

/**
 * @brief Manages a single DMA channel wired to LPSPI4 (SPI) for async
 *        byte-stream transmission.
 *
 * Usage:
 *   spi.transmitAsync(buffer, length);  // returns immediately
 *   // ... do other work ...
 *   while (!spi.isComplete()) {}        // or poll later
 */
class TeensySPIDMA {
public:
  /**
   * @brief Constructs the DMA SPI driver with configurable SPI settings.
   * @param clock  SPI clock frequency in Hz (default: 12 MHz for HD107S).
   * @param order  Bit order (default: MSBFIRST).
   * @param mode   SPI mode (default: SPI_MODE0).
   */
  TeensySPIDMA(uint32_t clock = 12000000, uint8_t order = MSBFIRST,
               uint8_t mode = SPI_MODE0)
      : transferComplete_(true), spiSettings_(clock, order, mode) {}

  /**
   * @brief Initializes SPI and DMA hardware. Must be called from setup(),
   *        not from global constructors (peripherals may not be ready yet).
   */
  void init() {
    instance_ = this;

    SPI.begin();
    SPI.beginTransaction(spiSettings_);

    // Configure LPSPI4 for DMA
    LPSPI4_CFGR1 |= LPSPI_CFGR1_NOSTALL;
    LPSPI4_DER   |= LPSPI_DER_TDDE;  // Enable TX DMA request
    LPSPI4_TCR    = (LPSPI4_TCR & ~LPSPI_TCR_FRAMESZ(31))
                  | LPSPI_TCR_FRAMESZ(7);  // 8-bit frames

    // Configure DMA channel
    dma_.begin(true);
    dma_.destination(reinterpret_cast<volatile uint32_t&>(LPSPI4_TDR));
    dma_.triggerAtHardwareEvent(DMAMUX_SOURCE_LPSPI4_TX);
    dma_.disableOnCompletion();
    dma_.interruptAtCompletion();
    dma_.attachInterrupt(dmaISR);
  }

  /**
   * @brief Starts an async DMA transfer. Returns immediately.
   * @param data Pointer to byte buffer (must remain valid until complete).
   * @param len  Number of bytes to transmit.
   */
  void transmitAsync(const uint8_t* data, size_t len) {
    if (!transferComplete_.load(std::memory_order_relaxed)) {
      waitComplete();
    }
    transferComplete_.store(false, std::memory_order_relaxed);
    dma_.sourceBuffer(data, len);
    dma_.enable();
  }

  bool isComplete() const {
    return transferComplete_.load(std::memory_order_relaxed);
  }

  void waitComplete() {
    while (!transferComplete_.load(std::memory_order_relaxed)) { /* spin */ }
  }

private:
  static void FASTRUN dmaISR() {
    if (instance_) {
      instance_->dma_.clearInterrupt();
      instance_->transferComplete_.store(true, std::memory_order_relaxed);
    }
  }

  DMAChannel dma_;
  std::atomic<bool> transferComplete_;
  SPISettings spiSettings_;

  /**
   * @brief Singleton pointer for ISR callback dispatch.
   *
   * Thread-safety: Safe under the single-core Cortex-M7 ISR preemption model.
   * Only written once from setup() (before interrupts that use it) and read
   * from the DMA completion ISR. No concurrent write-write or read-write
   * race is possible — ARM single-core guarantees that ISR preemption of
   * the main thread sees all prior stores.
   */
  static TeensySPIDMA* instance_;
};

TeensySPIDMA* TeensySPIDMA::instance_ = nullptr;

// ============================================================================
// DMALEDController — Double-buffered async LED controller
// ============================================================================

/**
 * @brief High-level double-buffered DMA LED controller for HD107S strips.
 * @tparam N Number of pixels.
 *
 * Typical ISR usage:
 *   controller.show(leds);    // triggers async DMA, returns immediately
 *   // ISR exits → DMA transfers in background → main loop gets more CPU
 */
template <int N>
class DMALEDController {
public:
  DMALEDController() : activeBuffer_(0), transferCount_(0), overrunCount_(0) {}

  /**
   * @brief One-time hardware initialization. Call from setup().
   */
  void begin() {
    spi_.init();
  }

  /**
   * @brief Applies corrections, loads the back buffer, and triggers DMA.
   * @param pixels Array of CRGB pixels.
   *
   * Non-blocking: returns as soon as DMA is triggered. The previous DMA
   * transfer is guaranteed complete before we start the next one.
   */
  void show(const CRGB* pixels) {
    if (!spi_.isComplete()) {
      overrunCount_++;
    }

    int back = 1 - activeBuffer_;
    frames_[back].load(pixels, N);
    spi_.transmitAsync(frames_[back].data(), frames_[back].size());
    activeBuffer_ = back;
    transferCount_++;
  }

  bool isReady() const { return spi_.isComplete(); }
  void waitForCompletion() { spi_.waitComplete(); }

  // --- Diagnostics ---
  uint32_t getTransferCount() const { return transferCount_; }
  uint32_t getOverrunCount() const { return overrunCount_; }

  // --- Configuration pass-throughs ---

  void setBrightness(uint8_t brightness) {
    HD107SFrame<N>::setBrightness(brightness);
  }

  void setTemperature(uint8_t r, uint8_t g, uint8_t b) {
    HD107SFrame<N>::setTemperature(r, g, b);
  }

  void setCorrection(uint8_t r, uint8_t g, uint8_t b) {
    HD107SFrame<N>::setCorrection(r, g, b);
  }

private:
  HD107SFrame<N> frames_[2];
  TeensySPIDMA spi_;
  int activeBuffer_;
  volatile uint32_t transferCount_;
  volatile uint32_t overrunCount_;
};

#endif // ARDUINO
