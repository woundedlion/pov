/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file dma_led.h
 * @brief Non-blocking DMA LED controller for HD107S (APA102-compatible) LEDs.
 *
 * Drives the strip via an async DMA pipeline so the main loop never blocks on
 * SPI:
 *   HD107SFrame  — pre-formatted protocol buffer with inline color correction
 *   TeensySPIDMA — low-level DMA+SPI hardware driver
 *   DMALEDController — double-buffered high-level controller
 *
 * Color correction uses the PROGMEM sRGB↔Linear LUTs from color_luts.h.
 * All corrections are applied in linear 16-bit space.
 *
 * Only compiles on Teensy 4.x (ARDUINO defined). On WASM/sim builds this
 * header is a no-op.
 */

#ifdef ARDUINO

#include <Arduino.h>
#include <SPI.h>
#include <DMAChannel.h>
#include <atomic>
#include "core/engine/platform.h" // HS_CHECK / hs::log used below; explicit, not via color.h
#include "core/color/color.h"

// HD107SFrame (protocol buffer + color correction) lives in its own header so
// its wire-format / correction math is host-unit-testable without the Teensy
// peripherals below. See tests/test_hd107s_frame.h.
#include "hd107s_frame.h"

// Pure framing / transfer-length / stale-transfer math, host-unit-tested
// without the Teensy peripherals below. See tests/test_dma_core.h.
#include "dma_led_core.h"

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
   *
   * Call exactly once; a second call traps (re-running SPI.begin() would leak
   * the open transaction).
   */
  void init() {
    HS_CHECK(!initialized_, "TeensySPIDMA::init() called twice");
    HS_CHECK(instance_ == nullptr,
             "TeensySPIDMA: only one instance supported (shared DMA-completion "
             "ISR singleton)");
    initialized_ = true;
    instance_ = this;

    SPI.begin();
    SPI.beginTransaction(spiSettings_);

    LPSPI4_CFGR1 |= LPSPI_CFGR1_NOSTALL;
    LPSPI4_DER   |= LPSPI_DER_TDDE;  // TX DMA request
    LPSPI4_TCR    = (LPSPI4_TCR & ~LPSPI_TCR_FRAMESZ(31))
                  | LPSPI_TCR_FRAMESZ(7);  // 8-bit frames

    dma_.begin(true);
    // The uint8_t cast is load-bearing: the 32-bit destination() overload sets
    // ATTR_DST to 32-bit transfers, conflicting with sourceBuffer()'s NBYTES=1
    // and raising an eDMA config error on enable — nothing transmits.
    dma_.destination(reinterpret_cast<volatile uint8_t&>(LPSPI4_TDR));
    dma_.triggerAtHardwareEvent(DMAMUX_SOURCE_LPSPI4_TX);
    dma_.disableOnCompletion();
    dma_.interruptAtCompletion();
    dma_.attachInterrupt(dmaISR);
  }

  /**
   * @brief Starts an async DMA transfer. Returns immediately.
   * @param data Pointer to byte buffer (must remain valid until complete).
   * @param len  Number of bytes to transmit.
   * @pre The caller has cleaned the buffer from cache (submitFrame() does this
   *      via frames_[back].flush()); this method only enables the DMA and does
   *      not flush.
   */
  void transmitAsync(const uint8_t* data, size_t len) {
    // Trap rather than spin on an in-flight transfer: this runs in the column
    // ISR, where spinning would deadlock — the DMA-completion ISR that clears
    // transferComplete_ cannot preempt an equal/lower-priority ISR.
    if (!transferComplete_.load(std::memory_order_relaxed)) {
      hs::log("FATAL: transmitAsync entered with a transfer still in flight — "
              "submitFrame() must guard with isComplete()");
      __builtin_trap();
    }
    transferComplete_.store(false, std::memory_order_relaxed);
    transferStartUs_ = micros();
    dma_.sourceBuffer(data, len);
    dma_.enable();
  }

  /**
   * @brief Reports whether the in-flight transfer has finished.
   * @return true once the in-flight transfer's completion ISR has fired.
   */
  bool isComplete() const {
    return transferComplete_.load(std::memory_order_relaxed);
  }

  /**
   * @brief Surfaces a permanently wedged DMA channel from the overrun-drop path.
   * @details submitFrame() drops on overrun rather than spinning, so a channel
   *          whose completion ISR never fires would otherwise stay masked
   *          forever. Traps once the in-flight transfer outlives the watchdog.
   *          Only fires while submitFrame() is being called.
   */
  void checkStaleTransfer() {
    if (transferComplete_.load(std::memory_order_relaxed)) return;
    if (dma::transfer_stale(transferStartUs_, micros(), kTransferWatchdogUs)) {
      hs::log("FATAL: DMA channel wedged — in-flight transfer outlived the "
              "watchdog on the overrun-drop path; completion ISR never fired");
      __builtin_trap();
    }
  }

private:
  /**
   * @brief Watchdog bound for checkStaleTransfer(), in µs. Far above any real
   *        full-frame DMA (single-digit ms), so only a wedged channel trips it.
   */
  static constexpr unsigned long kTransferWatchdogUs = 100000UL;

  /**
   * @brief DMA completion ISR: clears the interrupt and marks the transfer done.
   * @details Dispatched via the singleton instance_; runs in interrupt context.
   */
  static void FASTRUN dmaISR() {
    if (instance_) {
      instance_->dma_.clearInterrupt();
      instance_->transferComplete_.store(true, std::memory_order_relaxed);
    }
  }

  DMAChannel dma_; /**< eDMA channel wired to LPSPI4 TX for async transmission. */
  /**
   * @brief Completion flag handed between the DMA-completion ISR and the
   *        main/column thread.
   * @details Single-observer (ISR preempts the thread), so relaxed ordering
   *          suffices. Does NOT order the buffer for the DMA engine — the
   *          caller's arm_dcache_flush() before transmitAsync() does that.
   */
  std::atomic<bool> transferComplete_;
  /**
   * @brief micros() at which the in-flight transfer was enabled.
   * @details Touched only in column-ISR context (transmitAsync /
   *          checkStaleTransfer), never the completion ISR, so a plain scalar is
   *          correct. Only meaningful while transferComplete_ is false.
   */
  unsigned long transferStartUs_ = 0;
  SPISettings spiSettings_; /**< Cached SPI clock/bit-order/mode for this driver. */

  /**
   * @brief Single-init guard set once on the first (setup-time) init() call.
   * @details init() touches global SPI/DMA peripheral state that is not safe to
   *          re-run.
   */
  bool initialized_ = false;

  /**
   * @brief Singleton pointer for ISR callback dispatch. Exactly one
   *        TeensySPIDMA per image; init() traps on a second instance.
   * @details Written once from setup() before any ISR uses it, then read from
   *          the completion ISR — race-free under the single-observer model (see
   *          transferComplete_).
   */
  static TeensySPIDMA* instance_;
};

inline TeensySPIDMA* TeensySPIDMA::instance_ = nullptr;

// ============================================================================
// DMALEDController — Double-buffered async LED controller
// ============================================================================

/**
 * @brief High-level double-buffered DMA LED controller for HD107S strips.
 * @tparam N Number of pixels.
 * @note One instance per firmware image: it drives the singleton TeensySPIDMA
 *       backing the shared DMA-completion ISR, so a second begin() traps.
 *
 * Typical ISR usage (per column):
 *   auto& f = controller.backFrame();  // back buffer (not being DMA'd)
 *   // ... pack pixels into f via packPixel() ...
 *   controller.submitFrame();          // triggers async DMA, returns immediately
 *   // ISR exits → DMA transfers in background → main loop gets more CPU
 */
template <int N>
class DMALEDController {
public:
  /**
   * @brief Constructs the controller, optionally overriding the SPI clock.
   * @param clock SPI clock in Hz forwarded to TeensySPIDMA (default: 12 MHz).
   *              The Phantasm driver passes 24 MHz (see pov_segmented.h).
   */
  explicit DMALEDController(uint32_t clock = 12000000)
      : spi_(clock), activeBuffer_(0), transferCount_(0), overrunCount_(0) {}

  /**
   * @brief One-time hardware initialization. Call from setup().
   */
  void begin() {
    spi_.init();
  }

  /**
   * @brief Returns the back frame (not currently being DMA'd).
   * @return Reference to the back-buffer frame; pack pixels via packPixel(),
   *         then call submitFrame().
   */
  HD107SFrame<N>& backFrame() {
    return frames_[dma::next_buffer(activeBuffer_)];
  }

  /**
   * @brief Flushes the back frame and triggers async DMA transfer.
   * Call after all packPixel() calls on backFrame().
   * @param withBg If true, DMAs the composite buffer (image + trailing
   *               black frame) in a single transfer — zero gap, no spin.
   * @return true if the frame was handed to the DMA engine; false if dropped on
   *         overrun (prior transfer still in flight). The fail-dark latch gates
   *         on this; the steady-state column path ignores it (self-heals).
   */
  [[nodiscard]] bool submitFrame(bool withBg = false) {
    if (!spi_.isComplete()) {
      // Drop on overrun. A transfer that NEVER completes is a wedged channel,
      // not a transient, so surface it here — the drop path is where it shows.
      spi_.checkStaleTransfer();
      overrunCount_.fetch_add(1, std::memory_order_relaxed);
      return false;
    }
    int back = dma::next_buffer(activeBuffer_);
    frames_[back].flush();
    size_t len = dma::transfer_len(frames_[back].size(),
                                   frames_[back].sizeWithBg(), withBg);
    spi_.transmitAsync(frames_[back].data(), len);
    activeBuffer_ = back;
    transferCount_.fetch_add(1, std::memory_order_relaxed);
    return true;
  }

  // --- Diagnostics ---
  /**
   * @brief Returns the count of frames handed to the DMA engine since start.
   * @return Monotonic transfer counter (number of successful submitFrame()s).
   */
  uint32_t getTransferCount() const {
    return transferCount_.load(std::memory_order_relaxed);
  }
  /**
   * @brief Returns the count of frames dropped on overrun since start.
   * @return Monotonic overrun counter (frames dropped because a prior transfer
   *         was still in flight).
   */
  uint32_t getOverrunCount() const {
    return overrunCount_.load(std::memory_order_relaxed);
  }

  // --- Configuration pass-throughs ---
  // Write HD107SFrame<N>'s static color state, shared across all controllers of
  // the same N (one controller per image, so this is per-image in practice).

  /**
   * @brief Sets the global brightness applied to every packed pixel.
   * @param brightness Global brightness scale in [0, 255].
   */
  void setBrightness(uint8_t brightness) {
    HD107SFrame<N>::setBrightness(brightness);
  }

  /**
   * @brief Sets the white-point temperature gains applied per channel.
   * @param r Red temperature gain in [0, 255].
   * @param g Green temperature gain in [0, 255].
   * @param b Blue temperature gain in [0, 255].
   */
  void setTemperature(uint8_t r, uint8_t g, uint8_t b) {
    HD107SFrame<N>::setTemperature(r, g, b);
  }

  /**
   * @brief Sets the per-channel color-correction gains applied per pixel.
   * @param r Red correction gain in [0, 255].
   * @param g Green correction gain in [0, 255].
   * @param b Blue correction gain in [0, 255].
   */
  void setCorrection(uint8_t r, uint8_t g, uint8_t b) {
    HD107SFrame<N>::setCorrection(r, g, b);
  }

private:
  HD107SFrame<N> frames_[2]; /**< Double-buffered protocol frames (front/back). */
  TeensySPIDMA spi_; /**< Low-level async DMA+SPI driver for this strip. */
  /**
   * @brief Index (0/1) of the front buffer currently being DMA'd.
   * @details Plain int: every access is in the single column-ISR context; the
   *          completion ISR never touches it, so no barrier is needed.
   */
  int activeBuffer_;
  /**
   * @brief Monotonic count of frames successfully handed to the DMA engine.
   * @details Atomic (ISR RMW + cross-context read); relaxed — an independent
   *          counter, not a happens-before signal.
   */
  std::atomic<uint32_t> transferCount_;
  std::atomic<uint32_t> overrunCount_; /**< Monotonic count of frames dropped on overrun. */
};

#endif // ARDUINO
