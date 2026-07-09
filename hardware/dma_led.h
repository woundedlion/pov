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
    // ISR, where spinning would deadlock — the DMA-completion ISR that marks
    // transferComplete_ true cannot preempt an equal/lower-priority ISR.
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
    if (dma::transfer_stale(transferStartUs_, micros(), TRANSFER_WATCHDOG_US)) {
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
  static constexpr unsigned long TRANSFER_WATCHDOG_US = 100000UL;

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

// DMALEDController (double-buffered controller) lives in its own header so its
// overrun-drop / double-buffer / watchdog orchestration is host-testable against
// a mock transport. Included here, after TeensySPIDMA, so its default transport
// resolves. See tests/test_dma_controller.h.
#include "dma_led_controller.h"

#endif // ARDUINO
