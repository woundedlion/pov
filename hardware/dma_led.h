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
      : transfer_complete(true), spi_settings(clock, order, mode) {}

  /**
   * @brief Initializes SPI and DMA hardware. Must be called from setup(),
   *        not from global constructors (peripherals may not be ready yet).
   *
   * Call exactly once; a second call traps (re-running SPI.begin() would leak
   * the open transaction).
   */
  void init() {
    HS_CHECK(!initialized, "TeensySPIDMA::init() called twice");
    HS_CHECK(instance == nullptr,
             "TeensySPIDMA: only one instance supported (shared DMA-completion "
             "ISR singleton)");
    initialized = true;
    instance = this;

    SPI.begin();
    SPI.beginTransaction(spi_settings);

    LPSPI4_DER |= LPSPI_DER_TDDE; // TX DMA request
    // RXMSK discards received frames: nothing drains the RX FIFO, so without it
    // the FIFO fills after 16 bytes. With RX masked, NOSTALL stays clear so a TX
    // underrun pauses the clock rather than shifting out a stale byte — a pause
    // is invisible to the self-clocked HD107S protocol, a stale byte is not.
    LPSPI4_TCR = (LPSPI4_TCR & ~LPSPI_TCR_FRAMESZ(31)) | LPSPI_TCR_FRAMESZ(7) |
                 LPSPI_TCR_RXMSK; // 8-bit frames, TX-only
    LPSPI4_CFGR1 &= ~LPSPI_CFGR1_NOSTALL;

    dma_channel.begin(true);
    // The uint8_t cast is load-bearing: the 32-bit destination() overload sets
    // ATTR_DST to 32-bit transfers, conflicting with sourceBuffer()'s NBYTES=1
    // and raising an eDMA config error on enable — nothing transmits.
    dma_channel.destination(reinterpret_cast<volatile uint8_t &>(LPSPI4_TDR));
    dma_channel.triggerAtHardwareEvent(DMAMUX_SOURCE_LPSPI4_TX);
    dma_channel.disableOnCompletion();
    dma_channel.interruptAtCompletion();
    dma_channel.attachInterrupt(dma_isr);
  }

  /**
   * @brief Starts an async DMA transfer. Returns immediately.
   * @param data Pointer to byte buffer (must remain valid until complete).
   * @param len  Number of bytes to transmit.
   * @pre The caller has cleaned the buffer from cache (submitFrame() does this
   *      via frames[back].flush()); this method only enables the DMA and does
   *      not flush.
   */
  void transmitAsync(const uint8_t *data, size_t len) {
    // Trap rather than spin on an in-flight transfer: this runs in the column
    // ISR, where spinning would deadlock — the DMA-completion ISR that marks
    // transfer_complete true cannot preempt an equal/lower-priority ISR.
    if (!transfer_complete.load(std::memory_order_relaxed)) {
      hs::log("FATAL: transmitAsync entered with a transfer still in flight — "
              "submitFrame() must guard with isComplete()");
      __builtin_trap();
    }
    transfer_complete.store(false, std::memory_order_relaxed);
    transfer_start_us = micros();
    dma_channel.sourceBuffer(data, len);
    dma_channel.enable();
  }

  /**
   * @brief Reports whether the in-flight transfer has finished.
   * @return true once the in-flight transfer's completion ISR has fired.
   */
  bool isComplete() const {
    return transfer_complete.load(std::memory_order_relaxed);
  }

  /**
   * @brief Surfaces a permanently wedged DMA channel from the overrun-drop path.
   * @details submitFrame() drops on overrun rather than spinning, so a channel
   *          whose completion ISR never fires would otherwise stay masked
   *          forever. Traps once the in-flight transfer outlives the watchdog.
   *          Only fires while submitFrame() is being called.
   */
  void checkStaleTransfer() {
    if (transfer_complete.load(std::memory_order_relaxed))
      return;
    if (dma::transfer_stale(transfer_start_us, micros(),
                            TRANSFER_WATCHDOG_US)) {
      hs::log("FATAL: DMA channel wedged — in-flight transfer outlived the "
              "watchdog on the overrun-drop path; completion ISR never fired");
      __builtin_trap();
    }
  }

private:
  /**
   * @brief Watchdog bound for checkStaleTransfer(), in µs.
   * @details Covers both shipping configurations. Holosphere: 40 px → 336-byte
   *          composite at the 12 MHz default clock = 224 µs, column period
   *          1302 µs. Phantasm: 72 px → 600-byte composite at 24 MHz = 200 µs,
   *          column period 434 µs. 5 ms is over 20× either transfer, so only a
   *          wedged channel trips it, and under 12 column periods on either, so
   *          a wedge does not blank the strip for hundreds of columns before
   *          trapping.
   */
  static constexpr unsigned long TRANSFER_WATCHDOG_US = 5000UL;

  /**
   * @brief DMA completion ISR: clears the interrupt and marks the transfer done.
   * @details Dispatched via the singleton instance; runs in interrupt context.
   */
  static void FASTRUN dma_isr() {
    if (instance) {
      instance->dma_channel.clearInterrupt();
      instance->transfer_complete.store(true, std::memory_order_relaxed);
      // DMA_CINT is a posted write; without this the ISR can return before it
      // retires and the NVIC re-pends the interrupt it just cleared.
      asm volatile("dsb" ::: "memory");
    }
  }

  DMAChannel
      dma_channel; /**< eDMA channel wired to LPSPI4 TX for async transmission. */
  /**
   * @brief Completion flag handed between the DMA-completion ISR and the
   *        main/column thread.
   * @details Single-observer (ISR preempts the thread), so relaxed ordering
   *          suffices. Does NOT order the buffer for the DMA engine — the
   *          caller's arm_dcache_flush() before transmitAsync() does that.
   */
  std::atomic<bool> transfer_complete;
  /**
   * @brief micros() at which the in-flight transfer was enabled.
   * @details Touched only in column-ISR context (transmitAsync /
   *          checkStaleTransfer), never the completion ISR, so a plain scalar is
   *          correct. Only meaningful while transfer_complete is false.
   */
  unsigned long transfer_start_us = 0;
  SPISettings
      spi_settings; /**< Cached SPI clock/bit-order/mode for this driver. */

  /**
   * @brief Single-init guard set once on the first (setup-time) init() call.
   * @details init() touches global SPI/DMA peripheral state that is not safe to
   *          re-run.
   */
  bool initialized = false;

  /**
   * @brief Singleton pointer for ISR callback dispatch. Exactly one
   *        TeensySPIDMA per image; init() traps on a second instance.
   * @details Written once from setup() before any ISR uses it, then read from
   *          the completion ISR — race-free under the single-observer model (see
   *          transfer_complete).
   */
  static TeensySPIDMA *instance;
};

inline TeensySPIDMA *TeensySPIDMA::instance = nullptr;

// DMALEDController (double-buffered controller) lives in its own header so its
// overrun-drop / double-buffer / watchdog orchestration is host-testable against
// a mock transport. Included here, after TeensySPIDMA, so its default transport
// resolves. See tests/test_dma_controller.h.
#include "dma_led_controller.h"

#endif // ARDUINO
