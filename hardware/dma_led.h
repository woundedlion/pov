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
#include "color.h"

// HD107SFrame (protocol buffer + color correction) lives in its own header so
// its wire-format / correction math is host-unit-testable without the Teensy
// peripherals below. See tests/test_hd107s_frame.h.
#include "hd107s_frame.h"

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
   * Call exactly once. A second call would re-run SPI.begin() and leak the
   * already-open SPI transaction, so re-initialization is a programming error
   * that fail-fast traps rather than silently corrupting the peripheral state.
   */
  void init() {
    HS_CHECK(!initialized_, "TeensySPIDMA::init() called twice");
    initialized_ = true;
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
    if (transferComplete_.load(std::memory_order_relaxed)) return;
    // Watchdog: a healthy frame clocks out in well under a millisecond at
    // 12 MHz SPI (the largest strips are still single-digit ms). An unbounded
    // spin means the DMA completion ISR never fired — a wedged channel — which
    // a fail-fast headless device must surface, not freeze on. The bound is
    // ~2 orders of magnitude above the worst real transfer. On timeout, name
    // the site on Serial (the trap is a bare illegal instruction with no
    // message) then trap; a truncated Serial write is harmless since we halt
    // immediately, and that holds even in the column-ISR context submitFrame()
    // runs in (in normal operation submitFrame() overrun-drops instead of
    // entering this spin — see submitFrame()).
    const unsigned long wait_start = micros();
    while (!transferComplete_.load(std::memory_order_relaxed)) {
      if (micros() - wait_start >= kTransferWatchdogUs) {
        hs::log("FATAL: DMA transfer watchdog timeout — completion ISR never fired");
        __builtin_trap();
      }
    }
  }

private:
  // Watchdog bound for waitComplete() (µs). A full-frame DMA is single-digit
  // ms even for the largest strips at 12 MHz; 100 ms is far above any real
  // transfer, so only a wedged DMA channel trips it.
  static constexpr unsigned long kTransferWatchdogUs = 100000UL;

  static void FASTRUN dmaISR() {
    if (instance_) {
      instance_->dma_.clearInterrupt();
      instance_->transferComplete_.store(true, std::memory_order_relaxed);
    }
  }

  DMAChannel dma_;
  // Completion flag handed between the DMA-completion ISR (writes true) and the
  // main/column thread (reads, and writes false to start a transfer). relaxed
  // ordering is intentional and correct on this single-core Cortex-M7: the ISR
  // and the thread are the SAME observer (the ISR preempts), so there is no
  // multi-core reordering for acquire/release to constrain — same rationale as
  // the instance_ note below. Crucially, this flag does NOT order the buffer
  // for the DMA engine: that buffer→DMA coherence is provided by
  // arm_dcache_flush_delete() (cache flush + DSB) before dma_.enable(), which
  // the atomic's memory_order has no effect on. Do not "upgrade" to
  // acquire/release expecting a visibility fix — it only adds DMB cost here.
  std::atomic<bool> transferComplete_;
  SPISettings spiSettings_;

  // Single-init guard: init() touches global SPI/DMA peripheral state that is
  // not safe to re-run. Set once on the first (setup-time) init() call.
  bool initialized_ = false;

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

inline TeensySPIDMA* TeensySPIDMA::instance_ = nullptr;

// ============================================================================
// DMALEDController — Double-buffered async LED controller
// ============================================================================

/**
 * @brief High-level double-buffered DMA LED controller for HD107S strips.
 * @tparam N Number of pixels.
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
  DMALEDController() : activeBuffer_(0), transferCount_(0), overrunCount_(0) {}

  /**
   * @brief One-time hardware initialization. Call from setup().
   */
  void begin() {
    spi_.init();
  }

  /**
   * @brief Returns the back frame (not currently being DMA'd).
   * Pack pixels via packPixel(), then call submitFrame().
   */
  HD107SFrame<N>& backFrame() {
    return frames_[1 - activeBuffer_];
  }

  /**
   * @brief Flushes the back frame and triggers async DMA transfer.
   * Call after all packPixel() calls on backFrame().
   * @param withBg If true, DMAs the composite buffer (image + trailing
   *               black frame) in a single transfer — zero gap, no spin.
   */
  void submitFrame(bool withBg = false) {
    if (!spi_.isComplete()) {
      // Drop on overrun — never enter transmitAsync()'s waitComplete() spin from
      // the column ISR. That spin can deadlock: transferComplete_ is only set by
      // the DMA-completion ISR, which cannot preempt an equal/lower-priority ISR,
      // so at best it extends the ISR past the next column tick. A dropped frame
      // is a transient (platform.h), not an invariant violation: the in-flight
      // DMA keeps the previous column; the already-packed back buffer is simply
      // discarded and the next column repacks.
      overrunCount_.fetch_add(1, std::memory_order_relaxed);
      return;
    }
    int back = 1 - activeBuffer_;
    frames_[back].flush();
    size_t len = withBg ? frames_[back].sizeWithBg() : frames_[back].size();
    spi_.transmitAsync(frames_[back].data(), len);
    activeBuffer_ = back;
    transferCount_.fetch_add(1, std::memory_order_relaxed);
  }

  // --- Diagnostics ---
  // These counters are this driver's telemetry surface. They are intentionally
  // NOT mirrored into a read-only registered param (the MindSplatter "Particles"
  // pattern): the param system's only consumer is the daydream WASM GUI, which
  // is absent from this ARDUINO-only build, and Phantasm.ino has no param
  // reader. A caller that wants them (e.g. a future serial telemetry path) reads
  // them here.
  uint32_t getTransferCount() const {
    return transferCount_.load(std::memory_order_relaxed);
  }
  uint32_t getOverrunCount() const {
    return overrunCount_.load(std::memory_order_relaxed);
  }

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
  // Plain int, no atomic/volatile needed: every access lives in the single
  // column-ISR context — the reads in backFrame()/submitFrame() and the write in
  // submitFrame() all run from show_col(). The DMA-completion ISR
  // touches only transferComplete_, never this index, so there is no
  // cross-context reader to synchronize and a barrier here would be pure cost
  // (same single-observer model as the transferComplete_/instance_ notes). The
  // frames_[] buffer it selects is made coherent for the DMA engine by the
  // cache flush in transmitAsync() (see above), not by this index.
  int activeBuffer_;
  // Atomic, not volatile: incremented in the column-ISR context (show/
  // submitFrame) and read elsewhere; volatile makes neither the RMW nor the
  // cross-context read well-defined. Relaxed ordering suffices — these are
  // independent monotonic counters, not a happens-before signal.
  std::atomic<uint32_t> transferCount_;
  std::atomic<uint32_t> overrunCount_;
};

#endif // ARDUINO
