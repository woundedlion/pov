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
#include "core/platform.h" // HS_CHECK / hs::log used below; explicit, not via color.h
#include "core/color.h"

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
   *
   * This is the only place a wedged channel is surfaced. Its caller,
   * DMALEDController::submitFrame(), guards with isComplete() and *drops* on
   * overrun rather than spinning (never entering transmitAsync while a transfer
   * is in flight — transmitAsync traps if it ever is). So a channel whose
   * completion ISR never fires would otherwise be masked forever — isComplete()
   * stays false, every submitFrame() drops and bumps overrunCount, and the dark
   * strip the fail-fast doctrine exists to surface never traps. Called from
   * submitFrame's overrun branch, this traps once the in-flight transfer has
   * outlived the watchdog bound (healthy transfers finish in single-digit ms;
   * 100 ms means wedged). Cheap enough for the drop path: one micros() read and
   * a compare, and the drop path runs at most once per column. No-op while a
   * transfer is legitimately still completing.
   *
   * Coverage boundary, by design: this watchdog fires only while submitFrame()
   * is being called, i.e. while columns are being submitted. POVSegmented's
   * fail-dark path latches the strip (dark_latched_) after one accepted black
   * frame and then stops submitting for the whole dark window, so a channel that
   * wedges during a sustained dark window is NOT surfaced until submission
   * resumes. That gap is intentional and acceptable: the wedge symptom — a dark
   * strip — is identical to the intended dark output, so there is nothing to
   * fail-fast about while dark-latched. The wedge is surfaced the moment real
   * column work resumes (the next submitFrame() on overrun). Do not add an idle-
   * path poll to close this gap without weighing the per-wake cost on the
   * column ISR — the diagnostic value is nil while the output is meant to be dark.
   */
  void checkStaleTransfer() {
    if (transferComplete_.load(std::memory_order_relaxed)) return;
    if (micros() - transferStartUs_ >= kTransferWatchdogUs) {
      hs::log("FATAL: DMA channel wedged — in-flight transfer outlived the "
              "watchdog on the overrun-drop path; completion ISR never fired");
      __builtin_trap();
    }
  }

private:
  /**
   * @brief Watchdog bound for checkStaleTransfer(), in µs.
   * @details A full-frame DMA is single-digit ms even for the largest strips at
   *          12 MHz; 100 ms is far above any real transfer, so only a wedged DMA
   *          channel trips it.
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
   * @details The ISR writes true; the main/column thread reads, and writes
   *          false to start a transfer. relaxed ordering is intentional and
   *          correct on this single-core Cortex-M7: the ISR and the thread are
   *          the SAME observer (the ISR preempts), so there is no multi-core
   *          reordering for acquire/release to constrain — same rationale as
   *          the instance_ note below. Crucially, this flag does NOT order the
   *          buffer for the DMA engine: that buffer→DMA coherence is provided by
   *          arm_dcache_flush() (cache clean + DSB) before dma_.enable(),
   *          which the atomic's memory_order has no effect on. Do not "upgrade"
   *          to acquire/release expecting a visibility fix — it only adds DMB
   *          cost here.
   */
  std::atomic<bool> transferComplete_;
  /**
   * @brief micros() at which the in-flight transfer was enabled.
   * @details Written and read only in the column-ISR context (transmitAsync /
   *          checkStaleTransfer, both reached from submitFrame) — the
   *          DMA-completion ISR never touches it — so a plain scalar is correct,
   *          same single-observer model as DMALEDController's activeBuffer_.
   *          Only meaningful while transferComplete_ is false.
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
    return frames_[1 - activeBuffer_];
  }

  /**
   * @brief Flushes the back frame and triggers async DMA transfer.
   * Call after all packPixel() calls on backFrame().
   * @param withBg If true, DMAs the composite buffer (image + trailing
   *               black frame) in a single transfer — zero gap, no spin.
   * @return true if the frame was handed to the DMA engine; false if it was
   *         dropped on overrun (prior transfer still in flight). Callers that
   *         must KNOW a specific frame landed — e.g. the fail-dark latch, which
   *         may not hold a stale bright column — gate on this; the steady-state
   *         column path ignores it (a dropped image column self-heals next tick).
   */
  [[nodiscard]] bool submitFrame(bool withBg = false) {
    if (!spi_.isComplete()) {
      // Drop on overrun. A transfer that NEVER completes is a wedged channel,
      // not a transient, so surface it here — the drop path is where it shows.
      spi_.checkStaleTransfer();
      overrunCount_.fetch_add(1, std::memory_order_relaxed);
      return false;
    }
    int back = 1 - activeBuffer_;
    frames_[back].flush();
    size_t len = withBg ? frames_[back].sizeWithBg() : frames_[back].size();
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
   * @details Plain int, no atomic/volatile needed: every access lives in the
   *          single column-ISR context — the reads in backFrame()/submitFrame()
   *          and the write in submitFrame() all run from show_col(). The
   *          DMA-completion ISR touches only transferComplete_, never this
   *          index, so there is no cross-context reader to synchronize and a
   *          barrier here would be pure cost (same single-observer model as the
   *          transferComplete_/instance_ notes). The frames_[] buffer it selects
   *          is made coherent for the DMA engine by the cache flush in
   *          transmitAsync() (see above), not by this index.
   */
  int activeBuffer_;
  /**
   * @brief Monotonic count of frames successfully handed to the DMA engine.
   * @details Atomic, not volatile: incremented in the column-ISR context (show/
   *          submitFrame) and read elsewhere; volatile makes neither the RMW nor
   *          the cross-context read well-defined. Relaxed ordering suffices —
   *          these are independent monotonic counters, not a happens-before
   *          signal.
   */
  std::atomic<uint32_t> transferCount_;
  std::atomic<uint32_t> overrunCount_; /**< Monotonic count of frames dropped on overrun. */
};

#endif // ARDUINO
