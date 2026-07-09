/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file dma_led_controller.h
 * @brief Double-buffered high-level DMA LED controller, templated on its SPI/DMA
 *        transport.
 *
 * Split out of dma_led.h so the double-buffer / overrun-drop / watchdog
 * orchestration is host-unit-testable against a mock transport, without the
 * Teensy peripherals (see tests/test_dma_controller.h). On device the transport
 * defaults to TeensySPIDMA; a host test substitutes a recording mock.
 *
 * The framing/decision math this drives lives in dma_led_core.h and the wire
 * format in hd107s_frame.h, both already host-tested.
 *
 * On ARDUINO the default transport TeensySPIDMA must be visible here, so this
 * header is included from dma_led.h *after* that class is defined; host builds
 * (ARDUINO undefined) carry no default and name the transport explicitly.
 */

#include <atomic>
#include <cstddef>
#include <cstdint>

#include "dma_led_core.h"
#include "hd107s_frame.h"

/**
 * @brief High-level double-buffered DMA LED controller for HD107S strips.
 * @tparam N Number of pixels.
 * @tparam Transport SPI/DMA transport driving the wire; defaults to TeensySPIDMA
 *         on device. Must provide the transmit/completion/watchdog contract
 *         TeensySPIDMA exposes: Transport(uint32_t clock), init(), isComplete(),
 *         checkStaleTransfer(), and transmitAsync(const uint8_t*, size_t).
 * @note One instance per firmware image: it drives the singleton TeensySPIDMA
 *       backing the shared DMA-completion ISR, so a second begin() traps.
 *
 * Typical ISR usage (per column):
 *   auto& f = controller.backFrame();  // back buffer (not being DMA'd)
 *   // ... pack pixels into f via packPixel() ...
 *   controller.submitFrame();          // triggers async DMA, returns immediately
 *   // ISR exits → DMA transfers in background → main loop gets more CPU
 */
template <int N, class Transport
#ifdef ARDUINO
          = TeensySPIDMA
#endif
          >
class DMALEDController {
public:
  /**
   * @brief Constructs the controller, optionally overriding the SPI clock.
   * @param clock SPI clock in Hz forwarded to the transport (default: 12 MHz).
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
  Transport spi_; /**< Low-level async DMA+SPI transport for this strip. */
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
