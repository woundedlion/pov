/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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
#include <concepts>
#include <cstddef>
#include <cstdint>

#include "dma_led_core.h"
#include "hd107s_frame.h"

/**
 * @brief The transmit/completion/watchdog contract a controller transport must
 *        satisfy.
 * @tparam T Candidate transport type.
 * @details Mirrors what TeensySPIDMA exposes. Asserted on DMALEDController's
 *          template parameter, so a non-conforming transport fails at the
 *          instantiation site rather than inside submitFrame().
 */
template <class T>
concept LedTransport = std::constructible_from<T, uint32_t> &&
                       requires(T t, const uint8_t *data, size_t len) {
                         t.init();
                         { t.isComplete() } -> std::convertible_to<bool>;
                         t.checkStaleTransfer();
                         t.transmitAsync(data, len);
                       };

/**
 * @brief High-level double-buffered DMA LED controller for HD107S strips.
 * @tparam N Number of pixels.
 * @tparam Transport SPI/DMA transport driving the wire, satisfying
 *         LedTransport; defaults to TeensySPIDMA on device.
 * @note One instance per firmware image: it drives the singleton TeensySPIDMA
 *       backing the shared DMA-completion ISR, so a second begin() traps.
 * @note A static instance belongs in DMAMEM: the HD107SFrame buffers are the
 *       eDMA TX source, and only in cached OCRAM is their arm_dcache_flush()
 *       write-back meaningful (in DTCM it is a no-op). GCC silently drops the
 *       DMAMEM section attribute on a vague-linkage template static member, so
 *       such an instance must be defined as an explicit specialization, whose
 *       ordinary strong linkage keeps the attribute.
 *
 * Typical ISR usage (per column):
 *   auto& f = controller.backFrame();  // back buffer (not being DMA'd)
 *   // ... pack pixels into f via packPixel() ...
 *   controller.submitFrame();          // triggers async DMA, returns immediately
 *   // ISR exits → DMA transfers in background → main loop gets more CPU
 */
template <int N, LedTransport Transport
#ifdef ARDUINO
                 = TeensySPIDMA
#endif
          >
class DMALEDController {
public:
  /** @brief SPI clock a default-constructed controller runs at, in Hz. */
  static constexpr uint32_t DEFAULT_CLOCK_HZ = 12000000;

  /**
   * @brief Constructs the controller, optionally overriding the SPI clock.
   * @param clock SPI clock in Hz forwarded to the transport.
   *              The Phantasm driver passes 24 MHz (see pov_segmented.h).
   */
  explicit DMALEDController(uint32_t clock = DEFAULT_CLOCK_HZ)
      : spi(clock), activeBuffer(0), transferCount(0), overrunCount(0) {}

  /**
   * @brief One-time hardware initialization. Call from setup().
   */
  void begin() { spi.init(); }

  /**
   * @brief Returns the back frame (not currently being DMA'd).
   * @return Reference to the back-buffer frame; pack pixels via packPixel(),
   *         then call submitFrame().
   */
  HD107SFrame<N> &backFrame() { return frames[dma::next_buffer(activeBuffer)]; }

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
    if (!spi.isComplete()) {
      // Drop on overrun. A transfer that NEVER completes is a wedged channel,
      // not a transient, so surface it here — the drop path is where it shows.
      spi.checkStaleTransfer();
      overrunCount.fetch_add(1, std::memory_order_relaxed);
      return false;
    }
    int back = dma::next_buffer(activeBuffer);
    frames[back].flush();
    size_t len = dma::transfer_len(frames[back].size(),
                                   frames[back].sizeWithBg(), withBg);
    spi.transmitAsync(frames[back].data(), len);
    activeBuffer = back;
    transferCount.fetch_add(1, std::memory_order_relaxed);
    return true;
  }

  // --- Diagnostics ---
  /**
   * @brief Returns the count of frames handed to the DMA engine since start.
   * @return Monotonic transfer counter (number of successful submitFrame()s).
   */
  uint32_t getTransferCount() const {
    return transferCount.load(std::memory_order_relaxed);
  }
  /**
   * @brief Returns the count of frames dropped on overrun since start.
   * @return Monotonic overrun counter (frames dropped because a prior transfer
   *         was still in flight).
   */
  uint32_t getOverrunCount() const {
    return overrunCount.load(std::memory_order_relaxed);
  }

  // --- Configuration pass-throughs ---
  // Write HD107SFrame<N>'s static color state, shared across all controllers of
  // the same N (one controller per image, so this is per-image in practice).

  /**
   * @brief Sets the global brightness applied to every packed pixel.
   * @param brightness Global brightness scale in [0, 255].
   */
  static void setBrightness(uint8_t brightness) {
    HD107SFrame<N>::setBrightness(brightness);
  }

  /**
   * @brief Sets the white-point temperature gains applied per channel.
   * @param r Red temperature gain in [0, 255].
   * @param g Green temperature gain in [0, 255].
   * @param b Blue temperature gain in [0, 255].
   */
  static void setTemperature(uint8_t r, uint8_t g, uint8_t b) {
    HD107SFrame<N>::setTemperature(r, g, b);
  }

  /**
   * @brief Sets the per-channel color-correction gains applied per pixel.
   * @param r Red correction gain in [0, 255].
   * @param g Green correction gain in [0, 255].
   * @param b Blue correction gain in [0, 255].
   */
  static void setCorrection(uint8_t r, uint8_t g, uint8_t b) {
    HD107SFrame<N>::setCorrection(r, g, b);
  }

private:
  HD107SFrame<N>
      frames[2]; /**< Double-buffered protocol frames (front/back). */
  Transport spi; /**< Low-level async DMA+SPI transport for this strip. */
  /**
   * @brief Index (0/1) of the front buffer currently being DMA'd.
   * @details Plain int: every access is in the single column-ISR context; the
   *          completion ISR never touches it, so no barrier is needed.
   */
  int activeBuffer;
  /**
   * @brief Monotonic count of frames successfully handed to the DMA engine.
   * @details Atomic (ISR RMW + cross-context read); relaxed — an independent
   *          counter, not a happens-before signal.
   */
  std::atomic<uint32_t> transferCount;
  std::atomic<uint32_t>
      overrunCount; /**< Monotonic count of frames dropped on overrun. */
};
