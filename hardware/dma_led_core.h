/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file dma_led_core.h
 * @brief Pure double-buffer / transfer-length / transfer-duration /
 *        stale-transfer math for the DMA LED controller. No Teensy
 *        peripherals — split out of dma_led.h so the framing and watchdog
 *        decisions are host-unit-testable without a Teensy (see
 *        tests/test_dma_core.h). The Arduino-only TeensySPIDMA /
 *        DMALEDController in dma_led.h derive their behavior from these
 *        functions, so the host tests cover the real arithmetic.
 */

#include <cstddef>
#include <cstdint>

namespace dma {

inline constexpr unsigned long TRANSFER_WATCHDOG_US = 5000UL;
inline constexpr uint32_t DEFAULT_CLOCK_HZ = 12000000;
inline constexpr uint32_t SEGMENTED_CLOCK_HZ = 24000000;
inline constexpr uint32_t LPSPI_FUNCTIONAL_CLOCK_HZ = 240000000;

/**
 * @brief Toggles the double-buffer index between 0 and 1.
 * @param active Current front-buffer index, 0 or 1.
 * @return The other buffer index (1 - active).
 */
constexpr int next_buffer(int active) { return 1 - active; }

/**
 * @brief Selects the DMA transfer length for a frame.
 * @param base_size Image-frame size in bytes (HD107SFrame::size()).
 * @param composite_size Composite size including the trailing black frame
 *        (HD107SFrame::size_with_bg()).
 * @param with_bg If true, transmit the composite buffer; otherwise the image
 *        frame only.
 * @return composite_size when with_bg, else base_size.
 */
constexpr std::size_t transfer_len(std::size_t base_size,
                                   std::size_t composite_size, bool with_bg) {
  return with_bg ? composite_size : base_size;
}

/**
 * @brief Worst-case duration of one column's LED transfer, in µs.
 * @param bytes Bytes clocked out for the column (image plus any black strobe).
 * @param clock_hz Bit clock the transport runs at, in Hz.
 * @return Transfer duration including LPSPI byte framing, rounded up to µs.
 * @pre clock_hz > 0.
 * @details Rounded up so the driver's `column_interval_us > transfer_us` check
 *          never under-counts the transfer and admits a configuration that
 *          overruns the DMA every column.
 */
constexpr unsigned long transfer_us(unsigned long bytes,
                                    unsigned long clock_hz) {
  const uint64_t REQUESTED_DIVIDER =
      (LPSPI_FUNCTIONAL_CLOCK_HZ + static_cast<uint64_t>(clock_hz) - 1) /
      clock_hz;
  const uint64_t DIVIDER = REQUESTED_DIVIDER < 2 ? 2 : REQUESTED_DIVIDER;
  // Teensy SPI sets DBT and PCSSCK to (SCKDIV / 2), with SCKDIV = divider - 2.
  const uint64_t BYTE_CLOCKS = 8 * DIVIDER + 2 * ((DIVIDER - 2) / 2) + 4;
  return static_cast<unsigned long>(
      (static_cast<uint64_t>(bytes) * BYTE_CLOCKS * 1000000u +
       LPSPI_FUNCTIONAL_CLOCK_HZ - 1) /
      LPSPI_FUNCTIONAL_CLOCK_HZ);
}

/**
 * @brief Stale-transfer watchdog predicate.
 * @param start_us micros() timestamp when the in-flight transfer was enabled.
 * @param now_us Current micros() timestamp.
 * @param watchdog_us Watchdog bound in µs.
 * @return true once now_us - start_us reaches watchdog_us.
 * @details Uses unsigned wrap-safe subtraction (now_us - start_us), matching the
 *          device's `micros() - transfer_start_us`: the elapsed delta stays
 *          correct across an unsigned-long micros() rollover.
 */
constexpr bool transfer_stale(unsigned long start_us, unsigned long now_us,
                              unsigned long watchdog_us) {
  return now_us - start_us >= watchdog_us;
}

} // namespace dma
