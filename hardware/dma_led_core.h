/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file dma_led_core.h
 * @brief Pure double-buffer / transfer-length / stale-transfer math for the DMA
 *        LED controller. No Teensy peripherals — split out of dma_led.h so the
 *        framing and watchdog decisions are host-unit-testable without a Teensy
 *        (see tests/test_dma_core.h). The Arduino-only TeensySPIDMA /
 *        DMALEDController in dma_led.h derive their behavior from these
 *        functions, so the host tests cover the real arithmetic.
 */

#include <cstddef>

namespace dma {

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
 *        (HD107SFrame::sizeWithBg()).
 * @param with_bg If true, transmit the composite buffer; otherwise the image
 *        frame only.
 * @return composite_size when with_bg, else base_size.
 */
constexpr std::size_t transfer_len(std::size_t base_size,
                                   std::size_t composite_size, bool with_bg) {
  return with_bg ? composite_size : base_size;
}

/**
 * @brief Stale-transfer watchdog predicate.
 * @param start_us micros() timestamp when the in-flight transfer was enabled.
 * @param now_us Current micros() timestamp.
 * @param watchdog_us Watchdog bound in µs.
 * @return true once now_us - start_us reaches watchdog_us.
 * @details Uses unsigned wrap-safe subtraction (now_us - start_us), matching the
 *          device's `micros() - transferStartUs_`: the elapsed delta stays
 *          correct across an unsigned-long micros() rollover.
 */
constexpr bool transfer_stale(unsigned long start_us, unsigned long now_us,
                              unsigned long watchdog_us) {
  return now_us - start_us >= watchdog_us;
}

} // namespace dma
