/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host unit tests for the DMA LED controller's pure framing/decision math
 * (hardware/dma_led_core.h), which the Arduino-only TeensySPIDMA /
 * DMALEDController in dma_led.h derive their behavior from. Covers the
 * double-buffer toggle, the with_bg transfer-length select, and the
 * stale-transfer watchdog predicate at its boundaries including a micros()
 * unsigned-long rollover.
 */
#pragma once

#include "hardware/dma_led_core.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <limits>

namespace hs_test {
namespace dma_core {

// Compile-time proof the decision math folds (the device relies on these as
// constexpr/inline, off the column ISR hot path).
static_assert(dma::next_buffer(0) == 1);
static_assert(dma::next_buffer(1) == 0);
static_assert(dma::transfer_len(100, 200, false) == 100);
static_assert(dma::transfer_len(100, 200, true) == 200);
static_assert(!dma::transfer_stale(0, 0, 100));    // now == start
static_assert(!dma::transfer_stale(0, 99, 100));   // just below
static_assert(dma::transfer_stale(0, 100, 100));   // at bound
static_assert(dma::transfer_stale(0, 101, 100));   // above

/**
 * @brief Pin the double-buffer index toggle (0<->1).
 */
inline void test_next_buffer() {
  HS_EXPECT_EQ(dma::next_buffer(0), 1);
  HS_EXPECT_EQ(dma::next_buffer(1), 0);
}

/**
 * @brief Pin the transfer-length select for both with_bg values.
 */
inline void test_transfer_len() {
  const std::size_t base = 1444;
  const std::size_t composite = 2888;
  HS_EXPECT_EQ(dma::transfer_len(base, composite, false), base);
  HS_EXPECT_EQ(dma::transfer_len(base, composite, true), composite);
}

/**
 * @brief Pin the stale-transfer predicate at its watchdog boundaries.
 */
inline void test_transfer_stale_bounds() {
  const unsigned long wd = 100000UL;
  HS_EXPECT_FALSE(dma::transfer_stale(0, 0, wd));         // now == start
  HS_EXPECT_FALSE(dma::transfer_stale(5000, 5000, wd));
  HS_EXPECT_FALSE(dma::transfer_stale(0, wd - 1, wd));    // just below
  HS_EXPECT_TRUE(dma::transfer_stale(0, wd, wd));         // at bound
  HS_EXPECT_TRUE(dma::transfer_stale(0, wd + 1, wd));     // above
}

/**
 * @brief The predicate stays correct across an unsigned-long micros() rollover.
 * @details now_us < start_us when micros() has wrapped; the unsigned subtraction
 * still yields the true elapsed delta.
 */
inline void test_transfer_stale_wraparound() {
  const unsigned long wd = 100000UL;
  const unsigned long max = std::numeric_limits<unsigned long>::max();
  // start just before rollover, now just after: elapsed = 1 + (now+1), below wd.
  HS_EXPECT_FALSE(dma::transfer_stale(max - 10, 9, wd)); // elapsed 20
  // start before rollover, now far enough past it to exceed the watchdog.
  HS_EXPECT_TRUE(dma::transfer_stale(max - 10, wd, wd)); // elapsed wd + 11
  // exact boundary across the seam: elapsed == wd.
  HS_EXPECT_TRUE(dma::transfer_stale(max - 9, wd - 10, wd));
}

/**
 * @brief Run the DMA-core decision-math suite.
 * @return Number of failed expectations (0 on full pass).
 */
inline int run_dma_core_tests() {
  hs_test::ModuleFixture fixture("dma_core");

  test_next_buffer();
  test_transfer_len();
  test_transfer_stale_bounds();
  test_transfer_stale_wraparound();

  return fixture.result();
}

} // namespace dma_core
} // namespace hs_test
