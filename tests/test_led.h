/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/led.h — the correction-guard RAII happy path: a
 * normal scope exit must balance correction_guard_depth() back to 0 and leave
 * the FastLED selectors at the canonical baseline. The double-construct trap is
 * covered separately by the death module.
 *
 * Self-contained header — no external framework. run_led_tests() returns the
 * module failure count.
 */
#pragma once

#include "core/render/led.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace led_tests {

/**
 * @brief Verifies a NoColorCorrection guard raises the depth to 1 in scope and
 *        restores it to 0 on scope exit.
 */
inline void test_color_guard_balances() {
  HS_EXPECT_EQ(correction_guard_depth(), 0);
  {
    NoColorCorrection guard;
    HS_EXPECT_EQ(correction_guard_depth(), 1);
  }
  HS_EXPECT_EQ(correction_guard_depth(), 0);
}

/**
 * @brief Verifies a NoTempCorrection guard raises the depth to 1 in scope and
 *        restores it to 0 on scope exit.
 */
inline void test_temp_guard_balances() {
  HS_EXPECT_EQ(correction_guard_depth(), 0);
  {
    NoTempCorrection guard;
    HS_EXPECT_EQ(correction_guard_depth(), 1);
  }
  HS_EXPECT_EQ(correction_guard_depth(), 0);
}

/**
 * @brief Verifies sequential (non-overlapping) guards of either type each
 *        balance back to 0, so the shared counter never leaks across scopes.
 */
inline void test_sequential_guards_balance() {
  {
    NoColorCorrection a;
    HS_EXPECT_EQ(correction_guard_depth(), 1);
  }
  HS_EXPECT_EQ(correction_guard_depth(), 0);
  {
    NoTempCorrection b;
    HS_EXPECT_EQ(correction_guard_depth(), 1);
  }
  HS_EXPECT_EQ(correction_guard_depth(), 0);
  {
    NoColorCorrection c;
    HS_EXPECT_EQ(correction_guard_depth(), 1);
  }
  HS_EXPECT_EQ(correction_guard_depth(), 0);
}

/**
 * @brief Verifies NoColorCorrection disables both corrections for its scope and
 *        reinstates the TypicalLEDStrip/Candle baseline on scope exit.
 */
inline void test_color_guard_sets_selectors() {
  FastLED.setCorrection(TypicalLEDStrip);
  FastLED.setTemperature(Candle);
  {
    NoColorCorrection guard;
    HS_EXPECT_EQ(FastLED.last_correction, UncorrectedColor);
    HS_EXPECT_EQ(FastLED.last_temperature, UncorrectedTemperature);
  }
  HS_EXPECT_EQ(FastLED.last_correction, TypicalLEDStrip);
  HS_EXPECT_EQ(FastLED.last_temperature, Candle);
}

/**
 * @brief Verifies NoTempCorrection keeps TypicalLEDStrip color correction while
 *        disabling temperature, and reinstates the baseline on scope exit.
 */
inline void test_temp_guard_sets_selectors() {
  FastLED.setCorrection(UncorrectedColor);
  FastLED.setTemperature(Candle);
  {
    NoTempCorrection guard;
    HS_EXPECT_EQ(FastLED.last_correction, TypicalLEDStrip);
    HS_EXPECT_EQ(FastLED.last_temperature, UncorrectedTemperature);
  }
  HS_EXPECT_EQ(FastLED.last_correction, TypicalLEDStrip);
  HS_EXPECT_EQ(FastLED.last_temperature, Candle);
}

/**
 * @brief Verifies both destructors restore the canonical baseline rather than
 *        the correction that was active at construction.
 */
inline void test_guards_restore_baseline_not_prior_state() {
  FastLED.setCorrection(UncorrectedColor);
  FastLED.setTemperature(UncorrectedTemperature);
  { NoColorCorrection guard; }
  HS_EXPECT_EQ(FastLED.last_correction, TypicalLEDStrip);
  HS_EXPECT_EQ(FastLED.last_temperature, Candle);

  FastLED.setCorrection(UncorrectedColor);
  FastLED.setTemperature(UncorrectedTemperature);
  { NoTempCorrection guard; }
  HS_EXPECT_EQ(FastLED.last_correction, TypicalLEDStrip);
  HS_EXPECT_EQ(FastLED.last_temperature, Candle);
}

/**
 * @brief Runs every led test case in the module.
 * @return The module's failure count, as reported by end_module.
 */
inline int run_led_tests() {
  hs_test::ModuleFixture fixture("led");

  test_color_guard_balances();
  test_temp_guard_balances();
  test_sequential_guards_balance();
  test_color_guard_sets_selectors();
  test_temp_guard_sets_selectors();
  test_guards_restore_baseline_not_prior_state();

  return fixture.result();
}

} // namespace led_tests
} // namespace hs_test
