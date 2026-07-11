/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/render/led.h — the correction-guard RAII happy path:
 * a normal scope exit must balance correction_guard_depth() back to 0. The
 * double-construct trap is covered separately by the death module.
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
 * @brief Runs every led test case in the module.
 * @return The module's failure count, as reported by end_module.
 */
inline int run_led_tests() {
  hs_test::ModuleFixture fixture("led");

  test_color_guard_balances();
  test_temp_guard_balances();
  test_sequential_guards_balance();

  return fixture.result();
}

} // namespace led_tests
} // namespace hs_test
