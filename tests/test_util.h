/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/util.h — the circular-domain wrap/distance helpers and
 * the apply_if_changed live-slider idiom.
 *
 * Self-contained header — no external framework. run_util_tests() returns the
 * module failure count.
 */
#pragma once

#include "core/util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace util_tests {

// --- wrap(float, m) / wrap_t ------------------------------------------------

/**
 * @brief Verifies wrap(float, m) folds any real into [0, m) for an arbitrary base m.
 * @details Checks in-range identity, negative fold-up, above-m fold-down, multiple
 *          periods, a non-unit base, and that results stay in [0, m) across signs.
 */
inline void test_wrap_float() {
  HS_EXPECT_NEAR(wrap(0.25f, 1.0f), 0.25f, 1e-6f);
  HS_EXPECT_NEAR(wrap(-0.25f, 1.0f), 0.75f, 1e-6f);
  HS_EXPECT_NEAR(wrap(1.25f, 1.0f), 0.25f, 1e-6f);
  HS_EXPECT_NEAR(wrap(-3.25f, 1.0f), 0.75f, 1e-6f);
  HS_EXPECT_NEAR(wrap(7.0f, 5.0f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(wrap(-1.0f, 5.0f), 4.0f, 1e-6f);

  for (int i = -30; i <= 30; ++i) {
    float w = wrap(i * 0.37f, 2.0f);
    HS_EXPECT_TRUE(w >= 0.0f && w < 2.0f);
  }
}

/**
 * @brief Verifies wrap_t folds any real into the unit interval [0, 1) (the m==1 case).
 */
inline void test_wrap_t() {
  HS_EXPECT_NEAR(wrap_t(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(0.4f), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(1.4f), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(-0.25f), 0.75f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(-2.1f), 0.9f, 1e-5f);
  // Tiny negative t: t - floorf(t) rounds up to exactly 1.0f, which the guard
  // must fold back to 0 to keep the half-open [0, 1) contract.
  float tiny = wrap_t(-1e-8f);
  HS_EXPECT_TRUE(tiny >= 0.0f && tiny < 1.0f);
  HS_EXPECT_NEAR(tiny, 0.0f, 1e-6f);
  for (int i = -50; i <= 50; ++i) {
    float w = wrap_t(i * 0.13f);
    HS_EXPECT_TRUE(w >= 0.0f && w < 1.0f);
  }
}

// --- wrap(int, int) ---------------------------------------------------------

/**
 * @brief Verifies wrap(int, int) folds an integer into [0, m) using the integer overload.
 * @details The exact (int, int) match selects the integer overload, so no float math runs.
 */
inline void test_wrap_int() {
  HS_EXPECT_EQ(wrap(3, 5), 3);
  HS_EXPECT_EQ(wrap(5, 5), 0);
  HS_EXPECT_EQ(wrap(7, 5), 2);
  HS_EXPECT_EQ(wrap(-1, 5), 4);
  HS_EXPECT_EQ(wrap(-6, 5), 4);
  for (int i = -17; i <= 17; ++i) {
    int w = wrap(i, 6);
    HS_EXPECT_TRUE(w >= 0 && w < 6);
  }
}

// --- wrap(mixed int/float) --------------------------------------------------

/**
 * @brief Verifies a mixed (int, float) wrap call returns the common type (float).
 * @details The template yields std::common_type_t (float), so the fractional part
 *          survives rather than being truncated to an int; (float, int) also uses
 *          float math, which geometry.h relies on.
 */
inline void test_wrap_mixed_type() {
  static_assert(std::is_same_v<decltype(wrap(3, 2.5f)), float>,
                "wrap(int, float) returns the common type (float)");
  HS_EXPECT_NEAR(wrap(3, 2.5f), 0.5f, 1e-6f);  // fmod(3, 2.5) == 0.5
  HS_EXPECT_NEAR(wrap(7, 2.5f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(wrap(-1, 2.5f), 1.5f, 1e-6f);
  HS_EXPECT_NEAR(wrap(7.5f, 5), 2.5f, 1e-6f);
}

// --- fast_wrap --------------------------------------------------------------

/**
 * @brief Verifies fast_wrap folds x with a single add or subtract.
 * @details fast_wrap assumes x is at most one period out of range, valid only for
 *          x in [-W, 2W); the test covers in-range identity, one period high, one
 *          period low, and the whole legal window landing in [0, W).
 */
inline void test_fast_wrap() {
  constexpr int W = 8;
  for (int x = 0; x < W; ++x)
    HS_EXPECT_EQ(fast_wrap(x, W), x);
  HS_EXPECT_EQ(fast_wrap(W, W), 0);
  HS_EXPECT_EQ(fast_wrap(W + 3, W), 3);
  HS_EXPECT_EQ(fast_wrap(2 * W - 1, W), W - 1);
  HS_EXPECT_EQ(fast_wrap(-1, W), W - 1);
  HS_EXPECT_EQ(fast_wrap(-W, W), 0);
  for (int x = -W; x < 2 * W; ++x) {
    int w = fast_wrap(x, W);
    HS_EXPECT_TRUE(w >= 0 && w < W);
  }
}

// --- shortest_distance / fwd_distance ---------------------------------------

/**
 * @brief Verifies shortest_distance returns the smaller arc length on the circular domain.
 * @details The result is symmetric in its arguments, wraps across the seam when shorter,
 *          is exactly m/2 for antipodal points, and is always bounded by [0, m/2].
 */
inline void test_shortest_distance() {
  HS_EXPECT_NEAR(shortest_distance(0.0f, 1.0f, 10.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(shortest_distance(1.0f, 0.0f, 10.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(shortest_distance(1.0f, 9.0f, 10.0f), 2.0f, 1e-5f);
  HS_EXPECT_NEAR(shortest_distance(0.0f, 5.0f, 10.0f), 5.0f, 1e-5f);
  // Across the seam the wrap arc (2) must beat the direct arc (8).
  HS_EXPECT_TRUE(shortest_distance(1.0f, 9.0f, 10.0f) < 9.0f - 1.0f);
  for (int i = 0; i < 20; ++i) {
    float d = shortest_distance(i * 0.5f, 3.3f, 10.0f);
    HS_EXPECT_TRUE(d >= 0.0f && d <= 5.0f + 1e-5f);
  }
}

/**
 * @brief Verifies fwd_distance measures the arc from a to b in the positive direction only.
 * @details The distance wraps across the seam when b is "behind" a, yields 0 for equal
 *          points, lies in [0, m), and forward plus reverse forward distance equals m
 *          for distinct points.
 */
inline void test_fwd_distance() {
  HS_EXPECT_NEAR(fwd_distance(1.0f, 4.0f, 10.0f), 3.0f, 1e-5f);
  HS_EXPECT_NEAR(fwd_distance(9.0f, 1.0f, 10.0f), 2.0f, 1e-5f);
  HS_EXPECT_NEAR(fwd_distance(4.0f, 4.0f, 10.0f), 0.0f, 1e-5f);
  float f = fwd_distance(2.0f, 7.0f, 10.0f);
  float r = fwd_distance(7.0f, 2.0f, 10.0f);
  HS_EXPECT_NEAR(f + r, 10.0f, 1e-5f);
}

// --- apply_if_changed -------------------------------------------------------

/**
 * @brief Verifies apply_if_changed invokes the callable only when the value changes.
 * @details The callable fires only when the incoming value differs from the latched
 *          `last`, then `last` is updated — the live-slider debounce idiom. The test
 *          covers no-change (no call), change (one call, latched), and repeat (no call).
 */
inline void test_apply_if_changed() {
  int last = 5;
  int applied = -1;
  int call_count = 0;

  apply_if_changed(5, last, [&](int v) { applied = v; ++call_count; });
  HS_EXPECT_EQ(call_count, 0);
  HS_EXPECT_EQ(last, 5);

  apply_if_changed(8, last, [&](int v) { applied = v; ++call_count; });
  HS_EXPECT_EQ(call_count, 1);
  HS_EXPECT_EQ(applied, 8);
  HS_EXPECT_EQ(last, 8);

  apply_if_changed(8, last, [&](int v) { applied = v; ++call_count; });
  HS_EXPECT_EQ(call_count, 1);
  HS_EXPECT_EQ(last, 8);
}

/**
 * @brief Runs every util test case in the module.
 * @return The module's failure count, as reported by end_module.
 */
inline int run_util_tests() {
  hs_test::ModuleFixture fixture("util");

  test_wrap_float();
  test_wrap_t();
  test_wrap_int();
  test_wrap_mixed_type();
  test_fast_wrap();
  test_shortest_distance();
  test_fwd_distance();
  test_apply_if_changed();

  return fixture.result();
}

} // namespace util_tests
} // namespace hs_test
