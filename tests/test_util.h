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
#include "tests/test_harness.h"

namespace hs_test {
namespace util_tests {

// --- wrap(float, m) / wrap_t ------------------------------------------------

inline void test_wrap_float() {
  // In range is unchanged.
  HS_EXPECT_NEAR(wrap(0.25f, 1.0f), 0.25f, 1e-6f);
  // Negative folds back into [0, m).
  HS_EXPECT_NEAR(wrap(-0.25f, 1.0f), 0.75f, 1e-6f);
  // Above m folds down.
  HS_EXPECT_NEAR(wrap(1.25f, 1.0f), 0.25f, 1e-6f);
  // Multiple periods (fmod handles the bulk; one +m correction for the sign).
  HS_EXPECT_NEAR(wrap(-3.25f, 1.0f), 0.75f, 1e-6f);
  // A non-unit base.
  HS_EXPECT_NEAR(wrap(7.0f, 5.0f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(wrap(-1.0f, 5.0f), 4.0f, 1e-6f);

  // Result is always in [0, m) across several periods of both signs.
  for (int i = -30; i <= 30; ++i) {
    float w = wrap(i * 0.37f, 2.0f);
    HS_EXPECT_TRUE(w >= 0.0f && w < 2.0f);
  }
}

inline void test_wrap_t() {
  HS_EXPECT_NEAR(wrap_t(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(0.4f), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(1.4f), 0.4f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(-0.25f), 0.75f, 1e-6f);
  HS_EXPECT_NEAR(wrap_t(-2.1f), 0.9f, 1e-5f);
  for (int i = -50; i <= 50; ++i) {
    float w = wrap_t(i * 0.13f);
    HS_EXPECT_TRUE(w >= 0.0f && w < 1.0f);
  }
}

// --- wrap(int, int) ---------------------------------------------------------

inline void test_wrap_int() {
  HS_EXPECT_EQ(wrap(3, 5), 3);
  HS_EXPECT_EQ(wrap(5, 5), 0);
  HS_EXPECT_EQ(wrap(7, 5), 2);
  HS_EXPECT_EQ(wrap(-1, 5), 4);
  HS_EXPECT_EQ(wrap(-6, 5), 4);
  // Exact (int,int) match selects the integer overload (no float math).
  for (int i = -17; i <= 17; ++i) {
    int w = wrap(i, 6);
    HS_EXPECT_TRUE(w >= 0 && w < 6);
  }
}

// --- wrap(mixed int/float) --------------------------------------------------

// A mixed (int, float) call resolves to the template, which now returns
// std::common_type_t (float) instead of T=int — so the fractional part survives
// instead of being truncated to 0.
inline void test_wrap_mixed_type() {
  static_assert(std::is_same_v<decltype(wrap(3, 2.5f)), float>,
                "wrap(int, float) returns the common type (float)");
  HS_EXPECT_NEAR(wrap(3, 2.5f), 0.5f, 1e-6f);  // fmod(3,2.5)=0.5, not truncated
  HS_EXPECT_NEAR(wrap(7, 2.5f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(wrap(-1, 2.5f), 1.5f, 1e-6f); // negative folds up
  // (float, int) is unchanged: still float math (geometry.h relies on this).
  HS_EXPECT_NEAR(wrap(7.5f, 5), 2.5f, 1e-6f);
}

// --- fast_wrap --------------------------------------------------------------

inline void test_fast_wrap() {
  constexpr int W = 8;
  // In range: identity.
  for (int x = 0; x < W; ++x)
    HS_EXPECT_EQ(fast_wrap(x, W), x);
  // One period high: single subtraction.
  HS_EXPECT_EQ(fast_wrap(W, W), 0);
  HS_EXPECT_EQ(fast_wrap(W + 3, W), 3);
  HS_EXPECT_EQ(fast_wrap(2 * W - 1, W), W - 1);
  // One period low: single addition.
  HS_EXPECT_EQ(fast_wrap(-1, W), W - 1);
  HS_EXPECT_EQ(fast_wrap(-W, W), 0);
  // Whole legal [-W, 2W) window lands in [0, W).
  for (int x = -W; x < 2 * W; ++x) {
    int w = fast_wrap(x, W);
    HS_EXPECT_TRUE(w >= 0 && w < W);
  }
}

// --- shortest_distance / fwd_distance ---------------------------------------

inline void test_shortest_distance() {
  // Symmetric, in [0, m/2].
  HS_EXPECT_NEAR(shortest_distance(0.0f, 1.0f, 10.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(shortest_distance(1.0f, 0.0f, 10.0f), 1.0f, 1e-5f);
  // The short way wraps across the seam.
  HS_EXPECT_NEAR(shortest_distance(1.0f, 9.0f, 10.0f), 2.0f, 1e-5f);
  // Antipodal is exactly m/2.
  HS_EXPECT_NEAR(shortest_distance(0.0f, 5.0f, 10.0f), 5.0f, 1e-5f);
  // Always within [0, m/2].
  for (int i = 0; i < 20; ++i) {
    float d = shortest_distance(i * 0.5f, 3.3f, 10.0f);
    HS_EXPECT_TRUE(d >= 0.0f && d <= 5.0f + 1e-5f);
  }
}

inline void test_fwd_distance() {
  // Forward (positive direction), in [0, m).
  HS_EXPECT_NEAR(fwd_distance(1.0f, 4.0f, 10.0f), 3.0f, 1e-5f);
  // Wraps across the seam when b is "behind" a.
  HS_EXPECT_NEAR(fwd_distance(9.0f, 1.0f, 10.0f), 2.0f, 1e-5f);
  HS_EXPECT_NEAR(fwd_distance(4.0f, 4.0f, 10.0f), 0.0f, 1e-5f);
  // fwd + reverse forward distance == m for distinct points.
  float f = fwd_distance(2.0f, 7.0f, 10.0f);
  float r = fwd_distance(7.0f, 2.0f, 10.0f);
  HS_EXPECT_NEAR(f + r, 10.0f, 1e-5f);
}

// --- apply_if_changed -------------------------------------------------------

inline void test_apply_if_changed() {
  int last = 5;
  int applied = -1;
  int call_count = 0;

  // No change: callable not invoked, last untouched.
  apply_if_changed(5, last, [&](int v) { applied = v; ++call_count; });
  HS_EXPECT_EQ(call_count, 0);
  HS_EXPECT_EQ(last, 5);

  // Change: callable invoked once with the new value, last latched.
  apply_if_changed(8, last, [&](int v) { applied = v; ++call_count; });
  HS_EXPECT_EQ(call_count, 1);
  HS_EXPECT_EQ(applied, 8);
  HS_EXPECT_EQ(last, 8);

  // Same value again: no second invocation.
  apply_if_changed(8, last, [&](int v) { applied = v; ++call_count; });
  HS_EXPECT_EQ(call_count, 1);
  HS_EXPECT_EQ(last, 8);
}

inline int run_util_tests() {
  auto scope = hs_test::begin_module("util");

  test_wrap_float();
  test_wrap_t();
  test_wrap_int();
  test_wrap_mixed_type();
  test_fast_wrap();
  test_shortest_distance();
  test_fwd_distance();
  test_apply_if_changed();

  return hs_test::end_module(scope);
}

} // namespace util_tests
} // namespace hs_test
