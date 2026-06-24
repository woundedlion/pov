/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/easing.h and core/waves.h.
 *
 * Asserts the defining contract of each easing curve (anchored endpoints,
 * finiteness, monotonicity where the curve is monotone) and the bounded/shape
 * properties of the wave generators. t is sampled at exact rational fractions
 * of [0,1] to avoid float drift past 1.0 (where circular easings would go NaN).
 */
#pragma once

#include "core/easing.h"
#include "core/waves.h"
#include "tests/test_harness.h"

#include <cmath>

namespace hs_test {
namespace easing_waves_tests {

// Sample count over [0,1].
constexpr int N = 200;

/**
 * @brief Maps a sample index to its exact rational fraction of [0,1].
 * @param i Sample index in [0, N].
 * @return The fraction i/N as a float; exact at the rationals to avoid float
 *         drift past 1.0 (where circular easings would go NaN).
 */
static inline float frac(int i) { return static_cast<float>(i) / N; }

/**
 * @brief Asserts a curve is finite over [0,1] and optionally monotone.
 * @tparam Fn Callable taking a float t and returning a float.
 * @param f Curve to sample at frac(0)..frac(N).
 * @param monotone If true, also require each sample to be non-decreasing
 *        (within a 1e-4 tolerance) relative to the previous one.
 * @param name Label forwarded to HS_EXPECT for failure reporting.
 */
template <typename Fn>
static inline void check_curve(Fn f, bool monotone, const char *name) {
  float prev = f(0.0f);
  HS_EXPECT(std::isfinite(prev), name);
  for (int i = 1; i <= N; ++i) {
    float v = f(frac(i));
    HS_EXPECT(std::isfinite(v), name);
    if (monotone)
      HS_EXPECT(v >= prev - 1e-4f, name);
    prev = v;
  }
}

// ---------------------------------------------------------------------------
// easing.h — endpoints
// ---------------------------------------------------------------------------

/**
 * @brief Verifies each easing curve anchors its endpoints.
 * @details "in" curves satisfy f(0)=0, "out" curves f(1)=1, "in-out" curves
 *          both; ease_linear is the linear identity. Expo/elastic endpoints are
 *          special-cased to avoid powf domain issues at the boundary.
 */
inline void test_easing_endpoints() {
  HS_EXPECT_NEAR(ease_in_out_cubic(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_out_cubic(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_out_sin(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_out_sin(1.0f), 1.0f, 1e-5f);

  HS_EXPECT_NEAR(ease_in_sin(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_sin(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_cubic(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_cubic(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_circ(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_circ(1.0f), 1.0f, 1e-5f);

  HS_EXPECT_NEAR(ease_out_expo(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_expo(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_elastic(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_elastic(1.0f), 1.0f, 1e-5f);

  HS_EXPECT_NEAR(ease_linear(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(ease_linear(0.5f), 0.5f, 1e-6f);
  HS_EXPECT_NEAR(ease_linear(1.0f), 1.0f, 1e-6f);
}

/**
 * @brief Verifies every easing curve stays finite, and monotone where expected.
 * @details Non-overshooting curves are checked for monotone non-decreasing
 *          behavior; expo/elastic overshoot, so they get finiteness only.
 */
inline void test_easing_finite_and_monotone() {
  check_curve(ease_in_out_cubic, true, "ease_in_out_cubic");
  check_curve(ease_in_out_sin, true, "ease_in_out_sin");
  check_curve(ease_in_sin, true, "ease_in_sin");
  check_curve(ease_out_sin, true, "ease_out_sin");
  check_curve(ease_in_cubic, true, "ease_in_cubic");
  check_curve(ease_out_cubic, true, "ease_out_cubic");
  check_curve(ease_in_circ, true, "ease_in_circ");
  check_curve(ease_out_circ, true, "ease_out_circ");
  check_curve(ease_linear, true, "ease_linear");
  // Overshooting / non-monotone curves: finiteness only.
  check_curve(ease_out_expo, false, "ease_out_expo");
  check_curve(ease_out_elastic, false, "ease_out_elastic");
}

/**
 * @brief Pins reference interior values for the expo/elastic curves.
 * @details The finite-only interior check (test_easing_finite_and_monotone) would
 *          pass for any finite curve, so it cannot catch a shape regression in
 *          the two overshooting/non-monotone curves. Lock independently-derived
 *          analytic values at interior points:
 *            ease_out_expo(t)    = 1 - 2^(-10t)
 *              t=0.25 -> 1 - 2^-2.5  = 0.8232233
 *              t=0.5  -> 1 - 2^-5    = 0.96875
 *            ease_out_elastic(x) = 2^(-10x)·sin((10x-0.75)·2π/3) + 1
 *              x=0.25 -> 2^-2.5·sin(210°) + 1 = 0.9116117
 *              x=0.5  -> 2^-5·sin(150°)  + 1 = 1.015625  (overshoots past 1)
 *          The elastic x=0.5 value also pins the defining overshoot (> 1).
 */
inline void test_easing_expo_elastic_interior_reference() {
  HS_EXPECT_NEAR(ease_out_expo(0.25f), 0.8232233f, 1e-4f);
  HS_EXPECT_NEAR(ease_out_expo(0.5f), 0.96875f, 1e-5f);

  HS_EXPECT_NEAR(ease_out_elastic(0.25f), 0.9116117f, 1e-4f);
  HS_EXPECT_NEAR(ease_out_elastic(0.5f), 1.015625f, 1e-4f);
  HS_EXPECT_GT(ease_out_elastic(0.5f), 1.0f);
}

/**
 * @brief Verifies in-out curves are symmetric about the midpoint.
 * @details Each in-out curve passes through 0.5 at t=0.5.
 */
inline void test_easing_in_out_symmetry_midpoint() {
  HS_EXPECT_NEAR(ease_in_out_sin(0.5f), 0.5f, 1e-4f);
  HS_EXPECT_NEAR(ease_in_out_cubic(0.5f), 0.5f, 1e-4f);
}

// ---------------------------------------------------------------------------
// waves.h — sin / tri / square generators
// ---------------------------------------------------------------------------

/**
 * @brief Verifies sin_wave bounds and its forward phase convention.
 * @details The wave stays within [from,to] and starts at the `from` end at t=0.
 *          Phase advances forward (+phase) consistently with
 *          tri_wave/square_wave, so a quarter-cycle phase puts t=0 at the rising
 *          midpoint, matching tri_wave's direction at the same value.
 */
inline void test_sin_wave_bounds_and_phase() {
  auto w = sin_wave(2.0f, 5.0f, 1.0f, 0.0f);
  for (int i = 0; i <= N; ++i) {
    float v = w(frac(i) * 3.0f);
    HS_EXPECT(std::isfinite(v), "sin_wave finite");
    HS_EXPECT_GE(v, 2.0f - 1e-3f);
    HS_EXPECT_LE(v, 5.0f + 1e-3f);
  }
  HS_EXPECT_NEAR(sin_wave(0.0f, 1.0f, 1.0f, 0.0f)(0.0f), 0.0f, 1e-3f);
  // A quarter-cycle phase puts t=0 at the rising midpoint (forward phase).
  auto wp = sin_wave(0.0f, 1.0f, 1.0f, 0.25f);
  HS_EXPECT_NEAR(wp(0.0f), 0.5f, 1e-3f);
  HS_EXPECT(wp(0.05f) > 0.5f, "sin_wave phase 0.25 rises at t=0 like tri_wave");
  HS_EXPECT_NEAR(wp(0.0f), tri_wave(0.0f, 1.0f, 1.0f, 0.25f)(0.0f), 1e-3f);
}

/**
 * @brief Verifies tri_wave traces a symmetric triangle within bounds.
 * @details The wave runs trough to peak and back, staying within [from,to]
 *          across multiple periods.
 */
inline void test_tri_wave_shape() {
  auto w = tri_wave(0.0f, 1.0f, 1.0f, 0.0f);
  HS_EXPECT_NEAR(w(0.0f), 0.0f, 1e-4f);   // trough at start
  HS_EXPECT_NEAR(w(0.25f), 0.5f, 1e-4f);  // rising through mid
  HS_EXPECT_NEAR(w(0.5f), 1.0f, 1e-4f);   // peak
  HS_EXPECT_NEAR(w(0.75f), 0.5f, 1e-4f);  // falling through mid
  for (int i = 0; i <= N; ++i) {
    float v = w(frac(i) * 2.0f);
    HS_EXPECT_GE(v, -1e-3f);
    HS_EXPECT_LE(v, 1.0f + 1e-3f);
  }
}

/**
 * @brief Verifies square_wave emits only its two levels at the duty boundary.
 * @details The wave emits only `from` or `to`, switching at the duty-cycle
 *          boundary (here 0.5: first half of the period high, second half low).
 */
inline void test_square_wave_binary() {
  auto w = square_wave(0.0f, 1.0f, 1.0f, 0.5f, 0.0f);
  for (int i = 0; i <= N; ++i) {
    float v = w(frac(i) * 2.0f);
    bool is_from = hs_test::approx(v, 0.0f, 1e-5f);
    bool is_to = hs_test::approx(v, 1.0f, 1e-5f);
    HS_EXPECT_TRUE(is_from || is_to);
  }
  // Duty cycle: first half of the period is high, second half low.
  HS_EXPECT_NEAR(w(0.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(w(0.25f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(w(0.6f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(w(0.9f), 0.0f, 1e-5f);
}

/**
 * @brief Verifies square_wave stays periodic across the sign boundary.
 * @details A negative t*freq+phase must fold into [0,1) via wrap() so the wave
 *          stays periodic; a sign-preserving fold would leave a negative
 *          argument below the duty cycle and latch the wave permanently "on".
 *          Checks w(t-1) == w(t) over a full sweep.
 */
inline void test_square_wave_negative_phase() {
  auto w = square_wave(0.0f, 1.0f, 1.0f, 0.5f, 0.0f);
  HS_EXPECT_NEAR(w(-0.25f), 0.0f, 1e-5f); // wrap(-0.25)=0.75 -> low, not latched
  HS_EXPECT_NEAR(w(-0.9f), 1.0f, 1e-5f);  // wrap(-0.9)=0.1 -> high
  for (int i = 0; i < N; ++i) {
    float ft = frac(i); // [0,1)
    HS_EXPECT_NEAR(w(ft - 1.0f), w(ft), 1e-5f);
  }
}

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

/**
 * @brief Runs the full easing/waves suite under its own module scope.
 * @return Number of failures reported by the module.
 */
inline int run_easing_waves_tests() {
  auto scope = hs_test::begin_module("easing_waves");

  test_easing_endpoints();
  test_easing_finite_and_monotone();
  test_easing_expo_elastic_interior_reference();
  test_easing_in_out_symmetry_midpoint();

  test_sin_wave_bounds_and_phase();
  test_tri_wave_shape();
  test_square_wave_binary();
  test_square_wave_negative_phase();

  return hs_test::end_module(scope);
}

} // namespace easing_waves_tests
} // namespace hs_test

