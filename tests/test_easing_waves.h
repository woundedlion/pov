/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/easing.h and core/waves.h.
 *
 * Pure scalar functions: cheap to test, easy to regress silently. We assert the
 * defining contract of each easing curve (anchored endpoints, finiteness,
 * monotonicity where the curve is monotone) and the bounded/shape properties of
 * the wave generators. t is sampled at exact rational fractions of [0,1] to
 * avoid float drift past 1.0 (where circular easings would go NaN).
 */
#pragma once

#include "core/easing.h"
#include "core/waves.h"
#include "tests/test_harness.h"

#include <cmath>

namespace hs_test {
namespace easing_waves_tests {

constexpr int N = 64; // sample count over [0,1]

static inline float frac(int i) { return static_cast<float>(i) / N; }

// Assert a function is finite over [0,1] and (optionally) monotone non-decreasing.
template <typename Fn>
static inline void check_curve(Fn f, bool monotone, const char *name) {
  float prev = f(0.0f);
  HS_EXPECT(std::isfinite(prev), name);
  for (int i = 1; i <= N; ++i) {
    float v = f(frac(i));
    HS_EXPECT(std::isfinite(v), name);
    if (monotone)
      HS_EXPECT(v >= prev - 1e-4f, name); // non-decreasing within tolerance
    prev = v;
  }
}

// ---------------------------------------------------------------------------
// easing.h — endpoints
// ---------------------------------------------------------------------------

inline void test_easing_endpoints() {
  // "In" curves anchor at f(0)=0; "out" curves anchor at f(1)=1; "in-out" both.
  HS_EXPECT_NEAR(ease_in_out_bicubic(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_out_bicubic(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_out_sin(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_out_sin(1.0f), 1.0f, 1e-5f);

  HS_EXPECT_NEAR(ease_in_sin(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_sin(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_cubic(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_cubic(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_in_circ(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_circ(1.0f), 1.0f, 1e-5f);

  // Special-cased endpoints (avoid powf domain issues at the boundary).
  HS_EXPECT_NEAR(ease_out_expo(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_expo(1.0f), 1.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_elastic(0.0f), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(ease_out_elastic(1.0f), 1.0f, 1e-5f);

  // Linear is the identity.
  HS_EXPECT_NEAR(ease_mid(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(ease_mid(0.5f), 0.5f, 1e-6f);
  HS_EXPECT_NEAR(ease_mid(1.0f), 1.0f, 1e-6f);
}

inline void test_easing_finite_and_monotone() {
  // Monotone non-decreasing curves.
  check_curve(ease_in_out_bicubic, true, "ease_in_out_bicubic");
  check_curve(ease_in_out_sin, true, "ease_in_out_sin");
  check_curve(ease_in_sin, true, "ease_in_sin");
  check_curve(ease_out_sin, true, "ease_out_sin");
  check_curve(ease_in_cubic, true, "ease_in_cubic");
  check_curve(ease_out_cubic, true, "ease_out_cubic");
  check_curve(ease_in_circ, true, "ease_in_circ");
  check_curve(ease_out_circ, true, "ease_out_circ");
  check_curve(ease_mid, true, "ease_mid");
  // Overshooting / non-monotone curves: finiteness only.
  check_curve(ease_out_expo, false, "ease_out_expo");
  check_curve(ease_out_elastic, false, "ease_out_elastic");
}

inline void test_easing_in_out_symmetry_midpoint() {
  // In-out curves pass through 0.5 at the midpoint.
  HS_EXPECT_NEAR(ease_in_out_sin(0.5f), 0.5f, 1e-4f);
  HS_EXPECT_NEAR(ease_in_out_bicubic(0.5f), 0.5f, 1e-4f);
}

// ---------------------------------------------------------------------------
// waves.h — sin / tri / square generators
// ---------------------------------------------------------------------------

inline void test_sin_wave_bounds_and_phase() {
  auto w = sin_wave(2.0f, 5.0f, 1.0f, 0.0f);
  for (int i = 0; i <= N; ++i) {
    float v = w(frac(i) * 3.0f);
    HS_EXPECT(std::isfinite(v), "sin_wave finite");
    HS_EXPECT_GE(v, 2.0f - 1e-3f);
    HS_EXPECT_LE(v, 5.0f + 1e-3f);
  }
  // Phase convention: at t=0 the wave starts at the `from` end.
  HS_EXPECT_NEAR(sin_wave(0.0f, 1.0f, 1.0f, 0.0f)(0.0f), 0.0f, 1e-3f);
}

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

inline void test_square_wave_binary() {
  auto w = square_wave(0.0f, 1.0f, 1.0f, 0.5f, 0.0f);
  // Only ever emits `from` or `to`.
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

// ---------------------------------------------------------------------------
// Runner
// ---------------------------------------------------------------------------

inline int run_easing_waves_tests() {
  auto scope = hs_test::begin_module("easing_waves");

  test_easing_endpoints();
  test_easing_finite_and_monotone();
  test_easing_in_out_symmetry_midpoint();

  test_sin_wave_bounds_and_phase();
  test_tri_wave_shape();
  test_square_wave_binary();

  return hs_test::end_module(scope);
}

} // namespace easing_waves_tests
} // namespace hs_test

