/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cmath>
#include <limits>

#include "core/math/interpolate.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace interpolate_tests {

inline void test_clamp_progress() {
  static_assert(interp::clamp_progress(-1.0f) == 0.0f);
  static_assert(interp::clamp_progress(2.0f) == 1.0f);
  static_assert(interp::clamp_progress(0.25f) == 0.25f);

  HS_EXPECT_EQ(interp::clamp_progress(-1e30f), 0.0f);
  HS_EXPECT_EQ(interp::clamp_progress(0.0f), 0.0f);
  HS_EXPECT_EQ(interp::clamp_progress(0.5f), 0.5f);
  HS_EXPECT_EQ(interp::clamp_progress(1.0f), 1.0f);
  HS_EXPECT_EQ(interp::clamp_progress(1e30f), 1.0f);
  // Neither comparison holds for NaN, so it passes through rather than
  // collapsing onto an endpoint.
  HS_EXPECT_TRUE(std::isnan(
      interp::clamp_progress(std::numeric_limits<float>::quiet_NaN())));
}

inline void test_linear() {
  static_assert(interp::linear(2.0f, 6.0f, 0.5f) == 4.0f);

  HS_EXPECT_EQ(interp::linear(2.0f, 6.0f, 0.0f), 2.0f);
  HS_EXPECT_EQ(interp::linear(2.0f, 6.0f, 1.0f), 6.0f);
  HS_EXPECT_EQ(interp::linear(2.0f, 6.0f, 0.25f), 3.0f);
  HS_EXPECT_EQ(interp::linear(-4.0f, 4.0f, 0.5f), 0.0f);
  HS_EXPECT_EQ(interp::linear(3.0f, 3.0f, 0.5f), 3.0f);
  // Out of range snaps to the nearer endpoint instead of extrapolating.
  HS_EXPECT_EQ(interp::linear(2.0f, 6.0f, -1.0f), 2.0f);
  HS_EXPECT_EQ(interp::linear(2.0f, 6.0f, 4.0f), 6.0f);
  // A descending pair is ordinary, not a reversed domain.
  HS_EXPECT_EQ(interp::linear(6.0f, 2.0f, 0.5f), 4.0f);
}

inline void test_log_positive() {
  HS_EXPECT_EQ(interp::log_positive(1.0f, 4.0f, 0.0f), 1.0f);
  HS_EXPECT_EQ(interp::log_positive(1.0f, 4.0f, 1.0f), 4.0f);
  HS_EXPECT_NEAR(interp::log_positive(1.0f, 4.0f, 0.5f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(interp::log_positive(2.0f, 8.0f, 0.5f), 4.0f, 1e-6f);
  HS_EXPECT_NEAR(interp::log_positive(1.0f, 8.0f, 1.0f / 3.0f), 2.0f, 1e-5f);
  HS_EXPECT_EQ(interp::log_positive(1.0f, 4.0f, -1.0f), 1.0f);
  HS_EXPECT_EQ(interp::log_positive(1.0f, 4.0f, 2.0f), 4.0f);

  // Constant ratio: equal steps in progress multiply by a constant factor.
  const float third = interp::log_positive(1.0f, 27.0f, 1.0f / 3.0f);
  const float two_thirds = interp::log_positive(1.0f, 27.0f, 2.0f / 3.0f);
  HS_EXPECT_NEAR(third, 3.0f, 1e-4f);
  HS_EXPECT_NEAR(two_thirds / third, third, 1e-4f);

  // A non-positive endpoint has no logarithm; the fallback is linear, not NaN.
  HS_EXPECT_EQ(interp::log_positive(0.0f, 4.0f, 0.5f),
               interp::linear(0.0f, 4.0f, 0.5f));
  HS_EXPECT_EQ(interp::log_positive(1.0f, 0.0f, 0.5f),
               interp::linear(1.0f, 0.0f, 0.5f));
  HS_EXPECT_EQ(interp::log_positive(-2.0f, -8.0f, 0.5f),
               interp::linear(-2.0f, -8.0f, 0.5f));
  HS_EXPECT_TRUE(std::isfinite(interp::log_positive(0.0f, 4.0f, 0.5f)));
}

inline void test_shortest_periodic() {
  constexpr float TAU = 6.28318530717958647692f;

  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, 0.0f, 1.0f), 0.25f);
  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, 1.0f, 1.0f), 0.75f);
  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, -1.0f, 1.0f), 0.25f);
  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, 2.0f, 1.0f), 0.75f);

  // Interior arc, no wrap needed.
  HS_EXPECT_NEAR(interp::shortest_periodic(0.2f, 0.4f, 0.5f, 1.0f), 0.3f,
                 1e-6f);
  // The short arc runs backwards through 0, not forwards across 0.9 -> 0.1.
  HS_EXPECT_NEAR(interp::shortest_periodic(0.9f, 0.1f, 0.5f, 1.0f), 0.0f,
                 1e-6f);
  HS_EXPECT_NEAR(interp::shortest_periodic(0.9f, 0.1f, 0.25f, 1.0f), 0.95f,
                 1e-6f);
  // Exactly half a period: delta >= half takes the negative arc.
  HS_EXPECT_NEAR(interp::shortest_periodic(0.0f, 0.5f, 0.5f, 1.0f), 0.75f,
                 1e-6f);
  // An interpolant below zero is wrapped back up into range.
  HS_EXPECT_NEAR(interp::shortest_periodic(0.05f, 0.85f, 0.5f, 1.0f), 0.95f,
                 1e-6f);
  // Non-unit period.
  HS_EXPECT_NEAR(interp::shortest_periodic(0.0f, TAU * 0.25f, 0.5f, TAU),
                 TAU * 0.125f, 1e-5f);

  // Every result lands in [0, period) wherever the endpoints sit.
  for (int a_step = 0; a_step < 12; ++a_step) {
    for (int b_step = 0; b_step < 12; ++b_step) {
      const float a = static_cast<float>(a_step) / 12.0f;
      const float b = static_cast<float>(b_step) / 12.0f;
      for (int t_step = 1; t_step < 8; ++t_step) {
        const float value = interp::shortest_periodic(
            a, b, static_cast<float>(t_step) / 8.0f, 1.0f);
        HS_EXPECT_GE(value, 0.0f);
        HS_EXPECT_LT(value, 1.0f);
      }
    }
  }

  // No arc is defined without a period; the fallback is linear.
  HS_EXPECT_EQ(interp::shortest_periodic(0.9f, 0.1f, 0.5f, 0.0f),
               interp::linear(0.9f, 0.1f, 0.5f));
  HS_EXPECT_EQ(interp::shortest_periodic(0.9f, 0.1f, 0.5f, -1.0f),
               interp::linear(0.9f, 0.1f, 0.5f));
}

inline void test_normalized_linear() {
  // Endpoints come back as supplied — valid, and deliberately not
  // renormalized.
  const auto at_start =
      interp::normalized_linear<2>({3.0f, 4.0f}, {0.0f, 1.0f}, 0.0f, 1e-6f);
  HS_EXPECT_TRUE(at_start.valid);
  HS_EXPECT_EQ(at_start.values[0], 3.0f);
  HS_EXPECT_EQ(at_start.values[1], 4.0f);
  const auto at_end =
      interp::normalized_linear<2>({1.0f, 0.0f}, {0.0f, 5.0f}, 1.0f, 1e-6f);
  HS_EXPECT_TRUE(at_end.valid);
  HS_EXPECT_EQ(at_end.values[1], 5.0f);
  const auto past_end =
      interp::normalized_linear<2>({1.0f, 0.0f}, {0.0f, 5.0f}, 2.0f, 1e-6f);
  HS_EXPECT_TRUE(past_end.valid);
  HS_EXPECT_EQ(past_end.values[1], 5.0f);

  const auto midpoint =
      interp::normalized_linear<2>({1.0f, 0.0f}, {0.0f, 1.0f}, 0.5f, 1e-6f);
  HS_EXPECT_TRUE(midpoint.valid);
  HS_EXPECT_NEAR(midpoint.values[0], 0.70710678f, 1e-6f);
  HS_EXPECT_NEAR(midpoint.values[1], 0.70710678f, 1e-6f);

  // Renormalization means a non-unit input still yields a unit interpolant.
  const auto scaled = interp::normalized_linear<3>(
      {2.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 3.0f}, 0.5f, 1e-6f);
  HS_EXPECT_TRUE(scaled.valid);
  const float norm = std::sqrt(scaled.values[0] * scaled.values[0] +
                               scaled.values[1] * scaled.values[1] +
                               scaled.values[2] * scaled.values[2]);
  HS_EXPECT_NEAR(norm, 1.0f, 1e-6f);
  HS_EXPECT_EQ(scaled.values[1], 0.0f);

  // Antipodal endpoints pass through the origin: no direction to report, and
  // the un-normalized interpolant comes back rather than a NaN.
  const auto antipodal =
      interp::normalized_linear<2>({1.0f, 0.0f}, {-1.0f, 0.0f}, 0.5f, 1e-6f);
  HS_EXPECT_FALSE(antipodal.valid);
  HS_EXPECT_EQ(antipodal.values[0], 0.0f);
  HS_EXPECT_EQ(antipodal.values[1], 0.0f);

  // epsilon is an inclusive floor: a norm of exactly epsilon is invalid.
  const auto at_epsilon =
      interp::normalized_linear<2>({2.0f, 0.0f}, {0.0f, 0.0f}, 0.5f, 1.0f);
  HS_EXPECT_FALSE(at_epsilon.valid);
  HS_EXPECT_EQ(at_epsilon.values[0], 1.0f);
  const auto above_epsilon =
      interp::normalized_linear<2>({2.0f, 0.0f}, {0.0f, 0.0f}, 0.5f, 0.5f);
  HS_EXPECT_TRUE(above_epsilon.valid);
  HS_EXPECT_NEAR(above_epsilon.values[0], 1.0f, 1e-6f);

  // A single-component vector normalizes to a sign.
  const auto scalar = interp::normalized_linear<1>({4.0f}, {8.0f}, 0.5f, 1e-6f);
  HS_EXPECT_TRUE(scalar.valid);
  HS_EXPECT_NEAR(scalar.values[0], 1.0f, 1e-6f);
}

inline int run_interpolate_tests() {
  ModuleFixture fixture("interpolate");
  test_clamp_progress();
  test_linear();
  test_log_positive();
  test_shortest_periodic();
  test_normalized_linear();
  return fixture.result();
}

} // namespace interpolate_tests
} // namespace hs_test
