/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Interpolation endpoints, periodic wrapping, and normalized vector interpolation.
 */
#pragma once

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

#include "core/math/interpolate.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace interpolate_tests {

/** @brief Progress saturates to endpoints and passes NaN through. */
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

/** @brief Linear interpolation preserves endpoints and clamps out-of-range progress. */
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

/** @brief Positive endpoints interpolate geometrically; nonpositive endpoints fall back to linear. */
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

/** @brief Periodic interpolation takes the shortest arc with a negative half-period tie. */
inline void test_shortest_periodic() {
  constexpr float TAU = 6.28318530717958647692f;

  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, 0.0f, 1.0f), 0.25f);
  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, 1.0f, 1.0f), 0.75f);
  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, -1.0f, 1.0f), 0.25f);
  HS_EXPECT_EQ(interp::shortest_periodic(0.25f, 0.75f, 2.0f, 1.0f), 0.75f);

  // Interior arc, no wrap needed.
  HS_EXPECT_NEAR(interp::shortest_periodic(0.2f, 0.4f, 0.5f, 1.0f), 0.3f,
                 1e-6f);
  // The short arc increases through the wrap seam.
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

/** @brief Reference periodic wrapping through the standard remainder.
 * @tparam T Scalar type.
 * @param value Value to wrap.
 * @param period Domain period.
 * @return The half-open wrapped remainder. */
template <typename T> inline T fmod_wrap_reference(T value, T period) {
  T result = std::fmod(value, period);
  if (result < T{0})
    result += period;
  return result >= period ? T{0} : result;
}

/** @brief Reference shortest-arc interpolation through standard remainders.
 * @param a Start endpoint.
 * @param b End endpoint.
 * @param progress Interpolation progress.
 * @param period Domain period.
 * @return Wrapped interpolation or the supplied endpoint. */
inline float fmod_shortest_reference(float a, float b, float progress,
                                     float period) {
  if (progress <= 0.0f)
    return a;
  if (progress >= 1.0f)
    return b;
  if (period <= 0.0f)
    return a + (b - a) * progress;
  const float HALF = period * 0.5f;
  float delta = std::fmod(b - a, period);
  if (delta < -HALF)
    delta += period;
  if (delta >= HALF)
    delta -= period;
  return fmod_wrap_reference(a + delta * progress, period);
}

/** @brief Compares scalar bits, accepting any NaN representation.
 * @tparam T Floating-point scalar type.
 * @param actual Measured value.
 * @param expected Reference value. */
template <typename T> inline void expect_periodic_bits(T actual, T expected) {
  if (std::isnan(expected)) {
    HS_EXPECT_TRUE(std::isnan(actual));
  } else if constexpr (sizeof(T) == sizeof(uint32_t)) {
    HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual),
                 std::bit_cast<uint32_t>(expected));
  } else {
    HS_EXPECT_EQ(std::bit_cast<uint64_t>(actual),
                 std::bit_cast<uint64_t>(expected));
  }
}

/** @brief Wrapping matches fmod bits at boundaries and across finite inputs. */
inline void test_periodic_wrap_fmod_parity() {
  const auto check = []<typename T>() {
    const T INFINITY_VALUE = std::numeric_limits<T>::infinity();
    const T NAN_VALUE = std::numeric_limits<T>::quiet_NaN();
    const T MINIMUM = std::numeric_limits<T>::min();
    const T SUBNORMAL = std::numeric_limits<T>::denorm_min();
    const T MAXIMUM = std::numeric_limits<T>::max();
    for (T period : {T{0}, -T{0}, T{-1}, T{1}, T{2}, T{0.5}, MINIMUM, SUBNORMAL,
                     MAXIMUM, INFINITY_VALUE, -INFINITY_VALUE, NAN_VALUE}) {
      for (T value : {T{0},
                      -T{0},
                      T{0.25},
                      T{-0.25},
                      T{1},
                      T{-1},
                      MINIMUM,
                      -MINIMUM,
                      SUBNORMAL,
                      -SUBNORMAL,
                      MAXIMUM,
                      -MAXIMUM,
                      INFINITY_VALUE,
                      -INFINITY_VALUE,
                      NAN_VALUE,
                      period,
                      -period,
                      std::nextafter(period, T{0}),
                      std::nextafter(period, INFINITY_VALUE),
                      std::nextafter(-period, T{0}),
                      std::nextafter(-period, -INFINITY_VALUE)}) {
        expect_periodic_bits(math::wrap_positive(value, period),
                             fmod_wrap_reference(value, period));
      }
    }
  };
  check.template operator()<float>();
  check.template operator()<double>();
  uint32_t state = 0x72ec5217U;
  const auto next = [&] {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
  };
  for (int sample = 0; sample < 32768; ++sample) {
    const float VALUE = std::bit_cast<float>(next());
    const float PERIOD = std::bit_cast<float>(next() & 0x7fffffffU);
    expect_periodic_bits(math::wrap_positive(VALUE, PERIOD),
                         fmod_wrap_reference(VALUE, PERIOD));
  }
  for (uint64_t value :
       {uint64_t{0}, uint64_t{1}, uint64_t{1} << 54, (uint64_t{1} << 54) + 1,
        std::numeric_limits<uint64_t>::max() - 1}) {
    constexpr uint64_t PERIOD = uint64_t{1} << 63;
    HS_EXPECT_EQ(math::wrap_positive(value, PERIOD),
                 fmod_wrap_reference(value, PERIOD));
  }
}

/** @brief Shortest-arc interpolation preserves fmod bits over hue and extremes. */
inline void test_shortest_periodic_fmod_parity() {
  constexpr float INFINITY_VALUE = std::numeric_limits<float>::infinity();
  constexpr float NAN_VALUE = std::numeric_limits<float>::quiet_NaN();
  constexpr float SUBNORMAL = std::numeric_limits<float>::denorm_min();
  constexpr float MAXIMUM = std::numeric_limits<float>::max();
  for (float period : {0.0f, -0.0f, -1.0f, 1.0f, 2.0f, 0.5f, SUBNORMAL,
                       std::numeric_limits<float>::min(), MAXIMUM,
                       INFINITY_VALUE, -INFINITY_VALUE, NAN_VALUE}) {
    for (float a :
         {0.0f, -0.0f, SUBNORMAL, -SUBNORMAL, MAXIMUM, -MAXIMUM, INFINITY_VALUE,
          -INFINITY_VALUE, NAN_VALUE, period, -period,
          std::nextafter(period, 0.0f), std::nextafter(-period, 0.0f)}) {
      for (float b : {0.0f, -0.0f, SUBNORMAL, -SUBNORMAL, MAXIMUM, -MAXIMUM,
                      INFINITY_VALUE, -INFINITY_VALUE, NAN_VALUE, period,
                      -period, std::nextafter(period, INFINITY_VALUE),
                      std::nextafter(-period, -INFINITY_VALUE)}) {
        for (float progress :
             {-1.0f, -0.0f, 0.0f, 0.25f, 0.5f, 0.75f, 1.0f, 2.0f, NAN_VALUE}) {
          expect_periodic_bits(
              interp::shortest_periodic(a, b, progress, period),
              fmod_shortest_reference(a, b, progress, period));
        }
      }
    }
  }
  for (int a = 0; a < 65; ++a) {
    for (int b = 0; b < 65; ++b) {
      for (int progress = 0; progress <= 16; ++progress) {
        const float START = static_cast<float>(a) / 65.0f;
        const float END = static_cast<float>(b) / 65.0f;
        const float FRACTION = static_cast<float>(progress) / 16.0f;
        expect_periodic_bits(
            interp::shortest_periodic(START, END, FRACTION, 1.0f),
            fmod_shortest_reference(START, END, FRACTION, 1.0f));
      }
    }
  }
  uint32_t state = 0xb4572e1dU;
  const auto next = [&] {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
  };
  for (int sample = 0; sample < 32768; ++sample) {
    const float START = std::bit_cast<float>(next());
    const float END = std::bit_cast<float>(next());
    const float PERIOD = std::bit_cast<float>(next() & 0x7fffffffU);
    const float FRACTION = static_cast<float>(next() & 0xffffU) / 65535.0f;
    expect_periodic_bits(
        interp::shortest_periodic(START, END, FRACTION, PERIOD),
        fmod_shortest_reference(START, END, FRACTION, PERIOD));
  }
}

/** @brief Vector interpolation normalizes interior values and rejects degenerate directions. */
inline void test_normalized_linear() {
  // Endpoints come back as supplied, valid and not renormalized.
  const auto at_start =
      interp::normalized_linear<2>({3.0f, 4.0f}, {0.0f, 1.0f}, 0.0f, 1e-6f);
  HS_EXPECT_TRUE(at_start.valid);
  HS_EXPECT_EQ(at_start.values[0], 3.0f);
  HS_EXPECT_EQ(at_start.values[1], 4.0f);
  const auto at_end =
      interp::normalized_linear<2>({1.0f, 0.0f}, {0.0f, 5.0f}, 1.0f, 1e-6f);
  HS_EXPECT_TRUE(at_end.valid);
  HS_EXPECT_EQ(at_end.values[0], 0.0f);
  HS_EXPECT_EQ(at_end.values[1], 5.0f);
  const auto past_end =
      interp::normalized_linear<2>({1.0f, 0.0f}, {0.0f, 5.0f}, 2.0f, 1e-6f);
  HS_EXPECT_TRUE(past_end.valid);
  HS_EXPECT_EQ(past_end.values[0], 0.0f);
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
  HS_RUN(test_clamp_progress);
  HS_RUN(test_linear);
  HS_RUN(test_log_positive);
  HS_RUN(test_shortest_periodic);
  HS_RUN(test_periodic_wrap_fmod_parity);
  HS_RUN(test_shortest_periodic_fmod_parity);
  HS_RUN(test_normalized_linear);
  return fixture.result();
}

} // namespace interpolate_tests
} // namespace hs_test
