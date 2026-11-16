/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Scalar subject for Lerp cases and borrow contracts. */
#pragma once

namespace hs_test::animation_tests {

/**
 * @brief Minimal subject satisfying the lerp(start, target, t) interface
 * Animation::Lerp drives via type erasure.
 */
struct Lerpable {
  float value = 0.0f; /**< Interpolated scalar payload. */
  /**
   * @brief Linearly interpolates this subject's value between two endpoints.
   * @param a Start subject (source value).
   * @param b Target subject (destination value).
   * @param t Interpolation parameter in [0, 1].
   */
  void lerp(const Lerpable &a, const Lerpable &b, float t) {
    value = a.value + (b.value - a.value) * t;
  }
};

} // namespace hs_test::animation_tests
