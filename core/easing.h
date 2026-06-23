/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
 
 // Based on https://easings.net/

#pragma once

#include <cmath>
#include "3dmath.h"

// Convention for every easing below: t is the normalized time factor in [0, 1]
// and the functions are UNCLAMPED by design (matching the easings.net reference
// they mirror). They neither clamp the input nor bound the output, so an
// out-of-[0,1] t extrapolates the curve — the cubic/back/elastic variants in
// particular can return values outside [0, 1]. Callers that feed an unbounded
// drive (or index a fixed table with the result) must clamp t, or the output,
// themselves.

/**
 * @brief Easing function: Cubic Interpolation (In-Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_in_out_cubic(float t) {
  // No endpoint guard (unlike circ/expo/elastic below): this is a pure
  // polynomial already exact at the ends — 4t³ = 0 at t=0 and
  // 1 - (-2t+2)³/2 = 1 at t=1 — with no sqrt (NaN) or 2^x (unbounded) term that
  // would need pinning. The asymmetry is deliberate.
  return t < 0.5f ?
    4 * t * t * t :
    1 - (-2 * t + 2) * (-2 * t + 2) * (-2 * t + 2) / 2;
}

/**
 * @brief Easing function: Sine Interpolation (In-Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_in_out_sin(float t) {
  return -(cosf(PI_F * t) - 1) / 2;
}

/**
 * @brief Easing function: Sine Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_in_sin(float t) {
  return 1 - cosf((t * PI_F) / 2);
}

/**
 * @brief Easing function: Sine Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_sin(float t) {
  return sinf((t * PI_F) / 2);
}

/**
 * @brief Easing function: Cubic Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_in_cubic(float t) {
  return t * t * t;
}

/**
 * @brief Easing function: Circular Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_in_circ(float t) {
  // Clamp the radicand: callers may pass t slightly outside [0,1] at sequence
  // ends, and 1 - t*t goes negative there -> sqrtf(NaN). Matches the guarded
  // expo/elastic easings below.
  return 1 - sqrtf(fmaxf(0.0f, 1 - t * t));
}

/**
 * @brief Easing function: Linear interpolation (no easing).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_linear(float t) {
  return t;
}

/**
 * @brief Easing function: Exponential Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_expo(float t) {
  // Endpoint guards (callers pass exactly 0.0f / 1.0f at sequence ends). The
  // upper guard pins 1.0f, which the formula otherwise approaches as 1 - 2^-10
  // ~= 0.999. The lower guard floors at 0: 2^(-10t) grows without bound for
  // t < 0, so an out-of-range t would otherwise return a wildly negative value
  // (the radicand-clamped siblings -- ease_*_circ -- bound their out-of-range
  // ends the same way). Compares are float-typed so t is not promoted to double.
  return t <= 0.0f ? 0.0f : t == 1.0f ? 1.0f : 1.0f - powf(2.0f, -10.0f * t);
}

/**
 * @brief Easing function: Circular Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_circ(float t) {
  // Clamp the radicand; see ease_in_circ (NaN for t just outside [0,1]).
  return sqrtf(fmaxf(0.0f, 1 - (t - 1) * (t - 1)));
}

/**
 * @brief Easing function: Cubic easing out.
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_cubic(float t) {
  return 1 - powf(1 - t, 3);
}

/**
 * @brief Easing function: Elastic easing out.
 * @param x The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_elastic(float x) {
  const float c4 = (2 * PI_F) / 3;
  // Endpoint guards (callers pass exactly 0.0f / 1.0f at the ends). The lower
  // guard floors at 0 for all x <= 0, not just x == 0: 2^(-10x) grows without
  // bound for x < 0, so an out-of-range x would otherwise scale the sine into a
  // wildly out-of-range value (the radicand-clamped siblings bound their ends
  // the same way). Float-typed compares keep x in float.
  return x <= 0.0f ?
    0.0f : x == 1.0f ?
    1.0f : powf(2.0f, -10.0f * x) * sinf((x * 10.0f - 0.75f) * c4) + 1.0f;
}
