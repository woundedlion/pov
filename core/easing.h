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
  // Clamp the radicand: 1 - t*t goes negative for t outside [-1,1] -> sqrt NaN.
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
  // Endpoint guards: the upper pins exactly 1.0f (the formula only reaches
  // 1 - 2^-10), the lower floors at 0 (2^(-10t) explodes for t < 0).
  return t <= 0.0f ? 0.0f : t == 1.0f ? 1.0f : 1.0f - powf(2.0f, -10.0f * t);
}

/**
 * @brief Easing function: Circular Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_circ(float t) {
  // Clamp the radicand; see ease_in_circ.
  return sqrtf(fmaxf(0.0f, 1 - (t - 1) * (t - 1)));
}

/**
 * @brief Easing function: Cubic easing out.
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_cubic(float t) {
  float u = 1 - t;
  return 1 - u * u * u;
}

/**
 * @brief Easing function: Elastic easing out.
 * @param x The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
inline float ease_out_elastic(float x) {
  const float c4 = (2 * PI_F) / 3;
  // Endpoint guards: pin exactly 0/1 and floor at 0 (2^(-10x) explodes for x < 0).
  return x <= 0.0f ?
    0.0f : x == 1.0f ?
    1.0f : powf(2.0f, -10.0f * x) * sinf((x * 10.0f - 0.75f) * c4) + 1.0f;
}
