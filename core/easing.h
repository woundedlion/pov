/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
 
 // Based on https://easings.net/

#pragma once

#include <cmath>
#include "3dmath.h"

/**
 * @brief Easing function: Cubic Interpolation (In-Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_in_out_bicubic(float t) {
  return t < 0.5f ?
    4 * t * t * t :
    1 - (-2 * t + 2) * (-2 * t + 2) * (-2 * t + 2) / 2;
}

/**
 * @brief Easing function: Sine Interpolation (In-Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_in_out_sin(float t) {
  return -(cosf(PI_F * t) - 1) / 2;
}

/**
 * @brief Easing function: Sine Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_in_sin(float t) {
  return 1 - cosf((t * PI_F) / 2);
}

/**
 * @brief Easing function: Sine Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_out_sin(float t) {
  return sinf((t * PI_F) / 2);
}

/**
 * @brief Easing function: Cubic Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_in_cubic(float t) {
  return t * t * t;
}

/**
 * @brief Easing function: Circular Interpolation (In).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_in_circ(float t) {
  return 1 - sqrtf(1 - t * t);
}

/**
 * @brief Easing function: Linear interpolation (no easing).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_mid(float t) {
  return t;
}

/**
 * @brief Easing function: Exponential Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_out_expo(float t) {
  // Exact-endpoint guard (callers pass exactly 1.0f at sequence end): without it
  // the formula returns 1 - 2^-10 ~= 0.999 there. The compare is float-typed so
  // t is not promoted to double.
  return t == 1.0f ? 1.0f : 1.0f - powf(2.0f, -10.0f * t);
}

/**
 * @brief Easing function: Circular Interpolation (Out).
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_out_circ(float t) {
  return sqrtf(1 - (t - 1) * (t - 1));
}

/**
 * @brief Easing function: Cubic easing out.
 * @param t The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_out_cubic(float t) {
  return 1 - powf(1 - t, 3);
}

/**
 * @brief Easing function: Elastic easing out.
 * @param x The time factor (0.0 to 1.0).
 * @return The eased factor.
 */
static inline float ease_out_elastic(float x) {
  const float c4 = (2 * PI_F) / 3;
  // Exact-endpoint guards (callers pass exactly 0.0f / 1.0f at the ends); the
  // float-typed compares keep x in float.
  return x == 0.0f ?
    0.0f : x == 1.0f ?
    1.0f : powf(2.0f, -10.0f * x) * sinf((x * 10.0f - 0.75f) * c4) + 1.0f;
}
