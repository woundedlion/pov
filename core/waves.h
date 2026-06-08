/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "3dmath.h"
#include <cmath>
#include <functional>
#include "util.h"

/**
 * @brief Builds a sine oscillator over the range [from, to].
 * @param from The output value at the trough.
 * @param to The output value at the peak.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset, in cycles (scaled by 2 internally).
 * @return A lambda mapping time t to the wave value.
 */
inline auto sin_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    auto w = (sinf(freq * t * 2 * PI_F - (PI_F / 2) - (2 * phase)) + 1) / 2;
    return hs::lerp(from, to, w);
  };
}

/**
 * @brief Builds a triangle oscillator over the range [from, to].
 * @param from The output value at the trough.
 * @param to The output value at the peak.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset, in cycles.
 * @return A lambda mapping time t to the wave value.
 */
inline auto tri_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    float w = wrap(t * freq + phase, 1.0f);
    if (w < 0.5f) {
      w = 2.0f * w;
    } else {
      w = 2.0f * (1.0f - w);
    }
    return hs::lerp(from, to, w);
  };
}

/**
 * @brief Builds a square oscillator that toggles between from and to.
 * @param from The output value while "off" (low).
 * @param to The output value while "on" (high).
 * @param freq The frequency (cycles per unit time).
 * @param dutyCycle Fraction in [0, 1] of each cycle spent "on" (high).
 * @param phase The starting phase offset, in cycles.
 * @return A lambda mapping time t to the wave value.
 */
inline auto square_wave(float from, float to, float freq, float dutyCycle,
                        float phase) {
  return [=](float t) -> float {
    if (fmod(t * freq + phase, 1.0f) < dutyCycle) {
      return to;
    }
    return from;
  };
}
