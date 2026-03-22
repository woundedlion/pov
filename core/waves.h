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
 * @brief Generates a sine wave function with offset, amplitude, frequency, and
 * phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset.
 * @return A lambda that takes time (t) and returns a float.
 */
inline auto sin_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    auto w = (sinf(freq * t * 2 * PI_F - (PI_F / 2) - (2 * phase)) + 1) / 2;
    return lerp(from, to, w);
  };
}

/**
 * @brief Generates a triangle wave function with offset, amplitude, frequency,
 * and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset.
 * @return A lambda that takes time (t) and returns a float.
 */
inline auto tri_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    float w = wrap(t * freq + phase, 1.0f);
    if (w < 0.5f) {
      w = 2.0f * w;
    } else {
      w = 2.0f * (1.0f - w);
    }
    return lerp(from, to, w);
  };
}

/**
 * @brief Generates a square wave function with offset, amplitude, frequency,
 * duty cycle, and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param dutyCycle The percentage of time the wave is "on" (high).
 * @param phase The starting phase offset.
 * @return A lambda that takes time (t) and returns a float.
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
