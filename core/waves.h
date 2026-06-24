/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "3dmath.h"
#include <cmath>
#include "util.h"

/**
 * @brief Builds a sine oscillator over the range [from, to].
 * @param from The output value at the trough.
 * @param to The output value at the peak.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset, in cycles (phase = 1 is a full cycle),
 *        matching tri_wave/square_wave so the three are interchangeable.
 * @return A lambda mapping time t to the wave value.
 */
inline auto sin_wave(float from, float to, float freq, float phase) {
  // Hoist only 2π·phase: reassociating freq·t·2π could shift the last bit, and
  // the wave must stay bit-identical between sim and device.
  const float phase_term = 2 * PI_F * phase;
  return [=](float t) -> float {
    // −π/2 anchors t=0, phase=0 at the trough.
    auto w = (sinf(freq * t * 2 * PI_F - (PI_F / 2) + phase_term) + 1) / 2;
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
    // wrap_t folds a negative phase into [0,1), keeping the triangle continuous.
    float w = wrap_t(t * freq + phase);
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
  HS_CHECK(dutyCycle >= 0.0f && dutyCycle <= 1.0f,
           "square_wave: dutyCycle must be in [0,1]");
  return [=](float t) -> float {
    // wrap_t, not raw fmod: fmod keeps the dividend's sign, so a negative phase
    // would stay < dutyCycle and wrongly latch the wave "on".
    if (wrap_t(t * freq + phase) < dutyCycle) {
      return to;
    }
    return from;
  };
}
