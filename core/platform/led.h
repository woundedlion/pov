/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/platform.h"

/**
 * @file led.h
 * @brief LED pin constants and the color-correction RAII guards effects share.
 */

// USE_DMA_LEDS selects the HD107S DMA driver; undefined for WASM/sim and
// FastLED builds.

/**
 * @brief Analog pin used for seeding the random number generator.
 */
inline constexpr int PIN_RANDOM = 15;
/**
 * @brief Data pin for the LED strip
 */
inline constexpr int PIN_DATA = 11;
/**
 * @brief Clock pin for the LED strip.
 */
inline constexpr int PIN_CLOCK = 13;

// At most one NoColorCorrection or NoTempCorrection may be live at a time; a
// second live guard of either type traps.

/**
 * @brief Shared liveness flag for the correction guards.
 * @return Reference to the single process-wide flag (false = no guard active).
 * @note Non-atomic; construct/destroy correction guards only from the main
 * loop.
 */
inline bool &correction_guard_live() {
  static bool live = false;
  return live;
}

inline void acquire_correction_guard() {
  HS_CHECK(!correction_guard_live(),
           "NoColorCorrection and NoTempCorrection guards cannot overlap");
  correction_guard_live() = true;
}

// When using DMA LEDs, correction is done in the DMA pipeline — the guards carry
// only the liveness flag.
#ifdef USE_DMA_LEDS
/**
 * @brief Scope guard with no effect on the DMA driver's configured correction.
 */
struct NoColorCorrection {
  NoColorCorrection() { acquire_correction_guard(); }
  ~NoColorCorrection() { correction_guard_live() = false; }
  NoColorCorrection(const NoColorCorrection &) = delete;
  NoColorCorrection &operator=(const NoColorCorrection &) = delete;
};
/**
 * @brief Scope guard with no effect on the DMA driver's configured temperature.
 */
struct NoTempCorrection {
  NoTempCorrection() { acquire_correction_guard(); }
  ~NoTempCorrection() { correction_guard_live() = false; }
  NoTempCorrection(const NoTempCorrection &) = delete;
  NoTempCorrection &operator=(const NoTempCorrection &) = delete;
};
#else
// The destructors reinstate the baseline (TypicalLEDStrip color, Candle
// temperature), not the correction active at construction.

/**
 * @brief Reinstates the engine's canonical baseline (TypicalLEDStrip color,
 * Candle temperature) and clears the guard liveness flag.
 */
inline void restore_correction_baseline() {
  FastLED.setCorrection(TypicalLEDStrip);
  FastLED.setTemperature(Candle);
  correction_guard_live() = false;
}

/**
 * @brief RAII guard to disable both color and temperature correction for its
 * scope, restoring the TypicalLEDStrip/Candle baseline on destruction.
 */
struct NoColorCorrection {
  /**
   * @brief Disables both color and temperature correction for the guard's scope.
   */
  NoColorCorrection() {
    acquire_correction_guard();
    FastLED.setCorrection(UncorrectedColor);
    FastLED.setTemperature(UncorrectedTemperature);
  }
  /**
   * @brief Restores the TypicalLEDStrip/Candle baseline (restore-to-baseline,
   * not the correction active at construction).
   */
  ~NoColorCorrection() { restore_correction_baseline(); }
  NoColorCorrection(const NoColorCorrection &) = delete;
  NoColorCorrection &operator=(const NoColorCorrection &) = delete;
};

/**
 * @brief RAII guard to disable temperature correction (keeping TypicalLEDStrip
 * color correction) for its scope, restoring the Candle baseline on
 * destruction.
 */
struct NoTempCorrection {
  /**
   * @brief Disables temperature correction while keeping TypicalLEDStrip color
   * correction for the guard's scope.
   */
  NoTempCorrection() {
    acquire_correction_guard();
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(UncorrectedTemperature);
  }
  /**
   * @brief Restores the TypicalLEDStrip/Candle baseline (restore-to-baseline,
   * not the correction active at construction).
   */
  ~NoTempCorrection() { restore_correction_baseline(); }
  NoTempCorrection(const NoTempCorrection &) = delete;
  NoTempCorrection &operator=(const NoTempCorrection &) = delete;
};
#endif // !USE_DMA_LEDS
