/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file led.h
 * @brief LED configuration constants and color correction RAII guards.
 *
 * This file provides the shared constants and correction helpers
 * that effects depend on. The actual POVDisplay driver lives in
 * hardware/pov_single.h and is included directly by target .ino files.
 */
#pragma once
#include "engine/platform.h"

// USE_DMA_LEDS selects the DMA-based HD107S controller instead of the FastLED
// WS2801 path. Requires Teensy 4.x and is left undefined for WASM/sim and
// single-board builds. Targets that need it define it themselves before
// including the driver (e.g. targets/common/phantasm_target.h); this engine-wide
// header never defines it.

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

// CONTRACT — at most ONE guard may be live at a time: NoColorCorrection and
// NoTempCorrection share the single liveness flag below, so a second live guard
// of EITHER type traps. Enforced in both the FastLED and DMA-stub builds.

/**
 * @brief Shared liveness flag for the correction guards.
 * @return Reference to the single process-wide flag (false = no guard active).
 * @details A function-local static keeps one instance across translation units
 * without an out-of-line definition. Set after the liveness check in each
 * guard's ctor and cleared in its dtor.
 * @note Non-atomic and main-loop-only — never construct/destroy a correction
 * guard from an ISR or any preemptive context.
 */
inline bool &correction_guard_live() {
  static bool live = false;
  return live;
}

// When using DMA LEDs, correction is done in the DMA pipeline — the guards carry
// only the liveness flag.
#ifdef USE_DMA_LEDS
/**
 * @brief No-op stub: the DMA pipeline applies no color/temperature correction.
 */
struct NoColorCorrection {
  NoColorCorrection() {
    HS_CHECK(!correction_guard_live(),
             "at most one correction guard may be live at a time (see contract "
             "above)");
    correction_guard_live() = true;
  }
  ~NoColorCorrection() { correction_guard_live() = false; }
  NoColorCorrection(const NoColorCorrection &) = delete;
  NoColorCorrection &operator=(const NoColorCorrection &) = delete;
};
/**
 * @brief No-op stub: the DMA pipeline applies no temperature correction.
 * @details A distinct type, not an alias of NoColorCorrection: the FastLED
 * branch defines two, and overload/if-constexpr dispatch must resolve the same
 * way on both.
 */
struct NoTempCorrection {
  NoTempCorrection() {
    HS_CHECK(!correction_guard_live(),
             "at most one correction guard may be live at a time (see contract "
             "above)");
    correction_guard_live() = true;
  }
  ~NoTempCorrection() { correction_guard_live() = false; }
  NoTempCorrection(const NoTempCorrection &) = delete;
  NoTempCorrection &operator=(const NoTempCorrection &) = delete;
};
#else
// CONTRACT — restore-to-baseline, NOT save/restore: the destructors reinstate
// the engine's canonical baseline (TypicalLEDStrip color, Candle temperature),
// not the correction active at construction (FastLED exposes no getter).

/**
 * @brief Reinstates the engine's canonical baseline (TypicalLEDStrip color,
 * Candle temperature) and clears the guard liveness flag.
 * @details Shared by both correction guards' destructors (see the
 * restore-to-baseline contract above).
 */
inline void restore_correction_baseline() {
  FastLED.setCorrection(TypicalLEDStrip);
  FastLED.setTemperature(Candle);
  correction_guard_live() = false;
}

/**
 * @brief RAII guard to disable both color and temperature correction for its
 * scope, restoring the TypicalLEDStrip/Candle baseline on destruction (see the
 * restore-to-baseline contract above).
 */
struct NoColorCorrection {
  /**
   * @brief Disables both color and temperature correction for the guard's scope.
   */
  NoColorCorrection() {
    HS_CHECK(!correction_guard_live(),
             "at most one correction guard may be live at a time (see contract "
             "above)");
    correction_guard_live() = true;
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
 * color correction) for its scope, restoring the Candle baseline on destruction
 * (see the restore-to-baseline contract above).
 */
struct NoTempCorrection {
  /**
   * @brief Disables temperature correction while keeping TypicalLEDStrip color
   * correction for the guard's scope.
   */
  NoTempCorrection() {
    HS_CHECK(!correction_guard_live(),
             "at most one correction guard may be live at a time (see contract "
             "above)");
    correction_guard_live() = true;
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
