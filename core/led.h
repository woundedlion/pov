/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file led.h
 * @brief LED configuration constants and color correction RAII guards.
 *
 * This file provides the shared constants and correction helpers
 * that effects depend on. The actual POVDisplay driver lives in
 * hardware/pov_single.h and is included directly by target .ino files.
 */
#pragma once
#include "platform.h"

// USE_DMA_LEDS selects the DMA-based HD107S controller instead of the FastLED
// WS2801 path. It requires Teensy 4.x hardware and is left undefined for WASM/sim
// and single-board builds. Targets that need it define USE_DMA_LEDS themselves
// before including the driver (e.g. targets/Phantasm/Phantasm.ino); this
// engine-wide header deliberately never defines it, so no commented-out toggle is
// kept here to be uncommented by accident.

/**
 * @brief Analog pin used for seeding the random number generator.
 */
static constexpr int PIN_RANDOM = 15;
/**
 * @brief Data pin for the LED strip
 */
static constexpr int PIN_DATA = 11;
/**
 * @brief Clock pin for the LED strip.
 */
static constexpr int PIN_CLOCK = 13;

// When using DMA LEDs, correction is done in the DMA pipeline — stubs only.
#ifdef USE_DMA_LEDS
/**
 * @brief No-op stub: the DMA pipeline applies no color/temperature correction.
 */
struct NoColorCorrection {};
/**
 * @brief No-op stub: the DMA pipeline applies no temperature correction.
 */
struct NoTempCorrection {};
#else
// CONTRACT — restore-to-baseline, NOT save/restore. The destructors below
// reinstate the engine's canonical baseline (TypicalLEDStrip color, Candle
// temperature — the same values the drivers set once at init, see
// pov_single.h / pov_segmented.h), not whatever correction happened to be
// active when the guard was constructed. FastLED exposes no getter for the
// current global correction (and on device it is per-controller), so a true
// save/restore is not available here. This is correct for the one supported
// use — a single per-effect *member* (see effects_legacy.h), with effects
// constructed and destroyed sequentially so the ambient correction is always
// the baseline. Consequence: these guards do NOT nest — constructing one inside
// a non-baseline correction scope would restore the baseline on exit, not the
// outer scope's value.
//
// That non-nesting precondition is enforced by the shared depth counter below:
// each guard traps (HS_CHECK, always-on so it survives the NDEBUG device build)
// if another guard is already live, surfacing the misuse on the bench instead of
// silently dropping the outer correction at runtime. Single-threaded by
// construction (the per-effect render loop), so a plain int needs no atomicity.

/**
 * @brief Shared liveness counter for the correction guards.
 * @return Reference to the single process-wide depth (0 = no guard active).
 * @details A function-local static keeps one instance across translation units
 * without an out-of-line definition. Incremented after the nesting check in each
 * guard's ctor and decremented in its dtor.
 */
inline int &correction_guard_depth() {
  static int depth = 0;
  return depth;
}

/**
 * @brief Reinstates the engine's canonical baseline (TypicalLEDStrip color,
 * Candle temperature) and releases the guard liveness counter.
 * @details Shared by both correction guards' destructors (see the
 * restore-to-baseline contract above).
 */
inline void restore_correction_baseline() {
  FastLED.setCorrection(TypicalLEDStrip);
  FastLED.setTemperature(Candle);
  --correction_guard_depth();
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
    HS_CHECK(correction_guard_depth() == 0,
             "correction guards do not nest (see contract above)");
    ++correction_guard_depth();
    FastLED.setCorrection(UncorrectedColor);
    FastLED.setTemperature(UncorrectedTemperature);
  }
  /**
   * @brief Restores the TypicalLEDStrip/Candle baseline (restore-to-baseline,
   * not the correction active at construction).
   */
  ~NoColorCorrection() { restore_correction_baseline(); }
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
    HS_CHECK(correction_guard_depth() == 0,
             "correction guards do not nest (see contract above)");
    ++correction_guard_depth();
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(UncorrectedTemperature);
  }
  /**
   * @brief Restores the TypicalLEDStrip/Candle baseline (restore-to-baseline,
   * not the correction active at construction).
   */
  ~NoTempCorrection() { restore_correction_baseline(); }
};
#endif // !USE_DMA_LEDS
