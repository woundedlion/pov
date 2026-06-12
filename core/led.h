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

// Selects the DMA-based HD107S controller instead of the FastLED WS2801 path.
// Requires Teensy 4.x hardware; leave undefined for WASM/sim and single-board
// builds. Targets that need it define USE_DMA_LEDS themselves before including
// the driver (e.g. targets/Phantasm/Phantasm.ino) rather than toggling it here,
// so this engine-wide header stays neutral. Uncomment only for a one-off local
// experiment:
// #define USE_DMA_LEDS

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
/// No-op stub: the DMA pipeline applies no color/temperature correction.
struct NoColorCorrection {};
/// No-op stub: the DMA pipeline applies no temperature correction.
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

/**
 * @brief RAII guard to disable both color and temperature correction for its
 * scope, restoring the TypicalLEDStrip/Candle baseline on destruction (see the
 * restore-to-baseline contract above).
 */
struct NoColorCorrection {
  NoColorCorrection() {
    FastLED.setCorrection(UncorrectedColor);
    FastLED.setTemperature(UncorrectedTemperature);
  }
  ~NoColorCorrection() {
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
  }
};

/**
 * @brief RAII guard to disable temperature correction (keeping TypicalLEDStrip
 * color correction) for its scope, restoring the Candle baseline on destruction
 * (see the restore-to-baseline contract above).
 */
struct NoTempCorrection {
  NoTempCorrection() {
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(UncorrectedTemperature);
  }
  ~NoTempCorrection() {
    FastLED.setCorrection(TypicalLEDStrip);
    FastLED.setTemperature(Candle);
  }
};
#endif // !USE_DMA_LEDS
