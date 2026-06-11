/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file led.h
 * @brief LED configuration constants and color correction RAII guards.
 *
 * This file provides the shared constants and correction helpers
 * that effects (including legacy effects) depend on. The actual
 * POVDisplay driver lives in hardware/pov_single.h and is included
 * directly by target .ino files.
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
struct NoColorCorrection {};
struct NoTempCorrection {};
#else
/**
 * @brief RAII guard to disable both color and temperature correction for its
 * scope, restoring TypicalLEDStrip/Candle defaults on destruction.
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
 * color correction) for its scope, restoring the Candle default on destruction.
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
