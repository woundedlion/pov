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
#ifndef HOLOSPHERE_CORE_LED_H_
#define HOLOSPHERE_CORE_LED_H_
#include "platform.h"

// Uncomment to use DMA-based HD107S controller instead of FastLED WS2801.
// Requires Teensy 4.x hardware. Leave commented for WASM/sim builds.
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
 * @brief RAII guard to temporarily disable FastLED's color correction.
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
 * @brief RAII guard to temporarily disable FastLED's temperature correction.
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
#endif // HOLOSPHERE_CORE_LED_H_
