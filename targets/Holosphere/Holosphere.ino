/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file Holosphere.ino
 * @brief Holosphere — Single-Teensy POV display.
 *
 * Target: Teensy 4.0
 * Physical LEDs: 40 (20 per arm × 2 sides)
 * Virtual canvas: 96×20
 */

#include <FastLED.h>
#include <SPI.h>
#include <new>

#include "pov_single.h"
#include "targets/effects.h"
#include "core/engine/effects_legacy.h"

static constexpr int NUM_PIXELS = 40;    ///< Physical LEDs, both arm sides.
static constexpr unsigned int RPM = 480; ///< Nominal rotation speed.

#ifdef USE_DMA_LEDS
// Explicit specialization keeps ledController's DMAMEM section (pov_single.h).
HS_DEFINE_POV_SINGLE_LED_CONTROLLER(NUM_PIXELS, RPM);
#endif

/// Single-Teensy POV driver for this rig.
using POV = POVDisplay<NUM_PIXELS, RPM>;

/** @brief Brings up serial and debug telemetry, then starts the driver. */
void setup() {
  Serial.begin(9600);
  delay(1000);
  hs::configure_debug_telemetry();
  Serial.println("Hello");
  POV::begin();
}

/** @brief Shows RingSpin for 120 seconds. */
FLASHMEM static void run_show_sequence() {
  POV::show<RingSpin<CANVAS_W, CANVAS_H>>(120, false);
}

/** @brief Logs a heartbeat, then runs the show sequence. */
void loop() {
  Serial.println("Oh hi again");
  delay(1000);
  run_show_sequence();
}
