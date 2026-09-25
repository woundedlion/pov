/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Holosphere — Single-Teensy POV display
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

static constexpr int NUM_PIXELS = 40;
static constexpr unsigned int RPM = 480;

#ifdef USE_DMA_LEDS
// Explicit specialization keeps ledController's DMAMEM section (pov_single.h).
HS_DEFINE_POV_SINGLE_LED_CONTROLLER(NUM_PIXELS, RPM);
#endif

using POV = POVDisplay<NUM_PIXELS, RPM>;

void setup() {
  Serial.begin(9600);
  delay(1000);
  hs::configure_debug_telemetry();
  Serial.println("Hello");
  POV::begin();
}

FLASHMEM static void run_show_sequence() {
  POV::show<RingSpin<96, 20>>(120, false);
}

void loop() {
  Serial.println("Oh hi again");
  delay(1000);
  run_show_sequence();
}
