/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Phantasm — Multi-Teensy segmented POV display (288×144)
 *
 * Target: 4× Teensy 4.1
 * Total physical LEDs: 288
 * LEDs per segment: 72
 * Virtual canvas: 288×144
 *
 * Each Teensy reads its hardware ID at boot to determine which
 * segment of the LED strip it owns.
 *
 * TODO: Implement after pov_segmented.h driver is complete.
 */

// #include "effects.h"
// #include "../../src/hardware/pov_segmented.h"

static constexpr int TOTAL_PIXELS = 288;
static constexpr int NUM_SEGMENTS = 4;
static constexpr int PIXELS_PER_SEGMENT = TOTAL_PIXELS / NUM_SEGMENTS;
static constexpr unsigned int RPM = 1200;

void setup() {
  // TODO: Read hardware ID, initialize segmented POV driver
}

void loop() {
  // TODO: Run show sequence with segmented driver
}
