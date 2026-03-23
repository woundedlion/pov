/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Phantasm — Multi-Teensy segmented POV display (288×144)
 *
 * Target: 4× Teensy 4.1
 * Total physical LEDs: 288 (72 per segment, 144 per arm)
 * Virtual canvas: 288×144
 *
 * Each Teensy reads its hardware ID at boot (pins 21, 22) to determine
 * which segment of the LED strip it owns.  Segment 0 generates the
 * column clock; all boards fire their ISR on the shared clock line.
 *
 * Hardware ID assignment (active-low, ground to set):
 *   ID 0: both floating  — arm A outer (clock master)
 *   ID 1: pin 21 grounded — arm A inner
 *   ID 2: pin 22 grounded — arm B outer
 *   ID 3: both grounded   — arm B inner
 */
#include <FastLED.h>
#include <SPI.h>

#include "pov_segmented.h"
#include "effects.h"

static constexpr int TOTAL_PIXELS = 288;
static constexpr int NUM_SEGMENTS = 4;
static constexpr unsigned int RPM = 1200;

namespace {
POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM> *pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
  pov = new POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM>();
}

FLASHMEM static void run_show_sequence() {
  pov->show<Test<288, 144>>(120);
  pov->show<ChaoticStrings<288, 144>>(120);
  pov->show<Liquid2D<288, 144>>(120);
  pov->show<Flyby<288, 144>>(120);
  pov->show<Raymarch<288, 144>>(120);
  pov->show<GnomonicStars<288, 144>>(120);
  pov->show<DreamBalls<288, 144>>(120);
  pov->show<HopfFibration<288, 144>>(120);
  pov->show<HankinSolids<288, 144>>(120);
  pov->show<IslamicStars<288, 144>>(120);
  pov->show<MindSplatter<288, 144>>(120);
  pov->show<SphericalHarmonics<288, 144>>(120);
  pov->show<Voronoi<288, 144>>(120);
  pov->show<MobiusGrid<288, 144>>(120);
  pov->show<FlowField<288, 144>>(120);
  pov->show<MeshFeedback<288, 144>>(120);
  pov->show<PetalFlow<288, 144>>(120);
  pov->show<BZReactionDiffusion<288, 144>>(120);
  pov->show<GSReactionDiffusion<288, 144>>(120);
}

void loop() { run_show_sequence(); }
