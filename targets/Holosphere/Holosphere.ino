/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Holosphere — Single-Teensy POV display (96×20 physical strip)
 *
 * Target: Teensy 4.0
 * Physical LEDs: 40 (20 per arm × 2 sides)
 * Virtual canvas: 288×144 (or 96×20 for legacy effects)
 */
#include <FastLED.h>
#include <SPI.h>

#include "pov_single.h"
#include "effects.h"
#include "effects_legacy.h"

static constexpr int NUM_PIXELS = 40;
static constexpr unsigned int RPM = 480;

namespace {
POVDisplay<NUM_PIXELS, RPM> *pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
  pov = new POVDisplay<NUM_PIXELS, RPM>();
}

FLASHMEM static void run_show_sequence() {
  // Test/Debug Effects
  pov->show<Test<288, 144>>(120);
//  pov->show<TestShapes<288, 144>>(120);

  pov->show<ChaoticStrings<288, 144>>(120);
  pov->show<Liquid2D<288, 144>>(120);
  pov->show<Flyby<288, 144>>(120);
  pov->show<Raymarch<288, 144>>(120);
//  pov->show<SplineFlow<288, 144>>(120);
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
//  pov->show<Thrusters<288, 144>>(120);
}

void loop() { run_show_sequence(); }
