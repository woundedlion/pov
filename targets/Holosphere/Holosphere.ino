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

static constexpr int NUM_PIXELS = 288;
static constexpr unsigned int RPM = 480;

namespace {
POVDisplay<NUM_PIXELS, RPM> *pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
  Serial.println("Hello");
  pov = new POVDisplay<NUM_PIXELS, RPM>();
}

FLASHMEM static void run_show_sequence() {
  pov->show<HankinSolids<288, 144>>(120); // 145

  /*
  pov->show<RingSpin<288, 144>>(120); // 140
  pov->show<ChaoticStrings<288, 144>>(120);
  pov->show<Raymarch<288, 144>>(120);
  pov->show<SplineFlow<288, 144>>(120);
  pov->show<GnomonicStars<288, 144>>(120); //60
  pov->show<DreamBalls<288, 144>>(120); // 117
  pov->show<HankinSolids<288, 144>>(120);
  pov->show<IslamicStars<288, 144>>(120);
  pov->show<MindSplatter<288, 144>>(120);
  pov->show<Voronoi<288, 144>>(120);
  pov->show<MobiusGrid<288, 144>>(120);
  pov->show<FlowField<288, 144>>(120);
  pov->show<MeshFeedback<288, 144>>(120);
  pov->show<PetalFlow<288, 144>>(120); //35
  pov->show<BZReactionDiffusion<288, 144>>(120);
  pov->show<GSReactionDiffusion<288, 144>>(120);
  pov->show<Thrusters<288, 144>>(120);
  pov->show<TestShapes<288, 144>>(120); //140

  pov->show<Liquid2D<288, 144>>(120);  //85
  pov->show<Flyby<288, 144>>(120); // 85

  pov->show<Test<288, 144>>(120); // 60
  pov->show<SphericalHarmonics<288, 144>>(120);  //49
  pov->show<Comets<288, 144>>(120); //53

  pov->show<HopfFibration<288, 144>>(120); //72


  */
}

void loop() {
  Serial.println("Oh hi again");
  delay(1000);
  run_show_sequence();
}
