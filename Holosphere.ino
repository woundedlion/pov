/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#include <FastLED.h>
#include <SPI.h>

#include "led.h"
#include "effects.h"
#include "effects_legacy.h"

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
  pov->show<TestShapes<288, 144>>(120);

  pov->show<MeshFeedback<288, 144>>(120);
  pov->show<ChaoticStrings<288, 144>>(120);
  pov->show<Liquid2D<288, 144>>(120);
  pov->show<Flyby<288, 144>>(120);
  pov->show<Raymarch<288, 144>>(120);
  pov->show<SplineFlow<288, 144>>(120);
  pov->show<GnomonicStars<288, 144>>(120);

  // Core Effects
  pov->show<DreamBalls<288, 144>>(120);
  pov->show<SpinShapes<288, 144>>(120);
  pov->show<HopfFibration<288, 144>>(120);
  pov->show<HankinSolids<288, 144>>(120);
  pov->show<IslamicStars<288, 144>>(120);
  pov->show<MindSplatter<288, 144>>(120);
  pov->show<SphericalHarmonics<288, 144>>(120);
  pov->show<Voronoi<288, 144>>(120);
  pov->show<MobiusGrid<288, 144>>(120);
  pov->show<FlowField<288, 144>>(120);

  // Existing Effects
  pov->show<PetalFlow<288, 144>>(120);
  pov->show<BZReactionDiffusion<288, 144>>(120);
  pov->show<GSReactionDiffusion<288, 144>>(120);
  pov->show<Thrusters<288, 144>>(120);

  /*

pov->show<TheMatrix<40, 20, 135>>(45);
pov->show<ChainWiggle<96, 20>>(150);
pov->show<RingRotate<96, 20>>(120);
pov->show<RingShower<96, 20>>(120);
pov->show<Curves<96, 20>>(120);
pov->show<RingTwist<96, 20>>(150);
pov->show<Dynamo<96, 20>>(120);
pov->show<RingSpin<96, 20>>(120);
pov->show<Comets<96, 20>>(120);
pov->show<FlowField<96, 20>>(60);
pov->show<Kaleidoscope<96, 20>>(90);
pov->show<StarsFade<40, 20>>(60);
pov->show<DotTrails<96, 20>>(120);
pov->show<Burnout<40, 20, 0, 5>>(52);
pov->show<Spinner<48, 20>>(50);
pov->show<Spiral<4, 20, 0>>(10);
pov->show<Spiral<8, 20, 0>>(10);
pov->show<Spiral<16, 20, 0>>(10);
pov->show<Spiral<48, 20, 0>>(20);
pov->show<Spiral<48, 20, 1>>(10);

*/
}

void loop() { run_show_sequence(); }
