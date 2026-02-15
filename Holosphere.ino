/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#include <FastLED.h>
#include <SPI.h>
#include <vector>
#include "led.h"
#include "effects.h"
#include "effects_legacy.h"

namespace {
  POVDisplay<NUM_PIXELS, RPM>* pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
  pov = new POVDisplay<NUM_PIXELS, RPM>();

}

void loop() {
  // Test/Debug Effects
  pov->show<Test<96, 20>>(120);
  pov->show<TestShapes<96, 20>>(120);
  pov->show<TestTemporal<96, 20>>(120); // NEW
  pov->show<FlamingMesh<96, 20>>(120); // NEW

  // Core Effects
  pov->show<DreamBalls<96, 20>>(120); // NEW
  pov->show<SpinShapes<96, 20>>(120); // NEW
  pov->show<HopfFibration<96, 20>>(120); // NEW
  pov->show<HankinSolids<96, 20>>(120); // NEW
  pov->show<IslamicStars<96, 20>>(120); // NEW
  pov->show<LSystem<96, 20>>(120); // NEW
  pov->show<MetaballEffect<96, 20>>(120); // NEW
  pov->show<MindSplatter<96, 20>>(120); // NEW
  pov->show<SphericalHarmonics<96, 20>>(120); // NEW
  pov->show<Voronoi<96, 20>>(120); // NEW
  
  // Existing Effects
  pov->show<PetalFlow<96, 20>>(120);
  pov->show<BZReactionDiffusion<96, 20>>(120);
  pov->show<GSReactionDiffusion<96, 20>>(120);
  pov->show<Thrusters<96, 20>>(120);
  pov->show<TheMatrix<40, 20, 135>>(45);
  pov->show<ChainWiggle<96, 20>>(150);
  pov->show<RingRotate<96, 20>>(120);
  pov->show<RingShower<96, 20>>(120);
  pov->show<Curves<96, 20> >(120);
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
  pov->show<Spiral<4, 20, 0> >(10);
  pov->show<Spiral<8, 20, 0> >(10);
  pov->show<Spiral<16, 20, 0> >(10);
  pov->show<Spiral<48, 20, 0> >(20);
  pov->show<Spiral<48, 20, 1> >(10);
}
