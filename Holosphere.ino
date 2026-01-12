/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#include <FastLED.h>
#include <SPI.h>
#include <vector>
#include "led.h"
#include "effects_legacy.h"
#include "effects.h"

// TODO: FilterMask
// TODO: FilterDisplace

namespace {
  POVDisplay<NUM_PIXELS, RPM>* pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
  pov = new POVDisplay<NUM_PIXELS, RPM>();

}

void loop() {
  pov->show<Test<96>>(120);
  pov->show<TestScanPolygon<96>>(120);
  pov->show<TestPlotPolygon<96>>(120);
  pov->show<PetalFlow<96>>(120);
  pov->show<BZReactionDiffusion<96>>(120);
  pov->show<GSReactionDiffusion<96>>(120);
  pov->show<Thrusters<96>>(120);
  pov->show<TheMatrix<40, 135>>(45);
  pov->show<ChainWiggle<96>>(150);
  pov->show<RingRotate<96>>(120);
  pov->show<RingShower<96>>(120);
  pov->show<Curves<96> >(120);
  pov->show<RingTwist<96>>(150);
  pov->show<Dynamo<96>>(120);
  pov->show<RingSpin<96>>(120);
  pov->show<Comets<96>>(120);
  pov->show<FlowField<96>>(60);
  pov->show<Kaleidoscope<96>>(90);
  pov->show<StarsFade<40>>(60);
  pov->show<DotTrails<96>>(120);
  pov->show<Burnout<40, 0, 5>>(52);
  pov->show<Spinner<48>>(50);
  pov->show<Spiral<4, 0> >(10);
  pov->show<Spiral<8, 0> >(10);
  pov->show<Spiral<16, 0> >(10);
  pov->show<Spiral<48, 0> >(20);
  pov->show<Spiral<48, 1> >(10);
}


