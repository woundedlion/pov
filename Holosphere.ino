#include <FastLED.h>
#include <fixmath.h>
#include "led.h"
#include "images.h"

namespace {
	POVDisplay<NUM_PIXELS, RPM> pov;
	const int K = 1000;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
}

void loop() {
//	pov.show<Test<40> >(9999 * K);
//  pov.show<RingTrails<96> >(90 * K);

	pov.show<TheMatrix<40, 135> >(60 * K);
	pov.show<ChainWiggle<96> >(120 * K);
	pov.show<RingRotate<96> >(120 * K);
	pov.show<RingShower<96> >(90 * K);
	pov.show<Curves<96> >(90 * K);
	pov.show<RingTwist<96> >(120 * K);
	pov.show<Kaleidoscope<96> >(120 * K);
	pov.show<StarsFade<40> >(60 * K);
	pov.show<RingTrails<96> >(90 * K);
	pov.show<DotTrails<96> >(90 * K);
	pov.show<Burnout<40, 0, 5> >(52 * K);
	pov.show<WaveTrails<96> >(60 * K);
	pov.show<Spinner<48> >(50 * K);
//	pov.show<Fire<40, 150, 120> >(45 * K);

	pov.show<Spiral<4, 0> >(10 * K);
	pov.show<Spiral<8, 0> >(10 * K);
	pov.show<Spiral<16, 0> >(10 * K);
	pov.show<Spiral<48, 0> >(30 * K);
	pov.show<Spiral<48, 1> >(30 * K);
	

	/*
	pov.show<TheMatrix<40, 135> >(15 * K);
	pov.show<DotTrails<96> >(30 * K);
	pov.show<ChainWiggle<96> >(30 * K);
	pov.show<RingRotate<96> >(30 * K);
	pov.show<RingShower<96> >(15 * K);
	pov.show<Curves<96> >(30 * K);
	pov.show<RingTwist<96> >(30 * K);
	pov.show<Kaleidoscope<96> >(30 * K);
	pov.show<StarsFade<40> >(15 * K);
	pov.show<RingTrails<96> >(15 * K);
//	pov.show<Burnout<40, 0, 5> >(52 * K);
	pov.show<WaveTrails<96> >(60 * K);
	pov.show<Spinner<48> >(50 * K);
	pov.show<Fire<40, 150, 120> >(15 * K);

	pov.show<Spiral<4, 0> >(2 * K);
	pov.show<Spiral<8, 0> >(2 * K);
	pov.show<Spiral<16, 0> >(2 * K);
	pov.show<Spiral<48, 0> >(5 * K);
	pov.show<Spiral<48, 1> >(5 * K);
	pov.show<Spiral<48, 2> >(2 * K);
	*/
}


