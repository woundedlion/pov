#include <FastLED.h>
#include <SPI.h>
#include "led.h"

#include <vector>

namespace {
	POVDisplay<NUM_PIXELS, RPM> *pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
	pov = new POVDisplay<NUM_PIXELS, RPM>();

}

void loop() {
	pov->show<TheMatrix<40, 135> >(60);
	pov->show<ChainWiggle<96> >(150);
	pov->show<RingRotate<96> >(120);
	pov->show<RingShower<96> >(90);
	pov->show<Curves<96> >(120);
	pov->show<RingTwist<96> >(150);
	pov->show<Kaleidoscope<96> >(120);
	pov->show<StarsFade<40> >(60);
	pov->show<RingTrails<96> >(90);
	pov->show<DotTrails<96> >(120);
	pov->show<Burnout<40, 0, 5> >(52);
	pov->show<WaveTrails<96> >(60);
	pov->show<Spinner<48> >(50);

	pov->show<Spiral<4, 0> >(10);
	pov->show<Spiral<8, 0> >(10);
	pov->show<Spiral<16, 0> >(10);
	pov->show<Spiral<48, 0> >(30);
	pov->show<Spiral<48, 1> >(20);
	pov->show<PolyRot<96>>(45);


	//	pov->show<Fire<40, 150, 120> >(45);


	/*
	pov->show<TheMatrix<40, 135> >(15);
	pov->show<DotTrails<96> >(30);
	pov->show<ChainWiggle<96> >(30);
	pov->show<RingRotate<96> >(30);
	pov->show<RingShower<96> >(15);
	pov->show<Curves<96> >(30);
	pov->show<RingTwist<96> >(30);
	pov->show<Kaleidoscope<96> >(30);
	pov->show<StarsFade<40> >(15);
	pov->show<RingTrails<96> >(15);
//	pov->show<Burnout<40, 0, 5> >(52);
	pov->show<WaveTrails<96> >(60);
	pov->show<Spinner<48> >(50);
	pov->show<Fire<40, 150, 120> >(15);

	pov->show<Spiral<4, 0> >(2);
	pov->show<Spiral<8, 0> >(2);
	pov->show<Spiral<16, 0> >(2);
	pov->show<Spiral<48, 0> >(5);
	pov->show<Spiral<48, 1> >(5);
	pov->show<Spiral<48, 2> >(2);
	*/
}


