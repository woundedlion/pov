#include <FastLED.h>
#include <SPI.h>
#include <vector>
#include "led.h"
#include "effects_legacy.h"
#include "effects.h"

// TODO: Convert from float to float everywhere
// TODO: Compile-tie FilterChain
// TODO: FilterMask
// TODO: FilterDisplace

namespace {
	POVDisplay<NUM_PIXELS, RPM> *pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
	pov = new POVDisplay<NUM_PIXELS, RPM>();

}

void loop() {
//	pov->show<Dynamo<96>>(300);
// pov->show<RingSpin<96>>(300);
// pov->show<Comets<96>>(300);
//	pov->show<FlowField<96>>(300);

	{
		//	pov->show<RingShower<96>>(90);
			pov->show<Dynamo<96>>(300);

	}

/*
	pov->show<TheMatrix<40, 135>>(60);
	pov->show<ChainWiggle<96>>(150);
	pov->show<RingRotate<96>>(120);
	pov->show<RingShower<96>>(90);
	pov->show<Curves<96> >(120);
	pov->show<RingTwist<96>>(150);
	pov->show<Kaleidoscope<96>>(120);
	pov->show<StarsFade<40>>(60); 
	pov->show<RingTrails<96>>(90); // FIXME OR REPLACE WITH RingSpin?
	pov->show<DotTrails<96>>(120);
	pov->show<Burnout<40, 0, 5>>(52);
	pov->show<WaveTrails<96>>(60);
	pov->show<Spinner<48>>(50);

	pov->show<Spiral<4, 0> >(10);
	pov->show<Spiral<8, 0> >(10);
	pov->show<Spiral<16, 0> >(10);
	pov->show<Spiral<48, 0> >(30);
	pov->show<Spiral<48, 1> >(20);
*/

	//	pov->show<PolyRot<96>>(45);
	//	pov->show<Fire<40, 150, 120> >(45);
}


