#include <FastLED.h>
#include <fixmath.h>
#include "led.h"
#include "images.h"

namespace {
	const unsigned int RPM = 480;
	const uint8_t NUM_PIXELS = 40;
	POVDisplay<NUM_PIXELS, RPM> pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
}

const int K = 1000;
void loop() {
	while (true) { 
//    pov.show<RotateWave<64, 20> >(900000);
//    pov.show<Snake<40, 20, 1> >(30000);

		pov.show<RingTrails<80, 20> >(900 * K);

		pov.show<TheMatrix<40, 20, 135> >(90 * K);
		pov.show<DotTrails<96, 20> >(120 * K);
		pov.show<RingTwist<96, 20> >(120 * K);
		pov.show<Curves<96, 20> >(90 * K);
		pov.show<Kaleidoscope<100, 20> >(120 * K);
		pov.show<Rotate<40, 20> >(90 * K);
		pov.show<Spinner<48, 20> >(50 * K);
		pov.show<RingTrails<80, 20> >(90 * K);
		pov.show<WaveTrails<96, 20> >(120 * K);

		pov.show<Spiral<4, 20, 0> >(10 * K);
		pov.show<Spiral<8, 20, 0> >(10 * K);
		pov.show<Spiral<16, 20, 0> >(10 * K);
		pov.show<Spiral<48, 20, 0> >(60 * K);
		pov.show<Spiral<48, 20, 1> >(30 * K);
		pov.show<Spiral<48, 20, 2> >(15 * K);
		
		pov.show<StarsFade<40, 20> >(60 * K);
		pov.show<Burnout<40, 20, 0, 5> >(30 * K);
		pov.show<Fire<40, 20, 150, 120> >(90 * K);
	}
}


