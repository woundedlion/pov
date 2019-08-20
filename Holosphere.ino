#include <FastLED.h>
#include <fixmath.h>
#include "led.h"
#include "images.h"

namespace {
  const unsigned int RPM = 480;
  const uint8_t NUM_PIXELS = 40;
  
  CRGB leds[NUM_PIXELS];
  POVDisplay<NUM_PIXELS> pov(leds, RPM);
}

void setup() {
  Serial.begin(9600);
  delay(1000);
}

void loop() {
	while (true) { 
//    pov.show<RotateWave<64, 20> >(900000);
//    pov.show<Snake<40, 20, 1> >(30000);
	
		pov.show<TheMatrix<40, 20, 135> >(90000);
		pov.show<DotTrails<96, 20> >(90000);
		pov.show<RingTwist<80, 20> >(90000);
		pov.show<Curves<96, 20> >(90000);
		pov.show<Kaleidoscope<100, 20> >(90000);
		pov.show<Rotate<40, 20> >(90000);
		pov.show<Spinner<48, 20> >(50000);
		pov.show<RingTrails<80, 20> >(90000);
		pov.show<WaveTrails<96, 20> >(90000);

		pov.show<Spiral<4, 20, 0> >(10000);
		pov.show<Spiral<8, 20, 0> >(10000);
		pov.show<Spiral<16, 20, 0> >(10000);
		pov.show<Spiral<48, 20, 0> >(20000);
		pov.show<Spiral<48, 20, 1> >(20000);
		pov.show<Spiral<48, 20, 2> >(10000);
		pov.show<Spiral<48, 20, 3> >(10000);
		
		pov.show<StarsFade<40, 20> >(60000);
		pov.show<Burnout<40, 20, 0, 5> >(30000);
		pov.show<Fire<40, 20, 150, 120> >(60000);
	}
}


