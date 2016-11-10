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
}

void loop() {
  while (true) { 
//    pov.show<StarsFade<40, 20> >(900000);

//    pov.show<ChasingDots<40, 20> >(300000);
//    pov.show<RotateWave<64, 20> >(900000);

    pov.show<TheMatrix<40, 20, 135> >(60000);
    pov.show<Kaleidoscope<100, 20> >(120000);
    pov.show<Rotate<40, 20> >(90000);
    pov.show<Spiral<4, 20, 0> >(10000);
    pov.show<Spiral<8, 20, 0> >(10000);
    pov.show<Spiral<16, 20, 0> >(10000);
    pov.show<Spiral<48, 20, 0> >(20000);
    pov.show<Spiral<48, 20, 1> >(20000);
    pov.show<Spiral<48, 20, 2> >(10000);
    pov.show<Spiral<48, 20, 3> >(10000);
    pov.show<Spinner<48, 20, 20> >(10000);
    pov.show<Spinner<48, 20, 10> >(5000);
    pov.show<Spinner<48, 20, 4> >(5000);
    pov.show<Spinner<48, 20, 2> >(5000);
    pov.show<Spinner<48, 20, 1> >(5000);
    pov.show<Spinner<48, 20, 1, true> >(5000);
    pov.show<Spinner<48, 20, 2, true> >(5000);
    pov.show<Spinner<48, 20, 4, true> >(5000);
    pov.show<Spinner<48, 20, 10, true> >(5000);
    pov.show<Spinner<48, 20, 20, true> >(5000);
    pov.show<Plaid<40, 20, 2> >(30000);
    pov.show<Snake<40, 20, 1> >(30000);
    pov.show<Stars<40, 20> >(30000);
    pov.show<Fire<40, 20, 150, 120> >(60000);
    pov.show<Burnout<40, 20, 0, 5> >(30000);
    pov.show<PaletteFall<40, 20, 5000, true> >(30000);
    pov.show<PaletteFall<1, 20, 5000, false> >(30000);
  }
}


