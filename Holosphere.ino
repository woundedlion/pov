#include <FastLED.h>
#include "led.h"
#include "images.h"

namespace {
  const unsigned int RPM = 240;
  const uint8_t NUM_PIXELS = 40;
  
  CRGB leds[NUM_PIXELS];
  POVDisplay<NUM_PIXELS> pov(leds, RPM);
}

void setup() {
  Serial.begin(9600);
}

void loop() {
  while (true) { 
//    pov.show<Grid<40, 20, 4> >(5000);
//    pov.show<Image<26, 20, world> >(10000);
//    pov.show<Spiral<60, 20> >(10000);
    pov.show<Stars<40, 20> >(10000);
//    pov.show<TheMatrix<40, 20, 5, RPM> >(10000);

    pov.show<Spinner<48, 20, 20> >(3000);
    pov.show<Spinner<48, 20, 10> >(3000);
    pov.show<Spinner<48, 20, 4> >(3000);
    pov.show<Spinner<48, 20, 2> >(3000);
    pov.show<Spinner<48, 20, 1> >(3000);
    pov.show<Spinner<48, 20, 1, true> >(3000);
    pov.show<Spinner<48, 20, 2, true> >(3000);
    pov.show<Spinner<48, 20, 4, true> >(3000);
    pov.show<Spinner<48, 20, 10, true> >(3000);
    pov.show<Spinner<48, 20, 20, true> >(3000);

    pov.show<Fire<40, 20, 150, 120> >(10000);
//    pov.show<WaveBall<40, 20, 8, 4> >(10000);
  }
}


