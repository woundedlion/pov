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
//    pov.show<Grid<40, 20, 4> >(10000);
//    pov.show<Spiral<60, 20> >(10000);
//    pov.show<Stars<40, 20> >(10000);
//    pov.show<TheMatrix<40, 20, 5, RPM> >(10000);
//    pov.show<Spinner<40, 20, 16> >(10000);
    pov.show<Fire<40, 20, 100, 120> >(10000);
  }
}


