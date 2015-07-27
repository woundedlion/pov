#include <FastLED.h>
#include "led.h"
#include "images.h"

namespace {
  const uint32_t RPM = 240;
  const uint32_t NUM_PIXELS = 40;
 
  Spinner<40, 20> spinner;
  Grid<40, 20, 4> grid;
  
  const Effect* effects[] = {
    &spinner, 
    &grid,
  };
  
  CRGB leds[NUM_PIXELS];
  POVDisplay<NUM_PIXELS> pov(leds, RPM);
}

void setup() {
  Serial.begin(9600);
}

void loop() {
  unsigned long i = 0;
  while (true) { 
    pov.show(*effects[i], 5000);
    i = (i + 1) % (sizeof(effects) / sizeof(Effect *)); 
  }
}


