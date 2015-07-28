#include <FastLED.h>
#include "led.h"
#include "images.h"

namespace {
  const uint32_t RPM = 240;
  const uint32_t NUM_PIXELS = 40;
 
  Spiral<60, 20> spiral;
  Grid<40, 20, 4> grid;
  Stars<40, 20> stars;
  TheMatrix<40, 20> the_matrix;
  
  Effect* effects[] = {
    &grid,
    &spiral, 
    &stars,
//      &the_matrix,
  };
  
  CRGB leds[NUM_PIXELS];
  POVDisplay<NUM_PIXELS> pov(leds, RPM);
}

void setup() {
  Serial.begin(9600);
}

void loop() {
  unsigned int i = 0;
  while (true) { 
    pov.show(*effects[i], 5000);
    i = addmod8(i, 1, sizeof(effects) / sizeof(Effect *));
  }
}


