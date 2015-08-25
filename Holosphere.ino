#include <FastLED.h>
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
    pov.show<Grid<40, 20, 4> >(10000);

//    pov.show<PaletteFall<10, 20, 5000, true> >(30000);
//    pov.show<PaletteFall<20, 20, 5000, true> >(30000);
//    pov.show<PaletteFall<40, 20, 5000, true> >(30000);
//    pov.show<PaletteFall<1, 20, 5000, false> >(30000);

   pov.show<PaletteGrid<40, 20, 5000> >(30000);

    
//    pov.show<Water<40, 20, 1000> >(30000);

//    pov.show<Stars<40, 20> >(10000);
//    pov.show<Fire<40, 20, 150, 120> >(10000);
//  
//    pov.show<Spiral<4, 20, 0> >(30000);
//    pov.show<Spiral<8, 20, 0> >(30000);
//    pov.show<Spiral<16, 20, 0> >(30000);
//    pov.show<Spiral<48, 20, 0> >(30000);
//    pov.show<Spiral<48, 20, 1> >(30000);
//    pov.show<Spiral<48, 20, 2> >(30000);
//    pov.show<Spiral<48, 20, 3> >(30000);
//            
//    pov.show<Spinner<48, 20, 20> >(3000);
//    pov.show<Spinner<48, 20, 10> >(3000);
//    pov.show<Spinner<48, 20, 4> >(3000);
//    pov.show<Spinner<48, 20, 2> >(3000);
//    pov.show<Spinner<48, 20, 1> >(3000);
//    pov.show<Spinner<48, 20, 1, true> >(3000);
//    pov.show<Spinner<48, 20, 2, true> >(3000);
//    pov.show<Spinner<48, 20, 4, true> >(3000);
//    pov.show<Spinner<48, 20, 10, true> >(3000);
//    pov.show<Spinner<48, 20, 20, true> >(3000);

//    pov.show<Image<36, 20, world> >(10000);
   
  }
}


