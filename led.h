#include <Arduino.h>
#include "effects.h"

template <size_t S>
class POVDisplay
{
  public:
    POVDisplay(CRGB *leds, size_t rpm) :
    leds_(leds),
    rpm_(rpm)
    {
      FastLED.addLeds<WS2801, 11, 13, RGB, DATA_RATE_MHZ(4)>(leds, S);
    }

  template <typename Effect>
  void show(unsigned long duration)
  { 
    unsigned long show_start = millis();
    Effect effect;
    while (millis() - show_start < duration) {
      unsigned long elapsed_us = micros();
      for (int x = 0; x < effect.width(); ++x) {
        unsigned long col_us = micros();
        for (int y = 0; y < S / 2; ++y) {
          leds_[S / 2 - y - 1] = effect.get_pixel(x, y);
          leds_[S/2 + y] =  
            effect.get_pixel((x + (effect.width() / 2)) % effect.width(), y);
        }
        FastLED.show();
        if (effect.show_bg()) {
          FastLED.showColor(effect.bg_color());
        }
        effect.advance_col(x);
        unsigned long d_us = max(0, (1000000 / (effect.width() * static_cast<double>(rpm_ / 60))) - (micros() - col_us));
        delay(d_us / 1000);
        delayMicroseconds(d_us % 1000);
      }
      elapsed_us = micros() - elapsed_us;
//      Serial.println(elapsed_us);
    }    
  }

  private:
    CRGB *leds_;
    size_t rpm_;
};
