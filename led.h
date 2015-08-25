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
    unsigned long col_delay_us = 1000000 / (rpm_ / 60) / effect.width();
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
        col_us = micros() - col_us;
        if (col_us < col_delay_us) {
          delay((col_delay_us - col_us) / 1000);
          delayMicroseconds((col_delay_us - col_us) % 1000);
        }
      }
      effect.advance_frame();
      elapsed_us = micros() - elapsed_us;
      Serial.println(elapsed_us);
    }    
  }

  private:
    CRGB *leds_;
    size_t rpm_;
};
