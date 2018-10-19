#include <Arduino.h>
#include "effects.h"

template <int S>
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
      unsigned long frame_start = micros();
      for (int x = 0; x < effect.width() / 2; ++x) {
        show_col(effect, x, col_delay_us);
      }
      effect.advance_frame();
      for (int x = effect.width() / 2; x < effect.width(); ++x) {
        show_col(effect, x, col_delay_us);
      }
      effect.advance_frame();
      Serial.println(micros() - frame_start);
    }    
  }

  private:

    template <typename Effect>
    inline void show_col(Effect& effect, int x, unsigned long col_delay_us) {
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
          noInterrupts();
          delay((col_delay_us - col_us) / 1000);
          delayMicroseconds((col_delay_us - col_us) % 1000);
          interrupts();
        }
    }
  
    CRGB *leds_;
    size_t rpm_;
};
