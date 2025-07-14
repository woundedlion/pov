#pragma once
#include <Arduino.h>
#include <FastLED.h>

namespace {
	const unsigned int RPM = 480;
	const int NUM_PIXELS = 40;
	const int H = NUM_PIXELS / 2;
	const int PIN_RANDOM = 15;
}

namespace hs {
	double rand_dbl() {
		return static_cast<double>(::random(0, std::numeric_limits<int32_t>::max()))
			/ std::numeric_limits<int32_t>::max();
	}

	int rand_int(int min, int max) {
		return ::random(min, max);
	}
}

#include "effects.h"
#include "effects2.h"


template <int S, int RPM>
class POVDisplay
{
  public:

	POVDisplay()
	{
     randomSeed(analogRead(PIN_RANDOM));
		 FastLED.addLeds<WS2801, 11, 13, RGB, DATA_RATE_MHZ(6)>(leds_, S);
	}

	template <typename E>
	void show(unsigned long duration)
	{
		long start = millis();
		effect_ = new E();
		x_ = 0;
		IntervalTimer timer;
		timer.begin(show_col, 1000000 / (RPM / 60) / effect_->width());
		while (millis() - start < duration * 1000) {
			effect_->draw_frame();
		}
		timer.end();
		delete effect_;
		effect_ = nullptr;
	}

  private:

    static inline void show_col() {

			for (int y = 0; y < S / 2; ++y) {
          leds_[S / 2 - y - 1] = effect_->get_pixel(x_, y);
          leds_[S / 2 + y] =  
            effect_->get_pixel((x_ + (effect_->width() / 2)) % effect_->width(), y);
      }

			FastLED.show();
			if (effect_->show_bg()) {
				FastLED.showColor(CRGB::Black);
			}
			
			x_ = (x_ + 1) % effect_->width();
			if (x_ == 0 || x_ == effect_->width() / 2) {
				effect_->advance_display();  
			}
		}

		static CRGB leds_[S];
		static Effect* effect_;
		static int x_;
};

template<int S, int RPM>
int POVDisplay<S, RPM>::x_ = 0;

template<int S, int RPM>
Effect* POVDisplay<S, RPM>::effect_ = nullptr;

template<int S, int RPM>
CRGB POVDisplay<S, RPM>::leds_[S];
