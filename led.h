#include <Arduino.h>
#include "effects.h"

template <int S, int RPM>
class POVDisplay
{
  public:

	POVDisplay()
	{
		FastLED.addLeds<WS2801, 11, 13, RGB, DATA_RATE_MHZ(4)>(leds_, S);
	}

	template <typename E>
	void show(unsigned long duration)
	{
		effect_ = new E();
		x_ = 0;
		IntervalTimer timer;
		timer.begin(show_col, 1000000 / (RPM / 60) / effect_->width());
		delay(duration);
		timer.end();
		delete effect_;
		effect_ = NULL;
	}

  private:

    static inline void show_col() {
        for (int y = 0; y < S / 2; ++y) {
          leds_[S / 2 - y - 1] = effect_->get_pixel(x_, y);
          leds_[S/2 + y] =  
            effect_->get_pixel((x_ + (effect_->width() / 2)) % effect_->width(), y);
        }
        FastLED.show();
        if (effect_->show_bg()) {
          FastLED.showColor(effect_->bg_color());
        }
        effect_->advance_col(x_++);
		if (x_ == effect_->width() / 2) {
			effect_->advance_frame();
		}
		else if (x_ == effect_->width()) {
			x_ = 0;
			effect_->advance_frame();
		}
    }
  
    static CRGB leds_[S];
	static Effect* effect_;
	static int x_;
};

template<int NUM_PIXELS, int RPM>
int POVDisplay<NUM_PIXELS, RPM>::x_ = 0;

template<int NUM_PIXELS, int RPM>
Effect* POVDisplay<NUM_PIXELS, RPM>::effect_ = NULL;

template<int NUM_PIXELS, int RPM>
CRGB POVDisplay<NUM_PIXELS, RPM>::leds_[NUM_PIXELS];
