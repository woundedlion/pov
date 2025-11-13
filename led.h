#pragma once
#include <Arduino.h>
#include <FastLED.h>

static constexpr unsigned int RPM = 480;
static constexpr int NUM_PIXELS = 40;
static constexpr int H = NUM_PIXELS / 2;
static constexpr int H_VIRT = H + 3;
static constexpr int MAX_W = 96;
static constexpr int PIN_RANDOM = 15;

inline int XY(int x, int y) { return x * H + y; }
class Canvas;

class Effect {
  friend class Canvas;

public:
  Effect(int W) :
    persist_pixels(true),
    width_(W)
  {
    bufs_[0] = buffer_a;
    memset(bufs_[0], 0, sizeof(CRGB) * W * H);
    bufs_[1] = buffer_b;
    memset(bufs_[1], 0, sizeof(CRGB) * W * H);
  }

  virtual ~Effect() {
  };

  virtual void draw_frame() = 0;
  virtual bool show_bg() const = 0;

  virtual const CRGB& get_pixel(int x, int y) const {
    return bufs_[prev_][XY(x, y)];
  }

  inline int width() const { return width_; }
  inline bool buffer_free() const { return prev_ == next_; }
  inline void advance_display() { prev_ = next_; }
  inline void advance_buffer() {
    noInterrupts();
    cur_ = cur_ ? 0 : 1;
    interrupts();
    if (persist_pixels) {
      memcpy(bufs_[cur_], bufs_[prev_], sizeof(CRGB) * width_ * H);
    }
  }

  inline void queue_frame() {
    noInterrupts();
    next_ = cur_;
    interrupts();
  }

protected:
  bool persist_pixels;

private:
  volatile int prev_ = 0, cur_ = 0, next_ = 0;
  int width_;
  inline static CRGB buffer_a[MAX_W * H];
  inline static CRGB buffer_b[MAX_W * H];
  CRGB* bufs_[2];
};

class Canvas {
public:
  Canvas(Effect& effect) : effect_(effect) {
    while (!effect_.buffer_free()) {}
    effect_.advance_buffer();
    if (!effect_.persist_pixels) {
      clear_buffer();
    }
  }

  ~Canvas() { effect_.queue_frame(); }

  inline CRGB& operator()(int x, int y) {
    return effect_.bufs_[effect_.cur_][XY(x, y)];
  }

  inline CRGB& operator()(int xy) {
    return effect_.bufs_[effect_.cur_][xy];
  }

  const int width() { return effect_.width(); }

  void clear_buffer() {
    memset(effect_.bufs_[effect_.cur_], 0, sizeof(CRGB) * effect_.width_ * H);
  }

private:
  Effect& effect_;
};

template <int S, int RPM>
class POVDisplay
{
  public:

	POVDisplay()
	{
     randomSeed(analogRead(PIN_RANDOM));
		 FastLED.addLeds<WS2801, 11, 13, RGB, DATA_RATE_MHZ(6)>(leds_, S);
     FastLED.setCorrection(TypicalLEDStrip);
     FastLED.setTemperature(Candle);
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
