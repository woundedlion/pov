#include <Arduino.h>

typedef uint32_t Color;

class Effect
{
  public:
    virtual Color get_pixel(int x, int y) const = 0;
    virtual size_t width() const = 0;
    virtual size_t height() const = 0;
};

template <size_t W, size_t H>
class Image : public Effect
{
  public:
    Image(const unsigned char (*data)[H][3]) :
      data_(data)
    {}
    
    Color get_pixel(int x, int y) const {
      uint32_t c = pgm_read_dword(&data_[x][y][0]);
      return (c & 0x000000ff) << 16 | (c & 0x0000ff00) | (c & 0x00ff0000) >> 16;
    }
    
    size_t width() const { return W; }
    size_t height() const { return H; }

  private:
    const unsigned char (*data_)[H][3];  
};

template <size_t W, size_t H>
class Spinner : public Effect
{
  public:
    Spinner()
    {}
    
    Color get_pixel(int x, int y) const {
      unsigned long c = palette_[(x + y) % 4];
      return c;
    }
    
    size_t width() const { return W; }
    size_t height() const { return H; }

  private:
    
    const unsigned long palette_[4] = {0x00ff0000, 0x00ffff00, 0, 0};
};

template <size_t W, size_t H, size_t S>
class Grid : public Effect
{
  public:
    Grid()
    {}
    
    Color get_pixel(int x, int y) const {
      unsigned long c = 0;
      if (x % S == 0 || y % S == 0) {
        c = 0x0000ff00;
      } 
      return c;
    }
    
    size_t width() const { return W; }
    size_t height() const { return H; }
};

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

  void show(const Effect& effect, unsigned long duration)
  { 
    unsigned long show_start = millis();
    while (millis() - show_start < duration) {
      unsigned long elapsed_us = micros();
      for (int x = 0; x < effect.width(); ++x) {
        unsigned long col_us = micros();
        for (int y = 0; y < S / 2; ++y) {
          leds_[y] = effect.get_pixel(x, y);
          leds_[S - y - 1] =  
            effect.get_pixel((x + (effect.width() / 2)) % effect.width(), y);
        }
        FastLED.show();
        FastLED.showColor(CRGB::Black);
        unsigned long d_us = (1000000 / (effect.width() * static_cast<double>(rpm_ / 60))) - (micros() - col_us);
        unsigned long d_ms = d_us / 1000;
        delay(d_ms);
        delayMicroseconds(d_us % 1000);
      }
      elapsed_us = micros() - elapsed_us;
      Serial.println(elapsed_us);
    }    
  }

  private:
    CRGB *leds_;
    size_t rpm_;
};
