#include <FastLED.h>

class Effect
{
  public:
    virtual CRGB get_pixel(int x, int y) const = 0;
    virtual void reset() {}
    virtual void advance_frame() {};
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
    
    CRGB get_pixel(int x, int y) const {
      uint32_t c = pgm_read_dword(&data_[x][y][0]);
      return CRGB((c & 0xff) << 16, c & 0xff00, (c & 0xff0000) >> 16);
    }
    
    size_t width() const { return W; }
    size_t height() const { return H; }

  private:
    const unsigned char (*data_)[H][3];  
};

template <size_t W, size_t H, size_t S>
class Grid : public Effect
{
  public:
    Grid()
    {}
    
    CRGB get_pixel(int x, int y) const {
      unsigned long c = 0;
      if (x % S == 0 || y % S == 0) {
        c = CRGB::Green;
      } 
      return c;
    }
    
    size_t width() const { return W; }
    size_t height() const { return H; }
};

template <size_t W, size_t H>
class Spiral : public Effect
{
  public:
    Spiral()
    {
      fill_rainbow(palette_.entries, 16, 0, 256/16);
    }
    
    CRGB get_pixel(int x, int y) const {
      return palette_[(x + y) % 15];
    }
        
    size_t width() const { return W; }
    size_t height() const { return H; }

  private:
        
    CRGBPalette16 palette_;
};

template <size_t W, size_t H>
class Stars : public Effect
{
  public:
    Stars()
    {
    }
    
    CRGB get_pixel(int x, int y) const {
      return random8() > 250 ? CRGB::Goldenrod : CRGB::Black;
    }
        
    size_t width() const { return W; }
    size_t height() const { return H; }

  private:
};

template <size_t W, size_t H>
class TheMatrix : public Effect
{
  public:
    TheMatrix() :
      drop_timer_(10)
    {
      random16_add_entropy(micros());
      reset();
    }
    
    CRGB get_pixel(int x, int y) const {
      uint8_t dist = static_cast<int>(pos_[x]) - y;
      if (dist >= 0 && dist <= 5) {
        // trails fade to black
        return blend(CRGB::Green, CRGB::Black, 255 * dist / 5);
      } else {
        return CRGB::Black; 
      }
    }
   
    void reset() {
      memset8(&pos_[0], 0, sizeof(pos_));
      for (int i = 0; i < W; ++i) {
         // random fall rate from 0.1 to 1.0
         rate_[i] = static_cast<float>(random16(6554, 65535)) / 65536;
      }
    }
    
    void advance_frame() {
      if (drop_timer_) {
        fall();
      }
    } 
           
    size_t width() const { return W; }
    size_t height() const { return H; }

  private:
   
    void fall() {
      for (int i = 0; i < W; ++i) {
        pos_[i] += rate_[i];
      }  
    }
    
    CEveryNMillis drop_timer_; 
    float pos_[W];
    float rate_[W];
};

