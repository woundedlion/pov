#include <FastLED.h>

template <uint8_t W, uint8_t H, const unsigned char (*DATA)[20][3]>
class Image
{
  public:
    
    Image()
    {}
    
    CRGB get_pixel(int x, int y) const {
      uint32_t c = pgm_read_dword(&DATA[x][y][0]);
      return CRGB(c & 0xff, (c & 0xff00) >> 8, (c & 0xff0000) >> 16);
    }

    static bool show_bg() { return false; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {} 

    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

  private:

    const unsigned char (*data_)[H][3];  
};


template <uint8_t W, uint8_t H, uint8_t S>
class Grid
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
        
    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {} 
    
    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

};

template <uint8_t W, uint8_t H>
class Spiral
{
  public:
    Spiral()
    {
      fill_rainbow(palette_.entries, 16, 0, 256/16);
    }
    
    CRGB get_pixel(int x, int y) const {
      return palette_[(x + y) % 15];
    }
        
    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    
   
    void advance_col(uint8_t x) {} 
    
    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

  private:
        
    CRGBPalette16 palette_;
};

template <uint8_t W, uint8_t H>
class Stars
{
  public:
    Stars()
    {}
    
    CRGB get_pixel(int x, int y) const {
      return random8() > 250 ? CRGB::Goldenrod : CRGB::Black;
    }
        
    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {} 
    
    uint8_t width() const { return W; }
    uint8_t height() const { return H; }
};

template <uint8_t W, uint8_t H, uint8_t T, uint8_t RPM>
class TheMatrix
{
  public:
    TheMatrix()
    {
      random16_add_entropy(random());
      c_.setHue(random8());
      memset8(&pos_[0], 0, sizeof(pos_));
      for (int i = 0; i < W; ++i) {
         // random fall rate in pix/sec
         rate_[i] = random8(1,5);
      }
    }
    
    CRGB get_pixel(int x, int y) const {
      uint8_t dist = static_cast<uint8_t>(pos_[x]) - y;
      if (dist >= 0 && dist <= T) {
        // trails fade to black
        return CRGB(c_).fadeLightBy(255 * dist / T);
      } else {
        return CRGB::Black; 
      }
    }
               
    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {
      fall(x);    
    } 
    
    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

  private:
   
    inline void fall(uint8_t x) {
      pos_[x] += static_cast<float>(rate_[x]) * 60 / RPM;
    }
    
    float pos_[W];
    uint8_t rate_[W];
    CRGB c_;
};

template <uint8_t W, uint8_t H, uint8_t P, bool swap = false>
class Spinner
{
  public:
    Spinner() :
      spin_timer_(125),
      pos_(0)
    {
      for (int i = 0; i < 9; i += 8) {
        palette1_[i] = CRGB::Red;
        palette1_[i+1] = CRGB::Yellow;
        palette1_[i+2] = CRGB::Black;
        palette1_[i+3] = CRGB::Black;
        palette1_[i+4] = CRGB::Black;
        palette1_[i+5] = CRGB::Black;
        palette1_[i+6] = CRGB::Black;
        palette1_[i+7] = CRGB::Black;
        palette2_[i] = CRGB::Cyan;
        palette2_[i+1] = CRGB::Blue;
        palette2_[i+2] = CRGB::Black;
        palette2_[i+3] = CRGB::Black;
        palette2_[i+4] = CRGB::Black;
        palette2_[i+5] = CRGB::Black;
        palette2_[i+6] = CRGB::Black;
        palette2_[i+7] = CRGB::Black;
      } 
    }
    
    CRGB get_pixel(int x, int y) const {
      if (!swap) {
        if (y % (P * 2) < P ) {
          return palette1_[addmod8(x, pos_, 16)];
        } else {
          return palette2_[(uint8_t)(pos_ - x) % 16];
        }
      } else {
        if (y % (P * 2) < P ) {
          return palette2_[(uint8_t)(pos_ - x) % 16];
        } else {
          return palette1_[addmod8(x, pos_, 16)];
        }
      }
    }
   
    static bool show_bg() { return false; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {
      if (spin_timer_) {
        pos_ = addmod8(pos_, 1, 16);
      }
    } 
           
    uint8_t width() const { return W; }
    uint8_t height() const { return H; }
           
  private:
       
    CRGBPalette16 palette1_;
    CRGBPalette16 palette2_;
    CEveryNMillis spin_timer_;
    uint8_t pos_;
};

template <uint8_t W, uint8_t H, uint8_t COOL, uint8_t SPARK>
class Fire
{
  public:
  
    Fire() {
      random16_add_entropy(random());    
    }
  
    CRGB get_pixel(int x, int y) const {   
      return HeatColor(heat_[x][y]);
    }

    static bool show_bg() { return false; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {
      cool(x);
      rise(x);
      spark(x);    
    } 

    uint8_t width() const { return W; }
    uint8_t height() const { return H; }    
  
  private:
  
    inline void cool(uint8_t x) {
       for (uint8_t y = 0; y < H; ++y) {
         heat_[x][y] = qsub8(heat_[x][y], random(0, ((COOL * 10) / H) + 2));          
      }
    }

    inline void rise(uint8_t x) {
      for (uint8_t y = H - 1; y >= 2; --y) {
        heat_[x][H - 1 - y] = (heat_[x][H - 1 - y + 1] + heat_[x][H - 1 - y + 2] + heat_[x][H - 1 - y + 2]) / 3;
      }       
    }
 
    inline void spark(uint8_t x) {
      if (random8() < SPARK) {
        uint8_t y = random8(3);
        heat_[x][H - 1 - y] = qadd8(heat_[x][H - 1 - y], random8(160, 255));
      }
    }
    
    uint8_t heat_[W][H];
};

template <uint8_t W, uint8_t H, uint8_t BPM, uint8_t PERIOD>
class WaveBall
{
  public:
  
    WaveBall() : palette_(CRGB::Blue, CRGB::Yellow)
    {}
  
    CRGB get_pixel(int x, int y) const {   
      uint8_t brightness = scale16by8(
        beatsin8(BPM, 0, 255, 0, (x % PERIOD) * 256 / PERIOD) 
       +  beatsin8(BPM, 0, 255, 0, (y % PERIOD) * 256 / PERIOD), 128);
      return ColorFromPalette(palette_, brightness); 
    }

    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {} 

    uint8_t width() const { return W; }
    uint8_t height() const { return H; }
    
  private:
    
      CRGBPalette16 palette_;
};

//template <uint8_t W, uint8_t H>
//class Tadpoles
//{
//  public:
//  
//    Tadpoles()
//    {
//      fill_rainbow(palette_.entries, 16, 0, 256/16);
//    }
//  
//    CRGB get_pixel(int x, int y) const {   
//    }
//
//    static bool show_bg() { return true; }    
//    static CRGB bg_color() { return CRGB::Black; }    
//
//    void advance_col(uint8_t x) {} 
//
//    uint8_t width() const { return W; }
//    uint8_t height() const { return H; }
//    
//  private:
//    
//      CRGBPalette16 palette_;
//      uint8_t pixels_[W][H];
//};
