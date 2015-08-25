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

    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {} 
    void advance_frame() {} 

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
    void advance_frame() {} 
    
    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

};

template <uint8_t W, uint8_t H, uint8_t SPREAD>
class Spiral
{
  public:
    Spiral()
    {
      fill_rainbow(palette_.entries, 16, 0, 256 / 16);
    }
    
    CRGB get_pixel(int x, int y) const {
      if ((x + y) % (SPREAD + 1) == 0) {
        return palette_[(x + y) % 16];
      }
      return CRGB::Black;
    }
        
    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    
   
    void advance_col(uint8_t x) {} 
    void advance_frame() {} 
  
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
    void advance_frame() {} 
    
    uint8_t width() const { return W; }
    uint8_t height() const { return H; }
};

template <uint8_t W, uint8_t H, uint8_t P, bool swap = false>
class Spinner
{
  public:
    Spinner() :
      spin_timer_(75),
      pos_(0)
    {
      memset(palette1_, 0, sizeof(palette1_));
      memset(palette2_, 0, sizeof(palette1_));
      palette1_[0] = palette1_[8] = CRGB::Red;
      palette1_[1] = palette1_[9] = CRGB::Yellow;
      palette2_[0] = palette2_[8] = CRGB::Cyan;
      palette2_[1] = palette2_[9] = CRGB::Blue; 
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
        pos_ = addmod8(pos_, 1, W);
      }
    } 

    void advance_frame() {
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

    void advance_frame() {} 

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

template <uint8_t W, uint8_t H, uint8_t SPEED>
class Water
{
  public:
  
    Water() : 
      palette_(CRGB::Blue, CRGB::Black),
      frame_(0)
    {
      for (uint8_t i = 0; i < W; ++i) {
        if (i % 3 == 0) {
          water_[i] = 0;
          continue;
        } 
        uint8_t j = (i + 5) % 6;
        if (j == 0 || j == 1) {
           water_[i] = 1;
           continue;
        } 
        uint8_t k = (i + 2) % 6;
        if (k == 0 || k == 1) {
           water_[i] = -1;
           continue;
        } 
      }
    }
  
    CRGB get_pixel(int x, int y) const {
      uint8_t water_level = get_water_level(x);
      if (y >= water_level) {
        return palette_[max(0, y - (H - 16))];
      }
      return CRGB::Black;
    }

    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {
      frame_ = (millis() / SPEED) % 4;
      if (frame_ == 3) {
        frame_ = 1;
      } 
    } 

    void advance_frame() {} 

    uint8_t width() const { return W; }
    uint8_t height() const { return H; }
    
  private:
  
    uint8_t get_water_level(int x) const {
      switch (frame_) {
        case 0:
         return H/2 + water_[x]; 
        case 1:
         return H/2;
        case 2:
         return H/2 - water_[x]; 
      }
    }

    int8_t water_[W];
    uint8_t frame_;
    CRGBPalette16 palette_;
};

template <uint8_t W, uint8_t H, uint16_t DURATION, bool BG>
class PaletteFall
{
  public:

    PaletteFall() :
      timer_(DURATION),
      palette_idx_(0),
      palette_(RainbowColors_p),
      palette_offset_(0)
    {
    }

    CRGB get_pixel(int x, int y) const {
      return palette_[((H - 1 - y) + palette_offset_) % 16];
    }

    static bool show_bg() { return BG; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {
    }
    
    void advance_frame() {
      palette_offset_ = addmod8(palette_offset_, 1, 16);
      if (timer_) {
        switch_palette();
      }
    } 

    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

  private:

    void switch_palette() {
      switch (palette_idx_) {
        case 0:
          palette_ = RainbowColors_p;
          break;
        case 2:
          palette_ = PartyColors_p;
          break;
        case 3:
          palette_ = HeatColors_p;
          break;
        case 4:
          palette_ = RainbowStripeColors_p;
          break;
      }
      palette_idx_ = addmod8(palette_idx_, 1, 5);
    }

    CEveryNMillis timer_;
    uint8_t palette_idx_;
    CRGBPalette16 palette_;
    uint8_t palette_offset_;
};

template <uint8_t W, uint8_t H, uint16_t DURATION>
class PaletteGrid
{
  public:

    PaletteGrid() :
      timer_(DURATION),
      palette_idx_(0),
      palette_(RainbowColors_p),
      palette_offset_(0)
    {
    }

    CRGB get_pixel(int x, int y) const {
      x = x % 4;
      y = y % 4;
      if (x == 2 && y == 2) {
        return CRGB::Red;
//        return palette_[(palette_offset_ * 4 + 2) % 16];        
      };
      if (x == 0 || y == 0) {
        return CRGB::Blue;
//        return palette_[(palette_offset_ * 4 + 1) % 16];                
      }
      return CRGB::Red;
//      return palette_[palette_offset_ * 4];
    }

    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {
    }
    
    void advance_frame() {
//      palette_offset_ = addmod8(palette_offset_, 1, 16);
//      if (timer_) {
//        switch_palette();
//      }
    } 

    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

  private:

    void switch_palette() {
      switch (palette_idx_) {
        case 0:
          palette_ = RainbowColors_p;
          break;
        case 2:
          palette_ = PartyColors_p;
          break;
        case 3:
          palette_ = HeatColors_p;
          break;
        case 4:
          palette_ = RainbowStripeColors_p;
          break;
      }
      palette_idx_ = addmod8(palette_idx_, 1, 5);
    }

    CEveryNMillis timer_;
    uint8_t palette_idx_;
    CRGBPalette16 palette_;
    uint8_t palette_offset_;
};

template <uint8_t W, uint8_t H, uint8_t HUE>
class TheMatrix
{
  public:

    TheMatrix() :
      bgcolor_(CHSV(HUE, 255, 20)),
      fgcolor_(CHSV(HUE, 255, 255))
    {
      memset8(pixels_, 0, sizeof(pixels_));
      random16_add_entropy(random());
    }

    CRGB get_pixel(int x, int y) const {
      return blend(bgcolor_, fgcolor_, pixels_[x][y]);
    }

    static bool show_bg() { return true; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(uint8_t x) {
      fall(x);
      generate(x);
    }
    
    void advance_frame() {
    } 

    uint8_t width() const { return W; }
    uint8_t height() const { return H; }

  private:

    void fall(int x) {
      memmove8(&pixels_[x][1], &pixels_[x][0], H - 1);
      pixels_[x][0] = 0;
    }
    
    void generate(int x) {
      if (random8() < 20) {
        pixels_[x][2] = 255;
        pixels_[x][1] = 128;
        pixels_[x][0] = 50;
      }
    }
    
    uint8_t pixels_[W][H];
    CRGB bgcolor_;
    CRGB fgcolor_;
};


//
//template <uint8_t W, uint8_t H>
//class Plaid
//{
//  public:
//    
//    Plaid()
//    {}
//    
//    CRGB get_pixel(int x, int y) const {
//      
//    }
//
//    static bool show_bg() { return true; }    
//    static CRGB bg_color() { return CRGB::Black; }    
//
//    void advance_col(uint8_t x) {} 
//    void advance_frame() {} 
//
//    uint8_t width() const { return W; }
//    uint8_t height() const { return H; }
//
//  private:
//
//  const CRGB c1_ = CRGB::Red;
//  const CRGB c2_ = CRGB::Green;
//  const CRGB c3_ = CRGB::Blue;
//};


