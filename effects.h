#include <FastLED.h>
#include "rotate.h"

class Effect {
public:
	virtual ~Effect() {};
	virtual CRGB get_pixel(int x, int y) const = 0;
	virtual bool show_bg() { return true; };
	virtual CRGB bg_color() { return CRGB::Black; }
	virtual void advance_col(int x) = 0;
	virtual void advance_frame() = 0;
	virtual int width() const = 0;
	virtual int height() const = 0;
};

void trail_rainbow(CHSV& p, uint8_t hue_falloff = 32, uint8_t dim_falloff = 32) {
  CHSV& c = p;
  c.s = 255;
  c.h += hue_falloff;
  c.v = qsub8(c.v, dim_falloff);
}

void trail_rainbow_lin(CHSV& p, uint8_t hue_falloff = 32, uint8_t dim_falloff = 32) {
	CHSV& c = p;
	c.s = 255;
	if (c.v < 128) {
		c.h += hue_falloff / 2;
		c.v = qsub8(c.v, dim_falloff / 2);
	}
	else {
		c.h += hue_falloff;
		c.v = qsub8(c.v, dim_falloff);
	}
}

int wrap(int x, int m) {
	return (x >= 0 ? 
		x % m : 
		((x % m) + m) % m);
}

Fix16 wrap(const Fix16& x, int m) {
	Fix16 m_fix = fix16_from_int(m);
	return (x > int16_t(0) ? 
		fix16_mod(x, m_fix) : 
		fix16_mod(Fix16(fix16_mod(x, m_fix)) + m_fix, m_fix));
}

template <int W, int H>
void plot_aa(CHSV (&leds)[W][H], const Fix16& x, const Fix16& y, const CHSV& c) {
	int16_t x_i = fix16_to_int(fix16_floor(x));
	int16_t y_i = fix16_to_int(fix16_floor(y));
	Fix16 x_m = x - x_i;
	Fix16 y_m = y - y_i;
	const Fix16 one(fix16_one);
	const int16_t FULL = 255;
	int v = fix16_to_int((one - x_m) * (one - y_m) * FULL);
	leds[x_i][y_i] = CHSV(c.h, c.s, dim8_lin(qadd8(leds[x_i][y_i].v, v)));
	v = fix16_to_int(x_m * (one - y_m) * FULL);
	leds[(x_i + 1) % W][y_i] = CHSV(c.h, c.s, dim8_lin(qadd8(leds[(x_i + 1) % W][y_i].v, v)));
	if (y_i < H - 1) {
		v = fix16_to_int((one - x_m) * y_m * FULL);
		leds[x_i][y_i + 1] = CHSV(c.h, c.s, dim8_lin(qadd8(leds[x_i][y_i + 1].v, v)));
		v = fix16_to_int(x_m * y_m * FULL);
		leds[(x_i + 1) % W][y_i + 1] = CHSV(c.h, c.s, dim8_lin(qadd8(leds[(x_i + 1) % W][y_i + 1].v, v)));
	}
}

class Oscillator
{
public:
	
	Oscillator(int min, int max, int end_delay = 0) :
		t_(min),
		min_(min),
		max_(max),
		end_delay_(end_delay),
		end_count_(0),
		dir_(1)
	{}

	int get() {
		int r = t_;
		if ((t_ == max_ + 1  && dir_ == 1) || (t_ == 0 && dir_ == -1)) {
			if (end_count_++ >= end_delay_) {
				dir_ = dir_ * -1;
				end_count_ = 0;
			}
			else {
				return r;
			}
		}
		t_ += dir_;
		return r;
	}

	int dir() { return dir_; }
	void set_delay(int d) { end_delay_ = d; }

private:

	int t_;
	int min_;
	int max_;
	int end_delay_;
	int end_count_;
	int dir_;
};


template <int W, int H>
class RingTwist : public Effect
{
public:
	RingTwist() :
		seed_(104, 0, 255)
	{
		randomSeed(analogRead(PIN_RANDOM_));
		memset(leds_, 0, sizeof(leds_));
		for (int x = 0; x < W; x += W / COUNT_) {
			for (int y = 0; y < H; ++y) {
				leds_[x][y] = seed_;
			}
		}
	}

	CRGB get_pixel(int x, int y) const {
		return leds_[wrap(x - pos_[y], W)][y];
	}

	bool show_bg() { return false; }
	CRGB bg_color() { return CRGB::Black; }

	void advance_col(int x) {
		for (int y = 0; y < H; ++y) {
			if (leds_[x][y] != seed_) {
				decay(x, y);
			}
		}
	}

	void advance_frame() {
		pull(leader_);
		EVERY_N_MILLISECONDS(STOP_TIMER_) {
			if (stop_ < 0) {
				stop_ = pos_[leader_];
			}
		}
	}

	int width() const { return W; }
	int height() const { return H; }

private:

	void decay(int x, int y) {
		trail_rainbow(leds_[x][y], 8, 24);
	}

	void pull(int y) {
		if (stop_ != -1 && pos_[y] == stop_) {
			block(y);
			return;
		}
		move(y, pos_[y], wrap(pos_[y] + dir_, W), dir_);
		for (int i = y - 1; i >= 0; --i) {
			if (distance(pos_[i], pos_[i + 1], dir_) > lead_length_) {
				move(i, pos_[i], wrap(pos_[i + 1] - lead_length_ * dir_, W), dir_);
			}
		}
		for (int i = y + 1; i < H; ++i) {
			if (distance(pos_[i], pos_[i - 1], dir_) > lead_length_) {
				move(i, pos_[i], wrap(pos_[i - 1] - lead_length_ * dir_, W), dir_);
			}
		}
	}

	void block(int y) {
		bool all_stop = true;
		for (int i = 0; i < H; ++i) {
			if (distance(pos_[i], pos_[y], dir_) != 0) {
				move(i, pos_[i], wrap(pos_[i] + dir_, W), dir_);
				all_stop = false;
			}
		}
		if (all_stop) {
			stop_ = -1;
			dir_ = random(2) < 1 ? 1 : -1;
			leader_ = random(0, H);
			lead_length_ = lead_lengths_[random(1, sizeof(lead_lengths_) / sizeof(int))];
		}
	}

	int distance(int a, int b, int dir) {
		return dir > 0 ?
			wrap(b - a, W) :
			wrap(a - b, W);
	}

	void move (int y, int x0, int x1, int dir) {
		pos_[y] = x1;
		for (int x = 0; x < W; x += W / COUNT_) {
			for (int d = 0; d < distance(x0, x1, dir); ++d) {
				add_trail(y, wrap(x - dir, W), dir);
			}
		}
	}

	void add_trail(int y, int x, int dir) {
		if (leds_[wrap(x + dir, W)][y].v == 0 || leds_[x][y] == seed_) {
			return;
		}
		add_trail(y, wrap(x - dir, W), dir);
		leds_[x][y] = leds_[wrap(x + dir, W)][y];
		decay(x, y);
	}

	const CHSV seed_;
	const int COUNT_ = 2;
	const int STOP_TIMER_ = 20000;
	const int PIN_RANDOM_ = 15;
	CHSV leds_[W][H];
	int pos_[H] = { 0 };
	int leader_ = 0;
	int dir_ = 1;
	int stop_ = -1;
	int lead_length_ = 1;
	int lead_lengths_[7] = { 1, 2, 3, 4, 6, 12, 16 };
};

template <int W, int H>
class Curves : public Effect {
public:
	Curves() {
		memset(leds_, 0, sizeof(leds_));
		fill_gradient<CHSV>(palette_, sizeof(palette_) / sizeof(CHSV), 
			CHSV(0, 255, 255), 
			CHSV(48, 255, 255));
	}

	CRGB get_pixel(int x, int y) const {
		return leds_[x][y];
	}

	bool show_bg() { return false; }
	CRGB bg_color() { return CHSV(0, 0, 24); }

	void advance_col(int x) {
		for (int y = 0; y < H; ++y) {
			trail_rainbow_lin(leds_[x][y], 32, 48);			
		}
	}

	void advance_frame() {
		for (int i = 0; i < COUNT; ++i) {
			plot(num_, 255, (i * W / COUNT + offset_) % W, palette_[triwave8(color_offset_)]);
			plot(-num_, 255, (i * W / COUNT + W - 1 - offset_) % W, palette_[triwave8(color_offset_)]);
		}
		EVERY_N_MILLIS(150) {
			offset_ = (offset_ + 1) % W;
		}
		num_ = wrap(num_ - 1, 1024);
		color_offset_++;
	}

	int width() const { return W; }
	int height() const { return H; }

private:

	void plot(int n, int d, int c, const CHSV& color) {
		for (int y = 0; y < H; ++y) {
			leds_[wrap(y * n / d + c, W)][y] = color;
		}
	}

	int COUNT = 4;
	int offset_ = 0;
	uint8_t color_offset_ = 0;
	uint16_t num_ = 1024;
	uint8_t falloff_ = 48;
	CHSV leds_[W][H];
	CHSVPalette256 palette_;
};


template <int W, int H>
class WaveTrails : public Effect {
  public:
    WaveTrails() {
      random16_add_entropy(random());    
      memset(leds, 0, sizeof(leds));
    }
    
    CRGB get_pixel(int x, int y) const {
      return leds[x][y];
    }
        
    bool show_bg() { return false; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow(leds[x][y], 5, 24);
      }
      leds[x][beatsin16(5, H/15, H*14/15, 0, 
        map(x, 0, W-1, 0, 65535))] = CHSV(HUE_BLUE, 100, 255); 
      leds[x][beatsin16(5, 0, H-1, 0, 
        map(W - 1 - x, 0, W-1, 0, 65535))] = CHSV(HUE_GREEN, 100, 255); 
      leds[x][beatsin16(5, 0, H-1, 32767, 
        map(W - 1 - x, 0, W-1, 0, 65535))] = CHSV(HUE_AQUA, 100, 255); 
    } 

    void advance_frame() {
      
    }
    
    int width() const { return W; }
    int height() const { return H; }

private:

  CHSV leds[W][H];
};

template <uint8_t W, uint8_t H>
class DotTrails : public Effect {
  public:
    DotTrails() {
      random16_add_entropy(random());    
      memset(leds, 0, sizeof(leds));
      for (int i = 0; i < H - 1; ++i) {
        bool rev = i % 2 == 0;
        uint8_t x = random8() % W;
        dots[i] = Dot(x, rev, (random8() % 30) + 20);
      }
    }
    
    CRGB get_pixel(int x, int y) const {
      return leds[x][y];
    }
        
    bool show_bg() { return false; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      for (int y = 0; y < H; ++y) {
        if (dots[y].x == x) {
          if (dots[y].drawn) {
            dots[y].drawn = false;
          } else {
            leds[x][y] = CHSV(hue - 12, 0, 255);
            move(dots[y]);
          }
        } else {
          trail_rainbow(leds[x][y], 12, 12);
        }
      }
    } 

    void advance_frame() {
    }
    
    int width() const { return W; }
    int height() const { return H; }

private:

  struct Dot {
    Dot() : 
      x(0),
      rev(0),
      ttl(0),
      drawn(0){}

    Dot(uint8_t x, bool rev, int ttl) : 
      x(x),
      rev(rev),
      ttl(ttl),
      drawn(0){}
 
    uint8_t x;
    bool rev;
    int ttl;
    bool drawn;
  };

  void move(Dot& dot) {
    if (dot.ttl == 0) {
      dot.x = random8() % W;
      dot.ttl = (random8() % 30) + 20;
    } else if (dot.rev) {
      dot.x = (dot.x - 1 + W) % W;
      dot.ttl--;
    } else {
      dot.x = addmod8(dot.x, 1, W);
      dot.drawn = true;
      dot.ttl--;
    }
  }

  CHSV leds[W][H];
  Dot dots[H];
  uint8_t hue = HUE_RED;
};

template <uint8_t W, uint8_t H>
class RingTrails : public Effect {
public:
	RingTrails() {
		memset(leds, 0, sizeof(leds));
		memset(led_shown, 0, sizeof(led_shown));
	}

	CRGB get_pixel(int x, int y) const {
		led_shown[x][y]++;
		return leds[x][y];
	}

	bool show_bg() { return false; }
	CRGB bg_color() { return CRGB::Black; }

	void advance_col(int x) {
		Point p = projection.project(x, (H - 1) / 2);
		leds[p.x][p.y] = CHSV(hue - 32, 0, 255);
		led_shown[p.x][p.y] = 0;

		for (int y = 0; y < H; ++y) {
			if (led_shown[x][y] > 1) {
				trail_rainbow(leds[x][y], 32, 32);
			}
		}
	}

	void advance_frame() {
		uint8_t dl = beatsin8(5, 1, 2);
		uint8_t dg = beatsin8(7, 1, 2, 16384);
		uint8_t dp = beatsin8(11, 1, 2, 32768);
		projection.rotate(dl, dg, dp);
	}

	int width() const { return W; }
	int height() const { return H; }

private:

	typedef typename Projection<W, H>::Point Point;
	Projection<W, H> projection;

	CHSV leds[W][H];
	mutable unsigned int led_shown[W][H];
	uint8_t hue = HUE_RED;
};

template <uint8_t W, uint8_t H, const unsigned char (*DATA)[20][3]>
class Image : public Effect
{
  public:
    
    Image()
    {}
    
    CRGB get_pixel(int x, int y) const {
      uint32_t c = pgm_read_dword(&DATA[x][y][0]); 
      return CRGB(c & 0xff, (c & 0xff00) >> 8, (c & 0xff0000) >> 16);
    }

    bool show_bg() { return true; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {} 
    void advance_frame() {} 

    int width() const { return W; }
    int height() const { return H; }

  private:

    const unsigned char (*data_)[H][3];  
};

template <uint8_t W, uint8_t H, uint8_t S>
class Snake : public Effect
{
  public:
    Snake() : 
    head_x_(0),
    head_y_(0),
    head_x2_(W / 2),
    head_y2_(H / 2),
    x_timer_(125),
    y_timer_(125)
    {}
    
    CRGB get_pixel(int x, int y) const {
       if (is_snake(x, y, 0, 0)) {
          return CRGB::Green;
       }
       else if (is_snake(x, y, 5, 10)) {
          return CRGB::Lime;
       }
       else if (is_snake(x, y, 10, 5)) {
          return CRGB::DarkGreen;
       }
       else if (is_snake(x, y, 5, 15)) {
          return CRGB::Lime;
       }
       else if (is_snake(x, y, 15, 5)) {
          return CRGB::Green;
       }
       else if (is_snake(x, y, 10, 15)) {
          return CRGB::DarkGreen;
       }
       else if (is_snake(x, y, 15, 10)) {
          return CRGB::Violet;
       }
       return CRGB::Black;
    }
        
    bool show_bg() { return true; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      if (x_timer_) {
        head_x_ = addmod8(head_x_, 1, W);
        head_x2_ = addmod8(head_x2_, 1, W);
      }
      if (y_timer_) {
        head_y_ = addmod8(head_y_, 1, H);
        head_y2_ = addmod8(head_y2_, 1, H);
      }
    } 
    void advance_frame() {} 
    
    int width() const { return W; }
    int height() const { return H; }

private:

    bool is_snake(int x, int y, int x_offset, int y_offset) const {
       int head_x = (head_x_ + x_offset) % W;
       int head_y = (head_y_ + y_offset) % H;
       int head_x2 = (head_x2_ + x_offset) % W;
       int head_y2 = (head_y2_ + y_offset) % H;
       uint8_t x2 = W + 1 - x;
       uint8_t y2 = H + 1 - y;
       return (
              (x == head_x && y == head_y)
           || (x2 == head_x && y2 == head_y)
           || (x == head_x && y2 == head_y)
           || (x2 == head_x && y == head_y)
           || (x == head_x2 && y == head_y2)
           || (x2 == head_x2 && y2 == head_y2)
           || (x == head_x2 && y2 == head_y2)
           || (x2 == head_x2 && y == head_y2)
           
          );  
    }
  
    uint8_t head_x_;
    uint8_t head_y_;
    uint8_t head_x2_;
    uint8_t head_y2_;
    CEveryNMillis x_timer_;
    CEveryNMillis y_timer_;
};

template <uint8_t W, uint8_t H>
class Kaleidoscope : public Effect
{
  public:
    Kaleidoscope() {
      memset(leds_, 0, sizeof(leds_));
      fill_rainbow(palette_.entries, 16, 0, 256 / 16);
    }
    
    CRGB get_pixel(int x, int y) const {
      return leds_[x][y];
    }
        
    bool show_bg() { return false; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {  
      nscale8(&leds_[x][0], H, 60);      
    }
    
    void advance_frame() {
      for (uint8_t i = 0; i < counts_[count_]; ++i) {
        uint8_t x = addmod8(i * offset_, tx_, W); 
        uint8_t y = ty_;
        uint8_t x2 = W - 1 - x;
        uint8_t y2 =  H - 1 - y;
        leds_[x][y] = palette_[addmod8(palette_offset_, i * (16 / counts_[count_]), 16)];
        leds_[x2][y2] = leds_[x][y], palette_[addmod8(palette_offset_ , 16 - (i * (16 / counts_[count_])), 16)];
      }
      tx_ = addmod8(tx_, 2 , W); 
      ty_ = addmod8(ty_, inc_y_, H);
      if (ty_ == H - 1 || ty_ == 0) {
        inc_y_ = 0 - inc_y_;
        tx_ = addmod8(tx_,2, W); 
      }

      EVERY_N_MILLISECONDS(100) {
         palette_offset_ = addmod8(palette_offset_, 1, 16);
      }
      EVERY_N_MILLISECONDS(3000) {
        count_ = addmod8(count_, 1, sizeof(counts_));
        offset_ = W / counts_[count_];
      }
    }
    
    int width() const { return W; }
    int height() const { return H; }

private:
  
  CRGB leds_[W][H];
  uint8_t counts_[9] = {1, 2, 4, 5, 6, 8, 10, 12, 20};
  uint8_t count_ = 0;
  uint8_t offset_ = W / counts_[count_];
  uint8_t tx_ = 0;
  uint8_t ty_ = 0;
  int inc_y_ = 1;
  CRGBPalette16 palette_;
  uint8_t palette_offset_ = 0;
};

template <uint8_t W, uint8_t H, uint8_t SPREAD>
class Spiral : public Effect
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
        
    bool show_bg() { return true; }    
    CRGB bg_color() { return CRGB::Black; }    
   
    void advance_col(int x) {} 
    void advance_frame() {} 
  
    int width() const { return W; }
    int height() const { return H; }

  private:
        
    CRGBPalette16 palette_;
};

template <uint8_t W, uint8_t H>
class StarsFade : public Effect
{
  public:
    StarsFade() :
    hue_(0)
    {
      memset(leds_, 0, sizeof(leds_));
    }
    
    CRGB get_pixel(int x, int y) const {
        return leds_[x][y];
    }
        
    bool show_bg() { return true; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
        for (int y = 0; y < H; ++y) {
			if (random8() < 3) {
				leds_[x][y] = CHSV(hue_, 255, 255);
			}
			trail_rainbow(leds_[x][y], 16, 32);
		}
    }
    
    void advance_frame() {
      hue_++;
    } 
    
    int width() const { return W; }
    int height() const { return H; }
    
    private:

      CHSV leds_[W][H];
      uint8_t hue_;
};

template <uint8_t W, uint8_t H>
class Spinner : public Effect
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
        if (y % (p[i] * 2) < p[i] ) {
          return palette1_[addmod8(x, pos_, 16)];
        } else {
          return palette2_[(uint8_t)(pos_ - x) % 16];
        }
      } else {
        if (y % (p[i] * 2) < p[i] ) {
          return palette2_[(uint8_t)(pos_ - x) % 16];
        } else {
          return palette1_[addmod8(x, pos_, 16)];
        }
      }
    }
   
    bool show_bg() { return false; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      if (spin_timer_) {
        pos_ = addmod8(pos_, 1, W);
      }
    } 

    void advance_frame() {
      EVERY_N_SECONDS(5) {
        if ((!swap && (i == sizeof(p) - 1)) 
        || (swap && (i == 0)) ) {
          swap = !swap;
        } else {
          i += swap ? -1 : 1;
        }
      }
    } 
           
    int width() const { return W; }
    int height() const { return H; }
           
  private:
       
    CRGBPalette16 palette1_;
    CRGBPalette16 palette2_;
    CEveryNMillis spin_timer_;
    uint8_t pos_;
    bool swap = false;
    uint8_t p[5] = { 20, 10, 4, 2, 1 };
    uint8_t i = 0;
};

template <uint8_t W, uint8_t H, uint8_t COOL, uint8_t SPARK>
class Fire : public Effect
{
  public:
  
    Fire() {
      random16_add_entropy(random());    
    }
  
    CRGB get_pixel(int x, int y) const {   
      return HeatColor(heat_[x][y]);
    }

    bool show_bg() { return false; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      cool(x);
      rise(x);
      spark(x);    
    } 

    void advance_frame() {} 

    int width() const { return W; }
    int height() const { return H; }    
  
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

template <uint8_t W, uint8_t H, uint8_t HUE>
class TheMatrix : public Effect
{
  public:

    TheMatrix()
    {
      memset8(pixels_, 0, sizeof(pixels_));
      random16_add_entropy(random());
    }

    CRGB get_pixel(int x, int y) const {
      return blend(CHSV(HUE + (3 * y), 255, 40), 
                        CHSV(HUE + (3 * (H - 1 - y)), 255, 255),
                        pixels_[x][y]);
    }
    bool show_bg() { return true; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      fall(x);
      generate(x);
    }
    
    void advance_frame() {
    } 

    int width() const { return W; }
    int height() const { return H; }

  private:

    void fall(int x) {
      memmove8(&pixels_[x][1], &pixels_[x][0], H - 1);
      pixels_[x][0] = 0;
    }
    
    void generate(int x) {
      if (random8() < 15) {
        pixels_[x][2] = 255;
        pixels_[x][1] = 100;
        pixels_[x][0] = 30;
      }
    }
    
    uint8_t pixels_[W][H];
};

template <uint8_t W, uint8_t H, uint8_t HUE, uint8_t BURNRATE>
class Burnout : public Effect
{
  public:

    Burnout() :
    timer_(30),
    burn_idx_(0)
    {
      random16_add_entropy(random());
      memset8(pixels_, INIT, sizeof(pixels_));
      fill_seq(burn_, W * H);
      shuffle(burn_, W * H);
    }

    CRGB get_pixel(int x, int y) const {
      if (pixels_[x][y] & BURN) {
          return CHSV(HUE + (3 * y), 255, 255);
      }
      if (pixels_[x][y] == INIT) {
          return CHSV(HUE + (3 * y), 255, 40);
      }
      return CRGB::Black;
    }
    
    bool show_bg() { return true; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      if (timer_) {
        burnout();
      }
      fall(x);
    }
    
    void advance_frame() {
    } 

    int width() const { return W; }
    int height() const { return H; }

  private:

    enum State {
      NONE,
      INIT,
      BURN = 0x80,
    };
    
    void burnout() {
      if (burn_idx_ < W * H) {
        *((uint8_t *)pixels_ + burn_[burn_idx_++]) = BURN;
      }
    }

    void fall(int x) {
      for (int y = 0; y < H; ++y) {
        if (pixels_[x][y] & BURN) {
          pixels_[x][y] ^= BURN;
          if (y < (H - 1)) {
            pixels_[x][y + 1] |= BURN;
            ++y;          
          }
        }
      }
    }

    void fill_seq(short *arr, int len) {
      for (int i = 0; i < len; ++i) {
        arr[i] = i;
      } 
    }

    void shuffle(short *arr, int len) {
      for (int i = 0; i < len - 1; ++i) {
        short j = random(i, len);
        short temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
      }
    }
    
    CEveryNMillis timer_;
    
    uint8_t pixels_[W][H];
    short burn_[W * H];
    int burn_idx_;
};

template <uint8_t W, uint8_t H>
class Rotate : public Effect
{
  public:

    Rotate() {
      memset(buf, 0, sizeof(buf));
      fill_rainbow(pal.entries, 16, 0, 256 / 16);
    }

    CRGB get_pixel(int x, int y) const {
      return buf[x][y];
    }

    bool show_bg() { return bg; }    
    CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      if (x == 0 || x == W / 2) {
        projection.rotate(0, 5, 0);
      }
      for (int y = 0; y < H; ++y) {
        Point p = projection.project(x, y);
        if (p.y == 3 || p.y == 10 || p.y == 17) {
          buf[x][y] =  ColorFromPalette(pal, (p.x + c_off) * 255 / W);
        } else {
          buf[x][y] = CRGB::Black;
        }
      }      
    }
    
    void advance_frame() {
        ++c_off;
        EVERY_N_SECONDS(10) {
          bg = !bg;
        }
    }

    int width() const { return W; }
    int height() const { return H; }

  private:
  
    typedef typename Projection<W, H>::Point Point;
    Projection<W, H> projection;

    CRGB buf[W][H];
    CRGBPalette16 pal;
    uint8_t c_off = 0;
    bool bg = 0;
};

template <uint8_t W, uint8_t H>
class RotateWave : public Effect
{
  public:

    RotateWave() {
      memset(buf, 0, sizeof(buf));
      fill_rainbow(pal.entries, 16, 0, 256 / 16);
    }

    CRGB get_pixel(int x, int y) const {
      return buf[x][y];
    }

    bool show_bg() { return false; }    
    static CRGB bg_color() { return CRGB::Black; }    

    void advance_col(int x) {
      if (x == 0 || x == W / 2) {
        projection.rotate(0, 5, 0);
      }
      for (int y = 0; y < H; ++y) {
        Point p = projection.project(x, y);
        if (p.y == 2 || p.y == 18 || equals_wave(p)) {
          buf[x][y] =  ColorFromPalette(pal, (p.x + c_off) * 255 / W);
        } else {
          buf[x][y] = CRGB::Black;
        }
      }      
    }
    
    void advance_frame() {
        ++c_off;
    }

    int width() const { return W; }
    int height() const { return H; }

  private:
  
    typedef typename Projection<W, H>::Point Point;

    bool equals_wave(const Point& p) const {
      uint8_t freq = W / 3;
      int amp = beatsin8(5, 0, 3);
      uint8_t origin = H / 2;
      return p.y == map8(triwave8(map(mod8(p.x, freq), 0, freq - 1 , 0, 255)), origin - amp, origin + amp);
    }

    CRGB buf[W][H];
    CRGBPalette16 pal;
    uint8_t c_off = 0;
    Projection<W, H> projection;
};



