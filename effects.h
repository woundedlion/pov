#include <FastLED.h>
#include "particle.h"
#include "rotate.h"
#include "util.h"

///////////////////////////////////////////////////////////////////////////////

CRGB dim(const CRGB& c, double s) {
  if (s < 0.5) s = s * s;
  return CRGB(
    static_cast<uint8_t>(c.r * s),
    static_cast<uint8_t>(c.g * s),
    static_cast<uint8_t>(c.b * s));
}

void trail_rainbow(CRGB& c, uint8_t hue_falloff = 32,
                   uint8_t dim_falloff = 32) {
  auto p = rgb2hsv_approximate(c);
  p.s = 255;
  p.h += hue_falloff;
  p.v = qsub8(p.v, dim_falloff);
  c = p;
}

void trail_rainbow_lin(CRGB& c, uint8_t hue_falloff = 32,
                       uint8_t dim_falloff = 32) {
  auto p = rgb2hsv_approximate(c);
  p.s = 255;
  if (p.v < 128) {
    p.h += hue_falloff / 2;
    p.v = qsub8(p.v, dim_falloff / 2);
  } else {
    p.h += hue_falloff;
    p.v = qsub8(p.v, dim_falloff);
  }
  c = p;
}

void disintegrate(CHSV& p, int prob) {
  if (random8() <= prob) {
    p = CHSV(p.h, p.s, 0);
  }
}

/*
template <int W, int H>
void plot_aa(CHSV (&leds)[W][H], const Fix16& x, const Fix16& y,
             const CHSV& c) {
  int16_t x_i = fix16_to_int(fix16_floor(x));
  int16_t y_i = fix16_to_int(fix16_floor(y));
  Fix16 x_m = x - x_i;
  Fix16 y_m = y - y_i;
  const Fix16 one(fix16_one);
  const int16_t FULL = 255;
  int v = fix16_to_int((one - x_m) * (one - y_m) * FULL);
  leds[x_i][y_i] = CHSV(c.h, c.s, dim8_lin(qadd8(leds[x_i][y_i].v, v)));
  v = fix16_to_int(x_m * (one - y_m) * FULL);
  leds[(x_i + 1) % W][y_i] =
      CHSV(c.h, c.s, dim8_lin(qadd8(leds[(x_i + 1) % W][y_i].v, v)));
  if (y_i < H - 1) {
    v = fix16_to_int((one - x_m) * y_m * FULL);
    leds[x_i][y_i + 1] =
        CHSV(c.h, c.s, dim8_lin(qadd8(leds[x_i][y_i + 1].v, v)));
    v = fix16_to_int(x_m * y_m * FULL);
    leds[(x_i + 1) % W][y_i + 1] =
        CHSV(c.h, c.s, dim8_lin(qadd8(leds[(x_i + 1) % W][y_i + 1].v, v)));
  }
}
*/

void plot_aa(Canvas& cv, const float& x, const float& y, const CRGB& c) {
  int x_i = floor(x);
  int y_i = floor(y);
  double x_m = x - x_i;
  double y_m = y - y_i;

  double v = (1 - x_m) * (1 - y_m);
  cv(x_i, y_i) = dim(c, v);

  v = x_m * (1 - y_m);
  cv((x_i + 1) % cv.width(), y_i) = dim(c, v);

  if (y_i < H - 1) {
    v = (1 - x_m) * y_m;
    cv(x_i, y_i + 1) = dim(c, v);

    v = x_m * y_m;
    cv((x_i + 1) % cv.width(), y_i + 1) = dim(c, v);
  }
}

class PollingRandomTimer {
 public:
   PollingRandomTimer(uint32_t min_ms, uint32_t max_ms)
      : min_ms(min_ms), max_ms(max_ms) {
    randomSeed(analogRead(PIN_RANDOM));
    next = millis() + random(min_ms, max_ms);
  }

  bool elapsed() {
    if (millis() >= next) {
      next = millis() + random(min_ms, max_ms);
      return true;
    }
    return false;
  }

 private:
  uint32_t min_ms;
  uint32_t max_ms;
  uint32_t next;
};


/*
template <int W>
class Test : public Effect
{
public:

        Test() :
                Effect(W)
        {
                randomSeed(analogRead(PIN_RANDOM));
                int i = 0;
                for (int y = 0; y < H && i < N; ++y) {
                        for (int x = 0; x < W && i < N; ++x) {
                                dots.add_particle(Particle(x, y, Vector(0, 0)));
                                ++i;
                        }
                }
        }

        bool show_bg() const { return true; }

        void draw_frame() {
                Canvas c(*this);
                for (int x = 0; x < W; ++x) {
                        for (int y = 0; y < H; ++y) {
                                c(x, y) = CHSV(0, 0, 0);
                        }
                }
                for (int i = 0; i < dots.size(); ++i) {
                        Particle& p = dots[i];
                        plot_aa(c, p.x, p.y, CHSV(HUE_BLUE, 255, 255));
                }

                Serial.printf("frame %d\n", millis());
                dots.step(1);
        }

private:

        static const int N = 160;

        typedef typename ParticleSystem<W, H, N>::Particle Particle;
        typedef typename ParticleSystem<W, H, N>::Vector Vector;

        ParticleSystem<W, H, N> dots;
};
*/

class ComplementaryColorSequence {
public:
  ComplementaryColorSequence() : last(0, random(128, 255), 255) {}

  CHSV get() {
    CHSV r;
    if (random8() < 32) {
      last = CHSV(last.h + random(16, 32), random(128, 255), 255);
      return last;
    }
    else if (random8() < 64) {
      last = CHSV(last.h - random(16, 32), random(128, 255), 255);
      return last;
    }
    r = CHSV(last.h + 64, last.s, last.v);
    last = r;
    return r;
  }

private: 
  CHSV last;
};

///////////////////////////////////////////////////////////////////////////////////////////////////

/*
template <int W>
class Hoops : public Effect {
public:
  Hoops() : Effect(W), osc(1, 4) {
    for (size_t i = 0; i < (sizeof(hoops) / sizeof(Hoop)); ++i) {
      hoops[i] = Hoop(0, i, 0);
    }
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);

    EVERY_N_MILLIS(1000) { speed *= -1; }
    EVERY_N_MILLIS(4000) {
      gap = osc.get();
      speed = dir(speed) * random(1, 2);
    }

    pull(c, dots, 0, speed);
  };

private:
  struct Hoop {
    Hoop() : x(0), y(0), v(0), hue(HUE_RED) {}

    Hoop(int x, int y, int v) : x(x), y(y), v(v), hue(HUE_RED) {}

    Hoop(const Hoop& d) : x(d.x), y(d.y), v(d.v), hue(d.hue) {}

    int x;
    int y;
    int v;
    uint8_t hue;
  };

  int distance(int a, int b, int m) {
    int fwd = wrap(b - a, m);
    int rev = wrap(a - b, m);
    if (fwd <= rev) {
      return fwd;
    }
    return -rev;
  }

  int dir(int v) { return v < 0 ? -1 : 1; }

  void pull(Canvas& c, Hoop* hoops, int y, int v) {
    hoops[y].v = v;
    move(c, hoops[y]);
    for (int i = y - 1; i >= 0; --i) {
      drag(c, hoops[i + 1], hoops[i]);
    }
    for (int i = y + 1; i < H; ++i) {
      drag(c, hoops[i - 1], hoops[i]);
    }
  }

  void drag(Canvas& c, Hoop& leader, Hoop& follower) {
    int dest = wrap(follower.x + follower.v, W);
    if (abs(distance(dest, leader.x, W)) > gap) {
      // Move to furthest extent based on leader's prior path
      dest = wrap(leader.x - dir(leader.v) * (abs(leader.v) - 1 + gap), W);
      follower.v = distance(follower.x, dest, W);
      move(c, follower);
      // Move to gap's length behind leader
      dest = wrap(leader.x - dir(leader.v) * gap, W);
      follower.v = distance(follower.x, dest, W);
      move(c, follower);
      // Adjust speed to match leader
      follower.v = leader.v;
    }
    else {
      move(c, follower);
    }
  }

  void move(Canvas& c, Hoop& hoop) {
    hoop.x = wrap(dot.x + dot.v, W);
  }

  Hoop hoops[H];
  Oscillator osc;
  int gap = 1;
  int speed = 1;
};

*/

template <int W>
class RingShower : public Effect {
 public:
  RingShower() : Effect(W), new_ring_timer(62, 2000) {
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        c(x, y) = CHSV(0, 0, 0);
      }
    }

    if (new_ring_timer.elapsed()) {
      new_ring();
    }
    draw_rings(c);
    move_rings(c);
  }

 private:
  struct Ring {
    Ring() : position(H), speed(0), color(CHSV(0, 0, 0)) {}

    float position;
    float speed;
    Projection<W, H> rotation;
    CHSV color;
  };

  void new_ring() {
    for (size_t i = 0; i < sizeof(rings) / sizeof(Ring); ++i) {
      if (rings[i].position >= H) {
        rings[i].position = 0;
        rings[i].speed = random(30, 100) / 100.0;
        rings[i].color = colors.get();
        rings[i].rotation.rotate(random(360), random(360), random(360));
        break;
      }
    }
  }

  void draw_rings(Canvas& c) {
    for (size_t i = 0; i < sizeof(rings) / sizeof(Ring); ++i) {
      if (rings[i].position < H) {
        for (float x = 0; x < W; x += 1) {
          Point p = rings[i].rotation.project(Point(x, rings[i].position));
          plot_aa(c, p.x, p.y, rings[i].color);
        }
      }
    }
  }

  void move_rings(Canvas& c) {
    for (size_t i = 0; i < sizeof(rings) / sizeof(Ring); ++i) {
      if (rings[i].position < H) {
        rings[i].position += rings[i].speed;
      }
    }
  }

  typedef typename Projection<W, H>::Point Point;
  PollingRandomTimer new_ring_timer;
  Ring rings[8];
  unsigned long t = 0;
  ComplementaryColorSequence colors;
};

template <int W>
class ChainWiggle : public Effect {
 public:
  ChainWiggle() : Effect(W), osc(1, 4) {
    for (size_t i = 0; i < (sizeof(dots) / sizeof(Dot)); ++i) {
      dots[i] = Dot(0, i, 0);
    }
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow_lin(c(x, y), 24, 24);
      }
    }

    EVERY_N_MILLIS(1000) { speed *= -1; }

    EVERY_N_MILLIS(4000) {
      gap = osc.get();
      speed = dir(speed) * random(1, 2);
    }

    pull(c, dots, 0, speed);
  };

 private:
  struct Dot {
    Dot() : x(0), y(0), v(0), hue(HUE_RED) {}

    Dot(int x, int y, int v) : x(x), y(y), v(v), hue(HUE_RED) {}

    Dot(const Dot& d) : x(d.x), y(d.y), v(d.v), hue(d.hue) {}

    int x;
    int y;
    int v;
    uint8_t hue;
  };

  void replicate(Canvas& c, int x, int y, const CHSV& color, int n) {
    for (int i = 0; i < W; i += W / n) {
      c(wrap(x + i, W), y) = color;
    }
  }

  int distance(int a, int b, int m) {
    int fwd = wrap(b - a, m);
    int rev = wrap(a - b, m);
    if (fwd <= rev) {
      return fwd;
    }
    return -rev;
  }

  int dir(int v) { return v < 0 ? -1 : 1; }

  void reverse(Dot* dots, int y) { dots[y].v *= -1; }

  void pull(Canvas& c, Dot* dots, int y, int v) {
    dots[y].v = v;
    move(c, dots[y]);
    for (int i = y - 1; i >= 0; --i) {
      drag(c, dots[i + 1], dots[i]);
    }
    for (int i = y + 1; i < H; ++i) {
      drag(c, dots[i - 1], dots[i]);
    }
  }

  void drag(Canvas& c, Dot& leader, Dot& follower) {
    int dest = wrap(follower.x + follower.v, W);
    if (abs(distance(dest, leader.x, W)) > gap) {
      // Move to furthest extent based on leader's prior path
      dest = wrap(leader.x - dir(leader.v) * (abs(leader.v) - 1 + gap), W);
      follower.v = distance(follower.x, dest, W);
      move(c, follower);
      // Move to gap's length behind leader
      dest = wrap(leader.x - dir(leader.v) * gap, W);
      follower.v = distance(follower.x, dest, W);
      move(c, follower);
      // Adjust speed to match leader
      follower.v = leader.v;
    } else {
      move(c, follower);
    }
  }

  void move(Canvas& c, Dot& dot) {
    int dest = wrap(dot.x + dot.v, W);
    for (int i = dot.x;; i = wrap(i + dir(dot.v), W)) {
      replicate(c, i, dot.y, CHSV(dot.hue, 255, 255), replicas);
      if (i == dest) break;
    }
    dot.x = dest;
  }

  Dot dots[H];
  Oscillator osc;
  int gap = 1;
  int speed = 1;
  int replicas = 2;
};

template <int W>
class RingTwist : public Effect {
 public:
  RingTwist() : Effect(W), seed_(CHSV(104, 0, 255)) {
    randomSeed(analogRead(PIN_RANDOM));
    Canvas c(*this);
    for (int x = 0; x < W; x += W / COUNT_) {
      for (int y = 0; y < H; ++y) {
        c(x, y) = seed_;
      }
    }
  }

  bool show_bg() const { return false; }

  const CRGB& get_pixel(int x, int y) const {
    return Effect::get_pixel(wrap(x - pos_[y], W), y);
  }

  void draw_frame() {
    Canvas c(*this);
    canvas_ = &c;
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (c(x, y) != seed_) {
          decay(c(x, y));
        }
      }
    }

    pull(leader_);
    EVERY_N_MILLISECONDS(STOP_TIMER_) {
      if (stop_ < 0) {
        stop_ = pos_[leader_];
      }
    }
  }

 private:
  void decay(CRGB& p) { trail_rainbow(p, 8, 24); }

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
      lead_length_ =
          lead_lengths_[random(1, sizeof(lead_lengths_) / sizeof(int))];
    }
  }

  int distance(int a, int b, int dir) {
    return dir > 0 ? wrap(b - a, W) : wrap(a - b, W);
  }

  void move(int y, int x0, int x1, int dir) {
    pos_[y] = x1;
    for (int x = 0; x < W; x += W / COUNT_) {
      for (int d = 0; d < distance(x0, x1, dir); ++d) {
        add_trail(*canvas_, y, wrap(x - dir, W), dir);
      }
    }
  }

  void add_trail(Canvas& c, int y, int x, int dir) {
    if (c(wrap(x + dir, W), y) == CRGB(0, 0, 0) || c(x, y) == seed_) {
      return;
    }
    add_trail(c, y, wrap(x - dir, W), dir);
    c(x, y) = c(wrap(x + dir, W), y);
    decay(c(x, y));
  }

  const CRGB seed_;
  const int COUNT_ = 2;
  const int STOP_TIMER_ = 15000;
  int pos_[H] = {0};
  int leader_ = 0;
  int dir_ = 1;
  int stop_ = -1;
  int lead_length_ = 1;
  int lead_lengths_[7] = {1, 2, 3, 4, 6, 12, 16};
  Canvas* canvas_ = NULL;
};

/*
template <int W>
class Curves : public Effect {
 public:
  Curves() : Effect(W) {
    fill_gradient<CHSV>(palette_, sizeof(palette_) / sizeof(CHSV),
                        rgb2hsv_approximate(CRGB(6, 4, 47)),
                        rgb2hsv_approximate(CRGB(162, 84, 84)),
                        rgb2hsv_approximate(CRGB(252, 114, 0)), SHORTEST_HUES);
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow_lin(c(x, y), 0, 32);
        c(x, y).h = palette_[c(x, y).v].h;
      }
    }

    for (int i = 0; i < COUNT; ++i) {
      CHSV color = palette_[map8(triwave8(color_offset_), 100, 140)];
      plot(c, num_, 255, (i * W / COUNT + offset_) % W, color);
      plot(c, -num_, 255, (i * W / COUNT + W - 1 - offset_) % W, color);
    }

    num_ = wrap(num_ - 1, 1024);
    color_offset_++;
    offset_ = (offset_ + 1) % W;
  }

 private:
  void plot(Canvas& cv, int n, int d, int c, const CHSV& color) {
    for (int y = 0; y < H; ++y) {
      cv(wrap(y * n / d + c, W), y) = color;
    }
  }

  int COUNT = 4;
  int offset_ = 0;
  uint8_t color_offset_ = 0;
  uint16_t num_ = 1024;
  uint8_t falloff_ = 48;
  CHSVPalette256 palette_;
};
*/

template <uint8_t W, uint8_t HUE>
class TheMatrix : public Effect {
 public:
  TheMatrix() : Effect(W) {
    random16_add_entropy(random());
    memset(pixels_, 0, sizeof(pixels_));
  }

  bool show_bg() const { return true; }

  void draw_frame() {
    Canvas c(*this);
    EVERY_N_MILLIS(125) {
      for (int x = 0; x < W; ++x) {
        fall(x);
        generate(x);
        for (int y = 0; y < H; ++y) {
          c(x, y) =
              blend(CHSV(HUE + (3 * y), 255, 40),
                    CHSV(HUE + (3 * (H - 1 - y)), 255, 255), pixels_[x][y]);
        }
      }
    }
  }

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

template <int W>
class StarsFade : public Effect {
 private:
  uint8_t hue_;

 public:
  StarsFade() : Effect(W), hue_(0) {}

  bool show_bg() const { return true; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (c(x, y) != CRGB(0, 0, 0)) {
          trail_rainbow(c(x, y), 16, 32);
        } else if (random8() < 3) {
          c(x, y) = CRGB(CHSV(hue_, 255, 255));
        }
      }
    }
    hue_++;
  }
};

template <int W, int SPREAD>
class Spiral : public Effect {
 private:
  CHSVPalette16 palette_;

 public:
  Spiral() : Effect(W) {
    fill_rainbow(palette_.entries, 16, 0, 256 / 16);
    draw_frame();
  }

  bool show_bg() const { return true; }

  FASTRUN void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if ((x + y) % (SPREAD + 1) == 0) {
          c(x, y) = palette_[(x + y) % 16];
        } else {
          c(x, y) = CHSV(0, 0, 0);
        }
      }
    }
  }
};

template <int W>
class WaveTrails : public Effect {
 public:
  WaveTrails() : Effect(W) {}

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow(c(x, y), 5, 24);
      }
      c(x, beatsin16(5, H / 15, H * 14 / 15, 0, map(x, 0, W - 1, 0, 65535))) =
          CHSV(HUE_BLUE, 100, 255);
      c(x, beatsin16(5, 0, H - 1, 0, map(W - 1 - x, 0, W - 1, 0, 65535))) =
          CHSV(HUE_GREEN, 100, 255);
      c(x, beatsin16(5, 0, H - 1, 32767, map(W - 1 - x, 0, W - 1, 0, 65535))) =
          CHSV(HUE_AQUA, 100, 255);
    }
  }
};

template <uint8_t W>
class RingTrails : public Effect {
 public:
  RingTrails() : Effect(W), dot(0) {
    fill_gradient<CHSV>(palette, sizeof(palette) / sizeof(CHSV),
                        rgb2hsv_approximate(CRGB(6, 4, 47)),
                        rgb2hsv_approximate(CRGB(162, 84, 84)),
                        rgb2hsv_approximate(CRGB(252, 114, 0)), SHORTEST_HUES);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        ttl[x][y] = 0;
      }
    }
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow_lin(c(x, y), 10, 24);
      }
    }

    for (float x = 0; x < W; x += 1) {
      if (static_cast<int>(x) % 2 == 0) {
//        return;
      }
      Point p = projection.project(x, H / 2);
      CHSV color;
      color = CHSV(HUE_YELLOW, 0, 255);
      plot_aa(c, p.x, p.y, color);
    }

    uint8_t dl = beatsin8(1, 1, 3);
    uint8_t dg = beatsin8(2, 1, 3, 16384);
    uint8_t dp = beatsin8(3, 1, 3, 32768);
    projection.rotate(dl, dg, dp);
    dot = (dot + 1) % W;
  }

 private:
  typedef typename Projection<W, H>::Point Point;
  Projection<W, H> projection;
  uint8_t hue = HUE_RED;
  CHSVPalette256 palette;
  int ttl[W][H];
  int dot;
};

template <uint8_t W>
class Kaleidoscope : public Effect {
 public:
  Kaleidoscope() : Effect(W) {
    fill_rainbow(palette_.entries, 16, 0, 256 / 16);
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        trail_rainbow_lin(c(x, y), 0, 24);
      }
    }

    for (uint8_t i = 0; i < counts_[count_]; ++i) {
      uint8_t x = addmod8(i * offset_, tx_, W);
      uint8_t y = ty_;
      uint8_t x2 = W - 1 - x;
      uint8_t y2 = H - 1 - y;
      c(x, y) =
          palette_[addmod8(palette_offset_, i * (16 / counts_[count_]), 16)];
      c(x2, y2) = palette_[addmod8(palette_offset_,
                                   16 - (i * (16 / counts_[count_])), 16)];
    }

    tx_ = addmod8(tx_, 2, W);
    ty_ = addmod8(ty_, inc_y_, H);
    if (ty_ == H - 1 || ty_ == 0) {
      inc_y_ = 0 - inc_y_;
      tx_ = addmod8(tx_, 2, W);
      EVERY_N_MILLISECONDS(3000) {
        count_ = addmod8(count_, 1, sizeof(counts_));
        offset_ = W / counts_[count_];
      }
    }
    palette_offset_ = addmod8(palette_offset_, 1, 16);
  }

 private:
  uint8_t counts_[9] = {1, 2, 3, 4, 6, 8, 12, 16, 20};
  uint8_t count_ = 0;
  uint8_t offset_ = W / counts_[count_];
  uint8_t tx_ = 0;
  uint8_t ty_ = 0;
  int inc_y_ = 1;
  CHSVPalette16 palette_;
  uint8_t palette_offset_ = 0;
};

template <uint8_t W>
class RingRotate : public Effect {
 public:
  RingRotate() : Effect(W) { 
    fill_rainbow(pal.entries, 256, HUE_RED, 1); 
  }
  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        c(x, y) = CHSV(0, 0, 0);
      }
    }

    const float rings[] = {4.5, 9.5, 15.5};
    for (size_t i = 0; i < sizeof(rings) / sizeof(float); ++i) {
      float y = rings[i];
      for (float x = 0; x < W; ++x) {
        Point p(x, y);
        int dir = 1;
        if (i % 2 == 0) {
          dir = 0;
          p = projection.project(p);
        } else {
          p = projection2.project(p);
        }
        CHSV color = pal[(static_cast<uint8_t>(x) + c_off * dir) * 255 / W];
        plot_aa(c, p.x, p.y, color);
      }
    }

    EVERY_N_SECONDS(10) {
      r_off = (r_off + 1) % (sizeof(rotations) / sizeof(uint16_t) / 3 / 2);
    }

    projection.rotate(rotations[r_off][0][0], rotations[r_off][0][1],
                      rotations[r_off][0][2]);

    projection2.rotate(rotations[r_off][1][0], rotations[r_off][1][1],
                       rotations[r_off][1][2]);

    c_off += 2;
  }

 private:
  typedef typename Projection<W, H>::Point Point;
  Projection<W, H> projection;
  Projection<W, H> projection2;

  CHSVPalette256 pal;
  uint8_t c_off = 0;
  bool bg = 0;

  uint16_t rotations[4][2][3] = {
      {{0, 5, 0}, {0, 5, 0}},
      {{0, 1, 3}, {0, 2, 5}},
      {{0, 5, 0}, {0, 355, 0}},
      {{0, 355, 0}, {0, 3, 0}},
  };

  uint8_t r_off = 0;
};

template <uint8_t W, uint8_t HUE, uint8_t BURNRATE>
class Burnout : public Effect {
 public:
  Burnout() : Effect(W), timer_(30), burn_idx_(0) {
    random16_add_entropy(random());
    memset8(pixels_, INIT, sizeof(pixels_));
    fill_seq(burn_, W * H);
    shuffle(burn_, W * H);
  }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (pixels_[x][y] & BURN) {
          c(x, y) = CHSV(HUE + (3 * y), 255, 255);
        } else if (pixels_[x][y] == INIT) {
          c(x, y) = CHSV(HUE + (3 * y), 255, 40);
        } else {
          c(x, y) = CHSV(0, 0, 0);
        }
      }
      fall(x);
    }
    if (timer_) {
      burnout();
    }
  }

  bool show_bg() const { return true; }

 private:
  enum State {
    NONE,
    INIT,
    BURN = 0x80,
  };

  void burnout() {
    if (burn_idx_ < W * H) {
      *((uint8_t*)pixels_ + burn_[burn_idx_++]) = BURN;
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

  void fill_seq(short* arr, int len) {
    for (int i = 0; i < len; ++i) {
      arr[i] = i;
    }
  }

  void shuffle(short* arr, int len) {
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

template <uint8_t W, uint8_t COOL, uint8_t SPARK>
class Fire : public Effect {
 public:
  Fire() : Effect(W) { random16_add_entropy(random()); }

  bool show_bg() const { return false; }

  void draw_frame() {
    EVERY_N_MILLIS(125) {
      Canvas c(*this);
      for (int x = 0; x < W; ++x) {
        cool(x);
        rise(x);
        spark(x);
        for (int y = 0; y < H; ++y) {
          c(x, y) = rgb2hsv_approximate(HeatColor(heat_[x][y]));
        }
      }
    }
  }

 private:
  inline void cool(uint8_t x) {
    for (uint8_t y = 0; y < H; ++y) {
      heat_[x][y] = qsub8(heat_[x][y], random(0, ((COOL * 10) / H) + 2));
    }
  }

  inline void rise(uint8_t x) {
    for (uint8_t y = H - 1; y >= 2; --y) {
      heat_[x][H - 1 - y] = (heat_[x][H - 1 - y + 1] + heat_[x][H - 1 - y + 2] +
                             heat_[x][H - 1 - y + 2]) /
                            3;
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

template <uint8_t W>
class DotTrails : public Effect {
 public:
  DotTrails() : Effect(W) {
    random16_add_entropy(random());
    for (int i = 0; i < H - 1; ++i) {
      bool rev = i % 2 == 0;
      uint8_t x = random8() % W;
      dots[i] = Dot(x, rev, (random8() % 30) + 20);
    }
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (dots[y].x == x) {
          if (dots[y].drawn) {
            dots[y].drawn = false;
          } else {
            c(x, y) = CHSV(hue - 12, 0, 255);
            move(dots[y]);
          }
        } else {
          trail_rainbow(c(x, y), 12, 12);
        }
      }
    }
  }

 private:
  struct Dot {
    Dot() : x(0), rev(0), ttl(0), drawn(0) {}

    Dot(uint8_t x, bool rev, int ttl) : x(x), rev(rev), ttl(ttl), drawn(0) {}

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

  Dot dots[H];
  uint8_t hue = HUE_RED;
};

template <uint8_t W>
class Spinner : public Effect {
 public:
  Spinner() : Effect(W), spin_timer_(62), pos_(0) {
    memset(palette1_, 0, sizeof(palette1_));
    memset(palette2_, 0, sizeof(palette1_));
    palette1_[0] = palette1_[8] = CHSV(HUE_RED, 255, 255);
    palette1_[1] = palette1_[9] = CHSV(HUE_YELLOW, 255, 255);
    palette2_[0] = palette2_[8] = CHSV(HUE_AQUA, 255, 255);
    palette2_[1] = palette2_[9] = CHSV(HUE_BLUE, 255, 255);
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    for (int x = 0; x < W; ++x) {
      for (int y = 0; y < H; ++y) {
        if (!swap) {
          if (y % (p[i] * 2) < p[i]) {
            c(x, y) = palette1_[addmod8(x, pos_, 16)];
          } else {
            c(x, y) = palette2_[(uint8_t)(pos_ - x) % 16];
          }
        } else {
          if (y % (p[i] * 2) < p[i]) {
            c(x, y) = palette2_[(uint8_t)(pos_ - x) % 16];
          } else {
            c(x, y) = palette1_[addmod8(x, pos_, 16)];
          }
        }
      }
      if (spin_timer_) {
        pos_ = addmod8(pos_, 1, W);
      }
    }
    EVERY_N_SECONDS(5) {
      if ((!swap && (i == sizeof(p) - 1)) || (swap && (i == 0))) {
        swap = !swap;
      } else {
        i += swap ? -1 : 1;
      }
    }
  }

 private:
  CHSVPalette16 palette1_;
  CHSVPalette16 palette2_;
  CEveryNMillis spin_timer_;
  uint8_t pos_;
  bool swap = false;
  uint8_t p[5] = {20, 10, 4, 2, 1};
  uint8_t i = 0;
};
