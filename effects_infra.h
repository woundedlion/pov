#pragma once

#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "3dmath.h"
#include <FastLED.h>

inline int XY(int x, int y) { return x * H + y; }

template <uint8_t W, uint8_t H>
struct Dot {
  Dot(float x_cartesian, float y_cartesian, float z_cartesian, const CHSV& color) :
    lambda(atan2(y_cartesian, x_cartesian)),
    phi(asinf(z_cartesian)),
    color(color)
  {
    x = lambda_to_x(lambda);
    y = phi_to_y(phi);
  }

  Dot(int x, int y, const CHSV& color)
    : x(x), y(y), lambda(x_to_lambda(x)), phi(y_to_phi(y)), color(color) {}

  Dot(uint8_t x, uint8_t y, const CHSV& color)
    : x(x), y(y), lambda(x_to_lambda(x)), phi(y_to_phi(y)), color(color) {}

  Dot(float x, float y, const CHSV& color)
    : x(x), y(y), lambda(x_to_lambda(x)), phi(y_to_phi(y)), color(color) {}

  Dot(const Dot& p)
    : x(p.x), y(p.y), lambda(p.lambda), phi(p.phi), color(p.color) {}

  int xi() { return static_cast<int>(x + 0.5f) % W; }
  int yi() { return static_cast<int>(y + 0.5f); }

  float x_cartesian() const { return cosf(phi) * cosf(lambda); }
  float y_cartesian() const { return cosf(phi) * sinf(lambda); }
  float z_cartesian() const { return sinf(phi); }

  float x;
  float y;
  float lambda;  // longitude
  float phi;     // latitude
  CHSV color;

private:

  float x_to_lambda(float x) { return fmod(x * tau / W, tau) - pi; }
  float lambda_to_x(float lambda) { return fmod((lambda + pi) * W / tau, W); }
  float y_to_phi(float y) { return (H - 1 - y) * pi / (H - 1) - (pi / 2); }
  float phi_to_y(float phi) { return (H - 1) - ((phi + (pi / 2)) * (H - 1) / pi); }
};

template <uint8_t W, uint8_t H>
Dot<W, H> move(const Dot<W, H>& origin, float azimuth, float distance) {
  Dot<W, H> r(origin);
  float a = degrees_to_radians(azimuth);
  float d = distance / (W / 2); // angular distance d/R
  r.phi = asinf(sinf(origin.phi) * cosf(d) +
    cosf(origin.phi) * sinf(d) * cosf(a));
  r.lambda = origin.lambda + atan2f(sinf(a) * sinf(d) * cosf(origin.phi),
    cosf(d) - sinf(origin.phi) * sinf(r.phi));
  r.x = lambda_to_x(r.lambda);
  r.y = phi_to_y(r.phi);
  return r;
}

template<uint8_t W, uint8_t H>
Dot<W, H> rotate(const Dot<W, H>& d, const Quaternion& q) {
  Quaternion p(0, d.x_cartesian(), d.y_cartesian(), d.z_cartesian());
  auto r = q.inverse() * p * q;
  return Dot<W, H>(r.v.i, r.v.j, r.v.k, d.color);
}

template <uint8_t W, uint8_t H>
class Sprite {
public:

  Sprite() : orientation(1, 0, 0, 0) {}

  const Quaternion& get_orientation() const { return orientation; }

  void build(float x, float y, const CHSV& c) {
    dots.emplace_back(Dot<W, H>(x, y, c));
  }

  Sprite& rotate(const Vector& axis, float angle, float tween_step = 1) {
    auto rotation = make_rotation(axis, angle);
    Quaternion start = orientation;
    orientation *= rotation;
    for (float t = 0 + tween_step; t < 1; t += tween_step) {
      auto tween = slerp(start, orientation, t);
      for (const auto& d : dots) {
        painted.emplace_back(::rotate(d, tween));
      }
    }
    return *this;
  }

  void move(float azimuth, float distance) {
    for (auto& d : dots) {
      // TODO slerp and paint
      d = move(d, azimuth, distance);
    }
  }

  std::vector<Dot<W, H>> get_paint() {
    for (const auto& d : dots) {
      painted.emplace_back(::rotate(d, orientation));
    }
    std::vector<Dot<W, H>> r(painted);
    painted.clear();
    return r;
  }

private:

  std::vector<Dot<W, H>> dots;
  std::vector<Dot<W, H>> painted;
  Quaternion orientation;
};

///////////////////////////////////////////////////////////////////////////////

class Canvas;

class Effect {
  friend class Canvas;

public:
  Effect(int W) : width_(W) {
    bufs_[0] = new CHSV[W * H];
    memset(bufs_[0], 0, sizeof(CHSV) * W * H);
    bufs_[1] = new CHSV[W * H];
    memset(bufs_[1], 0, sizeof(CHSV) * W * H);
  }

  virtual ~Effect() {
    delete[] bufs_[0];
    delete[] bufs_[1];
  };

  virtual void draw_frame() = 0;
  virtual bool show_bg() const = 0;

  virtual const CHSV& get_pixel(int x, int y) const {
    return bufs_[prev_][XY(x, y)];
  }

  inline int width() const { return width_; }
  inline bool buffer_free() const { return prev_ == next_; }
  inline void advance_display() { prev_ = next_; }
  inline void advance_buffer() {
    cur_ = cur_ ? 0 : 1;
    memcpy(bufs_[cur_], bufs_[prev_], sizeof(CHSV) * width_ * H);
  }

  inline void queue_frame() { next_ = cur_; }

private:
  volatile int prev_ = 0, cur_ = 0, next_ = 0;
  int width_;
  CHSV* bufs_[2];
};

class Canvas {
public:
  Canvas(Effect& effect) : effect_(effect) {
    while (!effect_.buffer_free()) {}
    effect_.advance_buffer();
  }

  ~Canvas() { effect_.queue_frame(); }

  inline CHSV& operator()(int x, int y) {
    return operator()(XY(x, y));
  }

  inline CHSV& operator()(int xy) {
    return effect_.bufs_[effect_.cur_][xy];
  }

  const int width() { return effect_.width(); }

private:
  Effect& effect_;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Filter {
public:
  Filter() : next(nullptr) {}

  virtual Filter& chain(Filter& filter) {
    next = &filter;
    return *this;
  }

  virtual void plot(Canvas& canvas, float x, float y, const CHSV& c) = 0;

protected:

  void pass(Canvas& canvas, float x, float y, const CHSV& c) {
    if (next == nullptr) {
      Serial.printf("Canvas: (%f, %f) <- (%d, %d, %d)\n", x, y, c.h, c.s, c.v);
      canvas(x, y) = c;
    }
    else {
      next->plot(canvas, x, y, c);
    }
  }

  Filter* next;
};

template <uint8_t W, uint8_t H>
class FilterDecay : public Filter {
public:
  FilterDecay(int lifetime, const TProgmemRGBGradientPalette_bytes& palette_def) : lifetime(lifetime), palette(palette_def) {}

  void age(Canvas& canvas) {
    for (auto a = ages.begin(); a != ages.end(); ++a) {
      CHSV& pixel = canvas(a->first);
      int age = --(a->second);
      if (age == 0) {
        pixel = CHSV(0, 0, 0);
        a = ages.erase(a);
      }
      else {
        pixel = rgb2hsv_approximate(ColorFromPalette(palette, age * 255 / lifetime));
      }
    }
  }

  void plot(Canvas& canvas, float x, float y, const CHSV& c) {
    Serial.printf("Decay: (%f, %f) <- (%d, %d, %d)\n", x, y, c.h, c.s, c.v);
    auto a = std::find_if(ages.begin(), ages.end(),
      [=](const auto& a) -> bool {
        return a.first == XY(x, y);
      }
    );
    if (a == ages.end()) {
      ages.emplace_back(std::make_pair(XY(x, y), lifetime));
    }
    else {
      a->second = lifetime;
    }
    pass(canvas, x, y, c);
  }

private:
  typedef std::list<std::pair<int, int>> Ages;

  Ages ages;
  int lifetime;
  CRGBPalette256 palette;
};

template <uint8_t W, uint8_t H>
class FilterAntiAlias : public Filter {
public:
  FilterAntiAlias() {}

  void plot(Canvas& canvas, float x, float y, const CHSV& c) {
    Serial.printf("AA: (%f, %f) <- (%d, %d, %d)\n", x, y, c.h, c.s, c.v);

    float x_i = 0;
    float x_m = modf(x, &x_i);
    float y_i = 0;
    float y_m = modf(y, &y_i);
    constexpr int FULL = 255;

    uint8_t v = static_cast<uint8_t>(((1 - x_m) * (1 - y_m) * FULL) + 0.5f);
    pass(canvas, x_i, y_i,
      CHSV(c.h, c.s, dim8_lin(std::max(canvas(x_i, y_i).v, v))));

    v = static_cast<uint8_t>((x_m * (1 - y_m) * FULL) + 0.5f);
    pass(canvas, (static_cast<int>(x_i + 1)) % W, y_i,
      CHSV(c.h, c.s, dim8_lin(std::max(canvas(static_cast<int>(x_i + 1) % W, y_i).v, v))));

    if (y_i < H - 1) {
      v = static_cast<uint8_t>(((1 - x_m) * y_m * FULL) + 0.5f);
      pass(canvas, x_i, y_i + 1,
        CHSV(c.h, c.s, dim8_lin(std::max(canvas(x_i, y_i + 1).v, v))));

      v = static_cast<uint8_t>((x_m * y_m * FULL) + 0.5f);
      pass(canvas, static_cast<int>(x_i + 1) % W, y_i + 1,
        CHSV(c.h, c.s, dim8_lin(std::max(canvas(static_cast<int>(x_i + 1) % W, y_i + 1).v, v))));
    }
  }
};

template <typename Sprite>
void paint(Canvas& c, Filter& f, Sprite& sprite) {
  auto dots = sprite.get_paint();
  for (auto& d : dots) {
    f.plot(c, d.x, d.y, d.color);
  }
}
