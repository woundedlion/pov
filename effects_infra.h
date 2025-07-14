#pragma once

#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <cassert>
#include <array>
#include <memory>
#include <FastLED.h>
#include "3dmath.h"

typedef std::vector<Vector> VertexList;
typedef std::vector<std::vector<unsigned int>> AdjacencyList;
typedef std::function<CHSV(const Vector&, double)> ColorFn;
typedef std::function<double (double)> EasingFn;
typedef std::function<Vector (double)> PlotFn;
typedef std::function<double (double)> ShiftFn;
typedef std::function<CHSV(const CHSV&, const CHSV&)> BlendFn;
typedef std::function<void (double)> SpriteFn;
typedef std::function<void ()> TimerFn;
typedef std::function<double (double)> MutateFn;
typedef std::function<double(double)> WaveFn;


// inline int XY(int x, int y) { return x * H + y; }

static const int FPS = 16;
static const Vector X_AXIS(1, 0, 0);
static const Vector Y_AXIS(0, 1, 0);
static const Vector Z_AXIS(0, 0, 1);

struct Dodecahedron {
  VertexList vertices = {
    {1, 1, 1},       // 0
    {1, -1, 1},      // 1
    {1, 1, -1},      // 2
    {1, -1, -1},     // 3
    {-1, 1, 1},      // 4
    {-1, -1, 1},     // 5
    {-1, 1, -1},     // 6
    {-1, -1, -1},    // 7

    {0, g, 1 / g},   //8
    {0, -g, 1 / g},  // 9
    {0, g, -1 / g},  // 10
    {0, -g, -1 / g}, // 11

    {1 / g, 0, g},   // 12
    {1 / g, 0, -g},  // 13
    {-1 / g, 0, g},  // 14
    {-1 / g, 0, -g}, // 15

    {g, 1 / g, 0},   // 16
    {g, -1 / g, 0},  // 17
    {-g, 1 / g, 0},  // 18
    {-g, -1 / g, 0}, // 19


  };

  AdjacencyList edges = {
    {8, 12, 16},  // 0
    {9, 12, 17},  // 1
    {10, 13, 16}, // 2
    {11, 13, 17}, // 3
    {8, 14, 18},  // 4
    {9, 14, 19},  // 5
    {10, 15, 18}, // 6
    {11, 15, 19}, // 7
    {0, 4, 10},   // 8
    {1, 5, 11},   // 9
    {2, 6, 8},    // 10
    {3, 7, 9},    // 11
    {0, 1, 14},   // 12
    {2, 3, 15},   // 13
    {4, 5, 12},   // 14
    {6, 7, 13},   // 15
    {0, 2, 17},   // 16
    {1, 3, 16},   // 17
    {4, 6, 19},   // 18
    {5, 7, 18},   // 19
  };

  AdjacencyList euler_path = {
    {8, 12, 16},  // 0
    {9, 12, 17},  // 1
    {10, 13, 16}, // 2
    {11, 13, 17}, // 3
    {8, 14, 18},  // 4
    {9, 14, 19},  // 5
    {10, 15, 18}, // 6
    {11, 15, 19}, // 7
    {10},   // 8
    {11},   // 9
    {8},    // 10
    {9},    // 11
    {14},   // 12
    {15},   // 13
    {12},   // 14
    {13},   // 15
    {17},   // 16
    {16},   // 17
    {19},   // 18
    {18},   // 19
  };

  void dump() {
    Serial.printf("Dodecahedron\n");
    Serial.printf("vertices [%d] = \n", vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
      Serial.printf("%d: [%f, %f, %f]\n", i, vertices[i].i, vertices[i].j, vertices[i].k);
    }
    Serial.printf("euler path [%d] = \n", euler_path.size());
    for (size_t i = 0; i < euler_path.size(); ++i) {
      Serial.printf("%d: [", i);
      for (size_t j = 0; j < euler_path[i].size(); ++j) {
        Serial.printf("%d, ", euler_path[i][j]);
      }
      Serial.printf("]\n", i);
    }
  }
};

///////////////////////////////////////////////////////////////////////////////

struct Dot {
  Dot(const Vector& v, const CHSV& color) :
    position(v),
    color(color)
  {}

  Dot(const Dot& p)
    : position(p.position), color(p.color) {}

  Vector position;
  CHSV color;
};

struct Pixel {
  CHSV color = { 0, 0, 0 };
};

typedef std::vector<Dot> Dots;
typedef std::unordered_map<int, Pixel> Pixels;

CHSV dim(const CHSV& c, double s) {
  return CHSV(c.h, c.s, dim8_lin(s * c.v));
}

CHSV blend_over(const CHSV& c1, const CHSV& c2) {
  return CHSV(c2.h, c2.s, c2.v);
}

CHSV blend_over_max(const CHSV& c1, const CHSV& c2) {
  return CHSV(c2.h, c2.s, std::max(c1.v, c2.v));
}

CHSV blend_over_add(const CHSV& c1, const CHSV& c2) {
  return CHSV(c2.h, c2.s, qadd8(c1.v, c2.v));
}

//////////////////////////////////////////////////////////////////////////////////////////

double lerp(double from, double to, double t) {
  return ((to - from) * t) + from;
}

WaveFn sin_wave(double from, double to, double freq, double phase) {
  return [=](double t) -> double {
    auto w = (sin(freq * t * 2 * PI - (PI / 2) + PI - (2 * phase)) + 1) / 2;
    return lerp(from, to, w);
  };
}

WaveFn tri_wave(double from, double to, double freq, double phase) {
  return [=](double t) -> double {
    double w;
    if (t < 0.5) {
      w = 2 * t;
    } else {
      w = 2 - 2 * t;
    }
    return lerp(from, to, w);
  };
}

WaveFn square_wave(double from, double to, double freq, double dutyCycle, double phase) {
  return [=](double t) -> double {
    if (fmod(t * freq + phase, 1) < dutyCycle) {
      return to;
    }
    return from;
  };
}

///////////////////////////////////////////////////////////////////////////////

class Orientation {
public:
  Orientation() : orientations{ Quaternion() } {}
  Orientation(const Quaternion& q) {
    set(q);
  }

  int length() const { return orientations.size(); }
 
  Vector orient(const Vector& v) const {
    return rotate(Vector(v).normalize(), orientations.back());
  }
  
  Vector orient(const Vector& v, int i) const {
    return rotate(Vector(v).normalize(), orientations[i]);
  }
  
  Vector unorient(const Vector& v) const {
    return rotate(Vector(v).normalize(), orientations.back().inverse());
  }
  
  Vector unorient(const Vector& v, int i) const {
    return rotate(Vector(v).normalize(), orientations[i].inverse());
  }

  VertexList orient(const VertexList& vertices) const {
    VertexList r;
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(r), 
      [this](auto& v) {
        return orient(v);
      });
    return r;
  }
  VertexList orient(const VertexList& vertices, int i) const {
    VertexList r;
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(r),
      [=](auto& v) {
        return orient(v, i);
      });
    return r;
  }

  void clear() {
    orientations.clear();
  }

  const Quaternion& get() const {
    return orientations[orientations.size() - 1];
  }

  const Quaternion& get(int i) const {
    return orientations[i];
  }

  Orientation& set(const Quaternion& q) {
    orientations.assign({q});
    return *this;
  }

  Orientation& push(const Quaternion& q) {
    orientations.push_back(q);
    return *this;
  }

  Orientation& collapse() {
    if (orientations.size() > 1) {
      orientations.assign({ orientations.back() });
    }
    return *this;
  }

private:
  std::vector<Quaternion> orientations;
};

template <int W>
Dots draw_line(const Vector& v1, const Vector& v2, ColorFn color, bool long_way = false);

template <int W>
class Path {
public:
  Path() {}

  Path& append_line(const Vector& v1, Vector& v2, bool long_way = false) {
    if (points.size() > 0) {
      points.pop_back(); // Overlap previous segment
    }
    Dots seg = draw_line<W>(v1, v2, [](auto& , auto) { return CHSV(); }, long_way);
    std::transform(seg.begin(), seg.end(), std::back_inserter(points), 
      [](auto& d) { return d.position; });
    return *this;
  }

  Path& append_segment(PlotFn plot, double domain, double samples, EasingFn easing) {
    if (points.size() > 0) {
      points.pop_back(); // Overlap previous segment
    }
    for (double t = 0; t < samples; t++) {
      points.push_back(plot(easing(t / samples) * domain));
    }
    return *this;
  }

  Vector get_point(double t) const {
    assert(static_cast<size_t>(t * points.size()) < points.size());
    return points[static_cast<int>(t * points.size())];
  }

  size_t num_points() const { return points.size(); }

private:
  std::vector<Vector> points;
};

template <int W>
Dots draw_path(const Path<W>& path, ColorFn color) {
  Dots dots;
  size_t samples = path.num_points();
  for (size_t i = 0; i < samples; ++i) {
    auto v = path.get_point(static_cast<double>(i) / samples);
    dots.push_back(Dot(v, color(v, i / (samples - 1))));
  }
  return dots;
}

template <typename Poly>
void bisect(Poly& poly, const Orientation& orientation, const Vector& normal) {
  auto num_edges = poly.euler_path.size();
  for (size_t i = 0; i < num_edges; ++i) {
    for (auto j = poly.euler_path[i].begin(); j != poly.euler_path[i].end();) {
      Vector a(orientation.orient(poly.vertices[i]));
      Vector b(orientation.orient(poly.vertices[*j]));
      if (intersects_plane(a, b, normal)) {
        auto p = split_point(intersection(a, b, normal), normal);
        poly.vertices.push_back(orientation.unorient(p[0]));
        poly.vertices.push_back(orientation.unorient(p[1]));
        if (is_over(a, normal)) {
          poly.euler_path.push_back({ i });
          poly.euler_path.push_back({ *j });
        }
        else {
          poly.euler_path.push_back({ *j });
          poly.euler_path.push_back({ i });
        }
        j = poly.euler_path[i].erase(j);
      } else {
        ++j;
      }
    }
  }
}

CHSV distance_gradient(const Vector& v, const Vector& normal, CRGBPalette256 p1, CRGBPalette256 p2) {
  auto d = dot(v, normal);
  if (d > 0) {
    return rgb2hsv_approximate(p1[static_cast<int>(d * 255)]);
  } else {
    return rgb2hsv_approximate(p2[static_cast<int>(-d * 255)]);
  }
}
///////////////////////////////////////////////////////////////////////////////

template <int W>
Dots draw_line(const Vector& v1, const Vector& v2, ColorFn color, bool long_way /* = false*/) {
  Dots dots;
  Vector u(v1);
  Vector v(v2);
  u.normalize();
  v.normalize();
  double a = angle_between(u, v);
  Vector w = cross(v, u);
  if (long_way) {
    a = 2 * PI - a;
    w = (-w).normalize();
  }
  v = cross(u, w).normalize();

  double step = 2 * PI / W;
  for (double t = 0; t < a; t += step) {
    Vector vi(
      u.i * cos(t) + v.i * sin(t),
      u.j * cos(t) + v.j * sin(t),
      u.k * cos(t) + v.k * sin(t)
    );
    dots.emplace_back(Dot(vi, color(vi, t)));
  }
  return dots;
}

Dots draw_vertices(const VertexList& vertices, ColorFn color_fn) {
  Dots dots;
  for (Vector v : vertices) {
    dots.emplace_back(Dot(v.normalize(), color_fn(v, 0)));
  }
  return dots;
}

template <int W>
Dots draw_polyhedron(const VertexList& vertices, const AdjacencyList& edges, ColorFn color_fn) {
  Dots dots;

  for (size_t i = 0; i < edges.size(); ++i) {
    Vector a(vertices[i]);
    for (auto j : edges[i]) {
      Vector b(vertices[j]);
      auto seg = draw_line<W>(a, b, color_fn);
      dots.insert(dots.end(), seg.begin(), seg.end());
    }
  }

  return dots;
}

Vector calc_ring_point(double a, double radius, const Vector& u, const Vector& v, const Vector& w) {
  auto d = sqrt(pow(1 - radius, 2));
  return Vector(
    d * v.i + radius * u.i * cos(a) + radius * w.i * sin(a),
    d * v.j + radius * u.j * cos(a) + radius * w.j * sin(a),
    d * v.k + radius * u.k * cos(a) + radius * w.k * sin(a)
  ).normalize();
}

Vector fn_point(ShiftFn f, const Vector& normal, double radius, double angle) {
  Vector v(normal);
  if (radius > 1) {
    v = -v;
  }
  Vector u;
  if (v.i == 0 && v.j == 0) {
    u = cross(v, X_AXIS).normalize();
  }
  else {
    u = cross(v, Z_AXIS).normalize();
  }
  Vector w(cross(v, u));
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }

  auto vi = calc_ring_point(angle, radius, u, v, w);
  auto vp = calc_ring_point(angle, 1, u, v, w);
  auto axis = cross(v, vp).normalize();
  auto shift = make_rotation(axis, f(angle * PI / 2));
  return rotate(vi, shift);
};

template <int W>
Dots draw_fn(const Orientation& orientation, const Vector& normal, double radius, ShiftFn shift_fn, ColorFn color_fn) {
  Vector v(orientation.orient(normal));
  if (radius > 1) {
    v = -v;
  }
  Vector u;
  if (v.i == 0 && v.j == 0) {
    u = cross(v, X_AXIS).normalize();
  }
  else {
    u = cross(v, Z_AXIS).normalize();
  }
  Vector w(cross(v, u));
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }

  bool first = true;
  Vector start, from, to;
  double step = 1.0 / W;
  Dots dots;
  for (double t = 0; t < 1; t += step) {
    auto vi = calc_ring_point(t * 2 * PI, radius, u, v, w);
    auto vp = calc_ring_point(t * 2 * PI, 1, u, v, w);
    Vector axis = cross(v, vp).normalize();
    auto shift = make_rotation(axis, shift_fn(t));
    auto to = rotate(vi, shift);
    if (first) {
      dots.emplace_back(Dot(to, color_fn(to, t)));
      first = false;
      start = to;
    }
    else {
      auto seg = draw_line<W>(from, to, color_fn);
      dots.insert(dots.begin(), seg.begin(), seg.end());
    }
    from = to;
  }
  auto seg = draw_line<W>(from, start, color_fn);
  dots.insert(dots.begin(), seg.begin(), seg.end());

  return dots;
};


Vector ring_point(const Vector& normal, double radius, double angle, double phase = 0) {
  Vector v(normal);
  if (radius > 1) {
    v = -v;
  }
  Vector u;
  if (v.i == 0 && v.j == 0) {
    u = cross(v, X_AXIS).normalize();
  }
  else {
    u = cross(v, Z_AXIS).normalize();
  }
  Vector w(cross(v, u));
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }
  return calc_ring_point(angle + phase, radius, u, v, w);
};

template<int W>
Dots draw_ring(const Orientation& orientation, const Vector& normal, double radius, ColorFn color_fn, double phase = 0) {
  Dots dots;
  Vector v(orientation.orient(normal));
  if (radius > 1) {
    v = -v;
  }
  Vector u;
  if (v.i == 0 && v.j == 0) {
    u = cross(v, X_AXIS).normalize();
  }
  else {
    u = cross(v, Z_AXIS).normalize();
  }
  Vector w(cross(v, u));
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }

  double step = 2 * PI / W;
  for (double a = 0; a < 2 * PI; a += step) {
    auto vi = calc_ring_point(fmod((a + phase), (2 * PI)), radius, u, v, w);
    dots.emplace_back(Dot(vi, color_fn(vi, a / (2 * PI))));
  }

  return dots;
};

///////////////////////////////////////////////////////////////////////////////

template <int W>
class Filter {
public:
  Filter() : next(nullptr) {}

  virtual Filter& chain(Filter& filter) {
    next = &filter;
    return *this;
  }

  virtual void plot(Pixels& pixes, double x, double y, 
    const CHSV& c, double age, BlendFn blend_mode) = 0;

protected:
  
  void pass(Pixels& pixels, double x, double y,
    const CHSV& c, double age, BlendFn blend_mode)
  {
    if (next == nullptr) {
      pixels[XY(x, y)] = { blend_mode(pixels[XY(x, y)].color, c) };
    }
    else {
      next->plot(pixels, x, y, c, age, blend_mode);
    }
  }

  Filter* next;
};

template <int W>
class FilterRaw : public Filter<W> {
  void plot(Pixels& pixels, double x, double y, const CHSV& c, double age, BlendFn blend_mode) {
      this->pass(pixels, x, y, c, age, blend_mode);
    }
};

template <int W>
class FilterAntiAlias : public Filter<W> {
public:
  FilterAntiAlias() {}
  
  void plot(Pixels& pixels, double x, double y, const CHSV& c, double age, BlendFn blend_mode) {
    double x_i = 0;
    double x_m = modf(x, &x_i);
    double y_i = 0;
    double y_m = modf(y, &y_i);

    uint8_t v = static_cast<uint8_t>(((1 - x_m) * (1 - y_m) * c.v) + 0.5f);
    this->pass(pixels, x_i, y_i, CHSV(c.h, c.s, v), age, blend_mode);

    v = static_cast<uint8_t>((x_m * (1 - y_m) * c.v) + 0.5f);
    this->pass(pixels, (static_cast<int>(x_i + 1)) % W, y_i, CHSV(c.h, c.s, v), age, blend_mode);

    if (y_i < H - 1) {
      v = static_cast<uint8_t>(((1 - x_m) * y_m * c.v) + 0.5f);
      this->pass(pixels, x_i, y_i + 1, CHSV(c.h, c.s, v), age, blend_mode);

      v = static_cast<uint8_t>((x_m * y_m * c.v) + 0.5f);
      this->pass(pixels, static_cast<int>(x_i + 1) % W, y_i + 1, CHSV(c.h, c.s, v), age, blend_mode);
    }
  }
};

template <int W>
class FilterDecayMask : public Filter<W> {
public:
  FilterDecayMask(int lifetime) : lifetime(lifetime) {}

  void decay() {
    for (auto e = ttls.begin(); e != ttls.end();) {
      if (--e->second < std::numeric_limits<double>::epsilon()) {
        e = ttls.erase(e);
      }
      else {
        e++;
      }
    }
  }

  CHSV mask(int xy, const CHSV& c) {
    CHSV r(c);
    if (ttls.find(xy) != ttls.end()) {
      r.v = dim8_lin(r.v * ttls[xy] / lifetime);
    } else {
      r.v = 0;
    }
    return r;
  }

  void plot(Pixels& pixels, double x, double y, const CHSV& c, double age, BlendFn blend_mode) {
    ttls[XY(x, y)] = std::max(0.0, lifetime - age);;
    this->pass(pixels, x, y, c, age, blend_mode);
  }

private:
  std::unordered_map<int, double> ttls;
  int lifetime;
};


template <int W>
class FilterDecayTrails : public Filter<W> {
public:

  FilterDecayTrails(int lifetime, CRGBPalette256 palette) :
    lifetime(lifetime),
    palette(palette)
  {}

  void decay() {
    for (auto e = ttls.begin(); e != ttls.end();) {
      if (--e->second < std::numeric_limits<double>::epsilon()) {
        e = ttls.erase(e);
        trail_pixels.erase(e->first);
      }
      else {
        trail_pixels[e->first] = rgb2hsv_approximate(
          palette[static_cast<int>(ttls[e->first] * 255.0 / lifetime)]);
        e++;
      }
    }
  }

  void plot(Pixels& pixels, double x, double y, const CHSV& c, double age, BlendFn blend_mode) {
    int xy = XY(x, y);
    ttls[xy] = std::max(0.0, lifetime - age);
    trail_pixels[xy] = rgb2hsv_approximate(
      palette[static_cast<int>(ttls[xy] * 255.0 / lifetime)]);
    for (auto& [xy, color] : trail_pixels) {
      pixels[xy] = { color };
    }
    this->pass(pixels, x, y, c, age, blend_mode);
  }

private:

  std::unordered_map<int, double> ttls;
  std::unordered_map<int, CHSV> trail_pixels;
  int lifetime;
  CRGBPalette256 palette;
};

///////////////////////////////////////////////////////////////////////////////

class ProceduralPalette {
public:

  ProceduralPalette(
    std::array<double, 3> a, 
    std::array<double, 3> b,
    std::array<double, 3> c,
    std::array<double, 3> d) :
    a(a),
    b(b),
    c(c),
    d(d)
  {
  }

  CRGB get(double t) {
      return CRGB(
        255 * (a[0] + b[0] * cos(2 * PI * (c[0] * t + d[0]))),
        255 * (a[1] + b[1] * cos(2 * PI * (c[1] * t + d[1]))),
        255 * (a[2] + b[2] * cos(2 * PI * (c[2] * t + d[2])))
    );
  }

private:

  std::array<double, 3> a;
  std::array<double, 3> b;
  std::array<double, 3> c;
  std::array<double, 3> d;
};

CHSV dotted_brush(const CHSV& color, double freq, double duty_cycle, double phase, double t) {
  return dim(color, square_wave(0, 1, freq, duty_cycle, phase)(t));
}

///////////////////////////////////////////////////////////////////////////////

class Animation {
public:

  Animation(int duration, bool repeat) :
    duration(duration),
    repeat(repeat),
    canceled(false)
  {
  }

  void cancel() { canceled = true; }
  virtual bool done() { return canceled || (duration >= 0 && t >= duration); }

  virtual void step() {
    t++;
    if (done()) {
      if (repeat) {
        t = 0;
      }
    }
  }

protected:

  int duration;
  bool repeat;
  int t = 0;

private:

  bool canceled;
};

struct TimelineEvent {
  TimelineEvent(int start, std::shared_ptr<Animation> animation) :
    start(start),
    animation(animation) {}

  TimelineEvent(const TimelineEvent& e) :
  start(e.start),
  animation(e.animation) {}
  
  int start;
  std::shared_ptr<Animation> animation;
};

class Timeline {
public:

  Timeline() {}

  Timeline& add(double inSecs, std::shared_ptr<Animation> animation) {
    auto start = t + static_cast<int>(std::round(inSecs * FPS));
    for (auto e = events.begin(); e != events.end(); ++e) {
      if (e->start > start) {
        events.insert(e, TimelineEvent(start, animation));
        return *this;
      }
    }
    events.emplace_back(TimelineEvent(start, animation));
    return *this;
  }

  void step() {
    t++;
    Serial.printf("T %d: scanning %d events\n", t, events.size());
    events.erase(std::remove_if(events.begin(), events.end(), [](auto& e) { return e.animation->done(); }), events.end());
    auto snapshot(events);
    for (auto e = snapshot.begin(); e != snapshot.end(); ++e) {
      if (t >= e->start) {
        Serial.println("event step");
        e->animation->step();
      }
    }
  }

  int t = 0;

private:

  std::vector<TimelineEvent> events;
};

class RandomTimer : public Animation {
public:

  RandomTimer(int min, int max, TimerFn f, bool repeat = false) :
    Animation(-1, repeat),
    min(min),
    max(max),
    f(f),
    next(0)
  {
    reset();
  }

  void reset() {
    next = t + static_cast<int>(std::round(hs::rand_int(min, max)));
  }

  void step() {
    Serial.println("step RandomTimer");
    if (t >= next) {
      f();
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
    Animation::step();
  }
  
  private:

    int min;
    int max;
    TimerFn f;
    int next;
};

class PeriodicTimer : public Animation {
public:
  PeriodicTimer(int period, std::function<void()> f, bool repeat = false) :
    Animation(-1, repeat),
    period(period),
    f(f)
  {
    reset();
  }

  void reset() {
    next = t + period;
  }

  void step() {
    if (t >= next) {
      f();
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
    Animation::step();
  }

private:

  int period;
  std::function<void()> f;
  int next;
};

class Transition : public Animation {
public:
  Transition(double& mutant, double to, int duration, EasingFn easing_fn, bool quantized = false, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    to(to),
    easing_fn(easing_fn),
    quantized(quantized)
  {
  }

  void step() {
    Serial.println("step Transition");
    if (t == 0) {
      from = mutant;
    }
    Animation::step();
    auto t = std::min(1.0, static_cast<double>(this->t) / duration);
    auto n = easing_fn(t) * (to - from) + from;
    if (quantized) {
      n = std::floor(n);
    }
    mutant = n;
  }

private:

  double& mutant;
  double from;
  double to;
  EasingFn easing_fn;
  bool quantized;
};

class Mutation : public Animation {
public:

  Mutation(double& mutant, MutateFn f, int duration, EasingFn easing_fn, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    from(mutant),
    f(f),
    easing_fn(easing_fn)
  {}

  void step() {
    Serial.println("step Mutation");
    if (t == 0) {
      from = mutant;
    }
    auto t = std::min(1.0, static_cast<double>(this->t) / (duration - 1));
    mutant = f(easing_fn(t));
    Animation::step();
  }

private:

  double& mutant;
  double from;
  MutateFn f;
  EasingFn easing_fn;
};


class Sprite : public Animation {
public:

  Sprite(SpriteFn draw_fn, int duration,
    int fade_in_duration = 0, EasingFn fade_in_easing_fn = ease_mid,
    int fade_out_duration = 0, EasingFn fade_out_easing_fn = ease_mid) :
    Animation(duration, false),
    draw_fn(draw_fn),
    fader(fade_in_duration > 0 ? 0 : 1),
    fade_in_duration(fade_in_duration),
    fade_out_duration(fade_out_duration),
    fade_in(fader, 1, fade_in_duration, fade_in_easing_fn),
    fade_out(fader, 0, fade_out_duration, fade_out_easing_fn)
  {}

  void step() {
    Serial.println("step Sprite");
    if (!fade_in.done()) {
      fade_in.step();
    }
    else if (duration >= 0 && t >= (duration - fade_out_duration)) {
      fade_out.step();
    }
    draw_fn(fader);
    Animation::step();
  }

private:

  SpriteFn draw_fn;
  double fader;
  int fade_in_duration;
  int fade_out_duration;
  Transition fade_in;
  Transition fade_out;
};

template <int W>
class Motion : public Animation {
public:

  Motion(Orientation& orientation, std::shared_ptr<Path<W>> path, int duration, bool repeat = false) :
    Animation(duration, repeat),
    orientation(orientation),
    path(path),
    to(path->get_point(0))
  {}

  void step() {
    from = to;
    to = path->get_point(static_cast<double>(t) / duration);
    if (from != to) {
      auto axis = cross(from, to).normalize();
      auto angle = angle_between(from, to);
      auto origin = orientation.get();
      orientation.clear();
      for (auto a = MAX_ANGLE; angle - a > 0.0001; a += MAX_ANGLE) {
        orientation.push(make_rotation(axis, a) * origin);
      }
      orientation.push(make_rotation(axis, angle) * origin);
    }
    Animation::step();
  }

private:

  static constexpr double MAX_ANGLE = 2 * PI / W;
  Orientation& orientation;
  std::shared_ptr<Path<W>> path;
  Vector from;
  Vector to;
};

template <int W>
class Rotation : public Animation {
public:

  Rotation(Orientation& orientation, const Vector& axis, double angle, int duration, EasingFn easing_fn, bool repeat = false) :
    Animation(duration, repeat),
    orientation(orientation),
    axis(axis),
    total_angle(angle),
    easing_fn(easing_fn),
    from(0),
    to(0)
  {
  }

  void step() {
    Serial.println("step Rotation");
    from = to;
    to = easing_fn(static_cast<double>(t) / duration) * total_angle;
    auto angle = distance(from, to, total_angle);
    if (angle > 0.00001) {
      auto origin = orientation.get();
      for (auto a = MAX_ANGLE; angle - a > 0.00001; a += MAX_ANGLE) {
        orientation.push(make_rotation(axis, a) * origin);
      }
      orientation.push(make_rotation(axis, angle) * origin);
    }
    Animation::step();
  }

private:

  static constexpr double MAX_ANGLE = 2 * PI / W;
  Orientation& orientation;
  Vector axis;
  double total_angle;
  EasingFn easing_fn;
  double from;
  double to;
};

template <int W>
void rotate_between(Orientation& from, const Orientation& to) {
  auto diff = to.get() * from.get().inverse();
  auto angle = 2 * acos(diff.r);
  if (angle == 0) {
    return;
  }
  auto axis = Vector(diff.v.i, diff.v.j, diff.v.k).normalize();
  Rotation<W>(from, axis, angle, 1, ease_out_circ).step();
}

///////////////////////////////////////////////////////////////////////////////

template<int W>
void plot_dots(const Dots& dots, Filter<W>& filters, Pixels& pixels, double age = 0) {
  for (auto& dot : dots) {
    Spherical s(dot.position);
    double y = (s.phi * H) / PI;
    if (fabs(H - y) < 0.0001) {
      continue;
    }
    double x = fmod(((s.theta + PI) * W) / (2 * PI), W);
    filters.plot(pixels, x, y, dot.color, age, blend_over_add);
  }
}
