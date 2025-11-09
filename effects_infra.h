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
#include <variant>
#include <array>
#include "3dmath.h"
#include "FastNoiseLite.h"
#include "StaticCircularBuffer.h"

typedef CRGB Pixel;

struct PixelCoords {
  double x;
  double y;
};

struct Dot {
  Dot(const Vector& v, const Pixel& color) :
    position(v),
    color(color)
  {}

  Dot(const Dot& d)
    : position(d.position), color(d.color) {}

  Vector position;
  Pixel color;
};

using Dots = StaticCircularBuffer<Dot, 1024>;

typedef std::vector<Vector> VertexList;
typedef std::vector<std::vector<unsigned int>> AdjacencyList;
typedef std::function<Pixel(const Vector&, double)> ColorFn;
typedef std::function<Pixel(double x, double y, double t)> TrailFn;
typedef std::function<double(double)> EasingFn;
typedef std::function<Vector(double)> PlotFn;
typedef std::function<double(double)> ShiftFn;
typedef std::function<void(Canvas&, double)> SpriteFn;
typedef std::function<void(Canvas&)> TimerFn;
typedef std::function<double(double)> MutateFn;
typedef std::function<double(double)> WaveFn;

static const int FPS = 16;
static const Vector X_AXIS(1, 0, 0);
static const Vector Y_AXIS(0, 1, 0);
static const Vector Z_AXIS(0, 0, 1);


template <int W>
Vector pixel_to_vector(double x, double y) {
  return Vector(
    Spherical(
      (x * 2 * PI) / W,
      (y * PI) / (H - 1)
      )
  );
}

template <int W>
PixelCoords vector_to_pixel(const Vector& v) {
  auto s = Spherical(v);
  return PixelCoords({ wrap((s.theta * W) / (2 * PI), W), (s.phi * (H - 1)) / PI });
}

Vector random_vector() {
  return Vector(
    hs::rand_dbl(),
    hs::rand_dbl(),
    hs::rand_dbl()
  ).normalize();
}

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

    {0, PHI, 1 / PHI},   //8
    {0, -PHI, 1 / PHI},  // 9
    {0, PHI, -1 / PHI},  // 10
    {0, -PHI, -1 / PHI}, // 11

    {1 / PHI, 0, PHI},   // 12
    {1 / PHI, 0, -PHI},  // 13
    {-1 / PHI, 0, PHI},  // 14
    {-1 / PHI, 0, -PHI}, // 15

    {PHI, 1 / PHI, 0},   // 16
    {PHI, -1 / PHI, 0},  // 17
    {-PHI, 1 / PHI, 0},  // 18
    {-PHI, -1 / PHI, 0}, // 19


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

auto blend_add(const Pixel& c1, const Pixel& c2) {
  return Pixel(
    qadd8(c1.r, c2.r),
    qadd8(c1.g, c2.g),
    qadd8(c1.b, c2.b));
}


auto blend_alpha(double a) {
  return [a](const Pixel& c1, const Pixel& c2) {
    return Pixel(
      qadd8(c1.r * (1 - a), c2.r * a),
      qadd8(c1.g * (1 - a), c2.g * a),
      qadd8(c1.b * (1 - a), c2.b * a));
  };
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
  Orientation() : 
    num_frames(0)
  {
    set(Quaternion());
  }
  
  Orientation(const Quaternion& q) :
    num_frames(0)
  {
    set(q);
  }

  int length() const { return num_frames; }
 
  Vector orient(const Vector& v) const {
    return rotate(Vector(v).normalize(), orientations[num_frames - 1]);
  }
  
  Vector orient(const Vector& v, int i) const {
    return rotate(Vector(v).normalize(), orientations[i]);
  }
  
  Vector unorient(const Vector& v) const {
    return rotate(Vector(v).normalize(), orientations[num_frames - 1].inverse());
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

  const Quaternion& get() const {
    return orientations[num_frames - 1];
  }

  const Quaternion& get(int i) const {
    return orientations[i];
  }

  Orientation& set(const Quaternion& q) {
    orientations[0] = q;
    num_frames = 1;
    return *this;
  }

  Orientation& push(const Quaternion& q) {
    if (num_frames < MAX_FRAMES) {
      orientations[num_frames++] = q;
    }
    return *this;
  }

  Quaternion pop() {
    auto r = orientations.back();
    num_frames--;
    return r;
  }

  Orientation& collapse() {
    if (num_frames > 1) {
      orientations[0] = orientations[num_frames - 1];
      num_frames = 1;
    }
    return *this;
  }

private:
  static constexpr int MAX_FRAMES = MAX_W;
  std::array<Quaternion, MAX_FRAMES> orientations;
  size_t num_frames;
};

template <int W>
void draw_line(Dots& dots, const Vector& v1, const Vector& v2, ColorFn color, bool long_way = false);

template <int W>
class Path {
public:
  Path() {}

  Path& append_line(const Vector& v1, Vector& v2, bool long_way = false) {
    if (points.size() > 0) {
      points.pop_back(); // Overlap previous segment
    }
    Dots seg;
    draw_line<W>(seg, v1, v2, [](auto& , auto) { return Pixel(); }, long_way);
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

Pixel distance_gradient(const Vector& v, const Vector& normal, CRGBPalette256 p1, CRGBPalette256 p2) {
  auto d = dot(v, normal);
  if (d > 0) {
    return p1[static_cast<int>(d * 255)];
  } else {
    return p2[static_cast<int>(-d * 255)];
  }
}
///////////////////////////////////////////////////////////////////////////////

template <int W>
void draw_path(Dots& dots, const Path<W>& path, ColorFn color) {
  size_t samples = path.num_points();
  for (size_t i = 0; i < samples; ++i) {
    auto v = path.get_point(static_cast<double>(i) / samples);
    dots.push_back(Dot(v, color(v, i / (samples - 1))));
  }
}

void draw_vector(Dots& dots, const Vector& v, ColorFn color_fn) {
  Vector u(v);
  u.normalize();
  dots.emplace_back(Dot(u, color_fn(u, 0)));
}

template <int W>
void draw_line(Dots& dots, const Vector& v1, const Vector& v2, ColorFn color, bool long_way /* = false*/) {
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
}

void draw_vertices(Dots& dots, const VertexList& vertices, ColorFn color_fn) {
  for (Vector v : vertices) {
    dots.emplace_back(Dot(v.normalize(), color_fn(v, 0)));
  }
}

template <int W>
void draw_polyhedron(Dots& dots, const VertexList& vertices, const AdjacencyList& edges, ColorFn color_fn) {
  for (size_t i = 0; i < edges.size(); ++i) {
    Vector a(vertices[i]);
    for (auto j : edges[i]) {
      Vector b(vertices[j]);
      draw_line<W>(dots, a, b, color_fn);
    }
  }
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
  Vector axis = cross(v, vp).normalize();
  auto shift = make_rotation(axis, f(angle * PI / 2));
  return rotate(vi, shift);
};

template <int W>
void draw_fn(Dots& dots, const Orientation& orientation, const Vector& normal, double radius, ShiftFn shift_fn, ColorFn color_fn) {
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
      draw_line<W>(dots, from, to, color_fn);
    }
    from = to;
  }
  draw_line<W>(dots, from, start, color_fn);
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
void draw_ring(Dots& dots, const Vector& normal, double radius, ColorFn color_fn, double phase = 0) {
  Vector v(normal);
  if (radius > 1) {
    v = -v;
    phase = wrap(phase + PI, 2 * PI);
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

  virtual void plot(Canvas& canvas, double x, double y, 
    const Pixel& c, double age, double alpha) = 0;

protected:
  
  void pass(Canvas& canvas, double x, double y,
    const Pixel& c, double age, double alpha)
  {
    if (next == nullptr) {
      auto p = blend_alpha(alpha)(canvas(XY(x, y)), c);
      canvas(XY(x, y)) = p;
//      Serial.printf("plot_out (%f, %f,) (%d, %d, %d)\n", x, y, p.r, p.g, p.b);
    }
    else {
      next->plot(canvas, x, y, c, age, alpha);
    }
  }

  Filter* next;
};

template <int W>
class FilterRaw : public Filter<W> {
  void plot(Canvas& canvas, double x, double y, const Pixel& c, double age, double alpha) {
      this->pass(canvas, x, y, c, age, alpha);
    }
};

template <int W>
class FilterAntiAlias : public Filter<W> {
public:
  FilterAntiAlias() {}
  
  void plot(Canvas& canvas, double x, double y, const Pixel& c, double age, double alpha) {
    double x_i = 0;
    double x_m = modf(x, &x_i);
    double y_i = 0;
    double y_m = modf(y, &y_i);
//    Serial.printf("plot_in (%f, %f,) (%d, %d, %d) [%f]\n", x, y, c.r, c.g, c.b, alpha);

    double v = (1 - x_m) * (1 - y_m);
    this->pass(canvas, x_i, y_i, c, age, alpha * v);

    v = x_m * (1 - y_m);
    this->pass(canvas, (static_cast<int>(x_i + 1)) % W, y_i, c, age, alpha * v);

    if (y_i < H - 1) {
      v = (1 - x_m) * y_m;
      this->pass(canvas, x_i, y_i + 1, c, age, alpha * v);

      v = x_m * y_m;
      this->pass(canvas, static_cast<int>(x_i + 1) % W, y_i + 1, c, age, alpha * v);
    }
  }
};

template <int W, int MAX_PIXELS>
class FilterDecay : public Filter<W> {
public:

  FilterDecay(int lifetime) :
    lifetime(lifetime),
    num_pixels(0)
  {}

  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--ttls[i].ttl < 0.00001) {
        num_pixels--;
        if (i < num_pixels) {
          ttls[i] = std::move(ttls[num_pixels]);
          i--;
        }
      }
    }
  }

  void plot(Canvas& canvas, double x, double y, const Pixel& color, double age, double alpha) {
    if (age >= 0) {
      if (num_pixels < MAX_PIXELS) {
        ttls[num_pixels++] = { static_cast<float>(x), static_cast<float>(y), static_cast<float>(lifetime - age) };
      } else {
        Serial.println("FilterDecay full!");
      }
    }
    if (age <= 0) {
      this->pass(canvas, x, y, color, age, alpha);
    }
  }

  void trail(Canvas& canvas, TrailFn trailFn, double alpha) {
    for (int i = 0; i < num_pixels; ++i) {
      auto color = trailFn(ttls[i].x, ttls[i].y, 1 - (ttls[i].ttl / lifetime));
      this->pass(canvas, ttls[i].x, ttls[i].y, color, lifetime - ttls[i].ttl, alpha);
    }
  }

private:

  struct DecayPixel {
    float x;
    float y;
    float ttl;
  };

  int lifetime;
  std::array<DecayPixel, MAX_PIXELS> ttls;
  int num_pixels;
};

template<int W>
class FilterChromaticShift : public Filter<W> {
public:

  FilterChromaticShift()
  {
  }

  void plot(Canvas& canvas, double x, double y, const Pixel& color, double age, double alpha) {
    CRGB r(color.r, 0, 0);
    CRGB g(0, color.g, 0);
    CRGB b(0, 0, color.b);
    this->pass(canvas, x, y, color, age, alpha);
    this->pass(canvas, wrap(x + 1, W), y, r, age, alpha);
    this->pass(canvas, wrap(x + 2, W), y, g, age, alpha);
    this->pass(canvas, wrap(x + 3, W), y, b, age, alpha);
  }
};

template <int W>
class FilterOrient : public Filter<W> {
public:
  FilterOrient(Orientation& orientation) :
    orientation(orientation)
  {}

  void plot(Canvas& canvas, double x, double y, const Pixel& color, double age, double alpha) {
    auto v = pixel_to_vector<W>(x, y);
    auto r = vector_to_pixel<W>(orientation.orient(v));
    orientation.collapse();
    pass(canvas, r.x, r.y, color, age, alpha);
  }

private:

  Orientation& orientation;
};

template <int W>
class FilterReplicate : public Filter<W> {
public:

  FilterReplicate(int count) :
    count(std::clamp(count, 1, W))
  {}

  void plot(Canvas& canvas, double x, double y, const Pixel& color, double age, double alpha) {
    for (int i = 0; i < W; i += W / count) {
      this->pass(canvas, wrap(x + i, W), y, color, age, alpha);
    }
  }

private:

  int count;
};
///////////////////////////////////////////////////////////////////////////////

uint16_t to_short(double zero_to_one) {
  return std::clamp(std::round(zero_to_one * 65535), 0.0, 65535.0);
}

///////////////////////////////////////////////////////////////////////////////

class Palette {
public:
  virtual Pixel get(double t) const = 0;
};

enum class HarmonyType {
  TRIADIC,
  SPLIT_COMPLEMENTARY,
  COMPLEMENTARY,
  ANALOGOUS
};

enum class GradientShape {
  STRAIGHT,
  CIRCULAR,
  VIGNETTE,
  FALLOFF
};

class GenerativePalette : public Palette {
public:

  GenerativePalette(GradientShape gradient_shape, HarmonyType harmony_type) :
    gradient_shape(gradient_shape),
    harmony_type(harmony_type),
    seed_hue(static_cast<uint8_t>(hs::rand_int(0, 256)))
  {
    uint8_t h1 = seed_hue;
    uint8_t h2;
    uint8_t h3;

    seed_hue = static_cast<uint8_t>(
      (static_cast<uint32_t>(seed_hue) + static_cast<uint32_t>(G * 255.0)) % 256);
    calc_hues(h1, h2, h3, harmony_type);

    const uint8_t s1 = hs::rand_int(255 * 0.4, 255 * 0.8);
    const uint8_t s2 = hs::rand_int(255 * 0.4, 255 * 0.8);
    const uint8_t s3 = hs::rand_int(255 * 0.4, 255 * 0.8);

    const uint8_t v1 = hs::rand_int(255 * 0.1, 255 * 0.1);
    const uint8_t v2 = hs::rand_int(255 * 0.2, 255 * 0.5);
    const uint8_t v3 = hs::rand_int(255 * 0.6, 255 * 0.8);

    a = CHSV(h1, s1, v1);
    b = CHSV(h2, s2, v2);
    c = CHSV(h3, s3, v3);

    const Pixel vignette_color(0, 0, 0);
    switch (gradient_shape) {
    case GradientShape::VIGNETTE:
      shape = { 0, 0.1, 0.5, 0.9, 1 };
      colors = { vignette_color, a, b, c, vignette_color };
      size = 5;
      break;
    case GradientShape::STRAIGHT:
      shape = { 0, 0.5, 1 };
      colors = { a, b, c };
      size = 3;
      break;
    case GradientShape::CIRCULAR:
      shape = { 0, 0.33, 0.66, 1 };
      colors = { a, b, c, a };
      size = 4;
      break;
    case GradientShape::FALLOFF:
      shape = { 0, 0.33, 0.66, 0.9, 1 };
      colors = { a, b, c, vignette_color };
      size = 4;
      break;
    }
  }

  Pixel get(double t) const override {
    Serial.println(t);
    int seg = -1;
    for (int i = 0; i < size - 1; ++i) {
      if (t >= shape[i] && t < shape[i + 1]) {
        seg = i;
        break;
      }
    }
    if (seg < 0) {
      seg = size - 2;
    }

    auto start = shape[seg];
    auto end = shape[seg + 1];
    Pixel c1 = colors[seg];
    Pixel c2 = colors[seg + 1];

    auto r = c1.lerp16(c2, std::clamp((t - start) / (end - start) * 65535, 0.0, 65535.0));
    return r;
  }

  private:

    uint8_t wrap_hue(int hue) {
      return (hue % 256 + 256) % 256;
    }

    void calc_hues(uint8_t h1, uint8_t& h2, uint8_t& h3, HarmonyType harmony_type) {
      const int h1_int = h1;

      switch (harmony_type) {
        case HarmonyType::TRIADIC: {
          h2 = wrap_hue(h1_int + 256 / 3);
          h3 = wrap_hue(h1_int + (256 * 2) / 3);
          break;
        }

        case HarmonyType::SPLIT_COMPLEMENTARY: {
          const int complement = wrap_hue(h1_int + 256 / 2);
          const int offset = 256 / 12;
          h2 = wrap_hue(complement - offset);
          h3 = wrap_hue(complement + offset);
          break;
        }

        case HarmonyType::COMPLEMENTARY: {
          h2 = wrap_hue(h1_int + 256 / 2);
          const int offset = hs::rand_int(-256 / 36, 256 / 36);
          h3 = wrap_hue(h1_int + offset);
          break;
        }

        case HarmonyType::ANALOGOUS:
        default: {
          const int dir = (hs::rand_int(0, 1) == 0) ? 1 : -1;
          const int offset1 = dir * hs::rand_int(256 / 12, 256 / 6);
          h2 = wrap_hue(h1_int + offset1);
          const int offset2 = dir * hs::rand_int(256 / 12, 256 / 6);
          h3 = wrap_hue(h2 + offset2);
          break;
        }
      }
    }

    GradientShape gradient_shape;
    HarmonyType harmony_type;
    uint8_t seed_hue;
    Pixel a, b, c;
    std::array<double, 5> shape;
    std::array<Pixel, 5> colors;
    int size;
};

class ProceduralPalette : public Palette {
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

  Pixel get(double t) const override {
      return Pixel(
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

using PaletteVariant =
std::variant<
  GenerativePalette,
  ProceduralPalette
>;

auto vignette(const Palette& palette) {
  CRGB vignette_color(0, 0, 0);
  return [&](double t) {
    if (t < 0.2) {
      return vignette_color.lerp16(palette.get(0), to_short(t / 0.2));
    }
    else if (t >= 0.8) {
      return palette.get(1).lerp16(vignette_color, to_short((t - 0.8) / 0.2));
    }
    else {
      return palette.get((t - 0.2) / 0.6);
    }
  };
}


static const ProceduralPalette richSunset(
  { 0.309, 0.500, 0.500 }, // A
  { 1.000, 1.000, 0.500 }, // B
  { 0.149, 0.148, 0.149 }, // C
  { 0.132, 0.222, 0.521 }  // D
);

static const ProceduralPalette undersea(
  { 0.000, 0.000, 0.000 }, // A
  { 0.500, 0.276, 0.423 }, // B
  { 0.296, 0.296, 0.296 }, // C
  { 0.374, 0.941, 0.000 }  // D
);

static const ProceduralPalette lateSunset(
  { 0.337, 0.500, 0.096 }, // A
  { 0.500, 1.000, 0.176 }, // B
  { 0.261, 0.261, 0.261 }, // C
  { 0.153, 0.483, 0.773 }  // D
);

static const ProceduralPalette mangoPeel(
  { 0.500, 0.500, 0.500 }, // A
  { 0.500, 0.080, 0.500 }, // B
  { 0.431, 0.431, 0.431 }, // C
  { 0.566, 0.896, 0.236 }  // D
);

static const ProceduralPalette lemonLime(
  { 0.455, 0.455, 0.455 }, // A
  { 0.571, 0.151, 0.571 }, // B
  { 0.320, 0.320, 0.320 }, // C
  { 0.087, 0.979, 0.319 }  // D
);

static const ProceduralPalette algae(
  { 0.337, 0.500, 0.096 }, // A
  { 0.500, 1.000, 0.176 }, // B
  { 0.134, 0.134, 0.134 }, // C
  { 0.328, 0.658, 0.948 }  // D
);



Pixel dotted_brush(const Pixel& color, double freq, double duty_cycle, double phase, double t) {
  return dim(color, square_wave(0, 1, freq, duty_cycle, phase)(t));
}

///////////////////////////////////////////////////////////////////////////////

template <typename Derived>
class Animation {
public:

  void cancel() { canceled = true; }
  virtual bool done() { return canceled || (duration >= 0 && t >= duration); }

  virtual void step(Canvas& canvas) {
    t++;
    if (done()) {
      if (repeat) {
        t = 0;
      }
    }
  }

  Derived& then(std::function<void()> callback) {
    post = callback;
    return static_cast<Derived&>(*this);
  }

  void post_callback() {
    post();
  }

protected:

  Animation(int duration, bool repeat) :
    duration(duration),
    repeat(repeat),
    canceled(false),
    post([]() {})
  {
  }

  int duration;
  bool repeat;
  int t = 0;

private:

  bool canceled;
  std::function<void()> post;
};

class NullAnimation : public Animation<NullAnimation> {
public:
  NullAnimation() :
    Animation(0, false)
  {}
  void step(Canvas&) {}
  bool done() { return true; }
};

class RandomTimer : public Animation<RandomTimer> {
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
  
  void step(Canvas& canvas) {
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
    Animation::step(canvas);
  }
  
  private:

    int min;
    int max;
    TimerFn f;
    int next;
};

class PeriodicTimer : public Animation<PeriodicTimer> {
public:
  PeriodicTimer(int period, TimerFn f, bool repeat = false) :
    Animation(-1, repeat),
    period(period),
    f(f)
  {
    reset();
  }

  void reset() {
    next = t + period;
  }

  void step(Canvas& canvas) {
    if (t >= next) {
      f(canvas);
      if (repeat) {
        reset();
      }
      else {
        cancel();
      }
    }
    Animation::step(canvas);
  }

private:

  int period;
  TimerFn f;
  int next;
};

class Transition : public Animation<Transition> {
public:
  Transition(double& mutant, double to, int duration, EasingFn easing_fn, bool quantized = false, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    from(mutant),
    to(to),
    easing_fn(easing_fn),
    quantized(quantized)
  {
  }

  void step(Canvas& canvas) {
    if (t == 0) {
      from = mutant;
    }
    Animation::step(canvas);
    auto t = std::min(1.0, static_cast<double>(this->t) / duration);
    auto n = easing_fn(t) * (to - from) + from;
    if (quantized) {
      n = std::floor(n);
    }
    mutant.get() = n;
  }

  void rebind_mutant(double& new_mutant) {
    mutant = new_mutant;
  }

private:

  std::reference_wrapper<double> mutant;
  double from;
  double to;
  EasingFn easing_fn;
  bool quantized;
};

class Mutation : public Animation<Mutation> {
public:

  Mutation(double& mutant, MutateFn f, int duration, EasingFn easing_fn, bool repeat = false) :
    Animation(duration, repeat),
    mutant(mutant),
    from(mutant),
    f(f),
    easing_fn(easing_fn)
  {}

  void step(Canvas& canvas) {
    if (t == 0) {
      from = mutant;
    }
    auto t = std::min(1.0, static_cast<double>(this->t) / (duration - 1));
    mutant.get() = f(easing_fn(t));
    Animation::step(canvas);
  }


  void rebind_mutant(double& new_mutant) {
    mutant = new_mutant;
  }

private:

  std::reference_wrapper<double> mutant;
  double from;
  MutateFn f;
  EasingFn easing_fn;
};


class Sprite : public Animation<Sprite> {
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

  Sprite(Sprite&& other) noexcept
    : Animation(std::move(other)),
    draw_fn(std::move(other.draw_fn)),
    fader(other.fader),
    fade_in_duration(other.fade_in_duration),
    fade_out_duration(other.fade_out_duration),
    fade_in(std::move(other.fade_in)),
    fade_out(std::move(other.fade_out))
  {
    fade_in.rebind_mutant(this->fader);
    fade_out.rebind_mutant(this->fader);
  }

  Sprite& operator=(Sprite&& other) noexcept {
    if (this == &other) return *this;

    Animation::operator=(std::move(other));
    draw_fn = std::move(other.draw_fn);
    fader = other.fader;
    fade_in_duration = other.fade_in_duration;
    fade_out_duration = other.fade_out_duration;
    fade_in = std::move(other.fade_in);
    fade_in.rebind_mutant(this->fader);
    fade_out = std::move(other.fade_out);
    fade_out.rebind_mutant(this->fader);

    return *this;
  }

  void step(Canvas& canvas) {
    if (!fade_in.done()) {
      fade_in.step(canvas);
    }
    else if (duration >= 0 && t >= (duration - fade_out_duration)) {
      fade_out.step(canvas);
    }
    draw_fn(canvas, fader);
    Animation::step(canvas);
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
class Motion : public Animation<Motion<W>> {
public:

  Motion(Orientation& orientation, const Path<W>& path, int duration, bool repeat = false) :
    Animation<Motion<W>>(duration, repeat),
    orientation(orientation),
    path(path),
    to(this->path.get_point(0))
  {}

  void step(Canvas& canvas) {
    from = to;
    to = path.get().get_point(static_cast<double>(this->t) / this->duration);
    if (from != to) {
      Vector axis = cross(from, to).normalize();
      auto angle = angle_between(from, to);
      auto& origin = orientation.get().get();
      orientation.get().collapse();
      for (auto a = MAX_ANGLE; angle - a > 0.0001; a += MAX_ANGLE) {
        orientation.get().push(make_rotation(axis, a) * origin);
      }
      orientation.get().push(make_rotation(axis, angle) * origin);
    }
    Animation<Motion<W>>::step(canvas);
  }

private:

  static constexpr double MAX_ANGLE = 2 * PI / W;
  std::reference_wrapper<Orientation> orientation;
  std::reference_wrapper<const Path<W>> path;
  Vector from;
  Vector to;
};

template <int W>
class Rotation : public Animation<Rotation<W>> {
public:

  Rotation(Orientation& orientation, const Vector& axis, double angle, int duration, EasingFn easing_fn, bool repeat = false) :
    Animation<Rotation<W>>(duration, repeat),
    orientation(orientation),
    axis(axis),
    total_angle(angle),
    easing_fn(easing_fn),
    from(0),
    to(0)
  {
  }

  static void animate(Canvas& canvas, Orientation& orientation, const Vector& axis, double angle, EasingFn easing_fn) {
    Rotation<W> r(orientation, axis, angle, 1, easing_fn, false);
    r.step(canvas);
  }

  void step(Canvas& canvas) {
    Animation<Rotation<W>>::step(canvas);
    from = to;
    to = easing_fn(static_cast<double>(this->t) / this->duration) * total_angle;
    auto angle = distance(from, to, total_angle);
    if (angle > 0.00001) {
      auto origin = orientation.get().pop();
      for (auto a = MAX_ANGLE; angle - a > 0.00001; a += MAX_ANGLE) {
        orientation.get().push(make_rotation(axis, a) * origin);
      }
    }
  }

private:

  static constexpr double MAX_ANGLE = 2 * PI / W;
  std::reference_wrapper<Orientation> orientation;
  Vector axis;
  double total_angle;
  EasingFn easing_fn;
  double from;
  double to;
};

template<int W>
class RandomWalk : public Animation<RandomWalk<W>> {
public:
  RandomWalk(Orientation& orientation, const Vector& v_start) :
    Animation<RandomWalk<W>>(-1, false),
    orientation(orientation),
    v(Vector(v_start).normalize())
  {
    
    Vector u;
    if (std::abs(dot(v, X_AXIS)) > 0.9) {
      u = Y_AXIS;
    } else {
      u = X_AXIS;
    }
    direction = cross(v, u).normalize();
    noiseGenerator.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noiseGenerator.SetFrequency(NOISE_SCALE);
    noiseGenerator.SetSeed(hs::rand_int(0, 65535));
  }

  void step(Canvas& canvas) override {
    Animation<RandomWalk<W>>::step(canvas);
    double pivotAngle = noiseGenerator.GetNoise(this->t * NOISE_SCALE, 0.0) * PIVOT_STRENGTH;
    direction = rotate(direction, make_rotation(v, pivotAngle)).normalize();
    Vector walk_axis = cross(v, direction).normalize();
    v = rotate(v, make_rotation(walk_axis, WALK_SPEED)).normalize();
    direction = rotate(direction, make_rotation(walk_axis, WALK_SPEED)).normalize();
    Rotation<W>::animate(canvas, orientation, walk_axis, WALK_SPEED, ease_mid);
  }

private:

  static constexpr double WALK_SPEED = 0.2;
  static constexpr double PIVOT_STRENGTH = 1;
  static constexpr double NOISE_SCALE = 0.5;

  FastNoiseLite noiseGenerator;
  std::reference_wrapper<Orientation> orientation;
  Vector v; 
  Vector direction;
};

///////////////////////////////////////////////////////////////////////////////

using AnimationVariant = std::variant<
  NullAnimation,
  Sprite,
  Transition,
  Mutation,
  RandomTimer,
  PeriodicTimer,
  Rotation<MAX_W>,
  Motion<MAX_W>,
  RandomWalk<MAX_W>
>;

struct TimelineEvent {
  int start;
  AnimationVariant animation;
};

class Timeline {
public:

  Timeline() :
    num_events(0)
  {}

  template <typename A>
  Timeline& add(double in_frames, A animation) {
    if (num_events >= MAX_EVENTS) {
      Serial.println("Timeline full, failed to add animation!");
      return *this;
    }
    TimelineEvent& e = events[num_events++];
    e.start = t + in_frames;
    e.animation = std::move(animation);
    return *this;
  }

  void step(Canvas& canvas) {
    t++;
    for (int i = 0; i < num_events; ++i) {
      if (t < events[i].start) {
        continue;
      }
      std::visit([&](auto& a) { a.step(canvas); }, events[i].animation);
      bool done = std::visit([](auto& a) { return a.done(); }, events[i].animation);
      if (done) {
        std::visit([](auto& a) { a.post_callback(); }, events[i].animation);
        num_events--;
        if (i < num_events) {
          events[i] = std::move(events[num_events]);
          i--;
        }
      }
    }
  }

  int t = 0;

private:

  static constexpr int MAX_EVENTS = 32;
  std::array<TimelineEvent, MAX_EVENTS> events;
  int num_events;
};

///////////////////////////////////////////////////////////////////////////////

template<int W>
void plot_dots(const Dots& dots, Filter<W>& filters, Canvas& canvas, double age, double alpha) {
  for (auto& dot : dots) {
    auto p = vector_to_pixel<W>(dot.position);
    filters.plot(canvas, p.x, p.y, dot.color, age, alpha);
  }
}


///////////////////////////////////////////////////////////////////////////////

