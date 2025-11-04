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
#include "FastNoiseLite.h"

typedef CRGB Pixel;

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

struct PixelCoords {
  double x;
  double y;
};

typedef std::vector<Dot> Dots;
typedef std::vector<Vector> VertexList;
typedef std::vector<std::vector<unsigned int>> AdjacencyList;
typedef std::function<Pixel(const Vector&, double)> ColorFn;
typedef std::function<Pixel(double x, double y, double t)> TrailFn;
typedef std::function<double (double)> EasingFn;
typedef std::function<Vector (double)> PlotFn;
typedef std::function<double (double)> ShiftFn;
typedef std::function<void (Canvas&, double)> SpriteFn;
typedef std::function<void (Canvas&)> TimerFn;
typedef std::function<double (double)> MutateFn;
typedef std::function<double (double)> WaveFn;


inline int XY(int x, int y) { return x * H + y; }

static const int FPS = 16;
static const Vector X_AXIS(1, 0, 0);
static const Vector Y_AXIS(0, 1, 0);
static const Vector Z_AXIS(0, 0, 1);

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


Pixel dim(const Pixel c, double s) {
  if (s < 0.5) s = s * s;
  return Pixel(
    static_cast<uint8_t>(c.r * s), 
    static_cast<uint8_t>(c.g * s),
    static_cast<uint8_t>(c.b * s));
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
Dots draw_line(const Vector& v1, const Vector& v2, ColorFn color, bool long_way = false);

template <int W>
class Path {
public:
  Path() {}

  Path& append_line(const Vector& v1, Vector& v2, bool long_way = false) {
    if (points.size() > 0) {
      points.pop_back(); // Overlap previous segment
    }
    Dots seg = draw_line<W>(v1, v2, [](auto& , auto) { return Pixel(); }, long_way);
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

Pixel distance_gradient(const Vector& v, const Vector& normal, CRGBPalette256 p1, CRGBPalette256 p2) {
  auto d = dot(v, normal);
  if (d > 0) {
    return p1[static_cast<int>(d * 255)];
  } else {
    return p2[static_cast<int>(-d * 255)];
  }
}
///////////////////////////////////////////////////////////////////////////////

Dots draw_vector(const Vector& v, ColorFn color_fn) {
  auto u = v.normalize();
  return { Dot(u, color_fn(u, 0)) };
}

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

  virtual void plot(Canvas& canvas, double x, double y, 
    const Pixel& c, double age, double alpha) = 0;

protected:
  
  void pass(Canvas& canvas, double x, double y,
    const Pixel& c, double age, double alpha)
  {
    if (next == nullptr) {
      canvas(XY(x, y)) = { blend_alpha(alpha)(canvas(XY(x, y)), c) };
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

    double v = (1 - x_m) * (1 - y_m);
    this->pass(canvas, x_i, y_i, dim(c, v), age, alpha);

    v = x_m * (1 - y_m);
    this->pass(canvas, (static_cast<int>(x_i + 1)) % W, y_i, dim(c, v), age, alpha);

    if (y_i < H - 1) {
      v = (1 - x_m) * y_m;
      this->pass(canvas, x_i, y_i + 1, dim(c, v), age, alpha);

      v = x_m * y_m;
      this->pass(canvas, static_cast<int>(x_i + 1) % W, y_i + 1, dim(c, v), age, alpha);
    }
  }
};

template <int W>
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
        ttls[num_pixels++] = { x, y, lifetime - age };
      }
    }
    if (age <= 0) {
      pass(canvas, x, y, color, age, alpha);
    }
  }

  void trail(Canvas& canvas, TrailFn trailFn, double alpha) {
    for (int i = 0; i < num_pixels; ++i) {
      auto color = trailFn(ttls[i].x, ttls[i].y, 1 - (ttls[i].ttl / lifetime));
      pass(canvas, ttls[i].x, ttls[i].y, color, lifetime - ttls[i].ttl, alpha);
    }
  }

private:

  struct DecayPixel {
    double x;
    double y;
    float ttl;
  };

  int lifetime;
  static constexpr MAX_PIXELS = 6 * 1024;
  std::array<DecayPixel, MAX_PIXELS> ttls;
  size_t num_pixels;
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
    pass(canvas, x, y, color, age, alpha);
    pass(canvas, wrap(x + 1, W), y, r, age, alpha);
    pass(canvas, wrap(x + 2, W), y, g, age, alpha);
    pass(canvas, wrap(x + 3, W), y, b, age, alpha);
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
}

template <int W>
class FilterReplicate : public Filter<W> {
public:

  FilterReplicate(size_t count) :
    count(std::max(1, std::min(W, count)))
  {}

  void plot(Canvas& canvas, double x, double y, const Pixel& color, double age, double alpha) {
    for (int i = 0; i < W; i += W / count) {
      pass(canvas, wrap(x + i, W), y, color, age, alpha);
    }
  }

private:

  size_t count;
}
///////////////////////////////////////////////////////////////////////////////

uint16_t to_short(double zero_to_one) {
  return std::max(0, std::min(65535, std::round(zero_to_one * 65535)));
}

///////////////////////////////////////////////////////////////////////////////

auto vignette(const Palette& palette) {
  CRGB vignette_color(0, 0, 0);
  return [=](double t) {
    if (t < 0.2) {
      return vignette_color.lerp16(palette.get(0), to_short(t / 0.2));
    } else if (t >= 0.8) {
      return palette.get(1).lerp16(vignette_color, to_short((t - 0.8) / 0.2));
    } else {
      return palette.get((t - 0.2) / 0.6);
    }
  }
}

class Palette {
public:
  Pixel get(double t) const  = 0;
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

class GenerativePalette {
public:

  GenerativePalette(GradientShape shape, HarmonyType harmony_type) :
    gradient_shape(shape),
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
  }

  Pixel get(double t) const {
    std::array<double> shape;
    std::array<Pixel> colors;
    const Pixel vignette_color(0, 0, 0);

    switch (gradient_shape) {
    case GradientShape::VIGNETTE:
      shape = { 0, 0.1, 0.5, 0.9, 1 };
      colors = { vignette_color, a, b, c, vignette_color };
      break;
    case GradientShape::STRAIGHT:
      shape = { 0, 0.5, 1 };
      colors = { a, b, c };
      break;
    case GradientShape::CIRCULAR:
      shape = { 0, 0.33, 0.66, 1 };
      colors = { a, b, c, a };
      break;
    case GradientShape::FALLOFF:
      shape = { 0, 0.33, 0.66, 0.9, 1 };
      colors = { a, b, c, vignette_color };
      break;
    }

    int seg = -1;
    for (int i = 0; i < shape.size() - 1; ++i) {
      if (t >= shape[i] && t < shape[i + 1]) {
        seg = i;
        break;
      }
    }
    if (seg < 0) {
      seg = shape.size() - 1;
    }

    auto start = shape[seg];
    auto end = shape[seg + 1];
    auto c1 = colors[seg];
    auto c2 = colors[seg + 1];

    return c1.lerp16(c2, std::max(0, std::min(65535, (t - start) / (end - start) * 65535)));
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
}

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

  Pixel get(double t) const {
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

static const ProceduralPalette richSunset(
  { 0.309, 0.500, 0.500 }, // A
  { 1.000, 1.000, 0.500 }, // B
  { 0.149, 0.148, 0.149 }, // C
  { 0.132, 0.222, 0.521 }  // D
);

static const ProceduralPalette underSea(
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

class Animation {
public:

  Animation(int duration, bool repeat) :
    duration(duration),
    repeat(repeat),
    canceled(false),
    post([](){})
  {
  }

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

  Animation& then(std::function<void()> callback) {
    post = callback;
    return *this;
  }

  void post_callback() {
    post();
  }

protected:

  int duration;
  bool repeat;
  int t = 0;

private:

  bool canceled;
  std::function<void()> post;
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

  void step(Canvas& canvas) {
    Serial.println("step RandomTimer");
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

class PeriodicTimer : public Animation {
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

  void step(Canvas& canvas) {
    Serial.println("step Mutation");
    if (t == 0) {
      from = mutant;
    }
    auto t = std::min(1.0, static_cast<double>(this->t) / (duration - 1));
    mutant = f(easing_fn(t));
    Animation::step(canvas);
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

  void step(Canvas& canvas) {
    Serial.println("step Sprite");
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
class Motion : public Animation {
public:

  Motion(Orientation& orientation, std::unique_ptr<Path<W>> path, int duration, bool repeat = false) :
    Animation(duration, repeat),
    orientation(orientation),
    path(std::move(path)),
    to(this->path->get_point(0))
  {}

  void step(Canvas& canvas) {
    from = to;
    to = path->get_point(static_cast<double>(t) / duration);
    if (from != to) {
      auto axis = cross(from, to).normalize();
      auto angle = angle_between(from, to);
      auto origin = orientation.get();
      orientation.collapse();
      for (auto a = MAX_ANGLE; angle - a > 0.0001; a += MAX_ANGLE) {
        orientation.push(make_rotation(axis, a) * origin);
      }
      orientation.push(make_rotation(axis, angle) * origin);
    }
    Animation::step(canvas);
  }

private:

  static constexpr double MAX_ANGLE = 2 * PI / W;
  Orientation& orientation;
  std::unique_ptr<Path<W>> path;
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

  template <int W>
  static void animate(Canvas& canvas, Orientation& orientation, const Vector& axis, double angle, EasingFn easing_fn) {
    Rotation<W> r(orientation, axis, angle, 2, easing_fn, false);
    r.step(canvas);
    r.step(canvas);
  }

  void step(Canvas& canvas) {
    orientation.collapse();
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
    Animation::step(canvas);
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

template<int W>
class RandomWalk : public Animation {
public:
  RandomWalk(Orientation& orientation, const Vector& v_start) :
    Animation(-1, false),
    orientation(orientation),
    v(v_start.normalize())
  {
    
    Vector u;
    if (std::abs(this->v.dot(X_AXIS)) > 0.9) {
      u = Y_AXIS;
    }
    else {
      u = X_AXIS;
    }
    direction = cross(v, u).normalize();

    noiseGenerator.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
    noiseGenerator.SetFrequency(NOISE_SCALE);
  }

  void step(Canvas& canvas) override {
    Animation::step(canvas);
    double pivotAngle = noiseGenerator.GetNoise(t * NOISE_SCALE) * PIVOT_STRENGTH;
    direction = rotate(direction, make_rotation(v, pivotAngle)).normalize();
    auto walk_axis = cross(v, direction).normalize();
    v = rotate(v, make_rotation(walk_axis, WALK_SPEED)).normalize();
    direction = rotate(direction, make_rotation(walk_axis, WALK_SPEED)).normalize();
    Rotation::animate<W>(canvas, orientation, walk_axis, WALK_SPEED, ease_mid);
  }

private:
  static constexpr double WALK_SPEED = 0.2;
  static constexpr double PIVOT_STRENGTH = 0.1;
  static constexpr double NOISE_SCALE = 0.05;

  FastNoiseLite noiseGenerator;
  Orientation& orientation;
  Vector v; 
  Vector direction;
};

///////////////////////////////////////////////////////////////////////////////

using AnimationVariant = std::variant<
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
      Serial.println("Timeline full, failed to add animation!")
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
      std::visit([=](auto& a) { a.step(canvas); }, events[i].animation);
      bool done = std::visit([](auto& a) { return a.done(); }, events[i].animation);
      if (done) {
        std::visit([](auto& a) { a.post_callback(); }, events[i].animation)
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

  static constexpr MAX_EVENTS = 256;
  std::array<TimelineEvent, MAX_EVENTS> events;
  size_t num_events;
};

///////////////////////////////////////////////////////////////////////////////

template<int W>
void plot_dots(const Dots& dots, Filter<W>& filters, Canvas& canvas, double age, double alpha) {
  for (auto& dot : dots) {
    Spherical s(dot.position);
    double y = std::max(0, std::min(H - 1, (s.phi * H) / PI));
    double x = std::max(0, std::min(W - 1, fmod(((s.theta + PI) * W) / (2 * PI), W)));
    filters.plot(canvas, x, y, dot.color, age, alpha);
  }
}

template <int W>
Vector pixel_to_vector(double x, double y) {
  return Vector(
    Spherical(
      (x * 2 * PI) / W,
      (y * PI) / (H - 1),
    )
  );
}

template <int W>
PixelCoords vector_to_pixel(const Vector& v) {
  auto s = Spherical(v);
  return PixelCoords({ wrap((s.theta * W) / (2 * PI), W), (s.phi * (H - 1)) / PI });
}

///////////////////////////////////////////////////////////////////////////////

template <typename T, size_t N>
class StaticCircularBuffer {
public:

  StaticCircularBuffer() : 
    head(0), tail(0), count(0) 
  {}

  StaticCircularBuffer(std::initializer_list<T> items) 
    : head(0), tail(0), count(0) {
    assert(items.size() <= N && "Initializer list is larger than buffer capacity.");

    for (const T& item : items) {
      push(item);
    }
  }

  void push(const T& item) {
    assert(!is_full() && "Buffer overflow: cannot push to a full buffer.");

    buffer[tail] = item;
    tail = (tail + 1) % N;
    count++;
  }

  void push(T&& item) {
    assert(!is_full() && "Buffer overflow: cannot push to a full buffer.");
    buffer[tail] = std::move(item);
    tail = (tail + 1) % N;
    count++;
  }

  void pop() {
    assert(!is_empty() && "Buffer underflow: cannot pop from an empty buffer.");
    head = (head + 1) % N;
    count--;
  }

  T& front() {
    assert(!is_empty() && "Buffer is empty.");
    return buffer[head];
  }

  const T& front() const {
    assert(!is_empty() && "Buffer is empty.");
    return buffer[head];
  }

  T& operator[](size_t index) {
    assert(index < count && "Index out of bounds.");
    return buffer[(head + index) % N];
  }

  const T& operator[](size_t index) const {
    assert(index < count && "Index out of bounds.");
    return buffer[(head + index) % N];
  }

  constexpr bool is_empty() const {
    return count == 0;
  }

  constexpr bool is_full() const {
    return count == N;
  }

  constexpr size_t size() const {
    return count;
  }

  constexpr size_t capacity() const {
    return N;
  }

  void clear() {
    head = 0;
    tail = 0;
    count = 0;
  }

private:
  std::array<T, N> buffer;
  size_t head;
  size_t tail;
  size_t count;
};