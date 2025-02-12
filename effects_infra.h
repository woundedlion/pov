#pragma once

#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <cassert>
#include <array>
#include <FastLED.h>
#include "3dmath.h"

typedef std::vector<Vector> VertexList;
typedef std::vector<std::vector<unsigned int>> AdjacencyList;
typedef std::function<CHSV(const Vector&)> ColorFn;
typedef std::function<double (double)> EasingFn;
typedef std::function<Vector(double)> PlotFn;
typedef std::function<CHSV(const CHSV&, const CHSV&)> BlendFn;


// inline int XY(int x, int y) { return x * H + y; }

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

Dots draw_line(const Vector& v1, const Vector& v2, ColorFn color, bool long_way = false) {
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

  // TODO: Optimize angle step
  double step = 2.0 / 30;
  for (double t = 0; t < a; t += step) {
    Vector vi(
      u.i * cos(t) + v.i * sin(t),
      u.j * cos(t) + v.j * sin(t),
      u.k * cos(t) + v.k * sin(t)
    );
    dots.emplace_back(Dot(vi, color(vi)));
  }
  return dots;
}

Dots draw_vertices(const VertexList& vertices, ColorFn color) {
  Dots dots;
  for (Vector v : vertices) {
    dots.emplace_back(Dot(v.normalize(), color(v)));
  }
  return dots;
}

Dots draw_polyhedron(const VertexList& vertices, const AdjacencyList& edges, ColorFn color) {
  Dots dots;
  
  for (size_t i = 0; i < edges.size(); ++i) {
    Vector a(vertices[i]);
    for (auto j : edges[i]) {
      Vector b(vertices[j]);
      auto seg = draw_line(a, b, color);
      dots.insert(dots.end(), seg.begin(), seg.end());
    }
  }
  
  return dots;
}

Dots draw_ring(const Vector& normal, double radius, ColorFn color) {
  Dots dots;
  Vector u(0, 0, 0);
  Vector v(normal);
  Vector x_axis(1, 0, 0);
  Vector z_axis(0, 0, 1);
  if (radius > 1) {
    v = -v;
    radius = 2 - radius;
  }
  if (v.i == 0 && v.j == 0) {
    u = cross(v, x_axis);
  }
  else {
    u = cross(v, z_axis);
  }
  u.normalize();
  Vector w(cross(v, u));
  double d = sqrt(pow(1 - radius, 2));

  // TODO: optimize angle step
  double step = 2 * PI / 96;
  for (double t = 0; t < 2 * PI; t += step) {
    Vector vi(
      d * v.i + radius * u.i * cos(t) + radius * w.i * sin(t),
      d * v.j + radius * u.j * cos(t) + radius * w.j * sin(t),
      d * v.k + radius * u.k * cos(t) + radius * w.k * sin(t)
    );
    dots.emplace_back(Dot(vi, color(vi)));
  }

  return dots;
};

class Orientation {
public:
  Orientation() : orientations{ Quaternion() } {}

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

class Path {
public:
  Path() {}

  Path& append_line(const Vector& v1, Vector& v2, bool long_way = false) {
    if (points.size() > 0) {
      points.pop_back(); // Overlap previous segment
    }
    Dots seg = draw_line(v1, v2, [](auto& v) { return CHSV(); }, long_way);
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

class Motion {
public:
  Motion(const Path& path, double duration) :
    path(path),
    duration(duration),
    to(path.get_point(0)),
    from(to)
  {}

  bool done() const { return t >= duration; }

  void move(Orientation& orientation) {
    to = path.get_point(t / duration);
    if (from != to) {
      Vector axis = cross(from, to).normalize();
      double angle = angle_between(from, to);
      Quaternion origin = orientation.get();
      orientation.clear();
      for (double a = MAX_ANGLE; a < angle; a += MAX_ANGLE) {
        orientation.push(make_rotation(axis, a) * origin);
      }
      orientation.push(make_rotation(axis, angle) * origin);
    }
    from = to;
    ++t;
  }

private:
  // TODO: Optimize MAX_ANGLE
  const double MAX_ANGLE = 2 * PI / 96;
  const Path& path;
  double duration;
  Vector to;
  Vector from;
  double t = 0;
};

Dots draw_path(const Path& path, ColorFn color) {
  Dots dots;
  size_t samples = path.num_points();
  for (size_t i = 0; i < samples; ++i) {
    auto v = path.get_point(static_cast<double>(i) / samples);
    dots.push_back(Dot(v, color(v)));
  }
  return dots;
}

class Rotation {
public:
  Rotation(const Vector& axis, double angle, double duration) :
    axis(axis),
    total_angle(angle),
    duration(duration)
  {}

  bool done() const { return t >= duration; }

  void rotate(Orientation& orientation, EasingFn easing) {
    from = to;
    to = easing(t / duration) * total_angle;
    double angle = fabs(to - from);
    if (angle > std::numeric_limits<double>::epsilon()) {
      Quaternion origin = orientation.get();
      orientation.clear();
      for (double a = MAX_ANGLE; a < angle; a += MAX_ANGLE) {
        orientation.push(make_rotation(axis, a) * origin);
      }
      orientation.push(make_rotation(axis, angle) * origin);
    }
    ++t;
  }

  Rotation& set_axis(const Vector& v) {
    axis = v;
    return *this;
  }

private:
  // TODO: Optimize MAX_ANGLE
  const double MAX_ANGLE = 2 * PI / 96;
  Vector axis;
  double total_angle;
  double duration;
  double from = 0;
  double to = 0;
  double t = 0;
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

CHSV distance_gradient(const Vector& v, const Vector& normal, CRGBPalette256 p1, CRGBPalette256 p2) {
  auto d = dot(v, normal);
  if (d > 0) {
    return rgb2hsv_approximate(p1[static_cast<int>(d * 255)]);
  } else {
    return rgb2hsv_approximate(p2[static_cast<int>(-d * 255)]);
  }
}


///////////////////////////////////////////////////////////////////////////////
/*
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
*/

template <uint8_t W, uint8_t H>
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

template<uint8_t W, uint8_t H>
Pixels plot_dots(Filter<W, H>& filter, const Dots& dots, double age = 0) {
  Pixels pixels;
  for (auto& dot : dots) {
    Spherical s(dot.position);
    double y = (s.phi * H) / PI;
    if (fabs(H - y) < 0.0001) {
      continue;
    }
    double x = fmod(((s.theta + PI) * W) / (2 * PI), W);
    filter.plot(pixels, x, y, dot.color, age, blend_over_add);
  }
  return pixels;
}

template <uint8_t W, uint8_t H>
class FilterAntiAlias : public Filter<W, H> {
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

template <uint8_t W, uint8_t H>
class FilterDecayMask : public Filter<W, H> {
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


template <uint8_t W, uint8_t H>
class FilterDecayTrails : public Filter<W, H> {
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
