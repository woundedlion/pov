#pragma once
template <int W>

void draw_vector(Dots& dots, const Vector& v, ColorFn color_fn) {
  Vector u(v);
  u.normalize();
  dots.emplace_back(Dot(u, color_fn(u, 0)));
}

template <int W>
void draw_line(Dots& dots, const Vector& v1, const Vector& v2, ColorFn color, bool long_way) {
  Vector u(v1);
  Vector v(v2);
  u.normalize();
  v.normalize();
  float a = angle_between(u, v);
  Vector w = cross(v, u);
  if (long_way) {
    a = 2 * PI_F - a;
    w = (-w).normalize();
  }
  v = cross(u, w).normalize();

  float step = 2 * PI_F / W;
  for (float t = 0; t < a; t += step) {
    Vector vi(
      u.i * cosf(t) + v.i * sinf(t),
      u.j * cosf(t) + v.j * sinf(t),
      u.k * cosf(t) + v.k * sinf(t)
    );
    dots.emplace_back(Dot(vi, color(vi, t)));
  }
}

template <int W>
class Path {
public:
  Path() {}

  Path& append_line(const Vector& v1, const Vector& v2, bool long_way = false) {
    if (points.size() > 0) {
      points.pop_back(); // Overlap previous segment
    }
    Dots seg;
    draw_line<W>(seg, v1, v2, [](auto&, auto) { return Pixel(); }, long_way);
    std::transform(seg.begin(), seg.end(), std::back_inserter(points),
      [](auto& d) { return d.position; });
    return *this;
  }

  Path& append_segment(PlotFn plot, float domain, float samples, EasingFn easing) {
    if (points.size() > 0) {
      points.pop_back(); // Overlap previous segment
    }
    for (float t = 0; t < samples; t++) {
      points.push_back(plot(easing(t / samples) * domain));
    }
    return *this;
  }

  Vector get_point(float t) const {
    return points[static_cast<int>(t * (points.size() - 1))];
  }

  size_t num_points() const { return points.size(); }

  void collapse() {
    if (points.size() > 1) {
      points = { points.back() };
    }
  }

private:
  StaticCircularBuffer<Vector, 1024> points;
};

template <int W>
void draw_path(Dots& dots, const Path<W>& path, ColorFn color) {
  size_t samples = path.num_points();
  for (size_t i = 0; i < samples; ++i) {
    auto v = path.get_point(static_cast<float>(i) / samples);
    dots.push_back(Dot(v, color(v, i / (samples - 1))));
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

Vector calc_ring_point(float a, float radius, const Vector& u, const Vector& v, const Vector& w) {
  auto d = sqrtf((1 - radius) * (1 - radius));
  return Vector(
    d * v.i + radius * u.i * cosf(a) + radius * w.i * sinf(a),
    d * v.j + radius * u.j * cosf(a) + radius * w.j * sinf(a),
    d * v.k + radius * u.k * cosf(a) + radius * w.k * sinf(a)
  ).normalize();
}

Vector fn_point(ShiftFn f, const Vector& normal, float radius, float angle) {
  Vector v(normal);
  if (radius > 1) {
    v = -v;
  }
  Vector u;
  if (std::abs(dot(v, X_AXIS)) > (1 - TOLERANCE)) {
    u = cross(v, Y_AXIS);
  }
  else {
    u = cross(v, X_AXIS);
  }
  u.normalize();
  Vector w(cross(v, u));
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }

  auto vi = calc_ring_point(angle, radius, u, v, w);
  auto vp = calc_ring_point(angle, 1, u, v, w);
  Vector axis = cross(v, vp).normalize();
  auto shift = make_rotation(axis, f(angle * PI_F / 2));
  return rotate(vi, shift);
};

template <int W>
void draw_fn(Dots& dots, const Orientation& orientation, const Vector& normal, float radius, ShiftFn shift_fn, ColorFn color_fn) {
  Vector v(orientation.orient(normal));
  if (radius > 1) {
    v = -v;
  }
  Vector u;
  if (std::abs(dot(v, X_AXIS)) > (1 - TOLERANCE)) {
    u = cross(v, Y_AXIS);
  }
  else {
    u = cross(v, X_AXIS);
  }
  u.normalize();  Vector w(cross(v, u));
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }

  bool first = true;
  Vector start, from, to;
  float step = 1.0f / W;
  for (float t = 0; t < 1; t += step) {
    auto vi = calc_ring_point(t * 2 * PI_F, radius, u, v, w);
    auto vp = calc_ring_point(t * 2 * PI_F, 1, u, v, w);
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


Vector ring_point(const Vector& normal, float radius, float angle, float phase = 0) {
  Vector v(normal);
  if (radius > 1) {
    v = -v;
  }
  Vector u;
  if (std::abs(dot(v, X_AXIS)) > (1 - TOLERANCE)) {
    u = cross(v, Y_AXIS);
  }
  else {
    u = cross(v, X_AXIS);
  }
  u.normalize();
  Vector w(cross(v, u));
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }
  return calc_ring_point(angle + phase, radius, u, v, w);
};

template<int W>
void draw_ring(Dots& dots, const Vector& normal, float radius, ColorFn color_fn, float phase = 0) {
  Vector v(normal);
  if (radius > 1) {
    v = -v;
    phase = wrap(phase + PI_F, 2 * PI_F);
  }

  Vector u;
  if (std::abs(dot(v, X_AXIS)) > (1 - TOLERANCE)) {
    u = cross(v, Y_AXIS);
  } else {
    u = cross(v, X_AXIS);
  }
  u.normalize();

  Vector w(cross(v, u)); // already a normal vector
  if (radius > 1) {
    w = -w;
    radius = 2 - radius;
  }

  float step = 2 * PI_F / W;
  for (float a = 0; a < 2 * PI_F; a += step) {
    auto vi = calc_ring_point(fmod((a + phase), (2 * PI_F)), radius, u, v, w);
    dots.emplace_back(Dot(vi, color_fn(vi, a / (2 * PI_F))));
  }
};

Pixel dotted_brush(const Pixel& color, float freq, float duty_cycle, float phase, float t) {
  return dim(color, square_wave(0, 1, freq, duty_cycle, phase)(t));
}

template<int W>
void plot_dots(const Dots& dots, Filter<W>& filters, Canvas& canvas, float age, float alpha) {
  for (auto& dot : dots) {
    auto p = vector_to_pixel<W>(dot.position);
    filters.plot(canvas, p.x, p.y, gamma_correct(dot.color), age, alpha);
  }
}

typedef std::function<Vector(Vector)> OrientFn;
typedef std::function<void(OrientFn, float)> TweenFn;

void tween(const Orientation& orientation, TweenFn draw_fn) {
  size_t s = orientation.length();
  for (size_t i = 0; i < s; ++i) {
    draw_fn([orientation, i](const auto& v) { return orientation.orient(v, i); },
      static_cast<float>((s - 1 - i)) / s);
  }
}