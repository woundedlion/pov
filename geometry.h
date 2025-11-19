#pragma once

static constexpr Vector X_AXIS(1, 0, 0);
static constexpr Vector Y_AXIS(0, 1, 0);
static constexpr Vector Z_AXIS(0, 0, 1);

struct PixelCoords {
  float x;
  float y;
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


template <int W>
Vector pixel_to_vector(float x, float y) {
  return Vector(
    Spherical(
      (x * 2 * PI_F) / W,
      (y * PI_F) / (H_VIRT - 1)
    )
  );
}

template <int W>
PixelCoords vector_to_pixel(const Vector& v) {
  auto s = Spherical(v);
  assert(std::abs(v.length()) - 1 < TOLERANCE);
  PixelCoords p({ wrap((s.theta * W) / (2 * PI_F), W), (s.phi * (H_VIRT - 1)) / PI_F });
  return p;
}

constexpr float lerp(float from, float to, float t) {
  return ((to - from) * t) + from;
}

WaveFn sin_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    auto w = (sinf(freq * t * 2 * PI_F - (PI_F / 2) + PI_F - (2 * phase)) + 1) / 2;
    return lerp(from, to, w);
    };
}

WaveFn tri_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    float w = wrap(t * freq + phase, 1.0f);
    if (w < 0.5f) {
      w = 2.0f * w;
    } else {
      w = 2.0f * (1.0f - w);
    }
    return lerp(from, to, w);
    };
}

WaveFn square_wave(float from, float to, float freq, float dutyCycle, float phase) {
  return [=](float t) -> float {
    if (fmod(t * freq + phase, 1.0f) < dutyCycle) {
      return to;
    }
    return from;
    };
}

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
    return rotate(v, orientations[num_frames - 1]);
  }

  Vector orient(const Vector& v, int i) const {
    return rotate(v, orientations[i]);
  }

  Vector unorient(const Vector& v) const {
    return rotate(v, orientations[num_frames - 1].inverse());
  }

  Vector unorient(const Vector& v, int i) const {
    return rotate(v, orientations[i].inverse());
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
      [this, i](const auto& v) {
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
    } else {
      Serial.println("Orientation full, droping frame!");
      assert(false);
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
  static constexpr int MAX_FRAMES = MAX_W + 1;
  std::array<Quaternion, MAX_FRAMES> orientations;
  size_t num_frames;
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
      }
      else {
        ++j;
      }
    }
  }
}

Pixel distance_gradient(const Vector& v, const Vector& normal, CRGBPalette256 p1, CRGBPalette256 p2) {
  auto d = dot(v, normal);
  if (d > 0) {
    return p1[static_cast<int>(d * 255)];
  }
  else {
    return p2[static_cast<int>(-d * 255)];
  }
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

  Dodecahedron() {
    for (auto& v : vertices) {
      v.normalize();
    }
  }

  void dump() const {
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

Vector random_vector() {
  // Marsaglia's method
  float v1, v2, s;
  do {
    v1 = 2.0f * hs::rand_f() - 1.0f;
    v2 = 2.0f * hs::rand_f() - 1.0f;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1.0f || s == 0.0f);

  float sqrt_s = sqrtf(1.0f - s);
  return Vector(
    2.0f * v1 * sqrt_s,
    2.0f * v2 * sqrt_s,
    1.0f - 2.0f * s
  );
}

struct LissajousParams {
  float m1;
  float m2;
  float a;
  float domain;
};


Vector lissajous(float m1, float m2, float a, float t) {
  Vector v(
    sinf(m2 * t) * cosf(m1 * t - a * PI_F),
    cosf(m2 * t),
    sinf(m2 * t) * sinf(m1 * t - a * PI_F)
    );
  return v.normalize();
}