#pragma once

/**
 * @brief Unit vector along the Cartesian X-axis.
 */
static constexpr Vector X_AXIS(1, 0, 0);
/**
 * @brief Unit vector along the Cartesian Y-axis.
 */
static constexpr Vector Y_AXIS(0, 1, 0);
/**
 * @brief Unit vector along the Cartesian Z-axis.
 */
static constexpr Vector Z_AXIS(0, 0, 1);

/**
 * @brief Structure representing 2D floating-point pixel coordinates.
 */
struct PixelCoords {
  float x; /**< Horizontal coordinate. */
  float y; /**< Vertical coordinate. */
};

/**
 * @brief Structure representing a single rendered point (dot) in the scene.
 * @details Stores the 3D position and the color to be plotted.
 */
struct Dot {
  /**
   * @brief Constructs a Dot.
   * @param v The 3D position vector.
   * @param color The pixel color.
   */
  Dot(const Vector& v, const Pixel& color) :
    position(v),
    color(color)
  {}

  /**
   * @brief Copy constructor.
   * @param d The Dot to copy.
   */
  Dot(const Dot& d)
    : position(d.position), color(d.color) {}

  Vector position; /**< The 3D position (unit vector). */
  Pixel color; /**< The color of the dot. */
};

/**
 * @brief Type alias for a circular buffer used to store active dots/fragments.
 * @details Capacity is set to 1024.
 */
using Dots = StaticCircularBuffer<Dot, 1024>;


/**
 * @brief Converts 2D pixel coordinates to a 3D unit vector on the sphere.
 * @tparam W The width (number of columns) of the virtual LED display.
 * @param x The horizontal coordinate (0 to W).
 * @param y The vertical coordinate (0 to H_VIRT - 1).
 * @return The corresponding unit vector.
 */
template <int W>
Vector pixel_to_vector(float x, float y) {
  return Vector(
    Spherical(
      (x * 2 * PI_F) / W,
      (y * PI_F) / (H_VIRT - 1)
    )
  );
}

/**
 * @brief Converts a 3D unit vector back to 2D pixel coordinates.
 * @tparam W The width (number of columns) of the virtual LED display.
 * @param v The input unit vector.
 * @return The 2D PixelCoords.
 */
template <int W>
PixelCoords vector_to_pixel(const Vector& v) {
  auto s = Spherical(v);
  PixelCoords p({ wrap((s.theta * W) / (2 * PI_F), W), (s.phi * (H_VIRT - 1)) / PI_F });
  return p;
}

/**
 * @brief Linearly interpolates between two values.
 * @param from The starting value.
 * @param to The ending value.
 * @param t The interpolation factor (0.0 to 1.0).
 * @return The interpolated value.
 */
constexpr float lerp(float from, float to, float t) {
  return ((to - from) * t) + from;
}

/**
 * @brief Generates a sine wave function (WaveFn) with offset, amplitude, frequency, and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset.
 * @return A WaveFn functor that takes time (t) and returns a float.
 */
WaveFn sin_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    auto w = (sinf(freq * t * 2 * PI_F - (PI_F / 2) + PI_F - (2 * phase)) + 1) / 2;
    return lerp(from, to, w);
    };
}

/**
 * @brief Generates a triangle wave function (WaveFn) with offset, amplitude, frequency, and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset.
 * @return A WaveFn functor that takes time (t) and returns a float.
 */
WaveFn tri_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    float w = wrap(t * freq + phase, 1.0f);
    if (w < 0.5f) {
      w = 2.0f * w;
    }
    else {
      w = 2.0f * (1.0f - w);
    }
    return lerp(from, to, w);
    };
}

/**
 * @brief Generates a square wave function (WaveFn) with offset, amplitude, frequency, duty cycle, and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param dutyCycle The percentage of time the wave is "on" (high).
 * @param phase The starting phase offset.
 * @return A WaveFn functor that takes time (t) and returns a float.
 */
WaveFn square_wave(float from, float to, float freq, float dutyCycle, float phase) {
  return [=](float t) -> float {
    if (fmod(t * freq + phase, 1.0f) < dutyCycle) {
      return to;
    }
    return from;
    };
}

/**
 * @brief Class managing the current rotation state of an object, maintaining history for interpolation.
 * @details Stores a list of Quaternions (`orientations`) generated during the current frame step.
 */
class Orientation {
public:
  /**
   * @brief Default constructor (identity rotation).
   */
  Orientation() :
    num_frames(0)
  {
    set(Quaternion());
  }

  /**
   * @brief Constructs with a specific initial quaternion.
   * @param q The initial quaternion.
   */
  Orientation(const Quaternion& q) :
    num_frames(0)
  {
    set(q);
  }

  /**
   * @brief Gets the number of recorded orientation frames in the current step.
   * @return The number of frames.
   */
  int length() const { return num_frames; }

  /**
   * @brief Rotates a vector by the current (latest) orientation.
   * @param v The vector to orient.
   * @return The rotated vector.
   */
  Vector orient(const Vector& v) const {
    return rotate(v, orientations[num_frames - 1]);
  }

  /**
   * @brief Rotates a vector by an orientation at a specific historical frame index.
   * @param v The vector to orient.
   * @param i The frame index (0 being oldest, length-1 being current).
   * @return The rotated vector.
   */
  Vector orient(const Vector& v, int i) const {
    return rotate(v, orientations[i]);
  }

  /**
   * @brief Rotates a vector backward by the inverse of the current (latest) orientation.
   * @param v The vector to unorient.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector& v) const {
    return rotate(v, orientations[num_frames - 1].inverse());
  }

  /**
   * @brief Rotates a vector backward by the inverse of a specific historical orientation.
   * @param v The vector to unorient.
   * @param i The frame index.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector& v, int i) const {
    return rotate(v, orientations[i].inverse());
  }

  /**
   * @brief Rotates a list of vectors by the current (latest) orientation.
   * @param vertices The list of vectors to orient.
   * @return A new list of oriented vectors.
   */
  VertexList orient(const VertexList& vertices) const {
    VertexList r;
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(r),
      [this](auto& v) {
        return orient(v);
      });
    return r;
  }

  /**
   * @brief Rotates a list of vectors by a specific historical orientation.
   * @param vertices The list of vectors to orient.
   * @param i The frame index.
   * @return A new list of oriented vectors.
   */
  VertexList orient(const VertexList& vertices, int i) const {
    VertexList r;
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(r),
      [this, i](const auto& v) {
        return orient(v, i);
      });
    return r;
  }

  /**
   * @brief Gets the current (latest) quaternion.
   * @return The Quaternion reference.
   */
  const Quaternion& get() const {
    return orientations[num_frames - 1];
  }

  /**
   * @brief Gets the quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The Quaternion reference.
   */
  const Quaternion& get(int i) const {
    return orientations[i];
  }

  /**
   * @brief Sets the orientation, clearing all history.
   * @param q The new orientation quaternion.
   * @return Reference to the Orientation object.
   */
  Orientation& set(const Quaternion& q) {
    orientations[0] = q;
    num_frames = 1;
    return *this;
  }

  /**
   * @brief Pushes a new quaternion onto the history, tracking a motion step.
   * @param q The new rotation quaternion.
   * @return Reference to the Orientation object.
   */
  Orientation& push(const Quaternion& q) {
    if (num_frames < MAX_FRAMES) {
      orientations[num_frames++] = q;
    }
    else {
      Serial.println("Orientation full, droping frame!");
    }
    return *this;
  }

  /**
   * @brief Collapses the orientation history, retaining only the latest quaternion.
   * @details Used after rendering a motion step to reset the motion blur history.
   * @return Reference to the Orientation object.
   */
  Orientation& collapse() {
    if (num_frames > 1) {
      orientations[0] = orientations[num_frames - 1];
      num_frames = 1;
    }
    return *this;
  }

private:
  static constexpr int MAX_FRAMES = MAX_W + 1; /**< Max frames of history (related to display width). */
  std::array<Quaternion, MAX_FRAMES> orientations; /**< Storage for historical quaternions. */
  size_t num_frames; /**< The current number of active frames in history. */
};


/**
 * @brief Attempts to bisect the edges of a polyhedron that intersect a plane.
 * @details This function is used to dynamically add vertices and edges to a polyhedral wireframe
 * where it crosses a dividing plane (e.g., the horizon).
 * @tparam Poly The type of polyhedron (must expose `vertices` and `euler_path`).
 * @param poly The polyhedron to modify.
 * @param orientation The current rotation/orientation of the polyhedron.
 * @param normal The normal vector defining the dividing plane.
 */
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

/**
 * @brief Calculates a gradient color based on the vector's dot product with a normal.
 * @details This creates two color gradients extending from the dividing plane in opposite directions.
 * @param v The vector to color.
 * @param normal The plane normal.
 * @param p1 The palette for the positive side.
 * @param p2 The palette for the negative side.
 * @return The calculated gradient color (Pixel).
 */
Pixel distance_gradient(const Vector& v, const Vector& normal, CRGBPalette256 p1, CRGBPalette256 p2) {
  auto d = dot(v, normal);
  if (d > 0) {
    return p1[static_cast<int>(d * 255)];
  }
  else {
    return p2[static_cast<int>(-d * 255)];
  }
}

/**
 * @brief Structure representing a Dodecahedron (20-sided polyhedron).
 * @details Stores vertices as Vectors and edges as adjacency lists.
 */
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

  /**
   * @brief Constructor that normalizes all vertices to the unit sphere.
   */
  Dodecahedron() {
    for (auto& v : vertices) {
      v.normalize();
    }
  }

  /**
   * @brief Prints the polyhedron's data to the serial output for debugging.
   */
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

/**
 * @brief Generates a truly random 3D unit vector (direction) using Marsaglia's method.
 * @return A normalized random Vector.
 */
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

/**
 * @brief Parameters defining a Lissajous curve.
 */
struct LissajousParams {
  float m1; /**< Frequency coefficient for the axial components (X and Z). */
  float m2; /**< Frequency coefficient for the orbital component (Y). */
  float a; /**< Phase shift multiplier (multiplies PI_F). */
  float domain; /**< The total duration (t) over which the curve is drawn. */
};


/**
 * @brief Calculates a 3D point on the unit sphere corresponding to a spherical Lissajous curve.
 * @param m1 Frequency coefficient for XZ plane.
 * @param m2 Frequency coefficient for Y axis.
 * @param a Phase shift multiplier.
 * @param t Time variable (or position along the domain).
 * @return The calculated 3D point (unit vector).
 */
Vector lissajous(float m1, float m2, float a, float t) {
  Vector v(
    sinf(m2 * t) * cosf(m1 * t - a * PI_F),
    cosf(m2 * t),
    sinf(m2 * t) * sinf(m1 * t - a * PI_F)
  );
  return v.normalize();
}