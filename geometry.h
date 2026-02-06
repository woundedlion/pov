/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <FastLED.h>
#include "3dmath.h"
#include "color.h"
#include "static_circular_buffer.h"
#include "led.h"

 /**
  * @brief Represents a "Fragment" or a potential pixel/vertex with associated data registers.
  * Mirrors the JS Fragment structure for shader compatibility.
  */
struct Fragment {
  Vector pos;
  float v0 = 0.0f; /**< Register 0 (usually normalized progress t) */
  float v1 = 0.0f; /**< Register 1 (usually arc length/distance) */
  float v2 = 0.0f; /**< Register 2 (usually index/id) */
  float v3 = 0.0f; /**< Register 3 (auxiliary) */
  float age = 0.0f; /**< Age of the operation/trail */

  /**
   * @brief Linear interpolation between two fragments.
   * @param a Start fragment.
   * @param b End fragment.
   * @param t Interpolation factor (0.0 to 1.0).
   * @return The interpolated fragment.
   */
  static Fragment lerp(const Fragment& a, const Fragment& b, float t) {
    Fragment f;
    f.pos = a.pos + (b.pos - a.pos) * t; // Note: This is linear, not slerp. Slerp usually happens known externally.
    f.v0 = a.v0 + (b.v0 - a.v0) * t;
    f.v1 = a.v1 + (b.v1 - a.v1) * t;
    f.v2 = a.v2 + (b.v2 - a.v2) * t;
    f.v3 = a.v3 + (b.v3 - a.v3) * t;
    f.age = a.age + (b.age - a.age) * t;
    return f;
  }
};

/**
 * @brief Logic for handling Fragment shaders.
 */
struct ShaderResult {
  Pixel color;
  float alpha = 1.0f;
  uint8_t tag = 0;
};


/**
 * @brief A list of fragments, equivalent to 'Points' in the JS context but with full register support.
 */
using Fragments = std::vector<Fragment>;

/**
 * @brief Logic for no-op vertex shader.
 */
struct NullVertexShader {
    Fragment operator()(const Fragment& f) const { return f; }
};

/**
 * @brief Logic for no-op fragment shader.
 */
struct NullFragmentShader {
    Color4 operator()(const Vector&, const Fragment&) const { return Color4(0, 0, 0, 0); }
};

/**
 * @brief A lightweight view over a contiguous sequence of objects.
 * @tparam T The type of object (can be const).
 */



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
 * @brief Unit vector along the Cartesian Y-axis.
 */
static constexpr Vector UP = Y_AXIS;

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
   * @param color The pixel color with alpha.
   */
  Dot(const Vector& v, const Color4& color) :
    position(v),
    color(color)
  {
  }

  /**
   * @brief Copy constructor.
   * @param d The Dot to copy.
   */
  Dot(const Dot& d)
    : position(d.position), color(d.color) {
  }

  Vector position; /**< The 3D position (unit vector). */
  Color4 color; /**< The color of the dot. */
};

/**
 * @brief Type alias for a circular buffer used to store active dots/fragments.
 * @details Capacity is set to 1024.
 */
using Dots = StaticCircularBuffer<Dot, 1024>;

/**
 * @brief Type alias for a circular buffer used to store geometry points (Vectors).
 * @details Capacity is set to 1024.
 */
using Points = StaticCircularBuffer<Vector, 256>;

/**
 * @brief Struct to hold Log-Polar coordinates.
 */
struct LogPolar {
  float rho;
  float theta;
};

/**
 * @brief Converts a pixel y-coordinate to a spherical phi angle.
 * @param y The pixel y-coordinate [0, H_VIRT - 1].
 * @return The spherical phi angle in radians.
 */
float y_to_phi(float y) {
  return (y * PI_F) / (H_VIRT - 1);
}

/**
 * @brief Converts a spherical phi angle to a pixel y-coordinate.
 * @param phi The spherical phi angle in radians.
 * @return The pixel y-coordinate [0, H_VIRT - 1].
 */
float phi_to_y(float phi) {
  return (phi * (H_VIRT - 1)) / PI_F;
}

/**
 * @brief Converts 2D pixel coordinates to a 3D unit vector on the sphere.
 * @tparam W The width (number of columns) of the virtual LED display.
 * @param x The horizontal coordinate (0 to W).
 * @param y The vertical coordinate (0 to H_VIRT - 1).
 * @return The corresponding unit vector.
 */
template <int W>
struct PixelLUT {
  static std::array<Vector, W* H_VIRT> data;
  static bool initialized;
  static void init() {
    for (int y = 0; y < H_VIRT; y++) {
      for (int x = 0; x < W; x++) {
        data[x + y * W] = Vector(Spherical((x * 2 * PI_F) / W, y_to_phi(y)));
      }
    }
    initialized = true;
  }
};

template <int W> std::array<Vector, W* H_VIRT> PixelLUT<W>::data;
template <int W> bool PixelLUT<W>::initialized = false;

template <int W>
const Vector& pixel_to_vector(int x, int y) {
  if (!PixelLUT<W>::initialized) {
    PixelLUT<W>::init();
  }
  return PixelLUT<W>::data[x + y * W];
}

template <int W>
Vector pixel_to_vector(float x, float y) {
  if (std::abs(x - floor(x)) < TOLERANCE && std::abs(y - floor(y)) < TOLERANCE) {
    return pixel_to_vector<W>(static_cast<int>(x), static_cast<int>(y));
  }
  return Vector(
    Spherical(
      (x * 2 * PI_F) / W,
      y_to_phi(y)
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
  PixelCoords p({ wrap((s.theta * W) / (2 * PI_F), W), phi_to_y(s.phi) });
  return p;
}

/**
 * @brief Converts Log-Polar coordinates (rho, theta) to a vector on the unit sphere.
 * Maps: Log-Polar -> Complex Plane -> Inverse Stereographic -> Sphere
 * @param rho The log-radius (natural logarithm of the radius on the complex plane).
 * @param theta The angle in radians.
 * @return Normalized vector on the unit sphere.
 */
Vector logPolarToVector(float rho, float theta) {
  const float R = expf(rho);
  const float R2 = R * R;
  const float y = (R2 - 1.0f) / (R2 + 1.0f);
  const float r_xz = sqrtf(std::max(0.0f, 1.0f - y * y));
  return Vector(
    r_xz * cosf(theta),
    y,
    r_xz * sinf(theta)
  ).normalize();
}


/**
 * @brief Converts a vector on the unit sphere to Log-Polar coordinates.
 * Maps: Sphere -> Stereographic -> Complex Plane -> Log-Polar
 * @param v Normalized vector on the unit sphere.
 * @return Log-Polar coordinates.
 */
LogPolar vectorToLogPolar(const Vector& v) {
  const float denom = 1.0f - v.j;
  if (std::abs(denom) < 0.00001f) {
    return { 10.0f, 0.0f }; // Handle North Pole singularity
  }
  const float rho = 0.5f * logf((1.0f + v.j) / denom);
  const float theta = atan2f(v.k, v.i);
  return { rho, theta };
}

#if __cplusplus < 202002L
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
#endif

/**
 * @brief Calculates a point on the Fibonacci spiral on the unit sphere.
 * @param n The total number of points in the spiral.
 * @param eps The epsilon offset for the spiral.
 * @param i The index of the point to calculate.
 * @return The point on the unit sphere.
 */
Vector fib_spiral(int n, float eps, int i) {
  float phi = acosf(1.0f - (2.0f * (static_cast<float>(i) + eps)) / static_cast<float>(n));
  float theta = fmodf((2.0f * PI_F * static_cast<float>(i) * G), (2.0f * PI_F));
  // Y-up convention
  return Vector(
    sinf(phi) * cosf(theta),
    cosf(phi),
    sinf(phi) * sinf(theta)
  ).normalize();
}

/**
 * @brief Generates a sine wave function with offset, amplitude, frequency, and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset.
 * @return A lambda that takes time (t) and returns a float.
 */
auto sin_wave(float from, float to, float freq, float phase) {
  return [=](float t) -> float {
    auto w = (sinf(freq * t * 2 * PI_F - (PI_F / 2) + PI_F - (2 * phase)) + 1) / 2;
    return lerp(from, to, w);
    };
}

/**
 * @brief Generates a triangle wave function with offset, amplitude, frequency, and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param phase The starting phase offset.
 * @return A lambda that takes time (t) and returns a float.
 */
auto tri_wave(float from, float to, float freq, float phase) {
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
 * @brief Generates a square wave function with offset, amplitude, frequency, duty cycle, and phase.
 * @param from The minimum output value.
 * @param to The maximum output value.
 * @param freq The frequency (cycles per unit time).
 * @param dutyCycle The percentage of time the wave is "on" (high).
 * @param phase The starting phase offset.
 * @return A lambda that takes time (t) and returns a float.
 */
auto square_wave(float from, float to, float freq, float dutyCycle, float phase) {
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

  /**
   * @brief Access a mutable quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The Quaternion reference.
   */
  Quaternion& at(int i) {
    return orientations[i];
  }

  /**
   * @brief Increases the resolution of the history to 'count' steps, preserving shape via Slerp.
   * @param count The target number of steps in the history.
   */
  void upsample(int count) {
    if (num_frames >= count) return;
    if (count > MAX_FRAMES) count = MAX_FRAMES;

    std::array<Quaternion, MAX_FRAMES> old_orientations;
    std::copy(orientations.begin(), orientations.begin() + num_frames, old_orientations.begin());

    int old_num_frames = num_frames;
    num_frames = count;
    orientations[0] = old_orientations[0];
    orientations[count - 1] = old_orientations[old_num_frames - 1];

    for (int i = 1; i < count - 1; i++) {
      float t = static_cast<float>(i) / (count - 1);
      float old_val = t * (old_num_frames - 1);
      int idx_a = static_cast<int>(std::floor(old_val));
      int idx_b = static_cast<int>(std::ceil(old_val));
      float alpha = old_val - idx_a;
      orientations[i] = slerp(old_orientations[idx_a], old_orientations[idx_b], alpha);
    }
  }

private:
  static constexpr int MAX_FRAMES = MAX_W + 1; /**< Max frames of history (related to display width). */
  std::array<Quaternion, MAX_FRAMES> orientations; /**< Storage for historical quaternions. */
  int num_frames; /**< The current number of active frames in history. */
};


/**
 * @brief Helper to iterate over an Orientation's historical frames.
 * @param o The orientation to iterate.
 * @param callback The function to call for each frame: `void(const Quaternion&, float t)`.
 */
template <typename F>
void tween(const Orientation& o, F callback) {
    int len = o.length();
    if (len <= 1) {
        callback(o.get(), 1.0f);
        return;
    }
    for (int i = 0; i < len; ++i) {
        float t = static_cast<float>(i) / (len - 1);
        callback(o.get(i), t);
    }
}

/**
 * @brief Helper to iterate over any Tweenable object (Orientation or OrientationTrail).
 * @param o The object to iterate.
 * @param callback The function to call for each step: `void(const T&, float t)`.
 */
template <typename T, typename F>
void deep_tween(const T& o, F callback) {
    size_t len = o.length();
    if (len <= 1) {
        callback(o.get(0), 1.0f);
        return;
    }
    for (size_t i = 0; i < len; ++i) {
        float t = static_cast<float>(i) / (len - 1);
        callback(o.get(i), t);
    }
}


/**
 * @brief Calculates a gradient color based on the vector's dot product with a normal.
 * @details This creates two color gradients extending from the dividing plane in opposite directions.
 * @param v The vector to color.
 * @param normal The plane normal.
 * @param p1 The palette for the positive side.
 * @param p2 The palette for the negative side.
 * @return The calculated gradient color (Color4).
 */
Color4 distance_gradient(const Vector& v, const Vector& normal, CRGBPalette256 p1, CRGBPalette256 p2) {
  auto d = dot(v, normal);
  if (d > 0) {
    return Color4(p1[static_cast<int>(d * 255)], 1.0f);
  }
  else {
    return Color4(p2[static_cast<int>(-d * 255)], 1.0f);
  }
}











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

/**
 * @brief Creates a basis { u, v, w } from an orientation and normal.
 */
struct Basis {
  Vector u, v, w;
};

/**
 * @brief Creates a basis { u, v, w } from an orientation and normal.
 * @param orientation The orientation quaternion.
 * @param normal The normal vector (approximate 'v' axis).
 * @return The constructed Basis.
 */
static Basis make_basis(const Quaternion& orientation, const Vector& normal) {
  Vector ref_axis = (std::abs(dot(normal, X_AXIS)) > (1 - TOLERANCE)) ? Y_AXIS : X_AXIS;
  Vector v = rotate(normal, orientation).normalize();
  Vector ref = rotate(ref_axis, orientation).normalize();
  Vector u = cross(v, ref).normalize();
  Vector w = cross(v, u).normalize();
  return { u, v, w };
}

/**
 * @brief Represents a vertex in a Half-Edge mesh.
 */
struct HEVertex;

/**
 * @brief Represents a face in a Half-Edge mesh.
 */
struct HEFace;

/**
 * @brief Represents a half-edge in a Half-Edge mesh.
 */
struct HalfEdge {
  HEVertex* vertex = nullptr; /**< Vertex at the end of this half-edge. */
  HEFace* face = nullptr;     /**< Face this half-edge belongs to. */
  HalfEdge* next = nullptr;   /**< Next half-edge in the face loop. */
  HalfEdge* prev = nullptr;   /**< Previous half-edge in the face loop. */
  HalfEdge* pair = nullptr;   /**< Opposite half-edge. */
};

struct HEVertex {
  Vector position;
  HalfEdge* halfEdge = nullptr; /**< One of the half-edges pointing to this vertex. */
};

struct HEFace {
  HalfEdge* halfEdge = nullptr; /**< One of the half-edges bordering this face. */
};

/**
 * @brief Forward declaration of MeshState for HalfEdgeMesh.
 */
struct MeshState;

/**
 * @brief A simple dynamic mesh structure compatible with MeshOps templates.
 */
struct PolyMesh {
  std::vector<Vector> vertices;
  std::vector<std::vector<int>> faces;
};
class HalfEdgeMesh {
public:
  std::vector<HEVertex> vertices;
  std::vector<HEFace> faces;
  std::vector<HalfEdge> halfEdges; // Stored contiguously to allow pointers

  /**
   * @brief Constructs a HalfEdgeMesh from a standard mesh.
   * @param mesh The input mesh (vertices and faces).
   */
  HalfEdgeMesh(const MeshState& mesh);

  template <typename MeshT>
  HalfEdgeMesh(const MeshT& mesh) {
    // 1. Create Vertices
    vertices.resize(mesh.vertices.size());
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
      vertices[i].position = mesh.vertices[i];
    }

    // 2. Create Faces and HalfEdges
    size_t idx_offset = 0;
    std::map<std::pair<int, int>, HalfEdge*> edgeMap; // For pairing half-edges
    for (size_t f_idx = 0; f_idx < mesh.num_faces; ++f_idx) {
      faces.emplace_back();
      HEFace* currentFace = &faces.back();
      
      size_t count = mesh.face_counts[f_idx];
      size_t faceStartHeIdx = halfEdges.size();
      
      // Allocate edges for this face
      for (size_t i = 0; i < count; ++i) {
        halfEdges.emplace_back();
      }

      for (size_t i = 0; i < count; ++i) {
        int u = mesh.faces[idx_offset + i];
        int v = mesh.faces[idx_offset + (i + 1) % count];

        HalfEdge* he = &halfEdges[faceStartHeIdx + i];
        
        // Link basic geometry
        he->vertex = &vertices[v]; // Points TO v
        he->face = currentFace;
        
        // Circular links
        he->next = &halfEdges[faceStartHeIdx + (i + 1) % count];
        he->prev = &halfEdges[faceStartHeIdx + (i - 1 + count) % count];
        
        // Vertex ref (just needs one incoming edge)
        vertices[v].halfEdge = he;
        
        // Pair lookup
        if (edgeMap.count({v, u})) {
          HalfEdge* neighbor = edgeMap[{v, u}];
          he->pair = neighbor;
          neighbor->pair = he;
        } else {
          edgeMap[{u, v}] = he;
        }
      }
      currentFace->halfEdge = &halfEdges[faceStartHeIdx];
      idx_offset += count;
    }
  }
};


/**
 * @brief Represents the state of a mesh using static storage to avoid heap allocations.
 * @details Max vertices set to 64 to accommodate Truncated Icosahedron / Snub Dodecahedron (60 vertices).
 */
#include <memory>
#include <array>
#include <algorithm> // For std::swap
#include <cfloat>

// --- Spatial Structures ---
// (Moved to spatial.h)
#include "spatial.h"


/**
 * @brief Implementation of HalfEdgeMesh constructor for MeshState.
 */
inline HalfEdgeMesh::HalfEdgeMesh(const MeshState& mesh) {
  // 1. Create Vertices
  vertices.resize(mesh.num_vertices);
  for (size_t i = 0; i < mesh.num_vertices; ++i) {
    vertices[i].position = mesh.vertices[i];
  }

  // 2. Create Faces and HalfEdges
  faces.reserve(mesh.num_faces);
  std::map<std::pair<int, int>, HalfEdge*> edgeMap;

  const int* face_ptr = mesh.faces;

  for (size_t k = 0; k < mesh.num_faces; ++k) {
    size_t count = mesh.face_counts[k];
    
    faces.emplace_back();
    HEFace* currentFace = &faces.back();
    
    size_t faceStartHeIdx = halfEdges.size();
    
    // Allocate edges for this face
    for (size_t i = 0; i < count; ++i) {
      halfEdges.emplace_back();
    }

    for (size_t i = 0; i < count; ++i) {
      int u = face_ptr[i];
      int v = face_ptr[(i + 1) % count];

      HalfEdge* he = &halfEdges[faceStartHeIdx + i];
      
      // Link basic geometry
      he->vertex = &vertices[v]; // Points TO v
      he->face = currentFace;
      
      // Circular links
      he->next = &halfEdges[faceStartHeIdx + (i + 1) % count];
      he->prev = &halfEdges[faceStartHeIdx + (i - 1 + count) % count];
      
      // Vertex ref (just needs one incoming edge)
      vertices[v].halfEdge = he;
      
      // Pair lookup
      if (edgeMap.count({v, u})) {
        HalfEdge* neighbor = edgeMap[{v, u}];
        he->pair = neighbor;
        neighbor->pair = he;
      } else {
        edgeMap[{u, v}] = he;
      }
    }
    currentFace->halfEdge = &halfEdges[faceStartHeIdx];
    
    face_ptr += count;
  }
}

/**
 * @brief Pre-allocated buffer for morphing operations.
 * @details Stores the source and destination states and the interpolated paths.
 */
struct MorphBuffer {
  MeshState source;
  MeshState dest;

  struct Path {
    Vector start;
    Vector end;
    Vector axis;
    float angle;
  };
  std::array<Path, MeshState::MAX_VERTS> paths;

  // Storage for dynamic topology
  static constexpr size_t MAX_FACES = 1024;
  static constexpr size_t MAX_INDICES = 8192;
  
  std::array<uint8_t, MAX_FACES> dest_face_counts_storage;
  std::array<int, MAX_INDICES> dest_faces_storage;

  void load_dest(const PolyMesh& mesh) {
     dest.clear();
     if (mesh.vertices.size() > MeshState::MAX_VERTS) return; 
     
     dest.num_vertices = mesh.vertices.size();
     std::copy(mesh.vertices.begin(), mesh.vertices.end(), dest.vertices.begin());
     
     if (mesh.faces.size() > MAX_FACES) return;
     dest.num_faces = mesh.faces.size();
     dest.face_counts = dest_face_counts_storage.data();
     
     size_t idx_offset = 0;
     for(size_t i=0; i<mesh.faces.size(); ++i) {
         size_t n = mesh.faces[i].size();
         dest_face_counts_storage[i] = static_cast<uint8_t>(n);
         if (idx_offset + n > MAX_INDICES) return; // Overflow
         for(size_t k=0; k<n; ++k) {
             dest_faces_storage[idx_offset + k] = mesh.faces[i][k];
         }
         idx_offset += n;
     }
     dest.faces = dest_faces_storage.data();
  }

  void init_paths() {
      // Logic handled in MeshMorph
  }
};

/**
 * @brief Structure returned by compile_hankin.
 */
struct HankinInstruction {
  Vector pCorner;    /**< Corner vertex position. */
  Vector pPrev;      /**< Previous vertex position. */
  Vector pNext;      /**< Next vertex position. */
  int idxM1;         /**< Index of first midpoint (static vertex). */
  int idxM2;         /**< Index of second midpoint (static vertex). */
};

/**
 * @brief Compiled topological data for fast Hankin pattern updates.
 */
struct CompiledHankin {
  std::vector<Vector> staticVertices;         /**< Midpoints that don't move. */
  std::vector<Vector> dynamicVertices;        /**< Intersection points that move. */
  std::vector<HankinInstruction> dynamicInstructions; /**< Instructions to update dynamic vertices. */
  std::vector<std::vector<int>> faces;        /**< Resulting face topology. */
  int staticOffset;                           /**< Offset where dynamic vertices start. */
};

/**
 * @brief Operations on meshes (Dual, Hankin, etc.).
 */
struct MeshOps {
  /**
   * @brief Computes the dual of a mesh.
   */
  template <typename MeshT>
  static MeshT dual(const MeshT& mesh) {
    HalfEdgeMesh heMesh(mesh);
    MeshT dualMesh;
    std::unordered_map<HEFace*, int> faceToVertIdx;
    
    // New Vertices (Centroids of original faces)
    for (size_t i = 0; i < heMesh.faces.size(); ++i) {
      HEFace* face = &heMesh.faces[i];
      Vector c(0, 0, 0);
      int count = 0;
      HalfEdge* he = face->halfEdge;
      HalfEdge* start = he;
      do {
        c = c + he->vertex->position;
        count++;
        he = he->next;
      } while (he != start);
      
      c = c / static_cast<float>(count);
      dualMesh.vertices.push_back(c.normalize());
      faceToVertIdx[face] = static_cast<int>(i);
    }
    
    // New Faces (Cycles around original vertices)
    std::vector<HEVertex*> visitedVerts; // Naive set replacement
    
    // Helper to check if visited
    auto isVisited = [&](HEVertex* v) {
      return std::find(visitedVerts.begin(), visitedVerts.end(), v) != visitedVerts.end();
    };

    for (auto& heStart : heMesh.halfEdges) {
      if (!heStart.prev) continue; // Safety
      
      HEVertex* origin = heStart.prev->vertex; // The vertex 'start' comes FROM
      if (isVisited(origin)) continue;
      visitedVerts.push_back(origin);
      
      std::vector<int> faceIndices;
      HalfEdge* curr = &heStart;
      HalfEdge* startOrbit = curr;
      int safety = 0;
      
      do {
        if (!curr->face) break;
        faceIndices.push_back(faceToVertIdx[curr->face]);
        
        if (!curr->pair) break;
        curr = curr->pair->next;
        safety++;
      } while (curr != startOrbit && curr && safety < 100);
      
      if (faceIndices.size() > 2) {
        std::reverse(faceIndices.begin(), faceIndices.end()); // Maintain CCW
        dualMesh.faces.push_back(faceIndices);
      }
    }
    
    return dualMesh;
  }

  /**
   * @brief Compiles the topology for a Hankin pattern.
   */
  template <typename MeshT>
  static CompiledHankin compile_hankin(const MeshT& mesh) {
    HalfEdgeMesh heMesh(mesh);
    CompiledHankin compiled;
    
    std::map<HalfEdge*, int> heToMidpointIdx;
    std::map<HalfEdge*, int> heToDynamicIdx;
    
    // Helper to get/create midpoint index
    auto getMidpointIdx = [&](HalfEdge* he) {
      if (heToMidpointIdx.count(he)) return heToMidpointIdx[he];
      if (he->pair && heToMidpointIdx.count(he->pair)) return heToMidpointIdx[he->pair];
      
      Vector pA = he->prev ? he->prev->vertex->position : he->pair->vertex->position;
      Vector pB = he->vertex->position;
      Vector mid = (pA + pB) * 0.5f;
      mid.normalize();
      
      compiled.staticVertices.push_back(mid);
      int idx = static_cast<int>(compiled.staticVertices.size()) - 1;
      heToMidpointIdx[he] = idx;
      if (he->pair) heToMidpointIdx[he->pair] = idx;
      return idx;
    };
    
    // Ensure all midpoints
    for (auto& he : heMesh.halfEdges) {
      getMidpointIdx(&he);
    }
    
    compiled.staticOffset = static_cast<int>(compiled.staticVertices.size());
    
    // Star faces
    for (auto& face : heMesh.faces) {
      std::vector<int> starFaceIndices;
      HalfEdge* he = face.halfEdge;
      HalfEdge* startHe = he;
      
      do {
        HalfEdge* prev = he->prev;
        HalfEdge* curr = he;
        
        int idxM1 = getMidpointIdx(prev);
        int idxM2 = getMidpointIdx(curr);
        
        Vector pCorner = prev->vertex->position;
        Vector pPrev = (prev->prev ? prev->prev->vertex->position : prev->pair->vertex->position);
        Vector pNext = curr->vertex->position;
        
        compiled.dynamicInstructions.push_back({ pCorner, pPrev, pNext, idxM1, idxM2 });
        
        int dynIdx = static_cast<int>(compiled.dynamicVertices.size());
        heToDynamicIdx[curr] = dynIdx;
        compiled.dynamicVertices.emplace_back(); // Placeholder
        
        starFaceIndices.push_back(idxM1);
        starFaceIndices.push_back(compiled.staticOffset + dynIdx);
        
        he = he->next;
      } while (he != startHe);
      
      compiled.faces.push_back(starFaceIndices);
    }
    
    // Rosette faces
    std::vector<HEVertex*> visitedVerts;
    auto isVisited = [&](HEVertex* v) {
      return std::find(visitedVerts.begin(), visitedVerts.end(), v) != visitedVerts.end();
    };
    
    for (auto& heStart : heMesh.halfEdges) {
      if (!heStart.prev) continue;
      HEVertex* origin = heStart.prev->vertex;
      if (isVisited(origin)) continue;
      visitedVerts.push_back(origin);
      
      std::vector<int> rosetteIndices;
      HalfEdge* curr = &heStart;
      HalfEdge* startOrbit = curr;
      int safety = 0;
      
      do {
        rosetteIndices.push_back(heToMidpointIdx[curr]); // Static
        HalfEdge* nextEdge = curr->pair ? curr->pair->next : nullptr;
        if (!nextEdge) break;
        rosetteIndices.push_back(compiled.staticOffset + heToDynamicIdx[nextEdge]); // Dynamic
        curr = nextEdge;
        safety++;
      } while (curr != startOrbit && curr && safety < 100);
      
      if (rosetteIndices.size() > 2) {
        std::reverse(rosetteIndices.begin(), rosetteIndices.end());
        compiled.faces.push_back(rosetteIndices);
      }
    }
    
    return compiled;
  }
  
  /**
   * @brief Updates a compiled Hankin pattern.
   */
  template <typename MeshT> // Templated just to match style, though MeshT output is expected
  static MeshT update_hankin(CompiledHankin& compiled, float angle) {
    for (size_t i = 0; i < compiled.dynamicInstructions.size(); ++i) {
      const auto& instr = compiled.dynamicInstructions[i];
      Vector m1 = compiled.staticVertices[instr.idxM1];
      Vector m2 = compiled.staticVertices[instr.idxM2];
      
      Vector nEdge1 = cross(instr.pPrev, instr.pCorner).normalize();
      Quaternion q1 = make_rotation(m1, angle);
      Vector nHankin1 = rotate(nEdge1, q1);
      
      Vector nEdge2 = cross(instr.pCorner, instr.pNext).normalize();
      Quaternion q2 = make_rotation(m2, -angle);
      Vector nHankin2 = rotate(nEdge2, q2);
      
      Vector intersect = cross(nHankin1, nHankin2);
      if (dot(intersect, instr.pCorner) < 0) intersect = -intersect;
      
      compiled.dynamicVertices[i] = intersect.normalize();
    }
    
    MeshT result;
    result.vertices = compiled.staticVertices;
    result.vertices.insert(result.vertices.end(), compiled.dynamicVertices.begin(), compiled.dynamicVertices.end());
    result.faces = compiled.faces;
    return result;
  }
  
  /**
   * @brief Helper to do full hankin generation.
   */
  template <typename MeshT>
  static MeshT hankin(const MeshT& mesh, float angle) {
    auto compiled = compile_hankin(mesh);
    return update_hankin<MeshT>(compiled, angle);
  }

  // --- CONWAY OPERATORS ---

  /**
   * @brief Normalizes all vertices in the mesh to the unit sphere.
   */
  template <typename MeshT>
  static void normalize(MeshT& mesh) {
    for (auto& v : mesh.vertices) {
      v = v.normalize();
    }
  }

  /**
   * @brief Kis operator: Raises a pyramid on each face.
   */
  template <typename MeshT>
  static MeshT kis(const MeshT& mesh) {
    MeshT result;
    result.vertices = mesh.vertices; // Copy existing
    
    for (const auto& f : mesh.faces) {
      // Add centroid
      Vector centroid(0, 0, 0);
      for (int vi : f) {
        centroid = centroid + mesh.vertices[vi];
      }
      if (!f.empty()) centroid = centroid / static_cast<float>(f.size());
      
      result.vertices.push_back(centroid);
      int centerIdx = static_cast<int>(result.vertices.size()) - 1;

      // Create triangles
      for (size_t i = 0; i < f.size(); ++i) {
        int vi = f[i];
        int vj = f[(i + 1) % f.size()];
        result.faces.push_back({ vi, vj, centerIdx });
      }
    }

    normalize(result);
    return result;
  }

  /**
   * @brief Ambo operator: Truncates vertices to edge midpoints.
   */
  template <typename MeshT>
  static MeshT ambo(const MeshT& mesh) {
    MeshT result;
    std::map<std::pair<int, int>, int> edgeMap;

    // 1. Create vertices at edge midpoints
    for (const auto& f : mesh.faces) {
      for (size_t i = 0; i < f.size(); ++i) {
        int vi = f[i];
        int vj = f[(i + 1) % f.size()];
        int u = std::min(vi, vj);
        int v = std::max(vi, vj);
        
        if (edgeMap.find({u, v}) == edgeMap.end()) {
          Vector mid = (mesh.vertices[vi] + mesh.vertices[vj]) * 0.5f;
          result.vertices.push_back(mid);
          edgeMap[{u, v}] = static_cast<int>(result.vertices.size()) - 1;
        }
      }
    }

    // 2. Create faces
    // A. Shrink old faces
    for (const auto& f : mesh.faces) {
      std::vector<int> faceVerts;
      for (size_t i = 0; i < f.size(); ++i) {
        int vi = f[i];
        int vj = f[(i + 1) % f.size()];
        int u = std::min(vi, vj);
        int v = std::max(vi, vj);
        faceVerts.push_back(edgeMap[{u, v}]);
      }
      result.faces.push_back(faceVerts);
    }

    // B. Create new faces at old vertices
    std::map<std::pair<int, int>, std::vector<int>> edgeToFaces;
    for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
      const auto& f = mesh.faces[fi];
      for (size_t i = 0; i < f.size(); ++i) {
        int vi = f[i];
        int vj = f[(i + 1) % f.size()];
        int u = std::min(vi, vj);
        int v = std::max(vi, vj);
        edgeToFaces[{u, v}].push_back(static_cast<int>(fi));
      }
    }

    for (size_t vi = 0; vi < mesh.vertices.size(); ++vi) {
      std::vector<int> neighborMids;
      
      // Find a start face
      int startFaceIdx = -1;
      for(size_t fi = 0; fi < mesh.faces.size(); ++fi) {
        const auto& f = mesh.faces[fi];
        if (std::find(f.begin(), f.end(), static_cast<int>(vi)) != f.end()) {
          startFaceIdx = static_cast<int>(fi);
          break;
        }
      }
      if (startFaceIdx == -1) continue;

      int currFaceIdx = startFaceIdx;
      int safety = 0;
      do {
        const auto& face = mesh.faces[currFaceIdx];
        auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
        size_t idxInFace = std::distance(face.begin(), it);
        int nextVi = face[(idxInFace + 1) % face.size()];
        
        int u = std::min((int)vi, nextVi);
        int v = std::max((int)vi, nextVi);
        neighborMids.push_back(edgeMap[{u, v}]);

        // Find adjacent face
        const auto& adjFaces = edgeToFaces[{u, v}];
        auto nextFaceIt = std::find_if(adjFaces.begin(), adjFaces.end(), [&](int id){ return id != currFaceIdx; });
        if (nextFaceIt == adjFaces.end()) break;

        currFaceIdx = *nextFaceIt;
        safety++;
      } while (currFaceIdx != startFaceIdx && safety < 20);

      if (neighborMids.size() >= 3) {
        std::reverse(neighborMids.begin(), neighborMids.end());
        result.faces.push_back(neighborMids);
      }
    }

    normalize(result);
    return result;
  }

  /**
   * @brief Truncate operator: Cuts corners off the polyhedron.
   * @param t Truncation depth [0..0.5].
   */
  template <typename MeshT>
  static MeshT truncate(const MeshT& mesh, float t = 0.25f) {
    MeshT result;
    // Map edge (u,v) -> pair of new vertex indices {near_u, near_v}
    // Stored as key {min(u,v), max(u,v)} -> value {idx_near_key_first, idx_near_key_second}
    std::map<std::pair<int, int>, std::pair<int, int>> edgeMap;

    // 1. Create new vertices along edges
    for (const auto& f : mesh.faces) {
      for (size_t i = 0; i < f.size(); ++i) {
        int u = f[i];
        int v = f[(i + 1) % f.size()];
        int k1 = std::min(u, v);
        int k2 = std::max(u, v);
        
        if (edgeMap.find({k1, k2}) == edgeMap.end()) {
           Vector vU = mesh.vertices[u];
           Vector vV = mesh.vertices[v];
           
           // Near U
           // lerp implementation in 3dmath/geometry? Using manual: a + (b-a)*t
           Vector p1 = vU + (vV - vU) * t;
           result.vertices.push_back(p1);
           int idx1 = static_cast<int>(result.vertices.size()) - 1;
           
           // Near V
           Vector p2 = vU + (vV - vU) * (1.0f - t);
           result.vertices.push_back(p2);
           int idx2 = static_cast<int>(result.vertices.size()) - 1;
           
           if (u < v) {
             edgeMap[{k1, k2}] = {idx1, idx2};
           } else {
             edgeMap[{k1, k2}] = {idx2, idx1}; 
           }
           if (u < v) edgeMap[{u, v}] = {idx1, idx2};
           else       edgeMap[{v, u}] = {idx2, idx1};
        }
      }
    }

    // 2. Modified Faces (internal polygons)
    for (const auto& f : mesh.faces) {
      std::vector<int> faceVerts;
      for (size_t i = 0; i < f.size(); ++i) {
        int u = f[i];
        int v = f[(i + 1) % f.size()];
        int k1 = std::min(u, v);
        int k2 = std::max(u, v);
        
        std::pair<int, int> indices = edgeMap[{k1, k2}];

        if (u < v) {
          faceVerts.push_back(indices.first);
          faceVerts.push_back(indices.second);
        } else {
          faceVerts.push_back(indices.second);
          faceVerts.push_back(indices.first);
        }
      }
      result.faces.push_back(faceVerts);
    }
    
    // 3. Corner Faces
    std::map<std::pair<int, int>, std::vector<int>> edgeToFaces;
    for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
      const auto& f = mesh.faces[fi];
      for (size_t i = 0; i < f.size(); ++i) {
        int u = f[i];
        int v = f[(i + 1) % f.size()];
        int k1 = std::min(u, v);
        int k2 = std::max(u, v);
        edgeToFaces[{k1, k2}].push_back(static_cast<int>(fi));
      }
    }

    for (size_t vi = 0; vi < mesh.vertices.size(); ++vi) {
      int startFaceIdx = -1;
      for(size_t fi = 0; fi < mesh.faces.size(); ++fi) {
        const auto& f = mesh.faces[fi];
        if (std::find(f.begin(), f.end(), static_cast<int>(vi)) != f.end()) {
          startFaceIdx = static_cast<int>(fi);
          break;
        }
      }
      if (startFaceIdx == -1) continue;

      std::vector<int> polyVerts;
      int currFaceIdx = startFaceIdx;
      int safety = 0;
      do {
        const auto& face = mesh.faces[currFaceIdx];
        auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
        size_t idxInFace = std::distance(face.begin(), it);
        int nextVi = face[(idxInFace + 1) % face.size()];

        int k1 = std::min((int)vi, nextVi);
        int k2 = std::max((int)vi, nextVi);
        std::pair<int, int> indices = edgeMap[{k1, k2}];

        int idxNearVi = (static_cast<int>(vi) == k1) ? indices.first : indices.second;
        polyVerts.push_back(idxNearVi);

        // Find adjacent face
        const auto& adjFaces = edgeToFaces[{k1, k2}];
        auto nextFaceIt = std::find_if(adjFaces.begin(), adjFaces.end(), [&](int id){ return id != currFaceIdx; });
        if (nextFaceIt == adjFaces.end()) break;

        currFaceIdx = *nextFaceIt;
        safety++;
      } while (currFaceIdx != startFaceIdx && safety < 20);

      if (polyVerts.size() > 2) {
        std::reverse(polyVerts.begin(), polyVerts.end());
        result.faces.push_back(polyVerts);
      }
    }

    normalize(result);
    return result;
  }

  /**
   * @brief Snub operator: Creates a chiral semi-regular polyhedron.
   */
  template <typename MeshT>
  static MeshT snub(const MeshT& mesh) {
     MeshT result;
     // Structure: newVertsMap[faceIndex][vertIndexInFace] = globalIndex
     std::vector<std::vector<int>> newVertsMap(mesh.faces.size());
     const float SHRINK_FACTOR = 0.5f;

     // 1. Create new vertices (n per face)
     for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
         const auto& f = mesh.faces[fi];
         Vector centroid(0,0,0);
         for(int vi : f) centroid = centroid + mesh.vertices[vi];
         if (!f.empty()) centroid = centroid / static_cast<float>(f.size());

         newVertsMap[fi].resize(f.size());
         for(size_t i=0; i<f.size(); ++i) {
             Vector v = mesh.vertices[f[i]];
             // lerp(v, centroid, SHRINK_FACTOR)
             Vector newV = v + (centroid - v) * SHRINK_FACTOR;
             result.vertices.push_back(newV);
             newVertsMap[fi][i] = static_cast<int>(result.vertices.size()) - 1;
         }
     }

     // 2. Create "Face Faces" (shrunk originals)
     for(size_t fi = 0; fi < mesh.faces.size(); ++fi) {
         result.faces.push_back(newVertsMap[fi]);
     }

     // Helper: Build edge map
     std::map<std::pair<int, int>, std::vector<int>> edgeToFaces;
     for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
       const auto& f = mesh.faces[fi];
       for (size_t i = 0; i < f.size(); ++i) {
         int u = f[i];
         int v = f[(i + 1) % f.size()];
         int k1 = std::min(u, v);
         int k2 = std::max(u, v);
         edgeToFaces[{k1, k2}].push_back(static_cast<int>(fi));
       }
     }

     // 3. Create "Vertex Faces"
     for (size_t vi = 0; vi < mesh.vertices.size(); ++vi) {
         int startFaceIdx = -1;
         for(size_t fi = 0; fi < mesh.faces.size(); ++fi) {
           const auto& f = mesh.faces[fi];
           if (std::find(f.begin(), f.end(), static_cast<int>(vi)) != f.end()) {
             startFaceIdx = static_cast<int>(fi);
             break;
           }
         }
         if (startFaceIdx == -1) continue;

         std::vector<int> orderedFaces;
         int currFaceIdx = startFaceIdx;
         int safety = 0;
         do {
             orderedFaces.push_back(currFaceIdx);
             const auto& face = mesh.faces[currFaceIdx];
             auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
             size_t idxInFace = std::distance(face.begin(), it);
             // Previous vertex in cycle: (idx - 1)
             int prevVi = face[(idxInFace - 1 + face.size()) % face.size()];
             
             int k1 = std::min((int)vi, prevVi);
             int k2 = std::max((int)vi, prevVi);
             const auto& adjFaces = edgeToFaces[{k1, k2}];
             auto nextFaceIt = std::find_if(adjFaces.begin(), adjFaces.end(), [&](int id){ return id != currFaceIdx; });
             if (nextFaceIt == adjFaces.end()) break;

             currFaceIdx = *nextFaceIt;
             safety++;
         } while(currFaceIdx != startFaceIdx && safety < 20);

         std::vector<int> faceVerts;
         for(int fi : orderedFaces) {
             const auto& face = mesh.faces[fi];
             auto it = std::find(face.begin(), face.end(), static_cast<int>(vi));
             size_t idx = std::distance(face.begin(), it);
             faceVerts.push_back(newVertsMap[fi][idx]);
         }
         result.faces.push_back(faceVerts);
     }

     // 4. Create "Edge Triangles"
     std::set<std::pair<int, int>> processedEdges;
     for (size_t fi = 0; fi < mesh.faces.size(); ++fi) {
         const auto& f = mesh.faces[fi];
         for (size_t i = 0; i < f.size(); ++i) {
             int vi = f[i];
             int vj = f[(i + 1) % f.size()];
             int k1 = std::min(vi, vj);
             int k2 = std::max(vi, vj);
             
             if (processedEdges.count({k1, k2})) continue;
             processedEdges.insert({k1, k2});
             
             const auto& adj = edgeToFaces[{k1, k2}];
             if (adj.size() < 2) continue;
             
             int faceA = fi;
             int faceB = (adj[0] == static_cast<int>(fi)) ? adj[1] : adj[0];
             
             // Get local indices of u=vi, v=vj
             auto itA_u = std::find(mesh.faces[faceA].begin(), mesh.faces[faceA].end(), vi);
             int idxA_u = static_cast<int>(std::distance(mesh.faces[faceA].begin(), itA_u));
             
             auto itA_v = std::find(mesh.faces[faceA].begin(), mesh.faces[faceA].end(), vj);
             int idxA_v = static_cast<int>(std::distance(mesh.faces[faceA].begin(), itA_v));

             auto itB_u = std::find(mesh.faces[faceB].begin(), mesh.faces[faceB].end(), vi);
             int idxB_u = static_cast<int>(std::distance(mesh.faces[faceB].begin(), itB_u));

             auto itB_v = std::find(mesh.faces[faceB].begin(), mesh.faces[faceB].end(), vj);
             int idxB_v = static_cast<int>(std::distance(mesh.faces[faceB].begin(), itB_v));

             int A_u = newVertsMap[faceA][idxA_u];
             int A_v = newVertsMap[faceA][idxA_v];
             int B_u = newVertsMap[faceB][idxB_u];
             int B_v = newVertsMap[faceB][idxB_v];

             // Tri 1
             result.faces.push_back({A_v, A_u, B_v});
             // Tri 2
             result.faces.push_back({B_u, B_v, A_u});
         }
     }

     normalize(result);
     return result;
  }

  /**
   * @brief Gyro operator: dual(snub(mesh)).
   */
  template <typename MeshT>
  static MeshT gyro(const MeshT& mesh) {
    return dual(snub(mesh));
  }
  
  /**
   * @brief Finds the closest point on the mesh graph (BFS).
   * @param p The query point.
   * @param mesh The mesh to search.
   * @return The position of the closest vertex.
   */
  template <typename MeshT>
  static Vector closest_point_on_mesh_graph(const Vector& p, const MeshT& mesh) {
    if (mesh.vertices.empty()) return Vector(0,1,0);
    
    // Linear scan is probably faster for small meshes (<1000 verts) than building a graph
    // Holosphere meshes are typically small.
    Vector best = mesh.vertices[0];
    float minSq = distance_squared(p, best);
    
    for (const auto& v : mesh.vertices) {
      float sq = distance_squared(p, v);
      if (sq < minSq) {
        minSq = sq;
        best = v;
      }
    }
    return best;
  }


};
