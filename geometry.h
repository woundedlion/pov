/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>
#include <algorithm>
#include <array>
#include <map>
#include <set>
#include "platform.h"
#include "constants.h"
#include "3dmath.h"
#include "color.h"
#include "mesh.h"
#include "waves.h"
// transformers.h removed to break cycle

/**
 * @brief Represents a "Fragment" or a potential pixel/vertex with associated
 * data registers. Mirrors the JS Fragment structure for shader compatibility.
 */
struct Fragment {
  Vector pos;
  float v0 = 0.0f;   /**< Register 0 (usually normalized progress t) */
  float v1 = 0.0f;   /**< Register 1 (usually arc length/distance) */
  float v2 = 0.0f;   /**< Register 2 (usually index/id) */
  float v3 = 0.0f;   /**< Register 3 (auxiliary) */
  float size = 1.0f; /**< Size metric (e.g. radius/apothem) for normalization */
  float age = 0.0f;  /**< Age of the operation/trail */
  Color4 color = Color4(0, 0, 0, 0); /**< Output Color (RGBA) */

  /**
   * @brief Linear interpolation between two fragments.
   * @param a Start fragment.
   * @param b End fragment.
   * @param t Interpolation factor (0.0 to 1.0).
   * @return The interpolated fragment.
   */
  static Fragment lerp(const Fragment &a, const Fragment &b, float t) {
    Fragment f;
    f.pos =
        a.pos + (b.pos - a.pos) * t; // Note: This is linear, not slerp. Slerp
                                     // usually happens known externally.
    f.v0 = a.v0 + (b.v0 - a.v0) * t;
    f.v1 = a.v1 + (b.v1 - a.v1) * t;
    f.v2 = a.v2 + (b.v2 - a.v2) * t;
    f.v3 = a.v3 + (b.v3 - a.v3) * t;
    f.age = a.age + (b.age - a.age) * t;
    return f;
  }
};

/**
 * @brief A list of fragments, equivalent to 'Points' in the JS context but with
 * full register support.
 */
using Fragments = ArenaVector<Fragment>;

/**
 * @brief Logic for no-op vertex shader.
 */
struct NullVertexShader {
  void operator()(Fragment &f) const {}
};

/**
 * @brief Logic for no-op fragment shader.
 */
struct NullFragmentShader {
  Color4 operator()(const Vector &, const Fragment &) const {
    return Color4(0, 0, 0, 0);
  }
};


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
  Dot(const Vector &v, const Color4 &color) : position(v), color(color) {}

  /**
   * @brief Copy constructor.
   * @param d The Dot to copy.
   */
  Dot(const Dot &d) : position(d.position), color(d.color) {}

  Vector position; /**< The 3D position (unit vector). */
  Color4 color;    /**< The color of the dot. */
};

/**
 * @brief Type alias for a circular buffer used to store active dots/fragments.
 * @details Capacity is set to 1024.
 */
using Dots = StaticCircularBuffer<Dot, 1024>;

/**
 * @brief Type alias for a circular buffer used to store geometry points
 * (Vectors).
 * @details Capacity is set to 1024.
 */
using Points = StaticCircularBuffer<Vector, 1024>;

/**
 * @brief Struct to hold Log-Polar coordinates.
 */
struct LogPolar {
  float rho;
  float theta;
};

/**
 * @brief Converts a pixel y-coordinate to a spherical phi angle.
 * @param y The pixel y-coordinate [0, h_virt - 1].
 * @param h_virt The virtual height.
 * @return The spherical phi angle in radians.
 */
inline float y_to_phi(float y, int h_virt) { return (y * PI_F) / (h_virt - 1); }

/**
 * @brief Converts a spherical phi angle to a pixel y-coordinate.
 * @param phi The spherical phi angle in radians.
 * @param h_virt The virtual height.
 * @return The pixel y-coordinate [0, h_virt - 1].
 */
inline float phi_to_y(float phi, int h_virt) {
  return (phi * (h_virt - 1)) / PI_F;
}

template <int H> inline float phi_to_y(float phi) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  return (phi * (H_VIRT - 1)) / PI_F;
}

/**
 * @brief Precomputed lookup table for scanline phi angles.
 */
template <int H> struct PhiLUT {
  static constexpr int H_VIRT = H + hs::H_OFFSET;
  static std::array<float, H_VIRT> data;
  static bool initialized;
  static void init() {
    for (int y = 0; y < H_VIRT; y++) {
      data[y] = y_to_phi(static_cast<float>(y), H_VIRT);
    }
    initialized = true;
  }
};

template <int H> std::array<float, PhiLUT<H>::H_VIRT> PhiLUT<H>::data;
template <int H> bool PhiLUT<H>::initialized = false;

template <int H> inline float y_to_phi(int y) {
  if (!PhiLUT<H>::initialized) {
    PhiLUT<H>::init();
  }
  return PhiLUT<H>::data[y];
}

template <int H> inline float y_to_phi(float y) {
  if (std::abs(y - std::floor(y)) < TOLERANCE) {
    int iy = static_cast<int>(y);
    if (iy >= 0 && iy < PhiLUT<H>::H_VIRT) {
      return y_to_phi<H>(iy);
    }
  }
  constexpr int H_VIRT = H + hs::H_OFFSET;
  return (y * PI_F) / (H_VIRT - 1);
}


/**
 * @brief Split trig lookup tables for efficient vector reconstruction.
 * @details Caches sin/cos for theta (per column) and phi (per row) separately.
 * Reconstructs vectors with 3 multiplies instead of storing full Vectors.
 * Memory: ~(4*W + 4*H_VIRT) floats vs W*H_VIRT Vectors — a ~145x reduction.
 */
template <int W, int H> struct TrigLUT {
  static constexpr int H_VIRT = H + hs::H_OFFSET;
  static std::array<float, W> sin_theta;
  static std::array<float, W> cos_theta;
  static std::array<float, H_VIRT> sin_phi;
  static std::array<float, H_VIRT> cos_phi;
  static bool initialized;
  static void init() {
    if (!PhiLUT<H>::initialized) {
      PhiLUT<H>::init();
    }
    for (int x = 0; x < W; x++) {
      float theta = (x * 2 * PI_F) / W;
      sin_theta[x] = sinf(theta);
      cos_theta[x] = cosf(theta);
    }
    for (int y = 0; y < H_VIRT; y++) {
      float phi = PhiLUT<H>::data[y];
      sin_phi[y] = sinf(phi);
      cos_phi[y] = cosf(phi);
    }
    initialized = true;
  }
};

template <int W, int H> std::array<float, W> TrigLUT<W, H>::sin_theta;
template <int W, int H> std::array<float, W> TrigLUT<W, H>::cos_theta;
template <int W, int H>
std::array<float, TrigLUT<W, H>::H_VIRT> TrigLUT<W, H>::sin_phi;
template <int W, int H>
std::array<float, TrigLUT<W, H>::H_VIRT> TrigLUT<W, H>::cos_phi;
template <int W, int H> bool TrigLUT<W, H>::initialized = false;

/**
 * @brief Reconstruct a vector from pixel coordinates using split trig LUTs.
 * @tparam W Width.
 * @tparam H Height.
 * @param x X coordinate (column).
 * @param y Y coordinate (row).
 * @return Unit vector on the sphere.
 */
template <int W, int H> Vector pixel_to_vector(int x, int y) {
  if (!TrigLUT<W, H>::initialized) {
    TrigLUT<W, H>::init();
  }
  float sp = TrigLUT<W, H>::sin_phi[y];
  return Vector(sp * TrigLUT<W, H>::cos_theta[x], TrigLUT<W, H>::cos_phi[y],
                sp * TrigLUT<W, H>::sin_theta[x]);
}

template <int W, int H> Vector pixel_to_vector(float x, float y) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  if (std::abs(x - floor(x)) < TOLERANCE &&
      std::abs(y - floor(y)) < TOLERANCE) {
    return pixel_to_vector<W, H>(static_cast<int>(x), static_cast<int>(y));
  }
  return Vector(Spherical((x * 2 * PI_F) / W, y_to_phi<H_VIRT>(y)));
}

/**
 * @brief Converts a 3D unit vector back to 2D pixel coordinates.
 * @tparam W The width.
 * @tparam H The height.
 * @param v The input unit vector.
 * @return The 2D PixelCoords.
 */
template <int W, int H> PixelCoords vector_to_pixel(const Vector &v) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  auto s = Spherical(v);
  PixelCoords p({wrap((s.theta * W) / (2 * PI_F), W), phi_to_y(s.phi, H_VIRT)});
  return p;
}

/**
 * @brief Converts Log-Polar coordinates (rho, theta) to a vector on the unit
 * sphere. Maps: Log-Polar -> Complex Plane -> Inverse Stereographic -> Sphere
 * @param rho The log-radius (natural logarithm of the radius on the complex
 * plane).
 * @param theta The angle in radians.
 * @return Normalized vector on the unit sphere.
 */
inline Vector logPolarToVector(float rho, float theta) {
  const float R = expf(rho);
  const float R2 = R * R;
  const float y = (R2 - 1.0f) / (R2 + 1.0f);
  const float r_xz = sqrtf(std::max(0.0f, 1.0f - y * y));
  return Vector(r_xz * cosf(theta), y, r_xz * sinf(theta)).normalized();
}

/**
 * @brief Converts a vector on the unit sphere to Log-Polar coordinates.
 * Maps: Sphere -> Stereographic -> Complex Plane -> Log-Polar
 * @param v Normalized vector on the unit sphere.
 * @return Log-Polar coordinates.
 */
inline LogPolar vectorToLogPolar(const Vector &v) {
  const float denom = 1.0f - v.y;
  if (std::abs(denom) < 0.00001f) {
    return {10.0f, 0.0f}; // Handle North Pole singularity
  }
  const float rho = 0.5f * logf((1.0f + v.y) / denom);
  const float theta = atan2f(v.z, v.x);
  return {rho, theta};
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
inline Vector fib_spiral(int n, float eps, int i) {
  float phi = acosf(1.0f - (2.0f * (static_cast<float>(i) + eps)) /
                               static_cast<float>(n));
  float theta = fmodf((2.0f * PI_F * static_cast<float>(i) * G), (2.0f * PI_F));
  // Y-up convention
  return Vector(sinf(phi) * cosf(theta), cosf(phi), sinf(phi) * sinf(theta))
      .normalized();
}


/**
 * @brief Class managing the current rotation state of an object, maintaining
 * history for interpolation.
 * @details Stores a list of Quaternions (`orientations`) generated during the
 * current frame step.
 */
template <int W, int CAP> class Orientation {
public:
  static constexpr int CAPACITY = CAP;
  /**
   * @brief Default constructor (identity rotation).
   */
  Orientation() : num_frames(0) { set(Quaternion()); }

  /**
   * @brief Constructs with a specific initial quaternion.
   * @param q The initial quaternion.
   */
  Orientation(const Quaternion &q) : num_frames(0) { set(q); }

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
  Vector orient(const Vector &v) const {
    return rotate(v, orientations[num_frames - 1]);
  }

  /**
   * @brief Rotates a vector by an orientation at a specific historical frame
   * index.
   * @param v The vector to orient.
   * @param i The frame index (0 being oldest, length-1 being current).
   * @return The rotated vector.
   */
  Vector orient(const Vector &v, int i) const {
    return rotate(v, orientations[i]);
  }

  /**
   * @brief Rotates a vector backward by the inverse of the current (latest)
   * orientation.
   * @param v The vector to unorient.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v) const {
    return rotate(v, orientations[num_frames - 1].conjugate());
  }

  /**
   * @brief Rotates a vector backward by the inverse of a specific historical
   * orientation.
   * @param v The vector to unorient.
   * @param i The frame index.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v, int i) const {
    return rotate(v, orientations[i].conjugate());
  }

  /**
   * @brief Gets the current (latest) quaternion.
   * @return The Quaternion reference.
   */
  const Quaternion &get() const { return orientations[num_frames - 1]; }

  /**
   * @brief Gets the quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The Quaternion reference.
   */
  const Quaternion &get(int i) const { return orientations[i]; }

  /**
   * @brief Sets the orientation, clearing all history.
   * @param q The new orientation quaternion.
   * @return Reference to the Orientation object.
   */
  Orientation &set(const Quaternion &q) {
    orientations[0] = q;
    num_frames = 1;
    return *this;
  }

  /**
   * @brief Pushes a new quaternion onto the history, tracking a motion step.
   * @param q The new rotation quaternion.
   * @return Reference to the Orientation object.
   */
  Orientation &push(const Quaternion &q) {
    if (num_frames < CAPACITY) {
      orientations[num_frames++] = q;
    } else {
      hs::log("Orientation full, dropping frame!");
    }
    return *this;
  }

  /**
   * @brief Collapses the orientation history, retaining only the latest
   * quaternion.
   * @details Used after rendering a motion step to reset the motion blur
   * history.
   * @return Reference to the Orientation object.
   */
  Orientation &collapse() {
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
  Quaternion &at(int i) { return orientations[i]; }

  /**
   * @brief Increases the resolution of the history to 'count' steps, preserving
   * shape via Slerp.
   * @param count The target number of steps in the history.
   */
  void upsample(int count) {
    if (num_frames >= count)
      return;
    if (count > CAPACITY)
      count = CAPACITY;

    std::array<Quaternion, CAPACITY> old_orientations;
    std::copy(orientations.begin(), orientations.begin() + num_frames,
              old_orientations.begin());

    int old_num_frames = num_frames;

    for (int i = 0; i < count; ++i) {
      float t = static_cast<float>(i) / (count - 1);
      float source_float_index = t * (old_num_frames - 1);
      int idx = static_cast<int>(source_float_index);
      float frac = source_float_index - idx;

      orientations[i] = slerp(
          old_orientations[idx],
          old_orientations[std::min((int)old_num_frames - 1, idx + 1)], frac);
    }
    num_frames = count;
  }

private:
  std::array<Quaternion, CAPACITY>
      orientations; /**< Storage for historical quaternions. */
  int num_frames;   /**< The current number of active frames in history. */
};

/**
 * @brief Generates a truly random 3D unit vector (direction) using Marsaglia's
 * method.
 * @return A normalized random Vector.
 */
inline Vector random_vector() {
  // Marsaglia's method
  float v1, v2, s;
  do {
    v1 = 2.0f * hs::rand_f() - 1.0f;
    v2 = 2.0f * hs::rand_f() - 1.0f;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1.0f || s == 0.0f);

  float sqrt_s = sqrtf(1.0f - s);
  return Vector(2.0f * v1 * sqrt_s, 2.0f * v2 * sqrt_s, 1.0f - 2.0f * s);
}

/**
 * @brief Parameters defining a Lissajous curve.
 */
struct LissajousParams {
  float m1; /**< Frequency coefficient for the axial components (X and Z). */
  float m2; /**< Frequency coefficient for the orbital component (Y). */
  float a;  /**< Phase shift multiplier (multiplies PI_F). */
  float domain; /**< The total duration (t) over which the curve is drawn. */
};

/**
 * @brief Calculates a 3D point on the unit sphere corresponding to a spherical
 * Lissajous curve.
 * @param m1 Frequency coefficient for XZ plane.
 * @param m2 Frequency coefficient for Y axis.
 * @param a Phase shift multiplier.
 * @param t Time variable (or position along the domain).
 * @return The calculated 3D point (unit vector).
 */
inline Vector lissajous(float m1, float m2, float a, float t) {
  Vector v(sinf(m2 * t) * cosf(m1 * t - a * PI_F), cosf(m2 * t),
           sinf(m2 * t) * sinf(m1 * t - a * PI_F));
  return v.normalized();
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
inline Basis make_basis(const Quaternion &orientation, const Vector &normal) {
  Vector ref_axis =
      (std::abs(dot(normal, X_AXIS)) > (1 - TOLERANCE)) ? Y_AXIS : X_AXIS;
  Vector v = rotate(normal, orientation).normalized();
  Vector ref = rotate(ref_axis, orientation).normalized();
  Vector u = cross(v, ref).normalized();
  Vector w = cross(v, u).normalized();
  return {u, v, w};
}

/**
 * @brief Adjusted basis and radius for drawing on the opposite side of the
 * sphere.
 * @param basis The current basis {u, v, w}.
 * @param radius Angular radius (0-2).
 * @return A pair containing the adjusted Basis and radius.
 */
inline std::pair<Basis, float> get_antipode(const Basis &basis, float radius) {
  if (radius > 1.0f) {
    Basis new_basis;
    new_basis.u = -basis.u; // Flip U to maintain chirality
    new_basis.v = -basis.v; // Flip V (Antipode)
    new_basis.w = basis.w;  // W stays (Rotation axis)
    return {new_basis, 2.0f - radius};
  }
  return {basis, radius};
}