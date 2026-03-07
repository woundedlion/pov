/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <cfloat>
#include <limits>
#include <cassert>
#include <algorithm>
#include "platform.h"
#include "FastNoiseLite.h"

/**
 * @brief The Golden Ratio constant (Phi).
 */
static constexpr float PHI = 1.6180339887498948482045868343656f;
/**
 * @brief The inverse of the Golden Ratio (1/Phi).
 */
static constexpr float G = 1 / PHI;
/**
 * @brief Tolerance used for floating-point comparisons.
 */
static constexpr float TOLERANCE = 0.0001f;
/**
 * @brief Floating-point representation of PI.
 */
static constexpr float PI_F = static_cast<float>(PI);

/**
 * @brief Quintic interpolation kernel (smootherstep).
 * @param t Interpolation factor (clamped to 0.0 - 1.0).
 * @return The interpolated value: t^3 * (t * (t * 6 - 15) + 10).
 */
inline float quintic_kernel(float t) {
  t = hs::clamp(t, 0.0f, 1.0f);
  return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

struct Vector;

/**
 * @brief Represents a point in Spherical Coordinates.
 */
struct Spherical {
  /**
   * @brief Constructs a Spherical coordinate.
   * @param theta Azimuthal angle.
   * @param phi Polar angle.
   */
  Spherical(float theta, float phi) : theta(theta), phi(phi) {}
  /**
   * @brief Constructs a Spherical coordinate from a 3D Vector.
   * @param v The Cartesian vector.
   */
  Spherical(const Vector &v);

  float theta; /**< Azimuthal angle (usually longitude). */
  float phi;   /**< Polar angle (usually co-latitude). */
};

/**
 * @brief Represents a 3D Vector in Cartesian coordinates (x, y, z).
 */
struct Vector {
  /**
   * @brief Default constructor (initializes to zero vector).
   */
  constexpr Vector() {}
  /**
   * @brief Constructs a vector with explicit components.
   * @param x X-component.
   * @param y Y-component.
   * @param z Z-component.
   */
  constexpr Vector(float x, float y, float z) : x(x), y(y), z(z) {}
  /**
   * @brief Copy constructor.
   * @param v The vector to copy.
   */
  constexpr Vector(const Vector &v) : x(v.x), y(v.y), z(v.z) {}
  /**
   * @brief Constructs a vector from Spherical coordinates.
   * @param s The spherical coordinate.
   */
  Vector(const Spherical &s)
      : x(sinf(s.phi) * cosf(s.theta)), y(cosf(s.phi)),
        z(sinf(s.phi) * sinf(s.theta)) {}
  /**
   * @brief Constructs a vector from theta and phi angles.
   * @param theta Azimuthal angle.
   * @param phi Polar angle.
   */
  Vector(float theta, float phi) : Vector(Spherical(theta, phi)) {}

  /**
   * @brief Assignment operator.
   * @param v The vector to assign from.
   * @return Reference to the current object.
   */
  constexpr Vector &operator=(const Vector &v) {
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
  }

  /**
   * @brief Equality comparison using tolerance.
   * @param v The vector to compare.
   * @return True if components are within TOLERANCE.
   */
  bool operator==(const Vector &v) const {
    return std::abs(x - v.x) <= TOLERANCE && std::abs(y - v.y) <= TOLERANCE &&
           std::abs(z - v.z) <= TOLERANCE;
  }

  /**
   * @brief Inequality comparison.
   * @param v The vector to compare.
   * @return True if components are not equal within TOLERANCE.
   */
  bool operator!=(const Vector &v) const { return !(*this == v); }

  /**
   * @brief Unary negation operator.
   * @return A new vector with negated components.
   */
  constexpr Vector operator-() const { return Vector(-x, -y, -z); }

  Vector &operator+=(const Vector &v_other) {
    x += v_other.x;
    y += v_other.y;
    z += v_other.z;
    return *this;
  }
  Vector &operator-=(const Vector &v_other) {
    x -= v_other.x;
    y -= v_other.y;
    z -= v_other.z;
    return *this;
  }
  Vector &operator*=(float s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  Vector &operator/=(float s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  /**
   * @brief Calculates the magnitude (length) of the vector.
   * @return The magnitude.
   */
  constexpr float length() const { return sqrtf(x * x + y * y + z * z); }

  /**
   * @brief Alias for length().
   */
  constexpr float magnitude() const { return length(); }

  /**
   * @brief Normalizes the vector (scales to unit length).
   * @return Reference to the normalized vector.
   */
  Vector &normalize() {
    float m = length();
    if (m < std::numeric_limits<float>::epsilon()) {
      hs::log("Can't normalize a zero vector!");
      x = 1;
      y = z = 0;
    } else {
      x = x / m;
      y = y / m;
      z = z / m;
    }
    return *this;
  }

  float x = 0; /**< X-component. */
  float y = 0; /**< Y-component. */
  float z = 0; /**< Z-component. */
};

/**
 * @brief Constructs a Spherical coordinate from a 3D Vector.
 * @param v The Cartesian vector.
 */
inline Spherical::Spherical(const Vector &v) {
  Vector n(v);
  n.normalize();
  theta = atan2f(n.z, n.x);
  phi = acosf(hs::clamp(n.y, -1.0f, 1.0f));
}

struct Quaternion;
constexpr Quaternion operator*(const Quaternion &q1, const Quaternion &q2);
constexpr Quaternion operator/(const Quaternion &q, float s);
/**
 * @brief Scalar division decl.
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
constexpr Vector operator/(const Vector &v, float s);

/**
 * @brief Represents a Quaternion for 3D rotation (r + i*x + j*y + k*z).
 */
struct Quaternion {
  /**
   * @brief Default constructor (identity quaternion).
   */
  constexpr Quaternion() {}
  /**
   * @brief Copy constructor.
   * @param q The quaternion to copy.
   */
  constexpr Quaternion(const Quaternion &q) : r(q.r), v(q.v) {}
  /**
   * @brief Constructs a quaternion with explicit components.
   * @param a Real component (r).
   * @param b X imaginary component (i).
   * @param c Y imaginary component (j).
   * @param d Z imaginary component (k).
   */
  constexpr Quaternion(float a, float b, float c, float d) : r(a), v(b, c, d) {}
  /**
   * @brief Constructs a quaternion from a scalar and a Vector part.
   * @param a Real component (r).
   * @param v Vector part (x, y, z).
   */
  constexpr Quaternion(float a, const Vector &v) : r(a), v(v) {}

  /**
   * @brief Assignment operator.
   * @param q The quaternion to assign from.
   * @return Reference to the current object.
   */
  Quaternion &operator=(const Quaternion &q) {
    r = q.r;
    v = q.v;
    return *this;
  }

  /**
   * @brief Equality comparison using tolerance.
   * @param q The quaternion to compare.
   * @return True if components are within TOLERANCE.
   */
  bool operator==(const Quaternion &q) const {
    return std::abs(q.r - r) < TOLERANCE && q.v == v;
  }

  /**
   * @brief Quaternion multiplication compound assignment.
   * @param q The quaternion to multiply by.
   */
  Quaternion &operator*=(const Quaternion &q) {
    *this = *this * q;
    return *this;
  }
  Quaternion &operator+=(const Quaternion &q) {
    r += q.r;
    v += q.v;
    return *this;
  }
  Quaternion &operator-=(const Quaternion &q) {
    r -= q.r;
    v -= q.v;
    return *this;
  }
  Quaternion &operator*=(float s) {
    r *= s;
    v *= s;
    return *this;
  }

  /**
   * @brief Calculates the squared magnitude of the quaternion.
   * @return The squared magnitude.
   */
  constexpr float squared_magnitude() const {
    return r * r + v.x * v.x + v.y * v.y + v.z * v.z;
  }

  /**
   * @brief Calculates the inverse of the quaternion.
   * @return The inverse quaternion.
   */
  Quaternion inverse() const {
    float sq_mag = squared_magnitude();
    return Quaternion(r, -v) / sq_mag;
  }

  /**
   * @brief Calculates the conjugate of the quaternion.
   * @return The conjugate quaternion (r, -v).
   */
  Quaternion conjugate() const { return Quaternion(r, -v); }

  /**
   * @brief Unary negation operator.
   * @return A new quaternion with negated components.
   */
  constexpr Quaternion operator-() const { return Quaternion(-r, -v); }

  /**
   * @brief Calculates the magnitude (length) of the quaternion.
   * @return The magnitude.
   */
  constexpr float magnitude() const {
    return sqrtf(r * r + v.x * v.x + v.y * v.y + v.z * v.z);
  }

  /**
   * @brief Normalizes the quaternion (scales to unit magnitude).
   * @return Reference to the normalized quaternion.
   */
  Quaternion &normalize() {
    auto m = magnitude();
    if (m <= std::numeric_limits<float>::epsilon()) {
      hs::log("Can't normalize a zero Quaternion!");
      r = 1;
      v = Vector(0, 0, 0);
    } else {
      r = r / m;
      v = v / m;
    }
    return *this;
  }

  float r = 1; /**< Real component. */
  Vector v;    /**< Vector part (x, y, z). */
};

/**
 * @brief Represents a Complex number.
 */
struct Complex {
  float re; /**< Real component. */
  float im; /**< Imaginary component. */

  /**
   * @brief Default constructor.
   */
  constexpr Complex() : re(0), im(0) {}
  /**
   * @brief Constructs a complex number.
   * @param r Real part.
   * @param i Imaginary part.
   */
  constexpr Complex(float r, float i) : re(r), im(i) {}

  /**
   * @brief Complex addition.
   */
  constexpr Complex operator+(const Complex &b) const {
    return Complex(re + b.re, im + b.im);
  }
  /**
   * @brief Complex subtraction.
   */
  constexpr Complex operator-(const Complex &b) const {
    return Complex(re - b.re, im - b.im);
  }
  /**
   * @brief Complex multiplication.
   */
  constexpr Complex operator*(const Complex &b) const {
    return Complex(re * b.re - im * b.im, re * b.im + im * b.re);
  }

  /**
   * @brief Complex division.
   */
  Complex operator/(const Complex &b) const {
    float denom = b.re * b.re + b.im * b.im;
    if (std::abs(denom) < 0.000001f)
      return Complex(0, 0);
    return Complex((re * b.re + im * b.im) / denom,
                   (im * b.re - re * b.im) / denom);
  }
};

inline Vector fold_to_hemisphere(const Vector &v) {
  // Half-angle identities to calculate the new height (j') and the new
  // horizontal radius (r')
  float j_new = sqrtf((1.0f + v.y) * 0.5f);
  float r_new = sqrtf((1.0f - v.y) * 0.5f);

  // Get the length of the original horizontal components (x and z)
  float r_old = sqrtf(v.x * v.x + v.z * v.z);

  // Protect against division by zero at the exact poles
  if (r_old == 0.0f) {
    // The North Pole (j=1) stays the North Pole.
    // The South Pole (j=-1) is stretched into the entire equator ring.
    // We arbitrarily anchor the singularity to the X-axis (i=1).
    return v.y > 0.0f ? Vector(0.0f, 1.0f, 0.0f) : Vector(1.0f, 0.0f, 0.0f);
  }

  // Scale the original x and z coordinates to the new tightened radius
  float scale = r_new / r_old;

  return Vector(v.x * scale, j_new, v.z * scale);
}

/**
 * @brief Class to hold Mobius parameters with mutable components.
 */
struct MobiusParams {
  float aRe, /**< Real part of a. */
      aIm,   /**< Imaginary part of a. */
      bRe,   /**< Real part of b. */
      bIm,   /**< Imaginary part of b. */
      cRe,   /**< Real part of c. */
      cIm,   /**< Imaginary part of c. */
      dRe,   /**< Real part of d. */
      dIm;   /**< Imaginary part of d. */

  constexpr MobiusParams()
      : aRe(1), aIm(0), bRe(0), bIm(0), cRe(0), cIm(0), dRe(1), dIm(0) {}
  constexpr MobiusParams(float ar, float ai, float br, float bi, float cr,
                         float ci, float dr, float di)
      : aRe(ar), aIm(ai), bRe(br), bIm(bi), cRe(cr), cIm(ci), dRe(dr), dIm(di) {
  }

  Complex getA() const { return Complex(aRe, aIm); }
  Complex getB() const { return Complex(bRe, bIm); }
  Complex getC() const { return Complex(cRe, cIm); }
  Complex getD() const { return Complex(dRe, dIm); }
};

/**
 * @brief Stereographic Projection: Sphere -> Complex Plane.
 */
inline Complex stereo(const Vector &v) {
  float denom = 1.0f - v.y;
  if (std::abs(denom) < 0.0001f)
    return Complex(100, 100); // Infinity
  return Complex(v.x / denom, v.z / denom);
}

/**
 * @brief Inverse Stereographic Projection: Complex Plane -> Sphere.
 */
inline Vector inv_stereo(const Complex &z) {
  float r2 = z.re * z.re + z.im * z.im;
  return Vector(2 * z.re / (r2 + 1), (r2 - 1) / (r2 + 1), 2 * z.im / (r2 + 1));
}

/**
 * @brief Mobius Transformation: f(z) = (az + b) / (cz + d).
 */
inline Complex mobius(const Complex &z, const MobiusParams &params) {
  Complex num = (params.getA() * z) + params.getB();
  Complex den = (params.getC() * z) + params.getD();
  return num / den;
}

/**
 * @brief Gnomonic Projection: Sphere -> Plane (Equator at Infinity).
 * Projects from center (0,0,0) to plane z=1 (tangent at North Pole, i.e., k=1).
 */
inline Complex gnomonic(const Vector &v) {
  // Handle equator singularity with a large number instead of Infinity
  // Projects to plane at j=1.
  float div = (std::abs(v.y) < 1e-9f) ? 1e-9f * (v.y >= 0 ? 1.0f : -1.0f) : v.y;
  return Complex(v.x / div, v.z / div);
}

/**
 * @brief Inverse Gnomonic: Plane -> Sphere.
 * @param z Complex point on the plane.
 * @param original_sign The sign of the z-component (k) of the original vector.
 */
inline Vector inv_gnomonic(const Complex &z, float original_sign = 1.0f) {
  // Project (re, 1, im) back onto unit sphere
  float len = sqrtf(z.re * z.re + z.im * z.im + 1.0f);
  float inv_len = 1.0f / len;

  // Restore hemisphere sign (Upper or Lower)
  return Vector(z.re * inv_len * original_sign, // i
                inv_len * original_sign,        // j
                z.im * inv_len * original_sign  // k
  );
}

/**
 * @brief Rotates a vector by a unit quaternion.
 * @param v The vector to rotate.
 * @param q The unit rotation quaternion.
 * @return The rotated vector.
 */
inline Vector rotate(const Vector &v, const Quaternion &q) {
  Quaternion p(0, v);
  auto r = q * p * q.inverse();
  return r.v;
}

/**
 * @brief Vector addition.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The resulting vector.
 */
constexpr Vector operator+(const Vector &v1, const Vector &v2) {
  return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

/**
 * @brief Vector subtraction.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The resulting vector.
 */
constexpr Vector operator-(const Vector &v1, const Vector &v2) {
  return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

/**
 * @brief Scalar multiplication (vector * scalar).
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
constexpr Vector operator*(const Vector &v, float s) {
  return Vector(s * v.x, s * v.y, s * v.z);
}

/**
 * @brief Scalar multiplication (scalar * vector).
 * @param s Scalar.
 * @param v Vector.
 * @return The resulting vector.
 */
constexpr Vector operator*(float s, const Vector &v) { return v * s; }

/**
 * @brief Scalar division.
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
constexpr Vector operator/(const Vector &v, float s) {
  return Vector(v.x / s, v.y / s, v.z / s);
}

/**
 * @brief Calculates the dot product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The dot product (scalar).
 */
constexpr float dot(const Vector &v1, const Vector &v2) {
  return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

/**
 * @brief Calculates the cross product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The cross product vector.
 */
constexpr Vector cross(const Vector &v1, const Vector &v2) {
  return Vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                v1.x * v2.y - v1.y * v2.x);
}

/**
 * @brief Calculates the Euclidean distance between two vectors.
 * @param a First vector.
 * @param b Second vector.
 * @return The distance (scalar).
 */
constexpr float distance_between(const Vector &a, const Vector &b) {
  float dx = b.x - a.x;
  float dy = b.y - a.y;
  float dz = b.z - a.z;
  return sqrtf(dx * dx + dy * dy + dz * dz);
}

/**
 * @brief Calculates the square of the Euclidean distance between two vectors.
 * Useful for efficient comparisons without square root.
 * @param a First vector.
 * @param b Second vector.
 * @return The distance squared (scalar).
 */
constexpr float distance_squared(const Vector &a, const Vector &b) {
  float dx = b.x - a.x;
  float dy = b.y - a.y;
  float dz = b.z - a.z;
  return dx * dx + dy * dy + dz * dz;
}

/**
 * @brief Calculates the angle between two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The angle in radians.
 */
constexpr float angle_between(const Vector &v1, const Vector &v2) {
  float len_product = v1.length() * v2.length();
  if (len_product <= std::numeric_limits<float>::epsilon()) {
    hs::log("Cannot calculate angle between 0-vectors!");
    return 0;
  }
  float d = dot(v1, v2) / len_product;
  return acosf(hs::clamp(d, -1.0f, 1.0f));
}

/**
 * @brief Quaternion addition.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator+(const Quaternion &q1, const Quaternion &q2) {
  return Quaternion(q1.r + q2.r, q1.v + q2.v);
}

/**
 * @brief Quaternion subtraction.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator-(const Quaternion &q1, const Quaternion &q2) {
  return Quaternion(q1.r - q2.r, q1.v - q2.v);
}

/**
 * @brief Quaternion multiplication (Hamilton product).
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(const Quaternion &q1, const Quaternion &q2) {
  return Quaternion(q1.r * q2.r - dot(q1.v, q2.v),
                    q1.r * q2.v + q2.r * q1.v + cross(q1.v, q2.v));
}

/**
 * @brief Scalar multiplication (quaternion * scalar).
 * @param q Quaternion.
 * @param s Scalar.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(const Quaternion &q, float s) {
  return Quaternion(q.r * s, q.v * s);
}

/**
 * @brief Scalar multiplication (scalar * quaternion).
 * @param s Scalar.
 * @param q Quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(float s, const Quaternion &q) { return q * s; }

/**
 * @brief Scalar division.
 * @param q Quaternion.
 * @param s Scalar.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator/(const Quaternion &q, float s) {
  return Quaternion(q.r / s, q.v / s);
}

/**
 * @brief Calculates the dot product of two quaternions.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The dot product (scalar).
 */
constexpr float dot(const Quaternion &q1, const Quaternion &q2) {
  return (q1.r * q2.r) + (q1.v.x * q2.v.x) + (q1.v.y * q2.v.y) +
         (q1.v.z * q2.v.z);
}

/**
 * @brief Calculates the angle between two unit quaternions.
 * @param q1 First unit quaternion.
 * @param q2 Second unit quaternion.
 * @return The angle in radians.
 */
constexpr float angle_between(const Quaternion &q1, const Quaternion &q2) {
  return acosf(hs::clamp(dot(q1, q2), -1.0f, 1.0f));
}

/**
 * @brief Creates a rotation quaternion from an axis and an angle.
 * @param axis The rotation axis (must be unit vector).
 * @param theta The angle of rotation in radians.
 * @return The resulting unit rotation quaternion.
 */
inline Quaternion make_rotation(const Vector &axis, float theta) {
  return Quaternion(cosf(theta / 2), sinf(theta / 2) * axis).normalize();
}

/**
 * @brief Creates a rotation quaternion to rotate from one vector to another.
 * @param from The source vector (must be unit vector).
 * @param to The destination vector (must be unit vector).
 * @return The resulting unit rotation quaternion.
 */
inline Quaternion make_rotation(const Vector &from, const Vector &to) {
  float d = dot(from, to);

  // Handle antiparallel vectors (180 degrees apart)
  if (d < -1.0f + TOLERANCE) {
    // Choose an axis perpendicular to 'from'
    Vector axis = cross(Vector(1, 0, 0), from);
    if (axis.length() < TOLERANCE) {
      axis = cross(Vector(0, 1, 0), from);
    }
    axis.normalize();
    return make_rotation(axis, PI_F);
  }

  // Handle parallel vectors (0 degrees apart)
  if (d > 1.0f - TOLERANCE) {
    return Quaternion(1, 0, 0, 0); // Identity
  }

  auto axis = cross(from, to).normalize();
  auto angle = angle_between(from, to);
  return make_rotation(axis, angle);
}

/**
 * @brief Spherical Linear Interpolation (SLERP) between two Vectors.
 * @param v1 Starting vector.
 * @param v2 Ending vector.
 * @param t Interpolation factor (0.0 to 1.0).
 * @return The interpolated unit vector.
 */
inline Vector slerp(const Vector &v1, const Vector &v2, float t) {
  float d = hs::clamp(dot(v1, v2), -1.0f, 1.0f);
  // If vectors are extremely close, just lerp to avoid NaN
  if (d > 0.9999f) {
    return (v1 + (v2 - v1) * t).normalize();
  }
  float theta = acosf(d);
  float sin_theta = sinf(theta);
  if (sin_theta < 0.0001f) {
    return (v1 + (v2 - v1) * t).normalize();
  }
  float s1 = sinf((1 - t) * theta) / sin_theta;
  float s2 = sinf(t * theta) / sin_theta;
  return ((s1 * v1) + (s2 * v2)).normalize();
}

/**
 * @brief Spherical Linear Interpolation (SLERP) between two quaternions.
 * @param q1 Starting quaternion.
 * @param q2 Ending quaternion.
 * @param t Interpolation factor (0.0 to 1.0).
 * @param long_way If true, takes the longest path between quaternions.
 * @return The interpolated unit quaternion.
 */
inline Quaternion slerp(const Quaternion &q1, const Quaternion &q2, float t,
                        bool long_way = false) {
  float d = dot(q1, q2);
  Quaternion p(q1);
  Quaternion q(q2);

  if ((long_way && d > 0) || (!long_way && d < 0)) {
    p = -p;
    d = -d;
  }

  if (d > (1 - TOLERANCE)) {
    Quaternion r = p + t * (q - p);
    return r.normalize();
  }

  float theta = acosf(hs::clamp(d, -1.0f, 1.0f));
  float sin_theta = sinf(theta);
  if (sin_theta < 0.0001f) {
    Quaternion r = p + t * (q - p);
    return r.normalize();
  }
  float s1 = sinf((1 - t) * theta) / sin_theta;
  float s2 = sinf(t * theta) / sin_theta;
  return ((s1 * p) + (s2 * q)).normalize();
}
