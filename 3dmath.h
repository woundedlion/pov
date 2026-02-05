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

#include <vector>
#include <array>

/**
 * @brief A lightweight view over a contiguous sequence of objects.
 * @tparam T The type of object (can be const).
 */
template <typename T>
struct Span {
  T* ptr;
  size_t len;

  Span() : ptr(nullptr), len(0) {}
  Span(T* p, size_t l) : ptr(p), len(l) {}
  
  // Convert from vector
  template <typename U>
  Span(const std::vector<U>& v) : ptr(v.data()), len(v.size()) {}
  
  // Convert from array
  template <typename U, size_t N>
  Span(const std::array<U, N>& arr) : ptr(arr.data()), len(N) {}

  T& operator[](size_t i) const { return ptr[i]; }
  size_t size() const { return len; }
  T* begin() const { return ptr; }
  T* end() const { return ptr + len; }
  bool empty() const { return len == 0; }
};

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
  Spherical(const Vector& v);

  float theta; /**< Azimuthal angle (usually longitude). */
  float phi; /**< Polar angle (usually co-latitude). */
};

/**
 * @brief Represents a 3D Vector in Cartesian coordinates (i, j, k).
 */
struct Vector {
  /**
   * @brief Default constructor (initializes to zero vector).
   */
  constexpr Vector() {}
  /**
   * @brief Constructs a vector with explicit components.
   * @param i X-component.
   * @param j Y-component.
   * @param k Z-component.
   */
  constexpr Vector(float i, float j, float k) : i(i), j(j), k(k) {}
  /**
   * @brief Copy constructor.
   * @param v The vector to copy.
   */
  constexpr Vector(const Vector& v) : i(v.i), j(v.j), k(v.k) {}
  /**
   * @brief Constructs a vector from Spherical coordinates.
   * @param s The spherical coordinate.
   */
  Vector(const Spherical& s) :
    i(sinf(s.phi)* cosf(s.theta)),
    j(sinf(s.phi)* sinf(s.theta)),
    k(cosf(s.phi))
  {
  }
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
  constexpr Vector& operator=(const Vector& v) {
    i = v.i;
    j = v.j;
    k = v.k;
    return *this;
  }

  /**
   * @brief Equality comparison using tolerance.
   * @param v The vector to compare.
   * @return True if components are within TOLERANCE.
   */
  bool operator==(const Vector& v) const {
    return std::abs(i - v.i) <= TOLERANCE
      && std::abs(j - v.j) <= TOLERANCE
      && std::abs(k - v.k) <= TOLERANCE;
  }

  /**
   * @brief Inequality comparison.
   * @param v The vector to compare.
   * @return True if components are not equal within TOLERANCE.
   */
  bool operator!=(const Vector& v) const {
    return !(*this == v);
  }

  /**
   * @brief Unary negation operator.
   * @return A new vector with negated components.
   */
  constexpr Vector operator-() const {
    return Vector(-i, -j, -k);
  }

  Vector& operator+=(const Vector& v_other) {
    i += v_other.i; j += v_other.j; k += v_other.k;
    return *this;
  }
  Vector& operator-=(const Vector& v_other) {
    i -= v_other.i; j -= v_other.j; k -= v_other.k;
    return *this;
  }
  Vector& operator*=(float s) {
    i *= s; j *= s; k *= s;
    return *this;
  }
  Vector& operator/=(float s) {
    i /= s; j /= s; k /= s;
    return *this;
  }

  /**
   * @brief Calculates the magnitude (length) of the vector.
   * @return The magnitude.
   */
  constexpr float length() const { return sqrtf(i * i + j * j + k * k); }
  
  /**
   * @brief Alias for length().
   */
  constexpr float magnitude() const { return length(); }

  /**
   * @brief Normalizes the vector (scales to unit length).
   * @return Reference to the normalized vector.
   */
  Vector& normalize() {
    float m = length();
    if (m < std::numeric_limits<float>::epsilon()) {
      Serial.println("Can't normalize a zero vector!");
      i = 1;
      j = k = 0;
    }
    else {
      i = i / m;
      j = j / m;
      k = k / m;
    }
    return *this;
  }

  float i = 0; /**< X-component. */
  float j = 0; /**< Y-component. */
  float k = 0; /**< Z-component. */
};

/**
 * @brief Constructs a Spherical coordinate from a 3D Vector.
 * @param v The Cartesian vector.
 */
Spherical::Spherical(const Vector& v)
{
  Vector n(v);
  n.normalize();
  theta = atan2f(n.j, n.i);
  phi = acosf(std::clamp(n.k, -1.0f, 1.0f));
}

struct Quaternion;
constexpr Quaternion operator*(const Quaternion& q1, const Quaternion& q2);
constexpr Quaternion operator/(const Quaternion& q, float s);
constexpr Vector operator/(const Vector& v, float s);

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
  constexpr Quaternion(const Quaternion& q) : r(q.r), v(q.v) {}
  /**
   * @brief Constructs a quaternion with explicit components.
   * @param a Real component (r).
   * @param b X imaginary component (i).
   * @param c Y imaginary component (j).
   * @param d Z imaginary component (k).
   */
  constexpr  Quaternion(float a, float b, float c, float d) : r(a), v(b, c, d) {}
  /**
   * @brief Constructs a quaternion from a scalar and a Vector part.
   * @param a Real component (r).
   * @param v Vector part (i, j, k).
   */
  constexpr Quaternion(float a, const Vector& v) : r(a), v(v) {}

  /**
   * @brief Assignment operator.
   * @param q The quaternion to assign from.
   * @return Reference to the current object.
   */
  Quaternion& operator=(const Quaternion& q) {
    r = q.r;
    v = q.v;
    return *this;
  }

  /**
   * @brief Equality comparison using tolerance.
   * @param q The quaternion to compare.
   * @return True if components are within TOLERANCE.
   */
  bool operator==(const Quaternion& q) const {
    return std::abs(q.r - r) < TOLERANCE
      && q.v == v;
  }

  /**
   * @brief Quaternion multiplication compound assignment.
   * @param q The quaternion to multiply by.
   */
  Quaternion& operator*=(const Quaternion& q) {
    *this = *this * q;
    return *this;
  }
  Quaternion& operator+=(const Quaternion& q) {
    r += q.r; v += q.v;
    return *this;
  }
  Quaternion& operator-=(const Quaternion& q) {
    r -= q.r; v -= q.v;
    return *this;
  }
  Quaternion& operator*=(float s) {
    r *= s; v *= s;
    return *this;
  }

  /**
   * @brief Calculates the squared magnitude of the quaternion.
   * @return The squared magnitude.
   */
  constexpr float squared_magnitude() const {
    return r * r + v.i * v.i + v.j * v.j + v.k * v.k;
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
  Quaternion conjugate() const {
      return Quaternion(r, -v);
  }

  /**
   * @brief Unary negation operator.
   * @return A new quaternion with negated components.
   */
  constexpr Quaternion operator-() const {
    return Quaternion(-r, -v);
  }

  /**
   * @brief Calculates the magnitude (length) of the quaternion.
   * @return The magnitude.
   */
  constexpr float magnitude() const {
    return sqrtf(r * r + v.i * v.i + v.j * v.j + v.k * v.k);
  }

  /**
   * @brief Normalizes the quaternion (scales to unit magnitude).
   * @return Reference to the normalized quaternion.
   */
  Quaternion& normalize() {
    auto m = magnitude();
    if (m <= std::numeric_limits<float>::epsilon()) {
      Serial.println("Can't normalize a zaro Quaternion!");
      r = 1;
      v = Vector(0, 0, 0);
    }
    else {
      r = r / m;
      v = v / m;
    }
    return *this;
  }

  float r = 1; /**< Real component. */
  Vector v; /**< Vector part (i, j, k). */
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
  constexpr Complex operator+(const Complex& b) const { return Complex(re + b.re, im + b.im); }
  /**
   * @brief Complex subtraction.
   */
  constexpr Complex operator-(const Complex& b) const { return Complex(re - b.re, im - b.im); }
  /**
   * @brief Complex multiplication.
   */
  constexpr Complex operator*(const Complex& b) const { return Complex(re * b.re - im * b.im, re * b.im + im * b.re); }

  /**
   * @brief Complex division.
   */
  Complex operator/(const Complex& b) const {
    float denom = b.re * b.re + b.im * b.im;
    if (std::abs(denom) < 0.000001f) return Complex(0, 0);
    return Complex((re * b.re + im * b.im) / denom, (im * b.re - re * b.im) / denom);
  }
};

/**
 * @brief Class to hold Mobius parameters with mutable components.
 */
struct MobiusParams {
  float aRe, aIm, bRe, bIm, cRe, cIm, dRe, dIm;

  constexpr MobiusParams() : aRe(1), aIm(0), bRe(0), bIm(0), cRe(0), cIm(0), dRe(1), dIm(0) {}
  constexpr MobiusParams(float ar, float ai, float br, float bi, float cr, float ci, float dr, float di)
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
Complex stereo(const Vector& v) {
  float denom = 1.0f - v.k;
  if (std::abs(denom) < 0.0001f) return Complex(100, 100); // Infinity
  return Complex(v.i / denom, v.j / denom);
}

/**
 * @brief Inverse Stereographic Projection: Complex Plane -> Sphere.
 */
Vector inv_stereo(const Complex& z) {
  float r2 = z.re * z.re + z.im * z.im;
  return Vector(
    2 * z.re / (r2 + 1),
    2 * z.im / (r2 + 1),
    (r2 - 1) / (r2 + 1)
  );
}

/**
 * @brief Mobius Transformation: f(z) = (az + b) / (cz + d).
 */
Complex mobius(const Complex& z, const MobiusParams& params) {
  Complex num = (params.getA() * z) + params.getB();
  Complex den = (params.getC() * z) + params.getD();
  return num / den;
}

Vector mobius_transform(const Vector& v, const MobiusParams& params) {
    return inv_stereo(mobius(stereo(v), params));
}

/**
 * @brief Rotates a vector by a unit quaternion.
 * @param v The vector to rotate.
 * @param q The unit rotation quaternion.
 * @return The rotated vector.
 */
Vector rotate(const Vector& v, const Quaternion& q) {
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
constexpr Vector operator+(const Vector& v1, const Vector& v2) {
  return Vector(v1.i + v2.i, v1.j + v2.j, v1.k + v2.k);
}

/**
 * @brief Vector subtraction.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The resulting vector.
 */
constexpr Vector operator-(const Vector& v1, const Vector& v2) {
  return Vector(v1.i - v2.i, v1.j - v2.j, v1.k - v2.k);
}

/**
 * @brief Scalar multiplication (vector * scalar).
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
constexpr Vector operator*(const Vector& v, float s) {
  return Vector(s * v.i, s * v.j, s * v.k);
}

/**
 * @brief Scalar multiplication (scalar * vector).
 * @param s Scalar.
 * @param v Vector.
 * @return The resulting vector.
 */
constexpr Vector operator*(float s, const Vector& v) {
  return v * s;
}

/**
 * @brief Scalar division.
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
constexpr Vector operator/(const Vector& v, float s) {
  return Vector(v.i / s, v.j / s, v.k / s);
}

/**
 * @brief Calculates the dot product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The dot product (scalar).
 */
constexpr float dot(const Vector& v1, const Vector& v2) {
  return (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
}

/**
 * @brief Calculates the cross product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The cross product vector.
 */
constexpr Vector cross(const Vector& v1, const Vector& v2) {
  return Vector(
    v1.j * v2.k - v1.k * v2.j,
    v1.k * v2.i - v1.i * v2.k,
    v1.i * v2.j - v1.j * v2.i);
}

/**
 * @brief Calculates the Euclidean distance between two vectors.
 * @param a First vector.
 * @param b Second vector.
 * @return The distance (scalar).
 */
constexpr float distance_between(const Vector& a, const Vector& b) {
  float di = b.i - a.i;
  float dj = b.j - a.j;
  float dk = b.k - a.k;
  return sqrtf(di * di + dj * dj + dk * dk);
}

/**
 * @brief Calculates the square of the Euclidean distance between two vectors.
 * Usefule  for efficient comparisons without square root.
 * @param a First vector.
 * @param b Second vector.
 * @return The distance squared (scalar).
 */
constexpr float distance_squared(const Vector& a, const Vector& b) {
  float di = b.i - a.i;
  float dj = b.j - a.j;
  float dk = b.k - a.k;
  return di * di + dj * dj + dk * dk;
}

/**
 * @brief Calculates the angle between two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The angle in radians.
 */
constexpr float angle_between(const Vector& v1, const Vector& v2) {
  float len_product = v1.length() * v2.length();
  if (len_product <= std::numeric_limits<float>::epsilon()) {
    Serial.println("Cannot calculate angle between 0-vectors!");
    return 0;
  }
  float d = dot(v1, v2) / len_product;
  return acosf(std::clamp(d, -1.0f, 1.0f));
}

/**
 * @brief Quaternion addition.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(q1.r + q2.r, q1.v + q2.v);
}

/**
 * @brief Quaternion subtraction.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(q1.r - q2.r, q1.v - q2.v);
}

/**
 * @brief Quaternion multiplication (Hamilton product).
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(const Quaternion& q1, const Quaternion& q2) {
  return Quaternion(
    q1.r * q2.r - dot(q1.v, q2.v),
    q1.r * q2.v + q2.r * q1.v + cross(q1.v, q2.v));
}

/**
 * @brief Scalar multiplication (quaternion * scalar).
 * @param q Quaternion.
 * @param s Scalar.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(const Quaternion& q, float s) {
  return Quaternion(q.r * s, q.v * s);
}

/**
 * @brief Scalar multiplication (scalar * quaternion).
 * @param s Scalar.
 * @param q Quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(float s, const Quaternion& q) {
  return q * s;
}

/**
 * @brief Scalar division.
 * @param q Quaternion.
 * @param s Scalar.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator/(const Quaternion& q, float s) {
  return Quaternion(q.r / s, q.v / s);
}

/**
 * @brief Calculates the dot product of two quaternions.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The dot product (scalar).
 */
constexpr float dot(const Quaternion& q1, const Quaternion& q2) {
  return (q1.r * q2.r) + (q1.v.i * q2.v.i) + (q1.v.j * q2.v.j) + (q1.v.k * q2.v.k);
}

/**
 * @brief Calculates the angle between two unit quaternions.
 * @param q1 First unit quaternion.
 * @param q2 Second unit quaternion.
 * @return The angle in radians.
 */
constexpr float angle_between(const Quaternion& q1, const Quaternion& q2) {
  return acosf(std::clamp(dot(q1, q2), -1.0f, 1.0f));
}

/**
 * @brief Creates a rotation quaternion from an axis and an angle.
 * @param axis The rotation axis (must be unit vector).
 * @param theta The angle of rotation in radians.
 * @return The resulting unit rotation quaternion.
 */
Quaternion make_rotation(const Vector& axis, float theta) {
  return Quaternion(cosf(theta / 2), sinf(theta / 2) * axis).normalize();
}

/**
 * @brief Creates a rotation quaternion to rotate from one vector to another.
 * @param from The source vector (must be unit vector).
 * @param to The destination vector (must be unit vector).
 * @return The resulting unit rotation quaternion.
 */
Quaternion make_rotation(const Vector& from, const Vector& to) {
  auto axis = cross(from, to).normalize();
  auto angle = angle_between(from, to);
  return make_rotation(axis, angle);
}

/**
 * @brief Spherical Linear Interpolation (SLERP) between two quaternions.
 * @param q1 Starting quaternion.
 * @param q2 Ending quaternion.
 * @param t Interpolation factor (0.0 to 1.0).
 * @param long_way If true, takes the longest path between quaternions.
 * @return The interpolated unit quaternion.
 */
Quaternion slerp(const Quaternion& q1, const Quaternion& q2, float t, bool long_way = false) {
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

  float theta = acosf(std::clamp(d, -1.0f, 1.0f));
  float sin_theta = sinf(theta);
  float s1 = sinf((1 - t) * theta) / sin_theta;
  float s2 = sinf(t * theta) / sin_theta;
  return ((s1 * p) + (s2 * q)).normalize();
}

/**
 * @brief Checks if a vector is "over" (on the positive side of) a plane defined by a normal vector.
 * @param v The vector to check.
 * @param normal The normal vector defining the plane (must be unit vector).
 * @return True if the dot product is non-negative.
 */
bool is_over(const Vector& v, const Vector& normal) {
  return dot(normal, v) >= 0;
}

/**
 * @brief Checks if the segment between two vectors intersects the plane defined by a normal.
 * @param v1 The starting vector.
 * @param v2 The ending vector.
 * @param normal The plane normal (must be unit vector).
 * @return True if v1 and v2 are on opposite sides of the plane.
 */
bool intersects_plane(const Vector& v1, const Vector& v2, const Vector& normal) {
  return (is_over(v1, normal) && !is_over(v2, normal))
    || (!is_over(v1, normal) && is_over(v2, normal));
}

/**
 * @brief Calculates the intersection point of an arc (between u and v) and a plane (defined by normal).
 * @param u The starting vector (must be unit vector).
 * @param v The ending vector (must be unit vector).
 * @param normal The plane normal (must be unit vector).
 * @return The intersection vector (unit vector). Asserts on failure.
 */
Vector intersection(const Vector& u, const Vector& v, const Vector& normal) {
  Vector w = cross(v, u).normalize();
  Vector i1 = cross(w, normal).normalize();
  Vector i2 = cross(normal, w).normalize();

  float a1 = angle_between(u, v);
  float a2 = angle_between(i1, u);
  float a3 = angle_between(i1, v);
  if (std::abs(a2 + a3 - a1) < TOLERANCE) {
    return i1;
  }

  a1 = angle_between(u, v);
  a2 = angle_between(i2, u);
  a3 = angle_between(i2, v);
  if (std::abs(a2 + a3 - a1) < TOLERANCE) {
    return i2;
  }

  return Vector(0, 0, 0);
}

/**
 * @brief Constant derived from the screen width for splitting segments.
 */
static constexpr float SHIFT_DIV = 96;

/**
 * @brief Creates two slightly shifted vectors based on the intersection point and the normal.
 * @param v The intersection vector (must be unit vector).
 * @param normal The plane normal (must be unit vector).
 * @return A vector containing two normalized vectors, shifted slightly apart along the plane normal.
 */
std::vector<Vector> split_point(const Vector& v, const Vector& normal) {
  float shift = sinf(PI_F / SHIFT_DIV);
  return {
    {(v + (normal * shift)).normalize()},
    {(v + (-normal * shift)).normalize()},
  };
}