/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file quaternion_complex.h @brief Quaternion and complex number types. */

#include "core/math/fast_math.h"

namespace math {

struct Quaternion;
constexpr Quaternion operator*(const Quaternion &q1, const Quaternion &q2);
constexpr Quaternion operator/(const Quaternion &q, float s);
/**
 * @brief Scalar division of a vector (forward declaration).
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
  constexpr Quaternion(const Quaternion &q) = default;
  /**
   * @brief Constructs a quaternion with explicit components.
   * @param a Real component (r).
   * @param b X imaginary component (i).
   * @param c Y imaginary component (j).
   * @param d Z imaginary component (k).
   */
  constexpr explicit Quaternion(float a, float b, float c, float d)
      : r(a), v(b, c, d) {}
  /**
   * @brief Constructs a quaternion from a scalar and a Vector part.
   * @param a Real component (r).
   * @param v Vector part (x, y, z).
   */
  constexpr explicit Quaternion(float a, const Vector &v) : r(a), v(v) {}

  /**
   * @brief Assignment operator.
   * @param q The quaternion to assign from.
   * @return Reference to the current object.
   */
  Quaternion &operator=(const Quaternion &q) = default;

  /**
   * @brief Equality comparison using tolerance.
   * @param q The quaternion to compare.
   * @return True if components are within TOLERANCE.
   * @note Absolute per-component tolerance (TOLERANCE): scale-dependent and
   *   non-transitive (see Vector::operator==). Suitable for direct value
   *   compares, not for ordering or container keys. Component equality, not
   *   rotation equality: `q` and `-q` are the same rotation but compare unequal.
   */
  bool operator==(const Quaternion &q) const {
    return std::abs(q.r - r) <= math::TOLERANCE && q.v == v;
  }

  /**
   * @brief Inequality comparison using tolerance.
   * @param q The quaternion to compare.
   * @return True if any component differs by more than TOLERANCE.
   */
  bool operator!=(const Quaternion &q) const { return !(*this == q); }

  /**
   * @brief Quaternion multiplication compound assignment.
   * @param q The quaternion to multiply by.
   * @return Reference to this quaternion.
   */
  Quaternion &operator*=(const Quaternion &q) {
    *this = *this * q;
    return *this;
  }
  /**
   * @brief Component-wise addition compound assignment.
   * @param q The quaternion to add.
   * @return Reference to this quaternion.
   */
  Quaternion &operator+=(const Quaternion &q) {
    r += q.r;
    v += q.v;
    return *this;
  }
  /**
   * @brief Component-wise subtraction compound assignment.
   * @param q The quaternion to subtract.
   * @return Reference to this quaternion.
   */
  Quaternion &operator-=(const Quaternion &q) {
    r -= q.r;
    v -= q.v;
    return *this;
  }
  /**
   * @brief Scales this quaternion by a scalar in place.
   * @param s Scalar multiplier.
   * @return Reference to this quaternion.
   */
  Quaternion &operator*=(float s) {
    r *= s;
    v *= s;
    return *this;
  }

  /**
   * @brief Calculates the squared magnitude of the quaternion.
   * @return The squared magnitude.
   */
  [[nodiscard]] __attribute__((always_inline)) constexpr float
  squared_magnitude() const {
    return r * r + v.x * v.x + v.y * v.y + v.z * v.z;
  }

  /**
   * @brief Calculates the inverse of the quaternion.
   * @return The inverse quaternion.
   */
  [[nodiscard]] Quaternion inverse() const {
    float sq_mag = squared_magnitude();
    HS_CHECK(sq_mag >= math::EPS_NORMALIZE_SQ,
             "Quaternion::inverse: zero magnitude");
    return conjugate() / sq_mag;
  }

  /**
   * @brief Calculates the conjugate of the quaternion.
   * @return The conjugate quaternion (r, -v).
   */
  [[nodiscard]] Quaternion conjugate() const {
    return Quaternion(r, -v.x, -v.y, -v.z);
  }

  /**
   * @brief Unary negation operator.
   * @return A new quaternion with negated components.
   */
  constexpr Quaternion operator-() const {
    return Quaternion(-r, -v.x, -v.y, -v.z);
  }

  /**
   * @brief Calculates the magnitude (length) of the quaternion.
   * @return The magnitude.
   */
  [[nodiscard]] float magnitude() const { return sqrtf(squared_magnitude()); }

  /**
   * @brief Normalizes the quaternion in place (scales to unit magnitude).
   * @pre The squared magnitude is finite and at least EPS_NORMALIZE_SQ.
   */
  void normalize() {
    float m2 = squared_magnitude();
    HS_CHECK(m2 >= math::EPS_NORMALIZE_SQ,
             "Quaternion::normalize: zero magnitude");
    float m = sqrtf(m2);
    r = r / m;
    v = v / m;
  }

  /**
   * @brief Returns a unit-magnitude copy without mutating `this`.
   * @pre The squared magnitude is finite and at least EPS_NORMALIZE_SQ.
   * @return A unit-magnitude quaternion.
   * @details Traps on a zero-magnitude quaternion — a degenerate input is a
   * logic bug.
   */
  [[nodiscard]] __attribute__((always_inline)) Quaternion normalized() const {
    float m2 = squared_magnitude();
    HS_CHECK(m2 >= math::EPS_NORMALIZE_SQ, "Quaternion::normalized: zero");
    float m = sqrtf(m2);
    return Quaternion(r / m, v / m);
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
   * @brief Sum of the squared real and imaginary components.
   * @return re^2 + im^2.
   */
  constexpr float squared_magnitude() const { return re * re + im * im; }

  /**
   * @brief Euclidean magnitude, with ordinary float overflow behavior.
   * @return sqrt(re^2 + im^2).
   */
  float magnitude() const { return sqrtf(squared_magnitude()); }

  /**
   * @brief Complex conjugate, reflecting the imaginary component.
   * @return (re, -im).
   */
  constexpr Complex conjugate() const { return Complex(re, -im); }

  /**
   * @brief Exact component equality; NaNs compare unequal.
   * @return True when both components compare equal.
   */
  constexpr bool operator==(const Complex &) const = default;

  /**
   * @brief Complex addition.
   * @param b The right-hand operand.
   * @return The sum of this and `b`.
   */
  constexpr Complex operator+(const Complex &b) const {
    return Complex(re + b.re, im + b.im);
  }
  /**
   * @brief Complex subtraction.
   * @param b The right-hand operand.
   * @return The difference of this and `b`.
   */
  constexpr Complex operator-(const Complex &b) const {
    return Complex(re - b.re, im - b.im);
  }
  /**
   * @brief Complex multiplication.
   * @param b The right-hand operand.
   * @return The product of this and `b`.
   */
  constexpr Complex operator*(const Complex &b) const {
    return Complex(re * b.re - im * b.im, re * b.im + im * b.re);
  }

  /**
   * @brief Complex division.
   * @param b The divisor.
   * @return The quotient of this and `b`. This is ordinary complex division:
   * a zero divisor yields non-finite components. For the stereographic/Mobius
   * point-at-infinity conventions (magnitude clamp, 0/0 -> 0) use project_div().
   */
  Complex operator/(const Complex &b) const {
    float denom = b.squared_magnitude();
    if (denom < std::numeric_limits<float>::min() &&
        (b.re != 0.0f || b.im != 0.0f)) [[unlikely]] {
      const float den_re = b.re * COMPLEX_UNDERFLOW_LIFT;
      const float den_im = b.im * COMPLEX_UNDERFLOW_LIFT;
      denom = den_re * den_re + den_im * den_im;
      return Complex(
          ((re * den_re + im * den_im) / denom) * COMPLEX_UNDERFLOW_LIFT,
          ((im * den_re - re * den_im) / denom) * COMPLEX_UNDERFLOW_LIFT);
    }
    return Complex((re * b.re + im * b.im) / denom,
                   (im * b.re - re * b.im) / denom);
  }
};

/** @brief Principal complex square root.
 * @param z Argument.
 * @return Square root with nonnegative real part. */
inline Complex csqrt(const Complex &z) {
  const float radius = hypotf(z.re, z.im);
  if (radius == 0.0f)
    return Complex(0.0f, z.im);
  const float major = sqrtf(0.5f * radius + 0.5f * fabsf(z.re));
  if (z.re >= 0.0f)
    return Complex(major, z.im / (2.0f * major));
  return Complex(fabsf(z.im) / (2.0f * major), copysignf(major, z.im));
}

/** @brief Principal complex logarithm.
 * @param z Argument.
 * @return Logarithm with imaginary part in [-pi, pi]. */
inline Complex clog(const Complex &z) {
  return Complex(logf(hypotf(z.re, z.im)), atan2f(z.im, z.re));
}

/** @brief Complex hyperbolic cosine.
 * @param z Argument.
 * @return Hyperbolic cosine. */
inline Complex ccosh(const Complex &z) {
  const float exponential = expf(z.re);
  const float inverse = 1.0f / exponential;
  return Complex(0.5f * (exponential + inverse) * cosf(z.im),
                 0.5f * (exponential - inverse) * sinf(z.im));
}

/** @brief Complex hyperbolic sine.
 * @param z Argument.
 * @return Hyperbolic sine. */
inline Complex csinh(const Complex &z) {
  return Complex(sinhf(z.re) * cosf(z.im), coshf(z.re) * sinf(z.im));
}

/** @brief Principal inverse complex hyperbolic cosine.
 * @param z Argument.
 * @return Inverse hyperbolic cosine. */
inline Complex cacosh(const Complex &z) {
  return clog(z +
              csqrt(z - Complex(1.0f, 0.0f)) * csqrt(z + Complex(1.0f, 0.0f)));
}

} // namespace math
