/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file quaternion_complex.h @brief Quaternion and complex number types. */

#include <limits>
#include <utility>
#include "core/math/fast_math.h"

namespace math {

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
  Quaternion &operator*=(const Quaternion &q);
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
  [[nodiscard]] Quaternion inverse() const;

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
 * @details Per component, (r1 r2 - v1.v2, r1 v2 + r2 v1 + v1 x v2), with
 * no helper calls between the multiply-adds.
 * @param q1 First quaternion.
 * @param q2 Second quaternion.
 * @return The resulting quaternion.
 */
constexpr Quaternion operator*(const Quaternion &q1, const Quaternion &q2) {
  const float r1 = q1.r, x1 = q1.v.x, y1 = q1.v.y, z1 = q1.v.z;
  const float r2 = q2.r, x2 = q2.v.x, y2 = q2.v.y, z2 = q2.v.z;
  return Quaternion(r1 * r2 - ((x1 * x2 + y1 * y2) + z1 * z2),
                    (r1 * x2 + r2 * x1) + (y1 * z2 - z1 * y2),
                    (r1 * y2 + r2 * y1) + (z1 * x2 - x1 * z2),
                    (r1 * z2 + r2 * z1) + (x1 * y2 - y1 * x2));
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
 * @details A zero divisor is guarded only by an audit check.
 */
constexpr Quaternion operator/(const Quaternion &q, float s) {
  HS_AUDIT_CHECK(s != 0.0f);
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

inline Quaternion &Quaternion::operator*=(const Quaternion &q) {
  *this = *this * q;
  return *this;
}

inline Quaternion Quaternion::inverse() const {
  float sq_mag = squared_magnitude();
  HS_CHECK(sq_mag >= math::EPS_NORMALIZE_SQ,
           "Quaternion::inverse: zero magnitude");
  return conjugate() / sq_mag;
}

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
    if (fmaxf(fmaxf(fabsf(re), fabsf(im)), fmaxf(fabsf(b.re), fabsf(b.im))) >
        0x1p63f) [[unlikely]] {
      const double a = re, c = b.re, d = b.im;
      const double norm = c * c + d * d;
      return {static_cast<float>((a * c + im * d) / norm),
              static_cast<float>((im * c - a * d) / norm)};
    }
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
  if (radius < 0x1p-125f || radius > std::numeric_limits<float>::max())
      [[unlikely]] {
    const double real = z.re, imaginary = z.im;
    const double major =
        std::sqrt(0.5 * (std::hypot(real, imaginary) + std::abs(real)));
    const float minor = static_cast<float>(imaginary / (2.0 * major));
    if (z.re >= 0.0f)
      return {static_cast<float>(major), minor};
    return {fabsf(minor), copysignf(static_cast<float>(major), z.im)};
  }
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

/** @brief Computes both complex hyperbolic functions with shared factors.
 * @param z Argument.
 * @return Hyperbolic cosine followed by hyperbolic sine. */
HS_HOT_INLINE inline std::pair<Complex, Complex> ccosh_csinh(const Complex &z) {
  float real_cosh;
  float real_sinh;
  const float absolute = fabsf(z.re);
  if (absolute < 0x1p-5f || absolute > 88.0f) [[unlikely]] {
    real_cosh = coshf(z.re);
    real_sinh = sinhf(z.re);
  } else {
    const float exponential = expf(z.re);
    const float inverse = 1.0f / exponential;
    real_cosh = 0.5f * (exponential + inverse);
    real_sinh = 0.5f * (exponential - inverse);
  }
  const float cosine = cosf(z.im);
  const float sine = sinf(z.im);
  return {Complex(real_cosh * cosine, real_sinh * sine),
          Complex(real_sinh * cosine, real_cosh * sine)};
}

/** @brief Complex hyperbolic cosine.
 * @param z Argument.
 * @return Hyperbolic cosine. */
inline Complex ccosh(const Complex &z) { return ccosh_csinh(z).first; }

/** @brief Complex hyperbolic sine.
 * @param z Argument.
 * @return Hyperbolic sine. */
inline Complex csinh(const Complex &z) { return ccosh_csinh(z).second; }

/** @brief Principal inverse complex hyperbolic cosine.
 * @param z Argument.
 * @return Inverse hyperbolic cosine. */
inline Complex cacosh(const Complex &z) {
  return clog(z +
              csqrt(z - Complex(1.0f, 0.0f)) * csqrt(z + Complex(1.0f, 0.0f)));
}

} // namespace math
