/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file vector.h @brief Cartesian and spherical vector types and snorm16 direction packing. */

#include "core/math/scalar_noise.h"

namespace math {

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
   * @param v The Cartesian vector; its squared length must be finite and at
   *          least EPS_NORMALIZE_SQ. normalize() traps below that threshold.
   * @note Uses fast_atan2/fast_acos (peak errors ~3.8e-3 rad azimuth and
   *       ~5e-5 rad polar). The round trip through Vector(Spherical) is not exact.
   */
  explicit Spherical(const Vector &v);

  float theta; /**< Azimuthal angle from +x toward +z (atan2 convention; the
                    (theta, phi) ctor stores any value unchecked). */
  float phi;   /**< Co-latitude from +y, in [0, pi] (0 at +y, pi at -y); the
                    (theta, phi) ctor stores any value unchecked. */
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
   * @brief Constructs a zero vector from a null-pointer tag.
   */
  explicit constexpr Vector(std::nullptr_t) : x(0), y(0), z(0) {}
  /**
   * @brief Constructs a vector with explicit components.
   * @param x X-component.
   * @param y Y-component.
   * @param z Z-component.
   */
  constexpr Vector(float x, float y, float z) : x(x), y(y), z(z) {}
  /**
   * @brief Copy constructor (defaulted, trivial).
   * @param v The vector to copy.
   */
  constexpr Vector(const Vector &v) = default;
  /**
   * @brief Constructs a vector from Spherical coordinates.
   * @param s The spherical coordinate.
   */
  explicit Vector(const Spherical &s)
      : x(sinf(s.phi) * cosf(s.theta)), y(cosf(s.phi)),
        z(sinf(s.phi) * sinf(s.theta)) {}
  /**
   * @brief Builds a unit vector from spherical angles.
   * @param theta Azimuthal angle.
   * @param phi Polar angle.
   * @return Unit vector at (theta, phi).
   */
  static Vector from_spherical(float theta, float phi) {
    return Vector(Spherical(theta, phi));
  }

  /**
   * @brief Copy assignment (defaulted, trivial).
   * @param v The vector to copy from.
   * @return Reference to this vector.
   */
  constexpr Vector &operator=(const Vector &v) = default;

  /**
   * @brief Equality comparison using tolerance.
   * @param v The vector to compare.
   * @return True if components are within TOLERANCE.
   * @note Absolute per-component tolerance: scale-dependent and
   *   non-transitive, so unsuitable for ordering or container keys.
   */
  bool operator==(const Vector &v) const {
    return std::abs(x - v.x) <= math::TOLERANCE &&
           std::abs(y - v.y) <= math::TOLERANCE &&
           std::abs(z - v.z) <= math::TOLERANCE;
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

  /**
   * @brief Adds another vector component-wise into this vector.
   * @param v_other The vector to add.
   * @return Reference to this vector.
   */
  Vector &operator+=(const Vector &v_other) {
    x += v_other.x;
    y += v_other.y;
    z += v_other.z;
    return *this;
  }
  /**
   * @brief Subtracts another vector component-wise from this vector.
   * @param v_other The vector to subtract.
   * @return Reference to this vector.
   */
  Vector &operator-=(const Vector &v_other) {
    x -= v_other.x;
    y -= v_other.y;
    z -= v_other.z;
    return *this;
  }
  /**
   * @brief Scales this vector by a scalar in place.
   * @param s Scalar multiplier.
   * @return Reference to this vector.
   */
  Vector &operator*=(float s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  /**
   * @brief Divides this vector by a scalar in place.
   * @param s Scalar divisor; a zero yields the same +-Inf/NaN operator/(Vector,
   *          float) does, guarded only by an audit check.
   * @return Reference to this vector.
   */
  Vector &operator/=(float s) {
    HS_AUDIT_CHECK(s != 0.0f);
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  /**
   * @brief Calculates the magnitude (length) of the vector.
   * @return The magnitude.
   */
  [[nodiscard]] __attribute__((always_inline)) float length() const {
    return sqrtf(x * x + y * y + z * z);
  }

  /**
   * @brief Alias for length().
   * @return The magnitude (length) of the vector.
   */
  [[nodiscard]] float magnitude() const { return length(); }

  /**
   * @brief Normalizes the vector in place (scales to unit length).
   * @pre The squared length is finite and at least EPS_NORMALIZE_SQ.
   */
  void normalize() {
    float m2 = x * x + y * y + z * z;
    HS_CHECK(m2 >= math::EPS_NORMALIZE_SQ, "Vector::normalize: degenerate");
    float m = sqrtf(m2);
    x = x / m;
    y = y / m;
    z = z / m;
  }

  /**
   * @brief Returns a unit-length copy without mutating `this`.
   * @pre The squared length is finite and at least EPS_NORMALIZE_SQ.
   * @return A unit-length vector in the same direction.
   * @details Use normalized_or() for degenerate vectors.
   */
  [[nodiscard]] __attribute__((always_inline)) Vector normalized() const {
    float m2 = x * x + y * y + z * z;
    HS_CHECK(m2 >= math::EPS_NORMALIZE_SQ, "Vector::normalized: degenerate");
    float m = sqrtf(m2);
    return Vector(x / m, y / m, z / m);
  }

  float x = 0; /**< X-component. */
  float y = 0; /**< Y-component. */
  float z = 0; /**< Z-component. */
};

/**
 * @brief Normalizes, returning `fallback` when the vector is degenerate (zero
 * length) instead of trapping.
 * @param v The vector to normalize.
 * @pre The squared length of @p v is finite.
 * @param fallback Direction returned when `v` has near-zero length.
 * @return A unit-length copy of `v`, or `fallback` if `v` is degenerate.
 */
[[nodiscard]] __attribute__((always_inline)) inline Vector
normalized_or(const Vector &v, const Vector &fallback) {
  float m2 = v.x * v.x + v.y * v.y + v.z * v.z;
  if (m2 < math::EPS_NORMALIZE_SQ) {
    return fallback;
  }
  float m = sqrtf(m2);
  return Vector(v.x / m, v.y / m, v.z / m);
}

/**
 * @brief A unit-sphere direction packed into three snorm16 components.
 * @details Components are clamped to [-1, 1] and quantized on encode (max
 * chord ~2.6e-5), so a decoded value is near-unit rather than unit.
 */
struct Snorm3 {
  int16_t x; /**< X-component, snorm16. */
  int16_t y; /**< Y-component, snorm16. */
  int16_t z; /**< Z-component, snorm16. */

  static constexpr float SCALE = 32767.0f;         /**< Unit -> snorm16. */
  static constexpr float INV_SCALE = 1.0f / SCALE; /**< snorm16 -> unit. */

  /**
   * @brief Packs a direction, clamping each component to [-1, 1].
   * @param v The direction to pack.
   * @return The snorm16 triple.
   */
  static Snorm3 encode(const Vector &v) {
    return {quantize(v.x), quantize(v.y), quantize(v.z)};
  }

  /**
   * @brief Unpacks to a near-unit direction.
   * @return The decoded direction.
   */
  Vector decode() const {
    return Vector(x * INV_SCALE, y * INV_SCALE, z * INV_SCALE);
  }

private:
  __attribute__((always_inline)) static int16_t quantize(float c) {
    return static_cast<int16_t>(roundf(hs::clamp(c, -1.0f, 1.0f) * SCALE));
  }
};

} // namespace math
