/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file rotation_math.h @brief Vector and quaternion arithmetic and rotation construction. */

#include "core/math/quaternion_complex.h"

namespace math {

/**
 * @brief Rotates a vector by a unit quaternion.
 * @details Expanded q*v*conj(q): t = 2 (q.v x v), v + q.r t + q.v x t. The
 * doubling is applied to q.v, which is exact.
 * @param v The vector to rotate.
 * @param q The unit rotation quaternion.
 * @return The rotated vector.
 */
__attribute__((always_inline)) inline Vector rotate(const Vector &v,
                                                    const Quaternion &q) {
  const float qr = q.r, qx = q.v.x, qy = q.v.y, qz = q.v.z;
  const float x2 = qx + qx, y2 = qy + qy, z2 = qz + qz;
  const float tx = y2 * v.z - z2 * v.y;
  const float ty = z2 * v.x - x2 * v.z;
  const float tz = x2 * v.y - y2 * v.x;
  return Vector(v.x + qr * tx + (qy * tz - qz * ty),
                v.y + qr * ty + (qz * tx - qx * tz),
                v.z + qr * tz + (qx * ty - qy * tx));
}

/**
 * @brief Vector addition.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The resulting vector.
 */
__attribute__((always_inline)) constexpr Vector operator+(const Vector &v1,
                                                          const Vector &v2) {
  return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

/**
 * @brief Vector subtraction.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The resulting vector.
 */
__attribute__((always_inline)) constexpr Vector operator-(const Vector &v1,
                                                          const Vector &v2) {
  return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

/**
 * @brief Scalar multiplication (vector * scalar).
 * @param v Vector.
 * @param s Scalar.
 * @return The resulting vector.
 */
__attribute__((always_inline)) constexpr Vector operator*(const Vector &v,
                                                          float s) {
  return Vector(s * v.x, s * v.y, s * v.z);
}

/**
 * @brief Scalar multiplication (scalar * vector).
 * @param s Scalar.
 * @param v Vector.
 * @return The resulting vector.
 */
constexpr Vector operator*(float s, const Vector &v) { return v * s; }

// No release-build s == 0 trap (yields ±Inf/NaN); only an audit check.
constexpr Vector operator/(const Vector &v, float s) {
  HS_AUDIT_CHECK(s != 0.0f);
  return Vector(v.x / s, v.y / s, v.z / s);
}

/**
 * @brief Calculates the dot product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The dot product (scalar).
 */
__attribute__((always_inline)) constexpr float dot(const Vector &v1,
                                                   const Vector &v2) {
  return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

/**
 * @brief Calculates the cross product of two vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return The cross product vector.
 */
__attribute__((always_inline)) constexpr Vector cross(const Vector &v1,
                                                      const Vector &v2) {
  return Vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                v1.x * v2.y - v1.y * v2.x);
}

/**
 * @brief A unit quaternion expanded to its 3x3 rotation matrix.
 * @details Amortizes the quaternion sandwich over many vectors sharing one
 *          orientation: build once, then apply is three dot products.
 */
struct RotationMatrix {
  /// First matrix row; rotating v yields (dot(r0, v), dot(r1, v), dot(r2, v)).
  Vector r0;
  Vector r1; ///< Second matrix row.
  Vector r2; ///< Third matrix row.

  /** @brief The identity rotation. */
  RotationMatrix()
      : r0(1.0f, 0.0f, 0.0f), r1(0.0f, 1.0f, 0.0f), r2(0.0f, 0.0f, 1.0f) {}

  /**
   * @brief Expands a unit quaternion into its rows.
   * @param q The unit rotation quaternion.
   * @pre `q` is unit magnitude; a non-unit q distorts the transform and the error is
   *      baked into every vector the cached object goes on to rotate.
   */
  explicit RotationMatrix(const Quaternion &q) {
    HS_AUDIT_CHECK(
        !(std::fabs(q.squared_magnitude() - 1.0f) > math::EPS_UNIT_QUAT_SQ));
    const float qr = q.r, qx = q.v.x, qy = q.v.y, qz = q.v.z;
    const float x2 = qx + qx, y2 = qy + qy, z2 = qz + qz;
    const float xx = qx * x2, yy = qy * y2, zz = qz * z2;
    const float xy = qx * y2, xz = qx * z2, yz = qy * z2;
    const float rx = qr * x2, ry = qr * y2, rz = qr * z2;
    r0 = Vector(1.0f - (yy + zz), xy - rz, xz + ry);
    r1 = Vector(xy + rz, 1.0f - (xx + zz), yz - rx);
    r2 = Vector(xz - ry, yz + rx, 1.0f - (xx + yy));
  }

  /**
   * @brief Rotates a vector by the stored orientation.
   * @param v The vector to rotate.
   * @return The rotated vector.
   */
  __attribute__((always_inline)) Vector apply(const Vector &v) const {
    return Vector(dot(r0, v), dot(r1, v), dot(r2, v));
  }
};

/**
 * @brief Maps vectors from one world orientation into another.
 * @param from Orientation the input vectors have already passed through.
 * @param to Destination orientation.
 * @return `to` times the transpose of `from`.
 */
inline RotationMatrix rotation_between(const RotationMatrix &from,
                                       const RotationMatrix &to) {
  RotationMatrix result;
  result.r0 = {dot(to.r0, from.r0), dot(to.r0, from.r1), dot(to.r0, from.r2)};
  result.r1 = {dot(to.r1, from.r0), dot(to.r1, from.r1), dot(to.r1, from.r2)};
  result.r2 = {dot(to.r2, from.r0), dot(to.r2, from.r1), dot(to.r2, from.r2)};
  return result;
}

/**
 * @brief Calculates the Euclidean distance between two vectors.
 * @param a First vector.
 * @param b Second vector.
 * @return The distance (scalar).
 */
inline float distance_between(const Vector &a, const Vector &b) {
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
 * @note Traps unless both squared lengths are at least math::EPS_LEN_SQ, a
 * stricter floor than normalized()'s.
 */
inline float angle_between(const Vector &v1, const Vector &v2) {
  float m1 = dot(v1, v1);
  float m2 = dot(v2, v2);
  HS_CHECK(m1 >= math::EPS_LEN_SQ && m2 >= math::EPS_LEN_SQ,
           "angle_between: degenerate vector");
  float d = dot(v1, v2) / sqrtf(m1 * m2);
  return fast_acos(hs::clamp(d, -1.0f, 1.0f));
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

/**
 * @brief Creates a rotation quaternion from an axis and an angle.
 * @param axis The rotation axis (must be a non-zero unit vector).
 * @param theta The angle of rotation in radians.
 * @return The resulting unit rotation quaternion.
 * @pre `axis` is unit length. A zero axis traps at `theta = pi`; a non-unit
 *      axis does not trap, but the realized rotation angle is not `theta`.
 */
inline Quaternion make_rotation(const Vector &axis, float theta) {
  return Quaternion(cosf(theta / 2), sinf(theta / 2) * axis).normalized();
}

/**
 * @brief Returns +X, or +Y when @p v is near-parallel to +X.
 * @param v Vector to find a well-conditioned reference axis for.
 * @return +X unless @p v is near-parallel to it, in which case +Y.
 * @details Seed for building a cross-product frame around @p v.
 */
inline Vector least_parallel_axis(const Vector &v) {
  // Scale-invariant: compare v.x^2 against the threshold scaled by |v|^2, so a
  // non-unit v can't slip past the cosine test and seed a degenerate cross.
  constexpr float K = math::COS_AXIS_PARALLEL;
  return v.x * v.x > K * K * dot(v, v) ? Vector(0, 1, 0) : Vector(1, 0, 0);
}

/**
 * @brief Returns a unit vector perpendicular to @p v.
 * @param v Nonzero vector to find a perpendicular for; need not be unit.
 * @pre The cross product with least_parallel_axis(v) must have squared
 * magnitude at least EPS_NORMALIZE_SQ.
 * @return cross(v, least_parallel_axis(v)), normalized.
 * @details Which perpendicular is unspecified, only that it is unit and
 *          well-conditioned for every @p v.
 */
inline Vector perpendicular_axis(const Vector &v) {
  return cross(v, least_parallel_axis(v)).normalized();
}

/**
 * @brief Creates a rotation quaternion to rotate from one vector to another.
 * @param from The source vector (must be unit vector).
 * @param to The destination vector (must be unit vector).
 * @return The resulting unit rotation quaternion.
 * @details Half-vector form q = (1 + d, cross(from, to)) normalized, subject to
 * floating-point rounding and underflow. The antiparallel case (|cross|² below
 * EPS_ANTIPARALLEL_SQ) synthesizes a stable perpendicular axis and rotates by π.
 */
inline Quaternion make_rotation(const Vector &from, const Vector &to) {
  HS_CHECK(std::abs(dot(from, from) - 1.0f) < math::EPS_UNIT_VEC_SQ &&
               std::abs(dot(to, to) - 1.0f) < math::EPS_UNIT_VEC_SQ,
           "make_rotation(from, to): inputs must be unit vectors");
  float d = dot(from, to);
  Vector axis = cross(from, to);

  // Antiparallel (180°): seed the cross with least_parallel_axis(from) so the
  // axis we normalize through stays well-conditioned. Keyed on |cross|², not d:
  // d's ~1e-7 absolute error near -1 exceeds the crossover, so a raw-dot gate
  // misses exact antipodes whose float dot rounds above -1.
  if (d < 0.0f && dot(axis, axis) < math::EPS_ANTIPARALLEL_SQ) {
    Vector perp = cross(least_parallel_axis(from), from);
    perp.normalize();
    return make_rotation(perp, PI_F);
  }

  // |q|² = (1 + d)² + |cross|² >= EPS_ANTIPARALLEL_SQ, so normalized() is safe.
  return Quaternion(1.0f + d, axis).normalized();
}

/**
 * @brief Builds the unit rotation quaternion whose columns are the given
 * orthonormal axes.
 * @param cx Image of the body +X axis (must be unit).
 * @param cy Image of the body +Y axis (must be unit, orthogonal to cx).
 * @param cz Image of the body +Z axis (must be unit, = cross(cx, cy)).
 * @return The unit quaternion q with rotate(X_AXIS, q) == cx, etc.
 * @details Uses Shepperd's method with the largest quaternion component.
 * Orthonormality is an unchecked caller precondition; degenerate inputs need
 * not trap.
 */
inline Quaternion quaternion_from_basis(const Vector &cx, const Vector &cy,
                                        const Vector &cz) {
  // mIJ is row I, col J of the rotation matrix (column = image of a body axis).
  const float m00 = cx.x, m10 = cx.y, m20 = cx.z;
  const float m01 = cy.x, m11 = cy.y, m21 = cy.z;
  const float m02 = cz.x, m12 = cz.y, m22 = cz.z;
  const float trace = m00 + m11 + m22;
  Quaternion q;
  if (trace > 0.0f) {
    float s = sqrtf(trace + 1.0f) * 2.0f; // s = 4 * q.r
    q.r = 0.25f * s;
    q.v.x = (m21 - m12) / s;
    q.v.y = (m02 - m20) / s;
    q.v.z = (m10 - m01) / s;
  } else if (m00 > m11 && m00 > m22) {
    float s = sqrtf(1.0f + m00 - m11 - m22) * 2.0f; // s = 4 * q.v.x
    q.r = (m21 - m12) / s;
    q.v.x = 0.25f * s;
    q.v.y = (m01 + m10) / s;
    q.v.z = (m02 + m20) / s;
  } else if (m11 > m22) {
    float s = sqrtf(1.0f + m11 - m00 - m22) * 2.0f; // s = 4 * q.v.y
    q.r = (m02 - m20) / s;
    q.v.x = (m01 + m10) / s;
    q.v.y = 0.25f * s;
    q.v.z = (m12 + m21) / s;
  } else {
    float s = sqrtf(1.0f + m22 - m00 - m11) * 2.0f; // s = 4 * q.v.z
    q.r = (m10 - m01) / s;
    q.v.x = (m02 + m20) / s;
    q.v.y = (m12 + m21) / s;
    q.v.z = 0.25f * s;
  }
  return q.normalized();
}

} // namespace math
