/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spherical.h
 * @brief Spherical generators, bases, and transport helpers.
 */

#include "math/3dmath.h"
#include <utility>

namespace math {
/**
 * @brief Unit vector along the Cartesian X-axis.
 */
inline constexpr Vector X_AXIS(1, 0, 0);
/**
 * @brief Unit vector along the Cartesian Y-axis.
 */
inline constexpr Vector Y_AXIS(0, 1, 0);
/**
 * @brief Unit vector along the Cartesian Z-axis.
 */
inline constexpr Vector Z_AXIS(0, 0, 1);
/**
 * @brief Unit vector along the Cartesian Y-axis.
 */
inline constexpr Vector UP = Y_AXIS;

/**
 * @brief Calculates a point on the Fibonacci spiral on the unit sphere.
 * @param n The total number of points in the spiral (must be positive).
 * @param eps The epsilon offset for the spiral.
 * @param i The index of the point to calculate.
 * @return The point on the unit sphere.
 * @note Setup-time generator; exact trig is intentional (not a per-pixel path).
 */
inline Vector fib_spiral(int n, float eps, int i) {
  HS_CHECK(n > 0, "fib_spiral: n must be positive");
  // Clamp keeps |y| <= 1, so the ring radius never takes sqrtf of a negative.
  float y = hs::clamp(1.0f - (2.0f * (static_cast<float>(i) + eps)) /
                                 static_cast<float>(n),
                      -1.0f, 1.0f);
  float radius = sqrtf(1.0f - y * y);
  constexpr double INV_PHI_PRECISE = 0.6180339887498948482;
  const float theta = static_cast<float>(
      std::fmod(2.0 * PI * static_cast<double>(i) * INV_PHI_PRECISE, 2.0 * PI));
  // Y-up; unit by construction, so no normalize().
  return Vector(radius * cosf(theta), y, radius * sinf(theta));
}

/**
 * @brief Generates a uniformly distributed 3D unit vector (direction) from the
 * deterministic global RNG using Marsaglia's method.
 * @return A normalized random Vector.
 */
inline Vector random_vector() {
  float v1, v2, s;
  // Marsaglia rejection: accept when (v1,v2) lands in the open unit disk.
  do {
    v1 = 2.0f * hs::rand_f() - 1.0f;
    v2 = 2.0f * hs::rand_f() - 1.0f;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1.0f || s == 0.0f);

  float sqrt_one_minus_s = sqrtf(1.0f - s);
  return Vector(2.0f * v1 * sqrt_one_minus_s, 2.0f * v2 * sqrt_one_minus_s,
                1.0f - 2.0f * s);
}

/**
 * @brief Parameters defining a Lissajous curve.
 */
struct LissajousParams {
  float m1; /**< Frequency coefficient for the axial components (X and Z). */
  float m2; /**< Frequency coefficient for the orbital component (Y). */
  float a; /**< Phase shift in radians (matches the daydream lissajous tool). */
  float domain; /**< The total duration (t) over which the curve is drawn. */
};

/**
 * @brief Calculates a 3D point on the unit sphere corresponding to a spherical
 * Lissajous curve.
 * @param m1 Frequency coefficient for XZ plane.
 * @param m2 Frequency coefficient for Y axis.
 * @param a Phase shift in radians, used as-is (matches the daydream lissajous
 *          designer tools/lissajous.html).
 * @param t Time variable (or position along the domain).
 * @return The calculated 3D point (unit vector).
 * @note Setup-time generator; exact trig is intentional (not a per-pixel path).
 */
inline Vector lissajous(float m1, float m2, float a, float t) {
  // Unit by construction, so no normalize().
  return Vector(sinf(m2 * t) * cosf(m1 * t - a), cosf(m2 * t),
                sinf(m2 * t) * sinf(m1 * t - a));
}

/**
 * @brief An orthonormal basis { u, v, w } whose *middle* axis is the normal:
 * `v` is the normal, and `u`, `w` span the plane perpendicular to it.
 */
struct Basis {
  Vector u, v, w;
};

/**
 * @brief Rotates a basis by a unit quaternion; the axes stay orthonormal.
 * @param b Basis to rotate.
 * @param q Unit rotation quaternion.
 * @return The basis with every axis rotated by @p q.
 */
inline Basis rotate(const Basis &b, const Quaternion &q) {
  return {rotate(b.u, q), rotate(b.v, q), rotate(b.w, q)};
}

/**
 * @brief Creates a basis { u, v, w } from an orientation and normal.
 * @param orientation The orientation quaternion; MUST be unit length
 *   (HS_CHECK-trapped below — a non-unit quaternion would scale/shear the frame).
 * @param normal The normal vector; after rotation it becomes the 'v' axis. MUST
 *   be non-zero — a zero (or rotation-collapsed) normal traps in the
 *   `normalized()` of `v` below.
 * @return The constructed Basis.
 */
inline Basis make_basis(const Quaternion &orientation, const Vector &normal) {
  const float orientation_norm_sq = orientation.squared_magnitude();
  HS_CHECK(std::abs(orientation_norm_sq - 1.0f) < math::EPS_UNIT_QUAT_SQ,
           "make_basis: orientation |q|^2 is %d/1000, not unit",
           orientation_norm_sq < 1.0e6f
               ? static_cast<int>(orientation_norm_sq * 1000.0f)
               : static_cast<int>(INT32_MIN));
  Vector v = rotate(normal, orientation).normalized();
  // rotate preserves dot, so least_parallel_axis(normal) picks the same body
  // axis as the rotated frame; rotate it into the frame for the cross. Only its
  // direction matters, the cross below is normalized.
  Vector ref = rotate(least_parallel_axis(normal), orientation);
  Vector u = cross(v, ref).normalized();
  // v and u are orthonormal, so the cross is unit by construction.
  Vector w = cross(v, u);
  return {u, v, w};
}

/**
 * @brief First tangent-basis vector at a unit vertex.
 * @details The frame's second vector is cross(normal, u), unit because normal
 *          and u are orthonormal; callers derive it rather than store it.
 *          Crosses against +Y, swapping to +X within POLE_REFERENCE_COS of a
 *          pole where the +Y cross collapses. The band is far wider than the one
 *          perpendicular_axis() seeds with, which crosses against +X instead.
 * @param normal Unit vertex normal.
 * @return A unit tangent at @p normal.
 */
HS_FLASH_INLINE inline Vector tangent_axis(const Vector &normal) {
  constexpr float POLE_REFERENCE_COS = 0.99f;
  const Vector axis = std::abs(normal.y) > POLE_REFERENCE_COS ? X_AXIS : Y_AXIS;
  return cross(normal, axis).normalized();
}

/** Conditioning bound on |cross(from, to)|^2 for parallel_transport; below it
 *  the great circle is ill-determined. */
inline constexpr float MIN_TRANSPORT_CROSS_SQ = 1e-4f;

/**
 * @brief Parallel-transports a tangent along the great-circle arc between two
 *        unit vectors.
 * @param from Unit vector the tangent is attached to.
 * @param to Unit vector the tangent is carried to.
 * @param tangent Tangent at @p from.
 * @return The tangent at @p to.
 * @details Traps once the pair is inside MIN_TRANSPORT_CROSS_SQ of antipodal
 *   rather than only at the exact antipode: any non-tangent component of
 *   @p tangent is amplified by 2/|cross(from, to)| there, capped at 200x by the
 *   bound. Gated on |cross|^2, which stays accurate where 1 + dot cancels; the
 *   dot > 0 half-space short-circuits before the cross product.
 */
inline Vector parallel_transport(const Vector &from, const Vector &to,
                                 const Vector &tangent) {
  const float denominator = 1.0f + dot(from, to);
  HS_CHECK(denominator > 1.0f ||
               dot(cross(from, to), cross(from, to)) > MIN_TRANSPORT_CROSS_SQ,
           "parallel_transport: antipodal endpoints");
  return tangent - (from + to) * (dot(tangent, to) / denominator);
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

} // namespace math
