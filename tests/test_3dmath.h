/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/3dmath.h.
 *
 * Usage:
 *   #include "tests/test_3dmath.h"
 *   int main() { return hs_test::math3d::run_3dmath_tests(); }
 *
 * Self-contained header — no external test framework. All test functions
 * are inline; the runner returns the failure count for use as a process
 * exit code.
 */
#pragma once

#include "core/3dmath.h"
#include "core/rotate.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace math3d {

inline bool approx_vec(const Vector &a, const Vector &b, float tol) {
  return approx(a.x, b.x, tol) && approx(a.y, b.y, tol) &&
         approx(a.z, b.z, tol);
}
inline bool approx_quat(const Quaternion &a, const Quaternion &b, float tol) {
  return approx(a.r, b.r, tol) && approx_vec(a.v, b.v, tol);
}
inline bool approx_complex(const Complex &a, const Complex &b, float tol) {
  return approx(a.re, b.re, tol) && approx(a.im, b.im, tol);
}

#define HS_EXPECT_VEC(a, b, tol)                                               \
  HS_EXPECT(hs_test::math3d::approx_vec((a), (b), (tol)),                      \
            #a " ~= " #b " (tol=" #tol ")")
#define HS_EXPECT_QUAT(a, b, tol)                                              \
  HS_EXPECT(hs_test::math3d::approx_quat((a), (b), (tol)),                     \
            #a " ~= " #b " (tol=" #tol ")")
#define HS_EXPECT_COMPLEX(a, b, tol)                                           \
  HS_EXPECT(hs_test::math3d::approx_complex((a), (b), (tol)),                  \
            #a " ~= " #b " (tol=" #tol ")")

// ============================================================================
// Constants
// ============================================================================

inline void test_constants() {
  HS_EXPECT_NEAR(PHI, 1.61803398f, 1e-6f);
  HS_EXPECT_NEAR(G, 1.0f / PHI, 1e-6f);
  HS_EXPECT_NEAR(G * PHI, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::TOLERANCE, 0.0001f, 1e-9f);
  HS_EXPECT_NEAR(TOLERANCE, math::TOLERANCE, 1e-9f);
  HS_EXPECT_NEAR(PI_F, 3.14159265f, 1e-5f);
  HS_EXPECT_NEAR(STEREO_INF, 1e4f, 1e-3f);
}

// ============================================================================
// quintic_kernel
// ============================================================================

inline void test_quintic_kernel() {
  HS_EXPECT_NEAR(quintic_kernel(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(quintic_kernel(1.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(quintic_kernel(0.5f), 0.5f, 1e-6f);

  // Clamps below 0 and above 1
  HS_EXPECT_NEAR(quintic_kernel(-1.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(quintic_kernel(2.0f), 1.0f, 1e-6f);

  // Monotone increasing on [0, 1]
  float prev = -1.0f;
  for (int i = 0; i <= 20; ++i) {
    float v = quintic_kernel(i / 20.0f);
    HS_EXPECT_TRUE(v >= prev - 1e-6f);
    prev = v;
  }

  // Smootherstep is C2-continuous: derivative ~= 0 at endpoints
  float dl = quintic_kernel(0.01f) - quintic_kernel(0.0f);
  float dr = quintic_kernel(1.0f) - quintic_kernel(0.99f);
  HS_EXPECT_TRUE(std::abs(dl) < 1e-3f);
  HS_EXPECT_TRUE(std::abs(dr) < 1e-3f);
}

// ============================================================================
// fast_atan2 / fast_acos / fast_sinf / fast_cosf
// (measured peak errors: atan2 ~3.8e-3 rad, acos ~1.3e-4 rad, sin ~1.6e-3)
// ============================================================================

inline void test_fast_atan2() {
  HS_EXPECT_NEAR(fast_atan2(0.0f, 1.0f), 0.0f, 5e-3f);
  HS_EXPECT_NEAR(fast_atan2(1.0f, 0.0f), PI_F * 0.5f, 5e-3f);
  HS_EXPECT_NEAR(fast_atan2(0.0f, -1.0f), PI_F, 5e-3f);
  HS_EXPECT_NEAR(fast_atan2(-1.0f, 0.0f), -PI_F * 0.5f, 5e-3f);

  // Sweep against std::atan2
  for (int i = 0; i < 64; ++i) {
    float a = -PI_F + (i * 2.0f * PI_F) / 64.0f;
    float y = std::sin(a);
    float x = std::cos(a);
    HS_EXPECT_NEAR(fast_atan2(y, x), std::atan2(y, x), 5e-3f);
  }
}

inline void test_fast_acos() {
  HS_EXPECT_NEAR(fast_acos(1.0f), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(fast_acos(-1.0f), PI_F, 1e-3f);
  HS_EXPECT_NEAR(fast_acos(0.0f), PI_F * 0.5f, 1e-3f);

  // Out-of-range clamps
  HS_EXPECT_NEAR(fast_acos(1.5f), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(fast_acos(-1.5f), PI_F, 1e-3f);

  // Sweep against std::acos
  for (int i = 0; i <= 32; ++i) {
    float x = -1.0f + (i / 16.0f);
    if (x > 1.0f) x = 1.0f;
    HS_EXPECT_NEAR(fast_acos(x), std::acos(x), 1e-3f);
  }
}

inline void test_fast_sinf_cosf() {
  HS_EXPECT_NEAR(fast_sinf(0.0f), 0.0f, 2e-3f);
  HS_EXPECT_NEAR(fast_sinf(PI_F * 0.5f), 1.0f, 2e-3f);
  HS_EXPECT_NEAR(fast_sinf(PI_F), 0.0f, 2e-3f);
  HS_EXPECT_NEAR(fast_sinf(-PI_F * 0.5f), -1.0f, 2e-3f);

  HS_EXPECT_NEAR(fast_cosf(0.0f), 1.0f, 2e-3f);
  HS_EXPECT_NEAR(fast_cosf(PI_F * 0.5f), 0.0f, 2e-3f);
  HS_EXPECT_NEAR(fast_cosf(PI_F), -1.0f, 2e-3f);

  // Pythagorean identity holds approximately
  for (int i = 0; i < 32; ++i) {
    float a = -3.0f * PI_F + (i * 6.0f * PI_F) / 32.0f;
    float s = fast_sinf(a);
    float c = fast_cosf(a);
    HS_EXPECT_NEAR(s * s + c * c, 1.0f, 5e-3f);
  }

  // Range reduction: periodicity over ±2π
  for (int i = 0; i < 16; ++i) {
    float a = i * 0.3f;
    HS_EXPECT_NEAR(fast_sinf(a + 2.0f * PI_F), fast_sinf(a), 3e-3f);
    HS_EXPECT_NEAR(fast_sinf(a - 2.0f * PI_F), fast_sinf(a), 3e-3f);
  }
}

// ============================================================================
// Vector — construction
// ============================================================================

inline void test_vector_construction() {
  Vector v0;
  HS_EXPECT_NEAR(v0.x, 0.0f, 1e-7f);
  HS_EXPECT_NEAR(v0.y, 0.0f, 1e-7f);
  HS_EXPECT_NEAR(v0.z, 0.0f, 1e-7f);

  Vector v_int(0); // inplace_function compat constructor
  HS_EXPECT_VEC(v_int, Vector(0, 0, 0), 1e-7f);

  Vector v(1.0f, 2.0f, 3.0f);
  HS_EXPECT_NEAR(v.x, 1.0f, 0.0f);
  HS_EXPECT_NEAR(v.y, 2.0f, 0.0f);
  HS_EXPECT_NEAR(v.z, 3.0f, 0.0f);

  Vector vc(v);
  HS_EXPECT_VEC(vc, v, 0.0f);

  Vector vassigned;
  vassigned = v;
  HS_EXPECT_VEC(vassigned, v, 0.0f);
}

inline void test_vector_spherical_construction() {
  // theta=0, phi=π/2 → +X axis
  HS_EXPECT_VEC(Vector::from_spherical(0.0f, PI_F * 0.5f), Vector(1, 0, 0),
                2e-3f);
  // theta=π/2, phi=π/2 → +Z axis
  HS_EXPECT_VEC(Vector::from_spherical(PI_F * 0.5f, PI_F * 0.5f),
                Vector(0, 0, 1), 2e-3f);
  // phi=0 → +Y (north pole)
  HS_EXPECT_VEC(Vector::from_spherical(0.0f, 0.0f), Vector(0, 1, 0), 2e-3f);
  // phi=π → -Y (south pole)
  HS_EXPECT_VEC(Vector::from_spherical(0.0f, PI_F), Vector(0, -1, 0), 2e-3f);

  // Spherical-built unit vectors stay unit
  Vector u = Vector::from_spherical(0.8f, 1.1f);
  HS_EXPECT_NEAR(u.length(), 1.0f, 2e-3f);
}

inline void test_vector_equality() {
  Vector a(1, 2, 3), b(1, 2, 3), c(1.001f, 2, 3);
  HS_EXPECT_TRUE(a == b);
  HS_EXPECT_FALSE(a == c);
  HS_EXPECT_TRUE(a != c);
  HS_EXPECT_FALSE(a != b);

  // Tolerance behavior: 0.00005 difference equal, 0.001 difference not equal
  HS_EXPECT_TRUE(Vector(1, 0, 0) == Vector(1.00005f, 0, 0));
  HS_EXPECT_FALSE(Vector(1, 0, 0) == Vector(1.001f, 0, 0));
}

// ============================================================================
// Vector — arithmetic
// ============================================================================

inline void test_vector_arithmetic() {
  Vector a(1, 2, 3), b(4, 5, 6);
  HS_EXPECT_VEC(a + b, Vector(5, 7, 9), 1e-6f);
  HS_EXPECT_VEC(b - a, Vector(3, 3, 3), 1e-6f);
  HS_EXPECT_VEC(-a, Vector(-1, -2, -3), 1e-6f);
  HS_EXPECT_VEC(a * 2.0f, Vector(2, 4, 6), 1e-6f);
  HS_EXPECT_VEC(2.0f * a, Vector(2, 4, 6), 1e-6f);
  HS_EXPECT_VEC(a / 2.0f, Vector(0.5f, 1.0f, 1.5f), 1e-6f);

  Vector c(1, 2, 3);
  c += b;
  HS_EXPECT_VEC(c, Vector(5, 7, 9), 1e-6f);
  c -= b;
  HS_EXPECT_VEC(c, Vector(1, 2, 3), 1e-6f);
  c *= 3.0f;
  HS_EXPECT_VEC(c, Vector(3, 6, 9), 1e-6f);
  c /= 3.0f;
  HS_EXPECT_VEC(c, Vector(1, 2, 3), 1e-6f);
}

inline void test_vector_length() {
  HS_EXPECT_NEAR(Vector(3, 4, 0).length(), 5.0f, 1e-6f);
  HS_EXPECT_NEAR(Vector(0, 0, 0).length(), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(Vector(1, 2, 2).magnitude(), 3.0f, 1e-6f);
}

inline void test_vector_normalize() {
  Vector v(3, 0, 4);
  v.normalize();
  HS_EXPECT_NEAR(v.length(), 1.0f, 1e-6f);
  HS_EXPECT_VEC(v, Vector(0.6f, 0.0f, 0.8f), 1e-6f);

  // normalized() does not mutate the source
  Vector u(0, 5, 0);
  Vector n = u.normalized();
  HS_EXPECT_VEC(u, Vector(0, 5, 0), 1e-6f);
  HS_EXPECT_VEC(n, Vector(0, 1, 0), 1e-6f);

  // A zero-length vector now traps under the strict normalize()/normalized();
  // the guarded normalized_or() returns the supplied fallback, and a non-zero
  // vector normalizes as usual.
  HS_EXPECT_VEC(normalized_or(Vector(0, 0, 0), Vector(1, 0, 0)), Vector(1, 0, 0),
                1e-6f);
  HS_EXPECT_VEC(normalized_or(Vector(0, 6, 0), Vector(1, 0, 0)), Vector(0, 1, 0),
                1e-6f);
}

// ============================================================================
// Vector — free functions (dot, cross, distance, angle_between)
// ============================================================================

inline void test_dot_cross() {
  Vector x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);

  HS_EXPECT_NEAR(dot(x, y), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(dot(x, x), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(dot(x, -x), -1.0f, 1e-6f);
  HS_EXPECT_NEAR(dot(Vector(1, 2, 3), Vector(4, -5, 6)), 4 - 10 + 18, 1e-5f);

  // Right-handed cross products
  HS_EXPECT_VEC(cross(x, y), z, 1e-6f);
  HS_EXPECT_VEC(cross(y, z), x, 1e-6f);
  HS_EXPECT_VEC(cross(z, x), y, 1e-6f);
  // Anticommutative
  HS_EXPECT_VEC(cross(y, x), -z, 1e-6f);
  // a × a = 0
  HS_EXPECT_VEC(cross(Vector(2, 3, 5), Vector(2, 3, 5)), Vector(0, 0, 0),
                1e-6f);
  // a × b is perpendicular to both
  Vector a(1, 2, 3), b(4, -5, 6);
  Vector c = cross(a, b);
  HS_EXPECT_NEAR(dot(c, a), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(c, b), 0.0f, 1e-4f);
}

inline void test_distance() {
  Vector a(1, 2, 3), b(4, 6, 3);
  HS_EXPECT_NEAR(distance_between(a, b), 5.0f, 1e-6f);
  HS_EXPECT_NEAR(distance_squared(a, b), 25.0f, 1e-5f);
  HS_EXPECT_NEAR(distance_between(a, a), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(distance_squared(a, a), 0.0f, 1e-6f);
}

inline void test_angle_between_vectors() {
  Vector x(1, 0, 0), y(0, 1, 0);
  HS_EXPECT_NEAR(angle_between(x, y), PI_F * 0.5f, 1e-3f);
  HS_EXPECT_NEAR(angle_between(x, x), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(angle_between(x, -x), PI_F, 1e-3f);
  // Independent of magnitude
  HS_EXPECT_NEAR(angle_between(x * 5.0f, y * 0.3f), PI_F * 0.5f, 1e-3f);
  // A zero-length operand now traps (HS_CHECK) rather than returning 0, so it
  // is no longer exercised here — see normalized_or() for the guarded path.
}

// ============================================================================
// Spherical
// ============================================================================

inline void test_spherical() {
  Spherical s(0.5f, 1.2f);
  HS_EXPECT_NEAR(s.theta, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(s.phi, 1.2f, 1e-6f);

  // Roundtrip on a non-pole unit vector
  Vector v_orig(0.6f, 0.0f, 0.8f);
  Spherical s2(v_orig);
  Vector v2(s2);
  HS_EXPECT_VEC(v2, v_orig, 5e-3f);

  // Roundtrip on an off-equator unit vector
  Vector v3 = Vector(1.0f, 1.0f, 1.0f).normalized();
  Spherical s3(v3);
  Vector v3_back(s3);
  HS_EXPECT_VEC(v3_back, v3, 5e-3f);
}

// ============================================================================
// Quaternion
// ============================================================================

inline void test_quaternion_construction() {
  Quaternion id;
  HS_EXPECT_NEAR(id.r, 1.0f, 1e-7f);
  HS_EXPECT_VEC(id.v, Vector(0, 0, 0), 1e-7f);

  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);
  HS_EXPECT_NEAR(q.r, 0.5f, 0.0f);
  HS_EXPECT_VEC(q.v, Vector(0.5f, 0.5f, 0.5f), 0.0f);

  Quaternion qv(0.7f, Vector(0.1f, 0.2f, 0.3f));
  HS_EXPECT_NEAR(qv.r, 0.7f, 0.0f);
  HS_EXPECT_VEC(qv.v, Vector(0.1f, 0.2f, 0.3f), 0.0f);

  Quaternion qc(q);
  HS_EXPECT_QUAT(qc, q, 1e-7f);

  Quaternion qa;
  qa = q;
  HS_EXPECT_QUAT(qa, q, 1e-7f);
}

inline void test_quaternion_arithmetic() {
  Quaternion a(1, 2, 3, 4), b(0.5f, 1, 1.5f, 2);
  HS_EXPECT_QUAT(a + b, Quaternion(1.5f, 3, 4.5f, 6), 1e-6f);
  HS_EXPECT_QUAT(a - b, Quaternion(0.5f, 1, 1.5f, 2), 1e-6f);
  HS_EXPECT_QUAT(a * 2.0f, Quaternion(2, 4, 6, 8), 1e-6f);
  HS_EXPECT_QUAT(2.0f * a, Quaternion(2, 4, 6, 8), 1e-6f);
  HS_EXPECT_QUAT(a / 2.0f, Quaternion(0.5f, 1, 1.5f, 2), 1e-6f);
  HS_EXPECT_QUAT(-a, Quaternion(-1, -2, -3, -4), 1e-6f);

  Quaternion c(1, 2, 3, 4);
  c += b;
  HS_EXPECT_QUAT(c, Quaternion(1.5f, 3, 4.5f, 6), 1e-6f);
  c -= b;
  HS_EXPECT_QUAT(c, Quaternion(1, 2, 3, 4), 1e-6f);
  c *= 0.5f;
  HS_EXPECT_QUAT(c, Quaternion(0.5f, 1, 1.5f, 2), 1e-6f);
}

inline void test_quaternion_magnitude() {
  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);
  HS_EXPECT_NEAR(q.squared_magnitude(), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(q.magnitude(), 1.0f, 1e-6f);

  Quaternion p(1, 2, 2, 0);
  HS_EXPECT_NEAR(p.squared_magnitude(), 9.0f, 1e-5f);
  HS_EXPECT_NEAR(p.magnitude(), 3.0f, 1e-6f);
}

inline void test_quaternion_conjugate_inverse() {
  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f); // unit
  Quaternion conj = q.conjugate();
  HS_EXPECT_QUAT(conj, Quaternion(0.5f, -0.5f, -0.5f, -0.5f), 1e-6f);

  // For unit q, inverse and unit_inverse both equal conjugate
  HS_EXPECT_QUAT(q.inverse(), conj, 1e-6f);
  HS_EXPECT_QUAT(q.unit_inverse(), conj, 1e-6f);

  // q * q.inverse() == identity
  HS_EXPECT_QUAT(q * q.inverse(), Quaternion(1, 0, 0, 0), 1e-6f);

  // Non-unit quaternion: inverse() still gives a multiplicative inverse
  Quaternion p(2, 0, 0, 0); // |p|² = 4
  HS_EXPECT_QUAT(p * p.inverse(), Quaternion(1, 0, 0, 0), 1e-6f);
}

inline void test_quaternion_normalize() {
  Quaternion p(2, 0, 0, 0);
  p.normalize();
  HS_EXPECT_NEAR(p.magnitude(), 1.0f, 1e-6f);
  HS_EXPECT_QUAT(p, Quaternion(1, 0, 0, 0), 1e-6f);

  // A zero-magnitude quaternion now traps (HS_CHECK) rather than falling back.

  // normalized() does not mutate
  Quaternion u(3, 0, 0, 0);
  Quaternion n = u.normalized();
  HS_EXPECT_NEAR(u.r, 3.0f, 1e-6f);
  HS_EXPECT_QUAT(n, Quaternion(1, 0, 0, 0), 1e-6f);
}

inline void test_quaternion_multiplication() {
  Quaternion id;
  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);

  // Identity laws
  HS_EXPECT_QUAT(id * q, q, 1e-6f);
  HS_EXPECT_QUAT(q * id, q, 1e-6f);

  // Hamilton basis: i² = j² = k² = ijk = -1
  Quaternion i(0, 1, 0, 0), j(0, 0, 1, 0), k(0, 0, 0, 1);
  Quaternion neg_one(-1, 0, 0, 0);
  HS_EXPECT_QUAT(i * i, neg_one, 1e-6f);
  HS_EXPECT_QUAT(j * j, neg_one, 1e-6f);
  HS_EXPECT_QUAT(k * k, neg_one, 1e-6f);
  HS_EXPECT_QUAT(i * j * k, neg_one, 1e-6f);

  // Cyclic products
  HS_EXPECT_QUAT(i * j, k, 1e-6f);
  HS_EXPECT_QUAT(j * k, i, 1e-6f);
  HS_EXPECT_QUAT(k * i, j, 1e-6f);
  // Non-commutativity
  HS_EXPECT_QUAT(j * i, -k, 1e-6f);

  // *= consistency
  Quaternion qa(0.5f, 0.5f, 0.5f, 0.5f), qb(qa);
  qa *= qa;
  HS_EXPECT_QUAT(qa, qb * qb, 1e-6f);
}

inline void test_quaternion_equality() {
  Quaternion a(1, 2, 3, 4), b(1, 2, 3, 4);
  HS_EXPECT_TRUE(a == b);
  Quaternion c(1.001f, 2, 3, 4); // beyond TOLERANCE
  HS_EXPECT_FALSE(a == c);
  Quaternion d(1.00001f, 2, 3, 4); // within TOLERANCE
  HS_EXPECT_TRUE(a == d);
}

inline void test_dot_quaternion() {
  Quaternion a(1, 2, 3, 4), b(2, 3, 4, 5);
  HS_EXPECT_NEAR(dot(a, b), 1 * 2 + 2 * 3 + 3 * 4 + 4 * 5, 1e-5f);
  // Self-dot equals squared magnitude
  HS_EXPECT_NEAR(dot(a, a), a.squared_magnitude(), 1e-5f);
}

inline void test_angle_between_quaternions() {
  Quaternion id(1, 0, 0, 0);
  HS_EXPECT_NEAR(angle_between(id, id), 0.0f, 1e-4f);
  Quaternion qx(0, 1, 0, 0); // unit pure quaternion
  HS_EXPECT_NEAR(angle_between(id, qx), PI_F * 0.5f, 1e-4f);
}

// ============================================================================
// make_rotation / rotate
// ============================================================================

inline void test_make_rotation_axis_angle() {
  // Identity (angle = 0): real part = ±1, vector = 0
  Quaternion id = make_rotation(Vector(0, 1, 0), 0.0f);
  HS_EXPECT_NEAR(std::abs(id.r), 1.0f, 5e-3f);
  HS_EXPECT_VEC(id.v, Vector(0, 0, 0), 5e-3f);

  // Resulting quaternion is unit-length
  Quaternion qy90 = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  HS_EXPECT_NEAR(qy90.magnitude(), 1.0f, 1e-4f);

  // 90° around +Y rotates (1,0,0) to (0,0,-1) [right-handed]
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), qy90), Vector(0, 0, -1), 5e-3f);

  // 180° around +Z rotates (1,0,0) to (-1,0,0)
  Quaternion qz180 = make_rotation(Vector(0, 0, 1), PI_F);
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), qz180), Vector(-1, 0, 0), 5e-3f);
}

inline void test_make_rotation_from_to() {
  // Parallel → identity
  Quaternion id = make_rotation(Vector(1, 0, 0), Vector(1, 0, 0));
  HS_EXPECT_QUAT(id, Quaternion(1, 0, 0, 0), 1e-4f);

  // Perpendicular: x to y
  Quaternion q = make_rotation(Vector(1, 0, 0), Vector(0, 1, 0));
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), q), Vector(0, 1, 0), 5e-3f);

  // Antiparallel (degenerate): x to -x  → 180° rotation
  Quaternion qa = make_rotation(Vector(1, 0, 0), Vector(-1, 0, 0));
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), qa), Vector(-1, 0, 0), 5e-3f);

  // Generic case: rotation from one direction to another actually achieves it
  Vector from(1, 1, 0);
  from.normalize();
  Vector to(0, 1, 1);
  to.normalize();
  Quaternion qg = make_rotation(from, to);
  HS_EXPECT_VEC(rotate(from, qg), to, 5e-3f);
  HS_EXPECT_NEAR(qg.magnitude(), 1.0f, 1e-3f);
}

inline void test_rotate() {
  // Identity rotation leaves vectors unchanged
  Vector v(1, 2, 3);
  HS_EXPECT_VEC(rotate(v, Quaternion(1, 0, 0, 0)), v, 1e-6f);

  // Rotation preserves length
  Quaternion q =
      make_rotation(Vector(1, 2, 3).normalized(), 1.234f);
  Vector r = rotate(v, q);
  HS_EXPECT_NEAR(r.length(), v.length(), 1e-3f);

  // Composition: rotate(rotate(v, q1), q2) == rotate(v, q2 * q1)
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.3f);
  Quaternion q2 = make_rotation(Vector(1, 0, 0), 0.5f);
  Vector via_sequential = rotate(rotate(v, q1), q2);
  Vector via_composed = rotate(v, q2 * q1);
  HS_EXPECT_VEC(via_sequential, via_composed, 5e-3f);
}

// ============================================================================
// slerp
// ============================================================================

inline void test_vector_slerp() {
  Vector a(1, 0, 0), b(0, 1, 0);
  HS_EXPECT_VEC(slerp(a, b, 0.0f), a, 5e-3f);
  HS_EXPECT_VEC(slerp(a, b, 1.0f), b, 5e-3f);

  // Midpoint on unit sphere is (√2/2, √2/2, 0)
  Vector mid = slerp(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(mid.x, std::sqrt(2.0f) * 0.5f, 5e-3f);
  HS_EXPECT_NEAR(mid.y, std::sqrt(2.0f) * 0.5f, 5e-3f);
  HS_EXPECT_NEAR(mid.z, 0.0f, 5e-3f);

  // Nearly-identical vectors: falls back to lerp; result still unit-length
  Vector v1(1, 0, 0);
  Vector v2(0.99999f, 0.00001f, 0.0f);
  Vector lerp_result = slerp(v1, v2, 0.5f);
  HS_EXPECT_NEAR(lerp_result.length(), 1.0f, 1e-3f);

  // Antipodal endpoints: the lerp midpoint collapses to the zero vector at
  // t=0.5; slerp must degrade to a stable unit direction, not trap in
  // normalized(). Endpoints stay exact; the midpoint stays unit-length.
  Vector p(0, 1, 0), ap(0, -1, 0);
  HS_EXPECT_VEC(slerp(p, ap, 0.0f), p, 5e-3f);
  HS_EXPECT_VEC(slerp(p, ap, 1.0f), ap, 5e-3f);
  Vector anti_mid = slerp(p, ap, 0.5f);
  HS_EXPECT_NEAR(anti_mid.length(), 1.0f, 1e-3f);
}

inline void test_quaternion_slerp() {
  Quaternion id(1, 0, 0, 0);
  Quaternion q = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);

  // Endpoint at t=0 returns the starting quaternion
  HS_EXPECT_QUAT(slerp(id, q, 0.0f), id, 5e-3f);

  // Endpoint at t=1 returns q (possibly negated — same rotation).
  // Test via |dot| since q and -q represent the same orientation.
  Quaternion s1 = slerp(id, q, 1.0f);
  HS_EXPECT_NEAR(std::abs(dot(s1, q)), 1.0f, 5e-3f);

  // Halfway slerp is a unit quaternion
  Quaternion half = slerp(id, q, 0.5f);
  HS_EXPECT_NEAR(half.magnitude(), 1.0f, 1e-3f);

  // Composing half with itself recovers q's rotation (q^0.5 squared = q)
  Vector v(1, 0, 0);
  Vector r_twice = rotate(rotate(v, half), half);
  Vector r_full = rotate(v, q);
  HS_EXPECT_VEC(r_twice, r_full, 1e-2f);

  // long_way variant: the algorithm negates the start, so the long-arc path
  // begins at -q1. At t=0 it must therefore return -id (and the long-way
  // midpoint must differ from the short-way midpoint).
  Quaternion long_start = slerp(id, q, 0.0f, true);
  HS_EXPECT_QUAT(long_start, -id, 5e-3f);

  Quaternion mid_short = slerp(id, q, 0.5f, false);
  Quaternion mid_long = slerp(id, q, 0.5f, true);
  HS_EXPECT_FALSE(approx_quat(mid_short, mid_long, 1e-2f));
}

// ============================================================================
// fold_to_hemisphere
// ============================================================================

inline void test_fold_to_hemisphere() {
  // North pole maps to itself
  HS_EXPECT_VEC(fold_to_hemisphere(Vector(0, 1, 0)), Vector(0, 1, 0), 1e-4f);
  // South pole anchors to +X
  HS_EXPECT_VEC(fold_to_hemisphere(Vector(0, -1, 0)), Vector(1, 0, 0), 1e-4f);

  // Folded northern-hemisphere point lies on the unit sphere
  Vector vN = Vector(0.5f, 0.3f, 0.81f).normalized();
  Vector fN = fold_to_hemisphere(vN);
  HS_EXPECT_NEAR(fN.length(), 1.0f, 5e-3f);
  HS_EXPECT_TRUE(fN.y >= -1e-4f);

  // Folded southern-hemisphere point also stays unit & in the upper half
  Vector vS = Vector(0.5f, -0.5f, 0.7071f).normalized();
  Vector fS = fold_to_hemisphere(vS);
  HS_EXPECT_NEAR(fS.length(), 1.0f, 5e-3f);
  HS_EXPECT_TRUE(fS.y >= -1e-4f);
}

// ============================================================================
// Stereographic projection
// ============================================================================

inline void test_stereo_roundtrip() {
  Vector samples[] = {
      Vector(1, 0, 0),
      Vector(0, 0, 1),
      Vector(-1, 0, 0),
      Vector(0, -1, 0),
      Vector(0.6f, 0.0f, 0.8f),
      Vector(0.5f, 0.5f, 0.7071f).normalized(),
  };
  for (const Vector &v : samples) {
    Complex z = stereo(v);
    Vector back = inv_stereo(z);
    HS_EXPECT_VEC(back, v, 5e-3f);
  }

  // North pole maps to the infinity sentinel (azimuth undefined → +real axis)
  Complex zN = stereo(Vector(0, 1, 0));
  HS_EXPECT_NEAR(zN.re, STEREO_INF, 1.0f);
  HS_EXPECT_NEAR(zN.im, 0.0f, 1.0f);

  // Inside the pole cap (denom < STEREO_POLE_EPS) the sentinel preserves the
  // (x,z) azimuth at magnitude STEREO_INF rather than collapsing onto +real, so
  // a Mobius map can carry the cap's swirl to a finite point.
  Vector nearPole(0.006f, 0.99998f, 0.0021f); // |xz| small, denom ≈ 2e-5
  Complex zCap = stereo(nearPole.normalized());
  HS_EXPECT_NEAR(std::sqrt(zCap.re * zCap.re + zCap.im * zCap.im), STEREO_INF,
                 1.0f);
  HS_EXPECT_NEAR(std::atan2(zCap.im, zCap.re),
                 std::atan2(nearPole.z, nearPole.x), 1e-3f);

  // Large complex magnitude maps back to north pole
  Vector pole = inv_stereo(Complex(STEREO_INF, 0));
  HS_EXPECT_VEC(pole, Vector(0, 1, 0), 1e-3f);

  // Origin in plane maps to south pole
  HS_EXPECT_VEC(inv_stereo(Complex(0, 0)), Vector(0, -1, 0), 1e-3f);
}

// ============================================================================
// Complex
// ============================================================================

inline void test_complex_arithmetic() {
  Complex a(1, 2), b(3, 4);

  HS_EXPECT_COMPLEX(a + b, Complex(4, 6), 1e-6f);
  HS_EXPECT_COMPLEX(a - b, Complex(-2, -2), 1e-6f);
  // (1+2i)(3+4i) = 3 + 4i + 6i + 8i² = -5 + 10i
  HS_EXPECT_COMPLEX(a * b, Complex(-5, 10), 1e-5f);

  // Division round-trips: (a / b) * b ≈ a
  Complex q = a / b;
  HS_EXPECT_COMPLEX(q * b, a, 1e-4f);

  // Multiplicative identity
  HS_EXPECT_COMPLEX(a * Complex(1, 0), a, 1e-6f);

  // 0 / 0 → indeterminate, returns 0 by convention
  HS_EXPECT_COMPLEX(Complex(0, 0) / Complex(0, 0), Complex(0, 0), 1e-6f);

  // Nonzero / 0 → large magnitude in numerator direction
  Complex inf_dir = Complex(1, 0) / Complex(0, 0);
  HS_EXPECT_TRUE(std::abs(inf_dir.re) > 1e3f);
}

// ============================================================================
// Mobius
// ============================================================================

inline void test_mobius_params_accessors() {
  // The eight-float constructor fills the four Complex coefficients in order.
  MobiusParams p(1, 2, 3, 4, 5, 6, 7, 8);
  HS_EXPECT_COMPLEX(p.a, Complex(1, 2), 0.0f);
  HS_EXPECT_COMPLEX(p.b, Complex(3, 4), 0.0f);
  HS_EXPECT_COMPLEX(p.c, Complex(5, 6), 0.0f);
  HS_EXPECT_COMPLEX(p.d, Complex(7, 8), 0.0f);

  // The Complex-tuple constructor stores the coefficients as given.
  MobiusParams q(Complex(1, 2), Complex(3, 4), Complex(5, 6), Complex(7, 8));
  HS_EXPECT_COMPLEX(q.a, Complex(1, 2), 0.0f);
  HS_EXPECT_COMPLEX(q.d, Complex(7, 8), 0.0f);

  // Default constructor builds the identity Mobius (a=d=1, b=c=0)
  MobiusParams id;
  HS_EXPECT_COMPLEX(id.a, Complex(1, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.b, Complex(0, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.c, Complex(0, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.d, Complex(1, 0), 0.0f);
}

inline void test_mobius_transform() {
  Complex z(0.3f, 0.7f);

  // Identity Mobius leaves points unchanged
  MobiusParams id;
  HS_EXPECT_COMPLEX(mobius(z, id), z, 1e-4f);

  // Pure translation: (1·z + (2 - i)) / (0·z + 1) = z + (2 - i)
  MobiusParams trans(1, 0, 2.0f, -1.0f, 0, 0, 1, 0);
  HS_EXPECT_COMPLEX(mobius(z, trans), Complex(z.re + 2.0f, z.im - 1.0f), 1e-4f);

  // Pure scaling: (3·z) / 1 = 3z
  MobiusParams scl(3, 0, 0, 0, 0, 0, 1, 0);
  HS_EXPECT_COMPLEX(mobius(z, scl), Complex(z.re * 3.0f, z.im * 3.0f), 1e-4f);
}

// ============================================================================
// Gnomonic
// ============================================================================

inline void test_gnomonic_roundtrip() {
  // Generic upper-hemisphere point
  Vector vUp = Vector(0.3f, 0.8f, 0.4f).normalized();
  Complex zUp = gnomonic(vUp);
  Vector vUp_back = inv_gnomonic(zUp, 1.0f);
  HS_EXPECT_VEC(vUp_back, vUp, 5e-3f);

  // Generic lower-hemisphere point — sign tracked explicitly
  Vector vDn = Vector(0.3f, -0.8f, 0.4f).normalized();
  Complex zDn = gnomonic(vDn);
  Vector vDn_back = inv_gnomonic(zDn, -1.0f);
  HS_EXPECT_VEC(vDn_back, vDn, 5e-3f);

  // North-pole pre-image: gnomonic(0,1,0) = (0, 0); inv → (0, 1, 0)
  HS_EXPECT_COMPLEX(gnomonic(Vector(0, 1, 0)), Complex(0, 0), 1e-4f);
  HS_EXPECT_VEC(inv_gnomonic(Complex(0, 0), 1.0f), Vector(0, 1, 0), 1e-6f);

  // Saturated input → returns the pole, sign-dependent
  HS_EXPECT_VEC(inv_gnomonic(Complex(STEREO_INF, 0), 1.0f), Vector(0, 1, 0),
                1e-3f);
  HS_EXPECT_VEC(inv_gnomonic(Complex(STEREO_INF, 0), -1.0f), Vector(0, -1, 0),
                1e-3f);

  // Near-equator inputs get clamped to STEREO_INF
  Complex zEq = gnomonic(Vector(1.0f, 1e-10f, 0.0f));
  HS_EXPECT_TRUE(std::abs(zEq.re) >= STEREO_INF - 1.0f);
}

// ============================================================================
// Spline
// ============================================================================

inline void test_spline_cubic_endpoints() {
  const Vector p0(1, 0, 0);
  const Vector p1 = Vector(0.7f, 0.7f, 0).normalized();
  const Vector p2(0, 1, 0);
  const Vector p3 = Vector(-0.7f, 0.7f, 0).normalized();

  // cubic_fast hits the (normalized) endpoints
  HS_EXPECT_VEC(Spline::cubic_fast(p0, p1, p2, p3, 0.0f), p0, 1e-3f);
  HS_EXPECT_VEC(Spline::cubic_fast(p0, p1, p2, p3, 1.0f), p3, 1e-3f);

  // cubic_slerp hits the endpoints
  HS_EXPECT_VEC(Spline::cubic_slerp(p0, p1, p2, p3, 0.0f), p0, 5e-3f);
  HS_EXPECT_VEC(Spline::cubic_slerp(p0, p1, p2, p3, 1.0f), p3, 5e-3f);

  // Both variants yield unit-length results at midpoint
  Vector mFast = Spline::cubic_fast(p0, p1, p2, p3, 0.5f);
  Vector mSlerp = Spline::cubic_slerp(p0, p1, p2, p3, 0.5f);
  HS_EXPECT_NEAR(mFast.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(mSlerp.length(), 1.0f, 1e-3f);

  // Dispatch routes to the appropriate variant
  HS_EXPECT_VEC(Spline::cubic(p0, p1, p2, p3, 0.5f, SplineMode::Fast), mFast,
                1e-4f);
  HS_EXPECT_VEC(Spline::cubic(p0, p1, p2, p3, 0.5f, SplineMode::Geodesic),
                mSlerp, 1e-4f);
}

inline void test_spline_catmull_rom() {
  const Vector prev(1, 0, 0);
  Vector start(0.7f, 0.7f, 0);
  start.normalize();
  const Vector end(0, 1, 0);
  Vector next(-0.7f, 0.7f, 0);
  next.normalize();

  Vector cp1, cp2;

  // tension = 0 → geodesic segment: cp1 = start, cp2 = end
  Spline::catmull_rom_tangents(prev, start, end, next, 0.0f, cp1, cp2);
  HS_EXPECT_VEC(cp1, start, 1e-3f);
  HS_EXPECT_VEC(cp2, end, 1e-3f);

  // tension = 1 → full smoothing: cp1 = midpoint(prev, end),
  // cp2 = midpoint(start, next)
  Spline::catmull_rom_tangents(prev, start, end, next, 1.0f, cp1, cp2);
  HS_EXPECT_VEC(cp1, slerp(prev, end, 0.5f), 1e-3f);
  HS_EXPECT_VEC(cp2, slerp(start, next, 0.5f), 1e-3f);

  // Output control points are unit-length
  Spline::catmull_rom_tangents(prev, start, end, next, 0.5f, cp1, cp2);
  HS_EXPECT_NEAR(cp1.length(), 1.0f, 5e-3f);
  HS_EXPECT_NEAR(cp2.length(), 1.0f, 5e-3f);
}

// ============================================================================
// Runner
// ============================================================================

// ============================================================================
// wrap_index (core/rotate.h) — folds a float index into [0, m)
// ============================================================================

inline void test_wrap_index() {
  const int m = 288;

  // Non-negative inputs: integer + fractional parts preserved.
  HS_EXPECT_NEAR(wrap_index(0.0f, m), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(wrap_index(0.5f, m), 0.5f, 1e-5f);
  HS_EXPECT_NEAR(wrap_index(287.9f, m), 287.9f, 1e-3f);

  // Wrap at and above the period.
  HS_EXPECT_NEAR(wrap_index(static_cast<float>(m), m), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(wrap_index(m + 1.5f, m), 1.5f, 1e-4f);

  // Negative inputs must fold into [0, m). Regression: a truncating cast
  // returned them unchanged (e.g. -0.5 -> -0.5), yielding a negative pixel x.
  HS_EXPECT_NEAR(wrap_index(-0.5f, m), 287.5f, 1e-3f);
  HS_EXPECT_NEAR(wrap_index(-1.5f, m), 286.5f, 1e-3f);
  HS_EXPECT_NEAR(wrap_index(-static_cast<float>(m) + 0.25f, m), 0.25f, 1e-3f);

  // Result stays in [0, m) across several periods of both signs.
  for (int i = -3 * m; i <= 3 * m; ++i) {
    float w = wrap_index(i * 0.5f, m);
    HS_EXPECT_TRUE(w >= 0.0f && w < static_cast<float>(m));
  }
}

inline int run_3dmath_tests() {
  auto scope = hs_test::begin_module("3dmath");

  test_constants();
  test_quintic_kernel();

  test_fast_atan2();
  test_fast_acos();
  test_fast_sinf_cosf();

  test_vector_construction();
  test_vector_spherical_construction();
  test_vector_equality();
  test_vector_arithmetic();
  test_vector_length();
  test_vector_normalize();

  test_dot_cross();
  test_distance();
  test_angle_between_vectors();

  test_spherical();

  test_quaternion_construction();
  test_quaternion_arithmetic();
  test_quaternion_magnitude();
  test_quaternion_conjugate_inverse();
  test_quaternion_normalize();
  test_quaternion_multiplication();
  test_quaternion_equality();
  test_dot_quaternion();
  test_angle_between_quaternions();

  test_make_rotation_axis_angle();
  test_make_rotation_from_to();
  test_rotate();

  test_vector_slerp();
  test_quaternion_slerp();

  test_fold_to_hemisphere();
  test_stereo_roundtrip();
  test_complex_arithmetic();
  test_mobius_params_accessors();
  test_mobius_transform();
  test_gnomonic_roundtrip();

  test_spline_cubic_endpoints();
  test_spline_catmull_rom();

  test_wrap_index();

  return hs_test::end_module(scope);
}

} // namespace math3d
} // namespace hs_test

