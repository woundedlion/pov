/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Unit tests for core/math/3dmath.h.
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

#include "core/math/3dmath.h"
#include "core/math/rotate.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace math3d {

/**
 * @brief Tests whether two vectors agree componentwise within a tolerance.
 * @param a First vector operand.
 * @param b Second vector operand.
 * @param tol Per-component absolute tolerance.
 * @return True if every component of a and b agrees within tol.
 */
inline bool approx_vec(const Vector &a, const Vector &b, float tol) {
  return approx(a.x, b.x, tol) && approx(a.y, b.y, tol) &&
         approx(a.z, b.z, tol);
}
/**
 * @brief Tests whether two quaternions agree within a tolerance.
 * @param a First quaternion operand.
 * @param b Second quaternion operand.
 * @param tol Per-component absolute tolerance.
 * @return True if the scalar and vector parts of a and b agree within tol.
 */
inline bool approx_quat(const Quaternion &a, const Quaternion &b, float tol) {
  return approx(a.r, b.r, tol) && approx_vec(a.v, b.v, tol);
}
/**
 * @brief Tests whether two complex numbers agree within a tolerance.
 * @param a First complex operand.
 * @param b Second complex operand.
 * @param tol Per-component absolute tolerance.
 * @return True if the real and imaginary parts of a and b agree within tol.
 */
inline bool approx_complex(const Complex &a, const Complex &b, float tol) {
  return approx(a.re, b.re, tol) && approx(a.im, b.im, tol);
}

/**
 * @brief Tolerant equality assertions for vectors, quaternions, and complex
 *        values; the failure message stringizes the two compared expressions.
 */
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

/**
 * @brief Pins the math constants (golden ratio, tolerance, pi, stereo
 *        sentinel) to their expected values.
 */
inline void test_constants() {
  HS_EXPECT_NEAR(PHI, 1.61803398f, 1e-6f);
  HS_EXPECT_NEAR(INV_PHI, 1.0f / PHI, 1e-6f);
  HS_EXPECT_NEAR(INV_PHI * PHI, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::TOLERANCE, 0.0001f, 1e-9f);
  HS_EXPECT_NEAR(TOLERANCE, math::TOLERANCE, 1e-9f);
  HS_EXPECT_NEAR(PI_F, 3.14159265f, 1e-5f);
  HS_EXPECT_NEAR(STEREO_INF, 1e4f, 1e-3f);
}

// ============================================================================
// quintic_kernel
// ============================================================================

/**
 * @brief Verifies the smootherstep kernel: fixed points, clamping outside
 *        [0,1], monotonicity, and flat (C2) endpoints.
 */
inline void test_quintic_kernel() {
  HS_EXPECT_NEAR(quintic_kernel(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(quintic_kernel(1.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(quintic_kernel(0.5f), 0.5f, 1e-6f);

  HS_EXPECT_NEAR(quintic_kernel(-1.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(quintic_kernel(2.0f), 1.0f, 1e-6f);

  float prev = -1.0f;
  for (int i = 0; i <= 200; ++i) {
    float v = quintic_kernel(i / 200.0f);
    HS_EXPECT_TRUE(v >= prev);
    prev = v;
  }

  // C2-continuous: derivative ~= 0 at endpoints.
  float dl = quintic_kernel(0.01f) - quintic_kernel(0.0f);
  float dr = quintic_kernel(1.0f) - quintic_kernel(0.99f);
  HS_EXPECT_TRUE(std::abs(dl) < 1e-3f);
  HS_EXPECT_TRUE(std::abs(dr) < 1e-3f);
}

// ============================================================================
// fast_atan2 / fast_acos / fast_sinf / fast_cosf
// (measured peak errors: atan2 ~3.8e-3 rad, acos ~1.3e-4 rad, sin ~1.6e-3)
// ============================================================================

/**
 * @brief Verifies fast_atan2 at cardinal directions and across a full-circle
 *        sweep against std::atan2.
 */
inline void test_fast_atan2() {
  HS_EXPECT_NEAR(fast_atan2(0.0f, 1.0f), 0.0f, 4e-3f);
  HS_EXPECT_NEAR(fast_atan2(1.0f, 0.0f), PI_F * 0.5f, 4e-3f);
  HS_EXPECT_NEAR(fast_atan2(0.0f, -1.0f), PI_F, 4e-3f);
  HS_EXPECT_NEAR(fast_atan2(-1.0f, 0.0f), -PI_F * 0.5f, 4e-3f);

  for (int i = 0; i < 64; ++i) {
    float a = -PI_F + (i * 2.0f * PI_F) / 64.0f;
    float y = std::sin(a);
    float x = std::cos(a);
    HS_EXPECT_NEAR(fast_atan2(y, x), std::atan2(y, x), 4e-3f);
  }

  // Peak error (~3.76e-3) near a = -2.5702 rad, between the sweep's samples.
  {
    float a = -2.5702f, y = std::sin(a), x = std::cos(a);
    HS_EXPECT_NEAR(fast_atan2(y, x), std::atan2(y, x), 4e-3f);
  }
}

/**
 * @brief Verifies fast_acos at endpoints, out-of-range clamping to [0,pi], and
 *        across a sweep against std::acos.
 */
inline void test_fast_acos() {
  HS_EXPECT_NEAR(fast_acos(1.0f), 0.0f, 2e-4f);
  HS_EXPECT_NEAR(fast_acos(-1.0f), PI_F, 2e-4f);
  HS_EXPECT_NEAR(fast_acos(0.0f), PI_F * 0.5f, 2e-4f);

  HS_EXPECT_NEAR(fast_acos(1.5f), 0.0f, 2e-4f);
  HS_EXPECT_NEAR(fast_acos(-1.5f), PI_F, 2e-4f);

  for (int i = 0; i <= 32; ++i) {
    float x = -1.0f + (i / 16.0f);
    if (x > 1.0f) x = 1.0f;
    HS_EXPECT_NEAR(fast_acos(x), std::acos(x), 2e-4f);
  }

  // Peak error (~1.26e-4) near x = 0.1226, between the sweep's samples.
  HS_EXPECT_NEAR(fast_acos(0.1226f), std::acos(0.1226f), 2e-4f);
}

/**
 * @brief Verifies fast_cbrt anchors, the x<=0 -> 0 clamp, and the documented
 *        ~2.3e-5 peak relative error over [0,8] (plus a few values past 8).
 */
inline void test_fast_cbrt() {
  HS_EXPECT_NEAR(fast_cbrt(-1.0f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(fast_cbrt(0.0f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(fast_cbrt(1.0f), 1.0f, 2.3e-5f);
  HS_EXPECT_NEAR(fast_cbrt(8.0f), 2.0f, 2.3e-5f * 2.0f);

  for (int i = 1; i <= 256; ++i) {
    float x = (8.0f * i) / 256.0f;
    float rel = std::abs(fast_cbrt(x) - std::cbrt(x)) / std::cbrt(x);
    HS_EXPECT_TRUE(rel <= 2.3e-5f);
  }

  // Values past the documented [0,8] domain stay within the same rel error.
  for (float x : {27.0f, 100.0f, 1000.0f}) {
    float rel = std::abs(fast_cbrt(x) - std::cbrt(x)) / std::cbrt(x);
    HS_EXPECT_TRUE(rel <= 2.3e-5f);
  }
}

/**
 * @brief Verifies fast_sinf/fast_cosf at key angles, the Pythagorean identity,
 *        and periodicity.
 * @details The periodicity check exercises range reduction beyond ±2π.
 */
inline void test_fast_sinf_cosf() {
  HS_EXPECT_NEAR(fast_sinf(0.0f), 0.0f, 1.8e-3f);
  HS_EXPECT_NEAR(fast_sinf(PI_F * 0.5f), 1.0f, 1.8e-3f);
  HS_EXPECT_NEAR(fast_sinf(PI_F), 0.0f, 1.8e-3f);
  HS_EXPECT_NEAR(fast_sinf(-PI_F * 0.5f), -1.0f, 1.8e-3f);

  HS_EXPECT_NEAR(fast_cosf(0.0f), 1.0f, 1.8e-3f);
  HS_EXPECT_NEAR(fast_cosf(PI_F * 0.5f), 0.0f, 1.8e-3f);
  HS_EXPECT_NEAR(fast_cosf(PI_F), -1.0f, 1.8e-3f);

  // Sweep peak (~1.63e-3) sits near a = -9.22, in the range-reduction band.
  for (int i = 0; i <= 256; ++i) {
    float a = -3.0f * PI_F + (i * 6.0f * PI_F) / 256.0f;
    HS_EXPECT_NEAR(fast_sinf(a), std::sin(a), 1.8e-3f);
    HS_EXPECT_NEAR(fast_cosf(a), std::cos(a), 1.8e-3f);
  }

  for (int i = 0; i < 32; ++i) {
    float a = -3.0f * PI_F + (i * 6.0f * PI_F) / 32.0f;
    float s = fast_sinf(a);
    float c = fast_cosf(a);
    HS_EXPECT_NEAR(s * s + c * c, 1.0f, 5e-3f);
  }

  for (int i = 0; i < 16; ++i) {
    float a = i * 0.3f;
    HS_EXPECT_NEAR(fast_sinf(a + 2.0f * PI_F), fast_sinf(a), 3e-3f);
    HS_EXPECT_NEAR(fast_sinf(a - 2.0f * PI_F), fast_sinf(a), 3e-3f);
  }
}

// ============================================================================
// Vector — construction
// ============================================================================

/**
 * @brief Verifies Vector constructors (default-zero, scalar, component, copy)
 *        and assignment.
 */
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

/**
 * @brief Verifies Vector::from_spherical axis mappings, poles, and unit-length
 *        preservation.
 */
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

  Vector u = Vector::from_spherical(0.8f, 1.1f);
  HS_EXPECT_NEAR(u.length(), 1.0f, 2e-3f);
}

/**
 * @brief Verifies Vector ==/!= tolerant comparison (equal within TOLERANCE,
 *        unequal beyond).
 */
inline void test_vector_equality() {
  Vector a(1, 2, 3), b(1, 2, 3), c(1.001f, 2, 3);
  HS_EXPECT_TRUE(a == b);
  HS_EXPECT_FALSE(a == c);
  HS_EXPECT_TRUE(a != c);
  HS_EXPECT_FALSE(a != b);

  HS_EXPECT_TRUE(Vector(1, 0, 0) == Vector(1.00005f, 0, 0));
  HS_EXPECT_FALSE(Vector(1, 0, 0) == Vector(1.001f, 0, 0));
}

// ============================================================================
// Vector — arithmetic
// ============================================================================

/**
 * @brief Verifies Vector +, -, negate, scalar * and / (both orders), and
 *        compound assignment.
 */
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

/**
 * @brief Verifies Vector length()/magnitude() Euclidean norm.
 */
inline void test_vector_length() {
  HS_EXPECT_NEAR(Vector(3, 4, 0).length(), 5.0f, 1e-6f);
  HS_EXPECT_NEAR(Vector(0, 0, 0).length(), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(Vector(1, 2, 2).magnitude(), 3.0f, 1e-6f);
}

/**
 * @brief Verifies Vector normalize()/normalized() and the normalized_or()
 *        fallback.
 * @details The zero vector is rejected by the strict normalize() path, so
 *          normalized_or() supplies a fallback for the degenerate case.
 */
inline void test_vector_normalize() {
  Vector v(3, 0, 4);
  v.normalize();
  HS_EXPECT_NEAR(v.length(), 1.0f, 1e-6f);
  HS_EXPECT_VEC(v, Vector(0.6f, 0.0f, 0.8f), 1e-6f);

  Vector u(0, 5, 0);
  Vector n = u.normalized();
  HS_EXPECT_VEC(u, Vector(0, 5, 0), 1e-6f);
  HS_EXPECT_VEC(n, Vector(0, 1, 0), 1e-6f);

  // normalized_or() returns the fallback for a zero-length input (which traps
  // under strict normalize()), else normalizes as usual.
  HS_EXPECT_VEC(normalized_or(Vector(0, 0, 0), Vector(1, 0, 0)), Vector(1, 0, 0),
                1e-6f);
  HS_EXPECT_VEC(normalized_or(Vector(0, 6, 0), Vector(1, 0, 0)), Vector(0, 1, 0),
                1e-6f);
}

// ============================================================================
// Vector — free functions (dot, cross, distance, angle_between)
// ============================================================================

/**
 * @brief Verifies dot and cross: orthogonality, the right-handed basis,
 *        anticommutativity, and that a×b is perpendicular to both operands.
 */
inline void test_dot_cross() {
  Vector x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);

  HS_EXPECT_NEAR(dot(x, y), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(dot(x, x), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(dot(x, -x), -1.0f, 1e-6f);
  HS_EXPECT_NEAR(dot(Vector(1, 2, 3), Vector(4, -5, 6)), 4 - 10 + 18, 1e-5f);

  HS_EXPECT_VEC(cross(x, y), z, 1e-6f);
  HS_EXPECT_VEC(cross(y, z), x, 1e-6f);
  HS_EXPECT_VEC(cross(z, x), y, 1e-6f);
  HS_EXPECT_VEC(cross(y, x), -z, 1e-6f);
  HS_EXPECT_VEC(cross(Vector(2, 3, 5), Vector(2, 3, 5)), Vector(0, 0, 0),
                1e-6f);
  // a × b is perpendicular to both operands.
  Vector a(1, 2, 3), b(4, -5, 6);
  Vector c = cross(a, b);
  HS_EXPECT_NEAR(dot(c, a), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(dot(c, b), 0.0f, 1e-4f);
}

/**
 * @brief Verifies distance_between/distance_squared, including the
 *        coincident-point zero.
 */
inline void test_distance() {
  Vector a(1, 2, 3), b(4, 6, 3);
  HS_EXPECT_NEAR(distance_between(a, b), 5.0f, 1e-6f);
  HS_EXPECT_NEAR(distance_squared(a, b), 25.0f, 1e-5f);
  HS_EXPECT_NEAR(distance_between(a, a), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(distance_squared(a, a), 0.0f, 1e-6f);
}

/**
 * @brief Verifies angle_between for vectors: cardinal angles and
 *        magnitude-independence.
 */
inline void test_angle_between_vectors() {
  Vector x(1, 0, 0), y(0, 1, 0);
  HS_EXPECT_NEAR(angle_between(x, y), PI_F * 0.5f, 1e-3f);
  HS_EXPECT_NEAR(angle_between(x, x), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(angle_between(x, -x), PI_F, 1e-3f);
  // Independent of operand magnitude.
  HS_EXPECT_NEAR(angle_between(x * 5.0f, y * 0.3f), PI_F * 0.5f, 1e-3f);
}

// ============================================================================
// Spherical
// ============================================================================

/**
 * @brief Verifies Spherical accessors and Vector<->Spherical roundtrips (on-
 *        and off-equator).
 */
inline void test_spherical() {
  Spherical s(0.5f, 1.2f);
  HS_EXPECT_NEAR(s.theta, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(s.phi, 1.2f, 1e-6f);

  Vector v_orig(0.6f, 0.0f, 0.8f);
  Spherical s2(v_orig);
  Vector v2(s2);
  HS_EXPECT_VEC(v2, v_orig, 5e-3f);

  Vector v3 = Vector(1.0f, 1.0f, 1.0f).normalized();
  Spherical s3(v3);
  Vector v3_back(s3);
  HS_EXPECT_VEC(v3_back, v3, 5e-3f);
}

// ============================================================================
// Quaternion
// ============================================================================

/**
 * @brief Verifies Quaternion constructors (identity default, scalar+components,
 *        scalar+Vector, copy) and assignment.
 */
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

/**
 * @brief Verifies Quaternion +, -, scalar * and /, negate, and componentwise
 *        compound assignment.
 */
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

/**
 * @brief Verifies Quaternion squared_magnitude()/magnitude().
 */
inline void test_quaternion_magnitude() {
  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);
  HS_EXPECT_NEAR(q.squared_magnitude(), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(q.magnitude(), 1.0f, 1e-6f);

  Quaternion p(1, 2, 2, 0);
  HS_EXPECT_NEAR(p.squared_magnitude(), 9.0f, 1e-5f);
  HS_EXPECT_NEAR(p.magnitude(), 3.0f, 1e-6f);
}

/**
 * @brief Verifies conjugate(), inverse(), and unit_inverse().
 * @details For a unit q all three coincide; q*q^-1 = identity holds for unit
 *          and non-unit quaternions alike.
 */
inline void test_quaternion_conjugate_inverse() {
  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f); // unit
  Quaternion conj = q.conjugate();
  HS_EXPECT_QUAT(conj, Quaternion(0.5f, -0.5f, -0.5f, -0.5f), 1e-6f);

  // For unit q, inverse and unit_inverse both equal conjugate.
  HS_EXPECT_QUAT(q.inverse(), conj, 1e-6f);
  HS_EXPECT_QUAT(q.unit_inverse(), conj, 1e-6f);

  HS_EXPECT_QUAT(q * q.inverse(), Quaternion(1, 0, 0, 0), 1e-6f);

  // inverse() inverts a non-unit quaternion too.
  Quaternion p(2, 0, 0, 0);
  HS_EXPECT_QUAT(p * p.inverse(), Quaternion(1, 0, 0, 0), 1e-6f);
}

/**
 * @brief Verifies Quaternion normalize() (in place) and normalized()
 *        (non-mutating).
 */
inline void test_quaternion_normalize() {
  Quaternion p(2, 0, 0, 0);
  p.normalize();
  HS_EXPECT_NEAR(p.magnitude(), 1.0f, 1e-6f);
  HS_EXPECT_QUAT(p, Quaternion(1, 0, 0, 0), 1e-6f);

  Quaternion u(3, 0, 0, 0);
  Quaternion n = u.normalized();
  HS_EXPECT_NEAR(u.r, 3.0f, 1e-6f);
  HS_EXPECT_QUAT(n, Quaternion(1, 0, 0, 0), 1e-6f);
}

/**
 * @brief Verifies the Hamilton product: identity laws, basis relations
 *        (i²=j²=k²=ijk=-1, cyclic products), non-commutativity, and *=
 *        consistency.
 */
inline void test_quaternion_multiplication() {
  Quaternion id;
  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);

  HS_EXPECT_QUAT(id * q, q, 1e-6f);
  HS_EXPECT_QUAT(q * id, q, 1e-6f);

  // Hamilton basis: i² = j² = k² = ijk = -1
  Quaternion i(0, 1, 0, 0), j(0, 0, 1, 0), k(0, 0, 0, 1);
  Quaternion neg_one(-1, 0, 0, 0);
  HS_EXPECT_QUAT(i * i, neg_one, 1e-6f);
  HS_EXPECT_QUAT(j * j, neg_one, 1e-6f);
  HS_EXPECT_QUAT(k * k, neg_one, 1e-6f);
  HS_EXPECT_QUAT(i * j * k, neg_one, 1e-6f);

  HS_EXPECT_QUAT(i * j, k, 1e-6f);
  HS_EXPECT_QUAT(j * k, i, 1e-6f);
  HS_EXPECT_QUAT(k * i, j, 1e-6f);
  HS_EXPECT_QUAT(j * i, -k, 1e-6f);

  Quaternion qa(0.5f, 0.5f, 0.5f, 0.5f), qb(qa);
  qa *= qa;
  HS_EXPECT_QUAT(qa, qb * qb, 1e-6f);
}

/**
 * @brief Verifies Quaternion == tolerant comparison around TOLERANCE.
 */
inline void test_quaternion_equality() {
  Quaternion a(1, 2, 3, 4), b(1, 2, 3, 4);
  HS_EXPECT_TRUE(a == b);
  Quaternion c(1.001f, 2, 3, 4); // beyond TOLERANCE
  HS_EXPECT_FALSE(a == c);
  Quaternion d(1.00001f, 2, 3, 4); // within TOLERANCE
  HS_EXPECT_TRUE(a == d);
}

/**
 * @brief Verifies the 4-component dot product on quaternions; self-dot equals
 *        squared magnitude.
 */
inline void test_dot_quaternion() {
  Quaternion a(1, 2, 3, 4), b(2, 3, 4, 5);
  HS_EXPECT_NEAR(dot(a, b), 1 * 2 + 2 * 3 + 3 * 4 + 4 * 5, 1e-5f);
  HS_EXPECT_NEAR(dot(a, a), a.squared_magnitude(), 1e-5f);
}

/**
 * @brief Verifies angle_between for quaternions (the 4D angle between unit
 *        quaternions).
 */
inline void test_angle_between_quaternions() {
  Quaternion id(1, 0, 0, 0);
  HS_EXPECT_NEAR(angle_between(id, id), 0.0f, 1e-4f);
  Quaternion qx(0, 1, 0, 0); // unit pure quaternion
  HS_EXPECT_NEAR(angle_between(id, qx), PI_F * 0.5f, 1e-4f);
}

// ============================================================================
// make_rotation / rotate
// ============================================================================

/**
 * @brief Verifies make_rotation(axis, angle): identity at angle 0, unit-length
 *        result, and correct right-handed rotation of test vectors.
 */
inline void test_make_rotation_axis_angle() {
  Quaternion id = make_rotation(Vector(0, 1, 0), 0.0f);
  HS_EXPECT_NEAR(std::abs(id.r), 1.0f, 5e-3f);
  HS_EXPECT_VEC(id.v, Vector(0, 0, 0), 5e-3f);

  Quaternion qy90 = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);
  HS_EXPECT_NEAR(qy90.magnitude(), 1.0f, 1e-4f);

  // 90° around +Y rotates (1,0,0) to (0,0,-1) [right-handed]
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), qy90), Vector(0, 0, -1), 5e-3f);

  // 180° around +Z rotates (1,0,0) to (-1,0,0)
  Quaternion qz180 = make_rotation(Vector(0, 0, 1), PI_F);
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), qz180), Vector(-1, 0, 0), 5e-3f);
}

/**
 * @brief Verifies make_rotation(from, to): parallel (identity), perpendicular,
 *        the antiparallel degenerate case (180°), and a generic
 *        direction-to-direction rotation.
 */
inline void test_make_rotation_from_to() {
  // Parallel → identity
  Quaternion id = make_rotation(Vector(1, 0, 0), Vector(1, 0, 0));
  HS_EXPECT_QUAT(id, Quaternion(1, 0, 0, 0), 1e-4f);

  Quaternion q = make_rotation(Vector(1, 0, 0), Vector(0, 1, 0));
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), q), Vector(0, 1, 0), 5e-3f);

  // Antiparallel (degenerate): x to -x → 180° rotation.
  Quaternion qa = make_rotation(Vector(1, 0, 0), Vector(-1, 0, 0));
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), qa), Vector(-1, 0, 0), 5e-3f);

  Vector from(1, 1, 0);
  from.normalize();
  Vector to(0, 1, 1);
  to.normalize();
  Quaternion qg = make_rotation(from, to);
  HS_EXPECT_VEC(rotate(from, qg), to, 5e-3f);
  HS_EXPECT_NEAR(qg.magnitude(), 1.0f, 1e-3f);
}

/**
 * @brief Verifies quaternion_from_basis recovers the rotation whose columns are
 *        the given orthonormal axes, for an identity frame and a generic one.
 * @details Build an orthonormal frame by rotating the standard axes through a
 *        known quaternion, reconstruct a quaternion from that frame, and confirm
 *        it maps the body axes back onto the frame columns. Also checks the
 *        trace<=0 branch (a 180° frame) the Shepperd selection must handle.
 */
inline void test_quaternion_from_basis() {
  Quaternion id =
      quaternion_from_basis(Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1));
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), id), Vector(1, 0, 0), 1e-5f);
  HS_EXPECT_VEC(rotate(Vector(0, 0, 1), id), Vector(0, 0, 1), 1e-5f);

  Quaternion q0 = make_rotation(Vector(0.3f, 0.5f, 0.8f).normalized(), 1.1f);
  Vector cx = rotate(Vector(1, 0, 0), q0);
  Vector cy = rotate(Vector(0, 1, 0), q0);
  Vector cz = rotate(Vector(0, 0, 1), q0);
  Quaternion q = quaternion_from_basis(cx, cy, cz);
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), q), cx, 5e-3f);
  HS_EXPECT_VEC(rotate(Vector(0, 1, 0), q), cy, 5e-3f);
  HS_EXPECT_VEC(rotate(Vector(0, 0, 1), q), cz, 5e-3f);
  HS_EXPECT_NEAR(q.magnitude(), 1.0f, 1e-3f);

  // trace <= 0 branch: 180° rotation about Z (diagonal = (-1,-1,1)).
  Quaternion qz = quaternion_from_basis(Vector(-1, 0, 0), Vector(0, -1, 0),
                                        Vector(0, 0, 1));
  HS_EXPECT_VEC(rotate(Vector(1, 0, 0), qz), Vector(-1, 0, 0), 5e-3f);
  HS_EXPECT_VEC(rotate(Vector(0, 0, 1), qz), Vector(0, 0, 1), 5e-3f);
}

/**
 * @brief Verifies rotate(v, q): identity, length preservation, and the
 *        composition law rotate(rotate(v,q1),q2) == rotate(v, q2*q1).
 */
inline void test_rotate() {
  Vector v(1, 2, 3);
  HS_EXPECT_VEC(rotate(v, Quaternion(1, 0, 0, 0)), v, 1e-6f);

  Quaternion q =
      make_rotation(Vector(1, 2, 3).normalized(), 1.234f);
  Vector r = rotate(v, q);
  HS_EXPECT_NEAR(r.length(), v.length(), 1e-3f);

  // Composition: rotate(rotate(v, q1), q2) == rotate(v, q2 * q1).
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.3f);
  Quaternion q2 = make_rotation(Vector(1, 0, 0), 0.5f);
  Vector via_sequential = rotate(rotate(v, q1), q2);
  Vector via_composed = rotate(v, q2 * q1);
  HS_EXPECT_VEC(via_sequential, via_composed, 5e-3f);
}

// ============================================================================
// slerp
// ============================================================================

/**
 * @brief Verifies Vector slerp: endpoints, the unit-sphere midpoint, the
 *        near-identical lerp fallback, and the antipodal degenerate case.
 */
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

  // Nearly-identical vectors take the lerp fallback; result stays unit-length.
  Vector v1(1, 0, 0);
  Vector v2(0.99999f, 0.00001f, 0.0f);
  Vector lerp_result = slerp(v1, v2, 0.5f);
  HS_EXPECT_NEAR(lerp_result.length(), 1.0f, 1e-3f);

  // Antipodal endpoints: the great-circle direction is undefined, so slerp picks
  // a perpendicular axis and sweeps a monotone half-turn — the midpoint must NOT
  // collapse back onto p.
  Vector p(0, 1, 0), ap(0, -1, 0);
  HS_EXPECT_VEC(slerp(p, ap, 0.0f), p, 5e-3f);
  HS_EXPECT_VEC(slerp(p, ap, 1.0f), ap, 5e-3f);
  float a25 = dot(slerp(p, ap, 0.25f), p);
  float a50 = dot(slerp(p, ap, 0.50f), p);
  float a75 = dot(slerp(p, ap, 0.75f), p);
  HS_EXPECT_NEAR(slerp(p, ap, 0.5f).length(), 1.0f, 1e-3f);
  HS_EXPECT_GT(a25, a50);
  HS_EXPECT_GT(a50, a75);
  HS_EXPECT_NEAR(a50, 0.0f, 5e-3f);
}

/**
 * @brief Verifies Quaternion slerp: endpoints (q and -q are the same
 *        orientation), unit-length interpolants, the q^0.5-squared==q property,
 *        and the long_way (long-arc) variant.
 */
inline void test_quaternion_slerp() {
  Quaternion id(1, 0, 0, 0);
  Quaternion q = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);

  HS_EXPECT_QUAT(slerp(id, q, 0.0f), id, 5e-3f);

  // t=1 may return -q; |dot| == 1 since q and -q are the same orientation.
  Quaternion s1 = slerp(id, q, 1.0f);
  HS_EXPECT_NEAR(std::abs(dot(s1, q)), 1.0f, 5e-3f);

  Quaternion half = slerp(id, q, 0.5f);
  HS_EXPECT_NEAR(half.magnitude(), 1.0f, 1e-3f);

  // Composing half with itself recovers q (q^0.5 squared = q).
  Vector v(1, 0, 0);
  Vector r_twice = rotate(rotate(v, half), half);
  Vector r_full = rotate(v, q);
  HS_EXPECT_VEC(r_twice, r_full, 1e-2f);

  // long_way negates the start: t=0 returns -id, and its midpoint differs from
  // the short-arc one.
  Quaternion long_start = slerp(id, q, 0.0f, true);
  HS_EXPECT_QUAT(long_start, -id, 5e-3f);

  Quaternion mid_short = slerp(id, q, 0.5f, false);
  Quaternion mid_long = slerp(id, q, 0.5f, true);
  HS_EXPECT_FALSE(approx_quat(mid_short, mid_long, 1e-2f));
}

/**
 * @brief Verifies Vector slerp across antipodal endpoints sweeps monotonically:
 *        dot with the start strictly decreases as t increases, with no flip.
 */
inline void test_vector_slerp_antipodal_monotonic() {
  Vector p(0, 1, 0), ap(0, -1, 0);
  float prev = dot(slerp(p, ap, 0.0f), p);
  for (int i = 1; i <= 16; ++i) {
    float t = static_cast<float>(i) / 16.0f;
    Vector s = slerp(p, ap, t);
    HS_EXPECT_NEAR(s.length(), 1.0f, 1e-3f);
    float cur = dot(s, p);
    HS_EXPECT_GT(prev, cur);
    prev = cur;
  }
  HS_EXPECT_NEAR(prev, -1.0f, 5e-3f);
}

// ============================================================================
// Stereographic projection
// ============================================================================

/**
 * @brief Verifies stereo/inv_stereo roundtrips plus the pole handling.
 * @details The north pole maps to the STEREO_INF sentinel, the pole cap
 *          preserves azimuth at that magnitude, and the south pole corresponds
 *          to the plane origin.
 */
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

  // North pole maps to the infinity sentinel (azimuth undefined → +real axis).
  Complex zN = stereo(Vector(0, 1, 0));
  HS_EXPECT_NEAR(zN.re, STEREO_INF, 1.0f);
  HS_EXPECT_NEAR(zN.im, 0.0f, 1.0f);

  // Inside the pole cap (denom < STEREO_POLE_EPS) the sentinel preserves the
  // (x,z) azimuth at magnitude STEREO_INF rather than collapsing onto +real.
  Vector nearPole(0.006f, 0.99998f, 0.0021f); // |xz| small, denom ≈ 2e-5
  Complex zCap = stereo(nearPole.normalized());
  HS_EXPECT_NEAR(std::sqrt(zCap.re * zCap.re + zCap.im * zCap.im), STEREO_INF,
                 1.0f);
  HS_EXPECT_NEAR(std::atan2(zCap.im, zCap.re),
                 std::atan2(nearPole.z, nearPole.x), 1e-3f);

  // Large complex magnitude maps back to north pole.
  Vector pole = inv_stereo(Complex(STEREO_INF, 0));
  HS_EXPECT_VEC(pole, Vector(0, 1, 0), 1e-3f);

  // Plane origin maps to south pole.
  HS_EXPECT_VEC(inv_stereo(Complex(0, 0)), Vector(0, -1, 0), 1e-3f);
}

// ============================================================================
// Complex
// ============================================================================

/**
 * @brief Verifies Complex +, -, *, / including the division-by-zero
 *        conventions (0/0 -> 0, nonzero/0 -> large magnitude).
 */
inline void test_complex_arithmetic() {
  Complex a(1, 2), b(3, 4);

  HS_EXPECT_COMPLEX(a + b, Complex(4, 6), 1e-6f);
  HS_EXPECT_COMPLEX(a - b, Complex(-2, -2), 1e-6f);
  // (1+2i)(3+4i) = 3 + 4i + 6i + 8i² = -5 + 10i
  HS_EXPECT_COMPLEX(a * b, Complex(-5, 10), 1e-5f);

  // (a / b) * b ≈ a
  Complex q = a / b;
  HS_EXPECT_COMPLEX(q * b, a, 1e-4f);

  HS_EXPECT_COMPLEX(a * Complex(1, 0), a, 1e-6f);

  // 0 / 0 → 0 by convention.
  HS_EXPECT_COMPLEX(Complex(0, 0) / Complex(0, 0), Complex(0, 0), 1e-6f);

  // Nonzero / 0 → large magnitude in the numerator direction.
  Complex inf_dir = Complex(1, 0) / Complex(0, 0);
  HS_EXPECT_TRUE(std::abs(inf_dir.re) > 1e3f);
}

// ============================================================================
// Mobius
// ============================================================================

/**
 * @brief Verifies MobiusParams constructors (8-float, 4-Complex, identity
 *        default) and that the a,b,c,d coefficients land in order.
 */
inline void test_mobius_params_accessors() {
  MobiusParams p(1, 2, 3, 4, 5, 6, 7, 8);
  HS_EXPECT_COMPLEX(p.a, Complex(1, 2), 0.0f);
  HS_EXPECT_COMPLEX(p.b, Complex(3, 4), 0.0f);
  HS_EXPECT_COMPLEX(p.c, Complex(5, 6), 0.0f);
  HS_EXPECT_COMPLEX(p.d, Complex(7, 8), 0.0f);

  MobiusParams q(Complex(1, 2), Complex(3, 4), Complex(5, 6), Complex(7, 8));
  HS_EXPECT_COMPLEX(q.a, Complex(1, 2), 0.0f);
  HS_EXPECT_COMPLEX(q.d, Complex(7, 8), 0.0f);

  // Identity Mobius default: a=d=1, b=c=0.
  MobiusParams id;
  HS_EXPECT_COMPLEX(id.a, Complex(1, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.b, Complex(0, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.c, Complex(0, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.d, Complex(1, 0), 0.0f);
}

/**
 * @brief Verifies mobius(z, params): identity, pure translation, and pure
 *        scaling cases.
 */
inline void test_mobius_transform() {
  Complex z(0.3f, 0.7f);

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

/**
 * @brief Verifies gnomonic/inv_gnomonic roundtrips (hemisphere sign passed
 *        explicitly), the pole pre-image, saturated-input pole return, and
 *        near-equator clamping.
 */
inline void test_gnomonic_roundtrip() {
  Vector vUp = Vector(0.3f, 0.8f, 0.4f).normalized();
  Complex zUp = gnomonic(vUp);
  Vector vUp_back = inv_gnomonic(zUp, 1.0f);
  HS_EXPECT_VEC(vUp_back, vUp, 5e-3f);

  // Lower hemisphere: the hemisphere sign is passed to inv_gnomonic explicitly.
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

/**
 * @brief Verifies Spline::cubic_fast and cubic_slerp hit the endpoints, give
 *        unit-length midpoints, and that the cubic() dispatcher routes to the
 *        right variant by mode.
 */
inline void test_spline_cubic_endpoints() {
  const Vector p0(1, 0, 0);
  const Vector p1 = Vector(0.7f, 0.7f, 0).normalized();
  const Vector p2(0, 1, 0);
  const Vector p3 = Vector(-0.7f, 0.7f, 0).normalized();

  HS_EXPECT_VEC(Spline::cubic_fast(p0, p1, p2, p3, 0.0f), p0, 1e-3f);
  HS_EXPECT_VEC(Spline::cubic_fast(p0, p1, p2, p3, 1.0f), p3, 1e-3f);

  HS_EXPECT_VEC(Spline::cubic_slerp(p0, p1, p2, p3, 0.0f), p0, 5e-3f);
  HS_EXPECT_VEC(Spline::cubic_slerp(p0, p1, p2, p3, 1.0f), p3, 5e-3f);

  Vector mFast = Spline::cubic_fast(p0, p1, p2, p3, 0.5f);
  Vector mSlerp = Spline::cubic_slerp(p0, p1, p2, p3, 0.5f);
  HS_EXPECT_NEAR(mFast.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(mSlerp.length(), 1.0f, 1e-3f);

  // cubic() dispatches by mode.
  HS_EXPECT_VEC(Spline::cubic(p0, p1, p2, p3, 0.5f, SplineMode::Fast), mFast,
                1e-4f);
  HS_EXPECT_VEC(Spline::cubic(p0, p1, p2, p3, 0.5f, SplineMode::Geodesic),
                mSlerp, 1e-4f);
}

/**
 * @brief Verifies the cubic_fast degenerate case: control points whose t=0.5
 *        Bezier blend cancels to the zero vector must degrade gracefully to p1
 *        rather than trap in the strict normalized().
 */
inline void test_spline_cubic_fast_degenerate_fallback() {
  // Bezier weights at t=0.5 are 0.125/0.375/0.375/0.125, so for these antipodal
  // control points 0.125*p0 + 0.375*p1 + 0.375*p2 + 0.125*p3 cancels to (0,0,0).
  const Vector p0(1, 0, 0), p1(1, 0, 0), p2(-1, 0, 0), p3(-1, 0, 0);
  HS_EXPECT_VEC(Spline::cubic_fast(p0, p1, p2, p3, 0.5f), p1, 1e-6f);
}

/**
 * @brief Verifies Spline::catmull_rom_tangents: tension 0 gives a geodesic
 *        segment (tangents at the endpoints), tension 1 gives full smoothing
 *        (slerp midpoints), and the emitted control points stay unit-length.
 */
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

  Spline::catmull_rom_tangents(prev, start, end, next, 0.5f, cp1, cp2);
  HS_EXPECT_NEAR(cp1.length(), 1.0f, 5e-3f);
  HS_EXPECT_NEAR(cp2.length(), 1.0f, 5e-3f);
}

// ============================================================================
// wrap_index (core/math/rotate.h) — folds a float index into [0, m)
// ============================================================================

/**
 * @brief Verifies wrap_index folds a float index into [0, m): preserves
 *        in-range values, wraps at/above the period, folds negatives, and
 *        stays in range over many periods.
 */
inline void test_wrap_index() {
  const int m = 288;

  HS_EXPECT_NEAR(wrap_index(0.0f, m), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(wrap_index(0.5f, m), 0.5f, 1e-5f);
  HS_EXPECT_NEAR(wrap_index(287.9f, m), 287.9f, 1e-3f);

  HS_EXPECT_NEAR(wrap_index(static_cast<float>(m), m), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(wrap_index(m + 1.5f, m), 1.5f, 1e-4f);

  // Negatives fold into [0, m): -0.5 -> 287.5.
  HS_EXPECT_NEAR(wrap_index(-0.5f, m), 287.5f, 1e-3f);
  HS_EXPECT_NEAR(wrap_index(-1.5f, m), 286.5f, 1e-3f);
  HS_EXPECT_NEAR(wrap_index(-static_cast<float>(m) + 0.25f, m), 0.25f, 1e-3f);

  for (int i = -3 * m; i <= 3 * m; ++i) {
    float w = wrap_index(i * 0.5f, m);
    HS_EXPECT_TRUE(w >= 0.0f && w < static_cast<float>(m));
  }
}

// ============================================================================
// Runner
// ============================================================================

/**
 * @brief Runs every 3dmath test case.
 * @return Accumulated failure count (0 on success), suitable for use as a
 *         process exit code.
 */
inline int run_3dmath_tests() {
  hs_test::ModuleFixture fixture("3dmath");

  test_constants();
  test_quintic_kernel();

  test_fast_atan2();
  test_fast_acos();
  test_fast_sinf_cosf();
  test_fast_cbrt();

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
  test_quaternion_from_basis();
  test_rotate();

  test_vector_slerp();
  test_vector_slerp_antipodal_monotonic();
  test_quaternion_slerp();

  test_stereo_roundtrip();
  test_complex_arithmetic();
  test_mobius_params_accessors();
  test_mobius_transform();
  test_gnomonic_roundtrip();

  test_spline_cubic_endpoints();
  test_spline_cubic_fast_degenerate_fallback();
  test_spline_catmull_rom();

  test_wrap_index();

  return fixture.result();
}

} // namespace math3d
} // namespace hs_test

