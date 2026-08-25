/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/math/3dmath.h.
 *
 * Usage:
 *   #include "tests/test_3dmath.h"
 *   int main() { return hs_test::math3d_tests::run_3dmath_tests(); }
 *
 * Self-contained header — no external test framework. All test functions
 * are inline; the runner returns the failure count for use as a process
 * exit code.
 */
#pragma once

#include "core/math/3dmath.h"
#include "core/math/lenses.h"
#include "core/math/stereographic.h"
#include "core/math/rotate.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/vec_test_util.h"

#include <algorithm>

namespace hs_test {
namespace math3d_tests {

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
  HS_EXPECT_NEAR_REL(STEREO_INF, 1e4f, 1e-7f);
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
// (measured peak errors: atan2 ~3.8e-3 rad, acos ~5.0e-5 rad, sin ~1.6e-3)
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
 * @brief Verifies diamond_angle's [0,4) range, cardinal anchors, scale
 *        invariance, and strict monotonicity with std::atan2.
 * @details The sweep walks the circle counter-clockwise from +x, which is the
 *          order the pseudo-angle must reproduce for it to bin a direction. The
 *          tiny-negative-y probes cover the fourth-quadrant seam, where 4 + r
 *          rounds back up to exactly 4 and has to fold to 0.
 */
inline void test_diamond_angle() {
  HS_EXPECT_NEAR(diamond_angle(0.0f, 1.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(diamond_angle(1.0f, 0.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(diamond_angle(0.0f, -1.0f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(diamond_angle(-1.0f, 0.0f), 3.0f, 1e-6f);

  // Degenerate origin.
  HS_EXPECT_NEAR(diamond_angle(0.0f, 0.0f), 0.0f, 1e-6f);

  // Scale invariance: only the direction matters.
  for (int i = 0; i < 32; ++i) {
    float a = (i * 2.0f * PI_F) / 32.0f;
    float y = std::sin(a), x = std::cos(a);
    HS_EXPECT_NEAR(diamond_angle(y * 1e4f, x * 1e4f), diamond_angle(y, x),
                   1e-5f);
    HS_EXPECT_NEAR(diamond_angle(y * 1e-4f, x * 1e-4f), diamond_angle(y, x),
                   1e-5f);
  }

  // Strictly increasing over one counter-clockwise turn, always in [0, 4).
  float prev = -1.0f;
  for (int i = 0; i < 1024; ++i) {
    float a = (i * 2.0f * PI_F) / 1024.0f;
    float d = diamond_angle(std::sin(a), std::cos(a));
    HS_EXPECT_TRUE(d >= 0.0f && d < 4.0f);
    HS_EXPECT_TRUE(d > prev);
    prev = d;
  }

  // Fourth-quadrant seam: y -> 0- with x > 0 must stay inside [0, 4).
  for (float y : {-1e-3f, -1e-6f, -1e-9f, -1e-12f, -1e-20f, -1e-30f}) {
    float d = diamond_angle(y, 1.0f);
    HS_EXPECT_TRUE(d >= 0.0f && d < 4.0f);
  }
}

/**
 * @brief Verifies fast_reciprocal across normal inputs and outputs at its
 * documented peak relative error.
 */
inline void test_fast_reciprocal() {
  HS_EXPECT_NEAR(fast_reciprocal(1.0f), 1.0f, 7e-6f);
  HS_EXPECT_NEAR(fast_reciprocal(4.0f), 0.25f, 7e-6f * 0.25f);
  HS_EXPECT_NEAR(fast_reciprocal(0.25f), 4.0f, 7e-6f * 4.0f);

  for (int exponent = -120; exponent <= 120; ++exponent) {
    for (int mantissa = 0; mantissa < 8; ++mantissa) {
      const float x = std::ldexp(1.0f + 0.125f * mantissa, exponent);
      const float reference = 1.0f / x;
      HS_EXPECT_TRUE(std::abs(fast_reciprocal(x) - reference) / reference <=
                     7e-6f);
    }
  }
}

/**
 * @brief Verifies fast_rsqrt against 1/sqrt over a wide sweep at the documented
 *        ~5e-6 peak relative error.
 */
inline void test_fast_rsqrt() {
  HS_EXPECT_NEAR(fast_rsqrt(1.0f), 1.0f, 5e-6f);
  HS_EXPECT_NEAR(fast_rsqrt(4.0f), 0.5f, 5e-6f * 0.5f);
  HS_EXPECT_NEAR(fast_rsqrt(0.25f), 2.0f, 5e-6f * 2.0f);

  // Both exponent parities across ~12 decades: the bit-hack seed's quality
  // alternates with the low exponent bit, so a one-decade sweep would miss half
  // the error surface.
  for (int i = 0; i <= 512; ++i) {
    float x = std::pow(10.0f, -6.0f + (12.0f * i) / 512.0f);
    float ref = 1.0f / std::sqrt(x);
    HS_EXPECT_TRUE(std::abs(fast_rsqrt(x) - ref) / ref <= 5e-6f);
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
    HS_EXPECT_NEAR(fast_acos(x), std::acos(x), 2e-4f);
  }

  // Peak error (~5.0e-5) near x = 0.0807, between the sweep's samples.
  HS_EXPECT_NEAR(fast_acos(0.0807f), std::acos(0.0807f), 5.1e-5f);
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
 * @brief Verifies fast_cbrt3 tracks three separate fast_cbrt calls, clamps
 *        non-positive inputs, and stays accurate inside the documented
 *        ~1.3e11 / ~3e-13 window of the shared reciprocal.
 */
inline void test_fast_cbrt3() {
  for (int i = 0; i < 64; ++i) {
    float x[3], o[3];
    for (int j = 0; j < 3; ++j)
      x[j] = 1.0f + 65534.0f * ((i * 3 + j) % 61) / 60.0f;
    fast_cbrt3(x[0], x[1], x[2], o[0], o[1], o[2]);
    for (int j = 0; j < 3; ++j) {
      float ref = fast_cbrt(x[j]);
      HS_EXPECT_TRUE(std::abs(o[j] - ref) / ref <= 1e-6f);
      float rel = std::abs(o[j] - std::cbrt(x[j])) / std::cbrt(x[j]);
      HS_EXPECT_TRUE(rel <= 2.3e-5f);
    }
  }

  // Non-positive inputs clamp to 0 without poisoning the shared product.
  {
    float o[3];
    fast_cbrt3(-1.0f, 0.0f, 8.0f, o[0], o[1], o[2]);
    HS_EXPECT_NEAR(o[0], 0.0f, 1e-7f);
    HS_EXPECT_NEAR(o[1], 0.0f, 1e-7f);
    HS_EXPECT_NEAR(o[2], 2.0f, 2.3e-5f * 2.0f);
  }

  // Just inside each end of the documented window the result still tracks the
  // scalar helper; the tiny end is below fast_cbrt's own accurate domain, so
  // only the re-association is under test there.
  for (float v : {1.0e11f, 1.0e-12f}) {
    float o[3];
    fast_cbrt3(v, v, v, o[0], o[1], o[2]);
    float ref = fast_cbrt(v);
    for (int j = 0; j < 3; ++j) {
      HS_EXPECT_TRUE(std::isfinite(o[j]));
      HS_EXPECT_TRUE(std::abs(o[j] - ref) / ref <= 1e-5f);
    }
  }
}

/**
 * @brief Verifies fast_cbrt6 tracks six separate fast_cbrt calls and holds the
 *        documented ~2.3e-5 error against cbrtf across its usable domain.
 * @details The shared reciprocal re-associates the arithmetic, so agreement
 *          with fast_cbrt is ~4e-7 relative rather than exact. Also pins the
 *          x<=0 -> 0 clamp and the ~4.2e5 ceiling above which the
 *          six-denominator product overflows: a future caller widening the
 *          input range past the u16-magnitude domain trips this.
 */
inline void test_fast_cbrt6() {
  // Agreement with the scalar helper across the u16-magnitude LMS range.
  for (int i = 0; i < 64; ++i) {
    float x[6], o[6];
    for (int j = 0; j < 6; ++j)
      x[j] = 1.0f + 65534.0f * ((i * 6 + j) % 61) / 60.0f;
    fast_cbrt6(x, o);
    for (int j = 0; j < 6; ++j) {
      float ref = fast_cbrt(x[j]);
      HS_EXPECT_TRUE(std::abs(o[j] - ref) / ref <= 1e-6f);
      float rel = std::abs(o[j] - std::cbrt(x[j])) / std::cbrt(x[j]);
      HS_EXPECT_TRUE(rel <= 2.3e-5f);
    }
  }

  // Non-positive inputs clamp to 0 without poisoning the shared product.
  {
    const float x[6] = {-1.0f, 0.0f, -0.0f, 8.0f, 1.0f, 27.0f};
    float o[6];
    fast_cbrt6(x, o);
    HS_EXPECT_NEAR(o[0], 0.0f, 1e-7f);
    HS_EXPECT_NEAR(o[1], 0.0f, 1e-7f);
    HS_EXPECT_NEAR(o[2], 0.0f, 1e-7f);
    HS_EXPECT_NEAR(o[3], 2.0f, 2.3e-5f * 2.0f);
    HS_EXPECT_NEAR(o[4], 1.0f, 2.3e-5f);
    HS_EXPECT_NEAR(o[5], 3.0f, 2.3e-5f * 3.0f);
  }

  // Just under the documented overflow ceiling the result is still accurate.
  {
    const float v = 4.0e5f;
    const float x[6] = {v, v, v, v, v, v};
    float o[6];
    fast_cbrt6(x, o);
    for (int j = 0; j < 6; ++j) {
      HS_EXPECT_TRUE(std::isfinite(o[j]));
      float rel = std::abs(o[j] - std::cbrt(v)) / std::cbrt(v);
      HS_EXPECT_TRUE(rel <= 2.3e-5f);
    }
  }
}

/**
 * @brief Verifies fast_expf anchors, the large-magnitude saturation to 0, and
 *        the documented ~7.4e-4 peak relative error over the x<=0 domain.
 */
inline void test_fast_expf() {
  HS_EXPECT_NEAR(fast_expf(0.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(fast_expf(-1.0f), std::exp(-1.0f), 7.5e-4f);
  HS_EXPECT_NEAR(fast_expf(-10.0f), std::exp(-10.0f),
                 7.5e-4f * std::exp(-10.0f));

  for (int i = 0; i <= 512; ++i) {
    float x = -30.0f + (30.0f * i) / 512.0f;
    float ref = std::exp(x);
    float rel = std::abs(fast_expf(x) - ref) / ref;
    HS_EXPECT_TRUE(rel <= 7.5e-4f);
  }

  // Large-magnitude arguments saturate to 0.
  HS_EXPECT_NEAR(fast_expf(-100.0f), 0.0f, 1e-30f);
  HS_EXPECT_NEAR(fast_expf(-200.0f), 0.0f, 1e-30f);
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

/**
 * @brief Verifies fast_sincosf_0_pi reproduces fast_sinf/fast_cosf across its
 *        documented [0, pi] domain.
 */
inline void test_fast_sincosf_0_pi() {
  for (int i = 0; i <= 256; ++i) {
    float x = (i * PI_F) / 256.0f;
    float s, c;
    fast_sincosf_0_pi(x, s, c);
    HS_EXPECT_NEAR(s, fast_sinf(x), 1e-6f);
    HS_EXPECT_NEAR(c, fast_cosf(x), 1e-6f);
    HS_EXPECT_NEAR(s, std::sin(x), 1.8e-3f);
    HS_EXPECT_NEAR(c, std::cos(x), 1.8e-3f);
    HS_EXPECT_NEAR(s * s + c * c, 1.0f, 5e-3f);
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
  HS_EXPECT_EQ(a, b);
  HS_EXPECT_FALSE(a == c);
  HS_EXPECT_TRUE(a != c);
  HS_EXPECT_FALSE(a != b);

  HS_EXPECT_EQ(Vector(1, 0, 0), Vector(1.00005f, 0, 0));
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
  HS_EXPECT_VEC(normalized_or(Vector(0, 0, 0), Vector(1, 0, 0)),
                Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(normalized_or(Vector(0, 6, 0), Vector(1, 0, 0)),
                Vector(0, 1, 0), 1e-6f);
}

// ============================================================================
// Snorm3
// ============================================================================

/**
 * @brief Pins Snorm3's documented round-trip accuracy: per-component error
 *        within 1/65534, chord error within ~2.6e-5, endpoint codes exact, and
 *        out-of-domain components saturated rather than wrapped.
 */
inline void test_snorm3_roundtrip_bound() {
  // The +1e-7f absorbs the rounding of the decode multiply, which can carry a
  // worst-case quantization by up to one relative ulp past the exact 1/65534.
  constexpr float COMPONENT_BOUND = 1.0f / 65534.0f + 1e-7f;
  constexpr float CHORD_BOUND = 2.65e-5f;

  const Snorm3 endpoints = Snorm3::encode(Vector(1, 0, -1));
  HS_EXPECT_EQ(endpoints.x, 32767);
  HS_EXPECT_EQ(endpoints.y, 0);
  HS_EXPECT_EQ(endpoints.z, -32767);
  HS_EXPECT_VEC(endpoints.decode(), Vector(1, 0, -1), COMPONENT_BOUND);

  const Snorm3 saturated = Snorm3::encode(Vector(4.0f, -9.0f, 1.0f + 1e-3f));
  HS_EXPECT_EQ(saturated.x, 32767);
  HS_EXPECT_EQ(saturated.y, -32767);
  HS_EXPECT_EQ(saturated.z, 32767);

  hs::Pcg32 rng(20260803u);
  for (int i = 0; i < 256; ++i) {
    const Vector v = normalized_or(Vector(rand_uniform(rng, -1.0f, 1.0f),
                                          rand_uniform(rng, -1.0f, 1.0f),
                                          rand_uniform(rng, -1.0f, 1.0f)),
                                   Vector(1, 0, 0));
    const Vector decoded = Snorm3::encode(v).decode();
    HS_EXPECT_LE(std::abs(decoded.x - v.x), COMPONENT_BOUND);
    HS_EXPECT_LE(std::abs(decoded.y - v.y), COMPONENT_BOUND);
    HS_EXPECT_LE(std::abs(decoded.z - v.z), COMPONENT_BOUND);
    HS_EXPECT_LE(distance_between(decoded, v), CHORD_BOUND);
    // Near-unit, not unit: callers that need exact length must renormalize.
    HS_EXPECT_NEAR(decoded.length(), 1.0f, CHORD_BOUND);
  }
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
 * @brief Verifies conjugate() and inverse().
 * @details For a unit q the two coincide; q*q^-1 = identity holds for unit
 *          and non-unit quaternions alike.
 */
inline void test_quaternion_conjugate_inverse() {
  Quaternion q(0.5f, 0.5f, 0.5f, 0.5f); // unit
  Quaternion conj = q.conjugate();
  HS_EXPECT_QUAT(conj, Quaternion(0.5f, -0.5f, -0.5f, -0.5f), 1e-6f);

  HS_EXPECT_QUAT(q.inverse(), conj, 1e-6f);

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
  HS_EXPECT_EQ(a, b);
  Quaternion c(1.001f, 2, 3, 4); // beyond TOLERANCE
  HS_EXPECT_FALSE(a == c);
  Quaternion d(1.00001f, 2, 3, 4); // within TOLERANCE
  HS_EXPECT_EQ(a, d);
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
 * @brief Verifies least_parallel_axis picks +Y only near +/-X, is scale
 *        invariant, and always seeds a well-conditioned cross.
 */
inline void test_least_parallel_axis() {
  HS_EXPECT_VEC(least_parallel_axis(Vector(1, 0, 0)), Vector(0, 1, 0), 1e-6f);
  HS_EXPECT_VEC(least_parallel_axis(Vector(-1, 0, 0)), Vector(0, 1, 0), 1e-6f);
  HS_EXPECT_VEC(least_parallel_axis(Vector(0, 1, 0)), Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(least_parallel_axis(Vector(0, 0, 1)), Vector(1, 0, 0), 1e-6f);

  // Scale invariant: the cosine test is against |v|^2, not the raw component.
  HS_EXPECT_VEC(least_parallel_axis(Vector(50, 0, 0)), Vector(0, 1, 0), 1e-6f);
  HS_EXPECT_VEC(least_parallel_axis(Vector(0, 50, 0)), Vector(1, 0, 0), 1e-6f);

  // The whole point: cross(axis, v) never collapses. The worst case sits just
  // inside the COS_AXIS_PARALLEL switch, where sin^2 is ~2*TOLERANCE.
  for (int i = 0; i <= 128; ++i) {
    for (int j = 0; j <= 128; ++j) {
      float phi = (i * PI_F) / 128.0f;
      float theta = (j * 2.0f * PI_F) / 128.0f;
      Vector v(std::sin(phi) * std::cos(theta), std::cos(phi),
               std::sin(phi) * std::sin(theta));
      Vector c = cross(least_parallel_axis(v), v);
      HS_EXPECT_TRUE(dot(c, c) >= 1e-4f);
    }
  }
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

  Quaternion q = make_rotation(Vector(1, 2, 3).normalized(), 1.234f);
  Vector r = rotate(v, q);
  HS_EXPECT_NEAR(r.length(), v.length(), 1e-3f);

  // Composition: rotate(rotate(v, q1), q2) == rotate(v, q2 * q1).
  Quaternion q1 = make_rotation(Vector(0, 1, 0), 0.3f);
  Quaternion q2 = make_rotation(Vector(1, 0, 0), 0.5f);
  Vector via_sequential = rotate(rotate(v, q1), q2);
  Vector via_composed = rotate(v, q2 * q1);
  HS_EXPECT_VEC(via_sequential, via_composed, 5e-3f);
}

/**
 * @brief Pins RotationMatrix to rotate(): the expanded rows must reproduce the
 *        quaternion sandwich for every orientation, and the rows must stay
 *        orthonormal.
 */
inline void test_rotation_matrix_matches_rotate() {
  const Vector axes[] = {Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1),
                         Vector(1, 2, 3).normalized(),
                         Vector(-2, 0.5f, 1).normalized()};
  const Vector samples[] = {Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1),
                            Vector(1, 2, 3), Vector(-4, 0.25f, 2)};

  for (const Vector &axis : axes) {
    for (int i = 0; i <= 12; ++i) {
      float angle = (2.0f * PI_F * i) / 12.0f;
      Quaternion q = make_rotation(axis, angle);
      RotationMatrix m(q);
      for (const Vector &v : samples) {
        HS_EXPECT_VEC(m.apply(v), rotate(v, q), 1e-4f);
      }
      HS_EXPECT_NEAR(dot(m.r0, m.r0), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(dot(m.r1, m.r1), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(dot(m.r2, m.r2), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(dot(m.r0, m.r1), 0.0f, 1e-5f);
      HS_EXPECT_NEAR(dot(m.r0, m.r2), 0.0f, 1e-5f);
      HS_EXPECT_NEAR(dot(m.r1, m.r2), 0.0f, 1e-5f);
    }
  }
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
 * @brief Verifies nlerp_unit: endpoints, the unit-length midpoint, and the
 *        cancelling-blend fallback at math::EPS_BLEND_LEN_SQ.
 */
inline void test_vector_nlerp_unit() {
  const Vector a(1, 0, 0), b(0, 1, 0);
  HS_EXPECT_VEC(nlerp_unit(a, b, 0.0f), a, 1e-5f);
  HS_EXPECT_VEC(nlerp_unit(a, b, 1.0f), b, 1e-5f);

  const Vector mid = nlerp_unit(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.x, std::sqrt(2.0f) * 0.5f, 1e-5f);
  HS_EXPECT_NEAR(mid.y, std::sqrt(2.0f) * 0.5f, 1e-5f);
  HS_EXPECT_NEAR(mid.z, 0.0f, 1e-5f);

  for (int i = 0; i <= 10; ++i)
    HS_EXPECT_NEAR(nlerp_unit(a, b, static_cast<float>(i) / 10.0f).length(),
                   1.0f, 1e-5f);

  // Antipodal endpoints cancel: the blend carries no direction, so it holds a.
  const Vector p(0, 1, 0), ap(0, -1, 0);
  HS_EXPECT_VEC(nlerp_unit(p, ap, 0.5f), p, 1e-5f);
  HS_EXPECT_VEC(nlerp_unit(p, ap, 0.5f + 1e-5f), p, 1e-5f);
  HS_EXPECT_VEC(nlerp_unit(p, ap, 0.5f + 1e-3f), ap, 1e-5f);

  // Near-antipodal but non-cancelling: the residual still resolves to a unit
  // direction.
  const Vector near_ap = Vector(0.01f, -1.0f, 0.0f).normalized();
  HS_EXPECT_NEAR(nlerp_unit(p, near_ap, 0.5f).length(), 1.0f, 1e-4f);
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

  // Identical endpoints on the long arc drive d to -1 via the sign fixup, the
  // case that once collapsed to the zero quaternion; the fallback must stay unit.
  Quaternion long_degenerate = slerp(q, q, 0.5f, true);
  HS_EXPECT_NEAR(long_degenerate.magnitude(), 1.0f, 1e-3f);
}

inline void test_scaled_rotation_delta() {
  const Quaternion id;
  const Quaternion q = make_rotation(Vector(0, 1, 0), PI_F * 0.5f);

  // Both extremes are exact and skip the slerp entirely.
  HS_EXPECT_QUAT(scaled_rotation_delta(q, 1.0f), q, 1e-6f);
  HS_EXPECT_QUAT(scaled_rotation_delta(q, 0.0f), id, 1e-6f);
  // Identity in, identity out, at any fraction.
  HS_EXPECT_QUAT(scaled_rotation_delta(id, 0.5f), id, 5e-3f);

  // Half the arc, applied twice, recovers the full delta.
  const Quaternion half = scaled_rotation_delta(q, 0.5f);
  HS_EXPECT_NEAR(half.magnitude(), 1.0f, 1e-3f);
  const Vector v(1, 0, 0);
  HS_EXPECT_VEC(rotate(rotate(v, half), half), rotate(v, q), 1e-2f);

  // The scaled turn is monotone in amount: the angle away from the start grows.
  float previous = -1.0f;
  for (int step = 0; step <= 8; ++step) {
    const float amount = static_cast<float>(step) / 8.0f;
    const Quaternion scaled = scaled_rotation_delta(q, amount);
    HS_EXPECT_NEAR(scaled.magnitude(), 1.0f, 1e-3f);
    const float turned = 1.0f - dot(rotate(v, scaled), v);
    HS_EXPECT_GT(turned, previous);
    previous = turned;
  }
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
      Vector(1, 0, 0),          Vector(0, 0, 1),
      Vector(-1, 0, 0),         Vector(0, -1, 0),
      Vector(0.6f, 0.0f, 0.8f), Vector(0.5f, 0.5f, 0.7071f).normalized(),
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
  // At this scale the unit vector's y rounds to 1, so denom is exactly zero.
  Vector nearPole = Vector(6e-5f, 1.0f, 2.1e-5f).normalized();
  Complex zCap = stereo(nearPole);
  HS_EXPECT_NEAR(std::sqrt(zCap.re * zCap.re + zCap.im * zCap.im), STEREO_INF,
                 1.0f);
  HS_EXPECT_NEAR(std::atan2(zCap.im, zCap.re),
                 std::atan2(nearPole.z, nearPole.x), 1e-3f);

  // The cap boundary is the crossover where the raw quotient would reach the
  // sentinel, not a step: outside it the projection stays below STEREO_INF.
  Complex zOut = stereo(Vector(0.006f, 0.99998f, 0.0021f).normalized());
  HS_EXPECT_LT(std::sqrt(zOut.re * zOut.re + zOut.im * zOut.im), STEREO_INF);

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
 * @brief Verifies Complex +, -, *, / (ordinary complex division) and
 *        project_div's projection conventions (0/0 -> 0, nonzero/0 -> large
 *        magnitude in the numerator direction).
 */
inline void test_complex_arithmetic() {
  Complex a(1, 2), b(3, 4);

  HS_EXPECT_COMPLEX(a + b, Complex(4, 6), 1e-6f);
  HS_EXPECT_COMPLEX(a - b, Complex(-2, -2), 1e-6f);
  // (1+2i)(3+4i) = 3 + 4i + 6i + 8i² = -5 + 10i
  HS_EXPECT_COMPLEX(a * b, Complex(-5, 10), 1e-5f);

  // operator/ is ordinary complex division: (a / b) * b ≈ a.
  Complex q = a / b;
  HS_EXPECT_COMPLEX(q * b, a, 1e-4f);

  HS_EXPECT_COMPLEX(a * Complex(1, 0), a, 1e-6f);

  // project_div matches ordinary division away from the singularity.
  HS_EXPECT_COMPLEX(project_div(a, b), a / b, 1e-6f);

  // project_div convention: 0 / 0 → 0.
  HS_EXPECT_COMPLEX(project_div(Complex(0, 0), Complex(0, 0)), Complex(0, 0),
                    1e-6f);

  // project_div convention: nonzero / 0 → large magnitude in the numerator
  // direction.
  Complex inf_dir = project_div(Complex(1, 0), Complex(0, 0));
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
 * @brief Verifies mobius(z, params): identity, pure translation, pure scaling,
 *        and the inverting c != 0 branch including its pole.
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

  // Inversion: 1/z = conj(z) / |z|^2.
  MobiusParams inv(0, 0, 1, 0, 1, 0, 0, 0);
  const float r2 = z.re * z.re + z.im * z.im;
  HS_EXPECT_COMPLEX(mobius(z, inv), Complex(z.re / r2, -z.im / r2), 1e-4f);

  // General c != 0: (z + 1) / (z - 1).
  MobiusParams gen(1, 0, 1, 0, 1, 0, -1, 0);
  HS_EXPECT_COMPLEX(mobius(z, gen), Complex(-0.42857143f, -1.42857143f), 1e-4f);

  // Its pole z = -d/c = 1 vanishes the denominator: project_div substitutes the
  // point at infinity along the numerator's direction, which inv_stereo reads
  // back as the north pole.
  Complex at_pole = mobius(Complex(1, 0), gen);
  HS_EXPECT_NEAR(at_pole.re, STEREO_INF, 1.0f);
  HS_EXPECT_NEAR(at_pole.im, 0.0f, 1e-4f);
  HS_EXPECT_VEC(inv_stereo(at_pole), Vector(0, 1, 0), 1e-6f);

  // A pole of the degenerate map (ad - bc == 0) is the indeterminate 0/0 form.
  MobiusParams degenerate(1, 0, -1, 0, 1, 0, -1, 0);
  HS_EXPECT_COMPLEX(mobius(Complex(1, 0), degenerate), Complex(0, 0), 0.0f);
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

  // Saturated input → the equator point in the direction of z, sign-dependent
  // (gnomonic's singularity is the equator, not a pole).
  HS_EXPECT_VEC(inv_gnomonic(Complex(STEREO_INF, 0), 1.0f), Vector(1, 0, 0),
                1e-3f);
  HS_EXPECT_VEC(inv_gnomonic(Complex(STEREO_INF, 0), -1.0f), Vector(-1, 0, 0),
                1e-3f);
  HS_EXPECT_VEC(inv_gnomonic(Complex(0, -STEREO_INF), 1.0f), Vector(0, 0, -1),
                1e-3f);
  // A magnitude far past the sentinel must not square to infinity.
  HS_EXPECT_VEC(inv_gnomonic(Complex(3e30f, 4e30f), 1.0f),
                Vector(0.6f, 0.0f, 0.8f), 1e-3f);
  // The sentinel test is radial: a diagonal point past the recognition radius
  // snaps back even though neither component reaches it alone.
  HS_EXPECT_VEC(inv_gnomonic(Complex(4e3f, 4e3f), 1.0f),
                Vector(0.70710678f, 0.0f, 0.70710678f), 1e-3f);

  // Near-equator inputs get clamped to STEREO_INF
  Complex zEq = gnomonic(Vector(1.0f, 1e-10f, 0.0f));
  HS_EXPECT_TRUE(std::abs(zEq.re) >= STEREO_INF - 1.0f);

  // Round-trip identity through the singularity: |v.y| below ~2e-4 saturates
  // the projection, and the inverse must still land on the input.
  for (float y : {2e-4f, 1e-5f, 1e-9f, 0.0f, -1e-9f, -1e-5f, -2e-4f}) {
    for (float theta = 0.0f; theta < 2.0f * PI_F; theta += 0.37f) {
      const Vector v = Vector(cosf(theta), y, sinf(theta)).normalized();
      const Vector back = inv_gnomonic(gnomonic(v), copysignf(1.0f, y));
      HS_EXPECT_VEC(back, v, 1e-3f);
    }
  }

  // The floored divisor keys on the sign bit, so -0.0f projects like the tiny
  // negatives it is the limit of, not like +0.0f.
  Complex z_neg_zero = gnomonic(Vector(1.0f, -0.0f, 0.0f));
  Complex z_tiny_neg = gnomonic(Vector(1.0f, -1e-12f, 0.0f));
  Complex z_pos_zero = gnomonic(Vector(1.0f, 0.0f, 0.0f));
  HS_EXPECT_TRUE(std::signbit(z_neg_zero.re) == std::signbit(z_tiny_neg.re));
  HS_EXPECT_TRUE(std::signbit(z_neg_zero.re) != std::signbit(z_pos_zero.re));
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

/**
 * @brief Verifies hash01's frozen outputs, range, and seed independence.
 * @details The sorted-set check is the load-bearing one: a seed that only
 * permutes the lattice passes pointwise inequality but reproduces the same
 * multiset of values, so two seeds would be one stream re-indexed.
 */
inline void test_hash01() {
  // Frozen mixer output: integer-only, so exact on every host.
  HS_EXPECT_NEAR(hash01(42u, 7u), 0.0169150233f, 1e-9f);
  HS_EXPECT_NEAR(hash01(0u, 0u), 0.030199945f, 1e-9f);
  HS_EXPECT_NEAR(hash01(4294967295u, 123u), 0.722669423f, 1e-9f);
  for (uint32_t i = 0; i < 256; ++i) {
    float h = hash01(i, 0u);
    HS_EXPECT_GE(h, 0.0f);
    HS_EXPECT_LT(h, 1.0f);
  }
  int differing = 0;
  for (uint32_t i = 0; i < 32; ++i)
    if (hash01(i, 1u) != hash01(i, 2u))
      differing++;
  HS_EXPECT_GT(differing, 24);

  // Seeds select streams, not permutations: the value sets must differ.
  constexpr uint32_t N = 256;
  const uint32_t pairs[3][2] = {{0u, 1u}, {1u, 2u}, {7u, 8u}};
  for (const auto &pair : pairs) {
    float a[N], b[N];
    for (uint32_t i = 0; i < N; ++i) {
      a[i] = hash01(i, pair[0]);
      b[i] = hash01(i, pair[1]);
    }
    std::sort(a, a + N);
    std::sort(b, b + N);
    HS_EXPECT_TRUE(!std::equal(a, a + N, b));
  }
}

/**
 * @brief Verifies value noise hits lattice hashes at integers, stays in range,
 *        and is continuous — including across the x=0 cell boundary, where the
 *        negative-coordinate int cast could break adjacency.
 */
inline void test_value_noise() {
  // Integer coordinates sample the lattice hash exactly.
  HS_EXPECT_NEAR(value_noise_1d(3.0f, 5u), hash01(3u, 5u), 1e-6f);
  HS_EXPECT_NEAR(value_noise_1d(-2.0f, 5u),
                 hash01(static_cast<uint32_t>(-2), 5u), 1e-6f);

  // Range over positive and negative coordinates.
  for (int i = -50; i <= 50; ++i) {
    float x = i * 0.173f;
    float n = value_noise_1d(x, 9u);
    HS_EXPECT_GE(n, 0.0f);
    HS_EXPECT_LT(n, 1.0f);
    float n2 = value_noise_2d(x, x * 0.7f, 9u);
    HS_EXPECT_GE(n2, 0.0f);
    HS_EXPECT_LT(n2, 1.0f);
  }

  // Continuity: adjacent samples differ by a bounded step, including at 0.
  for (int i = -400; i < 400; ++i) {
    float x = i * 0.005f;
    HS_EXPECT_LT(
        std::fabs(value_noise_1d(x + 0.005f, 3u) - value_noise_1d(x, 3u)),
        0.05f);
    HS_EXPECT_LT(std::fabs(value_noise_2d(x + 0.005f, 0.4f, 3u) -
                           value_noise_2d(x, 0.4f, 3u)),
                 0.05f);
    HS_EXPECT_LT(std::fabs(value_noise_2d(0.4f, x + 0.005f, 3u) -
                           value_noise_2d(0.4f, x, 3u)),
                 0.05f);
  }

  // Not constant, and seeds decorrelate.
  HS_EXPECT_TRUE(value_noise_1d(0.5f, 0u) != value_noise_1d(7.5f, 0u) ||
                 value_noise_1d(2.5f, 0u) != value_noise_1d(9.5f, 0u));
  HS_EXPECT_TRUE(value_noise_2d(0.5f, 0.5f, 1u) !=
                 value_noise_2d(0.5f, 0.5f, 2u));
}

inline void test_twist_lens() {
  Vector input(0.6f, 0.5f, 0.6244998f);
  const Vector output = lenses::twist_lens(input);
  const float angle = lenses::TWIST_RATE * input.y;
  HS_EXPECT_NEAR(output.x, input.x * cosf(angle) - input.z * sinf(angle),
                 2e-3f);
  HS_EXPECT_NEAR(output.y, input.y, 1e-6f);
  HS_EXPECT_NEAR(output.z, input.x * sinf(angle) + input.z * cosf(angle),
                 2e-3f);
}

inline void test_kaleidoscope_lens() {
  const Vector input = Vector(-0.3f, 0.4f, -0.8660254f).normalized();
  const Vector output = lenses::kaleidoscope_lens(input);
  HS_EXPECT_TRUE(output.x >= 0.0f);
  HS_EXPECT_TRUE(output.z >= 0.0f);
  HS_EXPECT_TRUE(1.7320508075688772f * output.z <= output.x + 1e-6f);
  HS_EXPECT_NEAR(output.y, input.y, 1e-6f);
  HS_EXPECT_NEAR(output.magnitude(), input.magnitude(), 1e-5f);
}

/**
 * @brief Invokes @p visit with each direction of a latitude-longitude grid
 *        covering the sphere, poles included.
 * @param visit Callable taking one unit Vector.
 */
template <typename Visit> inline void for_each_sphere_direction(Visit visit) {
  constexpr int LATITUDE_STEPS = 8;
  constexpr int LONGITUDE_STEPS = 29;
  for (int latitude_step = -LATITUDE_STEPS; latitude_step <= LATITUDE_STEPS;
       ++latitude_step) {
    const float latitude =
        latitude_step * (0.5f * PI_F / static_cast<float>(LATITUDE_STEPS));
    const float radius = cosf(latitude);
    for (int longitude_step = 0; longitude_step < LONGITUDE_STEPS;
         ++longitude_step) {
      const float longitude =
          longitude_step * (TWO_PI_F / static_cast<float>(LONGITUDE_STEPS));
      visit(Vector(radius * cosf(longitude), sinf(latitude),
                   radius * sinf(longitude)));
    }
  }
}

/**
 * @brief Checks a chamber fold over a whole-sphere direction grid: the result
 *        satisfies every mirror half-space, the fold is an isometry, and a
 *        direction already inside is a fixed point.
 * @param mirrors Inward unit normals bounding the chamber.
 * @param fold Chamber fold under test.
 * @details A mistyped normal either leaves the fold non-convergent, which trips
 * its own reflection-limit check, or opens the chamber past a mirror, which the
 * half-space assertions catch.
 */
template <typename Fold>
inline void expect_chamber_fold(const std::array<Vector, 3> &mirrors,
                                Fold fold) {
  constexpr float WALL_TOLERANCE = 1e-5f;
  for (const Vector &normal : mirrors)
    HS_EXPECT_NEAR(normal.magnitude(), 1.0f, 1e-6f);
  for_each_sphere_direction([&](const Vector &input) {
    const Vector folded = fold(input);
    for (const Vector &normal : mirrors)
      HS_EXPECT_GE(dot(folded, normal), -WALL_TOLERANCE);
    HS_EXPECT_NEAR(folded.magnitude(), input.magnitude(), 1e-4f);
    HS_EXPECT_VEC(fold(folded), folded, 1e-6f);
  });
}

/**
 * @brief Sweeps every reflection-group table wired into the lens catalog
 *        through the generic fold, plus the dodecahedral specialization.
 */
inline void test_polyhedral_kaleidoscope_chambers() {
  const std::array<Vector, 3> tables[] = {
      lenses::TETRAHEDRAL_MIRRORS,     lenses::OCTAHEDRAL_MIRRORS,
      lenses::DODECAHEDRAL_MIRRORS,    lenses::TRIANGULAR_PRISM_MIRRORS,
      lenses::SQUARE_PRISM_MIRRORS,    lenses::PENTAGONAL_PRISM_MIRRORS,
      lenses::HEXAGONAL_PRISM_MIRRORS, lenses::OCTAGONAL_PRISM_MIRRORS,
  };
  for (const std::array<Vector, 3> &mirrors : tables)
    expect_chamber_fold(mirrors, [&mirrors](const Vector &v) {
      return lenses::polyhedral_kaleidoscope_lens(v, mirrors);
    });
  expect_chamber_fold(lenses::DODECAHEDRAL_MIRRORS, [](const Vector &v) {
    return lenses::dodecahedral_kaleidoscope_lens(v);
  });
}

inline void test_dodecahedral_kaleidoscope_specialization() {
  for_each_sphere_direction([](const Vector &input) {
    const Vector generic = lenses::polyhedral_kaleidoscope_lens(
        input, lenses::DODECAHEDRAL_MIRRORS);
    const Vector specialized = lenses::dodecahedral_kaleidoscope_lens(input);
    HS_EXPECT_VEC(specialized, generic, 1e-5f);
  });
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
  test_hash01();
  test_value_noise();
  test_twist_lens();
  test_kaleidoscope_lens();
  test_polyhedral_kaleidoscope_chambers();
  test_dodecahedral_kaleidoscope_specialization();

  test_fast_atan2();
  test_diamond_angle();
  test_fast_reciprocal();
  test_fast_rsqrt();
  test_fast_acos();
  test_fast_sinf_cosf();
  test_fast_sincosf_0_pi();
  test_fast_cbrt();
  test_fast_cbrt3();
  test_fast_cbrt6();
  test_fast_expf();

  test_vector_construction();
  test_vector_spherical_construction();
  test_vector_equality();
  test_vector_arithmetic();
  test_vector_length();
  test_vector_normalize();
  test_snorm3_roundtrip_bound();

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

  test_make_rotation_axis_angle();
  test_least_parallel_axis();
  test_make_rotation_from_to();
  test_quaternion_from_basis();
  test_rotate();
  test_rotation_matrix_matches_rotate();

  test_vector_slerp();
  test_vector_nlerp_unit();
  test_vector_slerp_antipodal_monotonic();
  test_quaternion_slerp();
  test_scaled_rotation_delta();

  test_stereo_roundtrip();
  test_complex_arithmetic();
  test_mobius_params_accessors();
  test_mobius_transform();
  test_gnomonic_roundtrip();

  test_wrap_index();

  return fixture.result();
}

} // namespace math3d_tests
} // namespace hs_test
