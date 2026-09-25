/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/math/3dmath.h and core/math/4dmath.h.
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
#include "core/math/4dmath.h"
#include "core/math/lenses.h"
#include "core/math/mobius.h"
#include "core/math/projection_patterns.h"
#include "core/math/rotate.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/vec_test_util.h"

#include <algorithm>
#include <array>

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
  HS_EXPECT_NEAR(math::PHI, 1.61803398f, 1e-6f);
  HS_EXPECT_NEAR(math::INV_PHI, 0.6180339887f, 1e-7f);
  HS_EXPECT_NEAR(math::TOLERANCE, 0.0001f, 1e-9f);
  HS_EXPECT_NEAR(math::PI_F, 3.14159265f, 1e-5f);
  HS_EXPECT_NEAR_REL(projections::STEREO_INF, 1e4f, 1e-7f);
}

// ============================================================================
// quintic_kernel
// ============================================================================

/**
 * @brief Verifies the smootherstep kernel: fixed points, clamping outside
 *        [0,1], monotonicity, and flat (C2) endpoints.
 */
inline void test_quintic_kernel() {
  HS_EXPECT_NEAR(math::quintic_kernel(0.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(math::quintic_kernel(1.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::quintic_kernel(0.5f), 0.5f, 1e-6f);

  HS_EXPECT_NEAR(math::quintic_kernel(-1.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(math::quintic_kernel(2.0f), 1.0f, 1e-6f);

  float prev = -1.0f;
  for (int i = 0; i <= 200; ++i) {
    float v = math::quintic_kernel(i / 200.0f);
    HS_EXPECT_TRUE(v >= prev);
    prev = v;
  }

  HS_EXPECT_NEAR(math::quintic_kernel(0.25f), 0.103515625f, 1e-7f);
  HS_EXPECT_NEAR(math::quintic_kernel(0.75f), 0.896484375f, 1e-7f);

  // Quintic endpoint increments are cubic in the step size.
  float dl = math::quintic_kernel(0.01f) - math::quintic_kernel(0.0f);
  float dr = math::quintic_kernel(1.0f) - math::quintic_kernel(0.99f);
  HS_EXPECT_TRUE(std::abs(dl) < 2e-5f);
  HS_EXPECT_TRUE(std::abs(dr) < 2e-5f);
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
  HS_EXPECT_NEAR(math::fast_atan2(0.0f, 1.0f), 0.0f, 4e-3f);
  HS_EXPECT_NEAR(math::fast_atan2(1.0f, 0.0f), math::PI_F * 0.5f, 4e-3f);
  HS_EXPECT_NEAR(math::fast_atan2(0.0f, -1.0f), math::PI_F, 4e-3f);
  HS_EXPECT_NEAR(math::fast_atan2(-1.0f, 0.0f), -math::PI_F * 0.5f, 4e-3f);

  for (int i = 0; i < 64; ++i) {
    float a = -math::PI_F + (i * 2.0f * math::PI_F) / 64.0f;
    float y = std::sin(a);
    float x = std::cos(a);
    HS_EXPECT_NEAR(math::fast_atan2(y, x), std::atan2(y, x), 4e-3f);
  }

  // Peak error (~3.76e-3) near a = -2.5702 rad, between the sweep's samples.
  {
    float a = -2.5702f, y = std::sin(a), x = std::cos(a);
    HS_EXPECT_NEAR(math::fast_atan2(y, x), std::atan2(y, x), 4e-3f);
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
  HS_EXPECT_NEAR(math::diamond_angle(0.0f, 1.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(math::diamond_angle(1.0f, 0.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::diamond_angle(0.0f, -1.0f), 2.0f, 1e-6f);
  HS_EXPECT_NEAR(math::diamond_angle(-1.0f, 0.0f), 3.0f, 1e-6f);

  // Degenerate origin.
  HS_EXPECT_NEAR(math::diamond_angle(0.0f, 0.0f), 0.0f, 1e-6f);

  // Scale invariance: only the direction matters.
  for (int i = 0; i < 32; ++i) {
    float a = (i * 2.0f * math::PI_F) / 32.0f;
    float y = std::sin(a), x = std::cos(a);
    HS_EXPECT_NEAR(math::diamond_angle(y * 1e4f, x * 1e4f),
                   math::diamond_angle(y, x), 1e-5f);
    HS_EXPECT_NEAR(math::diamond_angle(y * 1e-4f, x * 1e-4f),
                   math::diamond_angle(y, x), 1e-5f);
  }

  // Strictly increasing over one counter-clockwise turn, always in [0, 4).
  float prev = -1.0f;
  for (int i = 0; i < 1024; ++i) {
    float a = (i * 2.0f * math::PI_F) / 1024.0f;
    float d = math::diamond_angle(std::sin(a), std::cos(a));
    HS_EXPECT_TRUE(d >= 0.0f && d < 4.0f);
    HS_EXPECT_TRUE(d > prev);
    prev = d;
  }

  // Fourth-quadrant seam: y -> 0- with x > 0 must stay inside [0, 4).
  for (float y : {-1e-3f, -1e-6f, -1e-9f, -1e-12f, -1e-20f, -1e-30f}) {
    float d = math::diamond_angle(y, 1.0f);
    HS_EXPECT_TRUE(d >= 0.0f && d < 4.0f);
  }
}

/**
 * @brief Verifies fast_reciprocal across normal inputs and outputs at its
 * documented peak relative error.
 */
inline void test_fast_reciprocal() {
  HS_EXPECT_NEAR(math::fast_reciprocal(1.0f), 1.0f, 7e-6f);
  HS_EXPECT_NEAR(math::fast_reciprocal(4.0f), 0.25f, 7e-6f * 0.25f);
  HS_EXPECT_NEAR(math::fast_reciprocal(0.25f), 4.0f, 7e-6f * 4.0f);

  for (int exponent = -120; exponent <= 120; ++exponent) {
    for (int mantissa = 0; mantissa < 8; ++mantissa) {
      const float x = std::ldexp(1.0f + 0.125f * mantissa, exponent);
      const float reference = 1.0f / x;
      HS_EXPECT_TRUE(
          std::abs(math::fast_reciprocal(x) - reference) / reference <= 7e-6f);
    }
  }
}

/**
 * @brief Verifies fast_rsqrt against 1/sqrt over a wide sweep at the documented
 *        ~5e-6 peak relative error.
 */
inline void test_fast_rsqrt() {
  HS_EXPECT_NEAR(math::fast_rsqrt(1.0f), 1.0f, 5e-6f);
  HS_EXPECT_NEAR(math::fast_rsqrt(4.0f), 0.5f, 5e-6f * 0.5f);
  HS_EXPECT_NEAR(math::fast_rsqrt(0.25f), 2.0f, 5e-6f * 2.0f);

  // Both exponent parities across ~12 decades: the bit-hack seed's quality
  // alternates with the low exponent bit, so a one-decade sweep would miss half
  // the error surface.
  for (int i = 0; i <= 512; ++i) {
    float x = std::pow(10.0f, -6.0f + (12.0f * i) / 512.0f);
    float ref = 1.0f / std::sqrt(x);
    HS_EXPECT_TRUE(std::abs(math::fast_rsqrt(x) - ref) / ref <= 5e-6f);
  }
}

/**
 * @brief Verifies fast_acos at endpoints, out-of-range clamping to [0,pi], and
 *        across a sweep against std::acos.
 */
inline void test_fast_acos() {
  HS_EXPECT_NEAR(math::fast_acos(1.0f), 0.0f, 2e-4f);
  HS_EXPECT_NEAR(math::fast_acos(-1.0f), math::PI_F, 2e-4f);
  HS_EXPECT_NEAR(math::fast_acos(0.0f), math::PI_F * 0.5f, 2e-4f);

  HS_EXPECT_NEAR(math::fast_acos(1.5f), 0.0f, 2e-4f);
  HS_EXPECT_NEAR(math::fast_acos(-1.5f), math::PI_F, 2e-4f);

  for (int i = 0; i <= 32; ++i) {
    float x = -1.0f + (i / 16.0f);
    HS_EXPECT_NEAR(math::fast_acos(x), std::acos(x), 2e-4f);
  }

  // Peak error (~5.0e-5) near x = 0.0807, between the sweep's samples.
  HS_EXPECT_NEAR(math::fast_acos(0.0807f), std::acos(0.0807f), 5.1e-5f);
}

/**
 * @brief Verifies fast_cbrt anchors, the x<=0 -> 0 clamp, and the documented
 *        ~2.3e-5 peak relative error over [0,8] (plus a few values past 8).
 */
inline void test_fast_cbrt() {
  HS_EXPECT_NEAR(math::fast_cbrt(-1.0f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(math::fast_cbrt(0.0f), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(math::fast_cbrt(1.0f), 1.0f, 2.3e-5f);
  HS_EXPECT_NEAR(math::fast_cbrt(8.0f), 2.0f, 2.3e-5f * 2.0f);

  for (int i = 1; i <= 256; ++i) {
    float x = (8.0f * i) / 256.0f;
    float rel = std::abs(math::fast_cbrt(x) - std::cbrt(x)) / std::cbrt(x);
    HS_EXPECT_TRUE(rel <= 2.3e-5f);
  }

  // Values past the documented [0,8] domain stay within the same rel error.
  for (float x : {27.0f, 100.0f, 1000.0f}) {
    float rel = std::abs(math::fast_cbrt(x) - std::cbrt(x)) / std::cbrt(x);
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
    math::fast_cbrt3(x[0], x[1], x[2], o[0], o[1], o[2]);
    for (int j = 0; j < 3; ++j) {
      float ref = math::fast_cbrt(x[j]);
      HS_EXPECT_TRUE(std::abs(o[j] - ref) / ref <= 1e-6f);
      float rel = std::abs(o[j] - std::cbrt(x[j])) / std::cbrt(x[j]);
      HS_EXPECT_TRUE(rel <= 2.3e-5f);
    }
  }

  // Non-positive inputs clamp to 0 without poisoning the shared product.
  {
    float o[3];
    math::fast_cbrt3(-1.0f, 0.0f, 8.0f, o[0], o[1], o[2]);
    HS_EXPECT_NEAR(o[0], 0.0f, 1e-7f);
    HS_EXPECT_NEAR(o[1], 0.0f, 1e-7f);
    HS_EXPECT_NEAR(o[2], 2.0f, 2.3e-5f * 2.0f);
  }

  // Just inside each end of the documented window the result still tracks the
  // scalar helper; the tiny end is below fast_cbrt's own accurate domain, so
  // only the re-association is under test there.
  for (float v : {1.0e11f, 1.0e-12f}) {
    float o[3];
    math::fast_cbrt3(v, v, v, o[0], o[1], o[2]);
    float ref = math::fast_cbrt(v);
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
    math::fast_cbrt6(x, o);
    for (int j = 0; j < 6; ++j) {
      float ref = math::fast_cbrt(x[j]);
      HS_EXPECT_TRUE(std::abs(o[j] - ref) / ref <= 1e-6f);
      float rel = std::abs(o[j] - std::cbrt(x[j])) / std::cbrt(x[j]);
      HS_EXPECT_TRUE(rel <= 2.3e-5f);
    }
  }

  // Non-positive inputs clamp to 0 without poisoning the shared product.
  {
    const float x[6] = {-1.0f, 0.0f, -0.0f, 8.0f, 1.0f, 27.0f};
    float o[6];
    math::fast_cbrt6(x, o);
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
    math::fast_cbrt6(x, o);
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
  HS_EXPECT_NEAR(math::fast_expf(0.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::fast_expf(-1.0f), std::exp(-1.0f), 7.5e-4f);
  HS_EXPECT_NEAR(math::fast_expf(-10.0f), std::exp(-10.0f),
                 7.5e-4f * std::exp(-10.0f));

  for (int i = 0; i <= 512; ++i) {
    float x = -30.0f + (30.0f * i) / 512.0f;
    float ref = std::exp(x);
    float rel = std::abs(math::fast_expf(x) - ref) / ref;
    HS_EXPECT_TRUE(rel <= 7.5e-4f);
  }

  HS_EXPECT_GT(math::fast_expf(-87.0f), 0.0f);
  HS_EXPECT_EQ(math::fast_expf(-88.0f), 0.0f);
  HS_EXPECT_EQ(math::fast_expf(-100.0f), 0.0f);
  HS_EXPECT_EQ(math::fast_expf(-200.0f), 0.0f);
}

/**
 * @brief Verifies fast_sinf/fast_cosf at key angles, the Pythagorean identity,
 *        and periodicity.
 * @details The periodicity check exercises range reduction beyond ±2π.
 */
inline void test_fast_sinf_cosf() {
  HS_EXPECT_NEAR(math::fast_sinf(0.0f), 0.0f, 1.8e-3f);
  HS_EXPECT_NEAR(math::fast_sinf(math::PI_F * 0.5f), 1.0f, 1.8e-3f);
  HS_EXPECT_NEAR(math::fast_sinf(math::PI_F), 0.0f, 1.8e-3f);
  HS_EXPECT_NEAR(math::fast_sinf(-math::PI_F * 0.5f), -1.0f, 1.8e-3f);

  HS_EXPECT_NEAR(math::fast_cosf(0.0f), 1.0f, 1.8e-3f);
  HS_EXPECT_NEAR(math::fast_cosf(math::PI_F * 0.5f), 0.0f, 1.8e-3f);
  HS_EXPECT_NEAR(math::fast_cosf(math::PI_F), -1.0f, 1.8e-3f);

  // Sweep peak (~1.63e-3) sits near a = -9.22, in the range-reduction band.
  for (int i = 0; i <= 256; ++i) {
    float a = -3.0f * math::PI_F + (i * 6.0f * math::PI_F) / 256.0f;
    HS_EXPECT_NEAR(math::fast_sinf(a), std::sin(a), 1.8e-3f);
    HS_EXPECT_NEAR(math::fast_cosf(a), std::cos(a), 1.8e-3f);
  }

  for (int i = 0; i < 32; ++i) {
    float a = -3.0f * math::PI_F + (i * 6.0f * math::PI_F) / 32.0f;
    float s = math::fast_sinf(a);
    float c = math::fast_cosf(a);
    HS_EXPECT_NEAR(s * s + c * c, 1.0f, 5e-3f);
  }

  for (int i = 0; i < 16; ++i) {
    float a = i * 0.3f;
    HS_EXPECT_NEAR(math::fast_sinf(a + 2.0f * math::PI_F), math::fast_sinf(a),
                   3e-3f);
    HS_EXPECT_NEAR(math::fast_sinf(a - 2.0f * math::PI_F), math::fast_sinf(a),
                   3e-3f);
  }
}

/**
 * @brief Verifies fast_sincosf_0_pi reproduces fast_sinf/fast_cosf bit for bit
 *        across its documented [0, pi] domain.
 */
inline void test_fast_sincosf_0_pi() {
  for (int i = 0; i <= 256; ++i) {
    float x = (i * math::PI_F) / 256.0f;
    float s, c;
    math::fast_sincosf_0_pi(x, s, c);
    HS_EXPECT_EQ(s, math::fast_sinf(x));
    HS_EXPECT_EQ(c, math::fast_cosf(x));
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
  math::Vector v0;
  HS_EXPECT_NEAR(v0.x, 0.0f, 1e-7f);
  HS_EXPECT_NEAR(v0.y, 0.0f, 1e-7f);
  HS_EXPECT_NEAR(v0.z, 0.0f, 1e-7f);

  math::Vector v_int(0); // inplace_function compat constructor
  HS_EXPECT_VEC(v_int, math::Vector(0, 0, 0), 1e-7f);

  math::Vector v(1.0f, 2.0f, 3.0f);
  HS_EXPECT_NEAR(v.x, 1.0f, 0.0f);
  HS_EXPECT_NEAR(v.y, 2.0f, 0.0f);
  HS_EXPECT_NEAR(v.z, 3.0f, 0.0f);

  math::Vector vc(v);
  HS_EXPECT_VEC(vc, v, 0.0f);

  math::Vector vassigned;
  vassigned = v;
  HS_EXPECT_VEC(vassigned, v, 0.0f);
}

/**
 * @brief Verifies Vector::from_spherical axis mappings, poles, and unit-length
 *        preservation.
 */
inline void test_vector_spherical_construction() {
  // theta=0, phi=π/2 → +X axis
  HS_EXPECT_VEC(math::Vector::from_spherical(0.0f, math::PI_F * 0.5f),
                math::Vector(1, 0, 0), 2e-3f);
  // theta=π/2, phi=π/2 → +Z axis
  HS_EXPECT_VEC(
      math::Vector::from_spherical(math::PI_F * 0.5f, math::PI_F * 0.5f),
      math::Vector(0, 0, 1), 2e-3f);
  // phi=0 → +Y (north pole)
  HS_EXPECT_VEC(math::Vector::from_spherical(0.0f, 0.0f), math::Vector(0, 1, 0),
                2e-3f);
  // phi=π → -Y (south pole)
  HS_EXPECT_VEC(math::Vector::from_spherical(0.0f, math::PI_F),
                math::Vector(0, -1, 0), 2e-3f);

  math::Vector u = math::Vector::from_spherical(0.8f, 1.1f);
  HS_EXPECT_NEAR(u.length(), 1.0f, 2e-3f);
}

/**
 * @brief Verifies Vector ==/!= tolerant comparison (equal within TOLERANCE,
 *        unequal beyond).
 */
inline void test_vector_equality() {
  math::Vector a(1, 2, 3), b(1, 2, 3), c(1.001f, 2, 3);
  HS_EXPECT_EQ(a, b);
  HS_EXPECT_FALSE(a == c);
  HS_EXPECT_TRUE(a != c);
  HS_EXPECT_FALSE(a != b);

  HS_EXPECT_EQ(math::Vector(1, 0, 0), math::Vector(1.00005f, 0, 0));
  HS_EXPECT_FALSE(math::Vector(1, 0, 0) == math::Vector(1.001f, 0, 0));
}

// ============================================================================
// Vector — arithmetic
// ============================================================================

/**
 * @brief Verifies Vector +, -, negate, scalar * and / (both orders), and
 *        compound assignment.
 */
inline void test_vector_arithmetic() {
  math::Vector a(1, 2, 3), b(4, 5, 6);
  HS_EXPECT_VEC(a + b, math::Vector(5, 7, 9), 1e-6f);
  HS_EXPECT_VEC(b - a, math::Vector(3, 3, 3), 1e-6f);
  HS_EXPECT_VEC(-a, math::Vector(-1, -2, -3), 1e-6f);
  HS_EXPECT_VEC(a * 2.0f, math::Vector(2, 4, 6), 1e-6f);
  HS_EXPECT_VEC(2.0f * a, math::Vector(2, 4, 6), 1e-6f);
  HS_EXPECT_VEC(a / 2.0f, math::Vector(0.5f, 1.0f, 1.5f), 1e-6f);

  math::Vector c(1, 2, 3);
  c += b;
  HS_EXPECT_VEC(c, math::Vector(5, 7, 9), 1e-6f);
  c -= b;
  HS_EXPECT_VEC(c, math::Vector(1, 2, 3), 1e-6f);
  c *= 3.0f;
  HS_EXPECT_VEC(c, math::Vector(3, 6, 9), 1e-6f);
  c /= 3.0f;
  HS_EXPECT_VEC(c, math::Vector(1, 2, 3), 1e-6f);
}

/**
 * @brief Verifies Vector length()/magnitude() Euclidean norm.
 */
inline void test_vector_length() {
  HS_EXPECT_NEAR(math::Vector(3, 4, 0).length(), 5.0f, 1e-6f);
  HS_EXPECT_NEAR(math::Vector(0, 0, 0).length(), 0.0f, 1e-7f);
  HS_EXPECT_NEAR(math::Vector(1, 2, 2).magnitude(), 3.0f, 1e-6f);
}

/**
 * @brief Verifies Vector normalize()/normalized() and the normalized_or()
 *        fallback.
 * @details The zero vector is rejected by the strict normalize() path, so
 *          normalized_or() supplies a fallback for the degenerate case.
 */
inline void test_vector_normalize() {
  math::Vector v(3, 0, 4);
  v.normalize();
  HS_EXPECT_NEAR(v.length(), 1.0f, 1e-6f);
  HS_EXPECT_VEC(v, math::Vector(0.6f, 0.0f, 0.8f), 1e-6f);

  math::Vector u(0, 5, 0);
  math::Vector n = u.normalized();
  HS_EXPECT_VEC(u, math::Vector(0, 5, 0), 1e-6f);
  HS_EXPECT_VEC(n, math::Vector(0, 1, 0), 1e-6f);

  // normalized_or() returns the fallback for a zero-length input (which traps
  // under strict normalize()), else normalizes as usual.
  HS_EXPECT_VEC(
      math::normalized_or(math::Vector(0, 0, 0), math::Vector(1, 0, 0)),
      math::Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(
      math::normalized_or(math::Vector(0, 6, 0), math::Vector(1, 0, 0)),
      math::Vector(0, 1, 0), 1e-6f);
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

  const math::Snorm3 endpoints = math::Snorm3::encode(math::Vector(1, 0, -1));
  HS_EXPECT_EQ(endpoints.x, 32767);
  HS_EXPECT_EQ(endpoints.y, 0);
  HS_EXPECT_EQ(endpoints.z, -32767);
  HS_EXPECT_VEC(endpoints.decode(), math::Vector(1, 0, -1), COMPONENT_BOUND);

  const math::Snorm3 saturated =
      math::Snorm3::encode(math::Vector(4.0f, -9.0f, 1.0f + 1e-3f));
  HS_EXPECT_EQ(saturated.x, 32767);
  HS_EXPECT_EQ(saturated.y, -32767);
  HS_EXPECT_EQ(saturated.z, 32767);

  hs::Pcg32 rng(20260803u);
  for (int i = 0; i < 256; ++i) {
    const math::Vector v =
        math::normalized_or(math::Vector(rand_uniform(rng, -1.0f, 1.0f),
                                         rand_uniform(rng, -1.0f, 1.0f),
                                         rand_uniform(rng, -1.0f, 1.0f)),
                            math::Vector(1, 0, 0));
    const math::Vector decoded = math::Snorm3::encode(v).decode();
    HS_EXPECT_LE(std::abs(decoded.x - v.x), COMPONENT_BOUND);
    HS_EXPECT_LE(std::abs(decoded.y - v.y), COMPONENT_BOUND);
    HS_EXPECT_LE(std::abs(decoded.z - v.z), COMPONENT_BOUND);
    HS_EXPECT_LE(math::distance_between(decoded, v), CHORD_BOUND);
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
  math::Vector x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);

  HS_EXPECT_NEAR(math::dot(x, y), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(math::dot(x, x), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::dot(x, -x), -1.0f, 1e-6f);
  HS_EXPECT_NEAR(math::dot(math::Vector(1, 2, 3), math::Vector(4, -5, 6)),
                 4 - 10 + 18, 1e-5f);

  HS_EXPECT_VEC(math::cross(x, y), z, 1e-6f);
  HS_EXPECT_VEC(math::cross(y, z), x, 1e-6f);
  HS_EXPECT_VEC(math::cross(z, x), y, 1e-6f);
  HS_EXPECT_VEC(math::cross(y, x), -z, 1e-6f);
  HS_EXPECT_VEC(math::cross(math::Vector(2, 3, 5), math::Vector(2, 3, 5)),
                math::Vector(0, 0, 0), 1e-6f);
  // a × b is perpendicular to both operands.
  math::Vector a(1, 2, 3), b(4, -5, 6);
  math::Vector c = math::cross(a, b);
  HS_EXPECT_NEAR(math::dot(c, a), 0.0f, 1e-4f);
  HS_EXPECT_NEAR(math::dot(c, b), 0.0f, 1e-4f);
}

/**
 * @brief Verifies distance_between/distance_squared, including the
 *        coincident-point zero.
 */
inline void test_distance() {
  math::Vector a(1, 2, 3), b(4, 6, 3);
  HS_EXPECT_NEAR(math::distance_between(a, b), 5.0f, 1e-6f);
  HS_EXPECT_NEAR(math::distance_squared(a, b), 25.0f, 1e-5f);
  HS_EXPECT_NEAR(math::distance_between(a, a), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(math::distance_squared(a, a), 0.0f, 1e-6f);
}

/**
 * @brief Verifies angle_between for vectors: cardinal angles and
 *        magnitude-independence.
 */
inline void test_angle_between_vectors() {
  math::Vector x(1, 0, 0), y(0, 1, 0);
  HS_EXPECT_NEAR(math::angle_between(x, y), math::PI_F * 0.5f, 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(x, x), 0.0f, 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(x, -x), math::PI_F, 1e-3f);
  // Independent of operand magnitude.
  HS_EXPECT_NEAR(math::angle_between(x * 5.0f, y * 0.3f), math::PI_F * 0.5f,
                 1e-3f);
}

// ============================================================================
// Spherical
// ============================================================================

/**
 * @brief Verifies Spherical accessors and Vector<->Spherical roundtrips (on-
 *        and off-equator).
 */
inline void test_spherical() {
  math::Spherical s(0.5f, 1.2f);
  HS_EXPECT_NEAR(s.theta, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(s.phi, 1.2f, 1e-6f);

  math::Vector v_orig(0.6f, 0.0f, 0.8f);
  math::Spherical s2(v_orig);
  math::Vector v2(s2);
  HS_EXPECT_VEC(v2, v_orig, 5e-3f);

  math::Vector v3 = math::Vector(1.0f, 1.0f, 1.0f).normalized();
  math::Spherical s3(v3);
  math::Vector v3_back(s3);
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
  math::Quaternion id;
  HS_EXPECT_NEAR(id.r, 1.0f, 1e-7f);
  HS_EXPECT_VEC(id.v, math::Vector(0, 0, 0), 1e-7f);

  math::Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);
  HS_EXPECT_NEAR(q.r, 0.5f, 0.0f);
  HS_EXPECT_VEC(q.v, math::Vector(0.5f, 0.5f, 0.5f), 0.0f);

  math::Quaternion qv(0.7f, math::Vector(0.1f, 0.2f, 0.3f));
  HS_EXPECT_NEAR(qv.r, 0.7f, 0.0f);
  HS_EXPECT_VEC(qv.v, math::Vector(0.1f, 0.2f, 0.3f), 0.0f);

  math::Quaternion qc(q);
  HS_EXPECT_QUAT(qc, q, 1e-7f);

  math::Quaternion qa;
  qa = q;
  HS_EXPECT_QUAT(qa, q, 1e-7f);
}

/**
 * @brief Verifies Quaternion +, -, scalar * and /, negate, and componentwise
 *        compound assignment.
 */
inline void test_quaternion_arithmetic() {
  math::Quaternion a(1, 2, 3, 4), b(0.5f, 1, 1.5f, 2);
  HS_EXPECT_QUAT(a + b, math::Quaternion(1.5f, 3, 4.5f, 6), 1e-6f);
  HS_EXPECT_QUAT(a - b, math::Quaternion(0.5f, 1, 1.5f, 2), 1e-6f);
  HS_EXPECT_QUAT(a * 2.0f, math::Quaternion(2, 4, 6, 8), 1e-6f);
  HS_EXPECT_QUAT(2.0f * a, math::Quaternion(2, 4, 6, 8), 1e-6f);
  HS_EXPECT_QUAT(a / 2.0f, math::Quaternion(0.5f, 1, 1.5f, 2), 1e-6f);
  HS_EXPECT_QUAT(-a, math::Quaternion(-1, -2, -3, -4), 1e-6f);

  math::Quaternion c(1, 2, 3, 4);
  c += b;
  HS_EXPECT_QUAT(c, math::Quaternion(1.5f, 3, 4.5f, 6), 1e-6f);
  c -= b;
  HS_EXPECT_QUAT(c, math::Quaternion(1, 2, 3, 4), 1e-6f);
  c *= 0.5f;
  HS_EXPECT_QUAT(c, math::Quaternion(0.5f, 1, 1.5f, 2), 1e-6f);
}

/**
 * @brief Verifies Quaternion squared_magnitude()/magnitude().
 */
inline void test_quaternion_magnitude() {
  math::Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);
  HS_EXPECT_NEAR(q.squared_magnitude(), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(q.magnitude(), 1.0f, 1e-6f);

  math::Quaternion p(1, 2, 2, 0);
  HS_EXPECT_NEAR(p.squared_magnitude(), 9.0f, 1e-5f);
  HS_EXPECT_NEAR(p.magnitude(), 3.0f, 1e-6f);
}

/**
 * @brief Verifies conjugate() and inverse().
 * @details For a unit q the two coincide; q*q^-1 = identity holds for unit
 *          and non-unit quaternions alike.
 */
inline void test_quaternion_conjugate_inverse() {
  math::Quaternion q(0.5f, 0.5f, 0.5f, 0.5f); // unit
  math::Quaternion conj = q.conjugate();
  HS_EXPECT_QUAT(conj, math::Quaternion(0.5f, -0.5f, -0.5f, -0.5f), 1e-6f);

  HS_EXPECT_QUAT(q.inverse(), conj, 1e-6f);

  HS_EXPECT_QUAT(q * q.inverse(), math::Quaternion(1, 0, 0, 0), 1e-6f);

  // inverse() inverts a non-unit quaternion too.
  math::Quaternion p(2, 0, 0, 0);
  HS_EXPECT_QUAT(p * p.inverse(), math::Quaternion(1, 0, 0, 0), 1e-6f);
}

/**
 * @brief Verifies Quaternion normalize() (in place) and normalized()
 *        (non-mutating).
 */
inline void test_quaternion_normalize() {
  math::Quaternion p(2, 0, 0, 0);
  p.normalize();
  HS_EXPECT_NEAR(p.magnitude(), 1.0f, 1e-6f);
  HS_EXPECT_QUAT(p, math::Quaternion(1, 0, 0, 0), 1e-6f);

  math::Quaternion u(3, 0, 0, 0);
  math::Quaternion n = u.normalized();
  HS_EXPECT_NEAR(u.r, 3.0f, 1e-6f);
  HS_EXPECT_QUAT(n, math::Quaternion(1, 0, 0, 0), 1e-6f);
}

/**
 * @brief Verifies the Hamilton product: identity laws, basis relations
 *        (i²=j²=k²=ijk=-1, cyclic products), non-commutativity, and *=
 *        consistency.
 */
inline void test_quaternion_multiplication() {
  math::Quaternion id;
  math::Quaternion q(0.5f, 0.5f, 0.5f, 0.5f);

  HS_EXPECT_QUAT(id * q, q, 1e-6f);
  HS_EXPECT_QUAT(q * id, q, 1e-6f);

  // Hamilton basis: i² = j² = k² = ijk = -1
  math::Quaternion i(0, 1, 0, 0), j(0, 0, 1, 0), k(0, 0, 0, 1);
  math::Quaternion neg_one(-1, 0, 0, 0);
  HS_EXPECT_QUAT(i * i, neg_one, 1e-6f);
  HS_EXPECT_QUAT(j * j, neg_one, 1e-6f);
  HS_EXPECT_QUAT(k * k, neg_one, 1e-6f);
  HS_EXPECT_QUAT(i * j * k, neg_one, 1e-6f);

  HS_EXPECT_QUAT(i * j, k, 1e-6f);
  HS_EXPECT_QUAT(j * k, i, 1e-6f);
  HS_EXPECT_QUAT(k * i, j, 1e-6f);
  HS_EXPECT_QUAT(j * i, -k, 1e-6f);

  math::Quaternion qa(0.5f, 0.5f, 0.5f, 0.5f), qb(qa);
  qa *= qa;
  HS_EXPECT_QUAT(qa, qb * qb, 1e-6f);
}

/**
 * @brief Verifies Quaternion == tolerant comparison around TOLERANCE.
 */
inline void test_quaternion_equality() {
  math::Quaternion a(1, 2, 3, 4), b(1, 2, 3, 4);
  HS_EXPECT_EQ(a, b);
  math::Quaternion c(1.001f, 2, 3, 4); // beyond TOLERANCE
  HS_EXPECT_FALSE(a == c);
  math::Quaternion d(1.00001f, 2, 3, 4); // within TOLERANCE
  HS_EXPECT_EQ(a, d);
}

/**
 * @brief Verifies the 4-component dot product on quaternions; self-dot equals
 *        squared magnitude.
 */
inline void test_dot_quaternion() {
  math::Quaternion a(1, 2, 3, 4), b(2, 3, 4, 5);
  HS_EXPECT_NEAR(math::dot(a, b), 1 * 2 + 2 * 3 + 3 * 4 + 4 * 5, 1e-5f);
  HS_EXPECT_NEAR(math::dot(a, a), a.squared_magnitude(), 1e-5f);
}

// ============================================================================
// make_rotation / rotate
// ============================================================================

/**
 * @brief Verifies make_rotation(axis, angle): identity at angle 0, unit-length
 *        result, and correct right-handed rotation of test vectors.
 */
inline void test_make_rotation_axis_angle() {
  math::Quaternion id = math::make_rotation(math::Vector(0, 1, 0), 0.0f);
  HS_EXPECT_NEAR(std::abs(id.r), 1.0f, 5e-3f);
  HS_EXPECT_VEC(id.v, math::Vector(0, 0, 0), 5e-3f);

  math::Quaternion qy90 =
      math::make_rotation(math::Vector(0, 1, 0), math::PI_F * 0.5f);
  HS_EXPECT_NEAR(qy90.magnitude(), 1.0f, 1e-4f);

  // 90° around +Y rotates (1,0,0) to (0,0,-1) [right-handed]
  HS_EXPECT_VEC(math::rotate(math::Vector(1, 0, 0), qy90),
                math::Vector(0, 0, -1), 5e-3f);

  // 180° around +Z rotates (1,0,0) to (-1,0,0)
  math::Quaternion qz180 =
      math::make_rotation(math::Vector(0, 0, 1), math::PI_F);
  HS_EXPECT_VEC(math::rotate(math::Vector(1, 0, 0), qz180),
                math::Vector(-1, 0, 0), 5e-3f);
}

/**
 * @brief Verifies least_parallel_axis picks +Y only near +/-X, is scale
 *        invariant, and always seeds a well-conditioned cross, and that
 *        perpendicular_axis turns that seed into a unit tangent.
 */
inline void test_least_parallel_axis() {
  HS_EXPECT_VEC(math::least_parallel_axis(math::Vector(1, 0, 0)),
                math::Vector(0, 1, 0), 1e-6f);
  HS_EXPECT_VEC(math::least_parallel_axis(math::Vector(-1, 0, 0)),
                math::Vector(0, 1, 0), 1e-6f);
  HS_EXPECT_VEC(math::least_parallel_axis(math::Vector(0, 1, 0)),
                math::Vector(1, 0, 0), 1e-6f);
  HS_EXPECT_VEC(math::least_parallel_axis(math::Vector(0, 0, 1)),
                math::Vector(1, 0, 0), 1e-6f);

  // Scale invariant: the cosine test is against |v|^2, not the raw component.
  HS_EXPECT_VEC(math::least_parallel_axis(math::Vector(50, 0, 0)),
                math::Vector(0, 1, 0), 1e-6f);
  HS_EXPECT_VEC(math::least_parallel_axis(math::Vector(0, 50, 0)),
                math::Vector(1, 0, 0), 1e-6f);

  // The whole point: cross(axis, v) never collapses. The worst case sits just
  // inside the COS_AXIS_PARALLEL switch, where sin^2 is ~2*TOLERANCE.
  for (int i = 0; i <= 128; ++i) {
    for (int j = 0; j <= 128; ++j) {
      float phi = (i * math::PI_F) / 128.0f;
      float theta = (j * 2.0f * math::PI_F) / 128.0f;
      math::Vector v(std::sin(phi) * std::cos(theta), std::cos(phi),
                     std::sin(phi) * std::sin(theta));
      math::Vector c = math::cross(math::least_parallel_axis(v), v);
      HS_EXPECT_TRUE(math::dot(c, c) >= 1e-4f);
      math::Vector t = math::perpendicular_axis(v);
      HS_EXPECT_NEAR(math::dot(t, t), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(math::dot(t, v), 0.0f, 1e-6f);
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
  math::Quaternion id =
      math::make_rotation(math::Vector(1, 0, 0), math::Vector(1, 0, 0));
  HS_EXPECT_QUAT(id, math::Quaternion(1, 0, 0, 0), 1e-4f);

  math::Quaternion q =
      math::make_rotation(math::Vector(1, 0, 0), math::Vector(0, 1, 0));
  HS_EXPECT_VEC(math::rotate(math::Vector(1, 0, 0), q), math::Vector(0, 1, 0),
                5e-3f);

  // Antiparallel (degenerate): x to -x → 180° rotation.
  math::Quaternion qa =
      math::make_rotation(math::Vector(1, 0, 0), math::Vector(-1, 0, 0));
  HS_EXPECT_VEC(math::rotate(math::Vector(1, 0, 0), qa), math::Vector(-1, 0, 0),
                5e-3f);

  math::Vector from(1, 1, 0);
  from.normalize();
  math::Vector to(0, 1, 1);
  to.normalize();
  math::Quaternion qg = math::make_rotation(from, to);
  HS_EXPECT_VEC(math::rotate(from, qg), to, 5e-3f);
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
  math::Quaternion id = math::quaternion_from_basis(
      math::Vector(1, 0, 0), math::Vector(0, 1, 0), math::Vector(0, 0, 1));
  HS_EXPECT_VEC(math::rotate(math::Vector(1, 0, 0), id), math::Vector(1, 0, 0),
                1e-5f);
  HS_EXPECT_VEC(math::rotate(math::Vector(0, 0, 1), id), math::Vector(0, 0, 1),
                1e-5f);

  math::Quaternion q0 =
      math::make_rotation(math::Vector(0.3f, 0.5f, 0.8f).normalized(), 1.1f);
  math::Vector cx = math::rotate(math::Vector(1, 0, 0), q0);
  math::Vector cy = math::rotate(math::Vector(0, 1, 0), q0);
  math::Vector cz = math::rotate(math::Vector(0, 0, 1), q0);
  math::Quaternion q = math::quaternion_from_basis(cx, cy, cz);
  HS_EXPECT_VEC(math::rotate(math::Vector(1, 0, 0), q), cx, 5e-3f);
  HS_EXPECT_VEC(math::rotate(math::Vector(0, 1, 0), q), cy, 5e-3f);
  HS_EXPECT_VEC(math::rotate(math::Vector(0, 0, 1), q), cz, 5e-3f);
  HS_EXPECT_NEAR(q.magnitude(), 1.0f, 1e-3f);

  // trace <= 0 branch: 180° rotation about Z (diagonal = (-1,-1,1)).
  math::Quaternion qz = math::quaternion_from_basis(
      math::Vector(-1, 0, 0), math::Vector(0, -1, 0), math::Vector(0, 0, 1));
  HS_EXPECT_VEC(math::rotate(math::Vector(1, 0, 0), qz), math::Vector(-1, 0, 0),
                5e-3f);
  HS_EXPECT_VEC(math::rotate(math::Vector(0, 0, 1), qz), math::Vector(0, 0, 1),
                5e-3f);
}

/**
 * @brief Verifies rotate(v, q): identity, length preservation, and the
 *        composition law rotate(rotate(v,q1),q2) == rotate(v, q2*q1).
 */
inline void test_rotate() {
  math::Vector v(1, 2, 3);
  HS_EXPECT_VEC(math::rotate(v, math::Quaternion(1, 0, 0, 0)), v, 1e-6f);

  math::Quaternion q =
      math::make_rotation(math::Vector(1, 2, 3).normalized(), 1.234f);
  math::Vector r = math::rotate(v, q);
  HS_EXPECT_NEAR(r.length(), v.length(), 1e-3f);

  // Composition: rotate(rotate(v, q1), q2) == rotate(v, q2 * q1).
  math::Quaternion q1 = math::make_rotation(math::Vector(0, 1, 0), 0.3f);
  math::Quaternion q2 = math::make_rotation(math::Vector(1, 0, 0), 0.5f);
  math::Vector via_sequential = math::rotate(math::rotate(v, q1), q2);
  math::Vector via_composed = math::rotate(v, q2 * q1);
  HS_EXPECT_VEC(via_sequential, via_composed, 5e-3f);
}

/**
 * @brief Pins RotationMatrix to rotate(): the expanded rows must reproduce the
 *        quaternion sandwich for every orientation, and the rows must stay
 *        orthonormal.
 */
inline void test_rotation_matrix_matches_rotate() {
  const math::Vector axes[] = {math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                               math::Vector(0, 0, 1),
                               math::Vector(1, 2, 3).normalized(),
                               math::Vector(-2, 0.5f, 1).normalized()};
  const math::Vector samples[] = {math::Vector(1, 0, 0), math::Vector(0, 1, 0),
                                  math::Vector(0, 0, 1), math::Vector(1, 2, 3),
                                  math::Vector(-4, 0.25f, 2)};

  for (const math::Vector &axis : axes) {
    for (int i = 0; i <= 12; ++i) {
      float angle = (2.0f * math::PI_F * i) / 12.0f;
      math::Quaternion q = math::make_rotation(axis, angle);
      math::RotationMatrix m(q);
      for (const math::Vector &v : samples) {
        HS_EXPECT_VEC(m.apply(v), math::rotate(v, q), 1e-4f);
      }
      HS_EXPECT_NEAR(math::dot(m.r0, m.r0), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(math::dot(m.r1, m.r1), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(math::dot(m.r2, m.r2), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(math::dot(m.r0, m.r1), 0.0f, 1e-5f);
      HS_EXPECT_NEAR(math::dot(m.r0, m.r2), 0.0f, 1e-5f);
      HS_EXPECT_NEAR(math::dot(m.r1, m.r2), 0.0f, 1e-5f);
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
  math::Vector a(1, 0, 0), b(0, 1, 0);
  HS_EXPECT_VEC(math::slerp(a, b, 0.0f), a, 5e-3f);
  HS_EXPECT_VEC(math::slerp(a, b, 1.0f), b, 5e-3f);

  // Midpoint on unit sphere is (√2/2, √2/2, 0)
  math::Vector mid = math::slerp(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.length(), 1.0f, 1e-3f);
  HS_EXPECT_NEAR(mid.x, std::sqrt(2.0f) * 0.5f, 5e-3f);
  HS_EXPECT_NEAR(mid.y, std::sqrt(2.0f) * 0.5f, 5e-3f);
  HS_EXPECT_NEAR(mid.z, 0.0f, 5e-3f);

  // Nearly-identical vectors take the lerp fallback; result stays unit-length.
  math::Vector v1(1, 0, 0);
  math::Vector v2(0.99999f, 0.00001f, 0.0f);
  math::Vector lerp_result = math::slerp(v1, v2, 0.5f);
  HS_EXPECT_NEAR(lerp_result.length(), 1.0f, 1e-3f);

  // Antipodal endpoints: the great-circle direction is undefined, so slerp picks
  // a perpendicular axis and sweeps a monotone half-turn — the midpoint must NOT
  // collapse back onto p.
  math::Vector p(0, 1, 0), ap(0, -1, 0);
  HS_EXPECT_VEC(math::slerp(p, ap, 0.0f), p, 5e-3f);
  HS_EXPECT_VEC(math::slerp(p, ap, 1.0f), ap, 5e-3f);
  float a25 = math::dot(math::slerp(p, ap, 0.25f), p);
  float a50 = math::dot(math::slerp(p, ap, 0.50f), p);
  float a75 = math::dot(math::slerp(p, ap, 0.75f), p);
  HS_EXPECT_NEAR(math::slerp(p, ap, 0.5f).length(), 1.0f, 1e-3f);
  HS_EXPECT_GT(a25, a50);
  HS_EXPECT_GT(a50, a75);
  HS_EXPECT_NEAR(a50, 0.0f, 5e-3f);
}

/**
 * @brief Verifies nlerp_unit: endpoints, the unit-length midpoint, and the
 *        cancelling-blend fallback at math::EPS_BLEND_LEN_SQ.
 */
inline void test_vector_nlerp_unit() {
  const math::Vector a(1, 0, 0), b(0, 1, 0);
  HS_EXPECT_VEC(math::nlerp_unit(a, b, 0.0f), a, 1e-5f);
  HS_EXPECT_VEC(math::nlerp_unit(a, b, 1.0f), b, 1e-5f);

  const math::Vector mid = math::nlerp_unit(a, b, 0.5f);
  HS_EXPECT_NEAR(mid.x, std::sqrt(2.0f) * 0.5f, 1e-5f);
  HS_EXPECT_NEAR(mid.y, std::sqrt(2.0f) * 0.5f, 1e-5f);
  HS_EXPECT_NEAR(mid.z, 0.0f, 1e-5f);

  for (int i = 0; i <= 10; ++i)
    HS_EXPECT_NEAR(
        math::nlerp_unit(a, b, static_cast<float>(i) / 10.0f).length(), 1.0f,
        1e-5f);

  // Antipodal endpoints cancel: the blend carries no direction, so it holds a.
  const math::Vector p(0, 1, 0), ap(0, -1, 0);
  HS_EXPECT_VEC(math::nlerp_unit(p, ap, 0.5f), p, 1e-5f);
  HS_EXPECT_VEC(math::nlerp_unit(p, ap, 0.5f + 1e-5f), p, 1e-5f);
  HS_EXPECT_VEC(math::nlerp_unit(p, ap, 0.5f + 1e-3f), ap, 1e-5f);

  // Near-antipodal but non-cancelling: the residual still resolves to a unit
  // direction.
  const math::Vector near_ap = math::Vector(0.01f, -1.0f, 0.0f).normalized();
  HS_EXPECT_NEAR(math::nlerp_unit(p, near_ap, 0.5f).length(), 1.0f, 1e-4f);
}

/**
 * @brief Verifies Quaternion slerp: endpoints (q and -q are the same
 *        orientation), unit-length interpolants, the q^0.5-squared==q property,
 *        and the long_way (long-arc) variant.
 */
inline void test_quaternion_slerp() {
  math::Quaternion id(1, 0, 0, 0);
  math::Quaternion q =
      math::make_rotation(math::Vector(0, 1, 0), math::PI_F * 0.5f);

  HS_EXPECT_QUAT(math::slerp(id, q, 0.0f), id, 5e-3f);

  // t=1 may return -q; |dot| == 1 since q and -q are the same orientation.
  math::Quaternion s1 = math::slerp(id, q, 1.0f);
  HS_EXPECT_NEAR(std::abs(math::dot(s1, q)), 1.0f, 5e-3f);

  math::Quaternion half = math::slerp(id, q, 0.5f);
  HS_EXPECT_NEAR(half.magnitude(), 1.0f, 1e-3f);

  // Composing half with itself recovers q (q^0.5 squared = q).
  math::Vector v(1, 0, 0);
  math::Vector r_twice = math::rotate(math::rotate(v, half), half);
  math::Vector r_full = math::rotate(v, q);
  HS_EXPECT_VEC(r_twice, r_full, 1e-2f);

  // long_way negates the start: t=0 returns -id, and its midpoint differs from
  // the short-arc one.
  math::Quaternion long_start = math::slerp(id, q, 0.0f, true);
  HS_EXPECT_QUAT(long_start, -id, 5e-3f);

  math::Quaternion mid_short = math::slerp(id, q, 0.5f, false);
  math::Quaternion mid_long = math::slerp(id, q, 0.5f, true);
  HS_EXPECT_TRUE(std::isfinite(mid_short.magnitude()));
  HS_EXPECT_TRUE(std::isfinite(mid_long.magnitude()));
  HS_EXPECT_FALSE(approx_quat(mid_short, mid_long, 1e-2f));

  // Identical endpoints on the long arc drive d to -1 via the sign fixup; the
  // fallback must stay unit.
  math::Quaternion long_degenerate = math::slerp(q, q, 0.5f, true);
  HS_EXPECT_NEAR(long_degenerate.magnitude(), 1.0f, 1e-3f);
}

inline void test_scaled_rotation_delta() {
  const math::Quaternion id;
  const math::Quaternion q =
      math::make_rotation(math::Vector(0, 1, 0), math::PI_F * 0.5f);

  // Both extremes are exact and skip the slerp entirely.
  HS_EXPECT_QUAT(math::scaled_rotation_delta(q, 1.0f), q, 1e-6f);
  HS_EXPECT_QUAT(math::scaled_rotation_delta(q, 0.0f), id, 1e-6f);
  // Identity in, identity out, at any fraction.
  HS_EXPECT_QUAT(math::scaled_rotation_delta(id, 0.5f), id, 5e-3f);

  // Half the arc, applied twice, recovers the full delta.
  const math::Quaternion half = math::scaled_rotation_delta(q, 0.5f);
  HS_EXPECT_NEAR(half.magnitude(), 1.0f, 1e-3f);
  const math::Vector v(1, 0, 0);
  HS_EXPECT_VEC(math::rotate(math::rotate(v, half), half), math::rotate(v, q),
                1e-2f);

  // The scaled turn is monotone in amount: the angle away from the start grows.
  float previous = -1.0f;
  for (int step = 0; step <= 8; ++step) {
    const float amount = static_cast<float>(step) / 8.0f;
    const math::Quaternion scaled = math::scaled_rotation_delta(q, amount);
    HS_EXPECT_NEAR(scaled.magnitude(), 1.0f, 1e-3f);
    const float turned = 1.0f - math::dot(math::rotate(v, scaled), v);
    HS_EXPECT_GT(turned, previous);
    previous = turned;
  }
}

/**
 * @brief Verifies Vector slerp across antipodal endpoints sweeps monotonically:
 *        dot with the start strictly decreases as t increases, with no flip.
 */
inline void test_vector_slerp_antipodal_monotonic() {
  math::Vector p(0, 1, 0), ap(0, -1, 0);
  float prev = math::dot(math::slerp(p, ap, 0.0f), p);
  for (int i = 1; i <= 16; ++i) {
    float t = static_cast<float>(i) / 16.0f;
    math::Vector s = math::slerp(p, ap, t);
    HS_EXPECT_NEAR(s.length(), 1.0f, 1e-3f);
    float cur = math::dot(s, p);
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
  math::Vector samples[] = {
      math::Vector(1, 0, 0),
      math::Vector(0, 0, 1),
      math::Vector(-1, 0, 0),
      math::Vector(0, -1, 0),
      math::Vector(0.6f, 0.0f, 0.8f),
      math::Vector(0.5f, 0.5f, 0.7071f).normalized(),
  };
  for (const math::Vector &v : samples) {
    math::Complex z = projections::stereo(v);
    math::Vector back = projections::inv_stereo(z);
    HS_EXPECT_VEC(back, v, 5e-3f);
  }

  // North pole maps to the infinity sentinel (azimuth undefined → +real axis).
  math::Complex zN = projections::stereo(math::Vector(0, 1, 0));
  HS_EXPECT_NEAR(zN.re, projections::STEREO_INF, 1.0f);
  HS_EXPECT_NEAR(zN.im, 0.0f, 1.0f);

  // Inside the pole cap (denom < STEREO_POLE_EPS) the sentinel preserves the
  // (x,z) azimuth at magnitude STEREO_INF rather than collapsing onto +real.
  // At this scale the unit vector's y rounds to 1, so denom is exactly zero.
  math::Vector nearPole = math::Vector(6e-5f, 1.0f, 2.1e-5f).normalized();
  math::Complex zCap = projections::stereo(nearPole);
  HS_EXPECT_NEAR(std::sqrt(zCap.re * zCap.re + zCap.im * zCap.im),
                 projections::STEREO_INF, 1.0f);
  HS_EXPECT_NEAR(std::atan2(zCap.im, zCap.re),
                 std::atan2(nearPole.z, nearPole.x), 1e-3f);

  // The cap boundary is the crossover where the raw quotient would reach the
  // sentinel, not a step: outside it the projection stays below STEREO_INF.
  math::Complex zOut =
      projections::stereo(math::Vector(0.006f, 0.99998f, 0.0021f).normalized());
  HS_EXPECT_LT(std::sqrt(zOut.re * zOut.re + zOut.im * zOut.im),
               projections::STEREO_INF);

  // Large complex magnitude maps back to north pole.
  math::Vector pole =
      projections::inv_stereo(math::Complex(projections::STEREO_INF, 0));
  HS_EXPECT_VEC(pole, math::Vector(0, 1, 0), 1e-3f);

  // Plane origin maps to south pole.
  HS_EXPECT_VEC(projections::inv_stereo(math::Complex(0, 0)),
                math::Vector(0, -1, 0), 1e-3f);
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
  constexpr math::Complex VALUE(3.0f, -4.0f);
  static_assert(VALUE.squared_magnitude() == 25.0f);
  static_assert(VALUE.conjugate() == math::Complex(3.0f, 4.0f));
  static_assert(VALUE.conjugate().conjugate() == VALUE);
  static_assert(math::Complex(0.0f, -0.0f) == math::Complex(-0.0f, 0.0f));
  static_assert(math::Complex(1.0f, 0.0f) != math::Complex(1.00001f, 0.0f));
  HS_EXPECT_EQ(VALUE.magnitude(), 5.0f);
  HS_EXPECT_TRUE(VALUE * VALUE.conjugate() == math::Complex(25.0f, 0.0f));
  const math::Complex nonfinite(std::numeric_limits<float>::quiet_NaN(), 0);
  HS_EXPECT_FALSE(nonfinite == nonfinite.conjugate());
  for (float numerator : {0.0f, 1.0f, 1e15f}) {
    const math::Complex quotient =
        math::Complex(numerator, 0) / math::Complex(1e-23f, 0);
    HS_EXPECT_TRUE(std::isfinite(quotient.re));
    HS_EXPECT_EQ(quotient.im, 0.0f);
    HS_EXPECT_NEAR(quotient.re * 1e-23f / std::max(numerator, 1.0f),
                   numerator == 0.0f ? 0.0f : 1.0f, 2e-6f);
  }
  const math::Complex diagonal =
      math::Complex(1, 0) / math::Complex(1e-23f, 1e-23f);
  HS_EXPECT_NEAR(diagonal.re * 1e-23f, 0.5f, 2e-6f);
  HS_EXPECT_NEAR(diagonal.im * 1e-23f, -0.5f, 2e-6f);
  for (const float divisor : {3e-23f, 1e-21f}) {
    const math::Complex quotient =
        math::Complex(1, 0) / math::Complex(divisor, 0);
    HS_EXPECT_NEAR(quotient.re * divisor, 1.0f, 2e-6f);
    HS_EXPECT_EQ(quotient.im, 0.0f);
    HS_EXPECT_COMPLEX(math::project_div(math::Complex(divisor, divisor),
                                        math::Complex(divisor, 0)),
                      math::Complex(1, 1), 2e-6f);
  }
  math::Complex a(1, 2), b(3, 4);

  HS_EXPECT_COMPLEX(a + b, math::Complex(4, 6), 1e-6f);
  HS_EXPECT_COMPLEX(a - b, math::Complex(-2, -2), 1e-6f);
  // (1+2i)(3+4i) = 3 + 4i + 6i + 8i² = -5 + 10i
  HS_EXPECT_COMPLEX(a * b, math::Complex(-5, 10), 1e-5f);

  // operator/ is ordinary complex division: (a / b) * b ≈ a.
  math::Complex q = a / b;
  HS_EXPECT_COMPLEX(q * b, a, 1e-4f);

  HS_EXPECT_COMPLEX(a * math::Complex(1, 0), a, 1e-6f);

  // project_div matches ordinary division away from the singularity.
  HS_EXPECT_COMPLEX(math::project_div(a, b), a / b, 1e-6f);

  // project_div convention: 0 / 0 → 0.
  HS_EXPECT_COMPLEX(math::project_div(math::Complex(0, 0), math::Complex(0, 0)),
                    math::Complex(0, 0), 1e-6f);

  // project_div convention: nonzero / 0 → large magnitude in the numerator
  // direction.
  math::Complex inf_dir =
      math::project_div(math::Complex(1, 0), math::Complex(0, 0));
  HS_EXPECT_TRUE(std::abs(inf_dir.re) > 1e3f);

  // A divisor whose squared magnitude underflows is still a divisor, so a
  // tiny-but-equal homogeneous pair divides to 1 rather than to the sentinel.
  const math::Complex tiny(1e-30f, 0.0f);
  HS_EXPECT_COMPLEX(math::project_div(tiny, tiny), math::Complex(1, 0), 1e-6f);
  HS_EXPECT_COMPLEX(math::project_div(math::Complex(3e-30f, 4e-30f), tiny),
                    math::Complex(3, 4), 1e-5f);
  // A normal numerator over that same divisor is still the point at infinity.
  math::Complex underflow_inf = math::project_div(math::Complex(1, 0), tiny);
  HS_EXPECT_NEAR(underflow_inf.re, projections::STEREO_INF, 1e-1f);
  HS_EXPECT_NEAR(underflow_inf.im, 0.0f, 1e-6f);
}

// ============================================================================
// Mobius
// ============================================================================

/**
 * @brief Verifies MobiusParams constructors (8-float, 4-Complex, identity
 *        default) and that the a,b,c,d coefficients land in order.
 */
inline void test_mobius_params_accessors() {
  math::MobiusParams p(1, 2, 3, 4, 5, 6, 7, 8);
  HS_EXPECT_COMPLEX(p.a, math::Complex(1, 2), 0.0f);
  HS_EXPECT_COMPLEX(p.b, math::Complex(3, 4), 0.0f);
  HS_EXPECT_COMPLEX(p.c, math::Complex(5, 6), 0.0f);
  HS_EXPECT_COMPLEX(p.d, math::Complex(7, 8), 0.0f);

  math::MobiusParams q(math::Complex(1, 2), math::Complex(3, 4),
                       math::Complex(5, 6), math::Complex(7, 8));
  HS_EXPECT_COMPLEX(q.a, math::Complex(1, 2), 0.0f);
  HS_EXPECT_COMPLEX(q.d, math::Complex(7, 8), 0.0f);

  // Identity Mobius default: a=d=1, b=c=0.
  math::MobiusParams id;
  HS_EXPECT_COMPLEX(id.a, math::Complex(1, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.b, math::Complex(0, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.c, math::Complex(0, 0), 0.0f);
  HS_EXPECT_COMPLEX(id.d, math::Complex(1, 0), 0.0f);
}

/**
 * @brief Verifies mobius(z, params): identity, pure translation, pure scaling,
 *        and the inverting c != 0 branch including its pole.
 */
inline void test_mobius_transform() {
  math::Complex z(0.3f, 0.7f);

  math::MobiusParams id;
  HS_EXPECT_COMPLEX(math::mobius(z, id), z, 1e-4f);

  // Pure translation: (1·z + (2 - i)) / (0·z + 1) = z + (2 - i)
  math::MobiusParams trans(1, 0, 2.0f, -1.0f, 0, 0, 1, 0);
  HS_EXPECT_COMPLEX(math::mobius(z, trans),
                    math::Complex(z.re + 2.0f, z.im - 1.0f), 1e-4f);

  // Pure scaling: (3·z) / 1 = 3z
  math::MobiusParams scl(3, 0, 0, 0, 0, 0, 1, 0);
  HS_EXPECT_COMPLEX(math::mobius(z, scl),
                    math::Complex(z.re * 3.0f, z.im * 3.0f), 1e-4f);

  // Inversion: 1/z = conj(z) / |z|^2.
  math::MobiusParams inv(0, 0, 1, 0, 1, 0, 0, 0);
  const float r2 = z.re * z.re + z.im * z.im;
  HS_EXPECT_COMPLEX(math::mobius(z, inv), math::Complex(z.re / r2, -z.im / r2),
                    1e-4f);

  // General c != 0: (z + 1) / (z - 1).
  math::MobiusParams gen(1, 0, 1, 0, 1, 0, -1, 0);
  HS_EXPECT_COMPLEX(math::mobius(z, gen),
                    math::Complex(-0.42857143f, -1.42857143f), 1e-4f);

  // Its pole z = -d/c = 1 vanishes the denominator: project_div substitutes the
  // point at infinity along the numerator's direction, which inv_stereo reads
  // back as the north pole.
  math::Complex at_pole = math::mobius(math::Complex(1, 0), gen);
  HS_EXPECT_NEAR(at_pole.re, projections::STEREO_INF, 1.0f);
  HS_EXPECT_NEAR(at_pole.im, 0.0f, 1e-4f);
  HS_EXPECT_VEC(projections::inv_stereo(at_pole), math::Vector(0, 1, 0), 1e-6f);

  // A pole of the degenerate map (ad - bc == 0) is the indeterminate 0/0 form.
  math::MobiusParams degenerate(1, 0, -1, 0, 1, 0, -1, 0);
  HS_EXPECT_COMPLEX(math::mobius(math::Complex(1, 0), degenerate),
                    math::Complex(0, 0), 0.0f);
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
  math::Vector vUp = math::Vector(0.3f, 0.8f, 0.4f).normalized();
  math::Complex zUp = projections::gnomonic(vUp);
  math::Vector vUp_back = projections::inv_gnomonic(zUp, 1.0f);
  HS_EXPECT_VEC(vUp_back, vUp, 5e-3f);

  // Lower hemisphere: the hemisphere sign is passed to inv_gnomonic explicitly.
  math::Vector vDn = math::Vector(0.3f, -0.8f, 0.4f).normalized();
  math::Complex zDn = projections::gnomonic(vDn);
  math::Vector vDn_back = projections::inv_gnomonic(zDn, -1.0f);
  HS_EXPECT_VEC(vDn_back, vDn, 5e-3f);

  // North-pole pre-image: gnomonic(0,1,0) = (0, 0); inv → (0, 1, 0)
  HS_EXPECT_COMPLEX(projections::gnomonic(math::Vector(0, 1, 0)),
                    math::Complex(0, 0), 1e-4f);
  HS_EXPECT_VEC(projections::inv_gnomonic(math::Complex(0, 0), 1.0f),
                math::Vector(0, 1, 0), 1e-6f);

  // Saturated input → the equator point in the direction of z, sign-dependent
  // (gnomonic's singularity is the equator, not a pole).
  HS_EXPECT_VEC(projections::inv_gnomonic(
                    math::Complex(projections::STEREO_INF, 0), 1.0f),
                math::Vector(1, 0, 0), 1e-3f);
  HS_EXPECT_VEC(projections::inv_gnomonic(
                    math::Complex(projections::STEREO_INF, 0), -1.0f),
                math::Vector(-1, 0, 0), 1e-3f);
  HS_EXPECT_VEC(projections::inv_gnomonic(
                    math::Complex(0, -projections::STEREO_INF), 1.0f),
                math::Vector(0, 0, -1), 1e-3f);
  // A magnitude far past the sentinel must not square to infinity.
  HS_EXPECT_VEC(projections::inv_gnomonic(math::Complex(3e30f, 4e30f), 1.0f),
                math::Vector(0.6f, 0.0f, 0.8f), 1e-3f);
  // The sentinel test is radial: a diagonal point past the recognition radius
  // snaps back even though neither component reaches it alone.
  HS_EXPECT_VEC(projections::inv_gnomonic(math::Complex(4e3f, 4e3f), 1.0f),
                math::Vector(0.70710678f, 0.0f, 0.70710678f), 1e-3f);

  // Near-equator inputs get clamped to STEREO_INF
  math::Complex zEq = projections::gnomonic(math::Vector(1.0f, 1e-10f, 0.0f));
  HS_EXPECT_TRUE(std::abs(zEq.re) >= projections::STEREO_INF - 1.0f);

  // Round-trip identity through the singularity: |v.y| below ~2e-4 saturates
  // the projection, and the inverse must still land on the input.
  for (float y : {2e-4f, 1e-5f, 1e-9f, 0.0f, -1e-9f, -1e-5f, -2e-4f}) {
    for (float theta = 0.0f; theta < 2.0f * math::PI_F; theta += 0.37f) {
      const math::Vector v =
          math::Vector(cosf(theta), y, sinf(theta)).normalized();
      const math::Vector back = projections::inv_gnomonic(
          projections::gnomonic(v), copysignf(1.0f, y));
      HS_EXPECT_VEC(back, v, 1e-3f);
    }
  }

  // The floored divisor keys on the sign bit, so -0.0f projects like the tiny
  // negatives it is the limit of, not like +0.0f.
  math::Complex z_neg_zero =
      projections::gnomonic(math::Vector(1.0f, -0.0f, 0.0f));
  math::Complex z_tiny_neg =
      projections::gnomonic(math::Vector(1.0f, -1e-12f, 0.0f));
  math::Complex z_pos_zero =
      projections::gnomonic(math::Vector(1.0f, 0.0f, 0.0f));
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
  HS_EXPECT_NEAR(math::hash01(42u, 7u), 0.0169150233f, 1e-9f);
  HS_EXPECT_NEAR(math::hash01(0u, 0u), 0.030199945f, 1e-9f);
  HS_EXPECT_NEAR(math::hash01(4294967295u, 123u), 0.722669423f, 1e-9f);
  for (uint32_t i = 0; i < 256; ++i) {
    float h = math::hash01(i, 0u);
    HS_EXPECT_GE(h, 0.0f);
    HS_EXPECT_LT(h, 1.0f);
  }
  int differing = 0;
  for (uint32_t i = 0; i < 32; ++i)
    if (math::hash01(i, 1u) != math::hash01(i, 2u))
      differing++;
  HS_EXPECT_GT(differing, 24);

  // Seeds select streams, not permutations: the value sets must differ.
  constexpr uint32_t N = 256;
  const uint32_t pairs[3][2] = {{0u, 1u}, {1u, 2u}, {7u, 8u}};
  for (const auto &pair : pairs) {
    float a[N], b[N];
    for (uint32_t i = 0; i < N; ++i) {
      a[i] = math::hash01(i, pair[0]);
      b[i] = math::hash01(i, pair[1]);
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
  HS_EXPECT_NEAR(math::value_noise_1d(3.0f, 5u), math::hash01(3u, 5u), 1e-6f);
  HS_EXPECT_NEAR(math::value_noise_1d(-2.0f, 5u),
                 math::hash01(static_cast<uint32_t>(-2), 5u), 1e-6f);

  // Range over positive and negative coordinates.
  for (int i = -50; i <= 50; ++i) {
    float x = i * 0.173f;
    float n = math::value_noise_1d(x, 9u);
    HS_EXPECT_GE(n, 0.0f);
    HS_EXPECT_LT(n, 1.0f);
    float n2 = math::value_noise_2d(x, x * 0.7f, 9u);
    HS_EXPECT_GE(n2, 0.0f);
    HS_EXPECT_LT(n2, 1.0f);
  }

  // Continuity: adjacent samples differ by a bounded step, including at 0.
  for (int i = -400; i < 400; ++i) {
    float x = i * 0.005f;
    HS_EXPECT_LT(std::fabs(math::value_noise_1d(x + 0.005f, 3u) -
                           math::value_noise_1d(x, 3u)),
                 0.05f);
    HS_EXPECT_LT(std::fabs(math::value_noise_2d(x + 0.005f, 0.4f, 3u) -
                           math::value_noise_2d(x, 0.4f, 3u)),
                 0.05f);
    HS_EXPECT_LT(std::fabs(math::value_noise_2d(0.4f, x + 0.005f, 3u) -
                           math::value_noise_2d(0.4f, x, 3u)),
                 0.05f);
  }

  // Not constant, and seeds decorrelate.
  HS_EXPECT_TRUE(
      math::value_noise_1d(0.5f, 0u) != math::value_noise_1d(7.5f, 0u) ||
      math::value_noise_1d(2.5f, 0u) != math::value_noise_1d(9.5f, 0u));
  HS_EXPECT_TRUE(math::value_noise_2d(0.5f, 0.5f, 1u) !=
                 math::value_noise_2d(0.5f, 0.5f, 2u));
}

inline void test_twist_lens() {
  math::Vector input(0.6f, 0.5f, 0.6244998f);
  const math::Vector output = lenses::twist_lens(input);
  const float angle = lenses::TWIST_RATE * input.y;
  HS_EXPECT_NEAR(output.x, input.x * cosf(angle) - input.z * sinf(angle),
                 2e-3f);
  HS_EXPECT_NEAR(output.y, input.y, 1e-6f);
  HS_EXPECT_NEAR(output.z, input.x * sinf(angle) + input.z * cosf(angle),
                 2e-3f);
}

inline void test_kaleidoscope_lens() {
  const math::Vector input =
      math::Vector(-0.3f, 0.4f, -0.8660254f).normalized();
  const math::Vector output = lenses::kaleidoscope_lens(input);
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
    const float latitude = latitude_step * (0.5f * math::PI_F /
                                            static_cast<float>(LATITUDE_STEPS));
    const float radius = cosf(latitude);
    for (int longitude_step = 0; longitude_step < LONGITUDE_STEPS;
         ++longitude_step) {
      const float longitude =
          longitude_step *
          (math::TWO_PI_F / static_cast<float>(LONGITUDE_STEPS));
      visit(math::Vector(radius * cosf(longitude), sinf(latitude),
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
inline void expect_chamber_fold(const std::array<math::Vector, 3> &mirrors,
                                Fold fold) {
  constexpr float WALL_TOLERANCE = 1e-5f;
  for (const math::Vector &normal : mirrors)
    HS_EXPECT_NEAR(normal.magnitude(), 1.0f, 1e-6f);
  for_each_sphere_direction([&](const math::Vector &input) {
    const math::Vector folded = fold(input);
    for (const math::Vector &normal : mirrors)
      HS_EXPECT_GE(math::dot(folded, normal), -WALL_TOLERANCE);
    HS_EXPECT_NEAR(folded.magnitude(), input.magnitude(), 1e-4f);
    HS_EXPECT_VEC(fold(folded), folded, 1e-6f);
  });
}

/**
 * @brief Sweeps every reflection-group table wired into the lens catalog
 *        through the generic fold, plus the dodecahedral specialization.
 */
inline void test_polyhedral_kaleidoscope_chambers() {
  const std::array<math::Vector, 3> tables[] = {
      lenses::TETRAHEDRAL_MIRRORS,     lenses::OCTAHEDRAL_MIRRORS,
      lenses::DODECAHEDRAL_MIRRORS,    lenses::TRIANGULAR_PRISM_MIRRORS,
      lenses::SQUARE_PRISM_MIRRORS,    lenses::PENTAGONAL_PRISM_MIRRORS,
      lenses::HEXAGONAL_PRISM_MIRRORS, lenses::OCTAGONAL_PRISM_MIRRORS,
  };
  for (const std::array<math::Vector, 3> &mirrors : tables)
    expect_chamber_fold(mirrors, [&mirrors](const math::Vector &v) {
      return lenses::polyhedral_kaleidoscope_lens(v, mirrors);
    });
  expect_chamber_fold(lenses::DODECAHEDRAL_MIRRORS, [](const math::Vector &v) {
    return lenses::dodecahedral_kaleidoscope_lens(v);
  });
}

inline void test_dodecahedral_kaleidoscope_specialization() {
  for_each_sphere_direction([](const math::Vector &input) {
    const math::Vector generic = lenses::polyhedral_kaleidoscope_lens(
        input, lenses::DODECAHEDRAL_MIRRORS);
    const math::Vector specialized =
        lenses::dodecahedral_kaleidoscope_lens(input);
    HS_EXPECT_VEC(specialized, generic, 1e-5f);
  });
}

// ============================================================================
// 4D linear algebra
// ============================================================================

/** @brief The six coordinate planes a 4D orientation composes from. */
constexpr std::array<std::array<int, 2>, 6> PLANE_AXES = {
    {{{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}}};

/**
 * @brief Slack allowed on a single plane rotation's unit properties.
 * @details rotate_plane builds its rotor from cosf/sinf, so the rotor
 * identity holds to rounding; the measured bound is 8.4e-8.
 */
constexpr float PLANE_ROTATION_TOLERANCE = 1e-6f;

/**
 * @brief Slack allowed after all six planes have been composed.
 * @details Each factor contributes its own rounding, so the bound is the
 * single-rotation slack times the plane count; the measured bound is 4.0e-7.
 */
constexpr float COMPOSED_ROTATION_TOLERANCE = 4e-6f;

/**
 * @brief Squared length of a Vec4.
 * @param v Vector to measure.
 * @return The sum of squared components.
 */
inline float vec4_norm_squared(const math::Vec4 &v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
}

/**
 * @brief Dot product of two rows of a Mat4.
 * @param matrix Matrix to read.
 * @param i First row index.
 * @param j Second row index.
 * @return The dot product of rows i and j.
 */
inline float mat4_row_dot(const math::Mat4 &matrix, int i, int j) {
  float sum = 0.0f;
  for (int column = 0; column < math::VEC4_DIMENSIONS; ++column)
    sum += matrix.m[i][column] * matrix.m[j][column];
  return sum;
}

/**
 * @brief Draws a Vec4 with components in [-2, 2).
 * @param rng Generator supplying the components.
 * @return The sampled vector.
 */
inline math::Vec4 rand_vec4(hs::Pcg32 &rng) {
  math::Vec4 v{};
  for (int i = 0; i < math::VEC4_DIMENSIONS; ++i)
    v[i] = rand_uniform(rng, -2.0f, 2.0f);
  return v;
}

/**
 * @brief Pins Mat4::identity() and the row-major reading of Mat4::apply().
 */
inline void test_mat4_identity_and_apply() {
  const math::Mat4 id = math::Mat4::identity();
  for (int row = 0; row < math::VEC4_DIMENSIONS; ++row)
    for (int column = 0; column < math::VEC4_DIMENSIONS; ++column)
      HS_EXPECT_EQ(id.m[row][column], row == column ? 1.0f : 0.0f);

  const math::Vec4 v{{1.0f, -2.0f, 3.5f, 0.25f}};
  const math::Vec4 unchanged = id.apply(v);
  for (int i = 0; i < math::VEC4_DIMENSIONS; ++i)
    HS_EXPECT_EQ(unchanged[i], v[i]);

  // A transposed reading would still fix the identity, so apply() is also
  // scored against an asymmetric matrix.
  math::Mat4 ramp{};
  for (int row = 0; row < math::VEC4_DIMENSIONS; ++row)
    for (int column = 0; column < math::VEC4_DIMENSIONS; ++column)
      ramp.m[row][column] =
          static_cast<float>(row * math::VEC4_DIMENSIONS + column);
  const math::Vec4 mapped = ramp.apply(v);
  for (int row = 0; row < math::VEC4_DIMENSIONS; ++row) {
    float expected = 0.0f;
    for (int column = 0; column < math::VEC4_DIMENSIONS; ++column)
      expected += ramp.m[row][column] * v[column];
    HS_EXPECT_NEAR(mapped[row], expected, 1e-5f);
  }
}

/**
 * @brief Requires one plane rotation to be an isometry with orthonormal rows
 *        that leaves the two coordinates outside its plane untouched.
 * @param a First plane axis.
 * @param b Second plane axis.
 * @param angle Rotation angle in radians.
 * @param rng Generator supplying the sample points.
 */
inline void expect_plane_rotation_isometry(int a, int b, float angle,
                                           hs::Pcg32 &rng) {
  HS_CONTEXT(__func__, a, b);
  math::Mat4 rotation = math::Mat4::identity();
  math::rotate_plane(rotation, a, b, angle);

  for (int sample = 0; sample < 8; ++sample) {
    const math::Vec4 v = rand_vec4(rng);
    const math::Vec4 rotated = rotation.apply(v);
    HS_EXPECT_NEAR(vec4_norm_squared(rotated), vec4_norm_squared(v),
                   PLANE_ROTATION_TOLERANCE * vec4_norm_squared(v));
    // The rows outside the plane are still identity rows, so those two
    // coordinates survive bit for bit, not merely within tolerance.
    for (int i = 0; i < math::VEC4_DIMENSIONS; ++i)
      if (i != a && i != b)
        HS_EXPECT_EQ(rotated[i], v[i]);
  }

  for (int i = 0; i < math::VEC4_DIMENSIONS; ++i) {
    HS_EXPECT_NEAR(mat4_row_dot(rotation, i, i), 1.0f,
                   PLANE_ROTATION_TOLERANCE);
    for (int j = i + 1; j < math::VEC4_DIMENSIONS; ++j)
      HS_EXPECT_NEAR(mat4_row_dot(rotation, i, j), 0.0f,
                     PLANE_ROTATION_TOLERANCE);
  }
}

/**
 * @brief Sweeps every coordinate plane over a full turn of angles through the
 *        isometry properties.
 */
inline void test_rotate_plane_isometry() {
  hs::Pcg32 rng(20260825u);
  constexpr int ANGLE_STEPS = 12;
  for (const std::array<int, 2> &plane : PLANE_AXES)
    for (int step = 0; step < ANGLE_STEPS; ++step) {
      const float angle =
          -math::PI_F + (step * 2.0f * math::PI_F) / ANGLE_STEPS;
      expect_plane_rotation_isometry(plane[0], plane[1], angle, rng);
    }
}

/**
 * @brief Requires a quarter turn to carry the plane's first axis onto its
 *        second, so a rotation that degenerates to the identity cannot pass.
 */
inline void test_rotate_plane_quarter_turn() {
  for (const std::array<int, 2> &plane : PLANE_AXES) {
    HS_CONTEXT(__func__, plane[0], plane[1]);
    math::Mat4 rotation = math::Mat4::identity();
    math::rotate_plane(rotation, plane[0], plane[1], math::PI_F * 0.5f);
    math::Vec4 first_axis{};
    first_axis[plane[0]] = 1.0f;
    const math::Vec4 image = rotation.apply(first_axis);
    HS_EXPECT_NEAR(image[plane[0]], 0.0f, PLANE_ROTATION_TOLERANCE);
    HS_EXPECT_NEAR(image[plane[1]], 1.0f, PLANE_ROTATION_TOLERANCE);
  }
}

/**
 * @brief Requires +angle followed by -angle in the same plane to return the
 *        identity.
 */
inline void test_rotate_plane_inverse_composition() {
  for (const std::array<int, 2> &plane : PLANE_AXES)
    for (int step = 0; step < 8; ++step) {
      const float angle = 0.1f + step * 0.37f;
      HS_CONTEXT(__func__, plane[0], plane[1]);
      math::Mat4 rotation = math::Mat4::identity();
      math::rotate_plane(rotation, plane[0], plane[1], angle);
      math::rotate_plane(rotation, plane[0], plane[1], -angle);
      for (int row = 0; row < math::VEC4_DIMENSIONS; ++row)
        for (int column = 0; column < math::VEC4_DIMENSIONS; ++column)
          HS_EXPECT_NEAR(rotation.m[row][column], row == column ? 1.0f : 0.0f,
                         PLANE_ROTATION_TOLERANCE);
    }
}

inline void test_rotate_plane_composition_order() {
  math::Mat4 rotation = math::Mat4::identity();
  math::rotate_plane(rotation, 0, 1, math::PI_F * 0.5f);
  math::rotate_plane(rotation, 1, 2, math::PI_F * 0.5f);
  const math::Vec4 image = rotation.apply(math::Vec4{{1.0f, 0.0f, 0.0f, 0.0f}});
  for (int axis = 0; axis < math::VEC4_DIMENSIONS; ++axis)
    HS_EXPECT_NEAR(image[axis], axis == 2 ? 1.0f : 0.0f,
                   PLANE_ROTATION_TOLERANCE);
}

/**
 * @brief Requires a six-plane composition — the orientation shape HyperLattice
 *        builds — to stay an isometry with orthonormal rows.
 */
inline void test_rotate_plane_composition_stays_isometric() {
  math::Mat4 orientation = math::Mat4::identity();
  float phase = 0.37f;
  for (const std::array<int, 2> &plane : PLANE_AXES) {
    math::rotate_plane(orientation, plane[0], plane[1], phase);
    phase += 0.61f;
  }

  hs::Pcg32 rng(20260826u);
  for (int sample = 0; sample < 32; ++sample) {
    const math::Vec4 v = rand_vec4(rng);
    const math::Vec4 rotated = orientation.apply(v);
    HS_EXPECT_NEAR(vec4_norm_squared(rotated), vec4_norm_squared(v),
                   COMPOSED_ROTATION_TOLERANCE * vec4_norm_squared(v));
  }

  for (int i = 0; i < math::VEC4_DIMENSIONS; ++i) {
    HS_EXPECT_NEAR(mat4_row_dot(orientation, i, i), 1.0f,
                   COMPOSED_ROTATION_TOLERANCE);
    for (int j = i + 1; j < math::VEC4_DIMENSIONS; ++j)
      HS_EXPECT_NEAR(mat4_row_dot(orientation, i, j), 0.0f,
                     COMPOSED_ROTATION_TOLERANCE);
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

  test_mat4_identity_and_apply();
  test_rotate_plane_isometry();
  test_rotate_plane_quarter_turn();
  test_rotate_plane_inverse_composition();
  test_rotate_plane_composition_order();
  test_rotate_plane_composition_stays_isometric();

  return fixture.result();
}

} // namespace math3d_tests
} // namespace hs_test
