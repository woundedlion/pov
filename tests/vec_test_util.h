/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Shared tolerant-equality predicates and assertion macros for Vector,
 * Quaternion and Complex — the 3dmath, geometry, sdf, spatial, effects and
 * reaction_graph suites all assert on these types.
 */
#pragma once

#include "core/math/3dmath.h"
#include "tests/test_harness.h"

namespace hs_test {

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

} // namespace hs_test

/**
 * @brief Tolerant equality assertions for vectors, quaternions, and complex
 *        values; the failure message stringizes the two compared expressions.
 */
#define HS_EXPECT_VEC(a, b, tol)                                               \
  HS_EXPECT(hs_test::approx_vec((a), (b), (tol)),                              \
            #a " ~= " #b " (tol=" #tol ")")
#define HS_EXPECT_QUAT(a, b, tol)                                              \
  HS_EXPECT(hs_test::approx_quat((a), (b), (tol)),                             \
            #a " ~= " #b " (tol=" #tol ")")
#define HS_EXPECT_COMPLEX(a, b, tol)                                           \
  HS_EXPECT(hs_test::approx_complex((a), (b), (tol)),                          \
            #a " ~= " #b " (tol=" #tol ")")
