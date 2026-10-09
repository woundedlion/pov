/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Geodesic arc-distance oracle, tolerant-equality predicates, and assertions for Vector,
 * Quaternion and Complex.
 */
#pragma once

#include <algorithm>
#include <cmath>

#include "core/math/3dmath.h"
#include "tests/test_harness.h"

namespace hs_test {

/** @brief Snorm3 half-step quantization plus float encode/decode rounding. */
inline constexpr float SNORM3_COMPONENT_BOUND =
    0.5f / math::Snorm3::SCALE + 1e-7f;

/**
 * @brief Tests whether two vectors agree componentwise within a tolerance.
 * @param a First vector operand.
 * @param b Second vector operand.
 * @param tol Per-component absolute tolerance.
 * @return True if every component of a and b agrees within tol.
 */
inline bool approx_vec(const math::Vector &a, const math::Vector &b,
                       float tol) {
  return approx(a.x, b.x, tol) && approx(a.y, b.y, tol) &&
         approx(a.z, b.z, tol);
}

inline float max_component_delta(const math::Vector &a, const math::Vector &b) {
  return hs_test::fold_worst(
      std::abs(a.x - b.x),
      hs_test::fold_worst(std::abs(a.y - b.y), std::abs(a.z - b.z)));
}

/**
 * @brief Angle between two near-parallel vectors, below angle_between's floor.
 * @param a First vector (non-zero).
 * @param b Second vector (non-zero).
 * @return The angle in radians.
 * @details Differences the normalized endpoints in double.
 */
inline double small_angle_between(const math::Vector &a,
                                  const math::Vector &b) {
  const double ax = a.x, ay = a.y, az = a.z;
  const double bx = b.x, by = b.y, bz = b.z;
  const double na = std::sqrt(ax * ax + ay * ay + az * az);
  const double nb = std::sqrt(bx * bx + by * by + bz * bz);
  const double dx = ax / na - bx / nb;
  const double dy = ay / na - by / nb;
  const double dz = az / na - bz / nb;
  const double chord = std::sqrt(dx * dx + dy * dy + dz * dz);
  return 2.0 * std::asin(std::min(1.0, chord / 2.0));
}

/** @brief Returns angular distance from a unit direction to a geodesic arc. */
inline float arc_angular_distance(const math::Vector &p, const math::Vector &a,
                                  const math::Vector &b) {
  const auto angle = [](const math::Vector &u, const math::Vector &v) {
    return acosf(hs::clamp(math::dot(u, v), -1.0f, 1.0f));
  };
  const float endpoint_distance = std::min(angle(p, a), angle(p, b));
  math::Vector normal = math::cross(a, b);
  const float normal_length = normal.length();
  if (normal_length < 1e-6f)
    return endpoint_distance;
  normal = normal / normal_length;
  const float offset = math::dot(p, normal);
  math::Vector foot = p - normal * offset;
  const float foot_length = foot.length();
  if (foot_length > 1e-6f) {
    foot = foot / foot_length;
    const float span = angle(a, b);
    if (angle(a, foot) + angle(foot, b) <= span + 1e-4f)
      return asinf(hs::clamp(std::fabs(offset), 0.0f, 1.0f));
  }
  return endpoint_distance;
}
/**
 * @brief Tests whether two quaternions agree within a tolerance.
 * @param a First quaternion operand.
 * @param b Second quaternion operand.
 * @param tol Per-component absolute tolerance.
 * @return True if the scalar and vector parts of a and b agree within tol.
 */
inline bool approx_quat(const math::Quaternion &a, const math::Quaternion &b,
                        float tol) {
  return approx(a.r, b.r, tol) && approx_vec(a.v, b.v, tol);
}
/**
 * @brief Tests whether two complex numbers agree within a tolerance.
 * @param a First complex operand.
 * @param b Second complex operand.
 * @param tol Per-component absolute tolerance.
 * @return True if the real and imaginary parts of a and b agree within tol.
 */
inline bool approx_complex(const math::Complex &a, const math::Complex &b,
                           float tol) {
  return approx(a.re, b.re, tol) && approx(a.im, b.im, tol);
}

} // namespace hs_test

// Captures each operand once so loop-driven assertions don't re-evaluate side
// effects, then compares and (on failure) prints both values.
#define HS_EXPECT_APPROX(pred, a, b, tol)                                      \
  do {                                                                         \
    const auto hs_lhs = (a);                                                   \
    const auto hs_rhs = (b);                                                   \
    hs_test::report_cmp(hs_test::pred(hs_lhs, hs_rhs, (tol)), hs_lhs, hs_rhs,  \
                        #a " ~= " #b " (tol=" #tol ")", __func__, __FILE__,    \
                        __LINE__);                                             \
  } while (0)

/**
 * @brief Tolerant equality assertions for vectors, quaternions, and complex
 *        values; the failure message stringizes the two compared expressions
 *        and prints both operands componentwise.
 */
#define HS_EXPECT_VEC(a, b, tol) HS_EXPECT_APPROX(approx_vec, a, b, tol)
#define HS_EXPECT_QUAT(a, b, tol) HS_EXPECT_APPROX(approx_quat, a, b, tol)
#define HS_EXPECT_COMPLEX(a, b, tol) HS_EXPECT_APPROX(approx_complex, a, b, tol)
