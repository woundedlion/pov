/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cmath>

#include "core/math/spherical_harmonics.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace spherical_harmonics_tests {

constexpr double PI_D = 3.14159265358979323846;

/**
 * @brief Factorial in double, the reference SHMath::factorial approximates.
 * @param n Non-negative argument.
 * @return n! as a double.
 */
inline double factorial_reference(int n) {
  double result = 1.0;
  for (int i = 2; i <= n; ++i)
    result *= i;
  return result;
}

/**
 * @brief Associated Legendre P_l^m(x) in double, carrying the Condon-Shortley
 *        phase and the (1 - x^2)^(m/2) factor SHMath divides out.
 * @param l Degree.
 * @param m Order in [0, l].
 * @param x Argument in [-1, 1].
 * @return P_l^m(x).
 */
inline double legendre_reference(int l, int m, double x) {
  double pmm = 1.0;
  if (m > 0) {
    const double somx2 = std::sqrt((1.0 - x) * (1.0 + x));
    double fact = 1.0;
    for (int i = 1; i <= m; ++i) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  double pmmp1 = x * (2.0 * m + 1.0) * pmm;
  if (l == m + 1)
    return pmmp1;
  double pll = 0.0;
  for (int ll = m + 2; ll <= l; ++ll) {
    pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

/**
 * @brief Real spherical harmonic in double, in the kernel's frame.
 * @param l Degree.
 * @param m Order in [-l, l].
 * @param p Unit direction; +y is the polar axis and the azimuth is measured in
 *        the x-z plane from +x.
 * @return The harmonic value.
 */
inline double harmonic_reference(int l, int m, const Vector &p) {
  const int abs_m = std::abs(m);
  double norm = std::sqrt(
      ((2.0 * l + 1.0) / (4.0 * PI_D)) *
      (factorial_reference(l - abs_m) / factorial_reference(l + abs_m)));
  if (m != 0)
    norm *= std::sqrt(2.0);
  const double theta =
      std::atan2(static_cast<double>(p.z), static_cast<double>(p.x));
  const double azimuth =
      m < 0 ? std::sin(abs_m * theta) : std::cos(abs_m * theta);
  return norm * legendre_reference(l, abs_m, static_cast<double>(p.y)) *
         azimuth;
}

/** @brief Unit direction at polar angle @p phi and azimuth @p theta. */
inline Vector spherical_direction(double phi, double theta) {
  const double radius = std::sin(phi);
  return Vector(static_cast<float>(radius * std::cos(theta)),
                static_cast<float>(std::cos(phi)),
                static_cast<float>(radius * std::sin(theta)));
}

/** @brief SHMath::factorial reproduces the product for every argument that
 *         stays inside a float. */
inline void test_factorial_matches_the_product() {
  HS_EXPECT_EQ(SHMath::factorial(0), 1.0f);
  HS_EXPECT_EQ(SHMath::factorial(1), 1.0f);
  HS_EXPECT_EQ(SHMath::factorial(-3), 1.0f);
  for (int n = 2; n <= SHMath::MAX_FACTORIAL_ARGUMENT; ++n) {
    const double reference = factorial_reference(n);
    HS_EXPECT_NEAR_REL(static_cast<double>(SHMath::factorial(n)), reference,
                       1e-5);
  }
}

/** @brief SHMath::legendre_seed is (-1)^m (2m - 1)!!, the P_m^m factor the
 *         reduced polynomial leaves out. */
inline void test_legendre_seed_is_the_signed_double_factorial() {
  double expected = 1.0;
  for (int m = 0; m <= 10; ++m) {
    HS_EXPECT_NEAR_REL(static_cast<double>(SHMath::legendre_seed(m)), expected,
                       1e-6);
    expected *= -(2.0 * m + 1.0);
  }
}

/** @brief reduced_legendre times its omitted seed and sin(phi)^m factor is the
 *         full associated Legendre polynomial. */
inline void test_reduced_legendre_restores_the_full_polynomial() {
  for (int l = 0; l <= 6; ++l)
    for (int m = 0; m <= l; ++m)
      for (int step = 0; step <= 40; ++step) {
        const double x = -1.0 + step * (2.0 / 40.0);
        const double sine = std::sqrt(std::max(0.0, 1.0 - x * x));
        const double restored = static_cast<double>(SHMath::reduced_legendre(
                                    l, m, static_cast<float>(x))) *
                                static_cast<double>(SHMath::legendre_seed(m)) *
                                std::pow(sine, static_cast<double>(m));
        HS_EXPECT_NEAR(restored, legendre_reference(l, m, x), 2e-4);
      }
}

/** @brief The reduced polynomial matches the hand-written closed forms through
 *         l = 3, independently of the reference recurrence. */
inline void test_reduced_legendre_matches_closed_forms() {
  for (int step = 0; step <= 20; ++step) {
    const double x = -1.0 + step * (2.0 / 20.0);
    const double s = std::sqrt(std::max(0.0, 1.0 - x * x));
    const double closed[4][4] = {
        {1.0, 0.0, 0.0, 0.0},
        {x, -s, 0.0, 0.0},
        {0.5 * (3.0 * x * x - 1.0), -3.0 * x * s, 3.0 * (1.0 - x * x), 0.0},
        {0.5 * x * (5.0 * x * x - 3.0), -1.5 * (5.0 * x * x - 1.0) * s,
         15.0 * x * (1.0 - x * x), -15.0 * s * s * s}};
    for (int l = 0; l <= 3; ++l)
      for (int m = 0; m <= l; ++m) {
        const double restored = static_cast<double>(SHMath::reduced_legendre(
                                    l, m, static_cast<float>(x))) *
                                static_cast<double>(SHMath::legendre_seed(m)) *
                                std::pow(s, static_cast<double>(m));
        HS_EXPECT_NEAR(restored, closed[l][m], 2e-5);
      }
  }
}

/** @brief harmonic_scale is the per-mode normalization times the omitted seed.
 */
inline void test_harmonic_scale_folds_in_the_seed() {
  HS_EXPECT_NEAR(SHMath::harmonic_scale(0, 0), 0.2820947918f, 1e-7f);
  HS_EXPECT_NEAR(SHMath::harmonic_scale(1, 0), 0.4886025119f, 1e-7f);
  HS_EXPECT_NEAR(SHMath::harmonic_scale(1, 1), -0.4886025119f, 1e-7f);
  HS_EXPECT_NEAR(SHMath::harmonic_scale(1, -1), -0.4886025119f, 1e-7f);
  HS_EXPECT_NEAR(SHMath::harmonic_scale(2, 2), 0.5462742153f, 1e-7f);
  HS_EXPECT_NEAR(SHMath::harmonic_scale(3, 3), -0.5900435899f, 1e-7f);
}

/** @brief The Cartesian harmonic kernel reproduces the angular reference over
 *         the sphere for every mode through l = 5. */
inline void test_spherical_harmonic_matches_the_angular_reference() {
  for (int l = 0; l <= 5; ++l)
    for (int m = -l; m <= l; ++m) {
      const float scale = SHMath::harmonic_scale(l, m);
      for (int phi_step = 1; phi_step < 12; ++phi_step)
        for (int theta_step = 0; theta_step < 16; ++theta_step) {
          const double phi = phi_step * (PI_D / 12.0);
          const double theta = theta_step * (2.0 * PI_D / 16.0) - PI_D;
          const Vector p = spherical_direction(phi, theta);
          HS_EXPECT_NEAR(
              static_cast<double>(SHMath::spherical_harmonic(l, m, p, scale)),
              harmonic_reference(l, m, p), 3e-4);
        }
    }
}

/** @brief The poles keep only the m = 0 modes, where the azimuthal factor of
 *         every other mode vanishes. */
inline void test_spherical_harmonic_poles_keep_only_the_zonal_modes() {
  for (float pole : {1.0f, -1.0f}) {
    const Vector p(0.0f, pole, 0.0f);
    for (int l = 0; l <= 5; ++l)
      for (int m = -l; m <= l; ++m) {
        const float value =
            SHMath::spherical_harmonic(l, m, p, SHMath::harmonic_scale(l, m));
        if (m == 0)
          HS_EXPECT_NEAR(static_cast<double>(value),
                         harmonic_reference(l, m, p), 1e-5);
        else
          HS_EXPECT_EQ(value, 0.0f);
      }
  }
}

/** @brief The modes are orthonormal over the sphere, which pins the
 *         normalization constants no per-sample comparison can. */
inline void test_spherical_harmonics_are_orthonormal() {
  constexpr int MAX_DEGREE = 3;
  constexpr int MODES = (MAX_DEGREE + 1) * (MAX_DEGREE + 1);
  constexpr int PHI_STEPS = 400;
  constexpr int THETA_STEPS = 64;
  std::array<float, MODES> scale{};
  for (int idx = 0; idx < MODES; ++idx) {
    const std::pair<int, int> lm = SHMath::decode_lm(idx);
    scale[idx] = SHMath::harmonic_scale(lm.first, lm.second);
  }
  std::array<std::array<double, MODES>, MODES> gram{};
  const double d_phi = PI_D / PHI_STEPS;
  const double d_theta = 2.0 * PI_D / THETA_STEPS;
  std::array<double, MODES> value{};
  for (int phi_step = 0; phi_step < PHI_STEPS; ++phi_step) {
    const double phi = (phi_step + 0.5) * d_phi;
    const double weight = std::sin(phi) * d_phi * d_theta;
    for (int theta_step = 0; theta_step < THETA_STEPS; ++theta_step) {
      const Vector p = spherical_direction(phi, theta_step * d_theta - PI_D);
      for (int idx = 0; idx < MODES; ++idx) {
        const std::pair<int, int> lm = SHMath::decode_lm(idx);
        value[idx] =
            SHMath::spherical_harmonic(lm.first, lm.second, p, scale[idx]);
      }
      for (int row = 0; row < MODES; ++row)
        for (int column = 0; column < MODES; ++column)
          gram[row][column] += weight * value[row] * value[column];
    }
  }
  for (int row = 0; row < MODES; ++row)
    for (int column = 0; column < MODES; ++column)
      HS_EXPECT_NEAR(gram[row][column], row == column ? 1.0 : 0.0, 2e-3);
}

/** @brief decode_lm inverts the flat index for every level it can reach. */
inline void test_decode_lm_inverts_the_flat_index() {
  for (int idx = 0; idx < 4096; ++idx) {
    const std::pair<int, int> lm = SHMath::decode_lm(idx);
    const int l = lm.first;
    const int m = lm.second;
    HS_EXPECT_GE(l, 0);
    HS_EXPECT_LE(std::abs(m), l);
    HS_EXPECT_EQ(l * l + l + m, idx);
  }
}

inline int run_spherical_harmonics_tests() {
  hs_test::ModuleFixture fixture("spherical_harmonics");
  test_factorial_matches_the_product();
  test_legendre_seed_is_the_signed_double_factorial();
  test_reduced_legendre_restores_the_full_polynomial();
  test_reduced_legendre_matches_closed_forms();
  test_harmonic_scale_folds_in_the_seed();
  test_spherical_harmonic_matches_the_angular_reference();
  test_spherical_harmonic_poles_keep_only_the_zonal_modes();
  test_spherical_harmonics_are_orthonormal();
  test_decode_lm_inverts_the_flat_index();
  return fixture.result();
}

} // namespace spherical_harmonics_tests
} // namespace hs_test
