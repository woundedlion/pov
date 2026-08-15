/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file spherical_harmonics.h
 * @brief Real spherical harmonics evaluated in Cartesian form on the unit
 *        sphere.
 */

#include <cmath>
#include <cstdlib>
#include <utility>

#include "math/3dmath.h"

namespace SHMath {
/** @brief Largest n whose reciprocal factorial is still normal in a float;
 *  1/34! is subnormal and 1/35! is zero. */
inline constexpr int MAX_FACTORIAL_ARGUMENT = 33;

/**
 * @brief Factorial of n as a float.
 * @param n Non-negative integer whose factorial is computed; kept small.
 * @return n! as a float.
 * @details Float precision; large n would overflow float precision.
 */
inline float factorial(int n) {
  if (n <= 1)
    return 1.0f;
  float result = 1.0f;
  for (int i = 2; i <= n; i++)
    result *= i;
  return result;
}

/**
 * @brief Seed P_m^m of the associated Legendre recurrence, with its
 *        sin(phi)^m factor divided out.
 * @param m Order (m >= 0).
 * @return (-1)^m (2m - 1)!!, independent of the polynomial argument.
 */
inline float legendre_seed(int m) {
  float pmm = 1.0f;
  float fact = 1.0f;
  for (int i = 1; i <= m; i++) {
    pmm *= -fact;
    fact += 2.0f;
  }
  return pmm;
}

/**
 * @brief Associated Legendre polynomial with its sin(phi)^m factor and its
 *        P_m^m seed divided out.
 * @param l Degree (l >= 0).
 * @param m Order (0 <= m <= l).
 * @param x Argument, equal to cos(phi) with |x| <= 1.
 * @return P_l^m(x) / ((1 - x²)^(m/2) * P_m^m), a polynomial in x.
 * @details Standard upward recurrence in l, started from a unit seed. The
 * recurrence is linear and homogeneous in its two seeds, so the omitted P_m^m
 * scales the whole result and callers fold it into their per-mode constant
 * (harmonic_scale) instead of paying for it per sample. The sin(phi)^m factor
 * is restored by the caller's azimuthal term, which carries it in Cartesian
 * form.
 */
inline float reduced_legendre(int l, int m, float x) {
  float pmm = 1.0f;
  if (l == m)
    return pmm;

  float pmmp1 = x * (2.0f * m + 1.0f);
  if (l == m + 1)
    return pmmp1;

  float pll = 0;
  for (int ll = m + 2; ll <= l; ll++) {
    pll = ((2.0f * ll - 1.0f) * x * pmmp1 - (ll + m - 1.0f) * pmm) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

/**
 * @brief Precompute the spherical-harmonic normalization factor for (l, m).
 * @param l Degree (l >= 0).
 * @param m Order in [-l, l].
 * @return Normalization factor N, constant per shape.
 * @details Traps when l + |m| exceeds 33: the factorial ratio is subnormal at
 * 34 (1/34! = 3.4e-39, a few mantissa bits, zero under FTZ) and (l + |m|)!
 * overflows float from 35 on, so the ratio silently collapses.
 */
inline float normalization(int l, int m) {
  int abs_m = std::abs(m);
  HS_CHECK(l + abs_m <= MAX_FACTORIAL_ARGUMENT,
           "spherical harmonic normalization: l + |m| = %d collapses the float "
           "factorial ratio",
           l + abs_m);
  float N = sqrtf(((2.0f * l + 1.0f) / (4.0f * PI_F)) *
                  (factorial(l - abs_m) / factorial(l + abs_m)));
  return (m != 0) ? sqrtf(2.0f) * N : N;
}

/**
 * @brief Per-mode constant factor spherical_harmonic() expects as its N.
 * @param l Degree (l >= 0).
 * @param m Order in [-l, l].
 * @return normalization(l, m) times the Legendre seed P_|m|^|m|.
 */
inline float harmonic_scale(int l, int m) {
  return normalization(l, m) * legendre_seed(std::abs(m));
}

/**
 * @brief Evaluate a real spherical harmonic with a precomputed norm.
 * @param l Degree (l >= 0).
 * @param m Order in [-l, l]; sign selects cos/sin azimuthal factor.
 * @param p Unit direction in the shape's local frame; y is cos(phi).
 * @param N Precomputed per-mode factor from harmonic_scale(l, m).
 * @return The harmonic value.
 * @details sin(phi)^|m| * cos(|m| theta) is Re((x + iz)^|m|), and the sine
 * counterpart is Im, so the (1 - y²)^(|m|/2) factor of P_l^m is exactly what
 * the Cartesian power supplies. Evaluating the pair together leaves the whole
 * harmonic polynomial in (x, y, z): no angle, sine, square root, or division.
 */
inline float spherical_harmonic(int l, int m, const Vector &p, float N) {
  int abs_m = std::abs(m);
  float re = 1.0f, im = 0.0f;
  for (int i = 0; i < abs_m; i++) {
    float next_re = re * p.x - im * p.z;
    im = re * p.z + im * p.x;
    re = next_re;
  }
  float Q = reduced_legendre(l, abs_m, hs::clamp(p.y, -1.0f, 1.0f));
  return N * Q * ((m < 0) ? im : re);
}

/**
 * @brief Decode a flat harmonic index into its (l, m) pair.
 * @param idx Flat index where idx = l*l + l + m.
 * @return Pair {l, m} with l = floor(sqrt(idx)) and m in [-l, l].
 * @details Seeds l from a float sqrt, then snaps it to the exact integer floor:
 * sqrtf at (or just below) a perfect square can round to l-epsilon and truncate
 * to l-1, which would push m to +l — outside the level's valid [-l, l] band and
 * into the next level's order. The correction loops make l provably exact
 * (l*l <= idx < (l+1)*(l+1)), so the returned order is always valid. Cold path
 * (a few calls per frame).
 */
inline std::pair<int, int> decode_lm(int idx) {
  int l = static_cast<int>(sqrtf(static_cast<float>(idx)));
  while ((l + 1) * (l + 1) <= idx)
    ++l;
  while (l > 0 && l * l > idx)
    --l;
  return {l, idx - l * l - l};
}
} // namespace SHMath
