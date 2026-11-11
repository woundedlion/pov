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
#include <array>
#include <cstdlib>
#include <utility>

#include "core/math/3dmath.h"

namespace SHMath {
/** @brief Largest degree represented by the fixed polynomial table. */
inline constexpr int MAX_POLYNOMIAL_DEGREE = 8;

/** @brief Offset of one parity-compressed Legendre polynomial.
 * @param l Degree.
 * @param m Nonnegative order.
 * @return First coefficient index. */
constexpr int polynomial_offset(int l, int m) {
  int offset = 0;
  for (int degree = 0; degree <= l; ++degree)
    for (int order = 0; order <= degree; ++order) {
      if (degree == l && order == m)
        return offset;
      offset += (degree - order) / 2 + 1;
    }
  return offset;
}

/** @brief Constant-time offsets indexed by l*(l+1)/2+m. */
inline constexpr auto LEGENDRE_OFFSETS HS_PROGMEM_UNIQUE(legendre_offsets) =
    [] {
      std::array<uint8_t, 45> result{};
      for (int l = 0; l <= MAX_POLYNOMIAL_DEGREE; ++l)
        for (int m = 0; m <= l; ++m)
          result[l * (l + 1) / 2 + m] =
              static_cast<uint8_t>(polynomial_offset(l, m));
      return result;
    }();

/** @brief Builds reduced Legendre coefficients by upward recurrence.
 * @return Parity-compressed coefficients in ascending powers. */
constexpr auto make_legendre_coefficients() {
  std::array<float, 95> result{};
  for (int m = 0; m <= MAX_POLYNOMIAL_DEGREE; ++m) {
    std::array<double, 9> previous{}, current{};
    current[0] = 1.0;
    for (int l = m; l <= MAX_POLYNOMIAL_DEGREE; ++l) {
      const int parity = (l - m) % 2;
      const int offset = polynomial_offset(l, m);
      for (int i = parity; i <= l - m; i += 2)
        result[offset + i / 2] = static_cast<float>(current[i]);
      std::array<double, 9> next{};
      for (int i = 0; i < 9; ++i) {
        if (i > 0)
          next[i] = (2.0 * l + 1.0) * current[i - 1] / (l + 1 - m);
        next[i] -= (l + m) * previous[i] / (l + 1 - m);
      }
      previous = current;
      current = next;
    }
  }
  return result;
}

/** @brief Reduced Legendre polynomial coefficients through degree eight. */
inline constexpr auto LEGENDRE_COEFFICIENTS
    HS_PROGMEM_UNIQUE(legendre_coefficients) = make_legendre_coefficients();
static_assert(polynomial_offset(8, 8) + 1 == LEGENDRE_COEFFICIENTS.size());
/**
 * @brief Factorial of n as a float.
 * @param n Non-negative integer whose factorial is computed; kept small.
 * @return n! as a float.
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
 * @details Fixed polynomials through degree eight; larger degrees use recurrence.
 */
inline float reduced_legendre(int l, int m, float x) {
  if (l >= 0 && l <= MAX_POLYNOMIAL_DEGREE && m >= 0 && m <= l) {
    const int degree = l - m;
    const int offset = LEGENDRE_OFFSETS[l * (l + 1) / 2 + m];
    float value = LEGENDRE_COEFFICIENTS[offset + degree / 2];
    const float square = x * x;
    for (int i = degree / 2 - 1; i >= 0; --i)
      value = value * square + LEGENDRE_COEFFICIENTS[offset + i];
    return degree % 2 ? x * value : value;
  }
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

/** @brief Derivative of a reduced Legendre polynomial through degree eight.
 * @param l Degree in [0, 8].
 * @param m Order in [0, l].
 * @param x Polynomial argument.
 * @return Derivative with respect to x. */
HS_FLASH_INLINE inline float reduced_legendre_derivative(int l, int m,
                                                         float x) {
  HS_CHECK(l >= 0 && l <= MAX_POLYNOMIAL_DEGREE && m >= 0 && m <= l,
           "Legendre derivative: unsupported mode");
  const int degree = l - m;
  if (degree == 0)
    return 0.0f;
  const int offset = LEGENDRE_OFFSETS[l * (l + 1) / 2 + m];
  const int parity = degree % 2;
  float value = degree * LEGENDRE_COEFFICIENTS[offset + degree / 2];
  for (int i = degree / 2 - 1; i >= (parity == 0 ? 1 : 0); --i)
    value =
        value * x * x + (2 * i + parity) * LEGENDRE_COEFFICIENTS[offset + i];
  return parity ? value : x * value;
}

/**
 * @brief Precompute the spherical-harmonic normalization factor for (l, m).
 * @param l Degree (l >= 0).
 * @param m Order in [-l, l].
 * @return Normalization factor N, constant per shape.
 * @details Traps on |m| > l. Ratios smaller than the float range may
 * underflow to zero.
 */
inline float normalization(int l, int m) {
  int abs_m = std::abs(m);
  HS_CHECK(l >= 0 && abs_m <= l,
           "spherical harmonic normalization: order %d is outside [-%d, %d]", m,
           l, l);
  float ratio = 1.0f;
  for (int64_t i = static_cast<int64_t>(l) + abs_m; i > l - abs_m; --i)
    ratio /= static_cast<float>(i);
  float n = sqrtf(((2.0f * l + 1.0f) / (4.0f * math::PI_F)) * ratio);
  return (m != 0) ? sqrtf(2.0f) * n : n;
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
 * the Cartesian power supplies.
 */
inline float spherical_harmonic(int l, int m, const math::Vector &p, float N) {
  int abs_m = std::abs(m);
  float re = 1.0f, im = 0.0f;
  for (int i = 0; i < abs_m; i++) {
    float next_re = re * p.x - im * p.z;
    im = re * p.z + im * p.x;
    re = next_re;
  }
  float q = reduced_legendre(l, abs_m, hs::clamp(p.y, -1.0f, 1.0f));
  return N * q * ((m < 0) ? im : re);
}

/** @brief Evaluates the gradient for an admitted harmonic mode.
 * @param l Degree in [0, MAX_POLYNOMIAL_DEGREE].
 * @param m Nonnegative order in [0, l].
 * @param p Cartesian point.
 * @param N Prepared harmonic scale.
 * @param phase Complex phase multiplying the azimuthal power.
 * @return Ambient gradient of the polynomial extension. */
HS_HOT_FLASH_MEMBER inline math::Vector
spherical_harmonic_gradient_unchecked(int l, int m, const math::Vector &p,
                                      float N, math::Complex phase) {
  HS_AUDIT_CHECK(l >= 0 && l <= MAX_POLYNOMIAL_DEGREE && m >= 0 && m <= l);
  math::Complex power(1.0f, 0.0f), previous;
  const math::Complex w(p.x, p.z);
  for (int i = 0; i < m; ++i) {
    previous = power;
    power = power * w;
  }
  const math::Complex derivative =
      phase * previous * math::Complex(float(m), 0.0f);
  const int degree = l - m;
  const int offset = LEGENDRE_OFFSETS[l * (l + 1) / 2 + m];
  const float x = p.y;
  const float square = x * x;
  float series = LEGENDRE_COEFFICIENTS[offset + degree / 2];
  float series_slope = 0.0f;
  for (int i = degree / 2 - 1; i >= 0; --i) {
    series_slope = series_slope * square + series;
    series = series * square + LEGENDRE_COEFFICIENTS[offset + i];
  }
  const bool odd = degree % 2;
  const float q = odd ? x * series : series;
  const float dq =
      odd ? series + 2.0f * square * series_slope : 2.0f * x * series_slope;
  return {N * q * derivative.re, N * dq * (phase * power).re,
          -N * q * derivative.im};
}

/** @brief Cartesian gradient of a phased real harmonic.
 * @param l Degree in [0, 8].
 * @param m Nonnegative order in [0, l].
 * @param p Cartesian point.
 * @param N Prepared harmonic scale.
 * @param phase Complex phase multiplying the azimuthal power.
 * @return Ambient gradient of the polynomial extension. */
inline math::Vector
spherical_harmonic_gradient(int l, int m, const math::Vector &p, float N,
                            math::Complex phase = {1.0f, 0.0f}) {
  HS_CHECK(l >= 0 && l <= MAX_POLYNOMIAL_DEGREE && m >= 0 && m <= l,
           "harmonic gradient: unsupported mode");
  return spherical_harmonic_gradient_unchecked(l, m, p, N, phase);
}

/**
 * @brief Decode a flat harmonic index into its (l, m) pair.
 * @param idx Flat index where idx = l*l + l + m; must be non-negative.
 * @return Pair {l, m} with l = floor(sqrt(idx)) and m in [-l, l].
 * @details Corrects the rounded float-square-root seed to the exact integer
 * floor (l*l <= idx < (l+1)*(l+1)), keeping the returned order in [-l, l].
 */
inline std::pair<int, int> decode_lm(int idx) {
  HS_CHECK(idx >= 0, "decode_lm: flat index %d is negative", idx);
  int l = static_cast<int>(sqrtf(static_cast<float>(idx)));
  while (static_cast<int64_t>(l + 1) * (l + 1) <= idx)
    ++l;
  while (l > 0 && static_cast<int64_t>(l) * l > idx)
    --l;
  return {l, idx - l * l - l};
}
} // namespace SHMath
