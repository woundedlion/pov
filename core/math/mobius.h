/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file mobius.h
 * @brief Fractional-linear complex transforms and their sphere mappings.
 */

#include <cmath>
#include "math/stereographic.h"

namespace math {

/**
 * @brief mobius_transform's bound on squared |p|, linear |s|, and its
 * denominator floor near the homogeneous pole.
 */
inline constexpr float MOBIUS_POLE_EPS = 1e-12f;

/**
 * @brief Projection-domain complex division for the stereographic/Mobius maps.
 * @param num Numerator.
 * @param den Divisor.
 * @return num/den, except a quotient whose magnitude would reach STEREO_INF
 * clamps to the infinity sentinel along the numerator's direction, and an
 * exactly zero numerator returns (0,0).
 */
inline math::Complex project_div(const math::Complex &num,
                                 const math::Complex &den) {
  float den_re = den.re;
  float den_im = den.im;
  float num_re = num.re;
  float num_im = num.im;
  float denom = den_re * den_re + den_im * den_im;
  if (denom < std::numeric_limits<float>::min() &&
      (den_re != 0.0f || den_im != 0.0f)) {
    den_re *= COMPLEX_UNDERFLOW_LIFT;
    den_im *= COMPLEX_UNDERFLOW_LIFT;
    num_re *= COMPLEX_UNDERFLOW_LIFT;
    num_im *= COMPLEX_UNDERFLOW_LIFT;
    denom = den_re * den_re + den_im * den_im;
  }
  float num_mag = num_re * num_re + num_im * num_im;
  if (num_mag >= denom * (projections::STEREO_INF * projections::STEREO_INF)) {
    // Normalize by the larger component first: a numerator squared far above
    // the sentinel overflows to infinity and one far below it underflows to
    // zero, and either collapses the direction onto the origin.
    const float peak = fmaxf(std::abs(num.re), std::abs(num.im));
    if (peak == 0.0f)
      return math::Complex(0, 0);
    const float re = num.re / peak;
    const float im = num.im / peak;
    return projections::stereographic_detail::radial_scale(
        math::Complex(re, im), sqrtf(re * re + im * im),
        projections::STEREO_INF);
  }
  return math::Complex((num_re * den_re + num_im * den_im) / denom,
                       (num_im * den_re - num_re * den_im) / denom);
}

/** @brief Coefficients of a Mobius transform f(z) = (az + b) / (cz + d). */
struct MobiusParams {
  math::Complex a, b, c, d; /**< The four transform coefficients. */

  /**
   * @brief Default constructor producing the identity transform (a=d=1, b=c=0).
   */
  constexpr MobiusParams() : a(1, 0), b(0, 0), c(0, 0), d(1, 0) {}
  /**
   * @brief Constructs from four Complex coefficients.
   * @param coeff_a Coefficient a.
   * @param coeff_b Coefficient b.
   * @param coeff_c Coefficient c.
   * @param coeff_d Coefficient d.
   */
  constexpr MobiusParams(math::Complex coeff_a, math::Complex coeff_b,
                         math::Complex coeff_c, math::Complex coeff_d)
      : a(coeff_a), b(coeff_b), c(coeff_c), d(coeff_d) {}
  /**
   * @brief Constructs from eight floats (real/imaginary pairs per coefficient).
   * @param ar Real part of a.
   * @param ai Imaginary part of a.
   * @param br Real part of b.
   * @param bi Imaginary part of b.
   * @param cr Real part of c.
   * @param ci Imaginary part of c.
   * @param dr Real part of d.
   * @param di Imaginary part of d.
   */
  constexpr MobiusParams(float ar, float ai, float br, float bi, float cr,
                         float ci, float dr, float di)
      : a(ar, ai), b(br, bi), c(cr, ci), d(dr, di) {}
};

/**
 * @brief Mobius Transformation: f(z) = (az + b) / (cz + d).
 * @param z The complex input point.
 * @param params The four transform coefficients.
 * @return The transformed complex point.
 */
inline math::Complex mobius(const math::Complex &z,
                            const MobiusParams &params) {
  math::Complex num = (params.a * z) + params.b;
  math::Complex den = (params.c * z) + params.d;
  return project_div(num, den);
}

/**
 * @brief Applies a Mobius transformation to a vector.
 * @param v Unit vector to transform.
 * @param params Mobius transformation coefficients.
 * @return The transformed vector.
 * @details Fused stereographic projection, Mobius map and inverse projection
 * over the homogeneous plane coordinate (p : s); the pole is s = 0.
 */
inline math::Vector mobius_transform(const math::Vector &v,
                                     const MobiusParams &params) {
  float px = v.x, pz = v.z;
  float s = 1.0f - v.y;
  // The north pole and its epsilon neighborhood use the projective infinity
  // representative (1 : 0).
  if (px * px + pz * pz < MOBIUS_POLE_EPS && s < MOBIUS_POLE_EPS) {
    px = 1.0f;
    pz = 0.0f;
    s = 0.0f;
  }

  const float n_re = params.a.re * px - params.a.im * pz + params.b.re * s;
  const float n_im = params.a.re * pz + params.a.im * px + params.b.im * s;
  const float m_re = params.c.re * px - params.c.im * pz + params.d.re * s;
  const float m_im = params.c.re * pz + params.c.im * px + params.d.im * s;

  const float n2 = n_re * n_re + n_im * n_im;
  const float m2 = m_re * m_re + m_im * m_im;
  const float den = n2 + m2;
  // The absolute floor also catches nonsingular maps with small coefficients
  // near the pole; exact simultaneous zeros require a singular map.
  if (den < MOBIUS_POLE_EPS)
    return math::Vector(0.0f, 1.0f, 0.0f);

  const float inv = 1.0f / den;
  return math::Vector(2.0f * (n_re * m_re + n_im * m_im) * inv, (n2 - m2) * inv,
                      2.0f * (n_im * m_re - n_re * m_im) * inv);
}

/**
 * @brief Applies a gnomonic Mobius transformation to a vector.
 * @param v Unit vector to transform.
 * @param params Mobius transformation coefficients.
 * @return The transformed vector.
 * @details Projects to the gnomonic plane, applies the Mobius map, then
 * projects back to the hemisphere selected by the sign of v.y.
 */
inline math::Vector gnomonic_mobius_transform(const math::Vector &v,
                                              const MobiusParams &params) {
  math::Complex z = projections::gnomonic(v);
  math::Complex w = mobius(z, params);
  // copysignf keys on the sign bit, matching gnomonic's divisor floor for
  // v.y == -0.0f.
  return projections::inv_gnomonic(w, copysignf(1.0f, v.y));
}

} // namespace math
