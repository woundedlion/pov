/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file stereographic.h
 * @brief The stereographic, gnomonic and Mobius maps between the sphere and
 *        the complex plane, with the shading helpers for stereographic-space
 *        patterns: pole attenuation, pattern normalization and trig-argument
 *        clamping.
 */

#include "math/3dmath.h"

/**
 * @brief Conventional representation of the point at infinity on the complex
 * plane.
 * @details Single source of truth for the pole sentinel: every forward
 * projection that hits a singularity emits this magnitude, and every inverse
 * projection recognizes it (see the two thresholds below).
 */
inline constexpr float STEREO_INF = 1e4f;

/**
 * @brief |z| at/above which inv_stereo() treats its input as the infinity
 * sentinel.
 * @details Half of STEREO_INF: an intervening Mobius map can scale the emitted
 * sentinel toward (not past) zero, so the inverse needs margin below the emitted
 * magnitude to still snap it back to the pole. (Squared to avoid a sqrt.)
 */
inline constexpr float STEREO_INF_RECOGNIZE = STEREO_INF * 0.5f;

/**
 * @brief Soft-limit for a stereographic coordinate fed (times a pattern
 * frequency) into fast_sinf/fast_cosf, beyond which range reduction bands.
 * @details Unclamped, stereo()'s up-to-1e4 magnitudes times a pattern frequency
 * drive trig arguments toward ~2e5, where a float's ULP (~0.02 rad) makes
 * fast_sinf's range reduction band near the pole. Clamping to ±this bound holds
 * the reduction error to ~5e-4 rad in the pole cap; non-pole coordinates
 * (|z| ~ O(10)) never reach it.
 */
inline constexpr float STEREO_PATTERN_ARG_LIMIT = 4096.0f;

/**
 * @brief 1 - v.y below which stereo() is inside the north-pole cap and emits the
 * sentinel magnitude instead of the raw quotient.
 * @details Placed at the algebraic crossover: on the unit sphere
 * |stereo(v)| = sqrt((1 + v.y) / (1 - v.y)) reaches STEREO_INF exactly here.
 * Float spacing puts that crossover out of reach, though — 1 - v.y is exact for
 * v.y in [0.5, 1] and its smallest nonzero value (2^-24) already quotients to
 * ~5.8e3, so the cap is entered only at 1 - v.y == 0 and the sentinel steps up
 * ~1.7x from the largest magnitude the quotient can produce. Both magnitudes
 * clear STEREO_INF_RECOGNIZE, so inv_stereo returns either to the pole. Any
 * retune below 2^-24 is inert.
 */
inline constexpr float STEREO_POLE_EPS = 2.0f / (STEREO_INF * STEREO_INF);

/**
 * @brief Squared magnitude below which mobius_transform treats a homogeneous
 * pair (p : s) as the degenerate exact-pole form and substitutes the point at
 * infinity. Tighter than math::EPS_LEN_SQ so only a pair essentially exactly
 * zero loses its direction.
 */
inline constexpr float STEREO_DIV_NUM_EPS_SQ = 1e-12f;

/**
 * @brief Radius (a length, not squared) in the (x,z) plane below which the
 * north-pole azimuth is treated as undefined, so stereo() falls back to the
 * +real axis instead of scaling a near-zero direction.
 */
inline constexpr float STEREO_AZIMUTH_EPS = 1e-12f;

/**
 * @brief |v.y| (a length) at or below which gnomonic() floors the projection
 * divisor to avoid div-by-zero at the equator, clamping the result to the
 * sentinel. Unrelated to the coincidentally-equal math::EPS_NORMAL_SQ.
 */
inline constexpr float STEREO_EQUATOR_EPS = 1e-9f;

/**
 * @brief Projection-domain complex division for the stereographic/Mobius maps.
 * @param num Numerator.
 * @param den Divisor.
 * @return num/den, except a quotient whose magnitude would reach STEREO_INF
 * collapses to the infinity sentinel scaled along the numerator's direction (so
 * a near-singular divisor still yields a finite point). Only an exactly zero
 * numerator is the indeterminate 0/0 form, which returns (0,0); a nonzero
 * numerator keeps its direction however small it is.
 * @details NOT general complex division (see Complex::operator/): the STEREO_INF
 * clamp and the 0/0 -> 0 case are the point-at-infinity conventions the sphere
 * projections depend on.
 */
inline Complex project_div(const Complex &num, const Complex &den) {
  float denom = den.re * den.re + den.im * den.im;
  float num_mag = num.re * num.re + num.im * num.im;
  if (num_mag >= denom * (STEREO_INF * STEREO_INF)) {
    // Normalize by the larger component first: a numerator squared far above
    // the sentinel overflows to infinity and one far below it underflows to
    // zero, and either collapses the direction onto the origin.
    const float peak = std::max(std::abs(num.re), std::abs(num.im));
    if (peak == 0.0f)
      return Complex(0, 0);
    const float re = num.re / peak;
    const float im = num.im / peak;
    const float scale = STEREO_INF / sqrtf(re * re + im * im);
    return Complex(re * scale, im * scale);
  }
  return Complex((num.re * den.re + num.im * den.im) / denom,
                 (num.im * den.re - num.re * den.im) / denom);
}

/**
 * @brief Coefficients of a Mobius transform f(z) = (az + b) / (cz + d).
 * @details Stores the four coefficients as first-class `Complex` values;
 * animators mutate the `.re`/`.im` components in place. The eight-float
 * constructor is retained for terse literal initialization.
 */
struct MobiusParams {
  Complex a, b, c, d; /**< The four transform coefficients. */

  /** @brief No live-config hook for TransformerPool::prepare_frame(). */
  static constexpr bool NEEDS_REFRESH_FROM = false;
  /** @brief No derived state, so no prepare_frame() sync() hook. */
  static constexpr bool NEEDS_SYNC = false;

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
  constexpr MobiusParams(Complex coeff_a, Complex coeff_b, Complex coeff_c,
                         Complex coeff_d)
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
 * @brief Stereographic Projection: Sphere -> Complex Plane.
 * @param v Point on the unit sphere.
 * @return The projected complex-plane coordinate.
 * @details Inside the north-pole cap (v.y ≈ 1) the result is the infinity
 * sentinel magnitude carrying the (x,z) azimuth; only the exact pole, where the
 * azimuth is undefined, lands on the real axis.
 */
inline Complex stereo(const Vector &v) {
  float denom = 1.0f - v.y;
  if (denom < STEREO_POLE_EPS) {
    // North-pole cap: emit the sentinel but keep the (x,z) azimuth. At the exact
    // pole (x = z = 0) the azimuth is undefined → +real fallback.
    float r = sqrtf(v.x * v.x + v.z * v.z);
    if (r < STEREO_AZIMUTH_EPS)
      return Complex(STEREO_INF, 0.0f);
    float scale = STEREO_INF / r;
    return Complex(v.x * scale, v.z * scale);
  }
  return Complex(v.x / denom, v.z / denom);
}

/**
 * @brief Inverse Stereographic Projection: Complex Plane -> Sphere.
 * @param z Complex-plane coordinate (the infinity sentinel maps to the pole).
 * @return The corresponding point on the unit sphere.
 */
inline Vector inv_stereo(const Complex &z) {
  // |z| >= STEREO_INF_RECOGNIZE → North Pole (catches the sentinel and any point
  // within ~0.02° of the pole; squared compare avoids a sqrt).
  float r2 = z.re * z.re + z.im * z.im;
  if (r2 >= STEREO_INF_RECOGNIZE * STEREO_INF_RECOGNIZE)
    return Vector(0.0f, 1.0f, 0.0f);
  return Vector(2 * z.re / (r2 + 1), (r2 - 1) / (r2 + 1), 2 * z.im / (r2 + 1));
}

/**
 * @brief Mobius Transformation: f(z) = (az + b) / (cz + d).
 * @param z The complex input point.
 * @param params The four transform coefficients.
 * @return The transformed complex point.
 */
inline Complex mobius(const Complex &z, const MobiusParams &params) {
  Complex num = (params.a * z) + params.b;
  Complex den = (params.c * z) + params.d;
  return project_div(num, den);
}

/**
 * @brief Gnomonic Projection: Sphere -> Plane (Equator at Infinity).
 * @param v Point on the unit sphere.
 * @return The projected plane coordinate (equator points clamp to the sentinel).
 * @details Projects from center (0,0,0) to the plane y=1 (tangent at the North
 * Pole (0,1,0), i.e. j=1).
 * @note Identifies antipodes: `v` and `-v` project to the same plane
 * coordinate, so the hemisphere is lost. A caller that round-trips through
 * inv_gnomonic must track the sign of `v.y` and pass it back via
 * inv_gnomonic's `original_sign`.
 */
inline Complex gnomonic(const Vector &v) {
  // Floor the divisor to ±STEREO_EQUATOR_EPS to avoid div-by-zero at v.y == 0,
  // then clamp the magnitude to STEREO_INF. A near-equator point clamps to the
  // sentinel, which inv_gnomonic snaps back to the equator.
  // copysignf, not a >= 0 test: -0.0f must floor negative like the values it
  // is the limit of.
  float div = (std::abs(v.y) < STEREO_EQUATOR_EPS)
                  ? copysignf(STEREO_EQUATOR_EPS, v.y)
                  : v.y;
  float gx = v.x / div;
  float gz = v.z / div;
  // Radial clamp, matching project_div: clamping the components separately
  // would drag a saturated point towards the nearest diagonal, discarding the
  // azimuth the inverse reads back.
  const float magnitude_sq = gx * gx + gz * gz;
  if (magnitude_sq > STEREO_INF * STEREO_INF) {
    const float scale = STEREO_INF / sqrtf(magnitude_sq);
    gx *= scale;
    gz *= scale;
  }
  return Complex(gx, gz);
}

/**
 * @brief Inverse Gnomonic: Plane -> Sphere.
 * @param z Complex point on the plane.
 * @param original_sign Sign of the y-component (j) of the original vector, used
 * to restore the hemisphere the forward projection collapsed.
 * @return The corresponding point on the unit sphere; the plane's infinity is
 * the equator point in the direction of z, not a pole (the projection ray
 * flattens into y = 0 as |z| grows).
 */
inline Vector inv_gnomonic(const Complex &z, float original_sign) {
  // Clamped-to-infinity → equator, recognized from STEREO_INF_RECOGNIZE (margin
  // snaps a Mobius-shrunk sentinel back to the limit). Radial, matching the
  // forward clamp: a per-component test would make the snap-back radius
  // azimuth-dependent. Squared compare avoids a sqrt, and a magnitude past the
  // float range overflows to infinity, which still clears the bound.
  if (z.re * z.re + z.im * z.im >=
      STEREO_INF_RECOGNIZE * STEREO_INF_RECOGNIZE) {
    // Normalize by the larger component first: squaring a magnitude well past
    // the sentinel would overflow to infinity and yield a zero vector.
    const float scale = 1.0f / std::max(std::abs(z.re), std::abs(z.im));
    const float re = z.re * scale;
    const float im = z.im * scale;
    const float inv_len = original_sign / sqrtf(re * re + im * im);
    return Vector(re * inv_len, 0.0f, im * inv_len);
  }
  // Project (re, 1, im) back onto unit sphere
  float len = sqrtf(z.re * z.re + z.im * z.im + 1.0f);
  float inv_len = 1.0f / len;

  // Restore hemisphere sign (Upper or Lower)
  return Vector(z.re * inv_len * original_sign, // i
                inv_len * original_sign,        // j
                z.im * inv_len * original_sign  // k
  );
}

/**
 * @brief Smooth pole attenuation for stereographic-space effects.
 * @param r_sq Pre-computed |z|² (z.re² + z.im²).
 * @param pole_fade Attenuation radius (larger = wider fade zone).
 * @return Falloff factor 1/(1 + r²/pole_fade²) in (0, 1].
 * @details Stereographic projection sends the far pole to infinity, so |z|²
 * grows without bound near it; this falloff is 1 at the projection origin and
 * decays toward 0 with distance, taming that singularity.
 */
__attribute__((always_inline)) inline float pole_attenuation(float r_sq,
                                                             float pole_fade) {
  // Floor the radius so a 0 pole_fade can't divide by zero and poison the warp.
  const float pf = pole_fade > 1e-3f ? pole_fade : 1e-3f;
  return 1.0f / (1.0f + (r_sq / (pf * pf)));
}

/**
 * @brief Pole-attenuates a stereographic pattern value and maps it to [0, 1].
 * @param pattern Raw pattern value in [-1, 1].
 * @param r_sq Pre-computed |z|² driving the pole fade.
 * @param pole_fade Attenuation radius (larger = wider fade zone).
 * @return Pole-attenuated value normalized to [0, 1].
 */
inline float pole_normalize_pattern(float pattern, float r_sq,
                                    float pole_fade) {
  return (pattern * pole_attenuation(r_sq, pole_fade) + 1.0f) * 0.5f;
}

/**
 * @brief Scales a warped stereographic coordinate into bounded trig arguments.
 * @param w Warped stereographic coordinate.
 * @param pattern_freq Spatial frequency multiplier.
 * @return Frequency-scaled components clamped to ±STEREO_PATTERN_ARG_LIMIT.
 * @details Near the pole |w| -> STEREO_INF, so w*pattern_freq can reach ~2e5
 * where fast_sinf range reduction bands; the clamp keeps both components inside
 * the accurate range.
 */
inline Complex stereo_pattern_args(const Complex &w, float pattern_freq) {
  return Complex(hs::clamp(w.re * pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                           STEREO_PATTERN_ARG_LIMIT),
                 hs::clamp(w.im * pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                           STEREO_PATTERN_ARG_LIMIT));
}
