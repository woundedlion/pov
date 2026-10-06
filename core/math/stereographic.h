/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file stereographic.h
 * @brief Stereographic and gnomonic forward/inverse projection kernels.
 */

#include <cmath>
#include "math/3dmath.h"

namespace projections {

/**
 * @brief Conventional representation of the point at infinity on the complex
 * plane.
 * @details Forward projections emit this magnitude at a singularity; inverses
 * recognize it via STEREO_INF_RECOGNIZE.
 */
inline constexpr float STEREO_INF = 1e4f;

/**
 * @brief |z| at/above which inv_stereo() treats its input as the infinity
 * sentinel.
 * @details Half of STEREO_INF, so a sentinel shrunk by a Mobius map still
 * snaps back to the pole.
 */
inline constexpr float STEREO_INF_RECOGNIZE = STEREO_INF * 0.5f;

/**
 * @brief 1 - v.y below which stereo() is inside the north-pole cap and emits the
 * sentinel magnitude instead of the raw quotient.
 * @details |stereo(v)| reaches STEREO_INF exactly here; in float, the cap is
 * entered only at 1 - v.y == 0.
 */
inline constexpr float STEREO_POLE_EPS = 2.0f / (STEREO_INF * STEREO_INF);

/**
 * @brief Radius (a length, not squared) in the (x,z) plane below which the
 * north-pole azimuth is treated as undefined, so stereo() falls back to the
 * +real axis instead of scaling a near-zero direction.
 */
inline constexpr float STEREO_AZIMUTH_EPS = 1e-12f;

/**
 * @brief |v.y| (a length) at or below which gnomonic() floors the projection
 * divisor to avoid div-by-zero at the equator, clamping the result to the
 * sentinel.
 */
inline constexpr float STEREO_EQUATOR_EPS = 1e-9f;

namespace stereographic_detail {

/** @brief Scales a nonzero planar direction from its supplied length. */
inline math::Complex radial_scale(const math::Complex &direction, float length,
                                  float radius) {
  const float scale = radius / length;
  return math::Complex(direction.re * scale, direction.im * scale);
}

} // namespace stereographic_detail

/**
 * @brief Stereographic Projection: Sphere -> Complex Plane.
 * @param v Point on the unit sphere.
 * @return The projected complex-plane coordinate.
 * @details Inside the north-pole cap (v.y ≈ 1) the result is the infinity
 * sentinel magnitude carrying the (x,z) azimuth; a planar radius below
 * STEREO_AZIMUTH_EPS uses the +real axis.
 */
inline math::Complex stereo(const math::Vector &v) {
  float denom = 1.0f - v.y;
  if (denom < STEREO_POLE_EPS) {
    float r = sqrtf(v.x * v.x + v.z * v.z);
    if (r < STEREO_AZIMUTH_EPS)
      return math::Complex(STEREO_INF, 0.0f);
    return stereographic_detail::radial_scale(math::Complex(v.x, v.z), r,
                                              STEREO_INF);
  }
  return math::Complex(v.x / denom, v.z / denom);
}

/**
 * @brief Inverse Stereographic Projection: Complex Plane -> Sphere.
 * @param z Complex-plane coordinate (the infinity sentinel maps to the pole).
 * @return The corresponding point on the unit sphere.
 */
inline math::Vector inv_stereo(const math::Complex &z) {
  // Sentinel, or within ~0.02° of the pole → north pole.
  float r2 = z.squared_magnitude();
  if (r2 >= STEREO_INF_RECOGNIZE * STEREO_INF_RECOGNIZE)
    return math::Vector(0.0f, 1.0f, 0.0f);
  return math::Vector(2 * z.re / (r2 + 1), (r2 - 1) / (r2 + 1),
                      2 * z.im / (r2 + 1));
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
 * inv_gnomonic's `hemisphere_sign`.
 */
inline math::Complex gnomonic(const math::Vector &v) {
  // copysignf: -0.0f floors negative like the values it is the limit of.
  float div = (std::abs(v.y) < STEREO_EQUATOR_EPS)
                  ? copysignf(STEREO_EQUATOR_EPS, v.y)
                  : v.y;
  float gx = v.x / div;
  float gz = v.z / div;
  // Radial clamp preserves the azimuth inv_gnomonic reads back.
  const float magnitude_sq = gx * gx + gz * gz;
  if (magnitude_sq > STEREO_INF * STEREO_INF) {
    return stereographic_detail::radial_scale(math::Complex(gx, gz),
                                              sqrtf(magnitude_sq), STEREO_INF);
  }
  return math::Complex(gx, gz);
}

/**
 * @brief Inverse Gnomonic: Plane -> Sphere.
 * @param z Complex point on the plane.
 * @param hemisphere_sign +1 or -1, the sign of the y-component (j) of the
 * original vector, restoring the hemisphere the forward projection collapsed;
 * any other magnitude scales the result off the unit sphere.
 * @return The corresponding point on the unit sphere; the plane's infinity is
 * the equator point in the direction of hemisphere_sign * z, not a pole (the projection ray
 * flattens into y = 0 as |z| grows).
 */
inline math::Vector inv_gnomonic(const math::Complex &z,
                                 float hemisphere_sign) {
  // Sentinel → equator; a magnitude past the float range overflows to infinity,
  // which still clears the bound.
  if (z.squared_magnitude() >= STEREO_INF_RECOGNIZE * STEREO_INF_RECOGNIZE) {
    // Normalize by the larger component first: squaring a magnitude well past
    // the sentinel would overflow to infinity and yield a zero vector.
    const float scale = 1.0f / fmaxf(std::abs(z.re), std::abs(z.im));
    const float re = z.re * scale;
    const float im = z.im * scale;
    const math::Complex equator = stereographic_detail::radial_scale(
        math::Complex(re, im), sqrtf(re * re + im * im), hemisphere_sign);
    return math::Vector(equator.re, 0.0f, equator.im);
  }
  // Project (re, 1, im) back onto unit sphere
  float len = sqrtf(z.squared_magnitude() + 1.0f);
  float inv_len = 1.0f / len;

  // Restore hemisphere sign (Upper or Lower)
  return math::Vector(z.re * inv_len * hemisphere_sign, // i
                      inv_len * hemisphere_sign,        // j
                      z.im * inv_len * hemisphere_sign  // k
  );
}

} // namespace projections
