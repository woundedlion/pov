/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file stereographic.h
 * @brief Shading helpers for stereographic-space patterns: pole attenuation,
 *        pattern normalization, trig-argument clamping and the glitch lens.
 */

#include "math/3dmath.h"

/**
 * @brief Smooth pole attenuation for stereographic-space effects.
 * @param r_sq Pre-computed |z|² (z.re² + z.im²).
 * @param pole_fade Attenuation radius (larger = wider fade zone).
 * @return Falloff factor 1/(1 + r²/pole_fade²) in (0, 1].
 * @details Stereographic projection sends the far pole to infinity, so |z|²
 * grows without bound near it; this falloff is 1 at the projection origin and
 * decays toward 0 with distance, taming that singularity.
 */
inline float pole_attenuation(float r_sq, float pole_fade) {
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

/**
 * @brief Applies a trig-free glitch lens to a sphere direction.
 * @param v Unit direction vector on the sphere.
 * @return Direction after latitude doubling and azimuth tripling; returns
 * the up vector near the lens axis.
 */
inline Vector glitch_lens(const Vector &v) {
  const float x2 = v.x * v.x;
  const float z2 = v.z * v.z;
  const float radius2 = x2 + z2;
  constexpr float MIN_AXIS_RADIUS2 = 1e-6f;
  if (radius2 < MIN_AXIS_RADIUS2)
    return Vector(0.0f, 1.0f, 0.0f);

  const float inverse_radius2 = 1.0f / radius2;
  const float double_y = 2.0f * v.y;
  return Vector(double_y * v.x * (4.0f * x2 * inverse_radius2 - 3.0f),
                2.0f * v.y * v.y - 1.0f,
                double_y * v.z * (3.0f - 4.0f * z2 * inverse_radius2));
}
