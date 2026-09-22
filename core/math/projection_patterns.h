/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file projection_patterns.h
 * @brief Pole attenuation and bounded projection-space pattern inputs.
 */

#include "math/stereographic.h"

namespace projections {

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
 * @brief Smooth pole attenuation for stereographic-space effects.
 * @param r_sq Pre-computed |z|² (z.re² + z.im²).
 * @param singularity_fade Attenuation radius (larger = wider fade zone).
 * @return Falloff factor 1/(1 + r²/singularity_fade²) in (0, 1].
 * @details Stereographic projection sends the far pole to infinity, so |z|²
 * grows without bound near it; this falloff is 1 at the projection origin and
 * decays toward 0 with distance, taming that singularity.
 */
__attribute__((always_inline)) inline float
pole_attenuation(float r_sq, float singularity_fade) {
  // Floor the radius so a 0 singularity_fade can't divide by zero and poison the warp.
  const float pf = singularity_fade > 1e-3f ? singularity_fade : 1e-3f;
  return 1.0f / (1.0f + (r_sq / (pf * pf)));
}

/**
 * @brief Pole-attenuates a stereographic pattern value and maps it to [0, 1].
 * @param pattern Raw pattern value in [-1, 1].
 * @param r_sq Pre-computed |z|² driving the pole fade.
 * @param singularity_fade Attenuation radius (larger = wider fade zone).
 * @return Pole-attenuated value normalized to [0, 1].
 */
inline float pole_normalize_pattern(float pattern, float r_sq,
                                    float singularity_fade) {
  return (pattern * pole_attenuation(r_sq, singularity_fade) + 1.0f) * 0.5f;
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

} // namespace projections
