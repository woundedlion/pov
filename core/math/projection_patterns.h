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
 * @details Holds the range-reduction error to ~5e-4 rad in the pole cap.
 */
inline constexpr float STEREO_PATTERN_ARG_LIMIT = 4096.0f;

/**
 * @brief Smooth pole attenuation for stereographic-space effects.
 * @param r_sq Pre-computed |z|² (z.re² + z.im²).
 * @param singularity_fade Attenuation radius (larger = wider fade zone).
 * @return Falloff factor 1/(1 + r²/pf²), with pf = max(singularity_fade, 1e-3).
 * @details Stereographic projection sends the far pole to infinity, so |z|²
 * grows without bound near it; this falloff is 1 at the projection origin and
 * decays toward 0 with distance, taming that singularity.
 */
__attribute__((always_inline)) inline float
pole_attenuation(float r_sq, float singularity_fade) {
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
 */
inline math::Complex stereo_pattern_args(const math::Complex &w,
                                         float pattern_freq) {
  return math::Complex(hs::clamp(w.re * pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                                 STEREO_PATTERN_ARG_LIMIT),
                       hs::clamp(w.im * pattern_freq, -STEREO_PATTERN_ARG_LIMIT,
                                 STEREO_PATTERN_ARG_LIMIT));
}

} // namespace projections
