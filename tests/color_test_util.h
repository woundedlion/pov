/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Hue-difference wrap for OKLCH hue assertions, independent of
 * wrap_angle_pi so a hue assertion never compares the module under test
 * against itself.
 */
#pragma once

#include <cmath>

#include "core/math/3dmath.h"

namespace hs_test {

/**
 * @brief Per-channel budget of an identity pass through OKLab.
 * @details linear RGB -> fast_cbrt LMS -> OKLab and back is not exact; with the
 * rotation angle or chroma scale at identity, or an oklch_to_pixel round trip,
 * that conversion plus float -> uint16 rounding is the whole residue.
 */
inline constexpr float OKLAB_ROUND_TRIP_TOL = 16.0f;

/**
 * @brief Wraps a hue difference into [-PI, PI] for circular comparison.
 * @param dh Raw hue difference in radians.
 * @return The difference in [-PI, PI], or non-finite input unchanged.
 */
inline float wrap_hue_delta(float dh) {
  if (!std::isfinite(dh))
    return dh;
  return std::remainder(dh, 2.0f * math::PI_F);
}

} // namespace hs_test
