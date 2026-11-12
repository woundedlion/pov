/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file ring_knots.h
 * @brief Ring knot positions and periodic strided interpolation.
 */

#include "core/math/3dmath.h"
#include "core/platform/attributes.h"

namespace SDF {

/** @brief A ring's colatitude and knot-step trig, passed in VFP registers. */
struct KnotStep {
  float cos_t; /**< Cosine of the ring colatitude. */
  float sin_t; /**< Sine of the ring colatitude. */
  float cos_d; /**< Cosine of one knot cell's azimuth step. */
  float sin_d; /**< Sine of one knot cell's azimuth step. */
};

/**
 * @brief Steps a ring's knot positions around its azimuth by the
 * angle-addition recurrence, starting at knot 0.
 */
struct KnotWalk {
  const math::Basis &basis; /**< Ring frame; basis.v is the stack axis. */
  float cos_t;              /**< Cosine of the ring colatitude. */
  float sin_t;              /**< Sine of the ring colatitude. */
  float cos_d;              /**< Cosine of one knot cell's azimuth step. */
  float sin_d;              /**< Sine of one knot cell's azimuth step. */
  float cos_a = 1.0f;       /**< Cosine of the current knot's azimuth. */
  float sin_a = 0.0f;       /**< Sine of the current knot's azimuth. */

  /** @brief Unit position of the current knot. */
  HS_HOT_INLINE math::Vector position() const {
    return ((basis.v * cos_t) + ((basis.u * cos_a) + (basis.w * sin_a)) * sin_t)
        .normalized();
  }

  /** @brief Moves to the next knot. */
  HS_HOT_INLINE void advance() {
    const float next_cos = cos_a * cos_d - sin_a * sin_d;
    sin_a = sin_a * cos_d + cos_a * sin_d;
    cos_a = next_cos;
  }
};

/**
 * @brief Interpolates a periodic grid using uniform Catmull-Rom segments.
 * @tparam Value Callable returning a grid sample by index.
 * @param value Grid sample accessor.
 * @param d Grid stride; a positive divisor of n with 2*d <= n.
 * @param k Knot index in [0, n).
 * @param n Periodic knot count.
 * @return The interpolated value at k.
 */
template <typename Value>
HS_HOT_INLINE float catmull_rom_strided(Value value, int d, int k, int n) {
  const int r = k % d;
  if (r == 0)
    return value(k);
  auto wrap = [n](int index) HS_HOT_INLINE {
    return index < 0 ? index + n : index >= n ? index - n : index;
  };
  const int b = k - r;
  const float t = static_cast<float>(r) / d;
  const float p0 = value(wrap(b - d)), p1 = value(b), p2 = value(wrap(b + d)),
              p3 = value(wrap(b + 2 * d));
  return p1 + 0.5f * t *
                  (p2 - p0 +
                   t * (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3 +
                        t * (3.0f * (p1 - p2) + p3 - p0)));
}

} // namespace SDF
