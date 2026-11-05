/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
/** @file ripple_motion.h
 * @brief Frame-spanning ripple bands and vertex capture kernels. */
#include "core/animation/animation.h"
#include "core/math/ripple.h"
#include "core/memory.h"
#include <type_traits>

namespace RippleMotion {
/** @brief A band about a ripple's front at one end of its step. */
struct RippleBand {
  math::Vector center; ///< Ripple center.
  float cos_near;      ///< Cosine of the band's inner edge.
  float cos_far;       ///< Cosine of the band's outer edge.
  float step;          ///< The ripple's step arc.

  /**
     * @brief Whether a model vertex stands in the band.
     * @param v The vertex.
     * @return True between the band's edges.
     */
  HS_HOT_INLINE bool holds(const math::Vector &v) const {
    const float c = math::dot(v, center);
    return c <= cos_near && c >= cos_far;
  }
};

/**
   * @brief The bands about the active ripples that step this frame, one per
   *        end of each step, in pairs.
   * @tparam Pool Active ripple pool.
   * @param pool Pool exposing active_count() and active_params().
   * @param leg_arc The model's largest motion (radians).
   * @param arena Arena receiving the bands, for the caller's scope.
   * @param count Receives the number of bands.
   * @return The bands.
   */
template <typename Pool>
HS_FLASH_MEMBER RippleBand *ripple_bands(const Pool &pool, float leg_arc,
                                         Arena &arena, int &count) {
  RippleBand *bands = arena.allocate_n<RippleBand>(2 * pool.active_count() + 1);
  count = 0;
  for (int k = 0; k < pool.active_count(); ++k) {
    const Animation::RippleParams &p = pool.active_params(k);
    if (p.step_arc() == 0.0f)
      continue;
    for (const float phase : {p.previous_phase, p.phase}) {
      const float reach = 2.0f * p.half_width() + leg_arc;
      bands[count++] = {
          p.center, cosf(hs::clamp(phase - reach, 0.0f, math::PI_F)),
          cosf(hs::clamp(phase + reach, 0.0f, math::PI_F)), p.step_arc()};
    }
  }
  return bands;
}

/** @brief One ripple's frame constants and mutable capture sample.
   * @tparam Counters Optional test capture counters, or void. */
template <typename Counters = void> struct MotionRipple {
  Animation::RippleParams params; ///< Sample parameters and previous endpoint.
  float phase_delta = 0.0f;       ///< Phase change over the frame.
  float amplitude_delta = 0.0f;   ///< Amplitude change over the frame.

  /** @brief Captures the frame's initial wave and interpolation constants.
   * @param p Active wave with its retained endpoint. */
  HS_FLASH_INLINE void prepare(const Animation::RippleParams &p) {
    params = p.during_step(0.0f);
    phase_delta = p.phase - p.previous_phase;
    amplitude_delta = p.amplitude - p.previous_amplitude;
  }

  /** @brief Includes the final endpoint in the capture rejection band.
   * @param p Final wave endpoint. */
  HS_HOT_INLINE void widen(const Animation::RippleParams &p) {
    params.cos_threshold_min =
        fmaxf(params.cos_threshold_min, p.cos_threshold_min);
    params.cos_threshold_max =
        fminf(params.cos_threshold_max, p.cos_threshold_max);
  }

  /** @brief Applies the wave at a frame fraction after its band test.
     * @param v Vertex after preceding ripples.
     * @param u Capture fraction.
     * @return Displaced vertex. */
  HS_HOT_INLINE math::Vector transform(const math::Vector &v, float u) {
#if HS_ENABLE_TEST_HOOKS
    if constexpr (!std::is_void_v<Counters>)
      Counters::candidate();
#endif
    const float cosine = math::dot(v, params.center);
    if (cosine > params.cos_threshold_min || cosine < params.cos_threshold_max)
      return v;
#if HS_ENABLE_TEST_HOOKS
    if constexpr (!std::is_void_v<Counters>)
      Counters::interpolated();
#endif
    params.phase = params.previous_phase + phase_delta * u;
    params.amplitude = params.previous_amplitude + amplitude_delta * u;
    if (params.amplitude <= RippleKernel::NEGLIGIBLE_AMPLITUDE)
      return v;
    return RippleKernel::transform_at_distance(
        v, params, math::fast_acos(hs::clamp(cosine, -1.0f, 1.0f)),
        params.phase);
  }
};
} // namespace RippleMotion
