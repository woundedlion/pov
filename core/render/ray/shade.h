/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file shade.h
 * @brief Depth appearance, layer compositing and subray filtering. */

#include "color/baked_palette.h"
#include "color/layer_composite.h"
#include "render/ray/events.h"
#include "render/ray/march.h"

namespace Raycast {

/**
 * @brief Depth appearance for nonnegative ray distances.
 * @pre inv_far is finite and nonnegative, and t * inv_far <= 1 for every
 * shaded distance t.
 */
struct Appearance {
  float inv_far = .1f;     ///< Reciprocal of the fog distance.
  float near_start = 0;    ///< Distance where the near fade-in starts.
  float near_inv_span = 1; ///< Reciprocal of the near fade-in length.
  const BakedPalette *palette = nullptr; ///< Depth palette; required to shade.
  float gain = 1; /**< Brightness scale applied to every layer's color. */

  /**
   * @brief Whether an interval keeps depth palette lookups in [0, 1].
   * @param interval Ray-parameter range to shade.
   * @return True when the interval is valid and within the fog distance.
   */
  bool valid_for(Interval interval) const {
    return interval.valid() && finite(inv_far) && inv_far >= 0.0f &&
           interval.far * inv_far <= 1.0f;
  }

  /**
   * @brief Squared fog falloff times the near fade-in at distance t.
   * @param t Ray distance.
   * @return Opacity in [0, 1].
   */
  __attribute__((always_inline)) float opacity(float t) const {
    const float fog = fmaxf(0.0f, 1.0f - t * inv_far);
    return fog * fog * math::cubic_kernel((t - near_start) * near_inv_span);
  }
  /**
   * @brief Depth-graded palette color at distance t.
   * @param t Ray distance.
   * @return Palette color scaled by gain and nearness.
   */
  __attribute__((always_inline)) Pixel color(float t) const {
    const float nearness = 1 - t * inv_far;
    return palette->get_color_unit(nearness) *
           (gain * (.45f + .55f * nearness));
  }
  /**
   * @brief Adds the depth-graded color at distance t as one layer of
   *        coverage times opacity(t), without rounding the color first.
   * @param layers Front-to-back composite to add to.
   * @param t Ray distance.
   * @param coverage Layer coverage in [0, 1].
   */
  __attribute__((always_inline)) void composite(LayerComposite &layers, float t,
                                                float coverage) const {
    const float nearness = 1 - t * inv_far;
    float rgb[3];
    palette->get_color_unit_scaled(nearness, gain * (.45f + .55f * nearness),
                                   rgb);
    layers.add(rgb, coverage * opacity(t));
  }
};

/** @brief Shaded colour of a ray together with its trace result. */
struct ShadedTrace {
  Color4 color;      ///< Composited colour and alpha.
  TraceResult trace; ///< Status, counters and surface hit of the trace.
};

/**
 * @brief Composites every contribution of `trace_events` front to back.
 * @tparam Adapter Candidate stream adapter accepted by `trace_events`.
 * @param adapter Candidate streams to merge.
 * @param interval Ray-parameter range to trace.
 * @param limits Work budgets.
 * @param appearance Depth appearance applied to each layer.
 * @return Composited colour and the trace result; stops once saturated.
 */
template <typename Adapter>
__attribute__((always_inline)) inline ShadedTrace
shade_events(Adapter &adapter, Interval interval, const TraceLimits &limits,
             const Appearance &appearance) {
  LayerComposite composite;
  auto trace =
      trace_events(adapter, interval, limits,
                   [&](const Contribution &hit) __attribute__((always_inline)) {
                     HS_PROFILE_DEEP(hl_layer_composite);
                     appearance.composite(composite, hit.t, hit.coverage);
                     return !composite.saturated();
                   });
  return {composite.finish(), trace};
}

/** @brief Shades one verified boundary.
 * @details The footprint is validated and reserved for coverage-aware queries.
 * @tparam Query Distance query accepted by `surface_search`.
 * @param query Field to search.
 * @param ray Ray and parameter interval.
 * @param footprint Pixel cone footprint.
 * @param limits Step, refinement, and tolerance budgets.
 * @param appearance Depth appearance applied to the hit.
 * @return Shaded hit colour and the trace result.
 */
template <typename Query>
HS_HOT_FLASH_MEMBER ShadedTrace shade_surface(const Query &query,
                                              const Ray &ray,
                                              const Footprint &footprint,
                                              const TraceLimits &limits,
                                              const Appearance &appearance) {
  auto trace = surface_search(query, ray, footprint, limits);
  LayerComposite composite;
  if (trace.has_surface)
    appearance.composite(composite, trace.contribution.t,
                         trace.contribution.coverage);
  return {composite.finish(), trace};
}

/**
 * @brief Averages verified single-surface subrays; layered traces are excluded.
 * @tparam COUNT Subray count.
 * @tparam Trace Callable `ShadedTrace(const math::Vector&)`.
 * @param directions Subray directions.
 * @param trace Shades one subray.
 * @return Alpha-weighted mean colour, summed counters, and hit-fraction
 * coverage.
 */
template <size_t COUNT, typename Trace>
ShadedTrace verified_filter(const std::array<math::Vector, COUNT> &directions,
                            Trace trace) {
  static_assert(COUNT > 0);
  ShadedTrace result{};
  float red = 0, green = 0, blue = 0, alpha = 0;
  size_t hits = 0;
  for (const auto &direction : directions) {
    const auto sample = trace(direction);
    result.trace.counters.queries += sample.trace.counters.queries;
    result.trace.counters.steps += sample.trace.counters.steps;
    result.trace.counters.refinements += sample.trace.counters.refinements;
    result.trace.counters.candidates += sample.trace.counters.candidates;
    result.trace.counters.layers += sample.trace.counters.layers;
    if (sample.trace.status != TraceStatus::SURFACE &&
        sample.trace.status != TraceStatus::RANGE_COMPLETE &&
        result.trace.status != TraceStatus::INVALID_QUERY)
      result.trace.status = sample.trace.status;
    if (!sample.trace.has_surface || !sample.trace.contribution.verified)
      continue;
    if (!result.trace.has_surface)
      result.trace.contribution = sample.trace.contribution;
    result.trace.has_surface = true;
    ++hits;
    red += sample.color.color.r * sample.color.alpha;
    green += sample.color.color.g * sample.color.alpha;
    blue += sample.color.color.b * sample.color.alpha;
    alpha += sample.color.alpha;
  }
  result.trace.contribution.coverage = static_cast<float>(hits) / COUNT;
  if (result.trace.has_surface &&
      result.trace.status == TraceStatus::RANGE_COMPLETE)
    result.trace.status = TraceStatus::SURFACE;
  if (alpha > 0) {
    result.color = {{round_linear_channel(red / alpha),
                     round_linear_channel(green / alpha),
                     round_linear_channel(blue / alpha)},
                    alpha / COUNT};
  }
  return result;
}

} // namespace Raycast
