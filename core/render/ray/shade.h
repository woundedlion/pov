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
  float inv_far = .1f;
  float near_start = 0;
  float near_inv_span = 1;
  const BakedPalette *palette = nullptr;
  float gain = 1; /**< Brightness scale applied to every layer's color. */

  /** @brief Whether an interval keeps depth palette lookups in [0, 1]. */
  bool valid_for(Interval interval) const {
    return interval.valid() && finite(inv_far) && inv_far >= 0.0f &&
           interval.far * inv_far <= 1.0f;
  }

  __attribute__((always_inline)) float opacity(float t) const {
    const float fog = fmaxf(0.0f, 1.0f - t * inv_far);
    return fog * fog * math::cubic_kernel((t - near_start) * near_inv_span);
  }
  /** @brief Depth-graded palette color at distance t. */
  __attribute__((always_inline)) Pixel color(float t) const {
    const float nearness = 1 - t * inv_far;
    return palette->get_color_unit(nearness) *
           (gain * (.45f + .55f * nearness));
  }
  /**
   * @brief Adds the depth-graded color at distance t as one layer of
   *        coverage times opacity(t), without rounding the color first.
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

struct ShadedTrace {
  Color4 color;
  TraceResult trace;
};

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

/** @brief Averages verified single-surface subrays; layered traces are excluded. */
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
