/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cfloat>
#include "render/ray/contract.h"

namespace Raycast {

/** @brief Shared bounded progress driver; false ends the current policy. */
template <typename Step>
__attribute__((always_inline)) inline void march_steps(int limit, Step &&step) {
  for (int i = 0; i < limit; ++i)
    if (!step())
      break;
}

__attribute__((always_inline)) inline math::Vector
advance(const math::Vector &p, const math::Vector &direction, float step) {
  return math::Vector(p.x + direction.x * step, p.y + direction.y * step,
                      p.z + direction.z * step);
}

HS_O3_BEGIN
/** @brief Finds the closest sampled approach along a sphere-traced ray. */
template <typename Shape>
__attribute__((always_inline)) inline float
closest_approach(const Shape &shape, const math::Vector &origin,
                 const math::Vector &direction, float end_t, int max_steps,
                 float aa_width, math::Vector &closest_local,
                 float overrelaxation = 1.3f) {
  float t = 0.0f;
  math::Vector local_p = origin;
  closest_local = origin;
  float closest_d = FLT_MAX;
  float omega = overrelaxation;
  float prev_r = 0.0f;
  float step_len = 0.0f;
  march_steps(max_steps, [&]() {
    if (t > end_t)
      return false;
    float d = shape.distance(local_p);
    float r = d < 0.0f ? -d : d;
    if (omega > 1.0f && r + prev_r < step_len) {
      float back = prev_r - step_len;
      t += back;
      local_p = advance(local_p, direction, back);
      omega = 1.0f;
      prev_r = 0.0f;
      step_len = 0.0f;
      return true;
    }
    prev_r = r;
    if (d < closest_d) {
      closest_d = d;
      closest_local = local_p;
      if (closest_d <= aa_width * 0.02f)
        return false;
    } else if (closest_d < aa_width) {
      return false;
    }
    if (d < -aa_width)
      return false;
    step_len = std::max(d * 0.9f * omega, 1e-5f);
    t += step_len;
    local_p = advance(local_p, direction, step_len);
    return true;
  });
  return closest_d;
}
HS_O3_END

inline bool valid_sample(const QuerySample &sample) {
  return finite(sample.field) && finite(sample.clearance) &&
         sample.clearance >= 0.0f;
}

/** @brief Finds one verified boundary; proximity never creates coverage.
 * @param footprint Validated footprint reserved for coverage-aware queries.
 */
template <typename Query>
HS_HOT_FLASH_MEMBER TraceResult surface_search(const Query &query,
                                               const Ray &ray,
                                               const Footprint &footprint,
                                               const TraceLimits &limits = {}) {
  TraceResult result;
  if constexpr (requires { query.valid(); }) {
    if (!query.valid()) {
      result.status = TraceStatus::INVALID_QUERY;
      return result;
    }
  }
  const auto GUARANTEES = query.capabilities();
  if (!ray.valid() || !finite(footprint.angular_radius) ||
      footprint.angular_radius < 0.0f || !finite(footprint.radial_start) ||
      footprint.radial_start < 0.0f || !finite(limits.position_tolerance) ||
      limits.position_tolerance <= 0.0f || !finite(GUARANTEES.error) ||
      GUARANTEES.error < 0.0f ||
      GUARANTEES.error >= limits.position_tolerance) {
    result.status = TraceStatus::INVALID_QUERY;
    return result;
  }
  if constexpr (requires { query.check_trace_preconditions(); })
    query.check_trace_preconditions();
  result.status = TraceStatus::BUDGET_EXHAUSTED;
  if (limits.max_layers <= 0)
    return result;
  float t = ray.interval.near;
  bool first = true;
  bool inside = false;
  auto emit = [&](float position, const QuerySample &sample) {
    result.status = TraceStatus::SURFACE;
    result.has_surface = true;
    result.counters.layers = 1;
    result.contribution.t = position;
    result.contribution.material = sample.material;
    result.contribution.feature = sample.feature;
    result.contribution.verified = true;
  };
  auto evaluate = [&](float position, QuerySample &sample) HS_HOT_FLASH_MEMBER {
    if (result.counters.queries >= limits.max_queries) {
      result.status = TraceStatus::BUDGET_EXHAUSTED;
      return false;
    }
    sample = query.sample(ray.at(position));
    ++result.counters.queries;
    if (!valid_sample(sample)) {
      result.status = TraceStatus::INVALID_QUERY;
      return false;
    }
    return true;
  };
  march_steps(limits.max_steps, [&]() {
    QuerySample sample;
    if (!evaluate(t, sample))
      return false;
    ++result.counters.steps;
    if (sample.boundary && GUARANTEES.surface_verification) {
      emit(t, sample);
      return false;
    }
    if (first) {
      inside = sample.field < 0.0f;
      first = false;
      if ((inside && !GUARANTEES.interior_clearance) ||
          (!inside && !GUARANTEES.exterior_clearance) ||
          !GUARANTEES.surface_verification) {
        result.status = TraceStatus::UNSUPPORTED_START;
        return false;
      }
    } else if ((sample.field < 0.0f) != inside) {
      result.status = TraceStatus::INVALID_QUERY;
      return false;
    }
    const float CLEARANCE = std::max(sample.clearance - GUARANTEES.error, 0.0f);
    const float REMAINING = ray.interval.far - t;
    if (CLEARANCE > REMAINING) {
      result.status = TraceStatus::RANGE_COMPLETE;
      return false;
    }
    if (REMAINING <= 0.0f) {
      result.status = sample.field == 0.0f ? TraceStatus::UNRESOLVED
                                           : TraceStatus::RANGE_COMPLETE;
      return false;
    }
    if (CLEARANCE <= limits.position_tolerance) {
      result.status = TraceStatus::UNRESOLVED;
      if (result.counters.refinements >= limits.max_refinements)
        return false;
      const float PROBE_T =
          std::min(t + limits.position_tolerance, ray.interval.far);
      if (!(PROBE_T > t))
        return false;
      QuerySample probe;
      ++result.counters.refinements;
      if (!evaluate(PROBE_T, probe))
        return false;
      // Only this tolerance-sized interval may contain uncleared roots.
      if (probe.boundary || ((probe.field < 0.0f) != inside)) {
        emit((t + PROBE_T) * 0.5f, probe);
        return false;
      }
    }
    const float STEP = CLEARANCE * 0.9f;
    const float NEXT_T = t + STEP;
    if (!(NEXT_T > t) || NEXT_T > ray.interval.far) {
      result.status = TraceStatus::UNRESOLVED;
      return false;
    }
    t = NEXT_T;
    result.status = TraceStatus::BUDGET_EXHAUSTED;
    return true;
  });
  return result;
}

} // namespace Raycast
