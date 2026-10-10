/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file march.h
 * @brief Bounded verified surface searches. */

#include <cfloat>
#include "render/ray/contract.h"

namespace Raycast {

/** @brief Fraction of a safe distance bound taken per step. */
inline constexpr float STEP_SAFETY = 0.9f;

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

/** @brief True when two consecutive unbounding spheres leave a gap, so the
 * step between their centres crossed unbounded space. */
__attribute__((always_inline)) inline bool
spheres_disjoint(float radius, float prev_radius, float step) {
  return radius + prev_radius < step;
}

HS_O3_BEGIN
/** @brief Finds the closest sampled approach along a sphere-traced ray.
 * @details Steps are overrelaxed by `overrelaxation` until one overshoots;
 * the march then rewinds to the previous sample's bound and continues
 * unrelaxed. Stops on a hit within a small fraction of `aa_width`, or at the
 * first non-approaching sample once the closest approach is within
 * `aa_width`.
 */
template <typename Shape>
__attribute__((always_inline)) inline float
closest_approach(const Shape &shape, const math::Vector &origin,
                 const math::Vector &direction, float end_t, int max_steps,
                 float aa_width, math::Vector &closest_local,
                 float overrelaxation = 1.3f) {
  constexpr float HIT_FRACTION = 0.02f;
  constexpr float MIN_STEP = 1e-5f;
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
    if (omega > 1.0f && spheres_disjoint(r, prev_r, step_len)) {
      const float REWIND = prev_r - step_len;
      t += REWIND;
      local_p = advance(local_p, direction, REWIND);
      omega = 1.0f;
      prev_r = 0.0f;
      step_len = 0.0f;
      return true;
    }
    prev_r = r;
    if (d < closest_d) {
      closest_d = d;
      closest_local = local_p;
      const bool HIT = closest_d <= aa_width * HIT_FRACTION;
      if (HIT)
        return false;
    } else {
      const bool LEAVING_GRAZE = closest_d < aa_width;
      if (LEAVING_GRAZE)
        return false;
    }
    if (d < -aa_width)
      return false;
    step_len = fmaxf(d * STEP_SAFETY * omega, MIN_STEP);
    t += step_len;
    local_p = advance(local_p, direction, step_len);
    return true;
  });
  return closest_d;
}
HS_O3_END

__attribute__((always_inline)) inline bool finite_nonnegative(float value) {
  return finite(value) && value >= 0.0f;
}

inline bool valid_sample(const QuerySample &sample) {
  return finite(sample.field) && finite_nonnegative(sample.clearance);
}

__attribute__((always_inline)) inline bool
valid_footprint(const Footprint &footprint) {
  return finite_nonnegative(footprint.angular_radius) &&
         finite_nonnegative(footprint.radial_start);
}

/** @brief True when `limits` has a positive tolerance that `guarantees`
 * verifies boundaries within. */
__attribute__((always_inline)) inline bool
tolerance_verifiable(const TraceLimits &limits,
                     const QueryCapabilities &guarantees) {
  return finite(limits.position_tolerance) &&
         limits.position_tolerance > 0.0f &&
         finite_nonnegative(guarantees.error) &&
         guarantees.error < limits.position_tolerance;
}

/** @brief True when `guarantees` can sphere-trace and verify from a start on
 * the given side of the boundary. */
__attribute__((always_inline)) inline bool
can_start(const QueryCapabilities &guarantees, bool inside) {
  return (inside ? guarantees.interior_clearance
                 : guarantees.exterior_clearance) &&
         guarantees.surface_verification;
}

/** @brief Finds one verified boundary; proximity never creates coverage.
 * @details The footprint is validated and reserved for coverage-aware queries.
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
  if (!ray.valid() || !valid_footprint(footprint) ||
      !tolerance_verifiable(limits, GUARANTEES)) {
    result.status = TraceStatus::INVALID_QUERY;
    return result;
  }
  if constexpr (requires { query.check_trace_preconditions(); })
    query.check_trace_preconditions();
  result.status = TraceStatus::BUDGET_EXHAUSTED;
  if (limits.max_layers <= 0)
    return result;
  float t = ray.interval.near;
  QuerySample sample;
  bool sample_ready = false;
  bool first = true;
  bool inside = false;
  auto stop = [&](TraceStatus status) __attribute__((always_inline)) {
    result.status = status;
    return false;
  };
  auto emit = [&](float position, const QuerySample &hit)
                  __attribute__((always_inline)) {
                    result.has_surface = true;
                    result.counters.layers = 1;
                    result.contribution.t = position;
                    result.contribution.material = hit.material;
                    result.contribution.feature = hit.feature;
                    result.contribution.verified = true;
                    return stop(TraceStatus::SURFACE);
                  };
  auto evaluate = [&](float position, QuerySample &out) HS_HOT_FLASH_MEMBER {
    if (result.counters.queries >= limits.max_queries)
      return stop(TraceStatus::BUDGET_EXHAUSTED);
    out = query.sample(ray.at(position));
    ++result.counters.queries;
    if (!valid_sample(out))
      return stop(TraceStatus::INVALID_QUERY);
    return true;
  };
  auto keeps_start_side = [&]() __attribute__((always_inline)) {
    const bool SAMPLE_INSIDE = sample.field < 0.0f;
    if (first) {
      first = false;
      inside = SAMPLE_INSIDE;
      if (!can_start(GUARANTEES, inside))
        return stop(TraceStatus::UNSUPPORTED_START);
      return true;
    }
    if (SAMPLE_INSIDE != inside)
      return stop(TraceStatus::INVALID_QUERY);
    return true;
  };
  // Only this tolerance-sized interval may contain uncleared roots.
  auto refine = [&]() __attribute__((always_inline)) {
    if (result.counters.refinements >= limits.max_refinements)
      return stop(TraceStatus::UNRESOLVED);
    const float PROBE_T =
        fminf(t + limits.position_tolerance, ray.interval.far);
    if (!(PROBE_T > t))
      return stop(TraceStatus::UNRESOLVED);
    QuerySample probe;
    ++result.counters.refinements;
    if (!evaluate(PROBE_T, probe))
      return false;
    if (probe.boundary || ((probe.field < 0.0f) != inside))
      return emit((t + PROBE_T) * 0.5f, probe);
    t = PROBE_T;
    sample = probe;
    sample_ready = true;
    return true;
  };
  auto sphere_step = [&](float clearance) __attribute__((always_inline)) {
    const float NEXT_T = t + clearance * STEP_SAFETY;
    if (!(NEXT_T > t) || NEXT_T > ray.interval.far)
      return stop(TraceStatus::UNRESOLVED);
    t = NEXT_T;
    return true;
  };
  march_steps(limits.max_steps, [&]() __attribute__((always_inline)) {
    if (!sample_ready && !evaluate(t, sample))
      return false;
    sample_ready = false;
    ++result.counters.steps;
    if (sample.boundary && GUARANTEES.surface_verification)
      return emit(t, sample);
    if (!keeps_start_side())
      return false;
    const float CLEARANCE = fmaxf(sample.clearance - GUARANTEES.error, 0.0f);
    const float REMAINING = ray.interval.far - t;
    if (CLEARANCE > REMAINING)
      return stop(TraceStatus::RANGE_COMPLETE);
    if (REMAINING <= 0.0f)
      return stop(sample.field == 0.0f ? TraceStatus::UNRESOLVED
                                       : TraceStatus::RANGE_COMPLETE);
    if (CLEARANCE <= limits.position_tolerance)
      return refine();
    return sphere_step(CLEARANCE);
  });
  return result;
}

} // namespace Raycast
