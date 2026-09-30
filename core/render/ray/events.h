/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file events.h
 * @brief Ordered ray contribution streams and bounded tracing. */

#include <array>
#include <cmath>
#include "render/ray/contract.h"

namespace Raycast {

/**
 * @brief Merges bounded monotone candidate streams in deterministic order.
 * @details Adapter::GROUP_CAPACITY optionally bounds distinct merge identities
 * in one tolerance window; the default is STREAM_COUNT. Overflow retains the
 * buffered contributions and returns BUDGET_EXHAUSTED.
 */
template <typename Adapter, typename Consume>
__attribute__((always_inline)) inline TraceResult
trace_events(Adapter &adapter, Interval interval, const TraceLimits &limits,
             Consume consume, float relative_tolerance = 1.0e-4f) {
  TraceResult result;
  if (!interval.valid() || !finite(relative_tolerance) ||
      relative_tolerance < 0) {
    result.status = TraceStatus::INVALID_QUERY;
    return result;
  }
  constexpr size_t COUNT = Adapter::STREAM_COUNT;
  constexpr size_t CAPACITY = [] {
    if constexpr (requires { Adapter::GROUP_CAPACITY; })
      return Adapter::GROUP_CAPACITY;
    else
      return COUNT;
  }();
  static_assert(CAPACITY > 0 && CAPACITY <= COUNT);
  std::array<Contribution, CAPACITY> group;
  size_t grouped = 0;
  float group_end = 0;
  const auto flush = [&]() __attribute__((always_inline)) {
    for (size_t i = 0; i < grouped; ++i) {
      if (result.counters.layers >= limits.max_layers) {
        if (result.status != TraceStatus::INVALID_QUERY)
          result.status = TraceStatus::BUDGET_EXHAUSTED;
        return false;
      }
      ++result.counters.layers;
      if (!consume(group[i])) {
        if (result.status != TraceStatus::INVALID_QUERY)
          result.status = TraceStatus::SATURATED;
        return false;
      }
    }
    grouped = 0;
    return true;
  };
  while (true) {
    HS_PROFILE_DEEP(hl_event_step);
    float nearest = interval.far;
    size_t first = COUNT;
    for (size_t i = 0; i < COUNT; ++i) {
      if (!adapter.active(i))
        continue;
      const float t = adapter.distance(i);
      if (!finite(t)) {
        result.status = TraceStatus::INVALID_QUERY;
        first = COUNT;
        break;
      }
      if (t < nearest || (first == COUNT && t == interval.far)) {
        nearest = t;
        first = i;
      }
    }
    if (first == COUNT)
      break;
    if (grouped && nearest > group_end && !flush())
      return result;
    if (result.counters.candidates >= limits.max_candidates) {
      result.status = TraceStatus::BUDGET_EXHAUSTED;
      break;
    }
    ++result.counters.candidates;
    if (nearest >= interval.near) {
      Contribution candidate = adapter.candidate(first);
      if (!finite(candidate.t) || candidate.t != nearest ||
          !finite(candidate.coverage) || candidate.coverage < 0 ||
          candidate.coverage > 1 ||
          (candidate.has_normal && !finite(candidate.normal))) {
        result.status = TraceStatus::INVALID_QUERY;
        break;
      }
      if (candidate.coverage > 0) {
        if (!grouped)
          group_end = nearest + relative_tolerance * std::max(1.0f, nearest);
        size_t slot = 0;
        while (slot < grouped &&
               (group[slot].merge_identity != candidate.merge_identity ||
                group[slot].material != candidate.material ||
                group[slot].verified != candidate.verified))
          ++slot;
        if (slot == CAPACITY) {
          result.status = TraceStatus::BUDGET_EXHAUSTED;
          break;
        }
        if (slot == grouped)
          group[grouped++] = candidate;
        else if (candidate.coverage > group[slot].coverage) {
          const float first_t = group[slot].t;
          group[slot] = candidate;
          group[slot].t = first_t;
        }
      }
    }
    adapter.advance(first);
    if (adapter.active(first) && !(adapter.distance(first) > nearest)) {
      result.status = TraceStatus::INVALID_QUERY;
      break;
    }
  }
  flush();
  return result;
}

} // namespace Raycast
