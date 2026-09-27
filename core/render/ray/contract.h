/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "math/3dmath.h"

namespace Raycast {

inline bool finite(float value) {
  uint32_t bits;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return (bits & 0x7f800000u) != 0x7f800000u;
}

struct Interval {
  float near = 0.0f;
  float far = 10.0f;
  bool valid() const {
    return finite(near) && finite(far) && near >= 0.0f && far > near;
  }
};

inline bool finite(const math::Vector &p) {
  return finite(p.x) && finite(p.y) && finite(p.z);
}

struct Ray {
  math::Vector origin;
  math::Vector direction;
  Interval interval;
  math::Vector at(float t) const { return origin + direction * t; }
  bool valid() const {
    return interval.valid() && finite(origin) && finite(direction) &&
           fabsf(math::dot(direction, direction) - 1.0f) < 1e-4f;
  }
};

struct Footprint {
  float angular_radius = 0.0f;
  float radial_start = 0.0f;
  float at(float t) const { return (radial_start + t) * angular_radius; }
};

enum class TraceStatus {
  SURFACE,
  RANGE_COMPLETE,
  SATURATED,
  BUDGET_EXHAUSTED,
  UNRESOLVED,
  UNSUPPORTED_START,
  INVALID_QUERY
};

struct TraceLimits {
  int max_steps = 64;
  int max_refinements = 12;
  int max_candidates = 128;
  int max_layers = 32;
  int max_queries = 96;
  float position_tolerance = 1e-4f;
};

struct TraceCounters {
  int queries = 0;
  int steps = 0;
  int refinements = 0;
  int candidates = 0;
  int layers = 0;
};

struct Contribution {
  float t = 0.0f;
  float coverage = 1.0f;
  uint32_t material = 0;
  uint32_t feature = 0;
  uint64_t merge_identity = 0;
  bool verified = false;
  bool has_normal = false;
  math::Vector normal;
};

struct TraceResult {
  TraceStatus status = TraceStatus::RANGE_COMPLETE;
  TraceCounters counters;
  Contribution contribution;
  bool has_surface = false;
};

/** @brief Independent world-distance clearance and membership guarantees. */
struct QueryCapabilities {
  bool exterior_clearance = true;
  bool interior_clearance = false;
  bool surface_verification = true;
  float error = 0.0f;
};

struct QuerySample {
  float field = 0.0f;
  float clearance = 0.0f;
  bool boundary = false;
  uint32_t material = 0;
  uint32_t feature = 0;
};

} // namespace Raycast
