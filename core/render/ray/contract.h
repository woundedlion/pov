/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file contract.h
 * @brief Ray tracing query, budget and result contracts. */

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

/** @brief Completion reason for a bounded ray trace. */
enum class TraceStatus {
  SURFACE,          /**< Verified surface found. */
  RANGE_COMPLETE,   /**< Entire requested interval traversed. */
  SATURATED,        /**< Contribution consumer stopped accepting layers. */
  BUDGET_EXHAUSTED, /**< A tracing work limit was reached. */
  UNRESOLVED, /**< Search could not establish a safe next step or boundary. */
  UNSUPPORTED_START, /**< Query cannot trace from the initial membership. */
  INVALID_QUERY      /**< Invalid ray, query or sampling contract. */
};

/** @brief Per-trace work budgets and world-distance refinement tolerance. */
struct TraceLimits {
  int max_steps = 64;       /**< Maximum marching or traversal steps. */
  int max_refinements = 12; /**< Maximum boundary refinement iterations. */
  int max_candidates = 128; /**< Maximum candidate events inspected. */
  int max_layers = 32;      /**< Maximum composited contributions. */
  int max_queries = 96;     /**< Maximum distance query evaluations. */
  float position_tolerance = 1e-4f; /**< Boundary tolerance in world units. */
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

/** @brief Trace completion status and optional verified single-surface result. */
struct TraceResult {
  TraceStatus status = TraceStatus::RANGE_COMPLETE;
  TraceCounters counters;
  Contribution contribution;
  bool has_surface =
      false; /**< A verified single surface, not layered coverage. */
};

/** @brief Independent world-distance clearance and membership guarantees. */
struct QueryCapabilities {
  bool exterior_clearance =
      true; /**< Safe exterior distance bound available. */
  bool interior_clearance =
      false; /**< Safe interior distance bound available. */
  bool surface_verification = true; /**< Boundary membership can be verified. */
  float error =
      0.0f; /**< World-distance uncertainty of surface verification. */
};

/** @brief Signed field sample with an independent safe stepping distance. */
struct QuerySample {
  float field = 0.0f;     /**< Signed field value; negative denotes interior. */
  float clearance = 0.0f; /**< Nonnegative safe step in world units. */
  bool boundary = false;  /**< Query explicitly verifies a boundary here. */
  uint32_t material = 0;  /**< Material identity at the sample. */
  uint32_t feature = 0;   /**< Geometric feature identity at the sample. */
};

} // namespace Raycast
