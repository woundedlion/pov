/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file contract.h
 * @brief Ray tracing query, budget and result contracts. */

#include "math/3dmath.h"

namespace Raycast {

/**
 * @brief Bitwise finiteness test, immune to fast-math folding.
 * @param value Value to test.
 * @return False for infinities and NaNs.
 */
inline bool finite(float value) {
  uint32_t bits;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return (bits & 0x7f800000u) != 0x7f800000u;
}

/** @brief Ray-parameter range [near, far] in world units. */
struct Interval {
  float near = 0.0f; ///< Start distance; nonnegative.
  float far = 10.0f; ///< End distance; greater than `near`.
  /**
   * @brief Whether both ends are finite and 0 <= near < far.
   * @return True for a traceable interval.
   */
  bool valid() const {
    return finite(near) && finite(far) && near >= 0.0f && far > near;
  }
};

/**
 * @brief Whether every component of a vector is finite.
 * @param p Vector to test.
 * @return False if any component is infinite or NaN.
 */
inline bool finite(const math::Vector &p) {
  return finite(p.x) && finite(p.y) && finite(p.z);
}

/** @brief Parametric ray with a unit direction and a traced interval. */
struct Ray {
  math::Vector origin;    ///< Position at t = 0.
  math::Vector direction; ///< Unit-length direction.
  Interval interval;      ///< Parameter range to trace.
  /**
   * @brief Point at parameter t.
   * @param t Distance along `direction` from `origin`.
   * @return `origin + direction * t`.
   */
  math::Vector at(float t) const { return origin + direction * t; }
  /**
   * @brief Whether the ray is finite, unit-direction and has a valid interval.
   * @return True for a traceable ray.
   */
  bool valid() const {
    return interval.valid() && finite(origin) && finite(direction) &&
           fabsf(math::dot(direction, direction) - 1.0f) < 1e-4f;
  }
};

/** @brief Cone footprint of a pixel ray, growing linearly with distance. */
struct Footprint {
  float angular_radius = 0.0f; ///< Cone half-angle in radians; nonnegative.
  float radial_start = 0.0f;   ///< Distance from the cone apex to t = 0.
  /**
   * @brief Footprint radius at parameter t.
   * @param t Ray parameter.
   * @return World-unit radius `(radial_start + t) * angular_radius`.
   */
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

/** @brief Work spent by one trace, counted against `TraceLimits`. */
struct TraceCounters {
  int queries = 0;     ///< Distance query evaluations.
  int steps = 0;       ///< Marching or traversal steps.
  int refinements = 0; ///< Boundary refinement iterations.
  int candidates = 0;  ///< Candidate events inspected.
  int layers = 0;      ///< Contributions passed to the consumer.
};

/** @brief One surface or layer hit along a ray. */
struct Contribution {
  float t = 0.0f;              ///< Ray parameter of the hit.
  float coverage = 1.0f;       ///< Fractional coverage in [0, 1].
  uint32_t material = 0;       ///< Material identity at the hit.
  uint32_t feature = 0;        ///< Geometric feature identity at the hit.
  uint64_t merge_identity = 0; ///< Key under which near-coincident hits merge.
  bool verified = false;       ///< Boundary membership was verified.
  bool has_normal = false;     ///< Whether `normal` is set.
  math::Vector normal;         ///< Unit surface normal when `has_normal`.
};

/** @brief Trace completion status and optional verified single-surface result. */
struct TraceResult {
  TraceStatus status = TraceStatus::RANGE_COMPLETE; ///< Completion reason.
  TraceCounters counters;    ///< Work spent by the trace.
  Contribution contribution; ///< The surface hit when `has_surface`.
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
