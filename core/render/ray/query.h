/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file query.h
 * @brief Distance query adapters and placement transforms. */

#include "render/ray/camera.h"

namespace Raycast {

/** @brief Adapts a declared distance bound without changing the shape API. */
template <typename Shape> struct VolumeQuery {
  const Shape &shape;           ///< Borrowed shape; must outlive the query.
  QueryCapabilities guarantees; ///< Declared bounds of `shape`'s distance.

  /**
   * @brief Forwards the shape's own `valid()`, if it has one.
   * @return True when the shape is valid or has no check.
   */
  bool valid() const {
    if constexpr (requires { shape.valid(); })
      return shape.valid();
    return true;
  }
  /**
   * @brief Declared distance guarantees.
   * @return `guarantees`.
   */
  QueryCapabilities capabilities() const { return guarantees; }
  /** @brief Forwards the shape's trace precondition check, if it has one. */
  void check_trace_preconditions() const {
    if constexpr (requires { shape.check_trace_preconditions(); })
      shape.check_trace_preconditions();
  }
  /**
   * @brief Samples the shape's field and safe clearance at a point.
   * @param p Shape-space point.
   * @return Sample whose clearance honours `guarantees`; material and feature
   * are zero.
   */
  QuerySample sample(const math::Vector &p) const {
    const float DISTANCE = shape.distance(p);
    float field = DISTANCE;
    if constexpr (requires { shape.raw_distance(p); }) {
      // A cheap exterior bound cannot certify a boundary.
      if (DISTANCE <= 0.0f)
        field = shape.raw_distance(p);
    }
    return {field,
            field < 0.0f ? (guarantees.interior_clearance ? -DISTANCE : 0.0f)
                         : fmaxf(DISTANCE, 0.0f),
            field == 0.0f && guarantees.surface_verification, 0, 0};
  }
};

/** @brief Uniform placement preserving world-distance clearance units. */
template <typename Query> struct PlacedQuery {
  const Query &query;                ///< Borrowed query in local units.
  math::Vector center;               ///< World position of the local origin.
  math::Quaternion inverse_rotation; ///< Unit world-to-local rotation.
  float scale = 1.0f;                ///< Uniform local-to-world scale; > 0.

  /**
   * @brief Whether the placement is finite with a unit rotation and positive
   * scale, and the inner query is valid.
   * @return True for a usable placement.
   */
  bool valid() const {
    if constexpr (requires { query.valid(); })
      if (!query.valid())
        return false;
    return finite(center) && finite(scale) && scale > 0.0f &&
           finite(inverse_rotation.r) && finite(inverse_rotation.v) &&
           fabsf(inverse_rotation.r * inverse_rotation.r +
                 math::dot(inverse_rotation.v, inverse_rotation.v) - 1.0f) <
               1e-4f;
  }
  /** @brief Forwards the inner query's trace precondition check, if any. */
  void check_trace_preconditions() const {
    if constexpr (requires { query.check_trace_preconditions(); })
      query.check_trace_preconditions();
  }
  /**
   * @brief Inner guarantees with the error scaled to world units.
   * @return The inner capabilities, `error` multiplied by `scale`.
   */
  QueryCapabilities capabilities() const {
    auto result = query.capabilities();
    result.error *= scale;
    return result;
  }
  /** @brief Samples a placement for which valid() returned true. */
  QuerySample sample(const math::Vector &p) const {
    auto result = query.sample(
        math::rotate((p - center) * (1.0f / scale), inverse_rotation));
    result.field *= scale;
    result.clearance *= scale;
    return result;
  }
};

/**
 * @brief Samples a 3D or 4D query through a camera's ambient embedding.
 * @tparam Query Query sampled at ambient points.
 * @tparam FOUR_DIMENSIONAL Whether `Query` takes 4D points; must match the
 * camera's `SamplingDomain`.
 */
template <typename Query, bool FOUR_DIMENSIONAL> struct DomainQuery {
  const Query &query;           ///< Borrowed ambient-space query.
  const PreparedCamera &camera; ///< Borrowed camera embedding.

  /**
   * @brief Whether the camera is valid, matches `FOUR_DIMENSIONAL` and the
   * query is valid.
   * @return True for a usable domain query.
   */
  bool valid() const {
    if (!camera.valid() ||
        (camera.domain == SamplingDomain::SLICE_4D) != FOUR_DIMENSIONAL)
      return false;
    if constexpr (requires { query.valid(); })
      return query.valid();
    return true;
  }
  /**
   * @brief Inner query guarantees; the embedding is an isometry.
   * @return `query.capabilities()`.
   */
  QueryCapabilities capabilities() const { return query.capabilities(); }
  /** @brief Forwards the inner query's trace precondition check, if any. */
  void check_trace_preconditions() const {
    if constexpr (requires { query.check_trace_preconditions(); })
      query.check_trace_preconditions();
  }
  /** @brief Samples a domain for which valid() returned true. */
  QuerySample sample(const math::Vector &p) const {
    if constexpr (FOUR_DIMENSIONAL)
      return query.sample(camera.point4(p));
    else
      return query.sample(camera.point3(p));
  }
};

/// `DomainQuery` over a 3D ambient query.
template <typename Query> using DomainQuery3 = DomainQuery<Query, false>;
/// `DomainQuery` over a 4D ambient query.
template <typename Query> using DomainQuery4 = DomainQuery<Query, true>;

} // namespace Raycast
