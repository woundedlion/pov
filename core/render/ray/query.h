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
  const Shape &shape;
  QueryCapabilities guarantees;

  bool valid() const {
    if constexpr (requires { shape.valid(); })
      return shape.valid();
    return true;
  }
  QueryCapabilities capabilities() const { return guarantees; }
  void check_trace_preconditions() const {
    if constexpr (requires { shape.check_trace_preconditions(); })
      shape.check_trace_preconditions();
  }
  QuerySample sample(const math::Vector &p) const {
    const float DISTANCE = shape.distance(p);
    float field = DISTANCE;
    if constexpr (requires { shape.raw_distance(p); }) {
      // A cheap exterior bound cannot certify a boundary.
      if (DISTANCE <= 0.0f)
        field = shape.raw_distance(p);
    }
    return {field,
            field < 0.0f && !guarantees.interior_clearance ? 0.0f
                                                           : fabsf(DISTANCE),
            field == 0.0f && guarantees.surface_verification, 0, 0};
  }
};

/** @brief Uniform placement preserving world-distance clearance units. */
template <typename Query> struct PlacedQuery {
  const Query &query;
  math::Vector center;
  math::Quaternion inverse_rotation;
  float scale = 1.0f;

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
  void check_trace_preconditions() const {
    if constexpr (requires { query.check_trace_preconditions(); })
      query.check_trace_preconditions();
  }
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

template <typename Query, bool FOUR_DIMENSIONAL> struct DomainQuery {
  const Query &query;
  const PreparedCamera &camera;

  bool valid() const {
    if (!camera.valid() ||
        (camera.domain == SamplingDomain::SLICE_4D) != FOUR_DIMENSIONAL)
      return false;
    if constexpr (requires { query.valid(); })
      return query.valid();
    return true;
  }
  QueryCapabilities capabilities() const { return query.capabilities(); }
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

template <typename Query> using DomainQuery3 = DomainQuery<Query, false>;
template <typename Query> using DomainQuery4 = DomainQuery<Query, true>;

} // namespace Raycast
