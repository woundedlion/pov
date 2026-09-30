/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file periodic_surface.h
 * @brief Periodic implicit surface tracing. */

#include <cmath>
#include "math/3dmath.h"
#include "render/ray/contract.h"

namespace SDF {

enum class PeriodicSurfaceKind { COSINE, GYROID };

/** @brief Boundary of a periodic nodal field's nonpositive region in 3D. */
template <PeriodicSurfaceKind Kind> struct PeriodicSurface {
  float period = 1.0f;
  float iso = 0.0f;
  math::Vector origin{};

  bool valid() const {
    return Raycast::finite(period) && period > 0.0f && Raycast::finite(iso) &&
           Raycast::finite(origin.x) && Raycast::finite(origin.y) &&
           Raycast::finite(origin.z) && Raycast::finite(lipschitz());
  }

  float lipschitz() const {
    constexpr float FACTOR = Kind == PeriodicSurfaceKind::COSINE
                                 ? 1.7320508075688772f
                                 : 3.4641016151377544f;
    return FACTOR * (math::TWO_PI_F / period);
  }

  float field(const math::Vector &p) const {
    const math::Vector Q = (p - origin) * (math::TWO_PI_F / period);
    if constexpr (Kind == PeriodicSurfaceKind::COSINE)
      return cosf(Q.x) + cosf(Q.y) + cosf(Q.z) - iso;
    else
      return sinf(Q.x) * cosf(Q.y) + sinf(Q.y) * cosf(Q.z) +
             sinf(Q.z) * cosf(Q.x) - iso;
  }

  /** @brief Signed conservative clearance, not an exact signed distance. */
  float distance(const math::Vector &p) const { return field(p) / lipschitz(); }

  math::Vector gradient(const math::Vector &p) const {
    const float K = math::TWO_PI_F / period;
    const math::Vector Q = (p - origin) * K;
    if constexpr (Kind == PeriodicSurfaceKind::COSINE)
      return math::Vector(-sinf(Q.x), -sinf(Q.y), -sinf(Q.z)) * K;
    else
      return math::Vector(cosf(Q.x) * cosf(Q.y) - sinf(Q.z) * sinf(Q.x),
                          cosf(Q.y) * cosf(Q.z) - sinf(Q.x) * sinf(Q.y),
                          cosf(Q.z) * cosf(Q.x) - sinf(Q.y) * sinf(Q.z)) *
             K;
  }

  math::Vector normal(const math::Vector &p) const {
    const math::Vector G = gradient(p);
    const float LENGTH = G.magnitude();
    return LENGTH > 0.0f ? G * (1.0f / LENGTH) : math::Vector{};
  }

  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }

  Raycast::QuerySample sample(const math::Vector &p) const {
    const float VALUE = distance(p);
    return {VALUE, fabsf(VALUE), VALUE == 0.0f, 0, 0};
  }
};

using CosineSurface = PeriodicSurface<PeriodicSurfaceKind::COSINE>;
using GyroidSurface = PeriodicSurface<PeriodicSurfaceKind::GYROID>;

} // namespace SDF
