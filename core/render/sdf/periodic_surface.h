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

/** @brief Nodal field of a PeriodicSurface: Schwarz P cosine sum or gyroid. */
enum class PeriodicSurfaceKind { COSINE, GYROID };

/** @brief Boundary of a periodic nodal field's nonpositive region in 3D. */
template <PeriodicSurfaceKind Kind> struct PeriodicSurface {
  float period = 1.0f;   ///< World length of one field period; > 0.
  float iso = 0.0f;      ///< Field offset; the surface is field == 0.
  math::Vector origin{}; ///< World position of the field's phase origin.

  /** @brief True when every parameter is finite and `period` is positive.
   *  @return Whether the surface can be traced. */
  bool valid() const {
    return Raycast::finite(period) && period > 0.0f && Raycast::finite(iso) &&
           Raycast::finite(origin.x) && Raycast::finite(origin.y) &&
           Raycast::finite(origin.z) && Raycast::finite(lipschitz());
  }

  /** @brief Bound on the field gradient magnitude.
   *  @return sqrt(3) times the angular frequency 2*pi / `period`. */
  float lipschitz() const {
    constexpr float FACTOR = 1.7320508075688772f;
    return FACTOR * (math::TWO_PI_F / period);
  }

  /** @brief Nodal field value minus `iso`.
   *  @param p World point.
   *  @return Field value; nonpositive inside. */
  float field(const math::Vector &p) const {
    const math::Vector Q = (p - origin) * (math::TWO_PI_F / period);
    if constexpr (Kind == PeriodicSurfaceKind::COSINE)
      return cosf(Q.x) + cosf(Q.y) + cosf(Q.z) - iso;
    else
      return sinf(Q.x) * cosf(Q.y) + sinf(Q.y) * cosf(Q.z) +
             sinf(Q.z) * cosf(Q.x) - iso;
  }

  /** @brief Signed conservative clearance, not an exact signed distance.
   *  @param p World point.
   *  @return field() divided by its Lipschitz bound. */
  float distance(const math::Vector &p) const { return field(p) / lipschitz(); }

  /** @brief Analytic field gradient.
   *  @param p World point.
   *  @return Gradient of field() in world units. */
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

  /** @brief Unit outward normal.
   *  @param p World point.
   *  @return Normalized gradient(), or zero where the gradient vanishes. */
  math::Vector normal(const math::Vector &p) const {
    const math::Vector G = gradient(p);
    const float LENGTH = G.magnitude();
    return LENGTH > 0.0f ? G * (1.0f / LENGTH) : math::Vector{};
  }

  /** @brief Query capabilities: clearance on both sides, exact verification.
   *  @return Capabilities with zero verification error. */
  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }

  /** @brief Field sample with distance() as value and its magnitude as step.
   *  @param p World point.
   *  @return Sample with zero material and feature; boundary only at zero. */
  Raycast::QuerySample sample(const math::Vector &p) const {
    const float VALUE = distance(p);
    return {VALUE, fabsf(VALUE), VALUE == 0.0f, 0, 0};
  }
};

/// Schwarz P surface: cos x + cos y + cos z = iso.
using CosineSurface = PeriodicSurface<PeriodicSurfaceKind::COSINE>;
/// Gyroid surface: sin x cos y + sin y cos z + sin z cos x = iso.
using GyroidSurface = PeriodicSurface<PeriodicSurfaceKind::GYROID>;

} // namespace SDF
