/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file periodic_surface.h
 * @brief Periodic implicit surface tracing. */

#include <cmath>
#include "core/math/3dmath.h"
#include "core/render/ray/contract.h"

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
           Raycast::finite(origin) && Raycast::finite(lipschitz());
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
    const math::Vector q = (p - origin) * (math::TWO_PI_F / period);
    if constexpr (Kind == PeriodicSurfaceKind::COSINE)
      return cosf(q.x) + cosf(q.y) + cosf(q.z) - iso;
    else
      return sinf(q.x) * cosf(q.y) + sinf(q.y) * cosf(q.z) +
             sinf(q.z) * cosf(q.x) - iso;
  }

  /** @brief Signed conservative clearance, not an exact signed distance.
   *  @param p World point.
   *  @return field() divided by its Lipschitz bound. */
  float distance(const math::Vector &p) const { return field(p) / lipschitz(); }

  /** @brief Analytic field gradient.
   *  @param p World point.
   *  @return Gradient of field() in world units. */
  math::Vector gradient(const math::Vector &p) const {
    const float k = math::TWO_PI_F / period;
    const math::Vector q = (p - origin) * k;
    if constexpr (Kind == PeriodicSurfaceKind::COSINE)
      return math::Vector(-sinf(q.x), -sinf(q.y), -sinf(q.z)) * k;
    else
      return math::Vector(cosf(q.x) * cosf(q.y) - sinf(q.z) * sinf(q.x),
                          cosf(q.y) * cosf(q.z) - sinf(q.x) * sinf(q.y),
                          cosf(q.z) * cosf(q.x) - sinf(q.y) * sinf(q.z)) *
             k;
  }

  /** @brief Unit outward normal.
   *  @param p World point.
   *  @return Normalized gradient(), or zero where the gradient vanishes. */
  math::Vector normal(const math::Vector &p) const {
    const math::Vector g = gradient(p);
    const float length = g.magnitude();
    return length > 0.0f ? g * (1.0f / length) : math::Vector{};
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
    const float value = distance(p);
    return {value, fabsf(value), value == 0.0f, 0, 0};
  }
};

/// Schwarz P surface: cos x + cos y + cos z = iso.
using CosineSurface = PeriodicSurface<PeriodicSurfaceKind::COSINE>;
/// Gyroid surface: sin x cos y + sin y cos z + sin z cos x = iso.
using GyroidSurface = PeriodicSurface<PeriodicSurfaceKind::GYROID>;

} // namespace SDF
