/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Dense volume-ray reference model.
 */
#pragma once

#include <algorithm>
#include <cfloat>
#include "core/math/3dmath.h"
#include "core/render/scan/volume.h"

namespace hs_test::scan_tests {

/** @brief Dense fixed-step oracle, independent of sphere-tracing policy. */
struct VolumeReference {
  /** @brief A sampled position and signed distance. */
  struct Sample {
    math::Vector position; ///< Sample position.
    float distance;        ///< Signed distance at the position.
  };

  /** @brief Visits dense fixed-step samples until the visitor stops. */
  template <typename Shape, typename Visit>
  static void march(const Shape &shape, const math::Vector &origin,
                    const math::Vector &direction, float radius, Visit visit,
                    int steps = 16384) {
    const float end = radius - math::dot(origin, direction);
    for (int i = 0; i <= steps; ++i) {
      const float t = end * static_cast<float>(i) / steps;
      const math::Vector p(origin.x + direction.x * t,
                           origin.y + direction.y * t,
                           origin.z + direction.z * t);
      if (!visit(Sample{p, shape.distance(p)}))
        break;
    }
  }

  /** @brief Evaluates smooth coverage from signed distance. */
  static float coverage(float distance, float threshold, float width) {
    const double x =
        std::clamp((width - distance) / (width - threshold), 0.0f, 1.0f);
    return static_cast<float>(x * x * x * (10.0 + x * (-15.0 + 6.0 * x)));
  }

  /** @brief Finds the first near-surface distance minimum. */
  template <typename Shape>
  static Sample trace(const Shape &shape, const math::Vector &origin,
                      const math::Vector &direction, float radius, float aa) {
    Sample closest{origin, FLT_MAX};
    march(shape, origin, direction, radius, [&](Sample sample) {
      if (sample.distance < closest.distance)
        closest = sample;
      else if (closest.distance < aa)
        return false;
      return closest.distance > 0.02f * aa;
    });
    return closest;
  }

  /** @brief Finds a distance minimum after the first rising slope. */
  template <typename Shape>
  static Sample behind(const Shape &shape, const math::Vector &origin,
                       const math::Vector &direction, float radius) {
    Sample closest{origin, FLT_MAX};
    float previous = FLT_MAX;
    bool climbing = false;
    march(shape, origin, direction, radius, [&](Sample sample) {
      if (sample.distance > previous)
        climbing = true;
      if (climbing && sample.distance < previous &&
          sample.distance < closest.distance)
        closest = sample;
      previous = sample.distance;
      return true;
    });
    return closest;
  }
};

} // namespace hs_test::scan_tests
