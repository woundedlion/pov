/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "math/4dmath.h"
#include "render/ray/contract.h"

namespace Raycast {

enum class SamplingDomain { SPATIAL_3D, SLICE_4D };

/** @brief Camera embedding with three orthonormal ambient columns. */
struct PreparedCamera {
  SamplingDomain domain = SamplingDomain::SPATIAL_3D;
  math::Vec4 center{};
  math::Mat4 embedding = math::Mat4::identity();
  float radial_start = 0.0f;
  Interval interval;

  bool valid() const {
    if (domain != SamplingDomain::SPATIAL_3D &&
        domain != SamplingDomain::SLICE_4D)
      return false;
    if (!interval.valid() || !finite(radial_start) || radial_start < 0.0f)
      return false;
    for (int i = 0; i < 4; ++i) {
      if (!finite(center[i]))
        return false;
      for (int j = 0; j < 3; ++j)
        if (!finite(embedding.m[i][j]))
          return false;
    }
    if (domain == SamplingDomain::SPATIAL_3D &&
        (center[3] != 0.0f || embedding.m[3][0] != 0.0f ||
         embedding.m[3][1] != 0.0f || embedding.m[3][2] != 0.0f))
      return false;
    for (int i = 0; i < 3; ++i)
      for (int j = i; j < 3; ++j) {
        float dot = 0.0f;
        for (int k = 0; k < 4; ++k)
          dot += embedding.m[k][i] * embedding.m[k][j];
        if (fabsf(dot - (i == j ? 1.0f : 0.0f)) > 1e-4f)
          return false;
      }
    return true;
  }

  Ray ray(const math::Vector &direction) const {
    return {direction * radial_start, direction, interval};
  }

  math::Vec4 point4(const math::Vector &p) const {
    math::Vec4 result;
    for (int i = 0; i < 4; ++i)
      result[i] = center[i] + embedding.m[i][0] * p.x +
                  embedding.m[i][1] * p.y + embedding.m[i][2] * p.z;
    return result;
  }

  math::Vector point3(const math::Vector &p) const {
    const auto value = point4(p);
    return math::Vector(value[0], value[1], value[2]);
  }

  bool project_normal(const math::Vec4 &gradient, math::Vector &normal) const {
    float components[3]{};
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 4; ++i)
        components[j] += embedding.m[i][j] * gradient[i];
    const math::Vector projected(components[0], components[1], components[2]);
    const float LENGTH2 = math::dot(projected, projected);
    if (!finite(LENGTH2) || LENGTH2 <= 1e-20f) {
      normal = math::Vector();
      return false;
    }
    normal = projected * (1.0f / sqrtf(LENGTH2));
    return true;
  }
};

} // namespace Raycast
