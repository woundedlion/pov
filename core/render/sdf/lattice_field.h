/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <type_traits>
#include "math/4dmath.h"
#include "render/ray/contract.h"

namespace SDF {

/** @brief Thickened cubic or hypercubic edges in ambient world units. */
template <int Dimensions> struct WireLattice {
  static_assert(Dimensions == 3 || Dimensions == 4);
  using Point = std::conditional_t<Dimensions == 3, math::Vector, math::Vec4>;
  float cell_size = 1.0f;
  float wire_radius = 0.05f;
  math::Vec4 origin{};
  math::Mat4 rotation = math::Mat4::identity();

  bool valid() const {
    if (!Raycast::finite(cell_size) || cell_size <= 0.0f ||
        !Raycast::finite(wire_radius) || wire_radius <= 0.0f)
      return false;
    for (int i = 0; i < Dimensions; ++i) {
      if (!Raycast::finite(origin[i]))
        return false;
      for (int j = 0; j < Dimensions; ++j) {
        float dot = 0.0f;
        for (int k = 0; k < Dimensions; ++k) {
          if (!Raycast::finite(rotation.m[k][i]))
            return false;
          dot += rotation.m[k][i] * rotation.m[k][j];
        }
        if (fabsf(dot - (i == j ? 1.0f : 0.0f)) > 1e-5f)
          return false;
      }
    }
    return true;
  }

  math::Vec4 local(const Point &p) const {
    math::Vec4 input;
    if constexpr (Dimensions == 3)
      input = {{p.x, p.y, p.z, 0.0f}};
    else
      input = p;
    math::Vec4 result;
    for (int i = 0; i < Dimensions; ++i)
      for (int j = 0; j < Dimensions; ++j)
        result[i] += rotation.m[j][i] * (input[j] - origin[j]) / cell_size;
    return result;
  }

  math::Vec4 edge_offset(const Point &p, int &axis) const {
    math::Vec4 offset = local(p);
    axis = 0;
    for (int i = 0; i < Dimensions; ++i) {
      offset[i] -= roundf(offset[i]);
      if (fabsf(offset[i]) > fabsf(offset[axis]))
        axis = i;
    }
    offset[axis] = 0.0f;
    return offset;
  }

  float distance(const Point &p) const {
    int axis;
    const math::Vec4 OFFSET = edge_offset(p, axis);
    float sum = 0.0f;
    for (int i = 0; i < Dimensions; ++i)
      sum += OFFSET[i] * OFFSET[i];
    return cell_size * sqrtf(sum) - wire_radius;
  }

  Point gradient(const Point &p) const {
    int axis;
    const math::Vec4 OFFSET = edge_offset(p, axis);
    float sum = 0.0f;
    for (int i = 0; i < Dimensions; ++i)
      sum += OFFSET[i] * OFFSET[i];
    const float SCALE = sum > 0.0f ? 1.0f / sqrtf(sum) : 0.0f;
    math::Vec4 result;
    for (int i = 0; i < Dimensions; ++i)
      for (int j = 0; j < Dimensions; ++j)
        result[i] += rotation.m[i][j] * OFFSET[j] * SCALE;
    if constexpr (Dimensions == 3)
      return {result[0], result[1], result[2]};
    else
      return result;
  }

  Point normal(const Point &p) const { return gradient(p); }

  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }

  Raycast::QuerySample sample(const Point &p) const {
    int axis;
    const math::Vec4 OFFSET = edge_offset(p, axis);
    float sum = 0.0f;
    for (int i = 0; i < Dimensions; ++i)
      sum += OFFSET[i] * OFFSET[i];
    const float VALUE = cell_size * sqrtf(sum) - wire_radius;
    return {VALUE, fabsf(VALUE), VALUE == 0.0f, 0, static_cast<uint32_t>(axis)};
  }
};

} // namespace SDF
