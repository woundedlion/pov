/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file lattice_field.h
 * @brief Periodic lattice distance fields. */

#include <cmath>
#include <type_traits>
#include "math/4dmath.h"
#include "render/ray/contract.h"

namespace SDF {

/** @brief Thickened cubic or hypercubic edges in ambient world units. */
template <int Dimensions> struct WireLattice {
  static_assert(Dimensions == 3 || Dimensions == 4);
  /// World-space query point: math::Vector in 3D, math::Vec4 in 4D.
  using Point = std::conditional_t<Dimensions == 3, math::Vector, math::Vec4>;
  float cell_size = 1.0f;    ///< World-space edge length of one cell.
  float wire_radius = 0.05f; ///< Wire radius in world units.
  math::Vec4 origin{};       ///< World-space position of a lattice vertex.
  /// Lattice-to-world rotation; its leading block must be orthonormal.
  math::Mat4 rotation = math::Mat4::identity();

  /**
   * @brief Whether every parameter is finite and usable.
   * @return True for a positive cell size and wire radius and an orthonormal
   *   rotation block.
   */
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

  /**
   * @brief Maps a world point into lattice cells.
   * @param p World-space point.
   * @return Lattice coordinates; unused components are 0.
   */
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

  /**
   * @brief Lattice-space offset from the nearest edge to a point.
   * @param p World-space point.
   * @param axis Receives the axis the nearest edge runs along.
   * @return Offset in cells, zero on that axis.
   */
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

  /**
   * @brief Signed distance to the wire surface.
   * @param p World-space point.
   * @return World-space distance; negative inside a wire.
   */
  float distance(const Point &p) const {
    int axis;
    const math::Vec4 OFFSET = edge_offset(p, axis);
    float sum = 0.0f;
    for (int i = 0; i < Dimensions; ++i)
      sum += OFFSET[i] * OFFSET[i];
    return cell_size * sqrtf(sum) - wire_radius;
  }

  /**
   * @brief Unit gradient of distance().
   * @param p World-space point.
   * @return World-space unit vector away from the nearest edge; zero on it.
   */
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

  /**
   * @brief Surface normal; equal to gradient().
   * @param p World-space point.
   * @return World-space unit normal.
   */
  Point normal(const Point &p) const { return gradient(p); }

  /**
   * @brief Query guarantees of this field.
   * @return Exact exterior and interior clearance and exact surface tests.
   */
  Raycast::QueryCapabilities capabilities() const {
    return {true, true, true, 0.0f};
  }

  /**
   * @brief Field sample for ray queries.
   * @param p World-space point.
   * @return Signed distance, its magnitude as clearance, and the free axis
   *   as the feature.
   */
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
