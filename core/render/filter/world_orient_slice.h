/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <span>
#include <algorithm>
#include "render/filter/world_orient.h"
#include "math/geometry.h"
#include "color/color.h"

/**
 * @file world_orient_slice.h
 * @brief Filter::World::OrientSlice: picks an orientation from a borrowed list
 * by the point's projection onto a slicing axis.
 */

namespace Filter {
namespace World {

/**
 * @brief Selects an orientation from a list based on the point's projection
 * onto an axis. Useful for slicing objects with different rotations.
 */
class OrientSlice : public Is3D {
public:
  static constexpr bool requires_unit_world_input = true;
  /**
   * @brief Binds the slice selector to an orientation list and a slicing axis.
   * @param orientations Candidate orientations, indexed by axis projection.
   *        Borrowed, not copied: the backing array must outlive the filter.
   * @param axis Unit axis the point is projected onto to pick an orientation.
   */
  OrientSlice(std::span<const Orientation<>> orientations, const Vector &axis)
      : enabled(true), axis(axis.normalized()), orientations(orientations) {}

  /**
   * @brief Selects an orientation by axis projection and re-emits the point.
   * @param v World-space point to rotate.
   * @param color Source color, forwarded unchanged.
   * @param age Incoming age (frames); offset by fractional (1 - t) per tween step.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   * @details Passes through untouched when disabled or the orientation list is empty.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    if (!enabled || orientations.empty()) {
      pass(v, color, age, alpha);
      return;
    }

    float projection = v.x * axis.x + v.y * axis.y + v.z * axis.z;
    float dot_val = std::max(-1.0f, std::min(1.0f, projection));
    float t = hs::clamp(1.0f - fast_acos(dot_val) / PI_F, 0.0f, 1.0f);

    size_t count = orientations.size();
    size_t idx = static_cast<size_t>(floorf(t * count));
    if (idx >= count)
      idx = count - 1;

    const Orientation<> &q = orientations[idx];
    tween(q, [&](const Quaternion &rot, float tween_t) {
      pass(rotate(v, rot), color, age + (1.0f - tween_t), alpha);
    });
  }

  /**
   * @brief Re-emits a clip-cull edge under every candidate slice's rotation.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any candidate slice's tweened copy could intersect the band.
   * @details The endpoints may fall in different slices, so bound conservatively
   *          over all candidates rather than replicating the per-point selector.
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    if (!enabled || orientations.empty())
      return forward(a, b, pb);
    for (const Orientation<> &o : orientations) {
      bool hit = false;
      tween(o, [&](const Quaternion &q, float) {
        if (hit)
          return;
        hit = forward_rotated_edge(a, b, pb, q, forward);
      });
      if (hit)
        return true;
    }
    return false;
  }

  /**
   * @brief Sets the slicing axis, renormalizing to enforce the unit-length
   * contract that the projection bucket math assumes.
   * @param a New slicing axis (any non-zero length; renormalized internally).
   */
  void set_axis(const Vector &a) { axis = a.normalized(); }

  /**
   * @brief Accesses the current (unit-length) slicing axis.
   * @return The unit axis points are projected onto to select a slice.
   */
  const Vector &get_axis() const { return axis; }

  /**
   * @brief Enables or disables the slice selection.
   * @param value When false, points and edges pass through unrotated.
   */
  void set_enabled(bool value) { enabled = value; }

  /**
   * @brief Reports whether slice selection is active.
   * @return False when the filter passes points through unrotated.
   */
  bool is_enabled() const { return enabled; }

private:
  bool enabled; /**< When false, the filter passes points through unrotated. */
  Vector axis;  /**< Unit axis points are projected onto to select a slice. */
  std::span<const Orientation<>>
      orientations; /**< Candidate orientations indexed by projection; borrowed,
                         so the backing array must outlive the filter. */
};

} // namespace World
} // namespace Filter
