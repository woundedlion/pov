/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/filter/pipeline.h"
#include "math/geometry.h"
#include "color/color.h"

/**
 * @file world_orient.h
 * @brief Filter::World::Orient: rotates world points by a live Orientation,
 * sweeping its intra-frame SLERP history. Carries the rotated-edge cull helper
 * the rigid World stages share.
 */

namespace Filter {
namespace World {

/**
 * @brief Forwards a clip-cull edge rotated by @p q, carrying the planar basis.
 * @tparam FwdFn Downstream cull continuation
 *         `bool(const Vector&, const Vector&, const Basis*)`.
 * @param a,b Edge endpoints in world space (pre-rotation).
 * @param pb Optional planar basis, rotated alongside the endpoints.
 * @param q Rotation the stage applies to this copy of the edge.
 * @param forward Tail-of-pipeline cull continuation.
 * @return forward()'s verdict for the rotated edge.
 */
template <typename FwdFn>
__attribute__((always_inline)) inline bool
forward_rotated_edge(const Vector &a, const Vector &b, const Basis *pb,
                     const Quaternion &q, FwdFn &&forward) {
  if (pb) {
    Basis rb = rotate(*pb, q);
    return forward(rotate(a, q), rotate(b, q), &rb);
  }
  return forward(rotate(a, q), rotate(b, q), nullptr);
}

/**
 * @brief Rotates 3D points based on a dynamic Orientation.
 * @details Sweeps the orientation's intra-frame SLERP history and offsets `age`
 * by the fractional `(1 - t)`, producing temporal motion blur.
 */
class Orient : public Is3D {
public:
  /**
   * @brief Binds the filter to a live orientation source.
   * @param orientation Orientation whose SLERP history drives the rotation.
   */
  Orient(Orientation<> &orientation) : orientation(orientation) {}

  /**
   * @brief Rotates and re-emits the point across the orientation's tween sweep.
   * @param v World-space point to rotate.
   * @param color Source color, forwarded unchanged.
   * @param age Incoming age (frames); offset by the fractional (1 - t) per tween step.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    tween(orientation, [&](const Quaternion &q, float t) {
      pass(rotate(v, q), color, age + (1.0f - t), alpha);
    });
  }

  /**
   * @brief Re-emits a clip-cull edge under the rotation(s) applied at plot time.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any tweened copy of the edge could intersect the clip band.
   * @details Mirrors plot()'s tween so the cull spans the same motion-blur sweep
   *          the renderer draws. Without it the rasterizer would cull by the
   *          un-rotated latitude and drop geometry an off-axis orientation moves
   *          into a segment band (docs/segmented_stateful_effects_spec.md).
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    bool hit = false;
    tween(orientation, [&](const Quaternion &q, float) {
      if (hit)
        return;
      hit = forward_rotated_edge(a, b, pb, q, forward);
    });
    return hit;
  }

private:
  Orientation<>
      &orientation; /**< Live orientation source driving the rotation. */
};

} // namespace World
} // namespace Filter
