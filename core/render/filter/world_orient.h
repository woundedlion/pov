/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/animation/orientation.h"
#include "core/render/filter/pipeline.h"
#include "core/animation/animation.h"
#include "core/math/geometry.h"
#include "core/color/color.h"
#include "core/render/temporal.h"

/**
 * @file world_orient.h
 * @brief Filter::World::Orient: rotates world points by a live Orientation at
 * the time a display arm passes them, and the rotated-edge cull helper.
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
forward_rotated_edge(const math::Vector &a, const math::Vector &b,
                     const math::Basis *pb, const math::Quaternion &q,
                     FwdFn &&forward) {
  if (pb) {
    math::Basis rb = math::rotate(*pb, q);
    return forward(math::rotate(a, q), math::rotate(b, q), &rb);
  }
  return forward(math::rotate(a, q), math::rotate(b, q), nullptr);
}

/**
 * @brief Rotates 3D points based on a dynamic Orientation.
 * @details The orientation's history spans one frame, from the previous end
 * pose to the current pose. With temporal sampling enabled, each point is
 * emitted once per arm pass, at the pose of the frame fraction `u` whose
 * sweep column contains it, with `age` offset by `1 - u`; see
 * `Temporal::capture`. Otherwise it is emitted once at the current pose.
 */
class Orient : public Is3D {
public:
  /**
   * @brief Binds the filter to a live orientation source.
   * @param orientation Orientation whose frame history drives the rotation.
   */
  Orient(const math::Orientation<> &orientation) : orientation(orientation) {}
  Orient(const math::Orientation<> &&) = delete;

  /**
   * @brief Rotates the point and emits it at each arm pass in the frame.
   * @param v World-space point to rotate.
   * @param color Source color, forwarded unchanged.
   * @param age Incoming age (frames); offset by `1 - u` at capture fraction u.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const math::Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    if (orientation.length() == 1 || !Temporal::temporal_sampling_enabled()) {
      pass(math::rotate(v, orientation.get()), color, age, alpha);
      return;
    }
    const Temporal::Captures captures = Temporal::capture(
        [&](float u) { return math::rotate(v, orientation.get(u)); });
    for (int i = 0; i < captures.count; ++i)
      pass(captures.position[i], color, age + (1.0f - captures.fraction[i]),
           alpha);
  }

  /**
   * @brief Re-emits a clip-cull edge under the rotation(s) applied at plot time.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if the edge under any history pose could intersect the clip
   *         band.
   * @details Visits every history pose, including the previous end pose, so
   *          the cull spans the whole interval plot() can capture in.
   */
  template <typename FwdFn>
  bool cull_edge(const math::Vector &a, const math::Vector &b,
                 const math::Basis *pb, FwdFn &&forward) const {
    bool hit = false;
    for (int i = 0; i < orientation.length() && !hit; ++i)
      hit = forward_rotated_edge(a, b, pb, orientation.get(i), forward);
    return hit;
  }

private:
  const math::Orientation<>
      &orientation; /**< Live orientation source driving the rotation. */
};

} // namespace World
} // namespace Filter
