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
#include "core/memory.h"
#include "core/render/spacetime.h"
#include "core/render/spacetime/pose_table.h"

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

/** @brief Forwards the end pose to a motion-padded edge predicate.
 * @param orientation Frame orientation history.
 * @param a,b Source edge endpoints.
 * @param pb Optional source planar basis.
 * @param forward Predicate padded by the history's sweep arc.
 * @return The endpoint edge's visibility verdict. */
template <typename FwdFn>
HS_HOT_INLINE inline bool
cull_oriented_edge(const math::Orientation<> &orientation,
                   const math::Vector &a, const math::Vector &b,
                   const math::Basis *pb, FwdFn &&forward) {
  return forward_rotated_edge(a, b, pb, orientation.get(), forward);
}

/** @brief Angular bound on a frame orientation path.
 * @param orientation Frame orientation history.
 * @return Sum of sub-frame angles, radians; zero on a flat surface. */
HS_HOT_FLASH_MEMBER inline float
orientation_sweep_arc(const math::Orientation<> &orientation) {
  if (Spacetime::flat_surface())
    return 0.0f;
  constexpr float ARC_ERROR_PAD = 1.02e-4f;
  float bound = 0.0f;
  for (int i = 1; i < orientation.length(); ++i) {
    const float cosine =
        fabsf(math::dot(orientation.get(i - 1), orientation.get(i)));
    if (cosine < 1.0f)
      bound += 2.0f * math::fast_acos(cosine) + ARC_ERROR_PAD;
  }
  return bound;
}

/**
 * @brief Rotates 3D points based on a dynamic Orientation.
 * @details Explicit-time draws read the frame's pose table. Ordinary plots
 * rotate once at the current endpoint and preserve the age register.
 */
class Orient : public Is3D {
public:
  static constexpr bool rigid_world_transform = true; ///< Preserves geodesics.
  /// Evaluates orientation history over a frame.
  static constexpr bool has_temporal_motion = true;

  /**
   * @brief Binds the filter to a live orientation source.
   * @param orientation Orientation whose frame history drives the rotation.
   */
  Orient(const math::Orientation<> &orientation) : orientation(orientation) {}
  Orient(const math::Orientation<> &&) = delete;

  /** @brief Angular bound on the orientation's frame path.
   * @return Sum of sub-frame rotation angles, in radians. */
  HS_HOT_FLASH_MEMBER float sweep_arc() const {
    return orientation_sweep_arc(orientation);
  }

  /** @brief Forwards an edge at its endpoint pose for a padded cull.
   * @param a First endpoint.
   * @param b Second endpoint.
   * @param pb Optional planar basis.
   * @param forward Downstream cull continuation.
   * @return The downstream visibility verdict. */
  template <typename FwdFn>
  bool cull_edge_timed(const math::Vector &a, const math::Vector &b,
                       const math::Basis *pb, FwdFn &&forward) const {
    return cull_oriented_edge(orientation, a, b, pb, forward);
  }

  /**
   * @brief Persistent bytes init_storage() reserves.
   * @param width Canvas width the poses serve.
   * @return The pose table's entries.
   */
  static constexpr size_t storage_bytes(int width) {
    return Spacetime::PoseTable::storage_bytes(width);
  }

  /**
   * @brief Allocates the pose table explicit-time rotation reads.
   * @param arena Persistent arena.
   * @param width Canvas width the poses serve.
   * @details Call from effect init(), before any draw through at_time().
   */
  HS_COLD_MEMBER void init_storage(Arena &arena, int width) {
    poses.bind(arena, width);
    poses_ready = false;
  }

  /** @brief Rotates a point at an explicit frame fraction.
   * @param cv Canvas identifying the frame.
   * @param v Unrotated point.
   * @param t Frame fraction in [0, 1].
   * @return The point under the frame's pose table at @p t, built once per
   *         frame from the orientation's sub-frame history; at @p t = 1,
   *         without rotation over the frame, or with temporal sampling
   *         disabled, under the current pose. */
  HS_O3_FN HS_NOINLINE_NOCLONE math::Vector
  at_time(const Canvas &cv, const math::Vector &v, float t) {
    if (t == 1.0f || Spacetime::flat_surface() || orientation.length() == 1) {
      if (orientation.length() == 1)
        poses_ready = false;
      return math::rotate(v, orientation.get());
    }
    pose_table(cv);
    if (poses.stationary())
      return math::rotate(v, orientation.get());
    return poses.rotate(v, t);
  }

  /** @brief The frame's poses by column time.
   * @param cv Canvas identifying the frame.
   * @return The pose table, built once per frame from the orientation's
   *         sub-frame history; a single-pose history rebuilds it on each
   *         call, as its pose can change within the frame. */
  HS_HOT_FLASH_MEMBER const Spacetime::PoseTable *pose_table(const Canvas &cv) {
    if (!poses_ready || pose_frame != cv.frame_number()) {
      poses.build(orientation);
      pose_frame = cv.frame_number();
      poses_ready = orientation.length() != 1;
    }
    return &poses;
  }

  /** @brief Plots under the explicit pose, preserving the age register.
   * @param cv Canvas identifying the frame.
   * @param v Unrotated point.
   * @param color Source color.
   * @param age Age register.
   * @param alpha Source opacity.
   * @param t Frame fraction.
   * @param pass Downstream world continuation. */
  template <typename PassFnT>
  HS_HOT_INLINE void plot_at(Canvas &cv, const math::Vector &v,
                             const ::Pixel &color, float age, float alpha,
                             float t, PassFnT &&pass) {
    pass(at_time(cv, v, t), color, age, alpha);
  }

  /**
   * @brief Rotates the point once at the current endpoint.
   * @param v World-space point to rotate.
   * @param color Source color, forwarded unchanged.
   * @param age Incoming age (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const math::Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    pass(math::rotate(v, orientation.get()), color, age, alpha);
  }

  /**
   * @brief Forwards the end-pose edge to a motion-padded predicate.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation padded by sweep_arc().
   * @return The endpoint edge's visibility verdict.
   */
  template <typename FwdFn>
  bool cull_edge(const math::Vector &a, const math::Vector &b,
                 const math::Basis *pb, FwdFn &&forward) const {
    return cull_oriented_edge(orientation, a, b, pb, forward);
  }

  /** @brief Visits historical poses for still-raster sample spacing.
   * @param a,b Source position and tangent.
   * @param pb Optional source basis.
   * @param forward Sampling-pose visitor.
   * @return True if a visitor accepts a pose. */
  template <typename FwdFn>
  HS_HOT_FLASH_MEMBER bool
  sampling_edge(const math::Vector &a, const math::Vector &b,
                const math::Basis *pb, FwdFn &&forward) const {
    for (int i = 0; i < orientation.length(); ++i)
      if (forward_rotated_edge(a, b, pb, orientation.get(i), forward))
        return true;
    return false;
  }

private:
  const math::Orientation<>
      &orientation; /**< Live orientation source driving the rotation. */
  Spacetime::PoseTable poses; ///< Per-frame poses over arena entries.
  uint32_t pose_frame = 0;    ///< Frame the pose table was built for.
  bool poses_ready = false;   ///< The pose table holds pose_frame's poses.
};

} // namespace World
} // namespace Filter
