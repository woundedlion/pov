/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <bit>

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
  HS_HOT_FLASH_MEMBER float cull_bound() const {
    if (!Temporal::temporal_sampling_enabled())
      return 0.0f;
    float bound = 0.0f;
    for (int i = 1; i < orientation.length(); ++i) {
      const auto &a = orientation.get(i - 1);
      const auto &b = orientation.get(i);
      const float cosine = fabsf(a.r * b.r + math::dot(a.v, b.v));
      if (cosine < 1.0f)
        bound += 2.0f * math::fast_acos(cosine) + ARC_ERROR_PAD;
    }
    return bound;
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
    return forward_rotated_edge(a, b, pb, orientation.get(), forward);
  }

  /** @brief Rotates a point at an explicit frame fraction.
   * @param cv Canvas identifying the frame.
   * @param v Unrotated point.
   * @param t Frame fraction in [0, 1].
   * @return The point under the sampled pose. */
  HS_O3_FN HS_NOINLINE_NOCLONE math::Vector
  at_time(const Canvas &cv, const math::Vector &v, float t) {
    if (t == 1.0f || !Temporal::temporal_sampling_enabled() ||
        orientation.length() == 1) {
      if (orientation.length() == 1)
        poses_ready = false;
      return math::rotate(v, orientation.get());
    }
    if (!poses_ready || pose_frame != cv.frame_number())
      prepare_poses(cv.frame_number());
    if (!sample_ready || sample_time != t) {
      const float s = hs::clamp(t, 0.0f, 1.0f) * (POSE_COUNT - 1);
      const int i = std::min(static_cast<int>(s), POSE_COUNT - 2);
      const float f = s - i;
      sample_pose =
          math::Quaternion(poses[i][0] + (poses[i + 1][0] - poses[i][0]) * f,
                           poses[i][1] + (poses[i + 1][1] - poses[i][1]) * f,
                           poses[i][2] + (poses[i + 1][2] - poses[i][2]) * f,
                           poses[i][3] + (poses[i + 1][3] - poses[i][3]) * f);
      sample_pose *= 1.0f / sqrtf(sample_pose.squared_magnitude());
      sample_time = t;
      sample_ready = true;
      rotation_ready = false;
    }
    if (!rotation_ready ||
        std::bit_cast<uint32_t>(sample_input.x) !=
            std::bit_cast<uint32_t>(v.x) ||
        std::bit_cast<uint32_t>(sample_input.y) !=
            std::bit_cast<uint32_t>(v.y) ||
        std::bit_cast<uint32_t>(sample_input.z) !=
            std::bit_cast<uint32_t>(v.z)) {
      sample_input = v;
      sample_output = math::rotate(v, sample_pose);
      rotation_ready = true;
    }
    return sample_output;
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
    const Temporal::Captures captures = Temporal::capture(orientation, v);
    for (int i = 0; i < captures.count; ++i) {
      const float u = captures.fraction[i];
      pass(math::rotate(v, orientation.get(u)), color, age + (1.0f - u), alpha);
    }
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
  /** @brief Radian padding for twice the fast-acos approximation error. */
  static constexpr float ARC_ERROR_PAD = 1.02e-4f;
  const math::Orientation<>
      &orientation; /**< Live orientation source driving the rotation. */
  /** @brief Uniform pose lattice size. */
  static constexpr int POSE_COUNT = 16;
  float poses[POSE_COUNT][4]; /**< Sign-aligned quaternion components. */
  uint32_t pose_frame = 0;    /**< Frame whose sub-frame history is cached. */
  bool poses_ready = false;   /**< The pose lattice has been populated. */

  math::Quaternion sample_pose; /**< Normalized pose at the cached fraction. */
  math::Vector sample_input;   /**< Last point rotated under the cached pose. */
  math::Vector sample_output;  /**< Rotation of the cached input point. */
  float sample_time = 0.0f;    /**< Fraction of the cached normalized pose. */
  bool sample_ready = false;   /**< The fraction cache is populated. */
  bool rotation_ready = false; /**< The point cache is populated. */

  /** @brief Builds the pose lattice from the orientation history.
   * @param frame Frame identifier to cache. */
  HS_HOT_FLASH_MEMBER void prepare_poses(uint32_t frame) {
    for (int i = 0; i < POSE_COUNT; ++i) {
      math::Quaternion q =
          orientation.get(static_cast<float>(i) / (POSE_COUNT - 1));
      if (i > 0 && poses[i - 1][0] * q.r + poses[i - 1][1] * q.v.x +
                           poses[i - 1][2] * q.v.y + poses[i - 1][3] * q.v.z <
                       0.0f)
        q *= -1.0f;
      poses[i][0] = q.r;
      poses[i][1] = q.v.x;
      poses[i][2] = q.v.y;
      poses[i][3] = q.v.z;
    }
    pose_frame = frame;
    poses_ready = true;
    sample_ready = false;
  }
};

} // namespace World
} // namespace Filter
