/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Shared fixtures for transformers tests. */

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include "core/math/lenses.h"
#include "core/animation/transformer.h"
#include "core/render/canvas.h"
#include "core/render/filter/feedback_style.h"
#include "core/math/easing.h"
#include "tests/test_fixture.h"
#include "tests/vec_test_util.h"
#include "tests/test_harness.h"

#include <cstdint>

namespace hs_test::transformers_tests {

/** @brief Manhattan distance between two transformed points.
 * @param a First point.
 * @param b Second point.
 * @return Sum of absolute coordinate differences. */
inline float displacement(const math::Vector &a, const math::Vector &b) {
  return std::abs(a.x - b.x) + std::abs(a.y - b.y) + std::abs(a.z - b.z);
}

/// Ripple pool used by the transformer contract cases.
using RipplePool =
    TransformerPool<Animation::RippleParams, Animation::Ripple, 8>;
static_assert(std::is_constructible_v<RipplePool, Timeline &>);
static_assert(!std::is_convertible_v<Timeline &, RipplePool>);
static_assert(std::is_constructible_v<RippleTransformer<8>, Timeline &>);
static_assert(!std::is_convertible_v<Timeline &, RippleTransformer<8>>);
static_assert(std::is_constructible_v<BallDropTransformer<8>, Timeline &>);
static_assert(!std::is_convertible_v<Timeline &, BallDropTransformer<8>>);

// ============================================================================
// OrientTransformer
// ============================================================================

static_assert(std::is_constructible_v<OrientTransformer<4>,
                                      const math::Orientation<4> &>);
static_assert(!std::is_constructible_v<OrientTransformer<4>,
                                       const math::Orientation<4> &&>);

// ============================================================================
// Ripple step endpoints and substitute params
// ============================================================================

/**
 * @brief Angle between two transformed unit directions.
 * @param a First direction.
 * @param b Second direction.
 * @return Angle in radians; NaN when either direction is nonfinite.
 */
inline double transformed_angle(const math::Vector &a, const math::Vector &b) {
  if (!finite_vec(a) || !finite_vec(b))
    return std::numeric_limits<double>::quiet_NaN();
  return chord_angle_between(a, b);
}

// ============================================================================
// Transformer<> composition order follows spawn order across slot recycling
// ============================================================================

/**
 * @brief Params carrying only a spawn-order tag, for the composition-order test.
 */
struct OrderParams {
  int order = 0; /**< Composition tag written by TagAnim at spawn. */

  /** @brief No live-config hook for TransformerPool::prepare_frame(). */
  static constexpr bool NEEDS_REFRESH_FROM = false;
  /** @brief No derived state, so no prepare_frame() sync() hook. */
  static constexpr bool NEEDS_SYNC = false;
};

/**
 * @brief Shifts a base-10 digit of the tag into x on each application, so the
 *        final x reads back the exact order the warps composed in.
 * @param v Input point.
 * @param p Spawn-order tag.
 * @return Point carrying the appended tag in x.
 */
inline math::Vector order_transform(const math::Vector &v,
                                    const OrderParams &p) {
  return math::Vector(v.x * 10.0f + static_cast<float>(p.order), v.y, v.z);
}

/**
 * @brief Minimal animation that stamps its spawn-order tag into the params.
 * @details Finite and non-repeating so spawn() accepts it and a short duration
 *          frees its pool slot on completion.
 */
struct TagAnim : public Animation::AnimationBase<TagAnim> {
  /** @brief Stores an entity's composition-order tag.
   * @param params Destination parameters.
   * @param order Composition tag.
   * @param duration Animation duration in frames. */
  TagAnim(OrderParams &params, int order, int duration)
      : Animation::AnimationBase<TagAnim>(duration, false) {
    params.order = order;
  }
};

// ============================================================================
// FieldTransformer<> — scalar fields sum; bounds sum; slots recycle
// ============================================================================

/**
 * @brief Params carrying one scalar contribution, for the field tests.
 */
struct FieldTestParams {
  float value = 0.0f; /**< Contribution written by FieldTagAnim at spawn. */

  /** @brief No live-config hook for TransformerPool::prepare_frame(). */
  static constexpr bool NEEDS_REFRESH_FROM = false;
  /** @brief No derived state, so no prepare_frame() sync() hook. */
  static constexpr bool NEEDS_SYNC = false;

  /**
   * @brief Upper bound on |field_test_field| for this entity.
   * @return |value| (the field is value * v.x with |v.x| <= 1).
   */
  float field_bound() const { return std::fabs(value); }
};

/**
 * @brief Scales the sample's x by the entity's value, so the sum is readable
 *        and the point argument is proven to reach the field function.
 * @param v Sample direction.
 * @param p Scalar contribution parameters.
 * @return Product of x and the entity's value.
 */
inline float field_test_field(const math::Vector &v, const FieldTestParams &p) {
  return p.value * v.x;
}

/**
 * @brief Minimal animation that stamps its value into the params.
 * @details Finite and non-repeating so spawn() accepts it and a short duration
 *          frees its pool slot on completion.
 */
struct FieldTagAnim : public Animation::AnimationBase<FieldTagAnim> {
  /** @brief Stores an entity's scalar field contribution.
   * @param params Destination parameters.
   * @param value Scalar contribution.
   * @param duration Animation duration in frames. */
  FieldTagAnim(FieldTestParams &params, float value, int duration)
      : Animation::AnimationBase<FieldTagAnim>(duration, false) {
    params.value = value;
  }
};

/// Two-entity scalar field transformer used by the field cases.
using TestFieldTransformer =
    FieldTransformer<FieldTestParams, FieldTagAnim, field_test_field, 2>;

} // namespace hs_test::transformers_tests
