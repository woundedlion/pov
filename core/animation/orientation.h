/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file orientation.h
 * @brief Quaternion history for animation and motion blur.
 */

#include "math/3dmath.h"
#include <algorithm>
#include <array>

namespace math {
/**
 * @brief Class managing the current rotation state of an object, maintaining
 * history for interpolation.
 * @tparam CAP Maximum number of orientation frames retained in history.
 * @details Stores a list of Quaternions (`orientations`) generated during the
 * current frame step.
 */
template <int CAP = 4> class Orientation {
  static_assert(CAP >= 1, "Orientation requires CAP >= 1: the constructors "
                          "seed frame 0, so a zero-capacity history would "
                          "write past the storage array.");

public:
  static constexpr int CAPACITY = CAP;
  /**
   * @brief Default constructor (identity rotation).
   */
  Orientation() : num_frames(0) { set(Quaternion()); }

  /**
   * @brief Constructs with a specific initial quaternion.
   * @param q The initial quaternion.
   */
  explicit Orientation(const Quaternion &q) : num_frames(0) { set(q); }

  /**
   * @brief Gets the number of recorded orientation frames in the current step.
   * @return The number of frames.
   */
  int length() const { return num_frames; }

  /**
   * @brief Rotates a vector by the current (latest) orientation.
   * @param v The vector to orient.
   * @return The rotated vector.
   */
  Vector orient(const Vector &v) const {
    HS_CHECK(num_frames >= 1, "Orientation: no frames");
    return rotate(v, orientations[num_frames - 1]);
  }

  /**
   * @brief Rotates a vector by an orientation at a specific historical frame
   * index.
   * @param v The vector to orient.
   * @param i The frame index (0 being oldest, length-1 being current).
   * @return The rotated vector.
   */
  Vector orient(const Vector &v, int i) const {
    HS_CHECK(i >= 0 && i < num_frames, "Orientation: frame index out of range");
    return rotate(v, orientations[i]);
  }

  /**
   * @brief Rotates a vector backward by the inverse of the current (latest)
   * orientation.
   * @param v The vector to unorient.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v) const {
    HS_CHECK(num_frames >= 1, "Orientation: no frames");
    return rotate(v, orientations[num_frames - 1].conjugate());
  }

  /**
   * @brief Rotates a vector backward by the inverse of a specific historical
   * orientation.
   * @param v The vector to unorient.
   * @param i The frame index.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v, int i) const {
    HS_CHECK(i >= 0 && i < num_frames, "Orientation: frame index out of range");
    return rotate(v, orientations[i].conjugate());
  }

  /**
   * @brief Gets the current (latest) quaternion.
   * @return The Quaternion reference.
   */
  const Quaternion &get() const {
    HS_CHECK(num_frames >= 1, "Orientation: no frames");
    return orientations[num_frames - 1];
  }

  /**
   * @brief Gets the quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The Quaternion reference.
   */
  const Quaternion &get(int i) const {
    HS_CHECK(i >= 0 && i < num_frames, "Orientation: frame index out of range");
    return orientations[i];
  }

  /**
   * @brief Sets the orientation, clearing all history.
   * @param q The new orientation quaternion; MUST be unit length
   *   (HS_CHECK-trapped — a non-unit quaternion scales every rotated vector).
   * @return Reference to the Orientation object.
   */
  Orientation &set(const Quaternion &q) {
    HS_CHECK(std::abs(q.squared_magnitude() - 1.0f) < math::EPS_UNIT_QUAT_SQ,
             "Orientation: non-unit quaternion");
    orientations[0] = q;
    num_frames = 1;
    return *this;
  }

  /**
   * @brief Pushes a new quaternion onto the history, tracking a motion step.
   * @param q The new rotation quaternion; MUST be unit length
   *   (HS_CHECK-trapped — a non-unit quaternion scales every rotated vector).
   * @return Reference to the Orientation object.
   */
  Orientation &push(const Quaternion &q) {
    HS_CHECK(std::abs(q.squared_magnitude() - 1.0f) < math::EPS_UNIT_QUAT_SQ,
             "Orientation: non-unit quaternion");
    HS_CHECK(num_frames < CAPACITY, "Orientation: history full");
    orientations[num_frames++] = q;
    return *this;
  }

  /**
   * @brief Collapses the orientation history, retaining only the latest
   * quaternion.
   * @details Used after rendering a motion step to reset the motion blur
   * history.
   * @return Reference to the Orientation object.
   */
  Orientation &collapse() {
    if (num_frames > 1) {
      orientations[0] = orientations[num_frames - 1];
      num_frames = 1;
    }
    return *this;
  }

  /**
   * @brief Access a quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The const Quaternion reference.
   */
  const Quaternion &at(int i) const { return get(i); }

  /**
   * @brief Replaces a quaternion at a specific historical frame index.
   * @param i The frame index.
   * @param q Unit-length replacement quaternion.
   */
  void set_at(int i, const Quaternion &q) {
    HS_CHECK(i >= 0 && i < num_frames, "Orientation: frame index out of range");
    HS_CHECK(std::abs(q.squared_magnitude() - 1.0f) < math::EPS_UNIT_QUAT_SQ,
             "Orientation: non-unit quaternion");
    orientations[i] = q;
  }

  /**
   * @brief Normalizes and replaces a historical quaternion.
   * @param i The frame index.
   * @param q Replacement quaternion; must have nonzero magnitude.
   */
  void set_at_normalized(int i, Quaternion q) {
    HS_CHECK(i >= 0 && i < num_frames, "Orientation: frame index out of range");
    q.normalize();
    orientations[i] = q;
  }

  /**
   * @brief Increases the resolution of the history to 'count' steps by
   * Slerp-resampling the existing frames (uniform in source index, not arc
   * length).
   * @param count The target number of steps in the history.
   * @note When `count > CAPACITY` the trail is upsampled to `CAPACITY` instead:
   * graceful degradation (the write stays in-bounds, the current orientation
   * stays exact, only the motion-blur smear samples more coarsely), not a trap.
   * Raise `CAP` to trade RAM for a smoother fast smear.
   * @note A single-frame source (`num_frames == 1`, the common post-`set()`/
   * `collapse()` state) upsamples to a flat smear of that one frame; real motion
   * blur requires >=2 pushed frames.
   * @note Resampling is uniform in source index, so a trail whose frames were
   * pushed at a varying rate keeps that non-uniformity: the motion-blur streak
   * stays dense over the slow stretches and sparse over the fast ones rather
   * than spreading evenly along the arc.
   */
  void upsample(int count) {
    HS_CHECK(count >= 1, "Orientation: upsample count below 1");
    if (count > CAPACITY)
      count = CAPACITY;
    if (num_frames >= count)
      return;

    std::array<Quaternion, CAPACITY> old_orientations;
    std::copy(orientations.begin(), orientations.begin() + num_frames,
              old_orientations.begin());

    int old_num_frames = num_frames;

    for (int i = 0; i < count - 1; ++i) {
      float t = static_cast<float>(i) / (count - 1);
      float source_float_index = t * (old_num_frames - 1);
      int idx = static_cast<int>(source_float_index);
      float frac = source_float_index - idx;

      orientations[i] = slerp(
          old_orientations[idx],
          old_orientations[std::min((int)old_num_frames - 1, idx + 1)], frac);
    }
    // Endpoint maps exactly onto the last source frame (slerp(q, q, 0)).
    orientations[count - 1] = old_orientations[old_num_frames - 1].normalized();
    num_frames = count;
  }

private:
  std::array<Quaternion, CAPACITY>
      orientations; /**< Storage for historical quaternions. */
  int num_frames;   /**< The current number of active frames in history. */
};

} // namespace math
