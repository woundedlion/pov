/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

namespace Animation {

/**
 * @brief Manages a history of Orientation states.
 * @tparam OrientationType The orientation snapshot type stored in the trail.
 * @tparam CAPACITY The maximum number of snapshots to keep.
 */
template <typename OrientationType, int CAPACITY> class OrientationTrail {
public:
  /**
   * @brief Records a snapshot of the current orientation state.
   * @param source The orientation to copy.
   */
  void record(const OrientationType &source) { snapshots.push_back(source); }

  /**
   * @brief Gets the number of recorded snapshots.
   * @return Count of live snapshots in the trail.
   */
  size_t length() const { return snapshots.size(); }

  /**
   * @brief Gets a specific snapshot.
   * @param i Index into the history: 0 is the OLDEST snapshot, length()-1 the
   *          newest. (record() appends the newest at the end of the underlying
   *          ring buffer, whose operator[](0) is the oldest live element.)
   *          Matches the JS simulator's trail ordering.
   * @return Const reference to the requested snapshot.
   */
  const OrientationType &get(size_t i) const { return snapshots[i]; }

  /**
   * @brief Gets a specific snapshot (mutable). 0 is oldest, length()-1 newest.
   * @param i Index into the history (0 oldest, length()-1 newest).
   * @return Mutable reference to the requested snapshot.
   */
  OrientationType &get(size_t i) { return snapshots[i]; }

  /**
   * @brief Clears the history.
   */
  void clear() { snapshots.clear(); }

  /**
   * @brief Removes the oldest snapshot.
   */
  void expire() { snapshots.pop(); }

private:
  StaticCircularBuffer<OrientationType, CAPACITY> snapshots;
};

/**
 * @brief Manages a history of world-space Vector positions.
 * @tparam CAP The maximum number of snapshots to keep.
 */
template <int CAP> class VectorTrail {
public:
  static constexpr int CAPACITY = CAP; /**< Max retained snapshots. */

  /**
   * @brief Records a world-space position snapshot.
   * @param source The position to copy into the trail.
   */
  void record(const Vector &source) { snapshots.push_back(source); }

  /**
   * @brief Gets the number of recorded snapshots.
   * @return Count of live positions in the trail.
   */
  size_t length() const { return snapshots.size(); }

  /**
   * @brief Gets a specific position snapshot.
   * @param i Index into the history (0 oldest, length()-1 newest).
   * @return Const reference to the requested position.
   */
  const Vector &get(size_t i) const { return snapshots[i]; }

  /**
   * @brief Gets a specific position snapshot (mutable).
   * @param i Index into the history (0 oldest, length()-1 newest).
   * @return Mutable reference to the requested position.
   */
  Vector &get(size_t i) { return snapshots[i]; }

  /**
   * @brief Clears the history.
   */
  void clear() { snapshots.clear(); }

  /**
   * @brief Removes the oldest snapshot.
   */
  void expire() { snapshots.pop(); }

private:
  StaticCircularBuffer<Vector, CAP> snapshots;
};

} // namespace Animation

/**
 * @brief Helper to iterate over an Orientation's historical frames.
 * @tparam CAP Orientation sub-frame capacity.
 * @param o The orientation to iterate.
 * @param callback The function to call for each frame: `void(const Quaternion&,
 * float t)`.
 */
template <int CAP>
void tween(const Orientation<CAP> &o, TweenFn callback) {
  int len = o.length();
  int start = (len > 1) ? 1 : 0;
  for (int i = start; i < len; ++i) {
    // A lone snapshot is the newest sub-position → t = 1 (age-neutral).
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 1.0f;
    callback(o.get(i), t);
  }
}

/**
 * @brief Helper to iterate over a VectorTrail's historical frames.
 * @tparam CAPACITY Trail capacity.
 * @param trail The trail to iterate.
 * @param callback The function to call for each frame: `void(const Vector&,
 * float t)`.
 */
template <int CAPACITY>
void tween(const Animation::VectorTrail<CAPACITY> &trail,
           VectorTweenFn callback) {
  size_t len = trail.length();
  if (len == 0)
    return;

  for (size_t i = 0; i < len; ++i) {
    // A lone sample is the newest (head) position → t = 1 (age-neutral).
    float t = (len > 1) ? static_cast<float>(i) / (len - 1) : 1.0f;
    callback(trail.get(i), t);
  }
}

/**
 * @brief Helper to iterate over any Tweenable object (Orientation or
 * Animation::OrientationTrail).
 * @param trail The Tweenable object to iterate.
 * @param callback The function to call for each step: `void(const T&, float
 * t)`.
 */
void deep_tween(const Tweenable auto &trail, TweenFn callback) {
  size_t trail_len = trail.length();
  if (trail_len == 0)
    return;

  // The age ramp must end at 1.0 on the newest *plotted* orientation. A motionless
  // tail frame collapses to a single sub-frame that emits nothing, so normalizing
  // by trail_len would strand the ramp below 1.0; normalize against the newest
  // frame that actually contributes instead.
  size_t last = trail_len - 1;
  while (last > 0 && trail.get(last).length() <= 1)
    --last;

  // A fully motionless trail collapses to frame 0's lone sub-position, which must
  // read t = 1.0 (age-neutral) or quintic_kernel(0) renders a static head invisible.
  if (last == 0 && trail.get(0).length() == 1) {
    const auto &frame = trail.get(0);
    callback(frame.get(0), 1.0f);
    return;
  }

  // An interior length-1 frame holds only the boundary it shares with the prior
  // frame, so it emits nothing; counting it in the span would leave a hole in the
  // age ramp. Normalize against, and index by, the frames that actually contribute.
  size_t num_active = 0;
  for (size_t i = 0; i <= last; ++i)
    if (i == 0 || trail.get(i).length() > 1)
      ++num_active;
  float span = static_cast<float>(num_active);

  size_t active_idx = 0;
  for (size_t i = 0; i <= last; ++i) {
    const auto &frame = trail.get(i);
    size_t frame_size = frame.length();
    if (i != 0 && frame_size <= 1)
      continue;
    size_t start_j = (i == 0) ? 0 : 1;

    for (size_t j = start_j; j < frame_size; ++j) {
      const auto &q = frame.get(j);
      float sub_t =
          (frame_size > 1) ? static_cast<float>(j) / (frame_size - 1) : 0.0f;
      float global_t = (static_cast<float>(active_idx) + sub_t) / span;
      callback(q, global_t);
    }
    ++active_idx;
  }
}
