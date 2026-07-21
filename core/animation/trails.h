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
 * @brief Fixed-capacity history of snapshots over a StaticCircularBuffer.
 * @tparam T The snapshot type stored in the trail.
 * @tparam CAP The maximum number of snapshots to keep.
 */
template <typename T, int CAP> class Trail {
public:
  static constexpr int CAPACITY = CAP; /**< Max retained snapshots. */

  /**
   * @brief Records a snapshot.
   * @param source The value to copy into the trail.
   */
  void record(const T &source) { snapshots.push_back(source); }

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
  const T &get(size_t i) const { return snapshots[i]; }

  /**
   * @brief Gets a specific snapshot (mutable). 0 is oldest, length()-1 newest.
   * @param i Index into the history (0 oldest, length()-1 newest).
   * @return Mutable reference to the requested snapshot.
   */
  T &get(size_t i) { return snapshots[i]; }

  /**
   * @brief Clears the history.
   */
  void clear() { snapshots.clear(); }

  /**
   * @brief Removes the oldest snapshot.
   */
  void expire() { snapshots.pop_front(); }

private:
  StaticCircularBuffer<T, CAP> snapshots;
};

/** @brief History of Orientation states. */
template <typename OrientationType, int CAPACITY>
using OrientationTrail = Trail<OrientationType, CAPACITY>;

/** @brief History of world-space Vector positions. */
template <int CAP> using VectorTrail = Trail<Vector, CAP>;

/**
 * @brief Fixed-capacity history of unit-sphere positions, quantized to snorm16.
 * @tparam CAP The maximum number of snapshots to keep.
 * @details Stores 3x int16 per snapshot (6 bytes vs Vector's 12). Components
 * are clamped to [-1, 1] on record and round-trip with error <= 1/65534 per
 * component, so get() decodes and returns by value.
 */
template <int CAP> class QuantizedVectorTrail {
public:
  static constexpr int CAPACITY = CAP; /**< Max retained snapshots. */

  /**
   * @brief Records a snapshot, quantizing each component.
   * @param source The position to copy into the trail; components outside
   *        [-1, 1] are clamped.
   */
  void record(const Vector &source) {
    snapshots.push_back(
        {quantize(source.x), quantize(source.y), quantize(source.z)});
  }

  /**
   * @brief Gets the number of recorded snapshots.
   * @return Count of live snapshots in the trail.
   */
  size_t length() const { return snapshots.size(); }

  /**
   * @brief Decodes a specific snapshot.
   * @param i Index into the history: 0 is the OLDEST snapshot, length()-1 the
   *          newest (same ordering as Trail).
   * @return The decoded position, by value.
   */
  Vector get(size_t i) const { return decode(snapshots[i]); }

  /**
   * @brief Visits decoded snapshots with normalized oldest-to-newest progress.
   * @param callback Invoked as `void(const Vector &, float t)`.
   */
  void __attribute__((noinline)) tween(VectorTweenFn callback) const {
    const size_t len = snapshots.size();
    if (len == 0)
      return;
    if (len == 1) {
      callback(decode(snapshots.front()), 1.0f);
      return;
    }
    const float denominator = static_cast<float>(len - 1);
    snapshots.for_each([&](const Snorm3 &s, uint32_t i) {
      callback(decode(s), static_cast<float>(i) / denominator);
    });
  }

  /**
   * @brief Clears the history.
   */
  void clear() { snapshots.clear(); }

  /**
   * @brief Removes the oldest snapshot.
   */
  void expire() { snapshots.pop_front(); }

private:
  struct Snorm3 {
    int16_t x, y, z;
  };

  static constexpr float SCALE = 32767.0f;
  static constexpr float INV_SCALE = 1.0f / SCALE;

  static int16_t quantize(float c) {
    return static_cast<int16_t>(roundf(hs::clamp(c, -1.0f, 1.0f) * SCALE));
  }

  static Vector decode(const Snorm3 &s) {
    return Vector(s.x * INV_SCALE, s.y * INV_SCALE, s.z * INV_SCALE);
  }

  StaticCircularBuffer<Snorm3, CAP> snapshots;
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
 * @brief Helper to iterate over a QuantizedVectorTrail's historical frames.
 * @tparam CAPACITY Trail capacity.
 * @param trail The trail to iterate.
 * @param callback The function to call for each frame: `void(const Vector&,
 * float t)`.
 */
template <int CAPACITY>
void tween(const Animation::QuantizedVectorTrail<CAPACITY> &trail,
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
 * @brief deep_tween grouped by trail frame: one callback per contributing
 * frame with that frame's sub-positions and global t values.
 * @param trail The OrientationTrail to iterate.
 * @param callback Invoked per contributing frame as `void(const Quaternion
 * *qs, const float *ts, int count)`. Concatenated over frames, (qs, ts) is
 * exactly deep_tween's emission order. qs points into the frame's own
 * storage; ts into a buffer reused across callbacks — neither outlives the
 * call.
 */
template <typename FrameFn>
void deep_tween_frames(const Tweenable auto &trail, FrameFn &&callback) {
  using FrameT = std::remove_cvref_t<decltype(trail.get(0))>;
  float ts[FrameT::CAPACITY];

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
    ts[0] = 1.0f;
    callback(&trail.get(0).get(0), ts, 1);
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

    int count = 0;
    for (size_t j = start_j; j < frame_size; ++j) {
      float sub_t =
          (frame_size > 1) ? static_cast<float>(j) / (frame_size - 1) : 0.0f;
      ts[count++] = (static_cast<float>(active_idx) + sub_t) / span;
    }
    // Sub-positions are contiguous in the frame's storage, so one pointer
    // spans [start_j, frame_size).
    callback(&frame.get(static_cast<int>(start_j)), ts, count);
    ++active_idx;
  }
}

/**
 * @brief Helper to iterate over an Animation::OrientationTrail, flattening its
 * per-frame sub-positions. A bare Orientation has no per-frame structure to
 * flatten — use tween() for that.
 * @param trail The OrientationTrail to iterate.
 * @param callback The function to call for each step: `void(const T&, float
 * t)`.
 */
void deep_tween(const Tweenable auto &trail, TweenFn callback) {
  deep_tween_frames(trail,
                    [&](const Quaternion *qs, const float *ts, int count) {
                      for (int i = 0; i < count; ++i)
                        callback(qs[i], ts[i]);
                    });
}
