/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file pose_table.h
 * @brief An orientation's poses at the column times of one sweep, sampled
 * along its sub-frame history, as a generator and as a per-frame table.
 */

#include "core/animation/orientation.h"
#include "core/math/3dmath.h"
#include "core/memory.h"
#include "core/platform/attributes.h"
#include "core/render/spacetime/surface.h"
#include <algorithm>
#include <array>
#include <cmath>

namespace Spacetime {

/**
 * @brief Visits a pose history at each first-arm column's sweep time.
 * @tparam Fn Callable `(int column, const math::Quaternion &pose)`.
 * @param history Unit poses evenly spaced over the frame, oldest first: index
 *        0 is the previous endpoint and `count - 1` the current one.
 * @param count History length, at least 1.
 * @param width Even canvas width.
 * @param fn Receives each column `c` in [0, width / 2) with the pose at
 *        `surface_time(c, width)` as math::Orientation::get(float) samples
 *        it, up to rounding; consecutive poses share a hemisphere.
 * @details Each history interval is a slerp, over the columns whose times
 * fall in it. Sines of evenly spaced angles obey
 * `sin((k + 1) d) = 2 cos(d) sin(k d) - sin((k - 1) d)`, so an interval costs a
 * few sines rather than two per column; the recurrence runs in double, as its
 * rounding grows with the column count.
 */
template <typename Fn>
HS_HOT_INLINE void for_each_time_pose(const math::Quaternion *history,
                                      int count, int width, Fn &&fn) {
  const int half = width / 2;
  const int segments = count - 1;
  if (segments == 0) {
    for (int c = 0; c < half; ++c)
      fn(c, history[0]);
    return;
  }
  int c = 0;
  for (int i = 0; i < segments; ++i) {
    // Columns c with i <= c * segments / half < i + 1, the last interval
    // closed at the frame end.
    const int end =
        i + 1 == segments ? half : ((i + 1) * half + segments - 1) / segments;
    if (end == c)
      continue;
    const math::SlerpArc arc(history[i], history[i + 1]);
    if (arc.theta == 0.0f) {
      for (; c < end; ++c) {
        const float local = surface_time(c, width) * segments - i;
        fn(c, (arc.from + local * (arc.to - arc.from)).normalized());
      }
      continue;
    }
    // Local fraction of column c within interval i.
    const double first = static_cast<double>(c) * segments / half - i;
    const double theta = arc.theta;
    const double step = theta * segments / half;
    const double start = first * theta;
    const double twice_cos = 2.0 * math::series_cos(step);
    const double inverse_sine = 1.0 / static_cast<double>(arc.sine);
    // rising = sin(s * theta) and falling = sin((1 - s) * theta), scaled.
    double rising = math::series_sin(start) * inverse_sine;
    double rising_next = math::series_sin(start + step) * inverse_sine;
    double falling =
        first == 0.0 ? 1.0 : math::series_sin(theta - start) * inverse_sine;
    double falling_next = math::series_sin(theta - start - step) * inverse_sine;
    for (; c < end; ++c) {
      fn(c, (static_cast<float>(falling) * arc.from +
             static_cast<float>(rising) * arc.to)
                .normalized());
      const double up = twice_cos * rising_next - rising;
      rising = rising_next;
      rising_next = up;
      const double down = twice_cos * falling_next - falling;
      falling = falling_next;
      falling_next = down;
    }
  }
}

/** @brief Rotation bounds of a pose history's path over its frame. */
struct HistoryTurn {
  /// Bound on the rotation between any pose of the path and its last pose:
  /// the summed rotation angles of the history's intervals, radians.
  float arc;
  /// Bound on the rotation between poses of the path a frame fraction `du`
  /// apart, per unit `du`: the largest interval angle times the interval
  /// count, radians. Equals `arc` when every interval turns alike.
  float rate;
};

/**
 * @brief Rotation bounds of a pose history's path.
 * @param history Unit poses evenly spaced over the frame, oldest first.
 * @param count History length, at least 1.
 * @return The path's total and per-fraction rotation bounds.
 */
HS_FLASH_INLINE inline HistoryTurn history_turn(const math::Quaternion *history,
                                                int count) {
  float arc = 0.0f, fastest = 0.0f;
  for (int i = 1; i < count; ++i) {
    const math::Quaternion step = history[i] * history[i - 1].conjugate();
    const float angle =
        2.0f * asinf(fminf(sqrtf(math::dot(step.v, step.v)), 1.0f));
    arc += angle;
    fastest = fmaxf(fastest, angle);
  }
  return {arc, fastest * static_cast<float>(count - 1)};
}

/**
 * @brief Copies an orientation's sub-frame history.
 * @tparam CAP Orientation capacity.
 * @param orientation Orientation whose history spans the frame.
 * @param history Receives its poses, oldest first; at least CAP entries.
 * @return The history length.
 */
template <int CAP>
HS_HOT_INLINE int copy_history(const math::Orientation<CAP> &orientation,
                               math::Quaternion *history) {
  const int count = orientation.length();
  for (int i = 0; i < count; ++i)
    history[i] = orientation.get(i);
  return count;
}

/**
 * @brief An orientation's rotation matrices at the column times of one sweep.
 * @details Entry `i` in [0, times()) holds the pose at the sweep time of
 * first-arm column `i`, which opposing columns share, and entry times() the
 * frame's end. Entry 0 is the history's first pose and entry times() its last,
 * exactly, so a history that starts where the previous frame's ended starts
 * its table where that frame's ended. Arena storage holds the entries of the
 * one canvas width bind() names, each as the first two rows of its matrix;
 * a rotation's third row is the cross product of its first two.
 */
class PoseTable {
public:
  /**
   * @brief Arena bytes bind() takes.
   * @param width Even canvas width.
   * @return The entries' size, alignment slack included.
   */
  static constexpr size_t storage_bytes(int width) {
    return static_cast<size_t>(width / 2 + 1) * 2 * sizeof(math::Vector) +
           alignof(math::Vector);
  }

  /**
   * @brief Allocates the entries for a canvas width.
   * @param arena Arena that outlives every frame the table spans.
   * @param width Even canvas width the table serves.
   */
  HS_COLD_MEMBER void bind(Arena &arena, int width) {
    half = width / 2;
    rows = arena.allocate_n<math::Vector>(2 * (half + 1));
    stationary_frame = true;
  }

  /**
   * @brief Fills the table from a pose history.
   * @param history Unit poses evenly spaced over the frame, oldest first.
   * @param count History length, at least 1.
   */
  HS_FLASH_MEMBER void build(const math::Quaternion *history, int count) {
    const int width = 2 * half;
    for_each_time_pose(history, count, width,
                       [&](int c, const math::Quaternion &pose)
                           HS_HOT_INLINE { store(c, pose); });
    store(0, history[0]);
    store(half, history[count - 1]);
    const math::Quaternion &current = history[count - 1];
    stationary_frame = true;
    for (int i = 0; i + 1 < count; ++i)
      stationary_frame = stationary_frame && history[i].r == current.r &&
                         history[i].v.x == current.v.x &&
                         history[i].v.y == current.v.y &&
                         history[i].v.z == current.v.z;
  }

  /**
   * @brief Fills the table from an orientation's sub-frame history.
   * @tparam CAP Orientation capacity.
   * @param orientation Orientation whose history spans the frame.
   */
  template <int CAP>
  HS_FLASH_MEMBER void build(const math::Orientation<CAP> &orientation) {
    std::array<math::Quaternion, CAP> history;
    build(history.data(), copy_history(orientation, history.data()));
  }

  /** @brief Whether every pose of the built history equals its last.
   * @return True for a frame without rotation. */
  bool stationary() const { return stationary_frame; }

  /** @brief Distinct column times of the bound width.
   * @return Half the canvas width. */
  int times() const { return half; }

  /**
   * @brief The pose at a column time.
   * @param index Time index in [0, times()]; times() is the frame end.
   * @return The rotation at fraction `index / times()`.
   */
  HS_HOT_INLINE math::RotationMatrix at(int index) const {
    math::RotationMatrix pose;
    pose.r0 = rows[2 * index];
    pose.r1 = rows[2 * index + 1];
    pose.r2 = math::cross(pose.r0, pose.r1);
    return pose;
  }

  /**
   * @brief The pose a column is lit at.
   * @param x Canvas column in [0, 2 * times()).
   * @return The rotation at the column's sweep time.
   */
  HS_HOT_INLINE math::RotationMatrix column(int x) const {
    return at(x < half ? x : x - half);
  }

  /**
   * @brief Rotates a point by the pose at a frame fraction.
   * @param v Point to rotate.
   * @param t Frame fraction, clamped to [0, 1].
   * @return @p v under the entry at @p t, or blended linearly between the two
   *         entries bracketing it.
   */
  HS_HOT_INLINE math::Vector rotate(const math::Vector &v, float t) const {
    const float s = hs::clamp(t, 0.0f, 1.0f) * static_cast<float>(half);
    const int i = std::min(static_cast<int>(s), half - 1);
    const float f = s - static_cast<float>(i);
    const math::Vector a = at(i).apply(v);
    if (f == 0.0f)
      return a;
    return a + (at(i + 1).apply(v) - a) * f;
  }

private:
  /**
   * @brief Stores an entry's first two rows.
   * @param index Time index.
   * @param pose Unit pose.
   */
  HS_HOT_INLINE void store(int index, const math::Quaternion &pose) {
    const math::RotationMatrix matrix(pose);
    rows[2 * index] = matrix.r0;
    rows[2 * index + 1] = matrix.r1;
  }

  math::Vector *rows = nullptr; ///< Two arena rows per entry, by time index.
  int half = 1;                 ///< Column times of the bound width.
  bool stationary_frame =
      true; ///< Every pose of the built history equals its last.
};

} // namespace Spacetime
