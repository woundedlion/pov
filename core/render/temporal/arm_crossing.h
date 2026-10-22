/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file arm_crossing.h
 * @brief Closed-form arm passes of a moving point: where a longitude track
 * linear in time meets each sweeping arm, the analytic guard on that linear
 * model, and piecewise tracks.
 */

#include "core/math/3dmath.h"
#include "core/math/pixel_mapping.h"
#include "core/platform/attributes.h"
#include "core/render/temporal.h"
#include <cmath>
#include <utility>

namespace Temporal {

/** @brief One arm pass of a track. */
struct ArmPass {
  float time;  ///< Frame fraction of the pass.
  float along; ///< The pass's fraction of the track's span, in [0, 1].
  int arm;     ///< Branch: the pass lights unwrapped column `(arm + time) W/2`.
};

/**
 * @brief Arm passes of a column track linear over part of the frame.
 * @tparam W Canvas width.
 * @tparam Emit Callable `(const ArmPass &)`.
 * @param x0,x1 Unwrapped columns at fractions @p t0 and @p t1.
 * @param t0,t1 Frame fractions the track spans, `t0 < t1`.
 * @param frame_start Whether @p t0 starts the frame: a pass whose column
 *        lies within @p reach columns before the frame's first is emitted at
 *        @p t0.
 * @param frame_end Whether @p t1 ends the frame: a pass whose column lies
 *        within @p reach columns after the frame's last is emitted at @p t1.
 *        Otherwise a pass at @p t1 belongs to the following span.
 * @param emit Receives each pass: none, one, or two on either side of an arm
 *        cut.
 * @param reach Columns beyond the frame's ends whose passes are emitted at
 *        them.
 * @details Arm `n` lights unwrapped column `W / 2 * (n + t)` at fraction `t`,
 * so with `s = (t - t0) / (t1 - t0)` the pass is
 * `s = (x0 - (n + t0) W / 2) / ((t1 - t0) W / 2 - (x1 - x0))`, one division
 * per arm. A track keeping pace with the arm has no pass.
 */
template <int W, typename Emit>
HS_HOT_INLINE void arm_passes(float x0, float x1, float t0, float t1,
                              bool frame_start, bool frame_end, Emit &&emit,
                              float reach = 1.0f) {
  constexpr float HALF = W / 2.0f;
  const float dx = x1 - x0;
  const float span = t1 - t0;
  const float approach = span * HALF - dx;
  if (approach == 0.0f)
    return;
  const int first =
      static_cast<int>(floorf((fminf(x0, x1) - reach) / HALF + (1.0f - t1)));
  const int last =
      static_cast<int>(floorf((fmaxf(x0, x1) + reach) / HALF - t0));
  for (int branch = first; branch <= last; ++branch) {
    const float raw = (x0 - (branch + t0) * HALF) / approach;
    if (raw < 0.0f && (!frame_start || x0 < (branch + t0) * HALF - reach))
      continue;
    if ((raw >= 1.0f && !frame_end) ||
        (raw > 1.0f && x1 > (branch + t1) * HALF + reach))
      continue;
    const float s = hs::clamp(raw, 0.0f, 1.0f);
    emit(ArmPass{t0 + span * s, s, branch});
  }
}

/**
 * @brief Arm passes of a column track piecewise linear between knots evenly
 *        spaced over the frame.
 * @tparam W Canvas width.
 * @tparam Emit Callable `(const ArmPass &)`.
 * @param columns Unwrapped column at each knot `i / (knots - 1)`, each within
 *        a quarter turn of the one before.
 * @param knots Knot count, at least 2.
 * @param frame_start As arm_passes(), for the first piece.
 * @param frame_end As arm_passes(), for the last piece.
 * @param emit Receives each pass of each piece, its `along` relative to that
 *        piece; a pass at an inner knot is emitted once.
 */
template <int W, typename Emit>
HS_HOT_INLINE void piecewise_arm_passes(const float *columns, int knots,
                                        bool frame_start, bool frame_end,
                                        Emit &&emit) {
  const int pieces = knots - 1;
  for (int i = 0; i < pieces; ++i)
    arm_passes<W>(columns[i], columns[i + 1],
                  static_cast<float>(i) / static_cast<float>(pieces),
                  static_cast<float>(i + 1) / static_cast<float>(pieces),
                  frame_start && i == 0, frame_end && i + 1 == pieces, emit);
}

/** @brief One arm pass of a point moving linearly in screen space. */
struct LinearCapture {
  math::PixelCoords position; /**< Screen point at the pass. */
  float time;                 /**< Owning sweep fraction. */
};

/** @brief Arm passes of a screen trajectory linear over part of the frame.
 * @tparam W Canvas width.
 * @param a,b Projected positions at fractions @p t0 and @p t1, columns
 *        unwrapped against each other.
 * @param t0,t1 Frame fractions the trajectory spans, `t0 < t1`.
 * @param frame_start As arm_passes().
 * @param frame_end As arm_passes().
 * @param emit Receives each arm pass, its time and position on the
 *        trajectory. */
template <int W, typename Emit>
HS_HOT_INLINE void linear_captures(const math::PixelCoords &a,
                                   const math::PixelCoords &b, float t0,
                                   float t1, bool frame_start, bool frame_end,
                                   Emit &&emit) {
  const float dx = b.x - a.x, dy = b.y - a.y;
  arm_passes<W>(
      a.x, b.x, t0, t1, frame_start, frame_end,
      [&](const ArmPass &pass) HS_HOT_INLINE {
        emit(LinearCapture{{a.x + dx * pass.along, a.y + dy * pass.along},
                           pass.time});
      });
}

/** @brief Arm passes of a screen trajectory linear over the whole frame,
 *         including the passes whose splat reaches into the frame from
 *         outside it.
 * @tparam W Canvas width.
 * @param a,b Projected positions at t = 0 and t = 1, columns unwrapped
 *        against each other; `b.x - a.x` must stay below W / 2.
 * @param emit Receives each zero, one or two arm passes, its time clamped to
 *        [0, 1] and its position on the trajectory at that time. */
template <int W, typename Emit>
HS_HOT_INLINE void linear_captures(const math::PixelCoords &a,
                                   const math::PixelCoords &b, Emit &&emit) {
  linear_captures<W>(a, b, 0.0f, 1.0f, true, true, std::forward<Emit>(emit));
}

/** @brief A point's longitude over a span, modelled as linear in time. */
struct LongitudeTrack {
  float origin; ///< Column at the span's start, in [0, W).
  float travel; ///< Column change over the span, unwrapped.
  float rate;   ///< Pass columns per column of model error.
  float guard;  ///< Bound on a pass's column error under the model.
};

/**
 * @brief Models a point's longitude between two positions as linear in time.
 * @tparam W Canvas width.
 * @param from,to Positions at the span's start and end, unit length.
 * @param deviation Bound on the point's distance from the chord between
 *        them over the span.
 * @param span Frame fraction the positions span, in (0, 1].
 * @param track Receives the model.
 * @return False when the track can pass the polar axis, travels half the
 *         arm's sweep over the span or more, or its model errs by half a
 *         radian or more; the caller splits the span at a capture between or
 *         holds the point.
 * @details On the chord `p(s) = from + s d`, longitude has derivative
 * `c / |p|^2` with `c = from x d`, so its second derivative is at most
 * `2 |c| |d| / r^3`, `r` the chord's least distance from the axis, and the
 * linear model errs by an eighth of that. The point lies within the deviation
 * of the chord, which subtends at most `1.05 deviation / r` radians while the
 * deviation is at most `r / 2`. The guard adds a column for the column's own
 * width and the solved column's rounding, scaled by the rate at which the
 * sweep gains on the track.
 */
template <int W>
HS_HOT_INLINE bool model_track(const math::Vector &from, const math::Vector &to,
                               float deviation, float span,
                               LongitudeTrack &track) {
  const float sweep = span * (W / 2.0f);
  track.origin = math::vector_to_theta<W>(from);
  track.travel = centered_remainder(math::vector_to_theta<W>(to) - track.origin,
                                    static_cast<float>(W));
  const float cx = to.x - from.x, cz = to.z - from.z;
  const float chord2 = cx * cx + cz * cz;
  const float along =
      chord2 > 0.0f
          ? hs::clamp(-(from.x * cx + from.z * cz) / chord2, 0.0f, 1.0f)
          : 0.0f;
  const float nx = from.x + cx * along, nz = from.z + cz * along;
  const float near2 = nx * nx + nz * nz;
  if (!(4.0f * deviation * deviation <= near2 &&
        fabsf(track.travel) < 0.5f * sweep))
    return false;
  const float near = sqrtf(near2);
  const float radians =
      (0.25f * fabsf(from.x * cz - from.z * cx) * sqrtf(chord2) / near2 +
       1.05f * deviation) /
      near;
  if (!(radians < 0.5f))
    return false;
  track.rate = 1.0f / (1.0f - track.travel / sweep);
  track.guard =
      (radians * (1.0f / math::RADIANS_PER_COLUMN<W>)+1.0f) * track.rate;
  return true;
}

/**
 * @brief Arm passes of a modelled track.
 * @tparam W Canvas width.
 * @tparam Emit Callable `(const ArmPass &)`.
 * @param track Track over the span, as model_track() gives it.
 * @param t0,t1 Frame fractions the track spans.
 * @param frame_start As arm_passes().
 * @param frame_end As arm_passes().
 * @param emit Receives each pass; the point's true pass lies within
 *        `track.guard` columns of the emitted one.
 * @param reach As arm_passes(); `track.guard` also emits, at the frame's
 *        ends, every pass the model's error can carry into the frame.
 */
template <int W, typename Emit>
HS_HOT_INLINE void track_passes(const LongitudeTrack &track, float t0, float t1,
                                bool frame_start, bool frame_end, Emit &&emit,
                                float reach = 1.0f) {
  arm_passes<W>(track.origin, track.origin + track.travel, t0, t1, frame_start,
                frame_end, std::forward<Emit>(emit), reach);
}

} // namespace Temporal
