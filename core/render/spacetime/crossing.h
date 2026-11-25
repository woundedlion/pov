/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file crossing.h
 * @brief Closed-form arm passes of a moving point: where a longitude track
 * linear in time meets each sweeping arm, the analytic guard on that linear
 * model, and piecewise tracks.
 */

#include "core/math/3dmath.h"
#include "core/math/pixel_mapping.h"
#include "core/platform/attributes.h"
#include "core/platform/platform.h"
#include <cmath>
#include <utility>

namespace Spacetime {

/** @brief One arm pass of a track. */
struct Crossing {
  float time;  ///< Frame fraction of the pass.
  float along; ///< The pass's fraction of the track's span, in [0, 1].
  int arm;     ///< Branch: the pass lights unwrapped column `(arm + time) W/2`.
};

/** @brief Span ends at which crossings() emits the passes just beyond the
 *         span's columns. */
enum class FrameEnds : uint8_t {
  NEITHER = 0, ///< Neither end.
  START = 1,   ///< The span's start.
  END = 2,     ///< The span's end.
  BOTH = 3,    ///< Both ends.
};

/** @brief Whether ends include the span's start.
 * @param ends Span ends.
 * @return True for START and BOTH. */
constexpr bool includes_start(FrameEnds ends) { return int(ends) & 1; }

/** @brief Whether ends include the span's end.
 * @param ends Span ends.
 * @return True for END and BOTH. */
constexpr bool includes_end(FrameEnds ends) { return int(ends) & 2; }

/**
 * @brief Arm passes of a column track linear over part of the frame, on a
 *        canvas of a width known at run time.
 * @tparam Emit Callable `(const Crossing &)`.
 * @param width Canvas width.
 * @param x0,x1 Unwrapped columns at fractions @p t0 and @p t1.
 * @param t0,t1 Frame fractions the track spans, `t0 < t1`.
 * @param ends As crossings().
 * @param emit As crossings().
 * @param reach As crossings().
 */
template <typename Emit>
HS_HOT_INLINE void crossings(int width, float x0, float x1, float t0, float t1,
                             FrameEnds ends, Emit &&emit, float reach = 1.0f) {
  const float half = width / 2.0f;
  const float dx = x1 - x0;
  const float span = t1 - t0;
  const float approach = span * half - dx;
  if (approach == 0.0f)
    return;
  const int first =
      static_cast<int>(floorf((fminf(x0, x1) - reach) / half + (1.0f - t1)));
  const int last =
      static_cast<int>(floorf((fmaxf(x0, x1) + reach) / half - t0));
  for (int branch = first; branch <= last; ++branch) {
    const float raw = (x0 - (branch + t0) * half) / approach;
    if (raw < 0.0f &&
        (!includes_start(ends) || fabsf(x0 - (branch + t0) * half) > reach))
      continue;
    // The next span's start test reads this residual, so exactly one of the
    // two emits a pass at the shared column.
    const float end_residual = x1 - (branch + t1) * half;
    if ((end_residual * approach >= 0.0f && !includes_end(ends)) ||
        (raw > 1.0f && fabsf(end_residual) > reach))
      continue;
    const float s = hs::clamp(raw, 0.0f, 1.0f);
    emit(Crossing{t0 + span * s, s, branch});
  }
}

/**
 * @brief Arm passes of a column track linear over part of the frame.
 * @tparam W Canvas width.
 * @tparam Emit Callable `(const Crossing &)`.
 * @param x0,x1 Unwrapped columns at fractions @p t0 and @p t1.
 * @param t0,t1 Frame fractions the track spans, `t0 < t1`.
 * @param ends START emits at @p t0 a pass within @p reach columns before
 *        the span's first column; END emits at @p t1 a pass within @p reach
 *        columns after its last column. Include an end at the frame's start
 *        or end, or when the caller confines each pass to the span's columns;
 *        otherwise a pass at @p t1 belongs to the following span.
 * @param emit Receives each pass: none, one, or two on either side of an arm
 *        cut.
 * @param reach Extra columns whose passes can be emitted at the span ends.
 * @details Arm `n` lights unwrapped column `W / 2 * (n + t)` at fraction `t`,
 * so with `s = (t - t0) / (t1 - t0)` the pass is
 * `s = (x0 - (n + t0) W / 2) / ((t1 - t0) W / 2 - (x1 - x0))`, one division
 * per arm. A track keeping pace with the arm has no pass.
 */
template <int W, typename Emit>
HS_HOT_INLINE void crossings(float x0, float x1, float t0, float t1,
                             FrameEnds ends, Emit &&emit, float reach = 1.0f) {
  crossings(W, x0, x1, t0, t1, ends, std::forward<Emit>(emit), reach);
}

/**
 * @brief Arm passes of a column track piecewise linear between knots evenly
 *        spaced over the frame.
 * @tparam W Canvas width.
 * @tparam Emit Callable `(const Crossing &)`.
 * @param columns Unwrapped column at each knot `i / (knots - 1)`, each within
 *        a quarter turn of the one before.
 * @param knots Knot count, at least 2.
 * @param ends Span ends included by crossings(); START applies to the first
 *        piece and END to the last.
 * @param emit Receives each pass of each piece, its `along` relative to that
 *        piece; a pass at an inner knot is emitted once.
 */
template <int W, typename Emit>
HS_HOT_INLINE void piecewise_crossings(const float *columns, int knots,
                                       FrameEnds ends, Emit &&emit) {
  const int pieces = knots - 1;
  for (int i = 0; i < pieces; ++i)
    crossings<W>(
        columns[i], columns[i + 1],
        static_cast<float>(i) / static_cast<float>(pieces),
        static_cast<float>(i + 1) / static_cast<float>(pieces),
        static_cast<FrameEnds>(int(includes_start(ends) && i == 0) |
                               int(includes_end(ends) && i + 1 == pieces) << 1),
        emit);
}

/** @brief Solves a monotone quadratic screen trajectory's arm time.
 * @tparam W Display width.
 * @param origin Frame-start column.
 * @param velocity Frame-end minus frame-start column.
 * @param curve Coefficient of t(t-1).
 * @param arm_start Unwrapped first column of the arm pass.
 * @pre The trajectory stays inside this arm and its column velocity stays
 * below W/2 throughout the frame.
 * @return The crossing fraction. */
template <int W>
HS_HOT_INLINE float quadratic_screen_time(float origin, float velocity,
                                          float curve, float arm_start) {
  const float distance = origin - arm_start;
  const float approach = W / 2.0f - velocity + curve;
  const float discriminant = approach * approach - 4.0f * curve * distance;
  HS_AUDIT_CHECK(discriminant >= 0.0f,
                 "quadratic_screen_time: trajectory never meets the arm");
  const float time = 2.0f * distance / (approach + sqrtf(discriminant));
  HS_AUDIT_CHECK(time >= 0.0f && time <= 1.0f,
                 "quadratic_screen_time: crossing outside the frame");
  return time;
}

/** @brief One arm pass of a point moving linearly in screen space. */
struct ScreenCrossing {
  math::PixelCoords position; /**< Screen point at the pass. */
  float time;                 /**< Owning sweep fraction. */
};

/** @brief Arm passes of a screen trajectory linear over part of the frame.
 * @tparam W Canvas width.
 * @param a,b Projected positions at fractions @p t0 and @p t1, columns
 *        unwrapped against each other.
 * @param t0,t1 Frame fractions the trajectory spans, `t0 < t1`.
 * @param ends As crossings().
 * @param emit Receives each arm pass, its time and position on the
 *        trajectory. */
template <int W, typename Emit>
HS_HOT_INLINE void screen_crossings(const math::PixelCoords &a,
                                    const math::PixelCoords &b, float t0,
                                    float t1, FrameEnds ends, Emit &&emit) {
  const float dx = b.x - a.x, dy = b.y - a.y;
  crossings<W>(a.x, b.x, t0, t1, ends, [&](const Crossing &pass) HS_HOT_INLINE {
    emit(ScreenCrossing{{a.x + dx * pass.along, a.y + dy * pass.along},
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
HS_HOT_INLINE void screen_crossings(const math::PixelCoords &a,
                                    const math::PixelCoords &b, Emit &&emit) {
  screen_crossings<W>(a, b, 0.0f, 1.0f, FrameEnds::BOTH,
                      std::forward<Emit>(emit));
}

/** @brief A point's longitude over a span, modelled as linear in time. */
struct LongitudeTrack {
  float origin; ///< Column at the span's start, in [0, W).
  float travel; ///< Column change over the span, unwrapped.
  float rate;   ///< Pass columns per column of model error.
  float guard;  ///< Bound on a pass's column error under the model.
};

/** @brief Horizontal distance and chord terms at a segment's closest approach. */
struct HorizontalApproach {
  float chord2;    ///< Squared horizontal chord length.
  float dot0;      ///< Horizontal dot product of the origin and chord.
  float distance2; ///< Minimum squared horizontal distance from the axis.
};

/** @brief Finds the closest approach of a segment to the polar axis.
 * @param origin Segment start.
 * @param dx,dz Horizontal displacement from start to end.
 * @return Squared horizontal distance and the chord terms used to compute it. */
HS_HOT_INLINE inline HorizontalApproach
closest_horizontal_approach(const math::Vector &origin, float dx, float dz) {
  const float chord2 = dx * dx + dz * dz;
  const float dot0 = origin.x * dx + origin.z * dz;
  const float at = chord2 > 0.0f ? hs::clamp(-dot0 / chord2, 0.0f, 1.0f) : 0.0f;
  const float x = origin.x + dx * at, z = origin.z + dz * at;
  return {chord2, dot0, x * x + z * z};
}

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
  track.origin = math::vector_to_column<W>(from);
  track.travel = math::centered_remainder(
      math::vector_to_column<W>(to) - track.origin, static_cast<float>(W));
  const float cx = to.x - from.x, cz = to.z - from.z;
  const auto approach = closest_horizontal_approach(from, cx, cz);
  const float near2 = approach.distance2;
  if (!(4.0f * deviation * deviation <= near2 &&
        fabsf(track.travel) < 0.5f * sweep))
    return false;
  const float near = sqrtf(near2);
  const float radians = (0.25f * fabsf(from.x * cz - from.z * cx) *
                             sqrtf(approach.chord2) / near2 +
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
 * @tparam Emit Callable `(const Crossing &)`.
 * @param track Track over the span, as model_track() gives it.
 * @param t0,t1 Frame fractions the track spans.
 * @param ends As crossings().
 * @param emit Receives each pass; the point's true pass lies within
 *        `track.guard` columns of the emitted one.
 * @details Uses `track.guard` as the reach at the frame ends so model
 * error can carry a pass into the frame.
 */
template <int W, typename Emit>
HS_HOT_INLINE void track_crossings(const LongitudeTrack &track, float t0,
                                   float t1, FrameEnds ends, Emit &&emit) {
  crossings<W>(track.origin, track.origin + track.travel, t0, t1, ends,
               std::forward<Emit>(emit), track.guard);
}

} // namespace Spacetime
