/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <type_traits>
#include <cstdint>
#include <limits>

#include "core/color/color.h"
#include "core/math/3dmath.h"
#include "core/math/pixel_mapping.h"
#include "core/platform/attributes.h"
#include "core/render/canvas.h"
#include "core/render/filter/screen_anti_alias.h"
#include "core/render/temporal.h"
#include "core/render/temporal/arm_crossing.h"

/**
 * @file point_route.h
 * @brief Plot::PointRoute: clip gate and closed-form column-time splat for
 * single points that move along a chord over the frame.
 */

namespace Plot {

/**
 * @brief Plots single points, each at the frame fraction its column is lit.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @details A point's motion over the frame is its world positions at t = 0
 * and t = 1; between them it moves along their chord. Both endpoints are
 * projected once, and the projected column and row are taken as linear in t
 * between them, so each arm pass is one division
 * (Temporal::linear_captures()). The pass is splatted through the sink's own
 * anti-alias stage at its position then; where that would put a tap more
 * than SINGLE_TIME_PIXELS from where the point stood when the tap's column
 * was lit, or the taps straddle an arm boundary, each tap column takes the
 * point's position at that column's own fraction, read from the same line.
 *
 * Over a frame in which a point moves an arc `m`, a rigid rotation's path
 * leaves its chord by at most `m^2 / 8`, and the projected path leaves the
 * line between the projected endpoints by about `m^2 cos(phi) / (8 sin(phi))`
 * of arc at colatitude `phi`. Where that line spans more than
 * PIECE_COLUMNS columns, near a pole, the chord is projected at evenly spaced
 * fractions instead, so that each piece spans at most PIECE_COLUMNS columns,
 * and taps take their exact projected positions.
 */
template <int W, int H> class PointRoute {
public:
  /**
   * @brief Prepares the frame's clip gate.
   * @param canvas Canvas the frame draws into; its render clip is gated.
   * @param motion Bound on the chord between any point's two endpoints.
   * @details With temporal sampling disabled the frame is drawn at t = 1 and
   * @p motion is ignored.
   */
  HS_FLASH_MEMBER void prepare(const Canvas &canvas, float motion) {
    timed = Temporal::temporal_sampling_enabled();
    const float reach = timed ? motion : 0.0f;
    const ClipRegion &clip = canvas.clip();

    // Tap rows floor(y) and floor(y) + 1 reach [first, last) for y in
    // (first - 1, last).
    const int first_row = clip.render_y_start();
    const int last_row = clip.render_y_end();
    constexpr float INF = std::numeric_limits<float>::infinity();
    y_max =
        first_row > 0
            ? cosf(fmaxf(math::y_to_phi<H>(float(first_row - 1)) - ROW_SLACK,
                         0.0f)) +
                  reach
            : INF;
    y_min = last_row < H
                ? cosf(fminf(math::y_to_phi<H>(float(last_row)) + ROW_SLACK,
                             math::PI_F)) -
                      reach
                : -INF;

    // Tap columns floor(x) and floor(x) + 1 reach [first, last) for x in
    // [first - 1, last), widened by fast_atan2's error.
    const ClipRegion::XClip columns = clip.x_clip();
    column_test = ColumnTest::NONE;
    if (!columns.active)
      return;
    const int last_column = columns.wrap ? columns.re + W : columns.re;
    const float low =
        float(columns.rs - 1) * math::RADIANS_PER_COLUMN<W> - COLUMN_SLACK;
    const float high =
        float(last_column) * math::RADIANS_PER_COLUMN<W> + COLUMN_SLACK;
    if (high - low >= math::TWO_PI_F)
      return;
    low_x = cosf(low);
    low_z = sinf(low);
    high_x = cosf(high);
    high_z = sinf(high);
    column_reach = reach;
    column_test =
        high - low <= math::PI_F ? ColumnTest::BOTH : ColumnTest::EITHER;
  }

  /**
   * @brief Plots one point at each arm pass.
   * @tparam Sink Screen sink or pipeline whose float `plot()` splats through
   *         Filter::Screen::AntiAlias and whose integer `plot()` writes one
   *         pixel.
   * @tparam Start Callable `() -> math::Vector`.
   * @tparam Shade Callable `() -> Color4`, or `(float t) -> Color4` for a
   *         point whose color changes over the frame.
   * @param sink Prepared sink.
   * @param canvas Destination canvas.
   * @param end World position at t = 1, unit length.
   * @param start Returns the world position at t = 0, unit length and within
   *        the prepared motion of @p end; called only for timed points that
   *        may reach the clip.
   * @param shade Returns the point's color and opacity; pure, called only for
   *        points that are plotted: at most once, or at each pass's time when
   *        it takes the frame fraction.
   */
  template <typename Sink, typename Start, typename Shade>
  HS_HOT_FLASH_MEMBER void plot(Sink &sink, Canvas &canvas,
                                const math::Vector &end, Start &&start,
                                Shade &&shade) const {
    if (!reaches_clip(end))
      return;
    constexpr bool TIMED_SHADE = std::is_invocable_v<Shade &, float>;
    const math::PixelCoords b = math::vector_to_pixel<W, H>(end);
    if (!timed) {
      Color4 c;
      if constexpr (TIMED_SHADE)
        c = shade(1.0f);
      else
        c = shade();
      sink.plot(canvas, b.x, b.y, c.color, 0.0f, c.alpha);
      return;
    }
    const math::Vector origin = start();
    math::PixelCoords a = math::vector_to_pixel<W, H>(origin);
    a.x = b.x + Temporal::centered_remainder(a.x - b.x, float(W));
    if (a.x == b.x && a.y == b.y) {
      Color4 c;
      if constexpr (TIMED_SHADE)
        c = shade(Temporal::column_fraction<W>(b.x));
      else
        c = shade();
      sink.plot(canvas, b.x, b.y, c.color, 0.0f, c.alpha);
      return;
    }
    Color4 c;
    bool shaded = false;
    auto shade_once = [&](float t) HS_HOT_INLINE -> const Color4 & {
      if constexpr (TIMED_SHADE) {
        c = shade(t);
      } else if (!shaded) {
        c = shade();
        shaded = true;
      }
      return c;
    };
    const float dx = b.x - a.x, dy = b.y - a.y;
    if (bends(end, dx, dy)) {
      plot_pieces(sink, canvas, origin, end, a, b, shade_once);
      return;
    }
    auto at_time = [&](float t) HS_HOT_INLINE {
      return math::PixelCoords{a.x + dx * t, a.y + dy * t};
    };
    Temporal::linear_captures<W>(
        a, b, [&](const Temporal::LinearCapture &capture) HS_HOT_INLINE {
          emit(sink, canvas, capture, shade_once(capture.time), dx, dy,
               at_time);
        });
  }

private:
  /// Largest tap error a splat at one time keeps away from an arm boundary,
  /// pixels: the motion over the farther tap column's lag.
  static constexpr float SINGLE_TIME_PIXELS = 0.02f;
  /// Most columns a projected trajectory piece spans.
  static constexpr float PIECE_COLUMNS = 4.0f;
  /// Largest departure, in pixels, of a projected piece from its line.
  static constexpr float BEND_PIXELS = 0.05f;
  /// Most pieces a trajectory is projected in.
  static constexpr int MAX_PIECES = W / 8;

  /**
   * @brief Departure of a projected trajectory from the line between its
   *        projected endpoints.
   * @param end World position at t = 1.
   * @param dx,dy Column and row motion over the frame.
   * @param rho_sq Receives `sin(phi)^2` at @p end.
   * @param column_sq Receives `(4 sin(phi) e)^2` for a column departure `e`.
   * @param row_sq Receives `(8 sin(phi) e)^2` for a row departure `e`.
   * @details On the sphere `theta'' = -2 cot(phi) theta' phi'` and
   * `phi'' = sin(phi) cos(phi) theta'^2`, so over the frame the column leaves
   * the line by `|dx * dphi * cot(phi)| / 4` and the row by
   * `|sin(phi) cos(phi)| dtheta^2 / 8` radians.
   */
  HS_HOT_INLINE static void departure(const math::Vector &end, float dx,
                                      float dy, float &rho_sq, float &column_sq,
                                      float &row_sq) {
    const float lift = dy * math::RADIANS_PER_ROW<H>;
    const float turn = dx * math::RADIANS_PER_COLUMN<W>;
    const float y_sq = end.y * end.y;
    rho_sq = 1.0f - y_sq;
    const float column = dx * lift;
    const float row = turn * turn / math::RADIANS_PER_ROW<H>;
    column_sq = column * column * y_sq;
    row_sq = row * row * y_sq * rho_sq * rho_sq;
  }

  /**
   * @brief Whether a trajectory leaves the line between its projected
   *        endpoints by more than BEND_PIXELS or spans more than
   *        PIECE_COLUMNS columns.
   * @param end World position at t = 1.
   * @param dx,dy Column and row motion over the frame.
   * @return True when it is projected in pieces.
   */
  HS_HOT_INLINE static bool bends(const math::Vector &end, float dx, float dy) {
    if (fabsf(dx) > PIECE_COLUMNS)
      return true;
    float rho_sq, column_sq, row_sq;
    departure(end, dx, dy, rho_sq, column_sq, row_sq);
    constexpr float LIMIT_SQ = BEND_PIXELS * BEND_PIXELS;
    return column_sq > 16.0f * LIMIT_SQ * rho_sq ||
           row_sq > 64.0f * LIMIT_SQ * rho_sq;
  }

  /// Arc added to the gated rows for fast_acos's error.
  static constexpr float ROW_SLACK =
      2.0f * static_cast<float>(math::FAST_ACOS_MAX_ERROR);
  /// Arc added to the gated columns for fast_atan2's error.
  static constexpr float COLUMN_SLACK =
      2.0f * static_cast<float>(math::FAST_ATAN2_MAX_ERROR);

  /** @brief How the column gate combines its two meridian halfspaces. */
  enum class ColumnTest : uint8_t {
    NONE,   ///< Every column is drawn.
    BOTH,   ///< The drawn columns span at most half a turn.
    EITHER, ///< The drawn columns span more than half a turn.
  };

  /**
   * @brief Whether a point within the prepared motion of @p p may write a
   *        tap inside the clip.
   * @param p World position at t = 1.
   * @return False only when no such point reaches the clip.
   * @details A point within chord `d` of the drawn band lies within `d` of
   * each bounding meridian halfspace that contains the band and within `d`
   * of the band's unit-y range.
   */
  HS_HOT_INLINE bool reaches_clip(const math::Vector &p) const {
    if (p.y > y_max || p.y < y_min)
      return false;
    if (column_test == ColumnTest::NONE)
      return true;
    const bool after_low = low_x * p.z - low_z * p.x >= -column_reach;
    const bool before_high = p.x * high_z - p.z * high_x >= -column_reach;
    return column_test == ColumnTest::BOTH ? after_low && before_high
                                           : after_low || before_high;
  }

  /**
   * @brief Plots the passes of a trajectory that sweeps more than
   *        PIECE_COLUMNS columns, piece by piece.
   * @param sink Prepared sink.
   * @param canvas Destination canvas.
   * @param origin,end World positions at t = 0 and t = 1.
   * @param a,b Their projections, columns unwrapped against each other.
   * @param shade Returns the point's color and opacity at a frame fraction.
   */
  template <typename Sink, typename Shade>
  HS_FLASH_MEMBER static void
  plot_pieces(Sink &sink, Canvas &canvas, const math::Vector &origin,
              const math::Vector &end, const math::PixelCoords &a,
              const math::PixelCoords &b, Shade &shade) {
    const math::Vector chord = end - origin;
    auto at_time = [&](float t) HS_HOT_INLINE {
      return math::vector_to_pixel<W, H>((origin + chord * t).normalized());
    };
    float rho_sq, column_sq, row_sq;
    departure(end, b.x - a.x, b.y - a.y, rho_sq, column_sq, row_sq);
    // Departure scales as the square of a piece's share of the frame.
    const float bend = rho_sq > 0.0f ? fmaxf(column_sq / (16.0f * rho_sq),
                                             row_sq / (64.0f * rho_sq))
                                     : std::numeric_limits<float>::infinity();
    const float wanted = fmaxf(fabsf(b.x - a.x) / PIECE_COLUMNS,
                               sqrtf(sqrtf(bend) / BEND_PIXELS));
    const int pieces =
        static_cast<int>(ceilf(fminf(wanted, float(MAX_PIECES))));
    const float step = 1.0f / static_cast<float>(pieces);
    math::PixelCoords from = a;
    for (int k = 1; k <= pieces; ++k) {
      const float t0 = static_cast<float>(k - 1) * step;
      const float t1 = k == pieces ? 1.0f : static_cast<float>(k) * step;
      math::PixelCoords to = k == pieces ? b : at_time(t1);
      to.x = from.x + Temporal::centered_remainder(to.x - from.x, float(W));
      const float dx = (to.x - from.x) / (t1 - t0);
      const float dy = (to.y - from.y) / (t1 - t0);
      Temporal::linear_captures<W>(
          from, to, t0, t1, k == 1, k == pieces,
          [&](const Temporal::LinearCapture &capture) HS_HOT_INLINE {
            emit(sink, canvas, capture, shade(capture.time), dx, dy, at_time);
          });
      from = to;
    }
  }

  /**
   * @brief Splats one arm pass.
   * @param sink Prepared sink.
   * @param canvas Destination canvas.
   * @param capture The pass.
   * @param c Color and opacity.
   * @param dx,dy Column and row motion per frame at the pass.
   * @param at_time Returns the projected position at a frame fraction.
   */
  template <typename Sink, typename AtTime>
  HS_HOT_INLINE static void
  emit(Sink &sink, Canvas &canvas, const Temporal::LinearCapture &capture,
       const Color4 &c, float dx, float dy, AtTime &at_time) {
    math::PixelCoords position = capture.position;
    if (position.x < -float(W) || position.x >= 2.0f * W)
      position.x = Temporal::centered_remainder(position.x, float(W));
    if (fabsf(Temporal::centered_remainder(position.x, W / 2.0f)) > 1.0f &&
        Temporal::tap_motion_error<W>(position.x, dx, dy) <=
            SINGLE_TIME_PIXELS) {
      sink.plot(canvas, position.x, position.y, c.color, 0.0f, c.alpha);
      return;
    }
    using Splat = Filter::Screen::AntiAlias<W, H>;
    auto pass = [&](float x, float y, const Pixel &color, float age,
                    float alpha) HS_HOT_INLINE {
      sink.plot(canvas, static_cast<int>(x), static_cast<int>(y), color, age,
                alpha);
    };
    Splat::for_each_owned_column(
        position.x, capture.time, [&](int column, float at) HS_HOT_INLINE {
          math::PixelCoords p = at_time(at);
          if (p.x < -float(W) || p.x >= 2.0f * W)
            p.x = Temporal::centered_remainder(p.x, float(W));
          Splat::column_taps(p.x, p.y, c.color, 0.0f, c.alpha, column, pass);
        });
  }

  float y_max = 0.0f;        ///< Largest unit y that may reach the clip.
  float y_min = 0.0f;        ///< Smallest unit y that may reach the clip.
  float low_x = 0.0f;        ///< Unit x of the band's first bounding meridian.
  float low_z = 0.0f;        ///< Unit z of the band's first bounding meridian.
  float high_x = 0.0f;       ///< Unit x of the band's last bounding meridian.
  float high_z = 0.0f;       ///< Unit z of the band's last bounding meridian.
  float column_reach = 0.0f; ///< Motion admitted outside the meridians.
  ColumnTest column_test = ColumnTest::NONE; ///< Column gate combination.
  bool timed = false; ///< Whether the frame is drawn at column time.
};

} // namespace Plot
