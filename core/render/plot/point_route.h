/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <type_traits>
#include <cstdint>
#include <limits>
#include <utility>

#include "core/color/color.h"
#include "core/math/3dmath.h"
#include "core/math/pixel_mapping.h"
#include "core/platform/attributes.h"
#include "core/render/canvas.h"
#include "core/render/filter/screen_anti_alias.h"
#include "core/render/spacetime.h"
#include "core/render/spacetime/crossing.h"

/**
 * @file point_route.h
 * @brief Plot::PointRoute: clip gate and closed-form column-time splat for
 * single points with bounded polynomial motion over the frame.
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
 * (Spacetime::screen_crossings()). The pass is splatted through the sink's own
 * anti-alias stage at its position then. Taps share that position within the
 * physical tolerance; faster motion and arm cuts extrapolate projected velocity.
 *
 * Over a frame in which a point moves an arc `m`, a rigid rotation's path
 * leaves its chord by at most `m^2 / 8`, and the projected path leaves the
 * line between the projected endpoints by about `m^2 cos(phi) / (8 sin(phi))`
 * of arc at colatitude `phi`. Departure and tap error are measured in
 * physical row-pitch pixels: azimuth displacement is scaled by sin(phi).
 * Paths exceeding BEND_PIXELS are projected in analytically sized pieces.
 */
template <int W, int H> class PointRoute {
public:
  /** @brief Bounds projected speed in physical row-pitch pixels per frame.
   * @param dx,dy Column and row velocities.
   * @param row_lo,row_hi Latitude range of the trajectory.
   * @return Sum of the physical horizontal and vertical speed bounds. */
  HS_HOT_INLINE static float physical_speed(float dx, float dy, float row_lo,
                                            float row_hi) {
    const float EQUATOR = math::phi_to_y<H>(math::PI_F * 0.5f);
    const float nearest = hs::clamp(EQUATOR, row_lo, row_hi);
    const int row = hs::clamp(static_cast<int>(floorf(nearest)), 0, H - 1);
    const float sine =
        row <= EQUATOR && row + 1 >= EQUATOR
            ? 1.0f
            : fmaxf(math::TrigLUT<W, H>::sin_phi[row],
                    math::TrigLUT<W, H>::sin_phi[std::min(row + 1, H - 1)]);
    return fabsf(dy) + fabsf(dx) * sine * math::RADIANS_PER_COLUMN<W> *
                           math::ROWS_PER_RADIAN<H>;
  }

  /**
   * @brief Prepares the frame's clip gate.
   * @param canvas Canvas the frame draws into; its render clip is gated.
   * @param motion Bound on any point's displacement from its frame end.
   * @details With temporal sampling disabled the frame is drawn at t = 1 and
   * @p motion is ignored.
   */
  HS_FLASH_MEMBER void prepare(const Canvas &canvas, float motion) {
    clip_stamp = canvas.clip();
    prepared = true;
    timed = !Spacetime::flat_surface();
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const float reach = timed ? motion : 0.0f;
    const ClipRegion &clip = canvas.clip();

    // Tap rows floor(y) and floor(y) + 1 reach [first, last) for y in
    // (first - 1, last).
    const int first_row = clip.render_y_start();
    const int last_row = clip.render_y_end();
    north_limit =
        cosf(fmaxf(math::y_to_phi<H>(float(first_row - 1)) - ROW_SLACK, 0.0f));
    south_limit =
        cosf(fminf(math::y_to_phi<H>(float(last_row)) + ROW_SLACK, math::PI_F));
    projection_floor = fmaxf(BEND_PIXELS / math::ROWS_PER_RADIAN<H>,
                             fminf(math::TrigLUT<W, H>::sin_phi[0],
                                   math::TrigLUT<W, H>::sin_phi[H - 1]) -
                                 math::RADIANS_PER_ROW<H> - ROW_SLACK);
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
  HS_HOT_INLINE void plot(Sink &sink, Canvas &canvas, const math::Vector &end,
                          Start &&start, Shade &&shade) const {
    plot<Sink, Start, Shade, false>(sink, canvas, end,
                                    std::forward<Start>(start),
                                    std::forward<Shade>(shade), {});
  }

  /** @brief Plots a trajectory with a captured quadratic departure.
   * @tparam Sink,Start,Shade Destination, trajectory and shading callables.
   * @tparam CURVED Whether the midpoint departure may be nonzero.
   * @param sink,canvas Destination.
   * @param end World endpoint.
   * @param start World frame-start evaluator.
   * @param shade Sample color evaluator.
   * @param middle_offset Mid-frame departure from the endpoint chord. */
  template <typename Sink, typename Start, typename Shade, bool CURVED = true>
  HS_HOT_FLASH_MEMBER void
  plot(Sink &sink, Canvas &canvas, const math::Vector &end, Start &&start,
       Shade &&shade, const math::Vector &middle_offset) const {
    HS_CHECK(prepared && clip_stamp == canvas.clip(),
             "point route not prepared for this clip");
    const math::Vector offset = CURVED ? middle_offset : math::Vector{};
    if (!reaches_clip(end))
      return;
    constexpr bool TIMED_SHADE = std::is_invocable_v<Shade &, float>;
    if (!timed) {
      const math::PixelCoords b = math::vector_to_pixel<W, H>(end);
      Color4 c;
      if constexpr (TIMED_SHADE)
        c = shade(1.0f);
      else
        c = shade();
      sink.plot(canvas, b.x, b.y, c.color, 0.0f, c.alpha);
      return;
    }
    const math::Vector origin = start();
    const bool north = origin.y > north_limit && end.y > north_limit;
    const bool south = origin.y < south_limit && end.y < south_limit;
    if (north || south) {
      const float limit = north ? north_limit : -south_limit;
      const float y = north ? fminf(origin.y, end.y) + fminf(offset.y, 0.0f)
                            : -fmaxf(origin.y, end.y) - fmaxf(offset.y, 0.0f);
      const float horizontal =
          fmaxf(sqrtf(origin.x * origin.x + origin.z * origin.z),
                sqrtf(end.x * end.x + end.z * end.z)) +
          sqrtf(offset.x * offset.x + offset.z * offset.z);
      if (limit >= 0.0f && y > 0.0f &&
          horizontal * horizontal * limit * limit <
              y * y * (1.0f - limit * limit))
        return;
    }
    const math::PixelCoords b = math::vector_to_pixel<W, H>(end);
    math::PixelCoords a = math::vector_to_pixel<W, H>(origin);
    a.x = b.x + math::centered_remainder(a.x - b.x, float(W));
    if (a.x == b.x && a.y == b.y && math::dot(offset, offset) == 0.0f) {
      Color4 c;
      if constexpr (TIMED_SHADE)
        c = shade(Spacetime::SamplingSurface<W>::time(b.x));
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
#if HS_RUNTIME_DISPLAY_GEOMETRY
    const float PHYSICAL_GAIN = (1.0f + math::PI_F) * math::ROWS_PER_RADIAN<H>;
#else
    constexpr float PHYSICAL_GAIN =
        (1.0f + math::PI_F) * math::ROWS_PER_RADIAN<H>;
#endif
    const float norm2 = 1.0f - 0.25f * math::dot(end - origin, end - origin);
    const float offset2 = math::dot(offset, offset);
    bool chord_holds = offset2 == 0.0f && small_chord(origin, end);
    if (chord_holds) {
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::point_small_chords>();
    } else {
      const float acceleration = projected_acceleration(origin, end);
      const float remaining = BEND_PIXELS - acceleration * 0.125f;
      chord_holds =
          remaining >= 0.0f &&
          offset2 * (PHYSICAL_GAIN * PHYSICAL_GAIN + remaining * remaining) <=
              norm2 * remaining * remaining;
    }
    if (!chord_holds) {
      constexpr float CURVE_ERROR = 0.002f;
      const float ignored_error =
          offset2 > 0.0f && offset2 * (PHYSICAL_GAIN * PHYSICAL_GAIN +
                                       CURVE_ERROR * CURVE_ERROR) <=
                                norm2 * CURVE_ERROR * CURVE_ERROR
              ? CURVE_ERROR
              : 0.0f;
      plot_pieces(sink, canvas, origin, end, a, b, shade_once, offset,
                  ignored_error);
      return;
    }
    Spacetime::count_geometry<
        &Spacetime::GeometryCounters::point_capture_solves>();
    Spacetime::screen_crossings<W>(
        a, b, [&](const Spacetime::ScreenCrossing &capture) HS_HOT_INLINE {
          emit(sink, canvas, capture, shade_once(capture.time), dx, dy);
        });
  }

private:
  /// Largest departure, in pixels, of a projected piece from its line.
  static constexpr float BEND_PIXELS = 0.02f;

public:
  /** @brief Certifies the physical curvature bound of a short world chord.
   * @param origin,end Unit frame endpoints.
   * @return Whether projected_acceleration() is at most eight BEND_PIXELS.
   * @details For squared chord length m2 <= 0.01, both physical second
   * derivatives are bounded by 2 m2 / (norm2 rho). Squaring the certificate
   * avoids projection, roots and divisions on accepted chords. */
  HS_HOT_INLINE static bool small_chord(const math::Vector &origin,
                                        const math::Vector &end) {
    const math::Vector delta = end - origin;
    const float m2 = math::dot(delta, delta);
    const float norm2 = 1.0f - 0.25f * m2;
    const float horizontal2 = delta.x * delta.x + delta.z * delta.z;
    const float rho2 =
        fmaxf(0.0f, fminf(origin.x * origin.x + origin.z * origin.z,
                          end.x * end.x + end.z * end.z) -
                        0.25f * horizontal2);
    const float GAIN = math::ROWS_PER_RADIAN<H> / (4.0f * BEND_PIXELS);
    return m2 <= 0.01f && m2 * m2 * (GAIN * GAIN) <= rho2 * norm2 * norm2;
  }

  /**
   * @brief Bounds projected chord acceleration in physical row-pitch pixels.
   * @param origin,end Unit directions at the frame endpoints.
   * @param horizontal_floor Lower horizontal norm of the visible latitudes.
   * @return Bound on either projected second derivative, with longitude
   *         scaled by the direction's horizontal norm.
   * @details For a chord q, longitude acceleration is
   * `-2 cross(q_h, q_h') dot(q_h, q_h') / |q_h|^4`.
   * The closest horizontal approach bounds the denominator over the span.
   * A span of duration T departs from its projected chord by at most
   * acceleration times T squared divided by eight.
   */
  HS_HOT_INLINE static float
  projected_acceleration(const math::Vector &origin, const math::Vector &end,
                         float horizontal_floor = 0.0f) {
    const math::Vector v = end - origin;
    const float v2 = math::dot(v, v);
    if (!(v2 > 0.0f))
      return 0.0f;
    const float norm2 = fmaxf(0.0f, 1.0f - 0.25f * v2);
    const float horizontal_v2 = v.x * v.x + v.z * v.z;
    const float dot0 = origin.x * v.x + origin.z * v.z;
    const float at = horizontal_v2 > 0.0f
                         ? hs::clamp(-dot0 / horizontal_v2, 0.0f, 1.0f)
                         : 0.0f;
    const float x = origin.x + v.x * at, z = origin.z + v.z * at;
    const float rho2 =
        fmaxf(x * x + z * z, horizontal_floor * horizontal_floor * norm2);
    if (!(rho2 > 0.0f && norm2 > 0.0f))
      return std::numeric_limits<float>::infinity();
    const float cross = origin.x * v.z - origin.z * v.x;
    const float dot1 = dot0 + horizontal_v2;
    const float dot_max = fmaxf(fabsf(dot0), fabsf(dot1));
    const float dot_min =
        dot0 * dot1 <= 0.0f ? 0.0f : fminf(fabsf(dot0), fabsf(dot1));
    const float peak_dot =
        hs::clamp(fabsf(cross) * 0.7071067812f, dot_min, dot_max);
    const float peak_rho2 =
        horizontal_v2 > 0.0f
            ? (peak_dot * peak_dot + cross * cross) / horizontal_v2
            : rho2;
    const float rho3 = rho2 * sqrtf(rho2);
    const float longitude =
        horizontal_floor * horizontal_floor * norm2 > x * x + z * z
            ? 2.0f * fabsf(cross) * dot_max / (sqrtf(norm2) * rho3)
        : peak_rho2 > 0.0f ? 2.0f * fabsf(cross) * peak_dot /
                                 (sqrtf(norm2 * peak_rho2) * peak_rho2)
                           : 0.0f;
    const float latitude =
        cross * cross * fmaxf(fabsf(origin.y), fabsf(end.y)) / (norm2 * rho3) +
        v2 * sqrtf(v2) / (norm2 * norm2);
    return fmaxf(longitude, latitude) * math::ROWS_PER_RADIAN<H>;
  }

private:
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
   * @brief Splats crossings of physically bounded trajectory pieces.
   * @param sink Prepared sink.
   * @param canvas Destination canvas.
   * @param origin,end World positions at t = 0 and t = 1.
   * @param a,b Their projections, columns unwrapped against each other.
   * @param shade Returns the point's color and opacity at a frame fraction.
   * @param middle_offset Quadratic world departure from the endpoint chord.
   * @param ignored_error Physical error reserved for an omitted tiny departure.
   * @details The physical curvature bound selects the number of pieces.
   * Each piece solves its arm crossing once and emits one spatial AA splat.
   */
  template <typename Sink, typename Shade>
  HS_FLASH_MEMBER void
  plot_pieces(Sink &sink, Canvas &canvas, const math::Vector &origin,
              const math::Vector &end, const math::PixelCoords &a,
              const math::PixelCoords &b, Shade &shade,
              const math::Vector &middle_offset, float ignored_error) const {
    const float offset =
        ignored_error > 0.0f ? 0.0f : middle_offset.magnitude();
    const int curves =
        offset > 0.0f
            ? std::max(2, static_cast<int>(ceilf(
                              sqrtf(2.0f * offset * math::ROWS_PER_RADIAN<H> /
                                    BEND_PIXELS))))
            : 1;
    const math::Vector chord = end - origin;
    math::Vector from_world = origin;
    math::PixelCoords from = a;
    for (int curve = 0; curve < curves; ++curve) {
      const float begin = float(curve) / float(curves);
      const float finish = float(curve + 1) / float(curves);
      const math::Vector to_world =
          curve + 1 == curves
              ? end
              : (origin + chord * finish +
                 middle_offset * (4.0f * finish * (1.0f - finish)))
                    .normalized();
      const float acceleration =
          projected_acceleration(from_world, to_world, projection_floor);
      const int pieces = static_cast<int>(ceilf(fmaxf(
          1.0f, sqrtf(acceleration / (8.0f * (BEND_PIXELS - ignored_error))))));
      const math::Vector span = to_world - from_world;
      const math::Vector start_world = from_world;
      const float step = (finish - begin) / float(pieces);
      for (int k = 1; k <= pieces; ++k) {
        const float t0 = begin + float(k - 1) * step;
        const float t1 = k == pieces ? finish : begin + float(k) * step;
        Spacetime::count_geometry<
            &Spacetime::GeometryCounters::point_piece_evaluations>();
        math::PixelCoords to =
            curve + 1 == curves && k == pieces
                ? b
                : math::vector_to_pixel<W, H>(
                      (start_world + span * (float(k) / float(pieces)))
                          .normalized());
        to.x = from.x + math::centered_remainder(to.x - from.x, float(W));
        const float dx = (to.x - from.x) / (t1 - t0),
                    dy = (to.y - from.y) / (t1 - t0);
        Spacetime::count_geometry<
            &Spacetime::GeometryCounters::point_capture_solves>();
        Spacetime::screen_crossings<W>(
            from, to, t0, t1, true, true,
            [&](const Spacetime::ScreenCrossing &capture) HS_HOT_INLINE {
              emit(sink, canvas, capture, shade(capture.time), dx, dy, t0, t1);
            });
        from = to;
      }
      from_world = to_world;
    }
  }

public:
  /**
   * @brief Splats one arm pass.
   * @param sink Prepared sink.
   * @param canvas Destination canvas.
   * @param capture The pass.
   * @param c Color and opacity.
   * @param dx,dy Column and row motion per frame at the pass.
   * @param t0,t1 Time interval owned by this trajectory piece.
   * @details Taps within BEND_PIXELS physical error share the crossing position.
   * The complete support spans at most 1 / abs(W/2 - dx) frame fractions
   * from the crossing. Faster motion and arm cuts use each tap's column time.
   */
  template <typename Sink>
  HS_HOT_INLINE static void
  emit(Sink &sink, Canvas &canvas, const Spacetime::ScreenCrossing &capture,
       const Color4 &c, float dx, float dy, float t0 = 0.0f, float t1 = 1.0f) {
    math::PixelCoords position = capture.position;
    if (position.x < -float(W) || position.x >= 2.0f * W)
      position.x = math::centered_remainder(position.x, float(W));

    constexpr float HALF = W / 2.0f;
    const float approach = HALF - dx;
    const float tap_error =
        t0 == 0.0f && t1 == 1.0f ? BEND_PIXELS : 0.5f * BEND_PIXELS;
    const float allowance = tap_error * approach;
    const bool shared =
        fabsf(dy) + fabsf(dx) * math::RADIANS_PER_COLUMN<W> *
                        math::ROWS_PER_RADIAN<H> <=
            allowance ||
        physical_speed(dx, dy, position.y, position.y) <= allowance;
    const float fraction = position.x - floorf(position.x);
    const float lag = fmaxf(fraction, 1.0f - fraction);
    if (fabsf(math::centered_remainder(position.x, W / 2.0f)) > 1.0f &&
        approach > 0.0f && shared && capture.time - lag / approach >= t0 &&
        capture.time + lag / approach < t1) {
      sink.plot(canvas, position.x, position.y, c.color, 0.0f, c.alpha);
      return;
    }
    emit_columns(sink, canvas, capture, position, c, dx, dy, t0, t1);
  }

private:
  /** @brief Emits the complete column support of a moving splat.
   * @tparam Sink Prepared screen sink.
   * @param sink,canvas Destination.
   * @param capture Solved arm pass.
   * @param position Unwrapped screen position at the pass.
   * @param c Color and opacity.
   * @param dx,dy Projected frame velocity.
   * @param t0,t1 Time interval owned by the trajectory piece. */
  template <typename Sink>
  HS_HOT_FLASH_MEMBER static void
  emit_columns(Sink &sink, Canvas &canvas,
               const Spacetime::ScreenCrossing &capture,
               const math::PixelCoords &position, const Color4 &c, float dx,
               float dy, float t0, float t1) {
    constexpr float HALF = W / 2.0f;
    const float approach = HALF - dx;
    using Splat = Filter::Screen::AntiAlias<W, H>;
    auto pass = [&](float x, float y, const Pixel &color, float age,
                    float alpha) HS_HOT_INLINE {
      sink.plot(canvas, static_cast<int>(x), static_cast<int>(y), color, age,
                alpha);
    };
    if (fabsf(approach) < 1e-6f)
      return;
    // |x(t(column)) - column| < 1 is the complete linear splat support.
    const int arm =
        static_cast<int>(floorf(position.x / HALF - capture.time + 0.5f));
    const float center =
        arm * HALF +
        HALF * (position.x - dx * capture.time - arm * HALF) / approach;
    const float radius = HALF / fabsf(approach);
    const int first = std::max(static_cast<int>(ceilf(center - radius)),
                               static_cast<int>(ceilf((arm + t0) * HALF)));
    const int last = std::min(static_cast<int>(floorf(center + radius)),
                              static_cast<int>(ceilf((arm + t1) * HALF)) - 1);
    for (int unwrapped = first; unwrapped <= last; ++unwrapped) {
      if (unwrapped < arm * (W / 2) || unwrapped >= (arm + 1) * (W / 2))
        continue;
      const int column = math::fast_wrap(unwrapped, W);
      const float at = Spacetime::SamplingSurface<W>::time(column);
      const float lag = at - capture.time;
      float x = position.x + dx * lag;
      if (x < -float(W) || x >= 2.0f * W)
        x = math::centered_remainder(x, float(W));
      Splat::column_taps(x, position.y + dy * lag, c.color, 0.0f, c.alpha,
                         column, pass);
    }
  }

private:
  ClipRegion clip_stamp;    ///< Clip captured by prepare().
  bool prepared = false;    ///< Whether prepare() has run.
  float north_limit = 1.0f; ///< Northernmost splat latitude in the render rows.
  float south_limit =
      -1.0f; ///< Southernmost splat latitude in the render rows.
  float projection_floor = 0.0f; ///< Horizontal norm of the visible pole cap.
  float y_max = 0.0f;            ///< Largest unit y that may reach the clip.
  float y_min = 0.0f;            ///< Smallest unit y that may reach the clip.
  float low_x = 0.0f;        ///< Unit x of the band's first bounding meridian.
  float low_z = 0.0f;        ///< Unit z of the band's first bounding meridian.
  float high_x = 0.0f;       ///< Unit x of the band's last bounding meridian.
  float high_z = 0.0f;       ///< Unit z of the band's last bounding meridian.
  float column_reach = 0.0f; ///< Motion admitted outside the meridians.
  ColumnTest column_test = ColumnTest::NONE; ///< Column gate combination.
  bool timed = false; ///< Whether the frame is drawn at column time.
};

} // namespace Plot
