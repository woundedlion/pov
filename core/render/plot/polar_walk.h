/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file polar_walk.h
 * @brief Column-time strokes in a display pole's equidistant chart.
 */

#include "core/math/pixel_mapping.h"
#include "core/render/spacetime.h"
#include "core/render/canvas.h"
#include "core/render/plot/cull.h"
#include "core/render/plot/raster.h"

namespace Plot {

/** @brief Quadratic planar trajectory through three frame captures. */
struct PolarAnchor {
  math::PixelCoords end;      ///< Frame-end chart position in row pitches.
  math::PixelCoords velocity; ///< Linear coefficient about frame end.
  math::PixelCoords bend;     ///< Quadratic coefficient about frame end.

  /** @brief Fits frame-start, middle and end chart positions.
   * @param a,m,b Positions at times zero, one half and one.
   * @return Quadratic capture. */
  HS_HOT_INLINE static PolarAnchor
  capture(math::PixelCoords a, math::PixelCoords m, math::PixelCoords b) {
    const math::PixelCoords bend{2.0f * (a.x + b.x - 2.0f * m.x),
                                 2.0f * (a.y + b.y - 2.0f * m.y)};
    return {b, {b.x - a.x + bend.x, b.y - a.y + bend.y}, bend};
  }

  /** @brief Evaluates the captured trajectory.
   * @param t Frame fraction.
   * @return Chart position. */
  HS_HOT_INLINE math::PixelCoords at(float t) const {
    const float u = t - 1.0f;
    return {end.x + u * (velocity.x + u * bend.x),
            end.y + u * (velocity.y + u * bend.y)};
  }
};

/** @brief Meridian intersections and lattice taps of chart-straight pieces.
 * @tparam W,H Display dimensions. */
template <int W, int H> struct PolarWalk {
  /** @brief Projects a unit direction into a display pole's chart.
   * @param p Unit world direction.
   * @param north Selects the north pole.
   * @return Chart coordinates in row pitches. */
  HS_HOT_INLINE static math::PixelCoords project(const math::Vector &p,
                                                 bool north) {
    const float r = sqrtf(p.x * p.x + p.z * p.z);
    const float rho =
        math::precise_atan2(r, north ? p.y : -p.y) * math::ROWS_PER_RADIAN<H>;
    const float scale = r > 1e-9f ? rho / r : 0.0f;
    return {p.x * scale, p.z * scale};
  }

  /** @brief Intersects a planar segment with a diameter.
   * @param a,b Segment endpoints.
   * @param dx,dy Unit diameter direction.
   * @param s Receives the segment fraction.
   * @return Whether the half-open segment crosses the diameter. */
  HS_HOT_INLINE static bool crossing(math::PixelCoords a, math::PixelCoords b,
                                     float dx, float dy, float &s) {
    const float ca = a.x * dy - a.y * dx;
    const float cb = b.x * dy - b.y * dx;
    const float denominator = ca - cb;
    if (fabsf(denominator) < 1e-9f)
      return false;
    s = ca / denominator;
    return s >= 0.0f && s < 1.0f;
  }

  /** @brief Emits a moving chart segment at each supported column's time.
   * @tparam Sink Prepared screen sink with owned-column taps.
   * @param sink,canvas Destination sink and canvas.
   * @param a,b Captured segment endpoints.
   * @param north Selects the chart pole.
   * @param color Stroke registers.
   * @param step Screen lattice spacing.
   * @param balanced Compensate short samples with the balanced opacity law.
   * @param density Reference samples per unit lattice length.
   */
  template <typename Sink>
  HS_HOT_FLASH_MEMBER static void
  draw(Sink &sink, Canvas &canvas, const PolarAnchor &a, const PolarAnchor &b,
       bool north, const Color4 &color, float step, bool balanced = true,
       float density = 1.0f) {
    const ColumnWindow window = column_window(a, b);
    const int first = window.first, count = window.count;
    const float pole = north ? math::DisplayGeometry<H>::NORTH_POLE_ROW
                             : math::DisplayGeometry<H>::SOUTH_POLE_ROW;
    const float sign = north ? 1.0f : -1.0f;
    const float gain =
        density / step *
        balanced_sample_alpha(color.alpha, step / SCREEN_STEP_PX);
    for (int i = 0; i < count; ++i) {
      const int column = math::fast_wrap(first + i, W);
      if (canvas.render_x_clip().clipped(column))
        continue;
      const float t = Spacetime::SamplingSurface<W>::time(column);
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::polar_column_visits>();
      const auto from = a.at(t), to = b.at(t);
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::polar_geometry_samples>();
      const float dx = math::TrigLUT<W, H>::cos_theta(column);
      const float dy = math::TrigLUT<W, H>::sin_theta[column];
      const float ex = to.x - from.x, ey = to.y - from.y;
      const float along = from.x * dx + from.y * dy;
      const float across = from.y * dx - from.x * dy;
      const float ea = ex * dx + ey * dy;
      const float ec = ey * dx - ex * dy;
      float s = 0.0f;
      crossing(from, to, dx, dy, s);
      const float radius = along + s * ea;
      const float slope =
          fabsf(ec) > 1e-9f
              ? fabsf(radius * ea / ec) * math::RADIANS_PER_COLUMN<W>
              : std::numeric_limits<float>::infinity();
      auto emit = [&](float row, float x, float alpha) HS_HOT_INLINE {
        Spacetime::count_geometry<
            &Spacetime::GeometryCounters::polar_column_samples>();
        sink.plot_column_at(canvas, {x, pole + sign * row}, column, color.color,
                            0.0f, alpha, t);
      };
      if (!balanced && radius > 0.0f && slope <= 1.0f && along > 0.0f &&
          along + ea > 0.0f) {
        auto integral = [](float x) HS_HOT_INLINE {
          x = hs::clamp(x, -1.0f, 1.0f);
          return x < 0.0f ? 0.5f * (x + 1.0f) * (x + 1.0f)
                          : 1.0f - 0.5f * (1.0f - x) * (1.0f - x);
        };
        const float x0 = across / (along * math::RADIANS_PER_COLUMN<W>);
        const float x1 =
            (across + ec) / ((along + ea) * math::RADIANS_PER_COLUMN<W>);
        const float coverage = fabsf(integral(x1) - integral(x0));
        if (coverage > 0.0f) {
          auto moment = [](float x) HS_HOT_INLINE {
            x = hs::clamp(x, -1.0f, 1.0f);
            return x * x * (0.5f - fabsf(x) / 3.0f);
          };
          const float mass = integral(x1) - integral(x0);
          const float mean = (moment(x1) - moment(x0)) / mass;
          const float offset = mean * math::RADIANS_PER_COLUMN<W>;
          const float denominator = ec - offset * ea;
          const float u =
              hs::clamp((offset * along - across) / denominator, 0.0f, 1.0f);
          const float x = from.x + ex * u, y = from.y + ey * u;
          const float screen_rate = sqrtf(1.0f + slope * slope);
          const float horizontal = radius * math::RADIANS_PER_COLUMN<W>;
          const float physical_rate =
              sqrtf(horizontal * horizontal + slope * slope);
          const float min_pitch = math::RADIANS_PER_COLUMN<W> *
                                  math::ROWS_PER_RADIAN<H> * MIN_POLE_SCALE;
          const float pitch_scale = fminf(1.0f, SCREEN_STEP_PX * physical_rate /
                                                    (screen_rate * min_pitch));
          emit(sqrtf(x * x + y * y), float(column),
               gain * (SCREEN_STEP_PX / step) * screen_rate * coverage *
                   pitch_scale);
        }
        continue;
      }
      // Intersect the segment with the column's triangular AA support sector.
      float lo = 0.0f, hi = 1.0f;
      const float tangent = tanf(math::RADIANS_PER_COLUMN<W>);
      auto restrict = [&](float v, float delta) HS_HOT_INLINE {
        if (fabsf(delta) < 1e-9f) {
          if (v < 0.0f)
            hi = -1.0f;
        } else if (delta > 0.0f)
          lo = fmaxf(lo, -v / delta);
        else
          hi = fminf(hi, -v / delta);
      };
      restrict(along * tangent - across, ea * tangent - ec);
      restrict(along * tangent + across, ea * tangent + ec);
      if (hi <= lo)
        continue;
      if (radius > 0.0f && slope <= 1.0f && along > 0.0f && along + ea > 0.0f) {
        const float angle =
            math::centered_remainder(math::precise_atan2(to.y, to.x) -
                                         math::precise_atan2(from.y, from.x),
                                     math::TWO_PI_F) /
            math::RADIANS_PER_COLUMN<W>;
        const float radial = sqrtf(to.x * to.x + to.y * to.y) -
                             sqrtf(from.x * from.x + from.y * from.y);
        const float length = sqrtf(angle * angle + radial * radial);
        const int samples = std::max(1, static_cast<int>(ceilf(length / step)));
        const float inverse = 1.0f / samples;
        const float ratio = fminf(step, length * inverse) / SCREEN_STEP_PX;
        const float alpha =
            density * (!balanced && ratio < 1.0f
                           ? color.alpha * ratio
                           : balanced_sample_alpha(color.alpha, ratio));
        const int begin = std::max(0, static_cast<int>(ceilf(lo * samples)));
        const int end =
            std::min(samples - 1, static_cast<int>(floorf(hi * samples)));
        for (int sample = begin; sample <= end; ++sample) {
          const float u = sample * inverse;
          const float x = from.x + ex * u, y = from.y + ey * u;
          const float forward = along + ea * u;
          if (forward <= 0.0f)
            continue;
          const float column_at =
              float(column) + math::precise_atan2(across + ec * u, forward) /
                                  math::RADIANS_PER_COLUMN<W>;
          emit(sqrtf(x * x + y * y), column_at, alpha);
        }
        continue;
      }
      const float norm = ex * ex + ey * ey;
      if (norm < 1e-12f)
        continue;
      const float foot = -(from.x * ex + from.y * ey) / norm;
      auto radius_at = [&](float u) HS_HOT_INLINE {
        const float x = from.x + ex * u, y = from.y + ey * u;
        return sqrtf(x * x + y * y);
      };
      const float rlo = radius_at(lo), rhi = radius_at(hi);
      const float rmin = radius_at(hs::clamp(foot, lo, hi));
      const float rmax = fmaxf(rlo, rhi);
      const float distance2 =
          fmaxf(0.0f, from.x * from.x + from.y * from.y - foot * foot * norm);
      float radial_step = step, radial_phase = 0.0f, alpha = gain * step;
      const float first_radius = radius_at(0.0f);
      const float last_radius = radius_at(1.0f);
      if ((foot <= 0.0f || foot >= 1.0f) &&
          fabsf(last_radius - first_radius) > 1e-6f) {
        const float angle =
            math::centered_remainder(math::precise_atan2(to.y, to.x) -
                                         math::precise_atan2(from.y, from.x),
                                     math::TWO_PI_F) /
            math::RADIANS_PER_COLUMN<W>;
        const float dr = last_radius - first_radius;
        const float length = sqrtf(angle * angle + dr * dr);
        const int samples = std::max(1, static_cast<int>(ceilf(length / step)));
        const float inverse = 1.0f / samples;
        radial_step = fabsf(dr) * inverse;
        radial_phase = first_radius;
        const float ratio = fminf(step, length * inverse) / SCREEN_STEP_PX;
        alpha = density * (!balanced && ratio < 1.0f
                               ? color.alpha * ratio
                               : balanced_sample_alpha(color.alpha, ratio));
      }
      for (int row =
               static_cast<int>(ceilf((rmin - radial_phase) / radial_step));
           radial_phase + row * radial_step <= rmax; ++row) {
        const float r = radial_phase + row * radial_step;
        const float offset = sqrtf(fmaxf(0.0f, (r * r - distance2) / norm));
        for (int branch = 0; branch < 2; ++branch) {
          const float u = foot + (branch ? offset : -offset);
          if (u < lo || u >= hi || (branch && offset == 0.0f))
            continue;
          const float forward = along + ea * u;
          if (forward <= 0.0f)
            continue;
          const float x =
              float(column) + math::precise_atan2(across + ec * u, forward) /
                                  math::RADIANS_PER_COLUMN<W>;
          emit(r, x, alpha);
        }
      }
    }
  }

private:
  /** @brief Column range and radial hull of a moving segment. */
  struct ColumnWindow {
    int first; ///< First unwrapped column.
    int count; ///< Number of columns.
  };

  /** @brief Bounds the columns reached by the endpoints' quadratic hulls.
   * @param a,b Captured segment endpoints.
   * @return Column range and radial bound. */
  HS_FLASH_MEMBER static ColumnWindow column_window(const PolarAnchor &a,
                                                    const PolarAnchor &b) {
    const float origin = math::precise_atan2(a.end.y, a.end.x);
    float low = origin, high = origin;
    auto include = [&](math::PixelCoords p) HS_HOT_INLINE {
      const float angle =
          origin + math::centered_remainder(
                       math::precise_atan2(p.y, p.x) - origin, math::TWO_PI_F);
      low = fminf(low, angle);
      high = fmaxf(high, angle);
    };
    auto enclose = [&](const PolarAnchor &p) HS_HOT_INLINE {
      include(p.end);
      include(p.at(0.0f));
      include({p.end.x - 0.5f * p.velocity.x, p.end.y - 0.5f * p.velocity.y});
    };
    enclose(a);
    enclose(b);
    const int first =
        static_cast<int>(floorf(low / math::RADIANS_PER_COLUMN<W>)) - 1;
    const int count =
        high - low >= math::PI_F
            ? W
            : std::min(W, static_cast<int>(
                              ceilf(high / math::RADIANS_PER_COLUMN<W>)) +
                              2 - first);
    return {first, count};
  }
};

} // namespace Plot
