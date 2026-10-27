/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include "core/render/scan/radial_surface_math.h"

/** @file
 * @brief Column-time spherical rendering of paired radial surfaces.
 */

#include "core/render/scan/radial_surface_types.h"
#include "core/render/scan/radial_surface.h"
#include "core/render/scan/radial_surface_edge_aa.h"
#include "core/render/scan/shader.h"

namespace Scan {
/** @brief A planar chart's local pixel footprint. */
template <class Point> struct RadialChartFootprint {
  Point p, dx, dy;
  float radius, extent;
  bool far;
};

/** @brief Projects a tangent footprint into paired azimuthal charts.
 * @tparam Point Chart point type.
 * @param p Unit direction in {u,v,w} coordinates.
 * @param dx,dy Spherical pixel derivatives in the same coordinates.
 * @return Local chart footprint and hemisphere.
 */
template <class Point>
HS_HOT_INLINE RadialChartFootprint<Point>
radial_chart_footprint(math::Vector p, math::Vector dx, math::Vector dy) {
  const bool far = p.y < 0;
  if (far) {
    p.x = -p.x;
    p.y = -p.y;
    dx.x = -dx.x;
    dx.y = -dx.y;
    dy.x = -dy.x;
    dy.y = -dy.y;
  }
  const float sine = sqrtf(p.x * p.x + p.z * p.z);
  const float radius = radial_atan2(sine, p.y);
  const float scale = sine > 1e-5f ? radius / sine : 1;
  const float derivative =
      sine > 1e-3f ? (radius * p.y - sine) / (sine * sine * sine) : -1.0f / 3;
  auto project_derivative = [&](math::Vector d) HS_HOT_INLINE {
    return Point{scale * d.x + p.x * d.y * derivative,
                 scale * d.z + p.z * d.y * derivative};
  };
  const Point cx = project_derivative(dx), cy = project_derivative(dy);
  return {{p.x * scale, p.z * scale},
          cx,
          cy,
          radius,
          .5f * (sqrtf(cx.x * cx.x + cx.y * cx.y) +
                 sqrtf(cy.x * cy.x + cy.y * cy.y)),
          far};
}

/** @brief Draws two hemisphere charts using the common column-time scan.
 * @tparam W,H Canvas dimensions.
 * @tparam Surface Prepared similar-polygon family.
 * @tparam Area Pixel coverage integrator.
 * @tparam Prepare Per-column-time family preparation callable.
 * @tparam LowerBound Conservative radial interval lookup callable.
 * @tparam Shade First-arrival color callable.
 * @param canvas Destination.
 * @param near,far Families centered at the two chart poles.
 * @param area Reusable area workspace.
 * @param basis Frame-end chart basis.
 * @param columns Exact column-time pullback, or empty for a frozen frame.
 * @param prepare Prepares both families for an absolute first-arm column.
 * @param lower_bound Maps footprint radius to its first possible interval.
 * @param shade Maps a hit to a linear color.
 * @param stats Optional coverage census.
 */
template <int W, int H, class Surface, class Area, class Prepare, class Ensure,
          class LowerBound, class Shade>
HS_HOT_FLASH_MEMBER void draw_radial_surface(
    Canvas &canvas, Surface &near, Surface &far, Area &area,
    const math::Basis &basis, const ColumnBasis &columns, Prepare &&prepare,
    Ensure &&ensure, LowerBound &&lower_bound, Shade &&shade,
    SurfaceAreaStats *stats = nullptr, float row_displacement = 0) {
  using Point = typename Surface::Point;
  const float displacement_sine =
      sinf(-row_displacement * math::RADIANS_PER_ROW<H>);
  const float displacement_cosine =
      cosf(-row_displacement * math::RADIANS_PER_ROW<H>);
  bool near_ready = false, far_ready = false;
  Shader::draw_column_times<W, H>(
      canvas,
      [&](int column) {
        near_ready = far_ready = false;
        prepare(column);
        return column;
      },
      [&](int, int x) {
        math::Vector radial, up;
        if (columns.azimuth) {
          const int index = x % (W / 2);
          radial = x < W / 2 ? columns.azimuth[index] : -columns.azimuth[index];
          up = columns.up[index];
        } else {
          radial = {math::TrigLUT<W, H>::cos_theta(x), 0,
                    math::TrigLUT<W, H>::sin_theta[x]};
          up = math::Y_AXIS;
        }
        const auto tangent = math::cross(radial, up);
        auto local = [&](math::Vector v) {
          return math::Vector(math::dot(v, basis.u), math::dot(v, basis.v),
                              math::dot(v, basis.w));
        };
        const auto r = local(radial), v = local(up), t = local(tangent);
        return [&, r, v, t](float sine, float cosine) HS_HOT_INLINE {
          const float shifted_sine =
              sine * displacement_cosine + cosine * displacement_sine;
          cosine = cosine * displacement_cosine - sine * displacement_sine;
          sine = shifted_sine;
          const auto p = r * sine + v * cosine;
          const auto dx = t * (sine * math::TWO_PI_F / W);
          const auto dy = (r * cosine - v * sine) * math::RADIANS_PER_ROW<H>;
          const auto footprint = radial_chart_footprint<Point>(p, dx, dy);
          bool &ready = footprint.far ? far_ready : near_ready;
          if (!ready) {
            ensure(footprint.far);
            ready = true;
          }
          const auto &surface = footprint.far ? far : near;
          const float theta = radial_atan2(footprint.p.y, footprint.p.x);
          const int first =
              lower_bound(fmaxf(0, footprint.radius - footprint.extent));
          SurfaceColor result;
          result =
              area(surface, footprint.p, footprint.dx, footprint.dy, theta,
                   first, shade, stats, footprint.radius, footprint.extent);
          return Pixel(uint16_t(hs::clamp(result.r, 0.0f, 65535.0f)),
                       uint16_t(hs::clamp(result.g, 0.0f, 65535.0f)),
                       uint16_t(hs::clamp(result.b, 0.0f, 65535.0f)));
        };
      });
}
} // namespace Scan
