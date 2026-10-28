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
#include "core/render/scan/radial_surface_candidate_mask.h"
#include "core/render/scan/shader.h"

namespace Scan {
#if HS_ENABLE_TEST_HOOKS
/** @brief Deferred chart-work census. */
struct RadialChartWork {
  static inline uint64_t footprints = 0; ///< Constructed affine footprints.
};
#endif
/** @brief A planar chart's local pixel footprint. */
template <class Point> struct RadialChartFootprint {
  Point p;      ///< Chart position.
  Point dx;     ///< Chart derivative per pixel column.
  Point dy;     ///< Chart derivative per pixel row.
  float radius; ///< Distance from the chart pole.
  float extent; ///< Enclosing footprint radius.
  bool far;     ///< Whether the point uses the far chart.
};

/** @brief A chart center with reusable projection factors.
 * @tparam Point Planar chart point type. */
template <class Point> struct RadialChartPoint {
  Point p;             ///< Chart position.
  math::Vector sphere; ///< Hemisphere-local unit direction.
  float radius;        ///< Distance from the chart pole.
  float sine;          ///< Sine of the pole distance.
  float scale;         ///< Tangential chart scale.
  bool far;            ///< Whether the point uses the far chart.
};

/** @brief Projects a unit direction into paired azimuthal charts.
 * @tparam Point Chart point type.
 * @param p Unit direction in {u,v,w} coordinates.
 * @return Chart center and hemisphere. */
template <class Point>
HS_HOT_INLINE RadialChartPoint<Point> radial_chart_point(math::Vector p) {
  const bool far = p.y < 0;
  if (far) {
    p.x = -p.x;
    p.y = -p.y;
  }
  const float sine = sqrtf(p.x * p.x + p.z * p.z);
  const float radius = radial_atan2(sine, p.y);
  const float scale = sine > 1e-5f ? radius / sine : 1;
  return {{p.x * scale, p.z * scale}, p, radius, sine, scale, far};
}

/** @brief Projects tangent derivatives at a prepared chart center.
 * @tparam Point Chart point type.
 * @param center Prepared chart projection.
 * @param dx,dy Spherical pixel derivatives in {u,v,w} coordinates.
 * @return Local chart footprint and hemisphere. */
template <class Point>
HS_HOT_INLINE RadialChartFootprint<Point>
radial_chart_footprint(const RadialChartPoint<Point> &center, math::Vector dx,
                       math::Vector dy) {
#if HS_ENABLE_TEST_HOOKS
  ++RadialChartWork::footprints;
#endif
  if (center.far) {
    dx.x = -dx.x;
    dx.y = -dx.y;
    dy.x = -dy.x;
    dy.y = -dy.y;
  }
  const auto p = center.sphere;
  const float sine = center.sine, radius = center.radius, scale = center.scale;
  const float derivative =
      sine > 1e-3f ? (radius * p.y - sine) / (sine * sine * sine) : -1.0f / 3;
  auto project_derivative = [&](math::Vector d) HS_HOT_INLINE {
    return Point{scale * d.x + p.x * d.y * derivative,
                 scale * d.z + p.z * d.y * derivative};
  };
  const Point cx = project_derivative(dx), cy = project_derivative(dy);
  return {center.p,
          cx,
          cy,
          radius,
          .5f * (sqrtf(cx.x * cx.x + cx.y * cx.y) +
                 sqrtf(cy.x * cy.x + cy.y * cy.y)),
          center.far};
}

/** @brief Projects a tangent footprint into paired azimuthal charts.
 * @tparam Point Chart point type.
 * @param p Unit direction in {u,v,w} coordinates.
 * @param dx,dy Spherical pixel derivatives in the same coordinates.
 * @return Local chart footprint and hemisphere. */
template <class Point>
HS_HOT_INLINE RadialChartFootprint<Point>
radial_chart_footprint(math::Vector p, math::Vector dx, math::Vector dy) {
  return radial_chart_footprint(radial_chart_point<Point>(p), dx, dy);
}

/** @brief Routes certified interior pixels to point sampling and candidates to edge AA.
 * @tparam Surface Prepared similar-polygon radial family.
 * @tparam COUNT Enables optional work counters.
 * @tparam MASK Enables hemisphere candidate-mask routing.
 * @tparam Mask Conservative polar candidate-mask type.
 * @details Referenced families and masks must outlive the integrator. Masks must
 * cover the current families and pixel footprints when MASK is enabled. */
template <class Surface, bool COUNT = false, bool MASK = true,
          class Mask = RadialSurfaceCandidateMask<>>
class RadialSurfaceArea {
  RadialSurfaceEdgeAA<Surface> integrator; ///< Reusable edge-AA workspace.
  const Surface &near;                     ///< Near hemisphere identity.
  const Mask &near_mask;                   ///< Near hemisphere candidates.
  const Mask &far_mask;                    ///< Far hemisphere candidates.

public:
  /** @brief Binds the hemisphere identity and candidate masks.
   * @param near Near hemisphere family.
   * @param near_mask,far_mask Hemisphere candidate masks; unused when MASK is false. */
  HS_HOT_INLINE RadialSurfaceArea(const Surface &near, const Mask &near_mask,
                                  const Mask &far_mask)
      : near(near), near_mask(near_mask), far_mask(far_mask) {}

  /** @brief Integrates one chart pixel with a deferred affine footprint.
   * @tparam Footprint Callable returning chart derivatives and extent.
   * @tparam LowerBound Conservative radial interval lookup callable.
   * @tparam Shade First-arrival color callable.
   * @param surface Prepared near or far hemisphere family.
   * @param p Chart center.
   * @param theta Center argument, radians.
   * @param radius Center distance from the chart pole, radians.
   * @param footprint Builds the footprint only for edge-AA queries.
   * @param radial_first Maps chart radius to its first possible interval.
   * @param shade Maps a surface hit to linear color.
   * @param stats Optional coverage census.
   * @return Point color or area-weighted edge color. */
  template <class Footprint, class LowerBound, class Shade>
  HS_O3_FN HS_HOT_INLINE SurfaceColor
  operator()(const Surface &surface, typename Surface::Point p, float theta,
             float radius, Footprint &footprint, LowerBound &radial_first,
             Shade &shade, SurfaceAreaStats *stats) const {
    uint32_t candidate_bits = UINT32_MAX;
    int candidate_shift = 0;
    if constexpr (MASK) {
      if constexpr (COUNT)
        if (stats)
          ++stats->mask_queries;
      const auto &mask = &surface == &near ? near_mask : far_mask;
      candidate_bits = mask.interval_bits(radius, theta);
      candidate_shift = mask.interval_shift();
      if (!candidate_bits) {
        typename Surface::Stats work;
        const auto hit = surface.template sample<COUNT, true, false>(
            p, theta, radial_first(radius), COUNT ? &work : nullptr);
        if constexpr (COUNT)
          if (stats) {
            ++stats->pixels;
            stats->intervals += work.intervals;
            stats->point_intervals += work.intervals;
            stats->max_intervals =
                std::max(stats->max_intervals, int(work.intervals));
            stats->max_point_intervals =
                std::max(stats->max_point_intervals, int(work.intervals));
          }
        return hit.found() ? shade(hit) : SurfaceColor{};
      }
      if constexpr (COUNT)
        if (stats)
          ++stats->mask_hits;
    }
    const auto query = footprint();
    const int first = radial_first(fmaxf(0, radius - query.extent));
    return integrator.template sample<COUNT>(
        surface, p, query.dx, query.dy, theta, first, shade, stats, radius,
        query.extent, candidate_bits, candidate_shift, radial_first(radius));
  }
};

/** @brief Draws two hemisphere charts using the common column-time scan.
 * @tparam W,H Canvas dimensions.
 * @tparam Surface Prepared similar-polygon family.
 * @tparam Area Pixel integrator receiving a deferred footprint callable.
 * @tparam Prepare Per-column-time family preparation callable.
 * @tparam Ensure Lazy family preparation callable.
 * @tparam LowerBound Conservative radial interval lookup callable.
 * @tparam Shade First-arrival color callable.
 * @param canvas Destination.
 * @param near,far Families centered at the two chart poles.
 * @param area Reusable area workspace.
 * @param basis Frame-end chart basis.
 * @param columns Exact column-time pullback, or empty for a frozen frame.
 * @param prepare Prepares both families for an absolute first-arm column.
 * @param ensure Prepares the near or far family before its first query in a column.
 * @param lower_bound Maps footprint radius to its first possible interval.
 * @param shade Maps a hit to a linear color.
 * @param stats Optional coverage census.
 * @param row_displacement Vertical scan offset, rows.
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
          const auto center = radial_chart_point<Point>(p);
          bool &ready = center.far ? far_ready : near_ready;
          if (!ready) {
            ensure(center.far);
            ready = true;
          }
          const auto &surface = center.far ? far : near;
          const float theta = radial_atan2(center.p.y, center.p.x);
          auto footprint = [&]() HS_HOT_INLINE {
            const auto dx = t * (sine * math::TWO_PI_F / W);
            const auto dy = (r * cosine - v * sine) * math::RADIANS_PER_ROW<H>;
            return radial_chart_footprint(center, dx, dy);
          };
          const SurfaceColor result =
              area(surface, center.p, theta, center.radius, footprint,
                   lower_bound, shade, stats);
          return Pixel(uint16_t(hs::clamp(result.r, 0.0f, 65535.0f)),
                       uint16_t(hs::clamp(result.g, 0.0f, 65535.0f)),
                       uint16_t(hs::clamp(result.b, 0.0f, 65535.0f)));
        };
      });
}
} // namespace Scan
