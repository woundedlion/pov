/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file radial_surface_boundary.h
 * @brief Candidate first-arrival boundaries of similar-polygon families. */

#include <algorithm>
#include <cfloat>
#include <cmath>
#include "core/math/3dmath.h"
#include "core/render/scan/radial_surface_math.h"

namespace Scan {

/** @brief Rotationally periodic candidates of a bounded radial sweep. */
class RadialSurfaceBoundaryCandidates {
public:
  /** @brief Planar candidate-curve coordinate. */
  struct Point {
    float x; ///< Horizontal chart coordinate.
    float y; ///< Vertical chart coordinate.
  };

  /** @brief One unreplicated candidate first-arrival segment. */
  struct Candidate {
    Point a;                   ///< First segment endpoint.
    Point b;                   ///< Second segment endpoint.
    float start_angle = 0;     ///< Optional unwrapped first argument.
    float finish_angle = 0;    ///< Optional unwrapped second argument.
    bool angles_known = false; ///< Whether endpoint arguments are supplied.
    bool knot = false;         ///< Whether this is a knot-edge minimum.
    float center_angle =
        0; ///< Angular representative within half an edge sector.
    int origin_interval =
        0; ///< Source interval, or knot index for knot boundaries.
  };

  /** @brief Intersects an interval with a linear nonpositive half-line.
   * @param value,slope Coefficients of value+slope*u <= 0.
   * @param low,high Parameter interval to restrict.
   * @return Whether any parameter remains. */
  static HS_HOT_INLINE bool clip_parameter(float value, float slope, float &low,
                                           float &high) {
    if (slope > 0.0f)
      high = fminf(high, -value / slope);
    else if (slope < 0.0f)
      low = fmaxf(low, -value / slope);
    else if (value > 0.0f)
      return false;
    return low <= high;
  }

  /** @brief Restricts the outer-tip first-arrival range.
   * @param spin,radial Similarity motion determinant and radial product.
   * @param length2 Squared similarity step length.
   * @param h Unit edge determinant.
   * @param before_x,after_x Adjacent unit edge horizontal components.
   * @param slack Additional containment derivative tolerance.
   * @param low,high Parameter interval to restrict.
   * @return Whether any candidate range remains. */
  static HS_HOT_INLINE bool tip_range(float spin, float radial, float length2,
                                      float h, float before_x, float after_x,
                                      float slack, float &low, float &high) {
    return paired_range({spin * before_x - h * radial, -h * length2},
                        {spin * after_x - h * radial, -h * length2}, slack, low,
                        high);
  }

  /** @brief Restricts a knot-edge first-arrival minimum.
   * @param left_spin,left_radial Incoming similarity motion products.
   * @param right_spin,right_radial Outgoing similarity motion products.
   * @param shape_x,shape_y Unit edge vertex product and squared edge length.
   * @param h Unit edge determinant.
   * @param slack Additional containment derivative tolerance.
   * @param low,high Edge interval to restrict.
   * @return Whether any candidate range remains. */
  static HS_HOT_INLINE bool knot_range(float left_spin, float left_radial,
                                       float right_spin, float right_radial,
                                       float shape_x, float shape_y, float h,
                                       float slack, float &low, float &high) {
    return paired_range(
        {left_spin * shape_x - h * left_radial, left_spin * shape_y},
        {right_spin * shape_x - h * right_radial, right_spin * shape_y}, slack,
        low, high);
  }

  /** @brief Restricts a terminal edge's first-arrival range.
   * @param spin,radial Incoming similarity motion products.
   * @param shape_x,shape_y Unit edge vertex product and squared edge length.
   * @param h Unit edge determinant.
   * @param slack Additional containment derivative tolerance.
   * @param low,high Edge interval to restrict.
   * @return Whether any candidate range remains. */
  static HS_HOT_INLINE bool endpoint_range(float spin, float radial,
                                           float shape_x, float shape_y,
                                           float h, float slack, float &low,
                                           float &high) {
    const float value = spin * shape_x - h * radial, slope = spin * shape_y;
    const float tolerance =
        64 * FLT_EPSILON * (fabsf(value) + fabsf(slope) + 1e-8f) + slack;
    return clip_parameter(value - tolerance, slope, low, high);
  }

  /** @brief Visits first-arrival candidates in one repeated star period.
   * @tparam Surface Prepared similar-star surface type.
   * @tparam Visitor Callable receiving a Candidate.
   * @param surface Similarity knots and unit-star geometry.
   * @param visit Consumer for each unsplit candidate segment.
   * @param derivative_motion Uniform bound on containment derivative changes.
   * @param relative_derivative_motion Additional bound per local squared radius.
   * @details At most 3*N+2 candidates are emitted for N radial intervals.
   * Rotate segments by integer star periods to cover the full chart. */
  template <class Surface, class Visitor>
  static HS_HOT_FLASH_MEMBER void
  visit_candidates(const Surface &surface, Visitor &&visit,
                   float derivative_motion = 0,
                   float relative_derivative_motion = 0) {
    const int count = surface.interval_count();
    const float period = math::TWO_PI_F / (surface.edge_count() / 2);
    for (int i = 0; i < count; ++i) {
      const auto &a = surface.knot(i), &b = surface.knot(i + 1);
      const Point c{a.real, a.imag}, d{b.real - a.real, b.imag - a.imag};
      const float length2 = d.x * d.x + d.y * d.y;
      if (length2 == 0.0f && derivative_motion == 0.0f &&
          relative_derivative_motion == 0.0f)
        continue;
      const auto first = surface.edge(0);
      const float h = first.y;
      const bool positive = radial_cross(c, d) >= 0.0f;
      const auto before = surface.edge(positive ? 0 : surface.edge_count() - 1);
      const auto after = surface.edge(positive ? surface.edge_count() - 1 : 0);
      const float spin = radial_cross(c, d);
      const float radial = c.x * d.x + c.y * d.y;
      const float local_radius2 = fmaxf(a.real * a.real + a.imag * a.imag,
                                        b.real * b.real + b.imag * b.imag);
      float low = 0.0f, high = 1.0f;
      if (!tip_range(spin, radial, length2, h, before.x, after.x,
                     derivative_motion +
                         relative_derivative_motion * local_radius2,
                     low, high))
        continue;
      visit(Candidate{{c.x + low * d.x, c.y + low * d.y},
                      {c.x + high * d.x, c.y + high * d.y},
                      a.phase,
                      b.phase,
                      low == 0.0f && high == 1.0f,
                      false,
                      .5f * (a.phase + b.phase),
                      i});
    }
    for (int endpoint : {0, count})
      for (int edge = 0; edge < 2; ++edge) {
        const auto a = surface.vertex(edge), b = surface.vertex(edge + 1);
        const Point a_value{a.x, a.y}, e{b.x - a.x, b.y - a.y};
        float low = 0, high = 1;
        if (endpoint > 0) {
          const auto &k = surface.knot(endpoint),
                     &prior = surface.knot(endpoint - 1);
          const Point c{k.real, k.imag},
              d{k.real - prior.real, k.imag - prior.imag};
          const float spin = radial_cross(c, d);
          const auto &metrics = surface.edge_metrics(edge);
          const float local_radius2 =
              fmaxf(k.real * k.real + k.imag * k.imag,
                    prior.real * prior.real + prior.imag * prior.imag);
          if (!endpoint_range(spin, c.x * d.x + c.y * d.y, metrics.dot,
                              metrics.squared_length, metrics.determinant,
                              derivative_motion +
                                  relative_derivative_motion * local_radius2,
                              low, high))
            continue;
        }
        visit(Candidate{
            transform(surface.knot(endpoint),
                      {a_value.x + low * e.x, a_value.y + low * e.y}),
            transform(surface.knot(endpoint),
                      {a_value.x + high * e.x, a_value.y + high * e.y}),
            surface.knot(endpoint).phase + edge * period * .5f,
            surface.knot(endpoint).phase + (edge + 1) * period * .5f,
            low == 0 && high == 1, false,
            surface.knot(endpoint).phase + (edge + .5f) * period * .5f,
            endpoint});
      }

    const auto &m0 = surface.edge_metrics(0), &m1 = surface.edge_metrics(1);
    const Point shape0{m0.dot, m0.squared_length};
    const Point shape1{m1.dot, m1.squared_length};
    const float h0 = m0.determinant;
    const float h1 = m1.determinant;
    for (int i = 1; i < count; ++i) {
      const auto &c = surface.knot(i);
      const auto &previous = surface.knot(i - 1), &next = surface.knot(i + 1);
      const Point left{c.real - previous.real, c.imag - previous.imag};
      const Point right{next.real - c.real, next.imag - c.imag};
      const Point c_value{c.real, c.imag};
      const float left_spin = radial_cross(c_value, left),
                  right_spin = radial_cross(c_value, right);
      const float left_radial = c_value.x * left.x + c_value.y * left.y;
      const float right_radial = c_value.x * right.x + c_value.y * right.y;
      for (int edge = 0; edge < 2; ++edge) {
        const float h = edge == 0 ? h0 : h1;
        const Point shape = edge == 0 ? shape0 : shape1;
        const float local_radius2 = fmaxf(
            c_value.x * c_value.x + c_value.y * c_value.y,
            fmaxf(previous.real * previous.real + previous.imag * previous.imag,
                  next.real * next.real + next.imag * next.imag));
        float low = 0.0f, high = 1.0f;
        if (!knot_range(left_spin, left_radial, right_spin, right_radial,
                        shape.x, shape.y, h,
                        derivative_motion +
                            relative_derivative_motion * local_radius2,
                        low, high))
          continue;
        const auto unit_a = surface.vertex(edge), unit_e = surface.edge(edge);
        const Point a{unit_a.x, unit_a.y}, e{unit_e.x, unit_e.y};
        const Point world_a = transform(c, a), world_e = transform(c, e);
        visit(Candidate{
            {world_a.x + low * world_e.x, world_a.y + low * world_e.y},
            {world_a.x + high * world_e.x, world_a.y + high * world_e.y},
            0,
            0,
            false,
            true,
            c.phase + (edge + .5f) * period * .5f,
            i});
      }
    }
  }

private:
  /** @brief Restricts a range between incoming and outgoing half-lines.
   * @param left,right Incoming and outgoing derivative coefficients.
   * @param slack Additional containment derivative tolerance.
   * @param low,high Parameter interval to restrict.
   * @return Whether any candidate range remains. */
  static HS_HOT_INLINE bool paired_range(Point left, Point right, float slack,
                                         float &low, float &high) {
    const float tolerance = 64 * FLT_EPSILON *
                                (fabsf(left.x) + fabsf(left.y) +
                                 fabsf(right.x) + fabsf(right.y) + 1e-8f) +
                            slack;
    return fminf(left.x, left.x + left.y) <= tolerance &&
           fmaxf(right.x, right.x + right.y) >= -tolerance &&
           clip_parameter(left.x - tolerance, left.y, low, high) &&
           clip_parameter(-right.x - tolerance, -right.y, low, high);
  }
  /** @tparam Knot Similarity multiplier type.
   * @param c Similarity multiplier.
   * @param p Unit-star coordinate.
   * @return Transformed coordinate. */
  template <class Knot>
  static HS_HOT_INLINE Point transform(const Knot &c, Point p) {
    return radial_multiply(Point{c.real, c.imag}, p);
  }
};

} // namespace Scan
