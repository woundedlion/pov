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

  /** @brief Visits first-arrival candidates in one repeated star period.
   * @tparam Surface Prepared similar-star surface type.
   * @tparam Visitor Callable receiving a Candidate.
   * @param surface Similarity knots and unit-star geometry.
   * @param visit Consumer for each unsplit candidate segment.
   * @param derivative_motion Uniform bound on containment derivative changes.
   * @details At most 3*N+2 candidates are emitted for N radial intervals.
   * Rotate segments by integer star periods to cover the full chart. */
  template <class Surface, class Visitor>
  static HS_HOT_FLASH_MEMBER void
  visit_candidates_detailed(const Surface &surface, Visitor &&visit,
                            float derivative_motion = 0) {
    const int COUNT = surface.interval_count();
    const float PERIOD = math::TWO_PI_F / (surface.edge_count() / 2);
    for (int i = 0; i < COUNT; ++i) {
      const auto &a = surface.knot(i), &b = surface.knot(i + 1);
      const Point C{a.real, a.imag}, D{b.real - a.real, b.imag - a.imag};
      const float LENGTH2 = D.x * D.x + D.y * D.y;
      if (LENGTH2 == 0.0f && derivative_motion == 0.0f)
        continue;
      const auto first = surface.edge(0);
      const float H = first.y;
      const bool POSITIVE = cross(C, D) >= 0.0f;
      const auto before = surface.edge(POSITIVE ? 0 : surface.edge_count() - 1);
      const auto after = surface.edge(POSITIVE ? surface.edge_count() - 1 : 0);
      const float SPIN = cross(C, D);
      const float RADIAL = C.x * D.x + C.y * D.y;
      const float SLOPE = -H * LENGTH2;
      const Point L{SPIN * before.x - H * RADIAL, SLOPE};
      const Point R{SPIN * after.x - H * RADIAL, SLOPE};
      const float TOLERANCE =
          64.0f * FLT_EPSILON *
              (fabsf(L.x) + fabsf(L.y) + fabsf(R.x) + fabsf(R.y) + 1e-8f) +
          derivative_motion;
      float low = 0.0f, high = 1.0f;
      if (!clip_parameter(L.x - TOLERANCE, L.y, low, high) ||
          !clip_parameter(-R.x - TOLERANCE, -R.y, low, high))
        continue;
      visit(Candidate{{C.x + low * D.x, C.y + low * D.y},
                      {C.x + high * D.x, C.y + high * D.y},
                      a.phase,
                      b.phase,
                      low == 0.0f && high == 1.0f,
                      false,
                      .5f * (a.phase + b.phase),
                      i});
    }
    for (int endpoint : {0, COUNT})
      for (int edge = 0; edge < 2; ++edge) {
        const auto a = surface.vertex(edge), b = surface.vertex(edge + 1);
        const Point A{a.x, a.y}, E{b.x - a.x, b.y - a.y};
        float low = 0, high = 1;
        if (endpoint > 0) {
          const auto &k = surface.knot(endpoint),
                     &prior = surface.knot(endpoint - 1);
          const Point C{k.real, k.imag},
              D{k.real - prior.real, k.imag - prior.imag};
          const float SPIN = cross(C, D), H = cross(A, E);
          const float VALUE =
              SPIN * (A.x * E.x + A.y * E.y) - H * (C.x * D.x + C.y * D.y);
          const float SLOPE = SPIN * (E.x * E.x + E.y * E.y);
          const float TOLERANCE =
              64 * FLT_EPSILON * (fabsf(VALUE) + fabsf(SLOPE) + 1e-8f) +
              derivative_motion;
          if (!clip_parameter(VALUE - TOLERANCE, SLOPE, low, high))
            continue;
        }
        visit(Candidate{
            transform(surface.knot(endpoint),
                      {A.x + low * E.x, A.y + low * E.y}),
            transform(surface.knot(endpoint),
                      {A.x + high * E.x, A.y + high * E.y}),
            surface.knot(endpoint).phase + edge * PERIOD * .5f,
            surface.knot(endpoint).phase + (edge + 1) * PERIOD * .5f,
            low == 0 && high == 1, false,
            surface.knot(endpoint).phase + (edge + .5f) * PERIOD * .5f,
            endpoint});
      }

    const auto V0 = surface.vertex(0), V1 = surface.vertex(1);
    const auto E0 = surface.edge(0), E1 = surface.edge(1);
    const Point SHAPE0{V0.x * E0.x + V0.y * E0.y, E0.x * E0.x + E0.y * E0.y};
    const Point SHAPE1{V1.x * E1.x + V1.y * E1.y, E1.x * E1.x + E1.y * E1.y};
    const float H0 = V0.x * E0.y - V0.y * E0.x;
    const float H1 = V1.x * E1.y - V1.y * E1.x;
    for (int i = 1; i < COUNT; ++i) {
      const auto &c = surface.knot(i);
      const auto &previous = surface.knot(i - 1), &next = surface.knot(i + 1);
      const Point LEFT{c.real - previous.real, c.imag - previous.imag};
      const Point RIGHT{next.real - c.real, next.imag - c.imag};
      const Point C{c.real, c.imag};
      const float LEFT_SPIN = cross(C, LEFT), RIGHT_SPIN = cross(C, RIGHT);
      const float LEFT_RADIAL = C.x * LEFT.x + C.y * LEFT.y;
      const float RIGHT_RADIAL = C.x * RIGHT.x + C.y * RIGHT.y;
      for (int edge = 0; edge < 2; ++edge) {
        const float H = edge == 0 ? H0 : H1;
        const Point SHAPE = edge == 0 ? SHAPE0 : SHAPE1;
        const Point L{LEFT_SPIN * SHAPE.x - H * LEFT_RADIAL,
                      LEFT_SPIN * SHAPE.y};
        const Point R{RIGHT_SPIN * SHAPE.x - H * RIGHT_RADIAL,
                      RIGHT_SPIN * SHAPE.y};
        const float TOLERANCE =
            64.0f * FLT_EPSILON *
                (fabsf(L.x) + fabsf(L.y) + fabsf(R.x) + fabsf(R.y) + 1e-8f) +
            derivative_motion;
        float low = 0.0f, high = 1.0f;
        if (!clip_parameter(L.x - TOLERANCE, L.y, low, high) ||
            !clip_parameter(-R.x - TOLERANCE, -R.y, low, high))
          continue;
        const auto unit_a = surface.vertex(edge), unit_e = surface.edge(edge);
        const Point A{unit_a.x, unit_a.y}, E{unit_e.x, unit_e.y};
        const Point WORLD_A = transform(c, A), WORLD_E = transform(c, E);
        visit(Candidate{
            {WORLD_A.x + low * WORLD_E.x, WORLD_A.y + low * WORLD_E.y},
            {WORLD_A.x + high * WORLD_E.x, WORLD_A.y + high * WORLD_E.y},
            0,
            0,
            false,
            true,
            c.phase + (edge + .5f) * PERIOD * .5f,
            i});
      }
    }
  }

private:
  /** @param a,b Chart vectors.
   * @return Signed determinant. */
  static HS_HOT_INLINE float cross(Point a, Point b) {
    return a.x * b.y - a.y * b.x;
  }
  /** @param a,b Complex chart vectors.
   * @return Complex product. */
  static HS_HOT_INLINE Point multiply(Point a, Point b) {
    return {a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y};
  }
  /** @tparam Knot Similarity multiplier type.
   * @param c Similarity multiplier.
   * @param p Unit-star coordinate.
   * @return Transformed coordinate. */
  template <class Knot>
  static HS_HOT_INLINE Point transform(const Knot &c, Point p) {
    return multiply({c.real, c.imag}, p);
  }
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
};

} // namespace Scan
