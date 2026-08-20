/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <utility>
#include <type_traits>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <array>
#include <concepts>
#include "math/geometry.h"
#include "render/shading.h"
#include "color/color.h"
#include "platform/constants.h"
#include "render/clip.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include "containers/triangular_bitset.h"

/**
 * @file cull.h
 * @brief Edge samplers and the screen row/column span and clip-cull kernel the
 * curve rasterizer walks each segment through.
 */

namespace Plot {

/**
 * @brief Geodesic segment shorter than this (radians) collapses to a point.
 * @details 100× math::EPS_GEOMETRIC (1e-3 vs 1e-5): a slerp-axis stability
 * bound that picks the interpolation strategy, not a positional near-equality
 * test, so it does not track math::EPS_GEOMETRIC.
 */
inline constexpr float EPS_GEODESIC_SEGMENT = 0.001f;

/**
 * @brief Minimum |cross(a, b)|² for which the arc pole of a geodesic edge is
 *        taken from the cross product rather than a stable perpendicular.
 * @details The bound is on the quantity the pole normalization consumes rather
 * than on angle_between, whose derivative diverges as the normalized dot
 * approaches ±1: one ULP there moves the reported angle by ~3.5e-4 rad, so no
 * angular band narrow enough to be useful can also be wide enough to hold, and
 * an antipodal edge reaches a cross product it cannot normalize. |cross| =
 * sin(angle), so 1e-8 names the same geometric band an angular 1e-4 does,
 * without the amplification. Above it the cross components carry ~1e-7 of
 * absolute rounding, bounding the pole's direction error at ~2e-3 rad — a tenth
 * of a pixel at W=288.
 */
inline constexpr float EPS_ARC_POLE_SQ = 1e-8f;

/**
 * @brief Minimum |axis.y| for which a geodesic edge's endpoint columns bound
 *        its azimuth span.
 * @details Below this the great circle runs near the poles, where longitude is
 * ill-conditioned and the interior leaves the endpoint columns.
 */
inline constexpr float AXIS_Y_EPS = 1e-4f;

/**
 * @brief Minimum worst-case sin(φ) over a curve for which its plotted columns
 *        are trusted enough to cull by.
 * @details The azimuth Lipschitz bound scales as 1/sin(φ), and nearer the poles
 * the plotted column is float noise; below this the column cull bails.
 */
inline constexpr float MIN_SIN_PHI = 0.05f;

/**
 * @brief Floor on the adaptive sub-step length, as a fraction of base_step.
 * @details Caps sub-steps per segment so polar curves don't oversample: the
 * screen-velocity step sampler (screen_step) drives the step toward zero where
 * the azimuthal velocity diverges at the poles, and this is the lower clamp that
 * bounds it. A clamp, not a tolerance.
 */
inline constexpr float MIN_POLE_SCALE = 0.05f;

/**
 * @brief Target screen-space spacing (pixels) between adaptive sub-samples.
 * @details The rasterizer sizes each sub-step so consecutive samples land about
 * this far apart in SCREEN space. Slightly sub-pixel so the bilinear AntiAlias
 * splat of neighbouring samples overlaps and the rendered curve has no holes;
 * smaller = denser = smoother but costlier.
 */
inline constexpr float SCREEN_STEP_PX = 0.9f;

/**
 * @brief Columns of slack added on each side of a culled column span.
 * @details Absorbs plot rounding and the AntiAlias tap spread.
 */
inline constexpr int COL_PAD = 2;

/**
 * @brief Columns a padded span reaches past its fractional end.
 * @details The pad plus the boundary column ceil() adds.
 */
inline constexpr int COL_FOOTPRINT = COL_PAD + 1;

/**
 * @brief Columns outside the render band a clip cut is placed at.
 * @details A piece ending exactly on the band edge still overlaps it once
 * finish_col_span widens the span by COL_FOOTPRINT, so a cut there would leave
 * the outside piece visible and buy nothing. One column past that footprint
 * also absorbs the fast-trig error in the cut.
 */
inline constexpr int CLIP_CUT_COL_PAD = COL_FOOTPRINT + 1;

/**
 * @brief Rows outside the render band a clip cut is placed at.
 * @details could_intersect_y takes the band edge itself as an intersection, and
 * the row map's fast-trig round trip moves a cut by a fraction of a row.
 */
inline constexpr int CLIP_CUT_ROW_PAD = 1;

// --- Strategy Helpers ---

/**
 * @brief Forward azimuthal-equidistant projection of a sphere point to the plane.
 * @param p Unit sphere point to project.
 * @param basis Projection basis; center is basis.v, axes basis.u/basis.w.
 * @return Plane coordinates whose radius is the great-circle angle from the
 *         basis center (radians) and whose azimuth follows the basis u/w axes.
 * @details Shared by the planar rasterization strategy and the clip-cull
 * arc-extent sampler; azimuthal_unproject is the inverse map.
 */
static inline std::pair<float, float> azimuthal_project(const Vector &p,
                                                        const Basis &basis) {
  float R = angle_between(p, basis.v);
  if (R < math::EPS_GEOMETRIC)
    return {0.0f, 0.0f};
  float theta = fast_atan2(dot(p, basis.w), dot(p, basis.u));
  return {R * fast_cosf(theta), R * fast_sinf(theta)};
}

/**
 * @brief Inverse azimuthal-equidistant projection from the plane to the sphere.
 * @param Px Plane x-coordinate (azimuthal-equidistant).
 * @param Py Plane y-coordinate (azimuthal-equidistant).
 * @param basis Projection basis; center is basis.v, axes basis.u/basis.w.
 * @return Unit sphere point at great-circle angle sqrt(Px²+Py²) from basis.v.
 */
static inline Vector azimuthal_unproject(float Px, float Py,
                                         const Basis &basis) {
  HS_PLOT_COUNT(planar_unprojects);
  float R = sqrtf(Px * Px + Py * Py);
  if (R < math::EPS_GEOMETRIC)
    return basis.v;
  float sin_r;
  float cos_r;
  if (R <= PI_F) {
    fast_sincosf_0_pi(R, sin_r, cos_r);
  } else {
    sin_r = fast_sinf(R);
    cos_r = fast_cosf(R);
  }
  const Vector radial = (basis.u * Px) + (basis.w * Py);
  return (basis.v * cos_r) + (radial * (sin_r / R));
}

/**
 * @brief A rasterized sample: its unit-sphere position and unit tangent.
 * @details `tan` is the curve's unit tangent with respect to ARC LENGTH at the
 * sample, used by screen_step() to size the next sub-step. Zero for a degenerate
 * edge, where screen_step's speed floor maps it to a base_step (one-dot) step.
 */
struct SamplePT {
  Vector pos;
  Vector tan;
};

/**
 * @brief Renormalizes a sampled position with one Newton step.
 * @param v A sampler position, unit up to the fast sin/cos residual.
 * @return v scaled to unit length to second order.
 * @details The fast sin/cos kernels leave a sampled position up to 2e-3 off
 * unit, and vector_to_pixel's phi = acos(v.y) turns that into a near-pole row
 * offset. One Newton step leaves 5e-6 without the exact normalize's divides.
 */
static inline Vector newton_unit(const Vector &v) {
  const float norm2 = dot(v, v);
  return v * (1.5f - 0.5f * norm2);
}

constexpr int PLANAR_LEN_SAMPLES = 4;

/**
 * @brief Cumulative on-sphere arc length at PLANAR_LEN_SAMPLES+1 evenly-spaced
 *        PROJECTION samples of the azimuthal-equidistant straight edge whose
 *        projection starts at `proj` and spans (dx, dy).
 * @details arc_cumul[0] = 0; arc_cumul.back() is the PLANAR_LEN_SAMPLES-chord sum,
 * an underestimate of the rendered length (a few percent on a bowed edge), and
 * below the endpoints' angle_between on a radial edge, which does not bow.
 * Shared by the planar rasterizer (which inverts the table for arc-uniform
 * stepping) and rasterize()'s perimeter pre-pass (which takes the total), so
 * both sample identical points and sum identical lengths.
 */
static inline void
planar_arc_cumul(const std::pair<float, float> &proj, float dx, float dy,
                 const Basis &planar_basis,
                 std::array<float, PLANAR_LEN_SAMPLES + 1> &arc_cumul) {
  HS_PLOT_ADD(planar_arc_samples, PLANAR_LEN_SAMPLES + 1);
  arc_cumul[0] = 0.0f;
  Vector prev = azimuthal_unproject(proj.first, proj.second, planar_basis);
  for (int k = 1; k <= PLANAR_LEN_SAMPLES; ++k) {
    float p = static_cast<float>(k) / PLANAR_LEN_SAMPLES;
    Vector cur = azimuthal_unproject(proj.first + dx * p, proj.second + dy * p,
                                     planar_basis);
    arc_cumul[k] = arc_cumul[k - 1] + angle_between(prev, cur);
    prev = cur;
  }
}

/**
 * @brief Arc-uniform sampler for one azimuthal-equidistant straight edge.
 * @details Every edge sampler exposes both entry points the rasterizer needs:
 * `pos(t)` for the drawing phase, which wants the position alone, and
 * `operator()(t)` for the simulation phase, which also needs the tangent.
 */
struct PlanarEdgeSampler {
  std::pair<float, float> proj1; /**< Projection of the edge start. */
  float dx;                      /**< Projected chord x-component. */
  float dy;                      /**< Projected chord y-component. */
  const Basis *basis;            /**< Azimuthal-equidistant projection basis. */
  Vector chart_tangent;          /**< Constant chart-space edge tangent. */
  /** Cumulative on-sphere arc at evenly-spaced PROJECTION samples. */
  std::array<float, PLANAR_LEN_SAMPLES + 1> arc_cumul;
  float dist; /**< The edge's on-sphere length (radians). */

  /** @brief Maps normalized arc distance to normalized chart distance. */
  float projection_fraction(float s) const {
    if (dist < math::EPS_GEOMETRIC)
      return s;
    const float target = s * dist;
    int k = 0;
    while (k < PLANAR_LEN_SAMPLES - 1 && arc_cumul[k + 1] < target)
      ++k;
    const float seg = arc_cumul[k + 1] - arc_cumul[k];
    const float frac =
        seg > math::EPS_GEOMETRIC ? (target - arc_cumul[k]) / seg : 0.0f;
    return std::min(1.0f, std::max(0.0f, (static_cast<float>(k) + frac) /
                                             PLANAR_LEN_SAMPLES));
  }

  /** @brief Maps increasing arc fractions without rescanning prior intervals. */
  float projection_fraction_monotonic(float s, int &interval) const {
    if (dist < math::EPS_GEOMETRIC)
      return s;
    const float target = s * dist;
    while (interval < PLANAR_LEN_SAMPLES - 1 &&
           arc_cumul[interval + 1] < target)
      ++interval;
    const float seg = arc_cumul[interval + 1] - arc_cumul[interval];
    const float frac =
        seg > math::EPS_GEOMETRIC ? (target - arc_cumul[interval]) / seg : 0.0f;
    return std::min(1.0f, std::max(0.0f, (static_cast<float>(interval) + frac) /
                                             PLANAR_LEN_SAMPLES));
  }

  /**
   * @brief Unprojects the chart line at PROJECTION fraction p in [0,1].
   * @details Projection-uniform, so not arc-uniform under the anisotropic
   * metric; pos() maps an arc fraction onto it.
   */
  Vector unproject(float p) const {
    return azimuthal_unproject(proj1.first + dx * p, proj1.second + dy * p,
                               *basis);
  }

  /**
   * @brief Unprojects the chart line at PROJECTION fraction p, carrying the
   *        analytic tangent when `WithTangent`.
   * @tparam WithTangent Also derive and normalize the tangent; the position-only
   *         instantiation discards the rate terms before they are computed.
   * @details always_inline so the position-only callers pay for neither the
   * discarded tangent nor a call.
   */
  template <bool WithTangent>
  __attribute__((always_inline)) SamplePT sample_at(float p) const {
    const float x = proj1.first + dx * p;
    const float y = proj1.second + dy * p;
    const float r2 = x * x + y * y;
    if (r2 < math::EPS_GEOMETRIC * math::EPS_GEOMETRIC) {
      if constexpr (!WithTangent)
        return {basis->v, Vector()};
      HS_PLOT_COUNT(normalizations);
      return {basis->v, normalized_or(chart_tangent, Vector())};
    }

    const float radius = sqrtf(r2);
    const float inv_radius = 1.0f / radius;
    float sin_radius;
    float cos_radius;
    if (radius <= PI_F) {
      fast_sincosf_0_pi(radius, sin_radius, cos_radius);
    } else {
      sin_radius = fast_sinf(radius);
      cos_radius = fast_cosf(radius);
    }
    const Vector radial = (basis->u * x) + (basis->w * y);
    const float radial_scale = sin_radius * inv_radius;
    const Vector position = (basis->v * cos_radius) + (radial * radial_scale);
    if constexpr (!WithTangent)
      return {position, Vector()};

    const float radius_rate = (x * dx + y * dy) * inv_radius;
    const float scale_rate =
        (radius * cos_radius - sin_radius) * inv_radius * inv_radius;
    const Vector tangent = (basis->v * (-sin_radius * radius_rate)) +
                           (chart_tangent * radial_scale) +
                           (radial * (scale_rate * radius_rate));
    HS_PLOT_COUNT(normalizations);
    return {position, normalized_or(tangent, Vector())};
  }

  /**
   * @brief Position at arc fraction s in [0,1].
   * @details Inverts the piecewise-linear cumulative-arc table to a projection
   * parameter, then unprojects. A short scan over PLANAR_LEN_SAMPLES floats —
   * no trig.
   */
  Vector pos(float s) const { return unproject(projection_fraction(s)); }

  /** @brief Evaluates position and analytic tangent without a second unproject. */
  HS_FLASH_MEMBER SamplePT one_pass(float s) const {
    return sample_at<true>(projection_fraction(s));
  }

  /** @brief Evaluates an increasing sequence without rescanning arc intervals. */
  SamplePT one_pass_monotonic(float s, int &interval) const {
    return sample_at<true>(projection_fraction_monotonic(s, interval));
  }

  /** @brief Evaluates only position for an increasing sample sequence. */
  Vector position_monotonic(float s, int &interval) const {
    return sample_at<false>(projection_fraction_monotonic(s, interval)).pos;
  }

  /**
   * @brief Position and unit tangent at arc fraction s in [0,1].
   * @details Feeds the same screen-velocity sub-step sampler as the geodesic
   * path.
   */
  SamplePT operator()(float s) const { return one_pass(s); }
};

/** @brief Builds the reusable arc sampler for one planar edge. */
static inline PlanarEdgeSampler
make_planar_edge_sampler(const Vector &a, const Vector &b,
                         const Basis &planar_basis) {
  PlanarEdgeSampler sampler;
  sampler.proj1 = azimuthal_project(a, planar_basis);
  auto proj2 = azimuthal_project(b, planar_basis);
  sampler.dx = proj2.first - sampler.proj1.first;
  sampler.dy = proj2.second - sampler.proj1.second;
  sampler.basis = &planar_basis;
  sampler.chart_tangent =
      (planar_basis.u * sampler.dx) + (planar_basis.w * sampler.dy);
  planar_arc_cumul(sampler.proj1, sampler.dx, sampler.dy, planar_basis,
                   sampler.arc_cumul);
  sampler.dist = sampler.arc_cumul[PLANAR_LEN_SAMPLES];
  return sampler;
}

/**
 * @brief Planar interpolation strategy: builds an arc-uniform sampler for one edge.
 * @tparam ProcessSegmentFn Callable (sample, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param planar_basis Azimuthal-equidistant projection basis.
 * @param is_last_segment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length sampler, endpoints, on-sphere
 *                        length (radians), and the last-segment flag.
 * @details The path is a straight line in the azimuthal-equidistant projection.
 * Projection-uniform stepping is NOT arc-uniform under the anisotropic metric,
 * so a short cumulative-arc table is inverted to turn an arc-length fraction into
 * a projection parameter, making planar sampling arc-uniform with no new trig.
 */
template <typename ProcessSegmentFn>
static void
rasterize_planar_strategy(const Fragment &curr, const Fragment &next,
                          const Basis &planar_basis, bool is_last_segment,
                          ProcessSegmentFn &&process_segment) {
  PlanarEdgeSampler sampler =
      make_planar_edge_sampler(curr.pos, next.pos, planar_basis);

  process_segment(sampler, curr, next, sampler.dist, is_last_segment);
}

/**
 * @brief Chord-sum estimate (radians) of the on-sphere length of the
 *        azimuthal-equidistant straight edge a->b, the path planar
 *        interpolation actually renders.
 * @details Shares planar_arc_cumul with rasterize_planar_strategy, so
 * rasterize()'s perimeter pre-pass and per-segment arc accumulator sum exactly
 * the lengths the draw phase walks — the guarantee v1 relies on, not absolute
 * accuracy. An inscribed chord sum: short of the bowed arc it estimates, and
 * below angle_between(a, b) on a radial edge, whose chart line does not bow.
 */
static inline float planar_arc_length(const Vector &a, const Vector &b,
                                      const Basis &planar_basis) {
  auto p1 = azimuthal_project(a, planar_basis);
  auto p2 = azimuthal_project(b, planar_basis);
  std::array<float, PLANAR_LEN_SAMPLES + 1> arc_cumul;
  planar_arc_cumul(p1, p2.first - p1.first, p2.second - p1.second, planar_basis,
                   arc_cumul);
  return arc_cumul[PLANAR_LEN_SAMPLES];
}

/**
 * @brief Unit axis perpendicular to v, stable for an antipodal geodesic.
 * @param v Unit endpoint of a near-antipodal great-circle segment.
 * @return A unit axis perpendicular to v.
 * @details Antipodal endpoints leave cross(v1, v2) ~= 0 (infinitely many
 * connecting geodesics), so normalizing it yields a garbage axis. Cross v with
 * whichever world axis is least parallel to it and normalize for a well-defined
 * rotation axis.
 */
static inline Vector stable_perpendicular_axis(const Vector &v) {
  HS_PLOT_COUNT(normalizations);
  return cross(v, least_parallel_axis(v)).normalized();
}

/**
 * @brief Azimuthal-equidistant projection chart centered on a pole.
 * @param center Unit pole the planar chart is centered on (the 'v' axis).
 * @return A Basis {u, center, w} with u, w spanning the chart plane.
 */
static inline Basis planar_chart_basis(const Vector &center) {
  Vector ref = least_parallel_axis(center);
  HS_PLOT_ADD(normalizations, 2);
  Vector u = cross(center, ref).normalized();
  Vector w = cross(center, u).normalized();
  return {u, center, w};
}

/** @brief Constant sampler for a coincident-endpoint edge. */
struct DegenerateEdgeSampler {
  Vector p; /**< The collapsed edge's single position. */

  Vector pos(float) const { return p; }
  SamplePT operator()(float) const { return {p, Vector()}; }
};

/**
 * @brief Great-circle sampler for one edge.
 * @details v1 and v_perp are orthonormal, so pos is a unit slerp of the two and
 * tan = d(pos)/d(ang) is a unit vector out of the same sin/cos — the
 * screen-velocity sampler's tangent costs no extra trig.
 */
struct GeodesicEdgeSampler {
  /** @brief pos() is unit up to one newton_unit() correction. */
  static constexpr bool NEWTON_UNIT = true;

  Vector v1;        /**< Edge start (unit). */
  Vector v_perp;    /**< Unit vector perpendicular to v1 in the arc plane. */
  float total_dist; /**< The edge's on-sphere length (radians). */

  /** @brief Position at arc fraction t in [0,1]. */
  Vector pos(float t) const {
    float s, c;
    fast_sincosf_0_pi(total_dist * t, s, c);
    return (v1 * c) + (v_perp * s);
  }

  /** @brief Position and unit tangent at arc fraction t in [0,1]. */
  SamplePT operator()(float t) const {
    float s, c;
    fast_sincosf_0_pi(total_dist * t, s, c);
    return {(v1 * c) + (v_perp * s), (v_perp * c) - (v1 * s)};
  }
};

/**
 * @brief Shared per-edge geodesic setup: arc length and slerp axis.
 * @details have_axis is false on an edge the renderer collapses to a dot, which
 *          has no arc to pole.
 */
struct GeodesicEdgeSpan {
  float total;    /**< angle_between(a, b) in radians. */
  bool antipodal; /**< axis came from stable_perpendicular_axis, not cross. */
  bool have_axis; /**< axis holds a unit arc pole. */
  Vector axis;    /**< Unit arc pole (valid iff have_axis). */
};

/**
 * @brief Computes the shared geodesic edge setup once per edge.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @details The single slerp-axis decision every geodesic path builds on: the
 * renderer's sampler, Line::sample's presampled polyline, and the row/column
 * span bounds that cull them all resolve the same axis from it. always_inline
 * so the per-edge setup costs no call on the rasterizer's hot path.
 */
static __attribute__((always_inline)) inline GeodesicEdgeSpan
make_geodesic_edge_span(const Vector &a, const Vector &b) {
  GeodesicEdgeSpan es;
  es.total = angle_between(a, b);
  if (es.total < EPS_GEODESIC_SEGMENT) {
    es.antipodal = false;
    es.axis = Vector(0.0f, 0.0f, 0.0f);
    es.have_axis = false;
    return es;
  }
  Vector pole = cross(a, b);
  float pole_len_sq = dot(pole, pole);
  es.antipodal = pole_len_sq < EPS_ARC_POLE_SQ;
  if (es.antipodal) {
    es.axis = stable_perpendicular_axis(a);
  } else {
    HS_PLOT_COUNT(normalizations);
    es.axis = pole * (1.0f / sqrtf(pole_len_sq));
  }
  es.have_axis = true;
  return es;
}

/**
 * @brief Geodesic interpolation strategy: builds a great-circle sampler for one edge.
 * @tparam ProcessSegmentFn Callable (sample, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param is_last_segment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length sampler, endpoints, on-sphere
 *                        length (radians), and the last-segment flag.
 * @details Slerps about the axis make_geodesic_edge_span resolves; a
 * coincident-endpoint edge collapses to a constant sampler.
 */
HS_O3_BEGIN
template <typename ProcessSegmentFn>
static void rasterize_geodesic_strategy(const Fragment &curr,
                                        const Fragment &next,
                                        bool is_last_segment,
                                        ProcessSegmentFn &&process_segment) {
  HS_MSP_STALL_START(edge_setup_start);
  Vector v1 = curr.pos;
  const GeodesicEdgeSpan es = make_geodesic_edge_span(v1, next.pos);

  if (!es.have_axis) {
    HS_PLOT_COUNT(degenerate);
    HS_MSP_STALL_STOP(edge_setup, edge_setup_start);
    process_segment(DegenerateEdgeSampler{v1}, curr, next, es.total,
                    is_last_segment);
  } else {
    const GeodesicEdgeSampler sampler{v1, cross(es.axis, v1), es.total};
    HS_MSP_STALL_STOP(edge_setup, edge_setup_start);
    process_segment(sampler, curr, next, es.total, is_last_segment);
  }
}
HS_O3_END

constexpr int PLANAR_SPAN_SAMPLES = 8;

// make_planar_edge_sampler reads the arc-length table straight out of the span's
// interior samples at stride 2, so the span count must be exactly twice the
// arc-table count.
static_assert(PLANAR_SPAN_SAMPLES == 2 * PLANAR_LEN_SAMPLES);

/**
 * @brief Shared per-edge planar setup for the row/column span bounds.
 * @details Projects the edge and samples its chart line through the
 *          rasterizer's unprojection map once per edge. interior holds the
 *          k/PLANAR_SPAN_SAMPLES samples for k in [1, PLANAR_SPAN_SAMPLES);
 *          the endpoint sample is deferred to planar_col_span (the row span
 *          never needs it).
 */
struct PlanarEdgeSpan {
  std::pair<float, float> p1; /**< Projection of the edge start. */
  float dX;                   /**< Projected chord x-component. */
  float dY;                   /**< Projected chord y-component. */
  float gap_arc;              /**< Bound on each inter-sample arc length. */
  std::array<Vector, PLANAR_SPAN_SAMPLES - 1> interior;
};

/**
 * @brief Computes the shared planar edge setup once per edge.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param planar_basis Azimuthal-equidistant projection basis.
 */
static inline PlanarEdgeSpan make_planar_edge_span(const Vector &a,
                                                   const Vector &b,
                                                   const Basis &planar_basis) {
  PlanarEdgeSpan es;
  es.p1 = azimuthal_project(a, planar_basis);
  auto p2 = azimuthal_project(b, planar_basis);
  es.dX = p2.first - es.p1.first;
  es.dY = p2.second - es.p1.second;
  // The projected chord over-estimates the on-sphere arc, so gap_arc bounds
  // each inter-sample arc length.
  es.gap_arc = sqrtf(es.dX * es.dX + es.dY * es.dY) / PLANAR_SPAN_SAMPLES;
  for (int k = 1; k < PLANAR_SPAN_SAMPLES; ++k) {
    float p = static_cast<float>(k) / PLANAR_SPAN_SAMPLES;
    es.interior[k - 1] = azimuthal_unproject(
        es.p1.first + es.dX * p, es.p1.second + es.dY * p, planar_basis);
  }
  return es;
}

/**
 * @brief Builds a planar sampler from the exact quarter-point cull samples.
 * @param span Planar cull setup whose interior samples cover eighth points.
 * @param end Unprojected edge endpoint from planar_col_span().
 * @param planar_basis Azimuthal-equidistant projection basis.
 */
static inline PlanarEdgeSampler
make_planar_edge_sampler(const PlanarEdgeSpan &span, const Vector &end,
                         const Basis &planar_basis) {
  PlanarEdgeSampler sampler;
  sampler.proj1 = span.p1;
  sampler.dx = span.dX;
  sampler.dy = span.dY;
  sampler.basis = &planar_basis;
  sampler.chart_tangent =
      (planar_basis.u * sampler.dx) + (planar_basis.w * sampler.dy);
  HS_PLOT_ADD(planar_arc_samples, PLANAR_LEN_SAMPLES + 1);
  sampler.arc_cumul[0] = 0.0f;
  Vector prev =
      azimuthal_unproject(span.p1.first, span.p1.second, planar_basis);
  for (int k = 1; k < PLANAR_LEN_SAMPLES; ++k) {
    const Vector &cur = span.interior[k * 2 - 1];
    sampler.arc_cumul[k] = sampler.arc_cumul[k - 1] + angle_between(prev, cur);
    prev = cur;
  }
  sampler.arc_cumul[PLANAR_LEN_SAMPLES] =
      sampler.arc_cumul[PLANAR_LEN_SAMPLES - 1] + angle_between(prev, end);
  sampler.dist = sampler.arc_cumul[PLANAR_LEN_SAMPLES];
  return sampler;
}

/**
 * @brief Screen row of a unit-sphere y coordinate (the renderer's row map).
 * @tparam H Rasterization height (pixel grid).
 * @param y Unit-sphere y in [-1, 1] (clamped).
 */
template <int H> static inline float y_to_screen_row(float y) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  return phi_to_y_virtual(fast_acos(hs::clamp(y, -1.0f, 1.0f)), H_VIRT);
}

/**
 * @brief Geodesic screen-row span from precomputed endpoint rows.
 * @tparam H Rasterization height (pixel grid).
 * @param ra Precomputed y_to_screen_row<H>(a.y).
 * @param rb Precomputed y_to_screen_row<H>(b.y).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 * @details The arc y(t) has a turning point inside the span iff the forward
 * tangent's y-component flips sign between the endpoints; the extremal |y| is
 * the great circle's peak latitude sqrt(1 - n.y²) (n = arc pole). The span is
 * the exact closed-form y range, so no one-row epsilon. A degenerate setup
 * (no axis) keeps the endpoint rows.
 */
template <int H>
static __attribute__((always_inline)) inline void
geodesic_row_span_rows(float ra, float rb, const Vector &a, const Vector &b,
                       const GeodesicEdgeSpan &es, float &row_lo,
                       float &row_hi) {
  row_lo = std::min(ra, rb);
  row_hi = std::max(ra, rb);
  if (!es.have_axis)
    return;
  float t0 = cross(es.axis, a).y; // forward tangent y at a
  float t1 = cross(es.axis, b).y; // forward tangent y at b
  if ((t0 > 0.0f) != (t1 > 0.0f)) {
    // std::max(0, ...) absorbs the tiny negative that fast-math
    // renormalization of the axis can produce when |axis.y| ≈ 1 (a
    // near-polar arc pole), keeping the sqrt domain-safe.
    float peak = sqrtf(std::max(0.0f, 1.0f - es.axis.y * es.axis.y));
    float rp = y_to_screen_row<H>(t0 > 0.0f ? peak : -peak);
    row_lo = std::min(row_lo, rp);
    row_hi = std::max(row_hi, rp);
  }
}

/**
 * @brief Geodesic screen-row span from a precomputed edge setup.
 * @tparam H Rasterization height (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 */
template <int H>
static inline void geodesic_row_span(const Vector &a, const Vector &b,
                                     const GeodesicEdgeSpan &es, float &row_lo,
                                     float &row_hi) {
  geodesic_row_span_rows<H>(y_to_screen_row<H>(a.y), y_to_screen_row<H>(b.y), a,
                            b, es, row_lo, row_hi);
}

/**
 * @brief Planar screen-row span from a precomputed edge setup.
 * @tparam H Rasterization height (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_planar_edge_span(a, b, basis).
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 * @details No closed-form latitude extremum exists, so the endpoint rows are
 * extended over the shared samples and widened by the arc's Lipschitz bound.
 * The cull and renderer do NOT take bit-identical samples, so gap-freeness
 * comes from the Lipschitz + one-row margin: phi is 1-Lipschitz in angular
 * distance, so between samples |Δrow| ≤ (Δarc)·(H_VIRT−1)/π. The one-row
 * epsilon absorbs the sub-pixel difference between the unprojected sample
 * (≈unit to fast-math precision) and the renderer's normalized plot position.
 */
template <int H>
static inline void planar_row_span(const Vector &a, const Vector &b,
                                   const PlanarEdgeSpan &es, float &row_lo,
                                   float &row_hi) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  float ra = y_to_screen_row<H>(a.y);
  float rb = y_to_screen_row<H>(b.y);
  row_lo = std::min(ra, rb);
  row_hi = std::max(ra, rb);
  for (const Vector &s : es.interior) {
    float r = y_to_screen_row<H>(s.y);
    row_lo = std::min(row_lo, r);
    row_hi = std::max(row_hi, r);
  }
  float margin = es.gap_arc * (static_cast<float>(H_VIRT - 1) / PI_F) + 1.0f;
  row_lo -= margin;
  row_hi += margin;
}

/**
 * @brief Pads a fractional column interval and wraps it into a [0, W) arc.
 * @tparam W Rasterization width (pixel grid).
 * @param s_f Fractional start column.
 * @param len_f Fractional arc length in columns.
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 */
template <int W>
static __attribute__((always_inline)) inline void
finish_col_span(float s_f, float len_f, int &col_s, int &col_len) {
  const int lo = static_cast<int>(floorf(s_f)) - COL_PAD;
  const int hi = static_cast<int>(ceilf(s_f + len_f)) + COL_PAD;
  col_len = std::min(hi - lo + 1, W);
  col_s = ((lo % W) + W) % W;
}

/**
 * @brief Pads a column span whose start is already within one period.
 */
template <int W>
static inline void finish_col_span_one_period(float start, float length,
                                              int &col_s, int &col_len) {
  const int lo = static_cast<int>(floorf(start)) - COL_PAD;
  const int hi = static_cast<int>(ceilf(start + length)) + COL_PAD;
  col_len = std::min(hi - lo + 1, W);
  col_s = lo < 0 ? lo + W : lo;
}

/**
 * @brief Wraps a column delta into [0, W) with a single conditional add.
 * @tparam W Rasterization width (pixel grid).
 * @param d Delta to wrap; must lie in (-W, W).
 * @return The wrapped delta in [0, W), bit-identical to wrap(d, W) over that
 *   domain.
 * @details `d + W` can round up to exactly W for a tiny negative d; the upper
 * guard folds that back to 0, as wrap's own half-open guard does.
 */
template <int W>
static __attribute__((always_inline)) inline float wrap_one_period(float d) {
  assert(d > -W && d < W);
  if (d < 0.0f)
    d += W;
  return (d >= W) ? 0.0f : d;
}

/**
 * @brief Geodesic screen-column arc from precomputed endpoint columns.
 * @tparam W Rasterization width (pixel grid).
 * @param ca Precomputed vector_to_theta<W>(a).
 * @param cb Precomputed vector_to_theta<W>(b).
 * @param a Edge start (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 * @return False when no useful bound exists (degenerate cross on a
 *         non-collapsed edge, or a near-meridian axis whose y-component is
 *         float noise and the longitude can jump across a pole) — the caller
 *         must skip the horizontal cull.
 * @details Longitude is globally monotone along the rendered circle — with
 * pos(ang) = a·cos + cross(axis, a)·sin, the atan2(z, x) rate numerator
 * pos.x·tan.z - pos.z·tan.x folds to -axis.y, a constant. The arc therefore
 * sweeps from a's column toward the end column in the direction sign(-axis.y),
 * and one full revolution sweeps exactly W, so the directed modular difference
 * is the exact sweep. Antipodal symmetry (λ(-p) = λ(p) + π) makes every
 * half-circle sweep exactly W/2 and shorter arcs less, so the span is always
 * the endpoints' short-way separation — the direction only disambiguates the
 * near-antipodal boundary, where short-way is float noise. The column mapping
 * is the renderer's vector_to_theta.
 */
template <int W>
static inline bool geodesic_col_span_cols(float ca, float cb, const Vector &a,
                                          const GeodesicEdgeSpan &es,
                                          int &col_s, int &col_len) {
  float s_f, len_f;

  if (es.total < EPS_GEODESIC_SEGMENT) {
    // The renderer collapses the edge to a dot at a; span both endpoints the
    // short way around.
    const float d = wrap_one_period<W>(cb - ca);
    if (d <= W * 0.5f) {
      s_f = ca;
      len_f = d;
    } else {
      s_f = cb;
      len_f = W - d;
    }
  } else {
    if (!es.have_axis)
      return false;
    if (std::abs(es.axis.y) < AXIS_Y_EPS)
      return false;

    float ce;
    if (es.antipodal) {
      // The arbitrary-axis half-turn lands near, not on, b; take the column of
      // the point the renderer actually reaches.
      Vector v_perp = cross(es.axis, a);
      Vector end = a * fast_cosf(es.total) + v_perp * fast_sinf(es.total);
      ce = vector_to_theta<W>(end);
    } else {
      ce = cb;
    }

    if (es.axis.y < 0.0f) { // longitude increases from a
      s_f = ca;
      len_f = wrap_one_period<W>(ce - ca);
    } else {
      s_f = ce;
      len_f = wrap_one_period<W>(ca - ce);
    }
  }

  finish_col_span_one_period<W>(s_f, len_f, col_s, col_len);
  return true;
}

/**
 * @brief Geodesic screen-column arc from a precomputed edge setup.
 * @tparam W Rasterization width (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 * @return False when no useful bound exists — the caller must skip the
 *         horizontal cull.
 */
template <int W>
static inline bool geodesic_col_span(const Vector &a, const Vector &b,
                                     const GeodesicEdgeSpan &es, int &col_s,
                                     int &col_len) {
  return geodesic_col_span_cols<W>(vector_to_theta<W>(a), vector_to_theta<W>(b),
                                   a, es, col_s, col_len);
}

/**
 * @brief Most arc fractions a clip band can cut one geodesic edge at.
 * @details Two meridians met once each, two latitude circles met twice each.
 */
inline constexpr int GEODESIC_CLIP_MAX_SPLITS = 6;

/**
 * @brief The clip band's cut boundaries, in the terms the arc solve reads them.
 */
struct ClipCutBounds {
  float col_x[2]; /**< Boundary meridian directions, x component. */
  float col_z[2]; /**< Boundary meridian directions, z component. */
  float row_y[2]; /**< Boundary latitudes, as cos(phi). */
  bool cols;      /**< Column boundaries cut (x clipping is active). */
  bool rows;      /**< Row boundaries cut (band shorter than the canvas). */
};

/**
 * @brief Resolves the clip band's cut boundaries for a whole draw.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @return Boundary geometry to hand every geodesic_clip_splits call under
 *         @p cr.
 * @details The band is fixed for a draw, so its boundary directions and
 * latitudes resolve once rather than per edge.
 */
template <int W, int H>
static inline ClipCutBounds make_clip_cut_bounds(const ClipRegion &cr,
                                                 const ClipRegion::XClip &xc) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();

  ClipCutBounds cb{};
  cb.cols = xc.active;
  cb.rows = cr.render_y_start() > 0 || cr.render_y_end() < cr.h;

  if (cb.cols) {
    const int cols[2] = {xc.rs - CLIP_CUT_COL_PAD, xc.re + CLIP_CUT_COL_PAD};
    for (int i = 0; i < 2; ++i) {
      const int c = ((cols[i] % W) + W) % W;
      cb.col_x[i] = TrigLUT<W, H>::cos_theta(c);
      cb.col_z[i] = TrigLUT<W, H>::sin_theta[c];
    }
  }
  if (cb.rows) {
    const int rows[2] = {cr.render_y_start() - CLIP_CUT_ROW_PAD,
                         cr.render_y_end() + CLIP_CUT_ROW_PAD};
    for (int i = 0; i < 2; ++i)
      cb.row_y[i] = TrigLUT<W, H>::cos_phi[hs::clamp(rows[i], 0, H_VIRT - 1)];
  }
  return cb;
}

/**
 * @brief Arc fractions where a geodesic edge crosses the clip band's row and
 *        column boundaries.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b); must have an axis.
 * @param cb Clip-band boundaries from make_clip_cut_bounds.
 * @param ts Output, up to GEODESIC_CLIP_MAX_SPLITS fractions in (0, 1),
 *        ascending and separated enough that no piece is degenerate.
 * @return Number of fractions written.
 * @details Cutting the edge here and gating the pieces replaces a uniform chop:
 * every piece then lies wholly inside or wholly outside the band. The
 * boundaries are the RENDER band's, already widened by the clip margin to cover
 * filter reach; the only spacing added on top is the cull's own footprint
 * (CLIP_CUT_COL_PAD, CLIP_CUT_ROW_PAD), without which the outside piece lands
 * inside the span pad and is kept anyway.
 *
 * Both solves run against pos(ang) = a·cos(ang) + cross(axis, a)·sin(ang), the
 * arc the renderer walks, over the same boundary directions the projection
 * reads (TrigLUT). A boundary meridian's half-plane, direction d in the (x, z)
 * plane, is met where the cross with d vanishes and the dot is positive: one
 * atan2, resolved to the sign the half-plane wants. Longitude is globally
 * monotone along the circle and one edge sweeps at most half a revolution
 * (geodesic_col_span_cols), so a meridian is met at most once; the near-
 * horizontal arc pole that same bound refuses is skipped here too. A boundary
 * row is a latitude, and y(ang) folds to R·cos(ang − delta), so its two
 * crossings come from one acos.
 *
 * The solve's fast trig puts a cut within a small fraction of a pixel of the
 * boundary. That is harmless in either direction: each piece is gated by the
 * exact span of the arc between its own endpoints, so a mis-sided cut draws a
 * piece rather than dropping one.
 */
static inline int geodesic_clip_splits(const Vector &a, const Vector &b,
                                       const GeodesicEdgeSpan &es,
                                       const ClipCutBounds &cb, float *ts) {
  const Vector perp = cross(es.axis, a);
  float angs[GEODESIC_CLIP_MAX_SPLITS];
  int found = 0;

  // Roots repeat every turn and the edge sweeps at most half of one, so each
  // candidate has a single representative in [0, 2pi).
  const auto keep = [&](float ang) {
    if (ang < 0.0f)
      ang += 2.0f * PI_F;
    else if (ang >= 2.0f * PI_F)
      ang -= 2.0f * PI_F;
    if (ang > 0.0f && ang < es.total) {
      HS_CHECK(found < GEODESIC_CLIP_MAX_SPLITS);
      angs[found++] = ang;
    }
  };

  if (cb.cols && std::abs(es.axis.y) >= AXIS_Y_EPS) {
    for (int i = 0; i < 2; ++i) {
      const float dx = cb.col_x[i];
      const float dz = cb.col_z[i];
      const float cross_a = a.x * dz - a.z * dx;
      const float cross_b = b.x * dz - b.z * dx;
      // The cross runs sinusoidally in the arc angle, so over the at-most-half
      // turn the edge sweeps it has one root; equal end signs put that root
      // outside the arc, and the two roots an exactly antipodal edge could hold
      // sit on its endpoints, which keep() drops.
      if ((cross_a < 0.0f) == (cross_b < 0.0f))
        continue;
      const float cross_p = perp.x * dz - perp.z * dx;
      const float dot_a = a.x * dx + a.z * dz;
      const float dot_p = perp.x * dx + perp.z * dz;
      // sin, cos at the root are (-cross_a, cross_p) up to a positive scale, so
      // the half-plane's sign test needs no second trig call.
      const float ang = fast_atan2(-cross_a, cross_p);
      keep(dot_a * cross_p - dot_p * cross_a < 0.0f ? ang + PI_F : ang);
    }
  }

  if (cb.rows) {
    const float radius2 = a.y * a.y + perp.y * perp.y;
    if (radius2 > 0.0f) {
      const float radius = sqrtf(radius2);
      const float delta = fast_atan2(perp.y, a.y);
      // y folds to radius*cos(ang - delta), so the endpoints bound the arc
      // except where an extremum angle falls inside it.
      float y_lo = std::min(a.y, b.y);
      float y_hi = std::max(a.y, b.y);
      if (delta > 0.0f && delta < es.total)
        y_hi = radius;
      if (delta + PI_F < es.total)
        y_lo = -radius;
      for (int i = 0; i < 2; ++i) {
        const float y = cb.row_y[i];
        if (y < y_lo || y > y_hi)
          continue;
        const float half = fast_acos(hs::clamp(y / radius, -1.0f, 1.0f));
        keep(delta - half);
        keep(delta + half);
      }
    }
  }

  for (int i = 1; i < found; ++i) {
    const float key = angs[i];
    int j = i - 1;
    for (; j >= 0 && angs[j] > key; --j)
      angs[j + 1] = angs[j];
    angs[j + 1] = key;
  }

  // A piece shorter than the renderer's own collapse threshold would draw as a
  // dot, so fold coincident cuts (and cuts sitting on an endpoint) away.
  int n = 0;
  float last = 0.0f;
  for (int i = 0; i < found; ++i) {
    if (angs[i] - last < EPS_GEODESIC_SEGMENT ||
        es.total - angs[i] < EPS_GEODESIC_SEGMENT)
      continue;
    last = angs[i];
    ts[n++] = angs[i] / es.total;
  }
  return n;
}

/**
 * @brief Exact clip visibility of one geodesic edge.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam ColSpanFn bool(int &col_s, int &col_len), the edge's column arc;
 *         false when no useful bound exists.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param ra Precomputed y_to_screen_row<H>(a.y).
 * @param rb Precomputed y_to_screen_row<H>(b.y).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param col_span Column-arc source, evaluated only once the row span survives
 *        and x clipping is active.
 * @return True if the rendered edge could produce a pixel inside the clip.
 * @details The single definition of the geodesic segment cull: rasterize
 * evaluates it per edge through edge_visible_in_clip, and the hoisted trail
 * gates evaluate it per edge from shared per-point coordinates. Every caller
 * must agree exactly or the paths diverge.
 */
template <int W, int H, typename ColSpanFn>
static __attribute__((always_inline)) inline bool
exact_geodesic_edge_visible(const ClipRegion &cr, const ClipRegion::XClip &xc,
                            float ra, float rb, const Vector &a,
                            const Vector &b, const GeodesicEdgeSpan &es,
                            ColSpanFn &&col_span) {
  float row_lo, row_hi;
  geodesic_row_span_rows<H>(ra, rb, a, b, es, row_lo, row_hi);
  if (!cr.could_intersect_y(row_lo, row_hi))
    return false;
  if (!xc.active)
    return true;
  int col_s, col_len;
  if (!col_span(col_s, col_len))
    return true;
  return ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W);
}

/**
 * @brief Exact clip visibility of one geodesic edge, from a trail's hoisted
 *        per-point screen coordinates.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param rows Per-point screen rows from trail_gate_prologue.
 * @param cols Per-point screen columns from trail_gate_prologue; read only when
 *        x clipping is active, so it may be null otherwise.
 * @param e Edge index; the edge runs from point @p e to point @p e + 1.
 * @param a Edge start, trail point @p e.
 * @param b Edge end, trail point @p e + 1.
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @return True if the rendered edge could produce a pixel inside the clip.
 */
template <int W, int H>
static __attribute__((always_inline)) inline bool
exact_geodesic_edge_visible_hoisted(const ClipRegion &cr,
                                    const ClipRegion::XClip &xc,
                                    const float *rows, const float *cols,
                                    size_t e, const Vector &a, const Vector &b,
                                    const GeodesicEdgeSpan &es) {
  return exact_geodesic_edge_visible<W, H>(
      cr, xc, rows[e], rows[e + 1], a, b, es, [&](int &col_s, int &col_len) {
        return geodesic_col_span_cols<W>(cols[e], cols[e + 1], a, es, col_s,
                                         col_len);
      });
}

enum class RawGeodesicGateResult : uint8_t {
  CULLED,
  VISIBLE,
  EXACT_FALLBACK,
};

/**
 * @brief Gates a regular geodesic edge without angle or cross normalization.
 * @return Visibility, or EXACT_FALLBACK for numerically sensitive geometry.
 */
template <int W, int H>
static inline RawGeodesicGateResult
raw_geodesic_edge_gate(const ClipRegion &cr, const ClipRegion::XClip &xc,
                       float ra, float rb, float ca, float cb, const Vector &a,
                       const Vector &b) {
  constexpr float END_GUARD2 = 4.0e-6f;
  constexpr float AXIS_GUARD2 = 1.0e-4f;
  constexpr float TANGENT_GUARD2 = 1.0e-8f;
  constexpr float ROW_BOUNDARY_GUARD = 0.01f;
  const Vector c = cross(a, b);
  const float L2 = dot(c, c);
  const float d = dot(a, b);
  if (L2 <= END_GUARD2 || std::abs(d) >= 1.0f - END_GUARD2 * 0.5f)
    return RawGeodesicGateResult::EXACT_FALLBACK;

  const float cy2 = c.y * c.y;
  if (cy2 <= AXIS_GUARD2 * L2)
    return RawGeodesicGateResult::EXACT_FALLBACK;

  float row_lo = std::min(ra, rb);
  float row_hi = std::max(ra, rb);
  const float t0 = c.z * a.x - c.x * a.z;
  const float t1 = c.z * b.x - c.x * b.z;
  if (t0 * t0 <= TANGENT_GUARD2 * L2 || t1 * t1 <= TANGENT_GUARD2 * L2)
    return RawGeodesicGateResult::EXACT_FALLBACK;
  if ((t0 > 0.0f) != (t1 > 0.0f)) {
    const float peak = sqrtf(std::max(0.0f, (L2 - cy2) / L2));
    const float rp = y_to_screen_row<H>(t0 > 0.0f ? peak : -peak);
    row_lo = std::min(row_lo, rp);
    row_hi = std::max(row_hi, rp);
  }

  const float y_start = static_cast<float>(cr.render_y_start());
  const float y_end = static_cast<float>(cr.render_y_end());
  if (std::abs(row_hi - y_start) < ROW_BOUNDARY_GUARD ||
      std::abs(row_lo - y_end) < ROW_BOUNDARY_GUARD)
    return RawGeodesicGateResult::EXACT_FALLBACK;
  if (!cr.could_intersect_y(row_lo, row_hi))
    return RawGeodesicGateResult::CULLED;
  if (!xc.active)
    return RawGeodesicGateResult::VISIBLE;

  float col_start;
  float col_length;
  if (c.y < 0.0f) {
    col_start = ca;
    col_length = cb - ca;
  } else {
    col_start = cb;
    col_length = ca - cb;
  }
  if (col_length < 0.0f)
    col_length += W;

  int col_s, col_len;
  finish_col_span_one_period<W>(col_start, col_length, col_s, col_len);
  return ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W)
             ? RawGeodesicGateResult::VISIBLE
             : RawGeodesicGateResult::CULLED;
}

/**
 * @brief Planar screen-column arc from a precomputed edge setup.
 * @tparam W Rasterization width (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param planar_basis Azimuthal-equidistant projection basis.
 * @param es Shared setup from make_planar_edge_span(a, b, planar_basis); the
 *           end point enters through its projected chord.
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 * @param end_sample Optional output for the unprojected edge endpoint.
 * @return False when no useful bound exists — the edge nears a pole or the
 *         Lipschitz margin exceeds the short-way-delta proof — the caller
 *         must skip the horizontal cull.
 * @details Longitude is not monotone along the chart line, so accumulate
 * short-way column deltas over the shared samples and widen by the azimuth
 * Lipschitz bound (W/2π)/sin(φ) per inter-sample gap. sin(φ) over the whole
 * edge is bounded below from the samples minus the 1-Lipschitz gap drift; the
 * short-way delta reading is valid only while the per-gap bound stays under
 * W/4, and both that and the near-pole case fall back to no-cull (return
 * false).
 */
template <int W>
static inline bool planar_col_span(const Vector &a, const Basis &planar_basis,
                                   const PlanarEdgeSpan &es, int &col_s,
                                   int &col_len, Vector *end_sample = nullptr) {
  const float ca = vector_to_theta<W>(a);
  float s_f, len_f;

  {
    float min_sp2 =
        1.0f - a.y * a.y; // squared sin(phi), minimized over samples
    float cum = 0.0f, cum_lo = 0.0f, cum_hi = 0.0f;
    float prev = ca;
    auto step = [&](const Vector &s) {
      float c = vector_to_theta<W>(s);
      float d = c - prev;
      if (d > W * 0.5f)
        d -= W;
      else if (d < -W * 0.5f)
        d += W;
      cum += d;
      cum_lo = std::min(cum_lo, cum);
      cum_hi = std::max(cum_hi, cum);
      prev = c;
      min_sp2 = std::min(min_sp2, 1.0f - s.y * s.y);
    };
    for (const Vector &s : es.interior)
      step(s);
    Vector end = azimuthal_unproject(es.p1.first + es.dX, es.p1.second + es.dY,
                                     planar_basis);
    if (end_sample != nullptr)
      *end_sample = end;
    step(end);

    // Worst-case sin(phi) anywhere on the edge: phi is 1-Lipschitz in arc
    // length and sin is 1-Lipschitz in phi, so between samples it drifts by at
    // most gap_arc.
    const float sin_phi_worst = sqrtf(std::max(0.0f, min_sp2)) - es.gap_arc;
    if (sin_phi_worst < MIN_SIN_PHI)
      return false;
    // Column movement inside one gap; also the proof bound for reading each
    // sample-to-sample delta the short way (must stay well under W/2).
    const float margin =
        es.gap_arc * (static_cast<float>(W) / (2.0f * PI_F)) / sin_phi_worst +
        1.0f;
    if (margin >= W * 0.25f)
      return false;

    s_f = ca + cum_lo - margin;
    len_f = (cum_hi - cum_lo) + 2.0f * margin;
  }

  finish_col_span<W>(s_f, len_f, col_s, col_len);
  return true;
}

/** @brief Tests planar edge visibility from an existing cull sample set. */
template <int W, int H>
static inline bool planar_edge_visible_in_clip(const ClipRegion &cr,
                                               const ClipRegion::XClip &xc,
                                               const Vector &a, const Vector &b,
                                               const Basis &planar_basis,
                                               const PlanarEdgeSpan &span,
                                               Vector *end_sample = nullptr) {
  float row_lo, row_hi;
  int col_s, col_len;
  planar_row_span<H>(a, b, span, row_lo, row_hi);
  if (!cr.could_intersect_y(row_lo, row_hi))
    return false;
  if (!xc.active)
    return true;
  if (!planar_col_span<W>(a, planar_basis, span, col_s, col_len, end_sample))
    return true;
  return ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W);
}

/**
 * @brief One-refinement reciprocal square root for adaptive screen spacing.
 * @details The biased refinement bounds relative error to about 0.089%.
 */
static __attribute__((always_inline)) inline float screen_rsqrt(float x) {
  uint32_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  bits = 0x5f37642fu - (bits >> 1);
  float estimate;
  std::memcpy(&estimate, &bits, sizeof(estimate));
  return estimate * (1.500883f - 0.5f * x * estimate * estimate);
}

/**
 * @brief Adaptive sub-step length (radians of arc) for ~one-pixel screen steps.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param pos Current unit-sphere sample position.
 * @param tan Unit tangent at @p pos (with respect to arc length).
 * @param base_step Equatorial step 2π/W; also the maximum returned step.
 * @return Arc-length step that advances ~SCREEN_STEP_PX pixels on screen.
 * @details Converts the object-space tangent to a screen-space velocity (pixels
 * per radian of arc) under the canvas map x = θ·W/2π, y = φ·(H_VIRT-1)/π, then
 * returns step = SCREEN_STEP_PX / |v_screen|. With φ the colatitude and λ the
 * longitude, dφ/ds = -tan.y/sin(φ) and dλ/ds = (pos.x·tan.z - pos.z·tan.x)/sin²φ.
 * Tracking the full 2-D screen speed (not just longitudinal pole-crowding)
 * deposits ~one sample per pixel everywhere on the curve.
 *
 * Clamped to [base_step·MIN_POLE_SCALE, base_step]: the lower bound caps
 * oversampling at the poles (where dλ/ds diverges → speed → ∞ → step → 0); the
 * upper bound keeps the equator near one sample per column.
 */
template <int W, int H>
static inline float screen_step(const Vector &pos, const Vector &tan,
                                float base_step) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  const float KX = W / (2.0f * PI_F);   // columns per radian of longitude
  const float KY = (H_VIRT - 1) / PI_F; // rows per radian of colatitude
  // sin²φ = 1 - y²; floored so the pole (sin φ → 0) yields a finite, large
  // velocity (hence the min-clamped step) rather than a divide-by-zero.
  const float sin2 = std::max(1e-7f, 1.0f - pos.y * pos.y);
  const float vx_num = KX * (pos.x * tan.z - pos.z * tan.x);
  const float vy_num = KY * tan.y;
  // Factoring the common sin(phi) denominator avoids a separate reciprocal
  // square root while preserving the screen-speed floor.
  const float speed2_num =
      std::max(vx_num * vx_num + vy_num * vy_num * sin2, 1e-12f * sin2 * sin2);
  // Degenerate-speed floor: guards 1/speed when a zero/near-zero tangent stalls
  // the curve, yielding base_step rather than an unbounded step.
  const float step = SCREEN_STEP_PX * sin2 * screen_rsqrt(speed2_num);
  return std::max(base_step * MIN_POLE_SCALE, std::min(step, base_step));
}

#if HS_ENABLE_TEST_ORACLES
inline bool g_reference_screen_step = false;

/** @brief Caps rasterize()'s per-segment sub-step budget; 0 leaves it alone. */
inline size_t g_step_budget_override = 0;

template <int W, int H>
static inline float screen_step_reference(const Vector &pos, const Vector &tan,
                                          float base_step) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  const float KX = W / (2.0f * PI_F);
  const float KY = (H_VIRT - 1) / PI_F;
  const float sin2 = std::max(1e-7f, 1.0f - pos.y * pos.y);
  const float inv_sin = fast_rsqrt(sin2);
  const float dphi_ds = -tan.y * inv_sin;
  const float dlon_ds = (pos.x * tan.z - pos.z * tan.x) * inv_sin * inv_sin;
  const float vx = KX * dlon_ds;
  const float vy = KY * dphi_ds;
  const float speed2 = std::max(vx * vx + vy * vy, 1e-12f);
  const float step = SCREEN_STEP_PX * fast_rsqrt(speed2);
  return std::max(base_step * MIN_POLE_SCALE, std::min(step, base_step));
}
#endif

/**
 * @brief True when @p P statically declares it has no world cull stage, so a
 *        cull predicate may be evaluated against the raw geometry.
 * @tparam P Pipeline type; types without the has_world_cull member (e.g. the
 *           type-erased PipelineRef) are conservatively not hoistable.
 */
template <typename P> static consteval bool pipeline_hoistable_cull() {
  if constexpr (requires { P::has_world_cull; })
    return !P::has_world_cull;
  else
    return false;
}

/**
 * @brief True when @p P statically declares it has no world-space stage, so a
 *        caller may plot a point through precomputed screen coordinates.
 * @tparam P Pipeline type; types without the has_world_stage member (e.g. the
 *           type-erased PipelineRef) are conservatively not hoistable.
 */
template <typename P> static consteval bool pipeline_hoistable_projection() {
  if constexpr (requires { P::has_world_stage; })
    return !P::has_world_stage;
  else
    return false;
}

template <typename P> static consteval bool pipeline_direct_raster_path() {
  if constexpr (requires { P::direct_raster_path; })
    return P::direct_raster_path;
  else
    return false;
}

/**
 * @brief Conservative screen-length test: true only when the geodesic edge
 *        a->b provably spans at most SCREEN_STEP_PX on screen.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @details Tightened (never looser) form of the rasterizer fast-path test
 * `total_dist <= screen_step(sample(0))`, in multiplies only — no trig,
 * divides or square roots. theta and sin(theta) are eliminated via
 * sin(theta)*tangent = b - a*cos(theta) and theta/sin(theta) <= F on
 * theta <= base_step (enforced by the chord cap, which also keeps the edge
 * under screen_step's upper clamp). A false negative falls through to the
 * exact test; true also implies theta >= EPS_GEOMETRIC, so a routed edge can
 * never be one process_segment would have treated as degenerate.
 */
HS_O3_BEGIN
template <int W, int H>
static inline bool edge_fits_one_dot(const Vector &a, const Vector &b) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  constexpr float BASE = (2.0f * PI_F) / W;
  constexpr float B2 = BASE * BASE;
  static_assert(B2 < 1.0f, "chord/angle bounds assume base_step < 1 rad");
  constexpr float KX2 = (W / (2.0f * PI_F)) * (W / (2.0f * PI_F));
  constexpr float KY2 = ((H_VIRT - 1) / PI_F) * ((H_VIRT - 1) / PI_F);
  constexpr float SPX2 = SCREEN_STEP_PX * SCREEN_STEP_PX;
  // Preserve the fast-path implication under screen_rsqrt's <0.1% undershoot.
  constexpr float SCREEN_RSQRT_MIN2 = 0.999f * 0.999f;
  // chord^2 caps: (2 sin(BASE/2))^2 >= B2*(1 - B2/12) bounds theta <= BASE.
  // The lower cap keeps 1 - dot(a, b) orders of magnitude above float ULP:
  // below ~3.5e-4 rad the dot rounds to 1.0f, angle_between collapses to 0,
  // and the exact path treats the edge as degenerate (no interior dot).
  constexpr float CHORD2_MAX = B2 * (1.0f - B2 / 12.0f);
  constexpr float CHORD2_MIN = 4.0e-6f;
  // (theta/sin(theta))^2 <= F2 for theta <= BASE, plus float-rounding slack.
  constexpr float F2 = (1.0001f / ((1.0f - B2 / 6.0f) * (1.0f - B2 / 6.0f)));
  const Vector d = b - a;
  const float chord2 = dot(d, d);
  if (chord2 > CHORD2_MAX || chord2 < CHORD2_MIN)
    return false;
  const float sin2 = 1.0f - a.y * a.y;
  if (sin2 < 1e-7f)
    return false;
  const float c = dot(a, b);
  const float cx = a.x * b.z - a.z * b.x;
  const float ty = b.y - c * a.y;
  return F2 * (KX2 * cx * cx + KY2 * ty * ty * sin2) <=
         SPX2 * SCREEN_RSQRT_MIN2 * sin2 * sin2;
}
HS_O3_END

/**
 * @brief True when AntiAlias would emit any tap of a projected dot in @p cr.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param row Precomputed projected row.
 * @param col Precomputed projected column; unused when x clipping is inactive.
 * @details Mirrors Screen::AntiAlias's edge renormalization and 1e-8 tap
 * cutoff. The gate runs before shading, so it tests tap geometry only.
 */
template <int W, int H>
static inline bool antialiased_dot_visible_in_clip(const ClipRegion &cr,
                                                   const ClipRegion::XClip &xc,
                                                   float row, float col) {
  const float y_floor = floorf(row);
  const int y0 = static_cast<int>(y_floor);
  const int y1 = y0 + 1;
  const bool y0_ok = y0 >= 0 && y0 < H;
  const bool y1_ok = y1 >= 0 && y1 < H;
  if ((!y0_ok || !cr.contains_y(y0)) && (!y1_ok || !cr.contains_y(y1)))
    return false;

  const float x_floor = floorf(col);
  const int x0 = fast_wrap(static_cast<int>(x_floor), W);
  const int x1 = fast_wrap(x0 + 1, W);
  const float xs = quintic_kernel(col - x_floor);
  const float ys = quintic_kernel(row - y_floor);
  float wy0 = 1.0f - ys;
  float wy1 = ys;
  if (y0_ok && !y1_ok) {
    wy0 = 1.0f;
    wy1 = 0.0f;
  } else if (!y0_ok && y1_ok) {
    wy0 = 0.0f;
    wy1 = 1.0f;
  }
  auto visible = [&](int x, int y, float weight) {
    return weight > 1e-8f && cr.contains_y(y) && !xc.clipped(x);
  };
  return (y0_ok &&
          (visible(x0, y0, (1.0f - xs) * wy0) || visible(x1, y0, xs * wy0))) ||
         (y1_ok &&
          (visible(x0, y1, (1.0f - xs) * wy1) || visible(x1, y1, xs * wy1)));
}

/**
 * @brief Tier-3 clip visibility of one polyline edge, routed through the
 *        pipeline's world stages.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam PipelineT Pipeline type.
 * @param pipeline Render pipeline; world stages re-emit the edge under their
 *        plot-time rotations before the span test.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param pb Planar projection basis for the edge, or null for geodesic.
 * @return True if the rendered edge could produce a pixel inside the clip.
 * @details Geodesic edges route to exact_geodesic_edge_visible; the planar
 * branch is the single definition of the planar segment cull.
 */
template <int W, int H, typename PipelineT>
static inline bool
edge_visible_in_clip(PipelineT &pipeline, const ClipRegion &cr,
                     const ClipRegion::XClip &xc, const Vector &a,
                     const Vector &b, const Basis *pb) {
  auto pred = [&](const Vector &ea, const Vector &eb, const Basis *bp) {
    if (bp == nullptr) {
      // Geodesic: both span bounds share one edge setup (angle, arc pole).
      const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
      return exact_geodesic_edge_visible<W, H>(
          cr, xc, y_to_screen_row<H>(ea.y), y_to_screen_row<H>(eb.y), ea, eb,
          es, [&](int &col_s, int &col_len) {
            return geodesic_col_span<W>(ea, eb, es, col_s, col_len);
          });
    }
    // Planar: both span bounds share one projection + chart-line sample set.
    const PlanarEdgeSpan ps = make_planar_edge_span(ea, eb, *bp);
    return planar_edge_visible_in_clip<W, H>(cr, xc, ea, eb, *bp, ps);
  };
  if constexpr (requires { pipeline.could_intersect_clip(a, b, pb, pred); }) {
    return pipeline.could_intersect_clip(a, b, pb, pred);
  } else {
    // A filter pipeline (has any_crosses_segments) must answer the clip
    // query; a signature drift would else silently fall back to raw culling.
    static_assert(
        !requires { PipelineT::any_crosses_segments; },
        "pipeline exposes any_crosses_segments but not "
        "could_intersect_clip (signature drift)");
    return pred(a, b, pb);
  }
}

/**
 * @brief Conservative test: can a spherical cap reach a clip's render region?
 * @tparam H Canvas height in rows.
 * @param cr Clip region to test against.
 * @param dir Cap center direction (unit vector). A ring passes its axis with
 * half_angle = colatitude + displacement bound (the cap of that radius
 * contains the ring's band); a ring chunk passes its midpoint with
 * half_angle = chunk half-arc + displacement bound.
 * @param half_angle Cap angular radius including stroke/AA pad (radians).
 * @return False only when no fragment inside the cap can land in the clip's
 * render region; true is always safe.
 * @details Rows: the cap's polar range about the display's Y axis is
 * [beta - t2, beta + t2] (beta = center colatitude). Columns: a cap that
 * reaches either display pole spans all longitudes; otherwise its longitude
 * half-width about the center longitude is asin(sin t2 / sin beta), compared
 * against the clip's column wedge with the clip margin plus one pixel of
 * slack.
 */
template <int H>
inline bool cap_may_touch_clip(const ClipRegion &cr, const Vector &dir,
                               float half_angle) {
  float t2 = std::min(half_angle, PI_F);
  float beta = acosf(hs::clamp(dir.y, -1.0f, 1.0f));

  float phi_lo = std::max(beta - t2, 0.0f);
  float phi_hi = std::min(beta + t2, PI_F);
  if (phi_to_y<H>(phi_hi) < cr.render_y_start() ||
      phi_to_y<H>(phi_lo) >= cr.render_y_end())
    return false;

  if (cr.x_start == 0 && cr.x_end == cr.w)
    return true;
  if (beta <= t2 || PI_F - beta <= t2)
    return true;
  float dlam = asinf(hs::clamp(sinf(t2) / sinf(beta), 0.0f, 1.0f));
  float lam_v = atan2f(dir.z, dir.x);
  float width_px = static_cast<float>(cr.x_end - cr.x_start);
  float half_w = (width_px * 0.5f + cr.margin + 1.0f) * (2.0f * PI_F) / cr.w;
  float lam_c = (cr.x_start + width_px * 0.5f) * (2.0f * PI_F) / cr.w;
  float d = std::fabs(wrap_t((lam_v - lam_c) / (2.0f * PI_F) + 0.5f) - 0.5f) *
            (2.0f * PI_F);
  return d <= dlam + half_w;
}

enum class CartesianTrailGateResult : uint8_t {
  EXACT_FALLBACK,
  LATITUDE_REJECT,
  MERIDIAN_REJECT,
};

struct CartesianQuadrantClip {
  float latitude_sign = 0.0f;
  float latitude_threshold = 0.0f;
  float meridian_sign = 0.0f;
  float meridian_threshold = 0.0f;
  bool active = false;
};

/**
 * @brief Builds the Cartesian halfspace superset for a segmented quadrant.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @return Enabled thresholds only for the four hardware quadrant shapes.
 */
template <int W, int H>
static CartesianQuadrantClip
make_cartesian_quadrant_clip(const ClipRegion &cr) {
  CartesianQuadrantClip q;
  if (cr.w != W || cr.h != H || cr.margin < 0 || W % 2 != 0 || H % 2 != 0 ||
      cr.x_end - cr.x_start != W / 2 ||
      (cr.x_start != 0 && cr.x_start != W / 2) ||
      cr.y_end - cr.y_start != H / 2 ||
      (cr.y_start != 0 && cr.y_start != H / 2))
    return q;

  constexpr int H_VIRT = H + hs::H_OFFSET;
  if (cr.y_start == 0) {
    const float boundary = static_cast<float>(cr.render_y_end()) * PI_F /
                           static_cast<float>(H_VIRT - 1);
    q.latitude_sign = 1.0f;
    q.latitude_threshold = cosf(boundary);
  } else {
    const float boundary = static_cast<float>(cr.render_y_start()) * PI_F /
                           static_cast<float>(H_VIRT - 1);
    q.latitude_sign = -1.0f;
    q.latitude_threshold = -cosf(boundary);
  }

  const float half_width =
      PI_F * 0.5f + static_cast<float>(cr.margin + COL_FOOTPRINT) *
                        (2.0f * PI_F / static_cast<float>(W));
  if (half_width >= PI_F)
    return CartesianQuadrantClip{};
  q.meridian_sign = cr.x_start == 0 ? 1.0f : -1.0f;
  q.meridian_threshold = cosf(half_width);
  q.active = true;
  return q;
}

/**
 * @brief Conservatively rejects a geodesic trail outside a Cartesian quadrant.
 * @param clip Precomputed quadrant halfspace thresholds.
 * @param trail Unit-sphere geodesic polyline.
 * @return Rejecting halfspace, or exact fallback for every uncertain case.
 */
HS_O3_BEGIN
static inline CartesianTrailGateResult
cartesian_quadrant_trail_gate(const CartesianQuadrantClip &clip,
                              const Fragments &trail) {
  if (!clip.active || trail.size() < 2)
    return CartesianTrailGateResult::EXACT_FALLBACK;

#ifdef HS_PROFILE_MINDSPLATTER_STALLS
  hs::DwtStallBatch gate_batch(hs::g_mindsplatter_stalls.trail_gate);
#endif
  float latitude_max = -1.0f;
  float max_chord2 = 0.0f;
  for (size_t k = 0; k < trail.size(); ++k) {
    const Vector &p = trail[k].pos;
    latitude_max = std::max(latitude_max, clip.latitude_sign * p.y);
    if (k > 0) {
      const Vector d = p - trail[k - 1].pos;
      max_chord2 = std::max(max_chord2, dot(d, d));
    }
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
    gate_batch.step();
#endif
  }

  // Every point on a minor arc lies within half its arc of one endpoint, and
  // arc <= (pi/2)*chord. A unit-normal dot changes by at most angular distance.
  const float slack = (PI_F * 0.25f) * sqrtf(max_chord2);
  if (latitude_max + slack < clip.latitude_threshold - math::EPS_GEOMETRIC) {
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
    gate_batch.step();
    gate_batch.finish();
#endif
    return CartesianTrailGateResult::LATITUDE_REJECT;
  }

  float meridian_max = -1.0f;
  for (const Fragment &f : trail) {
    meridian_max = std::max(meridian_max, clip.meridian_sign * f.pos.z);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
    gate_batch.step();
#endif
  }
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
  gate_batch.step();
  gate_batch.finish();
#endif
  if (meridian_max + slack < clip.meridian_threshold - math::EPS_GEOMETRIC)
    return CartesianTrailGateResult::MERIDIAN_REJECT;
  return CartesianTrailGateResult::EXACT_FALLBACK;
}
HS_O3_END

static inline void
count_cartesian_trail_gate_result(CartesianTrailGateResult result) {
  if (result == CartesianTrailGateResult::LATITUDE_REJECT)
    HS_MSP_COUNT(cartesian_latitude_rejects);
  else if (result == CartesianTrailGateResult::MERIDIAN_REJECT)
    HS_MSP_COUNT(cartesian_meridian_rejects);
  else
    HS_MSP_COUNT(cartesian_fallbacks);
#if defined(HS_PROFILE_ENABLE) && defined(HS_PROFILE_CARTESIAN_COUNTS)
  static hs::CycleCounter latitude("plot_ps_cartesian_latitude_reject");
  static hs::CycleCounter meridian("plot_ps_cartesian_meridian_reject");
  static hs::CycleCounter fallback("plot_ps_cartesian_fallback");
  hs::CycleCounter *counter = &fallback;
  if (result == CartesianTrailGateResult::LATITUDE_REJECT)
    counter = &latitude;
  else if (result == CartesianTrailGateResult::MERIDIAN_REJECT)
    counter = &meridian;
  ++counter->count;
#else
  (void)result;
#endif
}

static inline void count_particle_edge_class(bool one_dot) {
  if (one_dot)
    HS_MSP_COUNT(one_dot_edges);
  else
    HS_MSP_COUNT(long_edges);
#if defined(HS_PROFILE_ENABLE) && defined(HS_PROFILE_EDGE_CLASS_COUNTS)
  static hs::CycleCounter one_dot_count("plot_ps_edge_one_dot");
  static hs::CycleCounter long_count("plot_ps_edge_long");
  ++(one_dot ? one_dot_count : long_count).count;
#else
  (void)one_dot;
#endif
}

static inline void count_particle_exact_gate_fallback() {
  HS_MSP_COUNT(exact_gate_fallbacks);
#if defined(HS_PROFILE_ENABLE) && defined(HS_PROFILE_EDGE_CLASS_COUNTS)
  static hs::CycleCounter exact_count("plot_ps_edge_exact_fallback");
  ++exact_count.count;
#endif
}

/**
 * @brief Hoisted per-point screen coordinates and whole-trail cull verdict.
 */
struct TrailGatePrologue {
  const float *rows; /**< Per-point screen rows, one per trail point. */
  const float
      *cols;     /**< Per-point screen columns, null when x-clip inactive. */
  bool rejected; /**< The whole trail is provably outside the clip. */
};

/**
 * @brief Computes one geodesic trail's per-point rows/columns and applies the
 *        conservative whole-trail row and column culls.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param trail Geodesic fragment polyline (>= 2 unit-position points).
 * @return The hoisted arrays plus whether the trail is culled whole.
 * @details rows and cols are bump-allocated from scratch_arena_a and the helper
 * opens no ScratchScope of its own: they stay valid until the CALLER's scope
 * unwinds, so a caller may keep them alive past the gate.
 */
template <int W, int H>
static __attribute__((always_inline)) inline TrailGatePrologue
trail_gate_prologue(const ClipRegion &cr, const ClipRegion::XClip &xc,
                    const Fragments &trail) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  const size_t n = trail.size();
  auto *rows = static_cast<float *>(
      scratch_arena_a.allocate(n * sizeof(float), alignof(float)));
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
  hs::DwtStallBatch gate_batch(hs::g_mindsplatter_stalls.trail_gate);
#endif
  float row_lo_t = 1e9f, row_hi_t = -1e9f;
  float min_sp2 = 1.0f;
  float max_chord2 = 0.0f;
  for (size_t k = 0; k < n; ++k) {
    const Vector &pt = trail[k].pos;
    rows[k] = y_to_screen_row<H>(pt.y);
    row_lo_t = std::min(row_lo_t, rows[k]);
    row_hi_t = std::max(row_hi_t, rows[k]);
    min_sp2 = std::min(min_sp2, 1.0f - pt.y * pt.y);
    if (k > 0) {
      const Vector d = pt - trail[k - 1].pos;
      max_chord2 = std::max(max_chord2, dot(d, d));
    }
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
    gate_batch.step();
#endif
  }
  // arc <= (pi/2)*chord on [0, pi]; an edge's interior latitude extremum lies
  // within arc/2 of an endpoint and phi is 1-Lipschitz in arc length, so this
  // margin covers every per-edge bulge peak.
  const float max_arc = (PI_F * 0.5f) * sqrtf(max_chord2);
  const float row_margin =
      (max_arc * 0.5f) * (static_cast<float>(H_VIRT - 1) / PI_F);
  if (!cr.could_intersect_y(row_lo_t - row_margin, row_hi_t + row_margin)) {
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
    gate_batch.step();
    gate_batch.finish();
#endif
    HS_MSP_COUNT(prologue_row_rejects);
    return {rows, nullptr, true};
  }

  float *cols = nullptr;
  if (xc.active) {
    cols = static_cast<float *>(
        scratch_arena_a.allocate(n * sizeof(float), alignof(float)));
    float cum = 0.0f, cum_lo = 0.0f, cum_hi = 0.0f;
    bool walk_safe = true;
    cols[0] = vector_to_theta<W>(trail[0].pos);
    for (size_t k = 1; k < n; ++k) {
      cols[k] = vector_to_theta<W>(trail[k].pos);
      // A geodesic edge's column sweep never exceeds W/2 (antipodal symmetry,
      // see geodesic_col_span_cols), so the short-way delta covers it
      // regardless of direction — except at ~exactly W/2, where the delta's
      // sign (which semicircle) is float noise.
      float d = cols[k] - cols[k - 1];
      if (d > W * 0.5f)
        d -= W;
      else if (d < -W * 0.5f)
        d += W;
      if (std::abs(d) >= W * 0.5f - 3.0f)
        walk_safe = false;
      // geodesic_col_span_cols refuses to bound an edge whose great-circle
      // axis is near-horizontal, and the per-edge tier then treats it as
      // visible; the endpoint columns walked here do not bound such an edge
      // either. |axis.y| = |cy| / |cross| and |cross| <= 1, so testing the
      // unnormalized cy covers every case it rejects.
      const Vector &ca_pos = trail[k - 1].pos;
      const Vector &cb_pos = trail[k].pos;
      if (std::abs(ca_pos.z * cb_pos.x - ca_pos.x * cb_pos.z) < AXIS_Y_EPS)
        walk_safe = false;
      cum += d;
      cum_lo = std::min(cum_lo, cum);
      cum_hi = std::max(cum_hi, cum);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
      gate_batch.step();
#endif
    }
    // Near a pole the plotted column is float noise (same caution as the
    // per-edge spans), so only cull by the column arc when the whole trail
    // provably stays clear.
    if (walk_safe && sqrtf(std::max(0.0f, min_sp2)) - max_arc >= MIN_SIN_PHI) {
      int col_s, col_len;
      finish_col_span<W>(cols[0] + cum_lo, cum_hi - cum_lo, col_s, col_len);
      if (!ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W)) {
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
        gate_batch.step();
        gate_batch.finish();
#endif
        HS_MSP_COUNT(prologue_column_rejects);
        return {rows, cols, true};
      }
    }
  }
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
  gate_batch.step();
  gate_batch.finish();
#endif
  return {rows, cols, false};
}

} // namespace Plot
