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
#include "engine/constants.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"

/**
 * @file plot.h
 * @brief The curve rasterizer: edge samplers and the Plot primitives that
 * stroke lines, rings, polygons, meshes and particle systems.
 */

namespace Plot {

/**
 * @brief Inner/outer radius ratio for star shapes.
 */
static constexpr float STAR_INNER_RATIO = ::STAR_INNER_RATIO;

/**
 * @brief Geodesic segment shorter than this (radians) collapses to a point.
 * @details 100× math::EPS_GEOMETRIC (1e-3 vs 1e-5): a slerp-axis stability
 * bound that picks the interpolation strategy, not a positional near-equality
 * test, so it does not track math::EPS_GEOMETRIC.
 */
static constexpr float EPS_GEODESIC_SEGMENT = 0.001f;

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
static constexpr float EPS_ARC_POLE_SQ = 1e-8f;

/**
 * @brief Minimum |axis.y| for which a geodesic edge's endpoint columns bound
 *        its azimuth span.
 * @details Below this the great circle runs near the poles, where longitude is
 * ill-conditioned and the interior leaves the endpoint columns.
 */
static constexpr float AXIS_Y_EPS = 1e-4f;

/**
 * @brief Minimum worst-case sin(φ) over a curve for which its plotted columns
 *        are trusted enough to cull by.
 * @details The azimuth Lipschitz bound scales as 1/sin(φ), and nearer the poles
 * the plotted column is float noise; below this the column cull bails.
 */
static constexpr float MIN_SIN_PHI = 0.05f;

/**
 * @brief Floor on the adaptive sub-step length, as a fraction of base_step.
 * @details Caps sub-steps per segment so polar curves don't oversample: the
 * screen-velocity step sampler (screen_step) drives the step toward zero where
 * the azimuthal velocity diverges at the poles, and this is the lower clamp that
 * bounds it. A clamp, not a tolerance.
 */
static constexpr float MIN_POLE_SCALE = 0.05f;

/**
 * @brief Target screen-space spacing (pixels) between adaptive sub-samples.
 * @details The rasterizer sizes each sub-step so consecutive samples land about
 * this far apart in SCREEN space. Slightly sub-pixel so the bilinear AntiAlias
 * splat of neighbouring samples overlaps and the rendered curve has no holes;
 * smaller = denser = smoother but costlier.
 */
static constexpr float SCREEN_STEP_PX = 0.9f;

/** @brief Adaptive raster sampling density. */
enum class RasterSamplingPolicy { DEFAULT, BALANCED, SELECTABLE };

/** @brief Balanced-policy target spacing in screen pixels. */
static constexpr float BALANCED_SCREEN_STEP_PX = 1.125f;

/** @brief Pole-floor multiple below which balanced sampling keeps exact cadence. */
static constexpr float BALANCED_POLE_GUARD_SCALE = 2.0f;

/** @brief Source-over-aware alpha gain for balanced sample spacing. */
static inline float balanced_sample_alpha(float alpha, float step_ratio) {
  const float gain = 1.0f + (step_ratio - 1.0f) * (0.88f - 0.20f * alpha);
  return std::min(1.0f, alpha * gain);
}

#ifdef HS_TEST_BUILD
inline uint32_t g_planar_full_samples = 0;
inline uint32_t g_planar_position_samples = 0;
#endif

/**
 * @brief Parameter delta for the planar strategy's finite-difference tangent.
 * @details The planar (azimuthal) map has no closed-form tangent, so its screen
 * velocity is taken from a short forward difference; small enough to track the
 * arc, large enough to stay clear of float cancellation.
 */
static constexpr float PLANAR_TAN_DT = 1.0f / 256.0f;

/**
 * @brief Antipode cutoff for the planar projection's stable-azimuth region.
 * @details The planar (azimuthal-equidistant) projection is singular at the
 * basis antipode (R→π: azimuth undefined). A control point whose dot with the
 * basis center is below this (≈ within 2.6° of the antipode) projects to an
 * unstable azimuth, so its segment falls back to a geodesic edge. cos(π − 0.045).
 */
static constexpr float COS_PLANAR_ANTIPODE = 0.999f;

/**
 * @brief Columns of slack added on each side of a culled column span.
 * @details Absorbs plot rounding and the AntiAlias tap spread.
 */
static constexpr int COL_PAD = 2;

/**
 * @brief Columns a padded span reaches past its fractional end.
 * @details The pad plus the boundary column ceil() adds.
 */
static constexpr int COL_FOOTPRINT = COL_PAD + 1;

/**
 * @brief Columns outside the render band a clip cut is placed at.
 * @details A piece ending exactly on the band edge still overlaps it once
 * finish_col_span widens the span by COL_FOOTPRINT, so a cut there would leave
 * the outside piece visible and buy nothing. One column past that footprint
 * also absorbs the fast-trig error in the cut.
 */
static constexpr int CLIP_CUT_COL_PAD = COL_FOOTPRINT + 1;

/**
 * @brief Rows outside the render band a clip cut is placed at.
 * @details could_intersect_y takes the band edge itself as an intersection, and
 * the row map's fast-trig round trip moves a cut by a fraction of a row.
 */
static constexpr int CLIP_CUT_ROW_PAD = 1;

/**
 * @brief Apply an optional per-control-point vertex shader to every fragment.
 * @tparam FragmentsT Fragment container type.
 * @param vertex_shader Vertex shader to run on each fragment; no-op if null.
 * @param pts Fragment container mutated in place.
 * @details Shared inline replacement for the per-primitive
 * `if (vertex_shader) for (auto &p : pts) vertex_shader(p);` block.
 */
template <typename FragmentsT>
inline void apply_vertex_shader(VertexShaderRef vertex_shader,
                                FragmentsT &pts) {
  if (vertex_shader) {
    for (auto &p : pts) {
      vertex_shader(p);
    }
  }
}

// --- Strategy Helpers ---
// Core rasterization logic for 3D lines and curves adapts step size based on
// screen-space density to avoid aliasing.

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
  SamplePT one_pass(float s) const {
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
   * @details The azimuthal map has no closed-form tangent, so take it from a
   * short forward difference (backward at the s=1 end); the difference
   * direction is the unit tangent regardless of magnitude. Feeds the same
   * screen-velocity sub-step sampler as the geodesic path.
   */
  SamplePT operator()(float s) const {
    Vector p = pos(s);
    bool fwd = (s + PLANAR_TAN_DT <= 1.0f);
    Vector q = pos(fwd ? s + PLANAR_TAN_DT : s - PLANAR_TAN_DT);
    Vector d = fwd ? (q - p) : (p - q);
    HS_PLOT_COUNT(normalizations);
    return {p, normalized_or(d, Vector())};
  }
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
    const float d = wrap(cb - ca, static_cast<float>(W));
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
      len_f = wrap(ce - ca, static_cast<float>(W));
    } else {
      s_f = ce;
      len_f = wrap(ca - ce, static_cast<float>(W));
    }
  }

  finish_col_span<W>(s_f, len_f, col_s, col_len);
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
static constexpr int GEODESIC_CLIP_MAX_SPLITS = 6;

/**
 * @brief Arc fractions where a geodesic edge crosses the clip band's row and
 *        column boundaries.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b); must have an axis.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
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
template <int W, int H>
static inline int geodesic_clip_splits(const Vector &a,
                                       const GeodesicEdgeSpan &es,
                                       const ClipRegion &cr,
                                       const ClipRegion::XClip &xc, float *ts) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();
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
    if (ang > 0.0f && ang < es.total)
      angs[found++] = ang;
  };

  if (xc.active && std::abs(es.axis.y) >= AXIS_Y_EPS) {
    for (const int col : {xc.rs - CLIP_CUT_COL_PAD, xc.re + CLIP_CUT_COL_PAD}) {
      const int c = ((col % W) + W) % W;
      const float dx = TrigLUT<W, H>::cos_theta(c);
      const float dz = TrigLUT<W, H>::sin_theta[c];
      const float cross_a = a.x * dz - a.z * dx;
      const float cross_p = perp.x * dz - perp.z * dx;
      const float dot_a = a.x * dx + a.z * dz;
      const float dot_p = perp.x * dx + perp.z * dz;
      // sin, cos at the root are (-cross_a, cross_p) up to a positive scale, so
      // the half-plane's sign test needs no second trig call.
      const float ang = fast_atan2(-cross_a, cross_p);
      keep(dot_a * cross_p - dot_p * cross_a < 0.0f ? ang + PI_F : ang);
    }
  }

  if (cr.render_y_start() > 0 || cr.render_y_end() < cr.h) {
    const float radius2 = a.y * a.y + perp.y * perp.y;
    if (radius2 > 0.0f) {
      const float radius = sqrtf(radius2);
      const float delta = fast_atan2(perp.y, a.y);
      for (const int row : {cr.render_y_start() - CLIP_CUT_ROW_PAD,
                            cr.render_y_end() + CLIP_CUT_ROW_PAD}) {
        const float target =
            TrigLUT<W, H>::cos_phi[hs::clamp(row, 0, H_VIRT - 1)] / radius;
        if (target < -1.0f || target > 1.0f)
          continue;
        const float half = fast_acos(target);
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

#ifdef HS_TEST_BUILD
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

/** Precomputed edge byte flags consumed by rasterize(). */
static constexpr uint8_t EDGE_VISIBLE = 1u << 0;
static constexpr uint8_t EDGE_ONE_DOT = 1u << 1;
static constexpr uint8_t EDGE_CLASSIFIED = 1u << 2;

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
 * @brief scratch_arena_a bytes rasterize binds for its own adaptive sub-step
 * cache, on top of the caller's Fragments buffer, which stays live across the
 * call.
 * @tparam W Rasterization width.
 * @return The cache size in bytes.
 * @details Effects sizing a custom scratch-A split must budget this alongside
 * their own buffers. A planar-basis draw binds count * (sizeof(float) + 1) more
 * for its per-segment arc/seam caches, where count is the segment count.
 */
template <int W> inline constexpr size_t rasterize_scratch_a_bytes() {
  constexpr size_t STEPS =
      2 * static_cast<size_t>(W) > 64 ? 2 * static_cast<size_t>(W) : 64;
  return STEPS * sizeof(float);
}

/**
 * @brief Optional rasterize() behaviors beyond the plain open geodesic
 * polyline; every field defaults to that common case.
 */
struct RasterOptions {
  /** Also draw the last→first edge. */
  bool close_loop = false;
  /**
   * Non-null selects azimuthal-equidistant interpolation (straight in the
   * projection); null uses geodesic edges.
   */
  const Basis *planar_basis = nullptr;
  /**
   * Open lines only: skip the final endpoint plot (each vertex is otherwise
   * plotted once by its outgoing segment), so abutting arcs tile a longer
   * curve without double-plotting the shared vertex.
   */
  bool omit_end = false;
  /**
   * Optional precomputed Tier-3 edge flags, one byte per segment-loop edge.
   * Bit 0 is visibility. Bits 1 and 2 optionally retain one-dot and
   * classification-known; legacy 0/1 visibility arrays remain valid.
   * Geodesic polylines only: a planar polyline's per-edge basis depends on
   * rasterize()'s seam pre-pass.
   */
  const uint8_t *edge_visible = nullptr;
  /**
   * Optional per-point screen rows, y_to_screen_row of each points[k].pos.
   * With point_cols, lets the single-dot shortcut skip the projection. Only
   * consumed when the pipeline has no world-space stage; both arrays or
   * neither.
   */
  const float *point_rows = nullptr;
  /**
   * Optional per-point screen columns, vector_to_theta of each
   * points[k].pos.
   */
  const float *point_cols = nullptr;
  /**
   * Optional last-to-first target fragment carrying seam registers for a
   * closed loop without an overlapping point.
   */
  const Fragment *loop_seam = nullptr;
  /** Enables balanced sampling for a SELECTABLE rasterizer. */
  bool balanced_sampling = false;
#ifdef HS_TEST_BUILD
  /** Rebuild a planar sampler after culling instead of reusing cull samples. */
  bool rebuild_planar_sampler = false;
#endif
};

/**
 * @brief Adaptively rasterize a fragment polyline onto the sphere.
 *
 * Walks consecutive fragment pairs, picks a geodesic or planar interpolation
 * strategy per segment, sub-steps each segment at ≈one-pixel SCREEN-space
 * density (screen_step, clamped near the poles), and plots through the pipeline.
 * Segments whose full screen-row span lies outside the active clip band are
 * culled.
 *
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam SinglePass Emit adaptive samples immediately instead of replaying a
 *         normalized step cache.
 * @tparam OpenGeodesic Compile out planar, closed-loop, seam and omit-end
 *         support for an open geodesic polyline.
 * @tparam DerivePlanarArcRegisters Recompute v0/v1 from the rendered planar
 *         perimeter.
 * @tparam InterpolateRegisters Interpolate source fragment registers at each
 *         adaptive sample.
 * @tparam SamplingPolicy Adaptive screen-space sample density.
 * @tparam PipelineT Pipeline type.
 * @tparam FragmentShaderT Fragment shader type for direct raster pipelines.
 * @param source_pipeline Render pipeline that plots fragments.
 * @param canvas Target canvas (supplies the active clip band).
 * @param points Fragment polyline to rasterize.
 * @param fragment_shader Per-fragment shader applied before plotting; must be
 *                        non-null (the per-pixel call sites below do not guard
 *                        it, and operator()'s null assert is stripped under
 *                        NDEBUG on-device).
 * @param opts Optional loop/projection/culling behaviors (see RasterOptions).
 */
HS_O3_BEGIN
template <int W, int H, bool SinglePass = false, bool OpenGeodesic = false,
          bool DerivePlanarArcRegisters = true,
          bool InterpolateRegisters = true,
          RasterSamplingPolicy SamplingPolicy = RasterSamplingPolicy::DEFAULT,
          typename PipelineT = PipelineRef,
          typename FragmentShaderT = FragmentShaderFn>
static void rasterize(PipelineT &source_pipeline, Canvas &canvas,
                      const Fragments &points, FragmentShaderT fragment_shader,
                      const RasterOptions &opts = {}) {
  if constexpr (OpenGeodesic)
    assert(!opts.close_loop && opts.planar_basis == nullptr && !opts.omit_end &&
           opts.loop_seam == nullptr);
  // A direct-raster sink writes through a cached framebuffer base; the canvas
  // double-buffers, so a stale base is the buffer the display is scanning out.
  if constexpr (requires { source_pipeline.prepared_for(canvas); })
    HS_CHECK(source_pipeline.prepared_for(canvas),
             "direct raster pipeline not prepared for this canvas");
  // Erasure collapses pipeline and shader; the erased call matches neither
  // clause so recursion ends.
  if constexpr (!pipeline_direct_raster_path<PipelineT>() &&
                (!std::same_as<std::decay_t<PipelineT>, PipelineRef> ||
                 !std::same_as<std::decay_t<FragmentShaderT>,
                               FragmentShaderFn>)) {
    PipelineRef erased(source_pipeline);
    FragmentShaderFn erased_shader(fragment_shader);
    rasterize<W, H, SinglePass, OpenGeodesic, DerivePlanarArcRegisters,
              InterpolateRegisters, SamplingPolicy>(erased, canvas, points,
                                                    erased_shader, opts);
    return;
  }
  HS_PLOT_COUNT(rings);
  const bool close_loop = OpenGeodesic ? false : opts.close_loop;
  const Basis *planar_basis = OpenGeodesic ? nullptr : opts.planar_basis;
  const bool omit_end = OpenGeodesic ? false : opts.omit_end;
  const uint8_t *edge_visible = opts.edge_visible;
  const float *point_rows = opts.point_rows;
  const float *point_cols = opts.point_cols;
  const Fragment *loop_seam = OpenGeodesic ? nullptr : opts.loop_seam;
  const bool balanced_sampling =
      SamplingPolicy == RasterSamplingPolicy::BALANCED ||
      (SamplingPolicy == RasterSamplingPolicy::SELECTABLE &&
       opts.balanced_sampling);
  auto &pipeline = source_pipeline;
  size_t len = points.size();
  // A degenerate path is not drawn — callers wanting a dot duplicate the vertex,
  // as Line::sample does.
  if (len < 2)
    return;
  // Trap a null shader once per polyline so the per-pixel fragment_shader()
  // calls below can't invoke a null thunk.
  if constexpr (std::same_as<std::decay_t<FragmentShaderT>, FragmentShaderFn>)
    HS_CHECK(fragment_shader, "rasterize requires a non-null fragment_shader");
  HS_CHECK(edge_visible == nullptr || planar_basis == nullptr,
           "precomputed edge visibility is geodesic-only");
  HS_CHECK(loop_seam == nullptr || close_loop,
           "a raster seam fragment requires a closed loop");
  HS_CHECK((point_rows == nullptr) == (point_cols == nullptr),
           "hoisted point projections take both rows and columns");

  size_t count = close_loop ? len : len - 1;
  HS_PLOT_ADD(edges, count);
  // SCRATCH ARENA CONTRACT (load-bearing): scratch_arena_a is a LIFO bump
  // allocator shared with Pixel::Feedback::flush; do not let a raw pointer into
  // it outlive the scope that produced it.
  ScratchScope sc_guard(scratch_arena_a);
  ArenaVector<float> steps_cache;
  // The cache holds ONE segment's adaptive sub-steps (cleared per segment).
  // Away from a pole, ≈SCREEN_STEP_PX spacing gives the usual screen-sweep
  // count. Near a pole, MIN_POLE_SCALE can instead create up to
  // 1/MIN_POLE_SCALE samples per base_step interval; those samples are not
  // covered by the planar W/H sweep derivation. At W=288 the measured
  // pole-crossing geodesic worst case is 546 of the 2·W=576 slots. A planar
  // chart line can bow farther, so the simulation loop retains its capacity
  // backstop. SinglePass emits as it goes and takes max_cache only as that
  // backstop, so it never binds the storage.
  size_t max_cache = rasterize_scratch_a_bytes<W>() / sizeof(float);
#ifdef HS_TEST_BUILD
  if (g_step_budget_override != 0 && g_step_budget_override < max_cache)
    max_cache = g_step_budget_override;
#endif
  if constexpr (!SinglePass)
    steps_cache.bind(scratch_arena_a, max_cache);

  // PLANAR ARC REGISTERS (v0/v1): under a planar basis the rendered edge bows
  // longer than the geodesic chord, so re-derive v0/v1 from the true rendered
  // arc (`cumul`/`seg_base` track it, `total_arc` normalizes v0). Skipped for
  // geodesic polylines or when DerivePlanarArcRegisters is false.
  const bool has_planar_basis = (planar_basis != nullptr);
  const bool override_uv = DerivePlanarArcRegisters && has_planar_basis;
  constexpr bool REUSE_PLANAR_CULL_SAMPLES =
      SinglePass && !DerivePlanarArcRegisters && !InterpolateRegisters &&
      pipeline_hoistable_cull<PipelineT>();
  auto segment_next = [&](size_t i) -> const Fragment & {
    if (loop_seam != nullptr && i + 1 == len)
      return *loop_seam;
    return points[(i + 1) % len];
  };
  float total_arc = 0.0f;
  // Per-segment rendered arc length and antipode-seam flag, reused by the draw
  // loop below so the seam decision is taken in exactly one place.
  ArenaVector<float> seg_arc_cache;
  ArenaVector<uint8_t> seg_seam_cache;
  if (override_uv) {
    seg_arc_cache.bind(scratch_arena_a, count);
    seg_seam_cache.bind(scratch_arena_a, count);
    const Vector &pcenter = planar_basis->v;
    for (size_t i = 0; i < count; i++) {
      const Vector &a = points[i].pos;
      const Vector &b = segment_next(i).pos;
      const bool seam = dot(a, pcenter) < -COS_PLANAR_ANTIPODE ||
                        dot(b, pcenter) < -COS_PLANAR_ANTIPODE;
      seg_seam_cache.push_back(seam ? 1 : 0);
      float seg =
          seam ? angle_between(a, b) : planar_arc_length(a, b, *planar_basis);
      seg_arc_cache.push_back(seg);
      total_arc += seg;
    }
  }
  float cumul = 0.0f;    // rendered arc reached so far (planar polylines only)
  float seg_base = 0.0f; // rendered arc at the in-flight segment's start

  auto shade_fragment = [&](const Vector &position, Fragment &fragment) {
    HS_MSP_STALL_START(shade_start);
    HS_MSP_COUNT(fragment_shader_calls);
    fragment_shader(position, fragment);
    HS_MSP_STALL_STOP(shade_palette, shade_start);
  };

  // Adaptively sub-step and plot one segment. `sample(t)` returns the sphere
  // point AND unit tangent at arc fraction t in [0,1] under the chosen strategy,
  // `sample.pos(t)` the point alone; `total_dist` is the segment's on-sphere
  // length (radians). Endpoints are omitted on interior / closed segments so a
  // shared vertex isn't plotted twice.
  auto process_segment = [&](auto &&sample, const Fragment &curr,
                             const Fragment &next, float total_dist,
                             bool is_last_segment) {
    // Rewrite the arc registers from the rendered arc when a planar basis is in
    // force (see the pre-pass above): `d` is the arc drawn so far within this
    // segment, `seg_base` the arc at its start. No-op for geodesic polylines.
    auto set_arc_uv = [&](Fragment &f, float d) {
      if constexpr (!DerivePlanarArcRegisters)
        return;
      if (!has_planar_basis)
        return;
      float arc = seg_base + d;
      f.v1 = arc;
      if (total_arc > math::EPS_GEOMETRIC)
        f.v0 = arc / total_arc;
    };
    // The degenerate and fast paths plot curr.pos/next.pos directly (original
    // sampled vertices), without the DRAWING PHASE renormalize that corrects
    // sample().pos's ~0.04% drift. Precondition: callers pass unit fragment
    // positions; the ~4e-6 an angle-addition vertex recurrence
    // (Star::sample_positions) leaves is two orders inside that drift and is
    // plotted as-is.
    // Degenerate (coincident endpoints): plot at most a single dot.
    if (total_dist < math::EPS_GEOMETRIC) {
      bool should_omit = close_loop || !is_last_segment || omit_end;
      if (!should_omit) {
        Fragment f_copy;
        if constexpr (InterpolateRegisters)
          f_copy = curr;
        f_copy.pos = curr.pos;
        f_copy.color = Color4(0, 0, 0, 0);
        set_arc_uv(f_copy, 0.0f);

        HS_PLOT_COUNT(shader_calls);
        shade_fragment(curr.pos, f_copy);
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, curr.pos, f_copy.color.color, f_copy.age,
                      f_copy.color.alpha);
      }
      return;
    }

    // Sub-step length at the segment start (also the first simulation step).
    const float base_step = (2.0f * PI_F) / W;
    auto balanced_step = [&](float default_step) {
      const float POLE_GUARD =
          base_step * MIN_POLE_SCALE * BALANCED_POLE_GUARD_SCALE;
      return default_step <= POLE_GUARD
                 ? default_step
                 : std::min(base_step, default_step * (BALANCED_SCREEN_STEP_PX /
                                                       SCREEN_STEP_PX));
    };
    auto adaptive_step = [&](const SamplePT &value) {
#ifdef HS_TEST_BUILD
      if (g_reference_screen_step)
        return screen_step_reference<W, H>(value.pos, value.tan, base_step);
#endif
      return screen_step<W, H>(value.pos, value.tan, base_step);
    };
    int planar_arc_interval = 0;
    auto adaptive_sample = [&](float t) -> SamplePT {
      HS_MSP_STALL_START(adaptive_start);
      SamplePT result;
      if constexpr (SinglePass && !DerivePlanarArcRegisters &&
                    !InterpolateRegisters && requires {
                      sample.one_pass_monotonic(t, planar_arc_interval);
                    }) {
        result = sample.one_pass_monotonic(t, planar_arc_interval);
#ifdef HS_TEST_BUILD
        ++g_planar_full_samples;
#endif
      } else if constexpr (SinglePass && requires { sample.one_pass(t); })
        result = sample.one_pass(t);
      else
        result = sample(t);
      HS_MSP_COUNT(adaptive_samples);
      HS_MSP_STALL_STOP(adaptive_sim, adaptive_start);
      return result;
    };
    HS_PLOT_COUNT(sim_samples);
    SamplePT smp = adaptive_sample(0.0f);
    float first_step = adaptive_step(smp);

    // FAST PATH: the whole segment spans ≤ one screen step, so a single dot
    // covers it. Keyed on SCREEN length, not arc length: a base_step arc can
    // still cross several pixels on a steep/near-polar segment, which an
    // arc-length test would undersample into a beaded line.
    if (total_dist <= first_step) {
      HS_PLOT_COUNT(one_dot);
      Fragment f;
      if constexpr (InterpolateRegisters)
        f = curr;
      f.pos = curr.pos;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, 0.0f);
      HS_PLOT_COUNT(shader_calls);
      shade_fragment(curr.pos, f);
      HS_PLOT_COUNT(plotted_samples);
      pipeline.plot(canvas, curr.pos, f.color.color, f.age, f.color.alpha);
      if (!close_loop && is_last_segment && !omit_end) {
        Fragment fl;
        if constexpr (InterpolateRegisters)
          fl = next;
        fl.pos = next.pos;
        fl.color = Color4(0, 0, 0, 0);
        set_arc_uv(fl, total_dist);
        HS_PLOT_COUNT(shader_calls);
        shade_fragment(next.pos, fl);
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, next.pos, fl.color.color, fl.age, fl.color.alpha);
      }
      return;
    }

    // Size each sub-step so consecutive samples land ~SCREEN_STEP_PX apart in
    // screen space. `smp`/`first_step` above seed the first iteration.
    if constexpr (SinglePass) {
      HS_PROFILE_DEEP(plot_seg_single_pass);
      float current_dist = 0.0f;
      float current_t = 0.0f;
      float desired_step = first_step;
      float default_desired_step = first_step;
      float previous_full_step = first_step;
      Vector previous_full_tangent = smp.tan;
      bool reuse_step = false;
      if constexpr (SamplingPolicy != RasterSamplingPolicy::DEFAULT) {
        if (balanced_sampling)
          desired_step = balanced_step(first_step);
      }
      size_t step_count = 0;
      float backstop_stretch = 1.0f;
      while (current_dist < total_dist) {
        Vector p;
        if constexpr (OpenGeodesic) {
          HS_PLOT_COUNT(normalizations);
#ifdef HS_TEST_BUILD
          if (g_reference_screen_step) {
            p = smp.pos.normalized();
          } else
#endif
          {
            // One Newton correction leaves only second-order length error.
            const float norm2 = dot(smp.pos, smp.pos);
            p = smp.pos * (1.5f - 0.5f * norm2);
          }
        } else if constexpr (SamplingPolicy != RasterSamplingPolicy::DEFAULT &&
                             requires { sample.one_pass(current_t); }) {
          HS_PLOT_COUNT(normalizations);
          if (balanced_sampling) {
            const float norm2 = dot(smp.pos, smp.pos);
            p = smp.pos * (1.5f - 0.5f * norm2);
          } else {
            p = smp.pos.normalized();
          }
        } else {
          HS_PLOT_COUNT(normalizations);
          p = smp.pos.normalized();
        }
        Fragment f;
        if constexpr (InterpolateRegisters)
          f = Fragment::lerp_registers(curr, next, current_t);
        f.pos = p;
        f.color = Color4(0, 0, 0, 0);
        set_arc_uv(f, current_dist);
        HS_PLOT_COUNT(shader_calls);
        shade_fragment(p, f);
        if constexpr (SamplingPolicy != RasterSamplingPolicy::DEFAULT) {
          if (balanced_sampling) {
            const float alpha_scale = desired_step / default_desired_step;
            f.color.alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
          }
        }
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, p, f.color.color, f.age, f.color.alpha);

        if (++step_count >= max_cache) {
          // Stretch factor matches the two-pass replay's; the hard stop below
          // bounds the extra steps it can cost.
          if (backstop_stretch == 1.0f) {
            HS_PLOT_COUNT(backstops);
            HS_SCAN_METRIC(hs::g_scan_metrics.plot_backstop_hits++);
            backstop_stretch = total_dist / current_dist;
          } else if (step_count >= 2 * max_cache) {
            break;
          }
        }
        HS_PLOT_MAX(steps_peak, step_count);
        float remaining = total_dist - current_dist;
        if (remaining <= desired_step) {
          current_dist = total_dist;
        } else if (remaining < 2.0f * desired_step) {
          current_dist += remaining * 0.5f;
        } else {
          current_dist += desired_step;
        }
        if (current_dist < total_dist) {
          current_t = current_dist / total_dist;
          HS_PLOT_COUNT(sim_samples);
          if constexpr (SamplingPolicy != RasterSamplingPolicy::DEFAULT &&
                        requires {
                          sample.position_monotonic(current_t,
                                                    planar_arc_interval);
                        }) {
            if (balanced_sampling && reuse_step) {
              HS_MSP_STALL_START(position_start);
              smp.pos =
                  sample.position_monotonic(current_t, planar_arc_interval);
              HS_MSP_COUNT(adaptive_samples);
              HS_MSP_STALL_STOP(adaptive_sim, position_start);
#ifdef HS_TEST_BUILD
              ++g_planar_position_samples;
#endif
              reuse_step = false;
            } else {
              smp = adaptive_sample(current_t);
              default_desired_step = adaptive_step(smp);
              if (balanced_sampling) {
                const float sin2 = 1.0f - smp.pos.y * smp.pos.y;
                reuse_step =
                    sin2 > 0.12f &&
                    default_desired_step > base_step * MIN_POLE_SCALE *
                                               BALANCED_POLE_GUARD_SCALE &&
                    default_desired_step < base_step * 0.9f &&
                    dot(smp.tan, previous_full_tangent) > 0.995f &&
                    fabsf(default_desired_step - previous_full_step) <
                        default_desired_step * 0.1f;
                previous_full_step = default_desired_step;
                previous_full_tangent = smp.tan;
              }
            }
          } else {
            smp = adaptive_sample(current_t);
            default_desired_step = adaptive_step(smp);
          }
          desired_step = default_desired_step;
          if constexpr (SamplingPolicy != RasterSamplingPolicy::DEFAULT) {
            if (balanced_sampling)
              desired_step = balanced_step(default_desired_step);
          }
          desired_step *= backstop_stretch;
        }
      }
      if (!close_loop && is_last_segment && !omit_end) {
        Fragment f;
        if constexpr (InterpolateRegisters)
          f = next;
        f.pos = next.pos;
        f.color = Color4(0, 0, 0, 0);
        set_arc_uv(f, total_dist);
        HS_PLOT_COUNT(shader_calls);
        shade_fragment(next.pos, f);
        if constexpr (SamplingPolicy != RasterSamplingPolicy::DEFAULT) {
          if (balanced_sampling) {
            const float alpha_scale = desired_step / default_desired_step;
            f.color.alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
          }
        }
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, next.pos, f.color.color, f.age, f.color.alpha);
      }
      return;
    }

    steps_cache.clear();
    float sim_dist = 0.0f;

    {
      HS_PROFILE_DEEP(plot_seg_sim);
      while (sim_dist < total_dist) {
        float step = steps_cache.is_empty() ? first_step : adaptive_step(smp);

        // Backstop: a pathological segment could exceed the 2*W cache. Stop
        // subdividing and let the normalized replay stretch the cached steps
        // over the rest of the segment (coarser sampling on an extreme arc is
        // fine).
        if (steps_cache.size() >= steps_cache.capacity()) {
          HS_PLOT_COUNT(backstops);
          HS_SCAN_METRIC(hs::g_scan_metrics.plot_backstop_hits++);
          break;
        }
        steps_cache.push_back(step);
        HS_PLOT_MAX(steps_peak, steps_cache.size());
        sim_dist += step;

        if (sim_dist < total_dist) {
          HS_PLOT_COUNT(sim_samples);
          smp = adaptive_sample(sim_dist / total_dist);
        }
      }
    }

    // The final step normally overshoots total_dist (scale <= 1) and the
    // normalized replay stretches the cached steps back to exactly total_dist.
    // On the backstop break path sim_dist can fall short (scale > 1) and the
    // replay stretches over the remaining segment instead.
    HS_CHECK(sim_dist > 0.0f);
    float scale = total_dist / sim_dist;
    bool omit_last = close_loop || !is_last_segment || omit_end;

    // DRAWING PHASE
    //
    // sample().pos is ~0.04% non-unit; vector_to_pixel's phi = acos(v.y) offsets
    // the row near the pole, so re-normalize the interpolated positions.
    HS_PROFILE_DEEP(plot_seg_draw);
    {
      HS_MSP_STALL_START(replay_start);
      HS_PLOT_COUNT(replay_samples);
      HS_PLOT_COUNT(normalizations);
      Vector start_pos = sample.pos(0.0f).normalized();
      Fragment f;
      if constexpr (InterpolateRegisters)
        f = Fragment::lerp_registers(curr, next, 0.0f);
      f.pos = start_pos;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, 0.0f);
      HS_MSP_STALL_STOP(normalized_replay, replay_start);

      HS_PLOT_COUNT(shader_calls);
      shade_fragment(start_pos, f);
      HS_PLOT_COUNT(plotted_samples);
      pipeline.plot(canvas, start_pos, f.color.color, f.age, f.color.alpha);
    }

    size_t loop_limit = omit_last ? steps_cache.size() - 1 : steps_cache.size();
    float current_dist = 0.0f;

    for (size_t j = 0; j < loop_limit; j++) {
      float step = steps_cache[j] * scale;
      current_dist += step;

      // total_dist > 0 here (HS_CHECK(sim_dist > 0) implies >=1 sim step).
      float t = current_dist / total_dist;

      // `t` (hence the drawn POSITION) is parameterized by the RENDERED arc
      // length. Registers are lerped from the control points; under a planar
      // basis set_arc_uv rewrites v0/v1 from the true rendered arc so a shader
      // keying off them as an arc-length proxy tracks the drawn position across
      // the planar bow. Geodesic edges keep the lerped registers.
      HS_MSP_STALL_START(replay_start);
      HS_PLOT_COUNT(replay_samples);
      HS_PLOT_COUNT(normalizations);
      Vector p = sample.pos(t).normalized();
      Fragment f;
      if constexpr (InterpolateRegisters)
        f = Fragment::lerp_registers(curr, next, t);
      f.pos = p;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, current_dist);
      HS_MSP_STALL_STOP(normalized_replay, replay_start);

      HS_PLOT_COUNT(shader_calls);
      shade_fragment(p, f);
      HS_PLOT_COUNT(plotted_samples);
      pipeline.plot(canvas, p, f.color.color, f.age, f.color.alpha);
    }
  };

  const auto &cr = canvas.clip();
  const bool clip_active = !cr.is_full();
  const auto xc = cr.x_clip();

  // Emits one shader-run dot for points[k]; the precomputed projection is
  // consumed only when no world stage would lift it back to a world vector.
  auto plot_dot = [&](const Fragment &src, size_t k) {
    Fragment f;
    if constexpr (InterpolateRegisters)
      f = src;
    f.pos = src.pos;
    f.color = Color4(0, 0, 0, 0);
    HS_PLOT_COUNT(shader_calls);
    shade_fragment(src.pos, f);
    if constexpr (pipeline_hoistable_projection<PipelineT>()) {
      if (point_rows != nullptr && point_cols != nullptr) {
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, point_cols[k], point_rows[k], f.color.color,
                      f.age, f.color.alpha);
        return;
      }
    }
    HS_PLOT_COUNT(plotted_samples);
    pipeline.plot(canvas, src.pos, f.color.color, f.age, f.color.alpha);
  };

  for (size_t i = 0; i < count; i++) {
    const Fragment &curr = points[i];
    const Fragment &next = segment_next(i);
    bool is_last_segment = (i == count - 1);
    PlanarEdgeSpan planar_cull_span;
    Vector planar_cull_end;
    bool reuse_planar_cull_samples = false;

    // --- Interpolation Strategy Selection ---
    // Branch-cut guard: the planar projection is singular at the basis antipode,
    // so a segment with an endpoint there falls back to a geodesic edge.
    bool antipodal_seam = false;
    if (has_planar_basis) {
      antipodal_seam =
          override_uv
              ? seg_seam_cache[i] != 0
              : dot(curr.pos, planar_basis->v) < -COS_PLANAR_ANTIPODE ||
                    dot(next.pos, planar_basis->v) < -COS_PLANAR_ANTIPODE;
    }
    const bool use_planar = planar_basis && !antipodal_seam;

    // Advance the rendered-arc accumulator for EVERY segment (drawn or culled) so
    // v0/v1 stay a true full-curve parameterization; seg_base snapshots the start
    // for the draw lambda. Skipped for geodesic polylines.
    if (override_uv) {
      seg_base = cumul;
      cumul += seg_arc_cache[i];
    }

    // Tier 3: Segment culling — skip if the edge's rendered row/column reach
    // (arc bulge included) lies outside the clip band; precomputed bits replace
    // the evaluation when the producer already ran the same predicate.
    if (clip_active) {
      HS_PLOT_COUNT(cull_tests);
      HS_PROFILE_DEEP(plot_seg_cull);
      bool visible;
      if constexpr (REUSE_PLANAR_CULL_SAMPLES) {
        bool rebuild_planar_sampler = false;
#ifdef HS_TEST_BUILD
        rebuild_planar_sampler = opts.rebuild_planar_sampler;
#endif
        if (edge_visible != nullptr) {
          visible = (edge_visible[i] & EDGE_VISIBLE) != 0;
        } else if (use_planar && xc.active && !rebuild_planar_sampler) {
          planar_cull_span =
              make_planar_edge_span(curr.pos, next.pos, *planar_basis);
          visible = planar_edge_visible_in_clip<W, H>(
              cr, xc, curr.pos, next.pos, *planar_basis, planar_cull_span,
              &planar_cull_end);
          reuse_planar_cull_samples = visible;
        } else {
          visible =
              edge_visible_in_clip<W, H>(pipeline, cr, xc, curr.pos, next.pos,
                                         use_planar ? planar_basis : nullptr);
        }
      } else {
        visible = edge_visible != nullptr
                      ? (edge_visible[i] & EDGE_VISIBLE) != 0
                      : edge_visible_in_clip<W, H>(
                            pipeline, cr, xc, curr.pos, next.pos,
                            use_planar ? planar_basis : nullptr);
      }
      if (!visible) {
        HS_PLOT_COUNT(culled);
        continue;
      }
    }

    // Single-dot shortcut: an edge proven to span <= one screen step renders
    // exactly as process_segment's fast path (set_arc_uv is a no-op without a
    // planar basis), so plot it without building the sampler. A predicate
    // false negative falls through and re-evaluates exactly.
    const bool one_dot =
        !has_planar_basis &&
        (edge_visible != nullptr && (edge_visible[i] & EDGE_CLASSIFIED) != 0
             ? (edge_visible[i] & EDGE_ONE_DOT) != 0
             : edge_fits_one_dot<W, H>(curr.pos, next.pos));
    if (one_dot) {
      HS_PLOT_COUNT(one_dot);
      plot_dot(curr, i);
      if (!close_loop && is_last_segment && !omit_end)
        plot_dot(next, i + 1);
      continue;
    }

    if (use_planar) {
      HS_PLOT_COUNT(planar);
      if constexpr (REUSE_PLANAR_CULL_SAMPLES) {
        PlanarEdgeSampler sampler;
        if (reuse_planar_cull_samples) {
          sampler = make_planar_edge_sampler(planar_cull_span, planar_cull_end,
                                             *planar_basis);
        } else {
          sampler = make_planar_edge_sampler(curr.pos, next.pos, *planar_basis);
        }
        process_segment(sampler, curr, next, sampler.dist, is_last_segment);
      } else {
        rasterize_planar_strategy(curr, next, *planar_basis, is_last_segment,
                                  process_segment);
      }
    } else {
      HS_PLOT_COUNT(geodesic);
      rasterize_geodesic_strategy(curr, next, is_last_segment, process_segment);
    }
  }
}
HS_O3_END

/**
 * @brief Per-primitive geometry/rasterization options for draw_fragments.
 *
 * `close_loop` and `planar_basis` default to the common (geodesic, open) case,
 * so most primitives only spell out `.capacity`.
 */
struct FragmentDrawParams {
  size_t capacity;         /**< Fragment buffer reservation (per-primitive). */
  bool close_loop = false; /**< Passed to rasterize (closes last→first edge). */
  bool omit_end =
      false; /**< Skip the final endpoint plot of an open line, so abutting arcs tile without a double-plot. */
  const Basis *planar_basis =
      nullptr; /**< Planar projection basis (null = geodesic). */
  Fragment *loop_seam =
      nullptr; /**< Optional closing target with seam-specific registers. */
};

/**
 * @brief Run the shared per-primitive draw ritual.
 *
 * Every Plot primitive opens a ScratchScope, binds a Fragments buffer, fills it,
 * applies the optional vertex shader, and rasterizes. The ScratchScope must
 * outlive the rasterize call (the arena backs the fragments).
 *
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam FillFn Callable (Fragments &) -> void supplying the primitive's
 *                sampling.
 * @param pipeline Render pipeline.
 * @param canvas Target canvas.
 * @param vertex_shader Optional per-vertex shader.
 * @param fragment_shader Per-fragment shader.
 * @param params Per-primitive capacity / close-loop / planar-basis options.
 * @param fill Fills the bound Fragments buffer with the primitive's samples.
 */
template <int W, int H, typename FillFn>
inline void draw_fragments(PipelineRef pipeline, Canvas &canvas,
                           VertexShaderRef vertex_shader,
                           FragmentShaderFn fragment_shader,
                           const FragmentDrawParams &params, FillFn &&fill) {
  ScratchScope frag_guard(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, params.capacity);
  fill(points);
  apply_vertex_shader(vertex_shader, points);
  if (params.loop_seam != nullptr && vertex_shader)
    vertex_shader(*params.loop_seam);
  rasterize<W, H>(pipeline, canvas, points, fragment_shader,
                  {.close_loop = params.close_loop,
                   .planar_basis = params.planar_basis,
                   .omit_end = params.omit_end,
                   .loop_seam = params.loop_seam});
}

/**
 * @brief Draws a geodesic line between two points.
 * Registers:
 *  v0: line progress (0..1)
 *  v1: Arc Length (radians)
 */
struct Line {
  /**
   * @brief One fragment at arc fraction t on a geodesic edge.
   * @param f1 Start fragment.
   * @param f2 End fragment.
   * @param es Shared setup from make_geodesic_edge_span(f1.pos, f2.pos); must
   *        have an axis.
   * @param perp cross(es.axis, f1.pos), the arc's start tangent direction.
   * @param t Arc fraction in [0, 1]; the endpoints reproduce f1/f2 exactly.
   * @return The interpolated fragment, registers included.
   */
  static Fragment sample_point(const Fragment &f1, const Fragment &f2,
                               const GeodesicEdgeSpan &es, const Vector &perp,
                               float t) {
    Fragment f = Fragment::lerp(f1, f2, t);
    if (t <= 0.0f)
      f.pos = f1.pos;
    else if (t >= 1.0f)
      f.pos = f2.pos;
    else {
      // fast trig's 0.17% error breaks c^2+s^2==1, so renormalize: callers
      // (vector_to_pixel's acos) require a unit position.
      float s, c;
      fast_sincosf_0_pi(es.total * t, s, c);
      Vector p = (f1.pos * c) + (perp * s);
      HS_PLOT_COUNT(normalizations);
      f.pos = p * fast_rsqrt(dot(p, p));
    }

    f.v0 = t;
    f.v1 = es.total * t;
    f.v2 = 0.0f;
    return f;
  }

  /**
   * @brief Samples a geodesic line between two points.
   * @param points Output fragment list; density+1 fragments are appended.
   * @param f1 Start fragment.
   * @param f2 End fragment.
   * @param density Number of sub-segments (>=1); the line is sampled at
   *                density+1 evenly-parameterized points.
   */
  static void sample(Fragments &points, const Fragment &f1, const Fragment &f2,
                     int density = 1) {
    if (density < 1)
      density = 1;

    const GeodesicEdgeSpan es = make_geodesic_edge_span(f1.pos, f2.pos);
    if (!es.have_axis) {
      Fragment f = f1;
      f.v0 = f.v1 = f.v2 = 0.0f;
      points.push_back(f);
      points.push_back(f); // draw at least a dot
      return;
    }

    const Vector perp = cross(es.axis, f1.pos);
    for (int i = 0; i <= density; ++i)
      points.push_back(
          sample_point(f1, f2, es, perp, static_cast<float>(i) / density));
  }

  /**
   * @brief Draws a geodesic line.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param f1 Start fragment.
   * @param f2 End fragment.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Fragment &f1,
                   const Fragment &f2, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = 4},
                         [&](Fragments &points) { sample(points, f1, f2); });
  }

  /**
   * @brief Draws a geodesic line without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param f1 Start fragment.
   * @param f2 End fragment.
   * @param fragment_shader Shader function.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Fragment &f1,
                   const Fragment &f2, FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, f1, f2, fragment_shader, {});
  }
};

/**
 * @brief Multiline primitive (Polyline).
 * Registers:
 *  v0: Path Progress (0.0 -> 1.0)
 *  v1: Cumulative Arc Length (radians) — geodesic chord-polygon length
 *  v2: Vertex Index
 * @note v0/v1 accumulate the GEODESIC distance between consecutive control
 *       points (the rendered arc under Multiline's geodesic edges). With a
 *       `planar_basis` the rasterizer re-derives v0/v1 from the longer
 *       azimuthal-equidistant arc, so the registers track the rendered position
 *       either way.
 */
struct Multiline {
  /**
   * @brief Samples a multiline path from a container of vertices.
   * @param points Output fragment list; one arc-length-parameterized fragment
   *               is appended per input vertex.
   * @param vertices Iterable container of Fragment.
   * @param closed If true, connects the last point to the first.
   * @return Closing target with continued seam registers when closed; a
   *         default Fragment otherwise.
   */
  static Fragment sample(Fragments &points, const auto &vertices,
                         bool closed = false) {
    auto it = std::begin(vertices);
    auto end = std::end(vertices);

    if (it == end)
      return {};

    float total_len = 0.0f;
    Fragment first = *it;
    Fragment prev = first;
    auto len_it = it;
    ++len_it;

    for (; len_it != end; ++len_it) {
      const Fragment &curr = *len_it;
      total_len += angle_between(prev.pos, curr.pos);
      prev = curr;
    }
    if (closed) {
      total_len += angle_between(prev.pos, first.pos);
    }

    if (total_len < math::EPS_GEOMETRIC) {
      // Avoid divide-by-zero on a degenerate path; v0 collapses toward 0, but the
      // path is geometrically a point so the lost progress parameter is moot.
      total_len = 1.0f;
    }

    float current_len = 0.0f;
    it = std::begin(vertices);
    prev = *it;

    Fragment f = prev;
    f.v0 = 0.0f;
    f.v1 = 0.0f;
    f.v2 = 0.0f;
    points.push_back(f);

    ++it;
    int idx = 1;
    for (; it != end; ++it) {
      const Fragment &curr = *it;
      float dist = angle_between(prev.pos, curr.pos);
      current_len += dist;

      f = curr;
      f.v0 = current_len / total_len;
      f.v1 = current_len;
      f.v2 = static_cast<float>(idx++);
      points.push_back(f);
      prev = curr;
    }

    if (closed) {
      float dist = angle_between(prev.pos, first.pos);
      current_len += dist;
      f = first;
      f.v0 = 1.0f;
      f.v1 = current_len;
      f.v2 = static_cast<float>(idx);
      return f;
    }
    return {};
  }

  /**
   * @brief Draws a multiline path.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param vertices Iterable container of Fragment.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param closed If true, connects the last point to the first.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &vertices,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, bool closed = false) {
    Fragment loop_seam;
    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = vertices.size() + 1,
                          .close_loop = closed,
                          .loop_seam = closed ? &loop_seam : nullptr},
                         [&](Fragments &points) {
                           loop_seam = sample(points, vertices, closed);
                         });
  }

  /**
   * @brief Draws a multiline path without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param vertices Iterable container of Fragment.
   * @param fragment_shader Shader function.
   * @param closed If true, connects the last point to the first.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &vertices,
                   FragmentShaderFn fragment_shader, bool closed = false) {
    draw<W, H>(pipeline, canvas, vertices, fragment_shader, {}, closed);
  }
};

/**
 * @brief Samples a closed parametric ring of `num_verts` vertices and appends
 * the overlap-close vertex.
 * @tparam PosFn Callable (int i) -> Vector giving vertex i's unit-sphere position.
 * @param points Output fragment list; num_verts+1 fragments are appended.
 * @param num_verts Number of ring vertices (i in [0, num_verts)).
 * @param pos_fn Returns the unit-sphere position of vertex i.
 * @details Each vertex carries the standard ring registers — v0: perimeter
 * progress (i / num_verts), v1: accumulated great-circle arc length from vertex
 * 0, v2: vertex index, age: 0. The trailing close vertex duplicates vertex 0's
 * position with v0 = 1 and the arc length continued across the wrap edge, so an
 * `omit_end` rasterize draws the wrap edge without a UV seam and without
 * replotting vertex 0. Shared skeleton for the accumulated-arc closed rings
 * (Star, Flower, DistortedRing); Ring keeps its own analytic-arc loop. For the
 * PLANAR callers the rasterizer overrides v0/v1, so these geodesic values seed
 * only the optional vertex shader.
 */
template <typename PosFn>
inline void sample_closed_ring(Fragments &points, int num_verts, PosFn pos_fn) {
  HS_CHECK(num_verts >= 1);
  float cumulative_len = 0.0f;
  size_t start_idx = points.size();
  for (int i = 0; i < num_verts; i++) {
    Fragment f;
    f.pos = pos_fn(i);
    if (i > 0)
      cumulative_len += angle_between(points.back().pos, f.pos);
    f.v0 = static_cast<float>(i) / num_verts;
    f.v1 = cumulative_len;
    f.v2 = static_cast<float>(i);
    f.age = 0;
    points.push_back(f);
  }

  // Manual close (overlap): duplicate vertex 0 with continued arc length.
  Fragment last = points[start_idx];
  cumulative_len += angle_between(points.back().pos, last.pos);
  last.v0 = 1.0f;
  last.v1 = cumulative_len;
  last.v2 = static_cast<float>(num_verts);
  points.push_back(last);
}

/**
 * @brief Tangent-plane vector at LUT index i rotated by phase.
 * @details Angle-addition identity: cos/sin(θ+φ) from the precomputed θ-grid,
 * then (u·cos_t + w·sin_t). Shared by the LUT-optimized Ring and DistortedRing
 * samplers so a future LUT-recovery correction stays in one place.
 */
template <int W, int H>
static inline Vector ring_tangent(int i, const Vector &u, const Vector &w,
                                  float cos_phase, float sin_phase) {
  float cos_t = TrigLUT<W, H>::cos_theta(i) * cos_phase -
                TrigLUT<W, H>::sin_theta[i] * sin_phase;
  float sin_t = TrigLUT<W, H>::sin_theta[i] * cos_phase +
                TrigLUT<W, H>::cos_theta(i) * sin_phase;
  return (u * cos_t) + (w * sin_t);
}

/**
 * @brief LUT index stride that samples a ring at ~one control point per screen
 *        pixel of its circumference.
 * @tparam W Rasterization width; the LUT grid and the equatorial column count.
 * @param r_val sin of the ring's colatitude — its circumference is 2π·r_val.
 * @details An azimuth step of stride·2π/W spans stride·r_val·2π/W radians of
 * arc, so stride <= 1/r_val holds consecutive control points within the
 * rasterizer's base_step (2π/W) — the density a great circle already gets at
 * stride 1. The rasterizer sub-steps every segment to SCREEN_STEP_PX, so a
 * coarser grid moves control points without thinning rendered coverage.
 * r_val is floored at 1/MAX_STRIDE, keeping a sub-pixel ring at
 * MIN_RING_SAMPLES vertices and the reciprocal finite.
 */
template <int W> static inline int ring_lut_stride(float r_val) {
  constexpr int MIN_RING_SAMPLES = 8;
  constexpr int MAX_STRIDE = std::max(1, W / MIN_RING_SAMPLES);
  const float inv =
      1.0f / std::max(r_val, 1.0f / static_cast<float>(MAX_STRIDE));
  return std::max(1, static_cast<int>(inv));
}

/** @brief Antipode-folded working basis and the ring's colatitude trig. */
struct RingFrame {
  Basis basis;
  float theta_eq;  ///< Ring colatitude (radians).
  float sin_theta; ///< Tangent-plane radius of the ring.
  float cos_theta; ///< Ring offset along the pole axis.
};

/**
 * @brief Resolves the frame of a radius-(0-2) ring drawn on `basis`.
 * @param basis Orientation basis.
 * @param radius Angular radius (0-2).
 * @details Folds radius > 1 to the antipode, then derives the colatitude and
 * its sine/cosine. Every ring-family sampler and DistortedRing::fn_point
 * resolve their frame here; a divergence would detach sampled points from the
 * visible ring.
 */
inline RingFrame ring_frame(const Basis &basis, float radius) {
  auto res = get_antipode(basis, radius);
  const float theta_eq = res.second * (PI_F / 2.0f);
  return {res.first, theta_eq, sinf(theta_eq), cosf(theta_eq)};
}

/**
 * @brief Ring primitives.
 * Registers:
 *  v0: Angular progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Index
 */
struct Ring {
  /**
   * @brief Samples a closed ring at `num_samples` evenly-spaced angles.
   * @param points Output fragment list; num_samples+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Ring radius (radians).
   * @param num_samples Number of evenly-spaced samples around the ring.
   * @param phase Rotation phase (radians).
   * @details Runtime sample count for the polygon samplers, whose vertex counts
   * do not match the TrigLUT grid; appends an overlap-close vertex. v1 is the
   * analytic arc length (theta·sin(theta_eq), theta_eq being the ring's
   * colatitude).
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_samples, float phase = 0) {
    HS_CHECK(num_samples >= 1);
    const RingFrame frame = ring_frame(basis, radius);
    const Vector &v = frame.basis.v;
    const Vector &u = frame.basis.u;
    const Vector &w = frame.basis.w;
    const float r_val = frame.sin_theta;
    const float d_val = frame.cos_theta;

    const float step = 2.0f * PI_F / num_samples;

    Vector first_pos; // Reused for the overlap-close vertex (i == 0 position).
    for (int i = 0; i < num_samples; i++) {
      float theta = i * step;
      float t = theta + phase;
      Vector u_temp = (u * cosf(t)) + (w * sinf(t));

      Fragment f;
      HS_PLOT_COUNT(normalizations);
      f.pos = ((v * d_val) + (u_temp * r_val)).normalized();
      if (i == 0)
        first_pos = f.pos;
      f.v0 = static_cast<float>(i) / num_samples;
      f.v1 = theta * r_val;
      f.v2 = static_cast<float>(i);
      f.age = 0;

      points.push_back(f);
    }

    // Manual Close (Overlap): the close vertex at theta == 2π has the same
    // position as the i == 0 vertex; reusing it also avoids the float error a
    // literal 2π+φ argument introduces.
    Fragment f;
    f.pos = first_pos;
    f.v0 = 1.0f;
    f.v1 = 2.0f * PI_F * r_val;
    f.v2 = static_cast<float>(num_samples);
    f.age = 0;
    points.push_back(f);
  }

  /**
   * @brief Closed ring on the LUT angle grid — LUT-optimized.
   * @tparam W,H Rasterization resolution; W is the LUT grid.
   * @param points Output fragment list; ceil(W/stride)+1 fragments are appended,
   *               at most W+1.
   * @param basis Orientation basis.
   * @param radius Ring radius (radians).
   * @param phase Rotation phase (radians).
   * @details The angle grid (i*2π/W) is exactly TrigLUT<W,H>::cos_theta and
   * sin_theta, so per-sample cosf(θ+φ)/sinf(θ+φ) becomes
   * the precomputed θ-grid plus one angle-addition against cos/sin(φ), saving
   * ~2*(W+1) libm trig calls per ring per frame. ring_lut_stride skips grid
   * indices a ring narrower than the equator cannot resolve; v0/v1/v2 stay keyed
   * to the grid index, so the analytic arc length and overlap close are
   * unchanged. The runtime int-num_samples overload stays for the polygon
   * samplers, whose vertex counts do not match the LUT grid.
   */
  template <int W, int H>
  static void sample(Fragments &points, const Basis &basis, float radius,
                     float phase = 0) {
    const RingFrame frame = ring_frame(basis, radius);
    const Vector &v = frame.basis.v;
    const Vector &u = frame.basis.u;
    const Vector &w = frame.basis.w;
    const float r_val = frame.sin_theta;
    const float d_val = frame.cos_theta;

    const float step = 2.0f * PI_F / W;
    const int stride = ring_lut_stride<W>(r_val);

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    const float cos_phase = cosf(phase);
    const float sin_phase = sinf(phase);

    for (int i = 0; i < W; i += stride) {
      Vector u_temp = ring_tangent<W, H>(i, u, w, cos_phase, sin_phase);

      Fragment f;
      HS_PLOT_COUNT(normalizations);
      f.pos = ((v * d_val) + (u_temp * r_val)).normalized();
      f.v0 = static_cast<float>(i) / W;
      f.v1 = (i * step) * r_val;
      f.v2 = static_cast<float>(i);
      f.age = 0;

      points.push_back(f);
    }

    // Manual Close (Overlap): θ = 2π folds to (cos φ, sin φ) by periodicity.
    Fragment f;
    Vector u_temp = (u * cos_phase) + (w * sin_phase);
    HS_PLOT_COUNT(normalizations);
    f.pos = ((v * d_val) + (u_temp * r_val)).normalized();
    f.v0 = 1.0f;
    f.v1 = (2.0f * PI_F) * r_val;
    f.v2 = static_cast<float>(W);
    f.age = 0;
    points.push_back(f);
  }

  /**
   * @brief Draws a ring.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Ring radius (radians).
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    draw_fragments<W, H>(
        pipeline, canvas, vertex_shader, fragment_shader,
        {.capacity = W + 2, .omit_end = true},
        [&](Fragments &points) { sample<W, H>(points, basis, radius, phase); });
  }

  /**
   * @brief Draws a ring without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Ring radius (radians).
   * @param fragment_shader Shader function.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, FragmentShaderFn fragment_shader,
                   float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, fragment_shader, {}, phase);
  }
};

/**
 * @brief Azimuthal-equidistant edge projection policy.
 * @note Always renders with PLANAR (azimuthal-equidistant) edges, which bow
 *       LONGER than the great-circle chord. The rasterizer re-derives v0/v1 from
 *       that true rendered arc, so both track the drawn position rather than the
 *       shorter chord polygon.
 */
struct PlanarProjection {
  static const Basis *edge_basis(const Basis &basis, float radius,
                                 Basis &storage) {
    storage = radius > 1.0f ? planar_chart_basis(-basis.v) : basis;
    return &storage;
  }

  static void finish_polygon_sample(Fragments &, size_t) {}
};

/** @brief Great-circle edge projection. */
struct GeodesicProjection {
  static const Basis *edge_basis(const Basis &, float, Basis &) {
    return nullptr;
  }

  static void finish_polygon_sample(Fragments &points, size_t start_idx) {
    float cumulative_length = 0.0f;
    for (size_t i = start_idx; i < points.size(); ++i) {
      if (i > start_idx)
        cumulative_length += angle_between(points[i - 1].pos, points[i].pos);
      points[i].v1 = cumulative_length;
    }
  }
};

/** @brief Polygon with projection-selected edges. */
template <typename Projection> struct Polygon {
  /**
   * @brief Samples a polygon.
   * @param points Output fragment list; num_sides+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Polygon radius.
   * @param num_sides Number of sides.
   * @param phase Rotation phase (radians).
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    HS_CHECK(num_sides >= 1);
    const size_t start_idx = points.size();
    Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);
    Projection::finish_polygon_sample(points, start_idx);
  }

  /**
   * @brief Draws a polygon.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Polygon radius.
   * @param num_sides Number of sides.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    Basis projection_basis;
    const Basis *edge_basis =
        Projection::edge_basis(basis, radius, projection_basis);

    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides + 2),
                          .omit_end = true,
                          .planar_basis = edge_basis},
                         [&](Fragments &points) {
                           sample(points, basis, radius, num_sides, phase);
                         });
  }

  /**
   * @brief Draws a polygon without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Polygon radius.
   * @param num_sides Number of sides.
   * @param fragment_shader Shader function.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, {},
               phase);
  }
};

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

/**
 * @brief Distorted Ring.
 * Registers:
 *  v0: Angular progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Index
 */
struct DistortedRing {
  /**
   * @brief Calculates a single point on a distorted ring.
   * @param shift_fn Radial distortion sampled at angle/(2π).
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param angle Angular position around the ring (radians).
   * @return Normalized unit sphere point on the distorted ring.
   * @details Same geometry as sample() at phase 0, so the returned point lands
   * on the drawn ring; any divergence would detach callers' sampled points from
   * the visible ring off Radius=1. Computed with direct cosf/sinf where
   * sample() uses TrigLUT angle addition, so the two agree in exact arithmetic
   * but are not bit-identical. There is no phase parameter: a ring drawn with a
   * non-zero phase is rotated away from the returned point.
   */
  static Vector fn_point(ScalarFn shift_fn, const Basis &basis, float radius,
                         float angle) {
    const RingFrame frame = ring_frame(basis, radius);
    const Vector &v = frame.basis.v;
    const Vector &u = frame.basis.u;
    const Vector &w = frame.basis.w;

    const float polar = frame.theta_eq + shift_fn(angle / (2.0f * PI_F));
    Vector u_temp = (u * cosf(angle)) + (w * sinf(angle));
    HS_PLOT_COUNT(normalizations);
    return ((v * cosf(polar)) + (u_temp * sinf(polar))).normalized();
  }

  template <int W, int H>
  /**
   * @brief Samples a distorted ring.
   * @tparam W,H Rasterization resolution (drives the W-sample count and LUT).
   * @param points Output fragment list; W+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param shift_fn Radial distortion sampled per vertex.
   * @param phase Rotation phase (radians).
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     ScalarFn shift_fn, float phase = 0) {
    const RingFrame frame = ring_frame(basis, radius);
    const Vector &v = frame.basis.v;
    const Vector &u = frame.basis.u;
    const Vector &w = frame.basis.w;
    const float r_val = frame.sin_theta;
    const float d_val = frame.cos_theta;

    const int num_samples = W;
    const float step = 2.0f * PI_F / num_samples;

    // Precompute phase for angle-addition: cos/sin(θ+φ) via TrigLUT
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    const float cos_phase = cosf(phase);
    const float sin_phase = sinf(phase);

    // Per-vertex point carries the shift-fn radial distortion; the loop, arc-
    // length accumulation, and overlap close are the shared closed-ring skeleton.
    sample_closed_ring(points, num_samples, [&](int i) {
      float theta = i * step;
      Vector u_temp = ring_tangent<W, H>(i, u, w, cos_phase, sin_phase);

      float shift = shift_fn(theta / (2.0f * PI_F));
      float cos_shift = cosf(shift);
      float sin_shift = sinf(shift);

      float v_scale = d_val * cos_shift - r_val * sin_shift;
      float u_scale = r_val * cos_shift + d_val * sin_shift;

      HS_PLOT_COUNT(normalizations);
      return ((v * v_scale) + (u_temp * u_scale)).normalized();
    });
  }

  /**
   * @brief Draws a distorted ring.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param shift_fn Distortion function.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, ScalarFn shift_fn,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = W + 2, .omit_end = true},
                         [&](Fragments &points) {
                           sample<W, H>(points, basis, radius, shift_fn, phase);
                         });
  }

  /**
   * @brief Draws a distorted ring without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param shift_fn Distortion function.
   * @param fragment_shader Shader function.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, ScalarFn shift_fn,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, shift_fn, fragment_shader, {},
               phase);
  }
};

/**
 * @brief Star shape with projection-selected edges.
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Vertex index
 */
template <typename Projection> struct Star {
public:
  /** @brief Sine/cosine values shared by every vertex at one radius. */
  struct RadiusTrig {
    float sine[2];
    float cosine[2];
  };

  /** @brief Sine/cosine values shared by every angular step. */
  struct StepTrig {
    float sine;
    float cosine;
  };

private:
  static void sample_positions_impl(Fragments &points, const Basis &work_basis,
                                    int num_sides, float phase,
                                    const RadiusTrig &radius_trig,
                                    const StepTrig &step_trig) {
    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;
    const size_t start_idx = points.size();
    float cos_theta = cosf(phase);
    float sin_theta = sinf(phase);
    for (int i = 0; i < num_sides * 2; ++i) {
      Fragment f;
      const float sin_r = radius_trig.sine[i & 1];
      const float cos_r = radius_trig.cosine[i & 1];
      f.pos =
          (v * cos_r) + (u * (cos_theta * sin_r)) + (w * (sin_theta * sin_r));
      points.push_back(f);
      const float next_cos =
          cos_theta * step_trig.cosine - sin_theta * step_trig.sine;
      sin_theta = sin_theta * step_trig.cosine + cos_theta * step_trig.sine;
      cos_theta = next_cos;
    }
    points.push_back(points[start_idx]);
  }

  template <bool EmitRegisters>
  static void sample_impl(Fragments &points, const Basis &basis, float radius,
                          int num_sides, float phase) {
    HS_CHECK(num_sides >= 1);
    const RingFrame frame = ring_frame(basis, radius);
    const Vector &v = frame.basis.v;
    const Vector &u = frame.basis.u;
    const Vector &w = frame.basis.w;

    float inner_radius = frame.theta_eq * STAR_INNER_RATIO;
    float angle_step = PI_F / num_sides;
    const float sin_radius[2] = {frame.sin_theta, sinf(inner_radius)};
    const float cos_radius[2] = {frame.cos_theta, cosf(inner_radius)};

    if constexpr (EmitRegisters) {
      auto position = [&](int i) {
        float theta = phase + i * angle_step;
        float sin_r = sin_radius[i & 1];
        float cos_r = cos_radius[i & 1];
        float cos_t = cosf(theta);
        float sin_t = sinf(theta);
        Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
        HS_PLOT_COUNT(normalizations);
        p.normalize();
        return p;
      };
      sample_closed_ring(points, num_sides * 2, position);
    } else {
      const RadiusTrig radius_trig = {
          {sin_radius[0], sin_radius[1]},
          {cos_radius[0], cos_radius[1]},
      };
      const StepTrig step_trig = {sinf(angle_step), cosf(angle_step)};
      sample_positions_impl(points, frame.basis, num_sides, phase, radius_trig,
                            step_trig);
    }
  }

  template <bool EmitRegisters>
  static void sample_continuous_impl(Fragments &points, const Basis &basis,
                                     float radius, int num_sides, float phase) {
    HS_CHECK(num_sides >= 1);
    HS_CHECK(radius >= 0.0f && radius <= 2.0f);

    const float outer_radius = radius * (PI_F / 2.0f);
    const float inner_radius =
        radius <= 1.0f
            ? outer_radius * STAR_INNER_RATIO
            : STAR_INNER_RATIO * (PI_F / 2.0f) +
                  (radius - 1.0f) * (PI_F - STAR_INNER_RATIO * (PI_F / 2.0f));
    const float angle_step = PI_F / num_sides;
    const float sin_radius[2] = {sinf(outer_radius), sinf(inner_radius)};
    const float cos_radius[2] = {cosf(outer_radius), cosf(inner_radius)};

    auto position = [&](int i) {
      const float theta = phase + i * angle_step;
      const float sin_r = sin_radius[i & 1];
      const float cos_r = cos_radius[i & 1];
      Vector p = (basis.v * cos_r) + (basis.u * (cosf(theta) * sin_r)) +
                 (basis.w * (sinf(theta) * sin_r));
      HS_PLOT_COUNT(normalizations);
      p.normalize();
      return p;
    };

    if constexpr (EmitRegisters) {
      sample_closed_ring(points, num_sides * 2, position);
    } else {
      const size_t start_idx = points.size();
      for (int i = 0; i < num_sides * 2; ++i) {
        Fragment f;
        f.pos = position(i);
        points.push_back(f);
      }
      points.push_back(points[start_idx]);
    }
  }

public:
  /** @brief Computes reusable radius trigonometry for position-only sampling. */
  static RadiusTrig radius_trig(float radius) {
    const float work_radius = radius > 1.0f ? 2.0f - radius : radius;
    const float outer_radius = work_radius * (PI_F / 2.0f);
    const float inner_radius = outer_radius * STAR_INNER_RATIO;
    return {{sinf(outer_radius), sinf(inner_radius)},
            {cosf(outer_radius), cosf(inner_radius)}};
  }

  /** @brief Computes reusable angular-step trigonometry for a side count. */
  static StepTrig step_trig(int num_sides) {
    const float angle_step = PI_F / num_sides;
    return {sinf(angle_step), cosf(angle_step)};
  }

  /**
   * @brief Samples a star shape.
   * @param points Output fragment list; num_sides*2+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of points.
   * @param phase Rotation phase (radians).
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    sample_impl<true>(points, basis, radius, num_sides, phase);
  }

  /**
   * @brief Samples only star positions, leaving fragment registers at defaults.
   * @param points Output fragment list; num_sides*2+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of points.
   * @param phase Rotation phase (radians).
   * @note Steps the vertex angle by an angle-addition recurrence and skips the
   *   per-vertex normalize sample() does; positions come back off-unit by
   *   ~4e-6, inside what rasterize()'s fast paths plot verbatim.
   */
  static void sample_positions(Fragments &points, const Basis &basis,
                               float radius, int num_sides, float phase = 0) {
    sample_impl<false>(points, basis, radius, num_sides, phase);
  }

  /** @brief Samples positions with caller-cached trigonometric values. */
  static void sample_positions(Fragments &points, const Basis &basis,
                               float radius, int num_sides, float phase,
                               const RadiusTrig &radius_trig,
                               const StepTrig &step_trig) {
    HS_CHECK(num_sides >= 1);
    const Basis work_basis = get_antipode(basis, radius).first;
    sample_positions_impl(points, work_basis, num_sides, phase, radius_trig,
                          step_trig);
  }

  /** @brief Samples Star levels that continue to the opposite pole. */
  static void sample_continuous(Fragments &points, const Basis &basis,
                                float radius, int num_sides, float phase = 0) {
    sample_continuous_impl<true>(points, basis, radius, num_sides, phase);
  }

  /** @brief Samples Star levels that remain continuous across the equator. */
  static void sample_continuous_positions(Fragments &points, const Basis &basis,
                                          float radius, int num_sides,
                                          float phase = 0) {
    sample_continuous_impl<false>(points, basis, radius, num_sides, phase);
  }

  /**
   * @brief Draws a star.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of points.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    Basis projection_basis;
    const Basis *edge_basis =
        Projection::edge_basis(basis, radius, projection_basis);

    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides * 2 + 2),
                          .omit_end = true,
                          .planar_basis = edge_basis},
                         [&](Fragments &points) {
                           sample(points, basis, radius, num_sides, phase);
                         });
  }

  /** @brief Draws a Star level that continues across the equator. */
  template <int W, int H>
  static void draw_continuous(PipelineRef pipeline, Canvas &canvas,
                              const Basis &basis, float radius, int num_sides,
                              FragmentShaderFn fragment_shader,
                              float phase = 0) {
    Basis projection_basis;
    const Basis *edge_basis =
        Projection::edge_basis(basis, radius, projection_basis);

    draw_fragments<W, H>(pipeline, canvas, {}, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides * 2 + 2),
                          .omit_end = true,
                          .planar_basis = edge_basis},
                         [&](Fragments &points) {
                           sample_continuous(points, basis, radius, num_sides,
                                             phase);
                         });
  }

  /**
   * @brief Draws a star without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of points.
   * @param fragment_shader Shader function.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, {},
               phase);
  }
};

/**
 * @brief Flower shape.
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians) — cumulative rendered planar arc
 *  v2: Vertex index
 * @note Always renders with PLANAR (azimuthal-equidistant) edges, which bow
 *       LONGER than the great-circle chord. The rasterizer re-derives v0/v1 from
 *       that true rendered arc, so both track the drawn position rather than the
 *       shorter chord polygon.
 */
struct Flower {
  /**
   * @brief Samples a flower shape.
   * @param points Output fragment list; num_sides*2+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of petals.
   * @param phase Rotation phase (radians).
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    HS_CHECK(num_sides >= 1);
    const RingFrame frame = ring_frame(basis, radius);
    const Vector &v = frame.basis.v;
    const Vector &u = frame.basis.u;
    const Vector &w = frame.basis.w;

    float apothem = PI_F - frame.theta_eq;
    float safe_apothem = std::min(apothem, PI_F - 1e-4f);
    float angle_step = PI_F / num_sides;
    const float sin_r = sinf(safe_apothem);
    const float cos_r = cosf(safe_apothem);

    // Constant polar radius per vertex; everything else is the shared closed-
    // ring skeleton.
    sample_closed_ring(points, num_sides * 2, [&](int i) {
      float theta = phase + i * angle_step;
      float cos_t = cosf(theta);
      float sin_t = sinf(theta);
      Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
      HS_PLOT_COUNT(normalizations);
      p.normalize();
      return p;
    });
  }

  /**
   * @brief Draws a flower.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of petals.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    // Center the chart on the antipode pole, opposite the petal ring: projecting
    // the constant-radius ring through the far-pole chart bows its straight edges
    // outward into petals.
    Basis planar_basis =
        planar_chart_basis(get_antipode(basis, radius).first.v);

    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides * 2 + 2),
                          .omit_end = true,
                          .planar_basis = &planar_basis},
                         [&](Fragments &points) {
                           sample(points, basis, radius, num_sides, phase);
                         });
  }

  /**
   * @brief Draws a flower without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of petals.
   * @param fragment_shader Shader function.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, {},
               phase);
  }
};

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

/**
 * @brief Gates one geodesic trail's edges against the clip in one hoisted pass.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam PipelineT Pipeline type; must have no world cull stage
 *         (pipeline_hoistable_cull), so the predicate sees the raw points.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param trail Geodesic fragment polyline (>= 2 unit-position points).
 * @param bits Output, one byte per edge (trail.size() - 1): 0 = culled, else
 *        1; valid as rasterize()'s edge_visible input.
 * @return False when no edge is visible; bits are then all zero.
 * @details The hoisted per-point coordinates and the whole-trail culls come
 * from trail_gate_prologue, shared with ParticleSystem::draw's gate.
 */
HS_O3_BEGIN
template <int W, int H, typename PipelineT>
static bool gate_trail_edges(const PipelineT &, const ClipRegion &cr,
                             const ClipRegion::XClip &xc,
                             const Fragments &trail, uint8_t *bits) {
  static_assert(pipeline_hoistable_cull<PipelineT>(),
                "gate_trail_edges requires a pipeline with no world cull "
                "stage; route others through edge_visible_in_clip");
  constexpr int H_VIRT = H + hs::H_OFFSET;
  const size_t n = trail.size();
  HS_CHECK(n >= 2);
  const size_t edges = n - 1;

  ScratchScope span_guard(scratch_arena_a);
  const TrailGatePrologue pro = trail_gate_prologue<W, H>(cr, xc, trail);
  if (pro.rejected) {
    std::fill_n(bits, edges, uint8_t{0});
    return false;
  }
  const float *rows = pro.rows;
  const float *cols = pro.cols;

  bool any = false;
  for (size_t e = 0; e < edges; ++e) {
    const Vector &ea = trail[e].pos;
    const Vector &eb = trail[e + 1].pos;

    // Cheap row tier: the exact span's interior extremum lies within arc/2 of
    // an endpoint and phi is 1-Lipschitz in arc length (arc <= (pi/2)*chord),
    // so the endpoint rows widened by chord*(H_VIRT-1)/4 contain the exact
    // span — a miss here implies the exact test below also misses, keeping
    // the bits identical while skipping the edge's cross/normalize/acos.
    {
      const Vector d = eb - ea;
      const float margin =
          sqrtf(dot(d, d)) * (static_cast<float>(H_VIRT - 1) * 0.25f);
      if (!cr.could_intersect_y(std::min(rows[e], rows[e + 1]) - margin,
                                std::max(rows[e], rows[e + 1]) + margin)) {
        bits[e] = 0;
        continue;
      }
    }

    const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
    const bool v = exact_geodesic_edge_visible_hoisted<W, H>(cr, xc, rows, cols,
                                                             e, ea, eb, es);
    bits[e] = v ? 1 : 0;
    any = any || v;
  }
  return any;
}
HS_O3_END

/**
 * @brief Mesh drawing.
 * Registers:
 *  v0: Edge Progress t (0.0 -> 1.0) per edge
 *  v1: Cumulative Arc Length (radians) per edge
 *  v2: Edge index
 */
struct Mesh {
  /**
   * @brief Max distinct vertices the edge-dedup bitset can track.
   * @details A mesh exceeding this traps while its faces are walked: at render
   * time for draw() (per frame for MeshFeedback), or at setup for
   * extract_edges(). Sized for a TriangularBitset of 128*127/2 bits = 1016
   * bytes.
   */
  static constexpr int DEDUP_CAPACITY = 128;

  /**
   * @brief Fragments one wireframe edge can produce: its two endpoints plus a
   *        cut at every clip-band boundary it crosses.
   */
  static constexpr int EDGE_MAX_POINTS = GEODESIC_CLIP_MAX_SPLITS + 2;

  /**
   * @brief Sample, shade, and rasterize one wireframe edge.
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @tparam PipelineT Pipeline type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh Mesh supplying vertex positions.
   * @param u,v Endpoint vertex indices (assumed in bounds — the callers run the
   *            cold OOB/capacity traps before delegating here).
   * @param edge_index Value written to each fragment's v2 register.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader; displaces the edge's two
   *        endpoints, which then bound one geodesic — it does not deform the
   *        edge's interior.
   * @details Shared body for both draw() overloads (face-walk and precomputed
   * edge list); keeping it in one place is why the two paths stay bit-identical.
   * The adaptive rasterizer already walks the true great-circle arc, so the edge
   * is handed over whole; under a clip it is cut at the band's boundaries
   * (geodesic_clip_splits) so the outside pieces cost nothing past their gate
   * verdict.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void draw_edge(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
                        int u, int v, int edge_index,
                        FragmentShaderFn fragment_shader,
                        VertexShaderRef vertex_shader) {
    Fragment fu;
    fu.pos = mesh.vertices[u];
    Fragment fv;
    fv.pos = mesh.vertices[v];

    const ClipRegion &cr = canvas.clip();
    const bool clip_active = !cr.is_full();
    const ClipRegion::XClip xc = cr.x_clip();

    // A vertex shader moves the endpoints after this test, so it opts out.
    if constexpr (pipeline_hoistable_cull<PipelineT>()) {
      if (clip_active && !vertex_shader &&
          !edge_visible_in_clip<W, H>(pipeline, cr, xc, fu.pos, fv.pos,
                                      nullptr))
        return;
    }

    ScratchScope edge_guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, EDGE_MAX_POINTS);

    const GeodesicEdgeSpan es = make_geodesic_edge_span(fu.pos, fv.pos);
    bool split = false;
    if constexpr (pipeline_hoistable_cull<PipelineT>())
      split = clip_active && !vertex_shader && es.have_axis;

    if (split) {
      float ts[GEODESIC_CLIP_MAX_SPLITS];
      const int cuts = geodesic_clip_splits<W, H>(fu.pos, es, cr, xc, ts);
      const Vector perp = cross(es.axis, fu.pos);
      points.push_back(Line::sample_point(fu, fv, es, perp, 0.0f));
      for (int i = 0; i < cuts; ++i)
        points.push_back(Line::sample_point(fu, fv, es, perp, ts[i]));
      points.push_back(Line::sample_point(fu, fv, es, perp, 1.0f));
    } else {
      Line::sample(points, fu, fv);
    }

    for (auto &p : points)
      p.v2 = static_cast<float>(edge_index); // Edge Index
    if (vertex_shader) {
      vertex_shader(points[0]);
      vertex_shader(points.back());
    }

    if constexpr (pipeline_hoistable_cull<PipelineT>()) {
      if (clip_active) {
        uint8_t bits[EDGE_MAX_POINTS - 1];
        HS_CHECK(points.size() >= 2 && points.size() <= EDGE_MAX_POINTS);
        if (points.size() == 2 && !vertex_shader) {
          // Uncut, unshaded: the whole-edge test above already ran on it.
          bits[0] = EDGE_VISIBLE;
        } else if (!gate_trail_edges<W, H>(pipeline, cr, xc, points, bits)) {
          return;
        }
        rasterize<W, H>(pipeline, canvas, points, fragment_shader,
                        {.edge_visible = bits});
        return;
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader);
  }

  /**
   * @brief Walk a mesh's faces and invoke fn(u, v) once per unique edge.
   * @tparam MeshT Mesh type.
   * @tparam Fn Per-edge callback type.
   * @param mesh Mesh whose faces are walked for edges.
   * @param visited Caller-owned dedup bitset; cleared before walking. Held by
   *                the caller so each path picks its own arena/scope.
   * @param fn Invoked as fn(u, v) for the first occurrence of each edge.
   * @details Shared face-walk/edge-dedup loop behind both draw() and
   * extract_edges().
   */
  template <typename MeshT, typename Fn>
  static void for_each_unique_edge(const MeshT &mesh,
                                   TriangularBitset<DEDUP_CAPACITY> &visited,
                                   Fn &&fn) {
    visited.clear();

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    size_t fi_size = mesh.get_faces_size();
    size_t offset = 0;

    for (size_t i = 0; i < num_f; ++i) {
      int count = fc[i];

      // Trap malformed mesh data: an offset/count pair disagreeing with the flat
      // index array yields out-of-bounds reads. Cold per-face check.
      HS_CHECK(offset + static_cast<size_t>(count) <= fi_size,
               "mesh face span exceeds face index array");

      for (int k = 0; k < count; ++k) {
        int u = fi[offset + k];
        int v = fi[offset + (k + 1) % count];
        int small = std::min(u, v);
        int large = std::max(u, v);

        // A vertex index past the dedup bitset's capacity is a mesh-sizing bug;
        // trap at the face-walk boundary rather than drop the edge.
        HS_CHECK(large < DEDUP_CAPACITY);

        if (!visited.test_and_set(small, large))
          fn(u, v);
      }
      offset += count;
    }
  }

  /**
   * @brief Draws a mesh (wireframe).
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh The mesh to draw.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    int edge_index = 0;

    // O(1) edge dedup in a 1016-byte triangular bit matrix, arena-allocated (deep
    // render chain, tight DTCM stack). Held in scratch_arena_b so the per-edge
    // scratch_arena_a scopes below keep their headroom.
    ScratchScope visited_guard(scratch_arena_b);
    auto &visited = *new (
        scratch_arena_b.allocate(sizeof(TriangularBitset<DEDUP_CAPACITY>),
                                 alignof(TriangularBitset<DEDUP_CAPACITY>)))
                        TriangularBitset<DEDUP_CAPACITY>();

    for_each_unique_edge(mesh, visited, [&](int u, int v) {
      // mesh.vertices[] only asserts in bounds (stripped on device), so guard the
      // per-edge setup boundary here. u,v come from uint16_t face data (non-
      // negative), so max(u,v) in bounds implies both endpoints are valid.
      HS_CHECK(static_cast<size_t>(std::max(u, v)) < mesh.vertices.size());

      draw_edge<W, H>(pipeline, canvas, mesh, u, v, edge_index, fragment_shader,
                      vertex_shader);

      edge_index++;
    });
  }

  /**
   * @brief Draws a mesh (wireframe) without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @tparam PipelineT Pipeline type (defaults to PipelineRef).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh The mesh to draw.
   * @param fragment_shader Shader function.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, mesh, fragment_shader, {});
  }

  /**
   * @brief Precomputed edge pair for static-topology meshes.
   */
  struct Edge {
    uint16_t u, v; /**< Endpoint vertex indices into the mesh's vertex array. */
  };

  /**
   * @brief Extract unique edges from a mesh (call once at setup time).
   * @tparam MeshT Mesh type.
   * @param mesh Mesh whose faces are walked for unique edges.
   * @param edges Output edge list; deduplicated edges are appended.
   */
  template <typename MeshT>
  static void extract_edges(const MeshT &mesh, ArenaVector<Edge> &edges) {
    // Dedup bitset (1016 B) in the arena, not the stack (deep setup chain). The
    // output `edges` lives in a separate persistent arena, so scratch_arena_b
    // cannot disturb it.
    ScratchScope visited_guard(scratch_arena_b);
    auto &visited = *new (
        scratch_arena_b.allocate(sizeof(TriangularBitset<DEDUP_CAPACITY>),
                                 alignof(TriangularBitset<DEDUP_CAPACITY>)))
                        TriangularBitset<DEDUP_CAPACITY>();

    for_each_unique_edge(mesh, visited, [&](int u, int v) {
      edges.push_back({(uint16_t)u, (uint16_t)v});
    });
  }

  /**
   * @brief Draw using a precomputed edge list (skips face walk + dedup).
   * @tparam W,H Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @tparam PipelineT Pipeline type (defaults to PipelineRef).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh Mesh supplying vertex positions.
   * @param edges Precomputed unique edge list.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param mask Optional dissolve ownership mask (see DissolveMask), keyed on
   *        the edge's endpoint indices; an unowned edge is skipped before any
   *        of its geometry work, so two complementary masks split one
   *        wireframe's cost across the two meshes of a transition.
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void
  draw(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
       const ArenaVector<Edge> &edges, FragmentShaderFn fragment_shader,
       VertexShaderRef vertex_shader = {}, const DissolveMask *mask = nullptr) {
    for (size_t ei = 0; ei < edges.size(); ++ei) {
      if (mask && !mask->owns(edges[ei].u, edges[ei].v))
        continue;
      // Setup-boundary OOB guard (see the face-walk overload above): the raw
      // edge list could outlive or mismatch its mesh, and mesh.vertices[] only
      // asserts (compiled out on device).
      HS_CHECK(edges[ei].u < mesh.vertices.size() &&
               edges[ei].v < mesh.vertices.size());

      draw_edge<W, H>(pipeline, canvas, mesh, edges[ei].u, edges[ei].v,
                      static_cast<int>(ei), fragment_shader, vertex_shader);
    }
  }
};

/**
 * @brief Particle System trails.
 * Registers:
 *  v0: Trail Progress (0.0=Tail -> 1.0=Head)
 *  v1: Reserved (always 0)
 *  v2: Particle ID, or the particle_v2 mapper's value
 *  v3: Normalized TTL
 */
struct ParticleSystem {
  /**
   * @brief Draws each active particle's history as a rasterized trail.
   * @tparam W,H Rasterization resolution.
   * @tparam HoistableCull True when raw point projections are valid for clip
   *         gating because the source pipeline has no world cull stage.
   * @tparam FuseVertex Apply the typed vertex shader as each point is emitted.
   * @tparam FragmentShaderT Fragment shader type.
   * @tparam VertexShaderFn Vertex shader type.
   * @tparam DeferredShaderT Deferred vertex shader type.
   * @tparam SystemT Particle-system type.
   * @tparam ParticleV2Fn Per-particle v2 mapper type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param system Particle system supplying the active pool and trail history.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader (position pass).
   * @param deferred_shader Optional second vertex pass, given each fragment and
   *        its original pre-shader position. Under an active clip it runs only
   *        for trails with at least one cull-surviving edge; a skipped trail
   *        renders nothing, so output is identical to an undeferred shader.
   *        Put per-point work that only affects shading registers here.
   * @param particle_v2 Optional mapper from particle and pool index to v2.
   *        Called once per materialized particle; the default stores the pool
   *        index.
   */
  template <int W, int H, bool HoistableCull, bool FuseVertex,
            bool SinglePassRaster, typename PipelineT, typename SystemT,
            typename FragmentShaderT, typename VertexShaderFn,
            typename DeferredShaderT, typename ParticleV2Fn>
  static void
  draw_impl(PipelineT &pipeline, Canvas &canvas, const SystemT &system,
            FragmentShaderT fragment_shader, VertexShaderFn vertex_shader,
            DeferredShaderT deferred_shader, ParticleV2Fn particle_v2) {
    int count = system.active();
    if (count == 0)
      return;

    const float max_life = static_cast<float>(system.max_life);
    HS_CHECK(std::isfinite(max_life) && max_life >= 1.0f &&
                 max_life <= 65535.0f,
             "ParticleSystem render max_life must be finite and in [1, 65535]");
    const float inv_max_life = 1.0f / max_life;
    const bool has_deferred_shader = [&] {
      if constexpr (std::is_same_v<std::decay_t<DeferredShaderT>,
                                   std::nullptr_t>)
        return false;
      else if constexpr (requires { static_cast<bool>(deferred_shader); })
        return static_cast<bool>(deferred_shader);
      else
        return true;
    }();

    // Segment-clip state for the trail-level deferred-shader gate below.
    const auto &cr = canvas.clip();
    const bool clip_active = !cr.is_full();
    const auto xc = cr.x_clip();
    CartesianQuadrantClip cartesian_clip;
    if constexpr (HoistableCull)
      cartesian_clip = make_cartesian_quadrant_clip<W, H>(cr);

    for (int i = 0; i < count; ++i) {
      const auto &p = system.pool[i];
      const size_t trail_len = p.history.length();
      HS_MSP_COUNT(resident_particles);
      if (p.life == 0)
        HS_MSP_COUNT(draining_histories);
      else {
        HS_MSP_COUNT(live_particles);
        if (trail_len == static_cast<size_t>(p.history.CAPACITY))
          HS_MSP_COUNT(full_histories);
        else
          HS_MSP_COUNT(partial_histories);
      }
      if (p.life == 0 || trail_len < 2)
        continue;
      const float v2 = [&] {
        if constexpr (std::is_same_v<ParticleV2Fn, std::nullptr_t>)
          return static_cast<float>(i);
        else
          return particle_v2(p, i);
      }();
      const float particle_life = static_cast<float>(p.life) * inv_max_life;
      ScratchScope trail_guard(scratch_arena_a);
      Fragments trail;
      trail.bind(scratch_arena_a, trail_len);
      // Original (pre-shader) positions, kept for the deferred pass.
      ArenaVector<Vector> orig;
      if (has_deferred_shader)
        orig.bind(scratch_arena_a, trail_len);
      {
        HS_PROFILE(plot_ps_tween);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
        hs::DwtStallBatch history_batch(
            hs::g_mindsplatter_stalls.history_vertex);
#endif
        p.history.tween([&](const Vector &v, float t) {
          trail.emplace_back(
              Fragment{.pos = v, .v0 = t, .v2 = v2, .v3 = particle_life});
          if constexpr (FuseVertex)
            vertex_shader(trail.back());
          if (has_deferred_shader)
            orig.push_back(v);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
          history_batch.step();
#endif
        });
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
        history_batch.finish();
#endif
      }

      if (trail.is_empty())
        continue;
      if constexpr (!FuseVertex) {
        HS_PROFILE(plot_ps_vertex);
        apply_vertex_shader(vertex_shader, trail);
      }

      // Trail-level gate: precompute each edge's cull verdict from the
      // position-shaded points. No visible edge means the trail renders
      // nothing, so the optional deferred pass and the rasterize call are
      // skipped whole; the bits feed rasterize so the cull is evaluated once.
      const uint8_t *vis = nullptr;
      const float *dot_rows = nullptr;
      const float *dot_cols = nullptr;
      if (clip_active && trail.size() >= 2) {
        HS_PROFILE(plot_ps_gate);
        const size_t edges = trail.size() - 1;
        auto *bits = static_cast<uint8_t *>(
            scratch_arena_a.allocate(edges, alignof(uint8_t)));
        bool any = false;
        if constexpr (HoistableCull) {
          CartesianTrailGateResult cartesian_result;
          {
            HS_PROFILE(plot_ps_cartesian_gate);
            cartesian_result =
                cartesian_quadrant_trail_gate(cartesian_clip, trail);
          }
          count_cartesian_trail_gate_result(cartesian_result);
          if (cartesian_result != CartesianTrailGateResult::EXACT_FALLBACK)
            continue;

          // No stage re-emits edges, so the predicate sees the raw points:
          // per-point rows/columns are computed once and shared by every edge,
          // and a conservative whole-trail bound rejects fully-invisible
          // trails before any per-edge work.
          const TrailGatePrologue pro =
              trail_gate_prologue<W, H>(cr, xc, trail);
          if (pro.rejected)
            continue;
          const float *rows = pro.rows;
          const float *cols = pro.cols;
          // rasterize's single-dot shortcut takes both projections or neither.
          if (cols != nullptr) {
            dot_rows = rows;
            dot_cols = cols;
          }

          for (size_t e = 0; e < edges; ++e) {
            HS_MSP_STALL_START(edge_gate_start);
            const Vector &ea = trail[e].pos;
            const Vector &eb = trail[e + 1].pos;
            const bool one_dot = edge_fits_one_dot<W, H>(ea, eb);
            count_particle_edge_class(one_dot);
            if (one_dot) {
              bool v = antialiased_dot_visible_in_clip<W, H>(
                  cr, xc, rows[e], cols != nullptr ? cols[e] : 0.0f);
              if (e + 1 == edges)
                v = v || antialiased_dot_visible_in_clip<W, H>(
                             cr, xc, rows[e + 1],
                             cols != nullptr ? cols[e + 1] : 0.0f);
              bits[e] = EDGE_CLASSIFIED | EDGE_ONE_DOT |
                        (v ? EDGE_VISIBLE : uint8_t{0});
              if (!v)
                HS_MSP_COUNT(edge_rejects);
              any = any || v;
              HS_MSP_STALL_STOP(trail_gate, edge_gate_start);
              continue;
            }
            bool v;
            const RawGeodesicGateResult raw = raw_geodesic_edge_gate<W, H>(
                cr, xc, rows[e], rows[e + 1], cols != nullptr ? cols[e] : 0.0f,
                cols != nullptr ? cols[e + 1] : 0.0f, ea, eb);
            if (raw != RawGeodesicGateResult::EXACT_FALLBACK) {
              v = raw == RawGeodesicGateResult::VISIBLE;
            } else {
              count_particle_exact_gate_fallback();
              const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
              v = exact_geodesic_edge_visible_hoisted<W, H>(cr, xc, rows, cols,
                                                            e, ea, eb, es);
            }
            bits[e] = EDGE_CLASSIFIED | (v ? EDGE_VISIBLE : uint8_t{0});
            if (!v)
              HS_MSP_COUNT(edge_rejects);
            any = any || v;
            HS_MSP_STALL_STOP(trail_gate, edge_gate_start);
          }
        } else {
          for (size_t e = 0; e < edges; ++e) {
            bits[e] = edge_visible_in_clip<W, H>(pipeline, cr, xc, trail[e].pos,
                                                 trail[e + 1].pos, nullptr)
                          ? 1
                          : 0;
            any = any || bits[e] != 0;
          }
        }
        if (!any)
          continue;
        HS_MSP_COUNT(visible_trails);
        vis = bits;
      }
      if (!clip_active)
        HS_MSP_COUNT(visible_trails);

      if (has_deferred_shader) {
        HS_PROFILE(plot_ps_deferred);
        for (size_t k = 0; k < trail.size(); ++k)
          deferred_shader(trail[k], orig[k]);
      }
      {
        HS_PROFILE(plot_ps_raster);
        rasterize<W, H, SinglePassRaster, true>(pipeline, canvas, trail,
                                                fragment_shader,
                                                {.edge_visible = vis,
                                                 .point_rows = dot_rows,
                                                 .point_cols = dot_cols});
      }
    }
  }

  /**
   * @brief Draws particle trails through the shared raster surface.
   * @details The source pipeline's static cull trait is retained separately so
   *          Cartesian, one-dot and raw-edge gates remain available after plot
   *          dispatch is erased. A pipeline explicitly declaring the direct
   *          raster path retains its compile-time plot calls.
   */
  template <int W, int H, typename PipelineT = PipelineRef,
            typename ParticleV2Fn = std::nullptr_t>
  static void draw(PipelineT &pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader,
                   DeferredShaderRef deferred_shader = {},
                   ParticleV2Fn particle_v2 = nullptr) {
    if constexpr (pipeline_direct_raster_path<PipelineT>()) {
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), false, false>(
          pipeline, canvas, system, fragment_shader, vertex_shader,
          deferred_shader, particle_v2);
    } else {
      PipelineRef erased(pipeline);
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), false, false>(
          erased, canvas, system, fragment_shader, vertex_shader,
          deferred_shader, particle_v2);
    }
  }

  /**
   * @brief Draws particle trails while applying a typed vertex shader during
   *        trail materialization.
   * @details Fuses point materialization and transformation into one traversal.
   * @tparam SinglePassRaster Emit adaptive geodesic samples without a replay
   *         cache. Applies only to a direct raster pipeline.
   * @tparam DeferredShaderT Deferred vertex shader type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param system Particle system supplying the active pool and trail history.
   * @param fragment_shader Shader function.
   * @param vertex_shader Typed vertex shader.
   * @param deferred_shader Optional deferred vertex shader.
   * @param particle_v2 Optional particle-to-v2 mapper.
   */
  template <int W, int H, bool SinglePassRaster = false, typename PipelineT,
            typename FragmentShaderT, typename VertexShaderFn,
            typename DeferredShaderT = DeferredShaderRef,
            typename ParticleV2Fn = std::nullptr_t>
  static void draw_fused_vertex(PipelineT &pipeline, Canvas &canvas,
                                const auto &system,
                                FragmentShaderT fragment_shader,
                                VertexShaderFn vertex_shader,
                                DeferredShaderT deferred_shader = {},
                                ParticleV2Fn particle_v2 = nullptr) {
    if constexpr (pipeline_direct_raster_path<PipelineT>()) {
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), true,
                SinglePassRaster>(pipeline, canvas, system, fragment_shader,
                                  vertex_shader, deferred_shader, particle_v2);
    } else {
      PipelineRef erased(pipeline);
      FragmentShaderFn erased_shader(fragment_shader);
      draw_impl<W, H, pipeline_hoistable_cull<PipelineT>(), true, false>(
          erased, canvas, system, erased_shader, vertex_shader, deferred_shader,
          particle_v2);
    }
  }

  /**
   * @brief Draws particle trails without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @tparam PipelineT Pipeline type (defaults to PipelineRef).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param system Particle system supplying the active pool and trail history.
   * @param fragment_shader Shader function.
   */
  template <int W, int H, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, system, fragment_shader, {});
  }
};

} // namespace Plot
