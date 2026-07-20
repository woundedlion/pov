/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <utility>
#include <type_traits>
#include <algorithm>
#include <cmath>
#include <array>
#include <concepts>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif
#include "math/geometry.h"
#include "render/shading.h"
#include "color/color.h"
#include "engine/constants.h"
#include "render/canvas.h"
#include "animation/animation.h"

namespace Plot {

/**
 * @brief Inner/outer radius ratio for star shapes.
 */
static constexpr float STAR_INNER_RATIO = ::STAR_INNER_RATIO;

/**
 * @brief Geodesic segment shorter than this (radians) collapses to a point.
 * @details Deliberately ~10× math::EPS_GEOMETRIC — a looser threshold for
 * slerp-axis stability, not positional near-equality.
 */
static constexpr float EPS_GEODESIC_SEGMENT = 0.001f;

/**
 * @brief Minimum |axis.y| for which a geodesic edge's endpoint columns bound
 *        its azimuth span.
 * @details Below this the great circle runs near the poles, where longitude is
 * ill-conditioned and the interior leaves the endpoint columns.
 */
static constexpr float AXIS_Y_EPS = 1e-4f;

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
 * @brief Apply an optional per-control-point vertex shader to every fragment.
 * @tparam Fragments Fragment container type.
 * @param vertex_shader Vertex shader to run on each fragment; no-op if null.
 * @param pts Fragment container mutated in place.
 * @details Shared inline replacement for the per-primitive
 * `if (vertex_shader) for (auto &p : pts) vertex_shader(p);` block.
 */
template <typename Fragments>
inline void apply_vertex_shader(VertexShaderRef vertex_shader, Fragments &pts) {
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
  float R = sqrtf(Px * Px + Py * Py);
  if (R < math::EPS_GEOMETRIC)
    return basis.v;
  float theta = fast_atan2(Py, Px);
  Vector axis = (basis.u * fast_cosf(theta)) + (basis.w * fast_sinf(theta));
  return (basis.v * fast_cosf(R)) + (axis * fast_sinf(R));
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
 * @details arc_cumul[0] = 0; arc_cumul.back() is the full rendered length. The
 * planar edge bows away from the great-circle chord, so the total exceeds the
 * chord's angle_between. Shared by the planar rasterizer (which inverts the table
 * for arc-uniform stepping) and rasterize()'s perimeter pre-pass (which takes the
 * total), so both sample identical points.
 */
static inline void
planar_arc_cumul(const std::pair<float, float> &proj, float dx, float dy,
                 const Basis &planar_basis,
                 std::array<float, PLANAR_LEN_SAMPLES + 1> &arc_cumul) {
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
 * @brief Planar interpolation strategy: builds an arc-uniform sampler for one edge.
 * @tparam ProcessSegmentFn Callable (sample, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param planar_basis Azimuthal-equidistant projection basis.
 * @param isLastSegment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length sampler (position + a
 *                        finite-difference unit tangent), endpoints, on-sphere
 *                        length (radians), and the last-segment flag.
 * @details The path is a straight line in the azimuthal-equidistant projection.
 * Projection-uniform stepping is NOT arc-uniform under the anisotropic metric,
 * so a short cumulative-arc table is inverted to turn an arc-length fraction into
 * a projection parameter, making planar sampling arc-uniform with no new trig.
 */
template <typename ProcessSegmentFn>
static void
rasterize_planar_strategy(const Fragment &curr, const Fragment &next,
                          const Basis &planar_basis, bool isLastSegment,
                          ProcessSegmentFn &&process_segment) {
  auto proj1 = azimuthal_project(curr.pos, planar_basis);
  auto proj2 = azimuthal_project(next.pos, planar_basis);
  float dx = proj2.first - proj1.first;
  float dy = proj2.second - proj1.second;

  // The path is a straight line in the azimuthal-equidistant projection,
  // parameterized here by the PROJECTION fraction p in [0,1] (not arc length).
  auto unproject = [=, &planar_basis](float p) -> Vector {
    return azimuthal_unproject(proj1.first + dx * p, proj1.second + dy * p,
                               planar_basis);
  };

  // Cumulative on-sphere arc length at evenly-spaced PROJECTION samples: the
  // path's true length, which drives the step count and lets map() invert an
  // arc-length fraction back to a projection parameter (projection-uniform
  // stepping is not arc-uniform under the anisotropic metric).
  std::array<float, PLANAR_LEN_SAMPLES + 1> arc_cumul;
  planar_arc_cumul(proj1, dx, dy, planar_basis, arc_cumul);
  const float dist = arc_cumul[PLANAR_LEN_SAMPLES];

  // Arc-length parameterization: s is the arc fraction in [0,1]. Invert the
  // piecewise-linear cumulative-arc table to a projection parameter, then
  // unproject. A short scan over PLANAR_LEN_SAMPLES floats — no trig.
  auto map_planar = [=](float s) -> Vector {
    if (dist < math::EPS_GEOMETRIC)
      return unproject(s);
    float target = s * dist;
    int k = 0;
    while (k < PLANAR_LEN_SAMPLES - 1 && arc_cumul[k + 1] < target)
      ++k;
    float seg = arc_cumul[k + 1] - arc_cumul[k];
    float frac =
        (seg > math::EPS_GEOMETRIC) ? (target - arc_cumul[k]) / seg : 0.0f;
    float p = (static_cast<float>(k) + frac) / PLANAR_LEN_SAMPLES;
    return unproject(std::min(1.0f, std::max(0.0f, p)));
  };

  // The azimuthal map has no closed-form tangent, so take it from a short forward
  // difference (backward at the s=1 end); the difference direction is the unit
  // tangent regardless of magnitude. Feeds the same screen-velocity sub-step
  // sampler as the geodesic path.
  auto sample_planar = [=](float s) -> SamplePT {
    Vector p = map_planar(s);
    bool fwd = (s + PLANAR_TAN_DT <= 1.0f);
    Vector q = map_planar(fwd ? s + PLANAR_TAN_DT : s - PLANAR_TAN_DT);
    Vector d = fwd ? (q - p) : (p - q);
    return {p, normalized_or(d, Vector())};
  };

  process_segment(sample_planar, curr, next, dist, isLastSegment);
}

/**
 * @brief On-sphere arc length (radians) of the azimuthal-equidistant straight
 *        edge a->b — the length actually rendered under planar interpolation.
 * @details Shares planar_arc_cumul with rasterize_planar_strategy, so
 * rasterize()'s perimeter pre-pass and per-segment arc accumulator sum exactly
 * the lengths the draw phase walks. The planar edge bows away from the
 * great-circle chord, so this exceeds angle_between(a, b).
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
 * rotation axis. Shared by rasterize_geodesic_strategy and Line::sample.
 */
static inline Vector stable_perpendicular_axis(const Vector &v) {
  return cross(v, least_parallel_axis(v)).normalized();
}

/**
 * @brief Azimuthal-equidistant projection chart centered on a pole.
 * @param center Unit pole the planar chart is centered on (the 'v' axis).
 * @return A Basis {u, center, w} with u, w spanning the chart plane.
 */
static inline Basis planar_chart_basis(const Vector &center) {
  Vector ref = least_parallel_axis(center);
  Vector u = cross(center, ref).normalized();
  Vector w = cross(center, u).normalized();
  return {u, center, w};
}

/**
 * @brief Geodesic interpolation strategy: builds a great-circle sampler for one edge.
 * @tparam ProcessSegmentFn Callable (sample, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param isLastSegment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length sampler (position + unit
 *                        tangent), endpoints, on-sphere length (radians), and
 *                        the last-segment flag.
 * @details Picks a stable perpendicular axis for antipodal/degenerate endpoints
 * and slerps along the great circle; a coincident-endpoint edge collapses to a
 * constant sampler.
 */
HS_O3_BEGIN
template <typename ProcessSegmentFn>
static void rasterize_geodesic_strategy(const Fragment &curr,
                                        const Fragment &next,
                                        bool isLastSegment,
                                        ProcessSegmentFn &&process_segment) {
  Vector v1 = curr.pos;
  Vector v2 = next.pos;
  float total_dist = angle_between(v1, v2);

  if (total_dist < EPS_GEODESIC_SEGMENT) {
    auto sample_degenerate = [=](float) -> SamplePT { return {v1, Vector()}; };
    process_segment(sample_degenerate, curr, next, total_dist, isLastSegment);
  } else {
    Vector axis;
    if (std::abs(PI_F - total_dist) < TOLERANCE) {
      axis = stable_perpendicular_axis(v1);
    } else {
      axis = cross(v1, v2).normalized();
    }

    Vector v_perp = cross(axis, v1);

    auto sample_geodesic = [=](float t) -> SamplePT {
      float ang = total_dist * t;
      float s = fast_sinf(ang);
      float c = fast_cosf(ang);
      // pos on the great circle; tan = d(pos)/d(ang) is a unit vector (v1,
      // v_perp orthonormal), so the screen-velocity sampler's tangent comes free
      // from the same sin/cos — no extra trig.
      return {(v1 * c) + (v_perp * s), (v_perp * c) - (v1 * s)};
    };
    process_segment(sample_geodesic, curr, next, total_dist, isLastSegment);
  }
}
HS_O3_END

/**
 * @brief Shared per-edge geodesic setup for the row/column span bounds.
 * @details axis mirrors rasterize_geodesic_strategy's slerp-axis selection;
 *          have_axis is false when cross(a, b) collapses on a non-antipodal
 *          edge (no stable arc pole exists).
 */
struct GeodesicEdgeSpan {
  float total;    /**< angle_between(a, b) in radians. */
  bool antipodal; /**< |π - total| < TOLERANCE. */
  bool have_axis; /**< axis holds a unit arc pole. */
  Vector axis;    /**< Unit arc pole (valid iff have_axis). */
};

/**
 * @brief Computes the shared geodesic edge setup once per edge.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 */
static inline GeodesicEdgeSpan make_geodesic_edge_span(const Vector &a,
                                                       const Vector &b) {
  GeodesicEdgeSpan es;
  es.total = angle_between(a, b);
  es.antipodal = std::abs(PI_F - es.total) < TOLERANCE;
  if (es.antipodal) {
    es.axis = stable_perpendicular_axis(a);
    es.have_axis = true;
    return es;
  }
  Vector c = cross(a, b);
  float L2 = dot(c, c);
  if (L2 <= math::EPS_GEOMETRIC * math::EPS_GEOMETRIC) {
    es.axis = Vector(0.0f, 0.0f, 0.0f);
    es.have_axis = false;
    return es;
  }
  es.axis = c * (1.0f / sqrtf(L2));
  es.have_axis = true;
  return es;
}

constexpr int PLANAR_SPAN_SAMPLES = 8;

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
 * @brief Screen row of a unit-sphere y coordinate (the renderer's row map).
 * @tparam H Rasterization height (pixel grid).
 * @param y Unit-sphere y in [-1, 1] (clamped).
 */
template <int H> static inline float y_to_screen_row(float y) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  return phi_to_y(fast_acos(hs::clamp(y, -1.0f, 1.0f)), H_VIRT);
}

/**
 * @brief Geodesic screen-row span from precomputed endpoint rows.
 * @tparam W,H Rasterization resolution (pixel grid).
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
template <int W, int H>
static __attribute__((always_inline)) inline void geodesic_row_span_rows(float ra, float rb, const Vector &a,
                                          const Vector &b,
                                          const GeodesicEdgeSpan &es,
                                          float &row_lo, float &row_hi) {
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
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 */
template <int W, int H>
static inline void geodesic_row_span(const Vector &a, const Vector &b,
                                     const GeodesicEdgeSpan &es, float &row_lo,
                                     float &row_hi) {
  geodesic_row_span_rows<W, H>(y_to_screen_row<H>(a.y), y_to_screen_row<H>(b.y),
                               a, b, es, row_lo, row_hi);
}

/**
 * @brief Planar screen-row span from a precomputed edge setup.
 * @tparam W,H Rasterization resolution (pixel grid).
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
template <int W, int H>
static inline void planar_row_span(const Vector &a, const Vector &b,
                                   const PlanarEdgeSpan &es, float &row_lo,
                                   float &row_hi) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  auto y_to_row = [](float y) {
    return phi_to_y(fast_acos(hs::clamp(y, -1.0f, 1.0f)), H_VIRT);
  };
  float ra = y_to_row(a.y);
  float rb = y_to_row(b.y);
  row_lo = std::min(ra, rb);
  row_hi = std::max(ra, rb);
  for (const Vector &s : es.interior) {
    float r = y_to_row(s.y);
    row_lo = std::min(row_lo, r);
    row_hi = std::max(row_hi, r);
  }
  float margin =
      es.gap_arc * (static_cast<float>(H_VIRT - 1) / PI_F) + 1.0f;
  row_lo -= margin;
  row_hi += margin;
}

/**
 * @brief Conservative screen-row span of a rendered edge, arc bulge included.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param planar_basis Non-null: edge is azimuthal-equidistant; null: geodesic.
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 * @details The clip cull must not test endpoints alone: an edge between two
 * points outside a clip band can still bulge through it. Rows come from
 * row = phi_to_y(acos(y)); the endpoint rows are extended by the arc's interior
 * latitude extremum (closed-form for geodesic, sampled with a Lipschitz margin
 * for planar). Runs once per coarse edge on the clip-only path.
 */
template <int W, int H>
static inline void edge_row_span(const Vector &a, const Vector &b,
                                 const Basis *planar_basis, float &row_lo,
                                 float &row_hi) {
  if (planar_basis == nullptr) {
    geodesic_row_span<W, H>(a, b, make_geodesic_edge_span(a, b), row_lo,
                            row_hi);
    return;
  }
  planar_row_span<W, H>(a, b, make_planar_edge_span(a, b, *planar_basis),
                        row_lo, row_hi);
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
static __attribute__((always_inline)) inline void finish_col_span(float s_f, float len_f, int &col_s,
                                   int &col_len) {
  constexpr int COL_PAD = 2;
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
 * near-antipodal boundary, where short-way is float noise. Axis selection
 * mirrors rasterize_geodesic_strategy and the column mapping is the renderer's
 * vector_to_theta; COL_PAD absorbs plot rounding and the AntiAlias tap spread.
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
 * @brief Planar screen-column arc from a precomputed edge setup.
 * @tparam W Rasterization width (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param planar_basis Azimuthal-equidistant projection basis.
 * @param es Shared setup from make_planar_edge_span(a, b, planar_basis); the
 *           end point enters through its projected chord.
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
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
                                   int &col_len) {
  constexpr float MIN_SIN_PHI = 0.05f;

  const float ca = vector_to_theta<W>(a);
  float s_f, len_f;

  {
    float min_sp2 = 1.0f - a.y * a.y; // squared sin(phi), minimized over samples
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
    step(azimuthal_unproject(es.p1.first + es.dX, es.p1.second + es.dY,
                             planar_basis));

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

/**
 * @brief Conservative screen-column arc of a rendered edge.
 * @tparam W Rasterization width (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param planar_basis Non-null: edge is azimuthal-equidistant; null: geodesic.
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 * @return False when no useful bound exists — the caller must skip the
 *         horizontal cull (see geodesic_col_span / planar_col_span).
 */
template <int W>
static inline bool edge_col_span(const Vector &a, const Vector &b,
                                 const Basis *planar_basis, int &col_s,
                                 int &col_len) {
  if (planar_basis == nullptr)
    return geodesic_col_span<W>(a, b, make_geodesic_edge_span(a, b), col_s,
                                col_len);
  return planar_col_span<W>(a, *planar_basis,
                            make_planar_edge_span(a, b, *planar_basis), col_s,
                            col_len);
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
  const float KX = W / (2.0f * PI_F);     // columns per radian of longitude
  const float KY = (H_VIRT - 1) / PI_F;   // rows per radian of colatitude
  // sin²φ = 1 - y²; floored so the pole (sin φ → 0) yields a finite, large
  // velocity (hence the min-clamped step) rather than a divide-by-zero.
  const float sin2 = std::max(1e-7f, 1.0f - pos.y * pos.y);
  const float inv_sin = 1.0f / sqrtf(sin2);
  const float dphi_ds = -tan.y * inv_sin;
  const float dlon_ds = (pos.x * tan.z - pos.z * tan.x) / sin2;
  const float vx = KX * dlon_ds;
  const float vy = KY * dphi_ds;
  const float speed = sqrtf(vx * vx + vy * vy);
  // Degenerate-speed floor: guards 1/speed when a zero/near-zero tangent stalls
  // the curve, yielding base_step rather than an unbounded step.
  const float step = SCREEN_STEP_PX / std::max(speed, 1e-6f);
  return std::max(base_step * MIN_POLE_SCALE, std::min(step, base_step));
}

/**
 * @brief True when @p P statically declares it has no world cull stage, so a
 *        caller may precompute per-point screen coordinates from raw geometry.
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
  return F2 * (KX2 * cx * cx + KY2 * ty * ty * sin2) <= SPX2 * sin2 * sin2;
}
HS_O3_END

/**
 * @brief Tier-3 clip visibility of one polyline edge, routed through the
 *        pipeline's world stages.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam PipelineT Pipeline type.
 * @param pipeline Render pipeline; world stages re-emit the edge under their
 *        plot-time rotations before the span test.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param band_len Clip band length in columns (seam-unwrapped).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param pb Planar projection basis for the edge, or null for geodesic.
 * @return True if the rendered edge could produce a pixel inside the clip.
 * @details The single definition of the segment cull: rasterize evaluates it
 * per edge, and Plot::ParticleSystem::draw precomputes it per trail to gate
 * the deferred shader — both must agree exactly or the two paths diverge.
 */
template <int W, int H, typename PipelineT>
static inline bool edge_visible_in_clip(PipelineT &pipeline,
                                        const ClipRegion &cr,
                                        const ClipRegion::XClip &xc,
                                        int band_len, const Vector &a,
                                        const Vector &b, const Basis *pb) {
  auto pred = [&](const Vector &ea, const Vector &eb, const Basis *bp) {
    float row_lo, row_hi;
    int col_s, col_len;
    if (bp == nullptr) {
      // Geodesic: both span bounds share one edge setup (angle, arc pole).
      const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
      geodesic_row_span<W, H>(ea, eb, es, row_lo, row_hi);
      if (!cr.could_intersect_y(row_lo, row_hi))
        return false;
      if (!xc.active)
        return true;
      if (!geodesic_col_span<W>(ea, eb, es, col_s, col_len))
        return true;
      return ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len, W);
    }
    // Planar: both span bounds share one projection + chart-line sample set.
    const PlanarEdgeSpan ps = make_planar_edge_span(ea, eb, *bp);
    planar_row_span<W, H>(ea, eb, ps, row_lo, row_hi);
    if (!cr.could_intersect_y(row_lo, row_hi))
      return false;
    if (!xc.active)
      return true;
    if (!planar_col_span<W>(ea, *bp, ps, col_s, col_len))
      return true;
    return ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len, W);
  };
  if constexpr (requires { pipeline.could_intersect_clip(a, b, pb, pred); }) {
    return pipeline.could_intersect_clip(a, b, pb, pred);
  } else {
    // A filter pipeline (has any_crosses_segments) must answer the clip
    // query; a signature drift would else silently fall back to raw culling.
    static_assert(!requires { PipelineT::any_crosses_segments; },
                  "pipeline exposes any_crosses_segments but not "
                  "could_intersect_clip (signature drift)");
    return pred(a, b, pb);
  }
}

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
 * @tparam PipelineT Pipeline type (defaults to PipelineRef).
 * @param pipeline Render pipeline that plots fragments.
 * @param canvas Target canvas (supplies the active clip band).
 * @param points Fragment polyline to rasterize.
 * @param fragment_shader Per-fragment shader applied before plotting; must be
 *                        non-null (the per-pixel call sites below do not guard
 *                        it, and operator()'s null assert is stripped under
 *                        NDEBUG on-device).
 * @param close_loop Also draw the last→first edge.
 * @param planar_basis Non-null selects azimuthal-equidistant interpolation
 *                     (straight in the projection); null uses geodesic edges.
 * @param omit_end Open lines only: skip the final endpoint plot (each vertex is
 *                 otherwise plotted once by its outgoing segment), so abutting
 *                 arcs tile a longer curve without double-plotting the shared
 *                 vertex.
 * @param edge_visible Optional precomputed Tier-3 visibility, one byte per
 *                     segment-loop edge. Producers must evaluate the same
 *                     edge_visible_in_clip predicate this function would.
 *                     Geodesic polylines only: a planar polyline's per-edge
 *                     basis depends on the seam pre-pass below.
 * @param point_rows Optional per-point screen rows, y_to_screen_row of each
 *                   points[k].pos. With point_cols, lets the single-dot
 *                   shortcut skip the projection. Only consumed when the
 *                   pipeline is hoistable (no world stage re-positions the
 *                   plot); both arrays or neither.
 * @param point_cols Optional per-point screen columns, vector_to_theta of
 *                   each points[k].pos.
 */
HS_O3_BEGIN
template <int W, int H, typename PipelineT = PipelineRef>
static void rasterize(PipelineT &pipeline, Canvas &canvas,
                      const Fragments &points, FragmentShaderFn fragment_shader,
                      bool close_loop = false,
                      const Basis *planar_basis = nullptr,
                      bool omit_end = false,
                      const uint8_t *edge_visible = nullptr,
                      const float *point_rows = nullptr,
                      const float *point_cols = nullptr) {
  size_t len = points.size();
  if (len < 2)
    return;
  // Trap a null shader once per polyline so the per-pixel fragment_shader()
  // calls below can't invoke a null thunk.
  HS_CHECK(fragment_shader, "rasterize requires a non-null fragment_shader");
  HS_CHECK(edge_visible == nullptr || planar_basis == nullptr,
           "precomputed edge visibility is geodesic-only");
  #ifdef __EMSCRIPTEN__
  double plot_t0 = emscripten_get_now();
  #endif

  size_t count = close_loop ? len : len - 1;
  // SCRATCH ARENA CONTRACT (load-bearing): scratch_arena_a is a LIFO bump
  // allocator shared with Pixel::Feedback::flush; do not let a raw pointer into
  // it outlive the scope that produced it.
  ScratchScope sc_guard(scratch_arena_a);
  ArenaVector<float> steps_cache;
  // The cache holds ONE segment's adaptive sub-steps (cleared per segment). Each
  // step advances ≈ one screen pixel, so a single segment needs ≲ W steps (the
  // step·≥MIN_POLE_SCALE lower clamp caps the per-segment count at the poles).
  // Size off W with 2× headroom; the simulation loop breaks at capacity as a
  // backstop.
  size_t max_cache = std::max((size_t)64, (size_t)(2 * W));
  steps_cache.bind(scratch_arena_a, max_cache);

  // PLANAR ARC REGISTERS (v0/v1): under a planar basis the rendered edge bows
  // longer than the geodesic chord, so re-derive v0/v1 from the true rendered
  // arc (`cumul`/`seg_base` track it, `total_arc` normalizes v0). Skipped for
  // geodesic polylines.
  const bool override_uv = (planar_basis != nullptr);
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
      const Vector &b = points[(i + 1) % len].pos;
      const bool seam = dot(a, pcenter) < -COS_PLANAR_ANTIPODE ||
                        dot(b, pcenter) < -COS_PLANAR_ANTIPODE;
      seg_seam_cache.push_back(seam ? 1 : 0);
      float seg = seam ? angle_between(a, b)
                       : planar_arc_length(a, b, *planar_basis);
      seg_arc_cache.push_back(seg);
      total_arc += seg;
    }
  }
  float cumul = 0.0f;    // rendered arc reached so far (planar polylines only)
  float seg_base = 0.0f; // rendered arc at the in-flight segment's start

  // Adaptively sub-step and plot one segment. `sample(t)` returns the sphere
  // point AND unit tangent at arc fraction t in [0,1] under the chosen strategy;
  // `total_dist` is the segment's on-sphere length (radians). Endpoints are
  // omitted on interior / closed segments so a shared vertex isn't plotted twice.
  auto process_segment = [&](auto &&sample, const Fragment &curr,
                             const Fragment &next, float total_dist,
                             bool isLastSegment) {
    // Rewrite the arc registers from the rendered arc when a planar basis is in
    // force (see the pre-pass above): `d` is the arc drawn so far within this
    // segment, `seg_base` the arc at its start. No-op for geodesic polylines.
    auto set_arc_uv = [&](Fragment &f, float d) {
      if (!override_uv)
        return;
      float arc = seg_base + d;
      f.v1 = arc;
      if (total_arc > math::EPS_GEOMETRIC)
        f.v0 = arc / total_arc;
    };
    // The degenerate and fast paths plot curr.pos/next.pos directly (original
    // sampled vertices, already unit), without the DRAWING PHASE renormalize that
    // corrects sample().pos's ~0.04% drift. Precondition: callers pass unit
    // fragment positions.
    // Degenerate (coincident endpoints): plot at most a single dot.
    if (total_dist < math::EPS_GEOMETRIC) {
      bool shouldOmit = close_loop || !isLastSegment || omit_end;
      if (!shouldOmit) {
        Fragment f_copy = curr;
        f_copy.color = Color4(0, 0, 0, 0);
        set_arc_uv(f_copy, 0.0f);

        fragment_shader(curr.pos, f_copy);
        pipeline.plot(canvas, curr.pos, f_copy.color.color, f_copy.age,
                      f_copy.color.alpha);
      }
      return;
    }

    // Sub-step length at the segment start (also the first simulation step).
    const float base_step = (2.0f * PI_F) / W;
    SamplePT smp = sample(0.0f);
    float first_step = screen_step<W, H>(smp.pos, smp.tan, base_step);

    // FAST PATH: the whole segment spans ≤ one screen step, so a single dot
    // covers it. Keyed on SCREEN length, not arc length: a base_step arc can
    // still cross several pixels on a steep/near-polar segment, which an
    // arc-length test would undersample into a beaded line.
    if (total_dist <= first_step) {
      Fragment f = curr;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, 0.0f);
      fragment_shader(curr.pos, f);
      pipeline.plot(canvas, curr.pos, f.color.color, f.age, f.color.alpha);
      if (!close_loop && isLastSegment && !omit_end) {
        Fragment fl = next;
        fl.color = Color4(0, 0, 0, 0);
        set_arc_uv(fl, total_dist);
        fragment_shader(next.pos, fl);
        pipeline.plot(canvas, next.pos, fl.color.color, fl.age,
                      fl.color.alpha);
      }
      return;
    }

    // SIMULATION PHASE — size each sub-step so consecutive samples land
    // ~SCREEN_STEP_PX apart in SCREEN space (screen_step). `smp`/`first_step`
    // above seed the first iteration.
    steps_cache.clear();
    float sim_dist = 0.0f;

    while (sim_dist < total_dist) {
      float step = steps_cache.is_empty()
                       ? first_step
                       : screen_step<W, H>(smp.pos, smp.tan, base_step);

      // Backstop: a pathological segment could exceed the 2*W cache. Stop
      // subdividing and let the normalized replay stretch the cached steps over
      // the rest of the segment (coarser sampling on an extreme arc is fine).
      if (steps_cache.size() >= steps_cache.capacity()) {
        HS_SCAN_METRIC(hs::g_scan_metrics.plot_backstop_hits++);
        break;
      }
      steps_cache.push_back(step);
      sim_dist += step;

      if (sim_dist < total_dist) {
        smp = sample(sim_dist / total_dist);
      }
    }

    // The final step normally overshoots total_dist (scale <= 1) and the
    // normalized replay stretches the cached steps back to exactly total_dist.
    // On the backstop break path sim_dist can fall short (scale > 1) and the
    // replay stretches over the remaining segment instead.
    HS_CHECK(sim_dist > 0.0f);
    float scale = total_dist / sim_dist;
    bool omitLast = close_loop || !isLastSegment || omit_end;

    // DRAWING PHASE
    //
    // sample().pos is ~0.04% non-unit; vector_to_pixel's phi = acos(v.y) offsets
    // the row near the pole, so re-normalize the interpolated positions.
    {
      Vector start_pos = sample(0.0f).pos.normalized();
      Fragment f = Fragment::lerp(curr, next, 0.0f);
      f.pos = start_pos;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, 0.0f);

      fragment_shader(start_pos, f);
      pipeline.plot(canvas, start_pos, f.color.color, f.age, f.color.alpha);
    }

    size_t loop_limit =
        omitLast ? steps_cache.size() - 1 : steps_cache.size();
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
      Vector p = sample(t).pos.normalized();
      Fragment f = Fragment::lerp(curr, next, t);
      f.pos = p;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, current_dist);

      fragment_shader(p, f);
      pipeline.plot(canvas, p, f.color.color, f.age, f.color.alpha);
    }
  };

  const auto &cr = canvas.clip();
  const bool clip_active = !cr.is_full();
  const auto xc = cr.x_clip();
  // Clip band as a cylindrical arc for the column cull.
  const int band_len = xc.wrap ? xc.re - xc.rs + W : xc.re - xc.rs;

  // Emits one shader-run dot for points[k]; the precomputed projection is
  // consumed only when no pipeline stage can re-position the plot.
  auto plot_dot = [&](const Fragment &src, size_t k) {
    Fragment f = src;
    f.color = Color4(0, 0, 0, 0);
    fragment_shader(src.pos, f);
    if constexpr (pipeline_hoistable_cull<PipelineT>()) {
      if (point_rows != nullptr && point_cols != nullptr) {
        pipeline.plot(canvas, point_cols[k], point_rows[k], f.color.color,
                      f.age, f.color.alpha);
        return;
      }
    }
    pipeline.plot(canvas, src.pos, f.color.color, f.age, f.color.alpha);
  };

  for (size_t i = 0; i < count; i++) {
    const Fragment &curr = points[i];
    const Fragment &next = points[(i + 1) % len];
    bool isLastSegment = (i == count - 1);

    // --- Interpolation Strategy Selection ---
    // Branch-cut guard: the planar projection is singular at the basis antipode,
    // so a segment with an endpoint there falls back to a geodesic edge. The seam
    // flag was decided once in the arc pre-pass (so the cached arc metric and the
    // rendered strategy cannot disagree) and is reused here, before the cull, so
    // the row-span bound matches the rendered arc shape.
    const bool antipodal_seam = override_uv && seg_seam_cache[i];
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
      const bool visible =
          edge_visible != nullptr
              ? edge_visible[i] != 0
              : edge_visible_in_clip<W, H>(
                    pipeline, cr, xc, band_len, curr.pos, next.pos,
                    use_planar ? planar_basis : nullptr);
      if (!visible)
        continue;
    }

    // Single-dot shortcut: an edge proven to span <= one screen step renders
    // exactly as process_segment's fast path (set_arc_uv is a no-op without a
    // planar basis), so plot it without building the sampler. A predicate
    // false negative falls through and re-evaluates exactly.
    if (!override_uv && edge_fits_one_dot<W, H>(curr.pos, next.pos)) {
      plot_dot(curr, i);
      if (!close_loop && isLastSegment && !omit_end)
        plot_dot(next, i + 1);
      continue;
    }

    if (use_planar) {
      rasterize_planar_strategy(curr, next, *planar_basis, isLastSegment,
                                process_segment);
    } else {
      rasterize_geodesic_strategy(curr, next, isLastSegment, process_segment);
    }
  }
  #ifdef __EMSCRIPTEN__
  canvas.add_render_us(emscripten_get_now() - plot_t0);
  #endif
}
HS_O3_END

/**
 * @brief Per-primitive geometry/rasterization options for draw_fragments.
 *
 * `close_loop` and `planar_basis` default to the common (geodesic, open) case,
 * so most primitives only spell out `.capacity`.
 */
struct FragmentDrawParams {
  size_t capacity;            /**< Fragment buffer reservation (per-primitive). */
  bool close_loop = false;    /**< Passed to rasterize (closes last→first edge). */
  bool omit_end = false;      /**< Skip the final endpoint plot of an open line, so abutting arcs tile without a double-plot. */
  const Basis *planar_basis = nullptr; /**< Planar projection basis (null = geodesic). */
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
  rasterize<W, H>(pipeline, canvas, points, fragment_shader, params.close_loop,
                  params.planar_basis, params.omit_end);
}

/**
 * @brief Draws a geodesic line between two points.
 * Registers:
 *  v0: line progress (0..1)
 *  v1: Arc Length (radians)
 */
struct Line {
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

    float angle = angle_between(f1.pos, f2.pos);
    // Same coincident-endpoint threshold as rasterize_geodesic_strategy: below
    // it the slerp axis is unstable, so collapse to a dot.
    if (angle < EPS_GEODESIC_SEGMENT) {
      Fragment f = f1;
      f.v0 = f.v1 = f.v2 = 0.0f;
      points.push_back(f);
      points.push_back(f); // draw at least a dot
      return;
    }

    Vector axis;
    if (std::abs(PI_F - angle) < TOLERANCE) {
      axis = stable_perpendicular_axis(f1.pos);
    } else {
      axis = cross(f1.pos, f2.pos).normalized();
    }

    for (int i = 0; i <= density; ++i) {
      float t = static_cast<float>(i) / density;

      Fragment f = Fragment::lerp(f1, f2, t);
      if (i == 0)
        f.pos = f1.pos;
      else if (i == density)
        f.pos = f2.pos;
      else {
        Quaternion q = make_rotation(axis, angle * t);
        f.pos = rotate(f1.pos, q);
      }

      f.v0 = t;
      f.v1 = angle * t;
      f.v2 = 0.0f;
      points.push_back(f);
    }
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
   * @param points Output fragment list; arc-length-parameterized fragments are
   *               appended.
   * @param vertices Iterable container of Fragment.
   * @param closed If true, connects the last point to the first.
   */
  static void sample(Fragments &points, const auto &vertices,
                     bool closed = false) {
    auto it = std::begin(vertices);
    auto end = std::end(vertices);

    if (it == end)
      return;

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
      points.push_back(f);
    }
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
    // sample() appends the wrap vertex carrying the v0=1 seam register;
    // close_loop drops the resulting degenerate closing edge.
    draw_fragments<W, H>(
        pipeline, canvas, vertex_shader, fragment_shader,
        {.capacity = vertices.size() + 2, .close_loop = closed},
        [&](Fragments &points) { sample(points, vertices, closed); });
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
 * position with v0 = 1 and the arc length continued across the wrap edge, so a
 * `close_loop` rasterize draws the final edge without a UV seam. Shared skeleton
 * for the accumulated-arc closed rings (Star, Flower, DistortedRing); Ring keeps
 * its own analytic-arc loop. For the PLANAR callers the rasterizer overrides
 * v0/v1, so these geodesic values seed only the optional vertex shader.
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
  if (points.size() > start_idx) {
    Fragment last = points[start_idx];
    cumulative_len += angle_between(points.back().pos, last.pos);
    last.v0 = 1.0f;
    last.v1 = cumulative_len;
    last.v2 = static_cast<float>(num_verts);
    points.push_back(last);
  }
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
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    const float theta_eq = work_radius * (PI_F / 2.0f);
    const float r_val = sinf(theta_eq);
    const float d_val = cosf(theta_eq);

    const float step = 2.0f * PI_F / num_samples;

    Vector first_pos; // Reused for the overlap-close vertex (i == 0 position).
    for (int i = 0; i < num_samples; i++) {
      float theta = i * step;
      float t = theta + phase;
      Vector u_temp = (u * cosf(t)) + (w * sinf(t));

      Fragment f;
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
    if (num_samples > 0) {
      Fragment f;
      f.pos = first_pos;
      f.v0 = 1.0f;
      f.v1 = 2.0f * PI_F * r_val;
      f.v2 = static_cast<float>(num_samples);
      f.age = 0;
      points.push_back(f);
    }
  }

  /**
   * @brief Full-resolution closed ring (W samples) — LUT-optimized.
   * @tparam W,H Rasterization resolution; W is the sample count and LUT grid.
   * @param points Output fragment list; W+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Ring radius (radians).
   * @param phase Rotation phase (radians).
   * @details For num_samples == W the angle grid (i*2π/W) is exactly
   * TrigLUT<W,H>::cos_theta/sin_theta, so per-sample cosf(θ+φ)/sinf(θ+φ) becomes
   * the precomputed θ-grid plus one angle-addition against cos/sin(φ), saving
   * ~2*(W+1) libm trig calls per ring per frame. Keeps Ring's analytic arc length
   * and overlap close. The runtime int-num_samples overload stays for the polygon
   * samplers, whose vertex counts do not match the LUT grid.
   */
  template <int W, int H>
  static void sample(Fragments &points, const Basis &basis, float radius,
                     float phase = 0) {
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    const float theta_eq = work_radius * (PI_F / 2.0f);
    const float r_val = sinf(theta_eq);
    const float d_val = cosf(theta_eq);

    const float step = 2.0f * PI_F / W;

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    const float cos_phase = cosf(phase);
    const float sin_phase = sinf(phase);

    for (int i = 0; i < W; i++) {
      Vector u_temp = ring_tangent<W, H>(i, u, w, cos_phase, sin_phase);

      Fragment f;
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
        pipeline, canvas, vertex_shader, fragment_shader, {.capacity = W + 2},
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
 * @brief Planar Polygon.
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians) — cumulative rendered planar arc
 *  v2: Vertex index
 * @note Always renders with PLANAR (azimuthal-equidistant) edges, which bow
 *       LONGER than the great-circle chord. The rasterizer re-derives v0/v1 from
 *       that true rendered arc, so both track the drawn position rather than the
 *       shorter chord polygon.
 */
struct PlanarPolygon {
  /**
   * @brief Samples a planar polygon.
   * @param points Output fragment list; num_sides+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Polygon radius.
   * @param num_sides Number of sides.
   * @param phase Rotation phase (radians).
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);
  }

  /**
   * @brief Draws a planar polygon.
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
    // Far-field radii (> 1) need the planar chart centered on the opposite pole,
    // else the azimuthal projection bows the polygon edges. radius <= 1 keeps the
    // supplied chart unchanged.
    Basis planar_basis = basis;
    if (radius > 1.0f) {
      planar_basis = planar_chart_basis(-basis.v);
    }

    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides + 2),
                          .close_loop = true, .planar_basis = &planar_basis},
                         [&](Fragments &points) {
                           sample(points, basis, radius, num_sides, phase);
                         });
  }

  /**
   * @brief Draws a planar polygon without a vertex shader.
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
 * @brief Spherical Polygon (Geodesic edges).
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Vertex index
 */
struct SphericalPolygon {
  /**
   * @brief Samples a spherical polygon.
   * @param points Output fragment list; num_sides+1 fragments are appended.
   * @param basis Orientation basis.
   * @param radius Polygon radius.
   * @param num_sides Number of sides.
   * @param phase Rotation phase (radians).
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    size_t start_idx = points.size();
    Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);

    // v1 = true geodesic chord length.
    float cumulative_length = 0.0f;
    for (size_t i = start_idx; i < points.size(); ++i) {
      points[i].v2 = static_cast<float>(i - start_idx);

      if (i > start_idx) {
        cumulative_length += angle_between(points[i - 1].pos, points[i].pos);
      }
      points[i].v1 = cumulative_length;
    }
  }

  /**
   * @brief Draws a spherical polygon.
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
    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides + 2),
                          .close_loop = true},
                         [&](Fragments &points) {
                           sample(points, basis, radius, num_sides, phase);
                         });
  }

  /**
   * @brief Draws a spherical polygon without a vertex shader.
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
   * @details Mirrors sample() exactly (at phase 0) so the returned point lands
   * on the drawn ring; any divergence would detach callers' sampled points from
   * the visible ring off Radius=1.
   */
  static Vector fn_point(ScalarFn shift_fn, const Basis &basis, float radius,
                         float angle) {
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    const float theta_eq = work_radius * (PI_F / 2.0f);
    const float polar = theta_eq + shift_fn(angle / (2.0f * PI_F));
    Vector u_temp = (u * cosf(angle)) + (w * sinf(angle));
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
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    const float theta_eq = work_radius * (PI_F / 2.0f);
    const float r_val = sinf(theta_eq);
    const float d_val = cosf(theta_eq);

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
                         {.capacity = W + 2, .close_loop = true},
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

  /**
   * @brief Samples a contiguous vertex arc [i_start, i_end] of the ring.
   * @tparam W,H Rasterization resolution (W is the full ring's vertex grid).
   * @param points Output fragment list; i_end - i_start + 1 fragments are
   * appended.
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param shift_fn Radial distortion sampled per vertex.
   * @param i_start First vertex index, in [0, W).
   * @param i_end Last vertex index (inclusive), in (i_start, W]: arcs never
   * cross the seam — a caller spanning it splits at vertex W and continues
   * from vertex 0.
   * @param phase Rotation phase (radians).
   * @details Every emitted vertex is bit-identical (position, v0, v2) to
   * sample()'s vertex at the same index — including i_end == W, which
   * reproduces the closed ring's overlap vertex (v0 = 1) — so an
   * arc-decomposed ring rasterizes the same pixels with the same interpolated
   * registers as the closed draw wherever both emit a segment. v1 accumulates
   * great-circle arc length from the arc's first vertex.
   */
  template <int W, int H>
  static void sample_arc(Fragments &points, const Basis &basis, float radius,
                         ScalarFn shift_fn, int i_start, int i_end,
                         float phase = 0) {
    HS_CHECK(i_start >= 0 && i_start < W && i_end > i_start && i_end <= W,
             "DistortedRing::sample_arc: bad arc range");
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    const float theta_eq = work_radius * (PI_F / 2.0f);
    const float r_val = sinf(theta_eq);
    const float d_val = cosf(theta_eq);

    const float step = 2.0f * PI_F / W;

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    const float cos_phase = cosf(phase);
    const float sin_phase = sinf(phase);

    float cumulative_len = 0.0f;
    for (int i = i_start; i <= i_end; ++i) {
      // idx folds only the i == W overlap vertex; theta and shift_fn's argument
      // must be the same floats the closed sample() produced, or the arc drifts
      // a ulp and rasterizes different pixels.
      int idx = i % W;
      float theta = idx * step;
      Vector u_temp = ring_tangent<W, H>(idx, u, w, cos_phase, sin_phase);

      float shift = shift_fn(theta / (2.0f * PI_F));
      float cos_shift = cosf(shift);
      float sin_shift = sinf(shift);

      float v_scale = d_val * cos_shift - r_val * sin_shift;
      float u_scale = r_val * cos_shift + d_val * sin_shift;

      Fragment f;
      f.pos = ((v * v_scale) + (u_temp * u_scale)).normalized();
      if (i > i_start)
        cumulative_len += angle_between(points.back().pos, f.pos);
      f.v0 = static_cast<float>(i) / W;
      f.v1 = cumulative_len;
      f.v2 = static_cast<float>(i);
      f.age = 0;
      points.push_back(f);
    }
  }

  /**
   * @brief Draws one open arc of a distorted ring.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param shift_fn Distortion function.
   * @param fragment_shader Shader function.
   * @param i_start First vertex index, in [0, W).
   * @param i_end Last vertex index (inclusive); see sample_arc's range contract.
   * @param phase Rotation phase.
   * @details Draws segments i_start .. i_end-1 and plots vertices
   * i_start .. i_end-1 (omit_end skips vertex i_end, which the closed draw
   * plots via its outgoing segment instead), so abutting arcs — or arcs plus
   * the closed draw's remaining segments — tile the ring with every vertex
   * plotted exactly once. A clip-culling caller must therefore only drop arcs
   * whose omitted end vertex is itself outside the clip.
   */
  template <int W, int H>
  static void draw_arc(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                       float radius, ScalarFn shift_fn,
                       FragmentShaderFn fragment_shader, int i_start, int i_end,
                       float phase = 0) {
    draw_fragments<W, H>(pipeline, canvas, {}, fragment_shader,
                         {.capacity = static_cast<size_t>(i_end - i_start + 2),
                          .close_loop = false,
                          .omit_end = true},
                         [&](Fragments &points) {
                           sample_arc<W, H>(points, basis, radius, shift_fn,
                                            i_start, i_end, phase);
                         });
  }

  /**
   * @brief Chunk-culled ring draw for partial clips: skips an unreachable ring
   * whole, otherwise cuts it into CHUNKS azimuth chunks whose bounding caps
   * gate the draw, emitting visible runs as open arcs.
   * @tparam W,H Rasterization resolution.
   * @tparam CHUNKS Azimuth chunks per ring; must divide W and fit a 32-bit mask.
   * @tparam BakeRun Callable as bake_run(int i0, int i1).
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param cr Clip region gating the cull (Effect::clip()).
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param shift_fn Distortion function.
   * @param fragment_shader Shader function.
   * @param reach Upper bound on |shift_fn| plus stroke/AA pad (radians); an
   * underestimate silently culls genuine arcs.
   * @param bake_run Called with the inclusive vertex range [i0, i1] whose
   * shift_fn samples the following draw reads (the full draw receives (0, W)),
   * so a LUT-backed shift_fn bakes only the ranges that render. A vertex's
   * parameter can float-round just below i0's LUT cell, so pad a bake one
   * cell below i0 (and one above i1 for the lerp neighbor), modulo the wrap.
   * @param phase Rotation phase.
   * @details Visible runs walk in ascending order — the same segment order the
   * closed draw rasterizes, which alpha blending's order dependence requires
   * wherever the ring's own taps overlap. A run touching the last chunk ends at
   * vertex W, so no arc crosses the seam, and arc vertices are bit-identical to
   * the closed ring's (see sample_arc), so clipped output matches full-frame
   * output exactly. A dropped chunk's caps clear the clip by reach, satisfying
   * draw_arc's contract that an omitted end vertex is itself outside the clip.
   */
  template <int W, int H, int CHUNKS = 24, typename BakeRun>
  static void draw_culled(PipelineRef pipeline, Canvas &canvas,
                          const ClipRegion &cr, const Basis &basis,
                          float radius, ScalarFn shift_fn,
                          FragmentShaderFn fragment_shader, float reach,
                          BakeRun bake_run, float phase = 0) {
    static_assert(W % CHUNKS == 0,
                  "chunk-to-vertex mapping requires W divisible by CHUNKS");
    static_assert(CHUNKS >= 1 && CHUNKS <= 31, "chunk mask is 32-bit");
    constexpr uint32_t FULL_MASK = (uint32_t{1} << CHUNKS) - 1;

    auto res = get_antipode(basis, radius);
    const Basis &work = res.first;
    float theta_eq = res.second * (PI_F / 2.0f);
    const bool try_cull = !cr.is_full();
    if (try_cull && !cap_may_touch_clip<H>(cr, work.v, theta_eq + reach))
      return;

    uint32_t vis = FULL_MASK;
    if (try_cull) {
      float cos_eq = cosf(theta_eq);
      float sin_eq = sinf(theta_eq);
      vis = 0;
      float chunk_reach = (PI_F / CHUNKS) * sin_eq + reach;
      for (int c = 0; c < CHUNKS; ++c) {
        // Midpoints carry the phase the drawn vertices apply via ring_tangent.
        float a = (c + 0.5f) * (2.0f * PI_F / CHUNKS) + phase;
        Vector mid = (work.v * cos_eq) +
                     ((work.u * cosf(a)) + (work.w * sinf(a))) * sin_eq;
        if (cap_may_touch_clip<H>(cr, mid, chunk_reach))
          vis |= 1u << c;
      }
      if (vis == 0)
        return;
    }

    if (vis == FULL_MASK) {
      bake_run(0, W);
      draw<W, H>(pipeline, canvas, basis, radius, shift_fn, fragment_shader,
                 phase);
      return;
    }

    constexpr int V = W / CHUNKS;
    int c = 0;
    while (c < CHUNKS) {
      while (c < CHUNKS && !(vis & (1u << c)))
        ++c;
      if (c == CHUNKS)
        break;
      int c0 = c;
      while (c < CHUNKS && (vis & (1u << c)))
        ++c;
      bake_run(c0 * V, c * V);
      draw_arc<W, H>(pipeline, canvas, basis, radius, shift_fn,
                     fragment_shader, c0 * V, c * V, phase);
    }
  }
};

/**
 * @brief Fibonacci Spiral primitive.
 */
struct Spiral {
  /**
   * @brief Samples a Fibonacci spiral.
   * @param fragments Output fragment list; n fragments are appended.
   * @param n Number of points.
   * @param eps Epsilon offset shifting points off the poles.
   */
  static void sample(Fragments &fragments, int n, float eps) {
    HS_CHECK(n >= 1);
    float cumulative_len = 0.0f;
    Vector last_pos;

    for (int i = 0; i < n; i++) {
      Vector pos = fib_spiral(n, eps, i);

      if (i > 0) {
        cumulative_len += angle_between(last_pos, pos);
      }
      last_pos = pos;

      Fragment f;
      f.pos = pos;
      f.v0 = (n > 1) ? static_cast<float>(i) / (n - 1) : 0.0f;
      f.v1 = cumulative_len;
      f.v2 = static_cast<float>(i);
      f.age = 0;
      fragments.push_back(f);
    }
  }

  /**
   * @brief Draws a Fibonacci spiral.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param n Number of points.
   * @param eps Epsilon offset shifting points off the poles.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, int n, float eps,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = n},
                         [&](Fragments &frags) { sample(frags, n, eps); });
  }

  /**
   * @brief Draws a Fibonacci spiral without a vertex shader.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param n Number of points.
   * @param eps Epsilon offset shifting points off the poles.
   * @param fragment_shader Shader function.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, int n, float eps,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, n, eps, fragment_shader, {});
  }
};

/**
 * @brief Star shape.
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians) — cumulative rendered planar arc
 *  v2: Vertex index
 * @note Always renders with PLANAR (azimuthal-equidistant) edges, which bow
 *       LONGER than the great-circle chord. The rasterizer re-derives v0/v1 from
 *       that true rendered arc, so both track the drawn position rather than the
 *       shorter chord polygon.
 */
struct Star {
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
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    float outer_radius = work_radius * (PI_F / 2.0f);
    float inner_radius = outer_radius * STAR_INNER_RATIO;
    float angle_step = PI_F / num_sides;

    // Alternating outer/inner radius per vertex; everything else is the shared
    // closed-ring skeleton.
    sample_closed_ring(points, num_sides * 2, [&](int i) {
      float theta = phase + i * angle_step;
      float r = (i % 2 == 0) ? outer_radius : inner_radius;
      float sin_r = sinf(r);
      float cos_r = cosf(r);
      float cos_t = cosf(theta);
      float sin_t = sinf(theta);
      Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
      p.normalize();
      return p;
    });
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
    Basis planar_basis = planar_chart_basis(get_antipode(basis, radius).first.v);

    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides * 2 + 2),
                          .close_loop = true, .planar_basis = &planar_basis},
                         [&](Fragments &points) {
                           sample(points, basis, radius, num_sides, phase);
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
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    float desired_outer_radius = work_radius * (PI_F / 2.0f);
    float apothem = PI_F - desired_outer_radius;
    float safe_apothem = std::min(apothem, PI_F - 1e-4f);
    float angle_step = PI_F / num_sides;

    // Constant polar radius per vertex; everything else is the shared closed-
    // ring skeleton.
    sample_closed_ring(points, num_sides * 2, [&](int i) {
      float theta = phase + i * angle_step;
      float sin_r = sinf(safe_apothem);
      float cos_r = cosf(safe_apothem);
      float cos_t = cosf(theta);
      float sin_t = sinf(theta);
      Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
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
    Basis planar_basis = planar_chart_basis(get_antipode(basis, radius).first.v);

    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = static_cast<size_t>(num_sides * 2 + 2),
                          .close_loop = true, .planar_basis = &planar_basis},
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
   * @details A mesh exceeding this traps on the cold setup path. Sized for a
   * TriangularBitset of 128*127/2 bits = 1016 bytes.
   */
  static constexpr int DEDUP_CAPACITY = 128;

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
   * @param vertex_shader Optional vertex shader.
   * @details Shared body for both draw() overloads (face-walk and precomputed
   * edge list); keeping it in one place is why the two paths stay bit-identical.
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

    ScratchScope edge_guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, 16);
    Line::sample(points, fu, fv, 10);

    if (vertex_shader) {
      for (auto &p : points) {
        p.v2 = static_cast<float>(edge_index); // Edge Index
        vertex_shader(p);
      }
    } else {
      for (auto &p : points) {
        p.v2 = static_cast<float>(edge_index);
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, false, nullptr);
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
        // trap on the cold setup path rather than drop the edge.
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
    auto &visited = *new (scratch_arena_b.allocate(
        sizeof(TriangularBitset<DEDUP_CAPACITY>),
        alignof(TriangularBitset<DEDUP_CAPACITY>))) TriangularBitset<DEDUP_CAPACITY>();

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
    auto &visited = *new (scratch_arena_b.allocate(
        sizeof(TriangularBitset<DEDUP_CAPACITY>),
        alignof(TriangularBitset<DEDUP_CAPACITY>))) TriangularBitset<DEDUP_CAPACITY>();

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
   */
  template <int W, int H, typename MeshT, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshT &mesh,
                   const ArenaVector<Edge> &edges,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader = {}) {
    for (size_t ei = 0; ei < edges.size(); ++ei) {
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
 * @brief Gates one geodesic trail's edges against the clip in one hoisted pass.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam PipelineT Pipeline type; must have no world cull stage
 *         (pipeline_hoistable_cull), so the predicate sees the raw points.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param band_len Clip band length in columns (seam-unwrapped).
 * @param trail Geodesic fragment polyline (>= 2 unit-position points).
 * @param bits Output, one byte per edge (trail.size() - 1): 0 = culled, else
 *        1; valid as rasterize()'s edge_visible input.
 * @return False when no edge is visible; bits are then all zero.
 * @details Per-edge verdicts are exactly edge_visible_in_clip's (rasterize
 * consumes them as its cull), with per-point rows/columns hoisted across
 * adjacent edges and a conservative whole-trail bound rejecting fully-
 * invisible trails first — the same scheme as ParticleSystem::draw's gate.
 */
HS_O3_BEGIN
template <int W, int H, typename PipelineT>
static bool gate_trail_edges(const PipelineT &, const ClipRegion &cr,
                             const ClipRegion::XClip &xc, int band_len,
                             const Fragments &trail, uint8_t *bits) {
  static_assert(pipeline_hoistable_cull<PipelineT>(),
                "gate_trail_edges requires a pipeline with no world cull "
                "stage; route others through edge_visible_in_clip");
  constexpr int H_VIRT = H + hs::H_OFFSET;
  constexpr float MIN_SIN_PHI = 0.05f;
  const size_t n = trail.size();
  HS_CHECK(n >= 2);
  const size_t edges = n - 1;

  ScratchScope span_guard(scratch_arena_a);
  auto *rows = static_cast<float *>(
      scratch_arena_a.allocate(n * sizeof(float), alignof(float)));
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
  }
  // arc <= (pi/2)*chord on [0, pi]; an edge's interior latitude extremum lies
  // within arc/2 of an endpoint and phi is 1-Lipschitz in arc length, so this
  // margin covers every per-edge bulge peak.
  const float max_arc = (PI_F * 0.5f) * sqrtf(max_chord2);
  const float row_margin =
      (max_arc * 0.5f) * (static_cast<float>(H_VIRT - 1) / PI_F);
  if (!cr.could_intersect_y(row_lo_t - row_margin, row_hi_t + row_margin)) {
    std::fill_n(bits, edges, uint8_t{0});
    return false;
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
    }
    // Near a pole the plotted column is float noise (same caution as the
    // per-edge spans), so only cull by the column arc when the whole trail
    // provably stays clear.
    if (walk_safe && sqrtf(std::max(0.0f, min_sp2)) - max_arc >= MIN_SIN_PHI) {
      int col_s, col_len;
      finish_col_span<W>(cols[0] + cum_lo, cum_hi - cum_lo, col_s, col_len);
      if (!ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len, W)) {
        std::fill_n(bits, edges, uint8_t{0});
        return false;
      }
    }
  }

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
    float row_lo, row_hi;
    geodesic_row_span_rows<W, H>(rows[e], rows[e + 1], ea, eb, es, row_lo,
                                 row_hi);
    bool v;
    if (!cr.could_intersect_y(row_lo, row_hi)) {
      v = false;
    } else if (!xc.active) {
      v = true;
    } else {
      int col_s, col_len;
      v = !geodesic_col_span_cols<W>(cols[e], cols[e + 1], ea, es, col_s,
                                     col_len) ||
          ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len, W);
    }
    bits[e] = v ? 1 : 0;
    any = any || v;
  }
  return any;
}
HS_O3_END

/**
 * @brief Particle System trails.
 * Registers:
 *  v0: Trail Progress (0.0=Head -> 1.0=Tail)
 *  v1: Reserved (always 0)
 *  v2: Particle ID
 *  v3: Normalized TTL
 */
struct ParticleSystem {
  /**
   * @brief Draws each active particle's history as a rasterized trail.
   * @tparam W,H Rasterization resolution.
   * @tparam PipelineT Pipeline type (defaults to PipelineRef).
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
   */
  template <int W, int H, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader,
                   DeferredShaderRef deferred_shader = {}) {
    int count = system.active();
    if (count == 0)
      return;

    const float max_life = static_cast<float>(system.max_life);
    HS_CHECK(std::isfinite(max_life) && max_life >= 1.0f &&
                 max_life <= 65535.0f,
             "ParticleSystem render max_life must be finite and in [1, 65535]");
    const float inv_max_life = 1.0f / max_life;

    // Segment-clip state for the trail-level deferred-shader gate below.
    const auto &cr = canvas.clip();
    const bool clip_active = !cr.is_full();
    const auto xc = cr.x_clip();
    const int band_len = xc.wrap ? xc.re - xc.rs + W : xc.re - xc.rs;

    for (int i = 0; i < count; ++i) {
      const auto &p = system.pool[i];
      ScratchScope trail_guard(scratch_arena_a);
      Fragments trail;
      // tween emits at most one fragment per retained trail position.
      trail.bind(scratch_arena_a,
                 std::remove_cvref_t<decltype(p.history)>::CAPACITY);
      // Original (pre-shader) positions, kept for the deferred pass.
      ArenaVector<Vector> orig;
      if (deferred_shader)
        orig.bind(scratch_arena_a,
                  std::remove_cvref_t<decltype(p.history)>::CAPACITY);
      {
        HS_PROFILE(plot_ps_tween);
        tween(p.history, [&](const Vector &v, float t) {
          Fragment f;
          f.pos = v;
          f.v0 = t;
          f.v1 = 0.0f;
          f.v2 = static_cast<float>(i);
          f.v3 = static_cast<float>(p.life) * inv_max_life;
          f.age = 0;
          f.color = Color4(0, 0, 0, 0);
          trail.push_back(f);
          if (deferred_shader)
            orig.push_back(v);
        });
      }

      if (trail.is_empty())
        continue;
      {
        HS_PROFILE(plot_ps_vertex);
        apply_vertex_shader(vertex_shader, trail);
      }

      if (!deferred_shader) {
        HS_PROFILE(plot_ps_raster);
        rasterize<W, H>(pipeline, canvas, trail, fragment_shader, false,
                        nullptr);
        continue;
      }

      // Trail-level gate: precompute each edge's cull verdict from the
      // position-shaded points. No visible edge means the trail renders
      // nothing, so the deferred pass and the rasterize call are skipped
      // whole; the bits feed rasterize so the cull is evaluated once.
      const uint8_t *vis = nullptr;
      const float *dot_rows = nullptr;
      const float *dot_cols = nullptr;
      if (clip_active && trail.size() >= 2) {
        HS_PROFILE(plot_ps_gate);
        const size_t edges = trail.size() - 1;
        auto *bits = static_cast<uint8_t *>(
            scratch_arena_a.allocate(edges, alignof(uint8_t)));
        bool any = false;
        if constexpr (pipeline_hoistable_cull<PipelineT>()) {
          // No stage re-emits edges, so the predicate sees the raw points:
          // per-point rows/columns are computed once and shared by every edge,
          // and a conservative whole-trail bound rejects fully-invisible
          // trails before any per-edge work. The per-edge pass must produce
          // exactly the bits edge_visible_in_clip would (rasterize consumes
          // them as its cull).
          constexpr int H_VIRT = H + hs::H_OFFSET;
          constexpr float MIN_SIN_PHI = 0.05f;
          const size_t n = trail.size();
          auto *rows = static_cast<float *>(
              scratch_arena_a.allocate(n * sizeof(float), alignof(float)));
          dot_rows = rows;
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
          }
          // arc <= (pi/2)*chord on [0, pi]; an edge's interior latitude
          // extremum lies within arc/2 of an endpoint and phi is 1-Lipschitz
          // in arc length, so this margin covers every per-edge bulge peak.
          const float max_arc = (PI_F * 0.5f) * sqrtf(max_chord2);
          const float row_margin =
              (max_arc * 0.5f) * (static_cast<float>(H_VIRT - 1) / PI_F);
          if (!cr.could_intersect_y(row_lo_t - row_margin,
                                    row_hi_t + row_margin))
            continue;

          float *cols = nullptr;
          if (xc.active) {
            cols = static_cast<float *>(
                scratch_arena_a.allocate(n * sizeof(float), alignof(float)));
            dot_cols = cols;
            float cum = 0.0f, cum_lo = 0.0f, cum_hi = 0.0f;
            bool walk_safe = true;
            cols[0] = vector_to_theta<W>(trail[0].pos);
            for (size_t k = 1; k < n; ++k) {
              cols[k] = vector_to_theta<W>(trail[k].pos);
              // A geodesic edge's column sweep never exceeds W/2 (antipodal
              // symmetry, see geodesic_col_span_cols), so the short-way delta
              // covers it regardless of direction — except at ~exactly W/2,
              // where the delta's sign (which semicircle) is float noise.
              float d = cols[k] - cols[k - 1];
              if (d > W * 0.5f)
                d -= W;
              else if (d < -W * 0.5f)
                d += W;
              if (std::abs(d) >= W * 0.5f - 3.0f)
                walk_safe = false;
              cum += d;
              cum_lo = std::min(cum_lo, cum);
              cum_hi = std::max(cum_hi, cum);
            }
            // Near a pole the plotted column is float noise (same caution as
            // the per-edge spans), so only cull by the column arc when the
            // whole trail provably stays clear.
            if (walk_safe &&
                sqrtf(std::max(0.0f, min_sp2)) - max_arc >= MIN_SIN_PHI) {
              int col_s, col_len;
              finish_col_span<W>(cols[0] + cum_lo, cum_hi - cum_lo, col_s,
                                 col_len);
              if (!ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len,
                                            W))
                continue;
            }
          }

          for (size_t e = 0; e < edges; ++e) {
            const Vector &ea = trail[e].pos;
            const Vector &eb = trail[e + 1].pos;
            const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
            float row_lo, row_hi;
            geodesic_row_span_rows<W, H>(rows[e], rows[e + 1], ea, eb, es,
                                         row_lo, row_hi);
            bool v;
            if (!cr.could_intersect_y(row_lo, row_hi)) {
              v = false;
            } else if (!xc.active) {
              v = true;
            } else {
              int col_s, col_len;
              v = !geodesic_col_span_cols<W>(cols[e], cols[e + 1], ea, es,
                                             col_s, col_len) ||
                  ClipRegion::arcs_overlap(xc.rs, band_len, col_s, col_len, W);
            }
            bits[e] = v ? 1 : 0;
            any = any || v;
          }
        } else {
          for (size_t e = 0; e < edges; ++e) {
            bits[e] = edge_visible_in_clip<W, H>(pipeline, cr, xc, band_len,
                                                 trail[e].pos, trail[e + 1].pos,
                                                 nullptr)
                          ? 1
                          : 0;
            any = any || bits[e] != 0;
          }
        }
        if (!any)
          continue;
        vis = bits;
      }

      {
        HS_PROFILE(plot_ps_deferred);
        for (size_t k = 0; k < trail.size(); ++k)
          deferred_shader(trail[k], orig[k]);
      }
      {
        HS_PROFILE(plot_ps_raster);
        rasterize<W, H>(pipeline, canvas, trail, fragment_shader, false,
                        nullptr, false, vis, dot_rows, dot_cols);
      }
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
