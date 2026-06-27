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
#include "geometry.h"
#include "color.h"
#include "constants.h"
#include "canvas.h"
#include "animation.h"

namespace Plot {

/**
 * @brief Inner/outer radius ratio for star shapes.
 * @details 1/φ² ≈ 0.382; mirrors SDF::STAR_INNER_RATIO.
 */
static constexpr float STAR_INNER_RATIO = 0.382f;

/**
 * @brief Geodesic segment shorter than this (radians) collapses to a point.
 * @details Deliberately ~10× math::EPS_GEOMETRIC — a looser threshold for
 * slerp-axis stability, not positional near-equality.
 */
static constexpr float EPS_GEODESIC_SEGMENT = 0.001f;

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
 * @details Zero-cost inline replacement for the identical
 * `if (vertex_shader) for (auto &p : pts) vertex_shader(p);` block repeated
 * across the primitives. The FunctionRef is passed by value (two pointers) and
 * the whole thing inlines away at -O3.
 */
template <typename Fragments>
inline void apply_vertex_shader(VertexShaderRef vertex_shader, Fragments &pts) {
  if (vertex_shader) {
    for (auto &p : pts) {
      vertex_shader(p);
    }
  }
}

/**
 * @brief Single point primitive.
 */
struct Point {
  /**
   * @brief Draws a single point.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param f Fragment to plot.
   * @param fragment_shader Shader function applied before plotting; must be
   *        non-null (a default-constructed FunctionRef would invoke a null
   *        thunk under NDEBUG, where operator()'s assert is stripped).
   */
  static void draw(PipelineRef pipeline, Canvas &canvas, const Fragment &f,
                   FragmentShaderFn fragment_shader) {
    // Cold (once per point), not per-pixel: trap a null shader here so it fails
    // deterministically instead of calling a null thunk under NDEBUG.
    HS_CHECK(fragment_shader, "Point::draw requires a non-null fragment_shader");
    Fragment f_copy = f;
    f_copy.color = Color4(0, 0, 0, 0);
    fragment_shader(f_copy.pos, f_copy);
    pipeline.plot(canvas, f_copy.pos, f_copy.color.color, f_copy.age,
                  f_copy.color.alpha);
  }
};

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
 * sample, used by screen_step() to size the next sub-step from screen velocity.
 * The geodesic strategy fills it in closed form (free from the slerp's sin/cos);
 * the planar strategy finite-differences it. Zero for a degenerate edge: a
 * coincident one (< EPS_GEOMETRIC) is plotted as a single dot before any
 * sampling, while a near-coincident geodesic one (< EPS_GEODESIC_SEGMENT, where
 * the slerp axis is unstable) does reach screen_step, whose speed floor maps the
 * zero tangent to a base_step (one-dot) step.
 */
struct SamplePT {
  Vector pos;
  Vector tan;
};

constexpr int kPlanarLenSamples = 4;

/**
 * @brief Cumulative on-sphere arc length at kPlanarLenSamples+1 evenly-spaced
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
                 std::array<float, kPlanarLenSamples + 1> &arc_cumul) {
  arc_cumul[0] = 0.0f;
  Vector prev = azimuthal_unproject(proj.first, proj.second, planar_basis);
  for (int k = 1; k <= kPlanarLenSamples; ++k) {
    float p = static_cast<float>(k) / kPlanarLenSamples;
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
 * Projection-uniform stepping is NOT arc-uniform under the anisotropic
 * azimuthal metric, so a short cumulative-arc table (LEN_SAMPLES samples) is
 * inverted to turn an arc-length fraction into a projection parameter, making
 * planar sampling arc-uniform (matching the geodesic strategy) with no new trig.
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

  // Cumulative on-sphere arc length at evenly-spaced PROJECTION samples.
  // total_dist drives the step count, so it must be the path's true on-sphere
  // length: the projected chord over-estimates it (the stretched tangential
  // metric) while the geodesic arc under-estimates it. Summing great-circle
  // arcs between a few samples gives the real length.
  //
  // The table additionally lets map() take an ARC-length fraction and invert it
  // back to a projection parameter. Projection-uniform stepping is NOT
  // arc-uniform under the anisotropic azimuthal metric, so feeding the
  // rasterizer's arc-fraction t straight into the projection-linear map would
  // cluster/gap samples. The inversion makes planar sampling arc-uniform —
  // matching the geodesic strategy — at the cost of a trig-free table scan per
  // sample (no new trig anywhere).
  std::array<float, kPlanarLenSamples + 1> arc_cumul;
  planar_arc_cumul(proj1, dx, dy, planar_basis, arc_cumul);
  const float dist = arc_cumul[kPlanarLenSamples];

  // Arc-length parameterization: s is the arc fraction in [0,1]. Invert the
  // piecewise-linear cumulative-arc table to a projection parameter, then
  // unproject. A short scan over kPlanarLenSamples floats — no trig.
  auto map_planar = [=](float s) -> Vector {
    if (dist < math::EPS_GEOMETRIC)
      return unproject(s);
    float target = s * dist;
    int k = 0;
    while (k < kPlanarLenSamples - 1 && arc_cumul[k + 1] < target)
      ++k;
    float seg = arc_cumul[k + 1] - arc_cumul[k];
    float frac =
        (seg > math::EPS_GEOMETRIC) ? (target - arc_cumul[k]) / seg : 0.0f;
    float p = (static_cast<float>(k) + frac) / kPlanarLenSamples;
    return unproject(std::min(1.0f, std::max(0.0f, p)));
  };

  // The azimuthal map has no closed-form tangent, so take it from a short forward
  // difference (backward at the s=1 end) — map_planar is arc-uniform, so the
  // difference direction is the unit tangent regardless of its magnitude. Feeds
  // the same screen-velocity sub-step sampler as the geodesic path.
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
  std::array<float, kPlanarLenSamples + 1> arc_cumul;
  planar_arc_cumul(p1, p2.first - p1.first, p2.second - p1.second, planar_basis,
                   arc_cumul);
  return arc_cumul[kPlanarLenSamples];
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
      axis = (std::abs(dot(v1, X_AXIS)) > math::COS_AXIS_PARALLEL)
                 ? cross(v1, Y_AXIS)
                 : cross(v1, X_AXIS);
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

/**
 * @brief Conservative screen-row span of a rendered edge, arc bulge included.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param planar_basis Non-null: edge is azimuthal-equidistant; null: geodesic.
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 * @details The clip cull must not test endpoints alone: a great-circle (or
 * planar-projected) edge between two points outside a clip band can still bulge
 * through it, and an endpoint-only test silently drops the arc — a gap at a
 * segment boundary on Phantasm hardware, where each board renders a Y-band.
 * Rows are computed via the canvas mapping row = phi_to_y(acos(y)); the endpoint
 * rows are extended by the arc's interior latitude extremum (closed-form for the
 * geodesic case, sampled with a Lipschitz margin for the planar case). Runs once
 * per coarse edge on the clip-only path.
 */
template <int W, int H>
static inline void edge_row_span(const Vector &a, const Vector &b,
                                 const Basis *planar_basis, float &row_lo,
                                 float &row_hi) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  auto y_to_row = [](float y) {
    return phi_to_y(fast_acos(hs::clamp(y, -1.0f, 1.0f)), H_VIRT);
  };
  float ra = y_to_row(a.y);
  float rb = y_to_row(b.y);
  row_lo = std::min(ra, rb);
  row_hi = std::max(ra, rb);

  if (planar_basis == nullptr) {
    // Geodesic edge: the arc y(t) = a.y·cos + v_perp.y·sin has its turning
    // point inside the span iff the forward tangent's y-component flips sign
    // between the endpoints. The extremal |y| is the great circle's peak
    // latitude sqrt(1 - n.y²) (n = arc pole), reached only when present.
    // No one-row epsilon here (unlike the planar branch below): the span is the
    // exact closed-form y range of the rendered arc — the endpoint y values are
    // the rendered endpoints, and y_to_row is monotonic in y — so there is no
    // project/unproject round-trip to leave the cull a sub-pixel short.
    Vector axis = cross(a, b);
    float L2 = dot(axis, axis);
    if (L2 > math::EPS_GEOMETRIC * math::EPS_GEOMETRIC) {
      Vector n = axis * (1.0f / sqrtf(L2));
      float t0 = cross(n, a).y; // forward tangent y at a
      float t1 = cross(n, b).y; // forward tangent y at b
      if ((t0 > 0.0f) != (t1 > 0.0f)) {
        // std::max(0, ...) absorbs the tiny negative that fast-math
        // renormalization of n can produce when |n.y| ≈ 1 (a near-polar
        // arc pole), keeping the sqrt domain-safe.
        float peak = sqrtf(std::max(0.0f, 1.0f - n.y * n.y));
        float rp = y_to_row(t0 > 0.0f ? peak : -peak);
        row_lo = std::min(row_lo, rp);
        row_hi = std::max(row_hi, rp);
      }
    }
  } else {
    // Planar edge: an azimuthal-equidistant straight line is not a great circle,
    // so there is no closed-form latitude extremum. Sample the arc through the
    // same unprojection MAP the rasterizer uses, then widen by the arc's
    // Lipschitz bound. The cull and the renderer do NOT take bit-identical
    // samples — the cull steps uniformly in PROJECTION space (p = k/SAMPLES)
    // while the renderer sub-steps adaptively in ARC length — so gap-freeness
    // comes from the Lipschitz + one-row margin below, not from matching sample
    // points. Working in row space keeps the margin uniform (it would blow up
    // near the poles in y): phi is 1-Lipschitz in angular distance, so between
    // samples |Δrow| ≤ (Δarc)·(H_VIRT−1)/π, and the projected chord length
    // over-estimates the on-sphere arc (the projection stretches tangential
    // distance). That margin makes the span provably gap-free without dense
    // sampling.
    auto p1 = azimuthal_project(a, *planar_basis);
    auto p2 = azimuthal_project(b, *planar_basis);
    float dX = p2.first - p1.first;
    float dY = p2.second - p1.second;
    constexpr int SAMPLES = 8;
    for (int k = 1; k < SAMPLES; ++k) {
      float p = static_cast<float>(k) / SAMPLES;
      float y = azimuthal_unproject(p1.first + dX * p, p1.second + dY * p,
                                    *planar_basis)
                    .y;
      float r = y_to_row(y);
      row_lo = std::min(row_lo, r);
      row_hi = std::max(row_hi, r);
    }
    // Lipschitz margin (covers the hump between samples) plus a one-row epsilon
    // that absorbs the sub-pixel difference between the unprojected sample (≈unit
    // to fast-math precision) and the renderer's normalized plot position.
    float margin = (sqrtf(dX * dX + dY * dY) / SAMPLES) *
                       (static_cast<float>(H_VIRT - 1) / PI_F) +
                   1.0f;
    row_lo -= margin;
    row_hi += margin;
  }
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
  const float step = SCREEN_STEP_PX / std::max(speed, 1e-6f);
  return std::max(base_step * MIN_POLE_SCALE, std::min(step, base_step));
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
 */
template <int W, int H, typename PipelineT = PipelineRef>
static void rasterize(PipelineT &pipeline, Canvas &canvas,
                      const Fragments &points, FragmentShaderFn fragment_shader,
                      bool close_loop = false,
                      const Basis *planar_basis = nullptr) {
  size_t len = points.size();
  if (len < 2)
    return;
  // Cold (once per polyline), not per-pixel: trap a null shader here so the
  // many per-pixel fragment_shader() calls below can't invoke a null thunk.
  HS_CHECK(fragment_shader, "rasterize requires a non-null fragment_shader");
  #ifdef __EMSCRIPTEN__
  double plot_t0 = emscripten_get_now();
  #endif

  size_t count = close_loop ? len : len - 1;
  // SCRATCH ARENA CONTRACT (load-bearing): this rasterizer and
  // Pixel::Feedback::flush (filter.h) both checkpoint scratch_arena_a. Sharing
  // it is safe by construction — scratch_arena_a is a LIFO bump allocator, so
  // stack-nested scopes (whether plot-inside-flush or flush-inside-plot) rewind
  // in order and never clobber a live allocation; ~ScratchScope HS_CHECKs that
  // LIFO discipline (memory.h). The residual rule the allocator can't enforce:
  // do not let a raw scratch_arena_a pointer outlive the scope that produced it
  // (e.g. across a plot/flush boundary) — read it only within its scope.
  ScratchScope sc_guard(scratch_arena_a);
  ArenaVector<float> steps_cache;
  // The cache holds ONE segment's adaptive sub-steps (cleared per segment). Each
  // step advances ≈ one screen pixel, so a single segment needs ≲ W steps (the
  // step·≥MIN_POLE_SCALE lower clamp caps the per-segment count at the poles).
  // Size off W with 2× headroom; the simulation loop breaks at capacity as a
  // backstop.
  size_t max_cache = std::max((size_t)64, (size_t)(2 * W));
  steps_cache.bind(scratch_arena_a, max_cache);

  // PLANAR ARC REGISTERS (v0/v1). The sample functions store v0 as a vertex
  // fraction and v1 as the GEODESIC chord cumulant. Under planar (azimuthal-
  // equidistant) interpolation the rendered edge bows LONGER than that chord, so
  // the stored values disagree with the drawn position. When a planar basis is in
  // force, re-derive v0/v1 from the TRUE rendered arc: `cumul` tracks the rendered
  // arc reached so far, `seg_base` snapshots a segment's start for the draw
  // lambda, and `total_arc` (one cold pre-pass over the segments, matching the
  // per-segment strategy the draw loop picks) normalizes v0 to [0,1]. Geodesic
  // polylines already carry the rendered arc, so this is skipped for them.
  const bool override_uv = (planar_basis != nullptr);
  float total_arc = 0.0f;
  // Per-segment rendered arc length, reused by the draw-loop cumulant below.
  ArenaVector<float> seg_arc_cache;
  if (override_uv) {
    seg_arc_cache.bind(scratch_arena_a, count);
    const Vector &pcenter = planar_basis->v;
    for (size_t i = 0; i < count; i++) {
      const Vector &a = points[i].pos;
      const Vector &b = points[(i + 1) % len].pos;
      const bool seam = dot(a, pcenter) < -COS_PLANAR_ANTIPODE ||
                        dot(b, pcenter) < -COS_PLANAR_ANTIPODE;
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
    // The degenerate and fast paths below plot curr.pos/next.pos directly, without
    // the renormalize the DRAWING PHASE applies further down: these are the
    // original sampled vertices (already unit), not the fast_sinf/cosf-interpolated
    // sample().pos outputs that drift ~0.04% off the unit sphere. Precondition:
    // callers pass unit fragment positions.
    // Degenerate (coincident endpoints): plot at most a single dot.
    if (total_dist < math::EPS_GEOMETRIC) {
      bool shouldOmit = (close_loop) ? true : !isLastSegment;
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

    // FAST PATH: the whole segment spans ≤ one screen step (~SCREEN_STEP_PX px),
    // so a single dot covers it — skip the simulation. Keyed on SCREEN length
    // (via screen_step), NOT arc length: a base_step arc can still cross several
    // pixels on a steep or near-polar segment, where an arc-length test (the old
    // `total_dist < base_step`) would skip the sub-stepping the curve needs and
    // render it as a beaded line.
    if (total_dist <= first_step) {
      Fragment f = curr;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, 0.0f);
      fragment_shader(curr.pos, f);
      pipeline.plot(canvas, curr.pos, f.color.color, f.age, f.color.alpha);
      if (!close_loop && isLastSegment) {
        Fragment fl = next;
        fl.color = Color4(0, 0, 0, 0);
        set_arc_uv(fl, total_dist);
        fragment_shader(next.pos, fl);
        pipeline.plot(canvas, next.pos, fl.color.color, fl.age,
                      fl.color.alpha);
      }
      return;
    }

    // 1. SIMULATION PHASE — size each sub-step so consecutive samples land
    // ~SCREEN_STEP_PX apart in SCREEN space (screen_step), using the strategy's
    // unit tangent at the step's start. This tracks the true 2-D screen speed
    // rather than the old sin(φ) longitudinal-only proxy, so the curve is
    // sampled ~one pixel per step everywhere instead of unevenly (see
    // screen_step's note). `smp`/`first_step` above seed the first iteration.
    steps_cache.clear();
    float sim_dist = 0.0f;

    while (sim_dist < total_dist) {
      float step = screen_step<W, H>(smp.pos, smp.tan, base_step);

      // Backstop: a pathological segment (e.g. a huge-radius shape wrapping the
      // sphere) could still exceed the 2*W cache. Stop subdividing and let the
      // normalized replay stretch the cached steps over the rest of the segment
      // — coarser sampling on an extreme arc is fine (and far better than
      // trapping a live show). Normal segments never reach this.
      if (steps_cache.size() >= steps_cache.capacity())
        break;
      steps_cache.push_back(step);
      sim_dist += step;

      if (sim_dist < total_dist) {
        smp = sample(sim_dist / total_dist);
      }
    }

    // The final step overshoots total_dist slightly, so sim_dist >= total_dist
    // and scale <= 1 — the normalized replay stretches the cached steps back to
    // exactly total_dist, correcting the overshoot.
    HS_CHECK(sim_dist > 0.0f);
    float scale = total_dist / sim_dist;
    bool omitLast = (close_loop) ? true : !isLastSegment;

    // 2. DRAWING PHASE
    //
    // sample()'s pos builds interpolated points with fast_sinf/fast_cosf,
    // which don't satisfy sin²+cos²=1 — so the result is ~0.04% non-unit. The
    // sink's vector_to_pixel skips normalization and takes phi = acos(v.y)
    // directly; near the pole acos has infinite slope, so a y of 0.9996 instead
    // of 1.0 lands ~1.3 rows below the pole. Re-normalize the interpolated
    // positions here (the sampled vertices are already unit) so polar samples
    // map to the correct row. Plot-path only; one sqrt per drawn sample.
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

      // total_dist > 0 here: HS_CHECK(sim_dist > 0) above implies >=1 sim step,
      // which the loop only takes while sim_dist < total_dist.
      float t = current_dist / total_dist;

      // `t` (hence the drawn POSITION) is parameterized by the RENDERED arc
      // length. The registers are lerped from the control points; under a planar
      // basis set_arc_uv then rewrites v0/v1 from the true rendered arc (current_dist
      // is the planar arc walked so far, since the planar sampler is arc-uniform),
      // so a shader keying off v1/v0 as an arc-length proxy tracks the drawn
      // position even across the planar edge's bow away from the great-circle
      // chord. Geodesic edges keep the lerped registers (already the rendered arc).
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

  for (size_t i = 0; i < count; i++) {
    const Fragment &curr = points[i];
    const Fragment &next = points[(i + 1) % len];
    bool isLastSegment = (i == count - 1);

    // --- Interpolation Strategy Selection ---
    // Branch-cut guard: the planar projection is singular at the basis antipode,
    // so a segment with an endpoint there would project to a garbage azimuth and
    // interpolate wildly. Fall back to a (well-defined) geodesic edge for that
    // segment — near the antipode the "straight in projection" intent is
    // meaningless anyway. Two dot products on the planar path only; no cost to
    // the geodesic callers. Decided before the cull so the row-span bound below
    // matches the arc shape that will actually be rendered.
    const Vector *pc = planar_basis ? &planar_basis->v : nullptr;
    const bool antipodal_seam =
        pc && (dot(curr.pos, *pc) < -COS_PLANAR_ANTIPODE ||
               dot(next.pos, *pc) < -COS_PLANAR_ANTIPODE);
    const bool use_planar = planar_basis && !antipodal_seam;

    // Advance the rendered-arc accumulator for EVERY segment (drawn or culled) so
    // v0/v1 stay a true full-curve parameterization; seg_base snapshots the start
    // for the draw lambda. Skipped for geodesic polylines.
    if (override_uv) {
      seg_base = cumul;
      cumul += seg_arc_cache[i];
    }

    // Tier 3: Segment culling — skip if the edge's full screen-row span (arc
    // latitude bulge included) lies outside the clip band.
    if (clip_active) {
      float row_lo, row_hi;
      edge_row_span<W, H>(curr.pos, next.pos, use_planar ? planar_basis : nullptr,
                          row_lo, row_hi);
      if (!cr.could_intersect_y(row_lo, row_hi))
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

/**
 * @brief Per-primitive geometry/rasterization options for draw_fragments.
 *
 * Groups the three settings that vary per primitive into a named aggregate so
 * call sites read as designated initializers (`{.capacity = …, .close_loop =
 * …}`) instead of an unlabeled `size_t, bool, const Basis*` trio. `close_loop`
 * and `planar_basis` default to the common (geodesic, open) case, so most
 * primitives only spell out `.capacity`.
 */
struct FragmentDrawParams {
  size_t capacity;            /**< Fragment buffer reservation (per-primitive). */
  bool close_loop = false;    /**< Passed to rasterize (closes last→first edge). */
  const Basis *planar_basis = nullptr; /**< Planar projection basis (null = geodesic). */
};

/**
 * @brief Run the shared per-primitive draw ritual.
 *
 * Every Plot primitive opens a ScratchScope, binds a Fragments buffer, fills it,
 * applies the optional vertex shader, and rasterizes. The ScratchScope must
 * outlive the rasterize call (the arena backs the fragments), so the whole
 * sequence lives in one helper rather than being split.
 *
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam FillFn Callable (Fragments &) -> void supplying the primitive's
 *                sampling; inlines at -O3 for a zero-cost per-primitive copy.
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
                  params.planar_basis);
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
    if (std::abs(angle) < TOLERANCE) {
      points.push_back(f1);
      points.push_back(f1); // draw at least a dot
      return;
    }

    Vector axis;
    if (std::abs(PI_F - angle) < TOLERANCE) {
      // Antipodal endpoints: cross() degenerates to ~0 (infinitely many
      // geodesics connect them), so normalizing it yields a garbage axis. Pick
      // a stable perpendicular axis instead, matching the geodesic strategy.
      axis = (std::abs(dot(f1.pos, X_AXIS)) > math::COS_AXIS_PARALLEL)
                 ? cross(f1.pos, Y_AXIS).normalized()
                 : cross(f1.pos, X_AXIS).normalized();
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
 * @brief Draws vertices without connection.
 */
struct Vertices {
  /**
   * @brief Draws a list of vertices as points.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param points List of points.
   * @param fragment_shader Shader function; must be non-null (the per-vertex
   *        call below is unguarded and operator()'s assert is NDEBUG-stripped).
   */
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &points,
                   FragmentShaderFn fragment_shader) {
    // Cold (once per call), not per-pixel: trap a null shader before the loop.
    HS_CHECK(fragment_shader,
             "Vertices::draw requires a non-null fragment_shader");
    for (const Fragment &p : points) {
      Fragment f = p;
      f.color = Color4(0, 0, 0, 0);

      fragment_shader(f.pos, f);
      pipeline.plot(canvas, f.pos, f.color.color, f.age, f.color.alpha);
    }
  }
};

/**
 * @brief Multiline primitive (Polyline).
 * Registers:
 *  v0: Path Progress (0.0 -> 1.0)
 *  v1: Cumulative Arc Length (radians) — geodesic chord-polygon length
 *  v2: Vertex Index
 * @note v0/v1 accumulate the GEODESIC (great-circle) distance between
 *       consecutive control points — the chord-polygon parameterization, which is
 *       the rendered arc under Multiline's geodesic edges. If a `planar_basis` is
 *       passed to the draw call, the rasterizer re-derives v0/v1 from the longer
 *       azimuthal-equidistant arc it actually draws, so the registers track the
 *       rendered position in either mode.
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
      // Avoid divide-by-zero on a degenerate (sub-epsilon) path. v0 then
      // collapses toward 0 for every vertex (current_len << total_len==1); the
      // path is geometrically a point, so the lost progress parameter is moot.
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
    // We manually close the loop in sample() if requested, so pass false to
    // rasterize (close_loop).
    draw_fragments<W, H>(
        pipeline, canvas, vertex_shader, fragment_shader,
        {.capacity = vertices.size() + 2},
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
 * position with v0 = 1, the arc length continued across the wrap edge, and
 * v2 = num_verts, so a `close_loop` rasterize draws the final edge back to the
 * start without a UV seam. This is the shared skeleton for the accumulated-arc
 * closed rings (Star, Flower, DistortedRing); Ring itself uses an analytic arc
 * length and keeps its own loop. For the PLANAR callers (Star, Flower) the
 * rasterizer overrides v0/v1 with the true rendered azimuthal arc, so these
 * stored geodesic values then seed only the optional vertex shader.
 */
template <typename PosFn>
inline void sample_closed_ring(Fragments &points, int num_verts, PosFn pos_fn) {
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
   * analytic arc length (theta·sin(radius)).
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
    const float arc_scale = sinf(work_radius);

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
      f.v1 = theta * arc_scale;
      f.v2 = static_cast<float>(i);
      f.age = 0;

      points.push_back(f);
    }

    // Manual Close (Overlap): the close vertex sits at theta == 2π, whose
    // position is identical to the i == 0 vertex (cos/sin(2π+φ) == cos/sin(φ));
    // reusing it also avoids the float error a literal 2π+φ argument introduces.
    if (num_samples > 0) {
      Fragment f;
      f.pos = first_pos;
      f.v0 = 1.0f;
      f.v1 = 2.0f * PI_F * arc_scale;
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
   * @details Ring::draw's only sampling configuration is num_samples == W, whose
   * angle grid (i*2π/W) is exactly TrigLUT<W,H>::cos_theta/sin_theta, so the
   * per-sample libm cosf(θ+φ)/sinf(θ+φ) becomes the precomputed θ-grid plus one
   * angle-addition against cos/sin(φ) — saving ~2*(W+1) libm trig calls per ring
   * per frame. Keeps Ring's analytic arc length (θ*arc_scale) and its own
   * overlap close; see sample_closed_ring's note on why Ring is not folded into
   * that helper. The runtime int-num_samples overload above stays for the
   * polygon samplers, whose vertex counts (num_sides, W/4, …) do not match the
   * LUT grid.
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
    const float arc_scale = sinf(work_radius);

    const float step = 2.0f * PI_F / W;

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    const float cos_phase = cosf(phase);
    const float sin_phase = sinf(phase);

    for (int i = 0; i < W; i++) {
      // Angle-addition identity: cos/sin(θ+φ) from the precomputed θ-grid.
      float cos_t = TrigLUT<W, H>::cos_theta(i) * cos_phase -
                    TrigLUT<W, H>::sin_theta[i] * sin_phase;
      float sin_t = TrigLUT<W, H>::sin_theta[i] * cos_phase +
                    TrigLUT<W, H>::cos_theta(i) * sin_phase;
      Vector u_temp = (u * cos_t) + (w * sin_t);

      Fragment f;
      f.pos = ((v * d_val) + (u_temp * r_val)).normalized();
      f.v0 = static_cast<float>(i) / W;
      f.v1 = (i * step) * arc_scale;
      f.v2 = static_cast<float>(i);
      f.age = 0;

      points.push_back(f);
    }

    // Manual Close (Overlap): θ = 2π folds to (cos φ, sin φ) by periodicity.
    Fragment f;
    Vector u_temp = (u * cos_phase) + (w * sin_phase);
    f.pos = ((v * d_val) + (u_temp * r_val)).normalized();
    f.v0 = 1.0f;
    f.v1 = (2.0f * PI_F) * arc_scale;
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
 * @note This shape always renders with PLANAR (azimuthal-equidistant) edges, which
 *       bow LONGER than the great-circle chord between vertices. The rasterizer
 *       re-derives v0/v1 from that true rendered arc (v0 the [0,1] fraction of the
 *       planar perimeter, v1 the cumulative planar arc in radians), so both track
 *       the drawn position rather than the shorter chord polygon.
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
    // else the azimuthal projection bows the polygon edges instead of keeping
    // them straight. radius <= 1 keeps the supplied chart unchanged.
    Basis planar_basis = basis;
    if (radius > 1.0f) {
      Vector v = -basis.v;
      Vector ref = (std::abs(dot(v, X_AXIS)) > math::COS_AXIS_PARALLEL) ? Y_AXIS
                                                                        : X_AXIS;
      Vector u = cross(v, ref).normalized();
      Vector w = cross(v, u).normalized();
      planar_basis = {u, v, w};
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
      // Angle-addition identity: cos/sin(θ+φ) from precomputed LUT
      float cos_t = TrigLUT<W, H>::cos_theta(i) * cos_phase -
                    TrigLUT<W, H>::sin_theta[i] * sin_phase;
      float sin_t = TrigLUT<W, H>::sin_theta[i] * cos_phase +
                    TrigLUT<W, H>::cos_theta(i) * sin_phase;
      Vector u_temp = (u * cos_t) + (w * sin_t);

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
    HS_CHECK(n >= 1); // matches the count guard the SplineChain samplers carry
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
 * @note This shape always renders with PLANAR (azimuthal-equidistant) edges, which
 *       bow LONGER than the great-circle chord between vertices. The rasterizer
 *       re-derives v0/v1 from that true rendered arc (v0 the [0,1] fraction of the
 *       planar perimeter, v1 the cumulative planar arc in radians), so both track
 *       the drawn position rather than the shorter chord polygon.
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
    Vector center = basis.v;
    if (radius > 1.0f)
      center = -center;

    Vector v = center;
    Vector ref = (std::abs(dot(v, X_AXIS)) > math::COS_AXIS_PARALLEL) ? Y_AXIS
                                                                      : X_AXIS;
    Vector u = cross(v, ref).normalized();
    Vector w = cross(v, u).normalized();
    Basis planar_basis = {u, v, w};

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
 * @note This shape always renders with PLANAR (azimuthal-equidistant) edges, which
 *       bow LONGER than the great-circle chord between vertices. The rasterizer
 *       re-derives v0/v1 from that true rendered arc (v0 the [0,1] fraction of the
 *       planar perimeter, v1 the cumulative planar arc in radians), so both track
 *       the drawn position rather than the shorter chord polygon.
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
    // The planar chart and sample()'s petal placement must derive from the same
    // antipode mapping, so both call get_antipode(basis, radius) on identical
    // inputs; keep them in lockstep if either mapping changes.
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    // Center the planar chart on work_basis.v, the pole *opposite* the petal ring
    // (which sits at colatitude apothem = PI - outer >= PI/2, near -work_basis.v):
    // the ring is a constant-radius polygon, and projecting it through the
    // far-pole chart (where R -> PI) bows the straight edges outward into petals.
    Vector center = work_basis.v;

    // Construct a planar basis aligned with the center
    Vector ref =
        (std::abs(dot(center, X_AXIS)) > math::COS_AXIS_PARALLEL) ? Y_AXIS
                                                                  : X_AXIS;
    Vector u_p = cross(center, ref).normalized();
    Vector w_p = cross(center, u_p).normalized();
    Basis planar_basis = {u_p, center, w_p};

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
   * @details A mesh exceeding this is a sizing bug (traps on the cold setup
   * path), not a recoverable case. Sized for a TriangularBitset of 128*127/2
   * bits = 1016 bytes.
   */
  static constexpr int kDedupCapacity = 128;

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
                                   TriangularBitset<kDedupCapacity> &visited,
                                   Fn &&fn) {
    visited.clear();

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    size_t offset = 0;

    for (size_t i = 0; i < num_f; ++i) {
      int count = fc[i];
      for (int k = 0; k < count; ++k) {
        int u = fi[offset + k];
        int v = fi[offset + (k + 1) % count];
        int small = std::min(u, v);
        int large = std::max(u, v);

        // A vertex index past the dedup bitset's capacity is a mesh-sizing bug
        // with no valid recovery: silently dropping the edge would mask it and
        // leave a wireframe with missing lines. Trap on the cold setup path
        // (platform.h).
        HS_CHECK(large < kDedupCapacity);

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

    // O(1) edge dedup in a 1016-byte triangular bit matrix over 128 vertices,
    // arena-allocated rather than stack-resident: this is on the deep mesh
    // render chain and Phantasm's DTCM stack is tight. Held in scratch_arena_b
    // for the whole call; the per-edge scratch_arena_a scopes below are
    // independent, so the heavily-used render scratch keeps its full headroom.
    ScratchScope visited_guard(scratch_arena_b);
    auto &visited = *new (scratch_arena_b.allocate(
        sizeof(TriangularBitset<kDedupCapacity>),
        alignof(TriangularBitset<kDedupCapacity>))) TriangularBitset<kDedupCapacity>();

    for_each_unique_edge(mesh, visited, [&](int u, int v) {
      // mesh.vertices[] only asserts in bounds (stripped under NDEBUG on the
      // device), so a malformed face index would read OOB silently on hardware.
      // This is the per-edge setup boundary, not a per-pixel path, so an
      // always-on HS_CHECK is contract-appropriate (platform.h). Checking only
      // the larger index covers both endpoints: u and v are read from uint16_t
      // face data, so they are non-negative — max(u,v) in bounds implies both
      // vertices[u] and vertices[v] are valid.
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
    // Dedup bitset (1016 B) in the arena, not on the stack: this runs at setup
    // on the deep mesh-build chain. The output `edges` lives in a separate arena
    // (persistent), so scratch_arena_b here cannot disturb it.
    ScratchScope visited_guard(scratch_arena_b);
    auto &visited = *new (scratch_arena_b.allocate(
        sizeof(TriangularBitset<kDedupCapacity>),
        alignof(TriangularBitset<kDedupCapacity>))) TriangularBitset<kDedupCapacity>();

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
 * @brief Particle System trails.
 * Registers:
 *  v0: Trail Progress (0.0=Head -> 1.0=Tail)
 *  v1: Trail Arc Length
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
   * @param vertex_shader Optional vertex shader.
   */
  template <int W, int H, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    int count = system.active_count;

    for (int i = 0; i < count; ++i) {
      const auto &p = system.pool[i];
      ScratchScope trail_guard(scratch_arena_a);
      Fragments trail;
      trail.bind(scratch_arena_a, 64);
      float cumulative_len = 0.0f;
      Vector last_pos;
      bool first = true;

      tween(p.history, [&](const Vector &v, float t) {
        if (!first) {
          cumulative_len += angle_between(last_pos, v);
        }
        last_pos = v;
        first = false;

        Fragment f;
        f.pos = v;
        f.v0 = t;
        f.v1 = cumulative_len;
        f.v2 = static_cast<float>(i);
        f.v3 = static_cast<float>(p.life) / static_cast<float>(system.max_life);
        f.age = 0;
        f.color = Color4(0, 0, 0, 0);
        trail.push_back(f);
      });

      if (!trail.is_empty()) {
        apply_vertex_shader(vertex_shader, trail);
        rasterize<W, H>(pipeline, canvas, trail, fragment_shader, false,
                        nullptr);
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

/**
 * @brief Single cubic Bézier curve on the sphere.
 * Registers:
 *  v0: curve progress (0→1)
 *  v1: cumulative arc length (radians)
 *  v2: 0 (single segment)
 */
struct Bezier {
  /**
   * @brief Samples num_samples+1 points along a cubic Bézier.
   * @param points Output fragment list; num_samples+1 fragments are appended.
   * @param p0 First control point.
   * @param p1 Second control point.
   * @param p2 Third control point.
   * @param p3 Fourth control point.
   * @param num_samples Sub-segment count (>0); t = i / num_samples.
   * @param mode Spline interpolation mode (geodesic by default).
   * @details v1 accumulates great-circle arc length (radians) between samples.
   */
  static void sample(Fragments &points, const Vector &p0, const Vector &p1,
                     const Vector &p2, const Vector &p3, int num_samples,
                     SplineMode mode = SplineMode::Geodesic) {
    HS_CHECK(num_samples > 0); // t = i / num_samples divides by the sample count
    float cumulative_len = 0.0f;
    Vector last_pos = p0;

    for (int i = 0; i <= num_samples; ++i) {
      float t = static_cast<float>(i) / num_samples;
      Vector pos = Spline::cubic(p0, p1, p2, p3, t, mode);

      if (i > 0)
        cumulative_len += angle_between(last_pos, pos);
      last_pos = pos;

      Fragment f;
      f.pos = pos;
      f.v0 = t;
      f.v1 = cumulative_len;
      f.v2 = 0.0f;
      f.age = 0;
      points.push_back(f);
    }
  }

  /**
   * @brief Draws a cubic Bézier from raw control vectors p0..p3.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param p0 First control point.
   * @param p1 Second control point.
   * @param p2 Third control point.
   * @param p3 Fourth control point.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param num_samples Sub-segment count (>0).
   * @param mode Spline interpolation mode (geodesic by default).
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &p0,
                   const Vector &p1, const Vector &p2, const Vector &p3,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader = {}, int num_samples = 32,
                   SplineMode mode = SplineMode::Geodesic) {
    draw_fragments<W, H>(pipeline, canvas, vertex_shader, fragment_shader,
                         {.capacity = num_samples + 1},
                         [&](Fragments &points) {
                           sample(points, p0, p1, p2, p3, num_samples, mode);
                         });
  }

  /**
   * @brief Convenience overload taking Fragment references (uses .pos).
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param f0 First control fragment.
   * @param f1 Second control fragment.
   * @param f2 Third control fragment.
   * @param f3 Fourth control fragment.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param num_samples Sub-segment count (>0).
   * @param mode Spline interpolation mode (geodesic by default).
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Fragment &f0,
                   const Fragment &f1, const Fragment &f2, const Fragment &f3,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader = {}, int num_samples = 32,
                   SplineMode mode = SplineMode::Geodesic) {
    draw<W, H>(pipeline, canvas, f0.pos, f1.pos, f2.pos, f3.pos,
               fragment_shader, vertex_shader, num_samples, mode);
  }
};

/**
 * @brief Catmull-Rom spline chain on the sphere.
 * Registers:
 *  v0: global progress (0→1 across full chain)
 *  v1: cumulative arc length (radians)
 *  v2: segment index
 */
struct SplineChain {
  /**
   * @brief Draws a Catmull-Rom spline through `control_points`.
   * @tparam W,H Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param control_points Spline control points (>=2).
   * @param tension Tangent shaping factor for the Catmull-Rom spans.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param closed If true, wraps the chain into a loop.
   * @param samples_per_segment Sub-samples per span (>0).
   * @param mode Spline interpolation mode (geodesic by default).
   * @details Converts each span to a cubic Bézier via Catmull-Rom tangents and
   * samples it; endpoints of a closed/open chain reuse phantom/clamped neighbors
   * for the end tangents.
   */
  template <int W, int H>
  static void
  draw(PipelineRef pipeline, Canvas &canvas, const Fragments &control_points,
       float tension, FragmentShaderFn fragment_shader,
       VertexShaderRef vertex_shader = {}, bool closed = false,
       int samples_per_segment = 16, SplineMode mode = SplineMode::Geodesic) {
    size_t n = control_points.size();
    if (n < 2)
      return;
    HS_CHECK(samples_per_segment > 0); // local_t = j / samples_per_segment

    // Segment 0 emits samples_per_segment+1 points (j=0..S); each later segment
    // emits S (j=1..S, sharing the previous segment's endpoint). Total is
    // seg_count*S + 1 for BOTH closed and open.
    const size_t frag_count =
        (closed ? n : (n - 1)) * samples_per_segment + 1;

    draw_fragments<W, H>(
        pipeline, canvas, vertex_shader, fragment_shader, {.capacity = frag_count},
        [&](Fragments &points) {
          float cumulative_len = 0.0f;
          Vector last_pos;
          bool first_point = true;

          size_t seg_count = closed ? n : n - 1;
          for (size_t i = 0; i < seg_count; ++i) {
            size_t i0 = (i == 0 && !closed) ? 0 : (i - 1 + n) % n;
            size_t i1 = i;
            size_t i2 = (i + 1) % n;
            size_t i3 = (i + 2 >= n && !closed) ? n - 1 : (i + 2) % n;

            Vector cp1, cp2;
            Spline::catmull_rom_tangents(
                control_points[i0].pos, control_points[i1].pos,
                control_points[i2].pos, control_points[i3].pos, tension, cp1,
                cp2);

            int start_j = (i == 0) ? 0 : 1;
            for (int j = start_j; j <= samples_per_segment; ++j) {
              float local_t = static_cast<float>(j) / samples_per_segment;
              Vector pos = Spline::cubic(control_points[i1].pos, cp1, cp2,
                                         control_points[i2].pos, local_t, mode);

              if (!first_point)
                cumulative_len += angle_between(last_pos, pos);
              last_pos = pos;
              first_point = false;

              Fragment f = Fragment::lerp(control_points[i1],
                                          control_points[i2], local_t);
              f.pos = pos;
              f.v0 = (static_cast<float>(i) + local_t) / seg_count;
              f.v1 = cumulative_len;
              f.v2 = static_cast<float>(i);
              points.push_back(f);
            }
          }
        });
  }
};

} // namespace Plot
