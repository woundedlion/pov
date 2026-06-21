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
 * @brief Floor on the sin(φ) step-density scale near the poles.
 * @details Caps sub-steps per segment so polar curves don't oversample. A
 * clamp, not a tolerance.
 */
static constexpr float MIN_POLE_SCALE = 0.05f;

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
 * @brief Planar interpolation strategy: builds an arc-uniform map for one edge.
 * @tparam ProcessSegmentFn Callable (map, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param planar_basis Azimuthal-equidistant projection basis.
 * @param isLastSegment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length map, endpoints, on-sphere
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
  constexpr int LEN_SAMPLES = 4;
  std::array<float, LEN_SAMPLES + 1> arc_cumul;
  arc_cumul[0] = 0.0f;
  Vector len_prev = unproject(0.0f);
  for (int k = 1; k <= LEN_SAMPLES; ++k) {
    Vector len_cur = unproject(static_cast<float>(k) / LEN_SAMPLES);
    arc_cumul[k] = arc_cumul[k - 1] + angle_between(len_prev, len_cur);
    len_prev = len_cur;
  }
  const float dist = arc_cumul[LEN_SAMPLES];

  // Arc-length parameterization: s is the arc fraction in [0,1]. Invert the
  // piecewise-linear cumulative-arc table to a projection parameter, then
  // unproject. A short scan over LEN_SAMPLES floats — no trig.
  auto map_planar = [=](float s) -> Vector {
    if (dist < math::EPS_GEOMETRIC)
      return unproject(s);
    float target = s * dist;
    int k = 0;
    while (k < LEN_SAMPLES - 1 && arc_cumul[k + 1] < target)
      ++k;
    float seg = arc_cumul[k + 1] - arc_cumul[k];
    float frac =
        (seg > math::EPS_GEOMETRIC) ? (target - arc_cumul[k]) / seg : 0.0f;
    float p = (static_cast<float>(k) + frac) / LEN_SAMPLES;
    return unproject(std::min(1.0f, std::max(0.0f, p)));
  };

  process_segment(map_planar, curr, next, dist, isLastSegment);
}

/**
 * @brief Geodesic interpolation strategy: builds a great-circle map for one edge.
 * @tparam ProcessSegmentFn Callable (map, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param isLastSegment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length map, endpoints, on-sphere
 *                        length (radians), and the last-segment flag.
 * @details Picks a stable perpendicular axis for antipodal/degenerate endpoints
 * and slerps along the great circle; a coincident-endpoint edge collapses to a
 * constant map.
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
    auto map_degenerate = [=](float) { return v1; };
    process_segment(map_degenerate, curr, next, total_dist, isLastSegment);
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

    auto map_geodesic = [=](float t) {
      float ang = total_dist * t;
      float s = fast_sinf(ang);
      float c = fast_cosf(ang);
      return (v1 * c) + (v_perp * s);
    };
    process_segment(map_geodesic, curr, next, total_dist, isLastSegment);
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
    // so there is no closed-form latitude extremum. Sample the arc with the same
    // unprojection the rasterizer uses (so the cull and the renderer agree to
    // the bit), then widen by the arc's Lipschitz bound. Working in row space
    // keeps the margin uniform (it would blow up near the poles in y): phi is
    // 1-Lipschitz in angular distance, so between samples |Δrow| ≤
    // (Δarc)·(H_VIRT−1)/π, and the projected chord length over-estimates the
    // on-sphere arc (the projection stretches tangential distance). That margin
    // makes the span provably gap-free without dense sampling.
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
 * @brief Adaptively rasterize a fragment polyline onto the sphere.
 *
 * Walks consecutive fragment pairs, picks a geodesic or planar interpolation
 * strategy per segment, sub-steps each segment at ≈one-pixel-column density
 * (with a sinφ floor near the poles), and plots through the pipeline. Segments
 * whose full screen-row span lies outside the active clip band are culled.
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
  double _plot_t0 = emscripten_get_now();
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
  ScratchScope _sc(scratch_arena_a);
  ArenaVector<float> _steps_cache;
  // The cache holds ONE segment's adaptive sub-steps (cleared per segment). Each
  // step advances ≈ one pixel column, so a single segment needs ≲ W steps (a
  // great-circle arc crosses each longitude ~once; the sinφ≥0.05 clamp caps the
  // pole case). Size off W with 2× headroom; the simulation loop breaks at
  // capacity as a backstop.
  size_t max_cache = std::max((size_t)64, (size_t)(2 * W));
  _steps_cache.bind(scratch_arena_a, max_cache);

  // Adaptively sub-step and plot one segment. `map(t)` returns the sphere point
  // at arc fraction t in [0,1] under the chosen strategy; `total_dist` is the
  // segment's on-sphere length (radians). Endpoints are omitted on interior /
  // closed segments so a shared vertex isn't plotted twice.
  auto process_segment = [&](auto &&map, const Fragment &curr,
                             const Fragment &next, float total_dist,
                             bool isLastSegment) {
    // The degenerate and fast paths below plot curr.pos/next.pos directly, without
    // the renormalize the DRAWING PHASE applies further down: these are the
    // original sampled vertices (already unit), not the fast_sinf/cosf-interpolated
    // map() outputs that drift ~0.04% off the unit sphere. Precondition: callers
    // pass unit fragment positions.
    // Degenerate (coincident endpoints): plot at most a single dot.
    if (total_dist < math::EPS_GEOMETRIC) {
      bool shouldOmit = (close_loop) ? true : !isLastSegment;
      if (!shouldOmit) {
        Fragment f_copy = curr;
        f_copy.color = Color4(0, 0, 0, 0);

        fragment_shader(curr.pos, f_copy);
        pipeline.plot(canvas, curr.pos, f_copy.color.color, f_copy.age,
                      f_copy.color.alpha);
      }
      return;
    }

    // FAST PATH: sub-pixel segment — skip simulation, just plot start
    const float base_step = (2.0f * PI_F) / W;
    if (total_dist < base_step) {
      Fragment f = curr;
      f.color = Color4(0, 0, 0, 0);
      fragment_shader(curr.pos, f);
      pipeline.plot(canvas, curr.pos, f.color.color, f.age, f.color.alpha);
      if (!close_loop && isLastSegment) {
        Fragment fl = next;
        fl.color = Color4(0, 0, 0, 0);
        fragment_shader(next.pos, fl);
        pipeline.plot(canvas, next.pos, fl.color.color, fl.age,
                      fl.color.alpha);
      }
      return;
    }

    // 1. SIMULATION PHASE
    _steps_cache.clear();
    float sim_dist = 0.0f;
    Vector p_temp = map(0.0f);

    while (sim_dist < total_dist) {
      float scale_factor =
          std::max(MIN_POLE_SCALE, sqrtf(std::max(0.0f, 1.0f - p_temp.y * p_temp.y)));
      float step = base_step * scale_factor;

      // Backstop: a pathological segment (e.g. a huge-radius shape wrapping the
      // sphere) could still exceed the 2*W cache. Stop subdividing and let the
      // normalized replay stretch the cached steps over the rest of the segment
      // — coarser sampling on an extreme arc is fine (and far better than
      // trapping a live show). Normal segments never reach this.
      if (_steps_cache.size() >= _steps_cache.capacity())
        break;
      _steps_cache.push_back(step);
      sim_dist += step;

      if (sim_dist < total_dist) {
        p_temp = map(sim_dist / total_dist);
      }
    }

    // sim_dist > 0 always holds here (the loop runs at least once past the
    // total_dist >= base_step fast-path guard above), so the fallback is dead;
    // use 1.0f (identity, no stretch) rather than 0.0f, which would have zeroed
    // every replayed step.
    float scale = (sim_dist > 0.0f) ? (total_dist / sim_dist) : 1.0f;
    bool omitLast = (close_loop) ? true : !isLastSegment;

    if (omitLast && _steps_cache.is_empty())
      return;

    // 2. DRAWING PHASE
    //
    // map_geodesic/map_planar build interpolated points with fast_sinf/fast_cosf,
    // which don't satisfy sin²+cos²=1 — so the result is ~0.04% non-unit. The
    // sink's vector_to_pixel skips normalization and takes phi = acos(v.y)
    // directly; near the pole acos has infinite slope, so a y of 0.9996 instead
    // of 1.0 lands ~1.3 rows below the pole. Re-normalize the interpolated
    // positions here (the sampled vertices are already unit) so polar samples
    // map to the correct row. Plot-path only; one sqrt per drawn sample.
    {
      Vector start_pos = map(0.0f).normalized();
      Fragment f = Fragment::lerp(curr, next, 0.0f);
      f.pos = start_pos;
      f.color = Color4(0, 0, 0, 0);

      fragment_shader(start_pos, f);
      pipeline.plot(canvas, start_pos, f.color.color, f.age, f.color.alpha);
    }

    size_t loop_limit =
        omitLast ? _steps_cache.size() - 1 : _steps_cache.size();
    float current_dist = 0.0f;

    for (size_t j = 0; j < loop_limit; j++) {
      float step = _steps_cache[j] * scale;
      current_dist += step;

      float t = (total_dist > 0.0f) ? (current_dist / total_dist) : 1.0f;

      // `t` (hence the drawn POSITION) is parameterized by the RENDERED arc
      // length — geodesic or planar, whichever this segment uses. The fragment
      // registers, however, are linearly interpolated from the control points'
      // values; for v1 (cumulative arc length) those were computed from geodesic
      // chords, so under planar interpolation the interpolated v1 follows the
      // chord polygon, not the longer rendered planar arc.
      Vector p = map(t).normalized();
      Fragment f = Fragment::lerp(curr, next, t);
      f.pos = p;
      f.color = Color4(0, 0, 0, 0);

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
  canvas.add_render_us(emscripten_get_now() - _plot_t0);
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
  ScratchScope _frag(scratch_arena_a);
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
      points.push_back(f1); // Draw at least a dot
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
      f.v2 = 0.0f; // Standard: Single lines have v2=0
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
 *  v1: Cumulative Arc Length (radians) — GEODESIC chord-polygon length
 *  v2: Vertex Index
 * @note v0/v1 accumulate the GEODESIC (great-circle) distance between
 *       consecutive control points — the chord-polygon parameterization. That
 *       equals the rendered arc length only under geodesic edge interpolation.
 *       Under PLANAR interpolation (a `planar_basis` passed to the draw call)
 *       the rendered edge is an azimuthal-equidistant straight line that is
 *       LONGER than the great-circle chord, so v1 under-counts the actual drawn
 *       arc length. A shader that needs an exact planar-arc parameterization must
 *       not rely on v1 in planar mode.
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

    // 1. Calculate total length to normalize v0
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

    // 2. Generate fragments
    float current_len = 0.0f;
    it = std::begin(vertices);
    prev = *it;

    // Add first point
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
 * length and keeps its own loop.
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
    // position is identical to the i == 0 vertex (cos/sin(2π+φ) == cos/sin(φ)).
    // Reuse it instead of recomputing the trig + normalize — and it avoids the
    // tiny float error a literal 2π+φ argument would introduce.
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
      float cos_t = TrigLUT<W, H>::cos_theta[i] * cos_phase -
                    TrigLUT<W, H>::sin_theta[i] * sin_phase;
      float sin_t = TrigLUT<W, H>::sin_theta[i] * cos_phase +
                    TrigLUT<W, H>::cos_theta[i] * sin_phase;
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
    // W samples keep circles smooth and avoid pinching at the poles.
    draw_fragments<W, H>(
        pipeline, canvas, vertex_shader, fragment_shader,
        {.capacity = W + 2, .close_loop = true},
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
 *  v1: Arc Length (radians) — GEODESIC chord length
 *  v2: Vertex index
 * @note v1 accumulates the GEODESIC (great-circle) distance between vertices, but
 *       this shape always renders with PLANAR (azimuthal-equidistant) edges, which
 *       are LONGER than the great-circle chord — so v1 under-counts the actual drawn
 *       arc length. A shader needing an exact planar-arc parameterization must not
 *       rely on v1. (Same hazard Multiline documents.)
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
    size_t start_idx = points.size();
    Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);
    float cumul = 0.0f;
    for (size_t i = start_idx; i < points.size(); i++) {
      points[i].v1 = cumul;
      if (i < points.size() - 1) {
        cumul += angle_between(points[i].pos, points[i + 1].pos);
      }
    }
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
    // exactly as Plot::Star does; without the flip the azimuthal projection bows
    // the polygon edges (the deliberate far-field petal-bowing Flower exploits)
    // instead of keeping them straight. radius <= 1 keeps the supplied chart
    // unchanged.
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

    // Re-calculate v1 to be true geodesic chord length
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
    // Mirror sample() exactly (at phase 0, as the renderer is invoked) so the
    // returned point lands on the *drawn* ring: same theta_eq mapping and
    // shift_fn(theta/2PI) domain. Any divergence here would detach callers'
    // sampled points (e.g. Thrusters' thrust points) from the visible ring off
    // Radius=1.
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
      float cos_t = TrigLUT<W, H>::cos_theta[i] * cos_phase -
                    TrigLUT<W, H>::sin_theta[i] * sin_phase;
      float sin_t = TrigLUT<W, H>::sin_theta[i] * cos_phase +
                    TrigLUT<W, H>::cos_theta[i] * sin_phase;
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
 *  v1: Arc Length (radians) — GEODESIC chord length
 *  v2: Vertex index
 * @note v1 accumulates the GEODESIC (great-circle) distance between vertices, but
 *       this shape always renders with PLANAR (azimuthal-equidistant) edges, which
 *       are LONGER than the great-circle chord — so v1 under-counts the actual drawn
 *       arc length. A shader needing an exact planar-arc parameterization must not
 *       rely on v1. (Same hazard Multiline documents.)
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
 *  v1: Arc Length (radians) — GEODESIC chord length
 *  v2: Vertex index
 * @note v1 accumulates the GEODESIC (great-circle) distance between vertices, but
 *       this shape always renders with PLANAR (azimuthal-equidistant) edges, which
 *       are LONGER than the great-circle chord — so v1 under-counts the actual drawn
 *       arc length. A shader needing an exact planar-arc parameterization must not
 *       rely on v1. (Same hazard Multiline documents.)
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
      // Unproject Polar -> Sphere
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
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    // Center the planar chart on work_basis.v, the pole *opposite* the petal ring
    // (which sits at colatitude apothem = PI - outer >= PI/2, near -work_basis.v).
    // This is intentional and load-bearing for the look, NOT a bug: the ring is a
    // constant-radius polygon, so it only reads as a flower because projecting it
    // through the far-pole chart (where R -> PI) bows the straight edges outward
    // into petals. Centering on -work_basis.v (the ring's own pole, as SDF/Scan
    // Flower do for their sector-folded geometry) collapses the edges to clean
    // tangent-plane chords -- i.e. a plain planar polygon, not a flower.
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

    // O(1) edge dedup in 1016 bytes (triangular bit matrix over 128 vertices)
    TriangularBitset<kDedupCapacity> visited;
    visited.clear();

    auto process_edge = [&](int u, int v) {
      int small = std::min(u, v);
      int large = std::max(u, v);

      // A vertex index past the dedup bitset's capacity is a mesh-sizing bug
      // with no valid recovery: silently dropping the edge would mask it and
      // leave a wireframe with missing lines. Trap on the cold per-edge setup
      // path, consistent with the OOB guard below (platform.h).
      HS_CHECK(large < kDedupCapacity);

      if (visited.test_and_set(small, large))
        return;

      // mesh.vertices[] only asserts in bounds (stripped under NDEBUG on the
      // device), so a malformed face index would read OOB silently on hardware.
      // This is the per-edge setup boundary, not a per-pixel path, so an
      // always-on HS_CHECK is contract-appropriate (platform.h). Checking only
      // `large` covers both endpoints: u and v are read from uint16_t face data
      // (call site below), so they are non-negative — large = max(u,v) in bounds
      // implies 0 <= small <= large < size(), hence vertices[u] and vertices[v]
      // are both valid.
      HS_CHECK(static_cast<size_t>(large) < mesh.vertices.size());

      Fragment fu;
      fu.pos = mesh.vertices[u];
      Fragment fv;
      fv.pos = mesh.vertices[v];

      ScratchScope _edge(scratch_arena_a);
      Fragments points;
      points.bind(scratch_arena_a, 16);
      // Fixed 10-segment pre-sample: the adaptive rasterize() below further
      // sub-steps each chord at ~1px, so this only sets the per-vertex shader
      // granularity and the great-circle seed points. Left fixed (not arc-length
      // adaptive) deliberately — lowering it would perturb the adaptive
      // sub-step/normalization pattern and thus the rendered edges, a change
      // that needs sign-off rather than a polish pass.
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

      edge_index++;
    };

    size_t offset = 0;
    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();

    for (size_t i = 0; i < num_f; ++i) {
      int count = fc[i];
      for (int k = 0; k < count; ++k) {
        int u = fi[offset + k];
        int v = fi[offset + (k + 1) % count];
        process_edge(u, v);
      }
      offset += count;
    }
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
    TriangularBitset<kDedupCapacity> visited;
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
        // A vertex index past the dedup bitset's capacity is a mesh-sizing
        // bug with no valid recovery: silently dropping the edge would mask
        // it and leave a wireframe with missing lines. Trap on this cold
        // setup path, matching the dynamic draw() overload above.
        HS_CHECK(large < kDedupCapacity);
        if (!visited.test_and_set(small, large)) {
          edges.push_back({(uint16_t)u, (uint16_t)v});
        }
      }
      offset += count;
    }
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

      Fragment fu;
      fu.pos = mesh.vertices[edges[ei].u];
      Fragment fv;
      fv.pos = mesh.vertices[edges[ei].v];

      ScratchScope _edge(scratch_arena_a);
      Fragments points;
      points.bind(scratch_arena_a, 16);
      // Fixed pre-sample; the adaptive rasterize() below refines to ~1px (see
      // the matching note on the other edge-draw path). Left fixed deliberately.
      Line::sample(points, fu, fv, 10);

      if (vertex_shader) {
        for (auto &p : points) {
          p.v2 = static_cast<float>(ei);
          vertex_shader(p);
        }
      } else {
        for (auto &p : points) {
          p.v2 = static_cast<float>(ei);
        }
      }
      rasterize<W, H>(pipeline, canvas, points, fragment_shader, false, nullptr);
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
      ScratchScope _trail(scratch_arena_a);
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
