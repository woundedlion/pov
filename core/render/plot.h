/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <array>
#include "math/geometry.h"
#include "render/shading.h"
#include "color/color.h"
#include "engine/constants.h"
#include "render/clip.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include "mesh/triangular_bitset.h"
#include "render/plot_raster.h"

/**
 * @file plot.h
 * @brief The curve rasterizer: the Plot primitives that stroke lines, rings,
 * polygons, meshes and particle systems, over plot_raster.h and plot_cull.h.
 */

namespace Plot {

/**
 * @brief Inner/outer radius ratio for star shapes.
 */
inline constexpr float STAR_INNER_RATIO = ::STAR_INNER_RATIO;

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
    HS_CHECK(density >= 1);

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
      // Avoid divide-by-zero on a degenerate path.
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
  assert(i >= 0 && i < W);
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
   * @param radius Ring radius as a fraction of the hemisphere.
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
   * @param radius Ring radius as a fraction of the hemisphere.
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
   * @param radius Ring radius as a fraction of the hemisphere.
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
   * @param radius Ring radius as a fraction of the hemisphere.
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
    HS_CHECK(num_sides >= 1);
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
    HS_CHECK(num_sides >= 1);
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
    HS_CHECK(num_sides >= 1);
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
    HS_CHECK(num_sides >= 1);
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
   * @param cb Clip-band cut boundaries, resolved once by the calling draw().
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
                        int u, int v, int edge_index, const ClipCutBounds &cb,
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
      const int cuts = geodesic_clip_splits(fu.pos, fv.pos, es, cb, ts);
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
          bits[0] = RasterOptions::EDGE_VISIBLE;
        } else if (!gate_trail_edges<W, H>(pipeline, cr, xc, points, bits)) {
          return;
        }
        rasterize<W, H>(pipeline, canvas, points, fragment_shader,
                        {.edge_flags = bits});
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
    auto &visited = *scratch_arena_b.make<TriangularBitset<DEDUP_CAPACITY>>();

    const ClipRegion &cr = canvas.clip();
    const ClipCutBounds cb = make_clip_cut_bounds<W, H>(cr, cr.x_clip());

    for_each_unique_edge(mesh, visited, [&](int u, int v) {
      // mesh.vertices[] only asserts in bounds (stripped on device), so guard the
      // per-edge setup boundary here. u,v come from uint16_t face data (non-
      // negative), so max(u,v) in bounds implies both endpoints are valid.
      HS_CHECK(static_cast<size_t>(std::max(u, v)) < mesh.vertices.size());

      draw_edge<W, H>(pipeline, canvas, mesh, u, v, edge_index, cb,
                      fragment_shader, vertex_shader);

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
    auto &visited = *scratch_arena_b.make<TriangularBitset<DEDUP_CAPACITY>>();

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
    const ClipRegion &cr = canvas.clip();
    const ClipCutBounds cb = make_clip_cut_bounds<W, H>(cr, cr.x_clip());

    for (size_t ei = 0; ei < edges.size(); ++ei) {
      if (mask && !mask->owns(edges[ei].u, edges[ei].v))
        continue;
      // Setup-boundary OOB guard (see the face-walk overload above): the raw
      // edge list could outlive or mismatch its mesh, and mesh.vertices[] only
      // asserts (compiled out on device).
      HS_CHECK(edges[ei].u < mesh.vertices.size() &&
               edges[ei].v < mesh.vertices.size());

      draw_edge<W, H>(pipeline, canvas, mesh, edges[ei].u, edges[ei].v,
                      static_cast<int>(ei), cb, fragment_shader, vertex_shader);
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
  template <typename SystemT> static consteval int trail_sample_stride() {
    if constexpr (requires { SystemT::TRAIL_SAMPLE_STRIDE; })
      return SystemT::TRAIL_SAMPLE_STRIDE;
    return 1;
  }

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
      constexpr int sample_stride = trail_sample_stride<SystemT>();
      const bool append_live_tip = [&] {
        if constexpr (sample_stride == 1) {
          return false;
        } else {
          // A max_life lowered below a live particle's remaining life would
          // underflow this uint16_t subtraction.
          const uint16_t age =
              p.life < system.max_life ? system.max_life - p.life : 0;
          return (age - 1) % sample_stride != 0;
        }
      }();
      const size_t point_count = trail_len + (append_live_tip ? 1 : 0);
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
      if (p.life == 0 || point_count < 2)
        continue;
      const float v2 = [&] {
        if constexpr (std::is_same_v<ParticleV2Fn, std::nullptr_t>)
          return static_cast<float>(i);
        else
          return particle_v2(p, i);
      }();
      const float particle_life =
          std::min(static_cast<float>(p.life) * inv_max_life, 1.0f);
      ScratchScope trail_guard(scratch_arena_a);
      Fragments trail;
      trail.bind(scratch_arena_a, point_count);
      // Original (pre-shader) positions, kept for the deferred pass.
      ArenaVector<Vector> orig;
      if (has_deferred_shader)
        orig.bind(scratch_arena_a, point_count);
      {
        HS_PROFILE(plot_ps_tween);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
        hs::DwtStallBatch history_batch(
            hs::g_mindsplatter_stalls.history_vertex);
#endif
        auto emit = [&](const Vector &v, float t) {
          trail.emplace_back(
              Fragment{.pos = v, .v0 = t, .v2 = v2, .v3 = particle_life});
          if constexpr (FuseVertex)
            vertex_shader(trail.back());
          if (has_deferred_shader)
            orig.push_back(v);
#ifdef HS_PROFILE_MINDSPLATTER_STALLS
          history_batch.step();
#endif
        };
        if constexpr (sample_stride == 1) {
          p.history.tween(emit);
        } else {
          const float denominator = static_cast<float>(point_count - 1);
          for (size_t j = 0; j < trail_len; ++j)
            emit(p.history.get(j), static_cast<float>(j) / denominator);
          if (append_live_tip)
            emit(p.position, 1.0f);
        }
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
              bits[e] = RasterOptions::EDGE_CLASSIFIED |
                        RasterOptions::EDGE_ONE_DOT |
                        (v ? RasterOptions::EDGE_VISIBLE : uint8_t{0});
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
            bits[e] = RasterOptions::EDGE_CLASSIFIED |
                      (v ? RasterOptions::EDGE_VISIBLE : uint8_t{0});
            if (!v)
              HS_MSP_COUNT(edge_rejects);
            any = any || v;
            HS_MSP_STALL_STOP(trail_gate, edge_gate_start);
          }
        } else {
          for (size_t e = 0; e < edges; ++e) {
            bits[e] = edge_visible_in_clip<W, H>(pipeline, cr, xc, trail[e].pos,
                                                 trail[e + 1].pos, nullptr)
                          ? RasterOptions::EDGE_VISIBLE
                          : uint8_t{0};
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
        rasterize<W, H,
                  RasterConfig{
                      .single_pass = SinglePassRaster && sample_stride == 1,
                      .open_geodesic = true}>(pipeline, canvas, trail,
                                              fragment_shader,
                                              {.edge_flags = vis,
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
