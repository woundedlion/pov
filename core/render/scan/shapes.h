/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include "render/sdf.h"
#include "render/shading.h"
#include "color/color.h"
#include "render/filter/pipeline.h"
#include "containers/static_circular_buffer.h"
#include "render/canvas.h"
#include "platform/platform.h"
#include "render/scan/raster.h"

/**
 * @file shapes.h
 * @brief The SDF-backed scan primitives: rings, polygons, circles, points, stars and flowers.
 */

namespace Scan {

/**
 * @brief Draws a ring whose radius is modulated around the circumference by
 *        shift_fn.
 */
struct DistortedRing {
  /**
   * @brief Rasterizes an undisplaced ring with exact polar centerline distance.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   * @param suppress_pole_fill Drop the degenerate exact-pole row instead of
   *        full-row filling it (dense ring stacks; see get_horizontal_intervals).
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw_flat(PipelineRef pipeline, Canvas &canvas,
                        const Basis &basis, float radius, float thickness,
                        FragmentShaderFn fragment_shader, float phase = 0,
                        bool debug_bb = false,
                        bool suppress_pole_fill = false) {
    SDF::FlatDistortedRing shape(basis, radius, thickness, phase);
    shape.suppress_pole_fill = suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a circumference-modulated ring.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param shift_fn Scalar modulation function over the circumference.
   * @param amplitude Modulation amplitude in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness, ScalarFn shift_fn,
                   float amplitude, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    SDF::DistortedRing shape(basis, radius, thickness, shift_fn, amplitude,
                             phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a ring whose centerline is a shift-knot polyline with
   *        exact stroke distance (see SDF::DistortedRing's knot overload).
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param knots lut_n + 1 centerline shifts, entry lut_n repeating entry 0;
   *        must outlive the call.
   * @param lut_n Number of knot cells.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   * @param suppress_pole_fill Drop the degenerate exact-pole row instead of
   *        full-row filling it (dense ring stacks; see get_horizontal_intervals).
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness, const float *knots, int lut_n,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false, bool suppress_pole_fill = false) {
    SDF::KnotPrefilter prefilter;
    SDF::DistortedRing shape(basis, radius, thickness, knots, lut_n, phase,
                             prefilter);
    shape.suppress_pole_fill = suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

/**
 * @brief Fused single-pass rasterizer for a stack of same-axis distorted
 *        rings.
 */
HS_O3_BEGIN
struct DistortedRingStack {
  /** @brief Candidate-window pad and the even-spacing tolerance it affords. */
  static constexpr float WINDOW_PAD = 1e-3f;

  /**
   * @brief Traps unless every occupied ring sits on its even-spacing slot with
   *        zero phase and slot 0's basis.
   * @param n_rings Stack size.
   * @param shapes Ring shapes indexed by slot.
   * @param slot_by_ring n_rings entries mapping ring index -> slot, -1 if
   *        culled.
   * @param n_slots Number of shapes; every occupied entry must index below it.
   * @details A spacing deviation within WINDOW_PAD still lands inside the
   * shared candidate window, so geometry survives; beyond it the ring is
   * windowed out and silently disappears. The fused scan derives one
   * phase-free azimuth per pixel for the whole stack, so a non-zero ring phase
   * would index the wrong knot cell and mislabel v0. That frame is read from
   * slot 0 alone, so a divergent basis would silently render its ring at slot
   * 0's orientation.
   */
  HS_COLD_MEMBER
  static void check_stack_preconditions(int n_rings,
                                        const SDF::DistortedRing *shapes,
                                        const int8_t *slot_by_ring,
                                        int n_slots) {
    const float delta = PI_F / (n_rings + 1);
    for (int i = 0; i < n_rings; ++i) {
      const int s = slot_by_ring[i];
      if (s < 0)
        continue;
      HS_CHECK(s < n_slots, "ring stack slot index out of range");
      HS_CHECK(std::abs(shapes[s].target_angle - delta * (i + 1)) <=
               WINDOW_PAD);
      HS_CHECK(shapes[s].phase == 0.0f);
      HS_CHECK(dot(shapes[s].normal, shapes[0].normal) >= 1.0f - TOLERANCE);
      HS_CHECK(dot(shapes[s].u, shapes[0].u) >= 1.0f - TOLERANCE);
      HS_CHECK(dot(shapes[s].w, shapes[0].w) >= 1.0f - TOLERANCE);
    }
  }

  /**
   * @brief Rasterizes every ring of an evenly spaced same-axis stack in one
   *        scan over the union band.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader: shader(int slot, const Vector &p,
   *         Fragment &f), with f populated as by process_pixel (v0 = azimuth
   *         t, v1 = raw distance, v2 = stroke coverage).
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param n_rings Stack size: ring i's centerline colatitude about the shared
   *        axis must be PI * (i + 1) / (n_rings + 1).
   * @param shapes n_slots knot-mode rings sharing one Basis and zero phase, in
   *        ascending ring order; culled rings are simply absent.
   * @param slot_by_ring n_rings entries mapping ring index -> slot in shapes,
   *        -1 for culled rings.
   * @param n_slots Number of shapes; at least 1.
   * @param shader Per-ring fragment shader (see RingShaderT).
   * @details The per-pixel frame shared by every ring at a pixel (axis dot,
   * fast_acos, fast_atan2) is computed once, the candidate rings fall out of
   * the frame's polar angle by arithmetic (the stack is evenly spaced), and
   * each candidate runs its own cos reject + exact polyline distance via
   * SDF::DistortedRing::distance_from_frame. Candidates evaluate in ascending
   * ring index, so the pixels plotted, their per-pixel blend order and their
   * colors match rasterizing the rings one by one at pole_lod_aggressiveness 0;
   * only the redundant per-ring frame recompute is elided. Blend weights match
   * to float rounding, not bit-exactly: the shared arithmetic is spelled once
   * but inlined into two different loops, which -ffast-math reassociates
   * independently. This scan shades every column, so a non-zero aggressiveness,
   * which decimates the per-ring scan_region walk, breaks that equivalence. The
   * aliased exact-pole rows are dropped, matching the per-ring path under
   * suppress_pole_fill. Takes no debug flag and does not read canvas.debug():
   * unlike RingGroup there is no per-ring fallback, so the bounding-box tint
   * never reaches a stack.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  static void draw(PipelineT &pipeline, Canvas &canvas, int n_rings,
                   const SDF::DistortedRing *shapes, const int8_t *slot_by_ring,
                   int n_slots, RingShaderT &&shader) {
    HS_CHECK(canvas.width() == W && canvas.height() == H);
    check_pipeline_prepared(pipeline, canvas);
    HS_CHECK(n_slots >= 1);
    check_stack_preconditions(n_rings, shapes, slot_by_ring, n_slots);
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4;
    const float *sin_theta = TrigLUT<W, H>::sin_theta.data();

    // Union band and the global candidate half-width.
    int y_lo = H, y_hi = -1;
    float b_max = 0.0f;
    for (int s = 0; s < n_slots; ++s) {
      auto b = shapes[s].template get_vertical_bounds<H>();
      y_lo = std::min(y_lo, b.y_min);
      y_hi = std::max(y_hi, b.y_max);
      b_max = std::max(b_max, shapes[s].max_thickness);
    }
    const auto &cr = canvas.clip();
    const auto xc = cr.x_clip();
    y_lo = std::max(y_lo, cr.render_y_start());
    y_hi = std::min(y_hi, cr.render_y_end() - 1);

    // Window pad: fast_acos error plus float theta/index inversion slop; a
    // ring wrongly windowed in is discarded by its own exact cos reject.
    const float b_win = b_max + WINDOW_PAD;
    const float inv_delta = (n_rings + 1) / PI_F;

    // The per-ring path suppresses the aliased exact-pole rows
    // (suppress_pole_fill); its full-scan fallback for a near-canvas-pole
    // axis (r_val below the projection floor) scans every row.
    SDF::AxisProjection ap = SDF::project_axis(shapes[0].normal);
    const bool skip_pole_rows = ap.r_val >= SDF::MIN_HORIZONTAL_PROJ;

    const Vector axis_v = shapes[0].normal;
    const Vector axis_u = shapes[0].u;
    const Vector axis_w = shapes[0].w;

    SDF::DistanceResult res;
    Fragment frag;
    for (int y = y_lo; y <= y_hi; ++y) {
      const float sp = TrigLUT<W, H>::sin_phi[y];
      const float cp = TrigLUT<W, H>::cos_phi[y];
      if (skip_pole_rows && std::abs(ap.r_val * sp) < SDF::INTERVAL_DENOM_EPS)
        continue;
      walk_clip_columns<W>(xc, [&](int x) {
        Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
        const float d = dot(p, axis_v);
        const float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
        int ilo = static_cast<int>(ceilf((polar - b_win) * inv_delta)) - 1;
        int ihi = static_cast<int>(floorf((polar + b_win) * inv_delta)) - 1;
        if (ilo < 0)
          ilo = 0;
        if (ihi > n_rings - 1)
          ihi = n_rings - 1;
        if (ilo > ihi)
          return;
        const float dot_u = dot(p, axis_u);
        const float dot_w = dot(p, axis_w);
        float azimuth = fast_atan2(dot_w, dot_u);
        if (azimuth < 0)
          azimuth += 2 * PI_F;
        const float t_norm = wrap_t(azimuth / (2 * PI_F));
        const float sin_polar =
            sqrtf(std::max(1.0f - d * d, SDF::DistortedRing::POLE_SIN2_FLOOR));
        for (int i = ilo; i <= ihi; ++i) {
          const int s = slot_by_ring[i];
          if (s < 0)
            continue;
          shapes[s].distance_from_frame(d, polar, sin_polar, t_norm, res);
          const float dd = res.dist;
          if (dd >= 0.0f)
            continue;
          // process_pixel's stroke epilogue with a slot-aware shader.
          const float aa = res.size;
          const float alpha = aa > 0.0f ? quintic_kernel(-dd / aa) : 0.0f;
          if (alpha <= MIN_ALPHA)
            continue;
          frag.color = Color4(0, 0, 0, 0);
          frag.pos = p;
          frag.v0 = res.t;
          frag.v1 = res.raw_dist;
          frag.v2 = alpha;
          frag.v3 = res.aux;
          frag.size = res.size;
          frag.age = 0;
          shader(s, p, frag);
          if (frag.color.alpha > MIN_ALPHA)
            pipeline.plot(canvas, x, y, frag.color.color, frag.age,
                          frag.color.alpha * alpha);
        }
      });
    }
  }
};
HS_O3_END

/**
 * @brief Draws a flat (tangent-plane) regular polygon projected onto the
 *        sphere.
 */
struct PlanarPolygon {
  /**
   * @brief Rasterizes a tangent-plane regular polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon plane.
   * @param radius Polygon circumradius in world units.
   * @param sides Number of polygon sides.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::PlanarPolygon shape(res.first, res.second, sides, phase,
                             radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color tangent-plane regular polygon.
   */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::PlanarPolygon shape(res.first, res.second, sides, phase,
                             radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a great-circle line segment of given thickness between two
 *        vectors.
 */
struct Line {
  /**
   * @brief Rasterizes a great-circle segment between two unit vectors.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param v1 First endpoint as a world-space unit vector.
   * @param v2 Second endpoint as a world-space unit vector.
   * @param thickness Stroke thickness in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &v1,
                   const Vector &v2, float thickness,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    SDF::Line shape(v1, v2, thickness);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

/**
 * @brief Draws a solid ring using SDF rasterization.
 */
struct Ring {
  /**
   * @brief Draws a solid ring from an orientation basis.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Ring shape(res.first, res.second, thickness, phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Draws a solid ring from a plane-normal vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param normal Plane normal as a world-space vector.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &normal,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    Basis basis = make_basis(Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, thickness,
                           fragment_shader, phase, debug_bb);
  }
};

/**
 * @brief Fused single-pass rasterizer for a small group of rings.
 */
HS_O3_BEGIN
struct RingGroup {
  /**
   * @brief Rasterizes every ring of a group in one scan over the union band.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader: shader(int slot, const Vector &p,
   *         Fragment &f). One Fragment serves the whole scan and only color,
   *         pos, v2 (stroke coverage), size and age are refreshed per pixel —
   *         no UVs, no raw distance. v0, v1 and v3 hold their struct defaults
   *         until the shader itself writes them, and from then on whatever the
   *         previous invocation left, so a shader must not read a register it
   *         did not set.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes; at least 1.
   * @param shader Per-ring fragment shader (see RingShaderT).
   * @param debug_bb When true, falls back to per-ring rasterizes so the
   *        bounding-box tint keeps per-shape scan bounds; canvas.debug() takes
   *        the same fallback. Each ring is then scanned against its own row
   *        intervals rather than the covering ring's, so even a conforming
   *        shader renders the AA-tail difference described above, and the
   *        fallback fills v0/v1/v3 per pixel on top of that.
   * @details Row intervals come from one covering ring — member 0 inflated by
   * the group's maximum plane/radius deviation plus thickness — which contains
   * every member's band, so the per-row interval math runs once, not per
   * member. Per pixel the members evaluate in ascending slot order via the
   * inline stroke_alpha eval, so blend order matches rasterizing the rings
   * one by one. At pole_lod_aggressiveness 0 the only output divergence is
   * AA-tail pixels (alpha barely above the 0.001 cutoff) that a member's own
   * interval clip drops but the covering scan paints. This scan shades every
   * column, so a non-zero aggressiveness, which decimates the per-ring
   * scan_region walk, widens the divergence beyond that dust.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  static void draw(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes,
                   int n, RingShaderT &&shader, bool debug_bb = false) {
    static constexpr int MAX_RINGS = 8;
    HS_CHECK(canvas.width() == W && canvas.height() == H);
    check_pipeline_prepared(pipeline, canvas);
    HS_CHECK(n >= 1 && n <= MAX_RINGS);
    if (debug_bb || canvas.debug()) {
      for (int s = 0; s < n; ++s) {
        auto slot_shader = [&](const Vector &p, Fragment &f) {
          shader(s, p, f);
        };
        Scan::rasterize<W, H, false>(pipeline, canvas, shapes[s], slot_shader,
                                     true);
      }
      return;
    }

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();

    int sy_min[MAX_RINGS], sy_max[MAX_RINGS];
    int y_lo = H, y_hi = -1;
    for (int s = 0; s < n; ++s) {
      auto b = shapes[s].template get_vertical_bounds<H>();
      sy_min[s] = b.y_min;
      sy_max[s] = b.y_max;
      y_lo = std::min(y_lo, b.y_min);
      y_hi = std::max(y_hi, b.y_max);
    }
    const auto &cr = canvas.clip();
    y_lo = std::max(y_lo, cr.render_y_start());
    y_hi = std::min(y_hi, cr.render_y_end() - 1);
    if (y_lo > y_hi)
      return;

    // Covering ring: every point of member s's centerline lies within its
    // plane/radius deviation of the middle member's, so the middle member
    // inflated by the worst deviation plus that member's thickness contains
    // the whole group's stroke band — and, centered mid-group, its width is
    // the true union's plus epsilon. The 1e-3 absorbs fast_acos error.
    const int mid = n / 2;
    float pad_th = shapes[mid].thickness;
    for (int s = 0; s < n; ++s) {
      if (s == mid)
        continue;
      float dev = fast_acos(hs::clamp(dot(shapes[mid].normal, shapes[s].normal),
                                      -1.0f, 1.0f)) +
                  std::abs(shapes[s].target_angle - shapes[mid].target_angle) +
                  1e-3f;
      pad_th = std::max(pad_th, shapes[s].thickness + dev);
    }
    SDF::Ring cover(shapes[mid].basis, shapes[mid].radius, pad_th);

    Fragment frag;

    // Row-local walk instead of scan_region: the covering ring emits at most
    // 2 arcs per row, so small stack buffers replace the arena-backed CSG
    // interval machinery — and the walk compiles inside this O3 region, where
    // scan_region (-Os) forces the O3 pixel body out of line and calls it per
    // pixel. Runs come from the shared emit_row_runs, so the walked columns are
    // scan_region's.
    const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4;
    const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
    const auto xc = cr.x_clip();
    StaticCircularBuffer<SDF::Interval, 4> intervals;
    StaticCircularBuffer<SDF::Interval, 8> norm;
    static_assert(decltype(intervals)::CAPACITY >=
                      SDF::sdf_max_spans<SDF::Ring>::value,
                  "intervals must hold the covering Ring's per-row emission");
    static_assert(decltype(norm)::CAPACITY == 2 * decltype(intervals)::CAPACITY,
                  "norm must hold 2 spans per input interval (seam split)");
    int active[MAX_RINGS];
    int n_active = 0;

    auto pixel_run = [&](int x1, int x2, int y, float sp, float cp) {
      for (int x = x1; x < x2; ++x) {
        Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
        for (int i = 0; i < n_active; ++i) {
          const int s = active[i];
          const float alpha = shapes[s].stroke_alpha(dot(p, shapes[s].normal));
          if (alpha <= MIN_ALPHA)
            continue;
          frag.color = Color4(0, 0, 0, 0);
          frag.pos = p;
          frag.v2 = alpha;
          frag.size = shapes[s].thickness;
          frag.age = 0;
          shader(s, p, frag);
          if (frag.color.alpha > MIN_ALPHA)
            pipeline.plot(canvas, x, y, frag.color.color, frag.age,
                          frag.color.alpha * alpha);
        }
      }
    };
    for (int y = y_lo; y <= y_hi; ++y) {
      // Band gate hoisted per row: bounds pad (0.95 * thickness) is narrower
      // than the SDF's own band, so an out-of-band row can still eval alpha >
      // 0.001; gating keeps each slot's domain identical to its solo
      // rasterize. Ascending slot order preserves blend order.
      n_active = 0;
      for (int s = 0; s < n; ++s)
        if (y >= sy_min[s] && y <= sy_max[s])
          active[n_active++] = s;
      if (n_active == 0)
        continue;
      float sp = TrigLUT<W, H>::sin_phi[y];
      float cp = TrigLUT<W, H>::cos_phi[y];
      auto emit = [&](int x1, int x2) { pixel_run(x1, x2, y, sp, cp); };

      if (cover.needs_full_row_scan(sp)) {
        clip_run(0, W, xc, emit);
        continue;
      }
      intervals.clear();
      cover.template get_horizontal_intervals<W, H>(y, [&](float t1, float t2) {
        SDF::push_interval(intervals, t1, t2);
      });
      // The full-row case is taken above, so the covering ring always answers.
      emit_row_runs<W>(true, intervals, norm, xc, emit);
    }
  }
};
HS_O3_END

/**
 * @brief Draws a disc as a zero-radius ring whose stroke spans it.
 * @details Stroke coverage falls off quintically from 1 at the center to 0 at
 * the rim; it modulates the plotted alpha and reaches the shader as the
 * fragment's register 2, so the shader picks the radial style on top of it.
 */
struct Circle {
  /**
   * @brief Draws a disc from an orientation basis.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the circle plane.
   * @param radius Circle radius in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    float th = radius * (PI_F / 2.0f);
    Ring::draw<W, H, ComputeUVs>(pipeline, canvas, basis, 0.0f, th,
                                 fragment_shader, 0, debug_bb);
  }

  /**
   * @brief Draws a disc from a plane-normal vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param normal Plane normal as a world-space vector.
   * @param radius Circle radius in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &normal,
                   float radius, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    Basis basis = make_basis(Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, fragment_shader,
                           debug_bb);
  }
};

/**
 * @brief Draws a dot as a zero-radius ring of the given thickness.
 * @details Stroke coverage falls off quintically from 1 at the center to 0 at
 * the rim, giving the soft glow effects shade through the fragment's
 * register 2.
 */
struct Point {
  /**
   * @brief Draws a dot centered on a unit vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param p Point center as a world-space unit vector.
   * @param thickness Point radius (ring thickness) in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &p,
                   float thickness, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    // Point is a Ring with radius 0.
    Basis basis = make_basis(Quaternion(), p);
    Ring::draw<W, H>(pipeline, canvas, basis, 0.0f, thickness, fragment_shader,
                     0.0f, debug_bb);
  }
};

/**
 * @brief Draws a solid star shape.
 */
struct Star {
  /**
   * @brief Rasterizes a solid star.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the star plane.
   * @param radius Star circumradius in world units.
   * @param sides Number of star points.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Star shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /** @brief Rasterizes a constant-color solid star. */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Star shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a solid flower shape.
 */
struct Flower {
  /**
   * @brief Rasterizes a solid flower.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the flower plane.
   * @param radius Flower circumradius in world units.
   * @param sides Number of flower petals.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /** @brief Rasterizes a constant-color solid flower. */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a solid spherical polygon.
 */
struct SphericalPolygon {
  /**
   * @brief Rasterizes a solid spherical polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon.
   * @param radius Polygon circumradius in world units.
   * @param sides Number of polygon sides.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    float offset = PI_F / sides;

    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color solid spherical polygon.
   * @tparam SineDistance Use edge-plane distance for the AA band.
   */
  template <int W, int H, bool SineDistance = false, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    float offset = PI_F / sides;
    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                radius > 1.0f);
    Scan::rasterize_solid<W, H, SineDistance>(pipeline, canvas, shape, color,
                                              debug_bb);
  }
};

} // namespace Scan
