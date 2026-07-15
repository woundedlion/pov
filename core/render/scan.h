/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <utility>
#include "render/sdf.h"
#include "mesh/mesh_classes.h"
#include "color/color.h"
#include "render/filter.h"
#include "engine/static_circular_buffer.h"
#include "render/canvas.h"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif
#include "engine/platform.h"

/**
 * @brief The Scan namespace contains volumetric (raster) drawing primitives.
 * @details General register mapping for Scan primitives:
 *  v0: Normalized parameter t (0-1) or angle
 *  v1: Raw distance or supplementary value
 *  v2: Stroke AA coverage (0-1), also applied by Scan at plot time; 0 for
 *      solid shapes; Mesh uses face index
 */
namespace Scan {

/**
 * @brief Processes a single pixel for rasterization: evaluates the shape SDF,
 *        computes anti-aliased coverage, runs the fragment shader, and plots.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
 * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
 * @param x Column index in [0, W).
 * @param y Row index in [0, H).
 * @param p World-space unit vector for the pixel center.
 * @param pipeline Plotting pipeline receiving the final color.
 * @param canvas Destination canvas.
 * @param shape SDF shape providing distance<ComputeUVs>().
 * @param fragment_shader Shader invoked to color the fragment.
 * @param debug_bb When true, forces plotting and tints the bounding box.
 * @param result_scratch Reused DistanceResult scratch (avoids per-pixel alloc).
 * @param frag_scratch Reused Fragment scratch (avoids per-pixel alloc).
 * @details The shader and pipeline are type-erased (FragmentShaderFn /
 * PipelineRef) so the scanline instantiates once per <W,H> rather than per
 * shader/pipeline; see the PipelineRef note in concepts.h.
 */
template <int W, int H, bool ComputeUVs = true,
          typename PipelineT = PipelineRef>
inline void process_pixel(int x, int y, const Vector &p, PipelineT &pipeline,
                          Canvas &canvas, const auto &shape,
                          FragmentShaderFn fragment_shader, bool debug_bb,
                          SDF::DistanceResult &result_scratch,
                          Fragment &frag_scratch) {
  shape.template distance<ComputeUVs>(p, result_scratch);

  float d = result_scratch.dist;
  float pixel_width = 2.0f * PI_F / W;
  constexpr bool solid = std::remove_cvref_t<decltype(shape)>::is_solid;
  float threshold = solid ? pixel_width : 0.0f;

  if (debug_bb || d < threshold) {
    float alpha = 1.0f;

    if constexpr (solid) {
      if (d <= -pixel_width) {
        // FAST PATH: Pixel is fully inside the shape. Skip AA
        alpha = 1.0f;
      } else {
        // SLOW PATH: Pixel is on the boundary. Calculate AA
        float t_aa = 0.5f - d / (2.0f * pixel_width);
        alpha = quintic_kernel(std::max(0.0f, std::min(1.0f, t_aa)));
      }
    } else {
      // Stroke falloff over the winning leaf's own half-width: result.size, not
      // shape.thickness (a CSG composite's wrapper carries a min/max thickness).
      // Inward-only ramp: d = centerline_dist - half_width, so d=0 is the tube
      // edge (alpha 0) and d=-size the centerline (alpha 1).
      float aa_thickness = result_scratch.size;
      if (aa_thickness > 0) {
        alpha = quintic_kernel(-d / aa_thickness);
      } else {
        alpha = 0.0f;
      }
    }

    if (!debug_bb && alpha <= 0.001f)
      return;

    // Scratch Fragment is reused across pixels; reset color each call so a
    // conditionally-writing shader starts from a clean color/alpha (matches
    // Plot::rasterize).
    frag_scratch.color = Color4(0, 0, 0, 0);
    frag_scratch.pos = p;
    frag_scratch.v0 = result_scratch.t;
    frag_scratch.v1 = result_scratch.raw_dist;
    frag_scratch.v2 = solid ? 0.0f : alpha;
    frag_scratch.v3 = result_scratch.aux;
    frag_scratch.size = result_scratch.size;
    frag_scratch.age = 0;

    fragment_shader(p, frag_scratch);

    if (debug_bb) {
      frag_scratch.color.color = frag_scratch.color.color.lerp16(
          Pixel(65535, 65535, 65535), 65535 / 2);
      frag_scratch.color.alpha = 1.0f;
      alpha = 1.0f;
    }

    if (frag_scratch.color.alpha > 0.001f) {
      pipeline.plot(canvas, x, y, frag_scratch.color.color, frag_scratch.age,
                    frag_scratch.color.alpha * alpha);
    }
  }
}

/**
 * @brief Shared pixel iteration utility for bounded spherical regions.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam IntervalFn Callable producing per-row longitude intervals.
 * @tparam PixelFn Callable invoked per visited pixel.
 * @param y_min First row to scan (inclusive).
 * @param y_max Last row to scan (inclusive).
 * @param get_intervals (int y, auto &out) -> bool. Pushes {float,float}
 *                      intervals via out(start, end). Returns true if intervals
 *                      were produced, false for full-row scan.
 * @param pixel_fn (int wx, int y, const Vector &p) -> void, called per pixel.
 * @param xc Column-arc clip; pixel runs are intersected with the arc before
 *           walking, so pixel_fn is never called for a clipped column.
 * @details Iterates y in [y_min, y_max], collects float intervals per row via
 * get_intervals, wraps x coordinates, and calls pixel_fn(wx, y, p) per pixel.
 *
 * Producer contract: emitted endpoints are in fractional column units and need
 * NOT lie in [0,W) — a start may be negative or a span may straddle θ=0 (this
 * is the single point that wraps and seam-splits them into range). Each interval
 * MUST have length <= W, though: a span longer than the full circle is a
 * full-row scan, and the seam-split norm buffer is sized for exactly one split
 * per input span. CSG shapes satisfy this by construction; BoundingSphere does
 * so via its min(W/2, ...) half-width clamp.
 */
template <int W, int H, typename IntervalFn, typename PixelFn>
inline void scan_region(int y_min, int y_max, IntervalFn &&get_intervals,
                        PixelFn &&pixel_fn, ClipRegion::XClip xc = {}) {
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();

  // Interval scratch (~1.5 KiB) lives in scratch_arena_b, not the stack:
  // Phantasm's DTCM stack is tight and scan_region is on the deepest render
  // chain. Per-call bump scope; norm is cleared per row below.
  //
  // intervals holds a top-level shape's full per-row emission. The widest top
  // node is Subtract/Intersection with both children at INTERVAL_SPAN_CAP:
  // |A|+|B|+2 == 2*INTERVAL_SPAN_CAP+2. norm holds one seam-split per span, so 2x.
  ScratchScope scratch(scratch_arena_b);
  static constexpr size_t TOP_SPAN_CAP = 2 * SDF::INTERVAL_SPAN_CAP + 2;
  using IntervalBuf =
      StaticCircularBuffer<std::pair<float, float>, TOP_SPAN_CAP>;
  using NormBuf =
      StaticCircularBuffer<std::pair<float, float>, 2 * TOP_SPAN_CAP>;
  static_assert(IntervalBuf::CAPACITY >= SDF::MergedIntervalBuffer::CAPACITY,
                "scan_region intervals must hold the largest top-level CSG "
                "emission (Subtract/Intersection: |A|+|B|+2)");
  static_assert(NormBuf::CAPACITY == 2 * IntervalBuf::CAPACITY,
                "norm must hold 2 spans per input interval (seam split)");
  auto &intervals = *new (scratch_arena_b.allocate(
      sizeof(IntervalBuf), alignof(IntervalBuf))) IntervalBuf();
  auto &norm = *new (scratch_arena_b.allocate(sizeof(NormBuf),
                                              alignof(NormBuf))) NormBuf();

  const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();

  auto walk = [&](int x1, int x2, int y, float sp, float cp) {
    for (int x = x1; x < x2; ++x)
      pixel_fn(x, y, Vector(sp * cos_theta[x], cp, sp * sin_theta[x]));
  };
  // Intersect a coalesced run [x1, x2) with the clip arc before walking. A
  // wrapping arc is [rs, W) ∪ [0, re); the two pieces are disjoint (re <= rs),
  // so no column is walked twice.
  auto walk_clipped = [&](int x1, int x2, int y, float sp, float cp) {
    if (!xc.active) {
      walk(x1, x2, y, sp, cp);
    } else if (xc.wrap) {
      walk(std::max(x1, xc.rs), x2, y, sp, cp);
      walk(x1, std::min(x2, xc.re), y, sp, cp);
    } else {
      walk(std::max(x1, xc.rs), std::min(x2, xc.re), y, sp, cp);
    }
  };

  // Inverted range (y_min > y_max) is a no-op: a disjoint CSG Intersection or a
  // fully-culled Face reports y_min=1, y_max=0, and the loop never runs.
  for (int y = y_min; y <= y_max; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];

    bool handled = get_intervals(
        y, [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });

    if (handled && !intervals.is_empty()) {
      // A single span covering the full circle (len >= W) paints every column;
      // detect it up front and skip the seam-split/sort/coalesce path. Coverage
      // assembled from multiple abutting spans is not caught here — it falls to
      // the slow path, which still paints every covered column.
      bool full_row = false;
      for (const auto &iv : intervals) {
        if (iv.second - iv.first >= static_cast<float>(W)) {
          full_row = true;
          break;
        }
      }

      if (full_row) {
        walk_clipped(0, W, y, sp, cp);
      } else {
        // Wrap each start into [0, W) and split any span crossing the x=0 seam so
        // the forward-sweep coalescer sees sorted, non-wrapping spans even when a
        // shape/CSG straddles θ=0; norm holds up to 2 spans per input interval.
        norm.clear();
        SDF::normalize_intervals_to_range<W>(intervals, norm);
        SDF::sort_intervals_by_start(norm);

        float current_end = -FLT_MAX;
        // Coalesce in integer pixel space too: last_x2 clamps each run's start
        // past the prior run's end so two spans sharing a fractional column don't
        // both paint it (double process_pixel / alpha).
        int last_x2 = 0;
        for (const auto &iv : norm) {
          if (iv.second <= current_end)
            continue;
          float start = std::max(iv.first, current_end);
          float end = iv.second;
          current_end = end;

          int x1 = static_cast<int>(floorf(start));
          int x2 = static_cast<int>(ceilf(end));
          // A zero-width interval (x1 == x2) still owns its pixel column; widen
          // to paint it. The last_x2 clamp then suppresses a following interval
          // mapping to the same column.
          if (x1 == x2)
            x2++;
          // Clamp to canvas columns. After the wrap/split above start>=0 and
          // end<=W, so x2>W only from the x1==x2 widen at the right edge; x1<0 is
          // a defensive floor. last_x2 keeps the integer sweep monotone.
          if (x1 < 0)
            x1 = 0;
          if (x2 > W)
            x2 = W;
          if (x1 < last_x2)
            x1 = last_x2;
          last_x2 = x2;
          walk_clipped(x1, x2, y, sp, cp);
        }
      }
      intervals.clear();
    } else if (!handled) {
      walk_clipped(0, W, y, sp, cp);
    }
  }
}

/**
 * @brief Computes bounding sphere y-range and per-row x intervals.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> struct BoundingSphere {
  int y_min, y_max;
  float center_theta; /**< Longitude of center in pixel units. */
  float angular_radius;
  float cos_rho;        /**< cos(angular_radius). */
  float cos_center_phi; /**< cos of the cap center's colatitude. */
  float sin_center_phi; /**< sin of the cap center's colatitude. */

  /**
   * @brief Constructs the bounding sphere clamped to the canvas rows.
   * @param center World-space unit vector at the sphere center.
   * @param bounds_radius Bounding radius in world units (sin of angular extent).
   */
  BoundingSphere(const Vector &center, float bounds_radius)
      : angular_radius(asinf(std::min(bounds_radius, 1.0f))) {
    center_theta = vector_to_theta<W>(center);
    float center_phi = acosf(hs::clamp(center.y, -1.0f, 1.0f));
    cos_rho = cosf(angular_radius);
    cos_center_phi = cosf(center_phi);
    sin_center_phi = sinf(center_phi);
    // Round the band outward (floor the top, ceil the bottom) so a fractional cap
    // edge keeps the fringe row it touches.
    y_min = std::max(
        0, static_cast<int>(floorf(phi_to_y<H>(center_phi - angular_radius))));
    y_max = std::min(H - 1, static_cast<int>(ceilf(
                                phi_to_y<H>(center_phi + angular_radius))));
  }

  /**
   * @brief Pushes a single longitude interval for row y based on the span at
   *        that latitude.
   * @tparam OutFn Callable accepting (float start, float end) in pixel units.
   * @param y Row index in [0, H).
   * @param out Sink receiving the interval for this row.
   * @return Always true (an interval is always produced).
   */
  template <typename OutFn> bool get_intervals(int y, OutFn &&out) const {
    // Phi trig from the static LUT (bit-identical to sinf(y_to_phi(y))), as on the
    // rest of the Volume hot path.
    float sin_phi = TrigLUT<W, H>::sin_phi[y];
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float theta_span;
    // Exact cap longitude half-width at this row from the spherical law of cosines:
    // cos(dtheta) = (cos rho - cos phi cos phi_c) / (sin phi sin phi_c). denom <= 0
    // only at a pole (this row or the cap center), where the cap spans the whole
    // row; that branch also dodges the 0/0 → NaN → int-cast UB.
    float denom = sin_phi * sin_center_phi;
    float cos_dtheta =
        denom > 0.0f ? (cos_rho - cos_phi * cos_center_phi) / denom : -1.0f;
    if (cos_dtheta <= -1.0f) {
      theta_span = static_cast<float>(W);
    } else {
      float dtheta = acosf(cos_dtheta < 1.0f ? cos_dtheta : 1.0f);
      theta_span = dtheta * W / (2.0f * PI_F);
    }
    // +1 absorbs ceil/round-off at the span edges; the downstream per-pixel
    // ray-sphere test rejects any extra column. The W/2 .. (W+1)/2 caps bound
    // the span length at W (scan_region's producer contract). Endpoints are
    // not clamped to [0,W); scan_region wraps them.
    int span = static_cast<int>(ceilf(theta_span)) + 1;
    int x_lo = std::min(W / 2, span);
    int x_hi = std::min((W + 1) / 2, span);
    out(center_theta - x_lo, center_theta + x_hi);
    return true;
  }
};

/**
 * @brief Scoped accumulator for per-draw render time (telemetry only).
 * @details Measures wall-clock over its lifetime and adds it to the canvas's
 * render-time counter on destruction. Off-Emscripten there is no JS perf clock,
 * so the type is empty and every use optimizes away to nothing.
 */
struct ScopedRenderTimer {
#ifdef __EMSCRIPTEN__
  Canvas &canvas_;
  double t0_;
  explicit ScopedRenderTimer(Canvas &canvas)
      : canvas_(canvas), t0_(emscripten_get_now()) {}
  ~ScopedRenderTimer() { canvas_.add_render_us(emscripten_get_now() - t0_); }
#else
  explicit ScopedRenderTimer(Canvas &) {}
#endif
  ScopedRenderTimer(const ScopedRenderTimer &) = delete;
  ScopedRenderTimer &operator=(const ScopedRenderTimer &) = delete;
};

/**
 * @brief Main rasterization routine for SDF shapes.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
 * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape SDF shape providing vertical bounds, horizontal intervals, and
 *              distance().
 * @param fragment_shader Shader invoked per covered pixel.
 * @param debug_bb When true, renders the bounding box for debugging.
 * @details Scans the bounding box, computes intervals, and executes the shader
 * for valid pixels.
 */
template <int W, int H, bool ComputeUVs = true,
          typename PipelineT = PipelineRef>
inline void rasterize(PipelineT &pipeline, Canvas &canvas, const auto &shape,
                      FragmentShaderFn fragment_shader, bool debug_bb = false) {
  bool effective_debug = debug_bb || canvas.debug();

  int y_lo, y_hi;
  const auto &cr = canvas.clip();
  const auto xc = cr.x_clip();
  auto bounds = shape.template get_vertical_bounds<H>();
  y_lo = bounds.y_min > cr.render_y_start() ? bounds.y_min : cr.render_y_start();
  y_hi = bounds.y_max < cr.render_y_end() - 1 ? bounds.y_max
                                               : cr.render_y_end() - 1;
  if (y_lo > y_hi)
    return;

  SDF::DistanceResult result_scratch;
  Fragment frag_scratch;

  ScopedRenderTimer timer_guard(canvas);
  scan_region<W, H>(
      y_lo, y_hi,
      [&](int y, auto &&out) {
        return shape.template get_horizontal_intervals<W, H>(y, out);
      },
      [&](int wx, int y, const Vector &p) {
        process_pixel<W, H, ComputeUVs>(wx, y, p, pipeline, canvas, shape,
                                        fragment_shader, effective_debug,
                                        result_scratch, frag_scratch);
      },
      xc);
}

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
   * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
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
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw_flat(PipelineT &pipeline, Canvas &canvas,
                        const Basis &basis, float radius, float thickness,
                        FragmentShaderFn fragment_shader, float phase = 0,
                        bool debug_bb = false, bool suppress_pole_fill = false) {
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
   * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
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
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const Basis &basis,
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
   * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param knots lut_n + 1 centerline shifts, entry lut_n repeating entry 0;
   *        must outlive the call.
   * @param lut_n Number of knot cells.
   * @param amplitude Modulation amplitude in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   * @param suppress_pole_fill Drop the degenerate exact-pole row instead of
   *        full-row filling it (dense ring stacks; see get_horizontal_intervals).
   */
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness, const float *knots,
                   int lut_n, float amplitude,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false, bool suppress_pole_fill = false) {
    SDF::DistortedRing shape(basis, radius, thickness, knots, lut_n, amplitude,
                             phase);
    shape.suppress_pole_fill = suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

/**
 * @brief Fused single-pass rasterizer for a stack of same-axis distorted
 *        rings.
 */
struct DistortedRingStack {
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
   * @param shapes n_slots knot-mode rings sharing one Basis, in ascending ring
   *        order; culled rings are simply absent.
   * @param slot_by_ring n_rings entries mapping ring index -> slot in shapes,
   *        -1 for culled rings.
   * @param n_slots Number of shapes; at least 1.
   * @param shader Per-ring fragment shader (see RingShaderT).
   * @details The per-pixel frame shared by every ring at a pixel (axis dot,
   * fast_acos, fast_atan2) is computed once, the candidate rings fall out of
   * the frame's polar angle by arithmetic (the stack is evenly spaced), and
   * each candidate runs its own cos reject + exact polyline distance via
   * SDF::DistortedRing::distance_from_frame. Candidates evaluate in ascending
   * ring index, so per-pixel blend order — and therefore output — is
   * bit-identical to rasterizing the rings one by one; only the redundant
   * per-ring frame recompute is elided. The aliased exact-pole rows are
   * dropped, matching the per-ring path under suppress_pole_fill.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  static void draw(PipelineT &pipeline, Canvas &canvas, int n_rings,
                   const SDF::DistortedRing *shapes,
                   const int8_t *slot_by_ring, int n_slots,
                   RingShaderT &&shader) {
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
    const float b_win = b_max + 1e-3f;
    const float inv_delta = (n_rings + 1) / PI_F;

    // The per-ring path suppresses the aliased exact-pole rows
    // (suppress_pole_fill); its full-scan fallback for a near-canvas-pole
    // axis (r_val below the projection floor) scans every row.
    SDF::AxisProjection ap = SDF::project_axis(shapes[0].normal);
    const bool skip_pole_rows = ap.R_val >= SDF::MIN_HORIZONTAL_PROJ;

    const Vector axis_v = shapes[0].normal;
    const Vector axis_u = shapes[0].u;
    const Vector axis_w = shapes[0].w;

    SDF::DistanceResult res;
    Fragment frag;
    for (int y = y_lo; y <= y_hi; ++y) {
      const float sp = TrigLUT<W, H>::sin_phi[y];
      const float cp = TrigLUT<W, H>::cos_phi[y];
      if (skip_pole_rows && std::abs(ap.R_val * sp) < SDF::INTERVAL_DENOM_EPS)
        continue;
      for (int x = 0; x < W; ++x) {
        if (xc.clipped(x))
          continue;
        Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
        const float d = dot(p, axis_v);
        const float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
        int ilo = static_cast<int>(ceilf((polar - b_win) * inv_delta)) - 2;
        int ihi = static_cast<int>(floorf((polar + b_win) * inv_delta));
        if (ilo < 0)
          ilo = 0;
        if (ihi > n_rings - 1)
          ihi = n_rings - 1;
        if (ilo > ihi)
          continue;
        const float dot_u = dot(p, axis_u);
        const float dot_w = dot(p, axis_w);
        float azimuth = fast_atan2(dot_w, dot_u);
        if (azimuth < 0)
          azimuth += 2 * PI_F;
        const float t_norm = wrap_t(azimuth / (2 * PI_F));
        const float sin_polar = sqrtf(
            std::max(1.0f - d * d, SDF::DistortedRing::POLE_SIN2_FLOOR));
        constexpr int PC = SDF::DistortedRing::PREFILTER_CHUNKS;
        const float chunk_u = (2.0f * PI_F / PC) * sin_polar;
        int pc = static_cast<int>(t_norm * PC);
        if (pc >= PC)
          pc = PC - 1;
        for (int i = ilo; i <= ihi; ++i) {
          const int s = slot_by_ring[i];
          if (s < 0)
            continue;
          const SDF::DistortedRing &shape = shapes[s];
          // Same 3-chunk gap test the polyline search opens with, applied
          // before the call: a gap-rejected ring cannot plot (its distance
          // would come back > thickness), so skipping it preserves output.
          bool prefiltered = false;
          if (chunk_u >= shape.thickness) {
            float lo, hi;
            shape.prefilter_band(pc, lo, hi);
            const float base = shape.target_angle - polar;
            if (std::max(base + lo, -(base + hi)) > shape.thickness)
              continue;
            prefiltered = true;
          }
          shape.distance_from_frame(d, polar, sin_polar, t_norm, res,
                                    prefiltered);
          const float dd = res.dist;
          if (dd >= 0.0f)
            continue;
          // process_pixel's stroke epilogue with a slot-aware shader.
          const float aa = res.size;
          const float alpha = aa > 0.0f ? quintic_kernel(-dd / aa) : 0.0f;
          if (alpha <= 0.001f)
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
          if (frag.color.alpha > 0.001f)
            pipeline.plot(canvas, x, y, frag.color.color, frag.age,
                          frag.color.alpha * alpha);
        }
      }
    }
  }
};

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
    float thickness = res.second * (PI_F / 2.0f);

    SDF::PlanarPolygon shape(res.first, thickness, sides, phase,
                             radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
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
   * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const Basis &basis,
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
   * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param normal Plane normal as a world-space vector.
   * @param radius Ring radius in world units.
   * @param thickness Ring stroke thickness in world units.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const Vector &normal,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    Basis basis = make_basis(Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, thickness,
                           fragment_shader, phase, debug_bb);
  }
};

/**
 * @brief Draws a solid circle (filled ring).
 */
struct Circle {
  /**
   * @brief Draws a solid circle from an orientation basis.
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
    // Circle is a Ring with inner radius 0.
    float th = radius * (PI_F / 2.0f);
    Ring::draw<W, H, ComputeUVs>(pipeline, canvas, basis, 0.0f, th,
                                 fragment_shader, 0, debug_bb);
  }

  /**
   * @brief Draws a solid circle from a plane-normal vector.
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
 * @brief Draws a solid point (thick circle).
 */
struct Point {
  /**
   * @brief Draws a solid point centered on a unit vector.
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
};

/**
 * @brief Draws a solid spherical polygon.
 */
struct SphericalPolygon {
  /**
   * @brief Rasterizes a solid spherical polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon.
   * @param radius Polygon circumradius in world units.
   * @param sides Number of polygon sides.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    float offset = PI_F / sides;

    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                radius > 1.0f);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

/**
 * @brief Rasterizes a polygonal mesh by drawing each face as an SDF::Face,
 *        threading the face index through register v2 so the shader can vary
 *        color per face.
 */
struct Mesh {
  /**
   * @brief Rasterizes every face of a mesh.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type (defaults to type-erased PipelineRef).
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @param fragment_shader Shader invoked per covered pixel; receives the face
   *                        index in register v2. v2 is a float, so the index is
   *                        exact only up to 2^24 faces.
   * @param scratch_arena Arena supplying per-face SDF::Face scratch storage.
   * @param debug_bb When true, renders the bounding box for debugging.
   * @param bake Optional congruence-class bake for this mesh (null = exact
   *        path everywhere, today's behavior). When present, each face is
   *        aligned to its canonical class shape after construction and the
   *        class distance LUT is bound for the probe loop.
   */
  template <int W, int H, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   bool debug_bb = false,
                   const MeshOps::MeshClassBake *bake = nullptr) {
    ScratchScope scope(scratch_arena);
    auto *scratch =
        static_cast<SDF::FaceScratchBuffer *>(scratch_arena.allocate(
            sizeof(SDF::FaceScratchBuffer), alignof(SDF::FaceScratchBuffer)));

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    const uint16_t *fo = mesh.get_face_offsets_data();
    size_t fi_size = mesh.get_faces_size();

    // An empty bake (build skipped) is equivalent to none; a populated one
    // must cover every face — records are indexed by face order.
    if (bake && !bake->face_recs.is_bound())
      bake = nullptr;
    HS_CHECK(!bake || bake->face_recs.size() == num_f,
             "mesh class bake face count disagrees with the mesh");

    for (size_t i = 0; i < num_f; ++i) {
      size_t count = fc[i];

      // Trap malformed mesh data: an offset/count pair disagreeing with the flat
      // index array yields an out-of-bounds span for SDF::Face. Cold per-face check.
      HS_CHECK(static_cast<size_t>(fo[i]) + count <= fi_size,
               "mesh face span exceeds face index array");

      std::span<const Vector> verts(mesh.vertices.data(), mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      // IIFE so the HS_PROFILE scope measures Face construction alone (prvalue
      // elided in place), not the rasterize below.
      SDF::Face shape = [&] {
        HS_PROFILE(scan_face_setup);
        return SDF::Face(verts, indices, 0.0f, *scratch, H + hs::H_OFFSET, H,
                         &canvas.clip());
      }();

      // Bind the face's congruence-class LUT: a vertex correlation aligns the
      // current projection to the canonical frame. Culled faces (y_min > y_max)
      // skip it; a missing LUT or degenerate alignment stays on the exact path.
      if (bake && shape.y_min <= shape.y_max) {
        const MeshOps::FaceClassRec &rec = bake->face_recs[i];
        if (rec.class_id != MeshOps::NO_CLASS) {
          const MeshOps::CongruenceClass &cls = bake->classes[rec.class_id];
          if (cls.lut.data && cls.n_verts == shape.count)
            shape.bind_class_lut(&cls.lut, cls.canon_xy, rec.vert_offset,
                                 rec.reflected != 0);
        }
      }

      auto wrapper = [&](const Vector &p, Fragment &f_in) {
        // Exact for i < 2^24 (float mantissa); meshes never approach that face count.
        f_in.v2 = static_cast<float>(i);
        fragment_shader(p, f_in);
      };

      { HS_PROFILE(scan_mesh_raster);
        Scan::rasterize<W, H, true>(pipeline, canvas, shape, wrapper, debug_bb);
      }
    }
  }
};

/**
 * @brief Full-screen per-pixel shader with SAMPLES× SSAA.
 *
 * Accepts a single callable ShaderFn(const Vector &v) -> Color4
 * that maps a world-space unit vector to a final color.
 * The utility calls it SAMPLES× per pixel at sub-pixel offsets and averages.
 */
struct Shader {
  // --- Shared SSAA helpers (used by both draw() overloads) -------------------
  /**
   * @brief Per-draw sub-pixel trig for the 2×2 SSAA sample grid, derived from
   *        the resident engine trig LUT.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @details Each sub-sample sits a constant ±0.25 px from an integer pixel,
   * i.e. a constant angular offset from that pixel's theta/phi under the same
   * parameterization as pixel_to_vector (theta = 2*pi*px/W,
   * phi = py*pi/(H_VIRT-1)). So sin/cos at a sub-sample follow from the
   * integer-pixel TrigLUT<W,H> tables by one angle-addition, keeping the
   * samples LUT-consistent with the non-SSAA path and needing no per-column
   * storage. The two rotation angles (d_theta, d_phi) are resolution constants;
   * their sin/cos are built once per draw (libm is not constexpr) into a
   * handful of floats.
   */
  template <int W, int H> struct SsaaGrid {
    float sin_phi[2];  /**< Current row's sin(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_phi[2];  /**< Current row's cos(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_dtheta;  /**< cos of the ±0.25 px column rotation. */
    float sin_dtheta;  /**< sin of the ±0.25 px column rotation. */
    float cos_dphi;    /**< cos of the ±0.25 px row rotation. */
    float sin_dphi;    /**< sin of the ±0.25 px row rotation. */

    SsaaGrid() {
      // d_theta = 0.25 px * (2*pi/W); d_phi = 0.25 px * (pi/(H_VIRT-1)).
      constexpr float d_theta = 0.5f * PI_F / static_cast<float>(W);
      constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
      constexpr float d_phi = 0.25f * PI_F / h_virt_minus_1;
      cos_dtheta = cosf(d_theta);
      sin_dtheta = sinf(d_theta);
      cos_dphi = cosf(d_phi);
      sin_dphi = sinf(d_phi);
    }

    /** @brief Loads the two phi trig pairs for pixel row y from the LUT. */
    void set_row(int y) {
      const float sy = TrigLUT<W, H>::sin_phi[y];
      const float cy = TrigLUT<W, H>::cos_phi[y];
      // Row 0 = y+0.25, row 1 = y-0.25 (the 2×2 grid's centered ±0.25 offsets).
      sin_phi[0] = sy * cos_dphi + cy * sin_dphi;
      cos_phi[0] = cy * cos_dphi - sy * sin_dphi;
      sin_phi[1] = sy * cos_dphi - cy * sin_dphi;
      cos_phi[1] = cy * cos_dphi + sy * sin_dphi;
    }

    /**
     * @brief World-space unit vector for sample i of pixel x in the current
     *        row (see set_row).
     * @param x Pixel column.
     * @param i Sample index in [0, 4); the low bit selects the column offset
     * (±0.25 px) and bit 1 the row offset (±0.25 px).
     */
    Vector at(int x, int i) const {
      const float st = TrigLUT<W, H>::sin_theta[x];
      const float ct = TrigLUT<W, H>::cos_theta(x);
      // Column 0 (i&1==0) = x+0.25, column 1 = x-0.25.
      const float s = (i & 1) ? -sin_dtheta : sin_dtheta;
      const float sin_theta = st * cos_dtheta + ct * s;
      const float cos_theta = ct * cos_dtheta - st * s;
      const float sp = sin_phi[(i >> 1) & 1];
      return Vector(sp * cos_theta, cos_phi[(i >> 1) & 1], sp * sin_theta);
    }
  };

  /**
   * @brief Validates the LUT-domain invariant shared by both draw overloads.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param cr Clip region whose bounds are checked against the LUT extents.
   * @details Checked once per draw, not per pixel: every (x,y) the loops feed to
   * pixel_to_vector indexes the trig LUTs within bounds.
   */
  template <int W, int H>
  static void check_lut_domain(const ClipRegion &cr) {
    HS_CHECK(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 &&
             cr.render_y_end() <= PhiLUT<H>::H_VIRT);
  }
  // --------------------------------------------------------------------------

  /**
   * @brief Full-screen per-pixel shader with SAMPLES× SSAA from a single
   *        callable.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam SAMPLES Number of sub-pixel samples per pixel (1 disables SSAA).
   * @tparam ShaderFn Callable ShaderFn(const Vector &v) -> Color4.
   * @param canvas Destination canvas.
   * @param shader Maps a world-space unit vector to a final color; invoked
   *               SAMPLES× per pixel at sub-pixel offsets and averaged.
   */
  template <int W, int H, int SAMPLES = 1, typename ShaderFn>
  static void draw(Canvas &canvas, ShaderFn &&shader) {
    // The sample-offset table has four distinct sub-pixel positions; only 1 and
    // the 2x2 grid (4) are valid.
    static_assert(SAMPLES == 1 || SAMPLES == 4,
                  "Scan::Shader SSAA supports only SAMPLES == 1 or 4");
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto xc = cr.x_clip();

    if constexpr (SAMPLES == 1) {
      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = 0; x < W; ++x) {
          if (xc.clipped(x))
            continue;
          Vector v = pixel_to_vector<W, H>(x, y);
          Color4 sample = shader(v);
          canvas(x, y) = sample.color * sample.alpha;
        }
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      if (!TrigLUT<W, H>::initialized)
        TrigLUT<W, H>::init();
      SsaaGrid<W, H> grid;

      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        grid.set_row(y);
        for (int x = 0; x < W; ++x) {
          if (xc.clipped(x))
            continue;
          // Premultiplied SSAA: accumulate each sample's coverage-weighted color
          // (color * alpha / N), matching the SAMPLES==1 path.
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            Color4 sample = shader(grid.at(x, i));
            accum += sample.color * (sample.alpha * inv_samples);
          }

          canvas(x, y) = accum;
        }
      }
    }
  }

  /**
   * @brief Full-screen per-pixel shader with SAMPLES× SSAA and a split
   *        vertex/fragment shader.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam SAMPLES Number of sub-pixel samples per pixel (1 disables SSAA).
   * @param canvas Destination canvas.
   * @param fragment_shader Per-sub-sample shader, called SAMPLES× per pixel.
   * @param vertex_shader Per-pixel shader, called once at the pixel center.
   * @details Splits expensive per-pixel work (vertex_shader, once at pixel
   * center) from per-sub-sample evaluation (fragment_shader, SAMPLES×).
   *
   * @note SAMPLES defaults to 1 (no SSAA), matching the single-callback overload.
   * Every call site passes SAMPLES explicitly.
   */
  template <int W, int H, int SAMPLES = 1>
  static void draw(Canvas &canvas, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    // Only 1 and 4 are supported (see the single-callback overload).
    static_assert(SAMPLES == 1 || SAMPLES == 4,
                  "Scan::Shader SSAA supports only SAMPLES == 1 or 4");
    // Cold (once per draw), not per-pixel: trap null shaders here so they fail
    // deterministically instead of calling a null thunk under NDEBUG.
    HS_CHECK(vertex_shader,
             "Scan::Shader::draw requires a non-null vertex_shader");
    HS_CHECK(fragment_shader,
             "Scan::Shader::draw requires a non-null fragment_shader");
    // frag_base is per pixel, not per draw: each pixel starts from a default
    // Fragment, so a vertex shader writing only some registers (v0-v3/size/age/
    // color) can't inherit the previous pixel's values.
    if constexpr (SAMPLES == 1) {
      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = 0; x < W; ++x) {
          if (xc.clipped(x))
            continue;
          Vector center_v = pixel_to_vector<W, H>(x, y);
          Fragment frag_base;
          frag_base.pos = center_v;
          vertex_shader(frag_base);
          fragment_shader(center_v, frag_base);
          // Premultiply by alpha, matching the single-callback overload and the
          // process_pixel/Volume contract.
          canvas(x, y) = frag_base.color.color * frag_base.color.alpha;
        }
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      if (!TrigLUT<W, H>::initialized)
        TrigLUT<W, H>::init();
      SsaaGrid<W, H> grid;

      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        grid.set_row(y);
        for (int x = 0; x < W; ++x) {
          if (xc.clipped(x))
            continue;
          Vector center_v = pixel_to_vector<W, H>(x, y);

          Fragment frag_base;
          frag_base.pos = center_v;
          vertex_shader(frag_base);

          // Premultiplied SSAA: accumulate coverage-weighted color directly (see
          // the single-callback overload).
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            Vector v = grid.at(x, i);

            Fragment sub_frag = frag_base;
            sub_frag.pos = v;

            fragment_shader(v, sub_frag);

            Color4 sample = sub_frag.color;
            accum += sample.color * (sample.alpha * inv_samples);
          }

          canvas(x, y) = accum;
        }
      }
    }
  }
};

/**
 * @brief Generic wrapper that places an SDF in world space via a center point
 *        and a rotation quaternion. Satisfies the Volume::draw shape concept.
 * @tparam SDF Underlying signed-distance shape type.
 * @details The quaternion q maps local→world: world_p = center +
 * rotate(local_p, q). ray_to_local uses q.inverse() to map world→local.
 */
template <typename SDF> struct TransformedVolume {
  const SDF &sdf;     /**< Underlying SDF evaluated in local space. */
  Vector center;      /**< World-space origin of the local frame. */
  Quaternion q_inv;   /**< Precomputed inverse rotation (world→local). */

  /**
   * @brief Constructs the transform from a center and a local→world rotation.
   * @param sdf Underlying SDF, stored by reference.
   * @param center World-space origin of the local frame.
   * @param q Local→world rotation; its inverse is precomputed.
   */
  TransformedVolume(const SDF &sdf, const Vector &center, const Quaternion &q)
      : sdf(sdf), center(center), q_inv(q.inverse()) {}

  /**
   * @brief Transforms a ray origin and direction from world to local space.
   * @param ro Ray origin in world space.
   * @param vd Ray direction in world space.
   * @return Pair of {local origin, local direction}.
   */
  std::pair<Vector, Vector> ray_to_local(const Vector &ro,
                                         const Vector &vd) const {
    return {rotate(ro - center, q_inv), rotate(vd, q_inv)};
  }

  /**
   * @brief Transforms only a ray origin from world to local space.
   * @param ro Ray origin in world space.
   * @return Local-space origin.
   * @details The local direction is constant across the draw, so the volume loop
   * precomputes it once and calls this per pixel to transform only the origin.
   */
  Vector origin_to_local(const Vector &ro) const {
    return rotate(ro - center, q_inv);
  }

  /**
   * @brief Evaluates the underlying SDF at a local-space point.
   * @param local_p Query point in local space.
   * @return Signed distance to the surface in local units.
   */
  float distance(const Vector &local_p) const { return sdf.distance(local_p); }
};

/**
 * @brief Raymarch volume renderer with orthographic projection.
 * @details The render loop is internal: callers provide a Shape with a
 * `float distance(const Vector&) const` method (evaluated per march step) and a
 * fragment shader (evaluated once per hit to populate Fragment registers for
 * shading).
 *
 * Coordinate-space contract:
 *   - `view_dir` is the normalized direction all rays travel (camera → scene).
 *   - Ray origins are computed via orthographic projection: each pixel's
 *     position is projected onto the plane perpendicular to `view_dir`,
 *     then offset backward along `view_dir`.
 *   - `bounds_center` and `view_dir` must both be in physical LED space.
 *   - Filter::World::Orient rotates the *output* position passed to the
 *     canvas, not the ray.
 */
struct Volume {
  /**
   * @brief Sphere-traces a ray in local space, recording the closest approach
   *        to the first surface reached.
   * @param shape Volume shape providing distance().
   * @param local_ro Ray origin in local space.
   * @param local_vd Unit ray direction in local space.
   * @param bounds_radius Bounding sphere radius (past-the-back early-out).
   * @param max_steps Maximum sphere-tracing steps.
   * @param aa_width Anti-aliasing band half-width (deep-hit early-out).
   * @param closest_local Output: local-space point of closest approach.
   * @return Signed distance at the closest approach (FLT_MAX if never sampled).
   * @details Inside the AA band, stops at the first rising local minimum (the
   * silhouette graze owning the pixel's coverage); marching past it would let a
   * deeper occluded surface steal the closest approach.
   */
  template <typename Shape>
  static __attribute__((always_inline)) float
  trace_closest(const Shape &shape, const Vector &local_ro,
                const Vector &local_vd, float bounds_radius, int max_steps,
                float aa_width, Vector &closest_local) {
    Vector local_p = local_ro;
    closest_local = local_ro;
    // Sentinel for "no surface seen yet": any real signed distance the
    // trace reports is smaller, so the first sample always wins.
    float closest_d = FLT_MAX;

    for (int i = 0; i < max_steps; ++i) {
      // Early out: ray has exited the back of the bounding sphere. The
      // local-space dot is compared against the world-space bounds_radius, valid
      // because ray_to_local is length-preserving (unit local_vd) and the caller
      // passes the shape center as bounds_center. Both are HS_CHECKed once per
      // draw at the top.
      if (local_p.x * local_vd.x + local_p.y * local_vd.y +
              local_p.z * local_vd.z >
          bounds_radius)
        break;

      float d = shape.distance(local_p);

      if (d < closest_d) {
        closest_d = d;
        closest_local = local_p;
      } else if (closest_d < aa_width) {
        // Rising past the first in-band local minimum: stop before a surface
        // behind the graze steals the closest approach.
        break;
      }

      if (d < -aa_width)
        break;

      // 1e-5 absolute stall-guard for the precision trace (fine steps near the
      // surface), bounded by max_steps and the early-out above. The probe loop
      // below instead uses a bounds_radius-relative floor for coarse punch-through.
      float advance = std::max(d * 0.9f, 1e-5f);
      local_p = Vector(local_p.x + local_vd.x * advance,
                       local_p.y + local_vd.y * advance,
                       local_p.z + local_vd.z * advance);
    }
    return closest_d;
  }

  /**
   * @brief Result of probing behind an AA halo for an occluded surface.
   */
  struct Occluder {
    bool solid;    /**< A solid surface sits behind the halo (behind is valid). */
    Vector behind; /**< Local-space hit point when solid, else the grazed
                        background edge's closest approach when soft > 0. */
    float soft;    /**< Coverage of a grazed background edge, for the corner fill. */
  };

  /**
   * @brief Marches behind an AA-halo pixel to find any surface the foreground edge
   *        occludes.
   * @param shape Volume shape providing distance().
   * @param closest_local Local-space closest approach (probe seed).
   * @param local_vd Unit ray direction in local space.
   * @param bounds_radius Bounding sphere radius (probe reach + step floor).
   * @param hit_threshold Solid-hit distance threshold.
   * @param aa_width Anti-aliasing band half-width (soft-occlusion falloff scale).
   * @return An Occluder: a solid hit point to antialias the edge over, or a
   *         grazed edge's closest approach and coverage for the corner where
   *         two edges meet.
   */
  template <typename Shape>
  static __attribute__((always_inline)) Occluder
  probe_occluder(const Shape &shape, const Vector &closest_local,
                 const Vector &local_vd, float bounds_radius, float hit_threshold,
                 float aa_width) {
    // March forward from the closest approach for a surface this halo occludes;
    // a solid hit is a self-occlusion edge (antialias over it). Step is floored
    // to punch past the stalled foreground; termination is the bounding sphere's
    // back face. With no solid hit, report a grazed background edge (local min of
    // pd) and its coverage for the corner fill.
    Vector probe = closest_local;
    float prev = FLT_MAX;  // previous step's distance
    bool climbing = false; // pd has risen off the foreground graze
    float min_behind = FLT_MAX;
    Vector min_pos = closest_local;
    // Bracket samples around the running minimum, as offsets along the ray from
    // min_pos, for the parabolic refinement below.
    float s = 0.0f, prev_s = 0.0f, min_s = 0.0f;
    float bef_s = 0.0f, bef_pd = FLT_MAX;
    float aft_s = 0.0f, aft_pd = FLT_MAX;
    bool need_aft = false;
    for (int i = 0; i < 24; ++i) {
      // Stop at the back of the bounding sphere: nothing left to occlude this halo.
      if (probe.x * local_vd.x + probe.y * local_vd.y + probe.z * local_vd.z >
          bounds_radius)
        break;
      float pd = shape.distance(probe);
      if (pd < hit_threshold)
        return {true, probe, 0.0f}; // solid surface behind the edge
      if (need_aft) {
        aft_s = s;
        aft_pd = pd;
        need_aft = false;
      }
      if (pd > prev)
        climbing = true; // moving away from the foreground graze
      else if (climbing && pd < min_behind) {
        min_behind = pd; // descending toward a surface behind
        min_pos = probe;
        min_s = s;
        bef_s = prev_s;
        bef_pd = prev;
        need_aft = true;
        aft_pd = FLT_MAX;
      }
      prev = pd;
      prev_s = s;
      float floor = bounds_radius * (i < 6 ? 0.04f : 0.12f);
      float step = std::max(pd * 0.9f, floor);
      probe = Vector(probe.x + local_vd.x * step, probe.y + local_vd.y * step,
                     probe.z + local_vd.z * step);
      s += step;
    }

    // The coarse floored stride quantizes the graze minimum, and the sampling
    // phase shifts every frame — corner coverage shimmers under motion. One
    // parabolic-interpolation step through the bracket tightens the minimum
    // (one extra distance eval, graze pixels only) and recovers a thin solid
    // chord the stride stepped over.
    if (min_behind < 2.0f * aa_width && bef_pd != FLT_MAX &&
        aft_pd != FLT_MAX) {
      float p = min_s - bef_s;
      float q = min_s - aft_s;
      float yb = bef_pd - min_behind;
      float ya = aft_pd - min_behind;
      float den = q * yb - p * ya;
      if (den < -1e-12f) {
        float ds = -0.5f * (q * q * yb - p * p * ya) / den;
        if (ds > -p && ds < -q) {
          Vector rp(min_pos.x + local_vd.x * ds, min_pos.y + local_vd.y * ds,
                    min_pos.z + local_vd.z * ds);
          float rpd = shape.distance(rp);
          if (rpd < min_behind) {
            min_behind = rpd;
            min_pos = rp;
          }
          if (min_behind < hit_threshold)
            return {true, min_pos, 0.0f};
        }
      }
    }

    float soft = (min_behind < aa_width)
                     ? quintic_kernel(1.0f - (min_behind - hit_threshold) /
                                                 (aa_width - hit_threshold))
                     : 0.0f;
    return {false, min_pos, soft};
  }

  /**
   * @brief Raymarches and shades a volume shape over its bounding sphere.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Shape Volume shape satisfying the concept below.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param bounds_center Bounding sphere center in physical LED space.
   * @param bounds_radius Bounding sphere radius in world units.
   * @param view_dir Normalized ray direction (camera → scene) in LED space.
   * @param shape Volume shape providing ray_to_local() and distance().
   * @param frag_fn Fragment shader invoked once per hit.
   * @param max_steps Maximum sphere-tracing steps per ray.
   * @param aa_width Anti-aliasing band half-width in world units.
   * @details Shape concept:
   *   std::pair<Vector, Vector> ray_to_local(const Vector &ro, const Vector
   *   &vd) const; Vector origin_to_local(const Vector &ro) const; float
   *   distance(const Vector &local_point) const;
   */
  template <int W, int H, typename Shape>
  static void
  draw(PipelineRef pipeline, Canvas &canvas, const Vector &bounds_center,
       float bounds_radius, const Vector &view_dir, const Shape &shape,
       FragmentShaderFn frag_fn, int max_steps = 15, float aa_width = 0.01f) {

    float vd_len = sqrtf(view_dir.x * view_dir.x + view_dir.y * view_dir.y +
                         view_dir.z * view_dir.z);
    float vd_inv = (vd_len > TOLERANCE) ? 1.0f / vd_len : 1.0f;
    Vector vd(view_dir.x * vd_inv, view_dir.y * vd_inv, view_dir.z * vd_inv);

    // Ray must start behind the farthest extent of the shape.
    float start_offset = 1.0f + bounds_radius;

    // bounds_center projected onto the view plane (⊥ vd).
    float bc_dot_vd = bounds_center.x * vd.x + bounds_center.y * vd.y +
                      bounds_center.z * vd.z;
    Vector bc_proj(bounds_center.x - bc_dot_vd * vd.x,
                   bounds_center.y - bc_dot_vd * vd.y,
                   bounds_center.z - bc_dot_vd * vd.z);
    float bounds_r2 = bounds_radius * bounds_radius;

    // Precompute the local-space view direction (shared across all pixels) and,
    // on the same cold transform, validate the volume preconditions. The per-step
    // early-out below compares a local-space dot against the world-space
    // bounds_radius, which holds only if ray_to_local is length-preserving
    // (|local_vd| == 1) and bounds_center maps to the shape's local origin (~0).
    // Trap a scaling shape or off-center bounds_center here, once per draw.
    auto [local_bc, local_vd] = shape.ray_to_local(bounds_center, vd);
    HS_CHECK(fabsf(local_vd.x * local_vd.x + local_vd.y * local_vd.y +
                   local_vd.z * local_vd.z - 1.0f) < TOLERANCE);
    HS_CHECK(local_bc.x * local_bc.x + local_bc.y * local_bc.y +
             local_bc.z * local_bc.z < TOLERANCE);
    // aa_width > 0 is the contract: the slow-path AA divides by (aa_width -
    // hit_threshold) == 0.9*aa_width, so a zero band-width gives 0/0 -> NaN.
    HS_CHECK(aa_width > 0.0f);

    BoundingSphere<W, H> bounds(bounds_center, bounds_radius);

    // Tier 2: Clamp volume bounds to clip region
    const auto &cr = canvas.clip();
    const auto vol_xc = cr.x_clip();
    int vol_y_lo =
        bounds.y_min > cr.render_y_start() ? bounds.y_min : cr.render_y_start();
    int vol_y_hi = bounds.y_max < cr.render_y_end() - 1 ? bounds.y_max
                                                        : cr.render_y_end() - 1;
    if (vol_y_lo > vol_y_hi)
      return;

    ScopedRenderTimer timer_guard(canvas);
    scan_region<W, H>(
        vol_y_lo, vol_y_hi,
        [&](int y, auto &&out) { return bounds.get_intervals(y, out); },
        [&](int, int, const Vector &p) {
          // Back-face cull
          float facing = p.x * vd.x + p.y * vd.y + p.z * vd.z;
          if (facing >= 0.0f)
            return;

          // Orthographic ray-sphere cull
          float pp_x = p.x - facing * vd.x;
          float pp_y = p.y - facing * vd.y;
          float pp_z = p.z - facing * vd.z;
          float dx = pp_x - bc_proj.x;
          float dy = pp_y - bc_proj.y;
          float dz = pp_z - bc_proj.z;
          if (dx * dx + dy * dy + dz * dz > bounds_r2)
            return;

          // Orthographic ray origin: outside the unit sphere
          Vector ro(pp_x - vd.x * start_offset, pp_y - vd.y * start_offset,
                    pp_z - vd.z * start_offset);

          // Transform the ray origin to local space once per pixel. The local
          // direction is constant across the draw (local_vd, computed above), so
          // only the origin is transformed here.
          Vector local_ro = shape.origin_to_local(ro);
          Vector closest_local;

          // --- Sphere tracing in local space ---
          float closest_d =
              trace_closest(shape, local_ro, local_vd, bounds_radius, max_steps,
                            aa_width, closest_local);

          if (closest_d >= aa_width)
            return;

          // --- Fragment shading ---
          Fragment frag;
          frag.pos = closest_local;
          frag.size = closest_d;
          frag_fn(closest_local, frag);

          // One-sided AA with quintic kernel
          float hit_threshold = aa_width * 0.1f;
          float edge_alpha;

          if (closest_d <= hit_threshold) {
            // FAST PATH: Solid hit. No probe needed.
            edge_alpha = 1.0f;
          } else {
            // SLOW PATH: fuzzy AA border. Standard one-sided AA coverage...
            edge_alpha = quintic_kernel(1.0f - (closest_d - hit_threshold) /
                                                   (aa_width - hit_threshold));

            // ...then probe behind the halo for a surface this edge occludes.
            Occluder occ = probe_occluder(shape, closest_local, local_vd,
                                          bounds_radius, hit_threshold, aa_width);
            if (occ.solid) {
              // Self-occlusion edge: antialias the foreground over the surface it
              // covers — lay the shaded background down, then blend the foreground
              // over it by the edge coverage. Smooth, vs. fading to black (fringe)
              // or snapping to opaque (jagged).
              Fragment bg;
              bg.pos = occ.behind;
              bg.size = 0.0f;
              frag_fn(occ.behind, bg);
              pipeline.plot(canvas, p, bg.color.color, 0.0f, bg.color.alpha);
              if (frag.color.alpha > 0.001f) {
                pipeline.plot(canvas, p, frag.color.color, 0.0f,
                              frag.color.alpha * edge_alpha);
              }
              return;
            }
            // No solid behind; a grazed background edge fills the corner,
            // shaded at its own point so the fill carries the background's
            // color, then the foreground blends over it.
            if (occ.soft > 0.001f) {
              Fragment bg;
              bg.pos = occ.behind;
              bg.size = 0.0f;
              frag_fn(occ.behind, bg);
              pipeline.plot(canvas, p, bg.color.color, 0.0f,
                            bg.color.alpha * occ.soft);
            }
          }

          if (frag.color.alpha * edge_alpha > 0.001f) {
            pipeline.plot(canvas, p, frag.color.color, 0.0f,
                          frag.color.alpha * edge_alpha);
          }
        },
        vol_xc);
  }
};

} // namespace Scan
