/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <utility>
#include "sdf.h"
#include "color.h"
#include "filter.h"
#include "static_circular_buffer.h"
#include "canvas.h"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif
#include "platform.h"

/**
 * @brief The Scan namespace contains volumetric (raster) drawing primitives.
 * @details General register mapping for Scan primitives:
 *  v0: Normalized parameter t (0-1) or angle
 *  v1: Raw distance or supplementary value
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
 * @details The shader and pipeline are intentionally type-erased
 * (FragmentShaderFn / PipelineRef), so the two per-pixel indirect calls below —
 * the shader at fragment_shader() and the plot at pipeline.plot() — do not
 * inline. This is a deliberate code-size tradeoff (one scanline instantiation
 * per <W,H> rather than per shader/pipeline combination); see the PipelineRef
 * note in concepts.h.
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
    frag_scratch.v2 = 0.0f;
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
                        PixelFn &&pixel_fn) {
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();

  // Interval scratch (~1.5 KiB) lives in scratch_arena_b, not the stack:
  // Phantasm's DTCM stack is tight and scan_region is on the deepest render
  // chain. Per-call bump scope; norm is cleared per row below.
  //
  // intervals holds a top-level shape's full per-row emission: a Union/SmoothUnion
  // merges both children (up to kIntervalSpanCap spans each) into one buffer, and
  // coalescing never grows the count, so 2*kIntervalSpanCap. norm holds one
  // seam-split per span, so 4*kIntervalSpanCap.
  ScratchScope scratch(scratch_arena_b);
  using IntervalBuf =
      StaticCircularBuffer<std::pair<float, float>, 2 * SDF::kIntervalSpanCap>;
  using NormBuf =
      StaticCircularBuffer<std::pair<float, float>, 4 * SDF::kIntervalSpanCap>;
  static_assert(IntervalBuf::kCapacity == SDF::MergedIntervalBuffer::kCapacity,
                "scan_region intervals must hold the largest top-level CSG "
                "emission (Subtract/Intersection: |A|+|B|+2)");
  static_assert(NormBuf::kCapacity == 2 * IntervalBuf::kCapacity,
                "norm must hold 2 spans per input interval (seam split)");
  auto &intervals = *new (scratch_arena_b.allocate(
      sizeof(IntervalBuf), alignof(IntervalBuf))) IntervalBuf();
  auto &norm = *new (scratch_arena_b.allocate(sizeof(NormBuf),
                                              alignof(NormBuf))) NormBuf();

  const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();

  // Inverted range (y_min > y_max) is a no-op: a disjoint CSG Intersection or a
  // fully-culled Face reports y_min=1, y_max=0, and the loop never runs.
  for (int y = y_min; y <= y_max; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];

    bool handled = get_intervals(
        y, [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });

    if (handled && !intervals.is_empty()) {
      // A shape spanning the full circle (len >= W) paints every column, so
      // detect it up front and skip the seam-split/sort/coalesce path.
      bool full_row = false;
      for (const auto &iv : intervals) {
        if (iv.second - iv.first >= static_cast<float>(W)) {
          full_row = true;
          break;
        }
      }

      if (full_row) {
        for (int x = 0; x < W; ++x)
          pixel_fn(x, y, Vector(sp * cos_theta[x], cp, sp * sin_theta[x]));
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
          for (int x = x1; x < x2; ++x)
            pixel_fn(x, y, Vector(sp * cos_theta[x], cp, sp * sin_theta[x]));
        }
      }
      intervals.clear();
    } else if (!handled) {
      for (int x = 0; x < W; ++x)
        pixel_fn(x, y, Vector(sp * cos_theta[x], cp, sp * sin_theta[x]));
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
    PixelCoords center_px = vector_to_pixel<W, H>(center);
    center_theta = center_px.x;
    // phi = acos(y) from the world y directly, not from center_px.y (avoids
    // compounding vector_to_pixel's phi_to_y rounding).
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
    // ray-sphere test rejects any extra column. min(W/2, ...) caps the span
    // length 2*x_half at W, scan_region's producer contract (interval length
    // <= W). Endpoints are not clamped to [0,W); scan_region wraps them.
    int x_half = std::min(W / 2, static_cast<int>(ceilf(theta_span)) + 1);
    out(center_theta - x_half, center_theta + x_half);
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
        if (xc.clipped(wx))
          return;
        process_pixel<W, H, ComputeUVs>(wx, y, p, pipeline, canvas, shape,
                                        fragment_shader, effective_debug,
                                        result_scratch, frag_scratch);
      });
}

/**
 * @brief Draws a ring whose radius is modulated around the circumference by
 *        shift_fn.
 */
struct DistortedRing {
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

    SDF::PlanarPolygon shape(res.first, thickness, sides, phase);
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
    SDF::Star shape(res.first, res.second, sides, phase);
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
    SDF::Flower shape(res.first, res.second, sides, phase);
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

    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset);
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
   */
  template <int W, int H, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   bool debug_bb = false) {
    ScratchScope scope(scratch_arena);
    auto *scratch =
        static_cast<SDF::FaceScratchBuffer *>(scratch_arena.allocate(
            sizeof(SDF::FaceScratchBuffer), alignof(SDF::FaceScratchBuffer)));

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    const uint16_t *fo = mesh.get_face_offsets_data();
    size_t fi_size = mesh.get_faces_size();

    for (size_t i = 0; i < num_f; ++i) {
      size_t count = fc[i];

      // Trap malformed mesh data: an offset/count pair disagreeing with the flat
      // index array yields an out-of-bounds span for SDF::Face. Cold per-face check.
      HS_CHECK(static_cast<size_t>(fo[i]) + count <= fi_size,
               "mesh face span exceeds face index array");

      std::span<const Vector> verts(mesh.vertices.data(), mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      // The IIFE lets the HS_PROFILE scope close right after the prvalue Face is
      // constructed in place (guaranteed elision), so it measures construction
      // alone, not the rasterize below.
      SDF::Face shape = [&] {
        HS_PROFILE(scan_face_setup);
        return SDF::Face(verts, indices, 0.0f, *scratch, H + hs::H_OFFSET, H);
      }();

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
   * @brief Compile-time table of sub-pixel sample offsets.
   * @tparam SAMPLES Number of sub-pixel samples per pixel.
   * @details Holds a centered 2×2 grid (±0.25, ±0.25), selected by the low two
   * bits of the sample index.
   */
  template <int SAMPLES> struct SampleOffsets {
    float x[SAMPLES]; /**< Per-sample x offsets in pixel units. */
    float y[SAMPLES]; /**< Per-sample y offsets in pixel units. */
  };
  /**
   * @brief Builds the sub-pixel sample offset table at compile time.
   * @tparam SAMPLES Number of sub-pixel samples per pixel.
   * @return Offset table with each sample on a centered 2×2 grid (±0.25, ±0.25).
   * @details Fully compile-time. The grid is centered at ±0.25, not the ±0.5
   * pixel corners, so the four samples lie strictly inside the pixel (corner
   * samples sit on shared boundaries, correlating adjacent pixels).
   */
  template <int SAMPLES>
  static constexpr SampleOffsets<SAMPLES> make_sample_offsets() {
    constexpr float eps = 0.25f;
    SampleOffsets<SAMPLES> o{};
    for (int i = 0; i < SAMPLES; ++i) {
      int qx = (i & 1) ? -1 : 1;
      int qy = (i & 2) ? -1 : 1;
      o.x[i] = eps * qx;
      o.y[i] = eps * qy;
    }
    return o;
  }

  /**
   * @brief Projects a sub-pixel coordinate to its world-space unit vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param px Sub-pixel column (pixel index + corner offset).
   * @param py Sub-pixel row (pixel index + corner offset).
   * @return World-space unit vector for the sub-pixel sample.
   * @details Uses explicit theta/phi trig (matching the non-SSAA
   * pixel_to_vector path so SSAA output stays byte-for-byte stable).
   */
  template <int W, int H>
  static inline Vector ssaa_sample_vector(float px, float py) {
    constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
    constexpr float w_float = static_cast<float>(W);
    float theta = (px * 2.0f * PI_F) / w_float;
    float phi = (py * PI_F) / h_virt_minus_1;
    float sin_phi = sinf(phi);
    return Vector(sin_phi * cosf(theta), cosf(phi), sin_phi * sinf(theta));
  }

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

    if constexpr (SAMPLES == 1) {
      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = cr.x_start; x < cr.x_end; ++x) {
          Vector v = pixel_to_vector<W, H>(x, y);
          Color4 sample = shader(v);
          canvas(x, y) = sample.color * sample.alpha;
        }
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      constexpr auto offsets = make_sample_offsets<SAMPLES>();

      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = cr.x_start; x < cr.x_end; ++x) {
          // Premultiplied SSAA: accumulate each sample's coverage-weighted color
          // (color * alpha * 1/N) and write it directly, matching the SAMPLES==1
          // path (sample.color * sample.alpha).
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            Vector v = ssaa_sample_vector<W, H>(
                static_cast<float>(x) + offsets.x[i],
                static_cast<float>(y) + offsets.y[i]);

            Color4 sample = shader(v);
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
   * @details Separates expensive per-pixel work (vertex_shader, called once at
   * pixel center) from per-sub-sample evaluation (fragment_shader, called
   * SAMPLES×).
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
    // frag_base is per pixel, not per draw: each pixel starts from a default
    // Fragment, so a vertex shader writing only some registers (v0-v3/size/age/
    // color) can't inherit the previous pixel's values.
    if constexpr (SAMPLES == 1) {
      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = cr.x_start; x < cr.x_end; ++x) {
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
      constexpr auto offsets = make_sample_offsets<SAMPLES>();

      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      ScopedRenderTimer timer_guard(canvas);
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = cr.x_start; x < cr.x_end; ++x) {
          Vector center_v = pixel_to_vector<W, H>(x, y);

          Fragment frag_base;
          frag_base.pos = center_v;
          vertex_shader(frag_base);

          // Premultiplied SSAA: accumulate coverage-weighted color directly (see
          // the single-callback overload).
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            Vector v = ssaa_sample_vector<W, H>(
                static_cast<float>(x) + offsets.x[i],
                static_cast<float>(y) + offsets.y[i]);

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
   * @details The per-pixel march needs a fresh local origin every pixel but the
   * local direction is constant across the draw (vd is fixed), so the volume
   * loop precomputes it once via ray_to_local() and calls this to skip the
   * redundant per-pixel direction rotation.
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
   * @brief Sphere-traces a ray in local space, recording the closest approach.
   * @param shape Volume shape providing distance().
   * @param local_ro Ray origin in local space.
   * @param local_vd Unit ray direction in local space.
   * @param bounds_radius Bounding sphere radius (past-the-back early-out).
   * @param max_steps Maximum sphere-tracing steps.
   * @param aa_width Anti-aliasing band half-width (deep-hit early-out).
   * @param closest_local Output: local-space point of closest approach.
   * @return Signed distance at the closest approach (FLT_MAX if never sampled).
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
    Vector behind; /**< Its local-space hit point (only when solid). */
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
   * @return An Occluder: a solid hit point to antialias the edge over, or a soft
   *         coverage for the corner where two edges meet.
   */
  template <typename Shape>
  static __attribute__((always_inline)) Occluder
  probe_occluder(const Shape &shape, const Vector &closest_local,
                 const Vector &local_vd, float bounds_radius, float hit_threshold,
                 float aa_width) {
    // Seed from the closest approach (not the march's terminal local_p — the loop
    // can exit having stepped past the bounding sphere) and march forward for a
    // surface this foreground halo occludes. A solid hit means the halo is a
    // self-occlusion edge, not a true silhouette, so the caller antialiases it over
    // that surface instead of fading it to black.
    //
    // closest_local sits near the FRONT of the volume, so a background feature can
    // lie up to a full diameter (2*bounds_radius) deeper; the probe must reach that
    // far or it falsely reports "no surface behind". The step is floored to punch
    // past the foreground surface the main trace stalled on (a pure sphere trace
    // would crawl its tangent halo): fine in the near field so a background just
    // behind the foreground (two features first crossing) isn't jumped over, coarse
    // beyond it so the budget still reaches a deep back face. Termination is the
    // back face of the bounding sphere (matching trace_closest's early-out).
    //
    // With no solid hit, a background edge can still graze in the AA band at the
    // corner where two edges meet (never pd < hit_threshold); detect it as a local
    // minimum of pd (rise then fall, robust when the gap is too thin for pd to
    // clear the band) and report its coverage for the corner fill.
    Vector probe = closest_local;
    float prev = FLT_MAX;  // previous step's distance
    bool climbing = false; // pd has risen off the foreground graze
    float min_behind = FLT_MAX;
    for (int i = 0; i < 24; ++i) {
      // Stop at the back of the bounding sphere: nothing left to occlude this halo.
      if (probe.x * local_vd.x + probe.y * local_vd.y + probe.z * local_vd.z >
          bounds_radius)
        break;
      float pd = shape.distance(probe);
      if (pd < hit_threshold)
        return {true, probe, 0.0f}; // solid surface behind the edge
      if (pd > prev)
        climbing = true; // moving away from the foreground graze
      else if (climbing)
        min_behind = std::min(min_behind, pd); // descending toward a surface behind
      prev = pd;
      float floor = bounds_radius * (i < 6 ? 0.04f : 0.12f);
      float step = std::max(pd * 0.9f, floor);
      probe = Vector(probe.x + local_vd.x * step, probe.y + local_vd.y * step,
                     probe.z + local_vd.z * step);
    }

    float soft = (min_behind < aa_width)
                     ? quintic_kernel(1.0f - (min_behind - hit_threshold) /
                                                 (aa_width - hit_threshold))
                     : 0.0f;
    return {false, closest_local, soft};
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
        [&](int wx, int, const Vector &p) {
          if (vol_xc.clipped(wx))
            return;
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
              if (frag.color.alpha > 0.001f) {
                Fragment bg;
                bg.pos = occ.behind;
                bg.size = 0.0f;
                frag_fn(occ.behind, bg);
                pipeline.plot(canvas, p, bg.color.color, 0.0f, bg.color.alpha);
                pipeline.plot(canvas, p, frag.color.color, 0.0f,
                              frag.color.alpha * edge_alpha);
              }
              return;
            }
            // No solid behind; a grazed background edge fills the corner ("over").
            edge_alpha += (1.0f - edge_alpha) * occ.soft;
          }

          if (frag.color.alpha * edge_alpha > 0.001f) {
            pipeline.plot(canvas, p, frag.color.color, 0.0f,
                          frag.color.alpha * edge_alpha);
          }
        });
  }
};

} // namespace Scan
