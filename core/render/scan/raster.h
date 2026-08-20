/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <type_traits>
#include <utility>
#include "render/sdf.h"
#include "render/shading.h"
#include "mesh/mesh_class_types.h"
#include "color/color.h"
#include "render/filter.h"
#include "containers/static_circular_buffer.h"
#include "render/canvas.h"
#include "platform/platform.h"

/**
 * @file raster.h
 * @brief The shared scanline kernel: pole LOD, the per-pixel body, interval and
 * row-span assembly, scan_region and the rasterize/rasterize_solid entry points.
 */
/**
 * @brief The Scan namespace contains area-filling scanline drawing primitives.
 * @details Pipeline typing follows the per-pixel body. A body that already
 * calls a type-erased FragmentShaderFn pays one indirect call per pixel
 * whatever the sink is, so those entry points take the sink as PipelineRef —
 * either by value or as a PipelineT defaulted to it — and keep the single
 * instantiation per <W,H> that erasure buys (see the PipelineRef note in
 * concepts.h). The bodies with no erased call — the constant-color
 * rasterize_solid and draw_solid family, the fused RingGroup and
 * DistortedRingStack walks, rasterize_face and Mesh — name PipelineT with no
 * default, so the caller's sink types the whole body.
 */
namespace Scan {

/** @brief Coverage/alpha at or below which a fragment is dropped unplotted. */
inline constexpr float MIN_ALPHA = 0.001f;

/**
 * @brief Columns one shade may cover at a row's colatitude.
 * @param sin_phi Sine of the row's colatitude; sign is ignored.
 * @return Run length in columns, 1 when decimation is off or unwarranted.
 * @details `pole_lod_aggressiveness / sin(phi)`, clamped to POLE_LOD_MAX_RUN.
 *          Returns 1 at aggressiveness 0, which makes every caller's walk
 *          bit-identical to an undecimated one.
 */
inline int pole_lod_run(float sin_phi) {
  const float lod = pole_lod_aggressiveness;
  if (lod <= 0.0f)
    return 1;
  const float a = sin_phi < 0.0f ? -sin_phi : sin_phi;
  if (a < 1e-6f)
    return POLE_LOD_MAX_RUN;
  const int run = static_cast<int>(lod / a);
  if (run <= 1)
    return 1;
  return run < POLE_LOD_MAX_RUN ? run : POLE_LOD_MAX_RUN;
}

/**
 * @brief Clearance beyond a walk's own threshold at which one probe vouches for
 *        a whole pole-LOD block.
 * @tparam W Canvas width in pixels.
 * @tparam ShapeT Shape the probe read distance() from.
 * @param run Columns in the block; 1 yields no slack.
 * @param sin_phi Sine of the row's colatitude.
 * @return Slack in the same units the walk's distance() reports.
 * @details Great-circle arc from a block's first column to its last -- a block
 *          of longitude, foreshortened by sin(phi) -- times the shape's own
 *          change-per-arc factor (SDF::arc_stretch). A probe farther than this
 *          from the surface cannot change side anywhere in the block. A walk
 *          whose reporting stretches further still (a gnomonic projection
 *          stretches by 1 + r^2) scales the result by its own factor.
 */
template <int W, typename ShapeT>
__attribute__((always_inline)) inline float pole_lod_slack(int run,
                                                           float sin_phi) {
  return static_cast<float>(run - 1) * (2.0f * PI_F / W) * sin_phi *
         SDF::arc_stretch<std::remove_cvref_t<ShapeT>>;
}

/**
 * @brief Factor from an angular step to the units a shape's distance() reports.
 * @return 1 for a shape reporting angular distance.
 * @details A block probe's slack is an arc, so a shape reporting in another
 * unit has to scale it. SDF::Face reports gnomonic-plane distance, which
 * stretches an angular step by up to 1 + r^2; max_dist bounds r over every
 * probe its cull admits.
 */
__attribute__((always_inline)) inline float report_stretch(const auto &) {
  return 1.0f;
}
__attribute__((always_inline)) inline float
report_stretch(const SDF::Face &shape) {
  return 1.0f + shape.max_dist_sq;
}

/**
 * @brief Whether a walk over a shape may settle a whole pole-LOD block from one
 *        probe of it.
 * @tparam ShapeT Shape the probe reads distance() from.
 * @details False with the knob compiled out, and false for a shape whose
 * distance() has no finite change-per-arc factor: no slack makes one probe
 * vouch for its neighbours, so the walk stays per column.
 */
template <typename ShapeT>
inline constexpr bool pole_lod_blocks =
    POLE_LOD_ENABLED &&
    SDF::arc_stretch<std::remove_cvref_t<ShapeT>> < SDF::ARC_STRETCH_UNBOUNDED;

/**
 * @brief Whether a probe's clearance report bounds a whole pole-LOD block.
 * @tparam ShapeT Shape the probe read distance() from.
 * @param threshold Clearance the walk tests the probe against.
 * @param block_slack Extra clearance demanded of a block probe.
 * @return True when a report at or above threshold + block_slack holds for
 *         every column in the block.
 * @details Past its reject band a shape may report FAR_SENTINEL instead of a
 * distance (SDF::reject_margin), and the sentinel bounds nothing. The block
 * test is trustworthy only where the widened threshold sits inside the margin,
 * so that a sentinel implies the surface is genuinely that far away.
 */
template <typename ShapeT>
__attribute__((always_inline)) inline bool
probe_bounds_block(float threshold, float block_slack) {
  return threshold + block_slack <=
         SDF::reject_margin<std::remove_cvref_t<ShapeT>>;
}

/**
 * @brief Extra clearance a block probe must carry, in the shape's report units.
 * @tparam W Canvas width in columns.
 * @tparam ShapeT Shape the probe reads distance() from.
 * @param run Columns in the block; 1 yields no slack.
 * @param p_y Probe's y coordinate on the unit sphere (cos of the colatitude).
 * @param shape Shape the probe reads distance() from.
 */
template <int W, typename ShapeT>
__attribute__((always_inline)) inline float
pole_lod_block_slack(int run, float p_y, const ShapeT &shape) {
  const float sin_phi = sqrtf(std::max(0.0f, 1.0f - p_y * p_y));
  return pole_lod_slack<W, ShapeT>(run, sin_phi) * report_stretch(shape);
}

/**
 * @brief Whether one probe settles a whole pole-LOD block on one side.
 * @tparam ShapeT Shape the probe read distance() from.
 * @param clearance Probe's distance from the surface on the side under test;
 *        negate a solid's report to test the interior side.
 * @param threshold Clearance the walk tests the probe against.
 * @param block_slack Extra clearance demanded of a block probe.
 * @return True when the probe holds for every column in the block.
 */
template <typename ShapeT>
__attribute__((always_inline)) inline bool
pole_lod_block_settles(float clearance, float threshold, float block_slack) {
  return clearance >= threshold + block_slack &&
         probe_bounds_block<ShapeT>(threshold, block_slack);
}

/**
 * @brief Validates that the canvas is the size the scan was instantiated for.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param canvas Destination canvas.
 * @details Checked once per draw, not per row or pixel. A narrower canvas
 * indexes the phi LUT past its rows and plots through a pipeline whose wrap
 * period disagrees with the framebuffer stride.
 */
template <int W, int H>
HS_NOINLINE_NOCLONE inline void check_canvas_dims(const Canvas &canvas) {
  HS_CHECK(canvas.width() == W && canvas.height() == H);
}

/**
 * @brief Validates that a direct-raster pipeline was prepared for this canvas.
 * @param pipeline Plotting pipeline the draw writes through.
 * @param canvas Destination canvas.
 * @details A direct-raster sink writes through a cached framebuffer base; the
 * canvas double-buffers, so a stale base is the buffer the display is scanning
 * out. Compiles away for a sink without the hook.
 */
template <typename PipelineT>
inline void check_pipeline_prepared(PipelineT &pipeline, Canvas &canvas) {
  if constexpr (requires { pipeline.prepared_for(canvas); })
    HS_CHECK(pipeline.prepared_for(canvas),
             "direct raster pipeline not prepared for this canvas");
}

/**
 * @brief Validates that a type-erased fragment shader refers to a callable.
 * @param fragment_shader Shader the draw invokes per pixel.
 * @details Checked once per draw, not per pixel: FunctionRef::operator() guards
 * an empty ref with assert alone, so an optimized build calls a null thunk.
 */
HS_NOINLINE_NOCLONE inline void
check_fragment_shader(FragmentShaderFn fragment_shader) {
  HS_CHECK(fragment_shader, "Scan: fragment_shader must be non-null");
}

/**
 * @brief Anti-aliased coverage of a solid shape at a signed distance.
 * @param d Signed distance to the surface, negative inside.
 * @param pixel_width Angular width of one column, the AA half-reach.
 * @return Coverage in [0, 1]; 1 at or inside a full pixel of depth.
 * @details The single spelling of the solid AA ramp, shared by every scan loop.
 */
__attribute__((always_inline)) inline float solid_coverage(float d,
                                                           float pixel_width) {
  if (d <= -pixel_width)
    return 1.0f;
  return quintic_kernel(0.5f - d / (2.0f * pixel_width));
}

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
 * @param max_run Columns [x, x+max_run) offered by the scan; see scan_region.
 * @return Columns consumed: max_run when the probe at x decided the whole
 *         offer (all clear, or a solid interior shaded once and splatted),
 *         else 1.
 * @details The shader is type-erased (FragmentShaderFn) and the pipeline
 * defaults to the erased PipelineRef, so the scanline instantiates once per
 * <W,H> rather than per shader/pipeline; a caller holding a typed sink names it
 * as PipelineT. See the PipelineRef note in concepts.h.
 */
template <int W, int H, bool ComputeUVs = true,
          typename PipelineT = PipelineRef>
inline int process_pixel(int x, int y, const Vector &p, PipelineT &pipeline,
                         Canvas &canvas, const auto &shape,
                         FragmentShaderFn fragment_shader, bool debug_bb,
                         SDF::DistanceResult &result_scratch,
                         Fragment &frag_scratch, int max_run = 1) {
  shape.template distance<ComputeUVs>(p, result_scratch);

  float d = result_scratch.dist;
  constexpr float pixel_width = 2.0f * PI_F / W;
  constexpr bool solid = std::remove_cvref_t<decltype(shape)>::is_solid;
  float threshold = solid ? pixel_width : 0.0f;

  // Strokes never splat: their coverage varies inside the band, so only the
  // all-clear case consumes the block.
  int span = 1;
  if constexpr (pole_lod_blocks<decltype(shape)>) {
    if (max_run > 1 && !debug_bb) {
      const float block_slack = pole_lod_block_slack<W>(max_run, p.y, shape);
      if (pole_lod_block_settles<decltype(shape)>(d, threshold, block_slack))
        return max_run;
      // A sentineled subtrahend loses Subtract's max, so an inside report is
      // bounded no further than a clear one.
      if (solid &&
          pole_lod_block_settles<decltype(shape)>(-d, pixel_width, block_slack))
        span = max_run;
    }
  }

  if (debug_bb || d < threshold) {
    float alpha = 1.0f;

    if constexpr (solid) {
      alpha = solid_coverage(d, pixel_width);
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

    if (!debug_bb && alpha <= MIN_ALPHA)
      return span;

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

    if (frag_scratch.color.alpha > MIN_ALPHA) {
      const float a = frag_scratch.color.alpha * alpha;
      for (int i = 0; i < span; ++i)
        pipeline.plot(canvas, x + i, y, frag_scratch.color.color,
                      frag_scratch.age, a);
    }
  }
  return span;
}

/**
 * @brief Wraps, sorts and coalesces a row's column spans into integer runs.
 * @tparam W Canvas width in pixels.
 * @tparam IntervalBufT Source span buffer type.
 * @tparam NormBufT Seam-split scratch buffer type; holds 2 spans per source span.
 * @param intervals Source spans in fractional column units, each of length <= W.
 * @param norm Scratch, cleared here, receiving coalesced integer column runs.
 * @details Wrapping each start into [0, W) and splitting spans that cross the
 * x=0 seam gives the forward sweep sorted, non-wrapping input even when a
 * shape/CSG straddles θ=0. Runs are emitted monotone and disjoint: last_x2
 * clamps a run's start past the previous run's end so two spans sharing a
 * fractional column do not both paint it (double shade / alpha). The clamp can
 * leave a run empty (x1 == x2); consumers iterate [x1, x2) and skip it.
 */
template <int W, typename IntervalBufT, typename NormBufT>
HS_NOINLINE_NOCLONE inline void coalesce_spans(const IntervalBufT &intervals,
                                               NormBufT &norm) {
  norm.clear();
  SDF::normalize_intervals_to_range<W>(intervals, norm);
  SDF::sort_intervals_by_start(norm);

  float current_end = -FLT_MAX;
  int last_x2 = 0;
  size_t write = 0;
  const size_t count = norm.size();
  for (size_t read = 0; read < count; ++read) {
    const auto iv = norm[read];
    if (iv.second <= current_end)
      continue;
    float start = std::max(iv.first, current_end);
    float end = iv.second;
    current_end = end;

    int x1 = static_cast<int>(floorf(start));
    int x2 = static_cast<int>(ceilf(end));
    // A zero-width interval (x1 == x2) still owns its pixel column; widen to
    // paint it.
    if (x1 == x2)
      x2++;
    // After the wrap/split start>=0 and end<=W, so x2>W only from the widen at
    // the right edge; x1<0 is a defensive floor.
    if (x1 < 0)
      x1 = 0;
    if (x2 > W)
      x2 = W;
    if (x1 < last_x2)
      x1 = last_x2;
    last_x2 = x2;
    norm[write++] = {static_cast<float>(x1), static_cast<float>(x2)};
  }
  while (norm.size() > write)
    norm.pop_back();
}

/**
 * @brief Splits a column run by the clip arc and emits the surviving pieces.
 * @tparam EmitFn Callable (int x1, int x2); an empty piece is its own no-op.
 * @param x1 First column of the run.
 * @param x2 One past the run's last column.
 * @param xc Column-arc clip.
 * @param emit Sink receiving each surviving piece.
 * @details A wrapping arc is [rs, W) ∪ [0, re); the two pieces are disjoint
 * (re <= rs) and emitted low piece first, so no column is emitted twice and
 * the pieces keep ascending column order.
 */
template <typename EmitFn>
__attribute__((always_inline)) inline void
clip_run(int x1, int x2, ClipRegion::XClip xc, EmitFn &&emit) {
  if (!xc.active) {
    emit(x1, x2);
  } else if (xc.wrap) {
    emit(x1, std::min(x2, xc.re));
    emit(std::max(x1, xc.rs), x2);
  } else {
    emit(std::max(x1, xc.rs), std::min(x2, xc.re));
  }
}

/**
 * @brief Walks the canvas columns a clip arc admits, skipping the rest whole.
 * @tparam W Canvas width in pixels.
 * @tparam BodyFn Callable (int x).
 * @param xc Column-arc clip.
 * @param body Sink receiving each surviving column, in the same ascending order
 *        a per-column `clipped()` test would leave.
 * @details The arc's pieces are collected before they are walked, so `body`
 * inlines once rather than once per clip_run emission.
 */
template <int W, typename BodyFn>
__attribute__((always_inline)) inline void
walk_clip_columns(ClipRegion::XClip xc, BodyFn &&body) {
  int lo[2];
  int hi[2];
  int pieces = 0;
  clip_run(0, W, xc, [&](int x1, int x2) {
    if (x1 >= x2)
      return;
    lo[pieces] = x1;
    hi[pieces] = x2;
    ++pieces;
  });
  for (int i = 0; i < pieces; ++i)
    for (int x = lo[i]; x < hi[i]; ++x)
      body(x);
}

/**
 * @brief Turns one row's emitted spans into clipped integer column runs.
 * @tparam W Canvas width in pixels.
 * @tparam IntervalBufT Source span buffer type.
 * @tparam NormBufT Seam-split scratch buffer; holds 2 spans per source span.
 * @tparam EmitFn Callable (int x1, int x2) receiving each run.
 * @param handled False when the producer declined the row, which paints it
 *                whole; the spans are then ignored.
 * @param intervals The row's emitted spans, in fractional column units.
 * @param norm Scratch for the seam split; clobbered.
 * @param xc Column-arc clip applied to every run before it is emitted.
 * @param emit Sink receiving the runs, ascending and disjoint.
 * @details Shared by scan_region and by the fused RingGroup walk, which cannot
 * call scan_region itself. Runs arrive in ascending column order, so a sink may
 * treat them as one left-to-right sweep.
 */
template <int W, typename IntervalBufT, typename NormBufT, typename EmitFn>
__attribute__((always_inline)) inline void
emit_row_runs(bool handled, const IntervalBufT &intervals, NormBufT &norm,
              ClipRegion::XClip xc, EmitFn &&emit) {
  if (!handled) {
    clip_run(0, W, xc, emit);
    return;
  }
  if (intervals.is_empty())
    return;

  // A single span covering the full circle (len >= W) paints every column;
  // detect it up front and skip the seam-split/sort/coalesce path. Coverage
  // assembled from multiple abutting spans is not caught here — it falls to the
  // slow path, which still paints every covered column.
  for (const auto &iv : intervals) {
    if (iv.second - iv.first >= static_cast<float>(W)) {
      clip_run(0, W, xc, emit);
      return;
    }
  }

  coalesce_spans<W>(intervals, norm);
  for (const auto &run : norm)
    clip_run(static_cast<int>(run.first), static_cast<int>(run.second), xc,
             emit);
}

/** Capacity of scan_region's per-row emission buffer, and so the compile-time
 *  ceiling on sdf_max_spans of any shape handed to a rasterizer. Covers a
 *  top-level Union/SmoothUnion (|A|+|B|) and Subtract (2·|A|); a top-level
 *  Intersection (2·|A| + 2·|B|) fits only while both children stay under
 *  INTERVAL_SPAN_CAP/2 spans. */
inline constexpr size_t TOP_SPAN_CAP = 2 * SDF::INTERVAL_SPAN_CAP + 2;

/** Traps at compile time when a top-level shape can emit more spans per row
 *  than scan_region's buffer holds. */
template <typename ShapeT>
inline constexpr bool fits_top_span_cap =
    SDF::sdf_max_spans<std::remove_cvref_t<ShapeT>>::value <= TOP_SPAN_CAP;

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
 * @param pixel_fn (int wx, int y, const Vector &p, int max_run) -> int, offered
 *                 the columns [wx, wx+max_run) with p at column wx; returns how
 *                 many of them it consumed — max_run when the wx probe decided
 *                 the whole offer, else 1 (never 0).
 * @param xc Column-arc clip; pixel runs are intersected with the arc before
 *           walking, so pixel_fn is never called for a clipped column.
 * @details Iterates y in [y_min, y_max], collects float intervals per row via
 * get_intervals, wraps x coordinates, and offers pixel_fn column runs.
 *
 * Near-pole rows offer whole blocks of `pole_lod_aggressiveness / sin(phi)`
 * columns (constants.h) so the sink can settle physically-overlapping columns
 * with one probe; the sink keeps per-column resolution wherever its probe
 * cannot vouch for the block. Only full canvas-aligned blocks are offered, so an
 * offer never straddles two blocks and a settled column always takes its shade
 * from its own block's anchor. A block a clip or span edge truncates goes per
 * column instead, so the columns beside a segment seam shade at full resolution
 * rather than from the anchor the neighbouring segment would have used. At
 * aggressiveness 0 every offer is one column and the walk is bit-identical to an
 * undecimated one.
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
  // intervals holds a top-level shape's full per-row emission; norm holds one
  // seam-split per span, so 2x.
  ScratchScope scratch(scratch_arena_b);
  using IntervalBuf = StaticCircularBuffer<SDF::Interval, TOP_SPAN_CAP>;
  using NormBuf = StaticCircularBuffer<SDF::Interval, 2 * TOP_SPAN_CAP>;
  static_assert(IntervalBuf::CAPACITY >= SDF::MergedIntervalBuffer::CAPACITY,
                "scan_region intervals must hold a top-level Union's merged "
                "emission (|A|+|B|)");
  static_assert(NormBuf::CAPACITY == 2 * IntervalBuf::CAPACITY,
                "norm must hold 2 spans per input interval (seam split)");
  auto &intervals = *scratch_arena_b.make<IntervalBuf>();
  auto &norm = *scratch_arena_b.make<NormBuf>();

  const float *cos_theta =
      TrigLUT<W, H>::sin_theta.data() + W / 4; // cos via +W/4
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();

  // A full canvas-aligned block of `stride` columns is offered to the sink as
  // one call probed at the block's first column; the sink returns how many
  // columns it consumed. A block the run truncates is walked per column, so a
  // settled column is always settled from its own block's anchor.
  auto walk = [&](int x1, int x2, int y, float sp, float cp, int stride) {
    for (int x = x1; x < x2; ++x) {
      if constexpr (POLE_LOD_ENABLED) {
        if (stride > 1 && x % stride == 0 && x + stride <= x2) {
          x += pixel_fn(x, y, Vector(sp * cos_theta[x], cp, sp * sin_theta[x]),
                        stride) -
               1;
          continue;
        }
      }
      pixel_fn(x, y, Vector(sp * cos_theta[x], cp, sp * sin_theta[x]), 1);
    }
  };

  // Inverted range (y_min > y_max) is a no-op: a disjoint CSG Intersection or a
  // fully-culled Face reports y_min=1, y_max=0, and the loop never runs.
  for (int y = y_min; y <= y_max; ++y) {
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    const int stride = pole_lod_run(sp);

    // Clear before the producer runs: a producer that emits and then returns
    // false must not leak spans into the next row.
    intervals.clear();
    bool handled = get_intervals(
        y, [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });

    emit_row_runs<W>(handled, intervals, norm, xc,
                     [&](int x1, int x2) { walk(x1, x2, y, sp, cp, stride); });
  }
}

/**
 * @brief Computes bounding sphere y-range and per-row x intervals.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> struct BoundingSphere {
  int y_min, y_max;
  float center_theta;   /**< Longitude of center in pixel units. */
  float cos_rho;        /**< cos of the cap's angular radius. */
  float cos_center_phi; /**< cos of the cap center's colatitude. */
  float sin_center_phi; /**< sin of the cap center's colatitude. */

  /**
   * @brief Constructs the bounding sphere clamped to the canvas rows.
   * @param center World-space unit vector at the sphere center.
   * @param bounds_radius Bounding radius in world units (sin of angular extent).
   */
  BoundingSphere(const Vector &center, float bounds_radius) {
    float angular_radius = asinf(std::min(bounds_radius, 1.0f));
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
 * @brief Intersects a shape's row band with the clip region's render rows.
 * @tparam BoundsT Bounds type exposing `y_min` and `y_max`.
 * @param bounds Shape vertical bounds.
 * @param cr Clip region of the destination canvas.
 * @param y_lo Receives the first row to scan, inclusive.
 * @param y_hi Receives the last row to scan, inclusive.
 * @return False when the intersection is empty and nothing should be scanned.
 */
template <typename BoundsT>
__attribute__((always_inline)) inline bool
clamp_rows_to_clip(const BoundsT &bounds, const ClipRegion &cr, int &y_lo,
                   int &y_hi) {
  y_lo =
      bounds.y_min > cr.render_y_start() ? bounds.y_min : cr.render_y_start();
  y_hi = bounds.y_max < cr.render_y_end() - 1 ? bounds.y_max
                                              : cr.render_y_end() - 1;
  return y_lo <= y_hi;
}

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
  static_assert(SDF::ScanShape<std::remove_cvref_t<decltype(shape)>, W, H>,
                "Scan::rasterize shape must expose is_solid, "
                "get_vertical_bounds<H>(), get_horizontal_intervals<W, H>() "
                "and distance<ComputeUVs>()");
  static_assert(
      fits_top_span_cap<decltype(shape)>,
      "top-level shape can emit more spans per row than scan_region's "
      "buffer holds; flatten the composition or raise "
      "INTERVAL_SPAN_CAP");

  check_canvas_dims<W, H>(canvas);
  check_pipeline_prepared(pipeline, canvas);
  check_fragment_shader(fragment_shader);
  bool effective_debug = debug_bb || canvas.debug();

  int y_lo, y_hi;
  const auto &cr = canvas.clip();
  const auto xc = cr.x_clip();
  auto bounds = shape.template get_vertical_bounds<H>();
  if (!clamp_rows_to_clip(bounds, cr, y_lo, y_hi))
    return;

  SDF::DistanceResult result_scratch;
  Fragment frag_scratch;

  scan_region<W, H>(
      y_lo, y_hi,
      [&](int y, auto &&out) {
        return shape.template get_horizontal_intervals<W, H>(y, out);
      },
      [&](int wx, int y, const Vector &p, int max_run) {
        return process_pixel<W, H, ComputeUVs>(
            wx, y, p, pipeline, canvas, shape, fragment_shader, effective_debug,
            result_scratch, frag_scratch, max_run);
      },
      xc);
}

/**
 * @brief Rasterizes a solid SDF with a constant color and typed pipeline.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam SineDistance Use the shape's sine-domain angular distance.
 * @tparam PipelineT Plotting pipeline type.
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Solid SDF shape providing bounds, intervals, and distance().
 * @param color Constant source color and alpha.
 * @param debug_bb When true, renders the bounding box for debugging.
 */
template <int W, int H, bool SineDistance = false, typename PipelineT>
HS_NOINLINE_NOCLONE inline void
rasterize_solid(PipelineT &pipeline, Canvas &canvas, const auto &shape,
                const Color4 &color, bool debug_bb = false) {
  static_assert(std::remove_cvref_t<decltype(shape)>::is_solid);
  static_assert(SDF::ScanShape<std::remove_cvref_t<decltype(shape)>, W, H>,
                "Scan::rasterize_solid shape must expose is_solid, "
                "get_vertical_bounds<H>(), get_horizontal_intervals<W, H>() "
                "and distance<ComputeUVs>()");
  static_assert(
      fits_top_span_cap<decltype(shape)>,
      "top-level shape can emit more spans per row than scan_region's "
      "buffer holds; flatten the composition or raise "
      "INTERVAL_SPAN_CAP");

  check_canvas_dims<W, H>(canvas);
  check_pipeline_prepared(pipeline, canvas);
  bool effective_debug = debug_bb || canvas.debug();

  int y_lo, y_hi;
  const auto &cr = canvas.clip();
  const auto xc = cr.x_clip();
  auto bounds = shape.template get_vertical_bounds<H>();
  if (!clamp_rows_to_clip(bounds, cr, y_lo, y_hi) ||
      (!effective_debug && color.alpha <= MIN_ALPHA))
    return;

  SDF::DistanceResult result;
  constexpr float PIXEL_WIDTH = 2.0f * PI_F / W;
  Pixel plot_color = color.color;
  if (effective_debug)
    plot_color = plot_color.lerp16(Pixel(65535, 65535, 65535), 65535 / 2);
  scan_region<W, H>(
      y_lo, y_hi,
      [&](int y, auto &&out) {
        return shape.template get_horizontal_intervals<W, H>(y, out);
      },
      [&](int x, int y, const Vector &p, int max_run) {
        if (effective_debug) {
          for (int i = 0; i < max_run; ++i)
            pipeline.plot(canvas, x + i, y, plot_color, 0, 1.0f);
          return max_run;
        }

        float d;
        if constexpr (SineDistance) {
          d = shape.sine_distance(p);
        } else {
          shape.template distance<false>(p, result);
          d = result.dist;
        }

        // The color is constant, so a block the surface cannot cross is exact
        // whether skipped or splatted at full coverage; edge blocks stay
        // per-column.
        int span = 1;
        if constexpr (pole_lod_blocks<decltype(shape)>) {
          if (max_run > 1) {
            const float block_slack =
                pole_lod_block_slack<W>(max_run, p.y, shape);
            if (pole_lod_block_settles<decltype(shape)>(d, PIXEL_WIDTH,
                                                        block_slack))
              return max_run;
            // A sentineled subtrahend loses Subtract's max, so an inside report
            // is bounded no further than a clear one.
            if (pole_lod_block_settles<decltype(shape)>(-d, PIXEL_WIDTH,
                                                        block_slack))
              span = max_run;
          }
        }
        if (d >= PIXEL_WIDTH)
          return span;

        const float coverage = solid_coverage(d, PIXEL_WIDTH);
        if (coverage <= MIN_ALPHA)
          return span;

        const float alpha = color.alpha * coverage;
        for (int i = 0; i < span; ++i)
          pipeline.plot(canvas, x + i, y, plot_color, 0, alpha);
        return span;
      },
      xc);
}

} // namespace Scan
