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
#include "engine/static_circular_buffer.h"
#include "render/canvas.h"
#include "engine/platform.h"

#ifdef HS_AA_AUDIT
#include "render/aa_audit.h"
#endif

/**
 * @file scan.h
 * @brief The scanline rasterizer: Scan::rasterize plus the SDF-backed Scan draw
 * primitives.
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
 * @details The single spelling of the solid AA ramp, shared by every scan loop
 * and by the HS_AA_AUDIT walk so an audit miss is never the two disagreeing.
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
      const float sin_phi = sqrtf(std::max(0.0f, 1.0f - p.y * p.y));
      const float block_slack =
          pole_lod_slack<W, decltype(shape)>(max_run, sin_phi);
      if (d >= threshold + block_slack &&
          probe_bounds_block<decltype(shape)>(threshold, block_slack))
        return max_run;
      // An inside probe is never sentineled, so the splat needs no margin.
      if (solid && d <= -pixel_width - block_slack)
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
 * fractional column do not both paint it (double shade / alpha).
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
            const float sin_phi = sqrtf(std::max(0.0f, 1.0f - p.y * p.y));
            const float block_slack =
                pole_lod_slack<W, decltype(shape)>(max_run, sin_phi);
            if (d >= PIXEL_WIDTH + block_slack &&
                probe_bounds_block<decltype(shape)>(PIXEL_WIDTH, block_slack))
              return max_run;
            if (d <= -PIXEL_WIDTH - block_slack)
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
  static void draw_flat(PipelineT &pipeline, Canvas &canvas, const Basis &basis,
                        float radius, float thickness,
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
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   * @param suppress_pole_fill Drop the degenerate exact-pole row instead of
   *        full-row filling it (dense ring stacks; see get_horizontal_intervals).
   */
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const Basis &basis,
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
   * suppress_pole_fill.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  static void draw(PipelineT &pipeline, Canvas &canvas, int n_rings,
                   const SDF::DistortedRing *shapes, const int8_t *slot_by_ring,
                   int n_slots, RingShaderT &&shader) {
    HS_CHECK(canvas.width() == W && canvas.height() == H);
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
    float circumradius = res.second * (PI_F / 2.0f);

    SDF::PlanarPolygon shape(res.first, circumradius, sides, phase,
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
    float circumradius = res.second * (PI_F / 2.0f);
    SDF::PlanarPolygon shape(res.first, circumradius, sides, phase,
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
   *        bounding-box tint keeps per-shape scan bounds. That fallback fills
   *        v0/v1/v3 per pixel, so a shader reading them renders differently
   *        under debug.
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
    // Circle is a Ring with inner radius 0.
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

HS_O3_BEGIN
/**
 * @brief Face-specialized scan loop for the mesh path, compiled at -O3 on the
 *        -Os device image.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam PipelineT Plotting pipeline type.
 * @tparam MinimalFragment When true, only v1 (raw distance) is refreshed per
 *         pixel. Nothing else ever writes the fragment: pos, v0, v2, v3, size,
 *         age and color hold their struct defaults until the shader itself
 *         writes them, and from then on whatever the previous pixel left. The
 *         shader must therefore write frag.color on every pixel it is called
 *         for, and must not read a register it did not set — with v2 at 0 and
 *         size at 1, mesh_face_index() reports face 0 and fragment_edge_dist()
 *         returns an unnormalized distance for every face.
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Face to rasterize.
 * @param fragment_shader Shader invoked per covered pixel.
 * @details Self-contained (the shared rasterize/scan_region/process_pixel
 * kernel stays -Os; GCC reuses that existing -Os instantiation rather than
 * re-optimizing it into a region caller). A Face's column intervals are
 * row-independent unless it touches a pole, so for every other face the
 * wrap/sort/coalesce pass — per-row in scan_region — and the clip-arc
 * intersection run once; the per-pixel body mirrors process_pixel's solid path.
 */
template <int W, int H, typename PipelineT, bool MinimalFragment = false,
          typename FragmentShaderT>
HS_NOINLINE_NOCLONE inline void
rasterize_face(PipelineT &pipeline, Canvas &canvas, const SDF::Face &shape,
               FragmentShaderT &fragment_shader) {
  const auto &cr = canvas.clip();
  const auto xc = cr.x_clip();

  int y_lo, y_hi;
  static constexpr size_t MAX_RUNS = 8;
  std::pair<int, int> runs[MAX_RUNS];
  size_t num_runs = 0;

  {
    HS_PROFILE_DEEP(raster_setup);
    auto bounds = shape.template get_vertical_bounds<H>();
    if (!clamp_rows_to_clip(bounds, cr, y_lo, y_hi))
      return;

    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
  }

  auto build_runs = [&](int y) {
    HS_PROFILE_DEEP(raster_setup);
    num_runs = 0;

    // Each span wraps into <= 2 seam pieces, each clip-split into <= 2 runs.
    StaticCircularBuffer<SDF::Interval, 4> intervals;
    static_assert(decltype(intervals)::CAPACITY >=
                      SDF::sdf_max_spans<SDF::Face>::value,
                  "intervals must hold a Face's per-row emission");
    static_assert(MAX_RUNS >= 4 * SDF::sdf_max_spans<SDF::Face>::value,
                  "runs must hold every span seam-split and clip-split");
    bool handled = shape.template get_horizontal_intervals<W, H>(
        y, [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });

    // Hand-rolled rather than the shared clip_run/emit_row_runs: this consumer
    // materializes runs into a cached array instead of walking them, and either
    // shared spelling costs 272 B of ITCM here.
    auto add_run = [&](int x1, int x2) {
      auto push = [&](int a, int b) {
        if (a >= b)
          return;
        HS_CHECK(num_runs < MAX_RUNS, "rasterize_face: run buffer overflow");
        runs[num_runs++] = {a, b};
      };
      if (!xc.active) {
        push(x1, x2);
      } else if (xc.wrap) {
        push(x1, std::min(x2, xc.re));
        push(std::max(x1, xc.rs), x2);
      } else {
        push(std::max(x1, xc.rs), std::min(x2, xc.re));
      }
    };

    if (!handled) {
      add_run(0, W);
    } else if (!intervals.is_empty()) {
      bool full_row = false;
      for (const auto &iv : intervals) {
        if (iv.second - iv.first >= static_cast<float>(W)) {
          full_row = true;
          break;
        }
      }
      if (full_row) {
        add_run(0, W);
      } else {
        StaticCircularBuffer<SDF::Interval, 8> norm;
        static_assert(decltype(norm)::CAPACITY ==
                          2 * decltype(intervals)::CAPACITY,
                      "norm must hold 2 spans per input interval (seam split)");
        coalesce_spans<W>(intervals, norm);
        for (const auto &run : norm)
          add_run(static_cast<int>(run.first), static_cast<int>(run.second));
      }
    }
  };

  // A pole-touching face widens its azimuth wedge toward the pole, so its runs
  // describe one row; every other face reuses one set for the whole band.
  const bool per_row = shape.pole_touch;
  if (!per_row) {
    build_runs(y_lo);
#ifndef HS_AA_AUDIT
    // The audit walk below visits every pixel regardless of the runs, so it
    // must not take this early out.
    if (num_runs == 0)
      return;
#endif
  }

  const float *cos_theta = TrigLUT<W, H>::sin_theta.data() + W / 4;
  const float *sin_theta = TrigLUT<W, H>::sin_theta.data();
  constexpr float pixel_width = 2.0f * PI_F / W;
  const uint32_t probe_flags = shape.probe_flags();

#ifdef HS_AA_AUDIT
  if (hs_aa::g_audit.enabled) {
    SDF::DistanceResult ares;
    for (int y = y_lo; y <= y_hi; ++y) {
      if (per_row)
        build_runs(y);
      float sp = TrigLUT<W, H>::sin_phi[y];
      float cp = TrigLUT<W, H>::cos_phi[y];
      for (int x = 0; x < W; ++x) {
        Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
        shape.template distance<true>(p, ares);
        float d = ares.dist;
        if (d >= pixel_width)
          continue;
        const float alpha = solid_coverage(d, pixel_width);
        if (alpha <= MIN_ALPHA)
          continue;
        int gap = W;
        bool in_run = false;
        for (size_t r = 0; r < num_runs; ++r) {
          if (x >= runs[r].first && x < runs[r].second) {
            in_run = true;
            break;
          }
          // Circular column distance to the nearest column of this run.
          const auto circ = [&](int a, int b) {
            int dd = ((a - b) % W + W) % W;
            return dd < W - dd ? dd : W - dd;
          };
          int g = std::min(circ(x, runs[r].first), circ(x, runs[r].second - 1));
          if (g < gap)
            gap = g;
        }
        if (in_run)
          hs_aa::g_audit.note_painted();
        else
          hs_aa::g_audit.note_missed(y, alpha, gap);
      }
    }
    for (int y = y_lo; y <= y_hi; ++y) {
      if (per_row)
        build_runs(y);
      for (size_t r = 0; r < num_runs; ++r) {
        const long long len = runs[r].second - runs[r].first;
        hs_aa::g_audit.note_probes(y, len);
      }
    }
  }
  if (!per_row && num_runs == 0)
    return;
#endif

  SDF::DistanceResult res;
  Fragment frag;

  // 1.0002 margin keeps the sqrt-free cull strictly conservative against the
  // widest threshold a probe is tested on below.
  const float base_reject_rad = pixel_width * 1.0002f;
  const float base_reject_dsq = base_reject_rad * base_reject_rad;

  // distance() reports gnomonic-plane units, which stretch an angular step by
  // up to 1 + r^2; max_dist bounds r over every probe the cull admits.
  [[maybe_unused]] const float plane_stretch = 1.0f + shape.max_dist_sq;

  HS_PROFILE_DEEP(raster_scan);
  for (int y = y_lo; y <= y_hi; ++y) {
    if (per_row) {
      build_runs(y);
      if (num_runs == 0)
        continue;
    }
    float sp = TrigLUT<W, H>::sin_phi[y];
    float cp = TrigLUT<W, H>::cos_phi[y];
    const int stride = pole_lod_run(sp);
    float block_slack = 0.0f;
    float reject_dsq = base_reject_dsq;
    if constexpr (pole_lod_blocks<SDF::Face>) {
      block_slack = pole_lod_slack<W, SDF::Face>(stride, sp) * plane_stretch;
      // the block test below reads d, so the cull must not sentinel a probe
      // inside that threshold
      const float reject_rad = (pixel_width + block_slack) * 1.0002f;
      reject_dsq = reject_rad * reject_rad;
    }

    for (size_t r = 0; r < num_runs; ++r) {
      const int rx2 = runs[r].second;
      for (int x = runs[r].first; x < rx2;) {
        Vector p(sp * cos_theta[x], cp, sp * sin_theta[x]);
        shape.template distance_with_flags<true>(p, res, reject_dsq,
                                                 probe_flags);
        const float d = res.dist;

        // Columns this one shade covers. Only a canvas-aligned block that fits
        // in the run and that the surface cannot cross qualifies; anything
        // holding an edge stays per-column, so coverage is identical either way
        // and only the shade's source column moves. The clear side additionally
        // needs the probe to carry a distance rather than the cull's sentinel;
        // an inside probe is never culled.
        int span = 1;
        if constexpr (pole_lod_blocks<SDF::Face>) {
          if (stride > 1 && x % stride == 0 && x + stride <= rx2 &&
              (d <= -pixel_width - block_slack ||
               (d >= pixel_width + block_slack &&
                probe_bounds_block<SDF::Face>(pixel_width, block_slack))))
            span = stride;
        }

        if (d >= pixel_width) {
          x += span;
          continue;
        }
        HS_SCAN_METRIC(hs::g_scan_metrics.shade_candidates++);
        HS_PROBE_COUNT(n_alpha);
        HS_PROBE_MARK(hs_ta);
        const float alpha = solid_coverage(d, pixel_width);
        HS_PROBE_SPAN(alpha, hs_ta);
        if (alpha <= MIN_ALPHA) {
          x += span;
          continue;
        }
        HS_PROFILE_DEEP(raster_shade);
        if constexpr (MinimalFragment) {
          frag.v1 = res.raw_dist;
        } else {
          frag.color = Color4(0, 0, 0, 0);
          frag.pos = p;
          frag.v0 = res.t;
          frag.v1 = res.raw_dist;
          frag.v2 = 0.0f;
          frag.v3 = res.aux;
          frag.size = res.size;
          frag.age = 0;
        }
        fragment_shader(p, frag);
        if (frag.color.alpha > MIN_ALPHA) {
          const float a = frag.color.alpha * alpha;
          for (int px = x; px < x + span; ++px) {
            if constexpr (requires {
                            pipeline.plot_in_bounds(
                                canvas, px, y, frag.color.color, frag.age, a);
                          }) {
              pipeline.plot_in_bounds(canvas, px, y, frag.color.color, frag.age,
                                      a);
            } else {
              pipeline.plot(canvas, px, y, frag.color.color, frag.age, a);
            }
          }
        }
        x += span;
      }
    }
  }
}
HS_O3_END

/**
 * @brief Arena-hosts one per-mesh Face scratch buffer.
 * @param arena Scratch arena supplying the storage.
 * @return The constructed buffer, valid until the caller's ScratchScope closes.
 * @details Default-init: every field is written by SDF::Face before it is read,
 * and value-init would zero ~7.7 KB per mesh draw.
 */
HS_NOINLINE_NOCLONE inline SDF::FaceScratchBuffer *
new_face_scratch(Arena &arena) {
  return arena.make_default<SDF::FaceScratchBuffer>();
}

HS_O3_BEGIN
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
   * @param bake Optional congruence-class bake for this mesh (null = exact
   *        path everywhere, today's behavior). When present, each face is
   *        aligned to its canonical class shape after construction and the
   *        class distance LUT is bound for the probe loop.
   * @param face_shader_setup Optional callback receiving the face index and
   *        size. When supplied, it runs once before rasterizing the face and
   *        the shader is invoked directly with only v1 refreshed per pixel —
   *        the callback never sees the fragment, so the index and size reach
   *        the shader only through what it hoists here, not through v2 and
   *        size. The fragment is per-face, so a shader on this path must write
   *        frag.color unconditionally: returning early leaves the previous
   *        pixel's color to be plotted again.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_NOINLINE_NOCLONE static void
  draw_impl(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
            FragmentShaderT &fragment_shader, Arena &scratch_arena,
            const MeshOps::MeshClassBake *bake,
            FaceShaderSetupT &face_shader_setup) {
    // Once per mesh, not per face: rasterize_face indexes the phi LUT by the
    // canvas' own rows and hands SDF::Face the template H.
    check_canvas_dims<W, H>(canvas);
    // The per-face wrapper below is itself always non-null, so the erased
    // shader it wraps has to be checked here or not at all.
    if constexpr (std::is_same_v<FragmentShaderT, FragmentShaderFn>)
      check_fragment_shader(fragment_shader);

    ScratchScope scope(scratch_arena);
    auto *scratch = new_face_scratch(scratch_arena);

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
        return SDF::Face(verts, indices, *scratch, H + hs::H_OFFSET, H,
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

      if constexpr (std::is_same_v<FaceShaderSetupT, std::nullptr_t>) {
        auto wrapper = [&](const Vector &p, Fragment &f_in) {
          // v2 carries the face index (decoded by mesh_face_index()).
          // Exact for i < 2^24 (float mantissa); meshes never approach that face count.
          f_in.v2 = static_cast<float>(i);
          fragment_shader(p, f_in);
        };
        FragmentShaderFn raster_shader = wrapper;
        {
          HS_PROFILE(scan_mesh_raster);
          rasterize_face<W, H, PipelineT, false>(pipeline, canvas, shape,
                                                 raster_shader);
        }
      } else {
        face_shader_setup(i, shape.size);
        {
          HS_PROFILE(scan_mesh_raster);
          rasterize_face<W, H, PipelineT, true>(pipeline, canvas, shape,
                                                fragment_shader);
        }
      }
    }
  }

  /** @brief Rasterizes a mesh through the type-erased fragment shader path. */
  template <int W, int H, typename PipelineT = PipelineRef,
            typename FaceShaderSetupT = std::nullptr_t>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake = nullptr,
                   FaceShaderSetupT face_shader_setup = nullptr) {
    draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                    bake, face_shader_setup);
  }

  /** @brief Rasterizes a mesh with an inlinable per-face fragment shader. */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  static void
  draw_specialized(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderT &fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake,
                   FaceShaderSetupT &face_shader_setup) {
    FunctionRef<void(size_t, float)> erased_setup = face_shader_setup;
    // Cold (once per draw): the erasure hides a null setup from draw_impl's
    // nullptr_t branch, which would then call a null thunk per face.
    HS_CHECK(erased_setup, "Scan::Mesh::draw_specialized requires a non-null "
                           "face_shader_setup");
    draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                    bake, erased_setup);
  }
};
HS_O3_END

/**
 * @brief Full-screen per-pixel shaders with SAMPLES× SSAA.
 *
 * Three entry points, in increasing order of caller control:
 * - draw(canvas, shader): one callable ShaderFn(const Vector &v) -> Color4,
 *   invoked SAMPLES× per pixel at sub-pixel offsets and averaged.
 * - draw(canvas, fragment_shader, vertex_shader): splits per-pixel setup
 *   (vertex_shader, once at the pixel center) from per-sub-sample evaluation.
 *   Both callables are required; a null one traps.
 * - draw_grid(canvas, vertex_shader, pixel_shader): hands the seeded fragment
 *   and the row's SsaaGrid to pixel_shader, which owns the sampling and returns
 *   the finished pixel.
 *
 * @details Every entry point assigns the finished premultiplied color to the
 * canvas rather than plotting it, so the destination is overwritten: alpha < 1
 * darkens the pixel instead of blending with what is under it, and no
 * plot-time filter stage (Filter::World / Filter::Pixel) sees it. These entry
 * points take no pipeline; an effect needing the filter chain must plot
 * through it itself.
 */
struct Shader {
  // --- Shared SSAA helpers (used by every entry point) -----------------------
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
    /** @brief Sub-pixel samples per pixel supplied by the 2×2 grid. */
    static constexpr int SAMPLES = 4;

    float sin_phi[2]; /**< Current row's sin(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_phi[2]; /**< Current row's cos(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_dtheta; /**< cos of the ±0.25 px column rotation. */
    float sin_dtheta; /**< sin of the ±0.25 px column rotation. */
    float cos_dphi;   /**< cos of the ±0.25 px row rotation. */
    float sin_dphi;   /**< sin of the ±0.25 px row rotation. */

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
     * @param i Sample index in [0, SAMPLES); the low bit selects the column
     * offset (±0.25 px) and bit 1 the row offset (±0.25 px).
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
   * @brief Validates the LUT-domain invariant shared by every entry point.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param cr Clip region whose bounds are checked against the LUT extents.
   * @details Checked once per draw, not per pixel: every (x,y) the loops feed to
   * pixel_to_vector indexes the trig LUTs within bounds. The row bound is H, not
   * the LUT's H_VIRT (H plus the pole offset): those rows also subscript the
   * canvas, which holds H of them, so the tighter bound is the binding one.
   */
  template <int W, int H> static void check_lut_domain(const ClipRegion &cr) {
    HS_CHECK(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 &&
             cr.render_y_end() <= H);
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
  HS_O3_FN static void draw(Canvas &canvas, ShaderFn &&shader) {
    // The sample-offset table has four distinct sub-pixel positions; only 1 and
    // the 2x2 grid (4) are valid.
    static_assert(SAMPLES == 1 || SAMPLES == 4,
                  "Scan::Shader SSAA supports only SAMPLES == 1 or 4");
    check_canvas_dims<W, H>(canvas);
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto xc = cr.x_clip();

    if constexpr (SAMPLES == 1) {
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        walk_clip_columns<W>(xc, [&](int x) {
          Vector v = pixel_to_vector<W, H>(x, y);
          Color4 sample = shader(v);
          canvas(x, y) = sample.color * sample.alpha;
        });
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      if (!TrigLUT<W, H>::initialized)
        TrigLUT<W, H>::init();
      SsaaGrid<W, H> grid;

      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        grid.set_row(y);
        walk_clip_columns<W>(xc, [&](int x) {
          // Premultiplied SSAA: accumulate each sample's coverage-weighted color
          // (color * alpha / N), matching the SAMPLES==1 path.
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            Color4 sample = shader(grid.at(x, i));
            accum += sample.color * (sample.alpha * inv_samples);
          }

          canvas(x, y) = accum;
        });
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
    check_canvas_dims<W, H>(canvas);
    // frag_base is per pixel, not per draw: each pixel starts from a default
    // Fragment, so a vertex shader writing only some registers (v0-v3/size/age/
    // color) can't inherit the previous pixel's values.
    if constexpr (SAMPLES == 1) {
      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        walk_clip_columns<W>(xc, [&](int x) {
          Vector center_v = pixel_to_vector<W, H>(x, y);
          Fragment frag_base;
          frag_base.pos = center_v;
          vertex_shader(frag_base);
          fragment_shader(center_v, frag_base);
          // Premultiply by alpha, matching the single-callback overload and the
          // process_pixel/Volume contract.
          canvas(x, y) = frag_base.color.color * frag_base.color.alpha;
        });
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      if (!TrigLUT<W, H>::initialized)
        TrigLUT<W, H>::init();
      SsaaGrid<W, H> grid;

      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        grid.set_row(y);
        walk_clip_columns<W>(xc, [&](int x) {
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
        });
      }
    }
  }

  /**
   * @brief Full-screen draw handing the effect the whole per-pixel body.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam VertexFn Callable VertexFn(Fragment&) — per-pixel seed.
   * @tparam PixelFn Callable PixelFn(Fragment&, const SsaaGrid<W,H>&, int x)
   *         -> Pixel — computes the finished (already premultiplied) pixel.
   * @param canvas Destination canvas.
   * @param vertex_shader Per-pixel shader, called once at the pixel center.
   * @param pixel_shader Owns the pixel: receives the seeded fragment, the
   *        sub-pixel SSAA grid for the current row, and the pixel column, and
   *        returns the final (premultiplied) pixel. A template, not a
   *        type-erased FunctionRef: the whole body inlines into this loop, and
   *        the effect can hoist work its sub-samples share.
   * @details Same outer scaffolding as the SSAA draw() overloads (clip,
   * LUT-domain check, trig-LUT init, per-row SsaaGrid); the
   * per-pixel work is delegated whole so the caller controls the sampling.
   */
  template <int W, int H, typename VertexFn, typename PixelFn>
  HS_O3_FN static void draw_grid(Canvas &canvas, VertexFn &&vertex_shader,
                                 PixelFn &&pixel_shader) {
    check_canvas_dims<W, H>(canvas);
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto xc = cr.x_clip();
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    SsaaGrid<W, H> grid;
    for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
      grid.set_row(y);
      walk_clip_columns<W>(xc, [&](int x) {
        Fragment frag_base;
        frag_base.pos = pixel_to_vector<W, H>(x, y);
        vertex_shader(frag_base);
        canvas(x, y) = pixel_shader(frag_base, grid, x);
      });
    }
  }
};

// Volume::draw has exactly one caller (effects/Raymarch.h), and
// TransformedVolume is instantiated only there, so no other effect pays ITCM
// for this region.
HS_O3_BEGIN
/**
 * @brief Generic wrapper that places an SDF in world space via a center point
 *        and a rotation quaternion. Satisfies the Volume::draw shape concept.
 * @tparam SDF Underlying signed-distance shape type.
 * @details The quaternion q maps local→world: world_p = center +
 * rotate(local_p, q). ray_to_local uses q.inverse() to map world→local.
 */
template <typename SDF> struct TransformedVolume {
  const SDF &sdf;   /**< Underlying SDF evaluated in local space. */
  Vector center;    /**< World-space origin of the local frame. */
  Quaternion q_inv; /**< Precomputed inverse rotation (world→local). */

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
  /** Sphere-trace step overrelaxation; 1 is the conservative march. */
  static constexpr float OVERRELAX_OMEGA = 1.3f;
  /** Occluder-probe march steps. */
  static constexpr int PROBE_STEPS = 24;
  /** Probe steps taken at the near step floor before the far one applies. */
  static constexpr int PROBE_NEAR_STEPS = 6;
  /** Near probe step floor, as a fraction of the bounding radius. */
  static constexpr float PROBE_FLOOR_NEAR = 0.04f;
  /** Far probe step floor, as a fraction of the bounding radius. */
  static constexpr float PROBE_FLOOR_FAR = 0.12f;

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
   *
   * Steps are overrelaxed by OVERRELAX_OMEGA (Keinert et al., "Enhanced Sphere
   * Tracing"): each sample's unbounding sphere must overlap its predecessor's,
   * and a step that breaks that overlap is rewound to the predecessor's surface
   * and the ray finished conservatively, so the trace still cannot cross a
   * surface undetected.
   */
  template <typename Shape>
  static __attribute__((always_inline)) float
  trace_closest(const Shape &shape, const Vector &local_ro,
                const Vector &local_vd, float bounds_radius, int max_steps,
                float aa_width, Vector &closest_local) {
    HS_PROFILE_DEEP(vol_trace);
    Vector local_p = local_ro;
    closest_local = local_ro;
    // Sentinel for "no surface seen yet": any real signed distance the
    // trace reports is smaller, so the first sample always wins.
    float closest_d = FLT_MAX;
    // Drops to 1 for the rest of the ray once a step fails the overlap test.
    float omega = OVERRELAX_OMEGA;
    float prev_r = 0.0f;
    float step_len = 0.0f;

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

      float r = d < 0.0f ? -d : d;
      if (omega > 1.0f && r + prev_r < step_len) {
        // Overrelaxed step left the previous unbounding sphere, so the interval
        // it skipped is unverified: rewind to that sphere's surface and finish
        // the ray conservatively. The rejected sample updates nothing.
        float back = prev_r - step_len;
        local_p =
            Vector(local_p.x + local_vd.x * back, local_p.y + local_vd.y * back,
                   local_p.z + local_vd.z * back);
        omega = 1.0f;
        prev_r = 0.0f;
        step_len = 0.0f;
        continue;
      }
      prev_r = r;

      if (d < closest_d) {
        closest_d = d;
        closest_local = local_p;
        // A frontal hit converges on the surface from outside, so the
        // d < -aa_width break below never fires for it.
        if (closest_d <= aa_width * 0.02f)
          break;
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
      step_len = std::max(d * 0.9f * omega, 1e-5f);
      local_p = Vector(local_p.x + local_vd.x * step_len,
                       local_p.y + local_vd.y * step_len,
                       local_p.z + local_vd.z * step_len);
    }
    return closest_d;
  }

  /**
   * @brief Result of probing behind an AA halo for an occluded surface.
   */
  struct Occluder {
    bool solid; /**< A solid surface sits behind the halo (behind is valid). */
    Vector behind; /**< Local-space hit point when solid, else the grazed
                        background edge's closest approach when soft > 0. */
    float
        soft; /**< Coverage of a grazed background edge, for the corner fill. */
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
  static __attribute__((always_inline)) Occluder probe_occluder(
      const Shape &shape, const Vector &closest_local, const Vector &local_vd,
      float bounds_radius, float hit_threshold, float aa_width) {
    HS_PROFILE_DEEP(vol_probe);
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
    for (int i = 0; i < PROBE_STEPS; ++i) {
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
      float floor = bounds_radius *
                    (i < PROBE_NEAR_STEPS ? PROBE_FLOOR_NEAR : PROBE_FLOOR_FAR);
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
   * @param bounds_center Bounding sphere center in physical LED space; must be
   *        a unit vector (on the canvas sphere).
   * @param bounds_radius Bounding sphere radius in world units.
   * @param view_dir Ray direction (camera → scene) in LED space; must point
   *        straight at bounds_center (a radial view). The scanned band is a cap
   *        around bounds_center, which is the orthographic footprint only under
   *        that view; a tilted one slides the footprint off the band and drops
   *        covered columns.
   * @param shape Volume shape providing ray_to_local() and distance().
   * @param frag_fn Fragment shader invoked once per hit.
   * @param max_steps Maximum sphere-tracing steps per ray.
   * @param aa_width Anti-aliasing band half-width in world units.
   * @details Shape concept:
   *   std::pair<Vector, Vector> ray_to_local(const Vector &ro, const Vector
   *   &vd) const; Vector origin_to_local(const Vector &ro) const; float
   *   distance(const Vector &local_point) const;
   *
   * Fragments are plotted at pixel centers by integer coordinates, so the
   * pipeline receives no sub-pixel positions from this draw.
   */
  template <int W, int H, typename Shape>
  static void
  draw(PipelineRef pipeline, Canvas &canvas, const Vector &bounds_center,
       float bounds_radius, const Vector &view_dir, const Shape &shape,
       FragmentShaderFn frag_fn, int max_steps = 15, float aa_width = 0.01f) {
    check_canvas_dims<W, H>(canvas);
    check_fragment_shader(frag_fn);

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
                 local_bc.z * local_bc.z <
             TOLERANCE);
    // The scan band below is a cap around bounds_center of angular radius
    // asin(bounds_radius), which equals the orthographic footprint only for a
    // radial view of a unit-length center: BoundingSphere reads center.y as
    // cos(phi), and a tilted view slides the footprint off the cap. Unit length
    // also backs the ray start offset above — farther out along the view axis a
    // ray can start in front of the shape.
    HS_CHECK(fabsf(dot(bounds_center, bounds_center) - 1.0f) < TOLERANCE);
    const Vector radial_err = cross(bounds_center, vd);
    HS_CHECK(bc_dot_vd < 0.0f && dot(radial_err, radial_err) < TOLERANCE);
    // aa_width > 0 is the contract: the slow-path AA divides by (aa_width -
    // hit_threshold) == 0.9*aa_width, so a zero band-width gives 0/0 -> NaN.
    HS_CHECK(aa_width > 0.0f);

    BoundingSphere<W, H> bounds(bounds_center, bounds_radius);

    // Tier 2: Clamp volume bounds to clip region
    const auto &cr = canvas.clip();
    const auto vol_xc = cr.x_clip();
    int vol_y_lo, vol_y_hi;
    if (!clamp_rows_to_clip(bounds, cr, vol_y_lo, vol_y_hi))
      return;

    scan_region<W, H>(
        vol_y_lo, vol_y_hi,
        [&](int y, auto &&out) { return bounds.get_intervals(y, out); },
        [&](int px, int py, const Vector &p, int max_run) {
          // Back-face cull
          float facing = p.x * vd.x + p.y * vd.y + p.z * vd.z;
          if (facing >= 0.0f)
            return 1;

          // Orthographic ray-sphere cull
          float pp_x = p.x - facing * vd.x;
          float pp_y = p.y - facing * vd.y;
          float pp_z = p.z - facing * vd.z;
          float dx = pp_x - bc_proj.x;
          float dy = pp_y - bc_proj.y;
          float dz = pp_z - bc_proj.z;
          if (dx * dx + dy * dy + dz * dz > bounds_r2)
            return 1;

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

          if (closest_d >= aa_width) {
            // Shifting the ray origin one column moves every sampled point by
            // at most that chord, and the SDF is 1-Lipschitz, so clearance
            // beyond the offer's own width holds for the whole block. Traced
            // hits stay per-column: shading varies along the surface, so a
            // splat would band.
            if constexpr (pole_lod_blocks<decltype(shape)>) {
              if (max_run > 1) {
                const float sin_phi = sqrtf(std::max(0.0f, 1.0f - p.y * p.y));
                const float block_slack =
                    pole_lod_slack<W, decltype(shape)>(max_run, sin_phi);
                if (closest_d >= aa_width + block_slack)
                  return max_run;
              }
            }
            return 1;
          }

          // --- Fragment shading ---
          Fragment frag;
          frag.pos = closest_local;
          frag.size = closest_d;
          {
            HS_PROFILE_DEEP(vol_shade);
            frag_fn(closest_local, frag);
          }

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
            Occluder occ =
                probe_occluder(shape, closest_local, local_vd, bounds_radius,
                               hit_threshold, aa_width);
            if (occ.solid) {
              // Self-occlusion edge: antialias the foreground over the surface it
              // covers — lay the shaded background down, then blend the foreground
              // over it by the edge coverage. Smooth, vs. fading to black (fringe)
              // or snapping to opaque (jagged).
              Fragment bg;
              bg.pos = occ.behind;
              bg.size = 0.0f;
              {
                HS_PROFILE_DEEP(vol_shade);
                frag_fn(occ.behind, bg);
              }
              if (bg.color.alpha > MIN_ALPHA) {
                HS_PROFILE_DEEP(vol_plot);
                pipeline.plot(canvas, px, py, bg.color.color, 0.0f,
                              bg.color.alpha);
              }
              if (frag.color.alpha * edge_alpha > MIN_ALPHA) {
                HS_PROFILE_DEEP(vol_plot);
                pipeline.plot(canvas, px, py, frag.color.color, 0.0f,
                              frag.color.alpha * edge_alpha);
              }
              return 1;
            }
            // No solid behind; a grazed background edge fills the corner,
            // shaded at its own point so the fill carries the background's
            // color, then the foreground blends over it.
            if (occ.soft > MIN_ALPHA) {
              Fragment bg;
              bg.pos = occ.behind;
              bg.size = 0.0f;
              {
                HS_PROFILE_DEEP(vol_shade);
                frag_fn(occ.behind, bg);
              }
              if (bg.color.alpha * occ.soft > MIN_ALPHA) {
                HS_PROFILE_DEEP(vol_plot);
                pipeline.plot(canvas, px, py, bg.color.color, 0.0f,
                              bg.color.alpha * occ.soft);
              }
            }
          }

          if (frag.color.alpha * edge_alpha > MIN_ALPHA) {
            HS_PROFILE_DEEP(vol_plot);
            pipeline.plot(canvas, px, py, frag.color.color, 0.0f,
                          frag.color.alpha * edge_alpha);
          }
          return 1;
        },
        vol_xc);
  }
};
HS_O3_END

} // namespace Scan
