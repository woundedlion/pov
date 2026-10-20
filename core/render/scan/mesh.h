/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cstddef>
#include <span>
#include <type_traits>
#include <utility>
#include "core/render/sdf.h"
#include "core/render/shading.h"
#include "core/render/sdf/face_classes.h"
#include "core/mesh/mesh_state.h"
#include "core/color/color.h"
#include "core/render/filter/pipeline.h"
#include "core/containers/static_circular_buffer.h"
#include "core/render/canvas.h"
#include "core/platform/platform.h"
#include "core/render/scan/raster.h"

/**
 * @file mesh.h
 * @brief rasterize_face, the -O3 face-specialized scan loop, and Scan::Mesh over it.
 */

namespace Scan {

/** @brief Half-open column run [first, second); left uninitialized. */
struct ColumnRun {
  int first;  ///< First column.
  int second; ///< One past the last column.
};

/**
 * @brief Whether a face's clipped row band can reach either arm of a window.
 * @param shape Constructed face.
 * @param window Render-x window.
 * @param clip Source clip defining the row band and display dimensions.
 * @return True when either arm can receive a pixel, or the window is inactive.
 */
HS_HOT_FLASH_MEMBER inline bool face_meets_window(const SDF::Face &shape,
                                                  const ColumnWindow &window,
                                                  const ClipRegion &clip) {
  if (!window.active)
    return true;
  ClipRegion side{clip.render_y_start(),
                  clip.render_y_end(),
                  window.begin,
                  window.end,
                  0,
                  clip.w,
                  clip.h};
  if (!shape.clip_rejects(side))
    return true;
  side.x_start += clip.w / 2;
  side.x_end += clip.w / 2;
  return !shape.clip_rejects(side);
}

/**
 * @brief Clips a row's column runs to a render-x window on both arms.
 * @tparam W Canvas width in pixels.
 * @tparam CAPACITY Run buffer capacity.
 * @param runs The row's runs; receives the clipped runs, each on one arm.
 * @param count Number of runs.
 * @param window Open render-x window.
 * @return Number of clipped runs.
 */
template <int W, size_t CAPACITY>
HS_HOT_INLINE inline size_t window_runs(ColumnRun *runs, size_t count,
                                        const ColumnWindow &window) {
  ColumnRun kept[CAPACITY + 1];
  size_t kept_count = 0;
  for (size_t r = 0; r < count; ++r)
    for (const int arm : {0, W / 2}) {
      const int a = std::max(runs[r].first, window.begin + arm);
      const int b = std::min(runs[r].second, window.end + arm);
      if (a < b)
        kept[kept_count++] = {a, b};
    }
  HS_CHECK(kept_count <= CAPACITY, "rasterize_face: run buffer overflow");
  std::copy_n(kept, kept_count, runs);
  return kept_count;
}

/**
 * @brief Whether a row's column intervals can leave a run inside a render-x
 *        window.
 * @tparam W Canvas width in pixels.
 * @tparam IntervalBufT Buffer of the row's column intervals.
 * @param intervals The row's intervals, in whole columns, unwrapped.
 * @param window Open render-x window.
 * @return False when every interval misses the window's columns on both arms,
 *         so that the row's clipped runs are empty.
 */
template <int W, typename IntervalBufT>
HS_HOT_INLINE inline bool intervals_meet_window(const IntervalBufT &intervals,
                                                const ColumnWindow &window) {
  const int length = window.end - window.begin;
  for (const auto &iv : intervals) {
    const int a = static_cast<int>(iv.start);
    // A zero-width interval still owns its column.
    const int len = std::max(static_cast<int>(iv.end) - a, 1);
    const int s = ((a % W) + W) % W;
    if (ClipRegion::arcs_overlap(window.begin, length, s, len, W) ||
        ClipRegion::arcs_overlap(window.begin + W / 2, length, s, len, W))
      return true;
  }
  return false;
}

HS_O3_BEGIN
/**
 * @brief Body of the face scan loop, for the fixed meridian or a rotated
 *        basis.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam PipelineT Plotting pipeline type.
 * @tparam MinimalFragment See rasterize_face.
 * @tparam COLUMN_TIME Whether pixels follow @p time's rotated basis.
 * @tparam INLINE_RUNS Whether the run builder inlines into the body, so a
 *         FLASH variant takes it along.
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Face to rasterize.
 * @param fragment_shader Shader invoked per covered pixel.
 * @param time Column-time state with a rotated basis; read only with
 *        COLUMN_TIME.
 * @details Under a rotated basis a pole-LOD block's slack grows by the
 * rotation between its first and last columns, and a block straddling an arm
 * boundary is walked per column. Only the columns an open RenderXWindow
 * admits are written.
 */
template <int W, int H, typename PipelineT, bool MinimalFragment,
          bool COLUMN_TIME, bool INLINE_RUNS, typename FragmentShaderT>
HS_HOT_INLINE inline void
rasterize_face_body(PipelineT &pipeline, Canvas &canvas, const SDF::Face &shape,
                    FragmentShaderT &fragment_shader,
                    const Temporal::ColumnTime *time) {
  const auto &cr = source_clip<W, H>(pipeline, canvas);
  const auto xc = cr.x_clip();

  int y_lo, y_hi;
  // A column-time or windowed run is split at the arm boundary, which one run
  // can cross.
  static constexpr size_t MAX_RUNS = 9;
  ColumnRun runs[MAX_RUNS];
  size_t num_runs = 0;
  const ColumnWindow &window = canvas.render_x_window();

  {
    HS_PROFILE_DEEP(raster_setup);
    auto bounds = shape.template get_vertical_bounds<H>();
    if (!clamp_rows_to_clip(bounds, cr, y_lo, y_hi))
      return;

    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
  }

  auto build_runs_body = [&](int y) HS_HOT_INLINE {
    HS_PROFILE_DEEP(raster_setup);
    num_runs = 0;

    // Each span wraps into <= 2 seam pieces, each clip-split into <= 2 runs.
    StaticCircularBuffer<SDF::Interval, 4> intervals;
    static_assert(decltype(intervals)::CAPACITY >=
                      SDF::sdf_max_spans<SDF::Face>::value,
                  "intervals must hold a Face's per-row emission");
    static_assert(MAX_RUNS >= 4 * SDF::sdf_max_spans<SDF::Face>::value + 1,
                  "runs must hold every span seam-split, clip-split and "
                  "arm-split");
    bool handled = shape.template get_horizontal_intervals<W, H>(
        y, [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });
    if (window.active && handled &&
        !intervals_meet_window<W>(intervals, window))
      return;

    auto append = [&](int a, int b) {
      HS_CHECK(num_runs < MAX_RUNS, "rasterize_face: run buffer overflow");
      runs[num_runs++] = {a, b};
    };
    auto push = [&](int a, int b) {
      if (a >= b)
        return;
      if constexpr (COLUMN_TIME) {
        if (a < W / 2 && b > W / 2) {
          append(a, W / 2);
          a = W / 2;
        }
      }
      append(a, b);
    };
    StaticCircularBuffer<SDF::Interval, 8> norm;
    static_assert(decltype(norm)::CAPACITY == 2 * decltype(intervals)::CAPACITY,
                  "norm must hold 2 spans per input interval (seam split)");
    emit_row_runs<W>(handled, intervals, norm, xc, push);
    if (window.active)
      num_runs = window_runs<W, MAX_RUNS>(runs, num_runs, window);
  };
  [[maybe_unused]] auto build_runs = [&](int y) { build_runs_body(y); };
  const auto rebuild_runs = [&](int y) HS_HOT_INLINE {
    if constexpr (INLINE_RUNS)
      build_runs_body(y);
    else
      build_runs(y);
  };

  const bool per_row =
      shape.template horizontal_intervals_vary_by_row<W, H>(y_lo, y_hi);
  int runs_y = y_lo;
  if (!per_row) {
    rebuild_runs(y_lo);
    if (num_runs == 0)
      return;
  }

  const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
  const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
  [[maybe_unused]] const math::Vector *azimuth = nullptr;
  [[maybe_unused]] const math::Vector *up = nullptr;
  // Arc the basis rotates per column of sweep time.
  [[maybe_unused]] float column_rotation = 0.0f;
  if constexpr (COLUMN_TIME) {
    azimuth = time->azimuth;
    up = time->up;
    column_rotation = time->sweep / (W / 2);
  }
  constexpr float PIXEL_WIDTH = 2.0f * math::PI_F / W;
  const uint32_t probe_flags = shape.probe_flags();
  constexpr bool TIMED_SHADER =
      COLUMN_TIME &&
      Temporal::reads_column_time<FragmentShaderT, const math::Vector &,
                                  Fragment &>;
  // A time-aware shader runs per column only while its state moves.
  [[maybe_unused]] bool per_column = TIMED_SHADER;
  if constexpr (TIMED_SHADER) {
    if constexpr (requires { fragment_shader.moves(); })
      per_column = fragment_shader.moves();
  }

  SDF::DistanceResult res;
  Fragment frag;

  // 1.0002 margin keeps the sqrt-free cull strictly conservative.
  const float base_reject_rad = PIXEL_WIDTH * 1.0002f;
  const float base_reject_dsq = base_reject_rad * base_reject_rad;
  [[maybe_unused]] const float plane_stretch = report_stretch(shape);
  float radial_cull = shape.max_dist;
  if (!(probe_flags & SDF::Face::PROBE_CONVEX)) {
    const float shade_reach = (probe_flags & SDF::Face::PROBE_LINEAR)
                                  ? PIXEL_WIDTH
                                  : tanf(PIXEL_WIDTH);
    radial_cull = shape.radius + shade_reach;
    if (probe_flags & SDF::Face::PROBE_HAS_LUT)
      radial_cull += (shape.lut_q_safe + 1) * shape.lut_dequant;
  }
  const float min_cos = 1.0f / sqrtf(1.0f + radial_cull * radial_cull);
  [[maybe_unused]] float block_min_cos = min_cos;
  if constexpr (pole_lod_blocks<SDF::Face>) {
    if (!(probe_flags & SDF::Face::PROBE_CONVEX) &&
        Render::pole_lod_aggressiveness > 0.0f) {
      float lod_reach = Render::pole_lod_aggressiveness * PIXEL_WIDTH;
      if constexpr (COLUMN_TIME)
        lod_reach += column_rotation * (Render::POLE_LOD_MAX_RUN - 1);
      const float block_reach = std::min(
          PIXEL_WIDTH + lod_reach * SDF::arc_stretch<SDF::Face> * plane_stretch,
          SDF::reject_margin<SDF::Face>);
      float block_radius =
          shape.radius + ((probe_flags & SDF::Face::PROBE_LINEAR)
                              ? block_reach
                              : tanf(block_reach));
      if (probe_flags & SDF::Face::PROBE_HAS_LUT)
        block_radius += (shape.lut_q_safe + 1) * shape.lut_dequant;
      block_min_cos = 1.0f / sqrtf(1.0f + block_radius * block_radius);
    }
  }

  HS_PROFILE_DEEP(raster_scan);
  for (int y = y_lo; y <= y_hi; ++y) {
    if (per_row) {
      if (y == y_lo ||
          !shape.template horizontal_intervals_equal_rows<W, H>(runs_y, y)) {
        rebuild_runs(y);
        runs_y = y;
      }
      if (num_runs == 0)
        continue;
    }
    float sp = math::TrigLUT<W, H>::sin_phi[y];
    float cp = math::TrigLUT<W, H>::cos_phi[y];
    const int stride = pole_lod_run(sp);
    float block_slack = 0.0f;
    float reject_dsq = base_reject_dsq;
    if constexpr (pole_lod_blocks<SDF::Face>) {
      if constexpr (COLUMN_TIME)
        block_slack = (pole_lod_slack<W, SDF::Face>(stride, sp) +
                       column_rotation * static_cast<float>(stride - 1) *
                           SDF::arc_stretch<SDF::Face>)*plane_stretch;
      else
        block_slack = pole_lod_slack<W, SDF::Face>(stride, sp) * plane_stretch;
      // the block test reads d, so the cull must not sentinel a probe inside
      // its threshold
      const float reject_rad = (PIXEL_WIDTH + block_slack) * 1.0002f;
      reject_dsq = reject_rad * reject_rad;
    }

    for (size_t r = 0; r < num_runs; ++r) {
      const int rx2 = runs[r].second;
      [[maybe_unused]] int next_block = runs[r].first;
      if constexpr (pole_lod_blocks<SDF::Face>)
        if (stride > 1)
          next_block = pole_lod_block_anchor(next_block, stride);
      // A column-time run lies on one arm: the far arm reads the near arm's
      // entry with sin(phi) negated.
      [[maybe_unused]] int arm_offset = 0;
      [[maybe_unused]] float arm_sp = sp;
      if constexpr (COLUMN_TIME) {
        if (runs[r].first >= W / 2) {
          arm_offset = W / 2;
          arm_sp = -sp;
        }
      }
      for (int x = runs[r].first; x < rx2;) {
        math::Vector p;
        if constexpr (COLUMN_TIME) {
          const math::Vector &e = azimuth[x - arm_offset];
          const math::Vector &v = up[x - arm_offset];
          p = math::Vector(arm_sp * e.x + cp * v.x, arm_sp * e.y + cp * v.y,
                           arm_sp * e.z + cp * v.z);
        } else {
          p = math::Vector(sp * cos_theta[x], cp, sp * sin_theta[x]);
        }
        float probe_min_cos = min_cos;
        [[maybe_unused]] bool block_fits = false;
        if constexpr (pole_lod_blocks<SDF::Face>) {
          block_fits = x + stride <= rx2;
          if (stride > 1 && x == next_block && block_fits &&
              probe_bounds_block<SDF::Face>(PIXEL_WIDTH, block_slack))
            probe_min_cos = block_min_cos;
        }
        shape.template distance_with_flags<true>(p, res, reject_dsq,
                                                 probe_flags, probe_min_cos);
        const float d = res.dist;

        // Columns this one shade covers: a canvas-aligned block within the run
        // that the surface cannot cross; edge-holding blocks stay per-column.
        int span = 1;
        if constexpr (pole_lod_blocks<SDF::Face>) {
          if (stride > 1 && x == next_block) {
            next_block += stride;
            if (block_fits && pole_lod_settle_with_slack<true>(
                                  d, PIXEL_WIDTH, block_slack, shape, &p) !=
                                  BlockSettle::PER_COLUMN)
              span = stride;
          }
        }

        if (d >= PIXEL_WIDTH) {
          x += span;
          continue;
        }
        HS_SCAN_METRIC(hs::g_scan_metrics.shade_candidates++);
        HS_PROBE_COUNT(n_alpha);
        HS_PROBE_MARK(hs_ta);
        const float alpha = solid_coverage(d, PIXEL_WIDTH);
        HS_PROBE_SPAN(alpha, hs_ta);
        if (alpha <= MIN_ALPHA) {
          x += span;
          continue;
        }
        HS_PROFILE_DEEP(raster_shade);
        if constexpr (MinimalFragment) {
#ifndef NDEBUG
          frag.poison_inputs();
#endif
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
        const auto plot_column = [&](int px, float a) HS_HOT_INLINE {
          if constexpr (requires {
                          pipeline.plot_in_bounds(
                              canvas, px, y, frag.color.color, frag.age, a);
                        }) {
            pipeline.plot_in_bounds(canvas, px, y, frag.color.color, frag.age,
                                    a);
          } else {
            pipeline.plot(canvas, px, y, frag.color.color, frag.age, a);
          }
        };
        const auto shade_once = [&]() HS_HOT_INLINE {
          fragment_shader(p, frag);
          if (frag.color.alpha > MIN_ALPHA) {
            const float a = frag.color.alpha * alpha;
            for (int px = x; px < x + span; ++px)
              plot_column(px, a);
          }
        };
        if constexpr (TIMED_SHADER) {
          if (per_column) {
            for (int px = x; px < x + span; ++px) {
              fragment_shader(p, frag,
                              static_cast<float>(px - arm_offset) / (W / 2));
              if (frag.color.alpha > MIN_ALPHA)
                plot_column(px, frag.color.alpha * alpha);
            }
          } else {
            shade_once();
          }
        } else {
          shade_once();
        }
        x += span;
      }
    }
  }
}

/**
 * @brief Face-specialized scan loop for the mesh path at the fixed meridian,
 *        compiled at -O3 on the -Os device image and placed in FLASH.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam PipelineT Plotting pipeline type.
 * @tparam MinimalFragment When true, only v1 (raw distance) is refreshed per
 *         pixel. Other input registers are NaN in debug builds and retain
 *         defaults or previous values in release builds. The shader must not
 *         read inputs it did not set. It must write color on every call; age
 *         retains its default or the value left by the shader.
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Face to rasterize.
 * @param fragment_shader Shader invoked per covered pixel.
 * @details Faces whose column intervals are row-independent build their runs
 * once. The per-pixel body mirrors process_pixel's solid path. Does not read
 * canvas.debug().
 */
template <int W, int H, typename PipelineT, bool MinimalFragment = false,
          typename FragmentShaderT>
HS_HOT_FLASH_MEMBER void rasterize_face(PipelineT &pipeline, Canvas &canvas,
                                        const SDF::Face &shape,
                                        FragmentShaderT &fragment_shader) {
  rasterize_face_body<W, H, PipelineT, MinimalFragment, false, true>(
      pipeline, canvas, shape, fragment_shader, nullptr);
}

/**
 * @brief rasterize_face along a rotated basis.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam PipelineT Plotting pipeline type.
 * @tparam MinimalFragment See rasterize_face.
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Face built at the current endpoint, with swept bounds.
 * @param fragment_shader Shader invoked per covered pixel.
 * @param time Column-time state with a rotated basis.
 */
template <int W, int H, typename PipelineT, bool MinimalFragment = false,
          typename FragmentShaderT>
HS_NOINLINE_NOCLONE inline void rasterize_face_column_time(
    PipelineT &pipeline, Canvas &canvas, const SDF::Face &shape,
    FragmentShaderT &fragment_shader, const Temporal::ColumnTime &time) {
  rasterize_face_body<W, H, PipelineT, MinimalFragment, true, false>(
      pipeline, canvas, shape, fragment_shader, &time);
}
HS_O3_END

/**
 * @brief Fills the per-row face azimuth pads for faces swept by an arc.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param pads H entries receiving the pads (radians).
 * @param sweep Largest arc any face point moves over the frame (radians).
 * @details The AA reach and the sweep add as arcs; a row closer to a pole
 * than their sum spans the whole row.
 */
template <int W, int H>
HS_FLASH_MEMBER void fill_swept_face_pads(float *pads, float sweep) {
  const float reach = SDF::face_azimuth_pad(W) + sweep;
  for (int y = 0; y < H; ++y) {
    const float sin_phi = math::TrigLUT<W, H>::sin_phi[y];
    pads[y] = sin_phi <= reach ? math::PI_F : asinf(reach / sin_phi);
  }
}

/**
 * @brief Fills the per-row face azimuth pads of a mesh scan.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param pads H entries receiving the pads (radians).
 * @param sweep Largest arc any face point moves over the frame (radians); 0
 *        for the fixed meridian.
 */
template <int W, int H> HS_O3_FN void fill_face_pads(float *pads, float sweep) {
  if (sweep > 0.0f) {
    fill_swept_face_pads<W, H>(pads, sweep);
  } else {
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    for (int y = 0; y < H; ++y)
      pads[y] = SDF::face_azimuth_pad(W, math::TrigLUT<W, H>::sin_phi[y]);
  }
}

/**
 * @brief Arena-hosts one per-mesh Face scratch buffer.
 * @param arena Scratch arena supplying the storage.
 * @return The constructed buffer, valid until the caller's ScratchScope closes.
 * @details Default-initialized; SDF::Face writes every field before reading it.
 */
HS_NOINLINE_NOCLONE inline SDF::FaceScratchBuffer *
new_face_scratch(Arena &arena) {
  return arena.make_default_initialized<SDF::FaceScratchBuffer>();
}

/// Slab face-column entry of a face that holds over the frame.
inline constexpr uint16_t HELD_FACE_COLUMN = 0xFFFF;

/**
 * @brief Canvas columns that can see a face over one frame.
 * @tparam W Canvas width in pixels.
 * @param vertices Mesh vertex pool, at the current endpoint.
 * @param indices The face's vertex indices.
 * @param count Number of face vertices.
 * @param reach Largest arc the face or a column's basis moves over the frame
 *        (radians).
 * @param columns Receives the first seeing column, in [0, W), and how many
 *        more follow it, wrapping past the last column.
 * @details Judged by the cap about the face's vertex mean that holds every
 * vertex, widened by @p reach; a cap past a hemisphere or holding a pole is
 * seen by every column.
 */
template <int W>
HS_FLASH_MEMBER void face_lit_columns(const math::Vector *vertices,
                                      const uint16_t *indices, size_t count,
                                      float reach, uint16_t *columns) {
  columns[0] = 0;
  columns[1] = W - 1;
  math::Vector sum(0.0f, 0.0f, 0.0f);
  for (size_t k = 0; k < count; ++k)
    sum = sum + vertices[indices[k]];
  const float length = sqrtf(math::dot(sum, sum));
  if (!(length > math::TOLERANCE))
    return;
  const math::Vector axis = sum * (1.0f / length);
  float min_cos = 1.0f;
  for (size_t k = 0; k < count; ++k) {
    const math::Vector &v = vertices[indices[k]];
    min_cos = fminf(min_cos, math::dot(axis, v) / sqrtf(math::dot(v, v)));
  }
  // A cap past a hemisphere does not bound its polygon.
  if (!(min_cos > 0.0f))
    return;
  int lo, hi;
  if (!lit_span<W>(axis, acosf(min_cos), reach, lo, hi) || hi - lo >= W - 1)
    return;
  columns[0] = static_cast<uint16_t>(((lo % W) + W) % W);
  columns[1] = static_cast<uint16_t>(hi - lo);
}

/**
 * @brief Whether a face's seeing columns meet a sweep slab's columns inside
 *        the render clip.
 * @tparam W Canvas width in pixels.
 * @param face The face's seeing columns, as face_lit_columns() gives them.
 * @param begin The slab's first half-turn column.
 * @param end One past the slab's last half-turn column.
 * @param xc Render clip.
 * @return True when a column of the slab, on either arm, sees the face and
 *         the clip admits it.
 */
template <int W>
HS_FLASH_MEMBER bool face_columns_meet(const uint16_t *face, int begin, int end,
                                       ClipRegion::XClip xc) {
  constexpr int HALF = W / 2;
  const ClipPieces clip = clip_pieces<W>(xc);
  const int start = face[0];
  const int stop = start + face[1] + 1;
  // The face's columns past the last one wrap to the first.
  for (const int shift : {0, W})
    for (const int arm : {0, HALF}) {
      const int lo = std::max(start - shift, begin + arm);
      const int hi = std::min(stop - shift, end + arm);
      for (int i = 0; i < clip.count; ++i)
        if (std::max(lo, clip.lo[i]) < std::min(hi, clip.hi[i]))
          return true;
    }
  return false;
}

/// Largest held vertex motion over a frame, in pixels of arc.
inline constexpr float HELD_VERTEX_PIXELS = 0.5f;

/// Arc a vertex moves over a frame, radians, up to which a deformation scan
/// holds it at the current endpoint.
template <int W>
inline constexpr float HELD_VERTEX_ARC =
    HELD_VERTEX_PIXELS * math::TWO_PI_F / W;

/**
 * @brief Flags the vertices that move between two endpoints.
 * @tparam W Canvas width in pixels.
 * @param previous Vertices at the previous endpoint.
 * @param current The same vertices at the current endpoint.
 * @param count Number of vertices.
 * @param moving Receives 1 for a vertex that moves more than
 *        HELD_VERTEX_ARC, else 0.
 */
template <int W>
HS_FLASH_MEMBER void flag_moving_vertices(const math::Vector *previous,
                                          const math::Vector *current,
                                          size_t count, uint8_t *moving) {
  const float held_cos = cosf(HELD_VERTEX_ARC<W>);
  for (size_t i = 0; i < count; ++i) {
    const float norms =
        math::dot(previous[i], previous[i]) * math::dot(current[i], current[i]);
    moving[i] = math::dot(previous[i], current[i]) < held_cos * sqrtf(norms);
  }
}

/**
 * @brief Scratch bytes a deformation scan's vertex and face tables take.
 * @param vertices Mesh vertex count.
 * @param faces Mesh face count.
 * @return The bytes of a per-vertex moving flag, pick_moving_faces() columns
 *         and deform_mesh() tables, with alignment slack.
 */
constexpr size_t mesh_motion_bytes(size_t vertices, size_t faces) {
  return 2 * vertices + 5 * faces + 64;
}

/// Vertex motion over one frame, in pixels of arc, past which the frame is a
/// jump and draws at its current endpoint.
inline constexpr float MAX_VERTEX_MOTION_PIXELS = 8.0f;

/**
 * @brief Largest arc any vertex moves between two endpoints.
 * @param previous Vertices at the previous endpoint.
 * @param current The same vertices at the current endpoint.
 * @param count Number of vertices.
 * @return The largest angle between a vertex's two positions (radians).
 */
HS_FLASH_MEMBER inline float max_vertex_arc(const math::Vector *previous,
                                            const math::Vector *current,
                                            size_t count) {
  float min_cos = 1.0f;
  for (size_t i = 0; i < count; ++i) {
    const float norms =
        math::dot(previous[i], previous[i]) * math::dot(current[i], current[i]);
    if (norms > 0.0f)
      min_cos =
          fminf(min_cos, math::dot(previous[i], current[i]) / sqrtf(norms));
  }
  return acosf(hs::clamp(min_cos, -1.0f, 1.0f));
}

/**
 * @brief Whether a frame's mesh motion uses capture-time deformation.
 * @tparam W Display width.
 * @param arc Largest vertex step over the frame, in radians.
 * @return True with temporal sampling enabled for a non-held, non-jump step.
 */
template <int W> HS_FLASH_MEMBER bool mesh_motion_active(float arc) {
  const float pixels = arc * (W / math::TWO_PI_F);
  return Temporal::temporal_sampling_enabled() && pixels > HELD_VERTEX_PIXELS &&
         pixels <= MAX_VERTEX_MOTION_PIXELS;
}

/**
 * @brief Recomputes the seeing columns of a deformation scan's moving faces.
 * @tparam W Canvas width in pixels.
 * @param canvas Canvas with the frame's column-time state installed.
 * @param mesh The scanned mesh, at the current endpoint.
 * @param reach Bound on how far a placed vertex strays from the current
 *        endpoint over the frame (radians).
 * @param columns Per face two entries, HELD_FACE_COLUMN first for a face
 *        that holds; the others receive the columns that can see the face,
 *        as face_lit_columns() gives them.
 */
template <int W>
HS_FLASH_MEMBER void moving_face_columns(const Canvas &canvas,
                                         const MeshState &mesh, float reach,
                                         uint16_t *columns) {
  const float sweep = canvas.column_time() ? canvas.column_time()->sweep : 0.0f;
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  for (size_t f = 0; f < mesh.num_faces(); ++f)
    if (columns[2 * f] != HELD_FACE_COLUMN)
      face_lit_columns<W>(mesh.vertices.data(), indices + offsets[f], counts[f],
                          sweep + reach, columns + 2 * f);
}

/**
 * @brief Picks the faces whose vertices move during a scan.
 * @tparam W Canvas width in pixels.
 * @param canvas Canvas with the frame's column-time state installed.
 * @param mesh The scanned mesh, at the current endpoint.
 * @param moving_vertices Nonzero per moving vertex; left nonzero only on the
 *        vertices of the picked faces.
 * @param reach Bound on how far a placed vertex strays from the current
 *        endpoint over the frame (radians).
 * @param columns Receives per face two entries: the columns that can see a
 *        picked face, as face_lit_columns() gives them, or HELD_FACE_COLUMN
 *        and 0 for the others.
 * @return The number of picked faces.
 * @details A face is picked when it holds a moving vertex and a column the
 * render clip admits can see it; the others draw once at the current
 * endpoint.
 */
template <int W>
HS_FLASH_MEMBER size_t pick_moving_faces(const Canvas &canvas,
                                         const MeshState &mesh,
                                         uint8_t *moving_vertices, float reach,
                                         uint16_t *columns) {
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  const size_t faces = mesh.num_faces();
  for (size_t f = 0; f < faces; ++f) {
    uint8_t moves = 0;
    for (size_t k = 0; k < counts[f]; ++k)
      moves |= moving_vertices[indices[offsets[f] + k]];
    columns[2 * f] = moves ? 0 : HELD_FACE_COLUMN;
    columns[2 * f + 1] = 0;
  }
  moving_face_columns<W>(canvas, mesh, reach, columns);
  std::fill_n(moving_vertices, mesh.vertices.size(), uint8_t(0));
  size_t picked = 0;
  for (size_t f = 0; f < faces; ++f) {
    uint16_t *face = columns + 2 * f;
    if (face[0] == HELD_FACE_COLUMN)
      continue;
    if (!face_columns_meet<W>(face, 0, W / 2, canvas.render_x_clip())) {
      face[0] = HELD_FACE_COLUMN;
      face[1] = 0;
      continue;
    }
    ++picked;
    for (size_t k = 0; k < counts[f]; ++k)
      moving_vertices[indices[offsets[f] + k]] = 1;
  }
  return picked;
}

/**
 * @brief Returns to holding the picked faces none of whose vertices moves.
 * @param mesh The scanned mesh.
 * @param moving Nonzero per vertex that moves.
 * @param columns The faces' seeing columns, as pick_moving_faces() gives them;
 *        a dropped face's become HELD_FACE_COLUMN and 0.
 * @return The number of faces still picked.
 */
HS_FLASH_MEMBER inline size_t drop_held_faces(const MeshState &mesh,
                                              const uint8_t *moving,
                                              uint16_t *columns) {
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  size_t picked = 0;
  for (size_t f = 0; f < mesh.num_faces(); ++f) {
    if (columns[2 * f] == HELD_FACE_COLUMN)
      continue;
    uint8_t moves = 0;
    for (size_t k = 0; k < counts[f]; ++k)
      moves |= moving[indices[offsets[f] + k]];
    if (moves) {
      ++picked;
    } else {
      columns[2 * f] = HELD_FACE_COLUMN;
      columns[2 * f + 1] = 0;
    }
  }
  return picked;
}

/**
 * @brief Unwraps a capture fraction around a window's centre, then clamps it.
 * @param u Periodic arm fraction.
 * @param centre Window centre in arm fractions.
 * @return Capture fraction in [0, 1].
 */
HS_HOT_INLINE inline float unwrap_capture(float u, float centre) {
  return hs::clamp(u + roundf(centre - u), 0.0f, 1.0f);
}

/**
 * @brief Projects a vertex into the frame's moving column basis.
 * @tparam W Display width.
 * @param p Vertex in endpoint camera space.
 * @param time Frame's column basis.
 * @param u Fraction selecting the local column basis.
 * @return Periodic arm fraction of the projected vertex.
 */
template <int W>
HS_HOT_INLINE float capture_fraction(const math::Vector &p,
                                     const Temporal::ColumnTime &time,
                                     float u) {
  const int x = hs::clamp(static_cast<int>(u * (W / 2)), 0, W / 2 - 1);
  const math::Vector tangent = math::cross(time.azimuth[x], time.up[x]);
  const float fraction =
      static_cast<float>(x) / (W / 2) +
      atan2f(math::dot(p, tangent), math::dot(p, time.azimuth[x])) / math::PI_F;
  return fraction - floorf(fraction);
}

/**
 * @brief Evaluates one moving vertex at its capture time with one correction.
 * @tparam W Display width.
 * @tparam VertexFn Callable (float u) returning the vertex at that fraction.
 * @param end Vertex at the current endpoint.
 * @param time Frame's column basis.
 * @param reach Bound on vertex motion in radians.
 * @param centre Window centre in arm fractions.
 * @param vertex Pure evaluation of the vertex path.
 * @return Captured vertex, or its endpoint inside the polar hold band.
 */
template <int W, typename VertexFn>
HS_HOT_INLINE math::Vector
capture_vertex(const math::Vector &end, const Temporal::ColumnTime &time,
               float reach, float centre, VertexFn &&vertex) {
  const float horizontal = sqrtf(end.x * end.x + end.z * end.z);
  const float bound = reach + time.sweep;
  if (horizontal <= sinf(bound) || bound / horizontal >= math::PI_F / 4.0f)
    return end;
  const float angle = atan2f(end.z, end.x) / math::PI_F;
  float u = unwrap_capture(angle - floorf(angle), centre);
  const math::Vector first = vertex(u);
  u = unwrap_capture(capture_fraction<W>(first, time, u), u);
  return vertex(u);
}

/**
 * @brief Scans shared capture-time vertices through two disjoint arm windows.
 * @tparam W Display width.
 * @tparam H Display height.
 * @tparam VertexFn Callable (size_t i, float u) evaluating an endpoint-space vertex.
 * @tparam ScanFn Callable () scanning the mesh into the canvas's open window.
 * @param canvas Canvas with column-time state when deformation is active.
 * @param mesh Endpoint mesh, overwritten by the captured vertices.
 * @param active Whether the mesh deforms during this frame.
 * @param columns Per-face seeing columns; HELD_FACE_COLUMN marks held faces.
 * @param moving Nonzero for vertices whose paths are captured; others hold.
 * @param scratch Arena for face selection and vertex placement flags.
 * @param reach Bound on vertex motion in radians.
 * @param vertex Pure evaluation of a vertex at a frame fraction.
 * @param scan Scans the mesh with its current vertices and selected faces.
 * @details Held faces draw once. Each window uses a common vertex pool for
 * all its moving faces; faces outside that window's seeing columns are skipped.
 */
template <int W, int H, typename VertexFn, typename ScanFn>
HS_HOT_FLASH_MEMBER void
deform_mesh(Canvas &canvas, MeshState &mesh, bool active,
            const uint16_t *columns, const uint8_t *moving, Arena &scratch,
            float reach, VertexFn &&vertex, ScanFn &&scan) {
  if (!active) {
    scan();
    return;
  }
  HS_CHECK(canvas.column_time(),
           "mesh deformation needs the frame's column-time state installed");
  ScratchScope guard(scratch);
  const size_t faces = mesh.num_faces(), vertices = mesh.vertices.size();
  uint8_t *selected = scratch.allocate_n<uint8_t>(faces);
  uint8_t *placed = scratch.allocate_n<uint8_t>(vertices);
  Temporal::ColumnTime time = *canvas.column_time();
  time.selected_faces = selected;
  float *face_pads = scratch.allocate_n<float>(H);
  fill_face_pads<W, H>(face_pads, time.azimuth ? time.sweep : 0.0f);
  time.face_pads = face_pads;
  ColumnTimeScope scope(canvas, time);
  const auto *counts = mesh.get_face_counts_data();
  const auto *indices = mesh.get_faces_data();
  const auto *offsets = mesh.get_face_offsets_data();
  bool held = false;
  for (size_t f = 0; f < faces; ++f) {
    const bool holds = columns[2 * f] == HELD_FACE_COLUMN;
    selected[f] =
        holds ? Temporal::MESH_FACE_HELD : Temporal::MESH_FACE_SKIPPED;
    held |= holds;
  }
  if (held)
    scan();
  for (int pass = 0; pass < 2; ++pass) {
    const int begin = pass * (W / 4), end = (pass + 1) * (W / 4);
    std::fill_n(placed, vertices, uint8_t(0));
    bool any = false;
    for (size_t f = 0; f < faces; ++f) {
      const bool draw = columns[2 * f] != HELD_FACE_COLUMN &&
                        face_columns_meet<W>(columns + 2 * f, begin, end,
                                             canvas.render_x_clip());
      selected[f] =
          draw ? Temporal::MESH_FACE_WINDOW : Temporal::MESH_FACE_SKIPPED;
      if (!draw)
        continue;
      any = true;
      for (size_t k = 0; k < counts[f]; ++k)
        placed[indices[offsets[f] + k]] = 1;
    }
    if (!any)
      continue;
    const float centre = pass == 0 ? 0.25f : 0.75f;
    for (size_t i = 0; i < vertices; ++i) {
      if (!placed[i] || !moving[i])
        continue;
      const math::Vector endpoint =
          pass == 0 ? mesh.vertices[i] : vertex(i, 1.0f);
      mesh.vertices[i] = capture_vertex<W>(
          endpoint, time, reach, centre,
          [&](float u) HS_HOT_INLINE { return vertex(i, u); });
    }
    RenderXWindow window(canvas, begin, end);
    scan();
  }
}

HS_O3_BEGIN
/**
 * @brief Rasterizes a polygonal mesh by drawing each face as an SDF::Face.
 * @details Without a face-shader setup callback the face index reaches the
 * shader in register v2; with one, it reaches the shader only through what
 * the callback hoists.
 * canvas.debug() does not tint a mesh. Under the canvas's column-time basis
 * each face's bounds are swept by its arc and every pixel is evaluated along
 * its column's rotated direction. An open RenderXWindow limits the scan to its
 * columns; in a deformation scan, Temporal::ColumnTime::selected_faces selects
 * the faces each pass draws.
 */
struct Mesh {
  /**
   * @brief Traps unless every flat face index addresses a vertex in the pool.
   * @param faces Flat per-face vertex index array.
   * @param num_indices Entries in @p faces.
   * @param num_verts Vertices the mesh carries.
   * @details SDF::Face indexes the vertex pool through an unchecked std::span.
   */
  HS_NOINLINE_NOCLONE
  static void check_face_index_domain(const uint16_t *faces, size_t num_indices,
                                      size_t num_verts) {
    for (size_t k = 0; k < num_indices; ++k)
      HS_CHECK(static_cast<size_t>(faces[k]) < num_verts,
               "mesh face index exceeds the vertex pool");
  }

  /**
   * @brief Rasterizes every face of a mesh.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type (deduced; pass PipelineRef to erase the sink).
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @param fragment_shader Shader invoked per covered pixel; without
   *                        @p face_shader_setup it receives the face index in
   *                        register v2. v2 is a float, so the index is exact
   *                        only up to 2^24 faces.
   * @param scratch_arena Arena supplying per-face SDF::Face scratch storage.
   * @param bake Optional congruence-class bake for this mesh (null = exact
   *        distance path for every face). When present, each face is
   *        aligned to its canonical class shape after construction and the
   *        class distance LUT is bound for the probe loop.
   * @param face_shader_setup Optional callback receiving the face index and
   *        size, run once before each face. When supplied, the shader is
   *        invoked directly with only v1 refreshed per pixel, so the index and
   *        size reach it only through what the callback hoists. The fragment
   *        persists across pixels, so the shader must write frag.color
   *        unconditionally.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_NOINLINE_NOCLONE static void
  draw_impl(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
            FragmentShaderT &fragment_shader, Arena &scratch_arena,
            const MeshOps::MeshClassBake *bake,
            FaceShaderSetupT &face_shader_setup) {
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    // the per-face wrapper is never null; check the erased shader it wraps
    if constexpr (std::is_same_v<FragmentShaderT, FragmentShaderFn>)
      check_fragment_shader(fragment_shader);

    ScratchScope scope(scratch_arena);
    auto *scratch = new_face_scratch(scratch_arena);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const Temporal::ColumnTime *time = column_time(pipeline, canvas);
    const float sweep = time && time->azimuth ? time->sweep : 0.0f;
    const float *azimuth_pads = time ? time->face_pads : nullptr;
    if (!azimuth_pads) {
      float *pads = scratch_arena.allocate_n<float>(H);
      fill_face_pads<W, H>(pads, sweep);
      azimuth_pads = pads;
    }
    const float bounds_margin =
        std::max(SDF::BOUNDS_MARGIN, math::TWO_PI_F / W) + sweep;

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    const uint16_t *fo = mesh.get_face_offsets_data();
    size_t fi_size = mesh.get_faces_size();

    HS_CHECK(mesh.get_face_offsets_size() == num_f,
             "solid mesh scan requires one face offset per face");

    check_face_index_domain(fi, fi_size, mesh.vertices.size());

    // An unbound bake is treated as none; records are indexed by face order.
    if (bake && !bake->face_recs.is_bound())
      bake = nullptr;
    HS_CHECK(!bake || bake->face_recs.size() == num_f,
             "mesh class bake face count disagrees with the mesh");

    const auto &cr = source_clip<W, H>(pipeline, canvas);
    const ColumnWindow &window = canvas.render_x_window();
    const uint8_t *selected_faces =
        canvas.column_time() ? canvas.column_time()->selected_faces : nullptr;
    const uint8_t drawn_face =
        window.active ? Temporal::MESH_FACE_WINDOW : Temporal::MESH_FACE_HELD;
    for (size_t i = 0; i < num_f; ++i) {
      size_t count = fc[i];

      HS_CHECK(static_cast<size_t>(fo[i]) + count <= fi_size,
               "mesh face span exceeds face index array");
      if (selected_faces && selected_faces[i] != drawn_face)
        continue;

      std::span<const math::Vector> verts(mesh.vertices.data(),
                                          mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      // IIFE scopes HS_PROFILE to Face construction (prvalue elided in place).
      SDF::Face shape = [&] {
        HS_PROFILE(scan_face_setup);
        return SDF::Face(verts, indices, *scratch, math::LatitudeGeometry(H), H,
                         &cr, azimuth_pads, bounds_margin);
      }();

      if (window.active &&
          (shape.y_min > shape.y_max || !face_meets_window(shape, window, cr)))
        continue;

      // Bind the face's congruence-class LUT; culled faces (y_min > y_max) skip
      // it.
      if (bake && shape.y_min <= shape.y_max) {
        const MeshOps::FaceClassRec &rec = bake->face_recs[i];
        if (rec.class_id != MeshOps::NO_CLASS) {
          HS_CHECK(rec.class_id < bake->classes.size(),
                   "mesh class bake face record names an unknown class");
          const MeshOps::CongruenceClass &cls = bake->classes[rec.class_id];
          if (cls.lut.data && cls.n_verts == shape.count)
            shape.bind_class_lut(&cls.lut, cls.canon_xy, rec.vert_offset,
                                 rec.reflected != 0);
        }
      }

      if constexpr (std::is_same_v<FaceShaderSetupT, std::nullptr_t>) {
        auto wrapper = [&](const math::Vector &p, Fragment &f_in) {
          // v2 carries the face index (decoded by mesh_face_index()).
          f_in.v2 = static_cast<float>(i);
          fragment_shader(p, f_in);
        };
        FragmentShaderFn raster_shader = wrapper;
        {
          HS_PROFILE(scan_mesh_raster);
          if (sweep > 0.0f)
            rasterize_face_column_time<W, H, PipelineT, false>(
                pipeline, canvas, shape, raster_shader, *time);
          else
            rasterize_face<W, H, PipelineT, false>(pipeline, canvas, shape,
                                                   raster_shader);
        }
      } else {
        face_shader_setup(i, shape.size);
        {
          HS_PROFILE(scan_mesh_raster);
          if (sweep > 0.0f)
            rasterize_face_column_time<W, H, PipelineT, true>(
                pipeline, canvas, shape, fragment_shader, *time);
          else
            rasterize_face<W, H, PipelineT, true>(pipeline, canvas, shape,
                                                  fragment_shader);
        }
      }
    }
  }

  /**
   * @brief Rasterizes a mesh through the type-erased fragment shader path.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @param fragment_shader Shader invoked per covered pixel; receives the
   *        face index in register v2.
   * @param scratch_arena Arena supplying per-face SDF::Face scratch storage.
   * @param bake Optional congruence-class bake for this mesh.
   */
  template <int W, int H, typename PipelineT>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake = nullptr) {
    std::nullptr_t face_shader_setup = nullptr;
    draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                    bake, face_shader_setup);
  }

  /**
   * @brief Rasterizes a mesh with an inlinable per-face fragment shader.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh Mesh providing vertices, face counts, indices, and offsets.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT Callable `void(size_t face, float size)`.
   * @param fragment_shader Shader invoked per covered pixel; must write
   *        frag.color unconditionally.
   * @param scratch_arena Arena supplying per-face SDF::Face scratch storage.
   * @param bake Optional congruence-class bake for this mesh.
   * @param face_shader_setup Non-null callback run once before each face.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  static void
  draw_specialized(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderT &fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake,
                   FaceShaderSetupT &face_shader_setup) {
    FunctionRef<void(size_t, float)> erased_setup = face_shader_setup;
    // the erasure hides a null setup from draw_impl's nullptr_t branch
    HS_CHECK(erased_setup, "Scan::Mesh::draw_specialized requires a non-null "
                           "face_shader_setup");
    draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                    bake, erased_setup);
  }

  /**
   * @brief Rasterizes a mesh shaded from a per-face blended-ramp table.
   * @tparam ShadingT Table exposing `const BakedPalette &ramp_for(size_t)`.
   * @param pipeline Filter pipeline applied on plot.
   * @param canvas Target canvas.
   * @param mesh Mesh to rasterize, already in camera space.
   * @param shading Per-face ramp table (an Animation::OpLeg::Shading).
   * @param gain Edge-distance shading gain, divided by the face size.
   * @param alpha Coverage applied to every fragment.
   * @param scratch_arena Arena backing the scan's per-face buffers.
   * @details The ramp and gradient scale resolve once per face. A table
   * carrying its crossfade endpoints (`bank` and `blend_face`) moves each
   * column-time column along the crossfade from the weight drawn a frame ago.
   */
  template <int W, int H, typename PipelineT, typename ShadingT>
  static void draw_opleg_shading(PipelineT &pipeline, Canvas &canvas,
                                 const MeshState &mesh, const ShadingT &shading,
                                 float gain, float alpha,
                                 Arena &scratch_arena) {
    FacePaletteShader fragment_shader;
    fragment_shader.alpha = alpha;
    auto select_face = [&](size_t fi, float size) {
      fragment_shader.set_palette(&shading.ramp_for(fi));
      if constexpr (requires { shading.blend_face(fragment_shader, fi); }) {
        if (shading.bank)
          shading.blend_face(fragment_shader, fi);
      }
      fragment_shader.scale = size > math::TOLERANCE ? gain / size : 0.0f;
    };
    draw_specialized<W, H>(pipeline, canvas, mesh, fragment_shader,
                           scratch_arena, nullptr, select_face);
  }
};
HS_O3_END

} // namespace Scan
