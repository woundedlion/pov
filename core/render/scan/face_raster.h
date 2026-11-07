/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstring>

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
#include "core/render/scan/swept_face.h"
#include "core/render/spacetime/crossing.h"

/**
 * @file face_raster.h
 * @brief Face-specialized scan loops and swept face bounds.
 */

namespace Scan {

#if HS_ENABLE_TEST_HOOKS
/** @brief Candidate observer for mesh replay verification. */
inline void (*face_audit_candidate)(size_t, int, int, float, bool) = nullptr;
/** @brief Face currently visited by the candidate observer. */
inline size_t face_audit_index = 0;
/** @brief Total faces, selected faces, total vertices and selected vertices, in order. */
inline size_t visibility_counts[4] = {};
/** @brief Deformation reach used by the last visibility pre-pass. */
inline float visibility_audit_reach = 0.0f;
/** @brief Disables structural culling for a replay reference. */
inline bool face_audit_reference = false;

#endif

/**
 * @brief Checks installed mesh tables against their consumer's dimensions.
 * @param pass Installed mesh pass.
 * @param faces Consumer face count.
 * @param vertices Consumer vertex count.
 */
inline void check_mesh_pass(const detail::MeshPass &pass, size_t faces,
                            size_t vertices) {
  HS_CHECK((!pass.visible_faces && !pass.visible_vertices) ||
               (pass.visibility_face_count == faces &&
                pass.visibility_vertex_count == vertices),
           "mesh visibility dimensions disagree with the mesh");
  HS_CHECK((!pass.selected_faces && !pass.face_copy && !pass.face_sweep) ||
               (pass.selection_face_count == faces &&
                pass.selection_vertex_count == vertices),
           "mesh deformation dimensions disagree with the mesh");
}

/** @brief Half-open column run [first, second); left uninitialized. */
struct ColumnRun {
  int first;  ///< First column.
  int second; ///< One past the last column.
};

/** @brief One face's projection of a column's rotating meridian basis. */
struct FaceProjectionColumn {
  float azimuth[3]; ///< Center, basis_u and basis_w dots with azimuth.
  float up[3];      ///< Center, basis_u and basis_w dots with up.
  uint32_t claim;   ///< Scratch claim identifying the projected face.
  int sector;       ///< Last angular sector visited down this column.

  /**
   * @brief Projects a column basis into a face's tangent frame.
   * @param face Face whose frame is projected.
   * @param e Column azimuth.
   * @param v Column up.
   */
  HS_HOT_INLINE void prepare(const SDF::PreparedFace &face,
                             const math::Vector &e, const math::Vector &v) {
    azimuth[0] = math::dot(e, face.center);
    azimuth[1] = math::dot(e, face.basis_u);
    azimuth[2] = math::dot(e, face.basis_w);
    up[0] = math::dot(v, face.center);
    up[1] = math::dot(v, face.basis_u);
    up[2] = math::dot(v, face.basis_w);
    claim = face.scratch_claim;
    sector = -1;
  }
};

/**
 * @brief Allocates column projection storage for one mesh draw.
 * @tparam W Canvas width in pixels.
 * @param arena Draw-scoped arena.
 * @return Storage with no face claiming a column.
 */
template <int W>
HS_HOT_FLASH_MEMBER FaceProjectionColumn *
new_face_projection_columns(Arena &arena) {
  auto *columns = arena.allocate_n<FaceProjectionColumn>(W / 2);
  for (int x = 0; x < W / 2; ++x)
    columns[x].claim = 0;
  return columns;
}

/**
 * @brief Whether a face's clipped row band can reach an arm of a window.
 * @param shape Constructed face.
 * @param window Render-x window.
 * @param clip Source clip defining the row band and display dimensions.
 * @return True when an admitted arm can receive a pixel, or the window is
 *         inactive.
 */
HS_HOT_FLASH_MEMBER inline bool
face_meets_window(const SDF::PreparedFace &shape,
                  const Spacetime::TimeWindow &window, const ClipRegion &clip) {
  if (!window.active)
    return true;
  const int arm = window.arm(clip.w);
  const int shift = arm < 0 ? 0 : (2 * arm - 1) * (clip.w / 2);
  ClipRegion side{clip.render_y_start(),
                  clip.render_y_end(),
                  window.begin - shift,
                  window.end - shift,
                  0,
                  clip.w,
                  clip.h};
  if (arm != 1 && !shape.clip_rejects(side))
    return true;
  if (arm == 0)
    return false;
  side.x_start += clip.w / 2;
  side.x_end += clip.w / 2;
  return !shape.clip_rejects(side);
}

/**
 * @brief Clips a row's column runs to a render-x window's arms.
 * @tparam W Canvas width in pixels.
 * @tparam CAPACITY Run buffer capacity.
 * @param runs The row's runs; receives the clipped runs, each on one arm.
 * @param count Number of runs.
 * @param window Open render-x window.
 * @return Number of clipped runs.
 */
template <int W, size_t CAPACITY>
HS_HOT_INLINE inline size_t window_runs(ColumnRun *runs, size_t count,
                                        const Spacetime::TimeWindow &window) {
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
  std::memcpy(runs, kept, kept_count * sizeof(ColumnRun));
  return kept_count;
}

/**
 * @brief Whether a row's column intervals can leave a run inside a render-x
 *        window.
 * @tparam W Canvas width in pixels.
 * @tparam IntervalBufT Buffer of the row's column intervals.
 * @param intervals The row's intervals, in whole columns, unwrapped.
 * @param window Open render-x window.
 * @return False when every interval misses the window's columns on both
 *         arms, so that the row's clipped runs are empty.
 * @pre The window admits both arms.
 */
template <int W, typename IntervalBufT>
HS_HOT_INLINE inline bool
intervals_meet_window(const IntervalBufT &intervals,
                      const Spacetime::TimeWindow &window) {
  const int length = window.end - window.begin;
  for (const auto &iv : intervals) {
    const int a = static_cast<int>(iv.start);
    // A zero-width interval still owns its column.
    const int len = std::max(static_cast<int>(iv.end) - a, 1);
    const int s = math::wrap(a, W);
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
 * @tparam FragmentShaderT Per-pixel shader type.
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Face to rasterize.
 * @param fragment_shader Shader invoked per covered pixel.
 * @param time Column-time state with a rotated basis; read only with
 *        COLUMN_TIME.
 * @param projections Per-column projection storage for COLUMN_TIME.
 * @details Under a rotated basis a pole-LOD block's slack grows by the
 * rotation between its first and last columns, and a block straddling an arm
 * boundary is walked per column. Only the columns an open TimeWindowScope
 * admits are written.
 */
template <int W, int H, typename PipelineT, bool MinimalFragment,
          bool COLUMN_TIME, bool INLINE_RUNS, typename FragmentShaderT>
HS_HOT_INLINE inline void rasterize_face_body(
    PipelineT &pipeline, Canvas &canvas, const SDF::PreparedFace &shape,
    FragmentShaderT &fragment_shader, const Scan::ColumnBasis *time,
    FaceProjectionColumn *projections) {
  const auto &cr = source_clip<W, H>(pipeline, canvas);
  const auto xc = cr.x_clip();

  int y_lo, y_hi;
  // A column-time or windowed run is split at the arm boundary, which one run
  // can cross.
  static constexpr size_t MAX_RUNS = 9;
  ColumnRun runs[MAX_RUNS];
  size_t num_runs = 0;
  const Spacetime::TimeWindow &window = canvas.time_window();

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
                      SDF::sdf_max_spans<SDF::PreparedFace>::value,
                  "intervals must hold a Face's per-row emission");
    static_assert(MAX_RUNS >=
                      4 * SDF::sdf_max_spans<SDF::PreparedFace>::value + 1,
                  "runs must hold every span seam-split, clip-split and "
                  "arm-split");
    bool handled = shape.template get_horizontal_intervals<W, H>(
        y, [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });
    // A one-arm window starts outside the half-turn.
    if (window.active && handled &&
        static_cast<unsigned>(2 * window.begin) < static_cast<unsigned>(W) &&
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
    column_rotation = time->sweep_rate / (W / 2);
  }
  constexpr float PIXEL_WIDTH = 2.0f * math::PI_F / W;
  const uint32_t probe_flags = shape.probe_flags();
  constexpr bool TIMED_SHADER =
      COLUMN_TIME &&
      Spacetime::reads_time<FragmentShaderT, const math::Vector &, Fragment &>;
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
  const float shade_reach = (probe_flags & SDF::PreparedFace::PROBE_LINEAR)
                                ? PIXEL_WIDTH
                                : tanf(PIXEL_WIDTH);
  float radial_cull = shape.radius + shade_reach;
  if (probe_flags & SDF::PreparedFace::PROBE_HAS_LUT)
    radial_cull += (shape.lut_q_safe + 1) * shape.lut_dequant;
  const float min_cos = 1.0f / sqrtf(1.0f + radial_cull * radial_cull);
  [[maybe_unused]] float block_min_cos = min_cos;
  if constexpr (pole_lod_blocks<SDF::PreparedFace>) {
    if (Render::pole_lod_aggressiveness > 0.0f) {
      float lod_reach = Render::pole_lod_aggressiveness * PIXEL_WIDTH;
      if constexpr (COLUMN_TIME)
        lod_reach += column_rotation * (Render::POLE_LOD_MAX_RUN - 1);
      const float block_reach = std::min(
          PIXEL_WIDTH +
              lod_reach * SDF::arc_stretch<SDF::PreparedFace> * plane_stretch,
          SDF::reject_margin<SDF::PreparedFace>);
      float block_radius =
          shape.radius + ((probe_flags & SDF::PreparedFace::PROBE_LINEAR)
                              ? block_reach
                              : tanf(block_reach));
      if (probe_flags & SDF::PreparedFace::PROBE_HAS_LUT)
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
    float reject_plane = base_reject_rad;
    if constexpr (pole_lod_blocks<SDF::PreparedFace>) {
      if constexpr (COLUMN_TIME)
        block_slack = (pole_lod_slack<W, SDF::PreparedFace>(stride, sp) +
                       column_rotation * static_cast<float>(stride - 1) *
                           SDF::arc_stretch<SDF::PreparedFace>)*plane_stretch;
      else
        block_slack =
            pole_lod_slack<W, SDF::PreparedFace>(stride, sp) * plane_stretch;
      // the block test reads d, so the cull must not sentinel a probe inside
      // its threshold
      const float reject_rad = (PIXEL_WIDTH + block_slack) * 1.0002f;
      reject_dsq = reject_rad * reject_rad;
      reject_plane = reject_rad;
    }

    for (size_t r = 0; r < num_runs; ++r) {
      const int rx2 = runs[r].second;
      [[maybe_unused]] int next_block = runs[r].first;
      if constexpr (pole_lod_blocks<SDF::PreparedFace>)
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
        if constexpr (pole_lod_blocks<SDF::PreparedFace>) {
          block_fits = x + stride <= rx2;
          if (stride > 1 && x == next_block && block_fits &&
              probe_bounds_block<SDF::PreparedFace>(PIXEL_WIDTH, block_slack))
            probe_min_cos = block_min_cos;
        }
        const float convex_reject =
            (probe_flags & SDF::PreparedFace::PROBE_LINEAR) ? reject_plane
                                                            : FLT_MAX;
        if constexpr (COLUMN_TIME) {
          const int column = x - arm_offset;
          FaceProjectionColumn &projection = projections[column];
          if (projection.claim != shape.scratch_claim)
            projection.prepare(shape, azimuth[column], up[column]);
          shape.template distance_projected<true>(
              arm_sp * projection.azimuth[0] + cp * projection.up[0],
              arm_sp * projection.azimuth[1] + cp * projection.up[1],
              arm_sp * projection.azimuth[2] + cp * projection.up[2], res,
              reject_dsq, probe_flags, probe_min_cos, convex_reject,
              &projection.sector);
        } else {
          shape.template distance_with_flags<true>(
              p, res, reject_dsq, probe_flags, probe_min_cos, convex_reject);
        }
        const float d = res.dist;
#if HS_ENABLE_TEST_HOOKS
        if (face_audit_candidate && d < PIXEL_WIDTH)
          face_audit_candidate(face_audit_index, x, y, d, shape.convex);
#endif

        // Columns this one shade covers: a canvas-aligned block within the run
        // that the surface cannot cross; edge-holding blocks stay per-column.
        int span = 1;
        if constexpr (pole_lod_blocks<SDF::PreparedFace>) {
          if (stride > 1 && x == next_block) {
            next_block += stride;
            if (block_fits && pole_lod_settle_with_slack<true>(
                                  d, PIXEL_WIDTH, block_slack, shape) !=
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
          load_fragment(frag, res, p, 0.0f);
        }
        const auto plot_column = [&](int px, float a) HS_HOT_INLINE {
          plot_admitted(pipeline, canvas, px, y, frag.color.color, frag.age, a);
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
              fragment_shader(p, frag, Spacetime::SamplingSurface<W>::time(px));
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
                                        const SDF::PreparedFace &shape,
                                        FragmentShaderT &fragment_shader) {
  rasterize_face_body<W, H, PipelineT, MinimalFragment, false, true>(
      pipeline, canvas, shape, fragment_shader, nullptr, nullptr);
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
 * @param projections Draw-scoped column projection storage.
 */
template <int W, int H, typename PipelineT, bool MinimalFragment = false,
          typename FragmentShaderT>
HS_NOINLINE_NOCLONE inline void rasterize_face_column_time(
    PipelineT &pipeline, Canvas &canvas, const SDF::PreparedFace &shape,
    FragmentShaderT &fragment_shader, const Scan::ColumnBasis &time,
    FaceProjectionColumn *projections) {
  rasterize_face_body<W, H, PipelineT, MinimalFragment, true, false>(
      pipeline, canvas, shape, fragment_shader, &time, projections);
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
  const float reach =
      math::TWO_PI_F / W + 2.0f * SDF::PreparedFace::AZIMUTH_ERROR + sweep;
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
 * @brief Arena-hosts one reusable face geometry and table slot.
 * @param arena Scratch arena supplying the storage.
 * @return Storage valid until the caller's ScratchScope closes.
 * @details Geometry is constructed separately for each admitted face.
 */
HS_NOINLINE_NOCLONE inline SDF::FaceStorage *new_face_storage(Arena &arena) {
  return arena.make_default_initialized<SDF::FaceStorage>();
}

/**
 * @brief Bound on the arc whose cosine is given.
 * @param c The cosine.
 * @return Radians, at least acos(c).
 */
HS_HOT_INLINE inline float cosine_arc(float c) {
  // Unit chord squared; below 1 the arc 2 asin(chord / 2) is at most
  // chord (1 + chord^2 / 16).
  const float chord2 = fmaxf(2.0f - 2.0f * c, 0.0f);
  if (chord2 < 1.0f)
    return sqrtf(chord2) * (1.0f + chord2 * (1.0f / 16.0f));
  return acosf(hs::clamp(c, -1.0f, 1.0f));
}

/// Maximum column span of a compact face pose.
inline constexpr int COMPACT_FACE_COLUMNS = 32;
/// Maximum column span of a wide face pose.
inline constexpr int WIDE_FACE_COLUMNS = 64;

/**
 * @brief Rotation bound across a compact face's seeing columns.
 * @tparam W Canvas width.
 * @param rate Frame rotation bound per frame fraction, as
 *        Scan::ColumnBasis::sweep_rate gives it.
 * @param columns Maximum seeing-column span.
 * @return Radians of residual rotation.
 */
template <int W>
HS_HOT_INLINE float compact_face_residual(float rate,
                                          int columns = COMPACT_FACE_COLUMNS) {
  return rate * (columns / 2 + 1.0f) / (W / 2) + rate * rate / math::PI_F;
}

/** @brief Shared row-padding tables for swept mesh faces. */
struct FacePadTables {
  const float *full;    ///< Full-motion pads.
  const float *compact; ///< Compact-pose pads.
  const float *wide;    ///< Wide-pose pads.

  /** @brief Reuses frame tables or prepares missing tables in the arena.
   * @tparam W,H Canvas dimensions.
   * @param arena Storage for missing tables.
   * @param time Column-time state, or null for a static frame.
   * @return Padding tables valid for the arena scope. */
  template <int W, int H>
  static FacePadTables prepare(Arena &arena, const Scan::ColumnBasis *time) {
    const float sweep = time && time->azimuth ? time->sweep_arc : 0.0f;
    FacePadTables result{time ? time->mesh.face_pads : nullptr,
                         time ? time->mesh.compact_face_pads : nullptr,
                         time ? time->mesh.wide_face_pads : nullptr};
    auto fill = [&](const float *&target, float motion) {
      if (!target) {
        float *filled = arena.allocate_n<float>(H);
        fill_face_pads<W, H>(filled, motion);
        target = filled;
      }
    };
    fill(result.full, sweep);
    if (sweep > 0.0f) {
      fill(result.compact, compact_face_residual<W>(time->sweep_rate));
      fill(result.wide,
           compact_face_residual<W>(time->sweep_rate, WIDE_FACE_COLUMNS));
    }
    return result;
  }
};

/** @brief Clip planes for a face's current vertex cap and rotation reach. */
template <int W, int H> struct FaceClipPlanes {
  math::Vector first; ///< First azimuth boundary normal.
  math::Vector last;  ///< Last azimuth boundary normal.
  float north;        ///< Northern admitted latitude cosine.
  float south;        ///< Southern admitted latitude cosine.
  float pad_cos;      ///< Cosine of the face's angular reach.
  float pad_sin;      ///< Sine of the face's angular reach.
  int columns;        ///< Admitted azimuth span.
  bool bounded;       ///< Reach leaves a bounding hemisphere.

  /**
   * @brief Prepares clip planes for one mesh scan pass.
   * @param clip Source render clip.
   * @param reach Rotation and AA reach in radians.
   */
  FaceClipPlanes(const ClipRegion &clip, float reach) {
    const auto xc = clip.x_clip();
    columns = xc.length(W);
    const int start = math::wrap(xc.rs, W);
    const int end = math::wrap(xc.rs + columns, W);
    first = {-math::TrigLUT<W, H>::sin_theta[start], 0.0f,
             math::TrigLUT<W, H>::cos_theta_data()[start]};
    last = {math::TrigLUT<W, H>::sin_theta[end], 0.0f,
            -math::TrigLUT<W, H>::cos_theta_data()[end]};
    const auto geometry = math::LatitudeGeometry(H);
    const float row = geometry.row_to_phi(1) - geometry.row_to_phi(0);
    north = cosf(fmaxf(geometry.row_to_phi(clip.render_y_start()) - row, 0.0f));
    south = cosf(
        fminf(geometry.row_to_phi(clip.render_y_end() - 1) + row, math::PI_F));
    bounded = reach < math::PI_F * 0.5f;
    pad_cos = cosf(reach);
    pad_sin = sinf(reach);
  }

  /**
   * @brief Rejects a cap missing every admitted pixel direction.
   * @param axis Unit cap axis.
   * @param cosine Cosine of the cap radius.
   * @return True when the padded cap misses the clip.
   */
  HS_HOT_INLINE bool rejects(const math::Vector &axis, float cosine) const {
    if (!bounded || cosine <= 0.0f)
      return false;
    const float sine = sqrtf(fmaxf(1.0f - cosine * cosine, 0.0f));
    const float expanded_cosine = cosine * pad_cos - sine * pad_sin;
    const float expanded_sine = sine * pad_cos + cosine * pad_sin;
    if (expanded_cosine <= 0.0f)
      return false;
    const float horizontal = sqrtf(fmaxf(1.0f - axis.y * axis.y, 0.0f));
    const float high =
        axis.y >= expanded_cosine
            ? 1.0f
            : axis.y * expanded_cosine + horizontal * expanded_sine;
    const float low =
        axis.y <= -expanded_cosine
            ? -1.0f
            : axis.y * expanded_cosine - horizontal * expanded_sine;
    if (high < south || low > north)
      return true;
    if (columns == W)
      return false;
    const bool misses_first = math::dot(axis, first) < -expanded_sine;
    const bool misses_last = math::dot(axis, last) < -expanded_sine;
    return columns <= W / 2 ? misses_first || misses_last
                            : misses_first && misses_last;
  }
};

/// Cap radius of a face that no cap about its vertex mean bounds.
inline constexpr float UNBOUNDED_FACE_CAP = -1.0f;

/**
 * @brief Axis of the cap about a face's vertex mean.
 * @param vertices Mesh vertex pool.
 * @param indices The face's vertex indices.
 * @param count Number of face vertices.
 * @param axis Receives the unit vertex mean.
 * @return False when the vertex mean is too short to give a direction.
 */
HS_HOT_INLINE inline bool face_cap_axis(const math::Vector *vertices,
                                        const uint16_t *indices, size_t count,
                                        math::Vector &axis) {
  math::Vector sum(0.0f, 0.0f, 0.0f);
  for (size_t k = 0; k < count; ++k)
    sum = sum + vertices[indices[k]];
  const float length = sqrtf(math::dot(sum, sum));
  if (!(length > math::TOLERANCE))
    return false;
  axis = sum * (1.0f / length);
  return true;
}

namespace detail {
/** @brief Cosine bound of a face about its vertex mean.
 * @param vertices Mesh vertex pool.
 * @param indices Face vertex indices.
 * @param count Number of face vertices.
 * @param axis Receives the unit cap axis.
 * @return Minimum vertex cosine, or UNBOUNDED_FACE_CAP. */
HS_HOT_INLINE inline float face_cap_cosine(const math::Vector *vertices,
                                           const uint16_t *indices,
                                           size_t count, math::Vector &axis) {
  if (!face_cap_axis(vertices, indices, count, axis))
    return UNBOUNDED_FACE_CAP;
  float min_cos = 1.0f;
  for (size_t k = 0; k < count; ++k) {
    const math::Vector &v = vertices[indices[k]];
    const float PROJECTION = math::dot(axis, v);
    if (!(PROJECTION > 0.0f))
      return UNBOUNDED_FACE_CAP;
    min_cos = fminf(min_cos, PROJECTION / sqrtf(math::dot(v, v)));
  }
  return min_cos;
}
} // namespace detail

/**
 * @brief Builds a compact face's seeing-time bounds frame.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param vertices Endpoint vertex pool.
 * @param indices Face vertex indices.
 * @param time Column basis, or null for a fixed meridian.
 * @param pose Receives the bounds frame and residual sweep.
 * @param pads Shared compact-face azimuth pads.
 * @param wide_pads Shared wide-face azimuth pads.
 * @param clip_planes Source clip expanded by this face pass's rotation reach.
 * @param culled Receives whether the face misses the expanded clip.
 * @return True when the seeing columns occupy a compact interval on one arm.
 */
template <int W, int H>
HS_HOT_FLASH_MEMBER bool
face_bounds_pose(std::span<const math::Vector> vertices,
                 std::span<const uint16_t> indices, const ColumnBasis *time,
                 SDF::FaceBoundsPose &pose, const float *pads,
                 const float *wide_pads,
                 const FaceClipPlanes<W, H> &clip_planes, bool &culled) {
  culled = false;
#if HS_ENABLE_TEST_HOOKS
  if (face_audit_reference)
    return false;
#endif
  math::Vector axis;
  const float cap_cosine = detail::face_cap_cosine(
      vertices.data(), indices.data(), indices.size(), axis);
  if (!(cap_cosine > 0.0f))
    return false;
  culled = clip_planes.rejects(axis, fmaxf(cap_cosine - 2e-6f, 0.0f));
  if (culled || !time || !time->azimuth || !(time->sweep_arc > 0.0f))
    return false;
  const float radius = cosine_arc(fminf(cap_cosine, 1.0f));
  int lo, hi;
  if (!lit_span<W>(axis, radius, time->sweep_arc, lo, hi) ||
      hi - lo > WIDE_FACE_COLUMNS)
    return false;
  constexpr int HALF = W / 2;
  const int arm = seeing_arm_start<W>(lo);
  if (hi >= arm + HALF)
    return false;
  const int m = (lo + hi) / 2 - arm;
  const auto FRAME = seeing_frame<W, H>(*time, m);
  pose.x = FRAME.equator * FRAME.cosine - FRAME.tangent * FRAME.sine;
  pose.y = FRAME.up;
  pose.z = FRAME.equator * FRAME.sine + FRAME.tangent * FRAME.cosine;
  const bool wide = hi - lo > COMPACT_FACE_COLUMNS;
  pose.residual = compact_face_residual<W>(
      time->sweep_rate, wide ? WIDE_FACE_COLUMNS : COMPACT_FACE_COLUMNS);
  pose.pads = wide ? wide_pads : pads;
  return true;
}

} // namespace Scan
