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
#include "core/render/spacetime/crossing.h"

/**
 * @file mesh.h
 * @brief rasterize_face, the -O3 face-specialized scan loop, and Scan::Mesh over it.
 */

namespace Scan {

#if HS_ENABLE_TEST_HOOKS
/** @brief Candidate observer for mesh replay verification. */
inline void (*face_audit_candidate)(size_t, int, int, float, bool) = nullptr;
/** @brief Face currently visited by the candidate observer. */
inline size_t face_audit_index = 0;
/** @brief Selected faces and vertices in the replayed frame. */
inline size_t visibility_counts[4] = {};
/** @brief Deformation reach used by the last visibility pre-pass. */
inline float visibility_audit_reach = 0.0f;
/** @brief Disables structural culling for a replay reference. */
inline bool face_audit_reference = false;

#endif

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
 * @param pipeline Plotting pipeline receiving the final colors.
 * @param canvas Destination canvas.
 * @param shape Face to rasterize.
 * @param fragment_shader Shader invoked per covered pixel.
 * @param time Column-time state with a rotated basis; read only with
 *        COLUMN_TIME.
 * @param projections Per-column projection storage for COLUMN_TIME.
 * @details Under a rotated basis a pole-LOD block's slack grows by the
 * rotation between its first and last columns, and a block straddling an arm
 * boundary is walked per column. Only the columns an open TimeWindow
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
  float radial_cull = shape.max_dist;
  {
    const float shade_reach = (probe_flags & SDF::PreparedFace::PROBE_LINEAR)
                                  ? PIXEL_WIDTH
                                  : tanf(PIXEL_WIDTH);
    radial_cull = shape.radius + shade_reach;
    if (probe_flags & SDF::PreparedFace::PROBE_HAS_LUT)
      radial_cull += (shape.lut_q_safe + 1) * shape.lut_dequant;
  }
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

/**
 * @brief Rotation bound across a compact face's seeing columns.
 * @tparam W Canvas width.
 * @param rate Frame rotation bound per frame fraction, as
 *        Scan::ColumnBasis::sweep_rate gives it.
 * @param columns Maximum seeing-column span.
 * @return Radians of residual rotation.
 */
template <int W>
HS_HOT_INLINE float compact_face_residual(float rate, int columns = 32) {
  return rate * (columns / 2 + 1.0f) / (W / 2) + rate * rate / math::PI_F;
}

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
  math::Vector axis(0, 0, 0);
  for (const uint16_t i : indices)
    axis = axis + vertices[i];
  const float length_sq = math::dot(axis, axis);
  if (!(length_sq > 1e-12f))
    return false;
  axis = axis * (1.0f / sqrtf(length_sq));
  float min_cos = 1.0f;
  for (const uint16_t i : indices)
    min_cos = fminf(min_cos, math::dot(axis, vertices[i]));
  if (!(min_cos > 0.0f))
    return false;
  float cap_cosine = 1.0f;
  for (const uint16_t i : indices) {
    const auto &v = vertices[i];
    cap_cosine = fminf(cap_cosine, math::dot(axis, v) / sqrtf(math::dot(v, v)));
  }
  culled = clip_planes.rejects(axis, fmaxf(cap_cosine - 2e-6f, 0.0f));
  if (culled || !time || !time->azimuth || !(time->sweep_arc > 0.0f))
    return false;
  const float radius = cosine_arc(fminf(min_cos, 1.0f));
  int lo, hi;
  if (!lit_span<W>(axis, radius, time->sweep_arc, lo, hi) || hi - lo > 64)
    return false;
  constexpr int HALF = W / 2;
  const int arm =
      static_cast<int>(floorf(static_cast<float>(lo) / HALF)) * HALF;
  if (hi >= arm + HALF)
    return false;
  const int m = (lo + hi) / 2 - arm;
  const math::Vector &e = time->azimuth[m];
  const math::Vector &v = time->up[m];
  const math::Vector t = math::cross(e, v);
  const float ct = math::TrigLUT<W, H>::cos_theta_data()[m];
  const float st = math::TrigLUT<W, H>::sin_theta[m];
  pose.x = e * ct - t * st;
  pose.y = v;
  pose.z = e * st + t * ct;
  const bool wide = hi - lo > 32;
  pose.residual = compact_face_residual<W>(time->sweep_rate, wide ? 64 : 32);
  pose.pads = wide ? wide_pads : pads;
  return true;
}

/**
 * @brief Whether an indexed element belongs to an optional visibility bitset.
 * @param bits Visibility bits, or null to admit every element.
 * @param i Element index.
 * @return True when the element is admitted.
 */
HS_HOT_INLINE inline bool mesh_visible(const uint8_t *bits, size_t i) {
  return !bits || (bits[i / 8] & (1u << (i % 8)));
}

/// Seeing-column entry of a face that holds over the frame.
inline constexpr uint16_t HELD_FACE_COLUMN = 0xFFFF;

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

/**
 * @brief Canvas columns that can see a cap over one frame.
 * @tparam W Canvas width in pixels.
 * @param axis Cap axis.
 * @param radius Cap radius (radians), or UNBOUNDED_FACE_CAP.
 * @param reach Largest arc the cap or a column's basis moves over the frame
 *        (radians).
 * @param columns Receives the columns as detail::face_lit_columns() gives them.
 */
template <int W>
HS_HOT_INLINE void cap_lit_columns(const math::Vector &axis, float radius,
                                   float reach, uint16_t *columns) {
  int lo, hi;
  if (radius == UNBOUNDED_FACE_CAP ||
      !lit_span<W>(axis, radius, reach, lo, hi) || hi - lo >= W - 1)
    return;
  columns[0] = static_cast<uint16_t>(math::wrap(lo, W));
  columns[1] = static_cast<uint16_t>(hi - lo);
}

namespace detail {
/**
 * @brief Cap about a face's vertex mean that holds every vertex.
 * @param vertices Mesh vertex pool.
 * @param indices The face's vertex indices.
 * @param count Number of face vertices.
 * @param axis Receives the cap axis, when the vertex mean gives a direction.
 * @return The cap radius (radians), or UNBOUNDED_FACE_CAP for a cap past a
 *         hemisphere or a vertex mean too short to give a direction.
 */
HS_HOT_INLINE inline float face_cap(const math::Vector *vertices,
                                    const uint16_t *indices, size_t count,
                                    math::Vector &axis) {
  if (!face_cap_axis(vertices, indices, count, axis))
    return UNBOUNDED_FACE_CAP;
  float min_cos = 1.0f;
  for (size_t k = 0; k < count; ++k) {
    const math::Vector &v = vertices[indices[k]];
    min_cos = fminf(min_cos, math::dot(axis, v) / sqrtf(math::dot(v, v)));
  }
  // A cap past a hemisphere does not bound its polygon.
  return min_cos > 0.0f ? cosine_arc(min_cos) : UNBOUNDED_FACE_CAP;
}
} // namespace detail

/**
 * @brief Selects the model-space faces and vertices a board can see this frame.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam Camera Model-to-endpoint rotation callable.
 * @param clip Board render clip.
 * @param mesh Undeformed model mesh.
 * @param camera Model-to-endpoint rotation.
 * @param reach Deformation displacement and capture-motion bound in radians.
 * @param arena Frame arena holding the selection bits.
 * @param time Column basis receiving the selection bits.
 */
template <int W, int H, typename Camera>
HS_HOT_FLASH_MEMBER void
select_mesh_visibility(const ClipRegion &clip, const MeshState &mesh,
                       const Camera &camera, float reach, Arena &arena,
                       ColumnBasis &time) {
  HS_PROFILE(mesh_visibility);
#if HS_ENABLE_TEST_HOOKS
  if (face_audit_reference)
    return;
  visibility_counts[0] += mesh.num_faces();
  visibility_counts[2] += mesh.vertices.size();
  visibility_audit_reach = reach;
#endif
  const size_t faces = mesh.num_faces(), vertices = mesh.vertices.size();
  auto *face_bits = arena.allocate_n<uint8_t>((faces + 7) / 8);
  auto *vertex_bits = arena.allocate_n<uint8_t>((vertices + 7) / 8);
  std::fill_n(face_bits, (faces + 7) / 8, uint8_t(0));
  std::fill_n(vertex_bits, (vertices + 7) / 8, uint8_t(0));
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  const auto xc = clip.x_clip();
  const float pad =
      reach + time.sweep_arc + SDF::BOUNDS_MARGIN_WIDE + SDF::MOVED_AXIS_SLACK;
  const float pad_cos = cosf(pad), pad_sin = sinf(pad);
  const auto geometry = math::LatitudeGeometry(H);
  const float row_step = geometry.row_to_phi(1) - geometry.row_to_phi(0);
  const float north =
      cosf(fmaxf(geometry.row_to_phi(clip.render_y_start()) - row_step, 0.0f));
  const float south = cosf(fminf(
      geometry.row_to_phi(clip.render_y_end() - 1) + row_step, math::PI_F));
  const int first_column = math::wrap(xc.rs, W);
  const int last_column = math::wrap(xc.rs + xc.length(W), W);
  const math::Vector first_normal(
      -math::TrigLUT<W, H>::sin_theta[first_column], 0.0f,
      math::TrigLUT<W, H>::cos_theta_data()[first_column]);
  const math::Vector last_normal(
      math::TrigLUT<W, H>::sin_theta[last_column], 0.0f,
      -math::TrigLUT<W, H>::cos_theta_data()[last_column]);
  float max_length_sq = 0.0f;
  for (const auto &v : mesh.vertices)
    max_length_sq = fmaxf(max_length_sq, math::dot(v, v));
  const float inverse_max_length =
      max_length_sq > 0.0f ? 1.0f / sqrtf(max_length_sq) : 0.0f;
  for (size_t f = 0; f < faces; ++f) {
    const size_t count = mesh.get_face_counts_data()[f];
    const uint16_t *indices =
        mesh.get_faces_data() + mesh.get_face_offsets_data()[f];
    math::Vector axis;
    if (pad < 0.5f * math::PI_F &&
        face_cap_axis(mesh.vertices.data(), indices, count, axis)) {
      float cosine = 1.0f;
      for (size_t k = 0; k < count; ++k)
        cosine = fminf(cosine, math::dot(axis, mesh.vertices[indices[k]]) *
                                   inverse_max_length);
      if (cosine > 0.0f) {
        const float sine = sqrtf(fmaxf(1.0f - cosine * cosine, 0.0f));
        const float expanded_cosine = cosine * pad_cos - sine * pad_sin;
        const float expanded_sine = sine * pad_cos + cosine * pad_sin;
        if (expanded_cosine > 0.0f) {
          axis = camera(axis);
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
            continue;
          const bool misses_first =
              math::dot(axis, first_normal) < -expanded_sine;
          const bool misses_last =
              math::dot(axis, last_normal) < -expanded_sine;
          if (xc.active &&
              (xc.length(W) <= W / 2 ? misses_first || misses_last
                                     : misses_first && misses_last))
            continue;
        }
      }
    }
#if HS_ENABLE_TEST_HOOKS
    ++visibility_counts[1];
#endif
    face_bits[f / 8] |= 1u << (f % 8);
    for (size_t k = 0; k < count; ++k)
      vertex_bits[indices[k] / 8] |= 1u << (indices[k] % 8);
  }
#if HS_ENABLE_TEST_HOOKS
  for (size_t i = 0; i < vertices; ++i)
    visibility_counts[3] += mesh_visible(vertex_bits, i);
#endif
  time.mesh.visible_faces = face_bits;
  time.mesh.visible_vertices = vertex_bits;
}

namespace detail {
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
 * @param radius Receives the face's cap radius for face_cap_columns(), when
 *        not null.
 * @details Judged by detail::face_cap() widened by @p reach; a cap past a hemisphere
 * or holding a pole is seen by every column.
 */
template <int W>
HS_FLASH_MEMBER void face_lit_columns(const math::Vector *vertices,
                                      const uint16_t *indices, size_t count,
                                      float reach, uint16_t *columns,
                                      float *radius = nullptr) {
  columns[0] = 0;
  columns[1] = W - 1;
  math::Vector axis;
  const float cap = detail::face_cap(vertices, indices, count, axis);
  if (radius)
    *radius = cap;
  cap_lit_columns<W>(axis, cap, reach, columns);
}
} // namespace detail

/**
 * @brief detail::face_lit_columns() of a face whose cap radius it recorded.
 * @tparam W Canvas width in pixels.
 * @param vertices Mesh vertex pool, as when the radius was recorded.
 * @param indices The face's vertex indices.
 * @param count Number of face vertices.
 * @param radius The recorded cap radius.
 * @param reach Largest arc the face or a column's basis moves over the frame
 *        (radians).
 * @param columns Receives the columns as detail::face_lit_columns() gives them.
 */
template <int W>
HS_FLASH_MEMBER void
face_cap_columns(const math::Vector *vertices, const uint16_t *indices,
                 size_t count, float radius, float reach, uint16_t *columns) {
  columns[0] = 0;
  columns[1] = W - 1;
  math::Vector axis;
  if (face_cap_axis(vertices, indices, count, axis))
    cap_lit_columns<W>(axis, radius, reach, columns);
}

namespace detail {
/**
 * @brief Whether a face's seeing columns meet a capture window's columns
 *        inside the render clip.
 * @param face The face's seeing columns, as detail::face_lit_columns() gives them.
 * @param begin The window's first half-turn column.
 * @param end One past the window's last half-turn column.
 * @param clip The render clip's column runs.
 * @param width Canvas width in pixels.
 * @param arm_only The one arm, 0 or 1, to test, or -1 for both.
 * @return True when a column of the window, on a tested arm, sees the face
 *         and the clip admits it.
 */
HS_FLASH_MEMBER inline bool face_columns_meet(const uint16_t *face, int begin,
                                              int end, const ClipPieces &clip,
                                              int width, int arm_only = -1) {
  const int half = width / 2;
  const int start = face[0];
  const int stop = start + face[1] + 1;
  // The face's columns past the last one wrap to the first.
  for (const int shift : {0, width})
    for (const int arm : {0, half}) {
      if (arm_only >= 0 && arm != arm_only * half)
        continue;
      const int lo = std::max(start - shift, begin + arm);
      const int hi = std::min(stop - shift, end + arm);
      for (int i = 0; i < clip.count; ++i)
        if (std::max(lo, clip.lo[i]) < std::min(hi, clip.hi[i]))
          return true;
    }
  return false;
}
} // namespace detail

namespace detail {
/**
 * @brief detail::face_columns_meet() under a render clip.
 * @tparam W Canvas width in pixels.
 * @param face The face's seeing columns, as detail::face_lit_columns() gives them.
 * @param begin The window's first half-turn column.
 * @param end One past the window's last half-turn column.
 * @param xc Render clip.
 * @return True when a column of the window, on either arm, sees the face and
 *         the clip admits it.
 */
template <int W>
HS_FLASH_MEMBER bool face_columns_meet(const uint16_t *face, int begin, int end,
                                       ClipRegion::XClip xc) {
  return detail::face_columns_meet(face, begin, end, clip_pieces<W>(xc), W);
}
} // namespace detail

/// Largest held vertex motion over a frame, in pixels of arc.
inline constexpr float HELD_VERTEX_PIXELS = 0.5f;

/// Arc a vertex moves over a frame, radians, up to which a deformation scan
/// holds it at the current endpoint.
template <int W>
inline constexpr float HELD_VERTEX_ARC =
    HELD_VERTEX_PIXELS * math::TWO_PI_F / W;

/// Moving-vertex flag: the vertex's path is captured.
inline constexpr uint8_t VERTEX_MOVES = 1;
/// Moving-vertex flags: shift of the vertex's reach level (vertex_flags()).
inline constexpr int VERTEX_REACH_SHIFT = 2;
/// Arc, radians, the lowest reach level stands for.
inline constexpr float VERTEX_REACH_UNIT = 1.0f / 1024.0f;
/// The arc each reach level stands for, 2^(1/3) apart, the top one
/// unbounded.
inline constexpr auto VERTEX_REACH_LEVELS = [] {
  std::array<float, 32> levels{};
  // 2^(1/3) to float precision.
  constexpr double STEP = 1.2599210498948732;
  double arc = VERTEX_REACH_UNIT;
  for (float &level : levels) {
    level = static_cast<float>(arc);
    arc *= STEP;
  }
  // The top level stands for any reach past the others.
  levels.back() = FLT_MAX;
  return levels;
}();

/**
 * @brief Flags of a moving vertex.
 * @param reach Bound on its distance from its endpoint over the frame
 *        (radians).
 * @return VERTEX_MOVES with the lowest reach level at least @p reach.
 */
HS_HOT_INLINE inline uint8_t vertex_flags(float reach) {
  // Binary search of the ascending levels; the top one bounds any reach.
  int level = 0;
  for (int step = 16; step > 0; step >>= 1)
    if (VERTEX_REACH_LEVELS[level + step - 1] < reach)
      level += step;
  return static_cast<uint8_t>(VERTEX_MOVES | (level << VERTEX_REACH_SHIFT));
}

/**
 * @brief A vertex's reach.
 * @param flags Its flags, as vertex_flags() gives them, or 0.
 * @return Radians: the bound its level stands for, 0 for a held vertex.
 */
HS_HOT_INLINE inline float vertex_reach(uint8_t flags) {
  return flags & VERTEX_MOVES
             ? VERTEX_REACH_LEVELS[(flags >> VERTEX_REACH_SHIFT) & 31]
             : 0.0f;
}

/**
 * @brief Bound on the arc between two positions of a vertex.
 * @param a,b The positions, any nonzero length.
 * @return Radians, at least the angle between them.
 */
HS_HOT_INLINE inline float vertex_arc(const math::Vector &a,
                                      const math::Vector &b) {
  const float norms = math::dot(a, a) * math::dot(b, b);
  return cosine_arc(norms > 0.0f ? math::dot(a, b) / sqrtf(norms) : 1.0f);
}

/**
 * @brief Flags the vertices that move between two endpoints.
 * @tparam W Canvas width in pixels.
 * @param previous Vertices at the previous endpoint.
 * @param current The same vertices at the current endpoint.
 * @param count Number of vertices.
 * @param moving Receives VERTEX_MOVES for a vertex that moves more than
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
    moving[i] = math::dot(previous[i], current[i]) < held_cos * sqrtf(norms)
                    ? VERTEX_MOVES
                    : 0;
  }
}

/**
 * @brief Records the reach of the moving vertices of a draw whose vertices
 *        move along the arcs between two endpoints.
 * @tparam W Display width.
 * @param previous Vertices at the previous endpoint.
 * @param current The same vertices at the current endpoint.
 * @param count Number of vertices.
 * @param moving Nonzero for the vertices to flag, as detail::pick_moving_faces()
 *        leaves them; each receives vertex_flags() of its arc when it moves
 *        more than HELD_VERTEX_ARC, 0 otherwise.
 */
template <int W>
HS_FLASH_MEMBER void flag_vertex_reaches(const math::Vector *previous,
                                         const math::Vector *current,
                                         size_t count, uint8_t *moving) {
  const float held_cos = cosf(HELD_VERTEX_ARC<W>);
  for (size_t i = 0; i < count; ++i) {
    if (!moving[i])
      continue;
    const float norms =
        math::dot(previous[i], previous[i]) * math::dot(current[i], current[i]);
    const float c =
        norms > 0.0f ? math::dot(previous[i], current[i]) / sqrtf(norms) : 1.0f;
    moving[i] = c < held_cos ? vertex_flags(cosine_arc(c)) : 0;
  }
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
  return !Spacetime::flat_surface() && pixels > HELD_VERTEX_PIXELS &&
         pixels <= MAX_VERTEX_MOTION_PIXELS;
}

/**
 * @brief Arc a mesh scan pads a face's rows by past its vertices.
 * @tparam W Canvas width in pixels.
 * @param sweep Largest arc a column's basis moves over the frame (radians).
 * @return Scan::Mesh's bounds margin with a column of slack.
 */
template <int W> HS_HOT_INLINE float scan_row_pad(float sweep) {
  return std::max(SDF::BOUNDS_MARGIN, math::TWO_PI_F / W) + sweep +
         math::TWO_PI_F / W;
}

/** @brief Vertex heights past which a moving face reaches no clip row. */
struct ClipRows {
  float south = -2.0f; ///< Reached by no face whose vertices all stand below.
  float north = 2.0f;  ///< Reached by no face whose vertices all stand above.

  /**
   * @brief The heights of a render clip's rows.
   * @tparam W Canvas width in pixels.
   * @param clip Render clip.
   * @param reach Bound on how far a placed vertex strays from the current
   *        endpoint over the frame (radians).
   * @param sweep Largest arc a column's basis moves over the frame (radians).
   * @return The heights; each lies past the equator by at least the padded
   *         reach, so a face beyond it keeps every placement in its
   *         hemisphere, where its edges bulge away from the clip.
   */
  template <int W>
  HS_FLASH_MEMBER static ClipRows of(const ClipRegion &clip, float reach,
                                     float sweep) {
    const math::LatitudeGeometry geometry(clip.h);
    const float pad = reach + scan_row_pad<W>(sweep);
    ClipRows rows;
    if (clip.render_y_end() < clip.h) {
      const float phi =
          std::max(geometry.row_to_phi(static_cast<float>(clip.render_y_end())),
                   0.5f * math::PI_F) +
          pad;
      if (phi < math::PI_F)
        rows.south = cosf(phi);
    }
    if (clip.render_y_start() > 0) {
      const float phi = std::min(geometry.row_to_phi(static_cast<float>(
                                     clip.render_y_start() - 1)),
                                 0.5f * math::PI_F) -
                        pad;
      if (phi > 0.0f)
        rows.north = cosf(phi);
    }
    return rows;
  }

  /**
   * @brief Whether a face can reach a clip row.
   * @param vertices Mesh vertex pool, at the current endpoint.
   * @param indices The face's vertex indices.
   * @param count Number of face vertices.
   * @return False when every vertex stands below `south` or above `north`.
   */
  HS_HOT_INLINE bool meet(const math::Vector *vertices, const uint16_t *indices,
                          size_t count) const {
    float low = 2.0f, high = -2.0f;
    for (size_t k = 0; k < count; ++k) {
      const float y = vertices[indices[k]].y;
      low = fminf(low, y);
      high = fmaxf(high, y);
    }
    return high >= south && low <= north;
  }
};

namespace detail {
/**
 * @brief Recomputes the seeing columns of a deformation scan's moving faces.
 * @tparam W Canvas width in pixels.
 * @param canvas Canvas with the frame's column-time state installed.
 * @param mesh The scanned mesh, at the current endpoint.
 * @param reach Bound on how far a placed vertex strays from the current
 *        endpoint over the frame (radians).
 * @param columns Per face two entries, HELD_FACE_COLUMN first for a face
 *        that holds; the others receive the columns that can see the face,
 *        as detail::face_lit_columns() gives them, or HELD_FACE_COLUMN and 0 when no
 *        placement of the face reaches a row the render clip admits.
 * @param radii Per-face cap radii detail::pick_moving_faces() recorded from the same
 *        vertices, or null to measure the caps.
 */
template <int W>
HS_FLASH_MEMBER void
moving_face_columns(const Canvas &canvas, const MeshState &mesh, float reach,
                    uint16_t *columns, const float *radii = nullptr) {
  const float sweep =
      canvas.column_basis() ? canvas.column_basis()->sweep_arc : 0.0f;
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  const ClipRows rows = ClipRows::of<W>(canvas.clip(), reach, sweep);
  for (size_t f = 0; f < mesh.num_faces(); ++f) {
    uint16_t *face = columns + 2 * f;
    if (face[0] == HELD_FACE_COLUMN)
      continue;
    if (!rows.meet(mesh.vertices.data(), indices + offsets[f], counts[f])) {
      face[0] = HELD_FACE_COLUMN;
      face[1] = 0;
      continue;
    }
    math::Vector axis;
    float radius = UNBOUNDED_FACE_CAP;
    if (!radii)
      radius = detail::face_cap(mesh.vertices.data(), indices + offsets[f],
                                counts[f], axis);
    else if (face_cap_axis(mesh.vertices.data(), indices + offsets[f],
                           counts[f], axis))
      radius = radii[f];
    face[0] = 0;
    face[1] = W - 1;
    cap_lit_columns<W>(axis, radius, sweep + reach, face);
  }
}
} // namespace detail

namespace detail {
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
 *        picked face, as detail::face_lit_columns() gives them, or HELD_FACE_COLUMN
 *        and 0 for the others.
 * @param radii When not null, receives the cap radius of each face holding a
 *        moving vertex that can reach a clip row, for detail::moving_face_columns().
 * @return The number of picked faces.
 * @details A face is picked when it holds a moving vertex and a column and a
 * row the render clip admits can see it; the others draw once at the current
 * endpoint.
 */
template <int W>
HS_FLASH_MEMBER size_t pick_moving_faces(const Canvas &canvas,
                                         const MeshState &mesh,
                                         uint8_t *moving_vertices, float reach,
                                         uint16_t *columns,
                                         float *radii = nullptr) {
  const uint8_t *counts = mesh.get_face_counts_data();
  const uint16_t *indices = mesh.get_faces_data();
  const uint16_t *offsets = mesh.get_face_offsets_data();
  const size_t faces = mesh.num_faces();
  const float sweep =
      canvas.column_basis() ? canvas.column_basis()->sweep_arc : 0.0f;
  const ClipRows rows = ClipRows::of<W>(canvas.clip(), reach, sweep);
  for (size_t f = 0; f < faces; ++f) {
    uint16_t *face = columns + 2 * f;
    uint8_t moves = 0;
    for (size_t k = 0; k < counts[f]; ++k)
      moves |= moving_vertices[indices[offsets[f] + k]];
    face[0] = HELD_FACE_COLUMN;
    face[1] = 0;
    if (!moves ||
        (canvas.column_basis() &&
         !mesh_visible(canvas.column_basis()->mesh.visible_faces, f)) ||
        !rows.meet(mesh.vertices.data(), indices + offsets[f], counts[f]))
      continue;
    math::Vector axis;
    const float cap = detail::face_cap(mesh.vertices.data(),
                                       indices + offsets[f], counts[f], axis);
    if (radii)
      radii[f] = cap;
    face[0] = 0;
    face[1] = W - 1;
    cap_lit_columns<W>(axis, cap, sweep + reach, face);
  }
  std::fill_n(moving_vertices, mesh.vertices.size(), uint8_t(0));
  size_t picked = 0;
  for (size_t f = 0; f < faces; ++f) {
    uint16_t *face = columns + 2 * f;
    if (face[0] == HELD_FACE_COLUMN)
      continue;
    if (!detail::face_columns_meet<W>(face, 0, W / 2, canvas.render_x_clip())) {
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
} // namespace detail

/**
 * @brief Returns to holding the picked faces none of whose vertices moves.
 * @param mesh The scanned mesh.
 * @param moving Nonzero per vertex that moves.
 * @param columns The faces' seeing columns, as detail::pick_moving_faces() gives them;
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
 * @brief Unwraps a capture fraction around a window's centre.
 * @param u Periodic arm fraction.
 * @param centre Window centre in arm fractions.
 * @return Capture fraction within half a frame of @p centre; past 1 or before
 *         0 for a vertex the next or previous arm's columns see.
 */
HS_HOT_INLINE inline float unwrap_capture(float u, float centre) {
  return u + roundf(centre - u);
}

/**
 * @brief Projects a vertex into the frame's moving column basis.
 * @param p Vertex in endpoint camera space.
 * @param time Frame's column basis.
 * @param u Fraction selecting the local column basis.
 * @param width Display width.
 * @return Periodic arm fraction of the projected vertex, to math::fast_atan2's
 *         accuracy.
 */
HS_HOT_INLINE inline float capture_fraction(const math::Vector &p,
                                            const Scan::ColumnBasis &time,
                                            float u, int width) {
  const int x = hs::clamp(static_cast<int>(u * (width / 2)), 0, width / 2 - 1);
  const math::Vector tangent = math::cross(time.azimuth[x], time.up[x]);
  const float fraction =
      Spacetime::surface_time(x, width) +
      math::fast_atan2(math::dot(p, tangent), math::dot(p, time.azimuth[x])) /
          math::PI_F;
  return math::wrap_t(fraction);
}

/**
 * @brief A vertex's column in the frame's moving column basis, packed.
 * @param p Vertex in endpoint camera space.
 * @param time Frame's column basis.
 * @param u Fraction selecting the local column basis.
 * @param width Display width.
 * @return capture_fraction() in 65536ths of an arm, wrapping at 1.
 */
HS_HOT_INLINE inline uint16_t vertex_column(const math::Vector &p,
                                            const Scan::ColumnBasis &time,
                                            float u, int width) {
  return static_cast<uint16_t>(static_cast<uint32_t>(
      capture_fraction(p, time, u, width) * 65536.0f + 0.5f));
}

#if HS_ENABLE_TEST_HOOKS
/** @brief Host work counts for moving mesh capture. */
struct MeshCaptureCounters {
  size_t deformations = 0;   ///< Deformation scans that classified faces.
  size_t vertices = 0;       ///< Moving vertices admitted by the segment clip.
  size_t evaluations = 0;    ///< Vertex callback evaluations.
  size_t captures = 0;       ///< Vertex capture solves.
  size_t shared_faces = 0;   ///< Moving faces scanned once in the shared pass.
  size_t seam_faces = 0;     ///< Seam-face copies in the arms' shared passes.
  size_t cross_faces = 0;    ///< Faces drawn in copies of their own per arm.
  size_t polar_faces = 0;    ///< Faces reaching a pole.
  size_t swept_faces = 0;    ///< Polar faces swept column pair by column pair.
  size_t sweep_columns = 0;  ///< Column pairs those sweeps scanned.
  size_t sweep_rows = 0;     ///< Rows those column pairs scanned, per arm.
  size_t sweep_tests = 0;    ///< Column pairs the sweeps tested.
  size_t sweep_pixels = 0;   ///< Probes of those rows within a pixel of a face.
  size_t sweep_probes = 0;   ///< Probes those sweeps made.
  size_t face_copies = 0;    ///< Faces scanned alone.
  size_t split_vertices = 0; ///< Vertices a face copy draws twice.
  size_t scans = 0;          ///< Mesh scan invocations.
  size_t face_setups = 0;    ///< Faces those scans constructed.
};
/** @brief Instrumentation for the current mesh measurement. */
inline MeshCaptureCounters mesh_capture_counters;
#endif

/**
 * @brief Evaluates one moving vertex at the time its column sees it.
 * @tparam VertexFn Callable (float u) returning the vertex at that fraction.
 * @param start The vertex's column at the frame's start, as vertex_column()
 *        gives it at u = 0.
 * @param end Its column at the frame's end, as vertex_column() gives it at
 *        u = 1.
 * @param centre Arm fraction the capture's branch lies within half a frame of.
 * @param width Display width.
 * @param vertex Closed-form vertex path, continuing the frame's motion past
 *        either endpoint.
 * @return Captured vertex.
 * @details The column a vertex stands in runs linearly between its columns
 * at the frame's two ends; the capture time is that track's arm pass
 * (Spacetime::crossings()) on the branch whose column time lies within half
 * a frame of @p centre, and the path is evaluated once, there. A vertex the
 * neighbouring arm's columns see is captured at that column's time, past 1
 * or before 0.
 */
template <typename VertexFn>
HS_HOT_INLINE math::Vector capture_vertex(uint16_t start, uint16_t end,
                                          float centre, int width,
                                          VertexFn &&vertex) {
#if HS_ENABLE_TEST_HOOKS
  ++mesh_capture_counters.captures;
  ++mesh_capture_counters.evaluations;
#endif
  constexpr float SCALE = 1.0f / 65536.0f;
  const float first = unwrap_capture(start * SCALE, centre);
  const float travel = unwrap_capture(end * SCALE, first) - first;
  // The track over the frames on either side, so a pass past either end of
  // this one is found too; its first-arm branch is the column's own time.
  const float half = width / 2.0f;
  float u = first;
  Spacetime::crossings(
      width, (first - travel) * half, (first + 2.0f * travel) * half, -1.0f,
      2.0f, false, false, [&](const Spacetime::Crossing &pass) HS_HOT_INLINE {
        if (pass.arm == 0)
          u = pass.time;
      });
  return vertex(u);
}

/** @brief Contiguous arm-time interval of a face. */
struct FaceCaptureInterval {
  float lo; ///< First vertex fraction; negative for polar faces.
  float hi; ///< Last unwrapped vertex fraction; a polar face's azimuth span
            ///< in half turns, 2 for one winding about the pole.
};

/**
 * @brief Finds the complement of a face's largest periodic vertex gap.
 * @param angles Signed vertex azimuths in half turns, one per face vertex.
 * @param count Number of face vertices.
 * @param reach Its fastest vertex's reach (vertex_reach()); 0 when none
 *        moves.
 * @param nearest Its nearest vertex's distance from the polar axis.
 * @param sweep The column basis's sweep over the frame.
 * @return Contiguous interval, or a negative lower bound for a face that
 *         winds about a pole, spans a half turn or reaches the pole: its
 *         fastest vertex's reach, with the basis's sweep, reaches a quarter
 *         of a right angle about the polar axis from one of its vertices.
 *         A polar face's upper bound is its azimuth span in half turns.
 */
HS_HOT_FLASH_MEMBER FaceCaptureInterval face_capture_interval(
    const float *angles, size_t count, float reach, float nearest, float sweep);

namespace detail {

/** @brief Borrowed face-capture callbacks with flash-resident dispatch. */
struct SweptMeshCallbacks {
  void *vertex_context; ///< Borrowed vertex evaluator.
  void *scan_context;   ///< Borrowed mesh scanner.
  math::Vector (*evaluate)(void *, size_t, float); ///< Vertex dispatch.
  void (*draw)(void *);                            ///< Scan dispatch.

  /**
   * @brief Borrows the callbacks for one deformation scan.
   * @tparam VertexFn Vertex path evaluator.
   * @tparam ScanFn Mesh scanner.
   * @param vertex Vertex path evaluator.
   * @param scan Mesh scanner.
   */
  template <typename VertexFn, typename ScanFn>
  SweptMeshCallbacks(VertexFn &vertex, ScanFn &scan)
      : vertex_context(const_cast<void *>(static_cast<const void *>(&vertex))),
        scan_context(const_cast<void *>(static_cast<const void *>(&scan))),
        evaluate(evaluate_vertex<VertexFn>), draw(scan_faces<ScanFn>) {}

  /**
   * @brief Evaluates a borrowed vertex path.
   * @tparam VertexFn Vertex path evaluator.
   * @param context Borrowed evaluator.
   * @param i Vertex index.
   * @param u Capture fraction.
   * @return Evaluated vertex.
   */
  template <typename VertexFn>
  HS_HOT_FLASH_MEMBER static math::Vector evaluate_vertex(void *context,
                                                          size_t i, float u) {
    return (*static_cast<VertexFn *>(context))(i, u);
  }

  /**
   * @brief Invokes a borrowed mesh scanner.
   * @tparam ScanFn Mesh scanner.
   * @param context Borrowed scanner.
   */
  template <typename ScanFn>
  HS_HOT_FLASH_MEMBER static void scan_faces(void *context) {
    (*static_cast<ScanFn *>(context))();
  }
};

/** @brief One deformation scan's frame state, shared by its face passes. */
struct SweptMesh {
  Canvas &canvas;          ///< Canvas carrying the scan's column-time state.
  MeshState &mesh;         ///< Vertex pool, overwritten by the captures.
  const uint16_t *columns; ///< Per-face seeing columns.
  const uint8_t *moving;   ///< vertex_flags() per vertex; 0 held.
  /// Moving vertices' columns at the frame's start (vertex_column()).
  const uint16_t *start;
  uint8_t *selected; ///< The scan's face selection table.
  Arena &scratch;    ///< Arena for the per-vertex and per-face tables.
  float pad; ///< Arc past its vertices that a moving face's scan can reach.
  ClipPieces clip;              ///< Render clip's column runs.
  SweptMeshCallbacks callbacks; ///< Vertex evaluator and mesh scanner.
};

/**
 * @brief Scans a moving mesh with each face's vertices at their column times.
 * @param deformation The scan's frame state, with the vertex pool at the
 *        current endpoint.
 * @details deform_mesh() describes the passes.
 */
HS_HOT_FLASH_MEMBER void scan_swept_faces(const SweptMesh &deformation);

} // namespace detail

/**
 * @brief Scans a moving mesh with every vertex at the time its column is lit.
 * @tparam W Display width.
 * @tparam H Display height.
 * @tparam VertexFn Callable (size_t i, float u) evaluating an endpoint-space vertex.
 * @tparam ScanFn Callable () scanning the mesh into the canvas's open window.
 * @param canvas Canvas with column-time state when deformation is active.
 * @param mesh Endpoint mesh, overwritten by the captured vertices.
 * @param active Whether the mesh deforms during this frame.
 * @param columns Per-face seeing columns; HELD_FACE_COLUMN marks held faces.
 * @param moving vertex_flags() for vertices whose paths are captured; 0 for
 *        those that hold.
 * @param start Each moving vertex's column at the frame's start, as
 *        vertex_column() gives it at u = 0 from the vertex evaluated there.
 * @param scratch Arena for face selection and vertex placement tables.
 * @param reach Bound on vertex motion in radians.
 * @param vertex Closed-form vertex path at a frame fraction, continuing the
 *        frame's motion past either endpoint.
 * @param scan Scans the mesh with its current vertices and selected faces.
 * @details Each moving vertex is captured at the sweep time of the column
 * that sees it (capture_vertex()), on its own arm or across the arm boundary
 * nearer it, once each, and shared by the faces around it. Held faces and
 * the moving faces clear of both arm boundaries draw in one scan. A face across an arm
 * boundary draws once per arm, clipped to that arm's columns, each edge's
 * vertex times continuing across the boundary the edge crosses, so frame
 * N's copy before a boundary and frame N+1's copy after it are one polygon
 * there. A face every column pair's meridian can reach (one around a pole,
 * spanning half a turn, or reaching the pole by its vertices' motion) draws
 * in one scan per half of the arms, each of its vertices at its time nearest
 * the half's middle; a vertex whose motion can carry it a quarter of a right
 * angle about the polar axis holds at its endpoint. A polar face whose single
 * placement would show a visible step draws one column pair at a time
 * instead, each pair's face standing as at the pair's sweep time
 * (Mesh::draw_face_sweep()).
 */
template <int W, int H, typename VertexFn, typename ScanFn>
HS_HOT_FLASH_MEMBER void
deform_mesh(Canvas &canvas, MeshState &mesh, bool active,
            const uint16_t *columns, const uint8_t *moving,
            const uint16_t *start, Arena &scratch, float reach,
            VertexFn &&vertex, ScanFn &&scan) {
  const detail::SweptMeshCallbacks callbacks(vertex, scan);
  if (!active) {
    callbacks.draw(callbacks.scan_context);
    return;
  }
  HS_CHECK(canvas.column_basis(),
           "mesh deformation needs the frame's column-time state installed");
  ScratchScope guard(scratch);
  uint8_t *selected = scratch.allocate_n<uint8_t>(mesh.num_faces());
  Scan::ColumnBasis time = *canvas.column_basis();
  time.mesh.selected_faces = selected;
  const float sweep = time.azimuth ? time.sweep_arc : 0.0f;
  float *face_pads = scratch.allocate_n<float>(H);
  fill_face_pads<W, H>(face_pads, sweep);
  time.mesh.face_pads = face_pads;
  if (sweep > 0.0f) {
    float *compact_pads = scratch.allocate_n<float>(H);
    float *wide_pads = scratch.allocate_n<float>(H);
    fill_face_pads<W, H>(compact_pads,
                         compact_face_residual<W>(time.sweep_rate));
    fill_face_pads<W, H>(wide_pads,
                         compact_face_residual<W>(time.sweep_rate, 64));
    time.mesh.compact_face_pads = compact_pads;
    time.mesh.wide_face_pads = wide_pads;
  }
  ColumnBasisScope scope(canvas, time);
  detail::scan_swept_faces({canvas, mesh, columns, moving, start, selected,
                            scratch, scan_row_pad<W>(sweep) + reach,
                            clip_pieces<W>(canvas.render_x_clip()), callbacks});
}

/// Angle past a swept face's edges within which its column pairs shade: a
/// column of antialiasing with slack.
template <int W>
inline constexpr float SWEEP_ROW_REACH = 1.25f * math::TWO_PI_F / W;

/**
 * @brief A face a deformation scan sweeps, standing at one column pair's
 *        sweep time at a time.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Per frame it projects the face's vertices, as SDF::PreparedFace projects
 * them, about the vertex centroid at the frame's middle, at the frame's
 * start, middle and end, and fits each projected vertex, a convex face's
 * half-planes and the face's size as quadratics in time through those
 * three. Per column pair it evaluates the fits by Horner's rule and finds
 * where the pair's meridian, a line in the projection, crosses the edges. A pixel is inside
 * where the meridian, walked from the pole away from the face, has crossed
 * an odd number of edges on the way to it; its distance is SDF::PreparedFace's: the
 * worst half-plane's on a convex face, else the nearest projected edge's, in
 * radians for a face whose size is. Its rows are those within
 * SWEEP_ROW_REACH of the crossings and of the vertices near the meridian, on
 * either arm and past either pole, and from the pole on when the face holds
 * it.
 */
template <int W, int H> class SweptFace {
public:
  /**
   * @brief Fits a sweep's face over the frame.
   * @param sweep The face and its vertex paths.
   * @param arena Arena holding the face's fits and per-column state for its
   *        scope.
   */
  HS_HOT_FLASH_MEMBER SweptFace(const Scan::PoleFaceSweep &sweep, Arena &arena)
      : count(static_cast<int>(sweep.count)), north(sweep.north) {
    edges = arena.allocate_n<Edge>(sweep.count);
    planes = arena.allocate_n<Plane>(sweep.count);
    motion = arena.allocate_n<Motion>(sweep.count);
    line2 = arena.allocate_n<float>(sweep.count);
    span = arena.allocate_n<Span>(sweep.count);
    order = arena.allocate_n<uint8_t>(sweep.count);
    for (int k = 0; k < count; ++k)
      order[k] = static_cast<uint8_t>(k);
    crossings[0] = arena.allocate_n<float>(sweep.count);
    crossings[1] = arena.allocate_n<float>(sweep.count);
    math::Vector middle(0.0f, 0.0f, 0.0f);
    for (int k = 0; k < count; ++k)
      middle = middle + sweep.paths[k].middle;
    axis = unit(middle);
    basis_u = math::perpendicular_axis(axis);
    basis_w = unit(math::cross(axis, basis_u));
    // The fitted quantities at the frame's start, middle and end.
    float sizes[3];
    bool convex_at[3];
    for (int j = 0; j < 3; ++j) {
      const float t = 0.5f * static_cast<float>(j);
      math::Vector sum(0.0f, 0.0f, 0.0f);
      for (int k = 0; k < count; ++k) {
        const math::Vector v = sweep.paths[k].at(t);
        sum = sum + v;
        float px, py;
        project(v, px, py);
        edges[k].x = px;
        edges[k].y = py;
        motion[k].x.sample(j, px);
        motion[k].y.sample(j, py);
      }
      float area2 = 0.0f, radius2 = 0.0f;
      stand_edges(area2, radius2);
      convex_at[j] = half_planes(area2);
      for (int k = 0; k < count; ++k) {
        motion[k].nx.sample(j, planes[k].nx);
        motion[k].ny.sample(j, planes[k].ny);
        motion[k].off.sample(j, planes[k].off);
      }
      float cx, cy;
      project(sum, cx, cy);
      sizes[j] = fmaxf(sqrtf(nearest_squared(cx, cy)),
                       sqrtf(radius2) * SDF::MIN_SIZE_RADIUS_RATIO);
    }
    convex = convex_at[0] && convex_at[1] && convex_at[2];
    linear = sizes[1] < LINEAR_SIZE;
    for (int j = 0; j < 3; ++j)
      if (!linear)
        sizes[j] = math::fast_atan2(sizes[j], 1.0f);
    for (int k = 0; k < count; ++k) {
      Motion &m = motion[k];
      m.x.fit();
      m.y.fit();
      m.nx.fit();
      m.ny.fit();
      m.off.fit();
    }
    for (int j = 0; j < 3; ++j)
      size_fit.sample(j, sizes[j]);
    size_fit.fit();
  }

  /** @brief The size shading normalizes by at the column pair column()
   *  stood the face at, as SDF::PreparedFace::size.
   *  @return Plane units or, for a large face, radians. */
  float size() const { return face_size; }

  /**
   * @brief Stands the face at a column pair's sweep time.
   * @param time Frame's column basis.
   * @param x The pair's first column, in [0, W / 2).
   * @param arms Arms to bound: bit 0 the first, bit 1 the second.
   * @param rows Receives each bounded arm's rows, culled where it misses the
   *        face.
   * @return Whether a bounded arm has rows; false where the face collapses,
   *         as a build culls.
   */
  HS_HOT_FLASH_MEMBER bool column(const Scan::ColumnBasis &time, int x,
                                  int arms, SDF::Bounds rows[2]) {
    rows[0] = rows[1] = SDF::BOUNDS_CULLED;
    const float t = Spacetime::SamplingSurface<W>::time(x);
    for (int k = 0; k < count; ++k) {
      edges[k].x = motion[k].x.at(t);
      edges[k].y = motion[k].y.at(t);
    }
    float area2 = 0.0f, radius2 = 0.0f;
    for (int k = 0; k < count; ++k) {
      Edge &e = edges[k];
      const Edge &next = edges[k + 1 == count ? 0 : k + 1];
      e.ex = next.x - e.x;
      e.ey = next.y - e.y;
      const float length2 = e.ex * e.ex + e.ey * e.ey;
      e.inv_length2 = length2 > 0.0f ? 1.0f / length2 : 0.0f;
      area2 += e.x * next.y - next.x * e.y;
      radius2 = fmaxf(radius2, e.x * e.x + e.y * e.y);
      if (convex) {
        planes[k].nx = motion[k].nx.at(t);
        planes[k].ny = motion[k].ny.at(t);
        planes[k].off = motion[k].off.at(t);
      }
    }
    if (fabsf(area2) <= SDF::COLLAPSED_AREA_RATIO * radius2)
      return false;
    face_size = size_fit.at(t);
    azimuth = time.azimuth[x];
    up = time.up[x];
    // The meridian's line in the projection, by its unit normal and offset.
    const math::Vector meridian = math::cross(azimuth, up);
    float la = math::dot(basis_u, meridian), lb = math::dot(basis_w, meridian);
    const float inv_norm = 1.0f / sqrtf(la * la + lb * lb);
    la *= inv_norm;
    lb *= inv_norm;
    const float le = math::dot(axis, meridian) * inv_norm;
    constexpr float REACH = SWEEP_ROW_REACH<W>;
    float lo[2] = {2.0f, 2.0f}, hi[2] = {-2.0f, -2.0f};
    // A projected point's arm and polar cosine along the meridian.
    const auto include = [&](float px, float py) HS_HOT_INLINE {
      const math::Vector q = axis + basis_u * px + basis_w * py;
      const int arm = math::dot(q, azimuth) >= 0.0f ? 0 : 1;
      const float c = math::dot(q, up) / sqrtf(math::dot(q, q));
      lo[arm] = fminf(lo[arm], c);
      hi[arm] = fmaxf(hi[arm], c);
      return std::pair<int, float>{arm, c};
    };
    crossed[0] = crossed[1] = 0;
    // Each edge's extent along the line, by the projection on its direction.
    along_x = -lb;
    along_y = la;
    float from_line = la * edges[0].x + lb * edges[0].y + le;
    float from_along = along_x * edges[0].x + along_y * edges[0].y;
    for (int k = 0; k < count; ++k) {
      Edge &e = edges[k];
      const Edge &next = edges[k + 1 == count ? 0 : k + 1];
      const float to_line = la * next.x + lb * next.y + le;
      const float to_along = along_x * next.x + along_y * next.y;
      span[k] = {fminf(from_along, to_along), fmaxf(from_along, to_along)};
      from_along = to_along;
      if (fabsf(from_line) <= REACH)
        include(e.x, e.y);
      // No pixel of the pair, on the meridian's line, comes nearer the edge
      // than its nearer end, unless the edge crosses the line.
      if ((from_line > 0.0f) == (to_line > 0.0f)) {
        const float apart = fminf(fabsf(from_line), fabsf(to_line));
        line2[k] = apart * apart;
      } else {
        line2[k] = 0.0f;
        const float s = from_line / (from_line - to_line);
        const auto [arm, c] = include(e.x + e.ex * s, e.y + e.ey * s);
        crossings[arm][crossed[arm]++] = c;
      }
      from_line = to_line;
    }
    // The edges nearest the line first; the order holds from the last pair,
    // so the insertion sort runs over nearly sorted keys.
    for (int i = 1; i < count; ++i) {
      const uint8_t k = order[i];
      const float key = line2[k];
      int j = i;
      for (; j > 0 && line2[order[j - 1]] > key; --j)
        order[j] = order[j - 1];
      order[j] = k;
    }
    // The pole nearer the face lies inside it when the walk from the other
    // crosses the first arm's edges an odd number of times.
    if (crossed[0] & 1) {
      for (int arm = 0; arm < 2; ++arm) {
        if (north)
          hi[arm] = 1.0f;
        else
          lo[arm] = -1.0f;
      }
    }
    // Each arm's polar angles, widened by the reach and math::fast_acos's
    // error, and past either pole onto the other arm.
    const float pad = REACH + static_cast<float>(math::FAST_ACOS_MAX_ERROR);
    float first[2], last[2];
    for (int arm = 0; arm < 2; ++arm) {
      const bool empty = lo[arm] > hi[arm];
      first[arm] = empty ? math::PI_F + pad : math::fast_acos(hi[arm]) - pad;
      last[arm] = empty ? -pad : math::fast_acos(lo[arm]) + pad;
    }
    bool any = false;
    for (int arm = 0; arm < 2; ++arm) {
      const int other = 1 - arm;
      float from = first[arm], to = last[arm];
      if (first[other] < 0.0f) {
        from = fminf(from, 0.0f);
        to = fmaxf(to, -first[other]);
      }
      if (last[other] > math::PI_F) {
        from = fminf(from, 2.0f * math::PI_F - last[other]);
        to = fmaxf(to, math::PI_F);
      }
      if (!(arms & (1 << arm)) || from > to)
        continue;
      rows[arm] =
          SDF::phi_bounds_to_rows<H>(fmaxf(from, 0.0f), fminf(to, math::PI_F));
      any = true;
    }
    return any;
  }

  /**
   * @brief The face's signed distance at a pixel of the column pair column()
   *        stood it at.
   * @param p Unit pixel direction.
   * @param res Receives the distance, negative inside, in the distance and
   *        raw-distance registers, and size().
   */
  HS_HOT_INLINE void probe(const math::Vector &p,
                           SDF::DistanceResult &res) const {
    const float depth = math::dot(p, axis);
    if (depth < math::TOLERANCE) {
      res = SDF::DistanceResult(SDF::FAR_SENTINEL, 0.0f, SDF::FAR_SENTINEL,
                                0.0f, face_size);
      return;
    }
    const float inv = 1.0f / depth;
    const float px = math::dot(p, basis_u) * inv,
                py = math::dot(p, basis_w) * inv;
    float signed_d;
    if (convex) {
      signed_d = -FLT_MAX;
      for (int k = 0; k < count; ++k) {
        const Plane &h = planes[k];
        signed_d = fmaxf(signed_d, h.nx * px + h.ny * py + h.off);
      }
    } else {
      const float c = math::dot(p, up);
      int passed;
      if (math::dot(p, azimuth) >= 0.0f) {
        passed = 0;
        for (int i = 0; i < crossed[0]; ++i)
          passed += north ? crossings[0][i] < c : crossings[0][i] > c;
      } else {
        passed = crossed[0];
        for (int i = 0; i < crossed[1]; ++i)
          passed += north ? crossings[1][i] > c : crossings[1][i] < c;
      }
      // From the edge nearest the last pixel, then the edges by their
      // distance from the meridian's line until that distance passes the
      // nearest found; an edge whose extent along the line keeps it farther
      // is passed over.
      float best = edge_squared(edges[nearest], px, py);
      const float along = along_x * px + along_y * py;
      for (int i = 0; i < count; ++i) {
        const int k = order[i];
        const float across = line2[k];
        if (across >= best)
          break;
        const float gap =
            fmaxf(fmaxf(span[k].lo - along, along - span[k].hi), 0.0f);
        if (across + gap * gap >= best)
          continue;
        const float d = edge_squared(edges[k], px, py);
        if (d < best) {
          best = d;
          nearest = k;
        }
      }
      const float d = sqrtf(best);
      signed_d = passed & 1 ? -d : d;
    }
    if (!linear)
      signed_d = math::fast_atan2(signed_d, 1.0f);
    res = SDF::DistanceResult(signed_d, 0.0f, signed_d, 0.0f, face_size);
  }

private:
  /// Size, in plane units, below which a face reports plane distances, as
  /// SDF::PreparedFace does.
  static constexpr float LINEAR_SIZE = 0.2f;

  /** @brief A quadratic in time through values at t = 0, 1/2 and 1. */
  struct Quadratic {
    float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f; ///< Coefficients, or samples.

    /// Records the sample at t = j / 2.
    HS_FLASH_INLINE void sample(int j, float value) {
      (j == 0 ? c0 : j == 1 ? c1 : c2) = value;
    }

    /// Turns the three samples into the quadratic through them.
    HS_FLASH_INLINE void fit() {
      const float s0 = c0, s1 = c1, s2 = c2;
      c1 = -3.0f * s0 + 4.0f * s1 - s2;
      c2 = 2.0f * s0 - 4.0f * s1 + 2.0f * s2;
    }

    /// The value at a frame fraction.
    HS_HOT_INLINE float at(float t) const { return c0 + t * (c1 + t * c2); }
  };

  /** @brief One projected vertex and the edge from it to the next. */
  struct Edge {
    float x, y;        ///< The vertex.
    float ex, ey;      ///< To the next.
    float inv_length2; ///< The edge's inverse squared length; 0 collapsed.
  };

  /** @brief An edge's extent along the meridian's line. */
  struct Span {
    float lo, hi; ///< Least and greatest projection on the line's direction.
  };

  /** @brief A convex face's outward half-plane along an edge. */
  struct Plane {
    float nx, ny, off; ///< Unit normal and offset.
  };

  /** @brief A vertex's and its edge's half-plane's fits over the frame. */
  struct Motion {
    Quadratic x, y;        ///< The projected vertex.
    Quadratic nx, ny, off; ///< The half-plane, on a convex face.
  };

  /// A direction scaled to unit length, or zero left zero.
  HS_HOT_INLINE static math::Vector unit(const math::Vector &v) {
    const float length2 = math::dot(v, v);
    return length2 > 0.0f ? v * (1.0f / sqrtf(length2)) : v;
  }

  /// Projects a direction, any length, about the axis as SDF::PreparedFace does.
  HS_FLASH_INLINE void project(const math::Vector &v, float &px,
                               float &py) const {
    float depth = math::dot(v, axis);
    if (fabsf(depth) < math::TOLERANCE)
      depth = copysignf(math::TOLERANCE, depth);
    const float inv = 1.0f / depth;
    px = math::dot(v, basis_u) * inv;
    py = math::dot(v, basis_w) * inv;
  }

  /// Fills the edges from the projected vertices, with twice the polygon's
  /// signed area and its squared circumradius about the axis.
  HS_FLASH_INLINE void stand_edges(float &area2, float &radius2) {
    for (int k = 0; k < count; ++k) {
      Edge &e = edges[k];
      const Edge &next = edges[k + 1 == count ? 0 : k + 1];
      e.ex = next.x - e.x;
      e.ey = next.y - e.y;
      const float length2 = e.ex * e.ex + e.ey * e.ey;
      e.inv_length2 = length2 > 0.0f ? 1.0f / length2 : 0.0f;
      area2 += e.x * next.y - next.x * e.y;
      radius2 = fmaxf(radius2, e.x * e.x + e.y * e.y);
      planes[k] = {0.0f, 0.0f, 0.0f};
    }
  }

  /**
   * @brief Finds whether the standing polygon is convex, as SDF::PreparedFace finds
   *        it, and its half-planes if so.
   * @param area2 Twice its signed area.
   * @return Whether it is.
   */
  HS_FLASH_INLINE bool half_planes(float area2) {
    bool pos = false, neg = false;
    int flips = 0;
    float dir = 0.0f;
    const Edge *previous = &edges[count - 1];
    for (int k = 0; k < count; ++k) {
      const Edge &e = edges[k];
      const float turn = previous->ex * e.ey - previous->ey * e.ex;
      const float scale =
          (previous->ex * previous->ex + previous->ey * previous->ey) *
          (e.ex * e.ex + e.ey * e.ey);
      previous = &e;
      if (turn * turn > SDF::TURN_EPS_SQ * scale) {
        pos |= turn > 0.0f;
        neg |= turn < 0.0f;
      }
      if (e.ey != 0.0f) {
        flips += dir * e.ey < 0.0f;
        dir = e.ey;
      }
    }
    if ((pos && neg) || flips < 1 || flips > 2)
      return false;
    const float sign = area2 >= 0.0f ? 1.0f : -1.0f;
    float centre_side = -FLT_MAX;
    for (int k = 0; k < count; ++k) {
      const Edge &e = edges[k];
      Plane &h = planes[k];
      const float length2 = e.ex * e.ex + e.ey * e.ey;
      if (length2 < 1e-12f)
        return false;
      const float inv = sign / sqrtf(length2);
      h.nx = e.ey * inv;
      h.ny = -e.ex * inv;
      h.off = -(h.nx * e.x + h.ny * e.y);
      centre_side = fmaxf(centre_side, h.off);
    }
    return centre_side < 0.0f;
  }

  /// Squared distance in the projection plane to an edge.
  HS_HOT_INLINE static float edge_squared(const Edge &e, float px, float py) {
    const float wx = px - e.x, wy = py - e.y;
    const float s =
        hs::clamp((wx * e.ex + wy * e.ey) * e.inv_length2, 0.0f, 1.0f);
    const float bx = wx - e.ex * s, by = wy - e.ey * s;
    return bx * bx + by * by;
  }

  /// Squared distance in the projection plane to the nearest edge.
  HS_FLASH_INLINE float nearest_squared(float px, float py) const {
    float best = FLT_MAX;
    for (int k = 0; k < count; ++k)
      best = fminf(best, edge_squared(edges[k], px, py));
    return best;
  }

  int count;               ///< Number of vertices.
  bool north;              ///< Whether it lies around the north pole.
  bool convex = false;     ///< Whether it probes by its half-planes.
  bool linear = true;      ///< Whether distances stay in plane units.
  Edge *edges;             ///< The vertices and edges.
  Plane *planes;           ///< A convex face's half-planes.
  Motion *motion;          ///< The fits.
  float *line2;            ///< Squared bound on each edge's distance from
                           ///< the column pair's pixels.
  uint8_t *order;          ///< The edges by line2, ascending.
  Span *span;              ///< Each edge's extent along the meridian's line.
  float along_x = 0.0f;    ///< The meridian's line's unit direction, x.
  float along_y = 0.0f;    ///< The meridian's line's unit direction, y.
  float *crossings[2];     ///< Cosines of the meridian's crossings, per arm.
  int crossed[2] = {0, 0}; ///< Crossings per arm.
  math::Vector axis;       ///< Projection axis.
  math::Vector basis_u;    ///< Projection x axis.
  math::Vector basis_w;    ///< Projection y axis.
  math::Vector azimuth;    ///< The pair's azimuth.
  math::Vector up;         ///< The pair's up.
  Quadratic size_fit;      ///< size()'s fit.
  float face_size = 0.0f;  ///< size().
  mutable int nearest = 0; ///< The edge nearest the last pixel probed.
};

HS_O3_BEGIN
/**
 * @brief Rasterizes a polygonal mesh by drawing each face as an SDF::PreparedFace.
 * @details Without a face-shader setup callback the face index reaches the
 * shader in register v2; with one, it reaches the shader only through what
 * the callback hoists.
 * canvas.debug() does not tint a mesh. Under the canvas's column-time basis
 * each face's bounds are swept by its arc and every pixel is evaluated along
 * its column's rotated direction. An open TimeWindow limits the scan to its
 * columns; in a deformation scan, Scan::ColumnBasis::selected_faces selects
 * the faces each pass draws.
 */
struct Mesh {
  /**
   * @brief Traps unless every flat face index addresses a vertex in the pool.
   * @param faces Flat per-face vertex index array.
   * @param num_indices Entries in @p faces.
   * @param num_verts Vertices the mesh carries.
   * @details SDF::PreparedFace indexes the vertex pool through an unchecked std::span.
   */
  HS_NOINLINE_NOCLONE
  static void check_face_index_domain(const uint16_t *faces, size_t num_indices,
                                      size_t num_verts) {
    for (size_t k = 0; k < num_indices; ++k)
      HS_CHECK(static_cast<size_t>(faces[k]) < num_verts,
               "mesh face index exceeds the vertex pool");
  }

  /**
   * @brief Shades and rasterizes one constructed face.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shape The face.
   * @param i The face's index in its mesh.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param face_shader_setup As draw_impl() takes it.
   * @param time The scan's column-time state; read when @p sweep is positive.
   * @param sweep Arc the scan's column basis turns over the frame, or 0 for
   *        the fixed meridian.
   * @param projections Draw-scoped column projection storage.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_INLINE static void raster_face(
      PipelineT &pipeline, Canvas &canvas, const SDF::PreparedFace &shape,
      size_t i, FragmentShaderT &fragment_shader,
      FaceShaderSetupT &face_shader_setup, const Scan::ColumnBasis *time,
      float sweep, FaceProjectionColumn *projections) {
#if HS_ENABLE_TEST_HOOKS
    face_audit_index = i;
#endif
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
              pipeline, canvas, shape, raster_shader, *time, projections);
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
              pipeline, canvas, shape, fragment_shader, *time, projections);
        else
          rasterize_face<W, H, PipelineT, true>(pipeline, canvas, shape,
                                                fragment_shader);
      }
    }
  }

  /**
   * @brief Rasterizes one face a deformation scan draws alone, from its
   *        copy's own vertices.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas, its render-x window open on the copy's
   *        arm.
   * @param copy The face copy.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param arena Arena for the face's scratch.
   * @param face_shader_setup As draw_impl() takes it.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_FLASH_MEMBER static void
  draw_face_copy(PipelineT &pipeline, Canvas &canvas,
                 const Scan::FaceArmCopy &copy,
                 FragmentShaderT &fragment_shader, Arena &arena,
                 FaceShaderSetupT &face_shader_setup) {
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    ScratchScope scope(arena);
    SDF::FaceStorage &storage = *new_face_storage(arena);
    const Scan::ColumnBasis *time = column_basis(pipeline, canvas);
    const float sweep = time && time->azimuth ? time->sweep_arc : 0.0f;
    const float *pads = time ? time->mesh.face_pads : nullptr;
    if (!pads) {
      float *filled = arena.allocate_n<float>(H);
      fill_face_pads<W, H>(filled, sweep);
      pads = filled;
    }
    const auto &cr = source_clip<W, H>(pipeline, canvas);
    const Spacetime::TimeWindow &window = canvas.time_window();
    const float *compact_pads = time ? time->mesh.compact_face_pads : nullptr;
    const float *wide_pads = time ? time->mesh.wide_face_pads : nullptr;
    if (sweep > 0.0f && !compact_pads) {
      float *filled = arena.allocate_n<float>(H);
      fill_face_pads<W, H>(filled, compact_face_residual<W>(time->sweep_rate));
      compact_pads = filled;
    }
    if (sweep > 0.0f && !wide_pads) {
      float *filled = arena.allocate_n<float>(H);
      fill_face_pads<W, H>(filled,
                           compact_face_residual<W>(time->sweep_rate, 64));
      wide_pads = filled;
    }
    SDF::FaceBoundsPose pose;
    const FaceClipPlanes<W, H> clip_planes(
        cr, std::max(SDF::BOUNDS_MARGIN, math::TWO_PI_F / W) + sweep +
                SDF::MOVED_AXIS_SLACK);
    bool culled;
    const bool seeing = face_bounds_pose<W, H>(
        {copy.vertices, copy.count}, {copy.indices, copy.count}, time, pose,
        compact_pads, wide_pads, clip_planes, culled);
    if (culled)
      return;
    SDF::Face face = [&]() HS_HOT_INLINE {
      HS_PROFILE(scan_face_setup);
      return SDF::Face({copy.vertices, copy.count}, {copy.indices, copy.count},
                       storage, math::LatitudeGeometry(H), H, &cr, pads,
                       std::max(SDF::BOUNDS_MARGIN, math::TWO_PI_F / W) + sweep,
                       seeing ? &pose : nullptr);
    }();
    SDF::PreparedFace &shape = face.prepared();
#if HS_ENABLE_TEST_HOOKS
    ++mesh_capture_counters.face_setups;
#endif
    if (shape.y_min > shape.y_max ||
        (window.active && !face_meets_window(shape, window, cr)))
      return;
    auto *projections =
        sweep > 0.0f ? new_face_projection_columns<W>(arena) : nullptr;
    raster_face<W, H>(pipeline, canvas, shape, copy.face, fragment_shader,
                      face_shader_setup, time, sweep, projections);
  }

  /**
   * @brief Rasterizes a mesh, or the one face a deformation scan draws
   *        alone or sweeps.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param mesh The mesh.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
   * @param bake Optional congruence-class bake for this mesh.
   * @param face_shader_setup As draw_impl() takes it.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_INLINE static void
  dispatch(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
           FragmentShaderT &fragment_shader, Arena &scratch_arena,
           const MeshOps::MeshClassBake *bake,
           FaceShaderSetupT &face_shader_setup) {
    const Scan::ColumnBasis *own = canvas.column_basis();
    if (own && own->mesh.face_sweep)
      draw_face_sweep<W, H>(pipeline, canvas, *own->mesh.face_sweep,
                            fragment_shader, scratch_arena, face_shader_setup);
    else if (own && own->mesh.face_copy)
      draw_face_copy<W, H>(pipeline, canvas, *own->mesh.face_copy,
                           fragment_shader, scratch_arena, face_shader_setup);
    else
      draw_impl<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena,
                      bake, face_shader_setup);
  }

  /**
   * @brief The azimuths whose meridian halves can come near a swept face.
   * @details A face clear of a pole lies within the azimuths of its vertices
   * at the frame's start, middle and end, widened by the reach its paths
   * bulge past their knots; a meridian half within a reach of it lies
   * within the asin of the reach over the sine of the face's least polar
   * angle over the frame. A face winding around the pole, reaching within the
   * reach of it or of its equator, spanning half a turn, or holding an edge
   * or a vertex path across nearly half a turn of azimuth, which an arc clear
   * of the pole turns through monotonically, is near every meridian half.
   */
  struct SweepWedge {
    float cx = 0.0f;        ///< Unit azimuth of the range's middle, x.
    float cz = 0.0f;        ///< Unit azimuth of the range's middle, z.
    float cos_half = -2.0f; ///< Cosine of its half-width; below -1 for all.

    /**
     * @brief The wedge of a sweep's face.
     * @param sweep The face and its vertex paths.
     * @param reach Arc within which a meridian half counts as near.
     * @return The wedge.
     */
    HS_HOT_FLASH_MEMBER static SweepWedge of(const Scan::PoleFaceSweep &sweep,
                                             float reach) {
      SweepWedge wedge;
      const size_t n = sweep.count;
      const math::Vector pole(0.0f, sweep.north ? 1.0f : -1.0f, 0.0f);
      const Spacetime::Worldline *paths = sweep.paths;
      const float reference = atan2f(paths[0].start.z, paths[0].start.x);
      float lo = 0.0f, hi = 0.0f, closest = math::PI_F;
      const auto knot = [&](size_t j, size_t k) HS_HOT_INLINE {
        return j == 0   ? paths[k].start
               : j == 1 ? paths[k].middle
                        : paths[k].end;
      };
      for (size_t j = 0; j < 3; ++j) {
        float winding = 0.0f;
        for (size_t k = 0; k < n; ++k) {
          const math::Vector a = knot(j, k);
          const math::Vector b = knot(j, k + 1 == n ? 0 : k + 1);
          const float turn = std::remainder(atan2f(b.z, b.x) - atan2f(a.z, a.x),
                                            math::TWO_PI_F);
          if (fabsf(turn) >= 0.9f * math::PI_F)
            return wedge;
          winding += turn;
          const float offset =
              std::remainder(atan2f(a.z, a.x) - reference, math::TWO_PI_F);
          lo = fminf(lo, offset);
          hi = fmaxf(hi, offset);
          closest = fminf(closest, arc_distance(pole, a, b));
        }
        if (fabsf(winding) > math::PI_F)
          return wedge;
      }
      // Each vertex keeps within its path's deviation of its chord.
      float bulge = 0.0f;
      for (size_t k = 0; k < n; ++k)
        bulge = fmaxf(bulge, paths[k].deviation());
      reach += bulge;
      for (size_t k = 0; k < n; ++k) {
        const math::Vector &a = paths[k].start, &b = paths[k].end;
        if (math::dot(a, pole) < sinf(reach) ||
            math::dot(b, pole) < sinf(reach))
          return wedge;
        if (fabsf(std::remainder(atan2f(b.z, b.x) - atan2f(a.z, a.x),
                                 math::TWO_PI_F)) >= 0.9f * math::PI_F)
          return wedge;
        closest = fminf(closest, arc_distance(pole, a, b));
      }
      if (closest <= reach || hi - lo >= math::PI_F)
        return wedge;
      const float spread = asinf(fminf(1.0f, sinf(reach) / sinf(closest)));
      const float half = 0.5f * (hi - lo) + spread;
      if (half >= math::PI_F)
        return wedge;
      const float middle = reference + 0.5f * (hi + lo);
      wedge.cx = cosf(middle);
      wedge.cz = sinf(middle);
      wedge.cos_half = cosf(half);
      return wedge;
    }

    /**
     * @brief Angle from a direction to the nearest point of an arc.
     * @param p Unit direction.
     * @param a The arc's start.
     * @param b The arc's end.
     * @return The angle, radians.
     */
    HS_HOT_FLASH_MEMBER static float arc_distance(const math::Vector &p,
                                                  const math::Vector &a,
                                                  const math::Vector &b) {
      const math::Vector ua = a.normalized(), ub = b.normalized();
      const math::Vector normal = math::cross(ua, ub);
      const float length = sqrtf(math::dot(normal, normal));
      const float to_a = acosf(hs::clamp(math::dot(p, ua), -1.0f, 1.0f));
      const float to_b = acosf(hs::clamp(math::dot(p, ub), -1.0f, 1.0f));
      if (length > 1e-6f) {
        const math::Vector unit = normal * (1.0f / length);
        const float height = math::dot(p, unit);
        const math::Vector foot = p - unit * height;
        if (math::dot(math::cross(ua, foot), unit) >= 0.0f &&
            math::dot(math::cross(foot, ub), unit) >= 0.0f)
          return asinf(hs::clamp(fabsf(height), 0.0f, 1.0f));
      }
      return fminf(to_a, to_b);
    }

    /**
     * @brief The arms of a column pair to scan.
     * @param azimuth The pair's first-arm meridian direction.
     * @param first Whether the clip admits the first arm's column.
     * @param second Whether the clip admits the second arm's column.
     * @return Bit 0 for the first arm, bit 1 for the second.
     */
    HS_HOT_INLINE int arms(const math::Vector &azimuth, bool first,
                           bool second) const {
      if (cos_half < -1.0f)
        return (first ? 1 : 0) | (second ? 2 : 0);
      const float along = (azimuth.x * cx + azimuth.z * cz) /
                          sqrtf(azimuth.x * azimuth.x + azimuth.z * azimuth.z);
      return (first && along >= cos_half ? 1 : 0) |
             (second && -along >= cos_half ? 2 : 0);
    }
  };

  /**
   * @brief Rasterizes a face a deformation scan sweeps, one column pair at a
   *        time, each pair's face standing as at its sweep time.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param sweep The face and its vertex paths.
   * @param fragment_shader Shader, as draw_impl() takes it.
   * @param arena Arena for the face's fits and per-column state.
   * @param face_shader_setup As draw_impl() takes it, run per column pair
   *        with the face's size there.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_FLASH_MEMBER static void
  draw_face_sweep(PipelineT &pipeline, Canvas &canvas,
                  const Scan::PoleFaceSweep &sweep,
                  FragmentShaderT &fragment_shader, Arena &arena,
                  FaceShaderSetupT &face_shader_setup) {
    HS_PROFILE(scan_face_sweep);
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    ScratchScope scope(arena);
    SweptFace<W, H> face(sweep, arena);
#if HS_ENABLE_TEST_HOOKS && defined(HS_SCAN_METRICS)
    const size_t near = hs::g_scan_metrics.shade_candidates;
#endif
    if constexpr (std::is_same_v<FaceShaderSetupT, std::nullptr_t>) {
      auto wrapper = [&](const math::Vector &p, Fragment &f_in) {
        // v2 carries the face index (decoded by mesh_face_index()).
        f_in.v2 = static_cast<float>(sweep.face);
        fragment_shader(p, f_in);
      };
      FragmentShaderFn raster_shader = wrapper;
      scan_swept_face<W, H>(pipeline, canvas, sweep, face, raster_shader,
                            face_shader_setup);
    } else {
      scan_swept_face<W, H>(pipeline, canvas, sweep, face, fragment_shader,
                            face_shader_setup);
    }
#if HS_ENABLE_TEST_HOOKS && defined(HS_SCAN_METRICS)
    mesh_capture_counters.sweep_pixels +=
        hs::g_scan_metrics.shade_candidates - near;
#endif
  }

  /**
   * @brief Scans a swept face column pair by column pair.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam FragmentShaderT Fragment shader callable.
   * @tparam FaceShaderSetupT As draw_impl() takes it.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param sweep The face and its vertex paths.
   * @param face The face, fitted over the frame.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param face_shader_setup As draw_impl() takes it, run per column pair
   *        with the face's size there.
   * @details Each column pair the clip admits on an arm the face can reach
   * (SweepWedge) stands the face at its sweep time (SweptFace::column()),
   * which bounds the pair's rows; each admitted arm's rows probe the
   * standing face along the pair's rotated basis and shade as the face
   * rasterizer shades, a time-aware shader at the pair's sweep time.
   */
  template <int W, int H, typename PipelineT, typename FragmentShaderT,
            typename FaceShaderSetupT>
  HS_HOT_FLASH_MEMBER static void
  scan_swept_face(PipelineT &pipeline, Canvas &canvas,
                  const Scan::PoleFaceSweep &sweep, SweptFace<W, H> &face,
                  FragmentShaderT &fragment_shader,
                  FaceShaderSetupT &face_shader_setup) {
    const Scan::ColumnBasis &time = *column_basis(pipeline, canvas);
    const auto &cr = source_clip<W, H>(pipeline, canvas);
    const ClipPieces clip = clip_pieces<W>(cr.x_clip());
    const auto admits = [&](int column) HS_HOT_INLINE {
      for (int r = 0; r < clip.count; ++r)
        if (column >= clip.lo[r] && column < clip.hi[r])
          return true;
      return false;
    };
    const SweepWedge wedge =
        SweepWedge::of(sweep, SWEEP_ROW_REACH<W> + time.sweep_arc);
    constexpr float PIXEL_WIDTH = 2.0f * math::PI_F / W;
    constexpr bool TIMED_SHADER =
        Spacetime::reads_time<FragmentShaderT, const math::Vector &,
                              Fragment &>;
    SDF::DistanceResult res;
    Fragment frag;
    for (int x = 0; x < W / 2; ++x) {
      const int arms =
          wedge.arms(time.azimuth[x], admits(x), admits(x + W / 2));
      if (!arms)
        continue;
#if HS_ENABLE_TEST_HOOKS
      ++mesh_capture_counters.sweep_tests;
#endif
      SDF::Bounds rows[2];
      bool stands;
      {
        HS_PROFILE(scan_sweep_column);
        stands = face.column(time, x, arms, rows);
      }
      if (!stands)
        continue;
      if constexpr (!std::is_same_v<FaceShaderSetupT, std::nullptr_t>)
        face_shader_setup(sweep.face, face.size());
      // A time-aware shader shades at the pair's time while its state moves.
      [[maybe_unused]] bool per_column = TIMED_SHADER;
      if constexpr (TIMED_SHADER) {
        if constexpr (requires { fragment_shader.moves(); })
          per_column = fragment_shader.moves();
      }
      const math::Vector &e = time.azimuth[x];
      const math::Vector &v = time.up[x];
      [[maybe_unused]] bool counted = false;
      for (int arm = 0; arm < 2; ++arm) {
        int y_lo, y_hi;
        if (rows[arm].y_min > rows[arm].y_max ||
            !clamp_rows_to_clip(rows[arm], cr, y_lo, y_hi))
          continue;
#if HS_ENABLE_TEST_HOOKS
        mesh_capture_counters.sweep_columns += !counted;
        counted = true;
        mesh_capture_counters.sweep_rows +=
            static_cast<size_t>(y_hi - y_lo + 1);
#endif
        const int column = x + arm * (W / 2);
        const float side = arm == 0 ? 1.0f : -1.0f;
        for (int y = y_lo; y <= y_hi; ++y) {
          const float sp = side * math::TrigLUT<W, H>::sin_phi[y];
          const float cp = math::TrigLUT<W, H>::cos_phi[y];
          const math::Vector p(sp * e.x + cp * v.x, sp * e.y + cp * v.y,
                               sp * e.z + cp * v.z);
#if HS_ENABLE_TEST_HOOKS
          ++mesh_capture_counters.sweep_probes;
#endif
          face.probe(p, res);
          const float d = res.dist;
          if (d >= PIXEL_WIDTH)
            continue;
          HS_SCAN_METRIC(hs::g_scan_metrics.shade_candidates++);
          const float alpha = solid_coverage(d, PIXEL_WIDTH);
          if (alpha <= MIN_ALPHA)
            continue;
          frag.color = Color4(0, 0, 0, 0);
          frag.pos = p;
          frag.v0 = res.t;
          frag.v1 = res.raw_dist;
          frag.v2 = 0.0f;
          frag.v3 = res.aux;
          frag.size = res.size;
          frag.age = 0;
          bool shaded = false;
          if constexpr (TIMED_SHADER) {
            if (per_column) {
              fragment_shader(p, frag, static_cast<float>(x) / (W / 2));
              shaded = true;
            }
          }
          if (!shaded)
            fragment_shader(p, frag);
          if (frag.color.alpha <= MIN_ALPHA)
            continue;
          const float a = frag.color.alpha * alpha;
          if constexpr (requires {
                          pipeline.plot_in_bounds(
                              canvas, column, y, frag.color.color, frag.age, a);
                        }) {
            pipeline.plot_in_bounds(canvas, column, y, frag.color.color,
                                    frag.age, a);
          } else {
            pipeline.plot(canvas, column, y, frag.color.color, frag.age, a);
          }
        }
      }
    }
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
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
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
    auto *storage = new_face_storage(scratch_arena);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const Scan::ColumnBasis *time = column_basis(pipeline, canvas);
    const float sweep = time && time->azimuth ? time->sweep_arc : 0.0f;
    auto *projections =
        sweep > 0.0f ? new_face_projection_columns<W>(scratch_arena) : nullptr;
    const float *azimuth_pads = time ? time->mesh.face_pads : nullptr;
    if (!azimuth_pads) {
      float *pads = scratch_arena.allocate_n<float>(H);
      fill_face_pads<W, H>(pads, sweep);
      azimuth_pads = pads;
    }
    const float *compact_pads = time ? time->mesh.compact_face_pads : nullptr;
    const float *wide_pads = time ? time->mesh.wide_face_pads : nullptr;
    if (sweep > 0.0f && !compact_pads) {
      float *filled = scratch_arena.allocate_n<float>(H);
      fill_face_pads<W, H>(filled, compact_face_residual<W>(time->sweep_rate));
      compact_pads = filled;
    }
    if (sweep > 0.0f && !wide_pads) {
      float *filled = scratch_arena.allocate_n<float>(H);
      fill_face_pads<W, H>(filled,
                           compact_face_residual<W>(time->sweep_rate, 64));
      wide_pads = filled;
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
    const Spacetime::TimeWindow &window = canvas.time_window();
    const FaceClipPlanes<W, H> clip_planes(cr, bounds_margin +
                                                   SDF::MOVED_AXIS_SLACK);
    const uint8_t *selected_faces =
        canvas.column_basis() ? canvas.column_basis()->mesh.selected_faces
                              : nullptr;
    const uint8_t drawn_face = window.active ? Scan::detail::FacePass::ARM_COPY
                                             : Scan::detail::FacePass::SHARED;
    for (size_t i = 0; i < num_f; ++i) {
      if (time && !mesh_visible(time->mesh.visible_faces, i))
        continue;
      size_t count = fc[i];

      HS_CHECK(static_cast<size_t>(fo[i]) + count <= fi_size,
               "mesh face span exceeds face index array");
      if (selected_faces && selected_faces[i] != drawn_face)
        continue;

      std::span<const math::Vector> verts(mesh.vertices.data(),
                                          mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      SDF::FaceBoundsPose pose;
      if (SDF::PreparedFace::phi_band_rejects(
              verts, indices, math::LatitudeGeometry(H), H, bounds_margin, &cr))
        continue;
      bool culled;
      const bool seeing =
          face_bounds_pose<W, H>(verts, indices, time, pose, compact_pads,
                                 wide_pads, clip_planes, culled);
      if (culled)
        continue;
      // IIFE scopes HS_PROFILE to Face construction (prvalue elided in place).
      SDF::Face face = [&] {
        HS_PROFILE(scan_face_setup);
        return SDF::Face(verts, indices, *storage, math::LatitudeGeometry(H), H,
                         &cr, azimuth_pads, bounds_margin,
                         seeing ? &pose : nullptr, true);
      }();
      SDF::PreparedFace &shape = face.prepared();
#if HS_ENABLE_TEST_HOOKS
      ++mesh_capture_counters.face_setups;
#endif

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

      raster_face<W, H>(pipeline, canvas, shape, i, fragment_shader,
                        face_shader_setup, time, sweep, projections);
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
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
   * @param bake Optional congruence-class bake for this mesh.
   */
  template <int W, int H, typename PipelineT>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   const MeshOps::MeshClassBake *bake = nullptr) {
    std::nullptr_t face_shader_setup = nullptr;
    dispatch<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena, bake,
                   face_shader_setup);
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
   * @param scratch_arena Arena supplying per-face SDF::PreparedFace scratch storage.
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
    dispatch<W, H>(pipeline, canvas, mesh, fragment_shader, scratch_arena, bake,
                   erased_setup);
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
