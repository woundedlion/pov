/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <type_traits>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <concepts>
#include <span>
#include <limits>
#include <optional>
#include "math/geometry.h"
#include "render/shading.h"
#include "render/clip.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "memory.h"
#include "render/plot/cull.h"

/**
 * @file raster.h
 * @brief rasterize(): the adaptive sub-stepping walk that turns a fragment
 * polyline into plotted samples.
 */

namespace Plot {

/**
 * @brief Applies an optional vertex shader to every control point.
 * @tparam FragmentsT Fragment container type.
 * @param vertex_shader Vertex shader to run; no-op if null.
 * @param pts Fragment container mutated in place.
 */
template <typename FragmentsT>
inline void apply_vertex_shader(VertexShaderRef vertex_shader,
                                FragmentsT &pts) {
  if (vertex_shader) {
    for (auto &p : pts)
      vertex_shader(p);
  }
}

/** @brief Adaptive raster sampling density. */
enum class RasterSamplingPolicy { DEFAULT, BALANCED, SELECTABLE };

/** @brief Balanced-policy target spacing in screen pixels. */
inline constexpr float BALANCED_SCREEN_STEP_PX = 1.125f;

/** @brief Pole-floor multiple below which balanced sampling keeps exact cadence. */
inline constexpr float BALANCED_POLE_GUARD_SCALE = 2.0f;

/** @brief Minimum sin²φ for balanced step reuse. */
inline constexpr float BALANCED_REUSE_MIN_SIN2 = 0.12f;

/** @brief Step-reuse ceiling, as a fraction of base_step. */
inline constexpr float BALANCED_REUSE_MAX_STEP_SCALE = 0.9f;

/** @brief Minimum tangent dot between consecutive full samples to reuse a step. */
inline constexpr float BALANCED_REUSE_MIN_TANGENT_DOT = 0.995f;

/** @brief Maximum relative step change, against the new step, to reuse a step. */
inline constexpr float BALANCED_REUSE_STEP_TOLERANCE = 0.1f;

/**
 * @brief Alpha gain compensating a stretched sample spacing.
 * @param alpha Per-sample coverage at the default spacing.
 * @param step_ratio Stretched step over the default step.
 * @return Gained coverage, clamped to 1.
 * @details Linear-in-alpha fit to source-over accumulation; it over-boosts at
 * high alpha, so a near-opaque stroke plots fully opaque.
 */
static inline float balanced_sample_alpha(float alpha, float step_ratio) {
  const float gain = 1.0f + (step_ratio - 1.0f) * (0.88f - 0.20f * alpha);
  return fminf(1.0f, alpha * gain);
}

#if HS_ENABLE_TEST_HOOKS
inline uint32_t g_planar_full_samples = 0;
inline uint32_t g_planar_position_samples = 0;
#endif

/**
 * @brief Antipode cutoff for the planar projection's stable-azimuth region.
 * @details The planar (azimuthal-equidistant) projection is singular at the
 * basis antipode (R→π: azimuth undefined). A control point whose dot with the
 * basis center is below -COS_PLANAR_ANTIPODE uses a geodesic edge.
 */
inline constexpr float COS_PLANAR_ANTIPODE = 0.999f;

/**
 * @brief Adaptive sub-step slots rasterize caches for one segment.
 * @tparam W Rasterization width.
 * @return The slot count: the 2*W screen sweep, floored at 64.
 */
template <int W> inline constexpr size_t rasterize_step_budget() {
  constexpr size_t SWEEP = 2 * static_cast<size_t>(W);
  return SWEEP > 64 ? SWEEP : 64;
}

/**
 * @brief Upper bound on the scratch_arena_a bytes rasterize binds for its own
 * caches, on top of the caller's Fragments buffer, which stays live across the
 * call.
 * @tparam W Rasterization width.
 * @param planar_segments Segment count of a planar-basis draw that derives arc
 *        registers; 0 for a geodesic polyline, which binds no per-segment
 *        cache.
 * @param trail_points Gated-trail point count, covering gate arrays live
 *        across the call or a preceding gate peak; 0 otherwise.
 * @return The cache size in bytes.
 * @details Covers the adaptive sub-step cache plus, under a planar basis, the
 * per-segment arc and seam caches. Deferred-shader position buffers are not
 * included.
 */
template <int W>
inline constexpr size_t rasterize_scratch_a_bytes(size_t planar_segments = 0,
                                                  size_t trail_points = 0) {
  return rasterize_step_budget<W>() * sizeof(float) +
         planar_segments * (sizeof(float) + sizeof(uint8_t)) +
         (trail_points == 0
              ? 0
              : trail_points * (2 * sizeof(float) + sizeof(uint8_t)) +
                    alignof(float));
}

/**
 * @brief Compile-time rasterize() configuration, passed as one NTTP.
 * @details Every field defaults to the plain cached-replay geodesic polyline.
 */
struct RasterConfig {
  /** Emit adaptive samples immediately instead of replaying a step cache. */
  bool single_pass = false;
  /**
   * Compile out planar, closed-loop, seam and omit-end support for an open
   * geodesic polyline.
   */
  bool open_geodesic = false;
  /** Recompute v0/v1 from the rendered planar perimeter. */
  bool derive_planar_arc_registers = true;
  /** Interpolate source fragment registers at each adaptive sample. */
  bool interpolate_registers = true;
  /** Adaptive screen-space sample density. */
  RasterSamplingPolicy sampling_policy = RasterSamplingPolicy::DEFAULT;
};

/** @brief Open or closed polyline, with seam registers only on a closed loop. */
class RasterLoop {
public:
  constexpr RasterLoop() = default;
  static constexpr RasterLoop closed(const Fragment *seam = nullptr) {
    return RasterLoop(seam);
  }
  constexpr bool is_closed() const { return closed_loop; }
  constexpr const Fragment *seam() const { return closing_fragment; }

private:
  const Fragment *closing_fragment = nullptr;
  bool closed_loop = false;
  constexpr explicit RasterLoop(const Fragment *seam)
      : closing_fragment(seam), closed_loop(true) {}
};

/** @brief Geodesic edges with optional visibility flags, or a planar chart. */
class RasterProjection {
public:
  constexpr RasterProjection() = default;
  /** @brief Supplies one visibility byte per rasterized edge, or none.
   * @details Flags combine EDGE_VISIBLE, EDGE_CLASSIFIED and EDGE_ONE_DOT.
   */
  static constexpr RasterProjection
  geodesic(std::span<const uint8_t> flags = {}) {
    return RasterProjection(nullptr, flags);
  }
  /** @brief Selects azimuthal-equidistant interpolation in the supplied chart.
   * @details Flags contain one byte per edge; planar rasterization reads only
   * EDGE_VISIBLE to replace the clip cull when the clip is active.
   */
  static constexpr RasterProjection
  planar(const math::Basis &basis, std::span<const uint8_t> flags = {}) {
    return RasterProjection(&basis, flags);
  }
  static RasterProjection planar(math::Basis &&,
                                 std::span<const uint8_t> = {}) = delete;
  static RasterProjection planar(const math::Basis &&,
                                 std::span<const uint8_t> = {}) = delete;
  constexpr const math::Basis *basis() const { return planar_basis; }
  constexpr std::span<const uint8_t> flags() const { return edge_flags; }

private:
  const math::Basis *planar_basis = nullptr;
  std::span<const uint8_t> edge_flags;
  constexpr RasterProjection(const math::Basis *basis,
                             std::span<const uint8_t> flags)
      : planar_basis(basis), edge_flags(flags) {}
};

/**
 * @brief Paired screen rows and columns for the same polyline points.
 * @details Rows use y_to_screen_row; columns use vector_to_theta. Valid only
 * for pipelines without a world-space stage.
 */
class PointProjections {
public:
  constexpr PointProjections() = default;
  template <size_t N>
  constexpr PointProjections(const float (&rows)[N], const float (&cols)[N])
      : screen_rows(rows), screen_cols(cols), count(N) {}
  static PointProjections paired(std::span<const float> rows,
                                 std::span<const float> cols) {
    HS_CHECK(rows.size() == cols.size(),
             "hoisted point projection rows and columns differ in length");
    return PointProjections(rows, cols);
  }
  constexpr const float *rows() const { return screen_rows; }
  constexpr const float *cols() const { return screen_cols; }
  constexpr size_t size() const { return count; }

private:
  const float *screen_rows = nullptr;
  const float *screen_cols = nullptr;
  size_t count = 0;
  constexpr PointProjections(std::span<const float> rows,
                             std::span<const float> cols)
      : screen_rows(rows.empty() ? nullptr : rows.data()),
        screen_cols(cols.empty() ? nullptr : cols.data()), count(rows.size()) {}
};

/**
 * @brief Optional rasterize() behaviors beyond the plain open geodesic
 * polyline; every field defaults to that common case.
 */
struct RasterOptions {
  /** edge_flags bit: the edge intersects the clip region. */
  static constexpr uint8_t EDGE_VISIBLE = Plot::EDGE_VISIBLE;
  /** edge_flags bit: the edge spans at most one screen step. */
  static constexpr uint8_t EDGE_ONE_DOT = 1u << 1;
  /** edge_flags bit: EDGE_ONE_DOT carries a verdict; else it is unclassified. */
  static constexpr uint8_t EDGE_CLASSIFIED = 1u << 2;

  RasterLoop loop{};
  RasterProjection projection{};
  /** Skip the final endpoint of an open line so adjoining arcs tile once. */
  bool omit_end = false;
  /**
   * Arc-fraction window outside which samples are not shaded or plotted; the
   * step schedule stays the whole segment's, so sample positions are
   * clip-independent. Single-segment polylines only.
   * @details Samples outside the window never reach pipeline stages. Widening
   * it may change history-stage state even when terminal clipping discards the
   * extra samples from the current framebuffer.
   */
  float plot_t_start = 0.0f;
  /** Upper bound of the plot_t_start window. */
  float plot_t_end = 1.0f;
  /** Hoisted projections used by the single-dot shortcut. */
  PointProjections point_projections{};
  /** Enables balanced sampling for a SELECTABLE rasterizer. */
  bool balanced_sampling = false;
#if HS_ENABLE_TEST_ORACLES
  /** Rebuild a planar sampler after culling instead of reusing cull samples. */
  bool rebuild_planar_sampler = false;
#endif
};

HS_O3_BEGIN
/** @brief Arc fraction of a replayed step: exactly 1 at the terminal step,
 * else @p dist over @p total_dist clamped to 1. */
HS_HOT_INLINE inline float replay_t(bool terminal, float dist,
                                    float total_dist) {
  return terminal ? 1.0f : fminf(dist / total_dist, 1.0f);
}

/**
 * @brief Moves an adaptive walk one step along its segment.
 * @param current_dist Arc walked so far, advanced in place.
 * @param total_dist Segment arc.
 * @param step Desired step.
 * @param endpoint_gap Set to the final step when it reaches the end.
 * @details A step within @p step of the end lands on it; under two steps the
 * walk halves the remainder so it never ends on a sliver.
 */
HS_HOT_INLINE inline void advance_distance(float &current_dist,
                                           float total_dist, float step,
                                           float &endpoint_gap) {
  float remaining = total_dist - current_dist;
  if (remaining <= step) {
    endpoint_gap = remaining;
    current_dist = total_dist;
  } else if (remaining < 2.0f * step) {
    current_dist += remaining * 0.5f;
  } else {
    current_dist += step;
  }
}

/**
 * @brief Balanced-policy step for a default-density step.
 * @param default_step Default-density screen step.
 * @param base_step Equatorial step 2π/W.
 * @return @p default_step near the pole floor; else the step stretched to
 *         BALANCED_SCREEN_STEP_PX, capped at @p base_step.
 */
HS_HOT_INLINE inline float balanced_step(float default_step, float base_step) {
  const float POLE_GUARD =
      base_step * MIN_POLE_SCALE * BALANCED_POLE_GUARD_SCALE;
  return default_step <= POLE_GUARD
             ? default_step
             : fminf(base_step,
                     default_step * (BALANCED_SCREEN_STEP_PX / SCREEN_STEP_PX));
}

/**
 * @brief Balanced-policy alpha scale of a segment's final endpoint.
 * @param endpoint_gap Arc from the last plotted sample to the endpoint.
 * @param default_step Default-density step at the last full sample.
 * @param step Step the walk was taking.
 * @return The gap the endpoint stands in for, clamped to
 *         [default_step, step], over default_step.
 */
HS_HOT_INLINE inline float
endpoint_alpha_scale(float endpoint_gap, float default_step, float step) {
  return hs::clamp(endpoint_gap, default_step, step) / default_step;
}

/** @brief True when @p p lies at the antipode of @p basis, where the planar
 * projection is singular. */
HS_HOT_INLINE inline bool at_planar_antipode(const math::Vector &p,
                                             const math::Basis &basis) {
  return math::dot(p, basis.v) < -COS_PLANAR_ANTIPODE;
}

/** @brief True when @p bit is set in edge @p i's visibility byte. */
HS_HOT_INLINE inline bool has_edge_flag(const uint8_t *edge_flags, size_t i,
                                        uint8_t bit) {
  return (edge_flags[i] & bit) != 0;
}

#include "render/plot/raster_walk.h"

/**
 * @brief True when a balanced walk may take its next step from a position-only
 * sample, reusing the current step.
 * @param sample Full sample just taken.
 * @param step Default-density step at @p sample.
 * @param previous_step Default-density step at the previous full sample.
 * @param previous_tangent Tangent at the previous full sample.
 * @param base_step Equatorial step 2π/W.
 * @details Requires a sample away from the poles, a step clear of the pole
 * floor and below the reuse ceiling, and a tangent and step that have barely
 * changed since the previous full sample.
 */
HS_HOT_INLINE inline bool can_reuse_step(const SamplePT &sample, float step,
                                         float previous_step,
                                         const math::Vector &previous_tangent,
                                         float base_step) {
  const float sin2 = 1.0f - sample.pos.y * sample.pos.y;
  return sin2 > BALANCED_REUSE_MIN_SIN2 &&
         step > base_step * MIN_POLE_SCALE * BALANCED_POLE_GUARD_SCALE &&
         step < base_step * BALANCED_REUSE_MAX_STEP_SCALE &&
         math::dot(sample.tan, previous_tangent) >
             BALANCED_REUSE_MIN_TANGENT_DOT &&
         fabsf(step - previous_step) < step * BALANCED_REUSE_STEP_TOLERANCE;
}

/**
 * @brief Adaptively rasterize a fragment polyline onto the sphere.
 *
 * Walks consecutive fragment pairs, picks a geodesic or planar interpolation
 * strategy per segment, sub-steps each segment at ≈one-pixel SCREEN-space
 * density (screen_step, clamped near the poles), and plots through the pipeline.
 * Segments whose row/column reach lies outside the active clip region are
 * culled.
 *
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam Cfg Compile-time behavior selection; see RasterConfig.
 * @tparam PipelineT Pipeline type.
 * @tparam FragmentShaderT Fragment shader type for direct raster pipelines.
 * @param source_pipeline Render pipeline that plots fragments.
 * @param canvas Target canvas (supplies the active clip band).
 * @param points Fragment polyline to rasterize.
 * @param fragment_shader Per-fragment shader applied before plotting; must be
 *                        non-null. An empty FragmentShaderFn traps once per
 *                        polyline; a typed shader cannot be empty.
 * @param opts Optional loop/projection/culling behaviors.
 */
template <int W, int H, RasterConfig Cfg = {}, typename PipelineT = PipelineRef,
          typename FragmentShaderT = FragmentShaderFn>
static void rasterize(PipelineT &source_pipeline, Canvas &canvas,
                      const Fragments &points, FragmentShaderT fragment_shader,
                      RasterOptions opts = {}) {
  static_assert(Cfg.single_pass ||
                    Cfg.sampling_policy == RasterSamplingPolicy::DEFAULT,
                "non-default raster sampling requires single_pass");
  constexpr bool SINGLE_PASS = Cfg.single_pass;
  constexpr bool OPEN_GEODESIC = Cfg.open_geodesic;
  constexpr bool DERIVE_PLANAR_ARC_REGISTERS = Cfg.derive_planar_arc_registers;
  constexpr bool INTERPOLATE_REGISTERS = Cfg.interpolate_registers;
  constexpr RasterSamplingPolicy SAMPLING_POLICY = Cfg.sampling_policy;
  if constexpr (OPEN_GEODESIC)
    HS_CHECK(!opts.loop.is_closed() && opts.projection.basis() == nullptr &&
                 !opts.omit_end && opts.loop.seam() == nullptr,
             "open_geodesic rasterize takes no loop, planar or omit-end "
             "options");
  // A canvas that is not W x H plots through a pipeline whose wrap period and
  // framebuffer stride disagree.
  HS_CHECK(canvas.width() == W && canvas.height() == H,
           "canvas size differs from the plot's W/H");
  // A direct-raster sink writes through a cached framebuffer base; the canvas
  // double-buffers, so a stale base is the buffer the display is scanning out.
  if constexpr (requires { source_pipeline.prepared_for(canvas); })
    HS_CHECK(source_pipeline.prepared_for(canvas),
             "direct raster pipeline not prepared for this canvas");
  // Erasure collapses pipeline and shader; the erased call matches neither
  // clause so recursion ends.
  if constexpr (!pipeline_direct_raster_path<PipelineT>() &&
                (!std::same_as<std::decay_t<PipelineT>, PipelineRef> ||
                 !std::same_as<std::decay_t<FragmentShaderT>,
                               FragmentShaderFn>)) {
    PipelineRef erased(source_pipeline);
    FragmentShaderFn erased_shader(fragment_shader);
    rasterize<W, H, Cfg>(erased, canvas, points, erased_shader, opts);
    return;
  }
  const bool world_identity = [&] {
    if constexpr (requires { source_pipeline.world_transform_is_identity; })
      return source_pipeline.world_transform_is_identity;
    else
      return true;
  }();
  HS_PLOT_COUNT(rings);
  const bool close_loop = OPEN_GEODESIC ? false : opts.loop.is_closed();
  const math::Basis *planar_basis =
      OPEN_GEODESIC ? nullptr : opts.projection.basis();
  const bool omit_end = OPEN_GEODESIC ? false : opts.omit_end;
  const uint8_t *edge_flags = opts.projection.flags().data();
  const float *point_rows = opts.point_projections.rows();
  const float *point_cols = opts.point_projections.cols();
  const Fragment *loop_seam = OPEN_GEODESIC ? nullptr : opts.loop.seam();
  const bool balanced_sampling =
      SAMPLING_POLICY == RasterSamplingPolicy::BALANCED ||
      (SAMPLING_POLICY == RasterSamplingPolicy::SELECTABLE &&
       opts.balanced_sampling);
  auto &pipeline = source_pipeline;
  const ScreenStepAxes step_axes =
      world_identity ? ScreenStepAxes{} : screen_step_axes(pipeline);
  size_t len = points.size();
  for (const Fragment &point : points)
    HS_CHECK(std::isfinite(point.pos.x) && std::isfinite(point.pos.y) &&
                 std::isfinite(point.pos.z),
             "rasterize control points must be finite");
  // Trap an empty shader once per polyline, even when every edge culls.
  if constexpr (std::same_as<std::decay_t<FragmentShaderT>, FragmentShaderFn>)
    HS_CHECK(fragment_shader, "rasterize requires a non-null fragment_shader");
  // A degenerate path is not drawn; a dot needs the vertex duplicated.
  if (len < 2)
    return;
  HS_CHECK(point_rows == nullptr || opts.point_projections.size() == len,
           "hoisted point projections need one entry per polyline point");

  size_t count = close_loop ? len : len - 1;
  HS_CHECK(edge_flags == nullptr || opts.projection.flags().size() == count,
           "edge_flags length must match the rasterized edge count");
  const float plot_t_start = opts.plot_t_start;
  const float plot_t_end = opts.plot_t_end;
  const bool plot_window = plot_t_start > 0.0f || plot_t_end < 1.0f;
  // Replay can overshoot t=1 by an ULP.
  const float plot_t_hi =
      plot_t_end < 1.0f ? plot_t_end : std::numeric_limits<float>::infinity();
  const bool PLOT_START =
      !plot_window || (plot_t_start <= 0.0f && plot_t_hi >= 0.0f);
  const bool PLOT_END =
      !plot_window || (plot_t_start <= 1.0f && plot_t_hi >= 1.0f);
  HS_CHECK(!plot_window || count == 1,
           "a plot window requires a single-segment polyline");
  HS_PLOT_ADD(edges, count);
  // scratch_arena_a is a LIFO bump allocator; a raw pointer into it must not
  // outlive the scope that produced it.
  ScratchScope sc_guard(scratch_arena_a);
  ArenaVector<float> steps_cache;
  // Cache one segment's adaptive steps; the simulation capacity is a backstop.
  // Single-pass emission uses the same bound without binding cache storage.
  size_t max_cache = rasterize_step_budget<W>();
#if HS_ENABLE_TEST_ORACLES
  if (g_step_budget_override != 0 && g_step_budget_override < max_cache)
    max_cache = g_step_budget_override;
#endif
  if constexpr (!SINGLE_PASS)
    steps_cache.bind(scratch_arena_a, max_cache);

  const bool has_planar_basis = (planar_basis != nullptr);
  constexpr bool REUSE_PLANAR_CULL_SAMPLES =
      SINGLE_PASS && !DERIVE_PLANAR_ARC_REGISTERS && !INTERPOLATE_REGISTERS &&
      pipeline_hoistable_cull<PipelineT>();
  auto segment_next = [&](size_t i) -> const Fragment & {
    if (loop_seam != nullptr && i + 1 == len)
      return *loop_seam;
    return points[(i + 1) % len];
  };
  PlanarArcTable<DERIVE_PLANAR_ARC_REGISTERS> arcs;
  PlanarArcStamp<DERIVE_PLANAR_ARC_REGISTERS> arc_uv =
      arcs.bind_and_measure(points, segment_next, count, planar_basis);

  // Adaptively sub-step and plot one segment. `sample(t)` returns the sphere
  // point and tangent estimate at arc fraction t in [0,1] under the chosen
  // strategy, `sample.pos(t)` the point alone; `total_dist` is the segment's
  // on-sphere length estimate (radians). Endpoints are omitted on interior /
  // closed segments so a shared vertex isn't plotted twice.
  auto process_segment = [&](auto &&sample, const Fragment &curr,
                             const Fragment &next, float total_dist,
                             bool is_last_segment) {
    constexpr bool NEWTON_UNIT_SAMPLER =
        requires { std::remove_cvref_t<decltype(sample)>::NEWTON_UNIT; };
    // Direct plotting requires unit fragment positions.
    // Coincident endpoints emit at most one dot.
    if (total_dist < math::EPS_GEOMETRIC) {
      bool should_omit = close_loop || !is_last_segment || omit_end;
      if (!should_omit && PLOT_START) {
        Fragment f;
        seed_fragment<INTERPOLATE_REGISTERS>(f, curr);
        arc_uv.stamp(f, 0.0f);
        shade_and_plot(pipeline, canvas, fragment_shader, curr.pos, f);
      }
      return;
    }

    // Equatorial step 2π/W: screen_step cap and balanced-threshold reference.
    const float base_step = (2.0f * math::PI_F) / W;
    auto adaptive_step = [&](const SamplePT &value) {
#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
      if (!g_reference_screen_step)
#endif
        if (!world_identity && step_axes.usable())
          return screen_step_from_axes<W, H>(value, step_axes);
      if constexpr (std::same_as<std::decay_t<PipelineT>, PipelineRef>) {
        return pipeline_screen_step<W, H>(pipeline, value, world_identity,
                                          &step_axes);
      } else {
        if (!world_identity)
          return pipeline_screen_step<W, H>(pipeline, value, false, &step_axes);
#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
        if (g_reference_screen_step)
          return screen_step_reference<W, H>(value.pos, value.tan, base_step);
#endif
        return screen_step<W, H>(value.pos, value.tan, base_step);
      }
    };
    int planar_arc_interval = 0;
    auto adaptive_sample = [&](float t) -> SamplePT {
      HS_PLOT_STALL_START(adaptive_start);
      SamplePT result;
      if constexpr (SINGLE_PASS && !DERIVE_PLANAR_ARC_REGISTERS &&
                    !INTERPOLATE_REGISTERS && requires {
                      sample.one_pass_monotonic(t, planar_arc_interval);
                    }) {
        result = sample.one_pass_monotonic(t, planar_arc_interval);
#if HS_ENABLE_TEST_HOOKS
        ++g_planar_full_samples;
#endif
      } else if constexpr (SINGLE_PASS && requires { sample.one_pass(t); })
        result = sample.one_pass(t);
      else
        result = sample(t);
      HS_PLOT_RENDER_COUNT(adaptive_samples);
      HS_PLOT_STALL_STOP(adaptive_sim, adaptive_start);
      return result;
    };
    HS_PLOT_COUNT(sim_samples);
    SamplePT smp = adaptive_sample(0.0f);
    float first_step = adaptive_step(smp);

    // FAST PATH: the whole segment spans ≤ one screen step, so a single dot
    // covers it. Keyed on screen length: a base_step arc can still cross
    // several pixels on a steep/near-polar segment.
    if (total_dist <= first_step) {
      HS_PLOT_COUNT(one_dot);
      if (PLOT_START) {
        Fragment f;
        seed_fragment<INTERPOLATE_REGISTERS>(f, curr);
        arc_uv.stamp(f, 0.0f);
        shade_and_plot(pipeline, canvas, fragment_shader, curr.pos, f);
      }
      if (!close_loop && is_last_segment && !omit_end && PLOT_END) {
        Fragment fl;
        seed_fragment<INTERPOLATE_REGISTERS>(fl, next);
        arc_uv.stamp(fl, total_dist);
        shade_and_plot(pipeline, canvas, fragment_shader, next.pos, fl);
      }
      return;
    }

    // Size each sub-step so consecutive samples land ~SCREEN_STEP_PX apart in
    // screen space. `smp`/`first_step` seed the first iteration.
    if constexpr (SINGLE_PASS) {
      HS_PROFILE_DEEP(plot_seg_single_pass);
      float current_dist = 0.0f;
      float current_t = 0.0f;
      float desired_step = first_step;
      float default_desired_step = first_step;
      float previous_full_step = first_step;
      math::Vector previous_full_tangent = smp.tan;
      bool reuse_step = false;
      if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT) {
        if (balanced_sampling)
          desired_step = balanced_step(first_step, base_step);
      }
      size_t step_count = 0;
      float backstop_stretch = 1.0f;
      // Arc from the last plotted sample to the segment end; the terminal step
      // is `remaining`, not the stretched `desired_step`.
      [[maybe_unused]] float endpoint_gap = total_dist;
      // Writes the unit position of `smp` to `p`.
      auto unit_position = [&](math::Vector &p) __attribute__((always_inline)) {
        if constexpr (OPEN_GEODESIC || NEWTON_UNIT_SAMPLER) {
          HS_PLOT_COUNT(normalizations);
#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
          if (g_reference_screen_step) {
            p = smp.pos.normalized();
          } else
#endif
            p = newton_unit(smp.pos);
        } else if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT &&
                             requires { sample.one_pass(current_t); }) {
          // Balanced walks use the Newton step, default-density walks the
          // exact normalize; the two can plot different rows.
          HS_PLOT_COUNT(normalizations);
          if (balanced_sampling) {
            p = newton_unit(smp.pos);
          } else {
            p = smp.pos.normalized();
          }
        } else {
          HS_PLOT_COUNT(normalizations);
          p = smp.pos.normalized();
        }
      };
      while (current_dist < total_dist) {
        math::Vector p;
        unit_position(p);
        if (!plot_window ||
            (current_t >= plot_t_start && current_t <= plot_t_hi)) {
          Fragment f;
          lerp_fragment<INTERPOLATE_REGISTERS>(f, curr, next, current_t, p);
          arc_uv.stamp(f, current_dist);
          std::optional<float> alpha_scale;
          if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT) {
            if (balanced_sampling)
              alpha_scale = desired_step / default_desired_step;
          }
          shade_and_plot(pipeline, canvas, fragment_shader, p, f, alpha_scale);
        }

        if (++step_count >= max_cache) {
          // Stretch factor matches the two-pass replay's; the hard stop bounds
          // the extra steps and can leave the segment's tail unplotted.
          if (backstop_stretch == 1.0f) {
            HS_PLOT_COUNT(backstops);
            HS_SCAN_METRIC(hs::g_scan_metrics.plot_backstop_hits++);
            backstop_stretch = total_dist / current_dist;
          } else if (step_count >= 2 * max_cache) {
            endpoint_gap = total_dist - current_dist;
            break;
          }
        }
        HS_PLOT_MAX(steps_peak, step_count);
        advance_distance(current_dist, total_dist, desired_step, endpoint_gap);
        if (current_dist < total_dist) {
          current_t = current_dist / total_dist;
          HS_PLOT_COUNT(sim_samples);
          if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT &&
                        requires {
                          sample.position_monotonic(current_t,
                                                    planar_arc_interval);
                        }) {
            if (balanced_sampling && reuse_step) {
              HS_PLOT_STALL_START(position_start);
              smp.pos =
                  sample.position_monotonic(current_t, planar_arc_interval);
              HS_PLOT_RENDER_COUNT(adaptive_samples);
              HS_PLOT_STALL_STOP(adaptive_sim, position_start);
#if HS_ENABLE_TEST_HOOKS
              ++g_planar_position_samples;
#endif
              reuse_step = false;
            } else {
              smp = adaptive_sample(current_t);
              default_desired_step = adaptive_step(smp);
              if (balanced_sampling) {
                reuse_step = world_identity &&
                             can_reuse_step(smp, default_desired_step,
                                            previous_full_step,
                                            previous_full_tangent, base_step);
                previous_full_step = default_desired_step;
                previous_full_tangent = smp.tan;
              }
            }
          } else {
            smp = adaptive_sample(current_t);
            default_desired_step = adaptive_step(smp);
          }
          desired_step = default_desired_step;
          if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT) {
            if (balanced_sampling)
              desired_step = balanced_step(default_desired_step, base_step);
          }
          desired_step *= backstop_stretch;
        }
      }
      if (!close_loop && is_last_segment && !omit_end &&
          (!plot_window || plot_t_hi >= 1.0f)) {
        Fragment f;
        seed_fragment<INTERPOLATE_REGISTERS>(f, next);
        arc_uv.stamp(f, total_dist);
        std::optional<float> alpha_scale;
        if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT) {
          // Gain the endpoint by the arc it actually stands in for, floored
          // at the default step so it never dims below the DEFAULT policy.
          if (balanced_sampling)
            alpha_scale = endpoint_alpha_scale(
                endpoint_gap, default_desired_step, desired_step);
        }
        shade_and_plot(pipeline, canvas, fragment_shader, next.pos, f,
                       alpha_scale);
      }
      return;
    }

    steps_cache.clear();
    float sim_dist = 0.0f;

    {
      HS_PROFILE_DEEP(plot_seg_sim);
      while (sim_dist < total_dist) {
        float step = steps_cache.is_empty() ? first_step : adaptive_step(smp);

        // Backstop on cache overflow: the normalized replay stretches the
        // cached steps over the rest of the segment.
        if (steps_cache.size() >= steps_cache.capacity()) {
          HS_PLOT_COUNT(backstops);
          HS_SCAN_METRIC(hs::g_scan_metrics.plot_backstop_hits++);
          break;
        }
        steps_cache.push_back(step);
        HS_PLOT_MAX(steps_peak, steps_cache.size());
        sim_dist += step;

        if (sim_dist < total_dist) {
          HS_PLOT_COUNT(sim_samples);
          smp = adaptive_sample(sim_dist / total_dist);
        }
      }
    }

    // scale <= 1 normally (the final step overshoots); > 1 after a backstop
    // break. Either way the replay spans exactly total_dist.
    HS_CHECK(sim_dist > 0.0f,
             "rasterize: simulated segment length is not positive");
    float scale = total_dist / sim_dist;
    bool omit_last = close_loop || !is_last_segment || omit_end;

    // Normalize interpolated positions before vector_to_pixel's acos(v.y).
    HS_PROFILE_DEEP(plot_seg_draw);
    if (!plot_window || plot_t_start <= 0.0f) {
      HS_PLOT_STALL_START(replay_start);
      HS_PLOT_COUNT(replay_samples);
      HS_PLOT_COUNT(normalizations);
      math::Vector start_pos = newton_unit(sample.pos(0.0f));
      Fragment f;
      lerp_fragment<INTERPOLATE_REGISTERS>(f, curr, next, 0.0f, start_pos);
      arc_uv.stamp(f, 0.0f);
      HS_PLOT_STALL_STOP(normalized_replay, replay_start);
      shade_and_plot(pipeline, canvas, fragment_shader, start_pos, f);
    }

    size_t loop_limit = omit_last ? steps_cache.size() - 1 : steps_cache.size();
    float current_dist = 0.0f;

    for (size_t j = 0; j < loop_limit; j++) {
      float step = steps_cache[j] * scale;
      current_dist += step;

      const bool terminal = !omit_last && j == loop_limit - 1;
      float t = replay_t(terminal, current_dist, total_dist);

      if (plot_window && (t < plot_t_start || t > plot_t_hi))
        continue;

      // `t` follows the rendered arc length. Under a planar basis
      // `PlanarArcStamp::stamp` rewrites the lerped v0/v1 from the sampled arc,
      // so arc-keyed shaders track the drawn position.
      HS_PLOT_STALL_START(replay_start);
      HS_PLOT_COUNT(replay_samples);
      HS_PLOT_COUNT(normalizations);
      math::Vector p = terminal ? next.pos : newton_unit(sample.pos(t));
      Fragment f;
      lerp_fragment<INTERPOLATE_REGISTERS>(f, curr, next, t, p);
      arc_uv.stamp(f, current_dist);
      HS_PLOT_STALL_STOP(normalized_replay, replay_start);
      shade_and_plot(pipeline, canvas, fragment_shader, p, f);
    }
  };

  const auto &cr = canvas.clip();
  const bool clip_active = !cr.is_full();
  const auto xc = cr.x_clip();

  // Emits one shader-run dot for points[k]; the precomputed projection is
  // consumed only when no world stage would lift it back to a world vector.
  auto plot_dot = [&](const Fragment &src, size_t k) {
    if (plot_window && ((k == 0 && !PLOT_START) || (k == 1 && !PLOT_END)))
      return;
    Fragment f;
    seed_fragment<INTERPOLATE_REGISTERS>(f, src);
    shade_fragment(fragment_shader, src.pos, f);
    if constexpr (pipeline_hoistable_projection<PipelineT>()) {
      if (point_rows != nullptr && point_cols != nullptr) {
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, point_cols[k], point_rows[k], f.color.color,
                      f.age, f.color.alpha);
        return;
      }
    }
    HS_PLOT_COUNT(plotted_samples);
    pipeline.plot(canvas, src.pos, f.color.color, f.age, f.color.alpha);
  };

  for (size_t i = 0; i < count; i++) {
    const Fragment &curr = points[i];
    const Fragment &next = segment_next(i);
    bool is_last_segment = (i == count - 1);
    PlanarEdgeSpan planar_cull_span;
    math::Vector planar_cull_end;
    bool reuse_planar_cull_samples = false;

    // --- Interpolation Strategy Selection ---
    // Branch-cut guard: the planar projection is singular at the basis antipode,
    // so a segment with an endpoint there falls back to a geodesic edge.
    bool antipodal_seam = false;
    if (has_planar_basis)
      antipodal_seam = arcs.is_active()
                           ? arcs.seam(i)
                           : at_planar_antipode(curr.pos, *planar_basis) ||
                                 at_planar_antipode(next.pos, *planar_basis);
    const bool use_planar = planar_basis && !antipodal_seam;

    arcs.advance(i, arc_uv);

    // Segment culling — skip if the edge's rendered row/column reach
    // (arc bulge included) lies outside the clip band; precomputed bits replace
    // the evaluation when the producer already ran the same predicate.
    if (clip_active) {
      HS_PLOT_COUNT(cull_tests);
      HS_PROFILE_DEEP(plot_seg_cull);
      bool visible;
      if constexpr (REUSE_PLANAR_CULL_SAMPLES) {
        bool rebuild_planar_sampler = false;
#if HS_ENABLE_TEST_ORACLES
        rebuild_planar_sampler = opts.rebuild_planar_sampler;
#endif
        if (edge_flags != nullptr) {
          visible = has_edge_flag(edge_flags, i, RasterOptions::EDGE_VISIBLE);
        } else if (use_planar && xc.active && !rebuild_planar_sampler) {
          planar_cull_span =
              make_planar_edge_span(curr.pos, next.pos, *planar_basis);
          visible = planar_edge_visible_in_clip<W, H>(
              cr, xc, curr.pos, next.pos, *planar_basis, planar_cull_span,
              &planar_cull_end);
          reuse_planar_cull_samples = visible;
        } else {
          visible =
              edge_visible_in_clip<W, H>(pipeline, cr, xc, curr.pos, next.pos,
                                         use_planar ? planar_basis : nullptr);
        }
      } else {
        visible =
            edge_flags != nullptr
                ? has_edge_flag(edge_flags, i, RasterOptions::EDGE_VISIBLE)
                : edge_visible_in_clip<W, H>(
                      pipeline, cr, xc, curr.pos, next.pos,
                      use_planar ? planar_basis : nullptr);
      }
      if (!visible) {
        HS_PLOT_COUNT(culled);
        continue;
      }
    }

    // Single-dot shortcut: an edge proven to span <= one screen step renders
    // exactly as process_segment's fast path (`PlanarArcStamp::stamp` is a
    // no-op without a planar basis), so plot it without building the sampler.
    // A predicate false negative falls through and re-evaluates exactly.
    const bool one_dot =
        world_identity && !has_planar_basis &&
        (edge_flags != nullptr &&
                 has_edge_flag(edge_flags, i, RasterOptions::EDGE_CLASSIFIED)
             ? has_edge_flag(edge_flags, i, RasterOptions::EDGE_ONE_DOT)
             : edge_fits_one_dot<W, H>(curr.pos, next.pos));
    if (one_dot) {
      HS_PLOT_COUNT(one_dot);
      plot_dot(curr, i);
      if (!close_loop && is_last_segment && !omit_end)
        plot_dot(next, i + 1);
      continue;
    }

    if (use_planar) {
      HS_PLOT_COUNT(planar);
      if constexpr (REUSE_PLANAR_CULL_SAMPLES) {
        PlanarEdgeSampler sampler;
        if (reuse_planar_cull_samples) {
          sampler = make_planar_edge_sampler(planar_cull_span, planar_cull_end,
                                             *planar_basis);
        } else {
          sampler = make_planar_edge_sampler(curr.pos, next.pos, *planar_basis);
        }
        process_segment(sampler, curr, next, sampler.dist, is_last_segment);
      } else {
        rasterize_planar_strategy(curr, next, *planar_basis, is_last_segment,
                                  process_segment);
      }
    } else {
      HS_PLOT_COUNT(geodesic);
      rasterize_geodesic_strategy(curr, next, is_last_segment, process_segment);
    }
  }
}
HS_O3_END

} // namespace Plot
