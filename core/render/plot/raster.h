/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <cassert>
#include <utility>
#include <type_traits>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <array>
#include <concepts>
#include "math/geometry.h"
#include "render/shading.h"
#include "platform/constants.h"
#include "render/clip.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include "render/plot/cull.h"

/**
 * @file raster.h
 * @brief rasterize(): the adaptive sub-stepping walk that turns a fragment
 * polyline into plotted samples, plus the trail gate that precomputes its
 * per-edge flags.
 */

namespace Plot {

/** @brief Adaptive raster sampling density. */
enum class RasterSamplingPolicy { DEFAULT, BALANCED, SELECTABLE };

/** @brief Balanced-policy target spacing in screen pixels. */
inline constexpr float BALANCED_SCREEN_STEP_PX = 1.125f;

/** @brief Pole-floor multiple below which balanced sampling keeps exact cadence. */
inline constexpr float BALANCED_POLE_GUARD_SCALE = 2.0f;

/** @brief Minimum sin²φ for balanced step reuse; √0.12 ≈ 7 × MIN_SIN_PHI. */
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
 * @details Linear-in-alpha fit to source-over accumulation. It tracks the exact
 * gain to a few percent up to alpha ~0.4 and over-boosts above that; at the
 * balanced ratio the clamp saturates past alpha ~0.85, so a near-opaque stroke
 * plots fully opaque and loses its soft edge.
 */
static inline float balanced_sample_alpha(float alpha, float step_ratio) {
  const float gain = 1.0f + (step_ratio - 1.0f) * (0.88f - 0.20f * alpha);
  return std::min(1.0f, alpha * gain);
}

#if HS_ENABLE_TEST_HOOKS
inline uint32_t g_planar_full_samples = 0;
inline uint32_t g_planar_position_samples = 0;
#endif

/**
 * @brief Antipode cutoff for the planar projection's stable-azimuth region.
 * @details The planar (azimuthal-equidistant) projection is singular at the
 * basis antipode (R→π: azimuth undefined). A control point whose dot with the
 * basis center is below this (≈ within 2.6° of the antipode) projects to an
 * unstable azimuth, so its segment falls back to a geodesic edge. cos(π − 0.045).
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
 * @param trail_points Point count of a ParticleSystem::draw trail; 0 for every
 *        other caller.
 * @return The cache size in bytes.
 * @details Effects sizing a custom scratch-A split must budget this alongside
 * their own buffers. Covers the adaptive sub-step cache plus, under a planar
 * basis, the per-segment arc and seam caches.
 *
 * ParticleSystem::draw allocates its trail-gate arrays — the per-edge
 * visibility bits and the hoisted per-point rows and columns — from the same
 * arena BEFORE the call and keeps them live across it, so @p trail_points folds
 * those three in (with the alignment slack between the byte and float blocks).
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
 * @details Named fields, so a call site reads
 * `rasterize<W, H, RasterConfig{.single_pass = true}>` instead of ordering four
 * bare booleans and a policy enum. Every field defaults to the plain
 * cached-replay geodesic polyline.
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

/**
 * @brief Optional rasterize() behaviors beyond the plain open geodesic
 * polyline; every field defaults to that common case.
 * @details Taken BY VALUE, never by const reference: a reference escapes the
 * aggregate's address, so every call site materializes it and no field reaches
 * the callee as a constant. Owned by the callee, IPA-SRA splits it back into
 * scalar arguments — worth 1,376 B of ITCM on the device image.
 */
struct RasterOptions {
  /** edge_flags bit: the edge intersects the clip region. */
  static constexpr uint8_t EDGE_VISIBLE = 1u << 0;
  /** edge_flags bit: the edge spans at most one screen step. */
  static constexpr uint8_t EDGE_ONE_DOT = 1u << 1;
  /** edge_flags bit: EDGE_ONE_DOT carries a verdict; else it is unclassified. */
  static constexpr uint8_t EDGE_CLASSIFIED = 1u << 2;

  /** Also draw the last→first edge. */
  bool close_loop = false;
  /**
   * Non-null selects azimuthal-equidistant interpolation (straight in the
   * projection); null uses geodesic edges.
   */
  const Basis *planar_basis = nullptr;
  /**
   * Open lines only: skip the final endpoint plot (each vertex is otherwise
   * plotted once by its outgoing segment), so abutting arcs tile a longer
   * curve without double-plotting the shared vertex.
   */
  bool omit_end = false;
  /**
   * Optional precomputed Tier-3 edge flags, one byte per rasterized edge:
   * points.size() - 1, or points.size() under close_loop. Each byte is
   * EDGE_VISIBLE, optionally OR'd with EDGE_CLASSIFIED | EDGE_ONE_DOT.
   * Geodesic polylines only: a planar polyline's per-edge basis depends on
   * rasterize()'s seam pre-pass.
   */
  const uint8_t *edge_flags = nullptr;
  /** Entries in edge_flags; asserted against the rasterized edge count. */
  size_t edge_flags_len = 0;
  /**
   * Optional per-point screen rows, y_to_screen_row of each points[k].pos.
   * With point_cols, lets the single-dot shortcut skip the projection. Only
   * consumed when the pipeline has no world-space stage; both arrays or
   * neither.
   */
  const float *point_rows = nullptr;
  /**
   * Optional per-point screen columns, vector_to_theta of each
   * points[k].pos.
   */
  const float *point_cols = nullptr;
  /**
   * Optional last-to-first target fragment carrying seam registers for a
   * closed loop without an overlapping point.
   */
  const Fragment *loop_seam = nullptr;
  /** Enables balanced sampling for a SELECTABLE rasterizer. */
  bool balanced_sampling = false;
#if HS_ENABLE_TEST_ORACLES
  /** Rebuild a planar sampler after culling instead of reusing cull samples. */
  bool rebuild_planar_sampler = false;
#endif
};

/**
 * @brief Gates one geodesic trail's edges against the clip in one hoisted pass.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam PipelineT Pipeline type; must have no world cull stage
 *         (pipeline_hoistable_cull), so the predicate sees the raw points.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param trail Geodesic fragment polyline (>= 2 unit-position points).
 * @param bits Output, one byte per edge (trail.size() - 1): 0 = culled, else
 *        EDGE_VISIBLE; valid as rasterize()'s edge_flags input for an open
 *        polyline.
 * @return False when no edge is visible; bits are then all zero.
 * @details The hoisted per-point coordinates and the whole-trail culls come
 * from trail_gate_prologue, shared with ParticleSystem::draw's gate.
 */
HS_O3_BEGIN
template <int W, int H, typename PipelineT>
static bool gate_trail_edges(const PipelineT &, const ClipRegion &cr,
                             const ClipRegion::XClip &xc,
                             const Fragments &trail, uint8_t *bits) {
  static_assert(pipeline_hoistable_cull<PipelineT>(),
                "gate_trail_edges requires a pipeline with no world cull "
                "stage; route others through edge_visible_in_clip");
  constexpr int H_VIRT = H + hs::H_OFFSET;
  const size_t n = trail.size();
  HS_CHECK(n >= 2);
  const size_t edges = n - 1;

  ScratchScope span_guard(scratch_arena_a);
  const TrailGatePrologue pro = trail_gate_prologue<W, H>(cr, xc, trail);
  if (pro.rejected) {
    std::fill_n(bits, edges, uint8_t{0});
    return false;
  }
  const float *rows = pro.rows;
  const float *cols = pro.cols;

  bool any = false;
  for (size_t e = 0; e < edges; ++e) {
    const Vector &ea = trail[e].pos;
    const Vector &eb = trail[e + 1].pos;

    // Cheap row tier: the exact span's interior extremum lies within arc/2 of
    // an endpoint and phi is 1-Lipschitz in arc length (arc <= (pi/2)*chord),
    // so the endpoint rows widened by chord*(H_VIRT-1)/4 contain the exact
    // span; the AA pad matches the exact test's own high end. A miss here
    // therefore implies the exact test below also misses, keeping the bits
    // identical while skipping the edge's cross/normalize/acos.
    {
      const Vector d = eb - ea;
      const float margin =
          sqrtf(dot(d, d)) * (static_cast<float>(H_VIRT - 1) * 0.25f);
      if (!cr.could_intersect_y(std::min(rows[e], rows[e + 1]) - margin,
                                std::max(rows[e], rows[e + 1]) + margin +
                                    GEODESIC_ROW_AA_PAD)) {
        bits[e] = 0;
        continue;
      }
    }

    const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
    const bool v = exact_geodesic_edge_visible_hoisted<W, H>(cr, xc, rows, cols,
                                                             e, ea, eb, es);
    bits[e] = v ? RasterOptions::EDGE_VISIBLE : uint8_t{0};
    any = any || v;
  }
  return any;
}
HS_O3_END

/**
 * @brief Adaptively rasterize a fragment polyline onto the sphere.
 *
 * Walks consecutive fragment pairs, picks a geodesic or planar interpolation
 * strategy per segment, sub-steps each segment at ≈one-pixel SCREEN-space
 * density (screen_step, clamped near the poles), and plots through the pipeline.
 * Segments whose full screen-row span lies outside the active clip band are
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
 *                        non-null (the per-pixel call sites below do not guard
 *                        it, so an empty ref traps once per pixel instead).
 * @param opts Optional loop/projection/culling behaviors; taken by value (see
 *             RasterOptions).
 */
HS_O3_BEGIN
template <int W, int H, RasterConfig Cfg = {}, typename PipelineT = PipelineRef,
          typename FragmentShaderT = FragmentShaderFn>
static void rasterize(PipelineT &source_pipeline, Canvas &canvas,
                      const Fragments &points, FragmentShaderT fragment_shader,
                      RasterOptions opts = {}) {
  constexpr bool SINGLE_PASS = Cfg.single_pass;
  constexpr bool OPEN_GEODESIC = Cfg.open_geodesic;
  constexpr bool DERIVE_PLANAR_ARC_REGISTERS = Cfg.derive_planar_arc_registers;
  constexpr bool INTERPOLATE_REGISTERS = Cfg.interpolate_registers;
  constexpr RasterSamplingPolicy SAMPLING_POLICY = Cfg.sampling_policy;
  if constexpr (OPEN_GEODESIC)
    HS_CHECK(!opts.close_loop && opts.planar_basis == nullptr &&
                 !opts.omit_end && opts.loop_seam == nullptr,
             "open_geodesic rasterize takes no loop, planar or omit-end "
             "options");
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
  HS_PLOT_COUNT(rings);
  const bool close_loop = OPEN_GEODESIC ? false : opts.close_loop;
  const Basis *planar_basis = OPEN_GEODESIC ? nullptr : opts.planar_basis;
  const bool omit_end = OPEN_GEODESIC ? false : opts.omit_end;
  const uint8_t *edge_flags = opts.edge_flags;
  const float *point_rows = opts.point_rows;
  const float *point_cols = opts.point_cols;
  const Fragment *loop_seam = OPEN_GEODESIC ? nullptr : opts.loop_seam;
  const bool balanced_sampling =
      SAMPLING_POLICY == RasterSamplingPolicy::BALANCED ||
      (SAMPLING_POLICY == RasterSamplingPolicy::SELECTABLE &&
       opts.balanced_sampling);
  auto &pipeline = source_pipeline;
  size_t len = points.size();
  // A degenerate path is not drawn — callers wanting a dot duplicate the vertex,
  // as Line::sample does.
  if (len < 2)
    return;
  // Trap a null shader once per polyline so the per-pixel fragment_shader()
  // calls below can't invoke a null thunk.
  if constexpr (std::same_as<std::decay_t<FragmentShaderT>, FragmentShaderFn>)
    HS_CHECK(fragment_shader, "rasterize requires a non-null fragment_shader");
  HS_CHECK(edge_flags == nullptr || planar_basis == nullptr,
           "precomputed edge visibility is geodesic-only");
  HS_CHECK(loop_seam == nullptr || close_loop,
           "a raster seam fragment requires a closed loop");
  HS_CHECK((point_rows == nullptr) == (point_cols == nullptr),
           "hoisted point projections take both rows and columns");

  size_t count = close_loop ? len : len - 1;
  HS_CHECK(edge_flags == nullptr || opts.edge_flags_len == count,
           "edge_flags length must match the rasterized edge count");
  HS_PLOT_ADD(edges, count);
  // SCRATCH ARENA CONTRACT (load-bearing): scratch_arena_a is a LIFO bump
  // allocator shared with Pixel::Feedback::flush; do not let a raw pointer into
  // it outlive the scope that produced it.
  ScratchScope sc_guard(scratch_arena_a);
  ArenaVector<float> steps_cache;
  // The cache holds ONE segment's adaptive sub-steps (cleared per segment).
  // Away from a pole, ≈SCREEN_STEP_PX spacing gives the usual screen-sweep
  // count. Near a pole, MIN_POLE_SCALE can instead create up to
  // 1/MIN_POLE_SCALE samples per base_step interval; those samples are not
  // covered by the planar W/H sweep derivation. At W=288 the measured
  // pole-crossing geodesic worst case is 546 of the 2·W=576 slots. A planar
  // chart line can bow farther, so the simulation loop retains its capacity
  // backstop. Single-pass emits as it goes and takes max_cache only as that
  // backstop, so it never binds the storage.
  size_t max_cache = rasterize_step_budget<W>();
#if HS_ENABLE_TEST_ORACLES
  if (g_step_budget_override != 0 && g_step_budget_override < max_cache)
    max_cache = g_step_budget_override;
#endif
  if constexpr (!SINGLE_PASS)
    steps_cache.bind(scratch_arena_a, max_cache);

  // PLANAR ARC REGISTERS (v0/v1): under a planar basis the rendered edge bows
  // longer than the geodesic chord, so re-derive v0/v1 from the true rendered
  // arc (`cumul`/`seg_base` track it, `total_arc` normalizes v0). Skipped for
  // geodesic polylines or when DERIVE_PLANAR_ARC_REGISTERS is false.
  const bool has_planar_basis = (planar_basis != nullptr);
  const bool override_uv = DERIVE_PLANAR_ARC_REGISTERS && has_planar_basis;
  constexpr bool REUSE_PLANAR_CULL_SAMPLES =
      SINGLE_PASS && !DERIVE_PLANAR_ARC_REGISTERS && !INTERPOLATE_REGISTERS &&
      pipeline_hoistable_cull<PipelineT>();
  auto segment_next = [&](size_t i) -> const Fragment & {
    if (loop_seam != nullptr && i + 1 == len)
      return *loop_seam;
    return points[(i + 1) % len];
  };
  float total_arc = 0.0f;
  // Per-segment rendered arc length and antipode-seam flag, reused by the draw
  // loop below so the seam decision is taken in exactly one place.
  ArenaVector<float> seg_arc_cache;
  ArenaVector<uint8_t> seg_seam_cache;
  if (override_uv) {
    seg_arc_cache.bind(scratch_arena_a, count);
    seg_seam_cache.bind(scratch_arena_a, count);
    const Vector &pcenter = planar_basis->v;
    for (size_t i = 0; i < count; i++) {
      const Vector &a = points[i].pos;
      const Vector &b = segment_next(i).pos;
      const bool seam = dot(a, pcenter) < -COS_PLANAR_ANTIPODE ||
                        dot(b, pcenter) < -COS_PLANAR_ANTIPODE;
      seg_seam_cache.push_back(seam ? 1 : 0);
      float seg =
          seam ? angle_between(a, b) : planar_arc_length(a, b, *planar_basis);
      seg_arc_cache.push_back(seg);
      total_arc += seg;
    }
  }
  float cumul = 0.0f;    // rendered arc reached so far (planar polylines only)
  float seg_base = 0.0f; // rendered arc at the in-flight segment's start

  auto shade_fragment = [&](const Vector &position, Fragment &fragment) {
    HS_MSP_STALL_START(shade_start);
    HS_MSP_COUNT(fragment_shader_calls);
    fragment_shader(position, fragment);
    HS_MSP_STALL_STOP(shade_palette, shade_start);
  };

  // Adaptively sub-step and plot one segment. `sample(t)` returns the sphere
  // point AND unit tangent at arc fraction t in [0,1] under the chosen strategy,
  // `sample.pos(t)` the point alone; `total_dist` is the segment's on-sphere
  // length (radians). Endpoints are omitted on interior / closed segments so a
  // shared vertex isn't plotted twice.
  auto process_segment = [&](auto &&sample, const Fragment &curr,
                             const Fragment &next, float total_dist,
                             bool is_last_segment) {
    constexpr bool NEWTON_UNIT_SAMPLER =
        requires { std::remove_cvref_t<decltype(sample)>::NEWTON_UNIT; };
    // Rewrite the arc registers from the rendered arc when a planar basis is in
    // force (see the pre-pass above): `d` is the arc drawn so far within this
    // segment, `seg_base` the arc at its start. No-op for geodesic polylines.
    auto set_arc_uv = [&](Fragment &f, float d) {
      if constexpr (!DERIVE_PLANAR_ARC_REGISTERS)
        return;
      if (!has_planar_basis)
        return;
      float arc = seg_base + d;
      f.v1 = arc;
      if (total_arc > math::EPS_GEOMETRIC)
        f.v0 = arc / total_arc;
    };
    // The degenerate and fast paths plot curr.pos/next.pos directly (original
    // sampled vertices), without the DRAWING PHASE renormalize that corrects
    // sample().pos's ~0.04% drift. Precondition: callers pass unit fragment
    // positions; the ~4e-6 an angle-addition vertex recurrence
    // (Star::sample_positions) leaves is two orders inside that drift and is
    // plotted as-is.
    // Degenerate (coincident endpoints): plot at most a single dot.
    if (total_dist < math::EPS_GEOMETRIC) {
      bool should_omit = close_loop || !is_last_segment || omit_end;
      if (!should_omit) {
        Fragment f_copy;
        if constexpr (INTERPOLATE_REGISTERS)
          f_copy = curr;
        f_copy.pos = curr.pos;
        f_copy.color = Color4(0, 0, 0, 0);
        set_arc_uv(f_copy, 0.0f);

        HS_PLOT_COUNT(shader_calls);
        shade_fragment(curr.pos, f_copy);
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, curr.pos, f_copy.color.color, f_copy.age,
                      f_copy.color.alpha);
      }
      return;
    }

    // Sub-step length at the segment start (also the first simulation step).
    const float base_step = (2.0f * PI_F) / W;
    auto balanced_step = [&](float default_step) {
      const float POLE_GUARD =
          base_step * MIN_POLE_SCALE * BALANCED_POLE_GUARD_SCALE;
      return default_step <= POLE_GUARD
                 ? default_step
                 : std::min(base_step, default_step * (BALANCED_SCREEN_STEP_PX /
                                                       SCREEN_STEP_PX));
    };
    auto adaptive_step = [&](const SamplePT &value) {
#if HS_ENABLE_TEST_ORACLES
      if (g_reference_screen_step)
        return screen_step_reference<W, H>(value.pos, value.tan, base_step);
#endif
      return screen_step<W, H>(value.pos, value.tan, base_step);
    };
    int planar_arc_interval = 0;
    auto adaptive_sample = [&](float t) -> SamplePT {
      HS_MSP_STALL_START(adaptive_start);
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
      HS_MSP_COUNT(adaptive_samples);
      HS_MSP_STALL_STOP(adaptive_sim, adaptive_start);
      return result;
    };
    HS_PLOT_COUNT(sim_samples);
    SamplePT smp = adaptive_sample(0.0f);
    float first_step = adaptive_step(smp);

    // FAST PATH: the whole segment spans ≤ one screen step, so a single dot
    // covers it. Keyed on SCREEN length, not arc length: a base_step arc can
    // still cross several pixels on a steep/near-polar segment, which an
    // arc-length test would undersample into a beaded line.
    if (total_dist <= first_step) {
      HS_PLOT_COUNT(one_dot);
      Fragment f;
      if constexpr (INTERPOLATE_REGISTERS)
        f = curr;
      f.pos = curr.pos;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, 0.0f);
      HS_PLOT_COUNT(shader_calls);
      shade_fragment(curr.pos, f);
      HS_PLOT_COUNT(plotted_samples);
      pipeline.plot(canvas, curr.pos, f.color.color, f.age, f.color.alpha);
      if (!close_loop && is_last_segment && !omit_end) {
        Fragment fl;
        if constexpr (INTERPOLATE_REGISTERS)
          fl = next;
        fl.pos = next.pos;
        fl.color = Color4(0, 0, 0, 0);
        set_arc_uv(fl, total_dist);
        HS_PLOT_COUNT(shader_calls);
        shade_fragment(next.pos, fl);
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, next.pos, fl.color.color, fl.age, fl.color.alpha);
      }
      return;
    }

    // Size each sub-step so consecutive samples land ~SCREEN_STEP_PX apart in
    // screen space. `smp`/`first_step` above seed the first iteration.
    if constexpr (SINGLE_PASS) {
      HS_PROFILE_DEEP(plot_seg_single_pass);
      float current_dist = 0.0f;
      float current_t = 0.0f;
      float desired_step = first_step;
      float default_desired_step = first_step;
      float previous_full_step = first_step;
      Vector previous_full_tangent = smp.tan;
      bool reuse_step = false;
      if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT) {
        if (balanced_sampling)
          desired_step = balanced_step(first_step);
      }
      size_t step_count = 0;
      float backstop_stretch = 1.0f;
      // Arc from the last plotted sample to the segment end; the terminal step
      // is `remaining`, not the stretched `desired_step`.
      [[maybe_unused]] float endpoint_gap = total_dist;
      while (current_dist < total_dist) {
        Vector p;
        if constexpr (OPEN_GEODESIC || NEWTON_UNIT_SAMPLER) {
          HS_PLOT_COUNT(normalizations);
#if HS_ENABLE_TEST_ORACLES
          if (g_reference_screen_step) {
            p = smp.pos.normalized();
          } else
#endif
            p = newton_unit(smp.pos);
        } else if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT &&
                             requires { sample.one_pass(current_t); }) {
          // Only the planar sampler reaches here. Its balanced walk corrects
          // with the Newton step, its default-density walk with the exact
          // normalize; the two land ~5e-6 of unit length apart and plot
          // different rows.
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
        Fragment f;
        if constexpr (INTERPOLATE_REGISTERS)
          f = Fragment::lerp_registers(curr, next, current_t);
        f.pos = p;
        f.color = Color4(0, 0, 0, 0);
        set_arc_uv(f, current_dist);
        HS_PLOT_COUNT(shader_calls);
        shade_fragment(p, f);
        if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT) {
          if (balanced_sampling) {
            const float alpha_scale = desired_step / default_desired_step;
            f.color.alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
          }
        }
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, p, f.color.color, f.age, f.color.alpha);

        if (++step_count >= max_cache) {
          // Stretch factor matches the two-pass replay's; the hard stop below
          // bounds the extra steps, and exits short of total_dist — the
          // segment's tail goes unplotted, where the two-pass replay always
          // stretches its cached steps over the whole segment.
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
        float remaining = total_dist - current_dist;
        if (remaining <= desired_step) {
          endpoint_gap = remaining;
          current_dist = total_dist;
        } else if (remaining < 2.0f * desired_step) {
          current_dist += remaining * 0.5f;
        } else {
          current_dist += desired_step;
        }
        if (current_dist < total_dist) {
          current_t = current_dist / total_dist;
          HS_PLOT_COUNT(sim_samples);
          if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT &&
                        requires {
                          sample.position_monotonic(current_t,
                                                    planar_arc_interval);
                        }) {
            if (balanced_sampling && reuse_step) {
              HS_MSP_STALL_START(position_start);
              smp.pos =
                  sample.position_monotonic(current_t, planar_arc_interval);
              HS_MSP_COUNT(adaptive_samples);
              HS_MSP_STALL_STOP(adaptive_sim, position_start);
#if HS_ENABLE_TEST_HOOKS
              ++g_planar_position_samples;
#endif
              reuse_step = false;
            } else {
              smp = adaptive_sample(current_t);
              default_desired_step = adaptive_step(smp);
              if (balanced_sampling) {
                const float sin2 = 1.0f - smp.pos.y * smp.pos.y;
                reuse_step =
                    sin2 > BALANCED_REUSE_MIN_SIN2 &&
                    default_desired_step > base_step * MIN_POLE_SCALE *
                                               BALANCED_POLE_GUARD_SCALE &&
                    default_desired_step <
                        base_step * BALANCED_REUSE_MAX_STEP_SCALE &&
                    dot(smp.tan, previous_full_tangent) >
                        BALANCED_REUSE_MIN_TANGENT_DOT &&
                    fabsf(default_desired_step - previous_full_step) <
                        default_desired_step * BALANCED_REUSE_STEP_TOLERANCE;
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
              desired_step = balanced_step(default_desired_step);
          }
          desired_step *= backstop_stretch;
        }
      }
      if (!close_loop && is_last_segment && !omit_end) {
        Fragment f;
        if constexpr (INTERPOLATE_REGISTERS)
          f = next;
        f.pos = next.pos;
        f.color = Color4(0, 0, 0, 0);
        set_arc_uv(f, total_dist);
        HS_PLOT_COUNT(shader_calls);
        shade_fragment(next.pos, f);
        if constexpr (SAMPLING_POLICY != RasterSamplingPolicy::DEFAULT) {
          if (balanced_sampling) {
            // Gain the endpoint by the arc it actually stands in for, floored
            // at the default step so it never dims below the DEFAULT policy.
            const float alpha_scale =
                hs::clamp(endpoint_gap, default_desired_step, desired_step) /
                default_desired_step;
            f.color.alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
          }
        }
        HS_PLOT_COUNT(plotted_samples);
        pipeline.plot(canvas, next.pos, f.color.color, f.age, f.color.alpha);
      }
      return;
    }

    steps_cache.clear();
    float sim_dist = 0.0f;

    {
      HS_PROFILE_DEEP(plot_seg_sim);
      while (sim_dist < total_dist) {
        float step = steps_cache.is_empty() ? first_step : adaptive_step(smp);

        // Backstop: a pathological segment could exceed the 2*W cache. Stop
        // subdividing and let the normalized replay stretch the cached steps
        // over the rest of the segment.
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

    // The final step normally overshoots total_dist (scale <= 1) and the
    // normalized replay stretches the cached steps back to exactly total_dist.
    // On the backstop break path sim_dist can fall short (scale > 1) and the
    // replay stretches over the remaining segment instead.
    HS_CHECK(sim_dist > 0.0f);
    float scale = total_dist / sim_dist;
    bool omit_last = close_loop || !is_last_segment || omit_end;

    // DRAWING PHASE
    //
    // sample().pos is non-unit by the fast sin/cos residual; vector_to_pixel's
    // phi = acos(v.y) offsets the row near the pole, so correct the interpolated
    // positions back to unit.
    HS_PROFILE_DEEP(plot_seg_draw);
    {
      HS_MSP_STALL_START(replay_start);
      HS_PLOT_COUNT(replay_samples);
      HS_PLOT_COUNT(normalizations);
      Vector start_pos = newton_unit(sample.pos(0.0f));
      Fragment f;
      if constexpr (INTERPOLATE_REGISTERS)
        f = Fragment::lerp_registers(curr, next, 0.0f);
      f.pos = start_pos;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, 0.0f);
      HS_MSP_STALL_STOP(normalized_replay, replay_start);

      HS_PLOT_COUNT(shader_calls);
      shade_fragment(start_pos, f);
      HS_PLOT_COUNT(plotted_samples);
      pipeline.plot(canvas, start_pos, f.color.color, f.age, f.color.alpha);
    }

    // The size() - 1 cannot underflow: HS_CHECK(sim_dist > 0) above implies the
    // simulation pushed at least one step.
    size_t loop_limit = omit_last ? steps_cache.size() - 1 : steps_cache.size();
    float current_dist = 0.0f;

    for (size_t j = 0; j < loop_limit; j++) {
      float step = steps_cache[j] * scale;
      current_dist += step;

      // total_dist > 0 here (HS_CHECK(sim_dist > 0) implies >=1 sim step).
      float t = current_dist / total_dist;

      // `t` (hence the drawn POSITION) is parameterized by the RENDERED arc
      // length. Registers are lerped from the control points; under a planar
      // basis set_arc_uv rewrites v0/v1 from the true rendered arc so a shader
      // keying off them as an arc-length proxy tracks the drawn position across
      // the planar bow. Geodesic edges keep the lerped registers.
      HS_MSP_STALL_START(replay_start);
      HS_PLOT_COUNT(replay_samples);
      HS_PLOT_COUNT(normalizations);
      Vector p = newton_unit(sample.pos(t));
      Fragment f;
      if constexpr (INTERPOLATE_REGISTERS)
        f = Fragment::lerp_registers(curr, next, t);
      f.pos = p;
      f.color = Color4(0, 0, 0, 0);
      set_arc_uv(f, current_dist);
      HS_MSP_STALL_STOP(normalized_replay, replay_start);

      HS_PLOT_COUNT(shader_calls);
      shade_fragment(p, f);
      HS_PLOT_COUNT(plotted_samples);
      pipeline.plot(canvas, p, f.color.color, f.age, f.color.alpha);
    }
  };

  const auto &cr = canvas.clip();
  const bool clip_active = !cr.is_full();
  const auto xc = cr.x_clip();

  // Emits one shader-run dot for points[k]; the precomputed projection is
  // consumed only when no world stage would lift it back to a world vector.
  auto plot_dot = [&](const Fragment &src, size_t k) {
    Fragment f;
    if constexpr (INTERPOLATE_REGISTERS)
      f = src;
    f.pos = src.pos;
    f.color = Color4(0, 0, 0, 0);
    HS_PLOT_COUNT(shader_calls);
    shade_fragment(src.pos, f);
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
    Vector planar_cull_end;
    bool reuse_planar_cull_samples = false;

    // --- Interpolation Strategy Selection ---
    // Branch-cut guard: the planar projection is singular at the basis antipode,
    // so a segment with an endpoint there falls back to a geodesic edge.
    bool antipodal_seam = false;
    if (has_planar_basis) {
      antipodal_seam =
          override_uv
              ? seg_seam_cache[i] != 0
              : dot(curr.pos, planar_basis->v) < -COS_PLANAR_ANTIPODE ||
                    dot(next.pos, planar_basis->v) < -COS_PLANAR_ANTIPODE;
    }
    const bool use_planar = planar_basis && !antipodal_seam;

    // Advance the rendered-arc accumulator for EVERY segment (drawn or culled) so
    // v0/v1 stay a true full-curve parameterization; seg_base snapshots the start
    // for the draw lambda. Skipped for geodesic polylines.
    if (override_uv) {
      seg_base = cumul;
      cumul += seg_arc_cache[i];
    }

    // Tier 3: Segment culling — skip if the edge's rendered row/column reach
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
          visible = (edge_flags[i] & RasterOptions::EDGE_VISIBLE) != 0;
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
        visible = edge_flags != nullptr
                      ? (edge_flags[i] & RasterOptions::EDGE_VISIBLE) != 0
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
    // exactly as process_segment's fast path (set_arc_uv is a no-op without a
    // planar basis), so plot it without building the sampler. A predicate
    // false negative falls through and re-evaluates exactly.
    const bool one_dot =
        !has_planar_basis &&
        (edge_flags != nullptr &&
                 (edge_flags[i] & RasterOptions::EDGE_CLASSIFIED) != 0
             ? (edge_flags[i] & RasterOptions::EDGE_ONE_DOT) != 0
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
