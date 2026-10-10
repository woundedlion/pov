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
/** @brief Single-pass planar samples taken with position and tangent. */
inline uint32_t g_planar_full_samples = 0;
/** @brief Planar samples taken for position only. */
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
  /**
   * @brief A closed loop whose last segment returns to the first point.
   * @param seam Registers for the closing endpoint, or null to reuse the
   *        first point's.
   * @return The closed loop.
   */
  static constexpr RasterLoop closed(const Fragment *seam = nullptr) {
    return RasterLoop(seam);
  }
  /** @return True for a closed loop. */
  constexpr bool is_closed() const { return closed_loop; }
  /** @return The closing-endpoint fragment, or null. */
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
   * @param flags One byte per edge, or empty for none.
   * @return The geodesic projection.
   */
  static constexpr RasterProjection
  geodesic(std::span<const uint8_t> flags = {}) {
    return RasterProjection(nullptr, flags);
  }
  /** @brief Selects azimuthal-equidistant interpolation in the supplied chart.
   * @details Flags contain one byte per edge; planar rasterization reads only
   * EDGE_VISIBLE to replace the clip cull when the clip is active.
   * @param basis Planar chart; must outlive the projection.
   * @param flags One byte per edge, or empty for none.
   * @return The planar projection.
   */
  static constexpr RasterProjection
  planar(const math::Basis &basis, std::span<const uint8_t> flags = {}) {
    return RasterProjection(&basis, flags);
  }
  static RasterProjection planar(math::Basis &&,
                                 std::span<const uint8_t> = {}) = delete;
  static RasterProjection planar(const math::Basis &&,
                                 std::span<const uint8_t> = {}) = delete;
  /** @return The planar chart, or null for geodesic edges. */
  constexpr const math::Basis *basis() const { return planar_basis; }
  /** @return Per-edge visibility bytes, or empty. */
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
  /**
   * @brief Views equal-length row and column arrays.
   * @tparam N Point count.
   * @param rows Screen row per point.
   * @param cols Screen column per point.
   */
  template <size_t N>
  constexpr PointProjections(const float (&rows)[N], const float (&cols)[N])
      : screen_rows(rows), screen_cols(cols), count(N) {}
  /**
   * @brief Views row and column spans; fails fast on a length mismatch.
   * @param rows Screen row per point.
   * @param cols Screen column per point.
   * @return The paired view; null pointers when empty.
   */
  static PointProjections paired(std::span<const float> rows,
                                 std::span<const float> cols) {
    HS_CHECK(rows.size() == cols.size(),
             "hoisted point projection rows and columns differ in length");
    return PointProjections(rows, cols);
  }
  /** @return Screen rows, or null. */
  constexpr const float *rows() const { return screen_rows; }
  /** @return Screen columns, or null. */
  constexpr const float *cols() const { return screen_cols; }
  /** @return Point count. */
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

  /** Open polyline or closed loop. */
  RasterLoop loop{};
  /** Geodesic or planar edges and per-edge flags. */
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
 * else @p dist over @p total_dist clamped to 1.
 * @param terminal The step reaches the segment end.
 * @param dist Arc walked at the step.
 * @param total_dist Segment arc.
 * @return Arc fraction in [0, 1]. */
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
 * projection is singular.
 * @param p Unit point.
 * @param basis Planar chart.
 * @return True past the COS_PLANAR_ANTIPODE cutoff. */
HS_HOT_INLINE inline bool at_planar_antipode(const math::Vector &p,
                                             const math::Basis &basis) {
  return math::dot(p, basis.v) < -COS_PLANAR_ANTIPODE;
}

/**
 * @brief True when @p bit is set in edge @p i's visibility byte.
 * @param edge_flags Per-edge visibility bytes; non-null.
 * @param i Edge index.
 * @param bit Flag mask.
 * @return Whether the bit is set.
 */
HS_HOT_INLINE inline bool has_edge_flag(const uint8_t *edge_flags, size_t i,
                                        uint8_t bit) {
  return (edge_flags[i] & bit) != 0;
}

#include "render/plot/raster_walk.h"

/**
 * @brief True when a balanced walk may take its next step from a position-only
 * sample, reusing the current step.
 * @return Whether the next sample may be position-only.
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
 * @brief Arc-fraction window inside which rasterize shades and plots samples.
 * @details A window that does not narrow [0, 1] admits every sample. A window
 * ending at or past t = 1 has no upper bound: replay can overshoot t = 1 by an
 * ULP.
 */
class PlotWindow {
public:
  /**
   * @param t_start Lower bound, RasterOptions::plot_t_start.
   * @param t_end Upper bound, RasterOptions::plot_t_end.
   */
  PlotWindow(float t_start, float t_end)
      : start(t_start),
        hi(t_end < 1.0f ? t_end : std::numeric_limits<float>::infinity()),
        active(t_start > 0.0f || t_end < 1.0f),
        start_vertex(!active || (start <= 0.0f && hi >= 0.0f)),
        end_vertex(!active || (start <= 1.0f && hi >= 1.0f)) {}

  /**
   * @brief True when the window narrows [0, 1].
   * @return Whether the window is active.
   */
  HS_HOT_INLINE bool is_active() const { return active; }
  /**
   * @brief True when a sample at @p t is plotted.
   * @param t Arc fraction.
   * @return Whether @p t lies in the window.
   */
  HS_HOT_INLINE bool contains(float t) const {
    return !active || (t >= start && t <= hi);
  }
  /**
   * @brief True when a sample at @p t lies strictly outside an active window.
   * @param t Arc fraction.
   * @return Whether @p t is excluded.
   */
  HS_HOT_INLINE bool excludes(float t) const {
    return active && (t < start || t > hi);
  }
  /**
   * @brief True when the segment's start vertex, t = 0, is plotted.
   * @return Whether t = 0 lies in the window.
   */
  HS_HOT_INLINE bool plots_start() const { return start_vertex; }
  /**
   * @brief True when the segment's end vertex, t = 1, is plotted.
   * @return Whether t = 1 lies in the window.
   */
  HS_HOT_INLINE bool plots_end() const { return end_vertex; }
  /**
   * @brief True when the window opens at or before t = 0.
   * @return Whether the window has no effective lower bound.
   */
  HS_HOT_INLINE bool opens_at_start() const { return !active || start <= 0.0f; }
  /**
   * @brief True when the window closes at or past t = 1.
   * @return Whether the window has no effective upper bound.
   */
  HS_HOT_INLINE bool closes_at_end() const { return !active || hi >= 1.0f; }
  /**
   * @brief True when vertex @p k of a one-dot edge falls outside the window.
   * @param k Vertex index: 0 for the start, 1 for the end.
   * @return Whether the vertex is skipped.
   */
  HS_HOT_INLINE bool skips_vertex(size_t k) const {
    return active && ((k == 0 && !plots_start()) || (k == 1 && !plots_end()));
  }

private:
  float start;
  float hi;
  bool active;
  bool start_vertex;
  bool end_vertex;
};

/**
 * @brief Shades and plots a segment endpoint.
 * @tparam INTERP Copy @p vertex's registers into the sample.
 * @param pipeline Pipeline that plots the endpoint.
 * @param canvas Target canvas.
 * @param fragment_shader Fragment shader.
 * @param arc_uv Planar v0/v1 stamp.
 * @param vertex Endpoint control point.
 * @param arc Arc from the segment start to @p vertex.
 * @param alpha_scale Balanced-sampling step ratio; none plots the shaded alpha.
 */
template <bool INTERP, typename PipelineT, typename FragmentShaderT,
          typename StampT>
HS_HOT_INLINE inline void
plot_vertex(PipelineT &pipeline, Canvas &canvas,
            FragmentShaderT &fragment_shader, const StampT &arc_uv,
            const Fragment &vertex, float arc,
            std::optional<float> alpha_scale = std::nullopt) {
  Fragment f;
  seed_fragment<INTERP>(f, vertex);
  arc_uv.stamp(f, arc);
  shade_and_plot(pipeline, canvas, fragment_shader, vertex.pos, f, alpha_scale);
}

/**
 * @brief State of a single-pass adaptive walk along one segment: each sample
 * is plotted as soon as it is taken.
 * @tparam Cfg Rasterizer configuration.
 */
template <RasterConfig Cfg> struct SinglePassWalk {
  /** Cfg selects a sampling policy other than DEFAULT. */
  static constexpr bool NON_DEFAULT_POLICY =
      Cfg.sampling_policy != RasterSamplingPolicy::DEFAULT;
  /** Samples interpolate the endpoints' registers. */
  static constexpr bool INTERPOLATE_REGISTERS = Cfg.interpolate_registers;

  /** Segment arc. */
  float total_dist;
  /** Equatorial step 2π/W. */
  float base_step;
  /** Balanced sampling is selected. */
  bool balanced;
  /** Sample at current_t; position-only after a reused step. */
  SamplePT smp;
  /** Arc walked so far. */
  float current_dist = 0.0f;
  /** Arc fraction of current_dist. */
  float current_t = 0.0f;
  /** Step to the next sample, after the balanced stretch and backstop. */
  float desired_step;
  /** Default-density step at the last full sample. */
  float default_desired_step;
  /** Default-density step at the previous full sample. */
  float previous_full_step;
  /** Tangent at the previous full sample. */
  math::Vector previous_full_tangent;
  /** The next step reuses desired_step from a position-only sample. */
  bool reuse_step = false;
  /** Samples taken. */
  size_t step_count = 0;
  /** Step multiplier set once the walk exhausts its step budget. */
  float backstop_stretch = 1.0f;
  /**
   * Arc from the last plotted sample to the segment end; the terminal step is
   * the remainder, not the stretched desired_step.
   */
  float endpoint_gap;

  /**
   * @param first Full sample at t = 0.
   * @param first_step Default-density step at @p first.
   * @param total_dist Segment arc.
   * @param base_step Equatorial step 2π/W.
   * @param balanced Balanced sampling is selected.
   */
  HS_HOT_INLINE SinglePassWalk(const SamplePT &first, float first_step,
                               float total_dist, float base_step, bool balanced)
      : total_dist(total_dist), base_step(base_step), balanced(balanced),
        smp(first), desired_step(first_step), default_desired_step(first_step),
        previous_full_step(first_step), previous_full_tangent(first.tan),
        endpoint_gap(total_dist) {
    if constexpr (NON_DEFAULT_POLICY) {
      if (balanced)
        desired_step = balanced_step(first_step, base_step);
    }
  }

  /**
   * @brief True while the walk has not reached the segment end.
   * @return Whether current_dist is short of total_dist.
   */
  HS_HOT_INLINE bool unfinished() const { return current_dist < total_dist; }

  /**
   * @brief Writes the unit position of smp to @p p; direct plotting requires
   * unit fragment positions.
   * @param p Output position.
   * @param sample Segment sampler.
   */
  template <typename SampleT>
  HS_HOT_INLINE void unit_position(math::Vector &p, SampleT &sample) const {
    constexpr bool NEWTON_UNIT_SAMPLER =
        requires { std::remove_cvref_t<SampleT>::NEWTON_UNIT; };
    if constexpr (Cfg.open_geodesic || NEWTON_UNIT_SAMPLER) {
      HS_PLOT_COUNT(normalizations);
#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
      if (g_reference_screen_step) {
        p = smp.pos.normalized();
      } else
#endif
        p = newton_unit(smp.pos);
    } else if constexpr (NON_DEFAULT_POLICY &&
                         requires { sample.one_pass(current_t); }) {
      // Balanced walks use the Newton step, default-density walks the
      // exact normalize; the two can plot different rows.
      HS_PLOT_COUNT(normalizations);
      if (balanced) {
        p = newton_unit(smp.pos);
      } else {
        p = smp.pos.normalized();
      }
    } else {
      HS_PLOT_COUNT(normalizations);
      p = smp.pos.normalized();
    }
  }

  /**
   * @brief Shades and plots the sample at @p p.
   * @param pipeline Pipeline that plots the sample.
   * @param canvas Target canvas.
   * @param fragment_shader Fragment shader.
   * @param arc_uv Planar v0/v1 stamp.
   * @param curr Segment start.
   * @param next Segment end.
   * @param p Unit sample position.
   */
  template <typename PipelineT, typename FragmentShaderT, typename StampT>
  HS_HOT_INLINE void plot_sample(PipelineT &pipeline, Canvas &canvas,
                                 FragmentShaderT &fragment_shader,
                                 const StampT &arc_uv, const Fragment &curr,
                                 const Fragment &next,
                                 const math::Vector &p) const {
    Fragment f;
    lerp_fragment<INTERPOLATE_REGISTERS>(f, curr, next, current_t, p);
    arc_uv.stamp(f, current_dist);
    std::optional<float> alpha_scale;
    if constexpr (NON_DEFAULT_POLICY) {
      if (balanced)
        alpha_scale = desired_step / default_desired_step;
    }
    shade_and_plot(pipeline, canvas, fragment_shader, p, f, alpha_scale);
  }

  /**
   * @brief Counts the sample just taken against the step budget.
   * @param max_steps Step budget.
   * @return True when the walk must stop, leaving the segment's tail
   *         unplotted.
   * @details Exhausting the budget once stretches every later step to fit the
   * rest of the segment, as the two-pass replay does; exhausting it twice
   * stops the walk.
   */
  HS_HOT_INLINE bool note_backstop(size_t max_steps) {
    if (++step_count >= max_steps) {
      if (backstop_stretch == 1.0f) {
        HS_PLOT_COUNT(backstops);
        HS_SCAN_METRIC(hs::g_scan_metrics.plot_backstop_hits++);
        backstop_stretch = total_dist / current_dist;
      } else if (step_count >= 2 * max_steps) {
        endpoint_gap = total_dist - current_dist;
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Moves the walk one desired_step along the segment.
   * @return True when the walk is still short of the segment end.
   */
  HS_HOT_INLINE bool advance() {
    advance_distance(current_dist, total_dist, desired_step, endpoint_gap);
    return current_dist < total_dist;
  }

  /**
   * @brief Samples the segment at the new current_dist and sets the next step.
   * @param sample Segment sampler.
   * @param planar_arc_interval Monotonic planar sampler cursor.
   * @param adaptive_sample Takes a full sample at an arc fraction.
   * @param adaptive_step Default-density step at a full sample.
   * @param world_identity The pipeline has no world-space transform.
   * @details A balanced walk whose step barely changed takes the next sample
   * position-only and reuses the step.
   */
  template <typename SampleT, typename AdaptiveSampleT, typename AdaptiveStepT>
  HS_HOT_INLINE void resample(SampleT &sample, int &planar_arc_interval,
                              AdaptiveSampleT &adaptive_sample,
                              AdaptiveStepT &adaptive_step,
                              bool world_identity) {
    current_t = current_dist / total_dist;
    HS_PLOT_COUNT(sim_samples);
    if constexpr (NON_DEFAULT_POLICY && requires {
                    sample.position_monotonic(current_t, planar_arc_interval);
                  }) {
      if (balanced && reuse_step) {
        HS_PLOT_STALL_START(position_start);
        smp.pos = sample.position_monotonic(current_t, planar_arc_interval);
        HS_PLOT_RENDER_COUNT(adaptive_samples);
        HS_PLOT_STALL_STOP(adaptive_sim, position_start);
#if HS_ENABLE_TEST_HOOKS
        ++g_planar_position_samples;
#endif
        reuse_step = false;
      } else {
        smp = adaptive_sample(current_t);
        default_desired_step = adaptive_step(smp);
        if (balanced) {
          reuse_step =
              world_identity &&
              can_reuse_step(smp, default_desired_step, previous_full_step,
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
    if constexpr (NON_DEFAULT_POLICY) {
      if (balanced)
        desired_step = balanced_step(default_desired_step, base_step);
    }
    desired_step *= backstop_stretch;
  }

  /**
   * @brief Shades and plots the segment's final endpoint.
   * @param pipeline Pipeline that plots the endpoint.
   * @param canvas Target canvas.
   * @param fragment_shader Fragment shader.
   * @param arc_uv Planar v0/v1 stamp.
   * @param next Segment end.
   */
  template <typename PipelineT, typename FragmentShaderT, typename StampT>
  HS_HOT_INLINE void plot_terminal(PipelineT &pipeline, Canvas &canvas,
                                   FragmentShaderT &fragment_shader,
                                   const StampT &arc_uv,
                                   const Fragment &next) const {
    std::optional<float> alpha_scale;
    if constexpr (NON_DEFAULT_POLICY) {
      // Gain the endpoint by the arc it actually stands in for, floored
      // at the default step so it never dims below the DEFAULT policy.
      if (balanced)
        alpha_scale = endpoint_alpha_scale(endpoint_gap, default_desired_step,
                                           desired_step);
    }
    plot_vertex<INTERPOLATE_REGISTERS>(pipeline, canvas, fragment_shader,
                                       arc_uv, next, total_dist, alpha_scale);
  }
};

/**
 * @brief True when edge @p i spans at most one screen step.
 * @tparam W,H Rasterization resolution.
 * @param edge_flags Per-edge visibility bytes, or null.
 * @param i Edge index.
 * @param a Edge start.
 * @param b Edge end.
 * @details A classified edge reads EDGE_ONE_DOT; any other edge evaluates
 * edge_fits_one_dot, which may return a false negative.
 * @return Whether the edge plots as a single dot.
 */
template <int W, int H>
HS_HOT_INLINE inline bool edge_spans_one_dot(const uint8_t *edge_flags,
                                             size_t i, const math::Vector &a,
                                             const math::Vector &b) {
  return edge_flags != nullptr &&
                 has_edge_flag(edge_flags, i, RasterOptions::EDGE_CLASSIFIED)
             ? has_edge_flag(edge_flags, i, RasterOptions::EDGE_ONE_DOT)
             : edge_fits_one_dot<W, H>(a, b);
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
  [[maybe_unused]] const bool balanced_sampling =
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
  const PlotWindow window(opts.plot_t_start, opts.plot_t_end);
  HS_CHECK(!window.is_active() || count == 1,
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
    // Coincident endpoints emit at most one dot.
    if (total_dist < math::EPS_GEOMETRIC) {
      bool should_omit = close_loop || !is_last_segment || omit_end;
      if (!should_omit && window.plots_start())
        plot_vertex<INTERPOLATE_REGISTERS>(pipeline, canvas, fragment_shader,
                                           arc_uv, curr, 0.0f);
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
      if (window.plots_start())
        plot_vertex<INTERPOLATE_REGISTERS>(pipeline, canvas, fragment_shader,
                                           arc_uv, curr, 0.0f);
      if (!close_loop && is_last_segment && !omit_end && window.plots_end())
        plot_vertex<INTERPOLATE_REGISTERS>(pipeline, canvas, fragment_shader,
                                           arc_uv, next, total_dist);
      return;
    }

    // Size each sub-step so consecutive samples land ~SCREEN_STEP_PX apart in
    // screen space. `smp`/`first_step` seed the first iteration.
    if constexpr (SINGLE_PASS) {
      HS_PROFILE_DEEP(plot_seg_single_pass);
      SinglePassWalk<Cfg> walk(smp, first_step, total_dist, base_step,
                               balanced_sampling);
      while (walk.unfinished()) {
        math::Vector p;
        walk.unit_position(p, sample);
        if (window.contains(walk.current_t))
          walk.plot_sample(pipeline, canvas, fragment_shader, arc_uv, curr,
                           next, p);
        if (walk.note_backstop(max_cache))
          break;
        HS_PLOT_MAX(steps_peak, walk.step_count);
        if (walk.advance())
          walk.resample(sample, planar_arc_interval, adaptive_sample,
                        adaptive_step, world_identity);
      }
      if (!close_loop && is_last_segment && !omit_end && window.plots_end())
        walk.plot_terminal(pipeline, canvas, fragment_shader, arc_uv, next);
      return;
    }

    // Walks the segment at the adaptive step, caching each step; returns the
    // simulated arc.
    auto simulate_steps = [&]() __attribute__((always_inline)) -> float {
      HS_PROFILE_DEEP(plot_seg_sim);
      float sim_dist = 0.0f;
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
      return sim_dist;
    };

    // Plots the cached steps, each scaled by `scale`; `omit_last` drops the
    // final endpoint.
    auto replay_steps = [&](float scale,
                            bool omit_last) __attribute__((always_inline)) {
      // Normalize interpolated positions before vector_to_pixel's acos(v.y).
      HS_PROFILE_DEEP(plot_seg_draw);
      if (window.plots_start()) {
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

      size_t loop_limit =
          omit_last ? steps_cache.size() - 1 : steps_cache.size();
      float current_dist = 0.0f;

      for (size_t j = 0; j < loop_limit; j++) {
        float step = steps_cache[j] * scale;
        current_dist += step;

        const bool terminal = !omit_last && j == loop_limit - 1;
        float t = replay_t(terminal, current_dist, total_dist);

        if (window.excludes(t))
          continue;

        // `t` follows the rendered arc length. Under a planar basis
        // `PlanarArcStamp::stamp` rewrites the lerped v0/v1 from the sampled
        // arc, so arc-keyed shaders track the drawn position.
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

    steps_cache.clear();
    const float sim_dist = simulate_steps();
    // scale <= 1 normally (the final step overshoots); > 1 after a backstop
    // break. Either way the replay spans exactly total_dist.
    HS_CHECK(sim_dist > 0.0f,
             "rasterize: simulated segment length is not positive");
    replay_steps(total_dist / sim_dist,
                 close_loop || !is_last_segment || omit_end);
  };

  const auto &cr = canvas.clip();
  const bool clip_active = !cr.is_full();
  const auto xc = cr.x_clip();

  // Emits one shader-run dot for points[k]; the precomputed projection is
  // consumed only when no world stage would lift it back to a world vector.
  auto plot_dot = [&](const Fragment &src, size_t k) {
    if (window.skips_vertex(k))
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
        edge_spans_one_dot<W, H>(edge_flags, i, curr.pos, next.pos);
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
