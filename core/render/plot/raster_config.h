/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include "core/math/geometry.h"
#include "core/render/shading.h"
#include "core/engine/concepts.h"
#include "core/memory.h"
#include "core/render/plot/cull.h"

/**
 * @file raster_config.h
 * @brief rasterize()'s configuration, options and scratch budgets, its
 * sampling-density constants, and the per-step helpers its walk is built from.
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

/**
 * @brief Density-scaled stroke coverage at a chosen sample spacing.
 * @param alpha Coverage at the default spacing.
 * @param step_ratio Sample spacing over the default spacing.
 * @param density Reference samples per unit lattice length.
 * @param balanced Apply balanced gain to short sample spacings.
 * @return Coverage clamped to 1.
 */
static inline float stroke_sample_alpha(float alpha, float step_ratio,
                                        float density, bool balanced) {
  const float coverage = !balanced && step_ratio < 1.0f
                             ? alpha * step_ratio
                             : balanced_sample_alpha(alpha, step_ratio);
  return fminf(1.0f, density * coverage);
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
  /** Marks a column-walk configuration; rasterize() rejects it. */
  bool timed = false;
  /** Recompute deforming edge lengths at each column's time and fold meridians
   *  by half a turn for exact antisymmetry. */
  bool column_edge_lengths = false;
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
 * @details Rows use y_to_screen_row; columns use vector_to_column. Valid only
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

static_assert(sizeof(void *) != 4 || sizeof(RasterOptions) <= 48,
              "RasterOptions exceeds the 32-bit hot-path size budget");

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
  const float pole_guard =
      base_step * MIN_POLE_SCALE * BALANCED_POLE_GUARD_SCALE;
  return default_step <= pole_guard
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

HS_O3_END

} // namespace Plot
