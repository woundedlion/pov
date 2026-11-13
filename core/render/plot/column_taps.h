/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/plot/edge_worldsheet.h"

/** @file column_taps.h @brief Column walk state and timed sample taps. */

namespace Plot {

HS_O3_BEGIN

/** @brief Limits a walk to a column arc, keeping its hull when the arc's gap
 * lies inside the walk.
 * @tparam W Display width.
 * @param walk Walk to limit.
 * @param direction Walk direction.
 * @param begin,length Column arc, ascending from begin; length W or more
 * leaves the walk unchanged. */
template <int W>
HS_HOT_INLINE void clamp_walk(detail::ColumnSegment &walk, int direction,
                              int begin, int length) {
  if (length >= W)
    return;
  const int head =
      direction > 0 ? begin : math::fast_wrap(begin + length - 1, W);
  const int into = math::fast_wrap(direction * (walk.start - head), W);
  if (into < length) {
    if (walk.count <= W - into)
      walk.count = std::min(walk.count, length - into);
    return;
  }
  const int gap = W - into;
  if (gap >= walk.count) {
    walk.count = 0;
    return;
  }
  walk.start = head;
  walk.count = std::min(walk.count - gap, length);
  walk.deadline -= gap;
}

/** @brief One edge's column walk: what its plain and split walks share.
 * @tparam W Display width. */
template <int W> struct ColumnWalk {
  PipelineRef &pipeline;                 /**< Destination pipeline. */
  Canvas &canvas;                        /**< Frame context. */
  const Fragments &points;               /**< The edge's endpoint registers. */
  TimedShaderRef shader;                 /**< Time-aware source shader. */
  const EdgeWorldsheet &frame;           /**< The edge's paths and bases. */
  const detail::ColumnSegment *segments; /**< Walks per arm pass. */
  int segment_count;                     /**< Number of walks. */
  const ClipRegion::XClip &band;         /**< Render band. */
  int direction;                         /**< Walk direction. */
  float radians;                         /**< Maximum world displacement. */
  bool single_time;      /**< End columns splat at their own time too. */
  const float *velocity; /**< Endpoint motion, pixels per frame. */
  bool source_positions; /**< Shaders read source-space positions. */
  bool chained;          /**< The end vertex is left unplotted. */
  bool turning;   /**< The arc can turn its longitude direction this frame. */
  float top_y;    /**< Highest emitted unit-sphere y. */
  float bottom_y; /**< Lowest emitted unit-sphere y. */
  float
      half_sine; /**< Sine of the window half-width: half a column, or a full column for local lattices. */
  float
      half_cosine; /**< Cosine of the window half-width: half a column, or a full column for local lattices. */
};

/** @brief A walked sample whose tap columns splat at their own times. */
struct ColumnTapSample {
  math::PixelCoords p; /**< Projected position at its column's time. */
  float t;             /**< Its column's sweep fraction. */
  float s;             /**< Arc fraction along the edge. */
};

/** @brief Plots a sample's taps in one column, moved to that column's time.
 * @tparam W Canvas width.
 * @param state The edge's walk.
 * @param sample The sample's position, time and motion.
 * @param column Wrapped tap column.
 * @param at The column's sweep fraction.
 * @param color,age,alpha Registers at @p at. */
template <int W>
HS_HOT_INLINE void plot_column_tap(const ColumnWalk<W> &state,
                                   const ColumnTapSample &sample, int column,
                                   float at, const Pixel &color, float age,
                                   float alpha) {
  if (at != sample.t)
    Spacetime::count_geometry<
        &Spacetime::GeometryCounters::column_boundary_samples>();
  const float *velocity = state.velocity;
  const float shift = at - sample.t;
  const float x =
      sample.p.x +
      (velocity[0] + (velocity[2] - velocity[0]) * sample.s) * shift;
  state.pipeline.plot_column_at(
      state.canvas,
      {x < 0.0f ? x + W
       : x >= W ? x - W
                : x,
       sample.p.y +
           (velocity[1] + (velocity[3] - velocity[1]) * sample.s) * shift},
      column, color, age, alpha, at);
}

/** @brief Prepares a sample's registers at the other side of an arm seam.
 * @tparam Cfg Register interpolation configuration.
 * @tparam UNIFORM_COLOR Reuse the prepared edge color.
 * @param points Edge endpoint registers.
 * @param sample Owner position, time and material fraction.
 * @param owner Registers shaded at the owner's time.
 * @param far Other tap's time.
 * @param shader Time-aware shader.
 * @return Registers for the other tap. */
template <RasterConfig Cfg, bool UNIFORM_COLOR>
HS_HOT_INLINE Fragment column_seam_registers(const Fragments &points,
                                             const ColumnTapSample &sample,
                                             const Fragment &owner, float far,
                                             const TimedShaderRef &shader) {
  Fragment other;
  if constexpr (Cfg.interpolate_registers)
    other = Fragment::lerp_registers(points[0], points[1], sample.s);
  if (far != sample.t) {
    if constexpr (UNIFORM_COLOR)
      other.color = owner.color;
    else
      shader(owner.pos, other, far);
  } else {
    other = owner;
  }
  return other;
}

/** @brief Splats a still sample's taps astride an arm seam.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam UNIFORM_COLOR Reuse the prepared edge color.
 * @param state The edge's walk.
 * @param sample The sample's position and time.
 * @param f The sample's registers, shaded at its time.
 * @details The taps are lit a frame apart: their registers interpolate in
 * time between the sample and one shade at the farthest tap's time. */
template <int W, int H, RasterConfig Cfg, bool UNIFORM_COLOR = false>
HS_HOT_FLASH_MEMBER void splat_seam_taps(const ColumnWalk<W> &state,
                                         const ColumnTapSample &sample,
                                         const Fragment &f) {
  using Splat = Filter::Screen::AntiAlias<W, H>;
  const float t = sample.t;
  float far = t;
  Splat::for_each_seam_column(sample.p.x, [&](int, float at) HS_HOT_INLINE {
    if (fabsf(at - t) > fabsf(far - t))
      far = at;
  });
  const Fragment other = column_seam_registers<Cfg, UNIFORM_COLOR>(
      state.points, sample, f, far, state.shader);
  const float span = far != t ? 1.0f / (far - t) : 0.0f;
  Splat::for_each_seam_column(
      sample.p.x, [&](int column, float at) HS_HOT_INLINE {
        const float w = (at - t) * span;
        const Color4 color = f.color.lerp(other.color, w);
        plot_column_tap(state, sample, column, at, color.color,
                        f.age + (other.age - f.age) * w, color.alpha);
      });
}

/** @brief Splats each of a moving sample's tap columns within half a frame
 *         of it at that column's time.
 * @tparam W,H Canvas dimensions.
 * @param state The edge's walk.
 * @param sample The sample's position, time and motion.
 * @param f The sample's registers. */
template <int W, int H>
HS_HOT_FLASH_MEMBER void splat_owned_taps(const ColumnWalk<W> &state,
                                          const ColumnTapSample &sample,
                                          const Fragment &f) {
  Filter::Screen::AntiAlias<W, H>::for_each_owned_column(
      sample.p.x, sample.t, [&](int column, float at) HS_HOT_INLINE {
        plot_column_tap(state, sample, column, at, f.color.color, f.age,
                        f.color.alpha);
      });
}

HS_O3_END

} // namespace Plot
