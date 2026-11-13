/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/plot/column_edge_walk.h"

/** @file column_edge_batch.h @brief Column-time mesh edge batch rasterization. */

namespace Plot {

HS_O3_BEGIN

/** @brief Rasterizes mesh edges one destination column at a time.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam RECORD_WORLD Record samples in source history.
 * @param pipeline Destination pipeline; its world stages are @p poses.
 * @param canvas Frame context.
 * @param points The walk's two endpoint fragments, written by
 * EdgeBatch::prepare.
 * @param shader Time-aware source shader.
 * @param batch Edges and their vertices' paths.
 * @param poses The world stages' poses by column time, or null without world
 * stages.
 * @param radians Maximum world displacement over the frame.
 * @param source_positions Whether shaders read source-space positions.
 * @param spin Rotation bound of the shared poses, radians.
 * @param setup The canvas's column-walk constants.
 * @details Each arm pass walks, in the edge's travel order, the render-band
 * columns its plan bounds, each with the edge's endpoints evaluated on their
 * paths at the column's time and rotated by its pose, and stops where the lit
 * run ends. Samples sit at arc fractions k/N, where N is the edge's length in
 * screen steps at that time, resolved locally when midpoint spacing exceeds
 * the bilinear footprint. A column plots the lattice points between
 * its two boundary meridians unweighted through the filter chain, so a moving
 * edge samples each column as the frozen edge at that column's time does. The
 * start vertex is lattice point zero; each edge also plots its end vertex
 * unless the batch is chained.
 * Columns strictly inside an edge splat at their own time; the end columns of
 * an edge moving faster than COLUMN_SINGLE_TAP_ERROR move each neighbor tap
 * to that column's time, and taps astride an arm seam interpolate their
 * registers toward the sample shaded at the far side's time. A stationary draw
 * with bounded midpoint spacing walks the lattice once. An
 * edge whose plan cannot bound it, such as one that can pass near a pole, is
 * walked as the two sides of its closest approach to the pole, out of line. */
template <int W, int H, RasterConfig Cfg, bool RECORD_WORLD>
HS_HOT_INLINE void rasterize_column_edges_impl(
    PipelineRef &pipeline, Canvas &canvas, Fragments &points,
    TimedShaderRef shader, const EdgeBatch &batch,
    const Spacetime::PoseTable *poses, float radians, bool source_positions,
    float spin, const ColumnWalkFrame<W, H> &setup) {
  const float half_sine = setup.half_sine;
  const float half_cosine = setup.half_cosine;
  const float full_sine = setup.full_sine;
  const float full_cosine = setup.full_cosine;

  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  const float top_y = setup.top_y;
  const float bottom_y = setup.bottom_y;
  const float physical_top_y = RECORD_WORLD ? top_y : setup.physical_top_y;
  Spacetime::count_geometry<&Spacetime::GeometryCounters::column_draw_setups>();
  const SweptEdgeCull<W> &cull = setup.cull;
  const auto &band = canvas.render_x_clip();
  ScratchScope draw_scope(scratch_arena_a);
  source_positions = source_positions && poses;
  auto &frame = *scratch_arena_a.make<EdgeWorldsheet>();
  for (size_t edge = 0; edge < batch.count; ++edge) {
    const MeshEdge &pair = batch.pairs[edge];
    HS_CHECK(pair.u < batch.motion.size() && pair.v < batch.motion.size(),
             "column walk: edge (%d, %d) outside vertex count %d", pair.u,
             pair.v, static_cast<int>(batch.motion.size()));
    if (!frame.prepare<W, H>(batch.motion, pair.u, pair.v, poses, spin,
                             source_positions, !batch.rigid_lattice, &cull)) {
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::column_edges_culled>();
      continue;
    }
    const bool local_lattice =
        batch.rigid_lattice || (pipeline.column_taps_ready() &&
                                (EdgeWorldsheet::needs_local_lattice<W, H>(
                                     frame.bases[0], frame.pieces[0]) ||
                                 EdgeWorldsheet::needs_local_lattice<W, H>(
                                     frame.bases[1], frame.pieces[1]) ||
                                 fmaxf(frame.pieces[0], frame.pieces[1]) >
                                     2 * rasterize_step_budget<W>()));
    const int direction =
        frame.bases[1].longitude_rate_numerator() >= 0.0f ? 1 : -1;
    detail::ColumnSegment segments[2];
    int segment_count = 1;
    float velocity[4] = {};
    bool single_time = true;
    if (radians == 0.0f && !local_lattice) {
      segments[0] = {0, 1, 1};
    } else {
      segment_count = frame.plan<W>(direction, segments, local_lattice);
      int kept = 0;
      for (int i = 0; i < segment_count; ++i) {
        detail::ColumnSegment walk = segments[i];
        clamp_walk<W>(walk, direction, band.rs, band.length(W));
        if (walk.count > 0)
          segments[kept++] = walk;
      }
      segment_count = kept;
      if (segment_count == 0) {
        Spacetime::count_geometry<
            &Spacetime::GeometryCounters::column_edges_culled>();
        continue;
      }
      if (!local_lattice) {
        frame.velocities<W, H>(velocity);
        single_time = fmaxf(fmaxf(fabsf(velocity[0]), fabsf(velocity[1])),
                            fmaxf(fabsf(velocity[2]), fabsf(velocity[3]))) <=
                      COLUMN_SINGLE_TAP_ERROR * (W / 2);
      }
    }
    float walk_top_y = top_y;
    if constexpr (!RECORD_WORLD) {
      if (local_lattice || single_time) {
        walk_top_y = fminf(top_y, physical_top_y);
        if (walk_top_y < top_y &&
            !frame.reaches_rows<W>(walk_top_y, bottom_y)) {
          Spacetime::count_geometry<
              &Spacetime::GeometryCounters::column_edges_culled>();
          continue;
        }
      }
    }
    Spacetime::count_geometry<
        &Spacetime::GeometryCounters::column_edges_prepared>();
    batch.prepare(points, batch.vertices[pair.u], batch.vertices[pair.v], edge,
                  frame.bases[1].total_dist);
    points[0].age += batch.age;
    points[1].age += batch.age;
    if (batch.uniform_color) {
      const math::Vector saved = points[0].pos;
      points[0].pos = source_positions ? EdgeWorldsheet::normalize(saved)
                                       : frame.bases[1].v1;
      shader(points[0].pos, points[0], 1.0f);
      points[0].pos = saved;
      points[1].color = points[0].color;
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::column_shader_calls>();
    }
    const auto sense = [&](const GeodesicEdgeSampler &basis) HS_HOT_INLINE {
      return basis.longitude_rate_numerator();
    };
    const bool split = radians != 0.0f && !frame.known;
    const ColumnWalk<W> walk{
        .pipeline = pipeline,
        .canvas = canvas,
        .points = points,
        .shader = shader,
        .frame = frame,
        .segments = segments,
        .segment_count = segment_count,
        .band = band,
        .direction = direction,
        .radians = radians,
        .single_time = single_time,
        .velocity = velocity,
        .source_positions = source_positions,
        .chained = batch.chained,
        .turning = local_lattice ||
                   (split && !(sense(frame.bases[0]) * sense(frame.bases[1]) >
                               4.0f * frame.deviation * frame.deviation)),
        .top_y = walk_top_y,
        .bottom_y = bottom_y,
        .half_sine = local_lattice ? full_sine : half_sine,
        .half_cosine = local_lattice ? full_cosine : half_cosine};
    if (batch.uniform_color) {
      if (local_lattice)
        walk_uniform_edge_columns<W, H, Cfg, false, true, RECORD_WORLD>(walk);
      else if (split)
        walk_uniform_edge_columns<W, H, Cfg, true, false, RECORD_WORLD>(walk);
      else
        walk_uniform_edge_columns<W, H, Cfg, false, false, RECORD_WORLD>(walk);
    } else if (local_lattice)
      walk_local_edge_columns<W, H, Cfg, RECORD_WORLD>(walk);
    else if (split)
      walk_split_edge_columns<W, H, Cfg, RECORD_WORLD>(walk);
    else
      walk_edge_columns<W, H, Cfg, false, false, false, RECORD_WORLD>(walk);
  }
}

/** @brief Walks edges without source-history recording.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @param pipeline,canvas,points Destination and endpoint registers.
 * @param shader Source shader.
 * @param batch Captured source edges.
 * @param poses World poses by column time.
 * @param radians Maximum world displacement.
 * @param source_positions Whether shaders read source positions.
 * @param spin World rotation bound.
 * @param setup The canvas's column-walk constants. */
template <int W, int H, RasterConfig Cfg>
FASTRUN HS_O3_FN HS_NOINLINE_NOCLONE void
rasterize_column_edges(PipelineRef &pipeline, Canvas &canvas, Fragments &points,
                       TimedShaderRef shader, const EdgeBatch &batch,
                       const Spacetime::PoseTable *poses, float radians,
                       bool source_positions, float spin,
                       const ColumnWalkFrame<W, H> &setup) {
  rasterize_column_edges_impl<W, H, Cfg, false>(pipeline, canvas, points,
                                                shader, batch, poses, radians,
                                                source_positions, spin, setup);
}

/** @brief Walks edges while recording source history.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @param pipeline,canvas,points Destination and endpoint registers.
 * @param shader Source shader.
 * @param batch Captured source edges.
 * @param poses World poses by column time.
 * @param radians Maximum world displacement.
 * @param source_positions Whether shaders read source positions.
 * @param spin World rotation bound.
 * @param setup The canvas's column-walk constants. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void rasterize_recording_column_edges(
    PipelineRef &pipeline, Canvas &canvas, Fragments &points,
    TimedShaderRef shader, const EdgeBatch &batch,
    const Spacetime::PoseTable *poses, float radians, bool source_positions,
    float spin, const ColumnWalkFrame<W, H> &setup) {
  rasterize_column_edges_impl<W, H, Cfg, true>(pipeline, canvas, points, shader,
                                               batch, poses, radians,
                                               source_positions, spin, setup);
}

HS_O3_END

} // namespace Plot
