/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>

#include "core/math/3dmath.h"
#include "core/memory.h"
#include "core/platform/attributes.h"
#include "core/render/canvas.h"
#include "core/render/plot/mesh.h"
#include "core/render/plot/point_route.h"

/**
 * @file column_chain.h
 * @brief Plot::ColumnChainDraw: chained edges drawn at their columns' times
 * through a point route and column-walk chunks.
 */

namespace Plot {

/**
 * @brief Draws chains of edges at their columns' times: steady short edges
 * through a point route, the rest batched into column-walk chunks.
 * @tparam W,H Canvas dimensions.
 * @tparam Capacity Chunk vertex capacity.
 * @tparam PipelineT Column-walk pipeline type.
 * @tparam Attributes Per-vertex chunk attribute types.
 * @details The chunk's storage is reclaimed by the caller's scratch scope;
 * the pipeline and canvas outlive the draw.
 */
template <int W, int H, size_t Capacity, typename PipelineT,
          typename... Attributes>
class ColumnChainDraw {
public:
  Mesh::ColumnChunk<Capacity, Attributes...>
      chunk;                 /**< Column-walk chunk of the chained edges. */
  GeometryCallback geometry; /**< Vertex motion callback; chained, handing the
                                shaders world positions. */

  /**
   * @brief Allocates the chunk and prepares the point route and column walk.
   * @param arena Chunk storage.
   * @param pipeline Column-walk pipeline.
   * @param canvas Destination canvas.
   * @param reach Bound on a point-routed sample's displacement from its frame
   * end, radians.
   * @param prepared_poses Shared column-time poses for a screen sink.
   * @param prepared_sweep Rotation bound of the shared poses, radians.
   */
  HS_HOT_INLINE
  ColumnChainDraw(Arena &arena, PipelineT &pipeline, Canvas &canvas,
                  float reach,
                  const Spacetime::PoseTable *prepared_poses = nullptr,
                  float prepared_sweep = 0.0f)
      : chunk(arena), setup(Mesh::column_draw_setup<W, H>(
                          pipeline, canvas, prepared_poses, prepared_sweep)),
        pipeline(pipeline), canvas(canvas) {
    geometry.source_positions = false;
    geometry.chained = true;
    route.prepare(canvas, reach);
  }

  /**
   * @brief Plots a short steady edge or appends it to the chunk.
   * @tparam RESIDENT Selects the ITCM point-route entry.
   * @tparam Sink Point-route pipeline type.
   * @tparam Shade Callable `(s, time)` or `(s)` returning Color4.
   * @tparam Store Callable `(size_t slot, size_t source)` filling vertex data.
   * @tparam Flush Callable that draws and empties the chunk.
   * @param sink Point-route pipeline.
   * @param source Edge indices in the caller's current chain.
   * @param starts,ends Unit world endpoints at frame start and end.
   * @param most Largest steady lattice plotted through the point route.
   * @param shade Edge sample color.
   * @param store Vertex data writer.
   * @param flush Chunk draw operation.
   * @param end_to_start Known rigid rotation from frame end to start.
   * @param frame_constant Allows differing endpoint counts, using their maximum.
   */
  template <bool RESIDENT = false, typename Sink, typename Shade,
            typename Store, typename Flush>
  HS_HOT_INLINE void
  plot_edge(Sink &sink, Mesh::Edge source, const math::Vector (&starts)[2],
            const math::Vector (&ends)[2], int most, const Shade &shade,
            const Store &store, const Flush &flush,
            const math::RotationMatrix *end_to_start = nullptr,
            bool frame_constant = false) {
    chunk.template plot_edge<W, H, RESIDENT>(
        route, sink, canvas, source, starts, ends, most, shade, store, flush,
        end_to_start, frame_constant);
  }

  /**
   * @brief Shades from interpolated vertex coordinates and draws the chunk.
   * @tparam Coordinate Callable `(size_t slot) -> float`.
   * @tparam Shade Callable `(size_t head, float coordinate, float time) ->
   * Color4`.
   * @param coordinate Vertex palette coordinate.
   * @param shade Interpolated edge sample color.
   * @param bound Callback motion bound in radians.
   */
  template <typename Coordinate, typename Shade>
  HS_HOT_INLINE void flush(const Coordinate &coordinate, const Shade &shade,
                           float bound) {
    chunk.template flush<W, H>(pipeline, canvas, setup, coordinate, shade,
                               geometry, bound);
  }

  /**
   * @brief Shades from edge endpoint slots and draws the chunk.
   * @tparam Shade Callable `(size_t head, size_t tail, float s, float time) ->
   * Color4`.
   * @param shade Edge sample color.
   * @param bound Callback motion bound in radians.
   */
  template <typename Shade>
  HS_HOT_INLINE void flush(const Shade &shade, float bound) {
    chunk.template flush<W, H>(pipeline, canvas, setup, shade, geometry, bound);
  }

private:
  const Mesh::ColumnDrawSetup<W, H>
      setup;              /**< Frame-fixed column-walk inputs. */
  PointRoute<W, H> route; /**< Prepared steady-edge point route. */
  PipelineT &pipeline;    /**< Column-walk pipeline. */
  Canvas &canvas;         /**< Destination canvas. */
};

} // namespace Plot
