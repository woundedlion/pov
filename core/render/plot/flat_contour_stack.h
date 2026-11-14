/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file flat_contour_stack.h
 * @brief Concentric contours drawn at the frame's end orientation.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>

#include "core/math/3dmath.h"
#include "core/memory.h"
#include "core/platform/attributes.h"
#include "core/platform/profiling.h"
#include "core/render/canvas.h"
#include "core/render/clip.h"
#include "core/render/filter/screen_direct_aa_sink.h"
#include "core/render/plot/chords.h"
#include "core/render/plot/contour_layout.h"
#include "core/render/plot/raster.h"
#include "core/render/plot/shapes.h"

namespace Plot {

/**
 * @brief Draws a ContourLayout's stack at one orientation, every contour
 *        sampled and rasterized as a Plot primitive.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details From DENSE_CONTOUR_COUNT contours up, planar stars stroke through
 * Plot::PlanarChords and spherical stars take balanced sampling. Flowers whose
 * petals cross the clip band's columns are split against it by
 * Plot::PlanarBandSplit.
 */
template <int W, int H> class FlatContourStack {
public:
  /** @brief Each of the 2 * sides flower edges takes FLOWER_PIECE_BUDGET /
   *  sides pieces, for at most 2 * FLOWER_PIECE_BUDGET pieces per flower. */
  static constexpr int FLOWER_PIECE_BUDGET = 24;
  /** @brief Most points a split flower takes. */
  static constexpr int FLOWER_MAX_POINTS =
      PlanarBandSplit<W, H>::max_points(2, FLOWER_PIECE_BUDGET);
  static_assert(CONTOUR_MAX_VERTICES <= 2 * FLOWER_PIECE_BUDGET);

  /** @brief Raster configuration of every sampled contour. */
  static constexpr RasterConfig SAMPLED_RASTER_CONFIG{
      .single_pass = true,
      .derive_planar_arc_registers = false,
      .interpolate_registers = false,
      .sampling_policy = RasterSamplingPolicy::SELECTABLE};
  static_assert(SAMPLED_RASTER_CONFIG.single_pass &&
                !SAMPLED_RASTER_CONFIG.derive_planar_arc_registers);

  /** @brief Arena bytes init_storage() takes.
   *  @return Byte count. */
  static constexpr size_t storage_bytes() {
    return PlanarBandSplit<W, H>::storage_bytes(FLOWER_MAX_POINTS);
  }

  /** @brief Peak scratch_arena_a bytes draw() takes, including nested raster
   *  storage.
   *  @return Byte count. */
  static constexpr size_t scratch_a_bytes() {
    return std::max(
        (CONTOUR_MAX_VERTICES + 2) * sizeof(Fragment) + alignof(Fragment) +
            PlanarChords<W, H>::scratch_a_bytes(),
        (CONTOUR_MAX_VERTICES + 1 + FLOWER_MAX_POINTS) * sizeof(Fragment) +
            2 * alignof(Fragment));
  }

  /** @brief Binds the flower split's flag storage.
   *  @param arena Arena the storage is allocated from. */
  HS_FLASH_MEMBER void init_storage(Arena &arena) {
    flower_split.init_storage(arena, FLOWER_MAX_POINTS);
  }

  /**
   * @brief Draws the contour stack at the frame's end orientation.
   * @tparam MAX Layout contour capacity.
   * @param canvas Target canvas.
   * @param sink Splat sink every stroke lands in.
   * @param chords Planar chord walker, bound for CONTOUR_MAX_VERTICES.
   * @param layout Prepared radii, waveform and draw order.
   * @param frame Shape, counts, twist, opacity, palette and orientation.
   */
  template <int MAX>
  HS_HOT_FLASH_MEMBER void
  draw(Canvas &canvas, Filter::Screen::DirectAntiAliasSink<W, H> &sink,
       PlanarChords<W, H> &chords, const ContourLayout<MAX> &layout,
       const ContourFrame &frame) {
    const int count = frame.count;
    const int sides = frame.sides;
    const ContourShape shape = frame.shape;
    const bool dense_contours =
        static_cast<float>(count) >= DENSE_CONTOUR_COUNT;
    const bool planar_star = shape == ContourShape::PLANAR_STAR;
    const math::Basis basis =
        math::make_basis(frame.orientation.get(), math::X_AXIS);
    const ClipRegion &clip = canvas.clip();
    CapCenter near_cap{}, far_cap{};
    if (planar_star) {
      near_cap = make_cap_center(clip, basis.v);
      far_cap = make_cap_center(clip, -basis.v);
    }

    if (planar_star && dense_contours)
      chords.prepare(clip);
    if (shape == ContourShape::FLOWER)
      flower_band = ClipBand<W, H>::of(clip);

    Color4 pair_color;
    const float global_alpha = frame.opacity;
    const bool continuous_star = shape == ContourShape::SPHERICAL_STAR;
    for (int ordinal = 0; ordinal < count; ++ordinal) {
      const ContourSlot slot = layout.slot(ordinal, count, continuous_star,
                                           frame.palette, pair_color);
      const int i = slot.i;
      const float radius = slot.radius;
      if (planar_star &&
          !planar_star_cap_visible<W, H>(clip, near_cap, far_cap, radius, 0.0f))
        continue;
      const ContourShading shading =
          layout.shading(i, radius, continuous_star, frame.opposite,
                         frame.amplitude, pair_color, global_alpha);
      const float contour_phase = shading.contour_phase;
      const Color4 shaded_color = shading.shaded_color;
      auto shader = [&](const math::Vector &, Fragment &fragment) {
        fragment.color = shaded_color;
      };
      dispatch_plot(sink, chords, layout, canvas, basis, shape, radius, sides,
                    shader, contour_phase, shaded_color, i, dense_contours);
    }
  }

private:
  /**
   * @brief Plot-rasterizes an open polyline sampled into scratch storage.
   * @tparam F Fragment-shader callable type.
   * @param sink Splat sink.
   * @param canvas Target canvas.
   * @param capacity Fragment slots to bind for the sampler.
   * @param planar_basis Azimuthal-equidistant chart for the edges, or nullptr
   * for geodesic edges.
   * @param balanced_sampling Selects the balanced raster sampling policy.
   * @param fragment_shader Per-fragment shader.
   * @param fill Callable that samples the primitive into the bound fragments.
   */
  template <typename F>
  HS_HOT_FLASH_MEMBER void
  draw_sampled(Filter::Screen::DirectAntiAliasSink<W, H> &sink, Canvas &canvas,
               size_t capacity, const math::Basis *planar_basis,
               bool balanced_sampling, const F &fragment_shader, auto &&fill) {
    ScratchScope guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, capacity);
    fill(points);
    rasterize<W, H, SAMPLED_RASTER_CONFIG>(
        sink, canvas, points, fragment_shader,
        {.projection = planar_basis ? RasterProjection::planar(*planar_basis)
                                    : RasterProjection{},
         .omit_end = true,
         .balanced_sampling = balanced_sampling});
  }

  /**
   * @brief Strokes a planar star with Plot::PlanarChords.
   * @tparam MAX Layout contour capacity.
   * @tparam F Fragment-shader callable type.
   * @param sink Splat sink.
   * @param chords Planar chord walker.
   * @param layout Layout holding the contour's radius trig.
   * @param canvas Target canvas.
   * @param basis Shared shape basis.
   * @param radius Shape radius in [0, 2].
   * @param sides Star point count.
   * @param fragment_shader Per-fragment shader, used by the pole runs.
   * @param color Contour color.
   * @param phase Star rotation in radians.
   * @param contour_index Contour slot in the layout's radius-trig table.
   */
  template <int MAX, typename F>
  HS_FLASH_MEMBER void
  draw_dense_planar_star(Filter::Screen::DirectAntiAliasSink<W, H> &sink,
                         PlanarChords<W, H> &chords,
                         const ContourLayout<MAX> &layout, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const F &fragment_shader, const Color4 &color,
                         float phase, int contour_index) {
    ScratchScope guard(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, static_cast<size_t>(sides * 2 + 2));
    math::Basis projection_basis;
    const math::Basis &planar_basis =
        *PlanarProjection::edge_basis(basis, radius, projection_basis);
    Star<PlanarProjection>::sample_chart_positions(
        points, chords.chart_x(), chords.chart_y(), basis, radius, sides, phase,
        layout.radius_trig(contour_index), layout.star_step_trig(),
        planar_basis);
    chords.draw_closed(sink, canvas, points, sides * 2, planar_basis, color,
                       fragment_shader);
  }

  /**
   * @brief Rasterizes a flower with its petal edges split against the clip
   *        band by Plot::PlanarBandSplit.
   * @tparam F Fragment-shader callable type.
   * @param sink Splat sink.
   * @param canvas Target canvas.
   * @param basis Shared shape basis.
   * @param planar_basis The flower's azimuthal-equidistant chart.
   * @param radius Shape radius in [0, 2].
   * @param sides Petal count.
   * @param fragment_shader Per-fragment shader.
   * @param phase Flower rotation in radians.
   */
  template <typename F>
  HS_FLASH_MEMBER void
  draw_banded_flower(Filter::Screen::DirectAntiAliasSink<W, H> &sink,
                     Canvas &canvas, const math::Basis &basis,
                     const math::Basis &planar_basis, float radius, int sides,
                     const F &fragment_shader, float phase) {
    const int edges = sides * 2;
    const int pieces = std::max(1, FLOWER_PIECE_BUDGET / sides);
    ScratchScope guard(scratch_arena_a);
    Fragments ring;
    ring.bind(scratch_arena_a, static_cast<size_t>(edges + 1));
    Flower::sample(ring, basis, radius, sides, phase);
    Fragments path;
    path.bind(
        scratch_arena_a,
        static_cast<size_t>(PlanarBandSplit<W, H>::max_points(edges, pieces)));
    const std::span<const uint8_t> flags = flower_split.split(
        path, ring, edges, pieces, planar_basis, flower_band);
    rasterize<W, H, SAMPLED_RASTER_CONFIG>(
        sink, canvas, path, fragment_shader,
        {.projection = RasterProjection::planar(planar_basis, flags),
         .omit_end = true});
  }

  /**
   * @brief Samples the selected primitive and draws it.
   * @tparam MAX Layout contour capacity.
   * @tparam F Fragment-shader callable type.
   * @param sink Splat sink.
   * @param chords Planar chord walker.
   * @param layout Layout holding the contour's radius trig.
   * @param canvas Target canvas.
   * @param basis Shared shape basis.
   * @param shape Primitive to draw.
   * @param radius Shape radius in [0, 2]; above 1 the shape wraps past the
   * pole and is charted about the antipode.
   * @param sides Polygon side, flower petal, or star point count.
   * @param fragment_shader Per-fragment shader.
   * @param shape_phase Primitive rotation in radians.
   * @param shape_color Color applied to this shape.
   * @param contour_index Index of this contour in the stack.
   * @param dense_contours Whether the stack is at or above DENSE_CONTOUR_COUNT
   * contours, selecting the screen-step-balanced star paths.
   * @details Cold (flash); the hot work is inside Plot::rasterize.
   */
  template <int MAX, typename F>
  HS_FLASH_MEMBER void
  dispatch_plot(Filter::Screen::DirectAntiAliasSink<W, H> &sink,
                PlanarChords<W, H> &chords, const ContourLayout<MAX> &layout,
                Canvas &canvas, const math::Basis &basis, ContourShape shape,
                float radius, int sides, const F &fragment_shader,
                float shape_phase, const Color4 &shape_color, int contour_index,
                bool dense_contours) {
    HS_PROFILE(ss_plot_dispatch);
    switch (shape) {
    case ContourShape::PLANAR_POLYGON: {
      math::Basis projection_basis;
      const math::Basis &planar_basis =
          *PlanarProjection::edge_basis(basis, radius, projection_basis);
      draw_sampled(sink, canvas, static_cast<size_t>(sides + 2), &planar_basis,
                   false, fragment_shader, [&](Fragments &points) {
                     Polygon<PlanarProjection>::sample(points, basis, radius,
                                                       sides, shape_phase);
                   });
      break;
    }
    case ContourShape::SPHERICAL_POLYGON:
      draw_sampled(sink, canvas, static_cast<size_t>(sides + 2), nullptr, false,
                   fragment_shader, [&](Fragments &points) {
                     Polygon<GeodesicProjection>::sample(points, basis, radius,
                                                         sides, shape_phase);
                   });
      break;
    case ContourShape::FLOWER: {
      math::Basis planar_basis =
          planar_chart_basis(math::get_antipode(basis, radius).first.v);
      if (flower_band.x_active) {
        draw_banded_flower(sink, canvas, basis, planar_basis, radius, sides,
                           fragment_shader, shape_phase);
        break;
      }
      draw_sampled(sink, canvas, static_cast<size_t>(sides * 2 + 2),
                   &planar_basis, false, fragment_shader,
                   [&](Fragments &points) {
                     Flower::sample(points, basis, radius, sides, shape_phase);
                   });
      break;
    }
    case ContourShape::PLANAR_STAR: {
      if (dense_contours) {
        draw_dense_planar_star(sink, chords, layout, canvas, basis, radius,
                               sides, fragment_shader, shape_color, shape_phase,
                               contour_index);
        break;
      }
      math::Basis projection_basis;
      const math::Basis &planar_basis =
          *PlanarProjection::edge_basis(basis, radius, projection_basis);
      draw_sampled(sink, canvas, static_cast<size_t>(sides * 2 + 2),
                   &planar_basis, false, fragment_shader,
                   [&](Fragments &points) {
                     Star<PlanarProjection>::sample_positions(
                         points, basis, radius, sides, shape_phase);
                   });
      break;
    }
    case ContourShape::SPHERICAL_STAR:
      draw_sampled(sink, canvas, static_cast<size_t>(sides * 2 + 2), nullptr,
                   dense_contours, fragment_shader, [&](Fragments &points) {
                     Star<GeodesicProjection>::sample_continuous_positions(
                         points, basis, radius, sides, shape_phase);
                   });
      break;
    }
  }

  PlanarBandSplit<W, H> flower_split; ///< Flower petal split against the band.
  ClipBand<W, H> flower_band;         ///< The frame's clip band.
};

} // namespace Plot
