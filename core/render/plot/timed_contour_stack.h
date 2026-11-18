/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file timed_contour_stack.h
 * @brief Concentric twisting contours drawn at their columns' times.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "core/math/3dmath.h"
#include "core/memory.h"
#include "core/platform/attributes.h"
#include "core/render/canvas.h"
#include "core/render/clip.h"
#include "core/render/filter/pipeline.h"
#include "core/render/filter/screen_anti_alias.h"
#include "core/render/filter/screen_direct_aa_sink.h"
#include "core/render/filter/world_orient.h"
#include "core/render/plot/chords.h"
#include "core/render/plot/column_bounds.h"
#include "core/render/plot/contour_layout.h"
#include "core/render/plot/mesh.h"
#include "core/render/plot/shapes.h"
#include "core/render/spacetime/pose_table.h"

namespace Plot {

/**
 * @brief Draws a ContourLayout's stack with every contour turned by its
 *        waveform and the stack turned by the orientation at each column's
 *        time.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details A sine contour turns rigidly about its axis with linearly advancing
 * twist phase; other waveforms hold their frame-end twist. Planar contours
 * stroke through Plot::PlanarChords; geodesic contours are closed chains the
 * column walk draws, in chunks, once their swept vertices can reach the
 * walk's rows and band.
 */
template <int W, int H> class TimedContourStack {
public:
  /** @brief Contour vertices one column walk draws. */
  static constexpr size_t COLUMN_CHUNK_VERTICES = 384;
  /** @brief Contours one column walk draws. */
  static constexpr size_t COLUMN_CHUNK_CONTOURS = 64;
  /** @brief Most anchor intervals a flower or planar polygon edge takes. */
  static constexpr int FLOWER_ANCHOR_INTERVALS =
      PlanarChords<W, H>::MAX_RADIAL_INTERVALS;

  /** @brief Binds the column-time filters to the stack's orientation.
   *  @param orientation Orientation the column-time world stage samples. */
  explicit TimedContourStack(math::Orientation<> &orientation)
      : filters{Filter::World::Orient(orientation),
                Filter::Screen::AntiAlias<W, H>()} {}

  /** @brief Arena bytes init_storage() takes.
   *  @return Byte count. */
  static constexpr size_t storage_bytes() {
    return Filter::World::Orient::storage_bytes(W);
  }

  /** @brief Binds the filters' column pose storage.
   *  @param arena Arena the storage is allocated from. */
  HS_FLASH_MEMBER void init_storage(Arena &arena) {
    filters.init_storage(arena);
  }

private:
  /** @brief A timed contour's twist over the frame and its shading. */
  struct TimedContour {
    math::Vector axis; /**< Axis its twist turns about. */
    float twist;       /**< Signed waveform amplitude, radians. */
    float wave;        /**< Waveform angle at the frame's end, radians. */
    float rate;        /**< Waveform angle change over the frame, radians. */
    float radius;      /**< Distance of its vertices from the axis, at most. */
    Color4 color;      /**< Stroke color and opacity. */
    math::Quaternion start_rotation;  ///< Analytic frame-start turn.
    math::Quaternion middle_rotation; ///< Analytic mid-frame turn.
    Spacetime::PathBound derivatives; ///< Cached trajectory derivative bounds.

    /** @brief Captures the contour's shared turns and derivative bounds. */
    HS_HOT_INLINE void prepare() {
      start_rotation = math::make_rotation(axis, turn_at(0.0f));
      middle_rotation = math::make_rotation(axis, turn_at(0.5f));
      derivatives = path();
    }

    /** @brief Evaluates the shared rigid turn.
     * @param position Frame-end vertex.
     * @param t Frame fraction.
     * @return Vertex at the requested fraction. */
    HS_HOT_INLINE math::Vector at(const math::Vector &position, float t) const {
      if (t == 1.0f)
        return position;
      if (t == 0.0f)
        return math::rotate(position, start_rotation);
      if (t == 0.5f)
        return math::rotate(position, middle_rotation);
      return math::rotate(position, math::make_rotation(axis, turn_at(t)));
    }

    /** @brief Twist from the frame's end back to a frame fraction.
     * @param t Frame fraction.
     * @return The angle to turn a frame-end vertex by, radians. */
    HS_HOT_INLINE float turn_at(float t) const {
      return twist * (sinf(wave - rate * (1.0f - t)) - sinf(wave));
    }

    /** @brief Bounds of a vertex's path over the frame.
     * @return The twist's derivative bounds at the contour's radius. */
    HS_HOT_INLINE Spacetime::PathBound path() const {
      const float amplitude = fabsf(twist) * rate;
      return Spacetime::PathBound::turning(amplitude, amplitude * rate,
                                           amplitude * rate * rate, radius);
    }
  };

  /** @brief A planar contour's edge chart and corner count. */
  struct PlanarChart {
    const math::Basis *basis; ///< Azimuthal-equidistant chart of the edges.
    int corners;              ///< Corners with chart coordinates set.
  };

public:
  /** @brief Peak scratch_arena_a bytes draw() takes, including nested raster
   *  storage.
   *  @return Byte count. */
  static constexpr size_t scratch_a_bytes() {
    return COLUMN_CHUNK_VERTICES *
               (sizeof(math::Vector) + sizeof(uint16_t) + sizeof(Mesh::Edge)) +
           COLUMN_CHUNK_CONTOURS * sizeof(TimedContour) +
           (CONTOUR_MAX_VERTICES + 2) *
               (sizeof(math::Vector) + sizeof(Fragment)) +
           std::max(Mesh::column_scratch_a_bytes(),
                    PlanarChords<W, H>::scratch_a_bytes() +
                        PlanarChords<W, H>::timed_scratch_a_bytes()) +
           128;
  }

  /**
   * @brief Draws the contour stack at its columns' times.
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
    const math::Basis basis =
        math::make_basis(math::Quaternion(), math::X_AXIS);
    const bool sine = frame.sine;
    const float rate = sine ? 2.0f * math::PI_F * frame.twist_step : 0.0f;
    const float global_alpha = frame.opacity;
    const bool continuous_star = shape == ContourShape::SPHERICAL_STAR;
    const auto cull = SweptEdgeCull<W>::template of<H>(canvas);
    const float spin = filters.sweep_arc();
    const Spacetime::PoseTable &poses = *filters.pose_table(canvas);
    const math::RotationMatrix end_pose = poses.at(poses.times());
    const math::Quaternion start_orientation = frame.orientation.sample(0.0f);
    const math::Quaternion middle_orientation = frame.orientation.sample(0.5f);
    const math::Quaternion end_orientation = frame.orientation.get();
    const math::Basis world_basis =
        math::make_basis(end_orientation, math::X_AXIS);
    const bool planar = shape == ContourShape::PLANAR_STAR ||
                        shape == ContourShape::PLANAR_POLYGON ||
                        shape == ContourShape::FLOWER;
    const ClipRegion &clip = canvas.clip();
    CapCenter near_cap{}, far_cap{};
    if (shape == ContourShape::PLANAR_STAR) {
      near_cap = make_cap_center(clip, world_basis.v);
      far_cap = make_cap_center(clip, -world_basis.v);
    }
    if (planar)
      chords.prepare(canvas, fabsf(frame.amplitude) * rate + spin +
                                 2.0f * math::PI_F / W);

    ScratchScope chunk_scope(scratch_arena_a);
    auto *vertices =
        scratch_arena_a.allocate_n<math::Vector>(COLUMN_CHUNK_VERTICES);
    auto *owners = scratch_arena_a.allocate_n<uint16_t>(COLUMN_CHUNK_VERTICES);
    auto *contours =
        scratch_arena_a.allocate_n<TimedContour>(COLUMN_CHUNK_CONTOURS);
    ArenaVector<Mesh::Edge> edges(scratch_arena_a, COLUMN_CHUNK_VERTICES);
    auto *world =
        scratch_arena_a.allocate_n<math::Vector>(CONTOUR_MAX_VERTICES + 2);
    Fragments ring;
    ring.bind(scratch_arena_a, static_cast<size_t>(CONTOUR_MAX_VERTICES + 2));
    size_t count_vertices = 0, count_contours = 0;
    float chunk_motion = 0.0f;

    auto path = [&](size_t i, float t) HS_HOT_FLASH_MEMBER {
      const TimedContour &contour = contours[owners[i]];
      return contour.at(vertices[i], t);
    };
    auto vertex_path = [&](size_t i) HS_HOT_FLASH_MEMBER {
      return contours[owners[i]].derivatives;
    };
    auto shader = [&](const math::Vector &, Fragment &f) HS_HOT_FLASH_MEMBER {
      f.color = contours[owners[edges[static_cast<size_t>(f.v2)].u]].color;
    };
    GeometryCallback geometry;
    geometry.evaluate = path;
    geometry.vertex_path = vertex_path;
    geometry.source_positions = false;
    geometry.chained = true;
    geometry.rigid_lattice = true;
    geometry.uniform_color = true;
    geometry.capture_error = 0.001f;
    auto flush = [&]() HS_HOT_INLINE {
      if (count_vertices == 0)
        return;
      Mesh::draw_column_edges<W, H>(
          sink, canvas, Mesh::VertexSpan{{vertices, count_vertices}}, edges,
          shader, geometry, chunk_motion, &poses, spin);
      count_vertices = 0;
      count_contours = 0;
      chunk_motion = 0.0f;
      edges.clear();
    };

    Color4 pair_color;
    for (int ordinal = 0; ordinal < count; ++ordinal) {
      const ContourSlot slot = layout.slot(ordinal, count, continuous_star,
                                           frame.palette, pair_color);
      const int i = slot.i;
      const float radius_t = slot.radius_t;
      const float radius = slot.radius;
      const ContourShading shading =
          layout.shading(i, radius, continuous_star, frame.opposite,
                         frame.amplitude, pair_color, global_alpha);
      const float direction = shading.direction;
      const float contour_phase = shading.contour_phase;
      const Color4 shaded_color = shading.shaded_color;
      const auto ring_frame_at = ring_frame(basis, radius);
      // A continuous star turns in the shared basis on both sides of the
      // equator; the other primitives are charted about the far pole past it.
      const math::Basis &turn_basis =
          continuous_star ? basis : ring_frame_at.basis;
      TimedContour contour{math::cross(turn_basis.u, turn_basis.w).normalized(),
                           sine ? direction * frame.amplitude : 0.0f,
                           2.0f * math::PI_F * (radius_t + frame.phase),
                           rate,
                           continuous_star ? 1.0f : ring_frame_at.sin_theta,
                           shaded_color,
                           {},
                           {},
                           {}};
      const float motion = fabsf(contour.twist) * rate * contour.radius;
      ring.clear();
      const bool geodesic_flower =
          shape == ContourShape::FLOWER &&
          ring_frame_at.cos_theta > COS_PLANAR_ANTIPODE;
      if (planar && !geodesic_flower) {
        if (shape == ContourShape::PLANAR_STAR) {
          if (!planar_star_cap_visible<W, H>(clip, near_cap, far_cap, radius,
                                             motion + spin))
            continue;
        }
        math::Basis projection_basis;
        const math::Basis *chart = &projection_basis;
        int corners = 2 * sides;
        if (shape == ContourShape::PLANAR_STAR) {
          chart = PlanarProjection::edge_basis(world_basis, radius,
                                               projection_basis);
          Star<PlanarProjection>::sample_chart_positions(
              ring, chords.chart_x(), chords.chart_y(), world_basis, radius,
              sides, contour_phase, layout.radius_trig(i),
              layout.star_step_trig(), *chart);
        } else {
          if (shape == ContourShape::FLOWER) {
            projection_basis = Flower::chart_basis(world_basis, radius);
            Flower::sample(ring, world_basis, radius, sides, contour_phase);
          } else {
            chart = PlanarProjection::edge_basis(world_basis, radius,
                                                 projection_basis);
            Polygon<PlanarProjection>::sample(ring, world_basis, radius, sides,
                                              contour_phase);
            corners = sides;
          }
          for (int k = 0; k < corners; ++k) {
            const auto at = azimuthal_project_exact(ring[k].pos, *chart);
            chords.chart_x()[k] = at.x;
            chords.chart_y()[k] = at.y;
          }
        }
        draw_planar_contour(sink, chords, canvas, ring, {chart, corners},
                            contour, shaded_color, start_orientation,
                            middle_orientation, end_orientation, shape, count,
                            spin, flush);
        continue;
      }
      if (geodesic_flower)
        Flower::sample(ring, basis, radius, sides, contour_phase);
      else if (continuous_star)
        Star<GeodesicProjection>::sample_continuous_positions(
            ring, basis, radius, sides, contour_phase);
      else
        Polygon<GeodesicProjection>::sample(ring, basis, radius, sides,
                                            contour_phase);
      const int points = static_cast<int>(ring.size()) - 1;
      for (int k = 0; k < points; ++k)
        world[k] = end_pose.apply(ring[k].pos);
      if (cull.rejects(world, points, motion + spin))
        continue;
      if (count_vertices + static_cast<size_t>(points) >
              COLUMN_CHUNK_VERTICES ||
          count_contours == COLUMN_CHUNK_CONTOURS)
        flush();
      contour.prepare();
      contours[count_contours] = contour;
      for (int k = 0; k < points; ++k) {
        vertices[count_vertices + k] = ring[k].pos;
        owners[count_vertices + k] = static_cast<uint16_t>(count_contours);
        edges.push_back(
            {static_cast<uint16_t>(count_vertices + k),
             static_cast<uint16_t>(count_vertices + (k + 1) % points)});
      }
      count_vertices += static_cast<size_t>(points);
      ++count_contours;
      chunk_motion = fmaxf(chunk_motion, motion);
    }
    flush();
  }

private:
  /**
   * @brief Strokes a sampled planar contour along its timed turn.
   * @tparam Flush Callable drawing the pending column-walk chunk.
   * @param sink Splat sink.
   * @param chords Planar chord walker holding the contour's chart.
   * @param canvas Target canvas.
   * @param ring Frame-end ring.
   * @param chart Edge chart and corner count.
   * @param contour Contour twist.
   * @param color Stroke color and opacity.
   * @param start_orientation Orientation at the frame's start.
   * @param middle_orientation Orientation at mid-frame.
   * @param end_orientation Orientation at the frame's end.
   * @param shape Planar primitive.
   * @param count Contours drawn.
   * @param spin Column sweep arc.
   * @param flush Called before the stroke so earlier contours draw first.
   */
  template <typename Flush>
  HS_HOT_INLINE void
  draw_planar_contour(Filter::Screen::DirectAntiAliasSink<W, H> &sink,
                      PlanarChords<W, H> &chords, Canvas &canvas,
                      const Fragments &ring, const PlanarChart &chart,
                      const TimedContour &contour, const Color4 &color,
                      const math::Quaternion &start_orientation,
                      const math::Quaternion &middle_orientation,
                      const math::Quaternion &end_orientation,
                      ContourShape shape, int count, float spin, Flush &flush) {
    const math::Quaternion start_turn =
        start_orientation *
        math::make_rotation(contour.axis, contour.turn_at(0.0f)) *
        end_orientation.conjugate();
    const float lever = PlanarChords<W, H>::rotation_radius(
        {chords.chart_x()[0], chords.chart_y()[0]},
        {chords.chart_x()[1], chords.chart_y()[1]});
    const float contour_speed =
        fabsf(contour.twist) * contour.rate * lever + spin;
    flush();
    chords.draw_closed(
        sink, canvas, ring, chart.corners, *chart.basis, color,
        {.start = math::RotationMatrix(start_turn),
         .middle = math::RotationMatrix(
             middle_orientation *
             math::make_rotation(contour.axis, contour.turn_at(0.5f)) *
             end_orientation.conjugate()),
         .max_intervals = shape == ContourShape::PLANAR_STAR
                              ? PlanarChords<W, H>::MAX_ANCHOR_INTERVALS
                              : FLOWER_ANCHOR_INTERVALS,
         .alternating_radii = shape == ContourShape::PLANAR_STAR,
         .reference_chords = shape == ContourShape::PLANAR_STAR &&
                             float(count) >= DENSE_CONTOUR_COUNT,
         .contour_speed = contour_speed});
  }

  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
      filters; ///< Orientation at column times, then anti-aliasing.
};

} // namespace Plot
