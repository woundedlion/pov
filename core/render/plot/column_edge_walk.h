/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/plot/column_taps.h"

/** @file column_edge_walk.h @brief Column-time edge lattice walking. */

namespace Plot {

HS_O3_BEGIN

/** @brief Walks one edge's planned columns and plots their lattice samples.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam SPLIT Walk the arc as its two sides of its closest approach to a
 * pole, each with its own column tests.
 * @tparam LOCAL Resolve the continuous lattice pitch from each column's arc.
 * @tparam UNIFORM_COLOR Reuse the prepared edge color.
 * @tparam RECORD_WORLD Record samples in source history.
 * @param state The edge's walk.
 * @details A turning arc walks every planned column; otherwise each pass stops
 * where its lit run ends. */
template <int W, int H, RasterConfig Cfg, bool SPLIT, bool LOCAL = false,
          bool UNIFORM_COLOR = false, bool RECORD_WORLD = false>
HS_HOT_INLINE void walk_edge_columns(const ColumnWalk<W> &state) {
  auto &[pipeline, canvas, points, shader, frame, segments, segment_count, band,
         direction, radians, single_time, velocity, source_positions, chained,
         turning, top_y, bottom_y, half_sine, half_cosine] = state;
  constexpr bool LOCAL_LATTICE = LOCAL;
  const auto arm_change = [&](int from) HS_HOT_INLINE {
    return (direction > 0 ? math::fast_wrap(from + 1, W) : from) % (W / 2) == 0;
  };
  size_t sample_count = 0;
  int lattice_count = 0;
  float lattice_step = 1.0f;
  for (int segment = 0; segment < segment_count; ++segment) {
    const detail::ColumnSegment &walk = segments[segment];
    int x = walk.start;
    bool seen = false;
    float carried = -1.0f;
    for (int column = 0; column < walk.count; ++column) {
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::columns_visited>();
      const int owner = x;
      if (radians != 0.0f && band.clipped(owner)) {
        const int skip = direction > 0
                             ? math::fast_wrap(band.rs - owner, W)
                             : math::fast_wrap(owner - band.re + 1, W);
        column += skip - 1;
        x = math::fast_wrap(owner + direction * skip, W);
        carried = -1.0f;
        continue;
      }
      x = math::fast_wrap(owner + direction, W);
      float t = radians == 0.0f && !LOCAL_LATTICE
                    ? 0.0f
                    : Spacetime::SamplingSurface<W>::time(float(owner));
      math::Vector a, b, source_a, source_b;
      GeodesicEdgeSampler sampler;
      float pieces, total;
      bool unit = false;
      if (radians != 0.0f || LOCAL_LATTICE) {
        frame.at(owner, t, a, b, source_positions, source_a, source_b);
        total = frame.bases[0].total_dist +
                (frame.bases[1].total_dist - frame.bases[0].total_dist) * t;
        pieces = frame.pieces[0] + (frame.pieces[1] - frame.pieces[0]) * t;
        if (LOCAL_LATTICE) {
          if constexpr (Cfg.column_edge_lengths) {
            sampler = EdgeWorldsheet::arc(EdgeWorldsheet::normalize(a),
                                          EdgeWorldsheet::normalize(b));
            total = sampler.total_dist;
          } else {
            sampler = EdgeWorldsheet::basis(a, b, total);
          }
          a = sampler.v1;
          unit = true;
        }
      } else {
        sampler = frame.bases[1];
        a = sampler.v1;
        b = frame.ends[1];
        total = sampler.total_dist;
        pieces = frame.pieces[1];
        unit = true;
        source_a = EdgeWorldsheet::at(frame.sources[0], 1.0f);
        source_b = EdgeWorldsheet::at(frame.sources[1], 1.0f);
      }
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::column_path_evaluations>();
      bool lit = radians == 0.0f, finished = false;
      float cut = 1.0f;
      math::Vector turn;
      if (SPLIT) {
        cut = EdgeWorldsheet::split_point(a, b, total, sampler, turn);
        unit = true;
      }
      for (int part = 0; part < (SPLIT ? 2 : 1); ++part) {
        float s_in = 0.0f, s_out = 1.0f;
        bool start_column = true, end_column = true;
        math::Vector part_a = a, part_b = b;
        float from = 0.0f, to = 1.0f, length = total;
        GeodesicEdgeSampler piece;
        bool pieced = false;
        if (SPLIT) {
          from = part ? cut : 0.0f;
          to = part ? 1.0f : cut;
          if (!(to > from))
            break;
          part_a = part ? turn : a;
          part_b = part ? b : turn;
          length = (to - from) * total;
        }
        if (radians != 0.0f || LOCAL_LATTICE) {
          const float inverse_length = length > 1e-7f ? 1.0f / length : 0.0f;
          const int meridian =
              Cfg.column_edge_lengths ? owner % (W / 2) : owner;
          const float arm =
              Cfg.column_edge_lengths && owner >= W / 2 ? -1.0f : 1.0f;
          const float column_cosine =
              arm * math::TrigLUT<W, H>::cos_theta(meridian);
          const float column_sine =
              arm * math::TrigLUT<W, H>::sin_theta[meridian];
          const int arc_direction =
              LOCAL_LATTICE ? (a.x * b.z - a.z * b.x >= 0.0f ? 1 : -1)
                            : direction;
          const float sine = arc_direction * half_sine;
          const math::Vector normal(
              -column_sine * half_cosine - column_cosine * sine, 0.0f,
              column_cosine * half_cosine - column_sine * sine);
          const math::Vector back(
              -column_sine * half_cosine + column_cosine * sine, 0.0f,
              column_cosine * half_cosine + column_sine * sine);
          const float end_side = math::dot(part_b, normal);
          const float back_side = math::dot(part_b, back);
          bool crossing = false;
          if (math::dot(part_a, normal) * end_side <= 0.0f) {
            if (!unit) {
              sampler = EdgeWorldsheet::basis(a, b, total);
              unit = true;
            }
            if (SPLIT) {
              piece = EdgeWorldsheet::basis(part_a, part_b, length);
              pieced = true;
            }
            const GeodesicEdgeSampler basis = SPLIT ? piece : sampler;
            const math::Vector along(
                column_cosine * half_cosine - column_sine * sine, 0.0f,
                column_sine * half_cosine + column_cosine * sine);
            crossing = column_basis_crossing(
                math::dot(basis.v1, normal), end_side,
                math::dot(basis.v_perp, normal), math::dot(basis.v1, along),
                math::dot(basis.v_perp, along), inverse_length, s_out);
          }
          end_column = !crossing && arc_direction * end_side <= 0.0f &&
                       arc_direction * back_side > 0.0f &&
                       part_b.x * column_cosine + part_b.z * column_sine > 0.0f;
          if (!crossing && !end_column)
            continue;
          if (!unit) {
            sampler = EdgeWorldsheet::basis(a, b, total);
            unit = true;
          }
          if (SPLIT && !pieced)
            piece = EdgeWorldsheet::basis(part_a, part_b, length);
          const GeodesicEdgeSampler basis = SPLIT ? piece : sampler;
          if (!LOCAL_LATTICE && carried >= 0.0f) {
            start_column = false;
            s_in = fminf(carried, s_out);
          } else {
            const float start_back = math::dot(basis.v1, back);
            const math::Vector behind(
                column_cosine * half_cosine + column_sine * sine, 0.0f,
                column_sine * half_cosine - column_cosine * sine);
            start_column =
                start_back * back_side > 0.0f ||
                !column_basis_crossing(
                    start_back, back_side, math::dot(basis.v_perp, back),
                    math::dot(basis.v1, behind),
                    math::dot(basis.v_perp, behind), inverse_length, s_in);
            s_in = start_column ? 0.0f : fminf(s_in, s_out);
          }
          carried = crossing && !SPLIT && !arm_change(owner) ? s_out : -1.0f;
          lit = true;
          if (crossing)
            Spacetime::count_geometry<
                &Spacetime::GeometryCounters::column_crossings>();
        }
        bool edge_start = start_column, edge_end = end_column;
        if (SPLIT) {
          s_in = start_column ? from : from + (to - from) * s_in;
          s_out = end_column ? to : from + (to - from) * s_out;
          edge_start = start_column && from == 0.0f;
          edge_end = end_column && to == 1.0f;
        }
        if (LOCAL_LATTICE) {
          const float at = 0.5f * (s_in + s_out);
          float sine, cosine;
          math::poly_sincosf_0_pi(total * at, sine, cosine);
          const float y = sampler.v1.y * cosine + sampler.v_perp.y * sine;
          const float tangent_y =
              sampler.v_perp.y * cosine - sampler.v1.y * sine;
          const float cross_y = sampler.longitude_rate_numerator();
          pieces =
              total / screen_step_components<W, H>(y, tangent_y, cross_y,
                                                   math::RADIANS_PER_COLUMN<W>);
        }
        const int count = LOCAL_LATTICE ? std::max(1, int(ceilf(pieces)))
                                        : EdgeWorldsheet::lattice(pieces);
        if (LOCAL_LATTICE || count != lattice_count) {
          lattice_count = count;
          lattice_step =
              1.0f / (LOCAL_LATTICE ? fmaxf(pieces, 1.0f) : float(count));
        }
        const float lattice_scale =
            LOCAL_LATTICE ? fmaxf(pieces, 1.0f) : float(count);
        const int first = edge_start ? 0 : int(ceilf(s_in * lattice_scale));
        const int last = edge_end ? count - int(chained)
                                  : int(ceilf(s_out * lattice_scale)) - 1;
        const bool plain = single_time || (!edge_start && !edge_end);
        GeodesicEdgeSampler source;
        if (source_positions) {
          if (frame.poses) {
            const math::RotationMatrix pose = frame.poses->column(owner);
            auto unrotate = [&](const math::Vector &v) HS_HOT_INLINE {
              return pose.r0 * v.x + pose.r1 * v.y + pose.r2 * v.z;
            };
            source = {unrotate(sampler.v1), unrotate(sampler.v_perp), total};
          } else {
            source = EdgeWorldsheet::basis(source_a, source_b, total);
          }
        }
        for (int k = first; k <= last; ++k) {
          Spacetime::count_geometry<
              &Spacetime::GeometryCounters::column_generated_samples>();
          if (!LOCAL_LATTICE &&
              ++sample_count > 2 * rasterize_step_budget<W>()) {
            HS_SCAN_METRIC(hs::g_scan_metrics.plot_backstop_hits++);
            break;
          }
          const float s = k >= count ? 1.0f : float(k) * lattice_step;
          float sine = 0.0f, cosine = 1.0f;
          if (k > 0 && k < count)
            math::poly_sincosf_0_pi(total * s, sine, cosine);
          const auto world =
              k == 0 ? sampler.v1
              : k >= count
                  ? EdgeWorldsheet::normalize(b)
                  : newton_unit(sampler.v1 * cosine + sampler.v_perp * sine);
          if (world.y > top_y || world.y < bottom_y)
            continue;
          Fragment f;
          if constexpr (Cfg.interpolate_registers)
            f = Fragment::lerp_registers(points[0], points[1], s);
          Spacetime::count_geometry<
              &Spacetime::GeometryCounters::column_projections>();
          const auto p = math::vector_to_pixel<W, H>(world);
          if (radians == 0.0f && !LOCAL_LATTICE)
            t = Spacetime::SamplingSurface<W>::time(p.x);
          if (plain)
            Spacetime::count_geometry<
                &Spacetime::GeometryCounters::column_single_samples>();
          else
            Spacetime::count_geometry<
                &Spacetime::GeometryCounters::column_tap_samples>();
          if (!source_positions)
            f.pos = world;
          else
            f.pos =
                k == 0 ? source.v1
                : k >= count
                    ? EdgeWorldsheet::normalize(source_b)
                    : newton_unit(source.v1 * cosine + source.v_perp * sine);
          if constexpr (UNIFORM_COLOR) {
            f.color = points[0].color;
          } else {
            shader(f.pos, f, t);
            Spacetime::count_geometry<
                &Spacetime::GeometryCounters::column_shader_calls>();
          }
          Spacetime::count_geometry<
              &Spacetime::GeometryCounters::column_samples>();
          if constexpr (RECORD_WORLD) {
            Spacetime::count_geometry<
                &Spacetime::GeometryCounters::column_history_visits>();
            pipeline.record_world_sample(canvas, f.pos, f.age, t,
                                         LOCAL_LATTICE ? owner : -1, p.x);
          }
          if (LOCAL_LATTICE) {
            pipeline.plot_column_at(canvas, p, owner, f.color.color, f.age,
                                    f.color.alpha, t);
            continue;
          }
          const bool arm_boundary =
              radians == 0.0f &&
              Spacetime::SamplingSurface<W>::near_arm_seam(p.x);
          if (plain && !arm_boundary) {
            pipeline.plot_projected_at(canvas, p, f.color.color, f.age,
                                       f.color.alpha, t);
            continue;
          }
          const ColumnTapSample sample{p, t, s};
          if (arm_boundary)
            splat_seam_taps<W, H, Cfg, UNIFORM_COLOR>(state, sample, f);
          else
            splat_owned_taps<W, H>(state, sample, f);
        }
        finished = finished || (end_column && to == 1.0f);
      }
      if (!lit) {
        carried = -1.0f;
        if (!turning && (seen || column + 1 >= walk.deadline))
          break;
        continue;
      }
      seen = true;
      if (!LOCAL_LATTICE && sample_count > 2 * rasterize_step_budget<W>())
        break;
      if (finished && !turning)
        break;
    }
  }
}

/** @brief walk_edge_columns() for split arcs, out of line.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam RECORD_WORLD Record samples in source history.
 * @param state The edge's walk. */
template <int W, int H, RasterConfig Cfg, bool RECORD_WORLD = false>
HS_HOT_FLASH_MEMBER void walk_split_edge_columns(const ColumnWalk<W> &state) {
  walk_edge_columns<W, H, Cfg, true, false, false, RECORD_WORLD>(state);
}

/** @brief Walks an edge with local column-time cadence.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam RECORD_WORLD Record samples in source history.
 * @param state Prepared column walk. */
template <int W, int H, RasterConfig Cfg, bool RECORD_WORLD = false>
HS_HOT_FLASH_MEMBER void walk_local_edge_columns(const ColumnWalk<W> &state) {
  walk_edge_columns<W, H, Cfg, false, true, false, RECORD_WORLD>(state);
}

/** @brief Walks an edge with its prepared uniform color.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam SPLIT,LOCAL Edge split and lattice policies.
 * @tparam RECORD_WORLD Record samples in source history.
 * @param state Prepared column walk. */
template <int W, int H, RasterConfig Cfg, bool SPLIT, bool LOCAL,
          bool RECORD_WORLD = false>
__attribute__((aligned(32))) HS_HOT_FLASH_MEMBER void
walk_uniform_edge_columns(const ColumnWalk<W> &state) {
  walk_edge_columns<W, H, Cfg, SPLIT, LOCAL, true, RECORD_WORLD>(state);
}

HS_O3_END

} // namespace Plot
