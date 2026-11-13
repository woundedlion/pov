/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/plot/column_bounds.h"

/** @file edge_worldsheet.h @brief Moving edge paths, world bases and column plans. */

namespace Plot {

HS_O3_BEGIN

/** @brief One moving edge's endpoint paths and frame-end world bases. */
struct EdgeWorldsheet {
  math::Vector curves[2][3];  /**< Each endpoint's drawn path as
                                c0 + c1 t + c2 t^2: its source path, or its
                                world path under poses fixed over the frame. */
  math::Vector sources[2][3]; /**< Each endpoint's source path, when shaders
                                 read source positions. */
  const Spacetime::PoseTable *poses =
      nullptr; /**< Poses by column time, or null when fixed. */
  GeodesicEdgeSampler bases[2]; /**< World bases at the frame's ends. */
  math::Vector ends[2];         /**< World end vertices at the frame's ends. */
  float pieces[2] = {}; /**< Edge length in screen steps at the frame's ends. */
  float deviation = 0.0f; /**< Bound on a drawn endpoint's distance from the
                             chord between its frame-end world positions. */
  Spacetime::LongitudeTrack tracks[2]; /**< Longitude models of both ends. */
  bool known = false; /**< plan() bounded both tracks and the arc's side. */

  /** @brief Normalizes a near-unit vector.
   * @param value Nonzero vector.
   * @return Unit vector with a refined reciprocal-square-root estimate. */
  HS_HOT_INLINE static math::Vector normalize(const math::Vector &value) {
    const float norm = math::dot(value, value);
    const float estimate = screen_rsqrt(norm);
    return value * (estimate * (1.5f - 0.5f * norm * estimate * estimate));
  }

  /** @brief Builds the orthonormal basis of the arc between two directions.
   * @param a,b Nonzero endpoint directions.
   * @param total Arc length, radians.
   * @return Unit start, unit tangent toward @p b, and @p total; a collapsed
   * edge keeps a zero tangent. */
  HS_HOT_INLINE static GeodesicEdgeSampler
  basis(const math::Vector &a, const math::Vector &b, float total) {
    const math::Vector v1 = normalize(a);
    const math::Vector tangent = b - v1 * math::dot(v1, b);
    const float norm = math::dot(tangent, tangent);
    return {v1, norm > 1e-12f ? normalize(tangent) : math::Vector{}, total};
  }

  /** @brief Finds where an arc comes closest to a pole.
   * @param a,b Endpoint directions.
   * @param total Arc length, radians.
   * @param basis Receives the arc's orthonormal basis.
   * @param turn Receives the arc's point at the returned fraction.
   * @return The arc fraction of its height's extreme, or one when the extreme
   * lies outside the arc.
   * @details Longitude along a great circle is monotone and turns by half a
   * turn between its highest and lowest points, so the arc's two sides of
   * the extreme each span less than half the longitudes. The height has an
   * extreme inside the arc exactly when its slopes at the two ends,
   * b.y - a.y (a.b) and -(a.y - b.y (a.b)), have opposite signs. */
  HS_HOT_FLASH_MEMBER static float
  split_point(const math::Vector &a, const math::Vector &b, float total,
              GeodesicEdgeSampler &basis, math::Vector &turn) {
    basis = EdgeWorldsheet::basis(a, b, total);
    const math::Vector end = normalize(b);
    const float cosine_ab = math::dot(basis.v1, end);
    float cut = 1.0f;
    if ((end.y - basis.v1.y * cosine_ab) * (basis.v1.y - end.y * cosine_ab) >
        0.0f) {
      float apex = atan2f(basis.v_perp.y, basis.v1.y);
      if (apex < 0.0f)
        apex += math::PI_F;
      if (apex > 0.0f && apex < total)
        cut = apex / total;
    }
    float sine, cosine;
    math::poly_sincosf_0_pi(cut * total, sine, cosine);
    turn = basis.v1 * cosine + basis.v_perp * sine;
    return cut;
  }

  /** @brief The geodesic basis of an arc between unit directions.
   * @param a,b Unit endpoints.
   * @return Start, unit tangent toward @p b and arc length; a collapsed arc
   * keeps a zero tangent. */
  HS_HOT_INLINE static GeodesicEdgeSampler arc(const math::Vector &a,
                                               const math::Vector &b) {
    const math::Vector chord = a - b;
    const float chord2 = math::dot(chord, chord);
    const float total =
        chord2 <= 0.0625f
            ? sqrtf(chord2) *
                  (1.0f + chord2 * (1.0f / 24.0f +
                                    chord2 * (3.0f / 640.0f +
                                              chord2 * (5.0f / 7168.0f))))
            : unit_arc_length(a, b);
    const auto span = make_geodesic_edge_span(a, b, total);
    return {a, span.have_axis ? math::cross(span.axis, a) : math::Vector{},
            span.total};
  }

  /** @brief Lattice intervals the walk spaces along an edge.
   * @param pieces Edge length in screen steps.
   * @return The length rounded, at least one. */
  HS_HOT_INLINE static int lattice(float pieces) {
    return std::max(1, int(pieces + 0.5f));
  }

  /** @brief Selects one lattice count for both ends of a frame.
   * @param start,end Endpoint lengths in screen steps.
   * @return The larger rounded count, at least one. */
  HS_HOT_INLINE static int frame_lattice(float start, float end) {
    return lattice(fmaxf(start, end));
  }

  /** @brief The lattice an edge keeps over the whole frame.
   * @tparam W,H Display dimensions.
   * @param starts Unit endpoints at the frame's start.
   * @param ends Unit endpoints at the frame's end.
   * @param most Largest count of interest.
   * @param end_to_start Rigid rotation from frame end to start, when known.
   * @param frame_constant Allows differing endpoint counts, using their maximum.
   * @param arcs Receive the edge's arcs at the frame's end and, for a count
   * of at most @p most, its start.
   * @return The admitted count at most @p most, or zero. */
  template <int W, int H>
  HS_HOT_FLASH_MEMBER static int
  steady_lattice(const math::Vector (&starts)[2], const math::Vector (&ends)[2],
                 int most, GeodesicEdgeSampler (&arcs)[2],
                 const math::RotationMatrix *end_to_start = nullptr,
                 bool frame_constant = false) {
    return steady_lattice_count<W, H>(starts, ends, most, arcs, end_to_start,
                                      frame_constant);
  }

  /** @brief Finds a steady lattice through the ITCM entry.
   * @tparam W,H Display dimensions.
   * @param starts,ends Unit endpoints at the frame's start and end.
   * @param most Largest count of interest.
   * @param arcs Receive the edge's arcs at the frame's end and start.
   * @param end_to_start Rigid rotation from frame end to start, when known.
   * @param frame_constant Allows differing endpoint counts, using their maximum.
   * @return The admitted count at most @p most, or zero. */
  template <int W, int H>
  HS_HOT_ITCM_MEMBER static int
  steady_lattice_resident(const math::Vector (&starts)[2],
                          const math::Vector (&ends)[2], int most,
                          GeodesicEdgeSampler (&arcs)[2],
                          const math::RotationMatrix *end_to_start = nullptr,
                          bool frame_constant = false) {
    return steady_lattice_count<W, H>(starts, ends, most, arcs, end_to_start,
                                      frame_constant);
  }

private:
  /** @brief Admits matching or frame-constant endpoint lattice counts.
   * @tparam W,H Display dimensions.
   * @param starts,ends Unit endpoints at the frame's start and end.
   * @param most Largest count of interest.
   * @param arcs Receive the edge's arcs at the frame's end and start.
   * @param end_to_start Rigid rotation from frame end to start, when known.
   * @param frame_constant Allows differing endpoint counts, using their maximum.
   * @return The admitted count at most @p most, or zero. */
  template <int W, int H>
  HS_HOT_INLINE static int steady_lattice_count(
      const math::Vector (&starts)[2], const math::Vector (&ends)[2], int most,
      GeodesicEdgeSampler (&arcs)[2], const math::RotationMatrix *end_to_start,
      bool frame_constant) {
    arcs[1] = arc(ends[0], ends[1]);
    const int count = lattice(pieces_of<W, H>(arcs[1]));
    if (count > most)
      return 0;
    arcs[0] = end_to_start
                  ? GeodesicEdgeSampler{starts[0],
                                        end_to_start->apply(arcs[1].v_perp),
                                        arcs[1].total_dist}
                  : arc(starts[0], starts[1]);
    const float first = pieces_of<W, H>(arcs[0]);
    if (!frame_constant)
      return lattice(first) == count ? count : 0;
    const int maximum = frame_lattice(first, float(count));
    return maximum <= most ? maximum : 0;
  }

public:
  /** @brief A lattice point of an arc, as the walk places it.
   * @param basis The arc.
   * @param s Arc fraction in (0, 1).
   * @return The unit point at @p s. */
  HS_HOT_INLINE static math::Vector
  lattice_point(const GeodesicEdgeSampler &basis, float s) {
    float sine, cosine;
    math::poly_sincosf_0_pi(basis.total_dist * s, sine, cosine);
    return newton_unit(basis.v1 * cosine + basis.v_perp * sine);
  }

  /** @brief Measures an edge's length in screen steps.
   * @tparam W,H Display dimensions.
   * @param basis Orthonormal geodesic basis.
   * @return Arc length over the screen step at the arc's middle. */
  template <int W, int H>
  HS_HOT_INLINE static float pieces_of(const GeodesicEdgeSampler &basis) {
    const auto middle = basis.accurate(0.5f);
    return basis.total_dist / screen_step<W, H>(middle.pos, middle.tan,
                                                math::RADIANS_PER_COLUMN<W>);
  }

  /** @brief Whether midpoint spacing exceeds the bilinear sample footprint.
   * @tparam W,H Display dimensions.
   * @tparam CERTIFIED Use a sufficient bound before the full spacing checks.
   * @param basis Orthonormal geodesic basis.
   * @param pieces Edge length in midpoint screen steps.
   * @return Whether the arc needs column-local sampling. */
  template <int W, int H, bool CERTIFIED = true>
  HS_HOT_INLINE static bool
  needs_local_lattice(const GeodesicEdgeSampler &basis, float pieces) {
#if HS_ENABLE_TEST_HOOKS
    if constexpr (CERTIFIED)
      ++g_column_lattice_checks;
#endif
    const float base_step = math::RADIANS_PER_COLUMN<W>;
    constexpr float SUPPORT = 2.0f;
    constexpr float KX = W / (2.0f * math::PI_F);
    const float ky = math::ROWS_PER_RADIAN<H>;
    const auto end = basis.accurate(1.0f);
    float peak_y2 = fmaxf(basis.v1.y * basis.v1.y, end.pos.y * end.pos.y);
    const bool turning = basis.v_perp.y * end.tan.y < 0.0f;
    if (turning)
      peak_y2 = fminf(1.0f, basis.v1.y * basis.v1.y +
                                basis.v_perp.y * basis.v_perp.y);
    const float count = float(lattice(pieces));
    const float sin2 = fmaxf(1e-7f, 1.0f - peak_y2);
    const float max_rate2 = fmaxf(KX * KX / sin2, ky * ky);
    if (basis.total_dist * basis.total_dist * max_rate2 <=
        SUPPORT * SUPPORT * count * count)
      return false;
    const float cross_y = basis.longitude_rate_numerator();
    if constexpr (CERTIFIED) {
      constexpr float ROUNDING_PAD = 8e-6f;
      constexpr float RSQRT_MIN = 0.9988f;
      const float safe_sin2 = fmaxf(1e-7f, sin2 - ROUNDING_PAD);
      const float row_slack =
          fmaxf(0.0f, basis.v1.y * basis.v1.y +
                          basis.v_perp.y * basis.v_perp.y - 1.0f) +
          ROUNDING_PAD;
      const float safe_cross = fabsf(cross_y) + ROUNDING_PAD;
      const float sin4 = safe_sin2 * safe_sin2;
      const float upper_speed = fmaxf(
          KX * KX * safe_cross * safe_cross +
              ky * ky * (safe_sin2 + row_slack) * safe_sin2,
          (SCREEN_STEP_PX * SCREEN_STEP_PX / (base_step * base_step)) * sin4);
      if (basis.total_dist * basis.total_dist * upper_speed <=
          SUPPORT * SUPPORT * count * count * sin4 * RSQRT_MIN * RSQRT_MIN) {
#if HS_ENABLE_TEST_HOOKS
        ++g_column_lattice_certificates;
#endif
        return false;
      }
    }
    float step = fminf(screen_step<W, H>(basis.v1, basis.v_perp, base_step),
                       screen_step<W, H>(end.pos, end.tan, base_step));
    if (turning) {
      const float peak_y = sqrtf(peak_y2);
      step = fminf(
          step, screen_step_components<W, H>(peak_y, 0.0f, cross_y, base_step));
    }
    if (basis.v1.y * end.pos.y < 0.0f) {
      const float tangent_y = sqrtf(fmaxf(0.0f, 1.0f - cross_y * cross_y));
      step = fminf(step, screen_step_components<W, H>(0.0f, tangent_y, cross_y,
                                                      base_step));
    }
    return basis.total_dist * SCREEN_STEP_PX > SUPPORT * count * step;
  }

  /** @brief Tests the edge's swept hull against a row interval.
   * @tparam W Display width.
   * @param top,bottom Highest and lowest emitted unit-sphere heights.
   * @return Whether the padded hull can reach the rows. */
  template <int W>
  HS_HOT_FLASH_MEMBER bool reaches_rows(float top, float bottom) const {
    const SweptEdgeCull<W> cull({}, top, bottom);
    return !cull.rejects({bases[0].v1, bases[1].v1, ends[0], ends[1]},
                         deviation + 1e-4f);
  }

  /** @brief Evaluates a path held as polynomial coefficients.
   * @param curve Coefficients c0, c1, c2.
   * @param t Frame fraction.
   * @return c0 + c1 t + c2 t^2. */
  HS_HOT_INLINE static math::Vector at(const math::Vector (&curve)[3],
                                       float t) {
    return curve[0] + (curve[1] + curve[2] * t) * t;
  }

  /** @brief Prepares one edge's paths, frame-end bases and deviation bound.
   * @tparam W,H Display dimensions.
   * @param motion Vertex captures and interpolants.
   * @param from,to Start and end vertex.
   * @param table The world stages' poses by column time, or null without
   * world stages.
   * @param spin Rotation bound of the poses over the frame, radians.
   * @param source Whether to keep the source paths.
   * @param measure_pitch Whether to measure the frame-end sampling pitch.
   * @param cull Optional swept clip gate before geodesic setup.
   * @return Whether the edge reaches the gate and its bases were prepared.
   * @details Each path is the parabola through its start, middle and end:
   * the middle capture for Spacetime::Interpolant::QUADRATIC, the arc's middle
   * for Spacetime::Interpolant::ARC, and the chord's for
   * Spacetime::Interpolant::CHORD, which makes it the chord. It lies within
   * |middle - (start + end) / 2| of its chord. Under poses fixed over the
   * frame the paths are rotated once. Moving poses turning at most spin
   * radians per frame about a point moving along a path that departs from its
   * chord by d move a drawn endpoint off its world chord by d plus an eighth
   * of the second-derivative bound spin (spin + 2 travel), travel bounding
   * the path's speed by its chord plus 4 d; captures add their rounding. */
  template <int W, int H>
  HS_HOT_FLASH_MEMBER bool
  prepare(const Spacetime::Worldlines<> &motion, uint16_t from, uint16_t to,
          const Spacetime::PoseTable *table, float spin, bool source,
          bool measure_pitch = true, const SweptEdgeCull<W> *cull = nullptr) {
    poses = table && !table->stationary() ? table : nullptr;
    const bool fixed = table && !poses;
    const math::RotationMatrix pose =
        fixed ? table->at(table->times()) : math::RotationMatrix();
    const math::RotationMatrix start_pose =
        poses ? poses->at(0) : math::RotationMatrix();
    const math::RotationMatrix end_pose =
        poses ? poses->at(poses->times()) : math::RotationMatrix();
    float bend = 0.0f, travel = 0.0f;
    math::Vector world[2][2];
    for (int k = 0; k < 2; ++k) {
      const uint16_t vertex = k ? to : from;
      math::Vector start = motion.start(vertex), end = motion.end(vertex);
      const math::Vector chord_middle = (start + end) * 0.5f;
      const auto kind = motion.interpolant(vertex);
      math::Vector middle = kind == Spacetime::Interpolant::CHORD ? chord_middle
                            : kind == Spacetime::Interpolant::QUADRATIC
                                ? motion.middle(vertex)
                                : motion.worldline(vertex).middle;
      const math::Vector offset = middle - chord_middle;
      bend = fmaxf(bend, math::dot(offset, offset));
      travel = fmaxf(travel, math::distance_squared(start, end));
      if (source) {
        sources[k][0] = start;
        sources[k][1] = middle * 4.0f - start * 3.0f - end;
        sources[k][2] = (start + end - middle * 2.0f) * 2.0f;
      }
      if (fixed) {
        start = pose.apply(start);
        middle = pose.apply(middle);
        end = pose.apply(end);
      }
      curves[k][0] = start;
      curves[k][1] = middle * 4.0f - start * 3.0f - end;
      curves[k][2] = (start + end - middle * 2.0f) * 2.0f;
      if (poses) {
        start = start_pose.apply(start);
        end = end_pose.apply(end);
      }
      world[0][k] = normalize(start);
      world[1][k] = normalize(end);
    }
    constexpr float CAPTURE_ROUNDING = 3e-5f;
    deviation = sqrtf(bend) + CAPTURE_ROUNDING;
    if (poses)
      deviation +=
          0.125f * spin * (spin + 2.0f * (sqrtf(travel) + 4.0f * deviation));
    if (cull &&
        cull->rejects({world[0][0], world[1][0], world[0][1], world[1][1]},
                      deviation + 1e-4f))
      return false;
    for (int k = 0; k < 2; ++k) {
      const math::Vector &a = world[k][0], &b = world[k][1];
      bases[k] = arc(a, b);
      ends[k] = b;
      pieces[k] = measure_pitch ? pieces_of<W, H>(bases[k]) : 0.0f;
    }
    return true;
  }

  /** @brief Evaluates both endpoints at a column's time.
   * @param owner Column.
   * @param t The column's frame fraction.
   * @param a,b Receive the drawn world endpoints, unnormalized.
   * @param source Whether to evaluate the source endpoints too.
   * @param source_a,source_b Receive the source endpoints when asked. */
  HS_HOT_INLINE void at(int owner, float t, math::Vector &a, math::Vector &b,
                        bool source, math::Vector &source_a,
                        math::Vector &source_b) const {
    a = at(curves[0], t);
    b = at(curves[1], t);
    if (source) {
      source_a = at(sources[0], t);
      source_b = at(sources[1], t);
    }
    if (poses) {
      const math::RotationMatrix pose = poses->column(owner);
      a = pose.apply(a);
      b = pose.apply(b);
    }
  }

  /** @brief Projected endpoint motion over the frame.
   * @tparam W,H Display dimensions.
   * @param motion Receives start x, start y, end x and end y velocities in
   * pixels per frame. */
  template <int W, int H> HS_HOT_INLINE void velocities(float motion[4]) const {
    motion[0] = tracks[0].travel;
    motion[1] =
        y_to_screen_row<H>(bases[1].v1.y) - y_to_screen_row<H>(bases[0].v1.y);
    motion[2] = tracks[1].travel;
    motion[3] = y_to_screen_row<H>(ends[1].y) - y_to_screen_row<H>(ends[0].y);
  }

  /** @brief Locates a column relative to the endpoint tracks' longitude arc.
   * @tparam W Display width.
   * @param column Absolute canvas column.
   * @param time Frame fraction.
   * @param initial_arc_span Signed frame-start endpoint separation, columns.
   * @param arc_span_delta Change in signed endpoint separation over the frame, columns.
   * @return Wrapped directed offset and absolute arc extent, columns. */
  template <int W>
  HS_HOT_INLINE detail::ColumnArcPosition
  arc_column_position(int column, float time, float initial_arc_span,
                      float arc_span_delta) const {
    const float arc_start = tracks[0].origin + tracks[0].travel * time;
    const float signed_extent = initial_arc_span + arc_span_delta * time;
    float offset = signed_extent < 0.0f ? arc_start - float(column)
                                        : float(column) - arc_start;
    offset -= float(W) * floorf(offset * (1.0f / W));
    return {offset, fabsf(signed_extent)};
  }

  /** @brief Bounds both endpoint tracks and their arc's longitude side.
   * @tparam W Display width.
   * @param initial_arc_span Receives signed frame-start endpoint separation, columns.
   * @param arc_span_delta Receives its change over the frame, columns.
   * @return Whether both tracks are bounded and their arc stays on the minor side. */
  template <int W>
  HS_HOT_INLINE bool prepare_longitude_tracks(float &initial_arc_span,
                                              float &arc_span_delta) {
    constexpr float HALF = W / 2.0f;
    const bool start_known = Spacetime::model_track<W>(
        bases[0].v1, bases[1].v1, deviation, 1.0f, tracks[0]);
    const bool end_known =
        Spacetime::model_track<W>(ends[0], ends[1], deviation, 1.0f, tracks[1]);
    if (!(start_known && end_known))
      return false;
    initial_arc_span =
        math::centered_remainder(tracks[1].origin - tracks[0].origin, float(W));
    arc_span_delta = tracks[1].travel - tracks[0].travel;
    return fmaxf(fabsf(initial_arc_span),
                 fabsf(initial_arc_span + arc_span_delta)) <
           HALF - tracks[0].guard - tracks[1].guard;
  }

  /** @brief Last possible lit column of a pass under bounded endpoint tracks.
   * @tparam W Display width.
   * @param pass,direction Arm pass and walk direction.
   * @param initial_arc_span,arc_span_delta Signed frame-start endpoint separation
   * and its change over the frame, columns.
   * @param last_crossing Last guarded endpoint crossing, in walk columns.
   * @param first_step First planned step relative to the pass head.
   * @return Exclusive end step relative to the pass head. */
  template <int W>
  HS_HOT_FLASH_MEMBER int
  planned_end(int pass, int direction, float initial_arc_span,
              float arc_span_delta, float last_crossing, int first_step) const {
    int end_step = W / 2;
    constexpr int PASS = W / 2;
    constexpr float HALF = W / 2.0f;
    const int final_column = pass * PASS + (direction > 0 ? PASS - 1 : 0);
    const float final_time = direction > 0 ? (HALF - 1.0f) / HALF : 0.0f;
    const auto position = arc_column_position<W>(
        final_column, final_time, initial_arc_span, arc_span_delta);
    const float start_margin = tracks[0].guard / tracks[0].rate;
    const float end_margin = tracks[1].guard / tracks[1].rate;
    if (position.offset > position.extent + end_margin &&
        position.offset < float(W) - start_margin) {
      if (!std::isfinite(last_crossing))
        return first_step;
      end_step = std::min(PASS, int(ceilf(last_crossing)) + 1);
    }
    return end_step;
  }

  /** @brief Bounds the columns each arm pass can light.
   * @tparam W Display width.
   * @param direction Walk direction.
   * @param segments Receives up to one walk per arm pass.
   * @param bounded Trim trailing columns certified outside AA support.
   * @return Number of walks.
   * @details Each column of an arm pass is drawn at its own time, so under the
   * linear track model of Spacetime::model_track() the sweep crosses each
   * endpoint's longitude at most once per pass, at a column known within the
   * track's guard. A geodesic edge is at most a half circle and so covers the
   * shorter side between its endpoints' longitudes. A pass's lit columns
   * begin at its first column when that column lies on the arc's side and
   * otherwise at the earliest guarded crossing, and an unlit walk gives up
   * after the latest one; a pass with neither lights nothing. Tracks the
   * model cannot bound, and endpoints whose arc can change sides over the
   * frame, walk the whole pass. */
  template <int W>
  HS_O3_FN HS_NOINLINE_NOCLONE int
  plan(int direction, detail::ColumnSegment segments[2], bool bounded = false) {
    constexpr int PASS = W / 2;
    constexpr float HALF = W / 2.0f;
    float initial_arc_span = 0.0f, arc_span_delta = 0.0f;
    known = prepare_longitude_tracks<W>(initial_arc_span, arc_span_delta);
    float first_crossing[2] = {std::numeric_limits<float>::infinity(),
                               std::numeric_limits<float>::infinity()};
    float last_crossing[2] = {-first_crossing[0], -first_crossing[0]};
    if (known)
      for (const auto &track : tracks)
        Spacetime::track_crossings<W>(
            track, 0.0f, 1.0f, true, true,
            [&](const Spacetime::Crossing &pass) HS_HOT_INLINE {
              const int arm = pass.arm & 1;
              const float sweep_column = pass.time * HALF;
              const float walk_step =
                  direction > 0 ? sweep_column : HALF - 1.0f - sweep_column;
              first_crossing[arm] =
                  fminf(first_crossing[arm], walk_step - track.guard);
              last_crossing[arm] =
                  fmaxf(last_crossing[arm], walk_step + track.guard);
            });
    int count = 0;
    for (int pass = 0; pass < 2; ++pass) {
      const int first_column = pass * PASS + (direction > 0 ? 0 : PASS - 1);
      if (!known) {
        segments[count++] = {first_column, PASS, PASS};
        continue;
      }
      int first_step = 0, deadline = PASS;
      if (first_crossing[pass] > 0.0f) {
        const float time = direction > 0 ? 0.0f : (HALF - 1.0f) / HALF;
        const auto position = arc_column_position<W>(
            first_column, time, initial_arc_span, arc_span_delta);
        if (position.offset <= position.extent) {
          deadline = 1;
        } else {
          if (!(first_crossing[pass] < HALF))
            continue;
          first_step = int(first_crossing[pass]);
          deadline = int(ceilf(last_crossing[pass])) + 1 - first_step;
        }
      } else if (last_crossing[pass] < HALF) {
        deadline = int(ceilf(last_crossing[pass])) + 1;
      }
      const int end_step =
          bounded
              ? planned_end<W>(pass, direction, initial_arc_span,
                               arc_span_delta, last_crossing[pass], first_step)
              : PASS;
      segments[count++] = {
          math::fast_wrap(first_column + direction * first_step, W),
          std::max(0, end_step - first_step), deadline};
    }
    return count;
  }
};

HS_O3_END

} // namespace Plot
