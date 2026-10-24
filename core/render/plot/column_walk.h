/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file column_walk.h
 * @brief Geometry callbacks and the column-time edge walk.
 */

/** @brief Borrowed geometric control-point evaluator and path bounds. */
struct GeometryCallback {
  /** Evaluates one original control point at a frame fraction. */
  FunctionRef<math::Vector(size_t, float)> evaluate;
  Spacetime::PathBound path{}; /**< Derivative bounds of every control point's
                                  path over the frame. */

  bool source_positions =
      true; /**< Shaders read source-space positions under world stages;
                 false hands them world positions. */
  /** Bounds each vertex's path, in place of `path`, when set. */
  FunctionRef<Spacetime::PathBound(size_t)> vertex_path{};
  bool chained = false; /**< Edges run along polylines, each vertex starting
                             at most one edge and ending at most one other:
                             the column walk plots each edge's start but not
                             its end, so a closed chain plots every vertex once
                             and an open one all but its last. */

  bool rigid_lattice =
      false; ///< Edges keep angular length; resolve their lattice at column time.
  bool uniform_color =
      false; ///< Shader writes only color, constant over each edge and frame.
  float capture_error =
      0.025f; ///< Maximum projected vertex capture error, pixels.

  /** @brief Constructs an exact evaluator.
   * @param fn Borrowed control-point callback.
   * @return The exact geometry descriptor. */
  static GeometryCallback exact(FunctionRef<math::Vector(size_t, float)> fn) {
    return {fn};
  }
};

HS_O3_BEGIN

/** @brief Angular bound between a meridian crossing's arc fraction and the
 * column walk sampler's parameter, radians. */
inline constexpr float COLUMN_ARC_PARAMETER_PAD = 1e-5f;

/** @brief Largest projected motion between adjacent column times, pixels, at
 * which a column holding an edge end splats its taps at its own time. */
inline constexpr float COLUMN_SINGLE_TAP_ERROR = 0.02f;

/** @brief Intersects a geodesic basis with a forward meridian.
 * @param start_side,end_side Endpoint distances from the meridian plane.
 * @param tangent_side Start tangent distance from the meridian plane.
 * @param start_along,tangent_along Basis projections along the meridian.
 * @param inverse_length Reciprocal arc length, or zero for a collapsed edge.
 * @param fraction Receives the material arc fraction.
 * @return Whether the minor arc crosses the forward meridian. */
HS_HOT_INLINE bool column_basis_crossing(float start_side, float end_side,
                                         float tangent_side, float start_along,
                                         float tangent_along,
                                         float inverse_length,
                                         float &fraction) {
  if (start_side * end_side > 0.0f || fabsf(start_side - end_side) < 1e-9f)
    return false;
  float angle = math::precise_atan2(-start_side, tangent_side);
  const float forward = tangent_side * start_along - start_side * tangent_along;
  if (angle < 0.0f) {
    angle += math::PI_F;
    if (forward > 0.0f)
      return false;
  } else if (forward < 0.0f)
    return false;
  fraction = hs::clamp(angle * inverse_length, 0.0f, 1.0f);
  return true;
}

/** @brief Intersects a minor arc with a destination meridian.
 * @param a,b Arc endpoint directions; a is unit length.
 * @param direction Horizontal meridian direction.
 * @param length Minor arc length.
 * @param tangent Unit tangent at the first endpoint.
 * @param fraction Resulting arc fraction.
 * @return Whether the meridian crosses the minor arc. */
HS_HOT_INLINE bool
column_edge_crossing(const math::Vector &a, const math::Vector &b,
                     const math::Vector &direction, float length,
                     const math::Vector &tangent, float &fraction) {
  const math::Vector normal(-direction.z, 0.0f, direction.x);
  return column_basis_crossing(
      math::dot(a, normal), math::dot(b, normal), math::dot(tangent, normal),
      math::dot(a, direction), math::dot(tangent, direction),
      length > 1e-7f ? 1.0f / length : 0.0f, fraction);
}

/** @brief A mesh edge as two vertex indices. */
struct MeshEdge {
  uint16_t u; ///< First endpoint vertex index.
  uint16_t v; ///< Second endpoint vertex index.
};

/** @brief Mesh edges sharing one column walk. */
struct EdgeBatch {
  size_t count; /**< Edge count. */
  void (*prepare)(Fragments &, const math::Vector &, const math::Vector &,
                  size_t,
                  float); /**< Writes an edge's endpoint registers into the
                                walk's fragments from its frame-end
                                vertices and index. */
  const math::Vector *vertices; /**< Frame-end vertex positions. */
  const Spacetime::Worldlines<>
      &motion;                /**< Source captures and interpolants of the
                              vertices this frame. */
  const MeshEdge *pairs;      /**< Vertex indices of each edge. */
  bool chained = false;       /**< Edges leave their end vertex unplotted. */
  bool rigid_lattice = false; ///< Resolve lattice and AA support per column.
  bool uniform_color = false; ///< Shader writes only one color per edge.
};

/** @brief Rejects edges whose swept arcs miss a draw's rows and columns.
 * @tparam W Display width.
 * @details Rows use the column walk's own sample cut; columns use the render
 * band widened by two columns on each side. */
template <int W> struct SweptEdgeCull {
  float top;           /**< Highest unit-sphere y the walk emits. */
  float bottom;        /**< Lowest unit-sphere y the walk emits. */
  math::Vector first;  /**< First longitude half-plane normal. */
  math::Vector second; /**< Second longitude half-plane normal. */
  int mode = 0; /**< 0: rows only; 1: hidden outside either half-plane; 2:
                   hidden inside both. */

  /** @brief Prepares the bounds for one draw.
   * @param columns Render-band column predicate.
   * @param top_y,bottom_y Unit-sphere y limits of emitted samples. */
  HS_HOT_FLASH_MEMBER SweptEdgeCull(const ClipRegion::XClip &columns,
                                    float top_y, float bottom_y)
      : top(top_y), bottom(bottom_y) {
    if (!columns.active)
      return;
    constexpr float MARGIN = 2.0f;
    const float start =
        (float(columns.rs) - 0.5f - MARGIN) * math::RADIANS_PER_COLUMN<W>;
    const float width = (float(columns.length(W)) + 2.0f * MARGIN) *
                        math::RADIANS_PER_COLUMN<W>;
    if (width >= 2.0f * math::PI_F)
      return;
    const float end = start + width;
    if (width <= math::PI_F) {
      mode = 1;
      first = {-sinf(start), 0.0f, cosf(start)};
      second = {sinf(end), 0.0f, -cosf(end)};
    } else {
      mode = 2;
      first = {-sinf(end), 0.0f, cosf(end)};
      second = {sinf(start), 0.0f, -cosf(start)};
    }
  }

  /** @brief Highest unit-sphere y a canvas's column walk emits.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows bound the walk.
   * @return The height of the row two above the clip. */
  template <int H> HS_HOT_INLINE static float top_of(const Canvas &canvas) {
    return cosf(math::DisplayGeometry<H>::row_to_phi(
        float(std::max(0, canvas.clip().y_start - 2))));
  }

  /** @brief Lowest unit-sphere y a canvas's column walk emits.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows bound the walk.
   * @return The height of the row two below the clip. */
  template <int H> HS_HOT_INLINE static float bottom_of(const Canvas &canvas) {
    return cosf(math::DisplayGeometry<H>::row_to_phi(
        float(std::min(H, canvas.clip().y_end + 2))));
  }

  /** @brief The bounds of a canvas's column walk.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows and render band bound the walk.
   * @return The cull. */
  template <int H> static SweptEdgeCull of(const Canvas &canvas) {
    return SweptEdgeCull(canvas.render_x_clip(), top_of<H>(canvas),
                         bottom_of<H>(canvas));
  }

  /** @brief Tests whether every geodesic between points a frame sweeps
   * stays outside the bounds.
   * @param points Positions whose normalized convex hull, padded, holds every
   * drawn sample.
   * @param count Number of points.
   * @param pad Bound on a drawn sample's distance from that hull.
   * @return True when no drawn sample can be emitted. */
  HS_HOT_FLASH_MEMBER bool rejects(const math::Vector *points, int count,
                                   float pad) const {
    return rejects_hull(points, count, pad);
  }

  /** @brief Tests whether every arc an edge sweeps over the frame stays
   * outside the bounds.
   * @param points Both endpoints at the frame's start and end.
   * @param pad Bound on a drawn endpoint's distance from its chord.
   * @return True when no drawn sample of the edge can be emitted. */
  HS_HOT_FLASH_MEMBER bool rejects(const math::Vector (&points)[4],
                                   float pad) const {
    return rejects_hull(points, 4, pad);
  }

private:
  /** @brief rejects() for any count, inlined so a fixed count unrolls.
   * @param points Positions whose padded hull holds every drawn sample.
   * @param count Number of points.
   * @param pad Bound on a drawn sample's distance from that hull.
   * @return True when no drawn sample can be emitted. */
  HS_HOT_INLINE bool rejects_hull(const math::Vector *points, int count,
                                  float pad) const {
    float low, high;
    if (Spacetime::swept_heights(points, count, pad, low, high) &&
        (high < bottom || low > top))
      return true;
    if (mode == 1)
      return Spacetime::swept_beyond(points, count, pad, -first) ||
             Spacetime::swept_beyond(points, count, pad, -second);
    return mode == 2 && Spacetime::swept_beyond(points, count, pad, first) &&
           Spacetime::swept_beyond(points, count, pad, second);
  }
};

namespace detail {
/** @brief One arm pass's walk over the columns an edge can light. */
struct ColumnSegment {
  int start = 0;    /**< First column in walk order. */
  int count = 0;    /**< Columns to walk. */
  int deadline = 0; /**< Columns walked unlit before the walk gives up. */
};

} // namespace detail

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
    const auto span = make_geodesic_edge_span(a, b);
    return {a, span.have_axis ? math::cross(span.axis, a) : math::Vector{},
            span.total};
  }

  /** @brief Lattice intervals the walk spaces along an edge.
   * @param pieces Edge length in screen steps.
   * @return The length rounded, at least one. */
  HS_HOT_INLINE static int lattice(float pieces) {
    return std::max(1, int(pieces + 0.5f));
  }

  /** @brief The lattice an edge keeps over the whole frame.
   * @tparam W,H Display dimensions.
   * @param starts Unit endpoints at the frame's start.
   * @param ends Unit endpoints at the frame's end.
   * @param most Largest count of interest.
   * @param end_to_start Rigid rotation from frame end to start, when known.
   * @param arcs Receive the edge's arcs at the frame's end and, for a count
   * of at most @p most, its start.
   * @return The walk's lattice intervals when the edge's length rounds to
   * the same count, at most @p most, at both ends of the frame, and so at
   * every time between; zero otherwise. */
  template <int W, int H>
  HS_HOT_FLASH_MEMBER static int
  steady_lattice(const math::Vector (&starts)[2], const math::Vector (&ends)[2],
                 int most, GeodesicEdgeSampler (&arcs)[2],
                 const math::RotationMatrix *end_to_start = nullptr) {
    arcs[1] = arc(ends[0], ends[1]);
    const int count = lattice(pieces_of<W, H>(arcs[1]));
    if (count > most)
      return 0;
    arcs[0] = end_to_start
                  ? GeodesicEdgeSampler{starts[0],
                                        end_to_start->apply(arcs[1].v_perp),
                                        arcs[1].total_dist}
                  : arc(starts[0], starts[1]);
    return lattice(pieces_of<W, H>(arcs[0])) == count ? count : 0;
  }

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
  HS_HOT_FLASH_MEMBER void
  prepare(const Spacetime::Worldlines<> &motion, uint16_t from, uint16_t to,
          const Spacetime::PoseTable *table, float spin, bool source,
          bool measure_pitch = true) {
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
    for (int k = 0; k < 2; ++k) {
      const math::Vector &a = world[k][0], &b = world[k][1];
      const auto span = make_geodesic_edge_span(a, b);
      bases[k] = {a,
                  span.have_axis ? math::cross(span.axis, a) : math::Vector{},
                  span.total};
      ends[k] = b;
      pieces[k] = measure_pitch ? pieces_of<W, H>(bases[k]) : 0.0f;
    }
    constexpr float CAPTURE_ROUNDING = 3e-5f;
    deviation = sqrtf(bend) + CAPTURE_ROUNDING;
    if (poses)
      deviation +=
          0.125f * spin * (spin + 2.0f * (sqrtf(travel) + 4.0f * deviation));
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

  /** @brief Last possible lit column of a pass under bounded endpoint tracks.
   * @tparam W Display width.
   * @param pass,direction Arm pass and walk direction.
   * @param span,widening Initial endpoint separation and its change, columns.
   * @param last Last guarded endpoint crossing, in walk order.
   * @param start First planned column in walk order.
   * @return Exclusive end in pass coordinates. */
  template <int W>
  HS_HOT_FLASH_MEMBER int planned_end(int pass, int direction, float span,
                                      float widening, float last,
                                      int start) const {
    int end = W / 2;
    constexpr int PASS = W / 2;
    constexpr float HALF = W / 2.0f;
    const int final_column = pass * PASS + (direction > 0 ? PASS - 1 : 0);
    const float final_time = direction > 0 ? (HALF - 1.0f) / HALF : 0.0f;
    const float a = tracks[0].origin + tracks[0].travel * final_time;
    const float width = span + widening * final_time;
    float into =
        width < 0.0f ? a - float(final_column) : float(final_column) - a;
    into -= float(W) * floorf(into * (1.0f / W));
    if (into > fabsf(width) + tracks[1].guard / tracks[1].rate &&
        into < float(W) - tracks[0].guard / tracks[0].rate) {
      if (!std::isfinite(last))
        return start;
      end = std::min(PASS, int(ceilf(last)) + 1);
    }
    return end;
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
    const bool start_known = Spacetime::model_track<W>(
        bases[0].v1, bases[1].v1, deviation, 1.0f, tracks[0]);
    const bool end_known =
        Spacetime::model_track<W>(ends[0], ends[1], deviation, 1.0f, tracks[1]);
    known = start_known && end_known;
    float span = 0.0f, widening = 0.0f;
    if (known) {
      span = math::centered_remainder(tracks[1].origin - tracks[0].origin,
                                      float(W));
      widening = tracks[1].travel - tracks[0].travel;
      known = fmaxf(fabsf(span), fabsf(span + widening)) <
              HALF - tracks[0].guard - tracks[1].guard;
    }
    float entry[2] = {std::numeric_limits<float>::infinity(),
                      std::numeric_limits<float>::infinity()};
    float last[2] = {-entry[0], -entry[0]};
    if (known)
      for (const auto &track : tracks)
        Spacetime::track_crossings<W>(
            track, 0.0f, 1.0f, true, true,
            [&](const Spacetime::Crossing &pass) HS_HOT_INLINE {
              const int arm = pass.arm & 1;
              const float at = pass.time * HALF;
              const float walk = direction > 0 ? at : HALF - 1.0f - at;
              entry[arm] = fminf(entry[arm], walk - track.guard);
              last[arm] = fmaxf(last[arm], walk + track.guard);
            },
            track.guard);
    int count = 0;
    for (int pass = 0; pass < 2; ++pass) {
      const int first = pass * PASS + (direction > 0 ? 0 : PASS - 1);
      if (!known) {
        segments[count++] = {first, PASS, PASS};
        continue;
      }
      int start = 0, deadline = PASS;
      if (entry[pass] > 0.0f) {
        const float time = direction > 0 ? 0.0f : (HALF - 1.0f) / HALF;
        const float a = tracks[0].origin + tracks[0].travel * time;
        const float width = span + widening * time;
        float into = width < 0.0f ? a - float(first) : float(first) - a;
        into -= float(W) * floorf(into * (1.0f / W));
        if (into <= fabsf(width)) {
          deadline = 1;
        } else {
          if (!(entry[pass] < HALF))
            continue;
          start = int(entry[pass]);
          deadline = int(ceilf(last[pass])) + 1 - start;
        }
      } else if (last[pass] < HALF) {
        deadline = int(ceilf(last[pass])) + 1;
      }
      const int end = bounded ? planned_end<W>(pass, direction, span, widening,
                                               last[pass], start)
                              : PASS;
      segments[count++] = {math::fast_wrap(first + direction * start, W),
                           std::max(0, end - start), deadline};
    }
    return count;
  }
};

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
  TimedShaderFn shader;                  /**< Time-aware source shader. */
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
  bool turning;    /**< The arc can turn its longitude direction this frame. */
  float top_y;     /**< Highest emitted unit-sphere y. */
  float bottom_y;  /**< Lowest emitted unit-sphere y. */
  float half_sine; /**< Sine of half a column. */
  float half_cosine; /**< Cosine of half a column. */
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
      column, color, age, alpha, sample.t);
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
                                             const TimedShaderFn &shader) {
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

/** @brief Walks one edge's planned columns and plots their lattice samples.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam SPLIT Walk the arc as its two sides of its closest approach to a
 * pole, each with its own column tests.
 * @tparam RIGID Resolve the continuous lattice pitch from each column's arc.
 * @tparam UNIFORM_COLOR Reuse the prepared edge color.
 * @param state The edge's walk.
 * @details A turning arc walks every planned column; otherwise each pass stops
 * where its lit run ends. */
template <int W, int H, RasterConfig Cfg, bool SPLIT, bool RIGID = false,
          bool UNIFORM_COLOR = false>
HS_HOT_INLINE void walk_edge_columns(const ColumnWalk<W> &state) {
  auto &[pipeline, canvas, points, shader, frame, segments, segment_count, band,
         direction, radians, single_time, velocity, source_positions, chained,
         turning, top_y, bottom_y, half_sine, half_cosine] = state;
  constexpr bool rigid_lattice = RIGID;
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
      float t = radians == 0.0f && !rigid_lattice
                    ? 0.0f
                    : Spacetime::SamplingSurface<W>::time(float(owner));
      math::Vector a, b, source_a, source_b;
      GeodesicEdgeSampler sampler;
      float pieces, total;
      bool unit = false;
      if (radians != 0.0f || rigid_lattice) {
        frame.at(owner, t, a, b, source_positions, source_a, source_b);
        total = frame.bases[0].total_dist +
                (frame.bases[1].total_dist - frame.bases[0].total_dist) * t;
        pieces = frame.pieces[0] + (frame.pieces[1] - frame.pieces[0]) * t;
        if (rigid_lattice) {
          sampler = EdgeWorldsheet::basis(a, b, total);
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
      size_t emitted = 0;
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
        if (radians != 0.0f || rigid_lattice) {
          const float inverse_length = length > 1e-7f ? 1.0f / length : 0.0f;
          const float column_cosine = math::TrigLUT<W, H>::cos_theta(owner);
          const float column_sine = math::TrigLUT<W, H>::sin_theta[owner];
          const int arc_direction =
              rigid_lattice ? (a.x * b.z - a.z * b.x >= 0.0f ? 1 : -1)
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
          if (!rigid_lattice && carried >= 0.0f) {
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
        if (rigid_lattice) {
          const float at = 0.5f * (s_in + s_out);
          float sine, cosine;
          math::poly_sincosf_0_pi(total * at, sine, cosine);
          const float y = sampler.v1.y * cosine + sampler.v_perp.y * sine;
          const float tangent_y =
              sampler.v_perp.y * cosine - sampler.v1.y * sine;
          const float cross_y =
              sampler.v1.x * sampler.v_perp.z - sampler.v1.z * sampler.v_perp.x;
          pieces =
              total / screen_step_components<W, H>(y, tangent_y, cross_y,
                                                   math::RADIANS_PER_COLUMN<W>);
        }
        const int count = rigid_lattice ? std::max(1, int(ceilf(pieces)))
                                        : EdgeWorldsheet::lattice(pieces);
        if (rigid_lattice || count != lattice_count) {
          lattice_count = count;
          lattice_step =
              1.0f / (rigid_lattice ? fmaxf(pieces, 1.0f) : float(count));
        }
        const float lattice_scale =
            rigid_lattice ? fmaxf(pieces, 1.0f) : float(count);
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
          if (!rigid_lattice && ++sample_count > 2 * rasterize_step_budget<W>())
            break;
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
          if (radians == 0.0f && !rigid_lattice)
            t = Spacetime::SamplingSurface<W>::time(p.x);
          ++emitted;
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
          if (rigid_lattice) {
            pipeline.plot_column_at(canvas, p, owner, f.color.color, f.age,
                                    f.color.alpha, t);
            continue;
          }
          const bool arm_boundary =
              radians == 0.0f &&
              fabsf(p.x - (W / 2.0f) * floorf(p.x * (2.0f / W) + 0.5f)) <= 1.0f;
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
      if (emitted > 1)
        Spacetime::count_geometry<
            &Spacetime::GeometryCounters::column_span_fills>();
      if (!rigid_lattice && sample_count > 2 * rasterize_step_budget<W>())
        break;
      if (finished && !turning)
        break;
    }
  }
}

/** @brief walk_edge_columns() for split arcs, out of line.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @param state The edge's walk. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void walk_split_edge_columns(const ColumnWalk<W> &state) {
  walk_edge_columns<W, H, Cfg, true>(state);
}

/** @brief Walks a rigid edge with local column-time cadence.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @param state Prepared column walk. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void walk_rigid_edge_columns(const ColumnWalk<W> &state) {
  walk_edge_columns<W, H, Cfg, false, true>(state);
}

/** @brief Walks an edge with its prepared uniform color.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam SPLIT,RIGID Edge split and lattice policies.
 * @param state Prepared column walk. */
template <int W, int H, RasterConfig Cfg, bool SPLIT, bool RIGID>
__attribute__((aligned(32))) HS_HOT_FLASH_MEMBER void
walk_uniform_edge_columns(const ColumnWalk<W> &state) {
  walk_edge_columns<W, H, Cfg, SPLIT, RIGID, true>(state);
}

/** @brief Rasterizes mesh edges one destination column at a time.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
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
 * @details Each arm pass walks, in the edge's travel order, the render-band
 * columns its plan bounds, each with the edge's endpoints evaluated on their
 * paths at the column's time and rotated by its pose, and stops where the lit
 * run ends. Samples sit at arc fractions k/N, where N is the edge's length in
 * screen steps at that time, and a column plots the lattice points between
 * its two boundary meridians unweighted through the filter chain, so a moving
 * edge samples each column as the frozen edge at that column's time does. The
 * start vertex is lattice point zero; each edge also plots its end vertex
 * unless the batch is chained.
 * Columns strictly inside an edge splat at their own time; the end columns of
 * an edge moving faster than COLUMN_SINGLE_TAP_ERROR move each neighbor tap
 * to that column's time, and taps astride an arm seam interpolate their
 * registers toward the sample shaded at the far side's time. A stationary draw walks the lattice once. An
 * edge whose plan cannot bound it, such as one that can pass near a pole, is
 * walked as the two sides of its closest approach to the pole, out of line. */
template <int W, int H, RasterConfig Cfg>
FASTRUN HS_O3_FN HS_NOINLINE_NOCLONE void
rasterize_column_edges(PipelineRef &pipeline, Canvas &canvas, Fragments &points,
                       TimedShaderFn shader, const EdgeBatch &batch,
                       const Spacetime::PoseTable *poses, float radians,
                       bool source_positions, float spin) {
  const float support = batch.rigid_lattice ? 1.0f : 0.5f;
  const float half_sine = sinf(support * math::RADIANS_PER_COLUMN<W>);
  const float half_cosine = cosf(support * math::RADIANS_PER_COLUMN<W>);

  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  const float top_y = SweptEdgeCull<W>::template top_of<H>(canvas);
  const float bottom_y = SweptEdgeCull<W>::template bottom_of<H>(canvas);
  Spacetime::count_geometry<&Spacetime::GeometryCounters::column_draw_setups>();
  const SweptEdgeCull<W> cull(canvas.render_x_clip(), top_y, bottom_y);
  const auto &band = canvas.render_x_clip();
  ScratchScope draw_scope(scratch_arena_a);
  source_positions = source_positions && poses;
  auto &frame = *scratch_arena_a.make<EdgeWorldsheet>();
  for (size_t edge = 0; edge < batch.count; ++edge) {
    const MeshEdge &pair = batch.pairs[edge];
    frame.prepare<W, H>(batch.motion, pair.u, pair.v, poses, spin,
                        source_positions, !batch.rigid_lattice);
    if (cull.rejects({frame.bases[0].v1, frame.bases[1].v1, frame.ends[0],
                      frame.ends[1]},
                     frame.deviation + 1e-4f)) {
      Spacetime::count_geometry<
          &Spacetime::GeometryCounters::column_edges_culled>();
      continue;
    }
    const int direction =
        frame.bases[1].v1.x * frame.bases[1].v_perp.z -
                    frame.bases[1].v1.z * frame.bases[1].v_perp.x >=
                0.0f
            ? 1
            : -1;
    detail::ColumnSegment segments[2];
    int segment_count = 1;
    float velocity[4] = {};
    bool single_time = true;
    if (radians == 0.0f && !batch.rigid_lattice) {
      segments[0] = {0, 1, 1};
    } else {
      segment_count = frame.plan<W>(direction, segments, batch.rigid_lattice);
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
      if (!batch.rigid_lattice) {
        frame.velocities<W, H>(velocity);
        single_time = fmaxf(fmaxf(fabsf(velocity[0]), fabsf(velocity[1])),
                            fmaxf(fabsf(velocity[2]), fabsf(velocity[3]))) <=
                      COLUMN_SINGLE_TAP_ERROR * (W / 2);
      }
    }
    Spacetime::count_geometry<
        &Spacetime::GeometryCounters::column_edges_prepared>();
    batch.prepare(points, batch.vertices[pair.u], batch.vertices[pair.v], edge,
                  frame.bases[1].total_dist);
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
      return basis.v1.x * basis.v_perp.z - basis.v1.z * basis.v_perp.x;
    };
    const bool split = radians != 0.0f && !frame.known;
    const ColumnWalk<W> walk{
        pipeline,
        canvas,
        points,
        shader,
        frame,
        segments,
        segment_count,
        band,
        direction,
        radians,
        single_time,
        velocity,
        source_positions,
        batch.chained,
        batch.rigid_lattice ||
            (split && !(sense(frame.bases[0]) * sense(frame.bases[1]) >
                        4.0f * frame.deviation * frame.deviation)),
        top_y,
        bottom_y,
        half_sine,
        half_cosine};
    if (batch.uniform_color) {
      if (batch.rigid_lattice)
        walk_uniform_edge_columns<W, H, Cfg, false, true>(walk);
      else if (split)
        walk_uniform_edge_columns<W, H, Cfg, true, false>(walk);
      else
        walk_uniform_edge_columns<W, H, Cfg, false, false>(walk);
    } else if (batch.rigid_lattice)
      walk_rigid_edge_columns<W, H, Cfg>(walk);
    else if (split)
      walk_split_edge_columns<W, H, Cfg>(walk);
    else
      walk_edge_columns<W, H, Cfg, false>(walk);
  }
}

HS_O3_END
