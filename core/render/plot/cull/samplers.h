/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/plot/cull.h.

// --- Strategy Helpers ---

/**
 * @brief Forward azimuthal-equidistant projection of a sphere point to the plane.
 * @param p Unit sphere point to project.
 * @param basis Projection basis; center is basis.v, axes basis.u/basis.w.
 * @return Plane coordinates whose radius is the great-circle angle from the
 *         basis center (radians) and whose azimuth follows the basis u/w axes.
 * @details azimuthal_unproject is the inverse map.
 */
static inline std::pair<float, float>
azimuthal_project(const math::Vector &p, const math::Basis &basis) {
  float R = math::angle_between(p, basis.v);
  if (R < math::EPS_GEOMETRIC)
    return {0.0f, 0.0f};
  float theta = math::fast_atan2(math::dot(p, basis.w), math::dot(p, basis.u));
  return {R * math::fast_cosf(theta), R * math::fast_sinf(theta)};
}

/**
 * @brief Inverse azimuthal-equidistant projection from the plane to the sphere.
 * @param Px Plane x-coordinate (azimuthal-equidistant).
 * @param Py Plane y-coordinate (azimuthal-equidistant).
 * @param basis Projection basis; center is basis.v, axes basis.u/basis.w.
 * @return Near-unit point from fast trig at chart radius sqrt(Px²+Py²).
 * @details Up to approximation error, the radius is the great-circle angle
 * only through PI; larger radii fold around the sphere.
 */
static inline math::Vector azimuthal_unproject(float Px, float Py,
                                               const math::Basis &basis) {
  HS_PLOT_COUNT(planar_unprojects);
  float R = sqrtf(Px * Px + Py * Py);
  if (R < math::EPS_GEOMETRIC)
    return basis.v;
  float sin_r;
  float cos_r;
  if (R <= math::PI_F) {
    math::fast_sincosf_0_pi(R, sin_r, cos_r);
  } else {
    sin_r = math::fast_sinf(R);
    cos_r = math::fast_cosf(R);
  }
  const math::Vector radial = (basis.u * Px) + (basis.w * Py);
  return (basis.v * cos_r) + (radial * (sin_r / R));
}

/**
 * @brief A rasterized sample: its sphere position and arc-length tangent.
 * @details `tan` estimates the curve's unit tangent with respect to ARC LENGTH
 * at the sample; zero for a degenerate edge.
 */
struct SamplePT {
  math::Vector pos;
  math::Vector tan;
};

/**
 * @brief Renormalizes a sampled position with one Newton step.
 * @param v A sampler position, unit up to the fast sin/cos residual.
 * @return v scaled to unit length to second order.
 * @details The fast sin/cos kernels leave a sampled position up to 2e-3 off
 * unit, which phi = acos(v.y) turns into a near-pole row offset. One Newton
 * step leaves 5e-6.
 */
static inline math::Vector newton_unit(const math::Vector &v) {
  const float norm2 = math::dot(v, v);
  return v * (1.5f - 0.5f * norm2);
}

constexpr int PLANAR_LEN_SAMPLES = 4;

/**
 * @brief Cumulative on-sphere length estimate at PLANAR_LEN_SAMPLES+1
 *        evenly-spaced PROJECTION samples of the azimuthal-equidistant straight edge whose
 *        projection starts at `proj` and spans (dx, dy).
 * @details arc_cumul[0] = 0; arc_cumul.back() is the PLANAR_LEN_SAMPLES-chord sum,
 * an arc-length estimate. Trig and angle approximations prevent a guaranteed
 * one-sided error bound, even on radial edges.
 */
static inline void
planar_arc_cumul(const std::pair<float, float> &proj, float dx, float dy,
                 const math::Basis &planar_basis,
                 std::array<float, PLANAR_LEN_SAMPLES + 1> &arc_cumul) {
  HS_PLOT_ADD(planar_arc_samples, PLANAR_LEN_SAMPLES + 1);
  arc_cumul[0] = 0.0f;
  math::Vector prev =
      azimuthal_unproject(proj.first, proj.second, planar_basis);
  for (int k = 1; k <= PLANAR_LEN_SAMPLES; ++k) {
    float p = static_cast<float>(k) / PLANAR_LEN_SAMPLES;
    math::Vector cur = azimuthal_unproject(proj.first + dx * p,
                                           proj.second + dy * p, planar_basis);
    arc_cumul[k] = arc_cumul[k - 1] + math::angle_between(prev, cur);
    prev = cur;
  }
}

/**
 * @brief Arc-uniform sampler for one azimuthal-equidistant straight edge.
 * @details `pos(t)` returns the position alone; `operator()(t)` adds the
 * tangent.
 */
struct PlanarEdgeSampler {
  std::pair<float, float> proj1; /**< Projection of the edge start. */
  float dx;                      /**< Projected chord x-component. */
  float dy;                      /**< Projected chord y-component. */
  const math::Basis *basis;      /**< Azimuthal-equidistant projection basis. */
  math::Vector chart_tangent;    /**< Constant chart-space edge tangent. */
  /** Cumulative on-sphere length estimate at evenly-spaced PROJECTION samples. */
  std::array<float, PLANAR_LEN_SAMPLES + 1> arc_cumul;
  float dist; /**< Sampled on-sphere edge length (radians). */

  /** @brief Maps normalized arc distance to normalized chart distance. */
  float projection_fraction(float s) const {
    int interval = 0;
    return projection_fraction_monotonic(s, interval);
  }

  /** @brief Maps increasing arc fractions without rescanning prior intervals. */
  float projection_fraction_monotonic(float s, int &interval) const {
    if (dist < math::EPS_GEOMETRIC)
      return s;
    const float target = s * dist;
    while (interval < PLANAR_LEN_SAMPLES - 1 &&
           arc_cumul[interval + 1] < target)
      ++interval;
    const float seg = arc_cumul[interval + 1] - arc_cumul[interval];
    const float frac =
        seg > math::EPS_GEOMETRIC ? (target - arc_cumul[interval]) / seg : 0.0f;
    return hs::clamp((static_cast<float>(interval) + frac) / PLANAR_LEN_SAMPLES,
                     0.0f, 1.0f);
  }

  /**
   * @brief Unprojects the chart line at PROJECTION fraction p in [0,1].
   * @details Projection-uniform, so not arc-uniform under the anisotropic
   * metric; pos() maps an arc fraction onto it.
   */
  math::Vector unproject(float p) const {
    return azimuthal_unproject(proj1.first + dx * p, proj1.second + dy * p,
                               *basis);
  }

  /**
   * @brief Unprojects the chart line at PROJECTION fraction p, carrying the
   *        analytic tangent when `WithTangent`.
   * @tparam WithTangent Also derive and normalize the tangent; the position-only
   *         instantiation discards the rate terms before they are computed.
   */
  template <bool WithTangent>
  __attribute__((always_inline)) SamplePT sample_at(float p) const {
    const float x = proj1.first + dx * p;
    const float y = proj1.second + dy * p;
    const float r2 = x * x + y * y;
    if (r2 < math::EPS_GEOMETRIC * math::EPS_GEOMETRIC) {
      if constexpr (!WithTangent)
        return {basis->v, math::Vector()};
      HS_PLOT_COUNT(normalizations);
      return {basis->v, math::normalized_or(chart_tangent, math::Vector())};
    }

    const float radius = sqrtf(r2);
    const float inv_radius = 1.0f / radius;
    float sin_radius;
    float cos_radius;
    if (radius <= math::PI_F) {
      math::fast_sincosf_0_pi(radius, sin_radius, cos_radius);
    } else {
      sin_radius = math::fast_sinf(radius);
      cos_radius = math::fast_cosf(radius);
    }
    const math::Vector radial = (basis->u * x) + (basis->w * y);
    const float radial_scale = sin_radius * inv_radius;
    const math::Vector position =
        (basis->v * cos_radius) + (radial * radial_scale);
    if constexpr (!WithTangent)
      return {position, math::Vector()};

    const float radius_rate = (x * dx + y * dy) * inv_radius;
    const float scale_rate =
        (radius * cos_radius - sin_radius) * inv_radius * inv_radius;
    const math::Vector tangent = (basis->v * (-sin_radius * radius_rate)) +
                                 (chart_tangent * radial_scale) +
                                 (radial * (scale_rate * radius_rate));
    HS_PLOT_COUNT(normalizations);
    return {position, math::normalized_or(tangent, math::Vector())};
  }

  /**
   * @brief Position at arc fraction s in [0,1].
   * @details Inverts the piecewise-linear cumulative-arc table to a projection
   * parameter, then unprojects.
   */
  math::Vector pos(float s) const { return unproject(projection_fraction(s)); }

  /** @brief Evaluates position and analytic tangent without a second unproject. */
  HS_FLASH_MEMBER SamplePT one_pass(float s) const {
    return sample_at<true>(projection_fraction(s));
  }

  /** @brief Evaluates an increasing sequence without rescanning arc intervals. */
  SamplePT one_pass_monotonic(float s, int &interval) const {
    return sample_at<true>(projection_fraction_monotonic(s, interval));
  }

  /** @brief Evaluates only position for an increasing sample sequence. */
  math::Vector position_monotonic(float s, int &interval) const {
    return sample_at<false>(projection_fraction_monotonic(s, interval)).pos;
  }

  /** @brief Position and unit tangent at arc fraction s in [0,1]. */
  SamplePT operator()(float s) const { return one_pass(s); }
};

/** @brief Builds the reusable arc sampler for one planar edge. */
static inline PlanarEdgeSampler
make_planar_edge_sampler(const math::Vector &a, const math::Vector &b,
                         const math::Basis &planar_basis) {
  PlanarEdgeSampler sampler;
  sampler.proj1 = azimuthal_project(a, planar_basis);
  auto proj2 = azimuthal_project(b, planar_basis);
  sampler.dx = proj2.first - sampler.proj1.first;
  sampler.dy = proj2.second - sampler.proj1.second;
  sampler.basis = &planar_basis;
  sampler.chart_tangent =
      (planar_basis.u * sampler.dx) + (planar_basis.w * sampler.dy);
  planar_arc_cumul(sampler.proj1, sampler.dx, sampler.dy, planar_basis,
                   sampler.arc_cumul);
  sampler.dist = sampler.arc_cumul[PLANAR_LEN_SAMPLES];
  return sampler;
}

/**
 * @brief Planar interpolation strategy: builds an arc-uniform sampler for one edge.
 * @tparam ProcessSegmentFn Callable (sample, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param planar_basis Azimuthal-equidistant projection basis.
 * @param is_last_segment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length sampler, endpoints, on-sphere
 *                        length (radians), and the last-segment flag.
 * @details The path is a straight line in the azimuthal-equidistant projection.
 * Projection-uniform stepping is not arc-uniform, so a short cumulative-arc
 * table maps an arc-length fraction to a projection parameter.
 */
template <typename ProcessSegmentFn>
static void
rasterize_planar_strategy(const Fragment &curr, const Fragment &next,
                          const math::Basis &planar_basis, bool is_last_segment,
                          ProcessSegmentFn &&process_segment) {
  PlanarEdgeSampler sampler =
      make_planar_edge_sampler(curr.pos, next.pos, planar_basis);

  process_segment(sampler, curr, next, sampler.dist, is_last_segment);
}

/**
 * @brief Chord-sum estimate (radians) of the on-sphere length of the
 *        azimuthal-equidistant straight edge a->b, the path planar
 *        interpolation actually renders.
 * @details Sums the same planar_arc_cumul lengths the planar sampler walks.
 * Trig and angle approximations prevent a guaranteed one-sided error bound,
 * even on radial edges.
 */
static inline float planar_arc_length(const math::Vector &a,
                                      const math::Vector &b,
                                      const math::Basis &planar_basis) {
  auto p1 = azimuthal_project(a, planar_basis);
  auto p2 = azimuthal_project(b, planar_basis);
  std::array<float, PLANAR_LEN_SAMPLES + 1> arc_cumul;
  planar_arc_cumul(p1, p2.first - p1.first, p2.second - p1.second, planar_basis,
                   arc_cumul);
  return arc_cumul[PLANAR_LEN_SAMPLES];
}

/**
 * @brief Unit axis perpendicular to v, stable for an antipodal geodesic.
 * @param v Unit endpoint of a near-antipodal great-circle segment.
 * @return A unit axis perpendicular to v.
 * @details Antipodal endpoints leave cross(v1, v2) ~= 0, so the axis is taken
 * from the world axis least parallel to v.
 */
static inline math::Vector stable_perpendicular_axis(const math::Vector &v) {
  HS_PLOT_COUNT(normalizations);
  return math::perpendicular_axis(v);
}

/**
 * @brief Azimuthal-equidistant projection chart centered on a pole.
 * @param center Unit pole the planar chart is centered on (the 'v' axis).
 * @return A Basis {u, center, w} with u, w spanning the chart plane.
 */
static inline math::Basis planar_chart_basis(const math::Vector &center) {
  HS_PLOT_ADD(normalizations, 2);
  math::Vector u = math::perpendicular_axis(center);
  math::Vector w = math::cross(center, u).normalized();
  return {u, center, w};
}

/** @brief Constant sampler for a coincident-endpoint edge. */
struct DegenerateEdgeSampler {
  math::Vector p; /**< The collapsed edge's single position. */

  math::Vector pos(float) const { return p; }
  SamplePT operator()(float) const { return {p, math::Vector()}; }
};

/**
 * @brief Great-circle sampler for one edge.
 * @details v1 and v_perp are orthonormal, so pos and tan are near-unit
 * combinations of the same approximate sin/cos — the
 * screen-velocity sampler's tangent costs no extra trig.
 */
struct GeodesicEdgeSampler {
  /** @brief pos() is unit up to one newton_unit() correction. */
  static constexpr bool NEWTON_UNIT = true;

  math::Vector v1;     /**< Edge start (unit). */
  math::Vector v_perp; /**< Unit vector perpendicular to v1 in the arc plane. */
  float total_dist;    /**< The edge's on-sphere length (radians). */

  /** @brief Position at arc fraction t in [0,1]. */
  math::Vector pos(float t) const {
    float s, c;
    math::fast_sincosf_0_pi(total_dist * t, s, c);
    return (v1 * c) + (v_perp * s);
  }

  /** @brief Near-unit position and tangent at arc fraction t in [0,1]. */
  SamplePT operator()(float t) const {
    float s, c;
    math::fast_sincosf_0_pi(total_dist * t, s, c);
    return {(v1 * c) + (v_perp * s), (v_perp * c) - (v1 * s)};
  }
};

/**
 * @brief Shared per-edge geodesic setup: arc length and slerp axis.
 * @details have_axis is false on an edge the renderer collapses to a dot, which
 *          has no arc to pole.
 */
struct GeodesicEdgeSpan {
  float
      total; /**< Arc-length estimate from unit_arc_length(a, b), in radians. */
  bool antipodal; /**< axis came from stable_perpendicular_axis, not cross. */
  bool have_axis; /**< axis holds a unit arc pole. */
  bool
      azimuth_bounded; /**< |axis.y| >= AXIS_Y_EPS: the walked axis resolves the longitude sweep direction. */
  math::Vector axis;   /**< Unit arc pole (valid iff have_axis). */
};

#if HS_ENABLE_TEST_ORACLES
/** @brief Number of geodesic edge spans built since the last test reset. */
inline uint32_t g_geodesic_edge_span_builds = 0;
#endif

/**
 * @brief Computes the shared geodesic edge setup once per edge.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 */
static __attribute__((always_inline)) inline GeodesicEdgeSpan
make_geodesic_edge_span(const math::Vector &a, const math::Vector &b) {
#if HS_ENABLE_TEST_ORACLES
  ++g_geodesic_edge_span_builds;
#endif
  GeodesicEdgeSpan es;
  es.total = unit_arc_length(a, b);
  if (es.total < EPS_GEODESIC_SEGMENT) {
    es.antipodal = false;
    es.axis = math::Vector(0.0f, 0.0f, 0.0f);
    es.have_axis = false;
    es.azimuth_bounded = false;
    return es;
  }
  math::Vector pole = math::cross(a, b);
  float pole_len_sq = math::dot(pole, pole);
  es.antipodal = pole_len_sq < EPS_ARC_POLE_SQ;
  if (es.antipodal) {
    es.axis = stable_perpendicular_axis(a);
    es.azimuth_bounded = std::abs(es.axis.y) >= AXIS_Y_EPS;
  } else {
    HS_PLOT_COUNT(normalizations);
    es.axis = pole * (1.0f / sqrtf(pole_len_sq));
    es.azimuth_bounded = std::abs(es.axis.y) >= AXIS_Y_EPS;
  }
  es.have_axis = true;
  return es;
}

/**
 * @brief Geodesic interpolation strategy: builds a great-circle sampler for one edge.
 * @tparam ProcessSegmentFn Callable (sample, curr, next, dist, isLast) -> void.
 * @param curr Start fragment of the edge.
 * @param next End fragment of the edge.
 * @param is_last_segment True if this is the final edge of the polyline.
 * @param process_segment Receives the arc-length sampler, endpoints, on-sphere
 *                        length (radians), and the last-segment flag.
 * @details Slerps about the axis make_geodesic_edge_span resolves; a
 * coincident-endpoint edge collapses to a constant sampler.
 */
HS_O3_BEGIN
template <typename ProcessSegmentFn>
static void rasterize_geodesic_strategy(const Fragment &curr,
                                        const Fragment &next,
                                        bool is_last_segment,
                                        ProcessSegmentFn &&process_segment) {
  HS_PLOT_STALL_START(edge_setup_start);
  math::Vector v1 = curr.pos;
  const GeodesicEdgeSpan es = make_geodesic_edge_span(v1, next.pos);

  if (!es.have_axis) {
    HS_PLOT_COUNT(degenerate);
    HS_PLOT_STALL_STOP(edge_setup, edge_setup_start);
    process_segment(DegenerateEdgeSampler{v1}, curr, next, es.total,
                    is_last_segment);
  } else {
    const GeodesicEdgeSampler sampler{v1, math::cross(es.axis, v1), es.total};
    HS_PLOT_STALL_STOP(edge_setup, edge_setup_start);
    process_segment(sampler, curr, next, es.total, is_last_segment);
  }
}
HS_O3_END

constexpr int PLANAR_SPAN_SAMPLES = 8;

// make_planar_edge_sampler(span) reads the arc table from the span's interior
// samples at stride 2.
static_assert(PLANAR_SPAN_SAMPLES == 2 * PLANAR_LEN_SAMPLES);

/**
 * @brief Shared per-edge planar setup for the row/column span bounds.
 * @details Projects the edge and samples its chart line through the
 *          rasterizer's unprojection map once per edge. interior holds the
 *          k/PLANAR_SPAN_SAMPLES samples for k in [1, PLANAR_SPAN_SAMPLES),
 *          excluding the endpoint.
 */
struct PlanarEdgeSpan {
  std::pair<float, float> p1; /**< Projection of the edge start. */
  float dX;                   /**< Projected chord x-component. */
  float dY;                   /**< Projected chord y-component. */
  float gap_arc;              /**< Bound on each inter-sample arc length. */
  std::array<math::Vector, PLANAR_SPAN_SAMPLES - 1> interior;
};

/**
 * @brief Computes the shared planar edge setup once per edge.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param planar_basis Azimuthal-equidistant projection basis.
 */
static inline PlanarEdgeSpan
make_planar_edge_span(const math::Vector &a, const math::Vector &b,
                      const math::Basis &planar_basis) {
  PlanarEdgeSpan es;
  es.p1 = azimuthal_project(a, planar_basis);
  auto p2 = azimuthal_project(b, planar_basis);
  es.dX = p2.first - es.p1.first;
  es.dY = p2.second - es.p1.second;
  // The projected chord over-estimates the on-sphere arc, so gap_arc bounds
  // each inter-sample arc length.
  es.gap_arc = sqrtf(es.dX * es.dX + es.dY * es.dY) / PLANAR_SPAN_SAMPLES;
  for (int k = 1; k < PLANAR_SPAN_SAMPLES; ++k) {
    float p = static_cast<float>(k) / PLANAR_SPAN_SAMPLES;
    es.interior[k - 1] = azimuthal_unproject(
        es.p1.first + es.dX * p, es.p1.second + es.dY * p, planar_basis);
  }
  return es;
}

/**
 * @brief Builds a planar sampler from the exact quarter-point cull samples.
 * @param span Planar cull setup whose interior samples cover eighth points.
 * @param end Unprojected edge endpoint from planar_col_span().
 * @param planar_basis Azimuthal-equidistant projection basis.
 */
static inline PlanarEdgeSampler
make_planar_edge_sampler(const PlanarEdgeSpan &span, const math::Vector &end,
                         const math::Basis &planar_basis) {
  PlanarEdgeSampler sampler;
  sampler.proj1 = span.p1;
  sampler.dx = span.dX;
  sampler.dy = span.dY;
  sampler.basis = &planar_basis;
  sampler.chart_tangent =
      (planar_basis.u * sampler.dx) + (planar_basis.w * sampler.dy);
  HS_PLOT_ADD(planar_arc_samples, PLANAR_LEN_SAMPLES + 1);
  sampler.arc_cumul[0] = 0.0f;
  math::Vector prev =
      azimuthal_unproject(span.p1.first, span.p1.second, planar_basis);
  for (int k = 1; k < PLANAR_LEN_SAMPLES; ++k) {
    const math::Vector &cur = span.interior[k * 2 - 1];
    sampler.arc_cumul[k] =
        sampler.arc_cumul[k - 1] + math::angle_between(prev, cur);
    prev = cur;
  }
  sampler.arc_cumul[PLANAR_LEN_SAMPLES] =
      sampler.arc_cumul[PLANAR_LEN_SAMPLES - 1] +
      math::angle_between(prev, end);
  sampler.dist = sampler.arc_cumul[PLANAR_LEN_SAMPLES];
  return sampler;
}
