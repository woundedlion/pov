/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/plot/cull.h.

/**
 * @brief Screen row of a unit-sphere y coordinate (the renderer's row map).
 * @tparam H Rasterization height (pixel grid).
 * @param y Unit-sphere y in [-1, 1] (clamped).
 */
template <int H> static inline float y_to_screen_row(float y) {
  return math::phi_to_y<H>(math::fast_acos(hs::clamp(y, -1.0f, 1.0f)));
}

/**
 * @brief Geodesic screen-row span from precomputed endpoint rows.
 * @tparam H Rasterization height (pixel grid).
 * @param ra Precomputed y_to_screen_row<H>(a.y).
 * @param rb Precomputed y_to_screen_row<H>(b.y).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 * @details The arc y(t) has a turning point inside the span iff the forward
 * tangent's y-component flips sign between the endpoints; the extremal |y| is
 * the great circle's peak latitude sqrt(1 - n.y²) (n = arc pole). The span is
 * the exact closed-form y range, so no one-row epsilon. A degenerate setup
 * (no axis) keeps the endpoint rows.
 */
template <int H>
static __attribute__((always_inline)) inline void
geodesic_row_span_rows(float ra, float rb, const math::Vector &a,
                       const math::Vector &b, const GeodesicEdgeSpan &es,
                       float &row_lo, float &row_hi) {
  row_lo = fminf(ra, rb);
  row_hi = fmaxf(ra, rb);
  if (!es.have_axis)
    return;
  float t0 = math::cross(es.axis, a).y; // forward tangent y at a
  float t1 = math::cross(es.axis, b).y; // forward tangent y at b
  if ((t0 > 0.0f) != (t1 > 0.0f)) {
    // fmaxf(0, ...) absorbs the tiny negative that fast-math
    // renormalization of the axis can produce when |axis.y| ≈ 1 (a
    // near-polar arc pole), keeping the sqrt domain-safe.
    float peak = sqrtf(fmaxf(0.0f, 1.0f - es.axis.y * es.axis.y));
    float rp = y_to_screen_row<H>(t0 > 0.0f ? peak : -peak);
    row_lo = fminf(row_lo, rp);
    row_hi = fmaxf(row_hi, rp);
  }
}

/**
 * @brief Geodesic screen-row span from a precomputed edge setup.
 * @tparam H Rasterization height (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 */
template <int H>
static inline void
geodesic_row_span(const math::Vector &a, const math::Vector &b,
                  const GeodesicEdgeSpan &es, float &row_lo, float &row_hi) {
  geodesic_row_span_rows<H>(y_to_screen_row<H>(a.y), y_to_screen_row<H>(b.y), a,
                            b, es, row_lo, row_hi);
}

/**
 * @brief Planar screen-row span from a precomputed edge setup.
 * @tparam H Rasterization height (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_planar_edge_span(a, b, basis).
 * @param row_lo Output: minimum screen row touched by the edge.
 * @param row_hi Output: maximum screen row touched by the edge.
 * @details No closed-form latitude extremum exists, so the endpoint rows are
 * extended over the shared samples and widened by the arc's Lipschitz bound.
 * The cull and renderer do NOT take bit-identical samples, so gap-freeness
 * comes from the Lipschitz + one-row margin: phi is 1-Lipschitz in angular
 * distance, so between samples |Δrow| ≤ (Δarc)·ROWS_PER_RADIAN<H>. The samples take
 * the renderer's newton_unit() correction first: phi = acos(y) amplifies the
 * fast-trig residual on the raw unprojection past the one-row epsilon once
 * sin(phi) falls under a few hundredths.
 */
template <int H>
static inline void planar_row_span(const math::Vector &a, const math::Vector &b,
                                   const PlanarEdgeSpan &es, float &row_lo,
                                   float &row_hi) {
  float ra = y_to_screen_row<H>(a.y);
  float rb = y_to_screen_row<H>(b.y);
  row_lo = fminf(ra, rb);
  row_hi = fmaxf(ra, rb);
  for (const math::Vector &s : es.interior) {
    float r = y_to_screen_row<H>(newton_unit(s).y);
    row_lo = fminf(row_lo, r);
    row_hi = fmaxf(row_hi, r);
  }
  float margin = es.gap_arc * math::ROWS_PER_RADIAN<H> + 1.0f;
  row_lo -= margin;
  row_hi += margin;
}

/**
 * @brief Pads a fractional column interval into an integer [lo, lo + len) run.
 * @tparam W Rasterization width (pixel grid).
 * @param start Fractional start column.
 * @param length Fractional arc length in columns.
 * @param lo Output: padded start column, unwrapped.
 * @param col_len Output: arc length in columns (may reach W = full width).
 */
template <int W>
static __attribute__((always_inline)) inline void
pad_col_span(float start, float length, int &lo, int &col_len) {
  lo = static_cast<int>(floorf(start)) - COL_PAD;
  const int hi = static_cast<int>(ceilf(start + length)) + COL_PAD;
  col_len = std::min(hi - lo + 1, W);
}

/**
 * @brief Pads a fractional column interval and wraps it into a [0, W) arc.
 * @tparam W Rasterization width (pixel grid).
 * @param s_f Fractional start column.
 * @param len_f Fractional arc length in columns.
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 */
template <int W>
static __attribute__((always_inline)) inline void
finish_col_span(float s_f, float len_f, int &col_s, int &col_len) {
  int lo;
  pad_col_span<W>(s_f, len_f, lo, col_len);
  col_s = ((lo % W) + W) % W;
}

/**
 * @brief Pads a column span whose start is already within one period.
 */
template <int W>
static inline void finish_col_span_one_period(float start, float length,
                                              int &col_s, int &col_len) {
  int lo;
  pad_col_span<W>(start, length, lo, col_len);
  col_s = lo < 0 ? lo + W : lo;
}

/**
 * @brief Wraps a column delta into [0, W) with a single conditional add.
 * @tparam W Rasterization width (pixel grid).
 * @param d Delta to wrap; must lie in (-W, W).
 * @return The wrapped delta in [0, W), bit-identical to wrap(d, W) over that
 *   domain.
 * @details `d + W` can round up to exactly W for a tiny negative d; the upper
 * guard folds that back to 0, as wrap's own half-open guard does.
 */
template <int W>
static __attribute__((always_inline)) inline float wrap_one_period(float d) {
  assert(d > -W && d < W);
  if (d < 0.0f)
    d += W;
  return (d >= W) ? 0.0f : d;
}

/**
 * @brief Geodesic screen-column arc from precomputed endpoint columns.
 * @tparam W Rasterization width (pixel grid).
 * @param ca Precomputed vector_to_theta<W>(a).
 * @param cb Precomputed vector_to_theta<W>(b).
 * @param a Edge start (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 * @return False when the arc pole has |y| below AXIS_Y_EPS, making longitude
 *         ill-conditioned; the caller must skip the horizontal cull.
 * @details Longitude is globally monotone along the rendered circle — with
 * pos(ang) = a·cos + cross(axis, a)·sin, the atan2(z, x) rate numerator
 * pos.x·tan.z - pos.z·tan.x folds to -axis.y, a constant. The arc therefore
 * sweeps from a's column toward the end column in the direction sign(-axis.y),
 * and one full revolution sweeps exactly W, so the directed modular difference
 * is the exact sweep. Antipodal symmetry (λ(-p) = λ(p) + π) makes every
 * half-circle sweep exactly W/2 and shorter arcs less, so the span is always
 * the endpoints' short-way separation — the direction only disambiguates the
 * near-antipodal boundary, where short-way is float noise. The column mapping
 * is the renderer's vector_to_theta.
 */
template <int W>
static inline bool
geodesic_col_span_cols(float ca, float cb, const math::Vector &a,
                       const GeodesicEdgeSpan &es, int &col_s, int &col_len) {
  float s_f, len_f;

  if (es.total < EPS_GEODESIC_SEGMENT) {
    // The renderer collapses the edge to a dot at a; span both endpoints the
    // short way around.
    const float d = wrap_one_period<W>(cb - ca);
    if (d <= W * 0.5f) {
      s_f = ca;
      len_f = d;
    } else {
      s_f = cb;
      len_f = W - d;
    }
  } else {
    if (!es.azimuth_bounded)
      return false;

    float ce;
    if (es.antipodal) {
      // The arbitrary-axis half-turn lands near, not on, b; take the column of
      // the point the renderer actually reaches.
      math::Vector v_perp = math::cross(es.axis, a);
      math::Vector end =
          a * math::fast_cosf(es.total) + v_perp * math::fast_sinf(es.total);
      ce = math::vector_to_theta<W>(end);
    } else {
      ce = cb;
    }

    if (es.axis.y < 0.0f) { // longitude increases from a
      s_f = ca;
      len_f = wrap_one_period<W>(ce - ca);
    } else {
      s_f = ce;
      len_f = wrap_one_period<W>(ca - ce);
    }
  }

  finish_col_span_one_period<W>(s_f, len_f, col_s, col_len);
  return true;
}

/**
 * @brief Geodesic screen-column arc from a precomputed edge setup.
 * @tparam W Rasterization width (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 * @return False when no useful bound exists — the caller must skip the
 *         horizontal cull.
 */
template <int W>
static inline bool
geodesic_col_span(const math::Vector &a, const math::Vector &b,
                  const GeodesicEdgeSpan &es, int &col_s, int &col_len) {
  return geodesic_col_span_cols<W>(math::vector_to_theta<W>(a),
                                   math::vector_to_theta<W>(b), a, es, col_s,
                                   col_len);
}

/**
 * @brief Most arc fractions a clip band can cut one geodesic edge at.
 * @details Two meridians met once each, two latitude circles met twice each.
 */
inline constexpr int GEODESIC_CLIP_MAX_SPLITS = 6;

/**
 * @brief The clip band's cut boundaries, in the terms the arc solve reads them.
 */
struct ClipCutBounds {
  float col_x[2]; /**< Boundary meridian directions, x component. */
  float col_z[2]; /**< Boundary meridian directions, z component. */
  float row_y[2]; /**< Boundary latitudes, as cos(phi). */
  bool cols;      /**< Column boundaries cut (x clipping is active). */
  bool rows;      /**< Row boundaries cut (band shorter than the canvas). */
};

/**
 * @brief Resolves the clip band's cut boundaries for a whole draw.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @return Boundary geometry to hand every geodesic_clip_splits call under
 *         @p cr.
 * @details The band is fixed for a draw, so its boundary directions and
 * latitudes resolve once rather than per edge.
 */
template <int W, int H>
static inline ClipCutBounds make_clip_cut_bounds(const ClipRegion &cr,
                                                 const ClipRegion::XClip &xc) {
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();

  ClipCutBounds cb{};
  cb.cols = xc.active;
  cb.rows = cr.render_y_start() > 0 || cr.render_y_end() < cr.h;

  if (cb.cols) {
    const int cols[2] = {xc.rs - CLIP_CUT_COL_PAD, xc.re + CLIP_CUT_COL_PAD};
    for (int i = 0; i < 2; ++i) {
      const int c = ((cols[i] % W) + W) % W;
      cb.col_x[i] = math::TrigLUT<W, H>::cos_theta(c);
      cb.col_z[i] = math::TrigLUT<W, H>::sin_theta[c];
    }
  }
  if (cb.rows) {
    const int rows[2] = {cr.render_y_start() - CLIP_CUT_ROW_PAD -
                             static_cast<int>(GEODESIC_ROW_AA_PAD),
                         cr.render_y_end() + CLIP_CUT_ROW_PAD};
    for (int i = 0; i < 2; ++i)
      cb.row_y[i] = math::TrigLUT<W, H>::cos_phi[hs::clamp(rows[i], 0, H - 1)];
  }
  return cb;
}

/**
 * @brief Arc fractions where a geodesic edge crosses the clip band's row and
 *        column boundaries.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b); must have an axis.
 * @param cb Clip-band boundaries from make_clip_cut_bounds.
 * @param ts Output, up to GEODESIC_CLIP_MAX_SPLITS fractions in (0, 1),
 *        ascending and separated enough that no piece is degenerate.
 * @return Number of fractions written.
 * @details Every resulting piece lies wholly inside or wholly outside the band. The
 * boundaries are the RENDER band's, already widened by the clip margin to cover
 * filter reach; the only spacing added on top is the cull's own footprint
 * (CLIP_CUT_COL_PAD, CLIP_CUT_ROW_PAD), without which the outside piece lands
 * inside the span pad and is kept anyway.
 *
 * Both solves run against pos(ang) = a·cos(ang) + cross(axis, a)·sin(ang), the
 * arc the renderer walks, over the same boundary directions the projection
 * reads (TrigLUT). A boundary meridian's half-plane, direction d in the (x, z)
 * plane, is met where the cross with d vanishes and the dot is positive: one
 * atan2, resolved to the sign the half-plane wants. Longitude is globally
 * monotone along the circle and one edge sweeps at most half a revolution
 * (geodesic_col_span_cols), so a meridian is met at most once; the near-
 * horizontal arc pole that same bound refuses is skipped here too. A boundary
 * row is a latitude, and y(ang) folds to R·cos(ang − delta), so its two
 * crossings come from one acos.
 *
 * The solve's fast trig puts a cut within a small fraction of a pixel of the
 * boundary. That is harmless in either direction: each piece is gated by the
 * exact span of the arc between its own endpoints, so a mis-sided cut draws a
 * piece rather than dropping one.
 */
static inline int geodesic_clip_splits(const math::Vector &a,
                                       const math::Vector &b,
                                       const GeodesicEdgeSpan &es,
                                       const ClipCutBounds &cb, float *ts) {
  const math::Vector perp = math::cross(es.axis, a);
  float angs[GEODESIC_CLIP_MAX_SPLITS];
  int found = 0;

  // Roots repeat every turn and the edge sweeps at most half of one, so each
  // candidate has a single representative in [0, 2pi).
  const auto keep = [&](float ang) {
    if (ang < 0.0f)
      ang += 2.0f * math::PI_F;
    else if (ang >= 2.0f * math::PI_F)
      ang -= 2.0f * math::PI_F;
    if (ang > 0.0f && ang < es.total) {
      HS_CHECK(found < GEODESIC_CLIP_MAX_SPLITS,
               "geodesic clip: more than %d split roots on one edge",
               GEODESIC_CLIP_MAX_SPLITS);
      angs[found++] = ang;
    }
  };

  if (cb.cols && es.azimuth_bounded) {
    for (int i = 0; i < 2; ++i) {
      const float dx = cb.col_x[i];
      const float dz = cb.col_z[i];
      const float cross_a = a.x * dz - a.z * dx;
      const float cross_b = b.x * dz - b.z * dx;
      // The cross runs sinusoidally in the arc angle, so over the at-most-half
      // turn the edge sweeps it has one root; equal end signs put that root
      // outside the arc, and the two roots an exactly antipodal edge could hold
      // sit on its endpoints, which keep() drops.
      if ((cross_a < 0.0f) == (cross_b < 0.0f))
        continue;
      const float cross_p = perp.x * dz - perp.z * dx;
      const float dot_a = a.x * dx + a.z * dz;
      const float dot_p = perp.x * dx + perp.z * dz;
      // sin, cos at the root are (-cross_a, cross_p) up to a positive scale, so
      // the half-plane's sign test needs no second trig call.
      const float ang = math::fast_atan2(-cross_a, cross_p);
      keep(dot_a * cross_p - dot_p * cross_a < 0.0f ? ang + math::PI_F : ang);
    }
  }

  if (cb.rows) {
    const float radius2 = a.y * a.y + perp.y * perp.y;
    if (radius2 > 0.0f) {
      const float radius = sqrtf(radius2);
      const float delta = math::fast_atan2(perp.y, a.y);
      // y folds to radius*cos(ang - delta), so the endpoints bound the arc
      // except where an extremum angle falls inside it.
      float y_lo = fminf(a.y, b.y);
      float y_hi = fmaxf(a.y, b.y);
      if (delta > 0.0f && delta < es.total)
        y_hi = radius;
      if (delta + math::PI_F < es.total)
        y_lo = -radius;
      for (int i = 0; i < 2; ++i) {
        const float y = cb.row_y[i];
        if (y < y_lo || y > y_hi)
          continue;
        const float half = math::fast_acos(hs::clamp(y / radius, -1.0f, 1.0f));
        keep(delta - half);
        keep(delta + half);
      }
    }
  }

  for (int i = 1; i < found; ++i) {
    const float key = angs[i];
    int j = i - 1;
    for (; j >= 0 && angs[j] > key; --j)
      angs[j + 1] = angs[j];
    angs[j + 1] = key;
  }

  // A piece shorter than the renderer's own collapse threshold would draw as a
  // dot, so fold coincident cuts (and cuts sitting on an endpoint) away.
  int n = 0;
  float last = 0.0f;
  for (int i = 0; i < found; ++i) {
    if (angs[i] - last < EPS_GEODESIC_SEGMENT ||
        es.total - angs[i] < EPS_GEODESIC_SEGMENT)
      continue;
    last = angs[i];
    ts[n++] = angs[i] / es.total;
  }
  return n;
}

/**
 * @brief Exact clip visibility of one geodesic edge.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam ColSpanFn bool(int &col_s, int &col_len), the edge's column arc;
 *         false when no useful bound exists.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param ra Precomputed y_to_screen_row<H>(a.y).
 * @param rb Precomputed y_to_screen_row<H>(b.y).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @param col_span Column-arc source, evaluated only once the row span survives
 *        and x clipping is active.
 * @return True if the rendered edge could produce a pixel inside the clip.
 * @details The exact geodesic segment cull used by rasterize through
 * edge_visible_in_clip. The parity-tested raw_geodesic_edge_gate is a fast
 * path for particle trails and defers here for sensitive geometry.
 */
template <int W, int H, typename ColSpanFn>
static __attribute__((always_inline)) inline bool
exact_geodesic_edge_visible(const ClipRegion &cr, const ClipRegion::XClip &xc,
                            float ra, float rb, const math::Vector &a,
                            const math::Vector &b, const GeodesicEdgeSpan &es,
                            ColSpanFn &&col_span) {
  float row_lo, row_hi;
  geodesic_row_span_rows<H>(ra, rb, a, b, es, row_lo, row_hi);
  if (!cr.could_intersect_y(row_lo, row_hi + GEODESIC_ROW_AA_PAD))
    return false;
  if (!xc.active)
    return true;
  int col_s, col_len;
  if (!col_span(col_s, col_len))
    return true;
  return ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W);
}

/**
 * @brief Exact clip visibility of one geodesic edge, from a trail's hoisted
 *        per-point screen coordinates.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param rows Per-point screen rows from trail_gate_prologue.
 * @param cols Per-point screen columns from trail_gate_prologue; read only when
 *        x clipping is active, so it may be null otherwise.
 * @param e Edge index; the edge runs from point @p e to point @p e + 1.
 * @param a Edge start, trail point @p e.
 * @param b Edge end, trail point @p e + 1.
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @return True if the rendered edge could produce a pixel inside the clip.
 */
template <int W, int H>
static __attribute__((always_inline)) inline bool
exact_geodesic_edge_visible_hoisted(const ClipRegion &cr,
                                    const ClipRegion::XClip &xc,
                                    const float *rows, const float *cols,
                                    size_t e, const math::Vector &a,
                                    const math::Vector &b,
                                    const GeodesicEdgeSpan &es) {
  return exact_geodesic_edge_visible<W, H>(
      cr, xc, rows[e], rows[e + 1], a, b, es, [&](int &col_s, int &col_len) {
        return geodesic_col_span_cols<W>(cols[e], cols[e + 1], a, es, col_s,
                                         col_len);
      });
}

enum class RawGeodesicGateResult : uint8_t {
  CULLED,
  VISIBLE,
  EXACT_FALLBACK,
};

/**
 * @brief Gates a regular geodesic edge without angle or cross normalization.
 * @return Visibility, or EXACT_FALLBACK for numerically sensitive geometry.
 */
template <int W, int H>
static inline RawGeodesicGateResult
raw_geodesic_edge_gate(const ClipRegion &cr, const ClipRegion::XClip &xc,
                       float ra, float rb, float ca, float cb,
                       const math::Vector &a, const math::Vector &b) {
  constexpr float END_GUARD2 = 4.0e-6f;
  constexpr float AXIS_GUARD2 = 1.0e-4f;
  constexpr float TANGENT_GUARD2 = 1.0e-8f;
  constexpr float ROW_BOUNDARY_GUARD = 0.01f;
  const math::Vector c = math::cross(a, b);
  const float L2 = math::dot(c, c);
  const float d = math::dot(a, b);
  if (L2 <= END_GUARD2 || std::abs(d) >= 1.0f - END_GUARD2 * 0.5f)
    return RawGeodesicGateResult::EXACT_FALLBACK;

  const float cy2 = c.y * c.y;
  if (cy2 <= AXIS_GUARD2 * L2)
    return RawGeodesicGateResult::EXACT_FALLBACK;

  float row_lo = fminf(ra, rb);
  float row_hi = fmaxf(ra, rb);
  const float t0 = c.z * a.x - c.x * a.z;
  const float t1 = c.z * b.x - c.x * b.z;
  if (t0 * t0 <= TANGENT_GUARD2 * L2 || t1 * t1 <= TANGENT_GUARD2 * L2)
    return RawGeodesicGateResult::EXACT_FALLBACK;
  if ((t0 > 0.0f) != (t1 > 0.0f)) {
    const float peak = sqrtf(fmaxf(0.0f, (L2 - cy2) / L2));
    const float rp = y_to_screen_row<H>(t0 > 0.0f ? peak : -peak);
    row_lo = fminf(row_lo, rp);
    row_hi = fmaxf(row_hi, rp);
  }

  const float y_start = static_cast<float>(cr.render_y_start());
  const float y_end = static_cast<float>(cr.render_y_end());
  const float padded_row_hi = row_hi + GEODESIC_ROW_AA_PAD;
  if (std::abs(padded_row_hi - y_start) < ROW_BOUNDARY_GUARD ||
      std::abs(row_lo - y_end) < ROW_BOUNDARY_GUARD)
    return RawGeodesicGateResult::EXACT_FALLBACK;
  if (!cr.could_intersect_y(row_lo, padded_row_hi))
    return RawGeodesicGateResult::CULLED;
  if (!xc.active)
    return RawGeodesicGateResult::VISIBLE;

  float col_start;
  float col_length;
  if (c.y < 0.0f) {
    col_start = ca;
    col_length = cb - ca;
  } else {
    col_start = cb;
    col_length = ca - cb;
  }
  if (col_length < 0.0f)
    col_length += W;

  int col_s, col_len;
  finish_col_span_one_period<W>(col_start, col_length, col_s, col_len);
  return ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W)
             ? RawGeodesicGateResult::VISIBLE
             : RawGeodesicGateResult::CULLED;
}

/**
 * @brief Planar screen-column arc from a precomputed edge setup.
 * @tparam W Rasterization width (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param planar_basis Azimuthal-equidistant projection basis.
 * @param es Shared setup from make_planar_edge_span(a, b, planar_basis); the
 *           end point enters through its projected chord.
 * @param col_s Output: arc start column, in [0, W).
 * @param col_len Output: arc length in columns (may reach W = full width).
 * @param end_sample Optional output for the unprojected edge endpoint.
 * @return False when no useful bound exists — the edge nears a pole or the
 *         Lipschitz margin exceeds the short-way-delta proof — the caller
 *         must skip the horizontal cull.
 * @details Longitude is not monotone along the chart line, so accumulate
 * short-way column deltas over the shared samples and widen by the azimuth
 * Lipschitz bound (W/2π)/sin(φ) per inter-sample gap. sin(φ) over the whole
 * edge is bounded below from the samples minus the 1-Lipschitz gap drift; the
 * short-way delta reading is valid only while the per-gap bound stays under
 * W/4, and both that and the near-pole case fall back to no-cull (return
 * false). The sin(φ) reads take the renderer's newton_unit() correction first:
 * the raw unprojection is off unit by the fast-trig residual, which near a pole
 * outweighs sin(φ) itself and would shrink the margin. The column mapping is
 * scale-invariant, so it reads the raw sample.
 */
template <int W>
static inline bool
planar_col_span(const math::Vector &a, const math::Basis &planar_basis,
                const PlanarEdgeSpan &es, int &col_s, int &col_len,
                math::Vector *end_sample = nullptr) {
  const float ca = math::vector_to_theta<W>(a);
  float s_f, len_f;

  {
    float min_sp2 =
        1.0f - a.y * a.y; // squared sin(phi), minimized over samples
    float cum = 0.0f, cum_lo = 0.0f, cum_hi = 0.0f;
    float prev = ca;
    auto step = [&](const math::Vector &s) {
      float c = math::vector_to_theta<W>(s);
      float d = c - prev;
      if (d > W * 0.5f)
        d -= W;
      else if (d < -W * 0.5f)
        d += W;
      cum += d;
      cum_lo = fminf(cum_lo, cum);
      cum_hi = fmaxf(cum_hi, cum);
      prev = c;
      const float sy = newton_unit(s).y;
      min_sp2 = fminf(min_sp2, 1.0f - sy * sy);
    };
    for (const math::Vector &s : es.interior)
      step(s);
    math::Vector end = azimuthal_unproject(es.p1.first + es.dX,
                                           es.p1.second + es.dY, planar_basis);
    if (end_sample != nullptr)
      *end_sample = end;
    step(end);

    // Worst-case sin(phi) anywhere on the edge: phi is 1-Lipschitz in arc
    // length and sin is 1-Lipschitz in phi, so between samples it drifts by at
    // most gap_arc.
    const float sin_phi_worst = sqrtf(fmaxf(0.0f, min_sp2)) - es.gap_arc;
    if (sin_phi_worst < MIN_SIN_PHI)
      return false;
    // Column movement inside one gap; also the proof bound for reading each
    // sample-to-sample delta the short way (must stay well under W/2).
    const float margin = es.gap_arc *
                             (static_cast<float>(W) / (2.0f * math::PI_F)) /
                             sin_phi_worst +
                         1.0f;
    if (margin >= W * 0.25f)
      return false;

    s_f = ca + cum_lo - margin;
    len_f = (cum_hi - cum_lo) + 2.0f * margin;
  }

  finish_col_span<W>(s_f, len_f, col_s, col_len);
  return true;
}

/**
 * @brief Tests planar edge visibility from an existing cull sample set.
 * @param cr Active raster clip region.
 * @param xc Precomputed horizontal clip interval.
 * @param a Unit-sphere edge start.
 * @param b Unit-sphere edge end.
 * @param planar_basis Basis used to unproject planar samples.
 * @param span Interior samples and arc-gap bound for the edge.
 * @param end_sample Optional output for the unprojected edge endpoint, written
 *        by the column-span pass. Only an active @p xc reaches that pass, so
 *        requesting the endpoint requires one.
 */
template <int W, int H>
static inline bool planar_edge_visible_in_clip(
    const ClipRegion &cr, const ClipRegion::XClip &xc, const math::Vector &a,
    const math::Vector &b, const math::Basis &planar_basis,
    const PlanarEdgeSpan &span, math::Vector *end_sample = nullptr) {
  assert(end_sample == nullptr || xc.active);
  float row_lo, row_hi;
  int col_s, col_len;
  planar_row_span<H>(a, b, span, row_lo, row_hi);
  if (!cr.could_intersect_y(row_lo, row_hi))
    return false;
  if (!xc.active)
    return true;
  if (!planar_col_span<W>(a, planar_basis, span, col_s, col_len, end_sample))
    return true;
  return ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W);
}

/**
 * @brief One-refinement reciprocal square root for adaptive screen spacing.
 * @details The biased refinement bounds relative error to about 0.089%.
 */
static __attribute__((always_inline)) inline float screen_rsqrt(float x) {
  uint32_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  bits = 0x5f37642fu - (bits >> 1);
  float estimate;
  std::memcpy(&estimate, &bits, sizeof(estimate));
  return estimate * (1.500883f - 0.5f * x * estimate * estimate);
}

/**
 * @brief Adaptive sub-step length (radians of arc) for ~one-pixel screen steps.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param pos_y y-component of the unit-sphere sample position.
 * @param tan_y y-component of the unit tangent (with respect to arc length).
 * @param cross_y pos.x·tan.z - pos.z·tan.x, the longitude-rate numerator.
 * @param base_step Equatorial step 2π/W; also the maximum returned step.
 * @return Arc-length step that advances ~SCREEN_STEP_PX pixels on screen.
 * @details Converts the object-space tangent to a screen-space velocity (pixels
 * per radian of arc) under the canvas map x = θ·W/2π, y = (φ - NORTH_PHI)·ROWS_PER_RADIAN<H>, then
 * returns step = SCREEN_STEP_PX / |v_screen|. With φ the colatitude and λ the
 * longitude, dφ/ds = -tan.y/sin(φ) and dλ/ds = (pos.x·tan.z - pos.z·tan.x)/sin²φ.
 * Tracking the full 2-D screen speed (not just longitudinal pole-crowding)
 * deposits ~one sample per pixel everywhere on the curve.
 *
 * Clamped to [base_step·MIN_POLE_SCALE, base_step]: the lower bound caps
 * oversampling at the poles (where dλ/ds diverges → speed → ∞ → step → 0); the
 * upper bound keeps the equator near one sample per column.
 */
template <int W, int H>
static inline float screen_step_components(float pos_y, float tan_y,
                                           float cross_y, float base_step) {
  const float KX = W / (2.0f * math::PI_F);  // columns per radian of longitude
  const float KY = math::ROWS_PER_RADIAN<H>; // rows per radian of colatitude
  // sin²φ = 1 - y²; floored so the pole (sin φ → 0) yields a finite, large
  // velocity (hence the min-clamped step) rather than a divide-by-zero.
  const float sin2 = fmaxf(1e-7f, 1.0f - pos_y * pos_y);
  const float vx_num = KX * cross_y;
  const float vy_num = KY * tan_y;
  // Factoring the common sin(phi) denominator avoids a separate reciprocal
  // square root while preserving the screen-speed floor.
  const float speed2_num =
      fmaxf(vx_num * vx_num + vy_num * vy_num * sin2, 1e-12f * sin2 * sin2);
  // Degenerate-speed floor: guards 1/speed when a zero/near-zero tangent stalls
  // the curve, yielding base_step rather than an unbounded step.
  const float step = SCREEN_STEP_PX * sin2 * screen_rsqrt(speed2_num);
  return hs::clamp(step, base_step * MIN_POLE_SCALE, base_step);
}

template <int W, int H>
static inline float screen_step(const math::Vector &pos,
                                const math::Vector &tan, float base_step) {
  return screen_step_components<W, H>(pos.y, tan.y,
                                      pos.x * tan.z - pos.z * tan.x, base_step);
}

#if HS_ENABLE_TEST_ORACLES
/** @brief Caps rasterize()'s per-segment sub-step budget; 0 leaves it alone. */
inline size_t g_step_budget_override = 0;
#endif

#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
inline bool g_reference_screen_step = false;

template <int W, int H>
static inline float screen_step_reference(const math::Vector &pos,
                                          const math::Vector &tan,
                                          float base_step) {
  const float KX = W / (2.0f * math::PI_F);
  const float KY = math::ROWS_PER_RADIAN<H>;
  const float sin2 = std::max(1e-7f, 1.0f - pos.y * pos.y);
  const float inv_sin = math::fast_rsqrt(sin2);
  const float dphi_ds = -tan.y * inv_sin;
  const float dlon_ds = (pos.x * tan.z - pos.z * tan.x) * inv_sin * inv_sin;
  const float vx = KX * dlon_ds;
  const float vy = KY * dphi_ds;
  const float speed2 = std::max(vx * vx + vy * vy, 1e-12f);
  const float step = SCREEN_STEP_PX * math::fast_rsqrt(speed2);
  return std::max(base_step * MIN_POLE_SCALE, std::min(step, base_step));
}
#endif
