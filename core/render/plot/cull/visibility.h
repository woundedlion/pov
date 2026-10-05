/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/plot/cull.h.

/**
 * @brief True when AntiAlias would emit any tap of a projected dot in @p cr.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param row Precomputed projected row.
 * @param col Precomputed projected column; unused when x clipping is inactive.
 * @details Tests the taps Screen::AntiAlias would emit, sharing splat_taps and
 * SPLAT_TAP_CUTOFF with it. The gate runs before shading, so it tests tap
 * geometry only.
 */
template <int W, int H>
static inline bool antialiased_dot_visible_in_clip(const ClipRegion &cr,
                                                   const ClipRegion::XClip &xc,
                                                   float row, float col) {
  const int y0 = static_cast<int>(floorf(row));
  const int y1 = y0 + 1;
  const bool y0_ok = y0 >= 0 && y0 < H;
  const bool y1_ok = y1 >= 0 && y1 < H;
  if ((!y0_ok || !cr.contains_y(y0)) && (!y1_ok || !cr.contains_y(y1)))
    return false;

  const Filter::Screen::SplatTaps t =
      Filter::Screen::splat_taps<W, H>(col, row);
  auto visible = [&](int x, int y, float weight) {
    return weight > Filter::Screen::SPLAT_TAP_CUTOFF && cr.contains_y(y) &&
           !xc.clipped(x);
  };
  return (t.y0_physical &&
          (visible(t.x0, t.y0, t.v00) || visible(t.x1, t.y0, t.v10))) ||
         (t.y1_physical &&
          (visible(t.x0, t.y1, t.v01) || visible(t.x1, t.y1, t.v11)));
}

template <typename PipelineT, typename Pred>
static inline bool
edge_visible_in_clip_dispatch(PipelineT &pipeline, const math::Vector &a,
                              const math::Vector &b, const math::Basis *pb,
                              Pred &&pred) {
  return pipeline_could_intersect_clip(pipeline, a, b, pb,
                                       std::forward<Pred>(pred));
}

/**
 * @brief World up axis, in object space, of every rigid world-stage copy.
 * @details screen_step reads only the y-components of the rotated position,
 * tangent and their cross product; for rotation R each equals a dot product of
 * the unrotated vector with Rᵀŷ. The copies' rotations do not depend on the
 * rasterized geometry, so one pass per stroke replaces a stage walk per sample.
 */
struct ScreenStepAxes {
  static constexpr int CAPACITY = 16;
  std::array<math::Vector, CAPACITY> up; /**< Rᵀŷ per copy. */
  int count = 0;                         /**< Live entries in up. */
  bool nonrigid = false; /**< A stage cannot report its copies. */
  bool overflow =
      false; /**< More copies than CAPACITY, or a stage forwarded a copy without its basis. */

  /** @brief True when the table replaces the per-sample stage walk. */
  bool usable() const { return !overflow; }
};

/**
 * @brief Minimum screen step over every world-stage copy, from their up axes.
 * @param sample Unrotated sample position and tangent.
 * @param axes Usable table from screen_step_axes().
 * @return The step pipeline_screen_step would return for @p sample.
 */
template <int W, int H>
__attribute__((always_inline)) inline float
screen_step_from_axes(const SamplePT &sample, const ScreenStepAxes &axes) {
  constexpr float BASE_STEP = (2.0f * math::PI_F) / W;
  if (math::dot(sample.tan, sample.tan) < math::EPS_NORMALIZE_SQ)
    return BASE_STEP;
  if (axes.nonrigid)
    return BASE_STEP * MIN_POLE_SCALE;
  const math::Vector n = math::cross(sample.pos, sample.tan);
  float step = BASE_STEP;
  for (int k = 0; k < axes.count; ++k) {
    const math::Vector &u = axes.up[k];
    step = fminf(step, screen_step_components<W, H>(
                           math::dot(sample.pos, u), math::dot(sample.tan, u),
                           -math::dot(n, u), BASE_STEP));
  }
  return step;
}

/**
 * @brief Collects each world-stage copy's up axis by sending the identity basis
 * through the pipeline's cull chain.
 * @param pipeline Pipeline whose world stages are walked.
 * @return The per-copy axes; nonrigid or overflow when the table is unusable.
 */
template <typename PipelineT>
HS_HOT_FLASH_MEMBER ScreenStepAxes screen_step_axes(PipelineT &pipeline) {
  ScreenStepAxes axes;
  const math::Basis identity{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  axes.nonrigid = edge_visible_in_clip_dispatch(
      pipeline, math::X_AXIS, math::Y_AXIS, &identity,
      [&](const math::Vector &, const math::Vector &, const math::Basis *rb) {
        if (!rb || axes.count == ScreenStepAxes::CAPACITY) {
          axes.overflow = true;
          return false;
        }
        axes.up[axes.count++] = math::Vector(rb->u.y, rb->v.y, rb->w.y);
        return false;
      });
  return axes;
}

/** @brief Screen step at the rendered latitude of every world-stage copy. */
template <int W, int H, typename PipelineT>
HS_HOT_FLASH_MEMBER float
pipeline_screen_step(PipelineT &pipeline, const SamplePT &sample,
                     bool world_identity,
                     const ScreenStepAxes *axes = nullptr) {
  constexpr float BASE_STEP = (2.0f * math::PI_F) / W;
  if (world_identity) {
#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
    if (g_reference_screen_step)
      return screen_step_reference<W, H>(sample.pos, sample.tan, BASE_STEP);
#endif
    return screen_step<W, H>(sample.pos, sample.tan, BASE_STEP);
  }
  if (math::dot(sample.tan, sample.tan) < math::EPS_NORMALIZE_SQ)
    return BASE_STEP;
  float step = BASE_STEP;
#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
  if (!g_reference_screen_step)
#endif
    if (axes && axes->usable())
      return screen_step_from_axes<W, H>(sample, *axes);
  // Rigid cull stages rotate both vectors; false visits every tween copy.
  const bool nonrigid = edge_visible_in_clip_dispatch(
      pipeline, sample.pos, sample.tan, nullptr,
      [&](const math::Vector &pos, const math::Vector &tan, const math::Basis *)
          HS_HOT_FLASH_MEMBER {
            float candidate;
#if HS_ENABLE_TEST_ORACLES || defined(HS_MINDSPLATTER_REPLAY)
            if (g_reference_screen_step)
              candidate = screen_step_reference<W, H>(pos, tan, BASE_STEP);
            else
#endif
              candidate = screen_step<W, H>(pos, tan, BASE_STEP);
            step = fminf(step, candidate);
            return false;
          });
  return nonrigid ? BASE_STEP * MIN_POLE_SCALE : step;
}

/**
 * @brief Clip visibility of one polyline edge, routed through the
 *        pipeline's world stages.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam PipelineT Pipeline type.
 * @param pipeline Render pipeline; world stages re-emit the edge under their
 *        plot-time rotations before the span test.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param pb Planar projection basis for the edge, or null for geodesic.
 * @return True if the rendered edge could produce a pixel inside the clip.
 * @details Geodesic edges route to exact_geodesic_edge_visible; the planar
 * branch is the single definition of the planar segment cull.
 */
template <int W, int H, typename PipelineT>
static inline bool
edge_visible_in_clip(PipelineT &pipeline, const ClipRegion &cr,
                     const ClipRegion::XClip &xc, const math::Vector &a,
                     const math::Vector &b, const math::Basis *pb) {
  auto pred = [&](const math::Vector &ea, const math::Vector &eb,
                  const math::Basis *bp) {
    if (bp == nullptr) {
      // Geodesic: both span bounds share one edge setup (angle, arc pole).
      const GeodesicEdgeSpan es = make_geodesic_edge_span(ea, eb);
      return exact_geodesic_edge_visible<W, H>(
          cr, xc, y_to_screen_row<H>(ea.y), y_to_screen_row<H>(eb.y), ea, eb,
          es, [&](int &col_s, int &col_len) {
            return geodesic_col_span<W>(ea, eb, es, col_s, col_len);
          });
    }
    // Planar: both span bounds share one projection + chart-line sample set.
    const PlanarEdgeSpan ps = make_planar_edge_span(ea, eb, *bp);
    return planar_edge_visible_in_clip<W, H>(cr, xc, ea, eb, *bp, ps);
  };
  return edge_visible_in_clip_dispatch(pipeline, a, b, pb, pred);
}

/**
 * @brief Clip visibility from a precomputed geodesic edge span.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @tparam PipelineT Pipeline type with no world cull stage.
 * @param pipeline Render pipeline.
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @param es Shared setup from make_geodesic_edge_span(a, b).
 * @return True if the rendered edge could produce a pixel inside the clip.
 */
template <int W, int H, typename PipelineT>
static inline bool
edge_visible_in_clip(PipelineT &pipeline, const ClipRegion &cr,
                     const ClipRegion::XClip &xc, const math::Vector &a,
                     const math::Vector &b, const GeodesicEdgeSpan &es) {
  static_assert(pipeline_hoistable_cull<PipelineT>(),
                "a precomputed geodesic span requires a hoistable cull");
  auto pred = [&](const math::Vector &ea, const math::Vector &eb,
                  const math::Basis *) {
    return exact_geodesic_edge_visible<W, H>(
        cr, xc, y_to_screen_row<H>(ea.y), y_to_screen_row<H>(eb.y), ea, eb, es,
        [&](int &col_s, int &col_len) {
          return geodesic_col_span<W>(ea, eb, es, col_s, col_len);
        });
  };
  return edge_visible_in_clip_dispatch(pipeline, a, b, nullptr, pred);
}

/** @brief Cap center and longitude wedge precomputed for a fixed clip. */
struct CapCenter {
  float beta;
  float sin_beta;
  float column_distance;
  float column_half_width;
  bool columns_active;
};

inline constexpr float CAP_ACOS_PAD = 1e-3f;
inline constexpr float CAP_ATAN2_PAD = 5e-3f;

inline CapCenter cap_column_center(const ClipRegion &cr,
                                   const math::Vector &dir, float beta) {
  CapCenter center{beta, 0, 0, 0, false};
  const auto columns = cr.x_clip();
  if (!columns.active)
    return center;
  center.columns_active = true;
  const float y = hs::clamp(dir.y, -1.0f, 1.0f);
  center.sin_beta = sqrtf(1.0f - y * y);
  const float width = static_cast<float>(columns.length(cr.w));
  center.column_half_width = (width * 0.5f + 1.0f) * math::TWO_PI_F / cr.w;
  const float longitude = math::fast_atan2(dir.z, dir.x);
  const float clip_longitude =
      (columns.rs + width * 0.5f) * math::TWO_PI_F / cr.w;
  center.column_distance =
      fabsf(math::wrap_t((longitude - clip_longitude) / math::TWO_PI_F + 0.5f) -
            0.5f) *
      math::TWO_PI_F;
  return center;
}

/** @brief Hoists a cap center's angles and column wedge for one clip. */
inline CapCenter make_cap_center(const ClipRegion &cr,
                                 const math::Vector &dir) {
  return cap_column_center(cr, dir,
                           math::fast_acos(hs::clamp(dir.y, -1.0f, 1.0f)));
}

template <int H, typename ColumnTest>
__attribute__((always_inline)) inline bool
cap_may_touch_clip_rows(const ClipRegion &cr, float beta, float half_angle,
                        ColumnTest &&columns) {
  const float t2 = fminf(half_angle, math::PI_F);
  const float phi_lo = fmaxf(beta - t2 - CAP_ACOS_PAD, 0.0f);
  const float phi_hi = fminf(beta + t2 + CAP_ACOS_PAD, math::PI_F);
  if (!cr.could_intersect_y(math::phi_to_y<H>(phi_lo),
                            math::phi_to_y<H>(phi_hi)))
    return false;
  return columns(t2);
}

inline bool cap_may_touch_clip_columns(const CapCenter &center, float t2,
                                       float sin_half_angle) {
  if (!center.columns_active || center.beta <= t2 + CAP_ACOS_PAD ||
      math::PI_F - center.beta <= t2 + CAP_ACOS_PAD)
    return true;
  const float dlam =
      math::PI_F / 2.0f -
      math::fast_acos(hs::clamp(sin_half_angle / center.sin_beta, 0.0f, 1.0f));
  return center.column_distance <=
         dlam + center.column_half_width + CAP_ACOS_PAD + CAP_ATAN2_PAD;
}

/** @brief Tests a cap using a center hoisted against the same clip. */
template <int H>
inline bool cap_may_touch_clip(const ClipRegion &cr, const CapCenter &center,
                               float half_angle, float sin_half_angle) {
  return cap_may_touch_clip_rows<H>(cr, center.beta, half_angle, [&](float t2) {
    return cap_may_touch_clip_columns(center, t2, sin_half_angle);
  });
}

/**
 * @brief Conservative test: can a spherical cap reach a clip's render region?
 * @tparam H Canvas height in rows.
 * @param cr Clip region to test against.
 * @param dir Cap center direction (unit vector). A ring passes its axis with
 * half_angle = colatitude + displacement bound (the cap of that radius
 * contains the ring's band); a ring chunk passes its midpoint with
 * half_angle = chunk half-arc + displacement bound.
 * @param half_angle Cap angular radius including stroke/AA pad (radians).
 * @param sin_half_angle sinf(min(half_angle, PI)), hoisted by callers testing
 * many caps of one radius.
 * @return False only when no fragment inside the cap can land in the clip's
 * render region; true is always safe.
 * @details Rows: the cap's polar range about the display's Y axis is
 * [beta - t2, beta + t2] (beta = center colatitude). Columns: a cap that
 * reaches either display pole spans all longitudes; otherwise its longitude
 * half-width about the center longitude is asin(sin t2 / sin beta), compared
 * against the clip's column wedge with the clip margin plus one pixel of
 * slack. The angles come from fast_acos and fast_atan2, each comparison
 * widened past their error bounds.
 */
template <int H>
inline bool cap_may_touch_clip(const ClipRegion &cr, const math::Vector &dir,
                               float half_angle, float sin_half_angle) {
  const float beta = math::fast_acos(hs::clamp(dir.y, -1.0f, 1.0f));
  return cap_may_touch_clip_rows<H>(cr, beta, half_angle, [&](float t2) {
    if (!cr.x_clip().active || beta <= t2 + CAP_ACOS_PAD ||
        math::PI_F - beta <= t2 + CAP_ACOS_PAD)
      return true;
    return cap_may_touch_clip_columns(cap_column_center(cr, dir, beta), t2,
                                      sin_half_angle);
  });
}

/** @brief Exclusive column boundary of an evenly divided, rounded-up chunk. */
template <int CHUNKS> constexpr int chunk_end(int c, int lut_n) {
  static_assert(CHUNKS > 0);
  return ((c + 1) * lut_n + CHUNKS - 1) / CHUNKS;
}

/** @brief Ring azimuth chunks reaching a clip, padded for stroke and column rounding. */
template <int H, int CHUNKS>
__attribute__((always_inline)) inline uint32_t
visible_chunk_mask(const ClipRegion &clip, const math::Basis &basis,
                   float theta, float cos_t, float sin_t, float band_r,
                   float thickness, const float *chunk_cos,
                   const float *chunk_sin) {
  static_assert(CHUNKS > 0 && CHUNKS < 32);
  constexpr uint32_t MASK = (1u << CHUNKS) - 1;
  const float chunk_reach = (math::PI_F / CHUNKS) * sin_t + band_r;
  const float sin_reach = sinf(fminf(chunk_reach, math::PI_F));
  uint32_t raw = 0u;
  for (int c = 0; c < CHUNKS; ++c) {
    math::Vector mid =
        (basis.v * cos_t) +
        ((basis.u * chunk_cos[c]) + (basis.w * chunk_sin[c])) * sin_t;
    if (cap_may_touch_clip<H>(clip, mid, chunk_reach, sin_reach))
      raw |= 1u << c;
  }
  if (!raw)
    return 0u;
  const float th_lo = theta - band_r;
  const float th_hi = theta + band_r;
  int pad_chunks = CHUNKS;
  if (th_lo > 0.0f && th_hi < math::PI_F) {
    float sin_lo = fminf(sinf(th_lo), sinf(th_hi));
    // A band hugging a pole drives sin_lo to zero; clamp before the cast.
    const float pad_f =
        ceilf(thickness * CHUNKS / (2.0f * math::PI_F * sin_lo));
    pad_chunks = 1 + static_cast<int>(
                         hs::clamp(pad_f, 0.0f, static_cast<float>(CHUNKS)));
  }
  if (2 * pad_chunks >= CHUNKS)
    return MASK;
  uint32_t visible = raw;
  for (int k = 1; k <= pad_chunks; ++k)
    visible |=
        (raw << k) | (raw >> (CHUNKS - k)) | (raw >> k) | (raw << (CHUNKS - k));
  return visible & MASK;
}

/** @brief cap_may_touch_clip() for a single cap, deriving sin(half_angle). */
template <int H>
inline bool cap_may_touch_clip(const ClipRegion &cr, const math::Vector &dir,
                               float half_angle) {
  return cap_may_touch_clip<H>(cr, dir, half_angle,
                               sinf(fminf(half_angle, math::PI_F)));
}

enum class CartesianTrailGateResult : uint8_t {
  EXACT_FALLBACK,
  LATITUDE_REJECT,
  MERIDIAN_REJECT,
};

struct CartesianQuadrantClip {
  float latitude_sign = 0.0f;
  float latitude_threshold = 0.0f;
  float meridian_sign = 0.0f;
  float meridian_threshold = 0.0f;
  bool active = false;
};

/**
 * @brief Builds the Cartesian halfspace superset for a segmented quadrant.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @return Enabled thresholds only for the four hardware quadrant shapes.
 */
template <int W, int H>
static CartesianQuadrantClip
make_cartesian_quadrant_clip(const ClipRegion &cr) {
  CartesianQuadrantClip q;
  if (cr.w != W || cr.h != H || cr.margin < 0 || W % 2 != 0 || H % 2 != 0 ||
      cr.x_end - cr.x_start != W / 2 ||
      (cr.x_start != 0 && cr.x_start != W / 2) ||
      cr.y_end - cr.y_start != H / 2 ||
      (cr.y_start != 0 && cr.y_start != H / 2))
    return q;

  if (cr.y_start == 0) {
    const float boundary = math::DisplayGeometry<H>::row_to_phi(
        static_cast<float>(cr.render_y_end()));
    q.latitude_sign = 1.0f;
    q.latitude_threshold = cosf(boundary);
  } else {
    const float boundary = math::DisplayGeometry<H>::row_to_phi(
        static_cast<float>(cr.render_y_start()) - GEODESIC_ROW_AA_PAD);
    q.latitude_sign = -1.0f;
    q.latitude_threshold = -cosf(boundary);
  }

  const float half_width =
      math::PI_F * 0.5f + static_cast<float>(cr.margin + COL_FOOTPRINT) *
                              (2.0f * math::PI_F / static_cast<float>(W));
  if (half_width >= math::PI_F)
    return CartesianQuadrantClip{};
  q.meridian_sign = cr.x_start == 0 ? 1.0f : -1.0f;
  q.meridian_threshold = cosf(half_width);
  q.active = true;
  return q;
}

/**
 * @brief Conservatively rejects a geodesic trail outside a Cartesian quadrant.
 * @param clip Precomputed quadrant halfspace thresholds.
 * @param trail Unit-sphere geodesic polyline.
 * @return Rejecting halfspace, or exact fallback for every uncertain case.
 */
HS_O3_BEGIN
static inline CartesianTrailGateResult
cartesian_quadrant_trail_gate(const CartesianQuadrantClip &clip,
                              const Fragments &trail) {
  if (!clip.active || trail.size() < 2)
    return CartesianTrailGateResult::EXACT_FALLBACK;

#ifdef HS_PROFILE_PLOT_STALLS
  hs::DwtStallBatch gate_batch(hs::g_plot_stalls.trail_gate);
#endif
  float latitude_max = -1.0f;
  float max_chord2 = 0.0f;
  for (size_t k = 0; k < trail.size(); ++k) {
    const math::Vector &p = trail[k].pos;
    latitude_max = fmaxf(latitude_max, clip.latitude_sign * p.y);
    if (k > 0) {
      const math::Vector d = p - trail[k - 1].pos;
      max_chord2 = fmaxf(max_chord2, math::dot(d, d));
    }
#ifdef HS_PROFILE_PLOT_STALLS
    gate_batch.step();
#endif
  }

  // Every point on a minor arc lies within half its arc of one endpoint, and
  // arc <= (pi/2)*chord. A unit-normal dot changes by at most angular distance.
  const float slack = (math::PI_F * 0.25f) * sqrtf(max_chord2);
  if (latitude_max + slack < clip.latitude_threshold - math::EPS_GEOMETRIC) {
#ifdef HS_PROFILE_PLOT_STALLS
    gate_batch.step();
    gate_batch.finish();
#endif
    return CartesianTrailGateResult::LATITUDE_REJECT;
  }

  float meridian_max = -1.0f;
  for (const Fragment &f : trail) {
    meridian_max = fmaxf(meridian_max, clip.meridian_sign * f.pos.z);
#ifdef HS_PROFILE_PLOT_STALLS
    gate_batch.step();
#endif
  }
#ifdef HS_PROFILE_PLOT_STALLS
  gate_batch.step();
  gate_batch.finish();
#endif
  if (meridian_max + slack < clip.meridian_threshold - math::EPS_GEOMETRIC)
    return CartesianTrailGateResult::MERIDIAN_REJECT;
  return CartesianTrailGateResult::EXACT_FALLBACK;
}
HS_O3_END

/**
 * @brief Hoisted per-point screen coordinates and whole-trail cull verdict.
 */
struct TrailGatePrologue {
  const float *rows; /**< Per-point screen rows, one per trail point. */
  const float
      *cols;     /**< Per-point screen columns, null when x-clip inactive. */
  bool rejected; /**< The whole trail is provably outside the clip. */
};

/**
 * @brief Computes one geodesic trail's per-point rows/columns and applies the
 *        conservative whole-trail row and column culls.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param cr Active clip region.
 * @param xc Precomputed x-clip predicate for @p cr.
 * @param trail Geodesic fragment polyline (>= 2 unit-position points).
 * @return The hoisted arrays plus whether the trail is culled whole.
 * @details rows and cols are bump-allocated from scratch_arena_a and the helper
 * opens no ScratchScope of its own: they stay valid until the CALLER's scope
 * unwinds, so a caller may keep them alive past the gate.
 */
template <int W, int H>
static __attribute__((always_inline)) inline TrailGatePrologue
trail_gate_prologue(const ClipRegion &cr, const ClipRegion::XClip &xc,
                    const Fragments &trail) {
  const size_t n = trail.size();
  auto *rows = scratch_arena_a.allocate_n<float>(n);
#ifdef HS_PROFILE_PLOT_STALLS
  hs::DwtStallBatch gate_batch(hs::g_plot_stalls.trail_gate);
#endif
  float row_lo_t = 1e9f, row_hi_t = -1e9f;
  float min_sp2 = 1.0f;
  float max_chord2 = 0.0f;
  for (size_t k = 0; k < n; ++k) {
    const math::Vector &pt = trail[k].pos;
    rows[k] = y_to_screen_row<H>(pt.y);
    row_lo_t = fminf(row_lo_t, rows[k]);
    row_hi_t = fmaxf(row_hi_t, rows[k]);
    min_sp2 = fminf(min_sp2, 1.0f - pt.y * pt.y);
    if (k > 0) {
      const math::Vector d = pt - trail[k - 1].pos;
      max_chord2 = fmaxf(max_chord2, math::dot(d, d));
    }
#ifdef HS_PROFILE_PLOT_STALLS
    gate_batch.step();
#endif
  }
  // arc <= (pi/2)*chord on [0, pi]; an edge's interior latitude extremum lies
  // within arc/2 of an endpoint and phi is 1-Lipschitz in arc length, so this
  // margin covers every per-edge bulge peak.
  const float max_arc = (math::PI_F * 0.5f) * sqrtf(max_chord2);
  const float row_margin = (max_arc * 0.5f) * math::ROWS_PER_RADIAN<H>;
  if (!cr.could_intersect_y(row_lo_t - row_margin,
                            row_hi_t + row_margin + GEODESIC_ROW_AA_PAD)) {
#ifdef HS_PROFILE_PLOT_STALLS
    gate_batch.step();
    gate_batch.finish();
#endif
    HS_PLOT_RENDER_COUNT(prologue_row_rejects);
    return {rows, nullptr, true};
  }

  float *cols = nullptr;
  if (xc.active) {
    cols = scratch_arena_a.allocate_n<float>(n);
    float cum = 0.0f, cum_lo = 0.0f, cum_hi = 0.0f;
    bool walk_safe = true;
    cols[0] = math::vector_to_theta<W>(trail[0].pos);
    for (size_t k = 1; k < n; ++k) {
      cols[k] = math::vector_to_theta<W>(trail[k].pos);
      // A geodesic edge's column sweep never exceeds W/2 (antipodal symmetry,
      // see geodesic_col_span_cols), so the short-way delta covers it
      // regardless of direction — except at ~exactly W/2, where the delta's
      // sign (which semicircle) is float noise.
      float d = cols[k] - cols[k - 1];
      if (d > W * 0.5f)
        d -= W;
      else if (d < -W * 0.5f)
        d += W;
      if (std::abs(d) >= W * 0.5f - 3.0f)
        walk_safe = false;
      // geodesic_col_span_cols refuses to bound an edge whose great-circle
      // axis is near-horizontal, and the per-edge tier then treats it as
      // visible; the endpoint columns walked here do not bound such an edge
      // either. |axis.y| = |cy| / |cross| and |cross| <= 1, so testing the
      // unnormalized cy covers every case it rejects.
      const math::Vector &ca_pos = trail[k - 1].pos;
      const math::Vector &cb_pos = trail[k].pos;
      if (std::abs(ca_pos.z * cb_pos.x - ca_pos.x * cb_pos.z) < AXIS_Y_EPS)
        walk_safe = false;
      cum += d;
      cum_lo = fminf(cum_lo, cum);
      cum_hi = fmaxf(cum_hi, cum);
#ifdef HS_PROFILE_PLOT_STALLS
      gate_batch.step();
#endif
    }
    // Near a pole the plotted column is float noise (same caution as the
    // per-edge spans), so only cull by the column arc when the whole trail
    // provably stays clear.
    if (walk_safe && sqrtf(fmaxf(0.0f, min_sp2)) - max_arc >= MIN_SIN_PHI) {
      int col_s, col_len;
      finish_col_span<W>(cols[0] + cum_lo, cum_hi - cum_lo, col_s, col_len);
      if (!ClipRegion::arcs_overlap(xc.rs, xc.length(W), col_s, col_len, W)) {
#ifdef HS_PROFILE_PLOT_STALLS
        gate_batch.step();
        gate_batch.finish();
#endif
        HS_PLOT_RENDER_COUNT(prologue_column_rejects);
        return {rows, cols, true};
      }
    }
  }
#ifdef HS_PROFILE_PLOT_STALLS
  gate_batch.step();
  gate_batch.finish();
#endif
  return {rows, cols, false};
}
