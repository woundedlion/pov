/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/scan/shapes.h"

/** @file ring_stack.h
 * @brief Same-axis distorted ring stacks. */

namespace Scan {

/**
 * @brief Fused single-pass rasterizer for a stack of same-axis distorted
 *        rings.
 */
HS_O3_BEGIN
struct DistortedRingStack {
  /** @brief Stroke-reach pad in the candidate table's bounds, absorbing
   *  fast_acos error and float slop. */
  static constexpr float REACH_PAD = 1e-3f;

  /**
   * @brief Caller-owned map from a pixel's (polar, azimuth) cell to the range
   *        of stack rings that can light it; rebuilt unless the range reports ready.
   * @tparam W Canvas width; sets the azimuth chunk count.
   * @tparam H Canvas height; sets the polar bin count.
   * @details Each cell holds the lowest and highest ring index whose stroke can
   * reach some pixel in the cell. The range is a superset: rings inside it are
   * still tested exactly, so the map only skips rings that cannot light.
   */
  template <int W, int H> struct CandidateTable {
    static constexpr int CHUNKS =
        W / 3 > 4 ? W / 3 : 4;     /**< Azimuth chunks over the full turn. */
    static constexpr int BINS = H; /**< Polar bins over [0, PI]. */
    static constexpr int MAX_RINGS = 255; /**< Ring indices fit a byte. */
    /** @brief Inclusive ring-index range; lo > hi when no ring reaches. */
    struct Range {
      uint8_t lo; ///< Lowest ring index that can reach the cell.
      uint8_t hi; ///< Highest ring index that can reach the cell.
    };
    Range cells[BINS * CHUNKS]; /**< Row-major [bin][chunk]. */
  };

  /**
   * @brief Traps unless every occupied ring is a knot ring with zero phase and
   *        slot 0's basis.
   * @param n_rings Stack size.
   * @param shapes Ring shapes indexed by slot.
   * @param slot_by_ring n_rings entries mapping ring index -> slot, -1 if
   *        culled.
   * @param n_slots Number of shapes; every occupied entry must index below it.
   * @details The fused scan derives one phase-free azimuth frame per pixel
   * from slot 0, so a ring with a non-zero phase or a divergent basis would
   * render wrong.
   */
  template <typename ShapeRange>
  HS_HOT_FLASH_MEMBER static void
  check_stack_preconditions(int n_rings, ShapeRange shapes,
                            const int8_t *slot_by_ring, int n_slots) {
    for (int i = 0; i < n_rings; ++i) {
      const int s = slot_by_ring[i];
      if (s < 0)
        continue;
      HS_CHECK(s < n_slots, "ring stack slot index out of range");
      HS_CHECK(shapes[s].knots != nullptr,
               "ring stack rings must be knot rings");
      HS_CHECK(shapes[s].phase == 0.0f, "ring stack rings must have no phase");
      HS_CHECK(math::dot(shapes[s].normal, shapes[0].normal) >=
                   1.0f - math::TOLERANCE,
               "ring stack rings must share slot 0's normal");
      HS_CHECK(math::dot(shapes[s].u, shapes[0].u) >= 1.0f - math::TOLERANCE,
               "ring stack rings must share slot 0's u axis");
      HS_CHECK(math::dot(shapes[s].w, shapes[0].w) >= 1.0f - math::TOLERANCE,
               "ring stack rings must share slot 0's w axis");
    }
  }

  /**
   * @brief Fills the candidate table for one frame's stack.
   * @param n_rings Stack size.
   * @param shapes Ring shapes indexed by slot.
   * @param slot_by_ring n_rings entries mapping ring index -> slot, -1 if
   *        culled.
   * @param table Table to fill.
   * @details A chunk's span is the knot range of the chunk and the k
   * neighbours on each side the stroke can reach at the ring band's narrowest
   * circle, widened by the thickness. A ring whose chart compresses past the
   * search budget claims its whole band.
   */
  template <int W, int H, typename ShapeRange>
  HS_HOT_FLASH_MEMBER static void
  build_candidate_table(int n_rings, ShapeRange shapes,
                        const int8_t *slot_by_ring,
                        CandidateTable<W, H> &table) {
    using Table = CandidateTable<W, H>;
    constexpr int C = Table::CHUNKS;
    constexpr float BIN_SCALE = Table::BINS / math::PI_F;
    for (auto &cell : table.cells)
      cell = typename Table::Range{255, 0};
    ScratchScope scratch(scratch_arena_b);
    float *clo = scratch_arena_b.allocate_n<float>(C);
    float *chi = scratch_arena_b.allocate_n<float>(C);
    for (int i = 0; i < n_rings; ++i) {
      const int s = slot_by_ring[i];
      if (s < 0)
        continue;
      const SDF::DistortedRing &ring = shapes[s];
      const float *kn = ring.knots;
      const int n = ring.lut_n;
      // Chunk c holds segments k with floor(k * C / n) <= c <=
      // floor((k + 1) * C / n), i.e. knots ceil(c * n / C) - 1 through
      // ceil((c + 1) * n / C), the last wrapping to knot 0.
      float gmin = 1e9f, gmax = -1e9f;
      int k_begin = 0;
      for (int c = 0; c < C; ++c) {
        const int k_next = ((c + 1) * n + C - 1) / C;
        float lo = kn[k_next == n ? 0 : k_next];
        float hi = lo;
        if constexpr (requires { shapes.knot_range(s, 0, lo, hi); })
          shapes.knot_range(s, k_next == n ? 0 : k_next, lo, hi);
        for (int k = k_begin; k < k_next; ++k) {
          float knot_lo = kn[k], knot_hi = kn[k];
          if constexpr (requires { shapes.knot_range(s, k, knot_lo, knot_hi); })
            shapes.knot_range(s, k, knot_lo, knot_hi);
          lo = __builtin_fminf(lo, knot_lo);
          hi = __builtin_fmaxf(hi, knot_hi);
        }
        clo[c] = lo;
        chi[c] = hi;
        gmin = __builtin_fminf(gmin, lo);
        gmax = __builtin_fmaxf(gmax, hi);
        k_begin = k_next - 1;
      }
      const float reach = ring.thickness + REACH_PAD;
      const float band_lo = fmaxf(0.0f, ring.target_angle + gmin - reach);
      const float band_hi = fminf(math::PI_F, ring.target_angle + gmax + reach);
      const float sin_min = fminf(sinf(band_lo), sinf(band_hi));
      bool whole = sin_min * SDF::DistortedRing::MAX_SEARCH_CELLS *
                       ring.knot_cell_angle <
                   reach;
      int k_chunks = C;
      if (!whole) {
        const float kf = reach / (math::TWO_PI_F * sin_min) * C;
        whole = kf >= C / 2;
        if (!whole)
          k_chunks = static_cast<int>(kf) + 1;
      }
      for (int c = 0; c < C; ++c) {
        float lo = gmin, hi = gmax;
        if (!whole) {
          lo = clo[c];
          hi = chi[c];
          for (int j = 1; j <= k_chunks; ++j) {
            const int cl = c - j < 0 ? c - j + C : c - j;
            const int cr = c + j >= C ? c + j - C : c + j;
            lo = fminf(lo, fminf(clo[cl], clo[cr]));
            hi = fmaxf(hi, fmaxf(chi[cl], chi[cr]));
          }
        }
        const float p0 = ring.target_angle + lo - reach;
        const float p1 = ring.target_angle + hi + reach;
        const int b0 = static_cast<int>(fmaxf(0.0f, p0) * BIN_SCALE);
        int b1 = static_cast<int>(fminf(math::PI_F, p1) * BIN_SCALE);
        if (b1 > Table::BINS - 1)
          b1 = Table::BINS - 1;
        if (b0 > b1)
          continue;
        // Rings arrive in ascending order: the first to reach a cell sets its
        // lo (255 until then), the latest its hi.
        typename Table::Range *cell = &table.cells[b0 * C + c];
        for (int b = b0; b <= b1; ++b, cell += C) {
          cell->lo = cell->lo < i ? cell->lo : static_cast<uint8_t>(i);
          cell->hi = static_cast<uint8_t>(i);
        }
      }
    }
  }

  /**
   * @brief Rasterizes every ring of a same-axis stack in one scan over the
   *        union band.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader: shader(int slot, const Vector &p,
   *         Fragment &f), with f populated as by process_pixel (v0 = azimuth
   *         t in [0, 1), v1 = raw distance, v2 = stroke coverage).
   * @tparam ShapeRange Indexable ring shapes; optional knot_range(slot, knot,
   *         lo, hi) widens knot bounds, and table_ready()/table_built() retain
   *         the candidate table across draws.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param n_rings Stack size; at most CandidateTable::MAX_RINGS.
   * @param shapes n_slots knot-mode rings sharing one Basis and zero phase, in
   *        ascending ring order; culled rings are simply absent. Their knot
   *        prefilters go unused and may be null.
   * @param slot_by_ring n_rings signed entries mapping ring index -> slot in
   *        shapes, -1 for culled rings.
   * @param n_slots Number of shapes, in [1, 127].
   * @param table Candidate map storage; rebuilt unless the range reports ready.
   * @param shader Per-ring fragment shader (see RingShaderT).
   * @details The per-pixel frame (axis dot, fast_acos, fast_atan2) is computed
   * once; the candidate table names the rings that can reach the pixel, each
   * evaluated in ascending ring index via
   * SDF::DistortedRing::distance_from_frame. At pole_lod_aggressiveness 0 the
   * output matches rasterizing the rings one by one under suppress_pole_fill,
   * to float rounding, when each plot writes one destination pixel.
   * Does not read canvas.debug(). Honours the canvas's column-time basis and
   * render-x window and ignores its birth fraction.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT,
            typename ShapeRange>
  static void draw(PipelineT &pipeline, Canvas &canvas, int n_rings,
                   ShapeRange shapes, const int8_t *slot_by_ring, int n_slots,
                   CandidateTable<W, H> &table, RingShaderT &&shader) {
    using Table = CandidateTable<W, H>;
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    HS_CHECK(n_slots >= 1, "ring stack needs at least one slot");
    HS_CHECK(n_slots <= INT8_MAX,
             "ring stack exceeds the signed slot index range");
    HS_CHECK(n_rings <= Table::MAX_RINGS,
             "ring stack exceeds the candidate table's ring index range");
    check_stack_preconditions(n_rings, shapes, slot_by_ring, n_slots);
    const Scan::ColumnBasis *time = column_basis(pipeline, canvas);
    if (time && time->azimuth) {
      draw_column_time<W, H>(pipeline, canvas, n_rings, shapes, slot_by_ring,
                             n_slots, table, shader, *time);
      return;
    }
    draw_fixed<W, H>(pipeline, canvas, n_rings, shapes, slot_by_ring, n_slots,
                     table, shader);
  }

  /**
   * @brief DistortedRingStack::draw at the fixed meridian, placed in FLASH.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @tparam ShapeRange Indexable ring shapes.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param n_rings Stack size.
   * @param shapes Ring shapes indexed by slot.
   * @param slot_by_ring Ring index to slot map.
   * @param n_slots Number of shapes.
   * @param table Candidate map storage.
   * @param shader Per-ring fragment shader.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT,
            typename ShapeRange>
  HS_HOT_FLASH_MEMBER static void
  draw_fixed(PipelineT &pipeline, Canvas &canvas, int n_rings,
             ShapeRange shapes, const int8_t *slot_by_ring, int n_slots,
             CandidateTable<W, H> &table, RingShaderT &shader) {
    draw_body<W, H, false>(pipeline, canvas, n_rings, shapes, slot_by_ring,
                           n_slots, table, shader, nullptr);
  }

  /**
   * @brief DistortedRingStack::draw along a rotated basis.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @tparam ShapeRange Indexable ring shapes.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param n_rings Stack size.
   * @param shapes Ring shapes indexed by slot.
   * @param slot_by_ring Ring index to slot map.
   * @param n_slots Number of shapes.
   * @param table Candidate map storage.
   * @param shader Per-ring fragment shader.
   * @param time Column-time state with a rotated basis.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT,
            typename ShapeRange>
  HS_NOINLINE_NOCLONE static void
  draw_column_time(PipelineT &pipeline, Canvas &canvas, int n_rings,
                   ShapeRange shapes, const int8_t *slot_by_ring, int n_slots,
                   CandidateTable<W, H> &table, RingShaderT &shader,
                   const Scan::ColumnBasis &time) {
    draw_body<W, H, true>(pipeline, canvas, n_rings, shapes, slot_by_ring,
                          n_slots, table, shader, &time);
  }

  /**
   * @brief Fused stack scan shared by the fixed and rotated-basis draws.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam COLUMN_TIME Whether pixels follow @p time's rotated basis.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @tparam ShapeRange Indexable ring shapes.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param n_rings Stack size.
   * @param shapes Ring shapes indexed by slot.
   * @param slot_by_ring Ring index to slot map.
   * @param n_slots Number of shapes.
   * @param table Candidate map storage.
   * @param shader Per-ring fragment shader.
   * @param time Column-time state; read only with COLUMN_TIME.
   */
  template <int W, int H, bool COLUMN_TIME, typename PipelineT,
            typename RingShaderT, typename ShapeRange>
  HS_HOT_INLINE static void
  draw_body(PipelineT &pipeline, Canvas &canvas, int n_rings, ShapeRange shapes,
            const int8_t *slot_by_ring, int n_slots,
            CandidateTable<W, H> &table, RingShaderT &shader,
            const Scan::ColumnBasis *time) {
    using Table = CandidateTable<W, H>;
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    [[maybe_unused]] const math::Vector *azimuth = nullptr;
    [[maybe_unused]] const math::Vector *up = nullptr;
    float sweep = 0.0f;
    if constexpr (COLUMN_TIME) {
      azimuth = time->azimuth;
      up = time->up;
      sweep = time->sweep_arc;
    }
    int y_lo = H, y_hi = -1;
    for (int s = 0; s < n_slots; ++s) {
      auto b =
          swept_rows<H>(shapes[s].template get_vertical_bounds<H>(), sweep);
      y_lo = std::min(y_lo, b.y_min);
      y_hi = std::max(y_hi, b.y_max);
    }
    const auto &cr = source_clip<W, H>(pipeline, canvas);
    const auto xc = cr.x_clip();
    const WindowPieces pieces = window_pieces<W>(xc, canvas.time_window());
    y_lo = std::max(y_lo, cr.render_y_start());
    y_hi = std::min(y_hi, cr.render_y_end() - 1);
    if (y_lo > y_hi)
      return;

    {
      HS_PROFILE(ring_stack_table);
      bool ready = false;
      if constexpr (requires { shapes.table_ready(); })
        ready = shapes.table_ready();
      if (!ready) {
        build_candidate_table<W, H>(n_rings, shapes, slot_by_ring, table);
        if constexpr (requires { shapes.table_built(); })
          shapes.table_built();
      }
    }
    constexpr float BIN_SCALE = Table::BINS / math::PI_F;

    // Aliased exact-pole rows are skipped unless the axis is near a canvas pole
    // (r_val below the projection floor).
    SDF::AxisProjection ap = SDF::project_axis(shapes[0].normal);
    const bool skip_pole_rows = ap.r_val >= SDF::MIN_HORIZONTAL_PROJ;

    const math::Vector axis_v = shapes[0].normal;
    const math::Vector axis_u = shapes[0].u;
    const math::Vector axis_w = shapes[0].w;

    for (int y = y_lo; y <= y_hi; ++y) {
      const float sp = math::TrigLUT<W, H>::sin_phi[y];
      const float cp = math::TrigLUT<W, H>::cos_phi[y];
      if (skip_pole_rows && std::abs(ap.r_val * sp) < SDF::INTERVAL_DENOM_EPS)
        continue;
      walk_pieces_once(pieces, [&](int x) __attribute__((always_inline)) {
        math::Vector p;
        if constexpr (COLUMN_TIME)
          p = rotated_direction<W>(azimuth, up, x, sp, cp);
        else
          p = math::Vector(sp * cos_theta[x], cp, sp * sin_theta[x]);
        const float d = math::dot(p, axis_v);
        const float polar = math::fast_acos(hs::clamp(d, -1.0f, 1.0f));
        const float point_azimuth = SDF::basis_azimuth(p, axis_u, axis_w, 0.0f);
        // Point azimuth is in [0, 2 PI], so wrap_t folds 1 onto 0.
        float t_norm = point_azimuth / (2 * math::PI_F);
        t_norm = t_norm >= 1.0f ? 0.0f : t_norm;
        int bin = static_cast<int>(polar * BIN_SCALE);
        if (bin > Table::BINS - 1)
          bin = Table::BINS - 1;
        int chunk = static_cast<int>(t_norm * Table::CHUNKS);
        if (chunk > Table::CHUNKS - 1)
          chunk = Table::CHUNKS - 1;
        const typename Table::Range cell =
            table.cells[bin * Table::CHUNKS + chunk];
        if (cell.lo > cell.hi)
          return;
        const float sin_polar =
            sqrtf(fmaxf(1.0f - d * d, SDF::DistortedRing::POLE_SIN2_FLOOR));
        for (int i = cell.lo; i <= cell.hi; ++i) {
          const int s = slot_by_ring[i];
          if (s < 0)
            continue;
          SDF::DistanceResult res;
          shapes[s].distance_from_frame(d, polar, sin_polar, t_norm, res);
          const float dd = res.dist;
          if (dd >= 0.0f)
            continue;
          // process_pixel's stroke epilogue with a slot-aware shader.
          const float aa = res.size;
          const float alpha = aa > 0.0f ? math::quintic_kernel(-dd / aa) : 0.0f;
          if (alpha <= MIN_ALPHA)
            continue;
          Fragment frag;
          frag.pos = p;
          frag.v0 = res.t;
          frag.v1 = res.raw_dist;
          frag.v2 = alpha;
          frag.v3 = res.aux;
          frag.size = res.size;
          frag.age = 0;
          shader(s, p, frag);
          if (frag.color.alpha <= MIN_ALPHA)
            continue;
          // The walk visits only clip-admitted columns of render rows.
          const float a = frag.color.alpha * alpha;
          plot_admitted(pipeline, canvas, x, y, frag.color.color, frag.age, a);
        }
      });
    }
  }
};
HS_O3_END

} // namespace Scan
