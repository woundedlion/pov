/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>
#include "core/render/sdf.h"
#include "core/render/shading.h"
#include "core/color/color.h"
#include "core/render/filter/pipeline.h"
#include "core/containers/static_circular_buffer.h"
#include "core/render/canvas.h"
#include "core/platform/platform.h"
#include "core/render/scan/raster.h"

/**
 * @file shapes.h
 * @brief The SDF-backed scan primitives: rings, polygons, circles, points, stars and flowers.
 */

namespace Scan {

/** @brief Optional DistortedRing draw settings; defaults draw a plain ring. */
struct RingDrawOptions {
  /** Angular phase offset in radians. */
  float phase = 0.0f;
  /** Render the bounding box for debugging. */
  bool debug_bb = false;
  /** Drop the degenerate exact-pole row instead of full-row filling it. */
  bool suppress_pole_fill = false;
};

/**
 * @brief Draws a ring whose radius is modulated around the circumference by
 *        shift_fn.
 */
struct DistortedRing {
  /**
   * @brief Rasterizes an undisplaced ring with exact polar centerline distance.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param options Phase, debug bounding box and pole-fill settings.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw_flat(PipelineRef pipeline, Canvas &canvas,
                        const math::Basis &basis, float radius, float thickness,
                        FragmentShaderFn fragment_shader,
                        const RingDrawOptions &options = {}) {
    SDF::FlatDistortedRing shape(basis, radius, thickness, options.phase);
    shape.suppress_pole_fill = options.suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      options.debug_bb);
  }

  /**
   * @brief Rasterizes a circumference-modulated ring.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param shift_fn Scalar modulation function over the circumference.
   * @param amplitude Modulation amplitude (radians); must upper-bound
   *        |shift_fn|.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Basis &basis, float radius, float thickness,
                   ScalarFn shift_fn, float amplitude,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    SDF::DistortedRing shape(basis, radius, thickness, shift_fn, amplitude,
                             phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a ring whose centerline is a shift-knot polyline with
   *        exact stroke distance (see SDF::DistortedRing's knot overload).
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param knots lut_n centerline shifts; closure wraps to entry 0.
   *        Storage must outlive the call.
   * @param lut_n Number of knot cells; at least 3.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param options Phase, debug bounding box and pole-fill settings.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void
  draw(PipelineRef pipeline, Canvas &canvas, const math::Basis &basis,
       float radius, float thickness, const float *knots, int lut_n,
       FragmentShaderFn fragment_shader, const RingDrawOptions &options = {}) {
    SDF::KnotPrefilter prefilter;
    SDF::DistortedRing shape(basis, radius, thickness, knots, lut_n,
                             options.phase, prefilter);
    shape.suppress_pole_fill = options.suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      options.debug_bb);
  }
};

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
   *        of stack rings that can light it, rebuilt by every draw.
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
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param n_rings Stack size; at most CandidateTable::MAX_RINGS.
   * @param shapes n_slots knot-mode rings sharing one Basis and zero phase, in
   *        ascending ring order; culled rings are simply absent. Their knot
   *        prefilters go unused and may be null.
   * @param slot_by_ring n_rings signed entries mapping ring index -> slot in
   *        shapes, -1 for culled rings.
   * @param n_slots Number of shapes, in [1, 127].
   * @param table Candidate map storage, rebuilt here.
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
    const Temporal::ColumnTime *time = column_time(pipeline, canvas);
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
                   const Temporal::ColumnTime &time) {
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
            const Temporal::ColumnTime *time) {
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
      sweep = time->sweep;
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
    const WindowPieces pieces = window_pieces<W>(xc, canvas.render_x_window());
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
        const float azimuth = SDF::basis_azimuth(p, axis_u, axis_w, 0.0f);
        // azimuth is in [0, 2 PI], so wrap_t reduces to folding 1 onto 0.
        float t_norm = azimuth / (2 * math::PI_F);
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
          if constexpr (requires {
                          pipeline.plot_in_bounds(
                              canvas, x, y, frag.color.color, frag.age, a);
                        })
            pipeline.plot_in_bounds(canvas, x, y, frag.color.color, frag.age,
                                    a);
          else
            pipeline.plot(canvas, x, y, frag.color.color, frag.age, a);
        }
      });
    }
  }
};
HS_O3_END

/**
 * @brief Draws a flat (tangent-plane) regular polygon projected onto the
 *        sphere.
 */
struct PlanarPolygon {
  /**
   * @brief Rasterizes a tangent-plane regular polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon plane.
   * @param radius Polygon circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Basis &basis, float radius, int sides,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::PlanarPolygon shape(res.first, res.second, sides, phase,
                             radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color tangent-plane regular polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon plane.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::PlanarPolygon shape(res.first, res.second, sides, phase,
                             radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a great-circle line segment of given thickness between two
 *        vectors.
 */
struct Line {
  /**
   * @brief Rasterizes a great-circle segment between two unit vectors.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param v1 First endpoint as a world-space unit vector.
   * @param v2 Second endpoint as a world-space unit vector.
   * @param thickness Stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const math::Vector &v1,
                   const math::Vector &v2, float thickness,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    SDF::Line shape(v1, v2, thickness);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

/**
 * @brief Draws a ring stroke using SDF rasterization.
 */
struct Ring {
  /**
   * @brief Draws a ring stroke from an orientation basis.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere, in [0, 2].
   * @param thickness Ring stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   * @details A radius past 1 keeps the caller's azimuth origin and handedness.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Basis &basis, float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    SDF::Ring shape(basis, radius, thickness, phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Draws a ring stroke from a plane-normal vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param normal Plane normal as a world-space vector.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Vector &normal, float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    math::Basis basis = math::make_basis(math::Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, thickness,
                           fragment_shader, phase, debug_bb);
  }
};

/**
 * @brief Fused single-pass rasterizer for a small group of rings.
 */
HS_O3_BEGIN
struct RingGroup {
  /** Maximum number of rings in one group. */
  static constexpr int MAX_RINGS = 8;

  /** Inclusive row spans where one column meets the covering band. */
  struct ColumnRows {
    int16_t lo[2]{}; /**< First rows. */
    int16_t hi[2]{}; /**< Last rows. */
    int count = 0;   /**< Number of disjoint spans. */
  };

  /** Angular solutions shared by the two arms at one sweep fraction. */
  struct ColumnBand {
    float center; /**< First-arm sinusoid phase. */
    float inner;  /**< Inner angular offset; greater than outer when empty. */
    float outer;  /**< Outer angular offset. */
  };

  /**
   * @brief Prepares a column's angular band intersections.
   * @param a First-arm coefficient of sin(phi).
   * @param b Coefficient of cos(phi).
   * @param low Lower cosine bound.
   * @param high Upper cosine bound.
   * @return Angular solutions about the first-arm phase.
   */
  HS_HOT_INLINE static ColumnBand column_band(float a, float b, float low,
                                              float high) {
    const float radius = sqrtf(a * a + b * b);
    if (low > radius || high < -radius)
      return {0.0f, 1.0f, 0.0f};
    if (radius < 1e-7f)
      return {0.0f, 0.0f, math::PI_F};
    const float inverse = 1.0f / radius;
    return {math::fast_atan2(a * inverse, b * inverse),
            math::fast_acos(hs::clamp(high * inverse, -1.0f, 1.0f)),
            math::fast_acos(hs::clamp(low * inverse, -1.0f, 1.0f))};
  }

  /**
   * @brief Fills the covering band's shared angular solutions once per draw.
   * @param bands Half-width table receiving angular solutions.
   * @param terms Per-column, per-slot axis terms.
   * @param half Number of distinct sweep fractions.
   * @param n Number of slots.
   * @param mid Covering slot index.
   * @param low Lower cosine bound.
   * @param high Upper cosine bound.
   */
  HS_HOT_FLASH_MEMBER static void fill_column_bands(ColumnBand *bands,
                                                    const float *terms,
                                                    int half, int n, int mid,
                                                    float low, float high) {
    for (int c = 0; c < half; ++c) {
      const float *axis = terms + 2 * (n * c + mid);
      bands[c] = column_band(axis[0], axis[1], low, high);
    }
  }

  /**
   * @brief Solves a column's sinusoidal axis dot against a cosine band.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @param a First-arm coefficient of sin(phi).
   * @param b Coefficient of cos(phi).
   * @param low Lower cosine bound.
   * @param high Upper cosine bound.
   * @param band Shared angular solutions for this sweep fraction.
   * @param far_arm Whether to reflect the first-arm phase.
   * @param y_lo First row the caller reads.
   * @param y_hi Last row the caller reads.
   * @return Up to two disjoint row spans with in-band endpoints. Within
   *         [@p y_lo, @p y_hi] their union is the band's whole-canvas spans';
   *         spans that cannot reach that window are dropped.
   * @pre The canvas trigonometric table is initialized.
   */
  template <int W, int H>
  HS_HOT_FLASH_MEMBER static ColumnRows
  column_rows(float a, float b, float low, float high, const ColumnBand &band,
              bool far_arm, int y_lo = 0, int y_hi = H - 1) {
    ColumnRows out;
    if (band.inner > band.outer)
      return out;
    const float center = far_arm ? -band.center : band.center;
    const float inner = band.inner, outer = band.outer;
    a = far_arm ? -a : a;
    constexpr float ANGLE_PAD = 4e-3f;
    const auto inside = [&](int y) HS_HOT_INLINE {
      const float d = a * math::TrigLUT<W, H>::sin_phi[y] +
                      b * math::TrigLUT<W, H>::cos_phi[y];
      return d >= low && d <= high;
    };
    SDF::Bounds candidates[6];
    int count = 0;
    for (int turn = -1; turn <= 1; ++turn) {
      const float origin = center + turn * math::TWO_PI_F;
      for (int side : {-1, 1}) {
        const float start = origin + (side < 0 ? -outer : inner) - ANGLE_PAD;
        const float end = origin + (side < 0 ? -inner : outer) + ANGLE_PAD;
        if (end < 0.0f || start > math::PI_F)
          continue;
        auto rows = SDF::phi_bounds_to_rows<H>(
            std::max(0.0f, start - 1e-5f), std::min(math::PI_F, end + 1e-5f));
        if (rows.y_max < y_lo || rows.y_min > y_hi)
          continue;
        while (rows.y_min <= rows.y_max && !inside(rows.y_min))
          ++rows.y_min;
        while (rows.y_min <= rows.y_max && !inside(rows.y_max))
          --rows.y_max;
        if (rows.y_min > rows.y_max)
          continue;
        int at = count++;
        while (at > 0 && candidates[at - 1].y_min > rows.y_min) {
          candidates[at] = candidates[at - 1];
          --at;
        }
        candidates[at] = rows;
      }
    }
    for (int i = 0; i < count; ++i) {
      const auto rows = candidates[i];
      if (out.count && rows.y_min <= out.hi[out.count - 1] + 1) {
        out.hi[out.count - 1] =
            std::max<int>(out.hi[out.count - 1], rows.y_max);
      } else {
        out.lo[out.count] = rows.y_min;
        out.hi[out.count++] = rows.y_max;
      }
    }
    return out;
  }

  /**
   * @brief Rasterizes every ring of a group in one scan over the union band.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader: shader(int slot, const Vector &p,
   *         Fragment &f). One Fragment serves the whole scan and only color,
   *         pos, v2 (stroke coverage), size and age are refreshed per pixel —
   *         no UVs, no raw distance. v0, v1 and v3 are NaN in debug builds and
   *         retain defaults or previous values in release builds. The shader
   *         must not read inputs it did not set.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes, in [1, MAX_RINGS].
   * @param shader Per-ring fragment shader (see RingShaderT).
   * @param debug_bb When true, or under canvas.debug(), falls back to per-ring
   *        rasterizes, which scan each ring's own row intervals and fill
   *        v0/v1/v3 per pixel.
   * @param previous_normals Optional per-slot axes at the previous frame
   *        endpoint. With temporal sampling on, column `x` evaluates each
   *        slot about the axis interpolated linearly from its previous axis
   *        to its current one at `sweep_fraction(x, W)` and normalised, and
   *        each slot's bounds grow by the angle between the two.
   * @details Row intervals come from one covering ring: the middle member (n/2)
   * inflated by the group's maximum plane/radius deviation plus thickness.
   * Members evaluate per pixel in ascending slot order. The canvas's
   * column-time basis is honoured and its birth fraction ignored. At
   * pole_lod_aggressiveness 0 and one destination pixel per plot, the only
   * divergence from rasterizing the rings one by one is AA-tail pixels a
   * member's own interval clip would drop.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  static void draw(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes,
                   int n, RingShaderT &&shader, bool debug_bb = false,
                   const math::Vector *previous_normals = nullptr) {
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    HS_CHECK(n >= 1 && n <= MAX_RINGS,
             "ring group size must be in [1, MAX_RINGS]");
    if (debug_bb || canvas.debug()) {
      for (int s = 0; s < n; ++s) {
        auto slot_shader = [&](const math::Vector &p, Fragment &f) {
          shader(s, p, f);
        };
        if constexpr (std::is_constructible_v<PipelineRef, PipelineT &>) {
          PipelineRef ref(pipeline);
          Scan::rasterize<W, H>(ref, canvas, shapes[s], slot_shader, true);
        } else {
          Scan::rasterize<W, H, true>(pipeline, canvas, shapes[s], slot_shader,
                                      true);
        }
      }
      return;
    }

    const Temporal::ColumnTime *time = column_time(pipeline, canvas);
    const bool rotated = time && time->azimuth;
    const bool moving_slots =
        previous_normals && Temporal::temporal_sampling_enabled();
    if (rotated || moving_slots) {
      draw_column_time<W, H>(pipeline, canvas, shapes, n, shader,
                             rotated ? time : nullptr,
                             moving_slots ? previous_normals : nullptr);
      return;
    }
    draw_fixed<W, H>(pipeline, canvas, shapes, n, shader);
  }

  /**
   * @brief RingGroup::draw at the fixed meridian with still slots, placed in
   *        FLASH.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param shader Per-ring fragment shader.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  HS_HOT_FLASH_MEMBER static void
  draw_fixed(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes,
             int n, RingShaderT &shader) {
    draw_body<W, H, false, true>(pipeline, canvas, shapes, n, shader, nullptr,
                                 nullptr);
  }

  /**
   * @brief RingGroup::draw under a rotated basis or moving slot axes.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param shader Per-ring fragment shader.
   * @param time Column-time state with a rotated basis, or null.
   * @param previous_normals Per-slot axes at the previous endpoint, or null.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  HS_NOINLINE_NOCLONE static void
  draw_column_time(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes,
                   int n, RingShaderT &shader, const Temporal::ColumnTime *time,
                   const math::Vector *previous_normals) {
    draw_body<W, H, true, false>(pipeline, canvas, shapes, n, shader, time,
                                 previous_normals);
  }

  /**
   * @brief Fills each column's ring-axis terms for the column-time scan.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param terms 2 * n * W / 2 entries receiving, per first-arm column and
   *        slot, the slot's axis at the column's sweep time dotted with the
   *        column's azimuth and with its up direction.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param azimuth Rotated first-arm azimuths, or null for the fixed meridian.
   * @param up Rotated up directions; read only with @p azimuth.
   * @param previous_normals Per-slot axes at the previous endpoint, or null.
   * @details A pixel at colatitude phi in column x has axis dot
   * `+-sin(phi) * a + cos(phi) * b`, negated on the far arm.
   */
  template <int W, int H>
  HS_NOINLINE_NOCLONE static void
  fill_column_axes(float *terms, const SDF::Ring *shapes, int n,
                   const math::Vector *azimuth, const math::Vector *up,
                   const math::Vector *previous_normals) {
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    // Slot s at fraction u has axis (start + u * motion) / |start + u * motion|.
    math::Vector start[MAX_RINGS], motion[MAX_RINGS];
    float length_terms[MAX_RINGS][3];
    for (int s = 0; s < n; ++s) {
      start[s] = previous_normals ? previous_normals[s] : shapes[s].normal;
      motion[s] = shapes[s].normal - start[s];
      length_terms[s][0] = math::dot(start[s], start[s]);
      length_terms[s][1] = 2.0f * math::dot(start[s], motion[s]);
      length_terms[s][2] = math::dot(motion[s], motion[s]);
    }
    for (int c = 0; c < W / 2; ++c) {
      const math::Vector e =
          azimuth ? azimuth[c] : math::Vector(cos_theta[c], 0.0f, sin_theta[c]);
      const math::Vector v = azimuth ? up[c] : math::Vector(0.0f, 1.0f, 0.0f);
      const float u = Temporal::sweep_fraction(c, W);
      float *out = terms + 2 * n * c;
      for (int s = 0; s < n; ++s) {
        const float length_sq =
            length_terms[s][0] +
            u * (length_terms[s][1] + u * length_terms[s][2]);
        float inverse = math::fast_rsqrt(length_sq);
        // a third Newton step: the ring's cosine band needs full precision
        inverse *= 1.5f - 0.5f * length_sq * inverse * inverse;
        out[2 * s] =
            (math::dot(e, start[s]) + u * math::dot(e, motion[s])) * inverse;
        out[2 * s + 1] =
            (math::dot(v, start[s]) + u * math::dot(v, motion[s])) * inverse;
      }
    }
  }

  /**
   * @brief Fused scan body shared by the fixed and column-time draws.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam COLUMN_TIME Whether @p time and @p previous_normals are read.
   * @tparam INLINE_WALK Whether the row walk inlines into the body, so a
   *         FLASH variant takes it along.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param shader Per-ring fragment shader.
   * @param time Column-time state with a rotated basis, or null.
   * @param previous_normals Per-slot axes at the previous endpoint, or null.
   */
  template <int W, int H, bool COLUMN_TIME, bool INLINE_WALK,
            typename PipelineT, typename RingShaderT>
  HS_HOT_INLINE static void
  draw_body(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes, int n,
            RingShaderT &shader, const Temporal::ColumnTime *time,
            const math::Vector *previous_normals) {
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();

    [[maybe_unused]] const math::Vector *azimuth = nullptr;
    [[maybe_unused]] const math::Vector *up = nullptr;
    [[maybe_unused]] float sweep = 0.0f;
    [[maybe_unused]] float slot_sweep[MAX_RINGS] = {};
    [[maybe_unused]] float max_slot_sweep = 0.0f;
    ScratchScope scratch(scratch_arena_b);
    [[maybe_unused]] float *axis_terms = nullptr;
    [[maybe_unused]] ColumnBand *column_bands = nullptr;
    if constexpr (COLUMN_TIME) {
      if (time) {
        azimuth = time->azimuth;
        up = time->up;
        sweep = time->sweep;
      }
      if (previous_normals) {
        for (int s = 0; s < n; ++s) {
          slot_sweep[s] = math::fast_acos(hs::clamp(
                              math::dot(shapes[s].normal, previous_normals[s]),
                              -1.0f, 1.0f)) +
                          1e-3f;
          max_slot_sweep = std::max(max_slot_sweep, slot_sweep[s]);
        }
      }
    }
    int sy_min[MAX_RINGS], sy_max[MAX_RINGS];
    int y_lo = H, y_hi = -1;
    for (int s = 0; s < n; ++s) {
      auto b = shapes[s].template get_vertical_bounds<H>();
      if constexpr (COLUMN_TIME)
        b = swept_rows<H>(b, sweep + slot_sweep[s]);
      sy_min[s] = b.y_min;
      sy_max[s] = b.y_max;
      y_lo = std::min(y_lo, b.y_min);
      y_hi = std::max(y_hi, b.y_max);
    }
    const auto &cr = source_clip<W, H>(pipeline, canvas);
    y_lo = std::max(y_lo, cr.render_y_start());
    y_hi = std::min(y_hi, cr.render_y_end() - 1);
    if (y_lo > y_hi)
      return;

    // Covering ring: the middle member inflated by the worst centerline
    // deviation plus thickness. The 1e-3 absorbs fast_acos error.
    const int mid = n / 2;
    float pad_th = shapes[mid].thickness;
    for (int s = 0; s < n; ++s) {
      if (s == mid)
        continue;
      float dev =
          math::fast_acos(hs::clamp(
              math::dot(shapes[mid].normal, shapes[s].normal), -1.0f, 1.0f)) +
          std::abs(shapes[s].target_angle - shapes[mid].target_angle) + 1e-3f;
      pad_th = std::max(pad_th, shapes[s].thickness + dev);
    }
    SDF::Ring cover(
        math::Basis{shapes[mid].u, shapes[mid].normal, shapes[mid].w},
        shapes[mid].radius, pad_th);
    // A column's lit pixels lie in the middle member's band, about its axis at
    // the column's time, widened by every member's deviation from it at
    // either endpoint; the lerped axes' deviations stay within the larger,
    // up to their normalisation.
    [[maybe_unused]] float cover_cos_min = -2.0f, cover_cos_max = 2.0f;
    if constexpr (COLUMN_TIME) {
      float pad = pad_th;
      if (previous_normals) {
        for (int s = 0; s < n; ++s) {
          if (s == mid)
            continue;
          const float dev =
              math::fast_acos(hs::clamp(
                  math::dot(previous_normals[mid], previous_normals[s]), -1.0f,
                  1.0f)) +
              std::abs(shapes[s].target_angle - shapes[mid].target_angle) +
              1e-3f;
          pad = std::max(pad, shapes[s].thickness + dev);
        }
        pad = pad * (1.0f + max_slot_sweep * max_slot_sweep) + 1e-3f;
      }
      SDF::band_edge_cosines(shapes[mid].target_angle - pad,
                             shapes[mid].target_angle + pad, cover_cos_max,
                             cover_cos_min);
    }
    if constexpr (COLUMN_TIME) {
      axis_terms = scratch_arena_b.allocate_n<float>(2 * n * (W / 2));
      fill_column_axes<W, H>(axis_terms, shapes, n, azimuth, up,
                             previous_normals);
      column_bands = scratch_arena_b.allocate_n<ColumnBand>(W / 2);
      fill_column_bands(column_bands, axis_terms, W / 2, n, mid, cover_cos_min,
                        cover_cos_max);
    }

    Fragment frag;

    // Walk the covering ring's row runs (shared emit_row_runs).
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    const auto xc = cr.x_clip();
    StaticCircularBuffer<SDF::Interval, 4> intervals;
    StaticCircularBuffer<SDF::Interval, 8> norm;
    static_assert(decltype(intervals)::CAPACITY >=
                      SDF::sdf_max_spans<SDF::Ring>::value,
                  "intervals must hold the covering Ring's per-row emission");
    static_assert(decltype(norm)::CAPACITY == 2 * decltype(intervals)::CAPACITY,
                  "norm must hold 2 spans per input interval (seam split)");
    int active[MAX_RINGS];
    int n_active = 0;

    auto pixel_run_body = [&](int x1, int x2, int y, float sp,
                              float cp) HS_HOT_INLINE {
      for (int x = x1; x < x2; ++x) {
        [[maybe_unused]] float s_phi = sp;
        [[maybe_unused]] const float *terms = nullptr;
        math::Vector p;
        if constexpr (COLUMN_TIME) {
          const bool far_arm = x >= W / 2;
          s_phi = far_arm ? -sp : sp;
          terms = axis_terms + 2 * n * (far_arm ? x - W / 2 : x);
          const float d_mid = s_phi * terms[2 * mid] + cp * terms[2 * mid + 1];
          if (d_mid < cover_cos_min || d_mid > cover_cos_max)
            continue;
        } else {
          p = math::Vector(sp * cos_theta[x], cp, sp * sin_theta[x]);
        }
        for (int i = 0; i < (COLUMN_TIME ? n : n_active); ++i) {
          const int s = COLUMN_TIME ? i : active[i];
          if constexpr (COLUMN_TIME)
            if (y < sy_min[s] || y > sy_max[s])
              continue;
          float d;
          if constexpr (COLUMN_TIME)
            d = s_phi * terms[2 * s] + cp * terms[2 * s + 1];
          else
            d = math::dot(p, shapes[s].normal);
          const float alpha = shapes[s].stroke_alpha(d);
          if (alpha <= MIN_ALPHA)
            continue;
          if constexpr (COLUMN_TIME)
            p = column_direction<W>(azimuth, up, x, sp, cp, cos_theta,
                                    sin_theta);
#ifndef NDEBUG
          frag.poison_inputs();
#endif
          frag.color = Color4(0, 0, 0, 0);
          frag.pos = p;
          frag.v2 = alpha;
          frag.size = shapes[s].thickness;
          frag.age = 0;
          shader(s, p, frag);
          if (frag.color.alpha > MIN_ALPHA)
            pipeline.plot(canvas, x, y, frag.color.color, frag.age,
                          frag.color.alpha * alpha);
        }
      }
    };
    // One out-of-line copy of the run body serves every emit site.
    [[maybe_unused]] auto pixel_run =
        [&](int x1, int x2, int y, float sp, float cp)
            HS_NOINLINE_NOCLONE { pixel_run_body(x1, x2, y, sp, cp); };
    if constexpr (COLUMN_TIME) {
      walk_clip_columns<W>(xc, [&](int x) HS_HOT_INLINE {
        const bool far_arm = x >= W / 2;
        const int c = far_arm ? x - W / 2 : x;
        const float *terms = axis_terms + 2 * n * c;
        const auto rows = column_rows<W, H>(
            terms[2 * mid], terms[2 * mid + 1], cover_cos_min, cover_cos_max,
            column_bands[c], far_arm, y_lo, y_hi);
        for (int span = 0; span < rows.count; ++span) {
          const int end = std::min<int>(y_hi, rows.hi[span]);
          for (int y = std::max<int>(y_lo, rows.lo[span]); y <= end; ++y)
            pixel_run_body(x, x + 1, y, math::TrigLUT<W, H>::sin_phi[y],
                           math::TrigLUT<W, H>::cos_phi[y]);
        }
      });
      return;
    }
    for (int y = y_lo; y <= y_hi; ++y) {
      // Per-slot domains match solo rasterization; slot order preserves blending.
      n_active = 0;
      for (int s = 0; s < n; ++s)
        if (y >= sy_min[s] && y <= sy_max[s])
          active[n_active++] = s;
      if (n_active == 0)
        continue;
      float sp = math::TrigLUT<W, H>::sin_phi[y];
      float cp = math::TrigLUT<W, H>::cos_phi[y];
      auto walk_row = [&](auto &emit) HS_HOT_INLINE {
        if (cover.needs_full_row_scan(sp)) {
          clip_run(0, W, xc, emit);
          return;
        }
        intervals.clear();
        cover.template get_horizontal_intervals<W, H>(
            y,
            [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });
        // The full-row case is taken above, so the covering ring always
        // answers.
        emit_row_runs<W>(true, intervals, norm, xc, emit);
      };
      if constexpr (INLINE_WALK) {
        auto emit = [&](int x1, int x2)
                        HS_HOT_INLINE { pixel_run_body(x1, x2, y, sp, cp); };
        walk_row(emit);
      } else {
        auto emit = [&](int x1, int x2) { pixel_run(x1, x2, y, sp, cp); };
        walk_row(emit);
      }
    }
  }
};
HS_O3_END

/**
 * @brief Draws a disc as a zero-radius ring whose stroke spans it.
 * @details Stroke coverage falls off quintically from 1 at the center to 0 at
 * the rim; it modulates the plotted alpha and reaches the shader as the
 * fragment's register 2, so the shader picks the radial style on top of it.
 */
struct Circle {
  /**
   * @brief Draws a disc from an orientation basis.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the circle plane.
   * @param radius Circle radius as a fraction of the hemisphere.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Basis &basis, float radius,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    float th = radius * (math::PI_F / 2.0f);
    Ring::draw<W, H, ComputeUVs>(pipeline, canvas, basis, 0.0f, th,
                                 fragment_shader, 0, debug_bb);
  }

  /**
   * @brief Draws a disc from a plane-normal vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param normal Plane normal as a world-space vector.
   * @param radius Circle radius as a fraction of the hemisphere.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Vector &normal, float radius,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    math::Basis basis = math::make_basis(math::Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, fragment_shader,
                           debug_bb);
  }
};

/**
 * @brief Draws a dot as a zero-radius ring of the given thickness.
 * @details Stroke coverage falls off quintically from 1 at the center to 0 at
 * the rim, giving the soft glow effects shade through the fragment's
 * register 2.
 */
struct Point {
  /**
   * @brief Draws a dot centered on a unit vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param p Point center as a world-space unit vector.
   * @param thickness Point radius as a stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const math::Vector &p,
                   float thickness, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    math::Basis basis = math::make_basis(math::Quaternion(), p);
    Ring::draw<W, H>(pipeline, canvas, basis, 0.0f, thickness, fragment_shader,
                     0.0f, debug_bb);
  }

  /**
   * @brief Draws a dot whose shader reads its column's sweep time.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param p Point center as a world-space unit vector.
   * @param thickness Point radius as a stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel with its column's
   *        sweep fraction.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const math::Vector &p,
                   float thickness, TimedShaderFn fragment_shader,
                   bool debug_bb = false) {
    const SDF::Ring shape(math::make_basis(math::Quaternion(), p), 0.0f,
                          thickness, 0.0f);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

/**
 * @brief Draws a solid star shape.
 */
struct Star {
  /**
   * @brief Rasterizes a solid star.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the star plane.
   * @param radius Star circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of star points.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Basis &basis, float radius, int sides,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Star shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color solid star.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the star plane.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of star points.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Star shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a solid flower shape.
 */
struct Flower {
  /**
   * @brief Rasterizes a solid flower.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the flower plane.
   * @param radius Flower circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of flower petals.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Basis &basis, float radius, int sides,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color solid flower.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the flower plane.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of flower petals.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a solid spherical polygon.
 * @details Both entry points add half a sector to the caller's phase, so phase
 * 0 puts a vertex on the basis u-axis.
 */
struct SphericalPolygon {
  /**
   * @brief Rasterizes a solid spherical polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon.
   * @param radius Polygon circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const math::Basis &basis, float radius, int sides,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    float offset = math::PI_F / sides;

    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color solid spherical polygon.
   * @tparam SineDistance Use edge-plane distance for the AA band.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool SineDistance = false, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    float offset = math::PI_F / sides;
    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                radius > 1.0f);
    Scan::rasterize_solid<W, H, SineDistance>(pipeline, canvas, shape, color,
                                              debug_bb);
  }
};

} // namespace Scan
