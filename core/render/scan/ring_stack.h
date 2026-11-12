/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/scan/shapes.h"
#include "core/render/sdf/ring_knots.h"
#include <bit>
#include <array>

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

/** @brief Packed second-window knot shifts.
 * @tparam KNOTS Maximum knots per ring. */
template <int KNOTS> class KnotAlternates {
public:
  static constexpr int WORDS = (KNOTS + 31) / 32; ///< Mask words per ring.
  float *values = nullptr;     ///< Packed late-window shifts.
  uint32_t *masks = nullptr;   ///< Per-ring alternate masks.
  uint16_t *starts = nullptr;  ///< First alternate per mask word.
  uint8_t *boundary = nullptr; ///< Current ring's boundary set.
  int ring = 0;                ///< Current authored ring.
  int count = 0;               ///< Retained alternate count.
  bool table_ready = false;    ///< Candidate table bounds both windows.
#if HS_ENABLE_TEST_HOOKS
  int duplicate_knots = 0; ///< Retained shifts scanned this frame.
  int stack_flushes = 0;   ///< Two-window scans this frame.
#endif

  /** @brief Allocates packed values and per-ring lookup tables.
   * @param values_arena Packed-value storage.
   * @param table_arena Mask, prefix and boundary storage.
   * @param slots Ring count.
   * @param capacity Packed-value capacity; at least KNOTS. */
  void allocate(Arena &values_arena, Arena &table_arena, int slots,
                int capacity) {
    allocate(values_arena, table_arena, slots, capacity, [](Arena &) {});
  }
  /** @brief Allocates lookup tables around caller-owned staging storage.
   * @tparam Storage Caller allocation callback type.
   * @param values_arena Packed-value storage.
   * @param table_arena Mask, prefix, caller and boundary storage.
   * @param slots Ring count.
   * @param capacity Packed-value capacity; at least KNOTS.
   * @param storage Allocates caller storage after prefixes and before boundary bytes. */
  template <typename Storage>
  void allocate(Arena &values_arena, Arena &table_arena, int slots,
                int capacity, Storage &&storage) {
    HS_CHECK(capacity <= UINT16_MAX,
             "knot alternates capacity exceeds prefix range");
    HS_CHECK(capacity >= KNOTS,
             "knot alternates capacity is smaller than one ring");
    slot_count = slots;
    value_capacity = capacity;
    values = values_arena.allocate_n<float>(capacity);
    masks = table_arena.make_n<uint32_t>(slots * WORDS);
    starts = table_arena.make_n<uint16_t>(slots * WORDS);
    storage(table_arena);
    boundary = table_arena.allocate_n<uint8_t>(KNOTS + 1);
    reset();
  }
  /** @brief Adds the split arena allocations to their footprints.
   * @param values Packed-value footprint.
   * @param table Mask, prefix and boundary footprint.
   * @param slots Ring count.
   * @param capacity Packed-value capacity. */
  static constexpr void add_footprint(ArenaFootprint &values,
                                      ArenaFootprint &table, int slots,
                                      int capacity) {
    add_footprint(values, table, slots, capacity, [](ArenaFootprint &) {});
  }
  /** @brief Adds lookup and interleaved caller storage to their footprints.
   * @tparam Storage Caller footprint callback type.
   * @param values Packed-value footprint.
   * @param table Mask, prefix, caller and boundary footprint.
   * @param slots Ring count.
   * @param capacity Packed-value capacity.
   * @param storage Adds caller allocations after prefixes and before boundary bytes. */
  template <typename Storage>
  static constexpr void add_footprint(ArenaFootprint &values,
                                      ArenaFootprint &table, int slots,
                                      int capacity, Storage &&storage) {
    values.add<float>(capacity);
    table.add<uint32_t>(slots * WORDS);
    table.add<uint16_t>(slots * WORDS);
    storage(table);
    table.add<uint8_t>(KNOTS + 1);
  }
  /** @return Packed-value capacity. */
  int capacity() const { return value_capacity; }
  /** @return Allocated ring count. */
  int slots() const { return slot_count; }
  /** @brief Starts a ring's packed range.
   * @param i Authored ring index. */
  HS_HOT_INLINE void begin_ring(int i) {
    ring = i;
    starts[i * WORDS] = static_cast<uint16_t>(count);
  }
  /** @return Whether another row fits.
   * @param lut_n Row knot count. */
  HS_HOT_INLINE bool fits(int lut_n) const {
    return count + lut_n <= capacity();
  }
  /** @brief Retains a differing late-window value and installs the early value.
   * @param k Knot index.
   * @param primary Value in the knot's window.
   * @param alternate Value in the opposite window.
   * @param time Knot sweep fraction.
   * @param row Receives early-window shifts. */
  HS_HOT_INLINE void retain(int k, float primary, float alternate, float time,
                            float *row) {
    row[k] = primary;
    if (primary == alternate)
      return;
    masks[ring * WORDS + k / 32] |= uint32_t{1} << (k % 32);
    row[k] = time < 0.5f ? primary : alternate;
    values[count++] = time < 0.5f ? alternate : primary;
  }
  /** @brief Builds a ring's per-word packed offsets.
   * @param i Authored ring index. */
  HS_HOT_INLINE void finish_ring(int i) {
    int index = starts[i * WORDS];
    for (int word = 0; word < WORDS; ++word) {
      starts[i * WORDS + word] = static_cast<uint16_t>(index);
      index += std::popcount(masks[i * WORDS + word]);
    }
  }
  /** @brief Widens primary bounds by the retained alternate, if present.
   * @param ring Authored ring index.
   * @param knot Knot index.
   * @param lo Primary lower bound, widened in place.
   * @param hi Primary upper bound, widened in place. */
  HS_HOT_INLINE void widen(int ring, int knot, float &lo, float &hi) const {
    const uint32_t *mask = masks + ring * WORDS;
    const uint32_t bit = uint32_t{1} << (knot % 32);
    if (!(mask[knot / 32] & bit))
      return;
    int index = starts[ring * WORDS + knot / 32];
    index += std::popcount(mask[knot / 32] & (bit - 1));
    lo = fminf(lo, values[index]);
    hi = fmaxf(hi, values[index]);
  }
  /** @brief Empties packed values and invalidates the candidate table. */
  void reset() {
    count = 0;
    table_ready = false;
  }

private:
  int slot_count = 0;     ///< Allocated ring count.
  int value_capacity = 0; ///< Packed-value capacity.
};
/**
   * @tparam SlotCount Slot capacity.
   * @tparam RowSize Entries per knot or hue row.
   * @brief The frame's drawn rings: placement-built DistortedRings with their
   * shift and hue rows, shading values and the ring-to-slot map.
   * @details A frame runs begin_frame(), then next() and commit() per drawn
   * ring, then release(). next() changes no state, so a ring abandoned after
   * its bake costs nothing; only commit() constructs a ring and counts it.
   */
template <int SlotCount, int RowSize> class RingStackPool {
  union RingSlot {
    char empty;
    SDF::DistortedRing ring;
    RingSlot() : empty{} {}
    ~RingSlot() {}
  };
  static_assert(sizeof(RingSlot) == sizeof(SDF::DistortedRing));

public:
  static constexpr int SLOTS = SlotCount; ///< Slot capacity.
  static constexpr int ROW = RowSize;     ///< Entries per knot or hue row.
  using Alternates =
      KnotAlternates<ROW - 1>;         ///< Alternate shifts for each ring.
  static constexpr int8_t CULLED = -1; /**< slot_of() for an undrawn ring. */
  static_assert(SLOTS <= INT8_MAX,
                "the ring-to-slot map is int8_t with CULLED = -1; a larger "
                "pool wraps slot indices negative");
  using ShapeStorage = std::array<RingSlot, SLOTS>; /**< Ring storage. */
  /** @brief Persistent bytes init_storage() takes, alignment included. */
  static constexpr size_t BYTES =
      SLOTS * ROW * (sizeof(float) + sizeof(Pixel)) +
      SLOTS * (2 * sizeof(float) + sizeof(int8_t)) + sizeof(ShapeStorage) +
      3 * alignof(float) + alignof(Pixel) + alignof(int8_t) +
      alignof(ShapeStorage);

  /** @brief The rows of the next free slot. */
  struct Pending {
    float *shift_row; /**< ROW centerline shifts; entry lut_n repeats entry
                           0 to close the polyline. */
    Pixel *hue_row;   /**< ROW hue-rotated ring colors. */
  };

  /** @brief Constructed rings indexed by slot. */
  struct ShapeView {
    ShapeStorage *storage;     ///< Pool slots; not owned.
    const RingStackPool *pool; ///< Slot map and alternate storage.
    /** @brief Bounds both sweep-window values of a knot.
       * @param slot Constructed ring slot.
       * @param knot Knot index.
       * @param lo Receives the smaller shift.
       * @param hi Receives the larger shift. */
    void knot_range(int slot, int knot, float &lo, float &hi) const {
      lo = hi = (*storage)[slot].ring.knots[knot];
      if (!pool->cache)
        return;
      pool->cache->widen(pool->ring_by_slot[slot], knot, lo, hi);
    }
    /** @return Whether the frame's shared candidate table is filled. */
    bool table_ready() const { return pool->cache && pool->cache->table_ready; }
    /** @brief Retains the candidate table across sweep windows. */
    void table_built() const {
      if (pool->cache)
        pool->cache->table_ready = true;
    }
    /**
       * @brief Ring constructed in a slot.
       * @param index Slot index in [0, SLOTS).
       * @return The slot's ring.
       */
    SDF::DistortedRing &operator[](size_t index) const {
      return (*storage)[index].ring;
    }
  };

  /** @brief Allocates zeroed hue rows and ring storage.
   * @param arena Persistent storage arena. */
  HS_FLASH_MEMBER void init_storage(Arena &arena) {
    shift_rows = arena.allocate_n<float>(SLOTS * ROW);
    // Zeroed, not raw: a culled azimuth chunk leaves its columns unbaked.
    hue_rows = arena.make_n<Pixel>(SLOTS * ROW);
    alphas = arena.allocate_n<float>(SLOTS);
    columns = arena.allocate_n<float>(SLOTS);
    slot_by_ring = arena.allocate_n<int8_t>(SLOTS);
    storage = arena.make<ShapeStorage>();
  }

  /** @brief Starts an empty frame with every ring culled.
   * @param frame_cache Optional alternate shifts. */
  HS_NOINLINE_NOCLONE void begin_frame(Alternates *frame_cache = nullptr) {
    cache = frame_cache;
    HS_CHECK(!cache || cache->slots() >= SLOTS,
             "ring stack pool alternate slots are too few");
    HS_CHECK(size() == 0, "ring stack pool begins a frame with live rings");
    for (int i = 0; i < SLOTS; ++i)
      slot_by_ring[i] = CULLED;
  }

  /** @return The next free slot's rows; repeated calls return the same rows. */
  HS_HOT_INLINE Pending next() const {
    HS_CHECK(count < SLOTS, "ring stack pool is full");
    return {shift_rows + count * ROW, hue_rows + count * ROW};
  }

  /**
     * @brief Constructs a ring in the slot next() returned and maps @p ring
     * to it; the pool must not be full.
     * @param ring Ring index in the stack.
     * @param frag_alpha Fragment alpha the shader applies to the ring.
     * @param lut_n Bake columns in the slot's rows.
     * @param ring_args DistortedRing constructor arguments.
     */
  template <typename... RingArgs>
  HS_HOT_FLASH_MEMBER void commit(int ring, float frag_alpha, int lut_n,
                                  RingArgs &&...ring_args) {
    const int s = count;
    ::new (static_cast<void *>(&(*storage)[s].ring))
        SDF::DistortedRing(std::forward<RingArgs>(ring_args)...);
    columns[s] = static_cast<float>(lut_n);
    alphas[s] = frag_alpha;
    slot_by_ring[ring] = static_cast<int8_t>(s);
    ring_by_slot[s] = ring;
    if (cache) {
      auto &shape = (*storage)[s].ring;
      float low = 0.0f, high = 0.0f;
      for (int k = 0; k < lut_n; ++k) {
        float lo, hi;
        shapes().knot_range(s, k, lo, hi);
        low = fminf(low, lo);
        high = fmaxf(high, hi);
      }
      shape.set_knot_band(low, high);
    }
    count = s + 1;
  }

  /**
     * @brief Destroys every committed ring and empties the pool.
     * @details The device ScalarFn member is not trivially destructible.
     */
  void release() {
    for (int s = 0; s < count; ++s) {
      using Ring = SDF::DistortedRing;
      (*storage)[s].ring.~Ring();
#if HS_ENABLE_TEST_HOOKS
      ++destroyed;
#endif
    }
    count = 0;
  }

  /** @return Committed rings this frame. */
  int size() const { return count; }
  /** @return Slot of the ring, or CULLED.
   * @param ring Authored ring index. */
  int8_t slot_of(int ring) const { return slot_by_ring[ring]; }
  /** @return The ring-to-slot map, SLOTS entries. */
  const int8_t *slot_map() const { return slot_by_ring; }
  /** @return The committed rings. */
  ShapeView shapes() const { return {storage, this}; }
  /** @return Hue row of the slot.
   * @param s Committed slot. */
  HS_HOT_INLINE const Pixel *hue_row(int s) const { return hue_rows + s * ROW; }
  /** @return Bake column count of the slot.
   * @param s Committed slot. */
  HS_HOT_INLINE float lut_columns(int s) const { return columns[s]; }
  /** @return Fragment alpha of the slot.
   * @param s Committed slot. */
  HS_HOT_INLINE float frag_alpha(int s) const { return alphas[s]; }

  /** @brief Replaces one slot's fragment alpha.
   * @param slot Committed slot.
   * @param alpha Fragment opacity. */
  void set_alpha(int slot, float alpha) { alphas[slot] = alpha; }
  /** @return Writable shift row of a committed slot.
   * @param slot Committed slot. */
  float *shift_row(int slot) { return shift_rows + slot * ROW; }
  /** @return Writable hue row of a committed slot.
   * @param slot Committed slot. */
  Pixel *mutable_hue_row(int slot) { return hue_rows + slot * ROW; }

private:
  Alternates *cache = nullptr;     ///< Current frame's alternate shifts.
  int ring_by_slot[SLOTS]{};       ///< Authored ring index of each slot.
  float *shift_rows = nullptr;     /**< SLOTS shift rows. */
  Pixel *hue_rows = nullptr;       /**< SLOTS hue rows. */
  float *alphas = nullptr;         /**< Per-slot fragment alpha. */
  float *columns = nullptr;        /**< Per-slot bake column count. */
  int8_t *slot_by_ring = nullptr;  /**< Ring index to slot, or CULLED. */
  ShapeStorage *storage = nullptr; /**< Placement storage for the rings. */
  int count = 0;                   /**< Committed rings. */
#if HS_ENABLE_TEST_HOOKS
public:
  int destroyed = 0; /**< Rings release() has destroyed. */
#endif
};
/**
   * @tparam Pool Ring pool type.
 * @brief Shades a stack fragment from its slot's hue row.
   * @details v2 is the stroke coverage the scan applies again on plot, so the
   * ring edge ramps as coverage squared. The stack hands v0 over in [0, 1).
   */
template <typename Pool> struct HueRowShader {
  const Pool *pool; ///< The frame's drawn rings.

  /**
     * @brief Shades one fragment.
     * @param s Slot of the ring lighting the fragment.
     * @param f Fragment to shade; v0 is the ring azimuth in [0, 1).
     */
  HS_HOT_INLINE void operator()(int s, const math::Vector &,
                                Fragment &f) const {
    const Pixel *hue = pool->hue_row(s);
    float x = f.v0 * pool->lut_columns(s);
    int j = static_cast<int>(x);
    f.color = Color4(
        hue[j].lerp16(hue[j + 1], frac_to_q16(math::quintic_kernel(x - j))),
        pool->frag_alpha(s) * f.v2);
  }
};
/** @brief Column-time sampling and two-window ring-stack scans. */
struct SweptRingStack {
  static constexpr float QUARTER = 0.25f; ///< Quarter-frame split.
  static constexpr int WINDOWS = 2;       ///< Half-turn sampling windows.
  /** @return Whether a sample needs its opposite endpoint.
   * @param t Frame fraction. */
  HS_HOT_INLINE static bool retains_endpoint(float t) {
    return t < QUARTER || t >= 1.0f - QUARTER;
  }
  /** @return Opposite endpoint for a sample.
   * @param t Frame fraction. */
  HS_HOT_INLINE static float other_endpoint(float t) {
    return t < 0.5f ? 1.0f : 0.0f;
  }
  /**
   * @brief Frame fraction at which a swept draw evaluates a ring point.
   * @param p The point where the stack holds it at the current endpoint.
   * @param middle Middle sweep fraction of the drawn window.
   * @return The sweep fraction of the column showing @p p, taken within half
   *         a frame of @p middle and clamped to the frame.
   * @details The window's columns read only points near them, so points
   * shown across an arm boundary from the window take the boundary's time.
   */
  HS_HOT_INLINE static float sweep_fraction(const math::Vector &p,
                                            float middle) {
    float u = math::fast_atan2(p.z, p.x) * (1.0f / math::PI_F);
    u -= floorf(u);
    if (u - middle >= 0.5f)
      u -= 1.0f;
    else if (middle - u > 0.5f)
      u += 1.0f;
    return hs::clamp(u, 0.0f, 1.0f);
  }

  /** @brief Marks knots whose supported geometry can cross an arm boundary.
   * @param basis Ring frame.
   * @param step Ring and knot-step trigonometry.
   * @param theta Ring colatitude.
   * @param n Knot count.
   * @param displacement Bound on polar displacement.
   * @param motion Frame rotation bound in radians.
   * @param thickness Stroke half-width.
   * @param support_cells Spline stencil reach in cells.
   * @param out Receives the boundary mask.
   * @details The supported azimuth interval covers the spline stencil and
   * bounded segment search. Its full displaced polar band is tested against
   * the arm plane widened by frame rotation and one column. */
  template <int W>
  HS_HOT_FLASH_MEMBER static void
  mark_arm_boundary(const math::Basis &basis, SDF::KnotStep step, float theta,
                    int n, float displacement, float thickness, float motion,
                    int support_cells, uint8_t *out) {
    const float low = fmaxf(0.0f, theta - displacement - thickness);
    const float high = fminf(math::PI_F, theta + displacement + thickness);
    const float min_sine = fmaxf(sqrtf(SDF::DistortedRing::POLE_SIN2_FLOOR),
                                 fminf(sinf(low), sinf(high)));
    const float cell = math::TWO_PI_F / n;
    const float search =
        SDF::DistortedRing::search_reach_cells(thickness, min_sine, n);
    const float delta = fminf(math::PI_F, (search + support_cells) * cell);
    const float epsilon =
        sinf(fminf(math::PI_F / 2, motion + math::TWO_PI_F / W)) + 1e-4f;
    const float axis = basis.v.z;
    const bool pole =
        (low == 0.0f || high == math::PI_F) && std::fabs(axis) <= epsilon;
    auto threshold = [&](float edge, bool minimum) HS_HOT_INLINE {
      auto at = [&](float angle) HS_HOT_INLINE {
        return (edge - axis * cosf(angle)) / fmaxf(1e-6f, sinf(angle));
      };
      float value =
          minimum ? fminf(at(low), at(high)) : fmaxf(at(low), at(high));
      const float cosine = axis / edge;
      if (std::fabs(cosine) < 1.0f) {
        const float angle = acosf(cosine);
        if (angle > low && angle < high)
          value = minimum ? fminf(value, at(angle)) : fmaxf(value, at(angle));
      }
      return value;
    };
    const float target_low = threshold(-epsilon, true);
    const float target_high = threshold(epsilon, false);
    const float radius = sqrtf(basis.u.z * basis.u.z + basis.w.z * basis.w.z);
    const float cosine = cosf(delta), sine = sinf(delta);
    SDF::KnotWalk walk{basis, step.cos_t, step.sin_t, step.cos_d, step.sin_d};
    for (int k = 0; k < n; ++k, walk.advance()) {
      const float center = basis.u.z * walk.cos_a + basis.w.z * walk.sin_a;
      const float tangent = basis.w.z * walk.cos_a - basis.u.z * walk.sin_a;
      float lo = center * cosine - std::fabs(tangent * sine);
      float hi = center * cosine + std::fabs(tangent * sine);
      if (center >= radius * cosine)
        hi = radius;
      if (center <= -radius * cosine)
        lo = -radius;
      out[k] = pole || (lo <= target_high + 1e-4f && hi >= target_low - 1e-4f);
    }
  }

  /** @brief Scans early and late knot windows, then releases frame storage.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @tparam PipelineT Scan pipeline.
   * @tparam Pool Ring pool type.
   * @tparam Alternates Packed alternate type.
   * @tparam Table Candidate table type.
   * @tparam Shader Fragment shader type.
   * @tparam Reshade Value-owned late-window shading visitor.
   * @param pipeline Scan pipeline.
   * @param canvas Destination canvas.
   * @param n_rings Authored ring count.
   * @param pool Committed rings.
   * @param alternates Packed late-window shifts.
   * @param table Candidate table.
   * @param shader Fragment shader.
   * @param reshade Calls begin(ring, slot), operator()(knot, shift) for each
   * installed late knot, then finish() after closing the shift row. */
  template <int W, int H, typename PipelineT, typename Pool,
            typename Alternates, typename Table, typename Shader,
            typename Reshade>
  HS_HOT_FLASH_MEMBER static void
  scan_two_windows(PipelineT &pipeline, Canvas &canvas, int n_rings, Pool &pool,
                   Alternates &alternates, Table &table, Shader &shader,
                   Reshade reshade) {
    static_assert(W % (2 * WINDOWS) == 0,
                  "sweep windows must tile the half-turn");
    if (pool.size() > 0) {
      for (int window_index = 0; window_index < WINDOWS; ++window_index) {
        if (window_index) {
          for (int ring = 0; ring < n_rings; ++ring) {
            const int slot = pool.slot_of(ring);
            if (slot < 0)
              continue;
            const int n = static_cast<int>(pool.lut_columns(slot));
            float *row = pool.shift_row(slot);
            const uint32_t *mask = alternates.masks + ring * Alternates::WORDS;
            int index = alternates.starts[ring * Alternates::WORDS];
            reshade.begin(ring, slot);
            for (int k = 0; k < n; ++k) {
              if (mask[k / 32] & (uint32_t{1} << (k % 32)))
                row[k] = alternates.values[index++];
              reshade(k, row[k]);
            }
            row[n] = row[0];
            reshade.finish();
          }
        }
        TimeWindowScope window(canvas, window_index * W / (2 * WINDOWS),
                               (window_index + 1) * W / (2 * WINDOWS));
        DistortedRingStack::draw<W, H>(pipeline, canvas, n_rings, pool.shapes(),
                                       pool.slot_map(), pool.size(), table,
                                       shader);
      }
    }
    pool.release();
#if HS_ENABLE_TEST_HOOKS
    alternates.duplicate_knots += alternates.count;
    ++alternates.stack_flushes;
#endif
    alternates.reset();
  }
};

} // namespace Scan
