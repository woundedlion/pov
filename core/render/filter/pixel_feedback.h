/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include "core/math/spherical_field.h"
#include "core/memory.h"
#include "core/render/filter/feedback_cap_plane.h"
#include "core/render/filter/feedback_column_time.h"
#include "core/render/filter/feedback_style.h"
#include "core/render/filter/feedback_warp_cache.h"
#include "core/render/filter/pipeline.h"
#include "core/render/spacetime.h"

/**
 * @file pixel_feedback.h
 * @brief Filter::Pixel::Feedback: the terminal filter that warps and
 * composites the previous frame.
 */

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
struct FeedbackWhiteBox;
}
#endif

namespace Filter {

namespace Pixel {

/**
 * @brief Style-aware terminal feedback filter that warps the previous frame.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The Style's spatial warp is computed on a spherical latitude-ring
 * field, then interpolated within and between rings. Terminal: it must be the
 * last Pipeline stage, and Pipeline::begin_frame() must run before the frame's
 * plot() calls, or Pipeline::begin_split_frame() before them and
 * Pipeline::finish_frame() after; flushing last blanks the frame at
 * alpha >= 1.
 * Rows under POLAR_TARGET_SINE interpolate offsets in their pole's cap plane,
 * where equirect longitude offsets blow up as 1/sin(phi). A pole row is one
 * point, composited once from its source a frame earlier.
 * With temporal sampling on, each pixel reads its source at the time elapsed
 * between their lit columns (FeedbackColumnTime). A source behind the target
 * arm's first column reads the stale frame the drawing buffer still holds;
 * one past the arm's last column reads the frame being drawn after its plots,
 * so the pixel is held: it keeps the part of its composite the previous
 * frames supply, the plots blend over it while its coverage tracks their
 * transmittance, and finish() adds the rest beneath that coverage. A frame
 * without the stale frame, or an unsplit one, reads such a source from the
 * previous frame at its age there.
 */
template <int W, int H> class Feedback : public Is2DWithHistory {
  using SphereField = hs::SphericalFieldLayout<W, H>;
  using CapPlane = FeedbackCapPlane<W, H>;
  using CapOffset = typename CapPlane::CapOffset;
  using CapPoint = typename CapPlane::CapPoint;
  using CapCell = typename CapPlane::CapCell;
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::FeedbackWhiteBox;
#endif
  using WarpCache = FeedbackWarpCache<W, H>;
  using WarpField = typename WarpCache::Buffers;
  using ColumnTime = FeedbackColumnTime<W>;
  using FadeTable = FeedbackFadeTable;
  static constexpr int HALF = ColumnTime::HALF;
  /** @brief A colour transform: one sample and its fade-table level to a
   *  pixel. */
  using PixelTransform = FunctionRef<::Pixel(float, float, float, int)>;

  /** @brief Constructs the shared spherical field layout.
   * @param ds Field downsample factor.
   * @return Matching north/south infill and rectangular column count. */
  static HS_HOT_INLINE constexpr SphereField field_for(int ds) {
    return SphereField(ds, ds, ds, W / ds);
  }

  /** @brief Coarse grid downsample the warp cache is sized for (the default
   *  Style's). Other values render uncached. */
  static constexpr int CACHE_DOWNSAMPLE = ::Feedback::Style{}.downsample;
  static constexpr int CACHE_COLUMNS = W / CACHE_DOWNSAMPLE;
  static constexpr SphereField CACHE_FIELD = field_for(CACHE_DOWNSAMPLE);
  /** @brief Cell count of the cached spherical warp field. */
  static constexpr int CACHE_CELLS = CACHE_COLUMNS * CACHE_FIELD.ring_count();

public:
  /// Pixel domain: runs after every screen stage.
  static constexpr int domain_rank = IsPixel::domain_rank;
  /** @brief Marks this as terminal: flush() writes the Canvas directly. */
  static constexpr bool is_terminal = true;
  /** @brief Latitude sine under which lattice rings carry cap-plane offsets
   *  and rows between two such rings convert each pixel's own target. */
  static constexpr float POLAR_TARGET_SINE = 0.5f;
  /** @brief Opaque store owns the frame: no history stage may precede it. */
  static constexpr bool terminal_replaces = true;

  // Covers only the default-constructed Style; a runtime-swapped style's
  // downsample is checked per flush.
  static_assert(
      ::Feedback::Style{}.downsample > 0 &&
          W % ::Feedback::Style{}.downsample == 0,
      "Feedback<W,H>: default style downsample must be > 0 and divide "
      "W");

  /**
   * @brief Binds the filter to a live feedback Style.
   * @param style Style supplying the spatial warp and color transforms.
   */
  explicit Feedback(::Feedback::Style &style) : feedback_style(&style) {}

  /**
   * @brief Passes a current-frame tap to the next stage, or blends it into a
   *        held pixel of a split frame with its coverage.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  HS_HOT_INLINE void plot(float x, float y, const ::Pixel &color, float age,
                          float alpha, PassFnT &&pass) {
    const int column = Filter::round_wrap_column<W>(x);
    const int past = column - ColumnTime::arm_start(column);
    if (past >= held_floor) [[unlikely]] {
      const int row = Filter::round_row<H>(y);
      if (row >= 0 && row < H && past >= row_floor[row]) {
        plot_near_arm_end(x, column, y, color, age, alpha, pass);
        return;
      }
    }
    pass(x, y, color, age, alpha);
  }

  /**
   * @brief Enables or disables feedback.
   * @param value When false, flush() only clears a stale drawing buffer.
   */
  void set_enabled(bool value) { enabled = value; }

  /**
   * @brief Scratch bytes for a full-width uncached flush at downsample ds.
   * @param ds Downsample factor.
   * @return Scratch byte count.
   */
  static size_t uncached_scratch_bytes(int ds) {
    HS_CHECK(ds > 0 && W % ds == 0,
             "feedback downsample %d must be > 0 and divide width %d", ds, W);
    const int columns = W / ds;
    const SphereField field = field_for(ds);
    const int rings = field.ring_count();
    const PolarRings polar = polar_rings(field, rings);
    return (2 * rings * columns + 2 * field.sample_count()) * sizeof(int16_t) +
           (polar.rows() * columns + polar.samples) * sizeof(CapOffset) +
           ((ds > 1 || !SphereField::HAS_NORTH_POLE ||
             !SphereField::HAS_SOUTH_POLE)
                ? W * sizeof(::Pixel)
                : 0);
  }

  /**
   * @brief Allocates the warp-field cache and fade table from the persistent
   * arena.
   * @param arena Persistent arena supplying STORAGE_BYTES bytes.
   * @details Call from effect init(), not the constructor, and again after
   * any arena reset. Without storage every flush needs
   * uncached_scratch_bytes(downsample) scratch bytes, and with temporal
   * sampling on FeedbackFadeTable::BYTES more.
   */
  HS_FLASH_MEMBER void init_storage(Arena &arena) {
    // Projected from the same incremental-rotation positions populate() hands
    // out, so the per-frame offset subtracts them exactly.
    const int cap_cells =
        polar_rings(CACHE_FIELD, CACHE_FIELD.ring_count()).rows() *
        CACHE_COLUMNS;
    // no death case: the reserve covers every display geometry the cap rows admit
    HS_CHECK(cap_cells <= CACHE_CAP_CELLS,
             "feedback: %d polar cap cells exceed the %d reserved", cap_cells,
             CACHE_CAP_CELLS);
    warp_cache.template init_storage<CACHE_CELLS>(
        arena, CACHE_FIELD, cap_cells,
        [](const math::Vector &position,
           const typename SphereField::Coordinates &point) {
          return lattice_origin(CACHE_FIELD, position, point);
        });
    fade_table.init_storage(arena);
    spans.reach = arena.allocate_n<int16_t>(H);
    spans.valid = false;
  }

  /**
   * @brief Accesses the bound Style.
   * @return Mutable reference to the bound feedback Style.
   */
  ::Feedback::Style &style() { return *feedback_style; }
  /**
   * @brief Accesses the bound Style.
   * @return Const reference to the bound feedback Style.
   */
  const ::Feedback::Style &style() const { return *feedback_style; }

private:
  template <int, int, typename...> friend struct ::Pipeline;

  /**
   * @brief Blends the distorted previous frame into the current frame.
   * @param cv Target canvas (reads cv.prev, writes the back (current-draw) buffer).
   * @param alpha Global blend alpha in [0, 1].
   * @details Computes a coarse warp field via the Style's space_fn, bilinearly
   * upsamples it, then composites the warped previous frame, honoring the
   * segment's cylindrical clip. Clears a stale drawing buffer unless column time reads it as the frame two back; a disabled filter only clears it.
   */
  HS_O3_BEGIN
  void flush(Canvas &cv, float alpha) {
    if (disabled(cv))
      return;
    ScratchScope scope(scratch_arena_a);
    FlushContext ctx = prepare_flush(cv, alpha, scope.get_arena());
    composite_frame(cv, ctx);
  }
  HS_O3_END

  /**
   * @brief Clears the frame a disabled filter leaves to the canvas.
   * @param cv Target canvas.
   * @return Whether the filter is disabled.
   */
  bool disabled(Canvas &cv) {
    require_unheld();
    if (enabled)
      return false;
    if (cv.keeps_stale_pixels())
      cv.clear_buffer();
    return true;
  }

  /**
   * @brief Composites a split frame short of the pixels whose sources lie
   *        past their arm's last column, and holds the frame for finish().
   * @param cv Target canvas.
   * @param alpha Global blend alpha in [0, 1].
   * @details The held frame keeps its scratch_arena_a allocations until
   * finish(); scopes opened meanwhile must close first. A frame no pixel of
   * which reads past a cut is released at once.
   */
  HS_FLASH_MEMBER void flush_split(Canvas &cv, float alpha) {
    if (disabled(cv))
      return;
    const size_t mark = scratch_arena_a.get_offset();
    FlushContext &ctx = *scratch_arena_a.make<FlushContext>(
        prepare_flush(cv, alpha, scratch_arena_a, true));
    held_floor = HALF;
    std::memset(row_floor, HALF, sizeof(row_floor));
    composite_frame(cv, ctx);
    if (held_floor == HALF) {
      scratch_arena_a.set_offset(mark);
      return;
    }
    held = &ctx;
    held_mark = mark;
  }

  /**
   * @brief Adds the part of each held pixel's composite that reads the frame
   *        as plotted, beneath the coverage its taps left there.
   * @details No-op without a held frame.
   */
  HS_FLASH_MEMBER void finish(Canvas &) {
    if (!held)
      return;
    FlushContext &ctx = *held;
    select_poles(ctx.current, ctx.current_poles);
    {
      HS_PROFILE(feedback_forward);
      composite_pass(ctx, true);
    }
    end_split();
  }
  HS_O3_BEGIN

private:
  /** @brief The lattice rings carrying 3D targets: rings [0, north_rings) and
   *  [south_ring, field_rows), with the sample index where each run ends or
   *  starts. */
  struct PolarRings {
    int north_rings;
    int north_samples;
    int south_ring;
    int south_sample;
    int samples;
    int ring_count;

    constexpr int rows() const { return north_rings + ring_count - south_ring; }
    /** @brief Compacted row of a polar ring. */
    int row(int field_y) const {
      return field_y < north_rings ? field_y
                                   : north_rings + field_y - south_ring;
    }
    /** @brief Compacted slot of a polar ring's lattice sample. */
    int slot(int sample) const {
      return sample < north_samples ? sample
                                    : north_samples + sample - south_sample;
    }
  };

  struct CoarseGrid {
    int downsample;
    SphereField field;
    int columns;
    int field_rows;
    PolarRings polar;
  };

  struct RenderBand {
    int y_begin;
    int y_end;
    int field_y_begin;
    int field_y_end;
    ClipRegion::XClip x_clip;
    std::bitset<W> coarse_columns_used;
  };

  struct WarpControl {
    int16_t x;
    int16_t y;
  };

  struct PixelAccumulator {
    uint32_t r = 0;
    uint32_t g = 0;
    uint32_t b = 0;

    void add(const ::Pixel &pixel) {
      r += pixel.r;
      g += pixel.g;
      b += pixel.b;
    }

    void remove(const ::Pixel &pixel) {
      r -= pixel.r;
      g -= pixel.g;
      b -= pixel.b;
    }

    ::Pixel average(int width) const {
      const uint32_t round = static_cast<uint32_t>(width / 2);
      return ::Pixel(static_cast<uint16_t>((r + round) / width),
                     static_cast<uint16_t>((g + round) / width),
                     static_cast<uint16_t>((b + round) / width));
    }
  };

  struct ColumnRun {
    int begin;
    int end;
  };

  struct ColumnRuns {
    ColumnRun items[6];
    int count;
  };

  enum class FlushMode : uint8_t { SKIP, HUE, PLAIN, GENERAL };

  /** @brief Where a lane's source lies relative to its target's arm. */
  enum class Side : uint8_t {
    WITHIN, ///< Within the arm: the previous frame.
    BEHIND, ///< Behind the arm's first column: the stale frame before it.
    AHEAD,  ///< From the arm's last column on: the frame being drawn.
  };

  /** @brief One bilinear sample on the [0, 65535] scale, unquantized. */
  struct Sample {
    float r;
    float g;
    float b;
  };

  /** @brief A lane's resolved read. */
  struct Lane {
    typename ColumnTime::Read read;
    Side side;
  };

  /** @brief A row's interpolation and reconstruction, resolved per frame. */
  struct RowGeometry {
    int16_t field_y0; ///< Ring at or above the row.
    int16_t field_y1; ///< Ring below it, or the same ring on the last.
    float wy1;        ///< Weight of field_y1.
    bool polar : 1;   ///< Both rings carry cap-plane offsets.
    bool north : 1;   ///< A polar row's pole is the north pole.
    bool pairs : 1;   ///< Composites column pairs.
    bool infill : 1;  ///< Reconstructs at its longitude footprint.
    bool pole : 1;    ///< Every column is the pole.
    bool wide : 1;    ///< Covers every column of each arm (cover_slot()).
    /// Lanes past a column that its reconstruction reads: the next pair, or
    /// the infill filter's radius.
    uint8_t reach;
    /// First of the row's finish samples (finish_window()).
    int16_t window_base;
    /// First of the row's coverage slots.
    uint16_t cover_base;
  };

  /** @brief Lane columns of a row whose sources read the frame as plotted. */
  struct AheadSpan {
    int16_t first; ///< First such column past each arm's first; HALF if none.
    int16_t last;  ///< Last such column, -1 if none; finish() extends it to
                   ///< the last lane whose step completes a held column.
  };

  /** @brief Per-frame composite inputs prepare_flush() resolves. */
  struct FlushContext {
    CoarseGrid grid;
    RenderBand band{};
    /// The warp field, each node at its column's lit time.
    WarpField warp{};
    /// The field at each arm's end, which the right nodes of the arm's last
    /// cells read; @c warp itself when the field holds still.
    WarpField end{};
    ::Pixel *filtered_row = nullptr;
    const ::Pixel *previous = nullptr;
    ::Pixel *current = nullptr;
    ::Pixel poles[SphereField::POLE_STORAGE_COUNT]{};
    /// The drawing buffer's stale frame, two frames back, which holds the
    /// content behind each arm's backward cut; null when it was cleared.
    const ::Pixel *back = nullptr;
    /// Shared pole values of the stale frame.
    ::Pixel back_poles[SphereField::POLE_STORAGE_COUNT]{};
    /// Per row, the even column count from each arm's first whose lanes may
    /// read behind its backward cut. Null without @c back.
    const int16_t *reach = nullptr;
    /// Per row, its interpolation and reconstruction.
    const RowGeometry *rows = nullptr;
    /// Per coverage slot, the transmittance of the taps plotted over the
    /// pixel in its upper 15 bits, zero where the pixel is final, and in its
    /// low bit whether the lane's source reads the frame as plotted. Null
    /// when unsplit.
    uint16_t *cover = nullptr;
    /// Per row, the lanes whose sources read the frame as plotted.
    AheadSpan *ahead_span = nullptr;
    /// Shared pole values of the frame being drawn, set by finish().
    ::Pixel current_poles[SphereField::POLE_STORAGE_COUNT]{};
    ColumnRuns runs{};
    float alpha = 0.0f;
    float fade = 0.0f;
    float pole_half_res = 0.0f;
    bool black_skips_color = false;
    /// Whether pixels read their sources at their elapsed column time.
    bool column_time = false;
    FlushMode mode = FlushMode::SKIP;
    float hue_k[9] = {};
    /// Per-elapsed-time transforms; read only with column_time.
    FadeTable::View fades{};
  };

  __attribute__((always_inline)) CoarseGrid
  make_coarse_grid(const Canvas &cv) const {
    const int downsample = feedback_style->downsample;
    HS_CHECK(downsample > 0 && W % downsample == 0,
             "feedback downsample %d must be > 0 and divide width %d",
             downsample, W);
    HS_CHECK(cv.width() == W,
             "feedback canvas width %d must equal template W %d", cv.width(),
             W);
    HS_CHECK(cv.height() == H,
             "feedback canvas height %d must equal template H %d", cv.height(),
             H);
    const int columns = W / downsample;
    const SphereField field = field_for(downsample);
    const int rings = field.ring_count();
    return {downsample, field, columns, rings, polar_rings(field, rings)};
  }

  static __attribute__((always_inline)) constexpr PolarRings
  polar_rings(const SphereField &field, int rings) {
    PolarRings polar{0, 0, rings, 0, 0, rings};
    bool leading = true;
    int total = 0;
    auto ring = field.ring(0);
    for (int i = 0; i < rings; ++i, ring = field.next_ring(ring)) {
      const bool is_polar =
          SphereField::latitude_sine(ring.y) < POLAR_TARGET_SINE;
      if (leading && is_polar) {
        polar.north_rings = i + 1;
        polar.north_samples = ring.offset + ring.samples;
      } else {
        leading = false;
        if (!is_polar) {
          polar.south_ring = rings;
        } else if (polar.south_ring == rings) {
          polar.south_ring = i;
          polar.south_sample = ring.offset;
        }
      }
      total = ring.offset + ring.samples;
    }
    if (polar.south_ring == rings)
      polar.south_sample = total;
    polar.samples = polar.north_samples + total - polar.south_sample;
    return polar;
  }

  static __attribute__((always_inline)) RenderBand
  make_render_band(const ClipRegion &clip, const CoarseGrid &grid) {
    RenderBand band{};
    band.y_begin = clip.render_y_start();
    band.y_end = clip.render_y_end();
    band.x_clip = clip.x_clip();

    band.field_y_begin = grid.field.ring_index_at_or_before(band.y_begin);
    band.field_y_end = grid.field.ring_index_at_or_after(band.y_end - 1);
    HS_CHECK(band.field_y_end >= band.field_y_begin,
             "feedback field band inverted: [%d,%d]", band.field_y_begin,
             band.field_y_end);

    if (band.x_clip.active) {
      for (int x = 0; x < W; ++x) {
        if (band.x_clip.clipped(x))
          continue;
        const int left = x / grid.downsample;
        const int right = (left + 1 < grid.columns) ? left + 1 : 0;
        band.coarse_columns_used[left] = true;
        band.coarse_columns_used[right] = true;
      }
    }
    return band;
  }

  static __attribute__((always_inline)) ColumnRuns
  make_column_runs(const ClipRegion::XClip &clip, bool split_arms = false) {
    ColumnRuns runs{};
    const auto append = [&](int begin, int end) {
      if (split_arms && begin < HALF && end > HALF) {
        runs.items[runs.count++] = {begin, HALF};
        begin = HALF;
      }
      runs.items[runs.count++] = {begin, end};
    };
    if (!clip.active) {
      append(0, W);
    } else if (clip.wrap) {
      if (clip.re > 0)
        append(0, clip.re);
      if (clip.rs < W)
        append(clip.rs, W);
    } else {
      append(clip.rs, clip.re);
    }
    return runs;
  }

  /** @brief Builds the cache identity for the current warp and render band.
   * @param band Rows and columns the frame renders.
   * @return Warp cache identity. */
  typename WarpCache::Key warp_cache_key(const RenderBand &band) const {
    const Animation::NoiseParams *noise = feedback_style->noise;
    return {feedback_style->space_fn,  noise,
            noise_config_key(noise),   feedback_style->amplitude,
            feedback_style->frequency, feedback_style->speed,
            feedback_style->scale,     noise ? noise->time : 0.0f,
            band.field_y_begin,        band.field_y_end,
            !Spacetime::flat_surface()};
  }

  /**
   * @brief The frame's populated warp field: the cache's when the Style and
   * band allow it, otherwise a scratch field.
   * @param reused Set when the frame reuses the cached field unchanged.
   * @param prior Receives the field one noise step before the returned one.
   */
  __attribute__((always_inline)) WarpField
  acquire_warp_field(Arena &scratch, const CoarseGrid &grid,
                     const RenderBand &band, bool &reused, WarpField &prior) {
    const bool stock_transform =
        feedback_style->space_fn == &::Feedback::noise_warp ||
        feedback_style->space_fn == &::Feedback::melt_warp;
    const bool cacheable = warp_cache.ready() && !band.x_clip.active &&
                           grid.downsample == CACHE_DOWNSAMPLE &&
                           stock_transform;
    if (cacheable)
      warp_cache.check_storage_alive();

    WarpField uncached{};
    if (!cacheable) {
      HS_CHECK(uncached_scratch_bytes(grid.downsample) <=
                   scratch.get_capacity() - scratch.get_offset(),
               "uncached feedback needs more scratch: missing cache, custom "
               "SpaceFn, nondefault downsample, or x clip");
      const int cells = grid.field_rows * grid.columns;
      const int polar_cells = grid.polar.rows() * grid.columns;
      uncached.x_offsets = scratch.allocate_n<int16_t>(cells);
      uncached.y_offsets = scratch.allocate_n<int16_t>(cells);
      uncached.cell_caps = polar_cells > 0
                               ? scratch.allocate_n<CapOffset>(polar_cells)
                               : nullptr;
    }

    const auto key = warp_cache_key(band);

    HS_PROFILE(feedback_populate);
    reused = cacheable;
    return warp_cache.acquire(
        cacheable ? &key : nullptr, uncached,
        [&](const WarpField &warp) __attribute__((always_inline)) {
          reused = false;
          ScratchScope populate_scope(scratch);
          WarpControl *controls =
              scratch.allocate_n<WarpControl>(grid.field.sample_count());
          CapOffset *sample_caps =
              grid.polar.samples > 0
                  ? scratch.allocate_n<CapOffset>(grid.polar.samples)
                  : nullptr;
          populate_warp_field(grid, band, warp, controls, sample_caps);
        },
        &prior);
  }

  /**
   * @brief Fills @p warp's cells for the band from the Style's warp.
   * @param grid Coarse layout of the frame.
   * @param band Rows and columns the frame renders.
   * @param warp Field to fill.
   * @param controls Scratch for every lattice sample's offset.
   * @param sample_caps Scratch for every polar lattice sample's cap offset.
   * @details Under column time a cell's offsets carry their elapsed scale,
   * which caps a forward column offset at a quarter turn; a pole ring's cap
   * offsets stay unscaled.
   */
  __attribute__((always_inline)) void
  populate_warp_field(const CoarseGrid &grid, const RenderBand &band,
                      const WarpField &warp, WarpControl *controls,
                      CapOffset *sample_caps) {
    hs::SphericalField<WarpControl, W, H> compact(controls, grid.field);
    // The cached origins belong to the cached layout.
    const typename SphereField::Coordinates *origins =
        grid.downsample == CACHE_DOWNSAMPLE ? warp_cache.origins() : nullptr;
    if (origins)
      warp_cache.check_storage_alive();
    // noinline keeps the per-sample warp out of flash-resident prepare_flush.
    compact.populate(
        band.field_y_begin, band.field_y_end,
        [&](const math::Vector &position,
            const typename SphereField::Coordinates &point,
            int index) __attribute__((noinline)) {
          math::Vector distorted;
          {
            HS_PROFILE_DEEP(fb_pop_warp);
            distorted = feedback_style->space_fn(position, *feedback_style);
          }
          HS_PROFILE_DEEP(fb_pop_project);
          const PolarRings &polar = grid.polar;
          if (index < polar.north_samples || index >= polar.south_sample) {
            const bool south = index >= polar.south_sample;
            const CapPoint from = CapPlane::cap_point(position, south);
            const CapPoint to = CapPlane::cap_point(distorted, south);
            sample_caps[polar.slot(index)] =
                CapPlane::encode_cap(to.u - from.u, to.v - from.v);
          }
          const auto projected = grid.field.project(distorted);
          const auto origin = origins
                                  ? origins[index]
                                  : lattice_origin(grid.field, position, point);
          float x_offset = projected.x - origin.x;
          const float y_offset = projected.y - origin.y;
          x_offset = unwrap_near(x_offset, 0.0f, static_cast<float>(W));
          return WarpControl{static_cast<int16_t>(hs::clamp(
                                 x_offset * WARP_SCALE, -32767.0f, 32767.0f)),
                             static_cast<int16_t>(hs::clamp(
                                 y_offset * WARP_SCALE, -32767.0f, 32767.0f))};
        });

    {
      HS_PROFILE_DEEP(fb_pop_expand);
      constexpr float INVERSE_WARP_SCALE = 1.0f / WARP_SCALE;
      const bool column_time = !Spacetime::flat_surface();
      const float forward_limit =
          column_time ? WRAP_PERIOD * 0.25f : WRAP_PERIOD * 0.5f;
      auto ring = grid.field.ring(band.field_y_begin);
      for (int field_y = band.field_y_begin; field_y <= band.field_y_end;
           ++field_y, ring = grid.field.next_ring(ring)) {
        const bool pole_ring = is_pole_row(ring.y);
        for (int coarse_x = 0; coarse_x < grid.columns; ++coarse_x) {
          if (band.x_clip.active && !band.coarse_columns_used[coarse_x])
            continue;
          const int x = coarse_x * grid.downsample;
          const auto longitude = grid.field.longitude_bounded(ring, x);
          const WarpControl a = controls[longitude.left];
          const WarpControl b = controls[longitude.right];
          const float bx = unwrap_near(b.x, a.x, WRAP_PERIOD);
          const int index = field_y * grid.columns + coarse_x;
          // Keep the stored offset canonical: the seam correction can lift the
          // interpolant a full period out, past both the int16_t range and the
          // single-step wrap sample_bilinear contracts for.
          const float offset_x =
              unwrap_near(hs::lerp(static_cast<float>(a.x), bx, longitude.mix),
                          0.0f, WRAP_PERIOD);
          const float offset_y = hs::lerp(
              static_cast<float>(a.y), static_cast<float>(b.y), longitude.mix);
          const float scale =
              column_time
                  ? ColumnTime::elapsed_scale(offset_x * INVERSE_WARP_SCALE)
                  : 1.0f;
          warp.x_offsets[index] = static_cast<int16_t>(
              hs::clamp(offset_x * scale, -WRAP_PERIOD * 0.5f, forward_limit));
          warp.y_offsets[index] = static_cast<int16_t>(
              hs::clamp(offset_y * scale, -32767.0f, 32767.0f));
          if (field_y < grid.polar.north_rings ||
              field_y >= grid.polar.south_ring) {
            const CapOffset da = sample_caps[grid.polar.slot(longitude.left)];
            const CapOffset db = sample_caps[grid.polar.slot(longitude.right)];
            CapOffset cap = CapPlane::lerp_offset(da, db, longitude.mix);
            if (column_time && !pole_ring)
              cap = CapPlane::scale_offset(cap, scale);
            warp.cell_caps[grid.polar.row(field_y) * grid.columns + coarse_x] =
                cap;
          }
        }
      }
    }
  }

  /**
   * @brief Resolves the frame's warp field and composite inputs.
   * @param cv Target canvas.
   * @param alpha Global blend alpha in [0, 1].
   * @param scratch Arena for the frame's uncached warp field and filter row.
   * @param split Whether the frame's held pixels composite after its plots.
   * @return The composite context; mode SKIP when the previous frame is black.
   */
  HS_HOT_FLASH_MEMBER FlushContext prepare_flush(Canvas &cv, float alpha,
                                                 Arena &scratch,
                                                 bool split = false) {
    const Animation::NoiseParams *bound = feedback_style->noise;
    HS_CHECK(!bound || (bound->amplitude == feedback_style->amplitude &&
                        bound->frequency == feedback_style->frequency &&
                        bound->speed == feedback_style->speed &&
                        bound->scale == feedback_style->scale),
             "feedback style scalars never reached the bound NoiseParams; call "
             "Style::sync_noise() after changing them");

    FlushContext ctx{make_coarse_grid(cv)};
    {
      HS_PROFILE(feedback_litscan);
      if (!any_pixel_lit(cv))
        return ctx;
    }

    const CoarseGrid &grid = ctx.grid;
    ctx.band = make_render_band(cv.clip(), grid);
    const RenderBand &band = ctx.band;
    bool reused = false;
    WarpField prior{};
    ctx.end = acquire_warp_field(scratch, grid, band, reused, prior);
    if (!reused)
      spans.valid = false;
    ctx.warp = ctx.end;
    ctx.filtered_row = (!band.x_clip.active && (split || grid.downsample > 1 ||
                                                !SphereField::HAS_NORTH_POLE ||
                                                !SphereField::HAS_SOUTH_POLE))
                           ? scratch.allocate_n<::Pixel>(W)
                           : nullptr;

    const float fade = feedback_style->fade;
    // Guard the float-to-int conversion independently of the check condition.
    HS_CHECK(
        fade >= 0.0f && fade <= 1.0f, "feedback fade %d/1000 must be in [0, 1]",
        (fade > -1.0e6f && fade < 1.0e6f) ? static_cast<int>(fade * 1000.0f)
                                          : static_cast<int>(INT32_MIN));
    ctx.fade = fade;
    ctx.alpha = alpha;
    feedback_style->sync_hue();
    ctx.column_time = !Spacetime::flat_surface();
    if (ctx.column_time && prior.x_offsets != ctx.end.x_offsets) {
      blend_lit_times(grid, band, ctx.end, prior);
      ctx.warp = prior;
    }
    if (ctx.column_time)
      ctx.fades = fade_table.acquire(*feedback_style, scratch);
    const bool hue_fade = feedback_style->color_fn == &::Feedback::hue_fade;
    ctx.black_skips_color = hue_fade;
    ctx.pole_half_res = feedback_style->pole_half_res;
    if (grid.polar.rows() > 0 && !math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    ctx.previous = cv.prev_data();
    ctx.current = cv.data();
    select_poles(ctx.previous, ctx.poles);
    ctx.runs = make_column_runs(band.x_clip, ctx.column_time);
    RowGeometry *rows = scratch.allocate_n<RowGeometry>(H);
    describe_rows(ctx, rows);
    ctx.rows = rows;
    const bool clipped =
        band.x_clip.active || band.y_begin > 0 || band.y_end < H;
    if (ctx.column_time && cv.keeps_stale_pixels() && alpha >= 1.0f &&
        !clipped) {
      ctx.back = ctx.current;
      select_poles(ctx.back, ctx.back_poles);
      // An unchanged cached field keeps its spans.
      const bool stored = reused && spans.reach;
      if (!(stored && spans.valid &&
            spans.pole_half_res == ctx.pole_half_res)) {
        int16_t *reach = stored ? spans.reach : scratch.allocate_n<int16_t>(H);
        measure_reach(grid, ctx.warp, ctx.end, rows, reach);
        spans.valid = stored;
        spans.pole_half_res = ctx.pole_half_res;
        ctx.reach = reach;
      } else {
        ctx.reach = spans.reach;
      }
      if (split) {
        ctx.cover = scratch.allocate_n<uint16_t>(COVER_SLOTS);
        std::memset(ctx.cover, 0, COVER_SLOTS * sizeof(uint16_t));
        ctx.ahead_span = scratch.allocate_n<AheadSpan>(H);
        for (int y = 0; y < H; ++y)
          ctx.ahead_span[y] = {static_cast<int16_t>(HALF), -1};
      }
    }

    const bool hue_identity =
        feedback_style->hue_ca == 1.0f && feedback_style->hue_sa == 0.0f;
    if (hue_fade && !hue_identity) {
      ::Feedback::scale_hue_matrix(feedback_style->hue_k, fade, ctx.hue_k);
      ctx.mode = FlushMode::HUE;
    } else {
      ctx.mode = hue_fade ? FlushMode::PLAIN : FlushMode::GENERAL;
    }
    return ctx;
  }

  /**
   * @brief Moves each node of a moving field to its column's lit time: an
   *        arm's first column reads the prior field, its end the current one.
   * @param grid Coarse layout of the frame.
   * @param band Rows and columns the frame renders.
   * @param now The current field.
   * @param prior The field one noise step before it; receives the blend.
   */
  HS_HOT_FLASH_MEMBER static void blend_lit_times(const CoarseGrid &grid,
                                                  const RenderBand &band,
                                                  const WarpField &now,
                                                  const WarpField &prior) {
    for (int field_y = band.field_y_begin; field_y <= band.field_y_end;
         ++field_y) {
      const bool polar =
          field_y < grid.polar.north_rings || field_y >= grid.polar.south_ring;
      for (int cx = 0; cx < grid.columns; ++cx) {
        const int column = cx * grid.downsample;
        const float lag = 1.0f - (column - ColumnTime::arm_start(column)) *
                                     ColumnTime::INVERSE_HALF;
        const int i = field_y * grid.columns + cx;
        const float x = now.x_offsets[i];
        const float from =
            unwrap_near(static_cast<float>(prior.x_offsets[i]), x, WRAP_PERIOD);
        prior.x_offsets[i] = static_cast<int16_t>(
            unwrap_near(x + lag * (from - x), 0.0f, WRAP_PERIOD));
        prior.y_offsets[i] = static_cast<int16_t>(
            now.y_offsets[i] + lag * (prior.y_offsets[i] - now.y_offsets[i]));
        if (polar) {
          const int cell = grid.polar.row(field_y) * grid.columns + cx;
          prior.cell_caps[cell] = CapPlane::lerp_offset(
              now.cell_caps[cell], prior.cell_caps[cell], lag);
        }
      }
    }
  }

  /**
   * @brief Resolves each row's ring pair, weights and reconstruction.
   * @param ctx Frame context with its grid, band and alpha set.
   * @param rows Receives H entries; rows outside the band stay zero.
   */
  HS_FLASH_MEMBER static void describe_rows(const FlushContext &ctx,
                                            RowGeometry *rows) {
    const CoarseGrid &grid = ctx.grid;
    const RenderBand &band = ctx.band;
    const PolarRings &polar = grid.polar;
    const bool opaque = ctx.alpha >= 1.0f;
    const float half_res_sine = half_res_limit(ctx.pole_half_res);
    std::memset(rows, 0, H * sizeof(RowGeometry));
    int field_y0 = band.field_y_begin;
    int field_y1 = field_y0 + (field_y0 < band.field_y_end ? 1 : 0);
    auto ring0 = grid.field.ring(field_y0);
    auto ring1 = field_y1 > field_y0 ? grid.field.next_ring(ring0) : ring0;
    int window_base = 0;
    int cover_base = 0;
    for (int y = band.y_begin; y < band.y_end; ++y) {
      while (y > ring1.y && field_y1 < band.field_y_end) {
        field_y0 = field_y1;
        ring0 = ring1;
        ++field_y1;
        ring1 = grid.field.next_ring(ring1);
      }
      // no death case: the ring chain ends on the band's last row
      HS_CHECK(y <= ring1.y, "feedback warp row %d past last control row %d", y,
               ring1.y);
      RowGeometry &g = rows[y];
      g.field_y0 = static_cast<int16_t>(field_y0);
      g.field_y1 = static_cast<int16_t>(field_y1);
      const int height = ring1.y - ring0.y;
      g.wy1 = height > 0 ? static_cast<float>(y - ring0.y) / height : 0.0f;
      g.polar = field_y1 < polar.north_rings || field_y0 >= polar.south_ring;
      g.north = field_y1 < polar.north_rings;
      g.pole = is_pole_row(y);
      g.infill = !band.x_clip.active && filtered_infill(grid, y);
      // The pair expansion needs plain stores, and a longitude-filtered row
      // is reconstructed at its footprint already.
      g.pairs = opaque && !g.infill && !g.pole &&
                SphereField::latitude_sine(y) < half_res_sine;
      const int width = g.infill ? grid.field.longitude_filter_width(y) : 1;
      // no death case: the bound widens longitude_filter_width()
      HS_CHECK(width <= filter_window_bound(y),
               "feedback infill row %d filter width %d exceeds its bound %d", y,
               width, filter_window_bound(y));
      g.reach = static_cast<uint8_t>(g.infill ? width / 2 : g.pairs ? 2 : 0);
      g.window_base = static_cast<int16_t>(window_base);
      window_base += 2 * finish_window(g);
      g.wide = SphereField::latitude_sine(y) < POLAR_TARGET_SINE;
      g.cover_base = static_cast<uint16_t>(cover_base);
      cover_base += 2 * cover_columns(g);
    }
    // no death case: COVER_ROWS bounds every display geometry's wide rows
    HS_CHECK(cover_base <= COVER_SLOTS,
             "feedback coverage needs %d slots of %d", cover_base, COVER_SLOTS);
  }

  /**
   * @brief Composites the warped previous frame through a color transform.
   * @tparam PAIR_PIXELS Whether to interleave two pixels per step through
   * @p transform_pair.
   * @param ctx Frame context from prepare_flush().
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param transform_pair Maps two sampled RGB triples sharing one level to
   *        two output pixels.
   * @param starts Whether to composite only the rows' lanes short of their
   *        reach, on rows reaching no further than a quarter turn, as
   *        between sweep_deep()'s halves; otherwise every row's lanes from its
   *        reach on.
   */
  template <bool PAIR_PIXELS, typename TransformPixelT, typename TransformPairT>
  void composite_pixels(const FlushContext &ctx,
                        TransformPixelT &&transform_pixel,
                        TransformPairT &&transform_pair, bool starts = false) {
    const CoarseGrid &grid = ctx.grid;
    const RenderBand &band = ctx.band;
    const int downsample = grid.downsample;
    const int coarse_columns = grid.columns;
    const int row_begin = band.y_begin;
    const int row_end = band.y_end;
    const float inverse_downsample = 1.0f / grid.downsample;
    const bool black_skips_color = ctx.black_skips_color;
    const auto blend = blend_alpha(ctx.alpha);
    const bool opaque = ctx.alpha >= 1.0f;
    const ::Pixel *previous = ctx.previous;
    ::Pixel *current = ctx.current;
    ::Pixel *filtered_row = ctx.filtered_row;
    const auto &poles = ctx.poles;
    const bool column_time = ctx.column_time;
    const int unit_level = column_time ? FadeTable::UNIT_LEVEL : 0;
    // Cells crossing an arm cut resolve per pixel.
    const bool straddled_cells = HALF % downsample != 0;
    constexpr float LEVELS_PER_COLUMN =
        FadeTable::STEPS_PER_FRAME * ColumnTime::INVERSE_HALF;
    for (int y = row_begin; y < row_end; ++y) {
      const RowGeometry &g = ctx.rows[y];
      const int row = y * W;
      if (g.pole) {
        if (!starts)
          composite_pole_row(ctx, y, g, transform_pixel);
        continue;
      }
      const bool defer_filter = g.infill && !opaque;
      ::Pixel *output = defer_filter ? filtered_row : current + row;
      const bool half_res = g.pairs;
      const float wy1 = g.wy1, wy0 = 1.0f - wy1;
      const int row0 = g.field_y0 * coarse_columns;
      const int row1 = g.field_y1 * coarse_columns;

      const ColumnRuns row_runs = row_lanes(ctx, g, y, starts);
      for (int r = 0; r < row_runs.count; ++r) {
        const int xs = row_runs.items[r].begin;
        const int xe = row_runs.items[r].end;
        // A half-resolution sample sits between its column pair, so it is the
        // pair's box average rather than the even column alone.
        const bool pairs = half_res && pair_lane(ctx, g, xs);
        const int stride = pairs ? 2 : 1;
        const float lane_offset = pairs ? 0.5f : 0.0f;
        if (g.polar) {
          composite_polar_row(ctx, xs, xe, y, g, output, stride, lane_offset,
                              defer_filter, transform_pixel);
          continue;
        }
        {
          int cx0 = xs / downsample;
          int sub = xs - cx0 * downsample;
          bool cell_stale = true;
          // Whether a pixel of the cell may read past its arm's columns.
          bool cell_resolved = false;
          // The buffer holding every source of an unresolved cell.
          const ::Pixel *source = previous;
          const auto *source_poles = &poles;
          int cell_level = unit_level;
          CellOffsets offsets{};
          auto cell = [&]() {
            HS_PROFILE_DEEP(fb_comp_cell);
            const int cx1 = (cx0 + 1 < coarse_columns) ? cx0 + 1 : 0;
            int reach = 0;
            offsets = cell_offsets(ctx, row0, row1, cx0, cx1, wy0, wy1, reach);
            if (column_time) {
              int column = cx0 * downsample;
              column -= column >= HALF ? HALF : 0;
              // No pixel's offset exceeds its cell's largest node offset, so
              // a cell at least that far inside its arm reads within it.
              const int room = std::min(column, HALF - 1 - downsample - column);
              cell_resolved = straddled_cells;
              source = previous;
              source_poles = &poles;
              if (reach > room * static_cast<int>(WARP_SCALE)) {
                // Otherwise the first pixel's source and the last's, past
                // the cell's far edge, bound the cell's.
                const float first =
                    column +
                    fminf(offsets.left_x, offsets.left_x + offsets.slope_x);
                const float last =
                    column + downsample +
                    fmaxf(offsets.left_x, offsets.left_x + offsets.slope_x);
                if (ctx.back && last <= -1.0f &&
                    first - column >= -static_cast<float>(HALF)) {
                  source = ctx.back;
                  source_poles = &ctx.back_poles;
                } else {
                  cell_resolved |= !(first >= 0.0f && last < HALF - 1.0f);
                }
              }
              cell_level =
                  static_cast<int>(FadeTable::UNIT_LEVEL + 0.5f -
                                   (offsets.left_x + 0.5f * offsets.slope_x) *
                                       LEVELS_PER_COLUMN);
              HS_AUDIT_CHECK(cell_resolved || (cell_level >= 0 &&
                                               cell_level < FadeTable::LEVELS),
                             "feedback cell fade level out of range");
            }
          };
          for (int x = xs; x < xe;) {
            if (cell_stale) {
              cell();
              cell_stale = false;
            }

            if (cell_resolved) {
              const int cell_end = std::min(xe, x + downsample - sub);
              x = composite_resolved_cell(ctx, transform_pixel, output, x,
                                          cell_end, y, sub, stride, lane_offset,
                                          offsets, opaque || defer_filter);
              sub = x - cx0 * downsample;
              while (sub >= downsample) {
                sub -= downsample;
                ++cx0;
                cell_stale = true;
              }
              continue;
            }

            if constexpr (PAIR_PIXELS) {
              if (downsample - sub >= 2 * stride && xe - x >= 2 * stride) {
                const float fx0 = (sub + lane_offset) * inverse_downsample;
                const float fx1 =
                    (sub + stride + lane_offset) * inverse_downsample;
                const float ddx0 = offsets.left_x + offsets.slope_x * fx0;
                const float ddy0 = offsets.left_y + offsets.slope_y * fx0;
                const float ddx1 = offsets.left_x + offsets.slope_x * fx1;
                const float ddy1 = offsets.left_y + offsets.slope_y * fx1;

                Sample s0, s1;
                {
                  HS_PROFILE_DEEP(fb_comp_sample);
                  s0 = sample_bilinear_prev(grid.field, source, *source_poles,
                                            x + lane_offset + ddx0, y + ddy0);
                  s1 = sample_bilinear_prev(grid.field, source, *source_poles,
                                            x + stride + lane_offset + ddx1,
                                            y + ddy1);
                }
                const bool black0 = skip_near_black(s0, black_skips_color);
                const bool black1 = skip_near_black(s1, black_skips_color);
                ::Pixel p0(0, 0, 0), p1(0, 0, 0);
                if (!(black0 && black1)) {
                  HS_PROFILE_DEEP(fb_comp_color);
                  transform_pair_call(transform_pair, s0.r, s0.g, s0.b, s1.r,
                                      s1.g, s1.b, cell_level, p0, p1);
                  p0 = black0 ? ::Pixel(0, 0, 0) : p0;
                  p1 = black1 ? ::Pixel(0, 0, 0) : p1;
                }

                HS_PROFILE_DEEP(fb_comp_write);
                ::Pixel &dst0 = output[x];
                dst0 = (opaque || defer_filter) ? p0 : blend(dst0, p0);
                ::Pixel &dst1 = output[x + stride];
                dst1 = (opaque || defer_filter) ? p1 : blend(dst1, p1);

                x += 2 * stride;
                sub += 2 * stride;
                if (sub >= downsample) {
                  sub -= downsample;
                  ++cx0;
                  cell_stale = true;
                }
                continue;
              }
            }

            const float fx = (sub + lane_offset) * inverse_downsample;
            const float ddx = offsets.left_x + offsets.slope_x * fx;
            const float ddy = offsets.left_y + offsets.slope_y * fx;

            Sample s;
            {
              HS_PROFILE_DEEP(fb_comp_sample);
              s = sample_bilinear_prev(grid.field, source, *source_poles,
                                       x + lane_offset + ddx, y + ddy);
            }
            ::Pixel p(0, 0, 0);
            if (!(skip_near_black(s, black_skips_color))) {
              HS_PROFILE_DEEP(fb_comp_color);
              p = transform_sample<PAIR_PIXELS>(transform_pixel, s, cell_level);
            }

            // Black must overwrite the stale double-buffer frame.
            HS_PROFILE_DEEP(fb_comp_write);
            ::Pixel &dst = output[x];
            dst = (opaque || defer_filter) ? p : blend(dst, p);

            x += stride;
            sub += stride;
            while (sub >= downsample) {
              sub -= downsample;
              ++cx0;
              cell_stale = true;
            }
          }
        }
      }
      if (starts)
        continue;
      if (half_res) {
        HS_PROFILE_DEEP(fb_comp_fill);
        expand_half_res_row(ctx, output, y);
      }
      if (g.infill) {
        HS_PROFILE_DEEP(fb_comp_filter);
        filter_infill_row(ctx, y, defer_filter);
      }
    }
  }

  /**
   * @brief Composites a pole row: one read of the pole's source a frame
   *        earlier, written to every column.
   * @param ctx Frame context from prepare_flush().
   * @param y The pole row.
   * @param g The row's geometry.
   * @param transform_pixel Maps a sample and its fade-table level to a pixel.
   */
  HS_FLASH_MEMBER void composite_pole_row(const FlushContext &ctx, int y,
                                          const RowGeometry &g,
                                          PixelTransform transform_pixel) {
    const CoarseGrid &grid = ctx.grid;
    // no death case: a pole's latitude sine is zero
    HS_CHECK(g.polar, "feedback pole row %d lies outside the polar rings", y);
    const CapCell cell = polar_cell(ctx, g, 0, 0);
    const auto at = CapPlane::polar_lane(cell, 0, 0.0f, g.north, false, 0.0f);
    const int level = ctx.column_time ? FadeTable::UNIT_LEVEL : 0;
    const Sample s =
        sample_bilinear_prev(grid.field, ctx.previous, ctx.poles, at.x, at.y);
    ::Pixel p(0, 0, 0);
    if (!skip_near_black(s, ctx.black_skips_color))
      p = transform_pixel(s.r, s.g, s.b, level);
    ::Pixel *output = ctx.current + y * W;
    const auto blend = blend_alpha(ctx.alpha);
    const bool plain_store = ctx.alpha >= 1.0f;
    for (int r = 0; r < ctx.runs.count; ++r)
      for (int x = ctx.runs.items[r].begin; x < ctx.runs.items[r].end; ++x)
        output[x] = plain_store ? p : blend(output[x], p);
  }

  /**
   * @brief Reconstructs a dense infill row at its spherical footprint, under
   *        column time within each arm.
   * @param ctx Frame context from prepare_flush().
   * @param y Row, composited into the current frame or, when @p deferred,
   *        into the filter row.
   * @param deferred Whether the row was composited into the filter row.
   */
  HS_HOT_FLASH_MEMBER void filter_infill_row(const FlushContext &ctx, int y,
                                             bool deferred) {
    ::Pixel *current = ctx.current + y * W;
    if (!deferred)
      std::memcpy(ctx.filtered_row, current, W * sizeof(::Pixel));
    const bool opaque = ctx.alpha >= 1.0f;
    const uint16_t alpha = frac_to_q16(ctx.alpha);
    auto store = [&](int x, const ::Pixel &pixel) HS_HOT_INLINE {
      ::Pixel &dst = current[x];
      dst = opaque ? pixel : dst.lerp16(pixel, alpha);
    };
    if (ctx.column_time)
      ctx.grid.field.template reconstruct_longitude_arms<PixelAccumulator>(
          ctx.filtered_row, y, store);
    else
      ctx.grid.field.template reconstruct_longitude_row<PixelAccumulator>(
          ctx.filtered_row, y, store);
  }

  /**
   * @brief Composites one polar row's lanes from their cap-plane offsets.
   * @details Each lane converts its own target back to field coordinates, one
   * lane per step.
   */
  HS_HOT_FLASH_MEMBER void
  composite_polar_row(const FlushContext &ctx, int xs, int xe, int y,
                      const RowGeometry &g, ::Pixel *output, int stride,
                      float lane_offset, bool defer_filter,
                      PixelTransform transform_pixel) {
    const CoarseGrid &grid = ctx.grid;
    const int downsample = grid.downsample;
    const int coarse_columns = grid.columns;
    const float inverse_downsample = 1.0f / downsample;
    const float cap_angle = cap_angle_of(y, g.north);
    const bool half_res = stride == 2;
    const bool black_skips_color = ctx.black_skips_color;
    const auto blend = blend_alpha(ctx.alpha);
    const bool plain_store = ctx.alpha >= 1.0f || defer_filter;
    {
      int cx0 = xs / downsample;
      int sub = xs - cx0 * downsample;
      bool cell_stale = true;
      CapCell cell{};
      for (int x = xs; x < xe;) {
        if (cell_stale) {
          const int cx1 = (cx0 + 1 < coarse_columns) ? cx0 + 1 : 0;
          cell = polar_cell(ctx, g, cx0, cx1);
          cell_stale = false;
        }

        const float fx = (sub + lane_offset) * inverse_downsample;
        const auto at =
            CapPlane::polar_lane(cell, x, cap_angle, g.north, half_res, fx);
        ::Pixel &dst = output[x];
        const float target = x + lane_offset;
        const int arm = ColumnTime::arm_start(x);
        Lane lane{ColumnTime::read(target, at.x, at.y), Side::WITHIN};
        if (ctx.column_time)
          lane = resolve_lane(target, arm, at.x, at.y);
        if (lane.side == Side::AHEAD && ctx.cover) {
          hold_lane(ctx, y, arm, x);
          dst = ::Pixel(0, 0, 0);
        } else {
          const Sample s = sample_lane(ctx, arm, lane);
          ::Pixel p(0, 0, 0);
          if (!skip_near_black(s, black_skips_color))
            p = transform_pixel(
                s.r, s.g, s.b,
                ctx.column_time ? FadeTable::level(lane.read.elapsed) : 0);
          dst = plain_store ? p : blend(dst, p);
        }

        x += stride;
        sub += stride;
        while (sub >= downsample) {
          sub -= downsample;
          ++cx0;
          cell_stale = true;
        }
      }
    }
  }

  /** @brief A polar ring's cap offsets in @p field, indexed by coarse
   *  column. */
  static HS_HOT_INLINE const CapOffset *
  polar_caps(const WarpField &field, const FlushContext &ctx, int field_y) {
    return field.cell_caps + ctx.grid.polar.row(field_y) * ctx.grid.columns;
  }

  /** @brief Whether a coarse column is an arm's first, whose node a cell
   *  ending there reads from the field at the arm's end. */
  static HS_HOT_INLINE bool arm_edge(const FlushContext &ctx, int cx) {
    return cx * ctx.grid.downsample % HALF == 0;
  }

  /** @brief A polar row's cell, its nodes at their lit times. */
  static HS_HOT_INLINE CapCell polar_cell(const FlushContext &ctx,
                                          const RowGeometry &g, int cx0,
                                          int cx1) {
    const CapOffset *caps0 = polar_caps(ctx.warp, ctx, g.field_y0);
    const CapOffset *caps1 = polar_caps(ctx.warp, ctx, g.field_y1);
    const WarpField &right = arm_edge(ctx, cx1) ? ctx.end : ctx.warp;
    return CapPlane::decode_cell(
        caps0[cx0], caps1[cx0], polar_caps(right, ctx, g.field_y0)[cx1],
        polar_caps(right, ctx, g.field_y1)[cx1], 1.0f - g.wy1, g.wy1);
  }

  /** @brief A row's angle from its pole, in radians. */
  static HS_HOT_INLINE float cap_angle_of(int y, bool north) {
    const float colatitude =
        SphereField::Geometry::row_to_phi(static_cast<float>(y));
    return north ? colatitude : math::PI_F - colatitude;
  }

  /** @brief A cell's offsets: at its left edge and their change across it,
   *  in pixels. */
  struct CellOffsets {
    float left_x;  ///< Column offset at the left edge.
    float slope_x; ///< Column offset change across the cell.
    float left_y;  ///< Row offset at the left edge.
    float slope_y; ///< Row offset change across the cell.
  };

  /**
   * @brief A cell's offsets blended between its two rings, its nodes at
   *        their lit times.
   * @param ctx Frame context.
   * @param row0,row1 First cell index of each ring.
   * @param cx0,cx1 The cell's left and right coarse columns.
   * @param wy0,wy1 The rings' weights.
   * @param reach Receives the largest node column offset, WARP_SCALE units.
   * @return The cell's offsets.
   */
  static HS_HOT_INLINE CellOffsets cell_offsets(const FlushContext &ctx,
                                                int row0, int row1, int cx0,
                                                int cx1, float wy0, float wy1,
                                                int &reach) {
    constexpr float INVERSE_WARP_SCALE = 1.0f / WARP_SCALE;
    const WarpField &right = arm_edge(ctx, cx1) ? ctx.end : ctx.warp;
    const int i00 = row0 + cx0, i10 = row0 + cx1;
    const int i01 = row1 + cx0, i11 = row1 + cx1;
    const int d00 = ctx.warp.x_offsets[i00];
    const int u10 = unwrap_near(right.x_offsets[i10], d00);
    const int u01 = unwrap_near(ctx.warp.x_offsets[i01], d00);
    const int u11 = unwrap_near(right.x_offsets[i11], d00);
    reach = std::max(std::max(std::abs(d00), std::abs(u10)),
                     std::max(std::abs(u01), std::abs(u11)));
    CellOffsets o;
    o.left_x = (static_cast<float>(d00) * wy0 + static_cast<float>(u01) * wy1) *
               INVERSE_WARP_SCALE;
    o.slope_x =
        (static_cast<float>(u10) * wy0 + static_cast<float>(u11) * wy1) *
            INVERSE_WARP_SCALE -
        o.left_x;
    o.left_y = (ctx.warp.y_offsets[i00] * wy0 + ctx.warp.y_offsets[i01] * wy1) *
               INVERSE_WARP_SCALE;
    o.slope_y = (right.y_offsets[i10] * wy0 + right.y_offsets[i11] * wy1) *
                    INVERSE_WARP_SCALE -
                o.left_y;
    return o;
  }

  /**
   * @brief Resolves a lane's read from its target and source.
   * @param target Target column, with its half-resolution lane offset.
   * @param arm First column of the target's arm.
   * @param source_x Source column, on any wrap branch.
   * @param source_y Source row.
   * @return The read and the side of the arm its source lies on.
   * @details The source unwraps to the shortest path, under a half turn ahead
   * and no further behind than keeps its taps past the other arm's lane of
   * the same column.
   */
  static HS_HOT_INLINE Lane resolve_lane(float target, int arm, float source_x,
                                         float source_y) {
    float delta = unwrap_near(source_x - target, 0.0f, static_cast<float>(W));
    if (delta >= static_cast<float>(HALF))
      delta -= static_cast<float>(W);
    delta = fmaxf(delta, -static_cast<float>(HALF - 1));
    const float column = (target - arm) + delta;
    const Side side = column < 0.0f           ? Side::BEHIND
                      : column >= HALF - 1.0f ? Side::AHEAD
                                              : Side::WITHIN;
    return {ColumnTime::read(target, target + delta, source_y), side};
  }

  /**
   * @brief Samples a resolved lane from the buffers holding its taps.
   * @param ctx Frame context.
   * @param arm First column of the target's arm.
   * @param lane The lane; AHEAD only in a frame without coverage. Its elapsed
   *        follows the buffer read.
   * @return The bilinear sample.
   * @details A source across a cut reads the frame holding it when the frame
   * has one, else the previous frame, a frame older past the cut and a frame
   * younger behind it.
   */
  HS_HOT_INLINE Sample sample_lane(const FlushContext &ctx, int arm,
                                   Lane &lane) const {
    const bool across = lane.side == Side::BEHIND  ? ctx.back != nullptr
                        : lane.side == Side::AHEAD ? ctx.cover != nullptr
                                                   : false;
    if (across)
      return sample_across(ctx, arm, lane.read.x, lane.read.y);
    if (lane.side != Side::WITHIN)
      lane.read.elapsed += lane.side == Side::AHEAD ? 1.0f : -1.0f;
    return sample_bilinear_prev(ctx.grid.field, ctx.previous, ctx.poles,
                                lane.read.x, lane.read.y);
  }

  /**
   * @brief Marks a lane whose source reads the frame as plotted, and opens
   *        the coverage of every column whose reconstruction reads it.
   * @param ctx Frame context with its coverage allocated.
   * @param y Row.
   * @param arm First column of the lane's arm.
   * @param x The lane's column.
   */
  HS_HOT_FLASH_MEMBER void hold_lane(const FlushContext &ctx, int y, int arm,
                                     int x) {
    const RowGeometry &g = ctx.rows[y];
    const int k = x - arm;
    auto open = [&](int column, uint16_t bits) {
      const int slot = cover_slot(g, arm, arm + column);
      // no death case: coverage spans every lane a capped offset can hold
      HS_CHECK(slot >= 0,
               "feedback row %d column %d is held outside its coverage", y,
               arm + column);
      ctx.cover[slot] |= bits;
      held_floor = std::min(held_floor, column);
      row_floor[y] = static_cast<uint8_t>(std::min<int>(row_floor[y], column));
    };
    open(k, AHEAD_BIT);
    AheadSpan &span = ctx.ahead_span[y];
    span.first = static_cast<int16_t>(std::min<int>(span.first, k));
    span.last = static_cast<int16_t>(std::max<int>(span.last, k));
    int lo = k, hi = k;
    if (g.pairs && k >= HALF - 2) {
      open(HALF - 3, TRANSMITTANCE_BITS);
    } else if (g.pairs) {
      lo = k - 1;
      hi = std::min(k + 2, HALF - 3);
    } else if (g.infill) {
      lo = std::max(k - g.reach, 0);
      hi = std::min(k + g.reach, HALF - 1);
    }
    for (int column = std::max(lo, 0); column <= hi; ++column)
      if (!g.infill ||
          std::abs(column - k) <= std::min(column, HALF - 1 - column))
        open(column, TRANSMITTANCE_BITS);
  }

  /** @brief Low coverage bit: the lane's source reads the frame as plotted. */
  static constexpr uint16_t AHEAD_BIT = 1;
  /** @brief The coverage bits holding a pixel's transmittance. */
  static constexpr uint16_t TRANSMITTANCE_BITS = 0xFFFE;

  /**
   * @brief Composites a cell's pixels whose sources may lie past their arm's
   * columns, each resolved by FeedbackColumnTime.
   * @param ctx Frame context from prepare_flush().
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param output Row being written.
   * @param x First column.
   * @param end One past the last column, within the cell.
   * @param y Row.
   * @param sub @p x's column within the cell.
   * @param stride Column step.
   * @param lane_offset Half-resolution lane offset of each column.
   * @param cell The cell's elapsed-scaled offsets.
   * @param plain_store Whether pixels overwrite rather than blend.
   * @return The first column past the cell's pixels.
   * @details The colour transforms stay out of line, in the composite's own
   * placement.
   */
  HS_HOT_FLASH_MEMBER int composite_resolved_cell(
      const FlushContext &ctx, PixelTransform transform_pixel, ::Pixel *output,
      int x, int end, int y, int sub, int stride, float lane_offset,
      const CellOffsets &cell, bool plain_store) {
    const float inverse_downsample = 1.0f / ctx.grid.downsample;
    const float row = static_cast<float>(y);
    const uint16_t alpha = frac_to_q16(ctx.alpha);
    while (x < end) {
      const float fx = (sub + lane_offset) * inverse_downsample;
      const float ddx = cell.left_x + cell.slope_x * fx;
      const float ddy = cell.left_y + cell.slope_y * fx;
      const float target = x + lane_offset;
      const int arm = ColumnTime::arm_start(x);
      Lane lane = resolve_lane(target, arm, target + ddx, row + ddy);
      ::Pixel p(0, 0, 0);
      if (lane.side == Side::AHEAD && ctx.cover) {
        hold_lane(ctx, y, arm, x);
        output[x] = p;
      } else {
        const Sample s = sample_lane(ctx, arm, lane);
        if (!skip_near_black(s, ctx.black_skips_color))
          p = transform_pixel(s.r, s.g, s.b,
                              FadeTable::level(lane.read.elapsed));
        ::Pixel &dst = output[x];
        dst = plain_store ? p : dst.lerp16(p, alpha);
      }
      x += stride;
      sub += stride;
    }
    return x;
  }

  /** @brief A row's state in a column-major sweep: its pole angle and each
   *  arm's current cell. */
  struct SweepRow {
    float cap_angle;    ///< A polar row's angle from its pole.
    int16_t cell_cx[2]; ///< Each arm's current cell's coarse column, or -1.
    union {
      CellOffsets offsets; ///< A cartesian row's offsets.
      CapCell cap;         ///< A polar row's cap-plane offsets.
    } cell[2];             ///< Each arm's current cell.
  };

  /**
   * @brief Each row's sweep state, no cell decoded yet.
   * @param ctx Frame context.
   * @param arena Arena for the H states.
   * @return The states, indexed by row.
   */
  static SweepRow *sweep_rows(const FlushContext &ctx, Arena &arena) {
    SweepRow *rows = arena.allocate_n<SweepRow>(H);
    for (int y = ctx.band.y_begin; y < ctx.band.y_end; ++y) {
      const RowGeometry &g = ctx.rows[y];
      rows[y].cap_angle = g.polar ? cap_angle_of(y, g.north) : 0.0f;
      rows[y].cell_cx[0] = rows[y].cell_cx[1] = -1;
    }
    return rows;
  }

  /**
   * @brief One lane's read in a column-major sweep.
   * @param ctx Frame context.
   * @param g The lane's row geometry.
   * @param y The lane's row.
   * @param x The lane's column.
   * @param row The row's sweep state; decodes the lane's cell when the lane
   *        enters a new one.
   * @return The resolved lane.
   */
  HS_HOT_INLINE Lane sweep_lane(const FlushContext &ctx, const RowGeometry &g,
                                int y, int x, SweepRow &row) const {
    const CoarseGrid &grid = ctx.grid;
    const int downsample = grid.downsample;
    const int cx0 = x / downsample;
    const int side = x >= HALF;
    auto &cell = row.cell[side];
    if (row.cell_cx[side] != cx0) {
      row.cell_cx[side] = static_cast<int16_t>(cx0);
      const int cx1 = cx0 + 1 < grid.columns ? cx0 + 1 : 0;
      int reach = 0;
      if (g.polar)
        cell.cap = polar_cell(ctx, g, cx0, cx1);
      else
        cell.offsets = cell_offsets(ctx, g.field_y0 * grid.columns,
                                    g.field_y1 * grid.columns, cx0, cx1,
                                    1.0f - g.wy1, g.wy1, reach);
    }
    const bool pair = pair_lane(ctx, g, x);
    const float lane_offset = pair ? 0.5f : 0.0f;
    const float fx = (x - cx0 * downsample + lane_offset) / downsample;
    const float target = x + lane_offset;
    const int arm = ColumnTime::arm_start(x);
    if (g.polar) {
      const auto at =
          CapPlane::polar_lane(cell.cap, x, row.cap_angle, g.north, pair, fx);
      return resolve_lane(target, arm, at.x, at.y);
    }
    const CellOffsets &o = cell.offsets;
    return resolve_lane(target, arm, target + o.left_x + o.slope_x * fx,
                        y + o.left_y + o.slope_y * fx);
  }

  /**
   * @brief Composites the lanes of the rows reaching past a quarter turn
   *        behind their arm's backward cut, column by column.
   * @param ctx Frame context with its stale frame.
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param begin,end Columns past each arm's first to sweep.
   * @details A read behind the cut lands past the other arm's lane of the
   * same column, so lane c of a row may be composited once no row's lanes
   * past c that it may read are. These rows' first quarter goes before and
   * their second after the other rows' lanes short of their reach, which
   * never write past a quarter turn and whose reads land past it.
   */
  HS_HOT_FLASH_MEMBER void sweep_deep(const FlushContext &ctx,
                                      PixelTransform transform_pixel, int begin,
                                      int end) {
    ScratchScope scope(scratch_arena_a);
    SweepRow *rows = sweep_rows(ctx, scope.get_arena());
    for (int c = begin; c < end; ++c)
      for (int y = ctx.band.y_begin; y < ctx.band.y_end; ++y)
        if (ctx.reach[y] > HALF / 2 && c < ctx.reach[y])
          behind_lane(ctx, transform_pixel, y, c, rows[y]);
  }

  /**
   * @brief Composites one column of a row in sweep_deep().
   * @param ctx Frame context with its stale frame.
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param y Row.
   * @param c Column past each arm's first.
   * @param row The row's sweep state.
   */
  HS_HOT_FLASH_MEMBER void behind_lane(const FlushContext &ctx,
                                       PixelTransform transform_pixel, int y,
                                       int c, SweepRow &row) {
    const RowGeometry &g = ctx.rows[y];
    if (g.pairs && (c & 1) && c < HALF - 2)
      return;
    for (int arm = 0; arm < W; arm += HALF) {
      const int x = arm + c;
      Lane lane = sweep_lane(ctx, g, y, x, row);
      ::Pixel &dst = ctx.current[y * W + x];
      if (lane.side == Side::AHEAD && ctx.cover) {
        hold_lane(ctx, y, arm, x);
        dst = ::Pixel(0, 0, 0);
        continue;
      }
      const Sample s = sample_lane(ctx, arm, lane);
      ::Pixel p(0, 0, 0);
      if (!skip_near_black(s, ctx.black_skips_color))
        p = transform_pixel(s.r, s.g, s.b, FadeTable::level(lane.read.elapsed));
      dst = p;
    }
  }

  /**
   * @brief Adds every held lane's composite from the frame as plotted beneath
   *        its coverage, in lit-time order so each read lands on columns
   *        already final.
   * @param ctx Held frame context.
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @details A read stops at the last final column of its source rows, at
   * most the reconstruction's reach and one column short of the antipode.
   * Each row keeps its last reconstruction window of samples, so a completed
   * column's addition rounds once.
   */
  HS_HOT_FLASH_MEMBER void sweep_ahead(const FlushContext &ctx,
                                       PixelTransform transform_pixel) {
    const RenderBand &band = ctx.band;
    ScratchScope scope(scratch_arena_a);
    int window_slots = 0;
    for (int y = band.y_begin; y < band.y_end; ++y)
      window_slots += 2 * finish_window(ctx.rows[y]);
    ::Pixel *windows = window_slots > 0
                           ? scope.get_arena().make_n<::Pixel>(window_slots)
                           : nullptr;
    // The rows holding a lane; each span's last becomes the last lane whose
    // step completes a held column.
    AheadSpan *spans = ctx.ahead_span;
    int16_t *visit_rows = scope.get_arena().allocate_n<int16_t>(H);
    SweepRow *rows = sweep_rows(ctx, scope.get_arena());
    int count = 0, begin = HALF, end = -1;
    for (int y = band.y_begin; y < band.y_end; ++y) {
      const AheadSpan span = ctx.ahead_span[y];
      if (span.first > span.last)
        continue;
      const RowGeometry &g = ctx.rows[y];
      const int last =
          g.infill  ? std::min(HALF - 1, span.last + 2 * g.reach)
          : g.pairs ? (span.last + 2 >= HALF - 2 ? HALF - 1 : span.last + 2)
                    : span.last;
      visit_rows[count++] = static_cast<int16_t>(y);
      spans[y].last = static_cast<int16_t>(last);
      begin = std::min<int>(begin, span.first);
      end = std::max(end, last);
    }
    for (int c = begin; c <= end; ++c) {
      for (int i = 0; i < count; ++i) {
        const int y = visit_rows[i];
        if (c < spans[y].first || c > spans[y].last)
          continue;
        const RowGeometry &g = ctx.rows[y];
        if (g.pairs && (c & 1) && c < HALF - 2)
          continue;
        for (int arm = 0, side = 0; arm < W; arm += HALF, ++side) {
          const int x = arm + c;
          ::Pixel p(0, 0, 0);
          const int slot = cover_slot(g, arm, x);
          const bool held = slot >= 0 && (ctx.cover[slot] & AHEAD_BIT);
          if (!held && !g.pairs && !g.infill)
            continue;
          if (held) {
            Lane lane = sweep_lane(ctx, g, y, x, rows[y]);
            // no death case: a lane resolves as it did before the plots
            HS_CHECK(lane.side == Side::AHEAD,
                     "feedback held lane %d of row %d no longer reads ahead", x,
                     y);
            // A column is final once every lane its reconstruction reads is
            // composited: the read stops at the last such column of its
            // source rows.
            const int row0 =
                hs::clamp(static_cast<int>(std::floor(lane.read.y)), 0, H - 1);
            const int reach =
                std::max<int>(ctx.rows[row0].reach,
                              ctx.rows[std::min(row0 + 1, H - 1)].reach);
            const float limit = static_cast<float>(arm + c + HALF - 2 - reach);
            if (lane.read.x > limit)
              lane.read = ColumnTime::read(
                  x + (pair_lane(ctx, g, x) ? 0.5f : 0.0f), limit, lane.read.y);
            const Sample s = sample_across(ctx, arm, lane.read.x, lane.read.y);
            if (!skip_near_black(s, ctx.black_skips_color))
              p = transform_pixel(s.r, s.g, s.b,
                                  FadeTable::level(lane.read.elapsed));
          }
          ::Pixel *window =
              windows ? windows + g.window_base + side * finish_window(g)
                      : nullptr;
          complete_columns(ctx, y, g, arm, c, p, window);
        }
      }
    }
  }

  /**
   * @brief Adds the columns a finish sample completes beneath their coverage.
   * @param ctx Held frame context.
   * @param y Row.
   * @param g The row's geometry.
   * @param arm First column of the lane's arm.
   * @param c The lane's column past its arm's first.
   * @param sample The lane's finish composite; black for an unheld lane.
   * @param window The row and arm's recent samples, finish_window() of them.
   */
  HS_HOT_FLASH_MEMBER static void
  complete_columns(const FlushContext &ctx, int y, const RowGeometry &g,
                   int arm, int c, const ::Pixel &sample, ::Pixel *window) {
    if (g.infill) {
      const int radius = g.reach;
      const int width = finish_window(g);
      window[c % width] = sample;
      auto complete = [&](int x) {
        const int reach = std::min(radius, std::min(x, HALF - 1 - x));
        const int count = 2 * reach + 1;
        PixelAccumulator sum;
        for (int k = x - reach; k <= x + reach; ++k)
          sum.add(window[k % width]);
        add_under(ctx, y, arm, arm + x, sum.average(count));
      };
      if (c >= radius)
        complete(c - radius);
      if (c == HALF - 1)
        for (int x = std::max(0, HALF - radius); x < HALF; ++x)
          complete(x);
      return;
    }
    if (g.pairs && c < HALF - 2) {
      const ::Pixel previous = window[0];
      if (c > 0)
        add_under(ctx, y, arm, arm + c - 1, blend_quarter(previous, sample));
      add_under(ctx, y, arm, arm + c, blend_quarter(sample, previous));
      window[0] = sample;
      return;
    }
    add_under(ctx, y, arm, arm + c, sample);
    if (!g.pairs)
      return;
    // The arm's last two columns are their own lanes; the pair column before
    // them blends toward their mean.
    if (c == HALF - 2)
      window[1] = sample;
    else
      add_under(ctx, y, arm, arm + c - 2,
                blend_quarter(window[0], pair_average(window[1], sample)));
  }

  /**
   * @brief Adds a column's finish composite beneath its coverage.
   * @param ctx Held frame context.
   * @param y Row.
   * @param arm First column of the column's arm.
   * @param x Column.
   * @param under The part of the column's composite the frame as plotted
   *        supplies.
   */
  static HS_HOT_INLINE void add_under(const FlushContext &ctx, int y, int arm,
                                      int x, const ::Pixel &under) {
    const int slot = cover_slot(ctx.rows[y], arm, x);
    if (slot < 0)
      return;
    const uint16_t t = ctx.cover[slot] & TRANSMITTANCE_BITS;
    if (t == 0)
      return;
    auto add = [t](uint16_t covered, uint16_t below) {
      return static_cast<uint16_t>(
          std::min<uint32_t>(covered + lerp_q16(0, below, t), 65535u));
    };
    ::Pixel &dst = ctx.current[y * W + x];
    dst =
        ::Pixel(add(dst.r, under.r), add(dst.g, under.g), add(dst.b, under.b));
  }

  /**
   * @brief Composites every pixel from the previous frames: the lanes that may
   *        read behind a cut first, in lit-time order, then each row.
   * @param cv Target canvas.
   * @param ctx Frame context from prepare_flush().
   */
  HS_HOT_INLINE void composite_frame(Canvas &cv, FlushContext &ctx) {
    if (cv.keeps_stale_pixels() && !ctx.back)
      cv.clear_buffer();
    if (ctx.mode == FlushMode::SKIP)
      return;
    HS_PROFILE(feedback_composite);
    composite_pass(ctx, false);
  }

  /**
   * @brief Composites the frame, or finishes its held pixels, through the
   *        frame's color transform.
   * @param ctx Frame context from prepare_flush().
   * @param finish Whether to add the held pixels' part from the frame as
   *        plotted instead of compositing from the previous frames.
   */
  void composite_pass(const FlushContext &ctx, bool finish) {
    switch (ctx.mode) {
    case FlushMode::HUE: {
      const float *k = ctx.column_time ? ctx.fades.hue_k : ctx.hue_k;
      composite_phase<true>(
          ctx, finish,
          [&](float r, float g, float b, int level) {
            return ::Feedback::hue_fade_apply(k + level * FadeTable::MATRIX, r,
                                              g, b);
          },
          [&](float r0, float g0, float b0, float r1, float g1, float b1,
              int level, ::Pixel &p0, ::Pixel &p1) {
            ::Feedback::hue_fade_apply2(k + level * FadeTable::MATRIX, r0, g0,
                                        b0, r1, g1, b1, p0, p1);
          });
      break;
    }
    case FlushMode::PLAIN:
      composite_plain(ctx, finish);
      break;
    case FlushMode::GENERAL:
      composite_general(ctx, finish);
      break;
    case FlushMode::SKIP:
      break;
    }
  }

  /** @brief composite_pass() through a colour transform. */
  template <bool PAIR_PIXELS, typename TransformPixelT, typename TransformPairT>
  HS_HOT_INLINE void composite_phase(const FlushContext &ctx, bool finish,
                                     TransformPixelT &&transform_pixel,
                                     TransformPairT &&transform_pair) {
    if (finish) {
      sweep_ahead(ctx, transform_pixel);
      return;
    }
    if (ctx.back) {
      sweep_deep(ctx, transform_pixel, 0, HALF / 2);
      composite_pixels<PAIR_PIXELS>(ctx, transform_pixel, transform_pair, true);
      sweep_deep(ctx, transform_pixel, HALF / 2, HALF);
    }
    composite_pixels<PAIR_PIXELS>(ctx, transform_pixel, transform_pair);
  }

  /** @brief Composites through the unrotated hue fade: a per-channel scale. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_plain(const FlushContext &ctx, bool finish) {
    const float *fades = ctx.column_time ? ctx.fades.fades : &ctx.fade;
    auto plain = [&](float r, float g, float b, int level) {
      const float fade = fades[level];
      return ::Pixel(quantize16(r * fade), quantize16(g * fade),
                     quantize16(b * fade));
    };
    composite_phase<false>(ctx, finish, plain, plain);
  }

  /** @brief Composites through the Style's type-erased color transform. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_general(const FlushContext &ctx, bool finish) {
    const float *fades = ctx.column_time ? ctx.fades.fades : &ctx.fade;
    auto general = [&](float r, float g, float b, int level) {
      return feedback_style->color_fn(
          ::Pixel(quantize16(r), quantize16(g), quantize16(b)), fades[level],
          *feedback_style);
    };
    composite_phase<false>(ctx, finish, general, general);
  }
  HS_O3_END

  // Round-to-nearest fade stalls dim channels short of black.
  static constexpr float NEAR_BLACK = 64.0f;
  static constexpr float WARP_SCALE = 128.0f;
  /** @brief Column offsets in WARP_SCALE units, one full turn apart. */
  static constexpr float WRAP_PERIOD = static_cast<float>(W) * WARP_SCALE;

  /**
   * @brief Whether a column of a half-resolution row stands for its pair.
   * @param ctx Frame context.
   * @param g The row's geometry.
   * @param x Column.
   * @return True except, under column time, for each arm's last two columns,
   *         which are their own lanes.
   */
  static HS_HOT_INLINE bool pair_lane(const FlushContext &ctx,
                                      const RowGeometry &g, int x) {
    return g.pairs &&
           !(ctx.column_time && x - ColumnTime::arm_start(x) >= HALF - 2);
  }

  /**
   * @brief A row's lanes for the row composite.
   * @param ctx Frame context.
   * @param g The row's geometry.
   * @param y Row.
   * @param starts Whether the lanes short of the row's reach are wanted, on a
   *        row reaching no further than a quarter turn.
   * @return The band's runs, or with a stale frame each arm's columns from
   *         the row's reach to its end; a half-resolution row's runs under
   *         column time end short of each arm's last two columns, which form
   *         runs of their own.
   */
  static ColumnRuns row_lanes(const FlushContext &ctx, const RowGeometry &g,
                              int y, bool starts) {
    ColumnRuns runs{};
    if (starts) {
      if (ctx.reach[y] <= HALF / 2)
        for (int arm = 0; arm < W; arm += HALF)
          runs.items[runs.count++] = {arm, arm + ctx.reach[y]};
    } else if (!ctx.back) {
      runs = ctx.runs;
    } else if (ctx.reach[y] < HALF) {
      for (int arm = 0; arm < W; arm += HALF)
        runs.items[runs.count++] = {arm + ctx.reach[y], arm + HALF};
    }
    if (!(ctx.column_time && g.pairs))
      return runs;
    ColumnRuns split{};
    for (int r = 0; r < runs.count; ++r) {
      const ColumnRun run = runs.items[r];
      const int tail = ColumnTime::arm_start(run.end - 1) + HALF - 2;
      if (run.begin < tail && run.end > tail) {
        split.items[split.count++] = {run.begin, tail};
        split.items[split.count++] = {tail, run.end};
      } else {
        split.items[split.count++] = run;
      }
    }
    return split;
  }

  /**
   * @brief Each row's reach behind the arms' backward cuts.
   * @param grid Coarse layout of the frame.
   * @param warp The frame's populated warp field.
   * @param end The field at each arm's end.
   * @param rows The frame's row geometry.
   * @param reach Receives, per row, the even column count from each arm's
   *        first whose lanes' sources may lie behind its backward cut; zero
   *        on a pole row.
   * @details A cell's pixels read no further than its nodes' extreme column
   * offsets in either field; a polar lane's source turns no further than its
   * cell's longest cap offset subtends, and a cap offset reaching the pole may
   * turn it anywhere.
   */
  HS_FLASH_MEMBER void measure_reach(const CoarseGrid &grid,
                                     const WarpField &warp,
                                     const WarpField &end,
                                     const RowGeometry *rows,
                                     int16_t *reach) const {
    constexpr float COLUMNS_PER_RADIAN = W / (2.0f * math::PI_F);
    // Columns between a cell's bound and a cut.
    constexpr float MARGIN = 1.0f;
    const PolarRings &polar = grid.polar;
    const int columns = grid.columns;
    const int downsample = grid.downsample;
    // Columns from each arm's first that the rows interpolating rings k and
    // k + 1 may read behind it: whole cells whose sources may cross the cut.
    auto pair_reach = [&](int k, int y0, int y1) HS_FLASH_MEMBER {
      const bool north = k + 1 < polar.north_rings;
      const bool polar_pair = north || k >= polar.south_ring;
      const float near =
          north ? SphereField::Geometry::row_to_phi(static_cast<float>(y0))
                : math::PI_F -
                      SphereField::Geometry::row_to_phi(static_cast<float>(y1));
      const int base = (polar_pair ? polar.row(k) : k) * columns;
      // The largest backward turn of `count` nodes; false when a polar node's
      // offset reaches the pole.
      auto extent = [&](int count, auto &&node, float &back) HS_HOT_INLINE {
        if (polar_pair) {
          int longest_sq = 0;
          for (const WarpField *field : {&warp, &end})
            for (int i = 0; i < count; ++i) {
              const CapOffset &o = field->cell_caps[node(i)];
              longest_sq = std::max(longest_sq, o.u * o.u + o.v * o.v);
            }
          const float longest =
              sqrtf(static_cast<float>(longest_sq)) / CapPlane::CAP_SCALE;
          // The turn's tangent bounds the turn.
          back = longest / sqrtf(near * near - longest * longest) *
                 COLUMNS_PER_RADIAN;
          return longest < near;
        }
        int low = 0;
        for (const WarpField *field : {&warp, &end})
          for (int i = 0; i < count; ++i)
            low = std::min<int>(low, field->x_offsets[node(i)]);
        back = -low / WARP_SCALE;
        return true;
      };
      float most_back = 0.0f;
      if (!extent(2 * columns, [&](int i) { return base + i; }, most_back))
        return HALF;
      int span = 0;
      for (int cx = 0; cx < columns; ++cx) {
        // Only cells within the pair's largest turn of a cut can join.
        const int column = cx * downsample % HALF;
        if (column - most_back >= MARGIN)
          continue;
        const int cx1 = cx + 1 < columns ? cx + 1 : 0;
        float back = 0.0f;
        extent(
            4,
            [&](int i) {
              return base + (i >> 1) * columns + (i & 1 ? cx1 : cx);
            },
            back);
        for (int arm = 0; arm < W; arm += HALF) {
          const int first = cx * downsample - arm;
          if (first + downsample <= 0 || first >= HALF)
            continue;
          if (first - back < MARGIN)
            span = std::max(span, first + downsample);
        }
      }
      return span;
    };
    auto even = [](int span) { return std::min(HALF, (span + 1) & ~1); };
    int pair_index = 0;
    auto lower_ring = grid.field.ring(0);
    auto upper_ring = grid.field.next_ring(lower_ring);
    int pair = pair_reach(0, lower_ring.y, upper_ring.y);
    for (int y = 0; y < H; ++y) {
      if (y > upper_ring.y && upper_ring.y < H - 1) {
        lower_ring = upper_ring;
        upper_ring = grid.field.next_ring(lower_ring);
        pair = pair_reach(++pair_index, lower_ring.y, upper_ring.y);
      }
      reach[y] = static_cast<int16_t>(rows[y].pole ? 0 : even(pair));
    }
  }

  /** @brief Latitude sine under which two columns subtend less than a
   *         number of row pitches.
   *  @param pole_half_res Row pitches (Feedback::Style::pole_half_res).
   *  @return The sine; rows below it composite column pairs. */
  static __attribute__((always_inline)) float
  half_res_limit(float pole_half_res) {
    return pole_half_res * W * SphereField::Geometry::RADIANS_PER_ROW *
           (1.0f / (4.0f * math::PI_F));
  }

  /** @brief Whether an unclipped row is an infill row reconstructed at its
   *         spherical longitude footprint.
   *  @param grid Coarse layout of the frame.
   *  @param y Row.
   *  @return True for a row between a pole and the first lattice ring whose
   *          longitude filter spans more than one column. */
  static __attribute__((always_inline)) bool
  filtered_infill(const CoarseGrid &grid, int y) {
    const int downsample = grid.downsample;
    const bool infill_band =
        (y < downsample && (!SphereField::HAS_NORTH_POLE || y > 0)) ||
        (y >= H - downsample && (!SphereField::HAS_SOUTH_POLE || y < H - 1));
    return infill_band && grid.field.longitude_filter_width(y) > 1;
  }

  /** @brief Whether every column of a row is its pole. */
  static __attribute__((always_inline)) bool is_pole_row(int y) {
    return (SphereField::HAS_NORTH_POLE && y == 0) ||
           (SphereField::HAS_SOUTH_POLE && y == H - 1);
  }

  /** @brief Three quarters of one pixel and a quarter of another, rounded. */
  static ::Pixel blend_quarter(const ::Pixel &own, const ::Pixel &other) {
    return ::Pixel(static_cast<uint16_t>((3u * own.r + other.r + 2u) >> 2),
                   static_cast<uint16_t>((3u * own.g + other.g + 2u) >> 2),
                   static_cast<uint16_t>((3u * own.b + other.b + 2u) >> 2));
  }

  /** @brief Rounded mean of two pixels. */
  static ::Pixel pair_average(const ::Pixel &a, const ::Pixel &b) {
    return ::Pixel(static_cast<uint16_t>((a.r + b.r + 1u) >> 1),
                   static_cast<uint16_t>((a.g + b.g + 1u) >> 1),
                   static_cast<uint16_t>((a.b + b.b + 1u) >> 1));
  }

  /**
   * @brief Expands a half-resolution run's pair samples back to every column.
   * @param output Row whose even columns in [begin, end) hold their pair's
   *        box average, sampled at the pair's midpoint.
   * @param begin First column of the run.
   * @param end One past the run's last column.
   * @param before The pair sample beside @p begin, or null.
   * @param after The pair sample beside @p end, or null.
   * @details Each column takes three quarters of its own pair's sample and a
   * quarter of the neighbouring pair's on its side. Given neighbours stand in
   * at the run's ends; otherwise a full row wraps across the seam and a
   * clipped run repeats its end pairs.
   */
  HS_HOT_FLASH_MEMBER static void
  reconstruct_half_res_run(::Pixel *output, int begin, int end,
                           const ::Pixel *before = nullptr,
                           const ::Pixel *after = nullptr) {
    const bool full_row = begin == 0 && end == W;
    auto blend = [](const ::Pixel &own, const ::Pixel &other) {
      return blend_quarter(own, other);
    };
    // The last pair's right neighbour wraps onto the first sample, which the
    // loop overwrites first.
    const ::Pixel first = output[begin];
    ::Pixel previous = before ? *before : full_row ? output[W - 2] : first;
    for (int x = begin; x < end; x += 2) {
      const ::Pixel own = output[x];
      const ::Pixel next = x + 2 < end ? output[x + 2]
                           : after     ? *after
                           : full_row  ? first
                                       : own;
      output[x] = blend(own, previous);
      if (x + 1 < end)
        output[x + 1] = blend(own, next);
      previous = own;
    }
  }

  /**
   * @brief Expands a half-resolution row's pair samples back to every column.
   * @param ctx Frame context.
   * @param output Row of pair samples, as for reconstruct_half_res_run().
   * @param y Row index.
   * @details Under column time a full row expands arm by arm short of each
   * arm's last two columns, which are their own lanes: the pair before an
   * arm's first column is the other arm's last two columns of the previous
   * frame, the pair after its last pair its own last two columns. Other rows
   * expand run by run.
   */
  HS_HOT_FLASH_MEMBER static void expand_half_res_row(const FlushContext &ctx,
                                                      ::Pixel *output, int y) {
    const ColumnRuns &runs = ctx.runs;
    if (!ctx.column_time || runs.count != 2 || runs.items[0].begin != 0 ||
        runs.items[1].end != W) {
      for (int r = 0; r < runs.count; ++r)
        reconstruct_half_res_run(output, runs.items[r].begin,
                                 runs.items[r].end);
      return;
    }
    const ::Pixel *previous_row = ctx.previous + y * W;
    for (int arm = 0; arm < W; arm += HALF) {
      const int other = arm == 0 ? HALF : 0;
      const ::Pixel before = pair_average(previous_row[other + HALF - 2],
                                          previous_row[other + HALF - 1]);
      const ::Pixel after =
          pair_average(output[arm + HALF - 2], output[arm + HALF - 1]);
      reconstruct_half_res_run(output, arm, arm + HALF - 2, &before, &after);
    }
  }

  /**
   * @brief Lifts @p v onto the wrap branch nearest @p reference.
   * @param v Column offset to unwrap.
   * @param reference Offset whose branch the result must share.
   * @param period Full turn in @p v's units.
   * @return @p v shifted by at most one period, within half a period of
   *         @p reference.
   */
  static __attribute__((always_inline)) float
  unwrap_near(float v, float reference, float period) {
    const float half = period * 0.5f;
    const float delta = v - reference;
    return delta > half ? v - period : (delta < -half ? v + period : v);
  }

  /**
   * @brief unwrap_near() for stored column offsets, in WARP_SCALE units.
   * @param v Stored column offset to unwrap.
   * @param reference Stored offset whose branch the result must share.
   * @return @p v shifted by at most one period, within half a period of
   *         @p reference; identical to the float form on these integers.
   */
  static __attribute__((always_inline)) int unwrap_near(int v, int reference) {
    constexpr int PERIOD = static_cast<int>(WRAP_PERIOD);
    constexpr int HALF_PERIOD = PERIOD / 2;
    const int delta = v - reference;
    return delta > HALF_PERIOD ? v - PERIOD
                               : (delta < -HALF_PERIOD ? v + PERIOD : v);
  }
  static_assert(static_cast<float>(static_cast<int>(WRAP_PERIOD)) ==
                        WRAP_PERIOD &&
                    static_cast<int>(WRAP_PERIOD) % 2 == 0,
                "Feedback<W,H>: the integer unwrap needs an even, exact "
                "period");

  /**
   * @brief Projected field coordinates of one lattice sample.
   * @param field Spherical layout the lattice belongs to.
   * @param position Lattice position as populate() hands it out.
   * @param point Exact field coordinates of the same sample.
   * @return @p point on a pole row, where the projection is degenerate,
   *         otherwise the projection of @p position, so the offset subtraction
   *         cancels the projection's approximation.
   */
  static typename SphereField::Coordinates
  lattice_origin(const SphereField &field, const math::Vector &position,
                 const typename SphereField::Coordinates &point) {
    const bool pole_row = (SphereField::HAS_NORTH_POLE && point.y == 0.0f) ||
                          (SphereField::HAS_SOUTH_POLE && point.y == H - 1);
    return pole_row ? point : field.project(position);
  }
  // populate_warp_field() canonicalizes the stored column offset onto
  // [-W*WARP_SCALE/2, W*WARP_SCALE/2] and clamps the elapsed-scaled offset
  // to that range before the int16_t cast.
  static_assert(W * WARP_SCALE * 0.5f <= 32767.0f,
                "Feedback<W,H>: canonical warp offset must fit int16_t");
  // Under runtime geometry, row offsets saturate at +/-32767 WARP_SCALE units.
#if !HS_RUNTIME_DISPLAY_GEOMETRY
  static_assert(
      math::PI_F * math::ROWS_PER_RADIAN<H> * WARP_SCALE <= 32767.0f,
      "Feedback<W,H>: fixed-geometry warp row offset must fit int16_t");
#endif

  /**
   * @brief Tests whether the previous frame has any non-black pixel.
   * @param cv Canvas whose previous-frame buffer is scanned.
   * @return True on the first lit pixel found, false if the frame is all black.
   * @details Scans the clip band, which covers the canvas for this full-frame stage.
   */
  static bool any_pixel_lit(const Canvas &cv) {
    const auto &clip = cv.clip();
    const ColumnRuns runs = make_column_runs(clip.x_clip());
    const ::Pixel *previous = cv.prev_data();
    for (int y = clip.render_y_start(); y < clip.render_y_end(); ++y) {
      const ::Pixel *row = previous + y * W;
      for (int i = 0; i < runs.count; ++i) {
        for (int x = runs.items[i].begin; x < runs.items[i].end; ++x) {
          const ::Pixel pixel = row[x];
          if (pixel.r | pixel.g | pixel.b)
            return true;
        }
      }
    }
    return false;
  }

  /**
   * @brief Chooses one color for a longitude-aliased pole row.
   * @param pole_row Base of the pole row in the previous frame.
   */
  HS_O3_FN static ::Pixel select_pole_sample(const ::Pixel *pole_row) {
    ::Pixel selected = pole_row[0];
    uint32_t selected_energy =
        static_cast<uint32_t>(selected.r) + selected.g + selected.b;
    for (int x = 1; x < W; ++x) {
      const ::Pixel candidate = pole_row[x];
      const uint32_t energy =
          static_cast<uint32_t>(candidate.r) + candidate.g + candidate.b;
      if (energy > selected_energy) {
        selected = candidate;
        selected_energy = energy;
      }
    }
    return selected;
  }

  /**
   * @brief Chooses each pole row's color in a frame.
   * @param frame Frame buffer.
   * @param poles Receives the shared value of each pole row present.
   */
  static void select_poles(const ::Pixel *frame,
                           ::Pixel (&poles)[SphereField::POLE_STORAGE_COUNT]) {
    if (SphereField::HAS_NORTH_POLE)
      poles[0] = select_pole_sample(frame);
    if (SphereField::HAS_SOUTH_POLE)
      poles[SphereField::HAS_NORTH_POLE ? 1 : 0] =
          select_pole_sample(frame + (H - 1) * W);
  }

  /**
   * @brief Samples a read across one of its target arm's cuts.
   * @param ctx Frame context.
   * @param arm First column of the target's arm.
   * @param x Source column, unwrapped near the target.
   * @param y Source row.
   * @return The bilinear sample, the taps on each side of the cut from the
   *         buffer holding that side: the stale frame behind the arm, the
   *         frame being drawn past it, the previous frame within it.
   */
  HS_HOT_FLASH_MEMBER Sample sample_across(const FlushContext &ctx, int arm,
                                           float x, float y) const {
    const SphereField &field = ctx.grid.field;
    const bool behind = x < arm;
    const ::Pixel *low = behind ? ctx.back : ctx.previous;
    const ::Pixel *high = behind ? ctx.previous : ctx.current;
    const auto &low_poles = behind ? ctx.back_poles : ctx.poles;
    const auto &high_poles = behind ? ctx.poles : ctx.current_poles;
    const float first = static_cast<float>(behind ? arm : arm + HALF);
    if (x >= first)
      return sample_bilinear_prev(field, high, high_poles, x, y);
    if (x <= first - 1.0f)
      return sample_bilinear_prev(field, low, low_poles, x, y);
    const Sample a =
        sample_bilinear_prev(field, low, low_poles, first - 1.0f, y);
    const Sample b = sample_bilinear_prev(field, high, high_poles, first, y);
    const float t = x - (first - 1.0f);
    return {a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t,
            a.b + (b.b - a.b) * t};
  }

  /** @brief Whether a near-black sample skips its color transform.
   * @param sample Unquantized linear sample.
   * @param enabled Whether near-black samples skip the color transform.
   * @return True when enabled and all channels are below NEAR_BLACK. */
  static HS_O3_FN HS_HOT_INLINE bool skip_near_black(const Sample &sample,
                                                     bool enabled) {
    return enabled && sample.r < NEAR_BLACK && sample.g < NEAR_BLACK &&
           sample.b < NEAR_BLACK;
  }

  /**
   * @brief Bilinearly samples the Canvas front buffer (previous frame).
   * @param field Spherical topology and interpolation policy.
   * @param prev Base of the previous-frame buffer, row-major with stride W.
   * @param poles Shared values for every aliased column of the pole rows.
   * @param bx Fractional column in [-W, 2W).
   * @param by Fractional row; north crossings reflect with a half-turn.
   * @return The interpolated channels, returned in registers.
   */
  HS_O3_FN __attribute__((noinline)) Sample
  sample_bilinear_prev(const SphereField &field, const ::Pixel *prev,
                       const ::Pixel (&poles)[SphereField::POLE_STORAGE_COUNT],
                       float bx, float by) const {
    Sample out;
    field.sample_bilinear_rgb(prev, poles, bx, by, out.r, out.g, out.b);
    return out;
  }

  /**
   * @brief Calls @p transform out of line.
   * @param transform Single-pixel colour transform.
   * @param r Sampled red.
   * @param g Sampled green.
   * @param b Sampled blue.
   * @param level Fade-table level.
   * @return The transformed pixel.
   */
  template <typename TransformPixelT>
  __attribute__((noinline)) static ::Pixel
  transform_pixel_call(TransformPixelT &transform, float r, float g, float b,
                       int level) {
    return transform(r, g, b, level);
  }

  /**
   * @brief Transforms a sample, out of line when the composite pairs pixels.
   * @tparam PAIR_PIXELS Whether the composite pairs pixels.
   * @param transform Single-pixel colour transform.
   * @param s Sample.
   * @param level Fade-table level.
   * @return The transformed pixel.
   */
  template <bool PAIR_PIXELS, typename TransformPixelT>
  static HS_HOT_INLINE ::Pixel transform_sample(TransformPixelT &transform,
                                                const Sample &s, int level) {
    if constexpr (PAIR_PIXELS)
      return transform_pixel_call(transform, s.r, s.g, s.b, level);
    else
      return transform(s.r, s.g, s.b, level);
  }

  /** @brief Calls @p transform out of line. */
  template <typename TransformPairT>
  __attribute__((noinline)) static void
  transform_pair_call(TransformPairT &transform, float r0, float g0, float b0,
                      float r1, float g1, float b1, int level, ::Pixel &p0,
                      ::Pixel &p1) {
    transform(r0, g0, b0, r1, g1, b1, level, p0, p1);
  }

  /** @brief Quantizes an unclamped [0, 65535]-scale channel to a ::Pixel
   *  component, round-to-nearest; NaN maps to the hi bound. */
  static HS_HOT_INLINE uint16_t quantize16(float v) {
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 65535.0f) + 0.5f);
  }

  /**
   * @brief The bound generator's own configuration, as a cache-key scalar.
   * @details Only the seed is keyed. Fractal and rotation settings on
   * NoiseParams::noise must remain at defaults, or the persistent arena must be
   * reset before re-running init_storage() after changes.
   */
  static uint32_t noise_config_key(const Animation::NoiseParams *noise) {
    return noise ? static_cast<uint32_t>(noise->seed) : 0u;
  }

  /**
   * @brief Blends a tap near an arm's end into a covered pixel, or passes it
   *        on.
   * @param x Column coordinate as plotted.
   * @param column @p x rounded and wrapped as the canvas writes it.
   * @param y Row coordinate.
   * @param color,age,alpha Shading registers.
   * @param pass Downstream 2D callback, the canvas write.
   * @details The color blends as the canvas write blends it, and the
   * transmittance takes the same blend toward zero, so finish() can add the
   * rest of the pixel's composite beneath both.
   */
  template <typename PassFnT>
  HS_HOT_FLASH_MEMBER void plot_near_arm_end(float x, int column, float y,
                                             const ::Pixel &color, float age,
                                             float alpha, PassFnT &pass) {
    const int row = Filter::round_row<H>(y);
    const int arm = ColumnTime::arm_start(column);
    const int slot =
        row < 0 || row >= H ? -1 : cover_slot(held->rows[row], arm, column);
    if (slot < 0) {
      pass(x, y, color, age, alpha);
      return;
    }
    ::Pixel &dst = held->current[row * W + column];
    uint16_t &cover = held->cover[slot];
    const uint16_t ahead = cover & AHEAD_BIT;
    if (alpha >= 1.0f) {
      dst = color;
      cover = ahead;
      return;
    }
    const uint16_t weight = frac_to_q16(alpha);
    dst = dst.lerp16(color, weight);
    cover =
        (lerp_q16(cover & TRANSMITTANCE_BITS, 0, weight) & TRANSMITTANCE_BITS) |
        ahead;
  }

  /** @brief Rows and lattice rings under POLAR_TARGET_SINE, at most, over
   *  the display geometries math::set_display_geometry() admits: its ends
   *  sampled, a row's angle taken a milliradian toward the cap and one row
   *  added at each end; describe_rows() and init_storage() check the bounds. */
  struct CapExtent {
    int rows;
    int rings;
  };
  static constexpr CapExtent CAP_EXTENT = [] {
#if HS_RUNTIME_DISPLAY_GEOMETRY
    constexpr float CAP = math::PI_F / 6.0f + 1.0e-3f;
    constexpr int STEPS = 24;
    CapExtent most{0, 0};
    for (int i = 0; i <= STEPS; ++i)
      for (int j = 0; j <= STEPS; ++j) {
        const float north = 0.25f * math::PI_F * i / STEPS;
        const float south = math::PI_F * (0.75f + 0.25f * j / STEPS);
        const float pitch = (south - north) / (H - 1);
        int top = 0, bottom = 0;
        for (int y = 0; y < H; ++y) {
          const float phi = north + y * pitch;
          top += phi < CAP;
          bottom += phi > math::PI_F - CAP;
        }
        top = std::min(top + 1, H);
        bottom = std::min(bottom + 1, H);
        const int rings = CACHE_FIELD.rings_within(top, bottom);
        most.rows = std::max(most.rows, std::min(top + bottom, H));
        most.rings = std::max(most.rings, rings);
      }
    return most;
#else
    int rows = 0;
    for (int y = 0; y < H; ++y)
      rows += SphereField::latitude_sine(y) < POLAR_TARGET_SINE;
    return CapExtent{rows,
                     polar_rings(CACHE_FIELD, CACHE_FIELD.ring_count()).rows()};
#endif
  }();

  /** @brief Coverage slots per arm of a row: every column of a row under
   *  POLAR_TARGET_SINE, where a cap-plane lane may turn any distance;
   *  otherwise the columns a quarter-turn forward offset, the elapsed
   *  scale's cap, carries past the arm's last column, and the pair column
   *  before them. */
  static constexpr int COVER_CARTESIAN = HALF / 2 + 1;
  /** @brief A row's coverage slots per arm. */
  static HS_HOT_INLINE int cover_columns(const RowGeometry &g) {
    return g.wide ? HALF : COVER_CARTESIAN;
  }
  /** @brief Coverage slots of a frame. */
  static constexpr int COVER_SLOTS =
      2 * (CAP_EXTENT.rows * HALF + (H - CAP_EXTENT.rows) * COVER_CARTESIAN);
  static_assert(COVER_SLOTS <= UINT16_MAX,
                "Feedback<W,H>: RowGeometry::cover_base holds a slot");
  /**
   * @brief Coverage slot of a column.
   * @param g The column's row geometry.
   * @param arm First column of the column's arm.
   * @param x Column.
   * @return Index into FlushContext::cover, or -1 for a column too far from
   *         its arm's end to have one.
   */
  static HS_HOT_INLINE int cover_slot(const RowGeometry &g, int arm, int x) {
    const int width = cover_columns(g);
    const int slot = arm + HALF - 1 - x;
    return slot < width ? g.cover_base + (arm != 0 ? width : 0) + slot : -1;
  }

  /**
   * @brief Bound on an infill row's longitude filter width.
   * @param y Row.
   * @return An odd width no smaller than the layout's.
   */
  static constexpr int filter_window_bound(int y) {
    constexpr int WIDEST = (W & 1) ? W : W - 1;
#if HS_RUNTIME_DISPLAY_GEOMETRY
    // math::set_display_geometry() puts the far end at least three quarters
    // of a turn from each pole, so a row d rows from an end lies at least
    // d * (3 pi / 4) / (H - 1) from its pole; a - a^3 / 6 bounds that angle's
    // sine from below.
    const int rows = std::min(y, H - 1 - y);
    const float angle = rows * (0.75f * math::PI_F / (H - 1));
    const float sine = angle - angle * angle * angle * (1.0f / 6.0f);
#else
    const float sine = SphereField::latitude_sine(y);
#endif
    if (sine < 2.0f / W)
      return WIDEST;
    int width = 1;
    while (width * sine < 1.0f)
      width += 2;
    return std::min(width + 2, WIDEST);
  }

  /** @brief Finish samples a row keeps per arm: the lanes its reconstruction
   *  reaches, at most the arm the infill filter narrows to. */
  static constexpr int finish_window(const RowGeometry &g) {
    return g.infill ? std::min(2 * g.reach + 1, HALF) : g.pairs ? 2 : 0;
  }

  /** @brief Finish window samples a frame may keep: a pair sample per arm on
   *  every row, and the infill rows' filter widths. */
  static constexpr int WINDOW_SLOTS = [] {
    int slots = 0;
    for (int y = 0; y < H; ++y) {
#if HS_RUNTIME_DISPLAY_GEOMETRY
      const bool infill = y < CACHE_DOWNSAMPLE || y >= H - CACHE_DOWNSAMPLE;
#else
      const bool infill =
          (y < CACHE_DOWNSAMPLE && (!SphereField::HAS_NORTH_POLE || y > 0)) ||
          (y >= H - CACHE_DOWNSAMPLE &&
           (!SphereField::HAS_SOUTH_POLE || y < H - 1));
#endif
      slots += 2 * (infill ? std::min(filter_window_bound(y), HALF) : 1);
    }
    return slots;
  }();

  /** @brief Traps while a split frame awaits finish(). */
  HS_FLASH_MEMBER void require_unheld() const {
    HS_CHECK(!held, "feedback: a split frame was never finished");
  }

  /** @brief Releases the held frame and its scratch. */
  HS_FLASH_MEMBER void end_split() {
    held = nullptr;
    held_floor = HALF;
    scratch_arena_a.set_offset(held_mark);
  }

  ::Feedback::Style *feedback_style; /**< Bound feedback Style (non-owning). */
  bool enabled = true;               /**< When false, flush() is skipped. */
  /** Split frame flush_split() holds for finish(), or null. */
  FlushContext *held = nullptr;
  /** scratch_arena_a offset finish() rewinds to. */
  size_t held_mark = 0;
  /** Least held column past an arm's first in the held frame; HALF holds no
   *  column. */
  int held_floor = HALF;
  /** The same per row. */
  uint8_t row_floor[H] = {};
  static_assert(HALF <= UINT8_MAX, "Feedback<W,H>: row_floor holds HALF");
  WarpCache warp_cache; /**< Persistent warp field and lattice origins. */
  /** @brief Persistent measure_reach() spans of the cached warp field. */
  struct Spans {
    int16_t *reach = nullptr;   ///< Per row, the reach behind each arm.
    float pole_half_res = 0.0f; ///< Style::pole_half_res they were measured at.
    bool valid = false;         ///< Whether they belong to the cached field.
  } spans;
  FadeTable fade_table; /**< Persistent per-elapsed-time transforms. */

public:
  /** @brief Cap-cell budget: the most polar rings a display geometry
   *  places (CAP_EXTENT). */
  static constexpr int CACHE_CAP_CELLS = CAP_EXTENT.rings * CACHE_COLUMNS;

  /** @brief Lattice-origin budget: the most samples a display geometry
   *  places. Under runtime geometry, row y lies between y * (3 pi / 4) /
   *  (H - 1) and pi / 4 more from the north pole (math::set_display_geometry()),
   *  and each ring takes the largest sine over that span. */
#if HS_RUNTIME_DISPLAY_GEOMETRY
  static constexpr int CACHE_ORIGIN_SAMPLES =
      CACHE_FIELD.sample_count([](int y) {
        constexpr float PITCH = 0.75f * math::PI_F / (H - 1);
        const float low = y * PITCH, high = 0.25f * math::PI_F + y * PITCH;
        if (low <= 0.5f * math::PI_F && high >= 0.5f * math::PI_F)
          return 1.0f;
        // Taylor sine through x^9 overestimates on [0, pi / 2].
        auto sine = [](float x) {
          x = x > 0.5f * math::PI_F ? math::PI_F - x : x;
          const float x2 = x * x;
          return x * (1.0f - x2 / 6.0f *
                                 (1.0f - x2 / 20.0f *
                                             (1.0f - x2 / 42.0f *
                                                         (1.0f - x2 / 72.0f))));
        };
        return std::min(1.0f, std::max(sine(low), sine(high)) + 1.0e-4f);
      });
#else
  static constexpr int CACHE_ORIGIN_SAMPLES = CACHE_FIELD.sample_count();
#endif

  /** @brief Persistent bytes of one warp field set over CACHE_CELLS: int16
   *  column and row offsets and polar cells' int16 cap-plane offsets. */
  static constexpr size_t FIELD_SET_BYTES =
      CACHE_CELLS * 2 * sizeof(int16_t) + CACHE_CAP_CELLS * 2 * sizeof(int16_t);
  /** @brief Persistent bytes init_storage() reserves over CACHE_CELLS: two
   *  alternating field sets, the lattice's projected origins, the fade table
   *  and the rows' reaches. */
  static constexpr size_t STORAGE_BYTES =
      2 * FIELD_SET_BYTES +
      CACHE_ORIGIN_SAMPLES * sizeof(typename SphereField::Coordinates) +
      FadeTable::BYTES + H * sizeof(int16_t);
  /** @brief Arena allocations init_storage() makes; each pads its share of
   *  STORAGE_BYTES by less than alignof(std::max_align_t). */
  static constexpr size_t STORAGE_ALLOCATIONS = 10;

  /**
   * @brief scratch_arena_a bytes a split frame holds from flush_split() to
   *        finish() with its warp field and fade table in storage.
   * @return The frame context, filter row, row geometry, per-row reaches,
   *         coverage and held-lane spans, each allocation padded.
   */
  static constexpr size_t split_scratch_bytes() {
    constexpr size_t PAD = alignof(std::max_align_t);
    return sizeof(FlushContext) + W * sizeof(::Pixel) +
           H * (sizeof(RowGeometry) + sizeof(int16_t) + sizeof(AheadSpan)) +
           COVER_SLOTS * sizeof(uint16_t) + 6 * PAD;
  }

  /**
   * @brief scratch_arena_a bytes a composite or finish() takes past the
   *        frame it holds, and releases before returning.
   * @return finish()'s sample windows, visited rows and sweep states, each
   *         allocation padded; sweep_deep()'s states are fewer.
   */
  static constexpr size_t transient_scratch_bytes() {
    constexpr size_t PAD = alignof(std::max_align_t);
    return WINDOW_SLOTS * sizeof(::Pixel) +
           H * (sizeof(int16_t) + sizeof(SweepRow)) + 3 * PAD;
  }
};

} // namespace Pixel

} // namespace Filter
