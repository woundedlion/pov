/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <algorithm>
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
 * where equirect longitude offsets blow up as 1/sin(phi).
 * With temporal sampling on, each pixel reads its source at the time elapsed
 * between their lit columns (FeedbackColumnTime): lattice offsets carry the
 * elapsed scale of their displacement, cells whose sources all lie in one
 * buffer fade by the cell's elapsed time, and the rest resolve per pixel.
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
  using ColumnSource = typename ColumnTime::Source;
  using FadeTable = FeedbackFadeTable;

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
   *        split frame's forward-band coverage when it lands in a band.
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
    if (column - ColumnTime::arm_start(column) >= held_floor) [[unlikely]] {
      plot_near_arm_end(x, column, y, color, age, alpha, pass);
      return;
    }
    pass(x, y, color, age, alpha);
  }

  /**
   * @brief Enables or disables feedback.
   * @param value When false, flush() is skipped entirely.
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
    const int COLUMNS = W / ds;
    const SphereField FIELD = field_for(ds);
    const int RINGS = FIELD.ring_count();
    const PolarRings POLAR = polar_rings(FIELD, RINGS);
    return (2 * RINGS * COLUMNS + 2 * FIELD.sample_count()) * sizeof(int16_t) +
           (POLAR.rows() * COLUMNS + POLAR.samples) * sizeof(CapOffset) +
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
    warp_cache.template init_storage<CACHE_CELLS>(
        arena, CACHE_FIELD,
        polar_rings(CACHE_FIELD, CACHE_FIELD.ring_count()).rows() *
            CACHE_COLUMNS,
        [](const math::Vector &position,
           const typename SphereField::Coordinates &point) {
          return lattice_origin(CACHE_FIELD, position, point);
        });
    fade_table.init_storage(arena);
    spans.reach = arena.allocate_n<int16_t>(H);
    spans.ahead = arena.allocate_n<Ahead>(H);
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
   * segment's cylindrical clip. No-op when disabled.
   */
  HS_O3_BEGIN
  void flush(Canvas &cv, float alpha) {
    require_unheld();
    if (!enabled) {
      if (cv.keeps_stale_pixels())
        cv.clear_buffer();
      return;
    }

    ScratchScope scope(scratch_arena_a);
    FlushContext ctx = prepare_flush(cv, alpha, scope.get_arena());
    if (cv.keeps_stale_pixels() && !ctx.back)
      cv.clear_buffer();
    if (ctx.mode == FlushMode::SKIP)
      return;

    HS_PROFILE(feedback_composite);
    if (ctx.back) {
      ctx.pass = Pass::ARM_STARTS;
      composite_pass(ctx);
      ctx.pass = Pass::REST;
    }
    composite_pass(ctx);
  }
  HS_O3_END

  /**
   * @brief Composites a split frame's columns outside its forward bands and
   *        holds the frame for finish().
   * @param cv Target canvas.
   * @param alpha Global blend alpha in [0, 1].
   * @details A frame that cannot split, or a row whose forward reach exceeds a
   * quarter turn, composites as flush() does. The held frame keeps its
   * scratch_arena_a allocations until finish(); scopes opened meanwhile must
   * close first.
   */
  HS_FLASH_MEMBER void flush_split(Canvas &cv, float alpha) {
    require_unheld();
    if (!enabled) {
      if (cv.keeps_stale_pixels())
        cv.clear_buffer();
      return;
    }
    const size_t mark = scratch_arena_a.get_offset();
    FlushContext &ctx = *scratch_arena_a.make<FlushContext>(
        prepare_flush(cv, alpha, scratch_arena_a, true));
    if (cv.keeps_stale_pixels() && !ctx.back)
      cv.clear_buffer();
    if (ctx.mode != FlushMode::SKIP) {
      HS_PROFILE(feedback_composite);
      if (ctx.back) {
        ctx.pass = Pass::ARM_STARTS;
        composite_pass(ctx);
        ctx.pass = Pass::REST;
      }
      composite_pass(ctx);
    }
    int floor = ColumnTime::HALF;
    if (ctx.ahead)
      for (int y = 0; y < H; ++y)
        floor = std::min(floor, static_cast<int>(ctx.ahead[y].held_from));
    if (floor == ColumnTime::HALF) {
      scratch_arena_a.set_offset(mark);
      return;
    }
    ctx.cover = scratch_arena_a.allocate_n<uint16_t>(2 * H * COVER_STRIDE);
    ctx.kept = scratch_arena_a.allocate_n<::Pixel>(2 * H);
    open_coverage(ctx);
    held = &ctx;
    held_mark = mark;
    held_floor = floor;
  }

  /**
   * @brief Composites a held frame's forward bands from the frame as plotted,
   *        under the coverage its taps left there.
   * @details A source past a band's forward cut reads the frame being drawn,
   * whose columns there the plots have finished. No-op without a held frame.
   */
  HS_FLASH_MEMBER void finish(Canvas &) {
    if (!held)
      return;
    FlushContext &ctx = *held;
    if (SphereField::HAS_NORTH_POLE)
      ctx.current_poles[0] = select_pole_sample(ctx.current);
    if (SphereField::HAS_SOUTH_POLE)
      ctx.current_poles[SphereField::HAS_NORTH_POLE ? 1 : 0] =
          select_pole_sample(ctx.current + (H - 1) * W);
    {
      HS_PROFILE(feedback_forward);
      ctx.pass = Pass::FORWARD;
      composite_pass(ctx);
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
    ColumnRun items[3];
    int count;
  };

  enum class FlushMode : uint8_t { SKIP, HUE, PLAIN, GENERAL };

  /** @brief Columns one composite pass writes. */
  enum class Pass : uint8_t {
    ALL,        ///< Every column of the band.
    ARM_STARTS, ///< Each arm's first reach columns, before any stale
                ///< column behind a cut is overwritten.
    REST,       ///< The columns ARM_STARTS and FORWARD leave.
    FORWARD,    ///< Each arm's last Ahead::band columns, after the plots.
  };

  /** @brief A row's forward band in a split frame. */
  struct Ahead {
    /// Columns before each arm's end whose sources may lie past its forward
    /// cut, even; zero where the row composites unsplit.
    int16_t band;
    /// Column past each arm's first from which taps are held for finish():
    /// the band, and on a half-resolution row the column before it, whose
    /// expansion reads the band's first pair; HALF without a band.
    int16_t held_from;
  };

  /** @brief Per-frame composite inputs prepare_flush() resolves. */
  struct FlushContext {
    CoarseGrid grid;
    RenderBand band{};
    WarpField warp{};
    ::Pixel *filtered_row = nullptr;
    const ::Pixel *previous = nullptr;
    ::Pixel *current = nullptr;
    ::Pixel poles[SphereField::POLE_STORAGE_COUNT]{};
    /// The drawing buffer's stale frame, two frames back, which holds the
    /// content behind each arm's backward cut; null when it was cleared.
    const ::Pixel *back = nullptr;
    /// Shared pole values of the stale frame.
    ::Pixel back_poles[SphereField::POLE_STORAGE_COUNT]{};
    /// Per row, the columns from each arm's first whose sources may lie
    /// behind its backward cut, even; negative where the reach exceeds a
    /// quarter turn and the row slides instead. Null without @c back.
    const int16_t *reach = nullptr;
    /// Per row, the forward band of a split frame; null when unsplit.
    const Ahead *ahead = nullptr;
    /// Per row and arm, the Q16 transmittance of the taps plotted over each
    /// held column (cover_index()); their colors accumulate in the drawing
    /// buffer. Null when unsplit.
    uint16_t *cover = nullptr;
    /// Per row and arm, the pair sample REST keeps for the held column before
    /// a half-resolution band.
    ::Pixel *kept = nullptr;
    /// Shared pole values of the frame being drawn, set by finish().
    ::Pixel current_poles[SphereField::POLE_STORAGE_COUNT]{};
    /// Columns the current pass writes.
    Pass pass = Pass::ALL;
    ColumnRuns runs{};
    typename SphereField::Ring control_ring1{};
    int field_y0 = 0;
    int field_y1 = 0;
    int control_y0 = 0;
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
      if (split_arms && begin < ColumnTime::HALF && end > ColumnTime::HALF) {
        runs.items[runs.count++] = {begin, ColumnTime::HALF};
        begin = ColumnTime::HALF;
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
   */
  __attribute__((always_inline)) WarpField
  acquire_warp_field(Arena &scratch, const CoarseGrid &grid,
                     const RenderBand &band, bool &reused) {
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
        });
  }

  /**
   * @brief Fills @p warp's cells for the band from the Style's warp.
   * @param grid Coarse layout of the frame.
   * @param band Rows and columns the frame renders.
   * @param warp Field to fill.
   * @param controls Scratch for every lattice sample's offset.
   * @param sample_caps Scratch for every polar lattice sample's cap offset.
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
      auto ring = grid.field.ring(band.field_y_begin);
      for (int field_y = band.field_y_begin; field_y <= band.field_y_end;
           ++field_y, ring = grid.field.next_ring(ring)) {
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
          warp.x_offsets[index] = static_cast<int16_t>(hs::clamp(
              offset_x * scale, -WRAP_PERIOD * 0.5f, WRAP_PERIOD * 0.5f));
          warp.y_offsets[index] = static_cast<int16_t>(
              hs::clamp(offset_y * scale, -32767.0f, 32767.0f));
          if (field_y < grid.polar.north_rings ||
              field_y >= grid.polar.south_ring) {
            const CapOffset da = sample_caps[grid.polar.slot(longitude.left)];
            const CapOffset db = sample_caps[grid.polar.slot(longitude.right)];
            CapOffset cap = CapPlane::lerp_offset(da, db, longitude.mix);
            if (column_time)
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
   * @param split Whether the frame's forward bands composite after its plots.
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
    ctx.warp = acquire_warp_field(scratch, grid, band, reused);
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
    if (ctx.column_time)
      ctx.fades = fade_table.acquire(*feedback_style, scratch);
    const bool hue_fade = feedback_style->color_fn == &::Feedback::hue_fade;
    ctx.black_skips_color = hue_fade;
    ctx.pole_half_res = feedback_style->pole_half_res;
    if (grid.polar.rows() > 0 && !math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    ctx.previous = cv.prev_data();
    ctx.current = cv.data();
    if (SphereField::HAS_NORTH_POLE)
      ctx.poles[0] = select_pole_sample(ctx.previous);
    if (SphereField::HAS_SOUTH_POLE)
      ctx.poles[SphereField::HAS_NORTH_POLE ? 1 : 0] =
          select_pole_sample(ctx.previous + (H - 1) * W);
    ctx.runs = make_column_runs(band.x_clip, ctx.column_time);
    const bool clipped =
        band.x_clip.active || band.y_begin > 0 || band.y_end < H;
    if (ctx.column_time && cv.keeps_stale_pixels() && alpha >= 1.0f &&
        !clipped) {
      ctx.back = ctx.current;
      if (SphereField::HAS_NORTH_POLE)
        ctx.back_poles[0] = select_pole_sample(ctx.back);
      if (SphereField::HAS_SOUTH_POLE)
        ctx.back_poles[SphereField::HAS_NORTH_POLE ? 1 : 0] =
            select_pole_sample(ctx.back + (H - 1) * W);
      // An unchanged cached field keeps its spans.
      const bool stored = reused && spans.reach;
      if (!(stored && spans.valid &&
            spans.pole_half_res == ctx.pole_half_res)) {
        int16_t *reach = stored ? spans.reach : scratch.allocate_n<int16_t>(H);
        Ahead *ahead = stored ? spans.ahead : scratch.allocate_n<Ahead>(H);
        measure_reach(grid, ctx.warp, half_res_limit(ctx.pole_half_res), reach,
                      ahead);
        spans.valid = stored;
        spans.pole_half_res = ctx.pole_half_res;
        ctx.reach = reach;
        ctx.ahead = ahead;
      } else {
        ctx.reach = spans.reach;
        ctx.ahead = spans.ahead;
      }
      if (!split)
        ctx.ahead = nullptr;
    }
    ctx.field_y0 = band.field_y_begin;
    ctx.field_y1 = ctx.field_y0 + (ctx.field_y0 < band.field_y_end ? 1 : 0);
    const auto control_ring0 = grid.field.ring(ctx.field_y0);
    ctx.control_ring1 = control_ring0;
    if (ctx.field_y1 > ctx.field_y0)
      ctx.control_ring1 = grid.field.next_ring(control_ring0);
    ctx.control_y0 = control_ring0.y;

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
   * @brief Composites the warped previous frame through a color transform.
   * @tparam PAIR_PIXELS Whether to interleave two pixels per step through
   * @p transform_pair.
   * @param ctx Frame context from prepare_flush().
   * @param transform_pixel Maps one sampled RGB triple and its fade-table
   *        level to the output pixel.
   * @param transform_pair Maps two sampled RGB triples sharing one level to
   *        two output pixels.
   */
  template <bool PAIR_PIXELS, typename TransformPixelT, typename TransformPairT>
  void composite_pixels(const FlushContext &ctx,
                        TransformPixelT &&transform_pixel,
                        TransformPairT &&transform_pair) {
    const CoarseGrid &grid = ctx.grid;
    const RenderBand &band = ctx.band;
    const int downsample = grid.downsample;
    const int coarse_columns = grid.columns;
    const int row_begin = band.y_begin;
    const int row_end = band.y_end;
    const int16_t *x_offsets = ctx.warp.x_offsets;
    const int16_t *y_offsets = ctx.warp.y_offsets;
    const PolarRings polar = grid.polar;
    constexpr float INVERSE_WARP_SCALE = 1.0f / WARP_SCALE;
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
    const bool straddled_cells = ColumnTime::HALF % downsample != 0;
    const bool forward = ctx.pass == Pass::FORWARD;
    constexpr float LEVELS_PER_COLUMN =
        FadeTable::STEPS_PER_FRAME * ColumnTime::INVERSE_HALF;
    int field_y0 = ctx.field_y0;
    int field_y1 = ctx.field_y1;
    auto control_ring1 = ctx.control_ring1;
    int control_y0 = ctx.control_y0;
    int control_y1 = control_ring1.y;
    const float half_res_sine = half_res_limit(ctx.pole_half_res);
    for (int y = row_begin; y < row_end; ++y) {
      const int row = y * W;
      const bool filter_output =
          !band.x_clip.active && filtered_infill(grid, y);
      const bool defer_filter = filter_output && !opaque;
      // A forward band composites beneath its coverage.
      ::Pixel *output = defer_filter || forward ? filtered_row : current + row;
      // The pair expansion needs plain stores, and a longitude-filtered row
      // is reconstructed at its footprint already.
      const bool half_res = opaque && !filter_output &&
                            SphereField::latitude_sine(y) < half_res_sine;
      while (y > control_y1 && field_y1 < band.field_y_end) {
        field_y0 = field_y1;
        control_y0 = control_y1;
        ++field_y1;
        control_ring1 = grid.field.next_ring(control_ring1);
        control_y1 = control_ring1.y;
      }
      // The last ring lands on the band's last row; short of it the weights
      // below extrapolate off a stale control pair.
      HS_CHECK(y <= control_y1, "feedback warp row %d past last control row %d",
               y, control_y1);
      // Interpolating outside the populated band silently corrupts pixels.
      HS_CHECK(field_y0 >= band.field_y_begin && field_y1 <= band.field_y_end,
               "feedback warp ring %d outside populated band [%d,%d]", field_y1,
               band.field_y_begin, band.field_y_end);
      const int control_height = control_y1 - control_y0;
      const float fy = control_height > 0
                           ? static_cast<float>(y - control_y0) / control_height
                           : 0.0f;
      const float wy0 = 1.0f - fy, wy1 = fy;
      const int row0 = field_y0 * coarse_columns;
      const int row1 = field_y1 * coarse_columns;

      const ColumnRuns pass_columns = pass_runs(ctx, y);
      if (pass_columns.count == 0)
        continue;
      const bool behind = ctx.back && ctx.reach[y] >= 0;
      // A forward band on a half-resolution row composites column pairs
      // short of each arm's last pair, then that pair column by column.
      const bool tail = forward && half_res;
      for (int part = 0; part <= int(tail); ++part) {
        const bool pairs = half_res && part == 0;
        const int stride = pairs ? 2 : 1;
        // A half-resolution sample sits between its column pair, so it is
        // the pair's box average rather than the even column alone.
        const float lane_offset = pairs ? 0.5f : 0.0f;
        const ColumnRuns row_runs =
            tail ? tail_runs(pass_columns, part) : pass_columns;
        if (field_y1 < polar.north_rings || field_y0 >= polar.south_ring) {
          if constexpr (PAIR_PIXELS)
            composite_polar_row(ctx, row_runs, y, field_y0, field_y1, wy0, wy1,
                                output, stride, lane_offset, pairs,
                                defer_filter, transform_pixel);
          else
            composite_polar_row_flash(ctx, row_runs, y, field_y0, field_y1, wy0,
                                      wy1, output, stride, lane_offset, pairs,
                                      defer_filter, transform_pixel);
        } else {
          for (int r = 0; r < row_runs.count; ++r) {
            const int xs = row_runs.items[r].begin;
            const int xe = row_runs.items[r].end;
            int cx0 = xs / downsample;
            int sub = xs - cx0 * downsample;
            bool cell_stale = true;
            // Whether a pixel of the cell may read past its arm's columns.
            bool cell_resolved = false;
            // The buffer holding every source of an unresolved cell.
            const ::Pixel *source = previous;
            const auto *source_poles = &poles;
            int cell_level = unit_level;
            float leftx = 0.0f, slopex = 0.0f;
            float lefty = 0.0f, slopey = 0.0f;
            auto cell = [&]() {
              HS_PROFILE_DEEP(fb_comp_cell);
              const int cx1 = (cx0 + 1 < coarse_columns) ? cx0 + 1 : 0;
              const int i00 = row0 + cx0, i10 = row0 + cx1;
              const int i01 = row1 + cx0, i11 = row1 + cx1;
              const int d00 = x_offsets[i00];
              const int u10 = unwrap_near(x_offsets[i10], d00);
              const int u01 = unwrap_near(x_offsets[i01], d00);
              const int u11 = unwrap_near(x_offsets[i11], d00);
              const float d10 = static_cast<float>(u10);
              const float d01 = static_cast<float>(u01);
              const float d11 = static_cast<float>(u11);
              leftx = (static_cast<float>(d00) * wy0 + d01 * wy1) *
                      INVERSE_WARP_SCALE;
              slopex = (d10 * wy0 + d11 * wy1) * INVERSE_WARP_SCALE - leftx;
              lefty = (y_offsets[i00] * wy0 + y_offsets[i01] * wy1) *
                      INVERSE_WARP_SCALE;
              slopey = (y_offsets[i10] * wy0 + y_offsets[i11] * wy1) *
                           INVERSE_WARP_SCALE -
                       lefty;
              if (column_time) {
                int column = cx0 * downsample;
                column -= column >= ColumnTime::HALF ? ColumnTime::HALF : 0;
                // No pixel's offset exceeds its cell's largest node offset, so
                // a cell at least that far inside its arm reads within it.
                const int reach =
                    std::max(std::max(std::abs(d00), std::abs(u10)),
                             std::max(std::abs(u01), std::abs(u11)));
                const int room = std::min(column, ColumnTime::HALF - 1 -
                                                      downsample - column);
                cell_resolved = straddled_cells;
                source = previous;
                source_poles = &poles;
                if (reach > room * static_cast<int>(WARP_SCALE)) {
                  // Otherwise the first pixel's source and the last's, past
                  // the cell's far edge, bound the cell's. A forward band reads
                  // the frame being drawn from the arm's last column on, and an
                  // arm's start the stale frame before its first.
                  const float first = column + fminf(leftx, leftx + slopex);
                  const float last =
                      column + downsample + fmaxf(leftx, leftx + slopex);
                  const bool within =
                      first >= 0.0f &&
                      (forward ? last < ColumnTime::HALF - 1.0f
                               : last <= ColumnTime::HALF - 0.5f);
                  if (!within && forward && first >= ColumnTime::HALF) {
                    source = current;
                    source_poles = &ctx.current_poles;
                  } else if (!within && behind && last <= -1.0f) {
                    source = ctx.back;
                    source_poles = &ctx.back_poles;
                  } else {
                    cell_resolved |= !within;
                  }
                }
                cell_level = static_cast<int>(FadeTable::UNIT_LEVEL + 0.5f -
                                              (leftx + 0.5f * slopex) *
                                                  LEVELS_PER_COLUMN);
              }
            };
            for (int x = xs; x < xe;) {
              if (cell_stale) {
                cell();
                cell_stale = false;
              }

              if (cell_resolved) {
                const int cell_end = std::min(xe, x + downsample - sub);
                x = composite_resolved_cell<PAIR_PIXELS>(
                    ctx, transform_pixel, output, x, cell_end, y, sub, stride,
                    lane_offset, {leftx, slopex, lefty, slopey},
                    opaque || defer_filter);
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
                  const float ddx0 = leftx + slopex * fx0;
                  const float ddy0 = lefty + slopey * fx0;
                  const float ddx1 = leftx + slopex * fx1;
                  const float ddy1 = lefty + slopey * fx1;

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
              const float ddx = leftx + slopex * fx;
              const float ddy = lefty + slopey * fx;

              Sample s;
              {
                HS_PROFILE_DEEP(fb_comp_sample);
                s = sample_bilinear_prev(grid.field, source, *source_poles,
                                         x + lane_offset + ddx, y + ddy);
              }
              ::Pixel p(0, 0, 0);
              if (!(skip_near_black(s, black_skips_color))) {
                HS_PROFILE_DEEP(fb_comp_color);
                if constexpr (PAIR_PIXELS)
                  p = transform_pixel_call(transform_pixel, s.r, s.g, s.b,
                                           cell_level);
                else
                  p = transform_pixel(s.r, s.g, s.b, cell_level);
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
      }
      if (ctx.pass == Pass::ARM_STARTS)
        continue;
      if (half_res) {
        HS_PROFILE_DEEP(fb_comp_fill);
        expand_half_res_row<PAIR_PIXELS>(ctx, output, y, transform_pixel);
      }
      if (forward)
        cover_row(ctx, y, output);
      if (filter_output) {
        HS_PROFILE_DEEP(fb_comp_filter);
        filter_infill_row(ctx, y, defer_filter);
      }
    }
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
   * @brief Composites one polar row from its cap-plane offsets, in the hue
   *        composite's placement.
   * @details Each lane converts its own target back to field coordinates, one
   * lane per step.
   */
  template <typename TransformPixelT>
  __attribute__((noinline)) void
  composite_polar_row(const FlushContext &ctx, const ColumnRuns &runs, int y,
                      int field_y0, int field_y1, float wy0, float wy1,
                      ::Pixel *output, int stride, float lane_offset,
                      bool half_res, bool defer_filter,
                      TransformPixelT &transform_pixel) {
    composite_polar_lanes(ctx, runs, y, field_y0, field_y1, wy0, wy1, output,
                          stride, lane_offset, half_res, defer_filter,
                          transform_pixel);
  }

  /** @brief composite_polar_row() in flash, for the composites that run
   *         from flash. */
  template <typename TransformPixelT>
  HS_FLASH_MEMBER void
  composite_polar_row_flash(const FlushContext &ctx, const ColumnRuns &runs,
                            int y, int field_y0, int field_y1, float wy0,
                            float wy1, ::Pixel *output, int stride,
                            float lane_offset, bool half_res, bool defer_filter,
                            TransformPixelT &transform_pixel) {
    composite_polar_lanes(ctx, runs, y, field_y0, field_y1, wy0, wy1, output,
                          stride, lane_offset, half_res, defer_filter,
                          transform_pixel);
  }

  /** @brief The body of composite_polar_row(). */
  template <typename TransformPixelT>
  __attribute__((always_inline)) void
  composite_polar_lanes(const FlushContext &ctx, const ColumnRuns &runs, int y,
                        int field_y0, int field_y1, float wy0, float wy1,
                        ::Pixel *output, int stride, float lane_offset,
                        bool half_res, bool defer_filter,
                        TransformPixelT &transform_pixel) {
    const CoarseGrid &grid = ctx.grid;
    const PolarRings &polar = grid.polar;
    const int downsample = grid.downsample;
    const int coarse_columns = grid.columns;
    const float inverse_downsample = 1.0f / downsample;
    const CapOffset *caps0 =
        ctx.warp.cell_caps + polar.row(field_y0) * coarse_columns;
    const CapOffset *caps1 =
        ctx.warp.cell_caps + polar.row(field_y1) * coarse_columns;
    const bool north = field_y1 < polar.north_rings;
    const float colatitude =
        SphereField::Geometry::row_to_phi(static_cast<float>(y));
    const float cap_angle = north ? colatitude : math::PI_F - colatitude;
    const bool black_skips_color = ctx.black_skips_color;
    const auto blend = blend_alpha(ctx.alpha);
    const bool plain_store = ctx.alpha >= 1.0f || defer_filter;
    for (int r = 0; r < runs.count; ++r) {
      const int xs = runs.items[r].begin;
      const int xe = runs.items[r].end;
      int cx0 = xs / downsample;
      int sub = xs - cx0 * downsample;
      bool cell_stale = true;
      CapCell cell{};
      for (int x = xs; x < xe;) {
        if (cell_stale) {
          const int cx1 = (cx0 + 1 < coarse_columns) ? cx0 + 1 : 0;
          cell = CapPlane::decode_cell(caps0, caps1, cx0, cx1, wy0, wy1);
          cell_stale = false;
        }

        const float fx = (sub + lane_offset) * inverse_downsample;
        const auto at =
            CapPlane::polar_lane(cell, x, cap_angle, north, half_res, fx);
        typename ColumnTime::Read read{at.x, at.y, 1.0f};
        int level = 0;
        bool across = false;
        int arm = 0;
        if (ctx.column_time) {
          const float target = x + lane_offset;
          arm = ColumnTime::arm_start(x);
          const float source =
              target + unwrap_near(at.x - target, 0.0f, static_cast<float>(W));
          across = across_cut(ctx, y, arm, source);
          if (across) {
            read = ColumnTime::read_across(target, source, at.y);
          } else {
            const float column = ColumnTime::arm_column(at.x, arm);
            read = column < ColumnTime::HALF - 0.5f
                       ? ColumnTime::read_within(target, arm, column, at.y)
                       : polar_read_past_cut(cell, x, cap_angle, north,
                                             half_res, fx, target, y, at);
          }
          level = FadeTable::level(read.elapsed);
        }
        const Sample s = sample_read(ctx, arm, across, read.x, read.y);
        ::Pixel p(0, 0, 0);
        if (!(skip_near_black(s, black_skips_color)))
          p = transform_pixel(s.r, s.g, s.b, level);
        ::Pixel &dst = output[x];
        dst = plain_store ? p : blend(dst, p);

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

  /**
   * @brief A polar lane's read for a source past one of its arm's cuts.
   * @param cell The lane's decoded cell.
   * @param x The lane's column.
   * @param cap_angle The row's angle from its pole, in radians.
   * @param north Whether the row's pole is the north pole.
   * @param half_res Whether the lane stands for its column pair.
   * @param fx The lane's fraction across @p cell.
   * @param target The lane's column with its half-resolution offset.
   * @param y The lane's row.
   * @param at The lane's source, past its arm.
   * @return The resolved read.
   */
  static HS_FLASH_MEMBER typename ColumnTime::Read
  polar_read_past_cut(const CapCell &cell, int x, float cap_angle, bool north,
                      bool half_res, float fx, float target, int y,
                      const typename SphereField::Coordinates &at) {
    return ColumnTime::read_past_cut(
        target, static_cast<float>(y), ColumnTime::arm_start(x),
        {target + unwrap_near(at.x - target, 0.0f, static_cast<float>(W)),
         at.y},
        [&]() HS_HOT_INLINE {
          const CapCell twice{{cell.left.u * 2.0f, cell.left.v * 2.0f},
                              {cell.slope.u * 2.0f, cell.slope.v * 2.0f}};
          const auto far =
              CapPlane::polar_lane(twice, x, cap_angle, north, half_res, fx);
          return ColumnSource{
              target + unwrap_near(far.x - target, 0.0f, static_cast<float>(W)),
              far.y};
        });
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
   * @brief Composites a cell's pixels whose sources may lie past their arm's
   * columns, each resolved by FeedbackColumnTime.
   * @tparam PAIR_PIXELS Whether the composite pairs pixels, and so keeps its
   * single-pixel transform out of line.
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
  template <bool PAIR_PIXELS, typename TransformPixelT>
  HS_HOT_FLASH_MEMBER int composite_resolved_cell(
      const FlushContext &ctx, TransformPixelT &transform_pixel,
      ::Pixel *output, int x, int end, int y, int sub, int stride,
      float lane_offset, const CellOffsets &cell, bool plain_store) {
    const float inverse_downsample = 1.0f / ctx.grid.downsample;
    const float row = static_cast<float>(y);
    const uint16_t alpha = frac_to_q16(ctx.alpha);
    auto resolve = [&](int column_x, int column_sub,
                       bool &across) HS_HOT_INLINE {
      const float fx = (column_sub + lane_offset) * inverse_downsample;
      const float ddx = cell.left_x + cell.slope_x * fx;
      const float ddy = cell.left_y + cell.slope_y * fx;
      const float target = column_x + lane_offset;
      const int arm = ColumnTime::arm_start(column_x);
      const float source = target + ddx;
      across = across_cut(ctx, y, arm, source);
      if (across)
        return ColumnTime::read_across(target, source, row + ddy);
      const float column = ColumnTime::arm_column(target + ddx, arm);
      return column < ColumnTime::HALF - 0.5f
                 ? ColumnTime::read_within(target, arm, column, row + ddy)
                 : ColumnTime::read_past_cut(
                       target, row, arm, {target + ddx, row + ddy},
                       [&]() {
                         return ColumnSource{target + 2.0f * ddx,
                                             row + 2.0f * ddy};
                       });
    };
    auto black = [&](const Sample &s) HS_HOT_INLINE {
      return skip_near_black(s, ctx.black_skips_color);
    };
    auto transform = [&](const Sample &s, int level) HS_HOT_INLINE {
      if constexpr (PAIR_PIXELS)
        return transform_pixel_call(transform_pixel, s.r, s.g, s.b, level);
      else
        return transform_pixel(s.r, s.g, s.b, level);
    };
    auto store = [&](int column_x, const ::Pixel &p) HS_HOT_INLINE {
      ::Pixel &dst = output[column_x];
      dst = plain_store ? p : dst.lerp16(p, alpha);
    };
    while (x < end) {
      bool across = false;
      const auto read0 = resolve(x, sub, across);
      const int level0 = FadeTable::level(read0.elapsed);
      const Sample s0 =
          sample_read(ctx, ColumnTime::arm_start(x), across, read0.x, read0.y);
      ::Pixel p(0, 0, 0);
      if (!black(s0))
        p = transform(s0, level0);
      store(x, p);
      x += stride;
      sub += stride;
    }
    return x;
  }

  /**
   * @brief Composites one pass's columns through the frame's color transform.
   * @param ctx Frame context from prepare_flush(), its pass selected.
   */
  void composite_pass(const FlushContext &ctx) {
    switch (ctx.mode) {
    case FlushMode::HUE: {
      const float *k = ctx.column_time ? ctx.fades.hue_k : ctx.hue_k;
      composite_pixels<true>(
          ctx,
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
      composite_plain(ctx);
      break;
    case FlushMode::GENERAL:
      composite_general(ctx);
      break;
    case FlushMode::SKIP:
      break;
    }
  }

  /** @brief Composites through the unrotated hue fade: a per-channel scale. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_plain(const FlushContext &ctx) {
    const float *fades = ctx.column_time ? ctx.fades.fades : &ctx.fade;
    auto plain = [&](float r, float g, float b, int level) {
      const float fade = fades[level];
      return ::Pixel(quantize16(r * fade), quantize16(g * fade),
                     quantize16(b * fade));
    };
    composite_pixels<false>(ctx, plain, plain);
  }

  /** @brief Composites through the Style's type-erased color transform. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_general(const FlushContext &ctx) {
    const float *fades = ctx.column_time ? ctx.fades.fades : &ctx.fade;
    auto general = [&](float r, float g, float b, int level) {
      return feedback_style->color_fn(
          ::Pixel(quantize16(r), quantize16(g), quantize16(b)), fades[level],
          *feedback_style);
    };
    composite_pixels<false>(ctx, general, general);
  }
  HS_O3_END

  // Round-to-nearest fade stalls dim channels short of black.
  static constexpr float NEAR_BLACK = 64.0f;
  static constexpr float WARP_SCALE = 128.0f;
  /** @brief Column offsets in WARP_SCALE units, one full turn apart. */
  static constexpr float WRAP_PERIOD = static_cast<float>(W) * WARP_SCALE;

  /**
   * @brief A row's columns for the current pass.
   * @param ctx Frame context, its pass selected.
   * @param y Row.
   * @return The band's runs for Pass::ALL, else each arm's columns the pass
   *         writes: its first reach columns, its last Ahead::band columns, or
   *         those between.
   */
  static ColumnRuns pass_runs(const FlushContext &ctx, int y) {
    if (ctx.pass == Pass::ALL)
      return ctx.runs;
    constexpr int HALF = ColumnTime::HALF;
    const int reach = ctx.reach[y] > 0 ? ctx.reach[y] : 0;
    const int band = ctx.ahead ? ctx.ahead[y].band : 0;
    int begin = reach, end = HALF - band;
    if (ctx.pass == Pass::ARM_STARTS) {
      begin = 0;
      end = reach;
    } else if (ctx.pass == Pass::FORWARD) {
      begin = HALF - band;
      end = HALF;
    }
    ColumnRuns runs{};
    if (begin < end)
      for (int arm = 0; arm < W; arm += HALF)
        runs.items[runs.count++] = {arm + begin, arm + end};
    return runs;
  }

  /**
   * @brief A half-resolution forward band's columns in two parts.
   * @param runs Each arm's forward band, as pass_runs() gives it.
   * @param part 0 for each band short of its last two columns, 1 for those
   *        two.
   * @return The part's runs.
   */
  static ColumnRuns tail_runs(const ColumnRuns &runs, int part) {
    ColumnRuns out{};
    for (int r = 0; r < runs.count; ++r) {
      const int split = std::max(runs.items[r].begin, runs.items[r].end - 2);
      const ColumnRun piece = part == 0 ? ColumnRun{runs.items[r].begin, split}
                                        : ColumnRun{split, runs.items[r].end};
      if (piece.begin < piece.end)
        out.items[out.count++] = piece;
    }
    return out;
  }

  /**
   * @brief Each row's reach across the arms' cuts.
   * @param grid Coarse layout of the frame.
   * @param warp The frame's populated warp field.
   * @param half_res_sine Latitude sine under which rows composite column
   *        pairs (half_res_limit()).
   * @param reach Receives, per row, the even column count from each arm's
   *        first whose sources may lie behind its backward cut, or -1 where
   *        that exceeds a quarter turn.
   * @param ahead Receives, per row, the forward band of a split frame, or
   *        null for an unsplit one. Pole rows, longitude-filtered infill rows
   *        and rows whose band reaches a quarter turn composite unsplit.
   * @details Each cell bounds its own pixels' sources, so a row's reach and
   * band cover the columns whose cells may read across a cut. A cell's pixels
   * read no further than its nodes' extreme column offsets; a polar lane's
   * source turns no further than its cell's longest cap offset subtends.
   */
  HS_FLASH_MEMBER void measure_reach(const CoarseGrid &grid,
                                     const WarpField &warp, float half_res_sine,
                                     int16_t *reach, Ahead *ahead) const {
    constexpr int HALF = ColumnTime::HALF;
    constexpr float COLUMNS_PER_RADIAN = W / (2.0f * math::PI_F);
    // Columns between a cell's bound and a cut.
    constexpr float MARGIN = 1.0f;
    const PolarRings &polar = grid.polar;
    const int columns = grid.columns;
    const int downsample = grid.downsample;
    struct Span {
      int back;    ///< Columns from each arm's first that may read behind it.
      int forward; ///< Columns before each arm's end that may read past it.
    };
    // The span of rows interpolating rings k and k + 1: each cell whose
    // pixels' sources may cross a cut joins whole, bounded as the composite
    // bounds it, from its first column less its backward offset to past its
    // last column plus its forward offset.
    auto pair_span = [&](int k, int y0, int y1) HS_FLASH_MEMBER {
      const bool north = k + 1 < polar.north_rings;
      const bool polar_pair = north || k >= polar.south_ring;
      // A polar lane's source turns by at most the angle its cap offset
      // subtends at the pair's cap angle nearest the pole.
      const float near =
          north ? SphereField::Geometry::row_to_phi(static_cast<float>(y0))
                : math::PI_F -
                      SphereField::Geometry::row_to_phi(static_cast<float>(y1));
      const CapOffset *caps = polar_pair ? warp.cell_caps : nullptr;
      const int16_t *offsets = warp.x_offsets;
      const int row0 = polar_pair ? polar.row(k) : k, row1 = row0 + 1;
      Span span{0, 0};
      for (int cx = 0; cx < columns; ++cx) {
        const int cx1 = cx + 1 < columns ? cx + 1 : 0;
        const int nodes[4] = {row0 * columns + cx, row0 * columns + cx1,
                              row1 * columns + cx, row1 * columns + cx1};
        float back = 0.0f, forward = 0.0f;
        if (polar_pair) {
          int longest_sq = 0;
          for (int node : nodes)
            longest_sq = std::max(longest_sq, caps[node].u * caps[node].u +
                                                  caps[node].v * caps[node].v);
          const float longest =
              sqrtf(static_cast<float>(longest_sq)) / CapPlane::CAP_SCALE;
          if (!(longest < near))
            return Span{HALF, HALF};
          // The turn's tangent bounds the turn.
          back = forward = longest / sqrtf(near * near - longest * longest) *
                           COLUMNS_PER_RADIAN;
        } else {
          for (int node : nodes) {
            back = fmaxf(back, -offsets[node] / WARP_SCALE);
            forward = fmaxf(forward, offsets[node] / WARP_SCALE);
          }
        }
        for (int arm = 0; arm < W; arm += HALF) {
          const int first = cx * downsample - arm;
          if (first + downsample <= 0 || first >= HALF)
            continue;
          if (first - back < MARGIN)
            span.back = std::max(span.back, first + downsample);
          if (first + downsample + forward > HALF - 1 - MARGIN)
            span.forward = std::max(span.forward, HALF - first);
        }
      }
      return span;
    };
    auto even = [](int span) { return (span + 1) & ~1; };
    int pair_index = 0;
    auto lower_ring = grid.field.ring(0);
    auto upper_ring = grid.field.next_ring(lower_ring);
    Span pair = pair_span(0, lower_ring.y, upper_ring.y);
    for (int y = 0; y < H; ++y) {
      if (y > upper_ring.y && upper_ring.y < H - 1) {
        lower_ring = upper_ring;
        upper_ring = grid.field.next_ring(lower_ring);
        pair = pair_span(++pair_index, lower_ring.y, upper_ring.y);
      }
      const int back_columns = even(pair.back);
      reach[y] = back_columns > HALF / 2 ? int16_t(-1)
                                         : static_cast<int16_t>(back_columns);
      if (!ahead)
        continue;
      const bool pole_row = (SphereField::HAS_NORTH_POLE && y == 0) ||
                            (SphereField::HAS_SOUTH_POLE && y == H - 1);
      const int band = even(pair.forward);
      const bool split =
          band > 0 && band < HALF / 2 && !pole_row && !filtered_infill(grid, y);
      const bool half_res = SphereField::latitude_sine(y) < half_res_sine;
      ahead[y] = split ? Ahead{static_cast<int16_t>(band),
                               static_cast<int16_t>(HALF - band - half_res)}
                       : Ahead{0, static_cast<int16_t>(HALF)};
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
   * @tparam PAIR_PIXELS Whether the composite pairs pixels, and so keeps its
   *         colour transform out of line.
   * @tparam TransformPixelT Callable (r, g, b, level) -> ::Pixel.
   * @param ctx Frame context, its pass selected.
   * @param output Row of pair samples, as for reconstruct_half_res_run().
   * @param y Row index.
   * @param transform_pixel Ages a sample by one frame at the unit level.
   * @details Under column time, a full non-pole row expands arm by arm.
   * Beside a cut, the pair lit just before an arm's first column is the other
   * arm's last pair of the previous frame; the pair lit just after its last
   * column is the other arm's first pair a frame on, aged one frame from this
   * one. A split row expands each arm short of its forward band, keeping the
   * last pair's sample in the column after it; the forward pass composites
   * the band in pairs up to its last pair, which it composites column by
   * column, and expands the pairs between the kept sample and that last pair.
   * Other rows expand run by run.
   */
  template <bool PAIR_PIXELS, typename TransformPixelT>
  HS_HOT_FLASH_MEMBER static void
  expand_half_res_row(const FlushContext &ctx, ::Pixel *output, int y,
                      TransformPixelT &transform_pixel) {
    const ColumnRuns &runs = ctx.runs;
    // Every column of a pole row is one point.
    const bool pole_row = (SphereField::HAS_NORTH_POLE && y == 0) ||
                          (SphereField::HAS_SOUTH_POLE && y == H - 1);
    if (!ctx.column_time || pole_row || runs.count != 2 ||
        runs.items[0].begin != 0 || runs.items[1].end != W) {
      for (int r = 0; r < runs.count; ++r)
        reconstruct_half_res_run(output, runs.items[r].begin,
                                 runs.items[r].end);
      return;
    }
    constexpr int HALF = ColumnTime::HALF;
    constexpr int LEVEL = FadeTable::UNIT_LEVEL;
    auto aged = [&](const ::Pixel &p) {
      if constexpr (PAIR_PIXELS)
        return transform_pixel_call(transform_pixel, p.r, p.g, p.b, LEVEL);
      else
        return transform_pixel(p.r, p.g, p.b, LEVEL);
    };
    const ::Pixel *previous_row = ctx.previous + y * W;
    const ::Pixel before[2] = {
        pair_average(previous_row[W - 2], previous_row[W - 1]),
        pair_average(previous_row[HALF - 2], previous_row[HALF - 1])};
    const int band = ctx.ahead ? ctx.ahead[y].band : 0;
    if (band == 0) {
      const ::Pixel after[2] = {aged(output[HALF]), aged(output[0])};
      reconstruct_half_res_run(output, 0, HALF, &before[0], &after[0]);
      reconstruct_half_res_run(output, HALF, W, &before[1], &after[1]);
      return;
    }
    for (int arm = 0, side = 0; arm < W; arm += HALF, ++side) {
      const int first = arm + HALF - band;
      if (ctx.pass == Pass::FORWARD) {
        const ::Pixel &last = ctx.kept[2 * y + side];
        const int tail = arm + HALF - 2;
        const ::Pixel after = pair_average(output[tail], output[tail + 1]);
        output[first - 1] =
            blend_quarter(last, first < tail ? output[first] : after);
        reconstruct_half_res_run(output, first, tail, &last, &after);
      } else {
        output[first - 1] = output[first - 2];
        reconstruct_half_res_run(output, arm, first - 1, &before[side]);
      }
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
    constexpr int HALF = PERIOD / 2;
    const int delta = v - reference;
    return delta > HALF ? v - PERIOD : (delta < -HALF ? v + PERIOD : v);
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

  /** @brief One bilinear sample on the [0, 65535] scale, unquantized. */
  struct Sample {
    float r;
    float g;
    float b;
  };

  /**
   * @brief Whether a read lies across one of its target arm's cuts in a
   *        buffer the current pass holds.
   * @param ctx Frame context, its pass selected.
   * @param y Target row.
   * @param arm First column of the target's arm.
   * @param x Source column, unwrapped near the target.
   * @return True behind the backward cut on a row whose reach the stale frame
   *         serves, and from the arm's last column on in the forward pass.
   */
  static HS_HOT_INLINE bool across_cut(const FlushContext &ctx, int y, int arm,
                                       float x) {
    return (ctx.back && ctx.reach[y] >= 0 && x < arm) ||
           (ctx.pass == Pass::FORWARD && x >= arm + ColumnTime::HALF - 1);
  }

  /**
   * @brief Samples a read from the buffers holding its taps.
   * @param ctx Frame context.
   * @param arm First column of the target's arm.
   * @param across Whether the read lies across a cut (across_cut()).
   * @param x Source column, unwrapped near the target.
   * @param y Source row.
   * @return The bilinear sample.
   */
  HS_HOT_INLINE Sample sample_read(const FlushContext &ctx, int arm,
                                   bool across, float x, float y) const {
    return across ? sample_across(ctx, arm, x, y)
                  : sample_bilinear_prev(ctx.grid.field, ctx.previous,
                                         ctx.poles, x, y);
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
  HS_FLASH_MEMBER Sample sample_across(const FlushContext &ctx, int arm,
                                       float x, float y) const {
    const SphereField &field = ctx.grid.field;
    const bool behind = x < arm;
    const ::Pixel *low = behind ? ctx.back : ctx.previous;
    const ::Pixel *high = behind ? ctx.previous : ctx.current;
    const auto &low_poles = behind ? ctx.back_poles : ctx.poles;
    const auto &high_poles = behind ? ctx.poles : ctx.current_poles;
    const float first =
        static_cast<float>(behind ? arm : arm + ColumnTime::HALF);
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
   * @brief Blends a tap near an arm's end into its row's forward-band
   *        coverage when it lands in the band, or passes it on.
   * @param x Column coordinate as plotted.
   * @param column @p x rounded and wrapped as the canvas writes it.
   * @param y Row coordinate.
   * @param color,age,alpha Shading registers.
   * @param pass Downstream 2D callback, the canvas write.
   * @details The color blends as the canvas write blends it, and the
   * transmittance takes the same blend toward zero, so finish() can composite
   * the band beneath both.
   */
  template <typename PassFnT>
  HS_FLASH_MEMBER void plot_near_arm_end(float x, int column, float y,
                                         const ::Pixel &color, float age,
                                         float alpha, PassFnT &pass) {
    const int row = Filter::round_row<H>(y);
    const int arm = ColumnTime::arm_start(column);
    if (row < 0 || row >= H || column - arm < held->ahead[row].held_from) {
      pass(x, y, color, age, alpha);
      return;
    }
    ::Pixel &dst = held->current[row * W + column];
    uint16_t &transmittance = held->cover[cover_index(row, arm, column)];
    if (alpha >= 1.0f) {
      dst = color;
      transmittance = 0;
      return;
    }
    const uint16_t weight = frac_to_q16(alpha);
    dst = dst.lerp16(color, weight);
    transmittance = lerp_q16(transmittance, 0, weight);
  }

  /** @brief Transmittance slots per row and arm: a band is even and under a
   *  quarter turn, and a half-resolution row holds one column more. */
  static constexpr int COVER_STRIDE = ColumnTime::HALF / 2;

  /**
   * @brief Transmittance slot of a held column.
   * @param row Row.
   * @param arm First column of the column's arm.
   * @param column Held column.
   * @return Index into FlushContext::cover.
   */
  static HS_HOT_INLINE int cover_index(int row, int arm, int column) {
    return (2 * row + (arm != 0)) * COVER_STRIDE +
           (arm + ColumnTime::HALF - 1 - column);
  }

  /**
   * @brief Clears a split frame's held columns to empty coverage, keeping
   *        each half-resolution row's last REST pair sample.
   * @param ctx Held frame context with its coverage allocated.
   */
  HS_FLASH_MEMBER static void open_coverage(FlushContext &ctx) {
    constexpr int HALF = ColumnTime::HALF;
    for (int y = 0; y < H; ++y) {
      const Ahead ahead = ctx.ahead[y];
      if (ahead.held_from >= HALF)
        continue;
      const bool half_res = ahead.held_from < HALF - ahead.band;
      ::Pixel *row = ctx.current + y * W;
      for (int arm = 0, side = 0; arm < W; arm += HALF, ++side) {
        const int first = arm + ahead.held_from;
        if (half_res)
          ctx.kept[2 * y + side] = row[first];
        for (int x = first; x < arm + HALF; ++x) {
          row[x] = ::Pixel(0, 0, 0);
          ctx.cover[cover_index(y, arm, x)] = 0xFFFF;
        }
      }
    }
  }

  /**
   * @brief Composites one row's held columns beneath their coverage.
   * @param ctx Held frame context.
   * @param y Row.
   * @param under Row composited for the held columns, indexed by column.
   */
  HS_FLASH_MEMBER static void cover_row(const FlushContext &ctx, int y,
                                        const ::Pixel *under) {
    constexpr int HALF = ColumnTime::HALF;
    const int from = ctx.ahead[y].held_from;
    if (from >= HALF)
      return;
    ::Pixel *row = ctx.current + y * W;
    auto add = [](uint16_t covered, uint16_t below, uint16_t transmittance) {
      return static_cast<uint16_t>(std::min<uint32_t>(
          covered + lerp_q16(0, below, transmittance), 65535u));
    };
    for (int arm = 0; arm < W; arm += HALF)
      for (int x = arm + from; x < arm + HALF; ++x) {
        const uint16_t t = ctx.cover[cover_index(y, arm, x)];
        ::Pixel &dst = row[x];
        // A column no tap reached holds no color.
        if (t == 0xFFFF)
          dst = under[x];
        else
          dst = ::Pixel(add(dst.r, under[x].r, t), add(dst.g, under[x].g, t),
                        add(dst.b, under[x].b, t));
      }
  }

  /** @brief Traps while a split frame awaits finish(). */
  HS_FLASH_MEMBER void require_unheld() const {
    HS_CHECK(!held, "feedback: a split frame was never finished");
  }

  /** @brief Releases the held frame and its scratch. */
  HS_FLASH_MEMBER void end_split() {
    held = nullptr;
    held_floor = ColumnTime::HALF;
    scratch_arena_a.set_offset(held_mark);
  }

  ::Feedback::Style *feedback_style; /**< Bound feedback Style (non-owning). */
  bool enabled = true;               /**< When false, flush() is skipped. */
  /** Split frame flush_split() holds for finish(), or null. */
  FlushContext *held = nullptr;
  /** scratch_arena_a offset finish() rewinds to. */
  size_t held_mark = 0;
  /** Least Ahead::held_from of the held frame; HALF holds no column. */
  int held_floor = ColumnTime::HALF;
  WarpCache warp_cache; /**< Persistent warp field and lattice origins. */
  /** @brief Persistent measure_reach() spans of the cached warp field. */
  struct Spans {
    int16_t *reach = nullptr;   ///< Per row, the reach behind each arm.
    Ahead *ahead = nullptr;     ///< Per row, the forward band.
    float pole_half_res = 0.0f; ///< Style::pole_half_res they were measured at.
    bool valid = false;         ///< Whether they belong to the cached field.
  } spans;
  FadeTable fade_table; /**< Persistent per-elapsed-time transforms. */

public:
  /** @brief Cap-cell budget; runtime geometry reserves a full-ring upper bound. */
#if HS_RUNTIME_DISPLAY_GEOMETRY
  static constexpr int CACHE_CAP_CELLS = CACHE_CELLS;
#else
  static constexpr int CACHE_CAP_CELLS =
      polar_rings(CACHE_FIELD, CACHE_FIELD.ring_count()).rows() * CACHE_COLUMNS;
#endif

  /** @brief Lattice-origin budget; runtime geometry reserves the cell upper bound. */
#if HS_RUNTIME_DISPLAY_GEOMETRY
  static constexpr int CACHE_ORIGIN_SAMPLES = CACHE_CELLS;
#else
  static constexpr int CACHE_ORIGIN_SAMPLES = CACHE_FIELD.sample_count();
#endif

  /** @brief Persistent bytes init_storage() reserves over CACHE_CELLS: two
   *  int16 warp fields, the lattice's projected origins, the polar cells'
   *  int16 cap-plane offsets, and the fade table. */
  static constexpr size_t STORAGE_BYTES =
      CACHE_CELLS * 2 * sizeof(int16_t) +
      CACHE_ORIGIN_SAMPLES * sizeof(typename SphereField::Coordinates) +
      CACHE_CAP_CELLS * 2 * sizeof(int16_t) + FadeTable::BYTES +
      H * (sizeof(int16_t) + sizeof(Ahead));
  /** @brief Arena allocations init_storage() makes; each pads its share of
   *  STORAGE_BYTES by less than alignof(std::max_align_t). */
  static constexpr size_t STORAGE_ALLOCATIONS = 8;

  /**
   * @brief scratch_arena_a bytes a split frame holds from flush_split() to
   *        finish() with its warp field and fade table in storage.
   * @return The frame context, filter row, per-row reaches and bands, and the
   *         forward bands' coverage and kept samples, each allocation padded.
   */
  static constexpr size_t split_scratch_bytes() {
    constexpr size_t PAD = alignof(std::max_align_t);
    return sizeof(FlushContext) + W * sizeof(::Pixel) +
           H * (sizeof(int16_t) + sizeof(Ahead)) +
           2 * H * COVER_STRIDE * sizeof(uint16_t) + 2 * H * sizeof(::Pixel) +
           6 * PAD;
  }
};

} // namespace Pixel

} // namespace Filter
