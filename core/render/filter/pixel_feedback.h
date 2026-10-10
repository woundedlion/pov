/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <algorithm>
#include <bitset>
#include <cstdint>
#include "math/spherical_field.h"
#include "memory.h"
#include "render/filter/feedback_cap_plane.h"
#include "render/filter/feedback_style.h"
#include "render/filter/feedback_warp_cache.h"
#include "render/filter/pipeline.h"

/**
 * @file pixel_feedback.h
 * @brief Filter::Pixel::Feedback: the terminal filter that warps and
 * composites the previous frame.
 */

namespace Filter {

namespace Pixel {

/**
 * @brief Style-aware terminal feedback filter that warps the previous frame.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The Style's spatial warp is computed on a spherical latitude-ring
 * field, then interpolated within and between rings. Terminal: it must be the
 * last Pipeline stage, and Pipeline::begin_frame() must run before the frame's
 * plot() calls; flushing last blanks the frame at alpha >= 1.
 * Rows under POLAR_TARGET_SINE interpolate offsets in their pole's cap plane,
 * where equirect longitude offsets blow up as 1/sin(phi).
 */
template <int W, int H> class Feedback : public Is2DWithHistory {
  using SphereField = hs::SphericalFieldLayout<W, H>;
  using CapPlane = FeedbackCapPlane<W, H>;
  using CapOffset = typename CapPlane::CapOffset;
  using CapPoint = typename CapPlane::CapPoint;
  using CapCell = typename CapPlane::CapCell;
  using WarpCache = FeedbackWarpCache<W, H>;
  using WarpField = typename WarpCache::Buffers;

  /** @brief Coarse grid downsample the warp cache is sized for (the default
   *  Style's). Other values render uncached. */
  static constexpr int CACHE_DOWNSAMPLE = ::Feedback::Style{}.downsample;
  static constexpr int CACHE_SOUTH_INFILL = CACHE_DOWNSAMPLE;
  static constexpr int CACHE_COLUMNS = W / CACHE_DOWNSAMPLE;
  static constexpr SphereField CACHE_FIELD{CACHE_DOWNSAMPLE, CACHE_DOWNSAMPLE,
                                           CACHE_SOUTH_INFILL, CACHE_COLUMNS};
  /** @brief Cell count of the cached spherical warp field; also bounds its
   *  lattice samples, since no ring carries more samples than the grid has
   *  columns. */
  static constexpr int CACHE_CELLS = CACHE_COLUMNS * CACHE_FIELD.ring_count();

public:
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
   * @brief Pass-through: current-frame pixels go straight to the next filter.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    pass(x, y, color, age, alpha);
  }

  /**
   * @brief Enables or disables feedback.
   * @param value When false, flush() is skipped entirely.
   */
  void set_enabled(bool value) { enabled = value; }

  /** @brief Scratch bytes for a full-width uncached flush at downsample ds. */
  static size_t uncached_scratch_bytes(int ds) {
    const int COLUMNS = W / ds;
    const SphereField FIELD(ds, ds, ds, COLUMNS);
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
   * @brief Allocates the warp-field cache from the persistent arena.
   * @param arena Persistent arena supplying STORAGE_BYTES bytes.
   * @details Call from effect init(), not the constructor, and again after
   * any arena reset. Without storage every flush needs
   * uncached_scratch_bytes(downsample) scratch bytes.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
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
    if (!enabled)
      return;

    ScratchScope scope(scratch_arena_a);
    const FlushContext ctx = prepare_flush(cv, alpha, scope.get_arena());
    if (ctx.mode == FlushMode::SKIP)
      return;

    HS_PROFILE(feedback_composite);
    switch (ctx.mode) {
    case FlushMode::HUE: {
      const auto &k = ctx.hue_k;
      composite_pixels<true>(
          ctx,
          [&](float r, float g, float b) {
            return ::Feedback::hue_fade_apply(k, r, g, b);
          },
          [&](float r0, float g0, float b0, float r1, float g1, float b1,
              ::Pixel &p0, ::Pixel &p1) {
            ::Feedback::hue_fade_apply2(k, r0, g0, b0, r1, g1, b1, p0, p1);
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
    ColumnRun items[2];
    int count;
  };

  enum class FlushMode : uint8_t { SKIP, HUE, PLAIN, GENERAL };

  /** @brief Per-frame composite inputs prepare_flush() resolves. */
  struct FlushContext {
    CoarseGrid grid;
    RenderBand band{};
    WarpField warp{};
    ::Pixel *filtered_row = nullptr;
    const ::Pixel *previous = nullptr;
    ::Pixel *current = nullptr;
    ::Pixel poles[SphereField::POLE_STORAGE_COUNT]{};
    ColumnRuns runs{};
    typename SphereField::Ring control_ring1{};
    int field_y0 = 0;
    int field_y1 = 0;
    int control_y0 = 0;
    float alpha = 0.0f;
    float fade = 0.0f;
    float pole_half_res = 0.0f;
    bool black_skips_color = false;
    FlushMode mode = FlushMode::SKIP;
    float hue_k[9] = {};
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
    const int south_infill = downsample;
    const SphereField field(downsample, downsample, south_infill, columns);
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
  make_column_runs(const ClipRegion::XClip &clip) {
    ColumnRuns runs{};
    if (!clip.active) {
      runs.items[runs.count++] = {0, W};
    } else if (clip.wrap) {
      if (clip.re > 0)
        runs.items[runs.count++] = {0, clip.re};
      if (clip.rs < W)
        runs.items[runs.count++] = {clip.rs, W};
    } else {
      runs.items[runs.count++] = {clip.rs, clip.re};
    }
    return runs;
  }

  /**
   * @brief The frame's populated warp field: the cache's when the Style and
   * band allow it, otherwise a scratch field.
   */
  __attribute__((always_inline)) WarpField acquire_warp_field(
      Arena &scratch, const CoarseGrid &grid, const RenderBand &band) {
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

    const Animation::NoiseParams *noise = feedback_style->noise;
    const typename WarpCache::Key key{
        feedback_style->space_fn,  noise,
        noise_config_key(noise),   feedback_style->amplitude,
        feedback_style->frequency, feedback_style->speed,
        feedback_style->scale,     noise ? noise->time : 0.0f,
        band.field_y_begin,        band.field_y_end};

    HS_PROFILE(feedback_populate);
    return warp_cache.acquire(
        cacheable ? &key : nullptr, uncached,
        [&](const WarpField &warp) __attribute__((always_inline)) {
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
              hs::lerp(static_cast<float>(a.x), bx, longitude.mix);
          warp.x_offsets[index] =
              static_cast<int16_t>(unwrap_near(offset_x, 0.0f, WRAP_PERIOD));
          warp.y_offsets[index] = static_cast<int16_t>(hs::lerp(
              static_cast<float>(a.y), static_cast<float>(b.y), longitude.mix));
          if (field_y < grid.polar.north_rings ||
              field_y >= grid.polar.south_ring) {
            const CapOffset da = sample_caps[grid.polar.slot(longitude.left)];
            const CapOffset db = sample_caps[grid.polar.slot(longitude.right)];
            warp.cell_caps[grid.polar.row(field_y) * grid.columns + coarse_x] =
                CapPlane::lerp_offset(da, db, longitude.mix);
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
   * @return The composite context; mode SKIP when the previous frame is black.
   */
  HS_HOT_FLASH_MEMBER FlushContext prepare_flush(Canvas &cv, float alpha,
                                                 Arena &scratch) {
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
    ctx.warp = acquire_warp_field(scratch, grid, band);
    ctx.filtered_row = (!band.x_clip.active &&
                        (grid.downsample > 1 || !SphereField::HAS_NORTH_POLE ||
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
    ctx.runs = make_column_runs(band.x_clip);
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
      const float sc = math::fast_cbrt(fade * (1.0f / 65535.0f));
      for (int i = 0; i < 9; ++i)
        ctx.hue_k[i] = feedback_style->hue_k[i] * sc;
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
   * @param transform_pixel Maps one sampled RGB triple to the output pixel.
   * @param transform_pair Maps two sampled RGB triples to two output pixels.
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
    const ColumnRuns runs = ctx.runs;
    int field_y0 = ctx.field_y0;
    int field_y1 = ctx.field_y1;
    auto control_ring1 = ctx.control_ring1;
    int control_y0 = ctx.control_y0;
    int control_y1 = control_ring1.y;
    // Latitude sine under which two columns subtend less than `pole_half_res`
    // row pitches.
    const float half_res_sine = ctx.pole_half_res * W *
                                SphereField::Geometry::RADIANS_PER_ROW *
                                (1.0f / (4.0f * math::PI_F));
    for (int y = row_begin; y < row_end; ++y) {
      const int row = y * W;
      const bool infill_band =
          (y < downsample && (!SphereField::HAS_NORTH_POLE || y > 0)) ||
          (y >= H - downsample && (!SphereField::HAS_SOUTH_POLE || y < H - 1));
      const bool filter_output = !band.x_clip.active && infill_band &&
                                 grid.field.longitude_filter_width(y) > 1;
      const bool defer_filter = filter_output && !opaque;
      ::Pixel *output = defer_filter ? filtered_row : current + row;
      // The pair expansion needs plain stores, and a longitude-filtered row
      // is reconstructed at its footprint already.
      const bool half_res = opaque && !filter_output &&
                            SphereField::latitude_sine(y) < half_res_sine;
      const int stride = half_res ? 2 : 1;
      // A half-resolution sample sits between its column pair, so it is the
      // pair's box average rather than the even column alone.
      const float lane_offset = half_res ? 0.5f : 0.0f;
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

      if (field_y1 < polar.north_rings || field_y0 >= polar.south_ring) {
        composite_polar_row(ctx, y, field_y0, field_y1, wy0, wy1, output,
                            stride, lane_offset, half_res, defer_filter,
                            transform_pixel);
      } else {
        for (int r = 0; r < runs.count; ++r) {
          const int xs = runs.items[r].begin;
          const int xe = runs.items[r].end;
          int cx0 = xs / downsample;
          int sub = xs - cx0 * downsample;
          bool cell_stale = true;
          float leftx = 0.0f, slopex = 0.0f;
          float lefty = 0.0f, slopey = 0.0f;
          auto cell = [&]() {
            HS_PROFILE_DEEP(fb_comp_cell);
            const int cx1 = (cx0 + 1 < coarse_columns) ? cx0 + 1 : 0;
            const int i00 = row0 + cx0, i10 = row0 + cx1;
            const int i01 = row1 + cx0, i11 = row1 + cx1;
            const int d00 = x_offsets[i00];
            const float d10 =
                static_cast<float>(unwrap_near(x_offsets[i10], d00));
            const float d01 =
                static_cast<float>(unwrap_near(x_offsets[i01], d00));
            const float d11 =
                static_cast<float>(unwrap_near(x_offsets[i11], d00));
            leftx = (static_cast<float>(d00) * wy0 + d01 * wy1) *
                    INVERSE_WARP_SCALE;
            slopex = (d10 * wy0 + d11 * wy1) * INVERSE_WARP_SCALE - leftx;
            lefty = (y_offsets[i00] * wy0 + y_offsets[i01] * wy1) *
                    INVERSE_WARP_SCALE;
            slopey = (y_offsets[i10] * wy0 + y_offsets[i11] * wy1) *
                         INVERSE_WARP_SCALE -
                     lefty;
          };
          for (int x = xs; x < xe;) {
            if (cell_stale) {
              cell();
              cell_stale = false;
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
                  s0 = sample_bilinear_prev(grid.field, previous, poles,
                                            x + lane_offset + ddx0, y + ddy0);
                  s1 = sample_bilinear_prev(grid.field, previous, poles,
                                            x + stride + lane_offset + ddx1,
                                            y + ddy1);
                }
                const bool black0 = black_skips_color && s0.r < NEAR_BLACK &&
                                    s0.g < NEAR_BLACK && s0.b < NEAR_BLACK;
                const bool black1 = black_skips_color && s1.r < NEAR_BLACK &&
                                    s1.g < NEAR_BLACK && s1.b < NEAR_BLACK;
                ::Pixel p0(0, 0, 0), p1(0, 0, 0);
                if (!(black0 && black1)) {
                  HS_PROFILE_DEEP(fb_comp_color);
                  transform_pair_call(transform_pair, s0.r, s0.g, s0.b, s1.r,
                                      s1.g, s1.b, p0, p1);
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
              s = sample_bilinear_prev(grid.field, previous, poles,
                                       x + lane_offset + ddx, y + ddy);
            }
            ::Pixel p(0, 0, 0);
            if (!(black_skips_color && s.r < NEAR_BLACK && s.g < NEAR_BLACK &&
                  s.b < NEAR_BLACK)) {
              HS_PROFILE_DEEP(fb_comp_color);
              p = transform_pixel(s.r, s.g, s.b);
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
      if (half_res) {
        HS_PROFILE_DEEP(fb_comp_fill);
        for (int r = 0; r < runs.count; ++r)
          reconstruct_half_res_run(output, runs.items[r].begin,
                                   runs.items[r].end);
      }
      if (filter_output) {
        HS_PROFILE_DEEP(fb_comp_filter);
        if (!defer_filter)
          std::copy_n(current + row, W, filtered_row);
        grid.field.template reconstruct_longitude_row<PixelAccumulator>(
            filtered_row, y, [&](int x, const ::Pixel &pixel) {
              ::Pixel &dst = current[row + x];
              dst = opaque ? pixel : blend(dst, pixel);
            });
      }
    }
  }

  /**
   * @brief Composites one polar row from its cap-plane offsets.
   * @details Each lane converts its own target back to field coordinates, one
   * lane per step.
   */
  template <typename TransformPixelT>
  __attribute__((noinline)) void
  composite_polar_row(const FlushContext &ctx, int y, int field_y0,
                      int field_y1, float wy0, float wy1, ::Pixel *output,
                      int stride, float lane_offset, bool half_res,
                      bool defer_filter, TransformPixelT &transform_pixel) {
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
    for (int r = 0; r < ctx.runs.count; ++r) {
      const int xs = ctx.runs.items[r].begin;
      const int xe = ctx.runs.items[r].end;
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

        const auto at =
            CapPlane::polar_lane(cell, x, cap_angle, north, half_res,
                                 (sub + lane_offset) * inverse_downsample);
        const Sample s = sample_bilinear_prev(grid.field, ctx.previous,
                                              ctx.poles, at.x, at.y);
        ::Pixel p(0, 0, 0);
        if (!(black_skips_color && s.r < NEAR_BLACK && s.g < NEAR_BLACK &&
              s.b < NEAR_BLACK))
          p = transform_pixel(s.r, s.g, s.b);
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

  /** @brief Composites through the unrotated hue fade: a per-channel scale. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_plain(const FlushContext &ctx) {
    const float fade = ctx.fade;
    auto plain = [&](float r, float g, float b) {
      return ::Pixel(quantize16(r * fade), quantize16(g * fade),
                     quantize16(b * fade));
    };
    composite_pixels<false>(ctx, plain, plain);
  }

  /** @brief Composites through the Style's type-erased color transform. */
  HS_FLASH_MEMBER __attribute__((flatten)) void
  composite_general(const FlushContext &ctx) {
    const float fade = ctx.fade;
    auto general = [&](float r, float g, float b) {
      return feedback_style->color_fn(
          ::Pixel(quantize16(r), quantize16(g), quantize16(b)), fade,
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
   * @brief Expands a half-resolution run's pair samples back to every column.
   * @param output Row whose even columns in [begin, end) hold their pair's
   *        box average, sampled at the pair's midpoint.
   * @param begin First column of the run.
   * @param end One past the run's last column.
   * @details Each column takes three quarters of its own pair's sample and a
   * quarter of the neighbouring pair's on its side. A full row wraps across
   * the seam; a clipped run repeats its end pairs.
   */
  static void reconstruct_half_res_run(::Pixel *output, int begin, int end) {
    const bool full_row = begin == 0 && end == W;
    auto blend = [](const ::Pixel &own, const ::Pixel &other) {
      return ::Pixel(static_cast<uint16_t>((3u * own.r + other.r + 2u) >> 2),
                     static_cast<uint16_t>((3u * own.g + other.g + 2u) >> 2),
                     static_cast<uint16_t>((3u * own.b + other.b + 2u) >> 2));
    };
    // The last pair's right neighbour wraps onto the first sample, which the
    // loop overwrites first.
    const ::Pixel first = output[begin];
    ::Pixel previous = full_row ? output[W - 2 + (W & 1)] : first;
    for (int x = begin; x < end; x += 2) {
      const ::Pixel own = output[x];
      const ::Pixel next = x + 2 < end ? output[x + 2] : full_row ? first : own;
      output[x] = blend(own, previous);
      if (x + 1 < end)
        output[x + 1] = blend(own, next);
      previous = own;
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
  // [-W*WARP_SCALE/2, W*WARP_SCALE/2] and casts it to int16_t unclamped.
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
   * @details Scans only this segment's clip band so another board's lit pixels
   * do not gate this board's flush.
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

  /** @brief Calls @p transform out of line. */
  template <typename TransformPairT>
  __attribute__((noinline)) static void
  transform_pair_call(TransformPairT &transform, float r0, float g0, float b0,
                      float r1, float g1, float b1, ::Pixel &p0, ::Pixel &p1) {
    transform(r0, g0, b0, r1, g1, b1, p0, p1);
  }

  /** @brief Quantizes an unclamped [0, 65535]-scale channel to a ::Pixel
   *  component, round-to-nearest; NaN maps to the hi bound. */
  static uint16_t quantize16(float v) {
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

  ::Feedback::Style *feedback_style; /**< Bound feedback Style (non-owning). */
  bool enabled = true;               /**< When false, flush() is skipped. */
  WarpCache warp_cache; /**< Persistent warp field and lattice origins. */

public:
  /** @brief Cap-cell budget; runtime geometry reserves a full-ring upper bound. */
#if HS_RUNTIME_DISPLAY_GEOMETRY
  static constexpr int CACHE_CAP_CELLS = CACHE_CELLS;
#else
  static constexpr int CACHE_CAP_CELLS =
      polar_rings(CACHE_FIELD, CACHE_FIELD.ring_count()).rows() * CACHE_COLUMNS;
#endif

  /** @brief Persistent bytes init_storage() reserves over CACHE_CELLS: two
   *  int16 warp fields, the lattice's projected origins, and the polar cells'
   *  int16 cap-plane offsets. */
  static constexpr size_t STORAGE_BYTES =
      CACHE_CELLS *
          (2 * sizeof(int16_t) + sizeof(typename SphereField::Coordinates)) +
      CACHE_CAP_CELLS * 2 * sizeof(int16_t);
};

} // namespace Pixel

} // namespace Filter
