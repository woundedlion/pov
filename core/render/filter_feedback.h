/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include "math/spherical_field.h"
#include "engine/styles.h"
#include "render/filter.h"

/**
 * @file filter_feedback.h
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
 * field, then interpolated within and between rings.
 * flush() iterates the full pixel grid within the active clip band. TERMINAL:
 * flush() composites directly into the Canvas rather than re-emitting downstream,
 * so it must be the last Pipeline stage. The effect must call flush() BEFORE
 * the frame's plot() calls (see `terminal_replaces`); flushing last, as a
 * non-replacing terminal permits, blanks the frame at alpha >= 1.
 *
 * The warp is stored as equirect pixel offsets, so it is anisotropic near the
 * poles: trails pinch and a low-frequency noise cap lenses over a few rows
 * (docs/feedback_pole_isotropy_plan.md).
 */
template <int W, int H> class Feedback : public Is2DWithHistory {
  using SphereField = hs::SphericalFieldLayout<W, H>;

  /** @brief Coarse grid downsample the warp cache is sized for (the default
   *  Style's; every preset keeps it). Other values render uncached. */
  static constexpr int CACHE_DOWNSAMPLE = ::Feedback::Style{}.downsample;
  static constexpr int CACHE_SOUTH_INFILL =
      CACHE_DOWNSAMPLE > hs::H_OFFSET ? CACHE_DOWNSAMPLE - hs::H_OFFSET : 0;
  static constexpr int CACHE_COLUMNS = W / CACHE_DOWNSAMPLE;
  /** @brief Cell count of the cached spherical warp field. */
  static constexpr int CACHE_CELLS =
      CACHE_COLUMNS * SphereField(CACHE_DOWNSAMPLE, CACHE_DOWNSAMPLE,
                                  CACHE_SOUTH_INFILL, CACHE_COLUMNS)
                          .ring_count();

public:
  static constexpr int domain_rank = 2;
  /** @brief Marks this as terminal: flush() writes the Canvas directly. */
  static constexpr bool is_terminal = true;
  /** @brief Opaque store owns the frame: no history stage may precede it. */
  static constexpr bool terminal_replaces = true;

  // Covers only the default-constructed Style; a runtime-swapped style's
  // downsample is validated in flush() (the HS_CHECK below).
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
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

  /** @brief Persistent bytes init_storage() reserves (two int16 warp fields). */
  static constexpr size_t STORAGE_BYTES = 2 * CACHE_CELLS * sizeof(int16_t);

  /**
   * @brief Allocates the warp-field cache from the persistent arena.
   * @param arena Persistent arena supplying 2 * CACHE_CELLS int16 slots.
   * @details Must be called from effect init(), not the constructor (arenas
   * aren't ready yet), and again after any compaction that resets the arena —
   * the cache is derived data, so it just re-populates on the next flush.
   * Without storage every flush renders uncached.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
    cached_warp_x = arena.allocate_n<int16_t>(CACHE_CELLS);
    cached_warp_y = arena.allocate_n<int16_t>(CACHE_CELLS);
    warp_cache_valid = false;
#ifndef NDEBUG
    stamp.record(arena);
#endif
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

    const Animation::NoiseParams *bound = feedback_style->noise;
    HS_CHECK(!bound || (bound->amplitude == feedback_style->amplitude &&
                        bound->frequency == feedback_style->frequency &&
                        bound->speed == feedback_style->speed &&
                        bound->scale == feedback_style->scale),
             "feedback style scalars never reached the bound NoiseParams; call "
             "Style::sync_noise() after changing them");

    const CoarseGrid grid = make_coarse_grid(cv);

    {
      HS_PROFILE(feedback_litscan);
      if (!any_pixel_lit(cv))
        return;
    }

    const RenderBand band = make_render_band(cv.clip(), grid);

    ScratchScope scope(scratch_arena_a);
    WarpField warp = select_warp_field(scope.get_arena(), grid, band);
    populate_warp_field(grid, band, warp);
    ::Pixel *filtered_row = scope.get_arena().allocate_n<::Pixel>(W);
    composite_previous_frame(cv, alpha, grid, band, warp, filtered_row);
  }

private:
  struct CoarseGrid {
    int downsample;
    SphereField field;
    int columns;
    int field_rows;
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

  struct WarpField {
    int16_t *x_offsets;
    int16_t *y_offsets;
    WarpControl *controls;
    bool needs_population;
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
    const int south_infill = std::max(downsample - hs::H_OFFSET, 0);
    const SphereField field(downsample, downsample, south_infill, columns);
    return {downsample, field, columns, field.ring_count()};
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

  __attribute__((always_inline)) WarpField select_warp_field(
      Arena &scratch, const CoarseGrid &grid, const RenderBand &band) {
    const bool stock_transform =
        feedback_style->space_fn == &::Feedback::noise_warp ||
        feedback_style->space_fn == &::Feedback::melt_warp;
    const bool cacheable = cached_warp_x && !band.x_clip.active &&
                           grid.downsample == CACHE_DOWNSAMPLE &&
                           stock_transform;
    if (cacheable)
      check_storage_alive();

    if (!cacheable) {
      const int cells = grid.field_rows * grid.columns;
      return {scratch.allocate_n<int16_t>(cells),
              scratch.allocate_n<int16_t>(cells),
              scratch.allocate_n<WarpControl>(grid.field.sample_count()), true};
    }

    const Animation::NoiseParams *noise = feedback_style->noise;
    const WarpKey key{feedback_style->space_fn,  noise,
                      noise_config_hash(noise),  feedback_style->amplitude,
                      feedback_style->frequency, feedback_style->speed,
                      feedback_style->scale,     noise ? noise->time : 0.0f,
                      band.field_y_begin,        band.field_y_end};

    const bool needs_population = !(warp_cache_valid && key == cached_warp_key);
    cached_warp_key = key;
    warp_cache_valid = true;
    return {cached_warp_x, cached_warp_y,
            needs_population
                ? scratch.allocate_n<WarpControl>(grid.field.sample_count())
                : nullptr,
            needs_population};
  }

  __attribute__((always_inline)) void
  populate_warp_field(const CoarseGrid &grid, const RenderBand &band,
                      const WarpField &warp) {
    HS_PROFILE(feedback_populate);
    if (!warp.needs_population)
      return;

    hs::SphericalField<WarpControl, W, H> compact(warp.controls, grid.field);
    compact.populate(
        band.field_y_begin, band.field_y_end,
        [&](const Vector &position,
            const typename SphereField::Coordinates &point) {
          Vector distorted;
          {
            HS_PROFILE_DEEP(fb_pop_warp);
            distorted = feedback_style->space_fn(position, *feedback_style);
          }
          HS_PROFILE_DEEP(fb_pop_project);
          // Both ends go through project() so its fast-trig error cancels; a
          // pole row's single backing vector has no azimuth to round-trip.
          // Rings stop at y = H - 1, which is the south pole only when
          // H_OFFSET == 0; the device's sub-pole rows carry no field samples.
          // h_offset_renorm_check compiles this file with the device H_OFFSET.
          bool pole_row = point.y == 0.0f;
          if constexpr (hs::H_OFFSET == 0) {
            constexpr float SOUTH_POLE_ROW = H + hs::H_OFFSET - 1;
            pole_row = pole_row || point.y == SOUTH_POLE_ROW;
          }
          const auto projected = grid.field.project(distorted);
          const auto origin = pole_row ? point : grid.field.project(position);
          float x_offset = projected.x - origin.x;
          const float y_offset = projected.y - origin.y;
          if (x_offset > W * 0.5f)
            x_offset -= W;
          else if (x_offset < -W * 0.5f)
            x_offset += W;
          return WarpControl{static_cast<int16_t>(hs::clamp(
                                 x_offset * WARP_SCALE, -32767.0f, 32767.0f)),
                             static_cast<int16_t>(hs::clamp(
                                 y_offset * WARP_SCALE, -32767.0f, 32767.0f))};
        });

    constexpr float WRAP_PERIOD = static_cast<float>(W) * WARP_SCALE;
    constexpr float HALF_WRAP_PERIOD = WRAP_PERIOD * 0.5f;
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
          const WarpControl a = warp.controls[longitude.left];
          WarpControl b = warp.controls[longitude.right];
          float bx = b.x;
          bx += (bx - a.x > HALF_WRAP_PERIOD)
                    ? -WRAP_PERIOD
                    : (bx - a.x < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
          const int index = field_y * grid.columns + coarse_x;
          // Keep the stored offset canonical: the seam correction can lift the
          // interpolant a full period out, past both the int16_t range and the
          // single-step wrap sample_bilinear contracts for.
          const float offset_x =
              hs::lerp(static_cast<float>(a.x), bx, longitude.mix);
          warp.x_offsets[index] = static_cast<int16_t>(
              offset_x > HALF_WRAP_PERIOD
                  ? offset_x - WRAP_PERIOD
                  : (offset_x < -HALF_WRAP_PERIOD ? offset_x + WRAP_PERIOD
                                                  : offset_x));
          warp.y_offsets[index] = static_cast<int16_t>(hs::lerp(
              static_cast<float>(a.y), static_cast<float>(b.y), longitude.mix));
        }
      }
    }
  }

  __attribute__((always_inline)) void
  composite_previous_frame(Canvas &cv, float alpha, const CoarseGrid &grid,
                           const RenderBand &band, const WarpField &warp,
                           ::Pixel *filtered_row) {
    const int downsample = grid.downsample;
    const int coarse_columns = grid.columns;
    const int row_begin = band.y_begin;
    const int row_end = band.y_end;
    const int16_t *x_offsets = warp.x_offsets;
    const int16_t *y_offsets = warp.y_offsets;
    constexpr float INVERSE_WARP_SCALE = 1.0f / WARP_SCALE;
    const float inverse_downsample = 1.0f / grid.downsample;
    const float fade = feedback_style->fade;
    HS_CHECK(std::isfinite(fade) && fade >= 0.0f,
             "feedback fade %f must be finite and non-negative", fade);
    feedback_style->sync_hue();
    const bool black_skips_color =
        feedback_style->color_fn == &::Feedback::hue_fade;
    // Round-to-nearest fade otherwise makes sub-50 channels immortal at
    // FADE_MAX.
    constexpr float NEAR_BLACK = 64.0f;
    const auto blend = blend_alpha(alpha);
    const bool opaque = alpha >= 1.0f;
    constexpr float WRAP_PERIOD = static_cast<float>(W) * WARP_SCALE;
    constexpr float HALF_WRAP_PERIOD = WRAP_PERIOD * 0.5f;
    const ::Pixel *previous = cv.prev_data();
    ::Pixel *current = cv.data();
    ::Pixel poles[SphereField::POLE_COUNT];
    poles[0] = select_pole_sample(previous);
    if constexpr (SphereField::POLE_COUNT == 2)
      poles[1] = select_pole_sample(previous + (H - 1) * W);
    const ColumnRuns runs = make_column_runs(band.x_clip);
    int field_y0 = band.field_y_begin;
    int field_y1 = field_y0 + (field_y0 < band.field_y_end ? 1 : 0);
    const auto control_ring0 = grid.field.ring(field_y0);
    auto control_ring1 = control_ring0;
    if (field_y1 > field_y0)
      control_ring1 = grid.field.next_ring(control_ring0);
    int control_y0 = control_ring0.y;
    int control_y1 = control_ring1.y;
    auto composite_pixels = [&](auto &&transform_pixel, auto &&transform_pair,
                                auto pair_pixels) {
      constexpr bool PAIR_PIXELS = decltype(pair_pixels)::value;
      for (int y = row_begin; y < row_end; ++y) {
        const int row = y * W;
        // The infill bands are dense on both sides, but the southern one backs
        // a pole only when H_OFFSET is 0; the device's last row is mid-latitude.
        bool infill_band = y > 0 && y < downsample;
        if constexpr (hs::H_OFFSET == 0)
          infill_band = infill_band || (y >= H - downsample && y < H - 1);
        const bool filter_output = !band.x_clip.active && infill_band &&
                                   grid.field.longitude_filter_width(y) > 1;
        const bool defer_filter = filter_output && !opaque;
        ::Pixel *output = defer_filter ? filtered_row : current + row;
        while (y > control_y1 && field_y1 < band.field_y_end) {
          field_y0 = field_y1;
          control_y0 = control_y1;
          ++field_y1;
          control_ring1 = grid.field.next_ring(control_ring1);
          control_y1 = control_ring1.y;
        }
        // The last ring lands on the band's last row; short of it the weights
        // below extrapolate off a stale control pair.
        HS_CHECK(y <= control_y1,
                 "feedback warp row %d past last control row %d", y,
                 control_y1);
        // Interpolating outside the populated band silently corrupts pixels.
        HS_CHECK(field_y0 >= band.field_y_begin && field_y1 <= band.field_y_end,
                 "feedback warp row %d outside populated band [%d,%d]",
                 field_y1, band.field_y_begin, band.field_y_end);
        const int control_height = control_y1 - control_y0;
        const float fy =
            control_height > 0
                ? static_cast<float>(y - control_y0) / control_height
                : 0.0f;
        const float wy0 = 1.0f - fy, wy1 = fy;
        const int row0 = field_y0 * coarse_columns;
        const int row1 = field_y1 * coarse_columns;

        for (int r = 0; r < runs.count; ++r) {
          const int xs = runs.items[r].begin;
          const int xe = runs.items[r].end;
          int cx0 = xs / downsample;
          int sub = xs - cx0 * downsample;
          float leftx = 0.0f, slopex = 0.0f;
          float lefty = 0.0f, slopey = 0.0f;
          auto cell = [&]() {
            HS_PROFILE_DEEP(fb_comp_cell);
            const int cx1 = (cx0 + 1 < coarse_columns) ? cx0 + 1 : 0;
            const int i00 = row0 + cx0, i10 = row0 + cx1;
            const int i01 = row1 + cx0, i11 = row1 + cx1;
            float d00 = x_offsets[i00], d10 = x_offsets[i10];
            float d01 = x_offsets[i01], d11 = x_offsets[i11];
            d10 += (d10 - d00 > HALF_WRAP_PERIOD)
                       ? -WRAP_PERIOD
                       : (d10 - d00 < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
            d01 += (d01 - d00 > HALF_WRAP_PERIOD)
                       ? -WRAP_PERIOD
                       : (d01 - d00 < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
            d11 += (d11 - d00 > HALF_WRAP_PERIOD)
                       ? -WRAP_PERIOD
                       : (d11 - d00 < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
            leftx = (d00 * wy0 + d01 * wy1) * INVERSE_WARP_SCALE;
            slopex = (d10 * wy0 + d11 * wy1) * INVERSE_WARP_SCALE - leftx;
            lefty = (y_offsets[i00] * wy0 + y_offsets[i01] * wy1) *
                    INVERSE_WARP_SCALE;
            slopey = (y_offsets[i10] * wy0 + y_offsets[i11] * wy1) *
                         INVERSE_WARP_SCALE -
                     lefty;
          };
          if (sub != 0)
            cell();

          for (int x = xs; x < xe;) {
            if (sub == 0)
              cell();

            if constexpr (PAIR_PIXELS) {
              if (downsample - sub >= 2 && xe - x >= 2) {
                const float fx0 = sub * inverse_downsample;
                const float fx1 = (sub + 1) * inverse_downsample;
                const float ddx0 = leftx + slopex * fx0;
                const float ddy0 = lefty + slopey * fx0;
                const float ddx1 = leftx + slopex * fx1;
                const float ddy1 = lefty + slopey * fx1;

                float sr0, sg0, sb0, sr1, sg1, sb1;
                {
                  HS_PROFILE_DEEP(fb_comp_sample);
                  sample_bilinear_prev(grid.field, previous, poles, x + ddx0,
                                       y + ddy0, sr0, sg0, sb0);
                  sample_bilinear_prev(grid.field, previous, poles,
                                       x + 1 + ddx1, y + ddy1, sr1, sg1, sb1);
                }
                ::Pixel p0(0, 0, 0), p1(0, 0, 0);
                {
                  HS_PROFILE_DEEP(fb_comp_color);
                  transform_pair(sr0, sg0, sb0, sr1, sg1, sb1, p0, p1);
                }
                // Keep both lanes on the paired transform path.
                const bool black0 = black_skips_color && sr0 < NEAR_BLACK &&
                                    sg0 < NEAR_BLACK && sb0 < NEAR_BLACK;
                const bool black1 = black_skips_color && sr1 < NEAR_BLACK &&
                                    sg1 < NEAR_BLACK && sb1 < NEAR_BLACK;
                p0 = black0 ? ::Pixel(0, 0, 0) : p0;
                p1 = black1 ? ::Pixel(0, 0, 0) : p1;

                HS_PROFILE_DEEP(fb_comp_write);
                ::Pixel &dst0 = output[x];
                dst0 = (opaque || defer_filter) ? p0 : blend(dst0, p0);
                ::Pixel &dst1 = output[x + 1];
                dst1 = (opaque || defer_filter) ? p1 : blend(dst1, p1);

                x += 2;
                sub += 2;
                if (sub == downsample) {
                  sub = 0;
                  ++cx0;
                }
                continue;
              }
            }

            const float fx = sub * inverse_downsample;
            const float ddx = leftx + slopex * fx;
            const float ddy = lefty + slopey * fx;

            float sr, sg, sb;
            {
              HS_PROFILE_DEEP(fb_comp_sample);
              sample_bilinear_prev(grid.field, previous, poles, x + ddx,
                                   y + ddy, sr, sg, sb);
            }
            ::Pixel p(0, 0, 0);
            if (!(black_skips_color && sr < NEAR_BLACK && sg < NEAR_BLACK &&
                  sb < NEAR_BLACK)) {
              HS_PROFILE_DEEP(fb_comp_color);
              p = transform_pixel(sr, sg, sb);
            }

            // Black must overwrite the stale double-buffer frame.
            HS_PROFILE_DEEP(fb_comp_write);
            ::Pixel &dst = output[x];
            dst = (opaque || defer_filter) ? p : blend(dst, p);

            ++x;
            if (++sub == downsample) {
              sub = 0;
              ++cx0;
            }
          }
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
    };
    dispatch_color_transform(fade, composite_pixels);
  }

  template <typename CompositeFnT>
  __attribute__((always_inline)) void
  dispatch_color_transform(float fade, CompositeFnT &&composite_pixels) {
    auto composite_scalar = [&](auto &&transform_pixel) {
      composite_pixels(transform_pixel, transform_pixel, std::false_type{});
    };

    HS_PROFILE(feedback_composite);
    const bool hue_identity =
        feedback_style->hue_ca == 1.0f && feedback_style->hue_sa == 0.0f;
    if (feedback_style->color_fn == &::Feedback::hue_fade && !hue_identity) {
      float k[9];
      const float sc = fast_cbrt(fade * (1.0f / 65535.0f));
      for (int i = 0; i < 9; ++i)
        k[i] = feedback_style->hue_k[i] * sc;
      composite_pixels(
          [&](float r, float g, float b) {
            return ::Feedback::hue_fade_apply(k, r, g, b);
          },
          [&](float r0, float g0, float b0, float r1, float g1, float b1,
              ::Pixel &p0, ::Pixel &p1) {
            ::Feedback::hue_fade_apply2(k, r0, g0, b0, r1, g1, b1, p0, p1);
          },
          std::true_type{});
    } else if (feedback_style->color_fn == &::Feedback::hue_fade) {
      auto plain = [&](float r, float g, float b) {
        return ::Pixel(quantize16(r * fade), quantize16(g * fade),
                       quantize16(b * fade));
      };
      composite_scalar(plain);
    } else {
      auto general = [&](float r, float g, float b) {
        return feedback_style->color_fn(
            ::Pixel(quantize16(r), quantize16(g), quantize16(b)), fade,
            *feedback_style);
      };
      composite_scalar(general);
    }
  }
  HS_O3_END

  static constexpr float WARP_SCALE = 128.0f;
  // populate_warp_field() canonicalizes the stored column offset onto
  // [-W*WARP_SCALE/2, W*WARP_SCALE/2] and casts it to int16_t unclamped.
  static_assert(W * WARP_SCALE * 0.5f <= 32767.0f,
                "Feedback<W,H>: canonical warp offset must fit int16_t");
  // The row offset spans the latitude range and is only runtime-clamped, so an
  // over-range displacement would saturate instead of trapping.
  static_assert((H + hs::H_OFFSET - 1) * WARP_SCALE <= 32767.0f,
                "Feedback<W,H>: warp row offset must fit int16_t");

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

  /**
   * @brief Bilinearly samples the Canvas front buffer (previous frame).
   * @param field Spherical topology and interpolation policy.
   * @param prev Base of the previous-frame buffer, row-major with stride W.
   * @param poles Shared values for every aliased column of the pole rows.
   * @param bx Fractional column in [-W, 2W).
   * @param by Fractional row; north crossings reflect with a half-turn.
   * @param r Out: interpolated red on the [0, 65535] scale, unquantized.
   * @param g Out: interpolated green.
   * @param b Out: interpolated blue.
   */
  HS_O3_FN
  void sample_bilinear_prev(const SphereField &field, const ::Pixel *prev,
                            const ::Pixel (&poles)[SphereField::POLE_COUNT],
                            float bx, float by, float &r, float &g,
                            float &b) const {
    field.sample_bilinear_rgb(prev, poles, bx, by, r, g, b);
  }

  /** @brief Quantizes an unclamped [0, 65535]-scale channel to a ::Pixel
   *  component, round-to-nearest; NaN maps to the hi bound. */
  static uint16_t quantize16(float v) {
    // clamp guards the cast against overshoot and maps NaN to the hi bound.
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 65535.0f) + 0.5f);
  }

  /**
   * @brief FNV-1a over the bound noise generator's configuration.
   * @details Seed, noise type and fractal settings feed the warp field but are
   * mirrored neither in Style nor behind an accessor, so the generator's object
   * representation is the only handle on them. Constant while the config is.
   */
  static uint32_t noise_config_hash(const Animation::NoiseParams *noise) {
    static_assert(std::is_trivially_copyable_v<FastNoiseLite>,
                  "hashing FastNoiseLite's bytes requires a trivial layout");
    if (!noise)
      return 0;
    const auto *bytes = reinterpret_cast<const unsigned char *>(&noise->noise);
    uint32_t hash = 2166136261u;
    for (size_t i = 0; i < sizeof(FastNoiseLite); ++i)
      hash = (hash ^ bytes[i]) * 16777619u;
    return hash;
  }

  /** @brief Inputs the coarse warp field is a pure function of (stock
   *  transforms only); equal keys make the cached field reusable. */
  struct WarpKey {
    ::Feedback::SpaceFn space_fn;
    const Animation::NoiseParams *noise;
    uint32_t noise_config;
    float amplitude;
    float frequency;
    float speed;
    float scale;
    float time;
    int field_y_begin;
    int field_y_end;
    bool operator==(const WarpKey &) const = default;
  };

  ::Feedback::Style *feedback_style; /**< Bound feedback Style (non-owning). */
  bool enabled = true;               /**< When false, flush() is skipped. */
  WarpKey cached_warp_key{};         /**< Key for the cached warp field. */
  bool warp_cache_valid = false;    /**< True when the cached field is valid. */
  int16_t *cached_warp_x = nullptr; /**< Arena-owned cached column deltas. */
  int16_t *cached_warp_y = nullptr; /**< Arena-owned cached row deltas. */
#ifndef NDEBUG
  ArenaBlockStamp
      stamp; /**< Arena state when the warp fields were allocated. */
#endif

  /**
   * @brief Debug-only use-after-free check on the arena-owned warp fields.
   * @details A compaction that resets or rewinds the persistent arena without a
   * fresh init_storage() leaves both pointers dangling while warp_cache_valid
   * still reads true, so flush() reads and writes CACHE_CELLS int16s through
   * them.
   */
  void check_storage_alive() const {
#ifndef NDEBUG
    constexpr size_t FIELD_BYTES = CACHE_CELLS * sizeof(int16_t);
    assert(!stamp.arena_reset() &&
           "Pixel::Feedback warp cache use-after-free!");
    assert(!stamp.block_uncovered(cached_warp_x, FIELD_BYTES) &&
           !stamp.block_uncovered(cached_warp_y, FIELD_BYTES) &&
           "Pixel::Feedback warp cache use-after-free (arena rewound below "
           "block)!");
    assert(!stamp.block_reissued(cached_warp_x, FIELD_BYTES) &&
           !stamp.block_reissued(cached_warp_y, FIELD_BYTES) &&
           "Pixel::Feedback warp cache use-after-free (block reclaimed by a "
           "rewind and reissued)!");
#endif
  }
};

} // namespace Pixel

} // namespace Filter
