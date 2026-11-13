/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file frame.h
 * @brief Feedback frame preparation and warp fields. */

#include "core/render/filter/pixel_feedback.h"

namespace Filter::Pixel {

template <int W, int H>
template <typename PassFnT>
inline HS_HOT_INLINE void Feedback<W, H>::plot(float x, float y,
                                               const ::Pixel &color, float age,
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

template <int W, int H>
inline size_t Feedback<W, H>::uncached_scratch_bytes(int ds) {
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

template <int W, int H>
inline HS_FLASH_MEMBER void Feedback<W, H>::init_storage(Arena &arena) {
  const int cap_cells =
      polar_rings(CACHE_FIELD, CACHE_FIELD.ring_count()).rows() * CACHE_COLUMNS;
  // no death case: the reserve covers every display geometry the cap rows admit
  HS_CHECK(cap_cells <= CACHE_CAP_CELLS,
           "feedback: %d polar cap cells exceed the %d reserved", cap_cells,
           CACHE_CAP_CELLS);
  // Projected from the same incremental-rotation positions populate() hands
  // out, so the per-frame offset subtracts them exactly.
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

HS_O3_BEGIN

template <int W, int H>
inline void Feedback<W, H>::flush(Canvas &cv, float alpha) {
  if (disabled(cv))
    return;
  ScratchScope scope(scratch_arena_a);
  FlushContext ctx = prepare_flush(cv, alpha, scope.get_arena());
  composite_frame(cv, ctx);
}

HS_O3_END

template <int W, int H> inline bool Feedback<W, H>::disabled(Canvas &cv) {
  require_unheld();
  if (enabled)
    return false;
  if (cv.keeps_stale_pixels())
    cv.clear_buffer();
  return true;
}

template <int W, int H>
inline HS_FLASH_MEMBER void Feedback<W, H>::flush_split(Canvas &cv,
                                                        float alpha) {
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

template <int W, int H>
inline HS_FLASH_MEMBER void Feedback<W, H>::finish(Canvas &) {
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

template <int W, int H>
inline __attribute__((always_inline)) typename Feedback<W, H>::CoarseGrid
Feedback<W, H>::make_coarse_grid(const Canvas &cv) const {
  const int downsample = feedback_style->downsample;
  HS_CHECK(downsample > 0 && W % downsample == 0,
           "feedback downsample %d must be > 0 and divide width %d", downsample,
           W);
  HS_CHECK(cv.width() == W, "feedback canvas width %d must equal template W %d",
           cv.width(), W);
  HS_CHECK(cv.height() == H,
           "feedback canvas height %d must equal template H %d", cv.height(),
           H);
  const int columns = W / downsample;
  const SphereField field = field_for(downsample);
  const int rings = field.ring_count();
  return {downsample, field, columns, rings, polar_rings(field, rings)};
}

template <int W, int H>
inline __attribute__((always_inline)) typename Feedback<W, H>::RenderBand
Feedback<W, H>::make_render_band(const ClipRegion &clip,
                                 const CoarseGrid &grid) {
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

template <int W, int H>
inline __attribute__((always_inline)) typename Feedback<W, H>::ColumnRuns
Feedback<W, H>::make_column_runs(const ClipRegion::XClip &clip,
                                 bool split_arms) {
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

template <int W, int H>
inline typename Feedback<W, H>::WarpCache::Key
Feedback<W, H>::warp_cache_key(const RenderBand &band) const {
  const Animation::NoiseParams *noise = feedback_style->noise;
  return {feedback_style->space_fn,  noise,
          noise_config_key(noise),   feedback_style->amplitude,
          feedback_style->frequency, feedback_style->speed,
          feedback_style->scale,     noise ? noise->time : 0.0f,
          band.field_y_begin,        band.field_y_end,
          !Spacetime::flat_surface()};
}

template <int W, int H>
inline __attribute__((always_inline)) typename Feedback<W, H>::WarpField
Feedback<W, H>::acquire_warp_field(Arena &scratch, const CoarseGrid &grid,
                                   const RenderBand &band, bool &reused,
                                   WarpField &prior) {
  const bool stock_transform =
      feedback_style->space_fn == &::Feedback::noise_warp ||
      feedback_style->space_fn == &::Feedback::melt_warp;
  const bool cacheable = warp_cache.ready() && !band.x_clip.active &&
                         grid.downsample == CACHE_DOWNSAMPLE && stock_transform;
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
    uncached.cell_caps =
        polar_cells > 0 ? scratch.allocate_n<CapOffset>(polar_cells) : nullptr;
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

template <int W, int H>
inline __attribute__((always_inline)) void Feedback<W, H>::populate_warp_field(
    const CoarseGrid &grid, const RenderBand &band, const WarpField &warp,
    WarpControl *controls, CapOffset *sample_caps) {
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
        return WarpControl{static_cast<int16_t>(hs::clamp(x_offset * WARP_SCALE,
                                                          -32767.0f, 32767.0f)),
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
        const float offset_y = hs::lerp(static_cast<float>(a.y),
                                        static_cast<float>(b.y), longitude.mix);
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER typename Feedback<W, H>::FlushContext
Feedback<W, H>::prepare_flush(Canvas &cv, float alpha, Arena &scratch,
                              bool split) {
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
  HS_CHECK(fade >= 0.0f && fade <= 1.0f,
           "feedback fade %d/1000 must be in [0, 1]",
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
  const bool clipped = band.x_clip.active || band.y_begin > 0 || band.y_end < H;
  if (ctx.column_time && cv.keeps_stale_pixels() && alpha >= 1.0f && !clipped) {
    ctx.back = ctx.current;
    select_poles(ctx.back, ctx.back_poles);
    // An unchanged cached field keeps its spans.
    const bool stored = reused && spans.reach;
    if (!(stored && spans.valid)) {
      int16_t *reach = stored ? spans.reach : scratch.allocate_n<int16_t>(H);
      measure_reach(grid, ctx.warp, ctx.end, rows, reach);
      spans.valid = stored;
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::blend_lit_times(const CoarseGrid &grid, const RenderBand &band,
                                const WarpField &now, const WarpField &prior) {
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

template <int W, int H>
inline HS_FLASH_MEMBER void
Feedback<W, H>::describe_rows(const FlushContext &ctx, RowGeometry *rows) {
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
  HS_CHECK(cover_base <= COVER_SLOTS, "feedback coverage needs %d slots of %d",
           cover_base, COVER_SLOTS);
}

HS_O3_END

} // namespace Filter::Pixel
