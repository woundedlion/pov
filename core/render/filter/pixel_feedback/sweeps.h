/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sweeps.h
 * @brief Feedback temporal sweeps and composite dispatch. */

#include "core/render/filter/pixel_feedback.h"
#include "core/render/filter/pixel_feedback/rows.h"

namespace Filter::Pixel {

HS_O3_BEGIN

template <int W, int H>
inline typename Feedback<W, H>::SweepRow *
Feedback<W, H>::sweep_rows(const FlushContext &ctx, Arena &arena) {
  SweepRow *rows = arena.allocate_n<SweepRow>(H);
  for (int y = ctx.band.y_begin; y < ctx.band.y_end; ++y) {
    const RowGeometry &g = ctx.rows[y];
    rows[y].cap_angle = g.polar ? cap_angle_of(y, g.north) : 0.0f;
    rows[y].cell_cx[0] = rows[y].cell_cx[1] = -1;
  }
  return rows;
}

template <int W, int H>
inline HS_HOT_INLINE typename Feedback<W, H>::Lane
Feedback<W, H>::sweep_lane(const FlushContext &ctx, const RowGeometry &g, int y,
                           int x, SweepRow &row) const {
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
  const float inverse_downsample = 1.0f / downsample;
  const float fx = (x - cx0 * downsample + lane_offset) * inverse_downsample;
  const float target = x + lane_offset;
  const int arm = ColumnTime::arm_start(x);
  if (g.polar) {
    const auto at =
        CapPlane::polar_lane(cell.cap, x, row.cap_angle, g.north, pair, fx);
    return resolve_lane(target, arm, at.x, at.y);
  }
  return cartesian_lane(x, y, x - cx0 * downsample, lane_offset,
                        inverse_downsample, cell.offsets);
}

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::sweep_deep(const FlushContext &ctx,
                           PixelTransform transform_pixel, int begin, int end) {
  ScratchScope scope(scratch_arena_a);
  SweepRow *rows = sweep_rows(ctx, scope.get_arena());
  for (int c = begin; c < end; ++c)
    for (int y = ctx.band.y_begin; y < ctx.band.y_end; ++y)
      if (ctx.reach[y] > HALF / 2 && c < ctx.reach[y])
        behind_lane(ctx, transform_pixel, y, c, rows[y]);
}

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::behind_lane(const FlushContext &ctx,
                            PixelTransform transform_pixel, int y, int c,
                            SweepRow &row) {
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::sweep_ahead(const FlushContext &ctx,
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
          // no death case: held sources stay at the cut within float roundoff
          HS_CHECK(lane.side == Side::AHEAD ||
                       lane.read.x >= arm + HALF - 1.0f - 1.0e-3f,
                   "feedback held lane %d of row %d no longer reads ahead", x,
                   y);
          // A column is final once every lane its reconstruction reads is
          // composited: the read stops at the last such column of its
          // source rows.
          const int row0 =
              hs::clamp(static_cast<int>(std::floor(lane.read.y)), 0, H - 1);
          const int reach = std::max<int>(
              ctx.rows[row0].reach, ctx.rows[std::min(row0 + 1, H - 1)].reach);
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::complete_columns(const FlushContext &ctx, int y,
                                 const RowGeometry &g, int arm, int c,
                                 const ::Pixel &sample, ::Pixel *window) {
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

template <int W, int H>
inline HS_HOT_INLINE void Feedback<W, H>::add_under(const FlushContext &ctx,
                                                    int y, int arm, int x,
                                                    const ::Pixel &under) {
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
  dst = ::Pixel(add(dst.r, under.r), add(dst.g, under.g), add(dst.b, under.b));
}

template <int W, int H>
inline HS_HOT_INLINE void Feedback<W, H>::composite_frame(Canvas &cv,
                                                          FlushContext &ctx) {
  if (cv.keeps_stale_pixels() && !ctx.back)
    cv.clear_buffer();
  if (ctx.mode == FlushMode::SKIP)
    return;
  HS_PROFILE(feedback_composite);
  composite_pass(ctx, false);
}

template <int W, int H>
template <bool PAIR_PIXELS, typename TransformPixelT, typename TransformPairT>
inline HS_HOT_INLINE void
Feedback<W, H>::composite_phase(const FlushContext &ctx, bool finish,
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

template <int W, int H>
inline HS_FLASH_MEMBER __attribute__((flatten)) void
Feedback<W, H>::composite_plain(const FlushContext &ctx, bool finish) {
  const float *fades = ctx.column_time ? ctx.fades.fades : &ctx.fade;
  auto plain = [&](float r, float g, float b, int level) {
    const float fade = fades[level];
    return ::Pixel(quantize16(r * fade), quantize16(g * fade),
                   quantize16(b * fade));
  };
  composite_phase<false>(ctx, finish, plain, plain);
}

template <int W, int H>
inline HS_FLASH_MEMBER __attribute__((flatten)) void
Feedback<W, H>::composite_general(const FlushContext &ctx, bool finish) {
  const float *fades = ctx.column_time ? ctx.fades.fades : &ctx.fade;
  auto general = [&](float r, float g, float b, int level) {
    return feedback_style->color_fn(
        ::Pixel(quantize16(r), quantize16(g), quantize16(b)), fades[level],
        *feedback_style);
  };
  composite_phase<false>(ctx, finish, general, general);
}

HS_O3_END

} // namespace Filter::Pixel
