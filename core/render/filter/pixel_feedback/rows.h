/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file rows.h
 * @brief Feedback row compositing and lane resolution. */

#include "core/render/filter/pixel_feedback.h"
#include "core/render/filter/pixel_feedback/frame.h"

namespace Filter::Pixel {

HS_O3_BEGIN

template <int W, int H>
template <bool PAIR_PIXELS, typename TransformPixelT, typename TransformPairT>
inline void Feedback<W, H>::composite_pixels(const FlushContext &ctx,
                                             TransformPixelT &&transform_pixel,
                                             TransformPairT &&transform_pair,
                                             bool starts) {
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
            cell_level = static_cast<int>(
                FadeTable::UNIT_LEVEL + 0.5f -
                (offsets.left_x + 0.5f * offsets.slope_x) * LEVELS_PER_COLUMN);
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

template <int W, int H>
inline HS_FLASH_MEMBER void
Feedback<W, H>::composite_pole_row(const FlushContext &ctx, int y,
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::filter_infill_row(const FlushContext &ctx, int y,
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void Feedback<W, H>::composite_polar_row(
    const FlushContext &ctx, int xs, int xe, int y, const RowGeometry &g,
    ::Pixel *output, int stride, float lane_offset, bool defer_filter,
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

template <int W, int H>
inline HS_HOT_INLINE typename Feedback<W, H>::CapCell
Feedback<W, H>::polar_cell(const FlushContext &ctx, const RowGeometry &g,
                           int cx0, int cx1) {
  const CapOffset *caps0 = polar_caps(ctx.warp, ctx, g.field_y0);
  const CapOffset *caps1 = polar_caps(ctx.warp, ctx, g.field_y1);
  const WarpField &right = arm_edge(ctx, cx1) ? ctx.end : ctx.warp;
  return CapPlane::decode_cell(
      caps0[cx0], caps1[cx0], polar_caps(right, ctx, g.field_y0)[cx1],
      polar_caps(right, ctx, g.field_y1)[cx1], 1.0f - g.wy1, g.wy1);
}

template <int W, int H>
inline HS_HOT_INLINE typename Feedback<W, H>::Lane
Feedback<W, H>::cartesian_lane(int x, int y, int sub, float lane_offset,
                               float inverse_downsample,
                               const CellOffsets &cell) {
  const float fx = (sub + lane_offset) * inverse_downsample;
  const float target = x + lane_offset;
  const float dx = cell.left_x + cell.slope_x * fx;
  const float dy = cell.left_y + cell.slope_y * fx;
  return resolve_lane(target, ColumnTime::arm_start(x), target + dx, y + dy);
}

template <int W, int H>
inline HS_HOT_INLINE typename Feedback<W, H>::CellOffsets
Feedback<W, H>::cell_offsets(const FlushContext &ctx, int row0, int row1,
                             int cx0, int cx1, float wy0, float wy1,
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
  o.slope_x = (static_cast<float>(u10) * wy0 + static_cast<float>(u11) * wy1) *
                  INVERSE_WARP_SCALE -
              o.left_x;
  o.left_y = (ctx.warp.y_offsets[i00] * wy0 + ctx.warp.y_offsets[i01] * wy1) *
             INVERSE_WARP_SCALE;
  o.slope_y = (right.y_offsets[i10] * wy0 + right.y_offsets[i11] * wy1) *
                  INVERSE_WARP_SCALE -
              o.left_y;
  return o;
}

template <int W, int H>
inline HS_HOT_INLINE typename Feedback<W, H>::Lane
Feedback<W, H>::resolve_lane(float target, int arm, float source_x,
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

template <int W, int H>
inline HS_HOT_INLINE typename Feedback<W, H>::Sample
Feedback<W, H>::sample_lane(const FlushContext &ctx, int arm,
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::hold_lane(const FlushContext &ctx, int y, int arm, int x) {
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER int Feedback<W, H>::composite_resolved_cell(
    const FlushContext &ctx, PixelTransform transform_pixel, ::Pixel *output,
    int x, int end, int y, int sub, int stride, float lane_offset,
    const CellOffsets &cell, bool plain_store) {
  const float inverse_downsample = 1.0f / ctx.grid.downsample;
  const uint16_t alpha = frac_to_q16(ctx.alpha);
  while (x < end) {
    const int arm = ColumnTime::arm_start(x);
    Lane lane =
        cartesian_lane(x, y, sub, lane_offset, inverse_downsample, cell);
    ::Pixel p(0, 0, 0);
    if (lane.side == Side::AHEAD && ctx.cover) {
      hold_lane(ctx, y, arm, x);
      output[x] = p;
    } else {
      const Sample s = sample_lane(ctx, arm, lane);
      if (!skip_near_black(s, ctx.black_skips_color))
        p = transform_pixel(s.r, s.g, s.b, FadeTable::level(lane.read.elapsed));
      ::Pixel &dst = output[x];
      dst = plain_store ? p : dst.lerp16(p, alpha);
    }
    x += stride;
    sub += stride;
  }
  return x;
}

HS_O3_END

} // namespace Filter::Pixel
