/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file reconstruction.h
 * @brief Feedback reach bounds and row reconstruction. */

#include "core/render/filter/pixel_feedback.h"
#include "core/render/filter/pixel_feedback/sweeps.h"

namespace Filter::Pixel {

template <int W, int H>
inline typename Feedback<W, H>::ColumnRuns
Feedback<W, H>::row_lanes(const FlushContext &ctx, const RowGeometry &g, int y,
                          bool starts) {
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

template <int W, int H>
inline HS_FLASH_MEMBER void
Feedback<W, H>::measure_reach(const CoarseGrid &grid, const WarpField &warp,
                              const WarpField &end, const RowGeometry *rows,
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
          [&](int i) { return base + (i >> 1) * columns + (i & 1 ? cx1 : cx); },
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

template <int W, int H>
inline __attribute__((always_inline)) bool
Feedback<W, H>::filtered_infill(const CoarseGrid &grid, int y) {
  const int downsample = grid.downsample;
  const bool infill_band =
      (y < downsample && (!SphereField::HAS_NORTH_POLE || y > 0)) ||
      (y >= H - downsample && (!SphereField::HAS_SOUTH_POLE || y < H - 1));
  return infill_band && grid.field.longitude_filter_width(y) > 1;
}

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::reconstruct_half_res_run(::Pixel *output, int begin, int end,
                                         const ::Pixel *before,
                                         const ::Pixel *after) {
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

template <int W, int H>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::expand_half_res_row(const FlushContext &ctx, ::Pixel *output,
                                    int y) {
  const ColumnRuns &runs = ctx.runs;
  if (!ctx.column_time || runs.count != 2 || runs.items[0].begin != 0 ||
      runs.items[1].end != W) {
    for (int r = 0; r < runs.count; ++r) {
      const ColumnRun run = runs.items[r];
      if (!ctx.column_time) {
        reconstruct_half_res_run(output, run.begin, run.end);
        continue;
      }
      const int arm = ColumnTime::arm_start(run.begin);
      const int end = std::min(run.end, arm + HALF - 2);
      if (run.begin >= end)
        continue;
      if (run.end == arm + HALF) {
        const ::Pixel after =
            pair_average(output[arm + HALF - 2], output[arm + HALF - 1]);
        reconstruct_half_res_run(output, run.begin, end, nullptr, &after);
      } else {
        reconstruct_half_res_run(output, run.begin, end);
      }
    }
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

template <int W, int H>
inline __attribute__((always_inline)) int
Feedback<W, H>::unwrap_near(int v, int reference) {
  constexpr int PERIOD = static_cast<int>(WRAP_PERIOD);
  constexpr int HALF_PERIOD = PERIOD / 2;
  const int delta = v - reference;
  return delta > HALF_PERIOD ? v - PERIOD
                             : (delta < -HALF_PERIOD ? v + PERIOD : v);
}

} // namespace Filter::Pixel
