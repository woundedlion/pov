/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sampling.h
 * @brief Feedback frame sampling and held pixel taps. */

#include "core/render/filter/pixel_feedback.h"
#include "core/render/filter/pixel_feedback/reconstruction.h"

namespace Filter::Pixel {

template <int W, int H>
inline bool Feedback<W, H>::any_pixel_lit(const Canvas &cv) {
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

template <int W, int H>
inline HS_O3_FN ::Pixel
Feedback<W, H>::select_pole_sample(const ::Pixel *pole_row) {
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

template <int W, int H>
inline void Feedback<W, H>::select_poles(
    const ::Pixel *frame, ::Pixel (&poles)[SphereField::POLE_STORAGE_COUNT]) {
  if (SphereField::HAS_NORTH_POLE)
    poles[0] = select_pole_sample(frame);
  if (SphereField::HAS_SOUTH_POLE)
    poles[SphereField::HAS_NORTH_POLE ? 1 : 0] =
        select_pole_sample(frame + (H - 1) * W);
}

template <int W, int H>
inline HS_HOT_FLASH_MEMBER typename Feedback<W, H>::Sample
Feedback<W, H>::sample_across(const FlushContext &ctx, int arm, float x,
                              float y) const {
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
  const Sample a = sample_bilinear_prev(field, low, low_poles, first - 1.0f, y);
  const Sample b = sample_bilinear_prev(field, high, high_poles, first, y);
  const float t = x - (first - 1.0f);
  return {a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t};
}

template <int W, int H>
inline HS_O3_FN __attribute__((noinline)) typename Feedback<W, H>::Sample
Feedback<W, H>::sample_bilinear_prev(
    const SphereField &field, const ::Pixel *prev,
    const ::Pixel (&poles)[SphereField::POLE_STORAGE_COUNT], float bx,
    float by) const {
  Sample out;
  field.sample_bilinear_rgb(prev, poles, bx, by, out.r, out.g, out.b);
  return out;
}

template <int W, int H>
template <typename PassFnT>
inline HS_HOT_FLASH_MEMBER void
Feedback<W, H>::plot_near_arm_end(float x, int column, float y,
                                  const ::Pixel &color, float age, float alpha,
                                  PassFnT &pass) {
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

} // namespace Filter::Pixel
