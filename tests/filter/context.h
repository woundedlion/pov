/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/filter.h.
 */
#pragma once

/** @file
 * @brief Shared filter test fixtures. */
#include <functional>

#include "tests/vec_test_util.h"

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <type_traits>

#include "core/render/filter.h"
#include "core/render/filter/feedback_warp_cache.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/render/canvas.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/temporal_oracle.h"
#include <vector>

namespace hs_test {

struct FeedbackWhiteBox {
  /** @brief Full uncached composites fit the advertised device-resolution budget. */
  static void test_uncached_begin_frame_budget() {
    constexpr int W = 288, H = 144;
    using Feedback = Filter::Pixel::Feedback<W, H>;
    const ScopedTemporalSampling sampling(true);
    for (int downsample : {8, 16}) {
      const size_t bytes = Feedback::uncached_scratch_bytes(downsample) +
                           Feedback::FadeTable::BYTES;
      const ScopedArenaSplit split(GLOBAL_ARENA_SIZE - bytes -
                                       DEFAULT_SCRATCH_B_SIZE,
                                   bytes, DEFAULT_SCRATCH_B_SIZE);
      ScratchScope persistent_scope(persistent_arena);
      StubEffect effect(W, H, {.replaces_frame = true});
      ::Feedback::Style style{};
      style.noise = nullptr;
      style.fade = 1.0f;
      style.downsample = downsample;
      Pipeline<W, H, Feedback> pipeline{Feedback(style)};
      {
        Canvas canvas(effect);
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x)
            canvas(x, y) = Pixel(12000, 24000, 36000);
      }
      effect.advance_display();
      Canvas canvas(effect);
      HS_EXPECT_TRUE(canvas.keeps_stale_pixels());
      pipeline.begin_frame(canvas, 1.0f);
      HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(), bytes);
      HS_EXPECT_TRUE(canvas(8, 8).r > 0);
    }
  }

  /** @brief Uncached coarse fields fit their advertised scratch budget. */
  static void test_uncached_scratch_budget() {
    constexpr int W = 32, H = 32;
    using Feedback = Filter::Pixel::Feedback<W, H>;
    for (bool temporal : {false, true}) {
      const ScopedTemporalSampling sampling(temporal);
      ScratchScope persistent_scope(persistent_arena);
      StubEffect effect(W, H, {.replaces_frame = true});
      ::Feedback::Style style{};
      style.noise = nullptr;
      style.fade = 1.0f;
      style.downsample = W / 2;
      Feedback feedback(style);
      feedback.fade_table.init_storage(persistent_arena);
      const size_t bytes = Feedback::uncached_scratch_bytes(style.downsample);
      std::vector<uint8_t> storage(bytes);
      Arena scratch(storage.data(), storage.size());
      {
        Canvas canvas(effect);
        canvas(8, 8) = Pixel(12000, 24000, 36000);
      }
      effect.advance_display();
      Canvas canvas(effect);
      const auto ctx = feedback.prepare_flush(canvas, 1.0f, scratch);
      HS_EXPECT_TRUE(ctx.geometry.rows != nullptr);
      if (temporal)
        HS_EXPECT_TRUE(ctx.frames.reach != nullptr);
      HS_EXPECT_LE(scratch.get_high_water_mark(), bytes);
    }
  }

  /** @brief Cartesian composite and finish agree at the arm cut. */
  static void test_cartesian_cut_roundoff() {
    constexpr int W = 288, H = 16, X = 143, Y = 8, DOWNSAMPLE = 4;
    using Feedback = Filter::Pixel::Feedback<W, H>;
    ::Feedback::Style style{};
    Feedback feedback(style);
    typename Feedback::FlushContext ctx{
        {{DOWNSAMPLE,
          typename Feedback::SphereField{DOWNSAMPLE},
          W / DOWNSAMPLE,
          0,
          {}}}};
    std::array<typename Feedback::RowGeometry, H> rows{};
    std::array<typename Feedback::AheadSpan, H> spans{};
    std::array<uint16_t, W * H> cover{};
    std::array<Pixel, W> output{};
    rows[Y].wide = true;
    ctx.geometry.rows = rows.data();
    ctx.coverage.ahead_span = spans.data();
    ctx.coverage.cover = cover.data();
    ctx.color.alpha = 1.0f;
    const typename Feedback::CellOffsets offsets{0.2f, -0.200005f / 0.75f, 0.0f,
                                                 0.0f};
    feedback.composite_resolved_cell(ctx, nullptr, output.data(), X, X + 1, Y,
                                     X % DOWNSAMPLE, 1, 0.0f, offsets, true);
    const int slot = feedback.cover_slot(rows[Y], 0, X);
    HS_EXPECT_TRUE((cover[slot] & Feedback::AHEAD_BIT) != 0);
    typename Feedback::SweepRow sweep{};
    sweep.cell_cx[0] = X / DOWNSAMPLE;
    sweep.cell[0].offsets = offsets;
    const auto lane = feedback.sweep_lane(ctx, rows[Y], Y, X, sweep);
    HS_EXPECT_TRUE(lane.side == Feedback::Side::AHEAD);
    HS_EXPECT_EQ(lane.read.x, static_cast<float>(X));
  }

  template <int W, int H>
  static bool holds(const Filter::Pixel::Feedback<W, H> &feedback,
                    const Canvas &canvas) {
    const auto grid = feedback.make_coarse_grid(canvas);
    const auto band = feedback.make_render_band(canvas.clip(), grid);
    return feedback.warp_cache.holds(feedback.warp_cache_key(band));
  }
  template <int W, int H>
  static size_t fills(const Filter::Pixel::Feedback<W, H> &feedback) {
    return feedback.warp_cache.cached_fills;
  }
  template <int W, int H>
  static void expand_clipped_row(Pixel *output, int begin, int end) {
    using Feedback = Filter::Pixel::Feedback<W, H>;
    typename Feedback::FlushContext ctx{
        {{1, typename Feedback::SphereField{1}, 0, 0, {}}}};
    ctx.color.column_time = true;
    ctx.geometry.runs.items[0] = {begin, end};
    ctx.geometry.runs.count = 1;
    Feedback::expand_half_res_row(ctx, output, 0);
  }
  /** @brief The coverage slot of a held pixel, or -1 for an unheld one. */
  template <int W, int H>
  static int held_slot(const Filter::Pixel::Feedback<W, H> &feedback, int y,
                       int x) {
    if (!feedback.split_hold.held)
      return -1;
    const int slot = feedback.cover_slot(
        feedback.split_hold.held->geometry.rows[y], x < W / 2 ? 0 : W / 2, x);
    return slot >= 0 && (feedback.split_hold.held->coverage.cover[slot] &
                         0xFFFE) != 0
               ? slot
               : -1;
  }
  /** @brief Held pixels a split frame's taps reached. */
  template <int W, int H>
  static int covered(const Filter::Pixel::Feedback<W, H> &feedback) {
    if (!feedback.split_hold.held)
      return 0;
    int count = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const int slot = held_slot(feedback, y, x);
        count += slot >= 0 && (feedback.split_hold.held->coverage.cover[slot] &
                               0xFFFE) != 0xFFFE;
      }
    return count;
  }
  /** @brief Replaces the part of each held pixel's composite the previous
   *  frames supplied with a given frame beneath its coverage, and releases
   *  the frame.
   *  @param under Frame beneath the held pixels, row-major.
   *  @param known The frame as flush_split() left it, row-major. */
  template <int W, int H>
  static void release(Filter::Pixel::Feedback<W, H> &feedback,
                      const Pixel *under, const Pixel *known) {
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const int slot = held_slot(feedback, y, x);
        if (slot < 0)
          continue;
        const uint16_t t =
            feedback.split_hold.held->coverage.cover[slot] & 0xFFFE;
        auto swap = [t](uint16_t covered, uint16_t was, uint16_t now) {
          const int32_t v = int32_t(covered) - int32_t(lerp_q16(0, was, t)) +
                            int32_t(lerp_q16(0, now, t));
          return static_cast<uint16_t>(std::min(std::max(v, 0), 65535));
        };
        Pixel &dst = feedback.split_hold.held->frames.current[y * W + x];
        const Pixel &k = known[y * W + x], &u = under[y * W + x];
        dst = Pixel(swap(dst.r, k.r, u.r), swap(dst.g, k.g, u.g),
                    swap(dst.b, k.b, u.b));
      }
    feedback.split_hold.release();
  }
};

} // namespace hs_test

namespace hs_test::filter_tests {

template <typename P>
concept RawFramePlotter = requires(P &pipeline, Canvas &canvas) {
  pipeline.plot(canvas, 0, 0, Pixel{}, 0.0f, 1.0f);
};

template <typename P>
concept RawTimedPlotter = requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_at(canvas, math::X_AXIS, Pixel{}, 0.0f, 1.0f, 0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_transformed_at(canvas, math::X_AXIS, Pixel{}, 0.0f, 1.0f, 0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_projected_at(canvas, math::PixelCoords{}, Pixel{}, 0.0f, 1.0f,
                             0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_column_at(canvas, math::PixelCoords{}, 0, Pixel{}, 0.0f, 1.0f,
                          0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_screen_at(canvas, 0.0f, 0.0f, Pixel{}, 0.0f, 1.0f, 0.5f);
};

template <typename P>
concept TerminalFlusher =
    requires(P &pipeline, Canvas &canvas) { pipeline.flush(canvas, 1.0f); };

template <typename P>
concept ReplacementFrameStarter = requires(P &pipeline, Canvas &canvas) {
  pipeline.begin_frame(canvas, 1.0f);
};
} // namespace hs_test::filter_tests
