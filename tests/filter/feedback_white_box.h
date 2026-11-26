/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief White-box access to Filter::Pixel::Feedback for filter tests. */

#include "tests/filter/context.h"

namespace hs_test {

/** @brief Test access to Filter::Pixel::Feedback internals, and the
 *         Feedback cases that need it. */
struct FeedbackWhiteBox {
  /** @brief Full uncached composites fit the advertised device-resolution budget. */
  static void test_uncached_begin_frame_budget() {
    constexpr int W = 288, H = 144;
    using Feedback = Filter::Pixel::Feedback<W, H>;
    const ScopedTemporalSampling sampling(true);
    for (int downsample : {8, 16}) {
      const size_t bytes = Feedback::uncached_scratch_bytes(downsample) +
                           Feedback::FadeTable::BYTES;
      constexpr size_t ALIGNMENT = alignof(std::max_align_t);
      const size_t PERSISTENT_BYTES =
          (GLOBAL_ARENA_SIZE - bytes - DEFAULT_SCRATCH_B_SIZE) / ALIGNMENT *
          ALIGNMENT;
      const ScopedArenaSplit split(PERSISTENT_BYTES, bytes,
                                   DEFAULT_SCRATCH_B_SIZE);
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

  /** @brief Column-time blending caps forward offsets after seam wrapping. */
  static void test_blended_warp_caps_forward_offsets() {
    constexpr int W = 288, H = 32, DOWNSAMPLE = 8, COLUMNS = W / DOWNSAMPLE;
    using Feedback = Filter::Pixel::Feedback<W, H>;
    constexpr int16_t FORWARD = Feedback::WRAP_PERIOD * 0.25f;
    constexpr int16_t BACKWARD = -Feedback::WRAP_PERIOD * 0.5f;
    const Feedback::CoarseGrid grid{DOWNSAMPLE,
                                    Feedback::SphereField{DOWNSAMPLE},
                                    COLUMNS,
                                    1,
                                    {0, 0, 1, 0, 0, 1}};
    const Feedback::RenderBand band{};
    std::array<int16_t, COLUMNS> prior_x, now_x, prior_y{}, now_y{};
    prior_x.fill(FORWARD);
    now_x.fill(BACKWARD);
    Feedback::blend_lit_times(grid, band, {now_x.data(), now_y.data(), nullptr},
                              {prior_x.data(), prior_y.data(), nullptr});
    for (int16_t x : prior_x)
      HS_EXPECT_EQ(x, FORWARD);
  }

  /** @brief Whether the warp cache holds the canvas's current band.
   *  @tparam W Feedback width.
   *  @tparam H Feedback height.
   *  @param feedback Filter under test.
   *  @param canvas Canvas whose clip selects the band.
   *  @return True on a cache hit. */
  template <int W, int H>
  static bool holds(const Filter::Pixel::Feedback<W, H> &feedback,
                    const Canvas &canvas) {
    const auto grid = feedback.make_coarse_grid(canvas);
    const auto band = feedback.make_render_band(canvas.clip(), grid);
    return feedback.warp_cache.holds(feedback.warp_cache_key(band));
  }
  /** @brief Warp cache fill count.
   *  @tparam W Feedback width.
   *  @tparam H Feedback height.
   *  @param feedback Filter under test.
   *  @return Cache population calls observed so far. */
  template <int W, int H>
  static size_t fills(const Filter::Pixel::Feedback<W, H> &feedback) {
    return feedback.warp_cache.cached_fills;
  }
  /** @brief Expands a half-resolution row in place over a clipped column span.
   *  @tparam W Feedback width.
   *  @tparam H Feedback height.
   *  @param output Row of @p W pixels expanded in place.
   *  @param begin First rendered column.
   *  @param end Column one past the last rendered one; less than @p begin for
   *         a span that wraps the seam. */
  template <int W, int H>
  static void expand_clipped_row(Pixel *output, int begin, int end) {
    using Feedback = Filter::Pixel::Feedback<W, H>;
    typename Feedback::FlushContext ctx{
        {{1, typename Feedback::SphereField{1}, 0, 0, {}}}};
    ctx.color.column_time = true;
    if (begin > end) {
      ctx.geometry.runs.items[0] = {0, end};
      ctx.geometry.runs.items[1] = {begin, W};
      ctx.geometry.runs.count = 2;
    } else {
      ctx.geometry.runs.items[0] = {begin, end};
      ctx.geometry.runs.count = 1;
    }
    ctx.frames.previous = output;
    Feedback::expand_half_res_row(ctx, output, 0);
  }
  /** @brief The coverage slot of a held pixel.
   *  @tparam W Feedback width.
   *  @tparam H Feedback height.
   *  @param feedback Filter under test.
   *  @param y Pixel row.
   *  @param x Pixel column.
   *  @return The slot, or -1 for an unheld pixel. */
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
  /** @brief Held pixels a split frame's taps reached.
   *  @tparam W Feedback width.
   *  @tparam H Feedback height.
   *  @param feedback Filter under test.
   *  @return Count of reached held pixels. */
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
   *  @tparam W Feedback width.
   *  @tparam H Feedback height.
   *  @param feedback Filter holding the split frame.
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
