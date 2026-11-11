/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Pixel::Feedback — Style binding and plot passthrough (no flush / no Canvas)
// ============================================================================

/**
 * @brief Verifies style() returns a live reference to the bound Style: reads,
 *        mutation, and the const overload all alias the original.
 */
inline void test_feedback_style_binding() {
  constexpr int W = 32, H = 32;
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  Filter::Pixel::Feedback<W, H> fb(style);

  HS_EXPECT_TRUE(&fb.style() == &style);

  const Filter::Pixel::Feedback<W, H> &cfb = fb;
  HS_EXPECT_TRUE(&cfb.style() == &style);

  HS_EXPECT_NEAR(fb.style().fade, ::Feedback::Style::Smoke().fade, 1e-6f);
  HS_EXPECT_EQ(fb.style().downsample, ::Feedback::Style::Smoke().downsample);

  // Mutating through the accessor is visible on the original Style.
  fb.style().fade = 0.5f;
  HS_EXPECT_NEAR(style.fade, 0.5f, 1e-6f);
}

/**
 * @brief Verifies Feedback::plot() forwards its input unchanged: one tap with
 *        the original coord, colour, age and alpha, no Canvas touched.
 * @details The feedback effect lives in flush, not plot.
 */
inline void test_feedback_plot_is_passthrough() {
  constexpr int W = 32, H = 32;
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  Filter::Pixel::Feedback<W, H> fb(style);

  int count = 0;
  float kx = -1, ky = -1, kage = -1, ka = -1;
  Pixel kc(0, 0, 0);
  Pixel src(7, 8, 9);
  fb.plot(3.0f, 4.0f, src, 1.5f, 0.75f,
          [&](float x, float y, const Pixel &c, float age, float a) {
            ++count;
            kx = x;
            ky = y;
            kc = c;
            kage = age;
            ka = a;
          });
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_NEAR(kx, 3.0f, 1e-6f);
  HS_EXPECT_NEAR(ky, 4.0f, 1e-6f);
  HS_EXPECT_PIXEL(kc, 7, 8, 9);
  HS_EXPECT_NEAR(kage, 1.5f, 1e-6f);
  HS_EXPECT_NEAR(ka, 0.75f, 1e-6f);
}

/**
 * @brief Verifies FeedbackWarpCache stores the key and marks itself valid
 *        only after the fill returns, keeping the field before a later noise
 *        time as its prior.
 */
inline void test_feedback_warp_cache_valid_only_after_fill() {
  constexpr int W = 32, H = 16;
  constexpr int DS = 4;
  using Cache = Filter::Pixel::FeedbackWarpCache<W, H>;
  constexpr hs::SphericalFieldLayout<W, H> FIELD(DS, DS, DS, W / DS);
  constexpr int CELLS = (W / DS) * FIELD.ring_count();

  ScratchScope persistent_scope(persistent_arena);
  Cache cache;
  HS_EXPECT_FALSE(cache.ready());
  cache.init_storage<CELLS>(
      persistent_arena, FIELD, W / DS,
      [](const math::Vector &, const Cache::Coordinates &point) {
        return point;
      });
  HS_EXPECT_TRUE(cache.ready());

  const Cache::Key first{
      &::Feedback::noise_warp, nullptr, 0u, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0,
      FIELD.ring_count() - 1,  true};
  Cache::Key second = first;
  second.time = 1.0f;
  HS_EXPECT_FALSE(cache.holds(first));

  int fills = 0;
  bool held_first_during_fill = true;
  bool held_second_during_fill = true;
  auto fill = [&](const Cache::Buffers &buffers) {
    ++fills;
    held_first_during_fill = cache.holds(first);
    held_second_during_fill = cache.holds(second);
    buffers.x_offsets[0] = static_cast<int16_t>(fills);
  };

  const Cache::Buffers a = cache.acquire(&first, Cache::Buffers{}, fill);
  HS_EXPECT_EQ(fills, 1);
  HS_EXPECT_FALSE(held_first_during_fill);
  HS_EXPECT_TRUE(cache.holds(first));
  HS_EXPECT_EQ(a.x_offsets[0], 1);

  const Cache::Buffers b = cache.acquire(&first, Cache::Buffers{}, fill);
  HS_EXPECT_EQ(fills, 1);
  HS_EXPECT_TRUE(b.x_offsets == a.x_offsets);

  Cache::Buffers prior{};
  held_first_during_fill = true;
  const Cache::Buffers later =
      cache.acquire(&second, Cache::Buffers{}, fill, &prior);
  HS_EXPECT_EQ(fills, 2);
  HS_EXPECT_FALSE(held_first_during_fill);
  HS_EXPECT_FALSE(held_second_during_fill);
  HS_EXPECT_TRUE(cache.holds(second));
  HS_EXPECT_FALSE(cache.holds(first));
  // A later noise time fills the other set and keeps the first as its prior.
  HS_EXPECT_TRUE(prior.x_offsets == a.x_offsets);
  HS_EXPECT_FALSE(later.x_offsets == a.x_offsets);
  HS_EXPECT_EQ(a.x_offsets[0], 1);
  HS_EXPECT_EQ(later.x_offsets[0], 2);
  cache.acquire(&second, Cache::Buffers{}, fill, &prior);
  HS_EXPECT_TRUE(prior.x_offsets == later.x_offsets);
  Cache::Key moved = second;
  moved.amplitude = 2.0f;
  const Cache::Buffers refilled =
      cache.acquire(&moved, Cache::Buffers{}, fill, &prior);
  HS_EXPECT_TRUE(prior.x_offsets == refilled.x_offsets);
  HS_EXPECT_TRUE(refilled.x_offsets == later.x_offsets);

  int16_t scratch_x[1] = {0};
  const Cache::Buffers scratch{scratch_x, nullptr, nullptr};
  const Cache::Buffers c = cache.acquire(nullptr, scratch, fill);
  HS_EXPECT_EQ(fills, 4);
  HS_EXPECT_TRUE(c.x_offsets == scratch_x);
  HS_EXPECT_EQ(scratch_x[0], 4);
  HS_EXPECT_EQ(refilled.x_offsets[0], 3);
  HS_EXPECT_TRUE(cache.holds(moved));
}

/** @brief Persistent feedback allocation follows lattice samples rather than rectangular cells. */
inline void test_feedback_cache_storage_matches_lattice_samples() {
  reset_globals();
  ScopedArenaSplit split(65536, 131072, 131072);
  using Feedback = Filter::Pixel::Feedback<288, 144>;
  using Layout = hs::SphericalFieldLayout<288, 144>;
  constexpr int DOWNSAMPLE = ::Feedback::Style{}.downsample;
  constexpr Layout FIELD(DOWNSAMPLE, DOWNSAMPLE, DOWNSAMPLE, 288 / DOWNSAMPLE);
  constexpr int CELLS = (288 / DOWNSAMPLE) * FIELD.ring_count();
  HS_EXPECT_LT(FIELD.sample_count(), CELLS);
  HS_EXPECT_EQ(Feedback::CACHE_ORIGIN_SAMPLES, FIELD.sample_count());
  ::Feedback::Style style{};
  Feedback feedback(style);
  const size_t before = persistent_arena.get_offset();
  feedback.init_storage(persistent_arena);
  HS_EXPECT_EQ(persistent_arena.get_offset() - before, Feedback::STORAGE_BYTES);
  std::printf("  feedback origins: cells=%d samples=%d saved=%zu bytes\n",
              CELLS, FIELD.sample_count(),
              (CELLS - FIELD.sample_count()) * sizeof(Layout::Coordinates));
}
