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

  HS_EXPECT_NEAR(fb.style().fade, 0.9f, 1e-6f);
  HS_EXPECT_EQ(fb.style().downsample, 4);

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
 *        only after the fill returns.
 */
inline void test_feedback_warp_cache_valid_only_after_fill() {
  constexpr int W = 32, H = 16;
  constexpr int DS = 4;
  using Cache = Filter::Pixel::FeedbackWarpCache<W, H>;
  constexpr hs::SphericalFieldLayout<W, H> field(DS, DS, DS, W / DS);
  constexpr int CELLS = (W / DS) * field.ring_count();

  ScratchScope persistent_scope(persistent_arena);
  Cache cache;
  HS_EXPECT_FALSE(cache.ready());
  cache.init_storage<CELLS>(
      persistent_arena, field, W / DS,
      [](const math::Vector &, const Cache::Coordinates &point) {
        return point;
      });
  HS_EXPECT_TRUE(cache.ready());

  const Cache::Key first{
      &::Feedback::noise_warp, nullptr, 0u, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0,
      field.ring_count() - 1};
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

  held_first_during_fill = true;
  cache.acquire(&second, Cache::Buffers{}, fill);
  HS_EXPECT_EQ(fills, 2);
  HS_EXPECT_FALSE(held_first_during_fill);
  HS_EXPECT_FALSE(held_second_during_fill);
  HS_EXPECT_TRUE(cache.holds(second));
  HS_EXPECT_FALSE(cache.holds(first));
  HS_EXPECT_EQ(a.x_offsets[0], 2);

  int16_t scratch_x[1] = {0};
  const Cache::Buffers scratch{scratch_x, nullptr, nullptr};
  const Cache::Buffers c = cache.acquire(nullptr, scratch, fill);
  HS_EXPECT_EQ(fills, 3);
  HS_EXPECT_TRUE(held_second_during_fill);
  HS_EXPECT_TRUE(c.x_offsets == scratch_x);
  HS_EXPECT_EQ(scratch_x[0], 3);
  HS_EXPECT_EQ(a.x_offsets[0], 2);
  HS_EXPECT_TRUE(cache.holds(second));
}
