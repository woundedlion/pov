/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// ============================================================================
// ColorWipe
// ----------------------------------------------------------------------------
// OKLCH-lerps between caller-owned snapshots. make_palette builds from three
// fixed keys via PaletteRecipes::from_colors without RNG draws.
// ============================================================================

/**
 * @brief Builds a deterministic palette with no RNG draws.
 * @param ka First key color.
 * @param kb Second key color.
 * @param kc Third key color.
 * @return A STRAIGHT-gradient palette over the supplied keys.
 */
inline GenerativePalette make_palette(CPixel ka, CPixel kb, CPixel kc) {
  return GenerativePalette(
      PaletteRecipes::from_colors(PaletteDomain::STRAIGHT, ka, kb, kc));
}

/**
 * @brief Verifies ColorWipe lands the source palette's keys exactly on the
 * target keys at completion and reports done() only on the final frame.
 */
inline void test_colorwipe_reaches_target_keys() {
  GenerativePalette from =
      make_palette(CPixel(10, 20, 30), CPixel(40, 50, 60), CPixel(70, 80, 90));
  GenerativePalette to =
      make_palette(CPixel(200, 0, 0), CPixel(0, 200, 0), CPixel(0, 0, 200));
  const GenerativePalette::Snapshot start = from.snapshot();
  GenerativePalette::Snapshot target = to.snapshot();

  const int duration = 6;
  Animation::ColorWipe wipe(from, start, target, duration, math::ease_linear);
  HS_EXPECT_FALSE(wipe.done());

  for (int i = 0; i < duration - 1; ++i)
    wipe.step(fake_canvas());
  HS_EXPECT_FALSE(wipe.done()); // not done until t == duration

  wipe.step(fake_canvas()); // t == duration: amount == 1 -> exact target keys
  HS_EXPECT_TRUE(wipe.done());
  const auto actual = from.snapshot();
  HS_EXPECT_EQ(actual.key_count, target.key_count);
  for (int i = 0; i < target.key_count; ++i) {
    const auto got = GenerativePalette::snapshot_key(actual, i);
    const auto want = GenerativePalette::snapshot_key(target, i);
    HS_EXPECT_EQ(got.L, want.L);
    HS_EXPECT_EQ(got.chroma, want.chroma);
    HS_EXPECT_EQ(got.h, want.h);
  }
}

/**
 * @brief Verifies the wipe reads its explicit start snapshot rather than the
 * subject's current value.
 */
inline void test_colorwipe_uses_owned_start_snapshot() {
  GenerativePalette from =
      make_palette(CPixel(10, 10, 10), CPixel(10, 10, 10), CPixel(10, 10, 10));
  GenerativePalette to = make_palette(
      CPixel(250, 250, 250), CPixel(250, 250, 250), CPixel(250, 250, 250));
  const GenerativePalette::Snapshot start = from.snapshot();
  const GenerativePalette::Snapshot target = to.snapshot();
  const int duration = 4;
  Animation::ColorWipe wipe(from, start, target, duration, math::ease_linear);

  from = make_palette(CPixel(200, 200, 200), CPixel(200, 200, 200),
                      CPixel(200, 200, 200));

  wipe.step(fake_canvas());
  const float modified_l = pixel_to_oklch(Pixel(CPixel(200, 200, 200))).L;
  HS_EXPECT_LT(GenerativePalette::snapshot_key(from.snapshot(), 0).L,
               modified_l);
}

/**
 * @brief Verifies a slow wipe resolves a new key level on essentially every
 *        frame.
 * @details The `advanced` assertion bounds distinct key updates during the fade.
 */
inline void test_colorwipe_slow_fade_resolves_every_frame() {
  GenerativePalette from =
      make_palette(CPixel(10, 10, 10), CPixel(10, 10, 10), CPixel(10, 10, 10));
  GenerativePalette to = make_palette(
      CPixel(250, 250, 250), CPixel(250, 250, 250), CPixel(250, 250, 250));
  const GenerativePalette::Snapshot start = from.snapshot();
  const GenerativePalette::Snapshot target = to.snapshot();
  const int duration = 600;
  Animation::ColorWipe wipe(from, start, target, duration, math::ease_linear);

  int advanced = 0;
  float prev = -1.0f;
  for (int i = 0; i < duration; ++i) {
    wipe.step(fake_canvas());
    float cur = GenerativePalette::snapshot_key(from.snapshot(), 0).L;
    if (i > 0 && cur != prev)
      ++advanced;
    prev = cur;
  }
  HS_EXPECT_GT(advanced, 400);
}

/**
 * @brief Verifies a wired pause flag freezes a ColorWipe: the source keys hold
 * and the timer does not advance while paused.
 */
inline void test_colorwipe_paused_holds_keys() {
  bool paused = true;
  GenerativePalette from =
      make_palette(CPixel(10, 10, 10), CPixel(10, 10, 10), CPixel(10, 10, 10));
  GenerativePalette to = make_palette(
      CPixel(250, 250, 250), CPixel(250, 250, 250), CPixel(250, 250, 250));
  const GenerativePalette::Snapshot start = from.snapshot();
  GenerativePalette::Snapshot target = to.snapshot();
  const float start_l = GenerativePalette::snapshot_key(from.snapshot(), 0).L;

  const int duration = 4;
  Animation::ColorWipe wipe(from, start, target, duration, math::ease_linear,
                            {.paused = &paused});
  for (int i = 0; i < duration; ++i)
    wipe.step(fake_canvas());
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(from.snapshot(), 0).L, start_l,
                 1e-6f);
  HS_EXPECT_FALSE(wipe.done());

  paused = false;
  for (int i = 0; i < duration; ++i)
    wipe.step(fake_canvas());
  HS_EXPECT_TRUE(wipe.done());
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(from.snapshot(), 0).L,
                 GenerativePalette::snapshot_key(target, 0).L, 1e-6f);
}
