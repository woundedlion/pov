/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// PaletteCycler
// ============================================================================

inline void test_palette_cycler_arena_components() {
  HS_EXPECT_EQ(PaletteCycler::display_arena_bytes(),
               BakedPalette::required_arena_bytes());
  HS_EXPECT_EQ(PaletteCycler::crossfade_arena_bytes(),
               2 * BakedPalette::required_arena_bytes());
  HS_EXPECT_EQ(PaletteCycler::morph_arena_bytes(),
               sizeof(GenerativePalette) + alignof(GenerativePalette));
}

inline void test_palette_cycler_key_morph_cycle() {
  const GenerativePalette first(PaletteRecipes::balanced_analogous(0.05f));
  const GenerativePalette second(PaletteRecipes::balanced_analogous(0.62f));
  const std::array<PaletteCycler::Entry, 2> entries = {{first, second}};

  alignas(std::max_align_t) static uint8_t
      cycler_buf[PaletteCycler::required_arena_bytes()];
  Arena cycler_arena(cycler_buf, sizeof(cycler_buf));
  alignas(std::max_align_t) static uint8_t
      ref_buf[2 * BakedPalette::required_arena_bytes()];
  Arena ref_arena(ref_buf, sizeof(ref_buf));

  static_assert(!std::is_copy_constructible_v<PaletteCycler>);
  static_assert(!std::is_copy_constructible_v<GeneratedPaletteBank>);
  PaletteCycler cycler;
  cycler.init(cycler_arena, entries.data(), entries.size(), 3, 4);

  BakedPaletteStorage ref;
  ref.bake(ref_arena, first);
  expect_baked_equal(cycler.palette(), ref);
  HS_EXPECT_EQ(cycler.current_index(), 0);

  cycler.step();
  cycler.step();
  expect_baked_equal(cycler.palette(), ref);
  HS_EXPECT_FALSE(cycler.fading());

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());
  expect_baked_equal(cycler.palette(), ref);

  cycler.step();
  GenerativePalette expected_morph;
  expected_morph.morph_palettes(first, second, 0.25f);
  BakedPaletteStorage expected;
  expected.bake(ref_arena, expected_morph);
  expect_baked_equal(cycler.palette(), expected);

  cycler.step();
  cycler.step();
  cycler.step();
  HS_EXPECT_FALSE(cycler.fading());
  HS_EXPECT_EQ(cycler.current_index(), 1);
  ref.rebake(second);
  expect_baked_equal(cycler.palette(), ref);

  for (int i = 0; i < 7; ++i)
    cycler.step();
  HS_EXPECT_FALSE(cycler.fading());
  HS_EXPECT_EQ(cycler.current_index(), 0);
  ref.rebake(first);
  expect_baked_equal(cycler.palette(), ref);
}

inline void test_palette_cycler_heterogeneous_crossfade() {
  Gradient ramp{{0.0f, CPixel(10u, 250u, 60u)}, {1.0f, CPixel(250u, 10u, 60u)}};
  SolidColorPalette solid(Color4(Pixel(9000, 21000, 43000), 1.0f));
  const GenerativePalette generative(PaletteRecipes::balanced_analogous(0.3f));

  alignas(std::max_align_t) static uint8_t
      buf[8 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPaletteStorage prebaked;
  prebaked.bake(arena, solid);

  const std::array<PaletteCycler::Entry, 3> entries = {
      {ramp, prebaked.view(), generative}};

  PaletteCycler cycler;
  cycler.init(arena, entries.data(), entries.size(), 1, 2);

  BakedPaletteStorage ramp_lut;
  ramp_lut.bake(arena, ramp);
  expect_baked_equal(cycler.palette(), ramp_lut);

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());
  cycler.step();
  BakedPaletteStorage expected;
  expected.bake(arena, ramp);
  expected.rebake_crossfade(ramp_lut, prebaked, 0.5f);
  expect_baked_equal(cycler.palette(), expected);

  cycler.step();
  HS_EXPECT_FALSE(cycler.fading());
  HS_EXPECT_EQ(cycler.current_index(), 1);
  expect_baked_equal(cycler.palette(), prebaked);

  cycler.step();
  cycler.step();
  cycler.step();
  HS_EXPECT_EQ(cycler.current_index(), 2);
  BakedPaletteStorage generative_lut;
  generative_lut.bake(arena, generative);
  expect_baked_equal(cycler.palette(), generative_lut);

  cycler.step();
  cycler.step();
  cycler.step();
  HS_EXPECT_EQ(cycler.current_index(), 0);
  expect_baked_equal(cycler.palette(), ramp_lut);
}

inline void test_palette_cycler_pause_and_static() {
  Gradient ramp{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};

  alignas(std::max_align_t) static uint8_t
      buf[6 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  const std::array<PaletteCycler::Entry, 1> single = {{ramp}};
  PaletteCycler static_cycler;
  static_cycler.init(arena, single.data(), single.size(), 1, 1);
  BakedPaletteStorage ref;
  ref.bake(arena, ramp);
  for (int i = 0; i < 5; ++i)
    static_cycler.step();
  HS_EXPECT_FALSE(static_cycler.fading());
  expect_baked_equal(static_cycler.palette(), ref);

  const GenerativePalette a(PaletteRecipes::balanced_analogous(0.1f));
  const GenerativePalette b(PaletteRecipes::balanced_analogous(0.6f));
  const std::array<PaletteCycler::Entry, 2> pair = {{a, b}};
  bool paused = true;
  PaletteCycler cycler;
  cycler.init(arena, pair.data(), pair.size(), 1, 2, nullptr, &paused);
  for (int i = 0; i < 5; ++i)
    cycler.step();
  HS_EXPECT_FALSE(cycler.fading());
  HS_EXPECT_EQ(cycler.current_index(), 0);
  paused = false;
  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());
}

inline void test_standalone_palette_rotations_morph_compatible() {
  constexpr float GOLDEN_STEP = 0.618034f;
  float rotation = 0.0f;
  GenerativePalette liquid_prev(
      EffectPaletteRecipes::standalone_liquid_at(0.0f));
  GenerativePalette flyby_prev(EffectPaletteRecipes::standalone_flyby_at(0.0f));
  for (int i = 0; i < 24; ++i) {
    rotation = math::wrap_t(rotation + GOLDEN_STEP);
    GenerativePalette liquid(
        EffectPaletteRecipes::standalone_liquid_at(rotation));
    GenerativePalette flyby(
        EffectPaletteRecipes::standalone_flyby_at(rotation));
    HS_EXPECT_TRUE(liquid_prev.morph_compatible(liquid));
    HS_EXPECT_TRUE(flyby_prev.morph_compatible(flyby));
    liquid_prev = liquid;
    flyby_prev = flyby;
  }
}

/** @brief Deterministic PaletteCycler provider walking analogous base hues. */
inline void scripted_next_palette(void *context, uint32_t sequence,
                                  GenerativePalette &out) {
  ++*static_cast<int *>(context);
  out = GenerativePalette(
      PaletteRecipes::balanced_analogous(0.1f + 0.2f * sequence));
}

/** @brief PaletteCycler provider whose recipe carries a BELL chroma curve. */
inline void scripted_tonal_palette(void *context, uint32_t sequence,
                                   GenerativePalette &out) {
  ++*static_cast<int *>(context);
  out = GenerativePalette(
      PaletteRecipes::tonal_monochrome(0.1f + 0.2f * sequence));
}

inline void test_palette_cycler_generated_cycle() {
  alignas(std::max_align_t) static uint8_t
      buf[PaletteCycler::generated_arena_bytes() +
          3 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  int provider_calls = 0;
  PaletteCycler cycler;
  cycler.init_generated(arena, scripted_next_palette, &provider_calls, 0, 2);
  HS_EXPECT_EQ(provider_calls, 2);

  BakedPaletteStorage ref;
  ref.bake(arena, GenerativePalette(PaletteRecipes::balanced_analogous(0.1f)));
  expect_baked_equal(cycler.palette(), ref);

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());

  cycler.step();
  GenerativePalette mid;
  mid.morph_palettes(
      GenerativePalette(PaletteRecipes::balanced_analogous(0.1f)),
      GenerativePalette(PaletteRecipes::balanced_analogous(0.3f)), 0.5f);
  BakedPaletteStorage expected;
  expected.bake(arena, mid);
  expect_baked_equal(cycler.palette(), expected);

  cycler.step();
  HS_EXPECT_FALSE(cycler.fading());
  HS_EXPECT_EQ(provider_calls, 3);
  ref.rebake(GenerativePalette(PaletteRecipes::balanced_analogous(0.3f)));
  expect_baked_equal(cycler.palette(), ref);

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());

  // Re-initializing in entry mode must drop the generated-cycle state.
  const GenerativePalette first(PaletteRecipes::balanced_analogous(0.4f));
  const GenerativePalette second(PaletteRecipes::balanced_analogous(0.9f));
  const std::array<PaletteCycler::Entry, 2> pair = {{first, second}};
  alignas(std::max_align_t) static uint8_t
      reinit_buf[PaletteCycler::required_arena_bytes() +
                 BakedPalette::required_arena_bytes()];
  Arena reinit_arena(reinit_buf, sizeof(reinit_buf));
  cycler.init(reinit_arena, pair.data(), pair.size(), 0, 2);

  cycler.step();
  cycler.step();
  GenerativePalette blended;
  blended.morph_palettes(first, second, 0.5f);
  BakedPaletteStorage blended_lut;
  blended_lut.bake(reinit_arena, blended);
  expect_baked_equal(cycler.palette(), blended_lut);
  HS_EXPECT_EQ(provider_calls, 3);
}

/** @brief A constant-chroma write keeps the axis curve, so a generated cycle
 *  whose recipe carries a non-constant chroma curve stays morphable. */
inline void test_palette_cycler_generated_chroma_keeps_morph() {
  GenerativePalette tonal(PaletteRecipes::tonal_monochrome(0.3f));
  const GenerativePalette untouched = tonal;
  tonal.set_constant_chroma(0.4f);
  HS_EXPECT_TRUE(tonal.morph_compatible(untouched));
  for (int i = 0; i <= 4; ++i)
    HS_EXPECT_NEAR(tonal.diagnose(i * 0.25f).q, 0.4f, 1e-3f);

  alignas(std::max_align_t) static uint8_t
      buf[PaletteCycler::generated_arena_bytes() +
          3 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  int provider_calls = 0;
  PaletteCycler cycler;
  cycler.init_generated(arena, scripted_tonal_palette, &provider_calls, 0, 2);
  cycler.set_generated_chroma(0.4f);
  // Each finish_fade re-checks the retired slot against a fresh provider
  // palette still carrying the recipe's curve.
  for (int frame = 0; frame < 12; ++frame)
    cycler.step();
  HS_EXPECT_GT(provider_calls, 4);
}

inline void test_palette_cycler_roster_hidden_advance_catches_up() {
  const GenerativePalette a(PaletteRecipes::balanced_analogous(0.2f));
  const GenerativePalette b(PaletteRecipes::balanced_analogous(0.7f));
  const GenerativePalette c(PaletteRecipes::tonal_monochrome(0.4f));
  HS_EXPECT_TRUE(a.morph_compatible(b));
  HS_EXPECT_FALSE(a.morph_compatible(c));
  for (const auto *target : {&b, &c}) {
    const std::array<PaletteCycler::Entry, 2> entries = {{a, *target}};
    for (int frames : {1, 3, 4, 5, 8, 17, 50}) {
      alignas(std::max_align_t)
          std::array<uint8_t, 2 * PaletteCycler::required_arena_bytes()>
              buffer{};
      Arena arena(buffer.data(), buffer.size());
      PaletteCycler full, hidden;
      full.init(arena, entries.data(), entries.size(), 3, 4);
      hidden.init(arena, entries.data(), entries.size(), 3, 4);
      for (int frame = 0; frame < frames; ++frame) {
        full.step();
        hidden.advance_without_display();
      }
      full.step();
      hidden.step();
      expect_baked_equal(full.palette(), hidden.palette());
    }
  }
}

inline void test_palette_cycler_hidden_advance_catches_up() {
  alignas(std::max_align_t)
      std::array<uint8_t, PaletteCycler::generated_arena_bytes() + 32>
          full_buf{};
  alignas(std::max_align_t)
      std::array<uint8_t, PaletteCycler::generated_arena_bytes() + 32>
          hidden_buf{};
  Arena full_arena(full_buf.data(), full_buf.size());
  Arena hidden_arena(hidden_buf.data(), hidden_buf.size());
  int full_provider_calls = 0;
  int hidden_provider_calls = 0;
  PaletteCycler full;
  PaletteCycler hidden;
  full.init_generated(full_arena, scripted_next_palette, &full_provider_calls,
                      0, 17);
  hidden.init_generated(hidden_arena, scripted_next_palette,
                        &hidden_provider_calls, 0, 17);

  for (int frame = 0; frame < 50; ++frame) {
    full.step();
    hidden.advance_without_display();
  }
  full.step();
  hidden.step();
  expect_baked_equal(full.palette(), hidden.palette());
  HS_EXPECT_EQ(full_provider_calls, hidden_provider_calls);
}

/** @brief The bake generation advances over a cycle and never holds still
 *         across a display-LUT change. */
inline void test_palette_cycler_bake_generation() {
  const GenerativePalette a(PaletteRecipes::balanced_analogous(0.2f));
  const GenerativePalette b(PaletteRecipes::balanced_analogous(0.7f));
  const std::array<PaletteCycler::Entry, 1> single = {{a}};
  const std::array<PaletteCycler::Entry, 2> pair = {{a, b}};

  alignas(std::max_align_t) static uint8_t
      buf[2 * PaletteCycler::required_arena_bytes() +
          PaletteCycler::generated_arena_bytes() +
          BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  PaletteCycler held;
  held.init(arena, single.data(), single.size(), 1, 2);
  const uint32_t after_init = held.bake_generation();
  HS_EXPECT_NE(after_init, 0u);
  for (int i = 0; i < 4; ++i)
    held.step();
  HS_EXPECT_EQ(held.bake_generation(), after_init);

  PaletteCycler cycler;
  cycler.init(arena, pair.data(), pair.size(), 2, 4);
  BakedPaletteStorage snapshot;
  snapshot.bake(arena, a);
  snapshot.rebake_copy(cycler.palette());
  uint32_t generation = cycler.bake_generation();
  uint32_t advances = 0;
  for (int i = 0; i < 24; ++i) {
    cycler.step();
    if (cycler.bake_generation() == generation) {
      expect_baked_equal(cycler.palette(), snapshot);
      continue;
    }
    ++advances;
    generation = cycler.bake_generation();
    snapshot.rebake_copy(cycler.palette());
  }
  HS_EXPECT_GT(advances, 0u);

  int provider_calls = 0;
  PaletteCycler generated;
  generated.init_generated(arena, scripted_next_palette, &provider_calls, 0, 3);
  const uint32_t before_chroma = generated.bake_generation();
  generated.set_generated_chroma(0.4f);
  HS_EXPECT_NE(generated.bake_generation(), before_chroma);
  generated.advance_without_display();
  generated.advance_without_display();
  const uint32_t hidden_generation = generated.bake_generation();
  generated.set_generated_chroma(0.8f);
  HS_EXPECT_EQ(generated.bake_generation(), hidden_generation);
  generated.step();
  HS_EXPECT_GT(generated.bake_generation(), hidden_generation);
}

inline void test_palette_cycler_zero_dwell_chains_fades() {
  const GenerativePalette a(PaletteRecipes::balanced_analogous(0.2f));
  const GenerativePalette b(PaletteRecipes::balanced_analogous(0.7f));
  const std::array<PaletteCycler::Entry, 2> pair = {{a, b}};

  alignas(std::max_align_t) static uint8_t
      buf[PaletteCycler::required_arena_bytes() +
          2 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  PaletteCycler cycler;
  cycler.init(arena, pair.data(), pair.size(), 0, 2);

  BakedPaletteStorage ref;
  ref.bake(arena, a);
  expect_baked_equal(cycler.palette(), ref);
  HS_EXPECT_FALSE(cycler.fading());

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());

  cycler.step();
  GenerativePalette mid;
  mid.morph_palettes(a, b, 0.5f);
  BakedPaletteStorage expected;
  expected.bake(arena, mid);
  expect_baked_equal(cycler.palette(), expected);

  cycler.step();
  HS_EXPECT_FALSE(cycler.fading());
  HS_EXPECT_EQ(cycler.current_index(), 1);
  ref.rebake(b);
  expect_baked_equal(cycler.palette(), ref);

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());
}

/** @brief Hidden generated harmonies catch up with visible morph and chroma state. */
inline void test_generated_palette_bank_routes_and_rechromas() {
  enum class Mode { TRIADIC, COMPLEMENTARY, ANALOGOUS };
  alignas(std::max_align_t) static uint8_t
      storage[GeneratedPaletteBank::required_arena_bytes() +
              BakedPalette::required_arena_bytes()];
  Arena arena(storage, sizeof(storage));
  GeneratedPaletteBank bank;
  bank.init(arena, 0.62f, nullptr);
  BakedPaletteStorage expected;
  expected.bake(arena,
                GenerativePalette(PaletteRecipes::balanced_analogous(0.0f)));
  auto compare = [&](Mode mode, PaletteHarmony harmony, int from_index,
                     float progress, float chroma, bool rechromed) {
    GenerativePalette from(PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, harmony, AxisCurve::ASCENDING,
        math::wrap_t(from_index * 159.0f / 256.0f), chroma));
    GenerativePalette to(PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, harmony, AxisCurve::ASCENDING,
        math::wrap_t((from_index + 1) * 159.0f / 256.0f), chroma));
    if (rechromed) {
      from.set_constant_chroma(chroma);
      to.set_constant_chroma(chroma);
    }
    HS_EXPECT_TRUE(from.morph_compatible(to));
    GenerativePalette morph;
    morph.morph_palettes(from, to, progress);
    expected.rebake(morph);
    hs_test::color_tests::expect_baked_equal(bank.palette(mode), expected);
  };
  for (int i = 0; i < 151; ++i)
    bank.step(Mode::TRIADIC);
  compare(Mode::TRIADIC, PaletteHarmony::TRIADIC, 0, 0.25f, 0.62f, false);
  bank.step(Mode::ANALOGOUS);
  compare(Mode::ANALOGOUS, PaletteHarmony::ANALOGOUS, 0, 151.0f / 600.0f, 0.62f,
          false);
  bank.set_chroma(0.4f);
  constexpr Mode modes[] = {Mode::TRIADIC, Mode::COMPLEMENTARY,
                            Mode::ANALOGOUS};
  constexpr PaletteHarmony harmonies[] = {PaletteHarmony::TRIADIC,
                                          PaletteHarmony::COMPLEMENTARY,
                                          PaletteHarmony::ANALOGOUS};
  for (int i = 0; i < 3; ++i) {
    bank.step(modes[i]);
    compare(modes[i], harmonies[i], 0, (152.0f + i) / 600.0f, 0.4f, true);
  }
  for (int i = 155; i < 602; ++i)
    bank.step(Mode::TRIADIC);
  bank.step(Mode::COMPLEMENTARY);
  compare(Mode::COMPLEMENTARY, PaletteHarmony::COMPLEMENTARY, 1, 1.0f / 600.0f,
          0.4f, true);
  uint32_t hue = 0;
  GenerativePalette previous;
  for (uint32_t sequence = 0; sequence < 3; ++sequence) {
    GenerativePalette generated;
    GeneratedPaletteBank::next_palette(hue, sequence, PaletteHarmony::ANALOGOUS,
                                       0.4f, generated);
    HS_EXPECT_EQ(hue, sequence * 159u);
    const GenerativePalette reference(PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, PaletteHarmony::ANALOGOUS,
        AxisCurve::ASCENDING, math::wrap_t(sequence * 159.0f / 256.0f), 0.4f));
    for (int i = 0; i <= 8; ++i)
      HS_EXPECT_EQ(generated.get(i / 8.0f).color,
                   reference.get(i / 8.0f).color);
    if (sequence > 0)
      HS_EXPECT_TRUE(previous.morph_compatible(generated));
    previous = generated;
  }
}
