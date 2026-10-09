/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Generative-palette tests and the EffectPaletteRecipes roster contract.
 */

/**
 * @brief Pins a compiled TRIADIC/BELL ramp to frozen colors at nine stops.
 * @details A change detector captured from `got` under the native clang
 *          toolchain, not an independent oracle; re-derive only after a
 *          deliberate retune. The tolerance is relative so the dark stops stay
 *          pinned.
 */
inline void test_generative_palette_deterministic() {
  const PaletteRecipe recipe = PaletteRecipes::profile(
      PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC, AxisCurve::BELL,
      PaletteRecipes::hue_turns(77), 0.86f);
  const GenerativePalette palette(recipe);
  const uint16_t expected[9][3] = {
      {298, 291, 32},       {1843, 3180, 346},     {1447, 12521, 7570},
      {3228, 28254, 35424}, {28456, 48750, 60539}, {4818, 24745, 55520},
      {9819, 2840, 51712},  {6742, 467, 6101},     {807, 48, 375}};
  const auto tolerance = [](uint16_t value) {
    return std::fmax(4.0f, 0.001f * static_cast<float>(value));
  };
  for (int i = 0; i < 9; ++i) {
    const Pixel got = palette.get(i / 8.0f).color;
    HS_EXPECT_NEAR(got.r, expected[i][0], tolerance(expected[i][0]));
    HS_EXPECT_NEAR(got.g, expected[i][1], tolerance(expected[i][1]));
    HS_EXPECT_NEAR(got.b, expected[i][2], tolerance(expected[i][2]));
  }
}

inline void test_effect_palette_recipe_roster() {
  namespace R = EffectPaletteRecipes;
  const float hue = PaletteRecipes::hue_turns(42);
  const R::Preset expected[] = {
      {"BZReactionDiffusion", false, R::bz_reaction_diffusion()},
      {"Comets", true, R::comets(hue)},
      {"DisplacementField / RingShower", true, R::displacement_field(hue)},
      {"Dynamo", true, R::dynamo(hue)},
      {"GSReactionDiffusion", true, R::gs_reaction_diffusion(hue)},
      {"MobiusRings", true, R::mobius_rings(hue)},
      {"Raymarch", false, R::raymarch()},
      {"Standalone Liquid", false, R::standalone_liquid()},
      {"Standalone Flyby", false, R::standalone_flyby()},
      {"HyperLattice", false, R::hyper_lattice()},
      {"MindSplatter", true, R::mind_splatter(hue)}};
  const auto presets = R::presets();
  static_assert(std::tuple_size_v<decltype(presets)> == std::size(expected));
  for (size_t i = 0; i < presets.size(); ++i) {
    HS_EXPECT_STREQ(presets[i].name, expected[i].name);
    HS_EXPECT_EQ(presets[i].random_hue, expected[i].random_hue);
    const GenerativePalette got(presets[i].recipe);
    const GenerativePalette want(expected[i].recipe);
    for (int k = 0; k <= 8; ++k) {
      const Pixel a = got.get(k / 8.0f).color;
      const Pixel b = want.get(k / 8.0f).color;
      HS_EXPECT_EQ(a.r, b.r);
      HS_EXPECT_EQ(a.g, b.g);
      HS_EXPECT_EQ(a.b, b.b);
    }
  }
  HS_EXPECT_NEAR(presets[1].recipe.hue.base_turns,
                 PaletteRecipes::hue_turns(42), 1e-6f);
  HS_EXPECT_EQ(presets[0].recipe.hue.mode, HueMode::CUSTOM);
  HS_EXPECT_EQ(presets[7].recipe.hue.mode, HueMode::CUSTOM);
  HS_EXPECT_EQ(presets[7].recipe.lightness.curve, AxisCurve::CUSTOM);
  HS_EXPECT_NEAR(presets[7].recipe.hue.custom_turns[1] -
                     presets[7].recipe.hue.custom_turns[0],
                 0.5f, 1e-6f);
}

inline void test_generative_palette_recipe_validation() {
  PaletteRecipe input;
  input.hue.base_turns = 1.25f;
  input.hue.spread_turns = 0.5f;
  input.lightness.range = -0.2f;
  input.falloff_start = 0.75f;
  input.input.offset = 0.8f;
  input.input.span = 0.5f;

  GenerativePalette output;
  PaletteRecipe canonical;
  PaletteCompileStatus status;
  HS_EXPECT_TRUE(
      GenerativePalette::try_compile(input, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::OK);
  HS_EXPECT_NEAR(canonical.hue.base_turns, 0.25f, 1e-6f);
  HS_EXPECT_NEAR(canonical.hue.spread_turns, 0.25f, 1e-6f);
  HS_EXPECT_NEAR(canonical.lightness.range, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(canonical.falloff_start, 0.90f, 1e-6f);
  HS_EXPECT_NEAR(canonical.input.offset, 0.8f, 1e-6f);
  HS_EXPECT_NEAR(canonical.input.span, 0.2f, 1e-6f);
  HS_EXPECT_TRUE(status.adjustments.wrapped_fields != 0);
  HS_EXPECT_TRUE(status.adjustments.clamped_fields != 0);
  HS_EXPECT_TRUE(status.adjustments.canonicalized_fields != 0);

  {
    PaletteRecipe loop;
    loop.domain = PaletteDomain::LOOP;
    loop.hue.mode = HueMode::SWEEP;
    loop.hue.sweep_turns = 1.0f + 5e-7f;
    GenerativePalette compiled;
    PaletteRecipe snapped;
    PaletteCompileStatus adjusted;
    HS_EXPECT_TRUE(
        GenerativePalette::try_compile(loop, compiled, snapped, adjusted));
    HS_EXPECT_EQ(snapped.hue.sweep_turns, 1.0f);
    HS_EXPECT_TRUE((adjusted.adjustments.canonicalized_fields &
                    (uint64_t{1} << static_cast<uint8_t>(
                         PaletteRecipeField::SWEEP_TURNS))) != 0);
  }
  const GenerativePalette::Snapshot before = output.snapshot();
  const PaletteRecipe canonical_before = canonical;
  input.lightness.center = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(
      GenerativePalette::try_compile(input, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::NON_FINITE);
  HS_EXPECT_EQ(status.field, PaletteRecipeField::LIGHTNESS_CENTER);
  const GenerativePalette::Snapshot after = output.snapshot();
  HS_EXPECT_EQ(std::memcmp(&before, &after, sizeof(before)), 0);
  HS_EXPECT_TRUE(canonical_before == canonical);

  struct RelationshipCount {
    PaletteHarmony harmony;
    uint8_t count;
  };
  constexpr RelationshipCount relationships[] = {
      {PaletteHarmony::MONOCHROMATIC, 2},
      {PaletteHarmony::ANALOGOUS, 3},
      {PaletteHarmony::ACCENTED_ANALOGOUS, 4},
      {PaletteHarmony::COMPLEMENTARY, 2},
      {PaletteHarmony::SPLIT_COMPLEMENTARY, 3},
      {PaletteHarmony::TRIADIC, 3},
      {PaletteHarmony::TETRADIC, 4},
      {PaletteHarmony::SQUARE, 4},
  };
  for (const auto relationship : relationships) {
    PaletteRecipe derived;
    derived.hue.harmony = relationship.harmony;
    HS_EXPECT_TRUE(
        GenerativePalette::try_compile(derived, output, canonical, status));
    HS_EXPECT_EQ(output.palette_key_count(), relationship.count);
  }
}

inline void test_generative_palette_canonical_ignores_inactive_fields() {
  PaletteRecipe first;
  first.hue.mode = HueMode::HARMONY;
  first.hue.harmony = PaletteHarmony::TRIADIC;
  first.lightness.curve = AxisCurve::CONSTANT;
  first.chroma.axis.curve = AxisCurve::CONSTANT;

  PaletteRecipe second = first;
  second.hue.spread_turns = std::numeric_limits<float>::quiet_NaN();
  second.hue.custom_turns.fill(std::numeric_limits<float>::quiet_NaN());
  second.lightness.custom.fill(std::numeric_limits<float>::quiet_NaN());
  second.chroma.axis.custom.fill(std::numeric_limits<float>::quiet_NaN());
  second.falloff_start = std::numeric_limits<float>::quiet_NaN();
  second.lightness.range = std::numeric_limits<float>::quiet_NaN();
  second.chroma.axis.range = std::numeric_limits<float>::quiet_NaN();

  GenerativePalette first_palette;
  GenerativePalette second_palette;
  PaletteRecipe first_canonical;
  PaletteRecipe second_canonical;
  PaletteCompileStatus status;
  HS_EXPECT_TRUE(GenerativePalette::try_compile(first, first_palette,
                                                first_canonical, status));
  HS_EXPECT_TRUE(GenerativePalette::try_compile(second, second_palette,
                                                second_canonical, status));
  HS_EXPECT_TRUE(first_canonical == second_canonical);
  const auto first_snapshot = first_palette.snapshot();
  const auto second_snapshot = second_palette.snapshot();
  HS_EXPECT_EQ(
      std::memcmp(&first_snapshot, &second_snapshot, sizeof(first_snapshot)),
      0);
}

inline void test_generative_palette_input_window() {
  PaletteRecipe recipe =
      PaletteRecipes::profile(PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC,
                              AxisCurve::CONSTANT, 0.17f, 0.86f);
  const GenerativePalette full(recipe);

  recipe.input.offset = 0.2f;
  recipe.input.span = 0.4f;
  const GenerativePalette windowed(recipe);
  HS_EXPECT_NEAR(windowed.palette_input_offset(), 0.2f, 1e-6f);
  HS_EXPECT_NEAR(windowed.palette_input_span(), 0.4f, 1e-6f);
  for (int i = 0; i <= 16; ++i) {
    const float t = i / 16.0f;
    const Pixel actual = windowed.get(t).color;
    const Pixel expected = full.get(0.2f + 0.4f * t).color;
    HS_EXPECT_EQ(actual.r, expected.r);
    HS_EXPECT_EQ(actual.g, expected.g);
    HS_EXPECT_EQ(actual.b, expected.b);
  }

  recipe.lightness.curve = AxisCurve::ASCENDING;
  recipe.lightness.center = 0.5f;
  recipe.lightness.range = 0.6f;
  PaletteRecipe full_envelope_recipe = recipe;
  full_envelope_recipe.input.offset = 0.0f;
  full_envelope_recipe.input.span = 1.0f;
  const GenerativePalette full_envelope(full_envelope_recipe);
  const GenerativePalette windowed_envelope(recipe);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(0.0f).L,
                 full_envelope.diagnose(0.2f).L, 1e-5f);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(1.0f).L,
                 full_envelope.diagnose(0.6f).L, 1e-5f);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(0.0f).h_path,
                 full.diagnose(0.2f).h_path, 1e-5f);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(1.0f).h_path,
                 full.diagnose(0.6f).h_path, 1e-5f);

  recipe.domain = PaletteDomain::MIRROR;
  const GenerativePalette cropped_mirror(recipe);
  HS_EXPECT_TRUE(cropped_mirror.mirrors_domain());
  for (int i = 0; i <= 16; ++i) {
    const Pixel left = cropped_mirror.get(i / 32.0f).color;
    const Pixel right = cropped_mirror.get(1.0f - i / 32.0f).color;
    HS_EXPECT_EQ(left.r, right.r);
    HS_EXPECT_EQ(left.g, right.g);
    HS_EXPECT_EQ(left.b, right.b);
  }
  const Pixel window_end = windowed_envelope.get(1.0f).color;
  const Pixel mirror_middle = cropped_mirror.get(0.5f).color;
  HS_EXPECT_EQ(mirror_middle.r, window_end.r);
  HS_EXPECT_EQ(mirror_middle.g, window_end.g);
  HS_EXPECT_EQ(mirror_middle.b, window_end.b);

  recipe.domain = PaletteDomain::VIGNETTE;
  const GenerativePalette cropped_vignette(recipe);
  HS_EXPECT_EQ(cropped_vignette.get(0.0f).color.r, 0);
  HS_EXPECT_EQ(cropped_vignette.get(0.0f).color.g, 0);
  HS_EXPECT_EQ(cropped_vignette.get(0.0f).color.b, 0);
  HS_EXPECT_EQ(cropped_vignette.get(1.0f).color.r, 0);
  HS_EXPECT_EQ(cropped_vignette.get(1.0f).color.g, 0);
  HS_EXPECT_EQ(cropped_vignette.get(1.0f).color.b, 0);
  const Pixel window_middle = windowed_envelope.get(0.5f).color;
  const Pixel vignette_middle = cropped_vignette.get(0.5f).color;
  HS_EXPECT_EQ(vignette_middle.r, window_middle.r);
  HS_EXPECT_EQ(vignette_middle.g, window_middle.g);
  HS_EXPECT_EQ(vignette_middle.b, window_middle.b);

  recipe.domain = PaletteDomain::FALLOFF;
  const GenerativePalette cropped_falloff(recipe);
  HS_EXPECT_EQ(cropped_falloff.get(1.0f).color.r, 0);
  HS_EXPECT_EQ(cropped_falloff.get(1.0f).color.g, 0);
  HS_EXPECT_EQ(cropped_falloff.get(1.0f).color.b, 0);
  const Pixel falloff_color_end = cropped_falloff.get(2.0f / 3.0f).color;
  HS_EXPECT_EQ(falloff_color_end.r, window_end.r);
  HS_EXPECT_EQ(falloff_color_end.g, window_end.g);
  HS_EXPECT_EQ(falloff_color_end.b, window_end.b);

  recipe.domain = PaletteDomain::LOOP;
  const GenerativePalette cropped_loop(recipe);
  HS_EXPECT_TRUE(cropped_loop.loops_domain());
  const Pixel loop_first = cropped_loop.get(0.0f).color;
  const Pixel loop_last = cropped_loop.get(1.0f).color;
  HS_EXPECT_EQ(loop_first.r, loop_last.r);
  HS_EXPECT_EQ(loop_first.g, loop_last.g);
  HS_EXPECT_EQ(loop_first.b, loop_last.b);
  const Pixel loop_color_end = cropped_loop.get(2.0f / 3.0f).color;
  HS_EXPECT_EQ(loop_color_end.r, window_end.r);
  HS_EXPECT_EQ(loop_color_end.g, window_end.g);
  HS_EXPECT_EQ(loop_color_end.b, window_end.b);
}

inline void test_generative_palette_resolves_axes_and_harmony() {
  PaletteRecipe recipe = PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                                 PaletteHarmony::COMPLEMENTARY,
                                                 AxisCurve::ASCENDING, 0.0f);
  recipe.lightness.center = 0.5f;
  recipe.lightness.range = 0.6f;
  const GenerativePalette palette(recipe);
  const auto keys = palette.snapshot();
  const auto key0 = GenerativePalette::snapshot_key(keys, 0);
  const auto key1 = GenerativePalette::snapshot_key(keys, 1);

  HS_EXPECT_EQ(keys.key_count, 2);
  HS_EXPECT_NEAR(key0.L, 0.2f, 3e-4f);
  HS_EXPECT_NEAR(key1.L, 0.8f, 3e-4f);
  HS_EXPECT_NEAR(key1.h - key0.h, math::PI_F, 1e-5f);

  recipe.lightness.curve = AxisCurve::BELL;
  const GenerativePalette bell(recipe);
  HS_EXPECT_NEAR(bell.diagnose(0.0f).L, 0.2f, 1e-5f);
  HS_EXPECT_NEAR(bell.diagnose(0.5f).L, 0.8f, 1e-5f);
  HS_EXPECT_NEAR(bell.diagnose(1.0f).L, 0.2f, 1e-5f);

  recipe.hue.direction = HueDirection::CLOCKWISE;
  const auto clockwise = GenerativePalette(recipe).snapshot();
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(clockwise, 1).h -
                     GenerativePalette::snapshot_key(clockwise, 0).h,
                 -math::PI_F, 1e-5f);

  recipe.hue.direction = HueDirection::COUNTERCLOCKWISE;
  recipe.hue.harmony = PaletteHarmony::TETRADIC;
  recipe.hue.spread_turns = 1.0f / 6.0f;
  const auto tetradic = GenerativePalette(recipe).snapshot();
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(tetradic, 1).h -
                     GenerativePalette::snapshot_key(tetradic, 0).h,
                 math::PI_F / 3.0f, 1e-5f);
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(tetradic, 2).h -
                     GenerativePalette::snapshot_key(tetradic, 0).h,
                 math::PI_F, 1e-5f);

  recipe.hue.harmony = PaletteHarmony::SQUARE;
  const auto square = GenerativePalette(recipe).snapshot();
  auto oversized = square;
  oversized.key_count = 255;
  const auto bounded_last = GenerativePalette::snapshot_key(square, 3);
  const auto oversized_last = GenerativePalette::snapshot_key(oversized, 254);
  HS_EXPECT_EQ(oversized_last.L, bounded_last.L);
  HS_EXPECT_EQ(oversized_last.chroma, bounded_last.chroma);
  HS_EXPECT_EQ(oversized_last.h, bounded_last.h);
  for (int i = 1; i < 4; ++i)
    HS_EXPECT_NEAR(GenerativePalette::snapshot_key(square, i).h -
                       GenerativePalette::snapshot_key(square, i - 1).h,
                   0.5f * math::PI_F, 1e-5f);
}

/**
 * @brief Torsion shears hue by lightness, is range-validated, and gates morphs.
 */
inline void test_generative_palette_hue_torsion() {
  constexpr float TORSION = 0.75f;
  PaletteRecipe recipe = PaletteRecipes::profile(
      PaletteDomain::STRAIGHT, PaletteHarmony::ANALOGOUS, AxisCurve::ASCENDING,
      0.1f, 0.5f);
  const GenerativePalette flat(recipe);
  recipe.hue_torsion = TORSION;
  const GenerativePalette sheared(recipe);
  HS_EXPECT_NEAR(flat.palette_hue_torsion(), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(sheared.palette_hue_torsion(), TORSION, 1e-6f);

  const int key_count = flat.palette_key_count();
  HS_EXPECT_EQ(key_count, 3);
  for (int i = 0; i < key_count; ++i) {
    const OKLCH base = flat.resolved_oklch_key(i);
    const OKLCH torsioned = sheared.resolved_oklch_key(i);
    HS_EXPECT_NEAR(torsioned.L, base.L, 1e-6f);
    HS_EXPECT_NEAR(torsioned.h, base.h + TORSION * (base.L - 0.5f), 1e-5f);
  }
  // The end keys sit off mid-lightness, so the shear is not a no-op.
  const float end_shear =
      sheared.resolved_oklch_key(0).h - flat.resolved_oklch_key(0).h;
  HS_EXPECT_GT(fabsf(end_shear), 0.1f);

  for (const float t : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
    const auto plain = flat.diagnose(t);
    HS_EXPECT_NEAR(plain.h_final, plain.h_path, 1e-6f);
    const auto twisted = sheared.diagnose(t);
    HS_EXPECT_NEAR(twisted.h_final - twisted.h_path,
                   TORSION * (twisted.L - 0.5f), 1e-5f);
  }

  GenerativePalette output;
  PaletteRecipe canonical;
  PaletteCompileStatus status;
  PaletteRecipe over = PaletteRecipes::balanced_analogous(0.0f);
  over.hue_torsion = 100.0f;
  HS_EXPECT_FALSE(
      GenerativePalette::try_compile(over, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::HUE_LIMIT);
  HS_EXPECT_EQ(status.field, PaletteRecipeField::HUE_TORSION);

  over.hue_torsion = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(
      GenerativePalette::try_compile(over, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::NON_FINITE);
  HS_EXPECT_EQ(status.field, PaletteRecipeField::HUE_TORSION);

  PaletteRecipe legal = PaletteRecipes::balanced_analogous(0.0f);
  legal.hue_torsion = -TORSION;
  HS_EXPECT_TRUE(
      GenerativePalette::try_compile(legal, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::OK);
  HS_EXPECT_NEAR(canonical.hue_torsion, -TORSION, 1e-6f);
  HS_EXPECT_NEAR(output.palette_hue_torsion(), -TORSION, 1e-6f);

  PaletteRecipe untwisted = PaletteRecipes::balanced_analogous(0.0f);
  PaletteRecipe twisted_recipe = untwisted;
  twisted_recipe.hue_torsion = TORSION;
  const GenerativePalette morph_a(untwisted);
  const GenerativePalette morph_b(twisted_recipe);
  HS_EXPECT_FALSE(morph_a.morph_compatible(morph_b));
  HS_EXPECT_FALSE(morph_b.morph_compatible(morph_a));
  untwisted.hue_torsion = TORSION;
  HS_EXPECT_TRUE(GenerativePalette(untwisted).morph_compatible(morph_b));
}

inline void test_generative_palette_blue_cusp_is_continuous() {
  PaletteRecipe recipe;
  recipe.hue.mode = HueMode::SWEEP;
  recipe.hue.base_turns = 98.0f / 256.0f;
  recipe.hue.sweep_turns = 1.0f;
  recipe.chroma.axis.center = 1.0f;
  recipe.chroma.headroom = 1.0f;
  recipe.lightness.curve = AxisCurve::ASCENDING;
  recipe.lightness.center = 0.495f;
  recipe.lightness.range = 0.41f;
  const GenerativePalette palette(recipe);

  auto previous = palette.diagnose(0.0f);
  float largest_chroma_step = 0.0f;
  for (int i = 1; i < 256; ++i) {
    const auto current = palette.diagnose(i / 255.0f);
    largest_chroma_step =
        hs_test::fold_worst(largest_chroma_step, fabsf(current.C - previous.C));
    HS_EXPECT_FALSE(current.fallback_mapped);
    previous = current;
  }
  HS_EXPECT_LT(largest_chroma_step, 0.04f);
}

inline void test_generative_palette_local_gamut_stays_in_gamut() {
  PaletteRecipe recipe =
      PaletteRecipes::profile(PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC,
                              AxisCurve::BELL, 0.17f, 0.86f);
  recipe.lightness.center = 0.52f;
  recipe.lightness.range = 0.72f;
  const GenerativePalette palette(recipe);
  for (int i = 0; i < 256; ++i) {
    const auto diagnostic = palette.diagnose(i / 255.0f);
    HS_EXPECT_FALSE(diagnostic.fallback_mapped);
    HS_EXPECT_TRUE(diagnostic.C <= diagnostic.C_max + 1e-5f);
  }
}

inline void test_generative_palette_domain_invariants() {
  const GenerativePalette mirror(
      PaletteRecipes::profile(PaletteDomain::MIRROR, PaletteHarmony::ANALOGOUS,
                              AxisCurve::BELL, 0.31f));
  alignas(std::max_align_t)
      uint8_t storage[BakedPalette::required_arena_bytes()];
  Arena arena(storage, sizeof(storage));
  BakedPaletteStorage baked;
  baked.bake(arena, mirror);
  for (int i = 0; i < 128; ++i) {
    const Pixel left = baked.get(i / 255.0f).color;
    const Pixel right = baked.get((255 - i) / 255.0f).color;
    HS_EXPECT_EQ(left.r, right.r);
    HS_EXPECT_EQ(left.g, right.g);
    HS_EXPECT_EQ(left.b, right.b);
  }

  const GenerativePalette loop(PaletteRecipes::isolight_spectral_loop(0.13f));
  arena.reset();
  baked.bake(arena, loop);
  const Color4 first = baked.get(0.0f);
  const Color4 last = baked.get(1.0f);
  HS_EXPECT_EQ(first.color.r, last.color.r);
  HS_EXPECT_EQ(first.color.g, last.color.g);
  HS_EXPECT_EQ(first.color.b, last.color.b);
  HS_EXPECT_EQ(first.alpha, last.alpha);
}

inline void test_generative_palette_morph_policy_contracts() {
  const GenerativePalette from(PaletteRecipes::balanced_analogous(0.1f));
  const GenerativePalette to(PaletteRecipes::balanced_analogous(0.3f));
  PaletteRecipe policy_recipe = PaletteRecipes::isolight_spectral_loop(0.5f);
  policy_recipe.color_path = ColorPath::OKLAB_CARTESIAN;
  policy_recipe.chroma.headroom = 0.8f;
  const GenerativePalette policy(policy_recipe);
  HS_EXPECT_TRUE(from.morph_compatible(to));
  HS_EXPECT_TRUE(policy.palette_domain() != from.palette_domain());
  HS_EXPECT_TRUE(policy.palette_color_path() != from.palette_color_path());
  HS_EXPECT_TRUE(policy.palette_headroom() != from.palette_headroom());
  GenerativePalette whole = policy;
  whole.morph_palettes(from, to, 0.5f);
  HS_EXPECT_EQ(whole.palette_domain(), from.palette_domain());
  HS_EXPECT_EQ(whole.palette_color_path(), from.palette_color_path());
  HS_EXPECT_EQ(whole.palette_headroom(), from.palette_headroom());
  GenerativePalette keys = policy;
  keys.morph_snapshots(from.snapshot(), to.snapshot(), 0.5f);
  HS_EXPECT_EQ(keys.palette_domain(), policy.palette_domain());
  HS_EXPECT_EQ(keys.palette_color_path(), policy.palette_color_path());
  HS_EXPECT_EQ(keys.palette_headroom(), policy.palette_headroom());
}

inline void test_generative_palette_snapshot_lerp() {
  GenerativePalette from(PaletteRecipes::balanced_analogous(0.0f));
  const GenerativePalette to(PaletteRecipes::balanced_analogous(0.75f));
  const auto first = from.snapshot();
  const auto last = to.snapshot();

  from.morph_snapshots(first, last, 0.5f);
  const auto middle = from.snapshot();
  const auto first_key = GenerativePalette::snapshot_key(first, 0);
  const auto last_key = GenerativePalette::snapshot_key(last, 0);
  const auto middle_key = GenerativePalette::snapshot_key(middle, 0);
  HS_EXPECT_NEAR(middle_key.L, 0.5f * (first_key.L + last_key.L), 3e-4f);
  HS_EXPECT_NEAR(middle_key.chroma, 0.5f * (first_key.chroma + last_key.chroma),
                 3e-4f);

  from.morph_snapshots(first, last, 1.0f);
  const GenerativePalette::Snapshot target = from.snapshot();
  HS_EXPECT_EQ(std::memcmp(&target, &last, sizeof(last)), 0);
}

/**
 * @brief Pins the key-morph hue path against a chain that accumulates past half
 *        a turn.
 * @details Each adjacent hue delta stays under half a turn, but the third key's
 * summed travel exceeds it; folding that into (-pi, pi] would sweep it backwards.
 */
inline void test_generative_palette_lerp_accumulates_segment_deltas() {
  PaletteRecipe recipe;
  recipe.hue.mode = HueMode::CUSTOM;
  const GenerativePalette from(recipe);
  recipe.hue.custom_turns[1] = 0.4f;
  recipe.hue.custom_turns[2] = 0.8f;
  const GenerativePalette to(recipe);
  HS_EXPECT_TRUE(from.morph_compatible(to));
  HS_EXPECT_EQ(from.palette_key_count(), (uint8_t)3);

  GenerativePalette morph;
  morph.morph_palettes(from, to, 0.5f);
  for (int i = 0; i < 3; ++i) {
    const float start = from.resolved_oklch_key(i).h;
    const float travel = to.resolved_oklch_key(i).h - start;
    HS_EXPECT_NEAR(morph.resolved_oklch_key(i).h, start + 0.5f * travel, 1e-3f);
  }
}

/**
 * @brief Pins the snapshot encode against dropping a key's hue.
 * @details The snapshot chroma quantum is coarser than the is_chromatic()
 * threshold, so the encode lifts chromatic keys to one quantum rather than
 * rounding them to gray.
 */
inline void test_generative_palette_snapshot_keeps_faint_chroma_chromatic() {
  PaletteRecipe recipe;
  recipe.chroma.axis.curve = AxisCurve::CUSTOM;
  for (int i = 0; i < PALETTE_MAX_KEYS; ++i)
    recipe.chroma.axis.custom[i] = 0.5f;
  recipe.chroma.axis.custom[1] = 1e-5f;
  const GenerativePalette palette(recipe);
  const GenerativePalette::Snapshot snapshot = palette.snapshot();
  const float chroma = GenerativePalette::snapshot_key(snapshot, 1).chroma;
  HS_EXPECT_GT(chroma, 0.0f);
  HS_EXPECT_GE(chroma, OKLCH_ACHROMATIC_C);
  HS_EXPECT_LT(chroma, 1e-3f);
}

inline void test_generative_palette_snapshot_keeps_absolute_gray_achromatic() {
  const GenerativePalette from(PaletteRecipes::from_colors(
      PaletteDomain::STRAIGHT, CPixel(255, 255, 255), CPixel(255, 0, 0),
      CPixel(0, 0, 255)));
  const GenerativePalette to(
      PaletteRecipes::from_colors(PaletteDomain::STRAIGHT, CPixel(0, 255, 0),
                                  CPixel(255, 0, 0), CPixel(0, 0, 255)));
  HS_EXPECT_TRUE(from.morph_compatible(to));
  HS_EXPECT_LT(GenerativePalette::snapshot_key(from.snapshot(), 0).chroma,
               OKLCH_ACHROMATIC_C);
  for (const float progress : {0.25f, 0.5f, 0.75f}) {
    GenerativePalette morph;
    morph.morph_palettes(from, to, progress);
    HS_EXPECT_NEAR(morph.resolved_oklch_key(0).h, to.resolved_oklch_key(0).h,
                   1e-4f);
  }
}

inline void test_generative_palette_lerp_target_aliases_this() {
  const GenerativePalette from(PaletteRecipes::balanced_analogous(0.0f));
  const GenerativePalette to(PaletteRecipes::balanced_analogous(0.25f));
  GenerativePalette expected;
  expected.morph_palettes(from, to, 0.25f);
  GenerativePalette aliased = to;
  aliased.morph_palettes(from, aliased, 0.25f);
  for (int i = 0; i <= 16; ++i) {
    const Pixel landed = aliased.get(i / 16.0f).color;
    const Pixel reference = expected.get(i / 16.0f).color;
    HS_EXPECT_EQ(landed.r, reference.r);
    HS_EXPECT_EQ(landed.g, reference.g);
    HS_EXPECT_EQ(landed.b, reference.b);
  }
}

inline void test_generative_palette_snapshot_lerp_closes_loop() {
  GenerativePalette morph(PaletteRecipes::isolight_spectral_loop(0.13f));
  const GenerativePalette target(PaletteRecipes::isolight_spectral_loop(0.37f));
  morph.morph_snapshots(morph.snapshot(), target.snapshot(), 1.0f);
  for (int i = 0; i <= 16; ++i) {
    const Pixel landed = morph.get(i / 16.0f).color;
    const Pixel expected = target.get(i / 16.0f).color;
    HS_EXPECT_EQ(landed.r, expected.r);
    HS_EXPECT_EQ(landed.g, expected.g);
    HS_EXPECT_EQ(landed.b, expected.b);
  }
}

inline void test_generative_palette_cartesian_path_neutralizes_midpoint() {
  PaletteRecipe arc_recipe = PaletteRecipes::profile(
      PaletteDomain::STRAIGHT, PaletteHarmony::COMPLEMENTARY,
      AxisCurve::CONSTANT, 0.0f, 0.72f);
  PaletteRecipe cartesian_recipe = arc_recipe;
  cartesian_recipe.color_path = ColorPath::OKLAB_CARTESIAN;
  const GenerativePalette arc(arc_recipe);
  const GenerativePalette cartesian(cartesian_recipe);
  HS_EXPECT_TRUE(cartesian.diagnose(0.25f).C < arc.diagnose(0.25f).C);
  HS_EXPECT_TRUE(cartesian.diagnose(0.25f).C <
                 cartesian.diagnose(0.0f).C * 0.8f);
  HS_EXPECT_TRUE(cartesian.diagnose(0.5f).C <
                 cartesian.diagnose(0.0f).C * 0.4f);
  HS_EXPECT_TRUE(cartesian.diagnose(0.75f).C <
                 cartesian.diagnose(1.0f).C * 0.8f);
}

inline void test_generative_palette_rejects_unavailable_path_minimum() {
  PaletteRecipe input;
  input.chroma.basis = ChromaBasis::PATH_MINIMUM;
  GenerativePalette output;
  PaletteRecipe canonical;
  PaletteCompileStatus status;
  HS_EXPECT_FALSE(
      GenerativePalette::try_compile(input, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::INVALID_ENUM);
  HS_EXPECT_EQ(status.field, PaletteRecipeField::CHROMA_BASIS);
}

inline void test_generative_palette_absolute_basis_canonicalizes_headroom() {
  const PaletteRecipe GRAY_BLUE =
      PaletteRecipes::from_colors(PaletteDomain::STRAIGHT, CPixel(0, 0, 0),
                                  CPixel(0, 0, 255), CPixel(255, 255, 255));
  HS_EXPECT_NEAR(GRAY_BLUE.hue.custom_turns[0], GRAY_BLUE.hue.custom_turns[1],
                 1e-6f);
  HS_EXPECT_NEAR(GRAY_BLUE.hue.custom_turns[2], GRAY_BLUE.hue.custom_turns[1],
                 1e-6f);
  PaletteRecipe input = PaletteRecipes::from_oklch_keys(
      PaletteDomain::STRAIGHT, OKLCH{0.5f, 0.10f, 0.0f},
      OKLCH{0.6f, 0.12f, 1.0f}, OKLCH{0.7f, 0.08f, 2.0f});
  HS_EXPECT_EQ(input.chroma.basis, ChromaBasis::ABSOLUTE);
  input.chroma.headroom = 0.8f;

  GenerativePalette output;
  PaletteRecipe canonical;
  PaletteCompileStatus status;
  HS_EXPECT_TRUE(
      GenerativePalette::try_compile(input, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::OK);
  HS_EXPECT_EQ(canonical.chroma.headroom, 1.0f);
  const uint64_t headroom_bit =
      uint64_t{1} << static_cast<uint8_t>(PaletteRecipeField::CHROMA_HEADROOM);
  HS_EXPECT_TRUE((status.adjustments.canonicalized_fields & headroom_bit) != 0);
}

inline void test_generative_palette_get_clamps_out_of_range() {
  const GenerativePalette palette(PaletteRecipes::balanced_analogous(0.2f));
  const Pixel low = palette.get(-1.0f).color;
  const Pixel first = palette.get(0.0f).color;
  const Pixel high = palette.get(2.0f).color;
  const Pixel last = palette.get(1.0f).color;
  HS_EXPECT_EQ(low.r, first.r);
  HS_EXPECT_EQ(low.g, first.g);
  HS_EXPECT_EQ(low.b, first.b);
  HS_EXPECT_EQ(high.r, last.r);
  HS_EXPECT_EQ(high.g, last.g);
  HS_EXPECT_EQ(high.b, last.b);
}

inline void test_generative_palette_get_nan_saturates_to_endpoint() {
  const GenerativePalette palette(
      PaletteRecipes::profile(PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC,
                              AxisCurve::CONSTANT, 0.0f, 0.86f));
  const Color4 endpoint = palette.get(1.0f);
  const Color4 singular = palette.get(std::numeric_limits<float>::quiet_NaN());
  HS_EXPECT_EQ(singular.color.r, endpoint.color.r);
  HS_EXPECT_EQ(singular.color.g, endpoint.color.g);
  HS_EXPECT_EQ(singular.color.b, endpoint.color.b);
}

inline void test_generative_palette_morph_compatible() {
  const GenerativePalette a(PaletteRecipes::balanced_analogous(0.0f));
  const GenerativePalette b(PaletteRecipes::balanced_analogous(0.75f));
  HS_EXPECT_TRUE(a.morph_compatible(b));
  HS_EXPECT_TRUE(b.morph_compatible(a));

  const GenerativePalette mirrored(PaletteRecipes::harmony(
      PaletteDomain::MIRROR, PaletteHarmony::ANALOGOUS, 0.0f));
  HS_EXPECT_FALSE(a.morph_compatible(mirrored));

  const GenerativePalette two_keys(PaletteRecipes::harmony(
      PaletteDomain::STRAIGHT, PaletteHarmony::COMPLEMENTARY, 0.0f));
  HS_EXPECT_FALSE(a.morph_compatible(two_keys));

  PaletteRecipe tight = PaletteRecipes::balanced_analogous(0.0f);
  tight.chroma.headroom = 0.8f;
  HS_EXPECT_FALSE(a.morph_compatible(GenerativePalette(tight)));

  PaletteRecipe mixed_lightness = PaletteRecipes::balanced_analogous(0.75f);
  mixed_lightness.lightness.curve = AxisCurve::ASCENDING;
  mixed_lightness.lightness.range = 0.4f;
  HS_EXPECT_FALSE(a.morph_compatible(GenerativePalette(mixed_lightness)));

  PaletteRecipe mixed_chroma = PaletteRecipes::balanced_analogous(0.75f);
  mixed_chroma.chroma.axis.curve = AxisCurve::DESCENDING;
  mixed_chroma.chroma.axis.range = 0.4f;
  HS_EXPECT_FALSE(a.morph_compatible(GenerativePalette(mixed_chroma)));

  const GenerativePalette loop_one(
      PaletteRecipes::isolight_spectral_loop(0.0f));
  const GenerativePalette loop_shifted(
      PaletteRecipes::isolight_spectral_loop(0.4f));
  HS_EXPECT_TRUE(loop_one.morph_compatible(loop_shifted));
  PaletteRecipe two_turn = PaletteRecipes::isolight_spectral_loop(0.0f);
  two_turn.hue.sweep_turns = 2.0f;
  HS_EXPECT_FALSE(loop_one.morph_compatible(GenerativePalette(two_turn)));
}

/** Worst channel gap measured at a 0.001 lerp step from either endpoint. */
constexpr int MEASURED_LERP_ENDPOINT_DELTA = 80;
/** Headroom over the measured gap. */
constexpr int LERP_ENDPOINT_HEADROOM = 2;
/** Endpoint-approach allowance for a 0.001 lerp step, 16-bit scale. */
constexpr int LERP_ENDPOINT_TOLERANCE =
    MEASURED_LERP_ENDPOINT_DELTA * LERP_ENDPOINT_HEADROOM;
/** Endpoint step for the mixed-curve arms, whose axes are further apart. */
constexpr float MIXED_CURVE_ENDPOINT_STEP = 0.0001f;

inline void test_generative_palette_lerp_mixed_curves_continuous() {
  const GenerativePalette from(PaletteRecipes::balanced_analogous(0.1f));
  PaletteRecipe custom_recipe = PaletteRecipes::balanced_analogous(0.1f);
  custom_recipe.lightness.curve = AxisCurve::CUSTOM;
  custom_recipe.lightness.custom[0] = 0.9f;
  custom_recipe.lightness.custom[1] = 0.2f;
  custom_recipe.lightness.custom[2] = 0.9f;
  const GenerativePalette to(custom_recipe);

  alignas(std::max_align_t) static uint8_t
      buf[4 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));
  BakedPaletteStorage from_lut;
  from_lut.bake(arena, from);
  BakedPaletteStorage to_lut;
  to_lut.bake(arena, to);
  BakedPaletteStorage morph_lut;
  morph_lut.bake(arena, from);

  GenerativePalette morph;
  morph.morph_palettes(from, to, 0.999f);
  morph_lut.rebake(morph);
  std::printf("  [lerp-endpoint] to worst=%d\n",
              expect_baked_near(morph_lut, to_lut, LERP_ENDPOINT_TOLERANCE));

  morph.morph_palettes(from, to, 0.001f);
  morph_lut.rebake(morph);
  std::printf("  [lerp-endpoint] from worst=%d\n",
              expect_baked_near(morph_lut, from_lut, LERP_ENDPOINT_TOLERANCE));

  // A snapshot morph between two analytic curves must not snap back to either
  // curve at the endpoints: the interior collapses the mixed pair to CUSTOM.
  const PaletteRecipe bell_recipe =
      PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                              PaletteHarmony::ANALOGOUS, AxisCurve::BELL, 0.1f);
  const GenerativePalette bell(bell_recipe);
  const GenerativePalette::Snapshot bell_keys = bell.snapshot();
  const GenerativePalette::Snapshot ascending_keys =
      GenerativePalette(PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                                PaletteHarmony::ANALOGOUS,
                                                AxisCurve::ASCENDING, 0.1f))
          .snapshot();
  BakedPaletteStorage endpoint_lut;
  endpoint_lut.bake(arena, bell);

  GenerativePalette snapshot_morph(bell_recipe);
  snapshot_morph.morph_snapshots(bell_keys, ascending_keys, 0.0f);
  endpoint_lut.rebake(snapshot_morph);
  snapshot_morph.morph_snapshots(bell_keys, ascending_keys,
                                 MIXED_CURVE_ENDPOINT_STEP);
  morph_lut.rebake(snapshot_morph);
  std::printf(
      "  [lerp-endpoint] mixed-curve from worst=%d\n",
      expect_baked_near(morph_lut, endpoint_lut, LERP_ENDPOINT_TOLERANCE));

  snapshot_morph.morph_snapshots(bell_keys, ascending_keys, 1.0f);
  endpoint_lut.rebake(snapshot_morph);
  snapshot_morph.morph_snapshots(bell_keys, ascending_keys,
                                 1.0f - MIXED_CURVE_ENDPOINT_STEP);
  morph_lut.rebake(snapshot_morph);
  std::printf(
      "  [lerp-endpoint] mixed-curve to worst=%d\n",
      expect_baked_near(morph_lut, endpoint_lut, LERP_ENDPOINT_TOLERANCE));
}

inline void test_generative_palette_lerp_interpolates_loop_seam() {
  const GenerativePalette a(PaletteRecipes::isolight_spectral_loop(0.0f));
  const GenerativePalette b(PaletteRecipes::isolight_spectral_loop(0.4f));
  GenerativePalette morph;
  for (const float amount : {0.25f, 0.5f, 0.75f}) {
    morph.morph_palettes(a, b, amount);
    const Pixel seam = morph.get(1.0f).color;
    const Pixel start = morph.get(0.0f).color;
    HS_EXPECT_TRUE(std::abs(int(seam.r) - int(start.r)) <= 220);
    HS_EXPECT_TRUE(std::abs(int(seam.g) - int(start.g)) <= 220);
    HS_EXPECT_TRUE(std::abs(int(seam.b) - int(start.b)) <= 220);
  }
}
