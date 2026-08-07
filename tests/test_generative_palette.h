/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

inline void test_generative_palette_deterministic() {
  const PaletteRecipe recipe = PaletteRecipes::profile(
      PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC, AxisCurve::BELL,
      PaletteRecipes::hue_turns(77), 0.86f);
  const GenerativePalette first(recipe);
  const GenerativePalette second(recipe);
  for (int i = 0; i < 256; ++i) {
    const Pixel a = first.get(i / 255.0f).color;
    const Pixel b = second.get(i / 255.0f).color;
    HS_EXPECT_EQ(a.r, b.r);
    HS_EXPECT_EQ(a.g, b.g);
    HS_EXPECT_EQ(a.b, b.b);
  }
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

  const GenerativePalette::Snapshot before = output.snapshot();
  const PaletteRecipe canonical_before = canonical;
  input.lightness.center = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(
      GenerativePalette::try_compile(input, output, canonical, status));
  HS_EXPECT_EQ(status.code, PaletteCompileCode::NON_FINITE);
  HS_EXPECT_EQ(status.field, PaletteRecipeField::LIGHTNESS_CENTER);
  const GenerativePalette::Snapshot after = output.snapshot();
  HS_EXPECT_EQ(std::memcmp(&before, &after, sizeof(before)), 0);
  HS_EXPECT_EQ(
      std::memcmp(&canonical_before, &canonical, sizeof(canonical_before)), 0);

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
  const GenerativePalette windowed_envelope(recipe);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(0.0f).L, 0.2f, 1e-5f);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(1.0f).L, 0.8f, 1e-5f);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(0.0f).h_path,
                 full.diagnose(0.2f).h_path, 1e-5f);
  HS_EXPECT_NEAR(windowed_envelope.diagnose(1.0f).h_path,
                 full.diagnose(0.6f).h_path, 1e-5f);

  recipe.domain = PaletteDomain::MIRROR;
  const GenerativePalette cropped_mirror(recipe);
  HS_EXPECT_FALSE(cropped_mirror.mirrors_domain());
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
  HS_EXPECT_NEAR(key1.h - key0.h, PI_F, 1e-5f);

  recipe.lightness.curve = AxisCurve::BELL;
  const GenerativePalette bell(recipe);
  HS_EXPECT_NEAR(bell.diagnose(0.0f).L, 0.2f, 1e-5f);
  HS_EXPECT_NEAR(bell.diagnose(0.5f).L, 0.8f, 1e-5f);
  HS_EXPECT_NEAR(bell.diagnose(1.0f).L, 0.2f, 1e-5f);

  recipe.hue.direction = HueDirection::CLOCKWISE;
  const auto clockwise = GenerativePalette(recipe).snapshot();
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(clockwise, 1).h -
                     GenerativePalette::snapshot_key(clockwise, 0).h,
                 -PI_F, 1e-5f);

  recipe.hue.direction = HueDirection::COUNTERCLOCKWISE;
  recipe.hue.harmony = PaletteHarmony::TETRADIC;
  recipe.hue.spread_turns = 1.0f / 6.0f;
  const auto tetradic = GenerativePalette(recipe).snapshot();
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(tetradic, 1).h -
                     GenerativePalette::snapshot_key(tetradic, 0).h,
                 PI_F / 3.0f, 1e-5f);
  HS_EXPECT_NEAR(GenerativePalette::snapshot_key(tetradic, 2).h -
                     GenerativePalette::snapshot_key(tetradic, 0).h,
                 PI_F, 1e-5f);

  recipe.hue.harmony = PaletteHarmony::SQUARE;
  const auto square = GenerativePalette(recipe).snapshot();
  for (int i = 1; i < 4; ++i)
    HS_EXPECT_NEAR(GenerativePalette::snapshot_key(square, i).h -
                       GenerativePalette::snapshot_key(square, i - 1).h,
                   0.5f * PI_F, 1e-5f);
}

inline void test_generative_palette_blue_cusp_is_continuous() {
  PaletteRecipe recipe;
  recipe.hue.mode = HueMode::SWEEP;
  recipe.hue.base_turns = 98.0f / 256.0f;
  recipe.hue.sweep_turns = 1.0f;
  recipe.chroma.center = 1.0f;
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
        std::max(largest_chroma_step, fabsf(current.C - previous.C));
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
  BakedPalette baked;
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

inline void test_generative_palette_snapshot_lerp() {
  GenerativePalette from(PaletteRecipes::balanced_analogous(0.0f));
  const GenerativePalette to(PaletteRecipes::balanced_analogous(0.75f));
  const auto first = from.snapshot();
  const auto last = to.snapshot();

  from.lerp(first, last, 0.5f);
  const auto middle = from.snapshot();
  const auto first_key = GenerativePalette::snapshot_key(first, 0);
  const auto last_key = GenerativePalette::snapshot_key(last, 0);
  const auto middle_key = GenerativePalette::snapshot_key(middle, 0);
  HS_EXPECT_NEAR(middle_key.L, 0.5f * (first_key.L + last_key.L), 3e-4f);
  HS_EXPECT_NEAR(middle_key.q, 0.5f * (first_key.q + last_key.q), 3e-4f);

  from.lerp(first, last, 1.0f);
  const GenerativePalette::Snapshot target = from.snapshot();
  HS_EXPECT_EQ(std::memcmp(&target, &last, sizeof(last)), 0);
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
  HS_EXPECT_EQ(status.code, PaletteCompileCode::INCOMPATIBLE_OPTIONS);
  HS_EXPECT_EQ(status.field, PaletteRecipeField::CHROMA_BASIS);
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

inline void test_mobius_longitude_singularity_saturates_to_endpoint() {
  const float z = 1.0f;
  const float R = std::sqrt((1.0f + z) / (1.0f - z));
  HS_EXPECT_TRUE(std::isinf(R));

  const float t = (std::log(R) + 2.5f) / 5.0f;
  const float wrapped = wrap(t, 1.0f);
  HS_EXPECT_TRUE(std::isnan(wrapped));

  const GenerativePalette palette(
      PaletteRecipes::profile(PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC,
                              AxisCurve::CONSTANT, 0.0f, 0.86f));
  const Color4 endpoint = palette.get(1.0f);
  const Color4 singular = palette.get(wrapped);
  HS_EXPECT_EQ(singular.color.r, endpoint.color.r);
  HS_EXPECT_EQ(singular.color.g, endpoint.color.g);
  HS_EXPECT_EQ(singular.color.b, endpoint.color.b);
}
