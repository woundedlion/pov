/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/color/generative_palette.h and
 * core/color/palette_cycler.h. Included by tests/test_color.h, whose
 * run_color_tests() calls these cases and whose roster row counts them.
 */
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>

#include "core/color/color.h"
#include "core/color/composition.h"
#include "core/color/effect_palette_recipes.h"
#include "core/color/generative_palette.h"
#include "core/color/palette_cycler.h"
#include "core/engine/memory.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace color_tests {

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

inline void test_effect_palette_recipe_roster() {
  const auto presets = EffectPaletteRecipes::presets();
  HS_EXPECT_EQ(presets.size(), size_t{9});
  HS_EXPECT_EQ(std::strcmp(presets[0].name, "BZReactionDiffusion"), 0);
  HS_EXPECT_EQ(std::strcmp(presets[2].name, "DisplacementField / RingShower"),
               0);
  HS_EXPECT_EQ(std::strcmp(presets[8].name, "ShaderBall Flyby"), 0);
  HS_EXPECT_FALSE(presets[0].random_hue);
  HS_EXPECT_TRUE(presets[1].random_hue);
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

inline void test_generative_palette_lerp_target_aliases_this() {
  const GenerativePalette from(PaletteRecipes::balanced_analogous(0.0f));
  const GenerativePalette to(PaletteRecipes::balanced_analogous(0.25f));
  GenerativePalette expected;
  expected.lerp(from, to, 0.25f);
  GenerativePalette aliased = to;
  aliased.lerp(from, aliased, 0.25f);
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
  morph.lerp(morph.snapshot(), target.snapshot(), 1.0f);
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

namespace {

/** @brief Asserts two baked LUTs agree at every sample point. */
inline void expect_baked_equal(const BakedPalette &a, const BakedPalette &b) {
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
    const float t = i / 255.0f;
    const Color4 ca = a.get(t);
    const Color4 cb = b.get(t);
    HS_EXPECT_EQ(ca.color.r, cb.color.r);
    HS_EXPECT_EQ(ca.color.g, cb.color.g);
    HS_EXPECT_EQ(ca.color.b, cb.color.b);
    HS_EXPECT_EQ(ca.alpha, cb.alpha);
  }
}

/** @brief Asserts two baked LUTs agree within a channel tolerance. */
inline void expect_baked_near(const BakedPalette &a, const BakedPalette &b,
                              int tolerance) {
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
    const float t = i / 255.0f;
    const Pixel pa = a.get(t).color;
    const Pixel pb = b.get(t).color;
    HS_EXPECT_TRUE(std::abs(int(pa.r) - int(pb.r)) <= tolerance);
    HS_EXPECT_TRUE(std::abs(int(pa.g) - int(pb.g)) <= tolerance);
    HS_EXPECT_TRUE(std::abs(int(pa.b) - int(pb.b)) <= tolerance);
  }
}

} // namespace

inline void test_baked_palette_rebake_crossfade() {
  Gradient ramp{{0.0f, CPixel(0u, 0u, 0u)}, {1.0f, CPixel(255u, 255u, 255u)}};
  Gradient warm{{0.0f, CPixel(200u, 40u, 10u)}, {1.0f, CPixel(20u, 80u, 220u)}};

  alignas(std::max_align_t) static uint8_t
      buf[4 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  BakedPalette from;
  from.bake(arena, ramp);
  BakedPalette to;
  to.bake(arena, warm);
  BakedPalette out;
  out.bake(arena, ramp);

  out.rebake_crossfade(from, to, 0.0f);
  expect_baked_equal(out, from);
  out.rebake_crossfade(from, to, 1.0f);
  expect_baked_equal(out, to);

  BakedPalette expected;
  expected.bake_blend(arena, from, to, 0.37f);
  out.rebake_crossfade(from, to, 0.37f);
  expect_baked_equal(out, expected);

  out.rebake_crossfade(from, to, std::numeric_limits<float>::quiet_NaN());
  expect_baked_equal(out, to);
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
  mixed_chroma.chroma.curve = AxisCurve::DESCENDING;
  mixed_chroma.chroma.range = 0.4f;
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
  BakedPalette from_lut;
  from_lut.bake(arena, from);
  BakedPalette to_lut;
  to_lut.bake(arena, to);
  BakedPalette morph_lut;
  morph_lut.bake(arena, from);

  GenerativePalette morph;
  morph.lerp(from, to, 0.999f);
  morph_lut.rebake(morph);
  expect_baked_near(morph_lut, to_lut, 2500);

  morph.lerp(from, to, 0.001f);
  morph_lut.rebake(morph);
  expect_baked_near(morph_lut, from_lut, 2500);
}

inline void test_generative_palette_lerp_interpolates_loop_seam() {
  const GenerativePalette a(PaletteRecipes::isolight_spectral_loop(0.0f));
  const GenerativePalette b(PaletteRecipes::isolight_spectral_loop(0.4f));
  GenerativePalette morph;
  for (const float amount : {0.25f, 0.5f, 0.75f}) {
    morph.lerp(a, b, amount);
    const Pixel seam = morph.get(1.0f).color;
    const Pixel start = morph.get(0.0f).color;
    HS_EXPECT_TRUE(std::abs(int(seam.r) - int(start.r)) <= 220);
    HS_EXPECT_TRUE(std::abs(int(seam.g) - int(start.g)) <= 220);
    HS_EXPECT_TRUE(std::abs(int(seam.b) - int(start.b)) <= 220);
  }
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

  PaletteCycler cycler;
  cycler.init(cycler_arena, entries.data(), entries.size(), 3, 4);

  BakedPalette ref;
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
  expected_morph.lerp(first, second, 0.25f);
  BakedPalette expected;
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

  BakedPalette prebaked;
  prebaked.bake(arena, solid);

  const std::array<PaletteCycler::Entry, 3> entries = {
      {ramp, prebaked, generative}};

  PaletteCycler cycler;
  cycler.init(arena, entries.data(), entries.size(), 1, 2);

  BakedPalette ramp_lut;
  ramp_lut.bake(arena, ramp);
  expect_baked_equal(cycler.palette(), ramp_lut);

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());
  cycler.step();
  BakedPalette expected;
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
  BakedPalette generative_lut;
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
  BakedPalette ref;
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

inline void test_shader_ball_palette_rotations_morph_compatible() {
  constexpr float GOLDEN_STEP = 0.618034f;
  float rotation = 0.0f;
  GenerativePalette liquid_prev(
      EffectPaletteRecipes::shader_ball_liquid_at(0.0f));
  GenerativePalette flyby_prev(
      EffectPaletteRecipes::shader_ball_flyby_at(0.0f));
  for (int i = 0; i < 24; ++i) {
    rotation = wrap_t(rotation + GOLDEN_STEP);
    GenerativePalette liquid(
        EffectPaletteRecipes::shader_ball_liquid_at(rotation));
    GenerativePalette flyby(
        EffectPaletteRecipes::shader_ball_flyby_at(rotation));
    HS_EXPECT_TRUE(liquid_prev.morph_compatible(liquid));
    HS_EXPECT_TRUE(flyby_prev.morph_compatible(flyby));
    liquid_prev = liquid;
    flyby_prev = flyby;
  }
}

namespace {

/** @brief Deterministic PaletteCycler provider walking analogous base hues. */
inline void scripted_next_palette(void *context, uint32_t sequence,
                                  GenerativePalette &out) {
  ++*static_cast<int *>(context);
  out = GenerativePalette(
      PaletteRecipes::balanced_analogous(0.1f + 0.2f * sequence));
}

} // namespace

inline void test_palette_cycler_generated_cycle() {
  alignas(std::max_align_t) static uint8_t
      buf[PaletteCycler::generated_arena_bytes() +
          3 * BakedPalette::required_arena_bytes()];
  Arena arena(buf, sizeof(buf));

  int provider_calls = 0;
  PaletteCycler cycler;
  cycler.init_generated(arena, scripted_next_palette, &provider_calls, 0, 2);
  HS_EXPECT_EQ(provider_calls, 2);

  BakedPalette ref;
  ref.bake(arena, GenerativePalette(PaletteRecipes::balanced_analogous(0.1f)));
  expect_baked_equal(cycler.palette(), ref);

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());

  cycler.step();
  GenerativePalette mid;
  mid.lerp(GenerativePalette(PaletteRecipes::balanced_analogous(0.1f)),
           GenerativePalette(PaletteRecipes::balanced_analogous(0.3f)), 0.5f);
  BakedPalette expected;
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
  blended.lerp(first, second, 0.5f);
  BakedPalette blended_lut;
  blended_lut.bake(reinit_arena, blended);
  expect_baked_equal(cycler.palette(), blended_lut);
  HS_EXPECT_EQ(provider_calls, 3);
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

  BakedPalette ref;
  ref.bake(arena, a);
  expect_baked_equal(cycler.palette(), ref);
  HS_EXPECT_FALSE(cycler.fading());

  cycler.step();
  HS_EXPECT_TRUE(cycler.fading());

  cycler.step();
  GenerativePalette mid;
  mid.lerp(a, b, 0.5f);
  BakedPalette expected;
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

} // namespace color_tests
} // namespace hs_test
