/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file effect_palette_recipes.h
 * @brief The per-effect PaletteRecipe builders, plus the preset roster the
 *        palette authoring tool opens with.
 */

#include <array>

#include "color/color.h"

/** @brief PaletteRecipe builders owned by the effects that render them. */
namespace EffectPaletteRecipes {

/** @brief The BZ reaction-diffusion ramp, pinned to three authored colors. */
HS_FLASH_MEMBER inline PaletteRecipe bz_reaction_diffusion() {
  return PaletteRecipes::from_oklch_keys(PaletteDomain::STRAIGHT,
                                         pixel_to_oklch(Pixel(36844, 10770, 3)),
                                         pixel_to_oklch(Pixel(0, 8112, 5753)),
                                         pixel_to_oklch(Pixel(2059, 0, 9668)));
}

/**
 * @brief The Comets ramp: a triadic harmony brightening along the domain.
 * @param base_turns Base hue in turns.
 * @return The recipe.
 */
HS_FLASH_MEMBER inline PaletteRecipe comets(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                 PaletteHarmony::TRIADIC, AxisCurve::ASCENDING,
                                 base_turns);
}

/**
 * @brief The DisplacementField / RingShower ramp: an analogous mirror at
 *        constant lightness.
 * @param base_turns Base hue in turns.
 * @return The recipe.
 */
HS_FLASH_MEMBER inline PaletteRecipe displacement_field(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::MIRROR,
                                 PaletteHarmony::ANALOGOUS, AxisCurve::CONSTANT,
                                 base_turns);
}

/**
 * @brief The Dynamo ramp: an analogous vignette brightening along the domain.
 * @param base_turns Base hue in turns.
 * @return The recipe.
 */
HS_FLASH_MEMBER inline PaletteRecipe dynamo(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::VIGNETTE,
                                 PaletteHarmony::ANALOGOUS,
                                 AxisCurve::ASCENDING, base_turns);
}

/**
 * @brief The Gray-Scott reaction-diffusion ramp: a split-complementary
 *        harmony brightening along the domain.
 * @param base_turns Base hue in turns; defaults to the authored hue.
 * @return The recipe.
 */
HS_FLASH_MEMBER inline PaletteRecipe
gs_reaction_diffusion(float base_turns = PaletteRecipes::hue_turns(160)) {
  return PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                 PaletteHarmony::SPLIT_COMPLEMENTARY,
                                 AxisCurve::ASCENDING, base_turns, 0.50f);
}

/**
 * @brief The MindSplatter trail ramp: a triadic harmony darkening along the
 *        domain at full local-gamut chroma.
 * @param base_turns Base hue in turns.
 * @return The recipe.
 * @details The recipe of record for the baked bank in
 *          core/color/triadic_palette_luts.h.
 */
HS_FLASH_MEMBER inline PaletteRecipe mind_splatter(float base_turns) {
  PaletteRecipe recipe;
  recipe.schema_version = PaletteRecipe::SCHEMA_VERSION;
  recipe.input.offset = 0.0f;
  recipe.input.span = 0.624f;
  recipe.domain = PaletteDomain::STRAIGHT;
  recipe.easing = SegmentEase::COSINE;
  recipe.color_path = ColorPath::OKLAB_CARTESIAN;
  recipe.hue.mode = HueMode::HARMONY;
  recipe.hue.harmony = PaletteHarmony::TRIADIC;
  recipe.hue.direction = HueDirection::CLOCKWISE;
  recipe.hue.base_turns = base_turns;
  recipe.hue.spread_turns = 0.07f;
  recipe.hue.sweep_turns = 1.0f;
  recipe.hue.custom_turns = {0.0f, 0.0f, 0.0f, 0.0f};
  recipe.lightness.curve = AxisCurve::DESCENDING;
  recipe.lightness.center = 0.37f;
  recipe.lightness.range = 0.5f;
  recipe.lightness.custom = {0.0f, 0.0f, 0.0f, 0.0f};
  recipe.chroma.curve = AxisCurve::CONSTANT;
  recipe.chroma.basis = ChromaBasis::LOCAL_GAMUT;
  recipe.chroma.center = 0.95f;
  recipe.chroma.range = 0.0f;
  recipe.chroma.headroom = 0.94f;
  recipe.chroma.custom = {0.0f, 0.0f, 0.0f, 0.0f};
  recipe.hue_torsion = 0.0f;
  recipe.falloff_start = 0.9f;
  return recipe;
}

/**
 * @brief Base hue of one entry of the baked MindSplatter bank.
 * @param index Entry index in [0, count).
 * @param count Bank size.
 * @return The base hue in turns.
 */
HS_FLASH_MEMBER inline float mind_splatter_bank_turns(int index, int count) {
  return 0.08611111f + static_cast<float>(index) / static_cast<float>(count);
}

/**
 * @brief The MobiusRings ramp: a split-complementary mirror at constant
 *        lightness.
 * @param base_turns Base hue in turns.
 * @return The recipe.
 */
HS_FLASH_MEMBER inline PaletteRecipe mobius_rings(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::MIRROR,
                                 PaletteHarmony::SPLIT_COMPLEMENTARY,
                                 AxisCurve::CONSTANT, base_turns);
}

/** @brief The Raymarch ramp: a complementary pair peaking mid-domain. */
HS_FLASH_MEMBER inline PaletteRecipe raymarch() {
  return PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                 PaletteHarmony::COMPLEMENTARY, AxisCurve::BELL,
                                 PaletteRecipes::hue_turns(219), 0.86f);
}

/** @brief The liquid recipe at an arbitrary hue rotation; every rotation is
 *  morph-compatible with every other
 *  (see test_shader_ball_palette_rotations_morph_compatible).
 *  @param rotation_turns Hue rotation in turns applied to every key.
 *  @return The recipe.
 *  @details Palindromic keys: hue travels a half turn out to the complement
 *  and back while lightness dives bright-dark-bright, so the shader's wrapped
 *  palette coordinate crosses the 1 -> 0 seam without a hard line. */
HS_FLASH_MEMBER inline PaletteRecipe
shader_ball_liquid_at(float rotation_turns) {
  constexpr float BASE_TURNS = 0.2933125f;
  PaletteRecipe recipe;
  recipe.domain = PaletteDomain::STRAIGHT;
  recipe.hue.mode = HueMode::CUSTOM;
  recipe.hue.custom_turns[0] = BASE_TURNS + rotation_turns;
  recipe.hue.custom_turns[1] = BASE_TURNS + rotation_turns + 0.5f;
  recipe.hue.custom_turns[2] = BASE_TURNS + rotation_turns;
  recipe.lightness.curve = AxisCurve::CUSTOM;
  recipe.lightness.custom[0] = 0.8798438f;
  recipe.lightness.custom[1] = 0.1623438f;
  recipe.lightness.custom[2] = 0.8798438f;
  recipe.chroma.center = 0.8871875f;
  return recipe;
}

/** @brief The liquid recipe at its authored hue. */
HS_FLASH_MEMBER inline PaletteRecipe shader_ball_liquid() {
  return shader_ball_liquid_at(0.0f);
}

/** @brief The flyby recipe at its authored hue. */
HS_FLASH_MEMBER inline PaletteRecipe shader_ball_flyby() {
  return PaletteRecipes::profile(
      PaletteDomain::STRAIGHT, PaletteHarmony::SPLIT_COMPLEMENTARY,
      AxisCurve::CONSTANT, PaletteRecipes::hue_turns(42));
}

/**
 * @brief The flyby split-complementary profile at an arbitrary base hue.
 * @param base_turns Base hue in turns.
 * @return The recipe.
 */
HS_FLASH_MEMBER inline PaletteRecipe shader_ball_flyby_at(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                 PaletteHarmony::SPLIT_COMPLEMENTARY,
                                 AxisCurve::CONSTANT, base_turns);
}

/** @brief One row of the authoring tool's preset roster. */
struct Preset {
  const char *name; /**< Display name. */
  /** @brief Whether the effect randomizes this recipe's base hue at run time. */
  bool random_hue;
  PaletteRecipe recipe; /**< The recipe at its preview hue. */
};

/** @brief A base hue drawn from the 256-step hue wheel. */
HS_FLASH_MEMBER inline float random_base_turns() {
  return PaletteRecipes::hue_turns(static_cast<uint32_t>(hs::rand_int(0, 256)));
}

/** @brief The preset roster, every recipe at a fixed preview hue. */
HS_FLASH_MEMBER inline std::array<Preset, 9> presets() {
  const float preview_hue = PaletteRecipes::hue_turns(42);
  return {{{"BZReactionDiffusion", false, bz_reaction_diffusion()},
           {"Comets", true, comets(preview_hue)},
           {"DisplacementField / RingShower", true,
            displacement_field(preview_hue)},
           {"Dynamo", true, dynamo(preview_hue)},
           {"GSReactionDiffusion", true, gs_reaction_diffusion(preview_hue)},
           {"MobiusRings", true, mobius_rings(preview_hue)},
           {"Raymarch", false, raymarch()},
           {"ShaderBall Liquid", false, shader_ball_liquid()},
           {"ShaderBall Flyby", false, shader_ball_flyby()}}};
}

} // namespace EffectPaletteRecipes
