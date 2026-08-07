/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>

#include "core/color/generative_palette.h"

namespace EffectPaletteRecipes {

HS_FLASH_MEMBER inline PaletteRecipe bz_reaction_diffusion() {
  return PaletteRecipes::from_oklch_keys(PaletteDomain::STRAIGHT,
                                         pixel_to_oklch(Pixel(36844, 10770, 3)),
                                         pixel_to_oklch(Pixel(0, 8112, 5753)),
                                         pixel_to_oklch(Pixel(2059, 0, 9668)));
}

HS_FLASH_MEMBER inline PaletteRecipe comets(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                 PaletteHarmony::TRIADIC, AxisCurve::ASCENDING,
                                 base_turns);
}

HS_FLASH_MEMBER inline PaletteRecipe displacement_field(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::MIRROR,
                                 PaletteHarmony::ANALOGOUS, AxisCurve::CONSTANT,
                                 base_turns);
}

HS_FLASH_MEMBER inline PaletteRecipe dynamo(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::VIGNETTE,
                                 PaletteHarmony::ANALOGOUS,
                                 AxisCurve::ASCENDING, base_turns);
}

HS_FLASH_MEMBER inline PaletteRecipe gs_reaction_diffusion() {
  return PaletteRecipes::profile(
      PaletteDomain::STRAIGHT, PaletteHarmony::SPLIT_COMPLEMENTARY,
      AxisCurve::ASCENDING, PaletteRecipes::hue_turns(160), 0.50f);
}

HS_FLASH_MEMBER inline PaletteRecipe mobius_grid(float base_turns) {
  return PaletteRecipes::profile(PaletteDomain::MIRROR,
                                 PaletteHarmony::SPLIT_COMPLEMENTARY,
                                 AxisCurve::CONSTANT, base_turns);
}

HS_FLASH_MEMBER inline PaletteRecipe raymarch() {
  return PaletteRecipes::profile(PaletteDomain::STRAIGHT,
                                 PaletteHarmony::COMPLEMENTARY, AxisCurve::BELL,
                                 PaletteRecipes::hue_turns(219), 0.86f);
}

HS_FLASH_MEMBER inline PaletteRecipe shader_ball_liquid() {
  constexpr float BASE_TURNS = 0.2933125f;
  PaletteRecipe recipe;
  recipe.domain = PaletteDomain::STRAIGHT;
  recipe.hue.mode = HueMode::CUSTOM;
  recipe.hue.custom_turns[0] = BASE_TURNS;
  recipe.hue.custom_turns[1] = BASE_TURNS + 0.5f;
  recipe.hue.custom_turns[2] = BASE_TURNS;
  recipe.lightness.curve = AxisCurve::CUSTOM;
  recipe.lightness.custom[0] = 0.8798438f;
  recipe.lightness.custom[1] = 0.1623438f;
  recipe.lightness.custom[2] = 0.8798438f;
  recipe.chroma.center = 0.8871875f;
  return recipe;
}

HS_FLASH_MEMBER inline PaletteRecipe shader_ball_flyby() {
  return PaletteRecipes::profile(
      PaletteDomain::STRAIGHT, PaletteHarmony::SPLIT_COMPLEMENTARY,
      AxisCurve::CONSTANT, PaletteRecipes::hue_turns(42));
}

struct Preset {
  const char *name;
  bool random_hue;
  PaletteRecipe recipe;
};

HS_FLASH_MEMBER inline float random_base_turns() {
  return PaletteRecipes::hue_turns(static_cast<uint32_t>(hs::rand_int(0, 256)));
}

inline std::array<Preset, 9> presets() {
  const float preview_hue = PaletteRecipes::hue_turns(42);
  return {{{"BZReactionDiffusion", false, bz_reaction_diffusion()},
           {"Comets", true, comets(preview_hue)},
           {"DisplacementField / RingShower", true,
            displacement_field(preview_hue)},
           {"Dynamo", true, dynamo(preview_hue)},
           {"GSReactionDiffusion", false, gs_reaction_diffusion()},
           {"MobiusGrid", true, mobius_grid(preview_hue)},
           {"Raymarch", false, raymarch()},
           {"ShaderBall Liquid", false, shader_ball_liquid()},
           {"ShaderBall Flyby", false, shader_ball_flyby()}}};
}

} // namespace EffectPaletteRecipes
