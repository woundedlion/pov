/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file authored_palette.h
 * @brief Authored source cache behavior. */

#include "core/color/authored_palette.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test::color_tests {

/** @brief Inactive source controls preserve every sample without rebaking. */
inline void test_authored_palette_inactive_fields_preserve_bake() {
  using Type = AuthoredPaletteType;
  for (Type type :
       {Type::GENERATED, Type::PROCEDURAL, Type::GRADIENT, Type::MUTATING}) {
    reset_globals();
    AuthoredPaletteResources resources;
    AuthoredPaletteSelection params;
    params.palette_type = static_cast<uint8_t>(type);
    const auto &palette =
        resources.prepare(persistent_arena, params, PALETTE_DEFAULT_CHROMA);
    std::array<Color4, 256> before;
    for (size_t index = 0; index < before.size(); ++index)
      before[index] = palette.get(float(index) / 255.0f);
    HS_EXPECT_EQ(resources.bake_generation(), 1u);
    for (int frame = 1; frame <= 60; ++frame) {
      const float t = float(frame) / 60.0f;
      if (type == Type::GENERATED || type == Type::GRADIENT) {
        params.palette_name = static_cast<uint8_t>(frame % 20);
        params.palette_name_b = static_cast<uint8_t>((frame + 7) % 20);
        params.palette_mutation = t;
      }
      if (type != Type::GENERATED) {
        params.palette_hue = t;
        params.palette_spread = t * 0.2f;
        params.palette_lightness = t;
        params.palette_lightness_range = t;
        params.palette_harmony = static_cast<uint8_t>(frame % 8);
        params.palette_lightness_curve = static_cast<uint8_t>(frame % 5);
      }
      if (type != Type::GRADIENT) {
        params.gradient_start_r = t;
        params.gradient_middle_g = t;
        params.gradient_end_b = t;
      }
      if (type == Type::PROCEDURAL) {
        params.palette_name_b = static_cast<uint8_t>(frame % 20);
        params.palette_mutation = t;
      }
      params.reserved[0] = static_cast<uint8_t>(frame);
      const auto &unchanged =
          resources.prepare(persistent_arena, params, PALETTE_DEFAULT_CHROMA);
      HS_EXPECT_EQ(resources.bake_generation(), 1u);
      for (size_t index = 0; index < before.size(); ++index) {
        const Color4 actual = unchanged.get(float(index) / 255.0f);
        HS_EXPECT_EQ(actual.color.r, before[index].color.r);
        HS_EXPECT_EQ(actual.color.g, before[index].color.g);
        HS_EXPECT_EQ(actual.color.b, before[index].color.b);
        HS_EXPECT_EQ(actual.alpha, before[index].alpha);
      }
    }
    const auto expect_rebake = [&](AuthoredPaletteSelection changed) {
      const auto generation = resources.bake_generation();
      (void)resources.prepare(persistent_arena, changed,
                              PALETTE_DEFAULT_CHROMA);
      HS_EXPECT_EQ(resources.bake_generation(), generation + 1);
      (void)resources.prepare(persistent_arena, params, PALETTE_DEFAULT_CHROMA);
      HS_EXPECT_EQ(resources.bake_generation(), generation + 2);
    };
    const auto check_float = [&](float AuthoredPaletteSelection::*field) {
      auto changed = params;
      changed.*field = params.*field < 0.5f ? 0.75f : 0.25f;
      expect_rebake(changed);
    };
    const auto check_index = [&](uint8_t AuthoredPaletteSelection::*field,
                                 size_t count) {
      auto changed = params;
      changed.*field = static_cast<uint8_t>((params.*field + 1) % count);
      expect_rebake(changed);
    };
    if (type == Type::GENERATED) {
      for (auto field : {&AuthoredPaletteSelection::palette_hue,
                         &AuthoredPaletteSelection::palette_spread,
                         &AuthoredPaletteSelection::palette_lightness,
                         &AuthoredPaletteSelection::palette_lightness_range})
        check_float(field);
      check_index(&AuthoredPaletteSelection::palette_harmony,
                  std::size(AUTHORED_PALETTE_HARMONY_IDS));
      check_index(&AuthoredPaletteSelection::palette_lightness_curve,
                  std::size(AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS));
    } else if (type == Type::GRADIENT) {
      for (auto field : {&AuthoredPaletteSelection::gradient_start_r,
                         &AuthoredPaletteSelection::gradient_start_g,
                         &AuthoredPaletteSelection::gradient_start_b,
                         &AuthoredPaletteSelection::gradient_middle_r,
                         &AuthoredPaletteSelection::gradient_middle_g,
                         &AuthoredPaletteSelection::gradient_middle_b,
                         &AuthoredPaletteSelection::gradient_end_r,
                         &AuthoredPaletteSelection::gradient_end_g,
                         &AuthoredPaletteSelection::gradient_end_b})
        check_float(field);
    } else {
      check_index(&AuthoredPaletteSelection::palette_name,
                  AuthoredPaletteResources::named_palette_count());
      if (type == Type::MUTATING) {
        check_index(&AuthoredPaletteSelection::palette_name_b,
                    AuthoredPaletteResources::named_palette_count());
        check_float(&AuthoredPaletteSelection::palette_mutation);
      }
    }
    const auto generation = resources.bake_generation();
    (void)resources.prepare(persistent_arena, params,
                            PALETTE_DEFAULT_CHROMA * 0.5f);
    HS_EXPECT_EQ(resources.bake_generation(), generation + 1);
  }
}

} // namespace hs_test::color_tests
