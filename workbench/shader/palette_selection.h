/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette_selection.h
 * @brief Cached workbench palette-source selection and baking.
 */

#include "core/color/palettes.h"
#include "core/render/pullback/operators.h"

/** @brief Baked authored palette with a cache keyed by its source controls. */
class ShaderPaletteResources {
public:
  /** @brief Resolves authored controls into an owned table.
   * @param arena Persistent table allocation arena.
   * @param params Current color operator parameters.
   * @return Palette view, valid until the next source change. */
  const BakedPalette &
  prepare(Arena &arena,
          const Pullback::Interp::Op::GeneratedPaletteParams &params) {
    const auto &selection =
        static_cast<const Pullback::Interp::Op::PaletteSelection &>(params);
    if (!initialized || !(selection == previous) ||
        params.palette_chroma != previous_chroma) {
      const auto bake = [&](const auto &source) {
        if (initialized)
          palette.rebake(source);
        else
          palette.bake(arena, source);
      };
      using Type = Pullback::Interp::Op::PaletteType;
      const auto type = static_cast<Type>(params.palette_type);
      if (type == Type::GENERATED) {
        PaletteRecipe recipe;
        recipe.hue.harmony =
            static_cast<PaletteHarmony>(params.palette_harmony);
        recipe.hue.base_turns = params.palette_hue;
        recipe.hue.spread_turns = params.palette_spread;
        recipe.lightness.curve =
            static_cast<AxisCurve>(params.palette_lightness_curve);
        recipe.lightness.center = params.palette_lightness;
        recipe.lightness.range = params.palette_lightness_range;
        recipe.chroma.axis.center = params.palette_chroma;
        const GenerativePalette source(recipe);
        bake(source);
      } else if (type == Type::GRADIENT) {
        const auto color = [](float r, float g, float b) {
          return CPixel(static_cast<uint8_t>(r * 255.0f + 0.5f),
                        static_cast<uint8_t>(g * 255.0f + 0.5f),
                        static_cast<uint8_t>(b * 255.0f + 0.5f));
        };
        const Gradient source{
            {0.0f, color(params.gradient_start_r, params.gradient_start_g,
                         params.gradient_start_b)},
            {0.5f, color(params.gradient_middle_r, params.gradient_middle_g,
                         params.gradient_middle_b)},
            {1.0f, color(params.gradient_end_r, params.gradient_end_g,
                         params.gradient_end_b)}};
        bake(ChromaSource{source, params.palette_chroma});
      } else {
        HS_CHECK(params.palette_name < std::size(COEFFICIENTS) &&
                     params.palette_name_b < std::size(COEFFICIENTS),
                 "invalid named palette index");
        const auto &a = COEFFICIENTS[params.palette_name];
        if (type == Type::MUTATING) {
          const auto &b = COEFFICIENTS[params.palette_name_b];
          MutatingPalette source(a.a, a.b, a.c, a.d, b.a, b.b, b.c, b.d);
          source.mutate(params.palette_mutation);
          bake(ChromaSource{source, params.palette_chroma});
        } else {
          HS_CHECK(type == Type::PROCEDURAL, "invalid authored palette type");
          const ProceduralPalette source(a.a, a.b, a.c, a.d);
          bake(ChromaSource{source, params.palette_chroma});
        }
      }
      previous = selection;
      previous_chroma = params.palette_chroma;
      initialized = true;
    }
    return palette.view();
  }

private:
  /** @brief Named cosine coefficients in roster order. */
  struct Coefficients {
    std::array<float, 3> a; ///< Bias.
    std::array<float, 3> b; ///< Amplitude.
    std::array<float, 3> c; ///< Frequency.
    std::array<float, 3> d; ///< Phase.
  };
  /** @brief Applies the common chroma control around the palette default. */
  struct ChromaSource {
    const Palette &source; ///< Source sampled before adjustment.
    float chroma;          ///< Common chroma value.
    /** @brief Samples with perceptual chroma scaling.
     * @param t Palette coordinate.
     * @return Source color with adjusted chroma. */
    Color4 get(float t) const {
      Color4 color = source.get(t);
      constexpr float DEFAULT_CHROMA =
          Pullback::Color::ColorControls{}.palette_chroma;
      if (chroma != DEFAULT_CHROMA) {
        OKLCH lch = pixel_to_oklch(color.color);
        lch.C *= chroma / DEFAULT_CHROMA;
        color.color = oklch_to_pixel(lch);
      }
      return color;
    }
  };
  /** @brief Emits one row of named cosine coefficients.
   * @param name Palette identifier.
   * @param A Bias triple.
   * @param B Amplitude triple.
   * @param C Frequency triple.
   * @param D Phase triple. */
#define HS_SHADER_PALETTE_COEFFICIENTS(name, A, B, C, D)                       \
  Coefficients{{HS_PALETTE_VEC3 A},                                            \
               {HS_PALETTE_VEC3 B},                                            \
               {HS_PALETTE_VEC3 C},                                            \
               {HS_PALETTE_VEC3 D}},
  inline static constexpr Coefficients COEFFICIENTS[] = {
      HS_PROCEDURAL_PALETTE_LIST(HS_SHADER_PALETTE_COEFFICIENTS)};
#undef HS_SHADER_PALETTE_COEFFICIENTS
  static_assert(std::size(COEFFICIENTS) ==
                std::size(Pullback::Interp::Op::PALETTE_NAME_IDS));

  BakedPaletteStorage palette;
  Pullback::Interp::Op::PaletteSelection previous;
  float previous_chroma = 0.0f;
  bool initialized = false;
};
