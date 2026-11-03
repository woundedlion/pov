/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file authored_palette.h
 * @brief Cached authored palette-source selection and baking.
 */

#include "core/color/palettes.h"
#include "core/color/palette_selection.h"

#include <utility>
#include <variant>

/** @brief Baked authored palette with a cache keyed by its source controls. */
class AuthoredPaletteResources {
public:
  /** @brief Resolves authored controls into an owned table.
   * @param arena Persistent table allocation arena.
   * @param params Current color operator parameters.
   * @param chroma Common perceptual chroma control.
   * @return Palette view, valid until the next source change. */
  HS_SIZE_FLASH_MEMBER const BakedPalette &
  prepare(Arena &arena, const AuthoredPaletteSelection &params, float chroma) {
    const auto &selection = params;
    if (!initialized || !(selection == previous) || chroma != previous_chroma) {
      ScratchScope scratch(scratch_arena_a);
      Source *source = scratch_arena_a.make<Source>();
      using Type = AuthoredPaletteType;
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
        recipe.chroma.axis.center = chroma;
        store(arena, source->emplace<GenerativePalette>(recipe));
      } else if (type == Type::GRADIENT) {
        const auto stop = [](float position, float r, float g, float b) {
          return std::pair<float, CPixel>{
              position, CPixel(static_cast<uint8_t>(r * 255.0f + 0.5f),
                               static_cast<uint8_t>(g * 255.0f + 0.5f),
                               static_cast<uint8_t>(b * 255.0f + 0.5f))};
        };
        const auto &gradient = source->emplace<Gradient>(
            {stop(0.0f, params.gradient_start_r, params.gradient_start_g,
                  params.gradient_start_b),
             stop(0.5f, params.gradient_middle_r, params.gradient_middle_g,
                  params.gradient_middle_b),
             stop(1.0f, params.gradient_end_r, params.gradient_end_g,
                  params.gradient_end_b)});
        store(arena, ChromaSource{gradient, chroma});
      } else {
        HS_CHECK(params.palette_name < std::size(COEFFICIENTS) &&
                     params.palette_name_b < std::size(COEFFICIENTS),
                 "invalid named palette index");
        const auto &a = COEFFICIENTS[params.palette_name];
        if (type == Type::MUTATING) {
          const auto &b = COEFFICIENTS[params.palette_name_b];
          auto &mutating = source->emplace<MutatingPalette>(a.a, a.b, a.c, a.d,
                                                            b.a, b.b, b.c, b.d);
          mutating.mutate(params.palette_mutation);
          store(arena, ChromaSource{mutating, chroma});
        } else {
          HS_CHECK(type == Type::PROCEDURAL, "invalid authored palette type");
          store(arena, ChromaSource{source->emplace<ProceduralPalette>(
                                        a.a, a.b, a.c, a.d),
                                    chroma});
        }
      }
      previous = selection;
      previous_chroma = chroma;
      initialized = true;
      ++bakes;
    }
    return palette.view();
  }

  /** @brief Number of table bakes. @return Monotonic bake generation. */
  uint32_t bake_generation() const { return bakes; }

  /** @brief Number of named palette recipes. @return Roster size. */
  static constexpr size_t named_palette_count() {
    return std::size(COEFFICIENTS);
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
      constexpr float DEFAULT_CHROMA = 0.62f;
      if (chroma != DEFAULT_CHROMA) {
        OKLCH lch = pixel_to_oklch(color.color);
        lch.C *= chroma / DEFAULT_CHROMA;
        color.color = oklch_to_pixel(lch);
      }
      return color;
    }
  };
  /** @brief Fills the table, allocating it on first use.
   * @tparam Sampled Palette source type.
   * @param arena Persistent table allocation arena.
   * @param sampled Source to sample. */
  template <typename Sampled>
  __attribute__((always_inline)) void store(Arena &arena,
                                            const Sampled &sampled) {
    if (initialized)
      palette.rebake(sampled);
    else
      palette.bake(arena, sampled);
  }
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

  /** @brief Scratch storage for whichever source the selection names. */
  using Source = std::variant<std::monostate, GenerativePalette, Gradient,
                              MutatingPalette, ProceduralPalette>;

  BakedPaletteStorage palette;
  AuthoredPaletteSelection previous;
  float previous_chroma = 0.0f;
  bool initialized = false;
  uint32_t bakes = 0;
};
