/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette_selection.h
 * @brief Authored palette source controls shared by fixed and interpreted pipelines. */

#include <array>
#include <cstdint>
#include "core/color/generative_palette.h"

/** @brief Palette source selected by authored color stages. */
enum class AuthoredPaletteType : uint8_t {
  CYCLING_GENERATED, ///< Animated generated harmony cycle.
  GENERATED,         ///< Authored perceptual palette.
  PROCEDURAL,        ///< Named cosine palette.
  GRADIENT,          ///< Three-stop perceptual gradient.
  MUTATING           ///< Interpolated cosine coefficients.
};
/// Authored palette source controls, persisted with the color operator.
struct AuthoredPaletteSelection {
  float palette_hue = 0.0f;     ///< Base hue in turns.
  float palette_spread = 0.07f; ///< Harmony spread in turns.
  float palette_lightness = PALETTE_DEFAULT_LIGHTNESS; ///< Lightness midpoint.
  float palette_lightness_range = 0.0f; ///< Lightness profile range.
  float palette_mutation = 0.0f;        ///< Coefficient interpolation weight.
  float gradient_start_r = 1.0f;        ///< First stop red, sRGB.
  float gradient_start_g = 0.15f;       ///< First stop green, sRGB.
  float gradient_start_b = 0.0f;        ///< First stop blue, sRGB.
  float gradient_middle_r = 0.1f;       ///< Middle stop red, sRGB.
  float gradient_middle_g = 1.0f;       ///< Middle stop green, sRGB.
  float gradient_middle_b = 0.4f;       ///< Middle stop blue, sRGB.
  float gradient_end_r = 0.1f;          ///< Last stop red, sRGB.
  float gradient_end_g = 0.2f;          ///< Last stop green, sRGB.
  float gradient_end_b = 1.0f;          ///< Last stop blue, sRGB.
  uint8_t palette_type = 0;             ///< An AuthoredPaletteType ordinal.
  /// Perceptual hue relationship.
  uint8_t palette_harmony = static_cast<uint8_t>(PaletteHarmony::TRIADIC);
  uint8_t palette_lightness_curve = 0; ///< Lightness axis profile.
  uint8_t palette_name = 0;            ///< First named cosine palette.
  uint8_t palette_name_b = 19;         ///< Second named cosine palette.
  std::array<uint8_t, 3> reserved{};   ///< Stable serialized base layout.
  constexpr bool operator==(const AuthoredPaletteSelection &) const = default;
};
static_assert(sizeof(AuthoredPaletteSelection) == 64);

/// Wire ids of authored palette source families.
inline constexpr const char *AUTHORED_PALETTE_TYPE_IDS[] = {
    "cycling-generated", "generated", "procedural", "gradient", "mutating"};
static_assert(std::size(AUTHORED_PALETTE_TYPE_IDS) ==
              static_cast<size_t>(AuthoredPaletteType::MUTATING) + 1);

/// Wire ids of perceptual hue relationships.
inline constexpr const char *AUTHORED_PALETTE_HARMONY_IDS[] = {
    "monochromatic",       "analogous", "accented-analogous", "complementary",
    "split-complementary", "triadic",   "tetradic",           "square"};
static_assert(std::size(AUTHORED_PALETTE_HARMONY_IDS) ==
              static_cast<size_t>(PaletteHarmony::SQUARE) + 1);

/// Wire ids of lightness axis curves.
inline constexpr const char *AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS[] = {
    "constant", "ascending", "descending", "bell", "cup"};
static_assert(std::size(AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS) ==
              static_cast<size_t>(AxisCurve::CUP) + 1);

/// Wire ids of named procedural palettes.
inline constexpr const char *AUTHORED_PALETTE_NAME_IDS[] = {
    "dark-rainbow", "blood-stream",   "vintage-sunset", "rich-sunset",
    "undersea",     "late-sunset",    "mango-peel",     "ice-melt",
    "lemon-lime",   "algae",          "embers",         "fire-glow",
    "dark-primary", "mauve-fade",     "lavender-lake",  "desert-rose",
    "bruised-moss", "bruised-banana", "bright-sunrise", "fire-and-ice",
    "peach-pop",    "tidal-jade",     "blue-lagoon",    "orange-crush",
    "plum-sunrise", "coral-blue",     "bruised-mango"};
