/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette_recipe.h
 * @brief Palette authoring recipes and compile diagnostics.
 */

#include <array>
#include "platform/platform.h"

/** @brief Unit cup curve over [0, 1]: 1 at both ends, 0 at the midpoint. */
__attribute__((always_inline)) inline float unit_cup(float t) {
  return fabsf(2.0f * t - 1.0f);
}

/** @brief Unit bell curve over [0, 1]: 0 at both ends, 1 at the midpoint. */
__attribute__((always_inline)) inline float unit_bell(float t) {
  return 1.0f - unit_cup(t);
}

enum class HueMode : uint8_t { HARMONY, SWEEP, CUSTOM };

enum class PaletteHarmony : uint8_t {
  MONOCHROMATIC,
  ANALOGOUS,
  ACCENTED_ANALOGOUS,
  COMPLEMENTARY,
  SPLIT_COMPLEMENTARY,
  TRIADIC,
  TETRADIC,
  SQUARE,
};

enum class HueDirection : uint8_t { SHORTEST, CLOCKWISE, COUNTERCLOCKWISE };

enum class AxisCurve : uint8_t {
  CONSTANT,
  ASCENDING,
  DESCENDING,
  BELL,
  CUP,
  CUSTOM,
};

/**
 * @brief The reference a chroma axis is measured against.
 * @details PATH_MINIMUM is reserved and unimplemented: a recipe naming it fails
 * compilation with INVALID_ENUM. Its ordinal is held rather than reclaimed
 * because the numbering is mirrored by the palette authoring tools and by
 * persisted recipes.
 */
enum class ChromaBasis : uint8_t { LOCAL_GAMUT, PATH_MINIMUM, ABSOLUTE };

enum class ColorPath : uint8_t { OKLCH_ARC, OKLAB_CARTESIAN };

enum class PaletteDomain : uint8_t {
  STRAIGHT,
  MIRROR,
  VIGNETTE,
  FALLOFF,
  LOOP,
};

enum class SegmentEase : uint8_t { LINEAR, COSINE, SMOOTHSTEP };

inline constexpr uint8_t PALETTE_MAX_KEYS = 4;

struct HueControls {
  HueMode mode = HueMode::HARMONY;
  PaletteHarmony harmony = PaletteHarmony::ANALOGOUS;
  HueDirection direction = HueDirection::SHORTEST;
  float base_turns = 0.0f;
  float spread_turns = 0.07f;
  float sweep_turns = 1.0f;
  /** CUSTOM uses exactly the first three keys; the fourth is canonicalized
   * to zero. Four-key runs are available through HARMONY. */
  std::array<float, PALETTE_MAX_KEYS> custom_turns{};
};

struct AxisControls {
  AxisCurve curve = AxisCurve::CONSTANT;
  float center = 0.62f;
  float range = 0.0f;
  std::array<float, PALETTE_MAX_KEYS> custom{};
};

struct ChromaControls {
  AxisControls axis;
  ChromaBasis basis = ChromaBasis::LOCAL_GAMUT;
  float headroom = 0.94f;
};
struct PaletteInputWindow {
  float offset = 0.0f;
  float span = 1.0f;
};

struct PaletteRecipe {
  static constexpr uint8_t SCHEMA_VERSION = 4;

  uint8_t schema_version = SCHEMA_VERSION;
  PaletteInputWindow input;
  PaletteDomain domain = PaletteDomain::STRAIGHT;
  SegmentEase easing = SegmentEase::COSINE;
  ColorPath color_path = ColorPath::OKLCH_ARC;
  HueControls hue;
  AxisControls lightness;
  ChromaControls chroma;
  float hue_torsion = 0.0f;
  float falloff_start = 0.90f;
};

/**
 * @brief The verdict a recipe compile returns.
 * @details INCOMPATIBLE_OPTIONS is reserved and never produced. Its ordinal is
 * held rather than reclaimed because the numbering is mirrored by the palette
 * authoring tools and by persisted status records.
 */
enum class PaletteCompileCode : uint8_t {
  OK,
  INVALID_SCHEMA,
  NON_FINITE,
  INVALID_ENUM,
  HUE_LIMIT,
  NON_INTEGER_LOOP_SWEEP,
  INVALID_FALLOFF_START,
  INCOMPATIBLE_OPTIONS,
};

enum class PaletteRecipeField : uint8_t {
  NONE = 0,
  PALETTE_DOMAIN = 1,
  EASING = 2,
  COLOR_PATH = 3,
  HUE_MODE = 4,
  HARMONY = 5,
  HUE_DIRECTION = 6,
  BASE_TURNS = 7,
  SPREAD_TURNS = 8,
  SWEEP_TURNS = 9,
  CUSTOM_TURNS_0 = 10,
  CUSTOM_TURNS_1 = 11,
  CUSTOM_TURNS_2 = 12,
  CUSTOM_TURNS_3 = 13,
  LIGHTNESS_CURVE = 14,
  LIGHTNESS_CENTER = 15,
  LIGHTNESS_RANGE = 16,
  LIGHTNESS_CUSTOM_0 = 17,
  LIGHTNESS_CUSTOM_1 = 18,
  LIGHTNESS_CUSTOM_2 = 19,
  LIGHTNESS_CUSTOM_3 = 20,
  CHROMA_CURVE = 21,
  CHROMA_BASIS = 22,
  CHROMA_CENTER = 23,
  CHROMA_RANGE = 24,
  CHROMA_CUSTOM_0 = 25,
  CHROMA_CUSTOM_1 = 26,
  CHROMA_CUSTOM_2 = 27,
  CHROMA_CUSTOM_3 = 28,
  CHROMA_HEADROOM = 29,
  HUE_TORSION = 30,
  FALLOFF_START = 31,
  SCHEMA_VERSION = 32,
  INPUT_OFFSET = 33,
  INPUT_SPAN = 34,
  COUNT,
};

struct PaletteAdjustments {
  uint64_t wrapped_fields = 0;
  uint64_t clamped_fields = 0;
  uint64_t canonicalized_fields = 0;
};

struct PaletteCompileStatus {
  PaletteCompileCode code = PaletteCompileCode::OK;
  PaletteRecipeField field = PaletteRecipeField::NONE;
  PaletteAdjustments adjustments;
};
