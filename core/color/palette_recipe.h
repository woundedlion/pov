/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette_recipe.h
 * @brief Palette authoring recipes and compile diagnostics.
 */

#include <array>
#include "math/3dmath.h"

/** @brief How key hues are generated: harmony offsets, an even sweep, or
 *  per-key custom turns. */
enum class HueMode : uint8_t { HARMONY, SWEEP, CUSTOM };

/** @brief Hue relationship among the keys in HARMONY mode. */
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

/** @brief Direction hue travels between successive keys. */
enum class HueDirection : uint8_t { SHORTEST, CLOCKWISE, COUNTERCLOCKWISE };

/**
 * @brief Profile of a lightness or chroma axis across the keys.
 * @details CONSTANT holds one value; BELL rises to the high bound
 * mid-run and CUP dips to the low bound; CUSTOM takes per-key values.
 */
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
 * compilation with INVALID_ENUM. The numbering is persisted; do not reclaim
 * the ordinal.
 */
enum class ChromaBasis : uint8_t { LOCAL_GAMUT, PATH_MINIMUM, ABSOLUTE };

/** @brief The space segments interpolate in. */
enum class ColorPath : uint8_t { OKLCH_ARC, OKLAB_CARTESIAN };

/** @brief How the lookup coordinate is folded across the key run. */
enum class PaletteDomain : uint8_t {
  STRAIGHT,
  MIRROR,
  VIGNETTE,
  FALLOFF,
  LOOP,
};

/** @brief Easing applied to progress within a segment between keys. */
enum class SegmentEase : uint8_t { LINEAR, COSINE, SMOOTHSTEP };

/// Maximum control keys in a palette.
inline constexpr uint8_t PALETTE_MAX_KEYS = 4;

/** @brief Hue authoring controls of a recipe. */
struct HueControls {
  bool operator==(const HueControls &) const = default;

  HueMode mode = HueMode::HARMONY; ///< Hue generation mode.
  /// Relationship used in HARMONY mode.
  PaletteHarmony harmony = PaletteHarmony::ANALOGOUS;
  /// Travel direction between keys.
  HueDirection direction = HueDirection::SHORTEST;
  float base_turns = 0.0f;    ///< First-key hue, turns; wrapped into [0, 1).
  float spread_turns = 0.07f; ///< Harmony spread, turns in [0, 0.25].
  /// SWEEP-mode hue travel across the run, turns; whole turns under LOOP.
  float sweep_turns = 1.0f;
  /** CUSTOM uses exactly the first three keys; the fourth is canonicalized
   * to zero. */
  std::array<float, PALETTE_MAX_KEYS> custom_turns{};
};

/** @brief Lightness or chroma axis controls of a recipe. */
struct AxisControls {
  bool operator==(const AxisControls &) const = default;

  AxisCurve curve = AxisCurve::CONSTANT; ///< Profile across the keys.
  float center = 0.62f; ///< Axis midpoint in [0, 1]; the value when CONSTANT.
  float range = 0.0f;   ///< Full low-to-high width in [0, 1] about `center`.
  /// Per-key values in [0, 1], read only when `curve` is CUSTOM.
  std::array<float, PALETTE_MAX_KEYS> custom{};
};

/** @brief Chroma axis controls and the basis they are measured in. */
struct ChromaControls {
  bool operator==(const ChromaControls &) const = default;

  AxisControls axis;                            ///< Chroma-control profile.
  ChromaBasis basis = ChromaBasis::LOCAL_GAMUT; ///< Chroma reference.
  /// Cap on gamut-relative chroma in [0, 1]; 1 under ABSOLUTE.
  float headroom = 0.94f;
};
/** @brief Sub-range of the key run that the lookup coordinate spans. */
struct PaletteInputWindow {
  bool operator==(const PaletteInputWindow &) const = default;

  float offset = 0.0f; ///< Window start along the key run, in [0, 1].
  float span = 1.0f;   ///< Window width, in [0, 1 - offset].
};

/** @brief Authoring description a GenerativePalette compiles from. */
struct PaletteRecipe {
  bool operator==(const PaletteRecipe &) const = default;

  static constexpr uint8_t SCHEMA_VERSION = 4; ///< Current recipe layout.

  /// Layout version; compile rejects any value but `SCHEMA_VERSION`.
  uint8_t schema_version = SCHEMA_VERSION;
  PaletteInputWindow input;                       ///< Sampled key-run window.
  PaletteDomain domain = PaletteDomain::STRAIGHT; ///< Coordinate folding.
  SegmentEase easing = SegmentEase::COSINE;       ///< Per-segment easing.
  ColorPath color_path = ColorPath::OKLCH_ARC;    ///< Interpolation space.
  HueControls hue;                                ///< Key hues.
  AxisControls lightness;                         ///< Key lightness profile.
  ChromaControls chroma;                          ///< Key chroma profile.
  /// Hue offset per unit of lightness away from 0.5, radians.
  float hue_torsion = 0.0f;
  /// FALLOFF-domain coordinate where visibility reaches 0, in (2/3, 1).
  float falloff_start = 0.90f;
};

/**
 * @brief The verdict a recipe compile returns.
 * @details INCOMPATIBLE_OPTIONS is reserved and never produced. The numbering
 * is persisted; do not reclaim the ordinal.
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

/** @brief A recipe field named by a compile diagnostic; the ordinal is its
 *  bit in the PaletteAdjustments masks. */
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

/** @brief Fields a successful compile rewrote, as PaletteRecipeField bits. */
struct PaletteAdjustments {
  uint64_t wrapped_fields = 0;       ///< Wrapped into their periodic range.
  uint64_t clamped_fields = 0;       ///< Clamped into their valid range.
  uint64_t canonicalized_fields = 0; ///< Reset to a canonical value.
};

/** @brief Result of compiling a PaletteRecipe. */
struct PaletteCompileStatus {
  PaletteCompileCode code = PaletteCompileCode::OK;    ///< Verdict.
  PaletteRecipeField field = PaletteRecipeField::NONE; ///< Failing field.
  PaletteAdjustments adjustments; ///< Rewrites; empty on failure.
};
