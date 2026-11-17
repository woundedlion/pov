/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Palette controls and mapping weights. */

#include "core/math/scalar_kernels.h"
#include "core/math/3dmath.h"
#include "core/color/palette_selection.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"
#include "core/color/hue_rotation_lut.h"

namespace Pullback::Color {

/** @brief Which cycling generated palette a colour stage reads; indexes
    `FrameContext::palettes`. */
enum class PaletteMode : uint8_t {
  TRIADIC = 0,       ///< Triadic-harmony palette.
  COMPLEMENTARY = 1, ///< Complementary-harmony palette.
  ANALOGOUS = 2      ///< Analogous-harmony palette.
};
/// Palette harmony mode wire values, by `PaletteMode` value.
inline constexpr const char *PALETTE_MODE_IDS[] = {"triadic", "complementary",
                                                   "analogous"};
static_assert(std::size(PALETTE_MODE_IDS) ==
              static_cast<size_t>(PaletteMode::ANALOGOUS) + 1);
/// Hue relationship of each `PaletteMode`, by value.
inline constexpr PaletteHarmony PALETTE_MODE_HARMONIES[] = {
    PaletteHarmony::TRIADIC, PaletteHarmony::COMPLEMENTARY,
    PaletteHarmony::ANALOGOUS};
static_assert(std::size(PALETTE_MODE_HARMONIES) == std::size(PALETTE_MODE_IDS));

/** @brief Applies sinusoidal oscillation to a mapping phase.
 *  @param phase Base mapping phase.
 *  @param depth Oscillation amplitude.
 *  @param oscillation_phase Oscillation phase in turns.
 *  @return Oscillated mapping offset. */
__attribute__((always_inline)) inline float
oscillated_mapping_offset(float phase, float depth, float oscillation_phase) {
  return phase + depth * math::fast_sinf(math::TWO_PI_F * oscillation_phase);
}

/** @brief Curve from the wrapped palette phase to the palette coordinate. */
enum class PaletteMapping : uint8_t {
  CUP = 0,    ///< `math::unit_cup` of the phase.
  BELL = 1,   ///< `math::unit_bell` of the phase.
  LINEAR = 2, ///< The phase itself.
  REVERSE = 3 ///< One minus the phase.
};

/** @brief Shape of the value-driven brightness gain. */
enum class BrightnessEnvelope : uint8_t {
  NONE = 0,      ///< Constant gain of 1.
  CUP = 1,       ///< `math::unit_cup` of the value.
  BELL = 2,      ///< `math::unit_bell` of the value.
  ASCENDING = 3, ///< Rises linearly with the value.
  DESCENDING = 4 ///< Falls linearly with the value.
};
/// Brightness envelope wire values, by `BrightnessEnvelope` value.
inline constexpr const char *BRIGHTNESS_ENVELOPE_IDS[] = {
    "none", "cup", "bell", "ascending", "descending"};
static_assert(std::size(BRIGHTNESS_ENVELOPE_IDS) ==
              static_cast<size_t>(BrightnessEnvelope::DESCENDING) + 1);

/// Gate of the brightness controls: live for every envelope but NONE.
inline constexpr TopologyGate BRIGHTNESS_ENVELOPE_GATE{
    "brightness-envelope",
    live_values(BrightnessEnvelope::CUP, BrightnessEnvelope::BELL,
                BrightnessEnvelope::ASCENDING, BrightnessEnvelope::DESCENDING)};

/** @brief Continuous palette and hue controls. */
struct ColorControls {
  float palette_chroma =
      PALETTE_DEFAULT_CHROMA; /**< Chroma the generated palettes are baked
                                      at. */
  /** Palette repeats across the value range. */
  float mapping_frequency = 1.0f;
  float mapping_phase = 0.0f;           /**< Offset into the palette. */
  float phase_oscillation_depth = 0.0f; /**< Amplitude of the sinusoidal wobble
                                             added to `mapping_phase`. */
  /** Per-frame advance of that wobble. */
  float phase_oscillation_speed = 0.0f;
  float brightness_bottom = 0.0f; /**< Gain at the envelope's low point. */
  float brightness_top = 1.0f;    /**< Gain at the envelope's high point. */
  float opacity_low =
      1.0f; /**< Alpha gain at field value 0 (after any transfer). */
  float opacity_high =
      1.0f; /**< Alpha gain at field value 1 (after any transfer). */

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      Field<ColorControls>{"palette-chroma", &ColorControls::palette_chroma,
                           nullptr, 0.0f, 1.0f, FieldCurve::LERP},
      Field<ColorControls>{"mapping-frequency",
                           &ColorControls::mapping_frequency, nullptr, 1.0f,
                           32.0f, FieldCurve::LOG_POSITIVE},
      Field<ColorControls>{"mapping-phase", &ColorControls::mapping_phase,
                           nullptr, -1.0f, 1.0f, FieldCurve::LERP},
      Field<ColorControls>{"phase-oscillation-depth",
                           &ColorControls::phase_oscillation_depth, nullptr,
                           0.0f, 1.0f, FieldCurve::LERP},
      Field<ColorControls>{"phase-oscillation-speed",
                           &ColorControls::phase_oscillation_speed, nullptr,
                           -0.01f, 0.01f, FieldCurve::LERP},
      Field<ColorControls>{
          "brightness-bottom", &ColorControls::brightness_bottom, nullptr, 0.0f,
          1.0f, FieldCurve::LERP, FieldGate::ALWAYS, BRIGHTNESS_ENVELOPE_GATE},
      Field<ColorControls>{"brightness-top", &ColorControls::brightness_top,
                           nullptr, 0.0f, 1.0f, FieldCurve::LERP,
                           FieldGate::ALWAYS, BRIGHTNESS_ENVELOPE_GATE},
      Field<ColorControls>{"value-opacity-low", &ColorControls::opacity_low,
                           nullptr, 0.0f, 1.0f, FieldCurve::LERP},
      Field<ColorControls>{"value-opacity-high", &ColorControls::opacity_high,
                           nullptr, 0.0f, 1.0f, FieldCurve::LERP},
  };

  /** @brief Copies the controls both palette endpoints read from the sampled
   *         column frame.
   *  @param column Controls supplying mapping, oscillation, brightness and
   *  opacity. */
  constexpr void assign_column_controls(const ColorControls &column) {
    mapping_frequency = column.mapping_frequency;
    mapping_phase = column.mapping_phase;
    phase_oscillation_depth = column.phase_oscillation_depth;
    phase_oscillation_speed = column.phase_oscillation_speed;
    brightness_bottom = column.brightness_bottom;
    brightness_top = column.brightness_top;
    opacity_low = column.opacity_low;
    opacity_high = column.opacity_high;
  }
  constexpr bool operator==(const ColorControls &) const = default;
};
static_assert(field_ids_unique<ColorControls>());
static_assert(field_defaults_in_range<ColorControls>());

/** @brief Composed-effect controls with a snapped palette mapping curve. */
struct ColorParams : ColorControls {
  PaletteMapping palette_mapping = PaletteMapping::LINEAR; ///< Mapping curve.
  /// The ColorControls fields; `palette_mapping` has no field entry.
  static constexpr auto FIELDS = concat_fields<ColorParams>(
      ColorControls::FIELDS, std::array<Field<ColorParams>, 0>{});
  /** @brief Copies the controls both palette endpoints read from the sampled
   *         column frame.
   *  @param column Controls supplying the column-sampled fields. */
  constexpr void assign_column_controls(const ColorParams &column) {
    ColorControls::assign_column_controls(column);
    palette_mapping = column.palette_mapping;
  }
  constexpr bool operator==(const ColorParams &) const = default;
};
static_assert(field_ids_unique<ColorParams>());
static_assert(field_defaults_in_range<ColorParams>());

/** @brief Builds the numeric palette-selection schema.
 * @tparam Owner Family inheriting AuthoredPaletteSelection.
 * @tparam Named Whether fields supply display labels.
 * @param generated Gate of generated palette controls.
 * @param mutating Gate of cosine mutation strength.
 * @param gradient Gate of RGB gradient stops.
 * @return Ordered palette selection fields. */
template <typename Owner, bool Named>
constexpr auto palette_selection_fields(TopologyGate generated = {},
                                        TopologyGate mutating = {},
                                        TopologyGate gradient = {}) {
  return std::array{
      Field<Owner>{"palette-hue", &Owner::palette_hue,
                   Named ? "Palette Hue" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, generated},
      Field<Owner>{"palette-spread", &Owner::palette_spread,
                   Named ? "Palette Spread" : nullptr, 0.0f, 0.25f,
                   FieldCurve::LERP, FieldGate::ALWAYS, generated},
      Field<Owner>{"palette-lightness", &Owner::palette_lightness,
                   Named ? "Palette Lightness" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, generated},
      Field<Owner>{"palette-lightness-range", &Owner::palette_lightness_range,
                   Named ? "Palette Lightness Range" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, generated},
      Field<Owner>{"palette-mutation", &Owner::palette_mutation,
                   Named ? "Palette Mutation" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, mutating},
      Field<Owner>{"gradient-start-r", &Owner::gradient_start_r,
                   Named ? "Gradient Start R" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-start-g", &Owner::gradient_start_g,
                   Named ? "Gradient Start G" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-start-b", &Owner::gradient_start_b,
                   Named ? "Gradient Start B" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-middle-r", &Owner::gradient_middle_r,
                   Named ? "Gradient Middle R" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-middle-g", &Owner::gradient_middle_g,
                   Named ? "Gradient Middle G" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-middle-b", &Owner::gradient_middle_b,
                   Named ? "Gradient Middle B" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-end-r", &Owner::gradient_end_r,
                   Named ? "Gradient End R" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-end-g", &Owner::gradient_end_g,
                   Named ? "Gradient End G" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
      Field<Owner>{"gradient-end-b", &Owner::gradient_end_b,
                   Named ? "Gradient End B" : nullptr, 0.0f, 1.0f,
                   FieldCurve::LERP, FieldGate::ALWAYS, gradient},
  };
}

/// Number of numeric palette-selection fields.
inline constexpr size_t PALETTE_SELECTION_FIELD_COUNT =
    palette_selection_fields<AuthoredPaletteSelection, false>().size();

/** @brief Full authored palette controls for composed pipelines. */
struct AuthoredColorParams : ColorParams, AuthoredPaletteSelection {
  uint8_t palette_mode =
      0; ///< Cycling harmony: triadic, complementary, analogous.
  BrightnessEnvelope brightness_envelope =
      BrightnessEnvelope::NONE; ///< Value envelope.
  /// Numeric palette selection and mapping controls.
  static const std::array<Field<AuthoredColorParams>,
                          ColorParams::FIELDS.size() +
                              PALETTE_SELECTION_FIELD_COUNT>
      FIELDS;
  /** @brief Copies the controls both palette endpoints read from the sampled
   *         column frame.
   *  @param column Controls supplying the column-sampled fields. */
  constexpr void assign_column_controls(const AuthoredColorParams &column) {
    ColorParams::assign_column_controls(column);
    brightness_envelope = column.brightness_envelope;
  }
  constexpr bool operator==(const AuthoredColorParams &) const = default;
};
// ARM GCC needs the complete class for secondary-base member offsets.
inline constexpr std::array<Field<AuthoredColorParams>,
                            ColorParams::FIELDS.size() +
                                PALETTE_SELECTION_FIELD_COUNT>
    AuthoredColorParams::FIELDS = concat_fields<AuthoredColorParams>(
        ColorParams::FIELDS,
        palette_selection_fields<AuthoredColorParams, true>());
static_assert(field_ids_unique<AuthoredColorParams>());
static_assert(field_defaults_in_range<AuthoredColorParams>());

/** @brief Whether two color families define the same palette endpoints.
 * @details Controls copied by `assign_column_controls` are excluded; every
 * other control is compared.
 * @tparam Controls Color family type.
 * @param a First family.
 * @param b Second family.
 * @return True when the families differ only in column-sampled controls. */
template <typename Controls>
HS_FLASH_INLINE constexpr bool same_palette_controls(const Controls &a,
                                                     const Controls &b) {
  Controls shared = a;
  shared.assign_column_controls(b);
  return shared == b;
}

/** @brief Blend weights over the PaletteMapping curves. */
struct PaletteMappingWeights {
  std::array<float, 4> values{}; ///< Weight per PaletteMapping, by value.
  /// PaletteMapping value when exactly one curve is selected, else 0xff.
  uint8_t exact = 0xff;

  /**
   * @brief Weights selecting exactly one curve.
   * @param mapping Curve to select.
   * @return Unit weight on `mapping`, with `exact` set.
   */
  static constexpr PaletteMappingWeights single(PaletteMapping mapping) {
    PaletteMappingWeights result;
    result.values[static_cast<size_t>(mapping)] = 1.0f;
    result.exact = static_cast<uint8_t>(mapping);
    return result;
  }

  /**
   * @brief Interpolates two weight sets.
   * @param a Weights at progress 0.
   * @param b Weights at progress 1.
   * @param progress Blend amount; clamped to the endpoints.
   * @return `a` when both select the same exact curve; otherwise the
   * per-curve linear blend, inexact unless progress hits an endpoint.
   */
  static constexpr PaletteMappingWeights lerp(const PaletteMappingWeights &a,
                                              const PaletteMappingWeights &b,
                                              float progress) {
    if (a.exact == b.exact && a.exact < a.values.size())
      return a;
    if (progress <= 0.0f)
      return a;
    if (progress >= 1.0f)
      return b;
    PaletteMappingWeights result;
    for (size_t index = 0; index < result.values.size(); ++index)
      result.values[index] =
          a.values[index] + (b.values[index] - a.values[index]) * progress;
    return result;
  }
};

using ::HueRotationLutView;
using ::UNIT_OPEN_MAX;
using ::prepare_hue_rotation_lut;
using ::sample_hue_rotation_lut;

/**
 * @brief Maps a field value to a palette coordinate through one curve.
 * @param value Field value in [0, 1].
 * @param mapping Curve applied to the wrapped phase.
 * @param frequency Palette repeats across the value range.
 * @param offset Phase offset, in palette cycles.
 * @return Palette coordinate in [0, 1]; `value` unchanged for the identity
 * LINEAR mapping.
 */
__attribute__((always_inline)) inline float
palette_mapping_coordinate(float value, PaletteMapping mapping, float frequency,
                           float offset) {
  if (mapping == PaletteMapping::LINEAR && frequency == 1.0f && offset == 0.0f)
    return value;
  const float phase =
      math::wrap_t(fminf(value, UNIT_OPEN_MAX) * frequency + offset);
  switch (mapping) {
  case PaletteMapping::CUP:
    return math::unit_cup(phase);
  case PaletteMapping::BELL:
    return math::unit_bell(phase);
  case PaletteMapping::LINEAR:
    return phase;
  case PaletteMapping::REVERSE:
    return 1.0f - phase;
  }
  return phase;
}

/**
 * @brief Maps a field value to a palette coordinate through blended curves.
 * @param value Field value in [0, 1].
 * @param weights Per-curve blend weights.
 * @param frequency Palette repeats across the value range.
 * @param offset Phase offset, in palette cycles.
 * @return Weighted sum of the curves at the wrapped phase.
 */
__attribute__((always_inline)) inline float
palette_mapping_coordinate(float value, const PaletteMappingWeights &weights,
                           float frequency, float offset) {
  if (weights.exact < weights.values.size())
    return palette_mapping_coordinate(
        value, static_cast<PaletteMapping>(weights.exact), frequency, offset);

  const float phase =
      math::wrap_t(fminf(value, UNIT_OPEN_MAX) * frequency + offset);
  const float cup = math::unit_cup(phase);
  const float bell = 1.0f - cup;
  return weights.values[static_cast<size_t>(PaletteMapping::CUP)] * cup +
         weights.values[static_cast<size_t>(PaletteMapping::BELL)] * bell +
         weights.values[static_cast<size_t>(PaletteMapping::LINEAR)] * phase +
         weights.values[static_cast<size_t>(PaletteMapping::REVERSE)] *
             (1.0f - phase);
}

/**
 * @brief Brightness gain of a field value under an envelope.
 * @param value Field value in [0, 1].
 * @param envelope Envelope shape.
 * @param bottom Gain where the shape is 0.
 * @param top Gain where the shape is 1.
 * @return 1 for BrightnessEnvelope::NONE, else lerp(bottom, top, shape).
 */
__attribute__((always_inline)) inline float
brightness_envelope_gain(float value, BrightnessEnvelope envelope, float bottom,
                         float top) {
  float shape = 1.0f;
  switch (envelope) {
  case BrightnessEnvelope::NONE:
    return 1.0f;
  case BrightnessEnvelope::CUP:
    shape = math::unit_cup(value);
    break;
  case BrightnessEnvelope::BELL:
    shape = math::unit_bell(value);
    break;
  case BrightnessEnvelope::ASCENDING:
    shape = value;
    break;
  case BrightnessEnvelope::DESCENDING:
    shape = 1.0f - value;
    break;
  }
  return hs::lerp(bottom, top, shape);
}

} // namespace Pullback::Color
