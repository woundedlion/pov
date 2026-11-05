/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/math/3dmath.h"
#include "core/color/palette_selection.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"
#include "core/color/hue_rotation_lut.h"
#include "core/color/palette_lab_lut.h"
#include "core/render/pullback/modulation_temporal.h"
#include "core/render/pullback/chamber_noise_cache.h"

/**
 * @file color.h
 * @brief Color and palette policies.
 */

namespace Pullback {

namespace Color {

/// Palette harmony mode wire values.
inline constexpr const char *PALETTE_MODE_IDS[] = {"triadic", "complementary",
                                                   "analogous"};
/// Brightness envelope wire values.
inline constexpr const char *BRIGHTNESS_ENVELOPE_IDS[] = {
    "none", "cup", "bell", "ascending", "descending"};

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

/// Gate of the brightness controls: live for every envelope but NONE.
inline constexpr TopologyGate BRIGHTNESS_ENVELOPE_GATE{
    "brightness-envelope",
    live_values(BrightnessEnvelope::CUP, BrightnessEnvelope::BELL,
                BrightnessEnvelope::ASCENDING, BrightnessEnvelope::DESCENDING)};

/** @brief Continuous palette and hue controls. */
struct ColorControls {
  float palette_chroma = 0.62f; /**< Chroma the generated palettes are baked
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

/** @brief Full authored palette controls for composed pipelines. */
struct AuthoredColorParams : ColorParams, AuthoredPaletteSelection {
  uint8_t palette_mode =
      0; ///< Cycling harmony: triadic, complementary, analogous.
  BrightnessEnvelope brightness_envelope =
      BrightnessEnvelope::NONE; ///< Value envelope.
  /// Numeric palette selection and mapping controls.
  static const std::array<Field<AuthoredColorParams>,
                          ColorParams::FIELDS.size() + 14>
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
                            ColorParams::FIELDS.size() + 14>
    AuthoredColorParams::FIELDS = concat_fields<AuthoredColorParams>(
        ColorParams::FIELDS,
        std::array{
            Field<AuthoredColorParams>{
                "palette-hue", &AuthoredColorParams::palette_hue, "Palette Hue",
                0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
                TopologyGate{}},
            Field<AuthoredColorParams>{
                "palette-spread", &AuthoredColorParams::palette_spread,
                "Palette Spread", 0.0f, 0.25f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "palette-lightness", &AuthoredColorParams::palette_lightness,
                "Palette Lightness", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "palette-lightness-range",
                &AuthoredColorParams::palette_lightness_range,
                "Palette Lightness Range", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "palette-mutation", &AuthoredColorParams::palette_mutation,
                "Palette Mutation", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-start-r", &AuthoredColorParams::gradient_start_r,
                "Gradient Start R", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-start-g", &AuthoredColorParams::gradient_start_g,
                "Gradient Start G", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-start-b", &AuthoredColorParams::gradient_start_b,
                "Gradient Start B", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-middle-r", &AuthoredColorParams::gradient_middle_r,
                "Gradient Middle R", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-middle-g", &AuthoredColorParams::gradient_middle_g,
                "Gradient Middle G", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-middle-b", &AuthoredColorParams::gradient_middle_b,
                "Gradient Middle B", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-end-r", &AuthoredColorParams::gradient_end_r,
                "Gradient End R", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-end-g", &AuthoredColorParams::gradient_end_g,
                "Gradient End G", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
            Field<AuthoredColorParams>{
                "gradient-end-b", &AuthoredColorParams::gradient_end_b,
                "Gradient End B", 0.0f, 1.0f, FieldCurve::LERP,
                FieldGate::ALWAYS, TopologyGate{}},
        });
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

/** @brief Per-frame generated-palette controls and immutable resources. */
struct GeneratedPaletteState {
  PaletteMappingWeights mapping;   ///< Mapping weights.
  float mapping_frequency;         ///< Palette repeats.
  float mapping_offset;            ///< Mapping phase.
  const BakedPalette *palette;     ///< Base palette.
  HueRotationLutView hue_rotation; ///< Optional base-palette hue table.
  BrightnessEnvelope brightness_envelope =
      BrightnessEnvelope::NONE;   ///< Gain shape.
  float brightness_bottom = 0.0f; ///< Low gain.
  float brightness_top = 1.0f;    ///< High gain.
  float opacity_low = 1.0f;       ///< Alpha at value zero.
  float opacity_high = 1.0f;      ///< Alpha at value one.
};
/** @brief Retained palette resources and column progress. */
struct TemporalPaletteResources {
  const BakedPalette *palette = nullptr; ///< Previous palette.
  const Pixel *hue_rotation = nullptr;   ///< Previous hue table.
  float fraction = 1.0f;                 ///< Current endpoint weight.
  bool same_palette = false; ///< Identical endpoint palette entries.
  bool phase_compatible =
      false; ///< Stable modulation recipes across the interval.
};
/** @brief Prepared generated-palette endpoint pair. */
struct TemporalGeneratedPaletteState {
  GeneratedPaletteState current;  ///< Current palette.
  GeneratedPaletteState previous; ///< Retained palette.
  float fraction = 1.0f;          ///< Current endpoint weight.
  PaletteLabTimeView lab{};       ///< Optional interval base cache.
};
/** @brief Samples an endpoint with shared palette controls.
 * @param sample Shared field coordinates.
 * @param state Endpoint resources.
 * @param value Prepared palette coordinate.
 * @param opacity Shared coverage and opacity.
 * @param gain Shared final brightness.
 * @return Canonical deferred endpoint. */
HS_HOT_INLINE inline ColorSample
apply_palette_endpoint_impl(const FieldSample &sample,
                            const GeneratedPaletteState &state, float value,
                            float opacity, float gain) {
  Color4 color = state.palette->get(value);
  color.alpha *= opacity;
  return ColorSample{
      color,
      sample,
      state.hue_rotation.active ? state.hue_rotation.data : nullptr,
      value,
      0.0f,
      gain};
}
/** @brief Samples a canonical endpoint from ITCM.
 * @param sample Shared field coordinates.
 * @param state Endpoint resources.
 * @param value Prepared palette coordinate.
 * @param opacity Shared coverage and opacity.
 * @param gain Shared final brightness.
 * @return Canonical deferred endpoint. */
HS_HOT_ITCM_MEMBER inline ColorSample
apply_palette_endpoint(const FieldSample &sample,
                       const GeneratedPaletteState &state, float value,
                       float opacity, float gain) {
  return apply_palette_endpoint_impl(sample, state, value, opacity, gain);
}
/** @brief Resolves palette opacity at the sampled field value.
 * @param sample Sampled field and coverage.
 * @param state Palette opacity controls.
 * @return Coverage-scaled opacity. */
HS_HOT_INLINE inline float palette_opacity(const FieldSample &sample,
                                           const GeneratedPaletteState &state) {
  return sample.coverage *
         hs::lerp(state.opacity_low, state.opacity_high, sample.value);
}
/** @brief Maps a sampled field into the color carrier.
 * @param sample Sampled field.
 * @param state Palette controls and resources.
 * @return Base color with deferred brightness. */
HS_HOT_ITCM_MEMBER inline ColorSample
apply_generated_palette(const FieldSample &sample,
                        const GeneratedPaletteState &state) {
  const float VALUE =
      palette_mapping_coordinate(sample.value, state.mapping,
                                 state.mapping_frequency, state.mapping_offset);
  Color4 color = state.palette->get(VALUE);
  color.alpha *= palette_opacity(sample, state);
  return {color,
          sample,
          state.hue_rotation.active ? state.hue_rotation.data : nullptr,
          VALUE,
          0.0f,
          brightness_envelope_gain(sample.value, state.brightness_envelope,
                                   state.brightness_bottom,
                                   state.brightness_top)};
}
/** @brief Completes deferred color adjustments.
 * @param sample Color carrier.
 * @return Final straight-alpha color. */
inline Color4 finish(const ColorSample &sample) {
  return finish_color_sample(sample);
}

/** @brief Blends normalized hue endpoints over the shortest unit-period arc.
 * @param previous Retained hue in [0, 1).
 * @param current Current hue in [0, 1).
 * @param fraction Current endpoint weight.
 * @return Wrapped hue, or the selected endpoint outside the interior. */
HS_HOT_INLINE inline float blend_unit_hue(float previous, float current,
                                          float fraction) {
  if (fraction <= 0.0f)
    return previous;
  if (fraction >= 1.0f)
    return current;
  assert(previous >= 0.0f && previous < 1.0f && current >= 0.0f &&
         current < 1.0f);
  float delta = current - previous;
  if (delta < -0.5f)
    delta += 1.0f;
  if (delta >= 0.5f)
    delta -= 1.0f;
  float result = previous + delta * fraction;
  if (result < 0.0f)
    result += 1.0f;
  return result >= 1.0f ? result - 1.0f : result;
}

/** @brief One color modulation's prepared endpoint pair. */
struct ModulationState {
  const Modulation::Prepared *current = nullptr;  ///< Current spatial source.
  const Modulation::Prepared *previous = nullptr; ///< Retained spatial source.
  Pattern::SourceKind source = Pattern::SourceKind::NOISE; ///< Scalar source.
  ColorModulation::Params mechanic; ///< Color operation and strength.
  bool enabled = false;             ///< Whether this layer contributes.
  bool paired_cube_ready =
      false; ///< Immutable endpoints admit shared cube taps.
};

/** @brief Eligible spatial source prepared at the column phase. */
struct PhaseModulationState : ModulationState {
  /// Interpolated axial phases; a directional sine uses the first wave.
  Pattern::SphereSmoothGrid::Prepared grid{};
  const Pattern::SphereNoise::Prepared *noise =
      nullptr;                 ///< Retained generator.
  math::Vector noise_offset{}; ///< Noise loop at the column phase.
  bool ready = false;          ///< Whether one shared color tail is admitted.
};

/** @brief Endpoint noise cubes blended by a baked-noise color tail. */
struct CubePhaseModulationState : PhaseModulationState {
  const Pattern::SphereLutView *cube_previous = nullptr; ///< Retained cube.
  const Pattern::SphereLutView *cube_current = nullptr;  ///< Current cube.
  float cube_fraction = 0.0f; ///< Current cube weight.
};
static_assert(sizeof(Pattern::SphereSmoothGrid::Prepared) == 40);

/** @brief Optional chamber patches shared by a direct-noise color tail. */
struct ChamberPhaseModulationState : PhaseModulationState {
  const Modulation::ChamberNoisePatch *noise_previous =
      nullptr; ///< Retained samples.
  const Modulation::ChamberNoisePatch *noise_current =
      nullptr;                       ///< Current samples.
  math::Quaternion noise_to_chamber; ///< Inverse column projection.
  float noise_fraction = 0.0f;       ///< Current endpoint weight.
  float noise_phase = 0.0f;          ///< Loop phase for uncached lookups.
  bool lazy_noise_offset = false;    ///< Resolve the loop offset on fallback.
};

/** @brief Samples one prepared phase, sharing per-pixel cube coordinates.
 * @tparam Phase Prepared modulation phase type.
 * @tparam SinePhase Whether to sample directional-sine phase sources.
 * @param field Source coordinates and displacement travel.
 * @param prepared Ready column phase.
 * @param scratch Lookups shared by this phase tail.
 * @return Source scalar before its mechanic. */
template <typename Phase, bool SinePhase = true>
HS_HOT_INLINE inline float sample_phase_scalar(const FieldSample &field,
                                               const Phase &prepared,
                                               PhaseScratch &scratch)
  requires std::is_base_of_v<PhaseModulationState, Phase>
{
  float scalar = 0.0f;
  if (prepared.source == Pattern::SourceKind::DISPLACEMENT)
    scalar = field.path_length;
  else if (prepared.source == Pattern::SourceKind::NOISE) {
    if constexpr (std::is_same_v<Phase, CubePhaseModulationState>)
      if (prepared.cube_current) {
        if (!scratch.cube_tap_ready) {
          scratch.cube_tap =
              Pattern::prepare_sphere_lut_sample_unchecked(field.sphere);
          scratch.cube_tap_ready = true;
        }
        const auto &TAP = scratch.cube_tap;
        return hs::lerp(
            Pattern::bounded_scalar(TAP.sample(*prepared.cube_previous)),
            Pattern::bounded_scalar(TAP.sample(*prepared.cube_current)),
            prepared.cube_fraction);
      }
    bool cached = false;
    if constexpr (std::is_same_v<Phase, ChamberPhaseModulationState>)
      if (prepared.noise_previous && prepared.noise_current)
        cached = Modulation::sample_chamber_noise(
            field.sphere, prepared.noise_to_chamber, *prepared.noise_previous,
            *prepared.noise_current, prepared.noise_fraction, scalar);
    if (!cached) {
      const auto &noise = *prepared.noise;
      auto noise_offset = prepared.noise_offset;
      if constexpr (std::is_same_v<Phase, ChamberPhaseModulationState>)
        if (prepared.lazy_noise_offset)
          noise_offset = math::noise_sphere_loop_offset(prepared.noise_phase);
      const auto q = math::noise_sphere_coordinate(
          Pattern::SphereDirection::admit(field.sphere).value(), noise.scale,
          noise_offset);
      scalar =
          Pattern::bounded_scalar(noise.noise.GetNoiseSingle(q.x, q.y, q.z));
    }
  } else {
    if constexpr (SinePhase)
      if (prepared.source == Pattern::SourceKind::DIRECTIONAL_SINE)
        return Pattern::DirectionalSine::sample(
            Pattern::SphereDirection::admit(field.sphere), prepared.grid.a);
    scalar = Pattern::SphereSmoothGrid::sample(
        Pattern::SphereDirection::admit(field.sphere), prepared.grid);
  }
  return scalar;
}

/** @brief Accumulates one deferred hue rotation.
 * @param hue Prior hue in turns.
 * @param scalar Source value.
 * @param mechanic Rotation controls.
 * @return Wrapped hue in turns. */
HS_HOT_INLINE inline float
accumulate_hue(float hue, float scalar,
               const ColorModulation::Params &mechanic) {
  return math::wrap_t(hue + math::wrap_t(mechanic.amount * scalar));
}

/** @brief Samples a cached RGB palette with the accumulated phase hue.
 * @param sample Field coordinates and coverage.
 * @param prepared Cached interval base and palette controls.
 * @param hue Accumulated hue in turns.
 * @return Final straight-alpha color. */
HS_HOT_INLINE inline Color4
cached_phase_color(const FieldSample &sample,
                   const TemporalGeneratedPaletteState &prepared, float hue) {
  const auto &state = prepared.current;
  const float VALUE =
      palette_mapping_coordinate(sample.value, state.mapping,
                                 state.mapping_frequency, state.mapping_offset);
  return Color4(PaletteRgbHueLut::sample(prepared.lab, VALUE, hue),
                palette_opacity(sample, state));
}

/** @brief Samples a cached RGB palette with the accumulated phase hue.
 * @param sample Field coordinates and coverage.
 * @param prepared Cached interval base and palette controls.
 * @param hue Accumulated hue in turns.
 * @return Final straight-alpha color. */
HS_HOT_FLASH_MEMBER inline Color4
finish_cached_phase(const FieldSample &sample,
                    const TemporalGeneratedPaletteState &prepared, float hue) {
  return cached_phase_color(sample, prepared, hue);
}

/** @brief Finishes a cached RGB palette through one prepared phase.
 * @tparam Phase Prepared modulation phase type.
 * @param sample Shared field coordinates.
 * @param prepared Prepared palette endpoints.
 * @param modulation Ready column phase.
 * @param output Receives the final straight-alpha color.
 * @return Whether this phase and RGB cache admit direct finishing. */
template <typename Phase>
HS_HOT_FLASH_MEMBER inline bool
shade_cached_phase(const FieldSample &sample,
                   const TemporalGeneratedPaletteState &prepared,
                   const Phase &modulation, Color4 &output)
  requires std::is_base_of_v<PhaseModulationState, Phase>
{
  if (!prepared.lab.low || !prepared.lab.rgb_low || !prepared.lab.rgb_high ||
      prepared.fraction <= 0.0f || prepared.fraction > 1.0f ||
      prepared.current.brightness_envelope != BrightnessEnvelope::NONE ||
      !modulation.ready ||
      (modulation.enabled && modulation.mechanic.amount != 0.0f &&
       modulation.mechanic.kind != ColorModulation::Mechanic::HUE_ROTATION))
    return false;
  PhaseScratch scratch;
  const float HUE =
      modulation.enabled && modulation.mechanic.amount != 0.0f
          ? accumulate_hue(
                0.0f, Color::sample_phase_scalar(sample, modulation, scratch),
                modulation.mechanic)
          : 0.0f;
  output = finish_cached_phase(sample, prepared, HUE);
  return true;
}

/** @brief Accumulates one deferred positive lightness lift.
 * @param lift Prior lift.
 * @param scalar Source value.
 * @param mechanic Lift controls.
 * @return Combined unit-range lift. */
HS_HOT_INLINE inline float
accumulate_lift(float lift, float scalar,
                const ColorModulation::Params &mechanic) {
  const float LIFT = mechanic.amount * hs::clamp(scalar, 0.0f, 1.0f);
  return lift + (1.0f - lift) * LIFT;
}

/** @brief Updates deferred perceptual state from one scalar.
 * @param sample Carrier to update.
 * @param scalar Source value.
 * @param mechanic Ordered perceptual operation. */
HS_HOT_INLINE inline void
apply_perceptual_scalar(ColorSample &sample, float scalar,
                        const ColorModulation::Params &mechanic) {
  if (mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION)
    sample.pending_hue = accumulate_hue(sample.pending_hue, scalar, mechanic);
  else
    sample.pending_lift =
        accumulate_lift(sample.pending_lift, scalar, mechanic);
}
/** @brief Applies one ordered modulation layer.
 * @param sample Color and source coordinates, updated in place.
 * @param state Layer controls and endpoint sources.
 * @param previous Whether to sample the retained endpoint.
 */
HS_HOT_ITCM_MEMBER inline void
apply_modulation_in_place(ColorSample &sample, const ModulationState &state,
                          bool previous = false) {
  if (!state.enabled || state.mechanic.amount == 0.0f)
    return;
  float scalar;
  if (state.source == Pattern::SourceKind::DISPLACEMENT)
    scalar = sample.field.path_length;
  else {
    const auto *source = previous ? state.previous : state.current;
    assert(source);
    scalar = Modulation::sample(
        *source, Pattern::SphereDirection::admit(sample.field.sphere));
  }
  if (state.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION)
    sample.pending_hue =
        accumulate_hue(sample.pending_hue, scalar, state.mechanic);
  else if (state.mechanic.kind ==
           ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT) {
    sample.pending_lift =
        accumulate_lift(sample.pending_lift, scalar, state.mechanic);
  } else {
    sample.color = ColorModulation::apply(sample.color, scalar, state.mechanic);
  }
}

/** @brief Applies one layer to retained and current color endpoints.
 * @param previous Retained color carrier.
 * @param current Current color carrier.
 * @param state Layer controls and immutable sources. */
HS_HOT_ITCM_MEMBER inline void
apply_modulation_pair_in_place(ColorSample &previous, ColorSample &current,
                               const ModulationState &state) {
  if (!state.enabled || state.mechanic.amount == 0.0f)
    return;
  if (state.source == Pattern::SourceKind::DISPLACEMENT ||
      previous.field.sphere != current.field.sphere) {
    apply_modulation_in_place(previous, state, true);
    apply_modulation_in_place(current, state, false);
    return;
  }
  assert(state.previous && state.current);
  const auto DIRECTION = Pattern::SphereDirection::admit(current.field.sphere);
  const auto SCALARS =
      Modulation::sample_pair(*state.previous, *state.current, DIRECTION);
  const auto apply = [&](ColorSample &sample, float scalar) {
    if (state.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION ||
        state.mechanic.kind ==
            ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT)
      apply_perceptual_scalar(sample, scalar, state.mechanic);
    else
      sample.color =
          ColorModulation::apply(sample.color, scalar, state.mechanic);
  };
  apply(previous, SCALARS.first);
  apply(current, SCALARS.second);
}

/** @brief Applies one ordered modulation layer to a copy.
 * @param input Color and source coordinates.
 * @param state Layer controls and endpoint sources.
 * @param previous Whether to sample the retained endpoint.
 * @return Modulated color carrier. */
inline ColorSample apply_modulation(const ColorSample &input,
                                    const ModulationState &state,
                                    bool previous = false) {
  ColorSample output = input;
  apply_modulation_in_place(output, state, previous);
  return output;
}

/// Linear-RGB distance from a gamut gate requiring the canonical palette base.
inline constexpr float PALETTE_LAB_GAMUT_GUARD_WIDTH = 0.001f;

/** @brief Resolves lightness and hue into rotated perceptual coordinates.
 * @param lab Base perceptual color.
 * @param hue Hue rotation in turns.
 * @param lift Lightness lift toward white.
 * @return Perceptual color before gamut mapping. */
HS_HOT_INLINE inline OKLab rotated_palette_lab(OKLab lab, float hue,
                                               float lift) {
  lab.L += (1.0f - lab.L) * lift;
  float cosine;
  float sine;
  hue_sincos(hue, cosine, sine);
  return rotate_chroma(lab, cosine, sine);
}

/** @brief Tests raw RGB proximity to a discontinuous gamut gate.
 * @param raw Hue-rotated linear RGB before gamut mapping.
 * @return Whether the canonical palette base is required. */
HS_HOT_INLINE inline bool near_palette_gamut_gate(const LinRGB &raw) {
  const float LEAST = __builtin_fminf(__builtin_fminf(raw.r, raw.g), raw.b);
  const float MOST = __builtin_fmaxf(__builtin_fmaxf(raw.r, raw.g), raw.b);
  return !(fabsf(LEAST - LINEAR_RGB_GAMUT_MIN) >
               PALETTE_LAB_GAMUT_GUARD_WIDTH &&
           fabsf(MOST - LINEAR_RGB_GAMUT_MAX) > PALETTE_LAB_GAMUT_GUARD_WIDTH);
}

/** @brief Restores a cached base whose raw RGB can cross the gamut gate.
 * @tparam Restore Canonical base callback.
 * @param sample Deferred color carrier, updated on fallback.
 * @param restore Replaces the carrier's RGB base and hue-table identity.
 * @return Whether the canonical callback ran. */
template <typename Restore>
HS_HOT_INLINE inline bool
restore_cached_palette_base_near_gamut_gate(ColorSample &sample,
                                            Restore &&restore) {
  if (!sample.has_cached_lab || sample.cached_rgb ||
      sample.brightness_gain == 0.0f)
    return false;
  if (!near_palette_gamut_gate(oklab_to_linear_rgb(rotated_palette_lab(
          sample.cached_lab, sample.pending_hue, sample.pending_lift))))
    return false;
  restore();
  sample.has_cached_lab = false;
  return true;
}

/** @brief Finishes perceptual color with canonical gamut-gate fallback.
 * @tparam Canonical Canonical finished-color callback.
 * @param lab Cached perceptual base.
 * @param hue Deferred hue rotation in turns.
 * @param lift Deferred lightness lift.
 * @param gain Final RGB brightness gain.
 * @param alpha Straight alpha.
 * @param canonical Produces canonical RGB-first finished color on fallback.
 * @return Final straight-alpha color. */
template <typename Canonical>
HS_HOT_INLINE inline Color4 finish_palette_lab(OKLab lab, float hue, float lift,
                                               float gain, float alpha,
                                               Canonical &&canonical) {
  if (gain == 0.0f)
    return Color4(Pixel(0, 0, 0), alpha);
  const auto LAB = rotated_palette_lab(lab, hue, lift);
  auto output = oklab_to_linear_rgb(LAB);
  if (near_palette_gamut_gate(output))
    return canonical();
  if (!linear_rgb_in_gamut(output.r, output.g, output.b)) {
    HS_PROFILE_DEEP(gamut_clip);
    output = oklab_to_linear_rgb(gamut_scale_to_boundary_lut(LAB));
  }
  if (gain != 1.0f) {
    output.r = hs::clamp(output.r, 0.0f, 1.0f) * gain;
    output.g = hs::clamp(output.g, 0.0f, 1.0f) * gain;
    output.b = hs::clamp(output.b, 0.0f, 1.0f) * gain;
  }
  return Color4(linrgb_to_pixel(output), alpha);
}

/** @brief Finishes a guarded cached carrier using one raw RGB evaluation.
 * @tparam Restore Canonical base callback.
 * @param sample Deferred color carrier, updated on fallback.
 * @param restore Replaces the carrier's RGB base and hue-table identity.
 * @return Final straight-alpha color. */
template <typename Restore>
HS_HOT_INLINE inline Color4 finish_fused_color(ColorSample &sample,
                                               Restore &&restore) {
  if (!sample.has_cached_lab || sample.cached_rgb)
    return finish_color_sample(sample);
  return finish_palette_lab(sample.cached_lab, sample.pending_hue,
                            sample.pending_lift, sample.brightness_gain,
                            sample.color.alpha, [&] {
                              restore();
                              sample.has_cached_lab = false;
                              return finish_color_sample(sample);
                            });
}

/// Qualified hue-rotation and noise-LUT modulation bounds.
inline constexpr std::array<ApproximationMetric, 3> MODULATION_LUT_METRICS{{
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MAXIMUM,
     7000.0f, "channel code"},
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN, 256.0f,
     "channel code"},
    {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
     3066.0f, "channel code"},
}};

/** @brief Generated-palette policy with retained immutable endpoints.
 * @tparam State Palette resource provider. */
template <typename State> struct GeneratedPalette : ApproximationDefaults {
  using SharedBaseOwner = GeneratedPalette<State>;  ///< Cache-contract owner.
  using CachedPhaseOwner = GeneratedPalette<State>; ///< Direct finish owner.
  static constexpr bool PHASE_BLEND_BASE = true;    ///< Common field controls.
  /** @brief Checks the generated-palette provider contract.
   * @tparam Binding Frame binding. */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        {
          State::mapping_weights(frame)
        } -> std::same_as<PaletteMappingWeights>;
        { State::mapping_frequency(frame) } -> std::same_as<float>;
        { State::mapping_phase(frame) } -> std::same_as<float>;
        { State::oscillation_depth(frame) } -> std::same_as<float>;
        { State::oscillation_phase(frame) } -> std::same_as<float>;
        { State::palette(frame) } -> std::same_as<const BakedPalette &>;
        { State::hue_rotation(frame) } -> std::same_as<HueRotationLutView>;
        {
          State::brightness_envelope(frame)
        } -> std::same_as<BrightnessEnvelope>;
        { State::brightness_bottom(frame) } -> std::same_as<float>;
        { State::brightness_top(frame) } -> std::same_as<float>;
        { State::opacity_low(frame) } -> std::same_as<float>;
        { State::opacity_high(frame) } -> std::same_as<float>;
      };
  static constexpr bool TEMPORAL =
      requires(const typename State::Binding::FrameState &frame) {
        State::temporal_resources(frame);
      }; ///< Whether previous palettes are retained.
  static constexpr bool LAB_CACHE =
      requires(const typename State::Binding::FrameState &frame) {
        State::palette_lab_lut(frame);
      }; ///< Whether an RGB-first base cache is available.
  static constexpr bool FUSED_TEMPORAL_COLOR = [] {
    if constexpr (requires { State::FUSED_TEMPORAL_COLOR; })
      return State::FUSED_TEMPORAL_COLOR;
    else
      return false;
  }(); ///< Whether the effect admits one-carrier endpoint color blending.
  static constexpr bool APPROXIMATE =
      LAB_CACHE; ///< Optional base approximation.
  static constexpr ApproximationOracleId ORACLE =
      LAB_CACHE ? ApproximationOracleId::PALETTE_RGB_BLEND_TO_LAB_LUT
                : ApproximationOracleId::NONE; ///< Exact RGB-first reference.
  static constexpr auto METRICS = [] {
    if constexpr (LAB_CACHE)
      return std::array<ApproximationMetric, 3>{
          {{ApproximationDomain::COLOR_CHANNEL,
            ApproximationAggregation::MAXIMUM, 7000.0f, "channel code"},
           {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN,
            256.0f, "channel code"},
           {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
            4000.0f, "channel code"}}};
    else
      return std::array<ApproximationMetric, 0>{};
  }(); ///< Base-cache qualification bounds.
  using Prepared =
      std::conditional_t<TEMPORAL, TemporalGeneratedPaletteState,
                         GeneratedPaletteState>; ///< Endpoint state.
  /** @brief Finishes an eligible cached phase without a deferred carrier.
   * @tparam Layer Canonical phase modulation policy.
   * @param sample Shared field coordinates.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared column phase.
   * @param output Receives the final straight-alpha color.
   * @return Whether the RGB cache supplied the final color. */
  template <typename Layer>
  HS_HOT_INLINE static bool
  shade_cached_phase(const FieldSample &sample, const Prepared &prepared,
                     const typename Layer::Prepared &modulation, Color4 &output)
    requires(TEMPORAL && LAB_CACHE &&
             std::is_base_of_v<PhaseModulationState, typename Layer::Prepared>)
  {
    return Color::shade_cached_phase(sample, prepared, modulation, output);
  }
  /** @brief Prepares immutable endpoint controls.
   * @tparam FrameState Frame context.
   * @param frame Current frame.
   * @return Prepared palette endpoints. */
  template <typename FrameState>
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    GeneratedPaletteState current{
        State::mapping_weights(frame),
        State::mapping_frequency(frame),
        oscillated_mapping_offset(State::mapping_phase(frame),
                                  State::oscillation_depth(frame),
                                  State::oscillation_phase(frame)),
        &State::palette(frame),
        State::hue_rotation(frame),
        State::brightness_envelope(frame),
        State::brightness_bottom(frame),
        State::brightness_top(frame),
        State::opacity_low(frame),
        State::opacity_high(frame)};
    if constexpr (TEMPORAL) {
      const auto resources = State::temporal_resources(frame);
      auto previous = current;
      previous.palette = resources.palette;
      previous.hue_rotation.data = resources.hue_rotation;
      TemporalGeneratedPaletteState result{current, previous,
                                           resources.fraction};
      if constexpr (LAB_CACHE)
        if (const auto *cache = State::palette_lab_lut(frame))
          result.lab = cache->view(resources.fraction);
      return result;
    } else
      return current;
  }
  /** @brief Produces one color endpoint.
   * @tparam FrameState Frame context.
   * @param sample Field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoints.
   * @param previous Selects the retained endpoint.
   * @return Color carrier. */
  template <typename FrameState>
  static ColorSample apply_endpoint(const FieldSample &sample,
                                    const FrameState &frame,
                                    const Prepared &prepared, bool previous) {
    (void)frame;
    if constexpr (TEMPORAL)
      return apply_generated_palette(sample, previous ? prepared.previous
                                                      : prepared.current);
    else {
      (void)previous;
      return apply_generated_palette(sample, prepared);
    }
  }
  /** @brief Samples both palettes with their shared mapping and gain controls.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoints.
   * @return Retained endpoint followed by current endpoint. */
  template <typename FrameState>
  static std::pair<ColorSample, ColorSample>
  apply_pair(const FieldSample &sample, const FrameState &frame,
             const Prepared &prepared)
    requires TEMPORAL
  {
    (void)frame;
    const auto &state = prepared.current;
    const float VALUE = palette_mapping_coordinate(sample.value, state.mapping,
                                                   state.mapping_frequency,
                                                   state.mapping_offset);
    const float OPACITY = palette_opacity(sample, state);
    const float GAIN =
        brightness_envelope_gain(sample.value, state.brightness_envelope,
                                 state.brightness_bottom, state.brightness_top);
    if constexpr (LAB_CACHE) {
      const auto endpoint =
          [&](const GeneratedPaletteState &endpoint_state) HS_HOT_INLINE {
            return apply_cache_fallback_endpoint(sample, endpoint_state, VALUE,
                                                 OPACITY, GAIN);
          };
      return {endpoint(prepared.previous), endpoint(state)};
    } else {
      return {apply_palette_endpoint(sample, prepared.previous, VALUE, OPACITY,
                                     GAIN),
              apply_palette_endpoint(sample, state, VALUE, OPACITY, GAIN)};
    }
  }
  /** @brief Samples a canonical endpoint for a cache-enabled fallback.
   * @param sample Shared field coordinates.
   * @param state Endpoint resources.
   * @param value Prepared palette coordinate.
   * @param opacity Shared coverage and opacity.
   * @param gain Shared final brightness.
   * @return Canonical deferred endpoint. */
  HS_HOT_FLASH_MEMBER static ColorSample
  apply_cache_fallback_endpoint(const FieldSample &sample,
                                const GeneratedPaletteState &state, float value,
                                float opacity, float gain)
    requires LAB_CACHE
  {
    return apply_palette_endpoint_impl(sample, state, value, opacity, gain);
  }
  /** @brief Tries one cached base without fetching either palette.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoint controls and cache view.
   * @param output Receives the cached color carrier.
   * @return Whether the cache supplied the base. */
  template <typename FrameState>
  HS_HOT_INLINE static bool
  apply_shared_base(const FieldSample &sample, const FrameState &frame,
                    const Prepared &prepared, ColorSample &output)
    requires(TEMPORAL && LAB_CACHE)
  {
    (void)frame;
    if (!prepared.lab.low || prepared.fraction <= 0.0f ||
        prepared.fraction > 1.0f)
      return false;
    apply_shared_base_ready(sample, frame, prepared, output);
    return true;
  }

  /** @brief Produces a cached base after column eligibility was checked.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared eligible cache view.
   * @param output Receives the cached color carrier.
   * @pre The base view is active and its fraction is inside (0, 1]. */
  template <typename FrameState>
  HS_HOT_INLINE static void
  apply_shared_base_ready(const FieldSample &sample, const FrameState &frame,
                          const Prepared &prepared, ColorSample &output)
    requires(TEMPORAL && LAB_CACHE)
  {
    (void)frame;
    const auto &state = prepared.current;
    const float VALUE = palette_mapping_coordinate(sample.value, state.mapping,
                                                   state.mapping_frequency,
                                                   state.mapping_offset);
    const float OPACITY = palette_opacity(sample, state);
    output.color = Color4(Pixel(0, 0, 0), OPACITY);
    output.field = sample;
    output.hue_rotation = nullptr;
    output.palette_coordinate = VALUE;
    output.pending_hue = 0.0f;
    output.brightness_gain =
        brightness_envelope_gain(sample.value, state.brightness_envelope,
                                 state.brightness_bottom, state.brightness_top);
    output.pending_lift = 0.0f;
    output.cached_rgb =
        prepared.lab.rgb_low && prepared.lab.rgb_high ? &prepared.lab : nullptr;
    output.cached_lab = output.cached_rgb
                            ? OKLab{}
                            : PaletteLabLut::sample(prepared.lab, VALUE);
    output.has_cached_lab = true;
  }

private:
  /** @brief Samples the original modulation endpoints at shared coordinates.
   * @param sample Shared field coordinates.
   * @param modulation Prepared modulation endpoints.
   * @return Previous and current endpoint scalars. */
  HS_HOT_INLINE static std::pair<float, float>
  sample_fused_scalars(const FieldSample &sample,
                       const ModulationState &modulation) {
    if (modulation.source == Pattern::SourceKind::DISPLACEMENT)
      return {sample.path_length, sample.path_length};
    assert(modulation.previous && modulation.current);
    const auto DIRECTION = Pattern::SphereDirection::admit(sample.sphere);
    if (modulation.paired_cube_ready) {
      const auto TAP = Pattern::prepare_sphere_lut_sample(DIRECTION);
      return {Pattern::bounded_scalar(TAP.sample(modulation.previous->cube)),
              Pattern::bounded_scalar(TAP.sample(modulation.current->cube))};
    }
    return Modulation::sample_pair(*modulation.previous, *modulation.current,
                                   DIRECTION);
  }

  /** @brief Resolves one canonical deferred hue endpoint.
   * @param scalar Source scalar.
   * @param amount Hue turns per source unit.
   * @return Deferred hue rotation in turns. */
  HS_HOT_INLINE static float fused_endpoint_hue(float scalar, float amount) {
    return accumulate_hue(0.0f, scalar, {.amount = amount});
  }

  /** @brief Blends wrapped endpoint hues along their shortest path.
   * @param previous Retained endpoint hue.
   * @param current Current endpoint hue.
   * @param fraction Current endpoint weight.
   * @return Blended hue in turns. */
  HS_HOT_INLINE static float blend_fused_hue(float previous, float current,
                                             float fraction) {
    return previous == current ? previous
                               : blend_unit_hue(previous, current, fraction);
  }

  /** @brief Resolves the canonical shortest-path endpoint hue blend.
   * @param scalars Previous and current source scalars.
   * @param amount Hue turns per source unit.
   * @param fraction Current endpoint weight.
   * @return Deferred hue rotation in turns. */
  HS_HOT_INLINE static float
  fused_pending_hue(const std::pair<float, float> &scalars, float amount,
                    float fraction) {
    const float PREVIOUS = fused_endpoint_hue(scalars.first, amount);
    const float CURRENT = fused_endpoint_hue(scalars.second, amount);
    return blend_fused_hue(PREVIOUS, CURRENT, fraction);
  }

  /** @brief Blends palette and modulation endpoints into one color carrier.
   * @tparam FrameState Frame context.
   * @param result Fresh color carrier to populate.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight. */
  template <typename FrameState>
  HS_HOT_INLINE static void
  prepare_fused_modulation(ColorSample &result, const FieldSample &sample,
                           const FrameState &frame, const Prepared &prepared,
                           const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    if (fraction <= 0.0f || fraction >= 1.0f) {
      result = apply_endpoint(sample, frame, prepared, fraction <= 0.0f);
      apply_modulation_in_place(result, modulation, fraction <= 0.0f);
      return;
    }
    const bool CACHED = [&] {
      if constexpr (LAB_CACHE)
        return apply_shared_base(sample, frame, prepared, result);
      else
        return false;
    }();
    Color4 previous;
    Color4 current;
    const Pixel *previous_table;
    const Pixel *current_table;
    if (CACHED) {
      previous = current = result.color;
      previous_table = current_table = nullptr;
    } else {
      const auto &state = prepared.current;
      const float VALUE = palette_mapping_coordinate(
          sample.value, state.mapping, state.mapping_frequency,
          state.mapping_offset);
      const float OPACITY = palette_opacity(sample, state);
      previous = prepared.previous.palette->get(VALUE);
      current = state.palette->get(VALUE);
      previous.alpha *= OPACITY;
      current.alpha *= OPACITY;
      previous_table = prepared.previous.hue_rotation.active
                           ? prepared.previous.hue_rotation.data
                           : nullptr;
      current_table =
          state.hue_rotation.active ? state.hue_rotation.data : nullptr;
      result.field = sample;
      result.hue_rotation = current_table;
      result.palette_coordinate = VALUE;
      result.brightness_gain = brightness_envelope_gain(
          sample.value, state.brightness_envelope, state.brightness_bottom,
          state.brightness_top);
    }
    float previous_hue = 0.0f;
    float current_hue = 0.0f;
    float previous_lift = 0.0f;
    float current_lift = 0.0f;
    if (modulation.enabled && modulation.mechanic.amount != 0.0f) {
      const auto SCALARS = sample_fused_scalars(sample, modulation);
      const auto apply = [&](Color4 &color, float &hue, float &lift,
                             float scalar) {
        if (modulation.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION)
          hue = accumulate_hue(hue, scalar, modulation.mechanic);
        else if (modulation.mechanic.kind ==
                 ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT) {
          lift = accumulate_lift(lift, scalar, modulation.mechanic);
        } else
          color = ColorModulation::apply(color, scalar, modulation.mechanic);
      };
      apply(previous, previous_hue, previous_lift, SCALARS.first);
      apply(current, current_hue, current_lift, SCALARS.second);
    }
    const bool SAME_BASE =
        previous.color == current.color && previous_table == current_table;
    result.color =
        Color4(previous.color.lerp16(current.color, frac_to_q16(fraction)),
               hs::lerp(previous.alpha, current.alpha, fraction));
    result.pending_hue = blend_fused_hue(previous_hue, current_hue, fraction);
    result.pending_lift = hs::lerp(previous_lift, current_lift, fraction);
    if (!SAME_BASE)
      result.hue_rotation = nullptr;
  }

  /** @brief Restores the canonical RGB-first palette base.
   * @param result Deferred carrier to restore.
   * @param prepared Prepared palette endpoints.
   * @param fraction Current endpoint weight. */
  HS_HOT_INLINE static void restore_fused_palette_base(ColorSample &result,
                                                       const Prepared &prepared,
                                                       float fraction) {
    const Pixel A =
        prepared.previous.palette->get_color(result.palette_coordinate);
    const Pixel B =
        prepared.current.palette->get_color(result.palette_coordinate);
    const Pixel *const A_TABLE = prepared.previous.hue_rotation.active
                                     ? prepared.previous.hue_rotation.data
                                     : nullptr;
    const Pixel *const B_TABLE = prepared.current.hue_rotation.active
                                     ? prepared.current.hue_rotation.data
                                     : nullptr;
    result.color.color = A.lerp16(B, frac_to_q16(fraction));
    result.hue_rotation = A == B && A_TABLE == B_TABLE ? B_TABLE : nullptr;
  }

public:
  /** @brief Blends palette and modulation endpoints into one guarded carrier.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight.
   * @return Deferred color after canonical endpoint blending. */
  template <typename FrameState>
  HS_HOT_ITCM_MEMBER static ColorSample
  apply_fused_modulation(const FieldSample &sample, const FrameState &frame,
                         const Prepared &prepared,
                         const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    ColorSample result;
    prepare_fused_modulation(result, sample, frame, prepared, modulation,
                             fraction);
    if constexpr (LAB_CACHE)
      restore_cached_palette_base_near_gamut_gate(result, [&] {
        restore_fused_palette_base(result, prepared, fraction);
      });
    return result;
  }

  /** @brief Shades fused endpoints without returning a deferred carrier.
   * @tparam FrameState Frame context.
   * @param sample Shared field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared palette endpoints.
   * @param modulation Prepared modulation endpoints.
   * @param fraction Current endpoint weight.
   * @return Final straight-alpha color. */
  template <typename FrameState>
  HS_HOT_ITCM_MEMBER static Color4
  shade_fused_modulation(const FieldSample &sample, const FrameState &frame,
                         const Prepared &prepared,
                         const ModulationState &modulation, float fraction)
    requires(TEMPORAL && FUSED_TEMPORAL_COLOR)
  {
    if constexpr (LAB_CACHE)
      if (fraction == 1.0f && prepared.fraction == 1.0f && prepared.lab.low &&
          !prepared.lab.rgb_low && !prepared.lab.rgb_high &&
          prepared.current.brightness_envelope == BrightnessEnvelope::NONE &&
          (!modulation.enabled || modulation.mechanic.amount == 0.0f ||
           modulation.mechanic.kind ==
               ColorModulation::Mechanic::HUE_ROTATION)) {
        const auto &state = prepared.current;
        const float VALUE = palette_mapping_coordinate(
            sample.value, state.mapping, state.mapping_frequency,
            state.mapping_offset);
        const float OPACITY = palette_opacity(sample, state);
        const float HUE = [&] {
          if (!modulation.enabled || modulation.mechanic.amount == 0.0f)
            return 0.0f;
          const float SCALAR = [&] {
            if (modulation.source == Pattern::SourceKind::DISPLACEMENT)
              return sample.path_length;
            assert(modulation.current);
            return Modulation::sample(
                *modulation.current,
                Pattern::SphereDirection::admit(sample.sphere));
          }();
          return fused_endpoint_hue(SCALAR, modulation.mechanic.amount);
        }();
        return finish_palette_lab(
            PaletteLabLut::sample(prepared.lab, VALUE), HUE, 0.0f, 1.0f,
            OPACITY, [&] {
              ColorSample canonical{
                  Color4(state.palette->get_color(VALUE), OPACITY), {}};
              canonical.palette_coordinate = VALUE;
              canonical.pending_hue = HUE;
              canonical.hue_rotation =
                  state.hue_rotation.active ? state.hue_rotation.data : nullptr;
              return finish_color_sample(canonical);
            });
      }
    if constexpr (LAB_CACHE)
      if (prepared.lab.low && fraction == prepared.fraction &&
          !prepared.lab.rgb_low && !prepared.lab.rgb_high && fraction > 0.0f &&
          fraction < 1.0f && modulation.enabled &&
          modulation.mechanic.amount != 0.0f &&
          modulation.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION &&
          prepared.current.brightness_envelope == BrightnessEnvelope::NONE) {
        const auto SCALARS = sample_fused_scalars(sample, modulation);
        const auto &state = prepared.current;
        const float VALUE = palette_mapping_coordinate(
            sample.value, state.mapping, state.mapping_frequency,
            state.mapping_offset);
        const float OPACITY = palette_opacity(sample, state);
        const float HUE =
            fused_pending_hue(SCALARS, modulation.mechanic.amount, fraction);
        return finish_palette_lab(
            PaletteLabLut::sample(prepared.lab, VALUE), HUE, 0.0f, 1.0f,
            OPACITY, [&] {
              ColorSample canonical{Color4(Pixel(0, 0, 0), OPACITY), {}};
              canonical.palette_coordinate = VALUE;
              canonical.pending_hue = HUE;
              restore_fused_palette_base(canonical, prepared, fraction);
              return finish_color_sample(canonical);
            });
      }
    ColorSample result;
    prepare_fused_modulation(result, sample, frame, prepared, modulation,
                             fraction);
    if constexpr (LAB_CACHE)
      return finish_fused_color(result, [&] {
        restore_fused_palette_base(result, prepared, fraction);
      });
    else
      return finish_color_sample(result);
  }

  /** @brief Produces the current color endpoint.
   * @tparam FrameState Frame context.
   * @param sample Field coordinates.
   * @param frame Frame context.
   * @param prepared Prepared endpoints.
   * @return Color carrier. */
  template <typename FrameState>
  static ColorSample apply(const FieldSample &sample, const FrameState &frame,
                           const Prepared &prepared) {
    return apply_endpoint(sample, frame, prepared, false);
  }
};

/** @brief Ordered color modulation policy.
 * @tparam State Layer resource provider. */
template <typename State> struct Modulate : ApproximationDefaults {
  using PairedOwner = Modulate; ///< Policy admitting paired endpoint updates.
  static constexpr bool APPROXIMATE =
      true; ///< Uses scalar and hue lookup tables.
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS; ///< Reference.
  static constexpr auto METRICS = MODULATION_LUT_METRICS; ///< Error bounds.

  using Prepared = decltype(State::prepare(
      std::declval<
          const typename State::Binding::FrameState &>())); ///< Layer state.
  /// Additional internal approximation qualified separately from the LUT metrics.
  static constexpr bool PHASE_NOISE_APPROXIMATION =
      std::is_same_v<Prepared, ChamberPhaseModulationState>;
  static constexpr bool FUSED_ENDPOINTS =
      std::same_as<Prepared, ModulationState>; ///< Canonical endpoint mode.
  static constexpr bool SINE_PHASE = [] {
    if constexpr (requires { State::SINE_PHASE; })
      return State::SINE_PHASE;
    else
      return false;
  }(); ///< Whether the column phase may carry a directional sine.
  /** @brief Checks the modulation provider contract.
   * @tparam Binding Frame binding. */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::prepare(frame) } -> std::same_as<Prepared>;
      };
  /** @brief Prepares source endpoints and continuous controls.
   * @tparam FrameState Frame context.
   * @param frame Frame context.
   * @return Layer endpoints. */
  template <typename FrameState>
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    if constexpr (std::same_as<Prepared, ModulationState>) {
      auto prepared = State::prepare(frame);
      prepared.paired_cube_ready =
          prepared.previous && prepared.current &&
          prepared.previous->evaluate == Modulation::sample_cube &&
          prepared.current->evaluate == Modulation::sample_cube &&
          prepared.previous->cube.active && prepared.current->cube.active;
      return prepared;
    } else
      return State::prepare(frame);
  }
  /** @brief Whether the prepared layer supports shared phase evaluation.
   * @param prepared Layer state.
   * @return Whether this column can use the shared tail. */
  static bool phase_ready(const Prepared &prepared)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    return prepared.ready;
  }
  /** @brief Applies the column's interpolated analytic source once.
   * @tparam FrameState Frame context.
   * @param sample Carrier to update.
   * @param frame Frame context.
   * @param prepared Layer state.
   * @param scratch Cached analytic samples for this phase. */
  template <typename FrameState>
  HS_HOT_FLASH_MEMBER static void
  apply_phase_in_place(ColorSample &sample, const FrameState &frame,
                       const Prepared &prepared, PhaseScratch &scratch)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    (void)frame;
    if (!prepared.enabled || prepared.mechanic.amount == 0.0f)
      return;
    apply_perceptual_scalar(
        sample, sample_phase_scalar(sample.field, prepared, scratch),
        prepared.mechanic);
  }
  /** @brief Applies the column's interpolated analytic source once.
   * @tparam FrameState Frame context.
   * @param sample Carrier to update.
   * @param frame Frame context.
   * @param prepared Layer state. */
  template <typename FrameState>
  HS_HOT_INLINE static void apply_phase_in_place(ColorSample &sample,
                                                 const FrameState &frame,
                                                 const Prepared &prepared)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    PhaseScratch scratch;
    apply_phase_in_place(sample, frame, prepared, scratch);
  }
  /** @brief Samples a prepared column's analytic phase source.
   * @param field Source coordinates and displacement travel.
   * @param prepared Ready column phase state.
   * @return Source scalar before the mechanic is applied. */
  HS_HOT_INLINE static float sample_phase_scalar(const FieldSample &field,
                                                 const Prepared &prepared)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    PhaseScratch scratch;
    return sample_phase_scalar(field, prepared, scratch);
  }
  /** @brief Samples a prepared column's analytic phase source.
   * @param field Source coordinates and displacement travel.
   * @param prepared Ready column phase state.
   * @param scratch Lookups shared with the tail's other layers.
   * @return Source scalar before the mechanic is applied. */
  HS_HOT_INLINE static float sample_phase_scalar(const FieldSample &field,
                                                 const Prepared &prepared,
                                                 PhaseScratch &scratch)
    requires std::is_base_of_v<PhaseModulationState, Prepared>
  {
    return Color::sample_phase_scalar<Prepared, SINE_PHASE>(field, prepared,
                                                            scratch);
  }
  /** @brief Modulates one endpoint.
   * @tparam FrameState Frame context.
   * @param input Color carrier.
   * @param frame Frame context.
   * @param prepared Layer endpoints.
   * @param previous Selects the retained endpoint.
   * @return Modulated carrier. */
  template <typename FrameState>
  static ColorSample apply_endpoint(const ColorSample &input,
                                    const FrameState &frame,
                                    const Prepared &prepared, bool previous) {
    (void)frame;
    return apply_modulation(input, prepared, previous);
  }
  /** @brief Modulates one endpoint without copying its carrier.
   * @tparam FrameState Frame context.
   * @param sample Color carrier to update.
   * @param frame Frame context.
   * @param prepared Layer endpoints.
   * @param previous Selects the retained endpoint. */
  template <typename FrameState>
  static void apply_endpoint_in_place(ColorSample &sample,
                                      const FrameState &frame,
                                      const Prepared &prepared, bool previous) {
    (void)frame;
    apply_modulation_in_place(sample, prepared, previous);
  }

  /** @brief Updates both endpoints through one prepared layer.
   * @tparam FrameState Frame context.
   * @param previous Retained carrier to update.
   * @param current Current carrier to update.
   * @param frame Frame context.
   * @param prepared Prepared layer state. */
  template <typename FrameState>
  static void
  apply_endpoints_in_place(ColorSample &previous, ColorSample &current,
                           const FrameState &frame, const Prepared &prepared) {
    (void)frame;
    apply_modulation_pair_in_place(previous, current, prepared);
  }
  /** @brief Modulates the current endpoint.
   * @tparam FrameState Frame context.
   * @param input Color carrier.
   * @param frame Frame context.
   * @param prepared Layer endpoints.
   * @return Modulated carrier. */
  template <typename FrameState>
  static ColorSample apply(const ColorSample &input, const FrameState &frame,
                           const Prepared &prepared) {
    return apply_endpoint(input, frame, prepared, false);
  }
};

} // namespace Color
} // namespace Pullback
