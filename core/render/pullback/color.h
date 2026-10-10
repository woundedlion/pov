/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "math/3dmath.h"
#include "render/pullback/contract.h"
#include "render/pullback/fields.h"
#include "color/noise_hue_palette.h"

/**
 * @file color.h
 * @brief Color and palette policies.
 */

namespace Pullback {

namespace Color {

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

/** @brief What drives the color stage's hue rotation, if anything. */
enum class HueMode : uint8_t {
  NONE = 0,  /**< No hue rotation; the palette color is used as sampled. */
  NOISE = 1, /**< Rotation amount read from a cube-face noise LUT. */
  PATH_LENGTH = 2,                /**< Accumulated path length. */
  WARP_DISPLACEMENT = PATH_LENGTH /**< Alias of PATH_LENGTH. */
};

/** Activation relations the chain's colorize operators export: the hue and
    brightness controls are read only under the matching topology value. */
inline constexpr TopologyGate HUE_ROTATION_GATE{
    "hue-shift-mode", live_values(HueMode::NOISE, HueMode::PATH_LENGTH)};
/// Gate of the hue-noise controls: live only in HueMode::NOISE.
inline constexpr TopologyGate HUE_NOISE_GATE{"hue-shift-mode",
                                             live_values(HueMode::NOISE)};
/// Gate of the brightness controls: live for every envelope but NONE.
inline constexpr TopologyGate BRIGHTNESS_ENVELOPE_GATE{
    "brightness-envelope",
    live_values(BrightnessEnvelope::CUP, BrightnessEnvelope::BELL,
                BrightnessEnvelope::ASCENDING, BrightnessEnvelope::DESCENDING)};

/** @brief Continuous palette and hue controls. */
struct ColorControls {
  float hue_shift_amount = 0.0f; /**< Hue rotation magnitude; 0 disables the
                                      rotation entirely. */
  float hue_noise_scale = 1.0f;  /**< Spatial scale of the hue-noise LUT. */
  float hue_noise_speed = 0.0f;  /**< Per-frame advance of the hue-noise loop
                                      phase; nonzero speed rebakes each frame. */
  float palette_chroma = 0.62f;  /**< Chroma the generated palettes are baked
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
  float opacity_low = 1.0f;       /**< Alpha gain at source value 0. */
  float opacity_high = 1.0f;      /**< Alpha gain at source value 1. */

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      Field<ColorControls>{"hue-shift-amount", &ColorControls::hue_shift_amount,
                           nullptr, -4.0f, 4.0f, FieldCurve::LERP,
                           FieldGate::ALWAYS, HUE_ROTATION_GATE},
      Field<ColorControls>{"hue-noise-scale", &ColorControls::hue_noise_scale,
                           nullptr, 1.0f / 64.0f, 8.0f,
                           FieldCurve::LOG_POSITIVE, FieldGate::ALWAYS,
                           HUE_NOISE_GATE},
      Field<ColorControls>{"hue-noise-speed", &ColorControls::hue_noise_speed,
                           nullptr, -0.001f, 0.001f, FieldCurve::LERP,
                           FieldGate::ALWAYS, HUE_NOISE_GATE},
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
  constexpr bool operator==(const ColorParams &) const = default;
};
static_assert(field_ids_unique<ColorParams>());
static_assert(field_defaults_in_range<ColorParams>());

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

using ::HueNoiseBakeCache;
using ::HueNoiseLutView;
using ::HueRotationLutView;
using ::UNIT_OPEN_MAX;
using ::hue_noise_face_direction;
using ::prepare_hue_noise_lut;
using ::prepare_hue_rotation_lut;
using ::sample_hue_noise_lut;
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

/** @brief Per-frame state of the generated-palette colorizer. */
struct GeneratedPaletteState {
  PaletteMappingWeights mapping; ///< Palette mapping curve weights.
  float mapping_frequency;       ///< Palette repeats across the value range.
  float mapping_offset;          ///< Phase offset including the oscillation.
  const BakedPalette *palette;   ///< Palette sampled; must be non-null.
  /** HueMode::NONE is carried as an inactive `hue_rotation` view. */
  HueMode hue_mode;
  float hue_shift_amount;          ///< Hue rotation magnitude.
  HueRotationLutView hue_rotation; ///< Hue rotation LUT; inactive disables.
  HueNoiseLutView hue_noise;       ///< Hue noise LUT read in HueMode::NOISE.
  BrightnessEnvelope brightness_envelope; ///< Value-driven brightness shape.
  float brightness_bottom;                ///< Gain at the envelope's low point.
  float brightness_top; ///< Gain at the envelope's high point.
  float opacity_low;    ///< Alpha gain at value 0.
  float opacity_high;   ///< Alpha gain at value 1.
};

/**
 * @brief Colours a field sample: palette lookup, hue rotation, brightness
 * envelope and value opacity.
 * @param sample Field sample; its coverage scales the alpha.
 * @param state Prepared per-frame colorizer state.
 * @return Colour with alpha = palette alpha * coverage * value opacity.
 */
HS_HOT_FLASH_MEMBER inline Color4
apply_generated_palette(const FieldSample &sample,
                        const GeneratedPaletteState &state) {
  const float palette_value =
      palette_mapping_coordinate(sample.value, state.mapping,
                                 state.mapping_frequency, state.mapping_offset);
  Color4 color;
  if (state.hue_rotation.active && state.hue_noise.active &&
      state.hue_mode == HueMode::NOISE) {
    color =
        Color4(sample_hue_rotation_lut(
                   state.hue_rotation, palette_value,
                   state.hue_shift_amount *
                       sample_hue_noise_lut(state.hue_noise, sample.sphere)),
               state.palette->get_alpha(palette_value));
  } else {
    color = state.palette->get(palette_value);
    if (state.hue_rotation.active && state.hue_mode == HueMode::PATH_LENGTH) {
      const float amount =
          math::wrap_t(state.hue_shift_amount * sample.path_length);
      if (amount != 0.0f)
        color.color =
            sample_hue_rotation_lut(state.hue_rotation, palette_value, amount);
    }
  }
  color.color =
      color.color *
      brightness_envelope_gain(sample.value, state.brightness_envelope,
                               state.brightness_bottom, state.brightness_top);
  color.alpha *= sample.coverage *
                 hs::lerp(state.opacity_low, state.opacity_high, sample.value);
  return color;
}

/** @brief Approximation bounds of the generated-palette colorizer. */
inline constexpr std::array<ApproximationMetric, 3> GENERATED_PALETTE_METRICS{{
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MAXIMUM,
     7000.0f, "channel code"},
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN, 256.0f,
     "channel code"},
    {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
     3066.0f, "channel code"},
}};

/**
 * @brief GeneratedPalette color policy.
 * @tparam State Provider with Binding and FrameState types and
 * mapping_weights(frame), mapping_frequency(frame), mapping_phase(frame),
 * oscillation_depth(frame), oscillation_phase(frame), palette(frame),
 * hue_mode(frame), hue_shift_amount(frame), hue_rotation(frame),
 * hue_noise(frame), brightness_envelope(frame), brightness_bottom(frame),
 * brightness_top(frame), opacity_low(frame), opacity_high(frame)
 * accessors.
 */
template <typename State> struct GeneratedPalette : ApproximationDefaults {
  static constexpr bool APPROXIMATE = true; ///< Output uses LUT approximations.
  /// Oracle the approximation is checked against.
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  /// Approximation bounds; GENERATED_PALETTE_METRICS.
  static constexpr auto METRICS = GENERATED_PALETTE_METRICS;

  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
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
        { State::hue_mode(frame) } -> std::same_as<HueMode>;
        { State::hue_shift_amount(frame) } -> std::same_as<float>;
        { State::hue_rotation(frame) } -> std::same_as<HueRotationLutView>;
        { State::hue_noise(frame) } -> std::same_as<HueNoiseLutView>;
        {
          State::brightness_envelope(frame)
        } -> std::same_as<BrightnessEnvelope>;
        { State::brightness_bottom(frame) } -> std::same_as<float>;
        { State::brightness_top(frame) } -> std::same_as<float>;
        { State::opacity_low(frame) } -> std::same_as<float>;
        { State::opacity_high(frame) } -> std::same_as<float>;
      };

  using Prepared = GeneratedPaletteState; ///< Per-frame colorizer state.

  /**
   * @brief Reads every accessor once and folds the phase oscillation in.
   * @tparam FrameState Frame state of the provider's binding.
   * @param frame Current frame state.
   * @return This frame's colorizer state.
   */
  template <typename FrameState>
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return {State::mapping_weights(frame),
            State::mapping_frequency(frame),
            State::mapping_phase(frame) +
                State::oscillation_depth(frame) *
                    math::fast_sinf(math::TWO_PI_F *
                                    State::oscillation_phase(frame)),
            &State::palette(frame),
            State::hue_mode(frame),
            State::hue_shift_amount(frame),
            State::hue_rotation(frame),
            State::hue_noise(frame),
            State::brightness_envelope(frame),
            State::brightness_bottom(frame),
            State::brightness_top(frame),
            State::opacity_low(frame),
            State::opacity_high(frame)};
  }

  /**
   * @brief Colours one field sample.
   * @tparam FrameState Frame state of the provider's binding.
   * @param sample Field sample to colour.
   * @param prepared This frame's colorizer state.
   * @return apply_generated_palette(sample, prepared).
   */
  template <typename FrameState>
  HS_O3_FN static Color4 apply(const FieldSample &sample, const FrameState &,
                               const Prepared &prepared) {
    return apply_generated_palette(sample, prepared);
  }
};

} // namespace Color

} // namespace Pullback
