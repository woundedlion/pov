/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"
#include "color/noise_hue_palette.h"

/**
 * @file color.h
 * @brief Color and palette policies.
 */

namespace Pullback {

namespace Color {

inline constexpr float UNIT_OPEN_MAX = 0x1.fffffep-1f;

enum class PaletteMapping : uint8_t {
  CUP = 0,
  BELL = 1,
  LINEAR = 2,
  REVERSE = 3
};

/** @brief Palette and hue parameters, shared by every composed effect. */
struct ColorParams {
  float hue_shift_amount = 0.0f; /**< Hue rotation magnitude; 0 disables the
                                      rotation entirely. */
  float hue_noise_scale = 1.0f;  /**< Spatial scale of the hue-noise LUT. */
  float hue_noise_speed = 0.0f;  /**< Per-frame advance of the hue-noise loop
                                      phase; a change rebuilds the LUT. */
  float palette_chroma = 0.62f;  /**< Chroma the generated palettes are baked
                                      at. */
  /** Palette repeats across the value range. */
  float mapping_frequency = 1.0f;
  float mapping_phase = 0.0f;           /**< Offset into the palette. */
  float phase_oscillation_depth = 0.0f; /**< Amplitude of the sinusoidal wobble
                                             added to `mapping_phase`. */
  /** Per-frame advance of that wobble. */
  float phase_oscillation_speed = 0.0f;
  float brightness_depth = 1.0f; /**< Depth of the brightness envelope; only
                                      registered when the envelope is not
                                      NONE. */
  float opacity_low = 1.0f;      /**< Alpha gain at source value 0. */
  float opacity_high = 1.0f;     /**< Alpha gain at source value 1. */
  /** Palette mapping curve; snapped, not blended, by interpolate(). */
  Pullback::Color::PaletteMapping palette_mapping =
      Pullback::Color::PaletteMapping::LINEAR;

  static constexpr auto FIELDS = std::array{
      Field<ColorParams>{"hue-shift-amount", &ColorParams::hue_shift_amount,
                         nullptr, -4.0f, 4.0f, FieldCurve::LERP},
      Field<ColorParams>{"hue-noise-scale", &ColorParams::hue_noise_scale,
                         nullptr, 1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<ColorParams>{"hue-noise-speed", &ColorParams::hue_noise_speed,
                         nullptr, -0.001f, 0.001f, FieldCurve::LERP},
      Field<ColorParams>{"palette-chroma", &ColorParams::palette_chroma,
                         nullptr, 0.0f, 1.0f, FieldCurve::LERP},
      Field<ColorParams>{"mapping-frequency", &ColorParams::mapping_frequency,
                         nullptr, 1.0f, 32.0f, FieldCurve::LOG_POSITIVE},
      Field<ColorParams>{"mapping-phase", &ColorParams::mapping_phase, nullptr,
                         -1.0f, 1.0f, FieldCurve::LERP},
      Field<ColorParams>{"phase-oscillation-depth",
                         &ColorParams::phase_oscillation_depth, nullptr, 0.0f,
                         1.0f, FieldCurve::LERP},
      Field<ColorParams>{"phase-oscillation-speed",
                         &ColorParams::phase_oscillation_speed, nullptr, -0.01f,
                         0.01f, FieldCurve::LERP},
      Field<ColorParams>{"brightness-depth", &ColorParams::brightness_depth,
                         nullptr, 0.0f, 1.0f, FieldCurve::LERP},
      Field<ColorParams>{"value-opacity-low", &ColorParams::opacity_low,
                         nullptr, 0.0f, 1.0f, FieldCurve::LERP},
      Field<ColorParams>{"value-opacity-high", &ColorParams::opacity_high,
                         nullptr, 0.0f, 1.0f, FieldCurve::LERP},
  };

  constexpr bool operator==(const ColorParams &) const = default;
};
static_assert(field_ids_unique<ColorParams>());

struct PaletteMappingWeights {
  std::array<float, 4> values{};
  uint8_t exact = 0xff;

  static constexpr PaletteMappingWeights single(PaletteMapping mapping) {
    PaletteMappingWeights result;
    result.values[static_cast<size_t>(mapping)] = 1.0f;
    result.exact = static_cast<uint8_t>(mapping);
    return result;
  }

  static constexpr PaletteMappingWeights lerp(const PaletteMappingWeights &a,
                                              const PaletteMappingWeights &b,
                                              float progress) {
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

enum class BrightnessEnvelope : uint8_t {
  NONE = 0,
  CUP = 1,
  BELL = 2,
  ASCENDING = 3,
  DESCENDING = 4
};

/** @brief What drives the color stage's hue rotation, if anything. */
enum class HueMode : uint8_t {
  NONE = 0,  /**< No hue rotation; the palette color is used as sampled. */
  NOISE = 1, /**< Rotation amount read from a cube-face noise LUT. */
  PATH_LENGTH = 2,                /**< Accumulated path length. */
  WARP_DISPLACEMENT = PATH_LENGTH /**< Workbench synonym. */
};

using ::HueNoiseBakeCache;
using ::HueNoiseLutView;
using ::HueRotationLutView;
using ::hue_noise_face_direction;
using ::prepare_hue_noise_lut;
using ::prepare_hue_rotation_lut;
using ::sample_hue_noise_lut;
using ::sample_hue_rotation_lut;

__attribute__((always_inline)) inline float
palette_mapping_coordinate(float value, PaletteMapping mapping, float frequency,
                           float offset) {
  if (mapping == PaletteMapping::LINEAR && frequency == 1.0f && offset == 0.0f)
    return value;
  const float phase =
      wrap_t(std::min(value, UNIT_OPEN_MAX) * frequency + offset);
  switch (mapping) {
  case PaletteMapping::CUP:
    return unit_cup(phase);
  case PaletteMapping::BELL:
    return unit_bell(phase);
  case PaletteMapping::LINEAR:
    return phase;
  case PaletteMapping::REVERSE:
    return 1.0f - phase;
  }
  return phase;
}

__attribute__((always_inline)) inline float
palette_mapping_coordinate(float value, const PaletteMappingWeights &weights,
                           float frequency, float offset) {
  if (weights.exact < weights.values.size())
    return palette_mapping_coordinate(
        value, static_cast<PaletteMapping>(weights.exact), frequency, offset);

  const float phase =
      wrap_t(std::min(value, UNIT_OPEN_MAX) * frequency + offset);
  const float cup = unit_cup(phase);
  const float bell = 1.0f - cup;
  return weights.values[static_cast<size_t>(PaletteMapping::CUP)] * cup +
         weights.values[static_cast<size_t>(PaletteMapping::BELL)] * bell +
         weights.values[static_cast<size_t>(PaletteMapping::LINEAR)] * phase +
         weights.values[static_cast<size_t>(PaletteMapping::REVERSE)] *
             (1.0f - phase);
}

__attribute__((always_inline)) inline float
brightness_envelope_gain(float value, BrightnessEnvelope envelope,
                         float depth) {
  float shape = 1.0f;
  switch (envelope) {
  case BrightnessEnvelope::NONE:
    return 1.0f;
  case BrightnessEnvelope::CUP:
    shape = unit_cup(value);
    break;
  case BrightnessEnvelope::BELL:
    shape = unit_bell(value);
    break;
  case BrightnessEnvelope::ASCENDING:
    shape = value;
    break;
  case BrightnessEnvelope::DESCENDING:
    shape = 1.0f - value;
    break;
  }
  return 1.0f - depth * (1.0f - shape);
}

struct GeneratedPaletteState {
  PaletteMappingWeights mapping;
  float mapping_frequency;
  float mapping_phase;
  float oscillation_depth;
  float oscillation_phase;
  const BakedPalette *palette;
  /** HueMode::NONE is carried as an inactive `hue_rotation` view. */
  HueMode hue_mode;
  float hue_shift_amount;
  HueRotationLutView hue_rotation;
  HueNoiseLutView hue_noise;
  BrightnessEnvelope brightness_envelope;
  float brightness_depth;
  float opacity_low;
  float opacity_high;
};

HS_HOT_FLASH_MEMBER inline Color4
apply_generated_palette(const FieldSample &sample,
                        const GeneratedPaletteState &state) {
  const float oscillation =
      state.oscillation_depth * fast_sinf(TWO_PI_F * state.oscillation_phase);
  const float palette_value = palette_mapping_coordinate(
      sample.value, state.mapping, state.mapping_frequency,
      state.mapping_phase + oscillation);
  Color4 color;
  if (state.hue_rotation.active && state.hue_noise.active &&
      state.hue_mode == HueMode::NOISE) {
    const NoiseHuePalette<BakedPalette> palette(
        state.palette, state.hue_rotation.data, state.hue_noise.data);
    color = palette.get(palette_value, sample.sphere, state.hue_shift_amount);
  } else {
    color = state.palette->get(palette_value);
    if (state.hue_rotation.active && state.hue_mode == HueMode::PATH_LENGTH) {
      const float amount = wrap_t(state.hue_shift_amount * sample.path_length);
      if (amount != 0.0f)
        color.color =
            sample_hue_rotation_lut(state.hue_rotation, palette_value, amount);
    }
  }
  color.color = color.color * brightness_envelope_gain(
                                  sample.value, state.brightness_envelope,
                                  state.brightness_depth);
  color.alpha *= sample.coverage *
                 hs::lerp(state.opacity_low, state.opacity_high, sample.value);
  return color;
}

/** @brief Approximation bounds of the generated-palette colorizer, shared by
    the policy and the chain operator. */
inline constexpr std::array<ApproximationMetric, 3> GENERATED_PALETTE_METRICS{{
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MAXIMUM,
     7000.0f, "channel code"},
    {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN, 256.0f,
     "channel code"},
    {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
     5400.0f, "channel code"},
}};

template <typename State> struct GeneratedPalette : ApproximationDefaults {
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  static constexpr auto METRICS = GENERATED_PALETTE_METRICS;

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
        { State::brightness_depth(frame) } -> std::same_as<float>;
        { State::opacity_low(frame) } -> std::same_as<float>;
        { State::opacity_high(frame) } -> std::same_as<float>;
      };

  template <typename FrameState>
  HS_FLASH_INLINE static Color4 apply(const FieldSample &sample,
                                      const FrameState &frame) {
    return apply_generated_palette(
        sample,
        {State::mapping_weights(frame), State::mapping_frequency(frame),
         State::mapping_phase(frame), State::oscillation_depth(frame),
         State::oscillation_phase(frame), &State::palette(frame),
         State::hue_mode(frame), State::hue_shift_amount(frame),
         State::hue_rotation(frame), State::hue_noise(frame),
         State::brightness_envelope(frame), State::brightness_depth(frame),
         State::opacity_low(frame), State::opacity_high(frame)});
  }
};

} // namespace Color

} // namespace Pullback
