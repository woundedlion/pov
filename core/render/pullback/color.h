/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

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

/** @brief Palette and hue parameters, shared by every look. */
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
      Field<ColorParams>{&ColorParams::hue_shift_amount, nullptr, -4.0f, 4.0f,
                         FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::hue_noise_scale, nullptr, 1.0f / 64.0f,
                         8.0f, FieldCurve::LOG_POSITIVE},
      Field<ColorParams>{&ColorParams::hue_noise_speed, nullptr, -0.001f,
                         0.001f, FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::palette_chroma, nullptr, 0.0f, 1.0f,
                         FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::mapping_frequency, nullptr, 1.0f, 32.0f,
                         FieldCurve::LOG_POSITIVE},
      Field<ColorParams>{&ColorParams::mapping_phase, nullptr, -1.0f, 1.0f,
                         FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::phase_oscillation_depth, nullptr, 0.0f,
                         1.0f, FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::phase_oscillation_speed, nullptr, -0.01f,
                         0.01f, FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::brightness_depth, nullptr, 0.0f, 1.0f,
                         FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::opacity_low, nullptr, 0.0f, 1.0f,
                         FieldCurve::LERP},
      Field<ColorParams>{&ColorParams::opacity_high, nullptr, 0.0f, 1.0f,
                         FieldCurve::LERP},
  };
};

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

enum class HueMode : uint8_t { NONE = 0, NOISE = 1, PATH_LENGTH = 2 };

struct HueRotationLutView {
  static constexpr int VALUE_STEPS = 64;
  static constexpr int HUE_STEPS = 16;
  static constexpr size_t SIZE = VALUE_STEPS * HUE_STEPS;

  const Pixel *data;
  bool active;
};

struct HueNoiseLutView {
  static constexpr int FACE_COUNT = 6;
  static constexpr int FACE_STEPS = 24;
  static constexpr int FACE_SIZE = FACE_STEPS * FACE_STEPS;
  static constexpr size_t SIZE = FACE_COUNT * FACE_SIZE;

  const int8_t *data;
  bool active;
};

template <typename Palette>
HS_FLASH_INLINE inline void
prepare_hue_rotation_lut(std::span<Pixel, HueRotationLutView::SIZE> output,
                         const Palette &palette) {
  for (int value_index = 0; value_index < HueRotationLutView::VALUE_STEPS;
       ++value_index) {
    const float value =
        UNIT_OPEN_MAX * value_index / (HueRotationLutView::VALUE_STEPS - 1);
    const HueRotateBase base = make_hue_rotate_base(palette.get(value));
    for (int hue_index = 0; hue_index < HueRotationLutView::HUE_STEPS;
         ++hue_index) {
      const float amount =
          static_cast<float>(hue_index) / HueRotationLutView::HUE_STEPS;
      output[value_index * HueRotationLutView::HUE_STEPS + hue_index] =
          hue_rotate_lut_gamut(base, amount).color;
    }
  }
}

__attribute__((always_inline)) inline Vector
hue_noise_face_direction(int face, float u, float v) {
  switch (face) {
  case 0:
    return Vector(1.0f, v, u).normalized();
  case 1:
    return Vector(-1.0f, v, -u).normalized();
  case 2:
    return Vector(u, 1.0f, v).normalized();
  case 3:
    return Vector(u, -1.0f, -v).normalized();
  case 4:
    return Vector(u, v, 1.0f).normalized();
  default:
    return Vector(-u, v, -1.0f).normalized();
  }
}

HS_FLASH_INLINE inline void
prepare_hue_noise_lut(std::span<int8_t, HueNoiseLutView::SIZE> output,
                      const FastNoiseLite &noise, float scale, float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  const Vector loop_offset(NOISE_LOOP_RADIUS * cosf(angle),
                           NOISE_LOOP_RADIUS * sinf(angle), 0.0f);
  constexpr float STEP = 2.0f / (HueNoiseLutView::FACE_STEPS - 1);
  for (int face = 0; face < HueNoiseLutView::FACE_COUNT; ++face) {
    const int face_offset = face * HueNoiseLutView::FACE_SIZE;
    for (int y = 0; y < HueNoiseLutView::FACE_STEPS; ++y) {
      const float v = -1.0f + STEP * y;
      for (int x = 0; x < HueNoiseLutView::FACE_STEPS; ++x) {
        const float u = -1.0f + STEP * x;
        const Vector direction = hue_noise_face_direction(face, u, v);
        const Vector q = scale * direction + loop_offset;
        const float sample =
            hs::clamp(noise.GetNoiseSingle(q.x, q.y, q.z), -1.0f, 1.0f);
        const int quantized =
            static_cast<int>(sample * 127.0f + (sample < 0.0f ? -0.5f : 0.5f));
        output[face_offset + y * HueNoiseLutView::FACE_STEPS + x] =
            static_cast<int8_t>(quantized);
      }
    }
  }
}

__attribute__((always_inline)) inline float
palette_mapping_coordinate(float value, PaletteMapping mapping, float frequency,
                           float offset) {
  if (mapping == PaletteMapping::LINEAR && frequency == 1.0f && offset == 0.0f)
    return value;
  const float phase =
      wrap_t(std::min(value, UNIT_OPEN_MAX) * frequency + offset);
  switch (mapping) {
  case PaletteMapping::CUP:
    return fabsf(2.0f * phase - 1.0f);
  case PaletteMapping::BELL:
    return 1.0f - fabsf(2.0f * phase - 1.0f);
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
  const float cup = fabsf(2.0f * phase - 1.0f);
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
    shape = fabsf(2.0f * value - 1.0f);
    break;
  case BrightnessEnvelope::BELL:
    shape = 1.0f - fabsf(2.0f * value - 1.0f);
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

__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view, float value,
                        float amount) {
  const float value_position = value * (HueRotationLutView::VALUE_STEPS - 1);
  const int value_low = static_cast<int>(value_position);
  const int value_high =
      std::min(value_low + 1, HueRotationLutView::VALUE_STEPS - 1);
  const uint16_t value_weight =
      frac_to_q16(value_position - static_cast<float>(value_low));

  const float hue_position = amount * HueRotationLutView::HUE_STEPS;
  int hue_low = static_cast<int>(hue_position);
  if (hue_position < static_cast<float>(hue_low))
    --hue_low;
  const uint16_t hue_weight =
      frac_to_q16(hue_position - static_cast<float>(hue_low));
  const int hue_index_low = hue_low & (HueRotationLutView::HUE_STEPS - 1);
  const int hue_index_high =
      (hue_index_low + 1) & (HueRotationLutView::HUE_STEPS - 1);
  const auto sample_row = [&](int row) {
    const int offset = row * HueRotationLutView::HUE_STEPS;
    return view.data[offset + hue_index_low].lerp16(
        view.data[offset + hue_index_high], hue_weight);
  };
  return sample_row(value_low).lerp16(sample_row(value_high), value_weight);
}

HS_FLASH_INLINE inline float sample_hue_noise_lut(const HueNoiseLutView &view,
                                                  const Vector &direction) {
  const float ax = fabsf(direction.x);
  const float ay = fabsf(direction.y);
  const float az = fabsf(direction.z);
  int face;
  float u;
  float v;
  if (ax >= ay && ax >= az) {
    const float inverse = 1.0f / ax;
    face = direction.x >= 0.0f ? 0 : 1;
    u = (direction.x >= 0.0f ? direction.z : -direction.z) * inverse;
    v = direction.y * inverse;
  } else if (ay >= az) {
    const float inverse = 1.0f / ay;
    face = direction.y >= 0.0f ? 2 : 3;
    u = direction.x * inverse;
    v = (direction.y >= 0.0f ? direction.z : -direction.z) * inverse;
  } else {
    const float inverse = 1.0f / az;
    face = direction.z >= 0.0f ? 4 : 5;
    u = (direction.z >= 0.0f ? direction.x : -direction.x) * inverse;
    v = direction.y * inverse;
  }

  constexpr float SCALE =
      0.5f * static_cast<float>(HueNoiseLutView::FACE_STEPS - 1);
  const float x_position = (u + 1.0f) * SCALE;
  const float y_position = (v + 1.0f) * SCALE;
  const int x_low =
      std::min(static_cast<int>(x_position), HueNoiseLutView::FACE_STEPS - 2);
  const int y_low =
      std::min(static_cast<int>(y_position), HueNoiseLutView::FACE_STEPS - 2);
  const float x_fraction = x_position - x_low;
  const float y_fraction = y_position - y_low;
  const int offset = face * HueNoiseLutView::FACE_SIZE +
                     y_low * HueNoiseLutView::FACE_STEPS + x_low;
  const float row_low =
      hs::lerp(static_cast<float>(view.data[offset]),
               static_cast<float>(view.data[offset + 1]), x_fraction);
  const float row_high = hs::lerp(
      static_cast<float>(view.data[offset + HueNoiseLutView::FACE_STEPS]),
      static_cast<float>(view.data[offset + HueNoiseLutView::FACE_STEPS + 1]),
      x_fraction);
  return hs::lerp(row_low, row_high, y_fraction) * (1.0f / 127.0f);
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
apply_generated_palette(const MaterialSample &sample,
                        const GeneratedPaletteState &state) {
  const float oscillation =
      state.oscillation_depth * fast_sinf(TWO_PI_F * state.oscillation_phase);
  const float palette_value = palette_mapping_coordinate(
      sample.value, state.mapping, state.mapping_frequency,
      state.mapping_phase + oscillation);
  Color4 color;
  if (state.hue_rotation.active && state.hue_mode == HueMode::NOISE) {
    const float amount = state.hue_shift_amount *
                         sample_hue_noise_lut(state.hue_noise, sample.sphere);
    color = Color4(
        sample_hue_rotation_lut(state.hue_rotation, palette_value, amount),
        1.0f);
  } else {
    color = state.palette->get(palette_value);
    if (state.hue_rotation.active) {
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

template <typename State> struct GeneratedPalette : ExactPolicy {
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  static constexpr std::array<ApproximationMetric, 3> METRICS{{
      {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MAXIMUM,
       7000.0f, "channel code"},
      {ApproximationDomain::COLOR_CHANNEL, ApproximationAggregation::MEAN,
       256.0f, "channel code"},
      {ApproximationDomain::FRAMEBUFFER, ApproximationAggregation::MAXIMUM,
       5400.0f, "channel code"},
  }};

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
  HS_FLASH_INLINE static Color4 apply(const MaterialSample &sample,
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
