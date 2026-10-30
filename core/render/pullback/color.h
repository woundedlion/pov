/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/math/3dmath.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"
#include "core/color/hue_rotation_lut.h"
#include "core/render/pullback/modulation_temporal.h"

/**
 * @file color.h
 * @brief Color and palette policies.
 */

namespace Pullback {

namespace Color {

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
};
/** @brief Prepared generated-palette endpoint pair. */
struct TemporalGeneratedPaletteState {
  GeneratedPaletteState current;  ///< Current palette.
  GeneratedPaletteState previous; ///< Retained palette.
  float fraction = 1.0f;          ///< Current endpoint weight.
};
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
  color.alpha *= sample.coverage *
                 hs::lerp(state.opacity_low, state.opacity_high, sample.value);
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

/** @brief One color modulation's prepared endpoint pair. */
struct ModulationState {
  const Modulation::Prepared *current = nullptr;  ///< Current spatial source.
  const Modulation::Prepared *previous = nullptr; ///< Retained spatial source.
  Pattern::SourceKind source = Pattern::SourceKind::NOISE; ///< Scalar source.
  ColorModulation::Params mechanic; ///< Color operation and strength.
  bool enabled = false;             ///< Whether this layer contributes.
};
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
    sample.pending_hue = math::wrap_t(
        sample.pending_hue + math::wrap_t(state.mechanic.amount * scalar));
  else if (state.mechanic.kind ==
           ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT) {
    const float LIFT = state.mechanic.amount * hs::clamp(scalar, 0.0f, 1.0f);
    sample.pending_lift += (1.0f - sample.pending_lift) * LIFT;
  } else {
    sample.color = ColorModulation::apply(sample.color, scalar, state.mechanic);
  }
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

/// Qualified generated-palette approximation bounds.
inline constexpr std::array<ApproximationMetric, 3> GENERATED_PALETTE_METRICS{{
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
  using Prepared =
      std::conditional_t<TEMPORAL, TemporalGeneratedPaletteState,
                         GeneratedPaletteState>; ///< Endpoint state.
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
      return {current, previous, resources.fraction};
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
  static constexpr bool APPROXIMATE =
      true; ///< Uses scalar and hue lookup tables.
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;    ///< Reference.
  static constexpr auto METRICS = GENERATED_PALETTE_METRICS; ///< Error bounds.

  using Prepared = ModulationState; ///< Immutable layer endpoints.
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::prepare(frame) } -> std::same_as<ModulationState>;
      };
  /** @brief Prepares source endpoints and continuous controls.
   * @tparam FrameState Frame context.
   * @param frame Frame context.
   * @return Layer endpoints. */
  template <typename FrameState>
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
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
