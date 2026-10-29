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

/** @brief What drives the color stage's hue rotation, if anything. */
enum class HueMode : uint8_t {
  NONE = 0, /**< No hue rotation; the palette color is used as sampled. */
  PATH_LENGTH = 1,                /**< Accumulated path length. */
  WARP_DISPLACEMENT = PATH_LENGTH /**< Alias of PATH_LENGTH. */
};

/** Activation relations the chain's colorize operators export: the hue and
    brightness controls are read only under the matching topology value. */
inline constexpr TopologyGate HUE_ROTATION_GATE{
    "hue-shift-mode", live_values(HueMode::PATH_LENGTH)};
/// Gate of the brightness controls: live for every envelope but NONE.
inline constexpr TopologyGate BRIGHTNESS_ENVELOPE_GATE{
    "brightness-envelope",
    live_values(BrightnessEnvelope::CUP, BrightnessEnvelope::BELL,
                BrightnessEnvelope::ASCENDING, BrightnessEnvelope::DESCENDING)};

/** @brief Continuous palette and hue controls. */
struct ColorControls {
  float hue_shift_amount = 0.0f; /**< Hue rotation magnitude; 0 disables the
                                      rotation entirely. */
  ColorModulation::Config modulation = [] {
    ColorModulation::Config config;
    config.seed = HUE_NOISE_SEED;
    return config;
  }(); ///< Independent spatial color modulation.
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
      Field<ColorControls>{"hue-shift-amount", &ColorControls::hue_shift_amount,
                           nullptr, -4.0f, 4.0f, FieldCurve::LERP,
                           FieldGate::ALWAYS, HUE_ROTATION_GATE},
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

/** @brief Per-frame state of the generated-palette colorizer. */
struct GeneratedPaletteState {
  PaletteMappingWeights mapping;   ///< Mapping weights.
  float mapping_frequency;         ///< Palette repeats.
  float mapping_offset;            ///< Mapping offset in turns.
  const BakedPalette *palette;     ///< Base palette.
  HueMode hue_mode;                ///< Renderer-owned path hue selection.
  float hue_shift_amount;          ///< Path hue rotation in turns.
  HueRotationLutView hue_rotation; ///< Optional path hue table.
  const Modulation::Prepared *modulation =
      nullptr;                      ///< Immutable source endpoint.
  ColorModulation::Params mechanic; ///< Shared column mechanic controls.
  BrightnessEnvelope brightness_envelope =
      BrightnessEnvelope::NONE;   ///< Brightness shape.
  float brightness_bottom = 0.0f; ///< Low gain.
  float brightness_top = 1.0f;    ///< High gain.
  float opacity_low = 1.0f;       ///< Alpha gain at value zero.
  float opacity_high = 1.0f;      ///< Alpha gain at value one.
  float palette_alpha = -1.0f;    ///< Constant palette alpha or sentinel.
  Color4 (*apply_mechanic)(Color4, float, const ColorModulation::Params &) =
      nullptr; ///< Resolved mechanic kernel.
};
/** @brief Previous resources and column progress. */
struct TemporalPaletteResources {
  const BakedPalette *palette = nullptr; ///< Retained palette.
  const Pixel *hue_rotation = nullptr;   ///< Retained path hue table.
  const Modulation::Prepared *modulation =
      nullptr;           ///< Retained source endpoint.
  float fraction = 1.0f; ///< Current endpoint weight.
  const math::Vector *hue_rotation_linear =
      nullptr;                       ///< Optional float path hue table.
  float palette_alpha = -1.0f;       ///< Previous constant alpha or sentinel.
  bool same_palette = false;         ///< Equal endpoint palette resources.
  bool continuous_modulation = true; ///< Modulation topology matches.
};
/** @brief Completed-color interpolation endpoints. */
struct TemporalGeneratedPaletteState {
  GeneratedPaletteState current;  ///< Current endpoint resources.
  GeneratedPaletteState previous; ///< Retained endpoint resources.
  float fraction;                 ///< Current endpoint weight.
  uint16_t weight;                ///< Q16 RGB interpolation weight.
  bool same_source = false; ///< Endpoints evaluate the same source inputs.
};
/** @brief Resolves the selected color operation outside scan loops.
 * @param kind Mechanic selection.
 * @return The selected immutable mechanic kernel. */
inline auto mechanic_kernel(ColorModulation::Mechanic kind)
    -> Color4 (*)(Color4, float, const ColorModulation::Params &) {
  using namespace ColorModulation;
  switch (kind) {
  case Mechanic::HUE_ROTATION:
    return +[](Color4 c, float v, const Params &p) {
      return HueRotation::apply(c, v, {p.amount});
    };
  case Mechanic::POSITIVE_LIGHTNESS_LIFT:
    return +[](Color4 c, float v, const Params &p) {
      return PositiveLightnessLift::apply(c, v, {p.amount});
    };
  case Mechanic::OPACITY:
    return +[](Color4 c, float v, const Params &p) {
      return Opacity::apply(c, v, {p.amount});
    };
  }
  HS_CHECK(false, "unknown modulation mechanic");
  return nullptr;
}
/** @brief Samples and modulates one endpoint at shared column coordinates.
 * @param sample Renderer field sample.
 * @param state Endpoint palette and shared continuous controls.
 * @param direction Admitted coordinate, shared by both endpoints.
 * @param shared_scalar Optional source sample shared by both endpoints.
 * @return Completed straight-alpha endpoint color. */
inline Color4
generated_palette_endpoint(const FieldSample &sample,
                           const GeneratedPaletteState &state,
                           const Pattern::SphereDirection *direction,
                           const float *shared_scalar = nullptr) {
  const float VALUE =
      palette_mapping_coordinate(sample.value, state.mapping,
                                 state.mapping_frequency, state.mapping_offset);
  Color4 color = state.palette->get(VALUE);
  if (state.modulation && state.modulation->config.enabled &&
      state.mechanic.amount != 0.0f) {
    HS_CHECK(direction, "modulation requires an admitted sphere direction");
    const float SCALAR =
        shared_scalar ? *shared_scalar
                      : Modulation::sample(*state.modulation, *direction);
    if (state.mechanic.kind == ColorModulation::Mechanic::HUE_ROTATION &&
        state.hue_rotation.active) {
      const float ROTATION = math::wrap_t(state.mechanic.amount * SCALAR);
      if (ROTATION != 0.0f)
        color.color =
            sample_hue_rotation_lut(state.hue_rotation, VALUE, ROTATION);
    } else {
      color = state.apply_mechanic
                  ? state.apply_mechanic(color, SCALAR, state.mechanic)
                  : ColorModulation::apply(color, SCALAR, state.mechanic);
    }
  }
  if (state.hue_mode == HueMode::PATH_LENGTH &&
      state.hue_shift_amount != 0.0f) {
    const float ROTATION =
        math::wrap_t(state.hue_shift_amount * sample.path_length);
    if (ROTATION != 0.0f)
      color.color =
          state.hue_rotation.active
              ? sample_hue_rotation_lut(state.hue_rotation, VALUE, ROTATION)
              : ColorModulation::HueRotation::apply(color, sample.path_length,
                                                    {state.hue_shift_amount})
                    .color;
  }
  const float GAIN =
      brightness_envelope_gain(sample.value, state.brightness_envelope,
                               state.brightness_bottom, state.brightness_top);
  if (GAIN != 1.0f)
    color.color = color.color * GAIN;
  color.alpha *= sample.coverage *
                 hs::lerp(state.opacity_low, state.opacity_high, sample.value);
  return color;
}
/** @brief Colors a field sample from one immutable source endpoint.
 * @param sample Renderer field sample.
 * @param state Prepared colorizer.
 * @return Completed straight-alpha color. */
HS_HOT_FLASH_MEMBER inline Color4
apply_generated_palette(const FieldSample &sample,
                        const GeneratedPaletteState &state) {
  if (state.modulation && state.modulation->config.enabled &&
      state.mechanic.amount != 0.0f) {
    const auto DIRECTION = Pattern::SphereDirection::from(sample.sphere);
    return generated_palette_endpoint(sample, state, &DIRECTION);
  }
  return generated_palette_endpoint(sample, state, nullptr);
}

/** @brief Interpolates completed straight-alpha endpoint colors.
 * @param sample Field sample with shared column controls.
 * @param prepared Endpoint resources and interpolation fraction.
 * @return RGB using Q16 rounding and alpha using the original float fraction. */
HS_HOT_FLASH_MEMBER inline Color4 apply_temporal_generated_palette(
    const FieldSample &sample, const TemporalGeneratedPaletteState &prepared) {
  if (prepared.fraction <= 0.0f)
    return apply_generated_palette(sample, prepared.previous);
  if (prepared.fraction >= 1.0f)
    return apply_generated_palette(sample, prepared.current);
  const bool ACTIVE = (prepared.current.modulation &&
                       prepared.current.modulation->config.enabled &&
                       prepared.current.mechanic.amount != 0.0f) ||
                      (prepared.previous.modulation &&
                       prepared.previous.modulation->config.enabled &&
                       prepared.previous.mechanic.amount != 0.0f);
  const auto endpoint = [&](const Pattern::SphereDirection *direction) {
    const float SCALAR =
        ACTIVE && prepared.same_source
            ? Modulation::sample(*prepared.current.modulation, *direction)
            : 0.0f;
    const float *shared = ACTIVE && prepared.same_source ? &SCALAR : nullptr;
    const Color4 A = generated_palette_endpoint(sample, prepared.previous,
                                                direction, shared);
    const Color4 B =
        generated_palette_endpoint(sample, prepared.current, direction, shared);
    return Modulation::interpolate_colors(A, B, prepared.fraction,
                                          prepared.weight);
  };
  if (ACTIVE) {
    const auto DIRECTION = Pattern::SphereDirection::from(sample.sphere);
    return endpoint(&DIRECTION);
  }
  return endpoint(nullptr);
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

/**
 * @brief GeneratedPalette color policy.
 * @tparam State Provider with Binding and FrameState types and
 * mapping_weights(frame), mapping_frequency(frame), mapping_phase(frame),
 * oscillation_depth(frame), oscillation_phase(frame), palette(frame),
 * hue_mode(frame), hue_shift_amount(frame), hue_rotation(frame),
 * modulation(frame), mechanic(frame), brightness_envelope(frame), brightness_bottom(frame),
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
        {
          State::modulation(frame)
        } -> std::same_as<const Modulation::Prepared *>;
        { State::mechanic(frame) } -> std::same_as<ColorModulation::Params>;
        {
          State::brightness_envelope(frame)
        } -> std::same_as<BrightnessEnvelope>;
        { State::brightness_bottom(frame) } -> std::same_as<float>;
        { State::brightness_top(frame) } -> std::same_as<float>;
        { State::opacity_low(frame) } -> std::same_as<float>;
        { State::opacity_high(frame) } -> std::same_as<float>;
      };

  /// Whether the provider supplies resources for temporal interpolation.
  static constexpr bool TEMPORAL =
      requires(const typename State::Binding::FrameState &frame) {
        {
          State::temporal_resources(frame)
        } -> std::same_as<TemporalPaletteResources>;
      };
  /// Per-frame colorizer state, with temporal endpoints when supplied.
  using Prepared = std::conditional_t<TEMPORAL, TemporalGeneratedPaletteState,
                                      GeneratedPaletteState>;

  /**
   * @brief Reads every accessor once and folds the phase oscillation in.
   * @tparam FrameState Frame state of the provider's binding.
   * @param frame Current frame state.
   * @return This frame's colorizer state.
   */
  template <typename FrameState>
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    const GeneratedPaletteState current{
        State::mapping_weights(frame),
        State::mapping_frequency(frame),
        oscillated_mapping_offset(State::mapping_phase(frame),
                                  State::oscillation_depth(frame),
                                  State::oscillation_phase(frame)),
        &State::palette(frame),
        State::hue_mode(frame),
        State::hue_shift_amount(frame),
        State::hue_rotation(frame),
        State::modulation(frame),
        State::mechanic(frame),
        State::brightness_envelope(frame),
        State::brightness_bottom(frame),
        State::brightness_top(frame),
        State::opacity_low(frame),
        State::opacity_high(frame),
        [&]() {
          if constexpr (requires { State::palette_alpha(frame); })
            return State::palette_alpha(frame);
          else
            return -1.0f;
        }(),
        mechanic_kernel(State::mechanic(frame).kind)};
    if constexpr (TEMPORAL) {
      const auto RESOURCES = State::temporal_resources(frame);
      GeneratedPaletteState previous = current;
      previous.palette = RESOURCES.palette;
      previous.hue_rotation.data = RESOURCES.hue_rotation;
      previous.hue_rotation.linear_data = RESOURCES.hue_rotation_linear;
      previous.palette_alpha = RESOURCES.palette_alpha;
      const bool SAME_TOPOLOGY =
          current.modulation && RESOURCES.modulation
              ? Modulation::same_topology(current.modulation->config,
                                          RESOURCES.modulation->config)
              : current.modulation == RESOURCES.modulation;
      previous.modulation = RESOURCES.continuous_modulation && SAME_TOPOLOGY
                                ? RESOURCES.modulation
                                : current.modulation;
      const bool SAME_SOURCE =
          current.modulation && previous.modulation &&
          Modulation::same_field(*current.modulation, *previous.modulation);
      const bool IDENTITY = !current.modulation ||
                            !current.modulation->config.enabled ||
                            current.mechanic.amount == 0.0f;
      const float FRACTION = RESOURCES.same_palette && (IDENTITY || SAME_SOURCE)
                                 ? 1.0f
                                 : RESOURCES.fraction;
      return {current, previous, FRACTION, frac_to_q16(FRACTION), SAME_SOURCE};
    } else
      return current;
  }

  /**
   * @brief Colours one field sample.
   * @tparam FrameState Frame state of the provider's binding.
   * @param sample Field sample to colour.
   * @param prepared This frame's colorizer state.
   * @return Colour with temporally interpolated resources when supplied.
   */
  template <typename FrameState>
  HS_O3_FN static Color4 apply(const FieldSample &sample, const FrameState &,
                               const Prepared &prepared) {
    if constexpr (TEMPORAL)
      return apply_temporal_generated_palette(sample, prepared);
    else
      return apply_generated_palette(sample, prepared);
  }
};

} // namespace Color

} // namespace Pullback
