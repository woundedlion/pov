/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include "core/render/pullback/operators/model.h"
#include "core/render/pullback/operators/common.h"
#include "core/render/pullback/operators/field.h"
#include "core/render/pullback/operators/project.h"
#include "core/render/pullback/operators/sample.h"
#include "core/render/pullback/operators/sphere.h"
#include "core/render/pullback/operators/warp.h"
#include "core/render/pullback/stage.h"

/**
 * @file operators.h
 * @brief Concrete chain-interpreter operator models: erased adapters over the
 *        shared carrier kernels, with param-block-backed policy math. The
 *        per-family headers carry the family batches; this header aggregates
 *        them and keeps the camera rotation and colorize exit operators.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Parameter family of sphere.rotate.v2. */
struct RotateChainParams {
  float wander = 0.0f;    /**< Fraction of the walk delta absorbed per frame. */
  float spin_rate = 0.0f; /**< Per-frame spin about Y, in radians. */

  static constexpr auto FIELDS = std::array{
      Field<RotateChainParams>{"wander", &RotateChainParams::wander, "Wander",
                               0.0f, 1.0f, FieldCurve::LERP},
      Field<RotateChainParams>{"spin-speed", &RotateChainParams::spin_rate,
                               "Spin Speed", 0.0f, 0.05f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<RotateChainParams>());
static_assert(field_defaults_in_range<RotateChainParams>());

/** @brief SPHERE endomorphism: the wandering, spinning camera. */
struct Rotate : ValueStateModel<SpatialWalkState> {
  static constexpr const char *ID = "sphere.rotate.v2";
  static constexpr const char *NAME = "Rotate";
  using Input = SphereSample;
  using Output = SphereSample;
  using Params = RotateChainParams;
  /** @brief Per-frame camera rotation. */
  struct Prepared {
    math::Quaternion conjugate; ///< Inverse of the spin-times-wander rotation.
  };

  static void init(State &state, InstanceId id) {
    init_walk(state, static_cast<int32_t>(id.stable_hash));
  }
  static void advance(State &state, const Params &params) {
    advance_walk(state, params.wander, params.spin_rate);
  }
  static Prepared prepare(const FrameContext &, const Params &,
                          const State &state) {
    return {(math::make_rotation(math::Y_AXIS, state.spin_phase) * state.wander)
                .conjugate()};
  }
  static SphereSample run(const SphereSample &input, const FrameContext &,
                          const Params &, const Prepared &prepared) {
    return Kernel::rotate_dir(input, prepared.conjugate);
  }
};

/** @brief Which frame palette the colorize operator reads; indexes
    `FrameContext::palettes`. */
enum class PaletteMode : uint8_t {
  TRIADIC = 0,       ///< Triadic-harmony palette.
  COMPLEMENTARY = 1, ///< Complementary-harmony palette.
  ANALOGOUS = 2      ///< Analogous-harmony palette.
};
/// Hue-rotation source of the colorize operator.
using HueShiftMode = Color::HueMode;
/// Brightness envelope of the colorize operator.
using EnvelopeMode = Color::BrightnessEnvelope;

/// Wire ids of `PaletteMode`, by enumerator value.
inline constexpr const char *PALETTE_MODE_IDS[] = {"triadic", "complementary",
                                                   "analogous"};
static_assert(std::size(PALETTE_MODE_IDS) ==
              static_cast<size_t>(PaletteMode::ANALOGOUS) + 1);
static_assert(std::size(PALETTE_MODE_IDS) == FrameContext{}.palettes.size(),
              "frame palette storage must cover every palette mode");
/// Wire ids of `Color::PaletteMapping`, by enumerator value.
inline constexpr const char *PALETTE_MAPPING_IDS[] = {"cup", "bell", "linear",
                                                      "reverse"};
static_assert(std::size(PALETTE_MAPPING_IDS) ==
              static_cast<size_t>(Color::PaletteMapping::REVERSE) + 1);
/// Wire ids of `HueShiftMode`, by enumerator value.
inline constexpr const char *HUE_SHIFT_MODE_IDS[] = {"none", "path-length"};
static_assert(std::size(HUE_SHIFT_MODE_IDS) ==
              static_cast<size_t>(HueShiftMode::PATH_LENGTH) + 1);
/// Wire ids of `EnvelopeMode`, by enumerator value.
inline constexpr const char *BRIGHTNESS_ENVELOPE_IDS[] = {
    "none", "cup", "bell", "ascending", "descending"};
static_assert(std::size(BRIGHTNESS_ENVELOPE_IDS) ==
              static_cast<size_t>(EnvelopeMode::DESCENDING) + 1);

/** @brief Parameter family of colorize.generated-palette.v3. */
struct GeneratedPaletteParams : Color::ColorControls {
  /// A `PaletteMode` value.
  uint8_t palette_mode = static_cast<uint8_t>(PaletteMode::TRIADIC);
  /// A `Color::PaletteMapping` value.
  uint8_t mapping_mode = static_cast<uint8_t>(Color::PaletteMapping::LINEAR);
  /// A `HueShiftMode` value.
  uint8_t hue_mode = static_cast<uint8_t>(HueShiftMode::NONE);
  /// An `EnvelopeMode` value.
  uint8_t envelope_mode = static_cast<uint8_t>(EnvelopeMode::NONE);

  /// Modulation enable values.
  static constexpr const char *ENABLED_IDS[] = {"disabled", "enabled"};
  /// Pattern source values.
  static constexpr const char *SOURCE_IDS[] = {
      "noise", "directional-sine", "sphere-interference", "sphere-smooth-grid"};
  /// Pattern evaluation values.
  static constexpr const char *EVALUATION_IDS[] = {"direct", "baked-cube"};
  /// Color mechanic values.
  static constexpr const char *MECHANIC_IDS[] = {
      "hue-rotation", "positive-lightness-lift", "opacity"};
  /// Noise octave formula values.
  static constexpr const char *BASIS_IDS[] = {"simplex", "fbm3", "ridged3"};
  /// Noise sampling policy values.
  static constexpr const char *VARIANT_IDS[] = {"configured-single",
                                                "octave-contour"};
  /// Modulation parameter schema appended to the palette controls.
  static constexpr auto EXTRA_SCHEMA = [] {
    std::array<ParamFieldInfo, 38> out{};
    const char *ids[] = {"modulation-amount",
                         "pattern-scale",
                         "pattern-speed",
                         "pattern-contrast",
                         "sine-axis-x",
                         "sine-axis-y",
                         "sine-axis-z",
                         "sine-frequency",
                         "sine-offset",
                         "sine-speed",
                         "interference-axis-x",
                         "interference-axis-y",
                         "interference-axis-z",
                         "interference-frequency",
                         "interference-offset",
                         "interference-speed",
                         "interference-axis-b-x",
                         "interference-axis-b-y",
                         "interference-axis-b-z",
                         "interference-frequency-b",
                         "interference-offset-b",
                         "grid-axis-x",
                         "grid-axis-y",
                         "grid-axis-z",
                         "grid-frequency",
                         "grid-offset",
                         "grid-speed",
                         "grid-axis-b-x",
                         "grid-axis-b-y",
                         "grid-axis-b-z",
                         "grid-frequency-b",
                         "grid-offset-b"};
    for (size_t i = 0; i < 32; ++i)
      out[i] = {ids[i], nullptr, -64.0f, 64.0f,   0.0f,    FieldCurve::LERP,
                false,  0,       0,      nullptr, nullptr, 0};
    out[0].min = -4.0f;
    out[0].max = 4.0f;
    out[1].min = 1.0f / 64.0f;
    out[1].max = 64.0f;
    out[1].def = 1.0f;
    out[1].curve = FieldCurve::LOG_POSITIVE;
    out[3].min = 0.0f;
    out[3].max = 8.0f;
    for (size_t i :
         {size_t{5}, size_t{7}, size_t{11}, size_t{13}, size_t{16}, size_t{19},
          size_t{22}, size_t{24}, size_t{27}, size_t{30}})
      out[i].def = 1.0f;
    for (size_t i : {size_t{7}, size_t{13}, size_t{19}, size_t{24}, size_t{30}})
      out[i].min = 0.0f;
    out[32] = {
        "modulation-enabled", nullptr, 0, 0, 0, FieldCurve::SNAP, true, 2, 1,
        ENABLED_IDS,          nullptr, 0};
    out[33] = {
        "modulation-source", nullptr, 0, 0, 0, FieldCurve::SNAP, true, 4, 0,
        SOURCE_IDS,          nullptr, 0};
    out[34] = {
        "modulation-evaluation", nullptr, 0, 0, 0, FieldCurve::SNAP, true, 2, 1,
        EVALUATION_IDS,          nullptr, 0};
    out[35] = {
        "modulation-mechanic", nullptr, 0, 0, 0, FieldCurve::SNAP, true, 3, 0,
        MECHANIC_IDS,          nullptr, 0};
    out[36] = {"pattern-noise-basis",
               nullptr,
               0,
               0,
               0,
               FieldCurve::SNAP,
               true,
               3,
               0,
               BASIS_IDS,
               nullptr,
               0};
    out[37] = {"pattern-noise-variant",
               nullptr,
               0,
               0,
               0,
               FieldCurve::SNAP,
               true,
               2,
               0,
               VARIANT_IDS,
               nullptr,
               0};
    return out;
  }();
  /** @brief Mutable canonical control addresses.
   * @param index Index in the appended schema.
   * @return Canonical configuration field address. */
  void *extra_address(size_t index) {
    auto &c = modulation;
    float *addresses[] = {&c.mechanic.amount,
                          &c.noise.scale,
                          &c.noise.speed,
                          &c.noise.contrast,
                          &c.sine.axis.x,
                          &c.sine.axis.y,
                          &c.sine.axis.z,
                          &c.sine.frequency,
                          &c.sine.offset,
                          &c.sine.speed,
                          &c.interference.axis.x,
                          &c.interference.axis.y,
                          &c.interference.axis.z,
                          &c.interference.frequency,
                          &c.interference.offset,
                          &c.interference.speed,
                          &c.interference.axis_b.x,
                          &c.interference.axis_b.y,
                          &c.interference.axis_b.z,
                          &c.interference.frequency_b,
                          &c.interference.offset_b,
                          &c.grid.axis.x,
                          &c.grid.axis.y,
                          &c.grid.axis.z,
                          &c.grid.frequency,
                          &c.grid.offset,
                          &c.grid.speed,
                          &c.grid.axis_b.x,
                          &c.grid.axis_b.y,
                          &c.grid.axis_b.z,
                          &c.grid.frequency_b,
                          &c.grid.offset_b};
    if (index < std::size(addresses))
      return addresses[index];
    switch (index) {
    case 32:
      return &c.enabled;
    case 33:
      return &c.source;
    case 34:
      return &c.evaluation;
    case 35:
      return &c.mechanic.kind;
    case 36:
      return &c.noise.basis;
    case 37:
      return &c.noise.variant;
    }
    HS_CHECK(false, "invalid modulation field index");
    return nullptr;
  }
  /** @brief Interpolates mechanic controls only across matching recipe topology.
   * @param out Interpolated parameters.
   * @param a Previous parameters.
   * @param b Current parameters.
   * @param t Current weight. */
  static void sample_extra(GeneratedPaletteParams &out,
                           const GeneratedPaletteParams &a,
                           const GeneratedPaletteParams &b, float t) {
    out.modulation = b.modulation;
    if (Modulation::same_topology(a.modulation, b.modulation))
      out.modulation.mechanic.amount = hs::lerp(
          a.modulation.mechanic.amount, b.modulation.mechanic.amount, t);
  }

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<GeneratedPaletteParams>{
          "palette-mode", &GeneratedPaletteParams::palette_mode,
          PALETTE_MODE_IDS, static_cast<uint8_t>(PaletteMode::TRIADIC)},
      TopologyField<GeneratedPaletteParams>{
          "palette-mapping", &GeneratedPaletteParams::mapping_mode,
          PALETTE_MAPPING_IDS,
          static_cast<uint8_t>(Color::PaletteMapping::LINEAR)},
      TopologyField<GeneratedPaletteParams>{
          "hue-shift-mode", &GeneratedPaletteParams::hue_mode,
          HUE_SHIFT_MODE_IDS, static_cast<uint8_t>(HueShiftMode::NONE)},
      TopologyField<GeneratedPaletteParams>{
          "brightness-envelope", &GeneratedPaletteParams::envelope_mode,
          BRIGHTNESS_ENVELOPE_IDS, static_cast<uint8_t>(EnvelopeMode::NONE)},
  };
};
static_assert(field_ids_unique<GeneratedPaletteParams>());
static_assert(field_defaults_in_range<GeneratedPaletteParams>());
static_assert(appended_block_size_matches<GeneratedPaletteParams,
                                          Color::ColorControls, 4>(),
              "appended parameter block must have the expected rounded size");

/** @brief Per-frame color phase clocks. */
struct ColorClockState {
  float oscillation_phase = 0.0f; ///< Mapping oscillation in turns.
  Modulation::State modulation;   ///< Selected source state.
};

/** @brief Prepares the palette mapping clock.
 * @param from Previous state.
 * @param to Current state.
 * @return Prepared oscillation interval. */
inline Spacetime::ClockInterval temporal_interval(const ColorClockState &from,
                                                  const ColorClockState &to) {
  return {from.oscillation_phase, to.oscillation_phase, 1.0f};
}

/**
 * @brief Samples the mapping oscillation clock at a frame fraction.
 * @tparam Params Operator param block; unused.
 * @param out Copy of the current state; receives the samples.
 * @param interval Prepared mapping clock.
 * @param fraction Frame fraction in [0, 1].
 */
template <typename Params>
inline void temporal_sample(ColorClockState &out,
                            const Spacetime::ClockInterval &interval,
                            const Params &, float fraction) {
  out.oscillation_phase = interval.sample(fraction);
}

/**
 * @brief Resolves a generated-palette param block into the per-frame colour
 *        state; mode fields out of range fail an `HS_CHECK`.
 * @param ctx Frame context supplying the palettes and hue LUTs.
 * @param params Operator parameters.
 * @param state Colour clocks.
 * @return The prepared colour state for this frame.
 */
inline Color::GeneratedPaletteState
prepare_generated_palette(const FrameContext &ctx,
                          const GeneratedPaletteParams &params,
                          const ColorClockState &state) {
  HS_CHECK(params.hue_mode <= static_cast<uint8_t>(HueShiftMode::PATH_LENGTH),
           "colorize.generated-palette: invalid hue shift mode");
  HS_CHECK(params.palette_mode < ctx.palettes.size(),
           "colorize.generated-palette: invalid palette mode");
  HS_CHECK(params.mapping_mode <=
               static_cast<uint8_t>(Color::PaletteMapping::REVERSE),
           "colorize.generated-palette: invalid palette mapping");
  HS_CHECK(params.envelope_mode <=
               static_cast<uint8_t>(EnvelopeMode::DESCENDING),
           "colorize.generated-palette: invalid brightness envelope");
  const Color::HueMode mode = static_cast<HueShiftMode>(params.hue_mode);
  const bool rotation_active = Color::hue_rotation_active(mode, params) &&
                               ctx.hue_rotation_lut != nullptr;
  const BakedPalette *palette = ctx.palettes[params.palette_mode];
  HS_CHECK(palette != nullptr,
           "colorize.generated-palette: frame context carries no palette");
  return {Color::PaletteMappingWeights::single(
              static_cast<Color::PaletteMapping>(params.mapping_mode)),
          params.mapping_frequency,
          Color::oscillated_mapping_offset(params.mapping_phase,
                                           params.phase_oscillation_depth,
                                           state.oscillation_phase),
          palette,
          mode,
          params.hue_shift_amount,
          {ctx.hue_rotation_lut, rotation_active},
          ctx.modulation,
          params.modulation.mechanic,
          static_cast<Color::BrightnessEnvelope>(params.envelope_mode),
          params.brightness_bottom,
          params.brightness_top,
          params.opacity_low,
          params.opacity_high,
          Color::mechanic_kernel(params.modulation.mechanic.kind)};
}

/** @brief FIELD-to-COLOR crossing with generated-palette clocks and endpoints. */
struct ColorizeGeneratedPaletteV3 : ValueStateModel<ColorClockState> {
  static constexpr const char *ID = "colorize.generated-palette.v3";
  static constexpr const char *NAME = "Generated Palette";
  using Input = FieldSample;
  using Output = Color4;
  using Params = GeneratedPaletteParams;
  using Prepared = Color::GeneratedPaletteState;

  /** @brief Validates source and mechanic recipes even when disabled.
   * @param params Canonical controls.
   * @return Error message, or null when admitted. */
  static const char *validate(const Params &params) {
    return Modulation::valid(params.modulation) &&
                   !(params.modulation.enabled &&
                     params.hue_mode ==
                         static_cast<uint8_t>(HueShiftMode::PATH_LENGTH))
               ? nullptr
               : Modulation::INVALID_CONFIGURATION;
  }
  /// Interpolation contract for source and mechanic controls.
  static constexpr const char *ADMISSIBILITY_CONVEXITY =
      "Source controls hold endpoint recipes; mechanic strength interpolates within matching recipe topology.";
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  static constexpr auto METRICS = Color::GENERATED_PALETTE_METRICS;

  static void advance(State &state, const Params &params) {
    state.oscillation_phase =
        math::wrap_t(state.oscillation_phase + params.phase_oscillation_speed);
    Modulation::advance(state.modulation, params.modulation);
  }
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    return prepare_generated_palette(ctx, params, state);
  }
  static Color4 run(const FieldSample &input, const FrameContext &,
                    const Params &, const Prepared &prepared) {
    return Color::apply_generated_palette(input, prepared);
  }
  /** @brief Interpolates completed colors using shared direction admission.
   * @param input Renderer field sample.
   * @param current Current prepared endpoint.
   * @param previous Retained prepared endpoint.
   * @param fraction Current endpoint interpolation weight.
   * @return Interpolated straight-alpha color. */
  static Color4 run_temporal(const FieldSample &input, const Prepared &current,
                             const Prepared &previous, float fraction) {
    const bool SAME_SOURCE =
        current.modulation && previous.modulation &&
        Modulation::same_field(*current.modulation, *previous.modulation);
    return Color::apply_temporal_generated_palette(
        input,
        {current, previous, fraction, frac_to_q16(fraction), SAME_SOURCE});
  }
  /** @brief Holds modulation topology at its current endpoint.
   * @param field Parameter field index.
   * @return Whether the field may change without snapping the column interval. */
  static constexpr bool temporal_hold_field(uint16_t field) {
    return field >= Params::FIELDS.size() + Params::TOPOLOGY.size();
  }
};

} // namespace Op

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
