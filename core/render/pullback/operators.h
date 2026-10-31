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
/// Wire ids of `EnvelopeMode`, by enumerator value.
inline constexpr const char *BRIGHTNESS_ENVELOPE_IDS[] = {
    "none", "cup", "bell", "ascending", "descending"};
static_assert(std::size(BRIGHTNESS_ENVELOPE_IDS) ==
              static_cast<size_t>(EnvelopeMode::DESCENDING) + 1);

/** @brief Palette source selected by the workbench color stage. */
enum class PaletteType : uint8_t {
  CYCLING_GENERATED, ///< Animated legacy harmony bank.
  GENERATED,         ///< Authored perceptual palette.
  PROCEDURAL,        ///< Named cosine palette.
  GRADIENT,          ///< Three-stop perceptual gradient.
  MUTATING           ///< Interpolated cosine coefficients.
};
/// Wire ids of palette source families.
inline constexpr const char *PALETTE_TYPE_IDS[] = {
    "cycling-generated", "generated", "procedural", "gradient", "mutating"};
/// Wire ids of perceptual harmony relationships.
inline constexpr const char *PALETTE_HARMONY_IDS[] = {
    "monochromatic",       "analogous", "accented-analogous", "complementary",
    "split-complementary", "triadic",   "tetradic",           "square"};
static_assert(std::size(PALETTE_HARMONY_IDS) ==
              static_cast<size_t>(PaletteHarmony::SQUARE) + 1);
/// Wire ids of lightness profiles supported by the palette controls.
inline constexpr const char *PALETTE_LIGHTNESS_CURVE_IDS[] = {
    "constant", "ascending", "descending", "bell", "cup"};
/// Wire ids of the named procedural palette roster.
inline constexpr const char *PALETTE_NAME_IDS[] = {
    "dark-rainbow", "blood-stream",   "vintage-sunset", "rich-sunset",
    "undersea",     "late-sunset",    "mango-peel",     "ice-melt",
    "lemon-lime",   "algae",          "embers",         "fire-glow",
    "dark-primary", "mauve-fade",     "lavender-lake",  "desert-rose",
    "bruised-moss", "bruised-banana", "bright-sunrise", "fire-and-ice",
    "peach-pop",    "tidal-jade",     "blue-lagoon",    "orange-crush",
    "plum-sunrise", "coral-blue",     "bruised-mango"};
/// Authored palette source controls, persisted with the color operator.
struct PaletteSelection {
  float palette_hue = 0.0f;             ///< Base hue in turns.
  float palette_spread = 0.07f;         ///< Harmony spread in turns.
  float palette_lightness = 0.62f;      ///< Lightness midpoint.
  float palette_lightness_range = 0.0f; ///< Lightness profile range.
  float palette_mutation = 0.0f;        ///< Coefficient interpolation weight.
  float gradient_start_r = 1.0f;        ///< First stop red, sRGB.
  float gradient_start_g = 0.15f;       ///< First stop green, sRGB.
  float gradient_start_b = 0.0f;        ///< First stop blue, sRGB.
  float gradient_middle_r = 0.1f;       ///< Middle stop red, sRGB.
  float gradient_middle_g = 1.0f;       ///< Middle stop green, sRGB.
  float gradient_middle_b = 0.4f;       ///< Middle stop blue, sRGB.
  float gradient_end_r = 0.1f;          ///< Last stop red, sRGB.
  float gradient_end_g = 0.2f;          ///< Last stop green, sRGB.
  float gradient_end_b = 1.0f;          ///< Last stop blue, sRGB.
  uint8_t palette_type = 0;             ///< A PaletteType ordinal.
  /// Perceptual hue relationship.
  uint8_t palette_harmony = static_cast<uint8_t>(PaletteHarmony::TRIADIC);
  uint8_t palette_lightness_curve = 0; ///< Lightness axis profile.
  uint8_t palette_name = 0;            ///< First named cosine palette.
  uint8_t palette_name_b = 19;         ///< Second named cosine palette.
  std::array<uint8_t, 3> reserved{};   ///< Stable serialized base layout.
  constexpr bool operator==(const PaletteSelection &) const = default;
};
static_assert(sizeof(PaletteSelection) == 64);
/// Activation gate for authored perceptual palette controls.
inline constexpr TopologyGate GENERATED_PALETTE_GATE{
    "palette-type", live_values(PaletteType::GENERATED)};
/// Activation gate for cosine mutation strength.
inline constexpr TopologyGate MUTATING_PALETTE_GATE{
    "palette-type", live_values(PaletteType::MUTATING)};
/// Activation gate for RGB gradient stops.
inline constexpr TopologyGate GRADIENT_PALETTE_GATE{
    "palette-type", live_values(PaletteType::GRADIENT)};

/** @brief Palette-only controls of colorize.generated-palette.v4. */
struct GeneratedPaletteParams : Color::ColorControls, PaletteSelection {
  uint8_t palette_mode =
      static_cast<uint8_t>(PaletteMode::TRIADIC); ///< Palette harmony.
  uint8_t mapping_mode =
      static_cast<uint8_t>(Color::PaletteMapping::LINEAR); ///< Mapping curve.
  uint8_t envelope_mode =
      static_cast<uint8_t>(EnvelopeMode::NONE); ///< Brightness shape.
  static constexpr auto FIELDS = concat_fields<GeneratedPaletteParams>(
      Color::ColorControls::FIELDS,
      std::array{
          Field<GeneratedPaletteParams>{
              "palette-hue", &GeneratedPaletteParams::palette_hue, nullptr,
              0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GENERATED_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "palette-spread", &GeneratedPaletteParams::palette_spread,
              nullptr, 0.0f, 0.25f, FieldCurve::LERP, FieldGate::ALWAYS,
              GENERATED_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "palette-lightness", &GeneratedPaletteParams::palette_lightness,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GENERATED_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "palette-lightness-range",
              &GeneratedPaletteParams::palette_lightness_range, nullptr, 0.0f,
              1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GENERATED_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "palette-mutation", &GeneratedPaletteParams::palette_mutation,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              MUTATING_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-start-r", &GeneratedPaletteParams::gradient_start_r,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-start-g", &GeneratedPaletteParams::gradient_start_g,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-start-b", &GeneratedPaletteParams::gradient_start_b,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-middle-r", &GeneratedPaletteParams::gradient_middle_r,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-middle-g", &GeneratedPaletteParams::gradient_middle_g,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-middle-b", &GeneratedPaletteParams::gradient_middle_b,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-end-r", &GeneratedPaletteParams::gradient_end_r,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-end-g", &GeneratedPaletteParams::gradient_end_g,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
          Field<GeneratedPaletteParams>{
              "gradient-end-b", &GeneratedPaletteParams::gradient_end_b,
              nullptr, 0.0f, 1.0f, FieldCurve::LERP, FieldGate::ALWAYS,
              GRADIENT_PALETTE_GATE},
      });

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<GeneratedPaletteParams>{
          "palette-mode", &GeneratedPaletteParams::palette_mode,
          PALETTE_MODE_IDS, 0},
      TopologyField<GeneratedPaletteParams>{
          "palette-mapping", &GeneratedPaletteParams::mapping_mode,
          PALETTE_MAPPING_IDS, 2},
      TopologyField<GeneratedPaletteParams>{
          "brightness-envelope", &GeneratedPaletteParams::envelope_mode,
          BRIGHTNESS_ENVELOPE_IDS, 0},
      TopologyField<GeneratedPaletteParams>{
          "palette-type", &GeneratedPaletteParams::palette_type,
          PALETTE_TYPE_IDS, 0},
      TopologyField<GeneratedPaletteParams>{
          "palette-harmony", &GeneratedPaletteParams::palette_harmony,
          PALETTE_HARMONY_IDS, static_cast<uint8_t>(PaletteHarmony::TRIADIC)},
      TopologyField<GeneratedPaletteParams>{
          "palette-lightness-curve",
          &GeneratedPaletteParams::palette_lightness_curve,
          PALETTE_LIGHTNESS_CURVE_IDS, 0},
      TopologyField<GeneratedPaletteParams>{
          "palette-name", &GeneratedPaletteParams::palette_name,
          PALETTE_NAME_IDS, 0},
      TopologyField<GeneratedPaletteParams>{
          "palette-name-b", &GeneratedPaletteParams::palette_name_b,
          PALETTE_NAME_IDS, 19},
  }; ///< Structural controls.
};
static_assert(field_ids_unique<GeneratedPaletteParams>());
static_assert(field_defaults_in_range<GeneratedPaletteParams>());

/** @brief Controls of one color.modulate.v1 operator. */
struct ModulationParams {
  ColorModulation::Config modulation = [] {
    ColorModulation::Config value;
    value.seed = HUE_NOISE_SEED;
    return value;
  }(); ///< Canonical source and mechanic recipe.
  static constexpr std::array<Field<ModulationParams>, 0>
      FIELDS{}; ///< Nested recipe controls.
  static constexpr std::array<TopologyField<ModulationParams>, 0>
      TOPOLOGY{}; ///< Nested recipe selectors.
  /// Modulation enable values.
  static constexpr const char *ENABLED_IDS[] = {"disabled", "enabled"};
  /// Pattern source values.
  static constexpr const char *SOURCE_IDS[] = {
      "noise", "directional-sine", "sphere-interference", "sphere-smooth-grid",
      "displacement"};
  static_assert(std::size(SOURCE_IDS) ==
                static_cast<size_t>(Pattern::SourceKind::DISPLACEMENT) + 1);
  /// Pattern evaluation values.
  static constexpr const char *EVALUATION_IDS[] = {"direct", "baked-cube"};
  static_assert(std::size(EVALUATION_IDS) ==
                static_cast<size_t>(Pattern::Evaluation::BAKED_CUBE) + 1);
  /// Color mechanic values.
  static constexpr const char *MECHANIC_IDS[] = {
      "hue-rotation", "positive-lightness-lift", "opacity"};
  static_assert(std::size(MECHANIC_IDS) ==
                static_cast<size_t>(ColorModulation::Mechanic::OPACITY) + 1);
  /// Noise octave formula values.
  static constexpr auto &BASIS_IDS = NOISE_BASIS_IDS;
  /// Noise sampling policy values.
  static constexpr const char *VARIANT_IDS[] = {"configured-single",
                                                "octave-contour"};
  static_assert(std::size(VARIANT_IDS) ==
                static_cast<size_t>(Pattern::NoiseVariant::OCTAVE_CONTOUR) + 1);
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
    out[2].min = -Pattern::MAX_NOISE_SPEED;
    out[2].max = Pattern::MAX_NOISE_SPEED;
    for (size_t i : {size_t{9}, size_t{15}, size_t{26}}) {
      out[i].min = -Pattern::MAX_ANALYTIC_SPEED;
      out[i].max = Pattern::MAX_ANALYTIC_SPEED;
    }
    out[3].min = 0.0f;
    out[3].max = 8.0f;
    for (size_t i :
         {size_t{5}, size_t{7}, size_t{11}, size_t{13}, size_t{16}, size_t{19},
          size_t{22}, size_t{24}, size_t{27}, size_t{30}})
      out[i].def = 1.0f;
    for (size_t i : {size_t{7}, size_t{13}, size_t{19}, size_t{24}, size_t{30}})
      out[i].min = 0.0f;
    out[32] = {"modulation-enabled",
               nullptr,
               0,
               0,
               0,
               FieldCurve::SNAP,
               true,
               std::size(ENABLED_IDS),
               1,
               ENABLED_IDS,
               nullptr,
               0};
    out[33] = {"modulation-source",
               nullptr,
               0,
               0,
               0,
               FieldCurve::SNAP,
               true,
               std::size(SOURCE_IDS),
               0,
               SOURCE_IDS,
               nullptr,
               0};
    out[34] = {"modulation-evaluation",
               nullptr,
               0,
               0,
               0,
               FieldCurve::SNAP,
               true,
               std::size(EVALUATION_IDS),
               1,
               EVALUATION_IDS,
               nullptr,
               0};
    out[35] = {"modulation-mechanic",
               nullptr,
               0,
               0,
               0,
               FieldCurve::SNAP,
               true,
               std::size(MECHANIC_IDS),
               0,
               MECHANIC_IDS,
               nullptr,
               0};
    out[36] = {"pattern-noise-basis",
               nullptr,
               0,
               0,
               0,
               FieldCurve::SNAP,
               true,
               std::size(BASIS_IDS),
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
               std::size(VARIANT_IDS),
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
  static void sample_extra(ModulationParams &out, const ModulationParams &a,
                           const ModulationParams &b, float t) {
    out.modulation = b.modulation;
    if (Modulation::same_topology(a.modulation, b.modulation))
      out.modulation.mechanic.amount = hs::lerp(
          a.modulation.mechanic.amount, b.modulation.mechanic.amount, t);
  }

  /** @brief Identifies modulation fields excluded from interval eligibility.
   * @param field Parameter schema index.
   * @return Whether the field belongs to modulation. */
  static constexpr bool temporal_hold_field(uint16_t field) {
    (void)field;
    return true;
  }
};
static_assert(field_ids_unique<ModulationParams>());

/** @brief Per-frame color phase clocks. */
struct ColorClockState {
  float oscillation_phase = 0.0f; ///< Mapping oscillation in turns.
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

/** @brief Prepares palette-only immutable state.
 * @param ctx Palette resources.
 * @param params Palette controls.
 * @param state Oscillation clock.
 * @return Palette endpoint state. */
inline Color::GeneratedPaletteState
prepare_generated_palette(const FrameContext &ctx,
                          const GeneratedPaletteParams &params,
                          const ColorClockState &state) {
  HS_CHECK(params.palette_mode < ctx.palettes.size(), "invalid palette mode");
  HS_CHECK(params.mapping_mode <=
               static_cast<uint8_t>(Color::PaletteMapping::REVERSE),
           "invalid palette mapping");
  HS_CHECK(params.envelope_mode <=
               static_cast<uint8_t>(EnvelopeMode::DESCENDING),
           "invalid brightness envelope");
  HS_CHECK(ctx.palettes[params.palette_mode] != nullptr,
           "missing palette endpoint");
  return {Color::PaletteMappingWeights::single(
              static_cast<Color::PaletteMapping>(params.mapping_mode)),
          params.mapping_frequency,
          Color::oscillated_mapping_offset(params.mapping_phase,
                                           params.phase_oscillation_depth,
                                           state.oscillation_phase),
          ctx.palettes[params.palette_mode],
          {ctx.hue_rotation_lut, ctx.hue_rotation_lut != nullptr},
          static_cast<EnvelopeMode>(params.envelope_mode),
          params.brightness_bottom,
          params.brightness_top,
          params.opacity_low,
          params.opacity_high};
}

/** @brief FIELD-to-COLOR generated-palette crossing. */
struct ColorizeGeneratedPaletteV4 : ValueStateModel<ColorClockState> {
  static constexpr const char *ID = "colorize.generated-palette.v4";
  static constexpr const char *NAME = "Generated Palette";
  using Input = FieldSample;
  using Output = ColorSample;
  using Params = GeneratedPaletteParams;
  using Prepared = Color::GeneratedPaletteState;
  static void advance(State &state, const Params &params) {
    state.oscillation_phase =
        math::wrap_t(state.oscillation_phase + params.phase_oscillation_speed);
  }
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    return prepare_generated_palette(ctx, params, state);
  }
  static ColorSample run(const FieldSample &input, const FrameContext &,
                         const Params &, const Prepared &prepared) {
    return Color::apply_generated_palette(input, prepared);
  }
};

/** @brief Clock and immutable resource key of one modulation node. */
struct ModulationClockState {
  std::string_view key;         ///< Exact instance label.
  Modulation::State modulation; ///< Source clock.
};

/** @brief COLOR-to-COLOR ordered modulation. */
struct ModulateColorV1 : ValueStateModel<ModulationClockState> {
  static constexpr const char *ID = "color.modulate.v1";
  static constexpr const char *NAME = "Modulation";
  using Input = ColorSample;
  using Output = ColorSample;
  using Params = ModulationParams;
  using Prepared = Color::ModulationState;
  static constexpr bool APPROXIMATE = true;
  static constexpr ApproximationOracleId ORACLE =
      ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  static constexpr auto METRICS = Color::GENERATED_PALETTE_METRICS;
  static void init(State &state, InstanceId id) { state.key = id.instance; }
  /** @brief Copies the source clock and binds its destination instance.
   * @param dst Destination state.
   * @param src Source state.
   * @param id Destination instance identity.
   * @return Successful migration status. */
  static Status migrate(State &dst, const State &src, InstanceId id) {
    dst = src;
    dst.key = id.instance;
    return Status::OK;
  }
  static void advance(State &state, const Params &params) {
    Modulation::advance(state.modulation, params.modulation);
  }
  /** @brief Checks modulation configuration admission.
   * @param params Modulation controls.
   * @return Null when admitted, otherwise a diagnostic. */
  static const char *validate(const Params &params) {
    return Modulation::admitted(params.modulation)
               ? nullptr
               : Modulation::INVALID_CONFIGURATION;
  }
  /// Admission rule for interpolation between source recipes.
  static constexpr const char *ADMISSIBILITY_CONVEXITY =
      "Source recipes hold endpoint topology; mechanic strength interpolates within matching topology.";
  /** @brief Tests whether a field holds its temporal endpoint.
   * @param field Parameter field index.
   * @return Whether the field is held. */
  static constexpr bool temporal_hold_field(uint16_t field) {
    return Params::temporal_hold_field(field);
  }
  static Prepared prepare(const FrameContext &ctx, const Params &params,
                          const State &state) {
    const Modulation::Prepared *source = nullptr;
    for (const auto &endpoint : ctx.modulations)
      if (endpoint.key == state.key) {
        source = endpoint.prepared;
        break;
      }
    HS_CHECK(source ||
                 params.modulation.source ==
                     Pattern::SourceKind::DISPLACEMENT ||
                 !params.modulation.enabled ||
                 params.modulation.mechanic.amount == 0.0f,
             "missing modulation source endpoint");
    return {source, source, params.modulation.source,
            params.modulation.mechanic, params.modulation.enabled};
  }
  static ColorSample run(const ColorSample &input, const FrameContext &,
                         const Params &, const Prepared &prepared) {
    return Color::apply_modulation(input, prepared);
  }
};

} // namespace Op
} // namespace Interp
} // namespace Pullback
#endif
