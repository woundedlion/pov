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
#include "core/render/pullback/operators/planar.h"
#include "core/render/pullback/operators/sample.h"
#include "core/render/pullback/operators/sphere.h"
#include "core/render/pullback/operators/tangent.h"
#include "core/render/pullback/operators/warp.h"
#include "core/render/pullback/stage.h"

/**
 * @file operators.h
 * @brief Concrete chain-interpreter operator models: erased adapters over the
 *        shared carrier kernels, with param-block-backed policy math. The
 *        per-family headers carry the family batches; this header aggregates
 *        them and keeps the camera rotation, colorize and modulation operators.
 */

namespace Pullback {

namespace Interp {

namespace Op {

/** @brief Whether a source morph shares the affine's effective lattice.
 * @param from Outgoing source parameters.
 * @param to Incoming source parameters.
 * @param affine_from Outgoing affine parameters.
 * @param affine_to Incoming affine parameters.
 * @return True for one kind, an unused lattice, or equal effective lattices. */
inline bool admits_planar_transition(const PlanarSampleParams &from,
                                     const PlanarSampleParams &to,
                                     const AffineWarpParams &affine_from,
                                     const AffineWarpParams &affine_to) {
  const auto consumes_lattice = [](const AffineWarpParams &affine) {
    return affine.path_length != 0 || std::round(affine.translation_x) != 0 ||
           std::round(affine.translation_y) != 0;
  };
  return from.source == to.source ||
         (!consumes_lattice(affine_from) && !consumes_lattice(affine_to)) ||
         same_effective_lattice(SamplePlanar::translation_lattice(from),
                                affine_from.lattice_period,
                                SamplePlanar::translation_lattice(to),
                                affine_to.lattice_period);
}

/** @brief Parameter family of sphere.rotate.v2. */
using RotateChainParams = SphereAnimation::RotationParams;

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
    advance_walk(state, params.wander, params.spin_speed);
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

/// Which frame palette the colorize operator reads.
using PaletteMode = Color::PaletteMode;
/// Brightness envelope of the colorize operator.
using EnvelopeMode = Color::BrightnessEnvelope;

/// Wire ids of `PaletteMode`, by enumerator value.
inline constexpr auto &PALETTE_MODE_IDS = Color::PALETTE_MODE_IDS;
static_assert(std::size(PALETTE_MODE_IDS) == FrameContext{}.palettes.size(),
              "frame palette storage must cover every palette mode");
/// Wire ids of `Color::PaletteMapping`, by enumerator value.
inline constexpr const char *PALETTE_MAPPING_IDS[] = {"cup", "bell", "linear",
                                                      "reverse"};
static_assert(std::size(PALETTE_MAPPING_IDS) ==
              static_cast<size_t>(Color::PaletteMapping::REVERSE) + 1);
/// Wire ids of `EnvelopeMode`, by enumerator value.
inline constexpr auto &BRIGHTNESS_ENVELOPE_IDS = Color::BRIGHTNESS_ENVELOPE_IDS;

/// Authored palette source family.
using PaletteType = ::AuthoredPaletteType;
/// Wire ids of palette source families.
inline constexpr auto &PALETTE_TYPE_IDS = ::AUTHORED_PALETTE_TYPE_IDS;
/// Wire ids of perceptual harmony relationships.
inline constexpr auto &PALETTE_HARMONY_IDS = ::AUTHORED_PALETTE_HARMONY_IDS;
static_assert(std::size(PALETTE_HARMONY_IDS) ==
              static_cast<size_t>(PaletteHarmony::SQUARE) + 1);
/// Wire ids of lightness profiles supported by the palette controls.
inline constexpr auto &PALETTE_LIGHTNESS_CURVE_IDS =
    ::AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS;
/// Wire ids of the named procedural palette roster.
inline constexpr auto &PALETTE_NAME_IDS = ::AUTHORED_PALETTE_NAME_IDS;
/// Authored palette source and parameters.
using PaletteSelection = ::AuthoredPaletteSelection;
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
      Color::palette_selection_fields<GeneratedPaletteParams, false>(
          GENERATED_PALETTE_GATE, MUTATING_PALETTE_GATE,
          GRADIENT_PALETTE_GATE));

  static constexpr auto TOPOLOGY = std::array{
      TopologyField<GeneratedPaletteParams>{
          "palette-mode", &GeneratedPaletteParams::palette_mode,
          PALETTE_MODE_IDS, 0},
      TopologyField<GeneratedPaletteParams>{
          "palette-mapping", &GeneratedPaletteParams::mapping_mode,
          PALETTE_MAPPING_IDS, 2, true},
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

  PreviewEnumWeights<std::size(PALETTE_MAPPING_IDS)> mapping_blend;
  ///< Transient palette-mapping weights.

  /** @brief Overrides or clears the preview's palette-mapping blend.
   * @param blend Option blend, or empty for the authored mapping.
   * @return True if an override was present. */
  bool set_enum_blend(const std::optional<EnumBlend> &blend) {
    return mapping_blend.set(blend);
  }

  /** @brief Samples the mapping weights at a column time.
   * @param value Sampled column parameters.
   * @param from Outgoing frame parameters.
   * @param to Incoming frame parameters.
   * @param fraction Column frame fraction. */
  static void sample_extra(GeneratedPaletteParams &value,
                           const GeneratedPaletteParams &from,
                           const GeneratedPaletteParams &to, float fraction) {
    value.mapping_blend.sample(from.mapping_blend.get(from.mapping_mode),
                               to.mapping_blend.get(to.mapping_mode), fraction);
  }
};
static_assert(field_ids_unique<GeneratedPaletteParams>());
static_assert(field_defaults_in_range<GeneratedPaletteParams>());

/** @brief Canonical controls of color.modulate.v1. */
using ModulationParams = Modulation::Controls;

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
  return {params.mapping_blend.get(params.mapping_mode),
          params.mapping_frequency,
          Color::oscillated_mapping_offset(params.mapping_phase,
                                           params.phase_oscillation_depth,
                                           state.oscillation_phase),
          ctx.palettes[params.palette_mode],
          {ctx.hue_rotation_lut, ctx.hue_rotation_lut != nullptr, nullptr,
           ctx.hue_rotation_gamut},
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
  /** @brief Allows preview changes in palette mapping to sample by column.
   * @param field Schema field index.
   * @param from Outgoing frame parameters.
   * @param to Incoming frame parameters.
   * @return True for an explicitly blended palette mapping. */
  static bool temporal_hold_field(uint16_t field, const Params &from,
                                  const Params &to) {
    return SCHEMA<ColorizeGeneratedPaletteV4>[field].blend &&
           (from.mapping_blend.active || to.mapping_blend.active);
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
  static constexpr auto METRICS = Color::MODULATION_LUT_METRICS;
  static void init(State &state, InstanceId id) { state.key = id.instance; }
  /** @brief Exact instance label of the modulation source.
   * @param state Instance state.
   * @return Source key. */
  static std::string_view modulation_key(const State &state) {
    return state.key;
  }
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
    const auto *source = Modulation::find_endpoint(ctx.modulations, state.key);
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
