/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Involution.h
 * @brief An affine lattice turns across an area-corrected cubic atlas beneath layered color.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>
#include "core/render/pullback/composed_effect.h"

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
template <typename FX> struct ComposedWhiteBox;
}
#endif

/** @brief Stage spec for `Involution`. */
struct InvolutionSpec : Pullback::Spec {
  /// Samples controls at each column's display time.
  static constexpr bool TEMPORAL_SAMPLING = true;
  /// Resolves stage topology from the authored presets.
  static constexpr bool AUTHORED_TOPOLOGY = true;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  /** @brief Authored sphere, projection, warp, source and color stages.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::RotationProvider<B, "camera">>,
      Pullback::Stage::Project<Pullback::AuthoredProjection<
          Pullback::ProjectionKind::CUBIC,
          Pullback::ProjectionProvider<B, Pullback::AtlasProjectionParams>>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::AffineParams, B, "outer_warp", true,
          Pullback::AffineForm::ROTATES |
              Pullback::AffineForm::TRANSLATES>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<
              Pullback::Source::AuthoredSquareCellsSampleParams, B>::Type,
          Pullback::Weight::Authored<Pullback::SourceProvider<
              B, Pullback::Source::AuthoredSquareCellsSampleParams>>,
          Pullback::ProjectionCoverage::Authored<Pullback::SourceProvider<
              B, Pullback::Source::AuthoredSquareCellsSampleParams>>>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE, true>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE, true>>>>;
};
/// Parameter storage derived from `InvolutionSpec`.
using InvolutionParams = Pullback::ParamsFor<InvolutionSpec>;

/** @brief An affine lattice turns across an area-corrected cubic atlas beneath layered color.
 * @tparam W Canvas width.
 * @tparam H Canvas height. */
template <int W, int H>
class Involution
    : public Pullback::ComposedEffect<W, H, Involution<W, H>, InvolutionSpec> {
#if HS_ENABLE_TEST_HOOKS
  friend struct ::hs_test::ComposedWhiteBox<Involution<W, H>>;
#endif
public:
  using Params = InvolutionParams; ///< Live parameter record.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::Involution::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "dda645c7d860a934f341f5a4dd8d2232092b934832331aa4a7fcc2a6193cd8a6";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "add6a794ad2c364eb058cbb38669b8988d1e1a4f52935ce833750a6b0e378573";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "involution"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  /** @brief Authored instance label of a native resource.
   * @param resource Native resource key.
   * @return Descriptor instance label. */
  static constexpr std::string_view authored_instance(std::string_view resource) {
    if (resource == "projection") return "project";
    if (resource == "outer_warp") return "warp";
    if (resource == "source") return "sample";
    if (resource == "color") return "colorize";
    return resource;
  }
  // clang-format on
  // End generated identity.
  /// Snapshot parameter layout version.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 5;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"camera">().wander = 1.0f;
    value.template get<"color">().brightness_envelope = Pullback::Color::BrightnessEnvelope::NONE;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.16041666f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_hue = 0.30701753f;
    value.template get<"color">().palette_lightness = 0.36549708f;
    value.template get<"color">().palette_lightness_curve = 4;
    value.template get<"color">().palette_lightness_range = 0.6754386f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().palette_mode = 2;
    value.template get<"color">().palette_spread = 0.059210528f;
    value.template get<"color">().palette_type = 1;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"modulate">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"modulate">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"modulate">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"modulate">().modulation.grid.axis.x = 0.0f;
    value.template get<"modulate">().modulation.grid.axis.y = 1.0f;
    value.template get<"modulate">().modulation.grid.axis.z = 0.0f;
    value.template get<"modulate">().modulation.grid.frequency = 1.0f;
    value.template get<"modulate">().modulation.grid.frequency_b = 1.0f;
    value.template get<"modulate">().modulation.grid.offset = 0.0f;
    value.template get<"modulate">().modulation.grid.offset_b = 0.0f;
    value.template get<"modulate">().modulation.grid.speed = 0.0f;
    value.template get<"modulate">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"modulate">().modulation.interference.axis_b.y = 0.0f;
    value.template get<"modulate">().modulation.interference.axis_b.z = 0.0f;
    value.template get<"modulate">().modulation.interference.axis.x = 0.0f;
    value.template get<"modulate">().modulation.interference.axis.y = 1.0f;
    value.template get<"modulate">().modulation.interference.axis.z = 0.0f;
    value.template get<"modulate">().modulation.interference.frequency = 1.0f;
    value.template get<"modulate">().modulation.interference.frequency_b = 1.0f;
    value.template get<"modulate">().modulation.interference.offset = 0.0f;
    value.template get<"modulate">().modulation.interference.offset_b = 0.0f;
    value.template get<"modulate">().modulation.interference.speed = 0.0f;
    value.template get<"modulate">().modulation.mechanic.amount = 0.4970588f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate">().modulation.noise.scale = 1.3329504f;
    value.template get<"modulate">().modulation.noise.speed = 0.00011176471f;
    value.template get<"modulate">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate">().modulation.sine.frequency = 13.741177f;
    value.template get<"modulate">().modulation.sine.offset = 0.0f;
    value.template get<"modulate">().modulation.sine.speed = 0.0f;
    value.template get<"modulation">().modulation.mechanic.amount = 0.83882546f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
    value.template get<"modulation">().modulation.noise.scale = 0.14453939f;
    value.template get<"modulation">().modulation.noise.speed = -0.0000010416667f;
    value.template get<"projection">().area_correction = 0.7f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().coordinate_scale = 1.3574561f;
    value.template get<"projection">().frame = 1;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"source">().cell_radius = 0.2777796f;
    value.template get<"source">().cell_scale = 3.9611204f;
    value.template get<"source">().cell_shape = 1.0f;
    value.template get<"source">().cell_softness = 0.5f;
    value.template get<"source">().coverage_mode = 3;
    value.template get<"source">().edge_width = 0.14035088f;
    value.template get<"source">().weight_mode = 1;
    value.template get<"outer_warp">().lattice_period = 1.0f;
    value.template get<"outer_warp">().repeat_cells = 8.0f;
    value.template get<"outer_warp">().rotation_rate = 1.0005351f;
    value.template get<"outer_warp">().scale_x = 1.0f;
    value.template get<"outer_warp">().scale_y = 1.0f;
    value.template get<"outer_warp">().shear = 0.0f;
    value.template get<"outer_warp">().speed = 0.0065036383f;
    value.template get<"outer_warp">().translation_x = 1.0f;
    value.template get<"outer_warp">().translation_y = 1.0f;
    return value;
  }
  // clang-format on
  // End generated params.

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, Involution<W, H>,
                                    InvolutionSpec>::choreography_hooks()
        .template extend<
            Involution,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &Involution::initial_params>>();
  }
};
