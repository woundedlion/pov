/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Iridescence.h
 * @brief Composed pullback effect:
 *        a dodecahedral potential flow folding a grid under procedural color.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `Iridescence`. */
struct IridescenceSpec : Pullback::Spec {
  /// Samples animation at each column time.
  static constexpr bool TEMPORAL_SAMPLING = true;
  /// Generated palette harmony.
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  /// Identity projection frame.
  static constexpr bool ANIMATED_PROJECTION = false;
  /** @brief Potential flow, equirectangular grid, procedural palette, baked
   *         noise lightness and baked noise hue.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
          Pullback::PotentialParams, B, false, "displace">::Type>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::EQUIRECTANGULAR, B>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::GridSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              Pullback::ProjectionCoverageMode::WEIGHT, B, void>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate2", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>>;
};
/// Parameter block derived from `IridescenceSpec`.
using IridescenceParams = Pullback::ParamsFor<IridescenceSpec>;

/**
 * @brief A dodecahedral potential flow folding a grid under procedural color.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class Iridescence : public Pullback::ComposedEffect<W, H, Iridescence<W, H>,
                                                    IridescenceSpec> {

public:
  using Params = IridescenceParams; ///< Live parameter struct.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::Iridescence::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "cf0ea806a8eaf497509191e74cb73f5faee04d8f3411b99c2a042d91b50d05b6";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "9775be782118aa161753da4555e78adf1f9c8d9a47b5b63046a7f9b62a68c99c";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "iridescence"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /// Params layout version; snapshots tagged with another are rejected.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"color">().brightness_bottom = 0.0f;
    value.template get<"color">().brightness_envelope = Pullback::Color::BrightnessEnvelope::NONE;
    value.template get<"color">().brightness_top = 1.0f;
    value.template get<"color">().gradient_end_b = 1.0f;
    value.template get<"color">().gradient_end_g = 0.2f;
    value.template get<"color">().gradient_end_r = 0.1f;
    value.template get<"color">().gradient_middle_b = 0.4f;
    value.template get<"color">().gradient_middle_g = 1.0f;
    value.template get<"color">().gradient_middle_r = 0.1f;
    value.template get<"color">().gradient_start_b = 0.0f;
    value.template get<"color">().gradient_start_g = 0.15f;
    value.template get<"color">().gradient_start_r = 1.0f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_harmony = 5;
    value.template get<"color">().palette_hue = 0.0f;
    value.template get<"color">().palette_lightness = 0.62f;
    value.template get<"color">().palette_lightness_curve = 0;
    value.template get<"color">().palette_lightness_range = 0.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().palette_mode = 2;
    value.template get<"color">().palette_mutation = 0.0f;
    value.template get<"color">().palette_name = 10;
    value.template get<"color">().palette_name_b = 19;
    value.template get<"color">().palette_spread = 0.07f;
    value.template get<"color">().palette_type = 2;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"displace">().drift = 0.03f;
    value.template get<"displace">().integrator = 0;
    value.template get<"displace">().orientation_azimuth = 0.0f;
    value.template get<"displace">().orientation_polar = 0.98867774f;
    value.template get<"displace">().strength = 1.0163976f;
    value.template get<"displace">().terms[0].amplitude = 1.0117648f;
    value.template get<"displace">().terms[0].kind = 1;
    value.template get<"displace">().terms[0].l = 1;
    value.template get<"displace">().terms[0].m = 1;
    value.template get<"displace">().terms[0].phase = 0.9424778f;
    value.template get<"displace">().terms[0].type = 3;
    value.template get<"displace">().terms[1].amplitude = 0.0f;
    value.template get<"displace">().terms[1].kind = 0;
    value.template get<"displace">().terms[1].l = 1;
    value.template get<"displace">().terms[1].m = 0;
    value.template get<"displace">().terms[1].phase = 0.0f;
    value.template get<"displace">().terms[1].type = 0;
    value.template get<"displace">().terms[2].amplitude = 0.0f;
    value.template get<"displace">().terms[2].kind = 0;
    value.template get<"displace">().terms[2].l = 1;
    value.template get<"displace">().terms[2].m = 0;
    value.template get<"displace">().terms[2].phase = 0.0f;
    value.template get<"displace">().terms[2].type = 0;
    value.template get<"displace">().terms[3].amplitude = 0.0f;
    value.template get<"displace">().terms[3].kind = 0;
    value.template get<"displace">().terms[3].l = 1;
    value.template get<"displace">().terms[3].m = 0;
    value.template get<"displace">().terms[3].phase = 0.0f;
    value.template get<"displace">().terms[3].type = 0;
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
    value.template get<"modulate">().modulation.mechanic.amount = 0.8f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate">().modulation.noise.scale = 2.0f;
    value.template get<"modulate">().modulation.noise.speed = 0.00037058824f;
    value.template get<"modulate">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate">().modulation.sine.offset = 0.0f;
    value.template get<"modulate">().modulation.sine.speed = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"modulate2">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis.x = 0.0f;
    value.template get<"modulate2">().modulation.grid.axis.y = 1.0f;
    value.template get<"modulate2">().modulation.grid.axis.z = 0.0f;
    value.template get<"modulate2">().modulation.grid.frequency = 1.0f;
    value.template get<"modulate2">().modulation.grid.frequency_b = 1.0f;
    value.template get<"modulate2">().modulation.grid.offset = 0.0f;
    value.template get<"modulate2">().modulation.grid.offset_b = 0.0f;
    value.template get<"modulate2">().modulation.grid.speed = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"modulate2">().modulation.interference.axis_b.y = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis_b.z = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis.x = 0.0f;
    value.template get<"modulate2">().modulation.interference.axis.y = 1.0f;
    value.template get<"modulate2">().modulation.interference.axis.z = 0.0f;
    value.template get<"modulate2">().modulation.interference.frequency = 1.0f;
    value.template get<"modulate2">().modulation.interference.frequency_b = 1.0f;
    value.template get<"modulate2">().modulation.interference.offset = 0.0f;
    value.template get<"modulate2">().modulation.interference.offset_b = 0.0f;
    value.template get<"modulate2">().modulation.interference.speed = 0.0f;
    value.template get<"modulate2">().modulation.mechanic.amount = 0.7764706f;
    value.template get<"modulate2">().modulation.enabled = true;
    value.template get<"modulate2">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate2">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate2">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate2">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate2">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate2">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate2">().modulation.noise.scale = 1.0f;
    value.template get<"modulate2">().modulation.noise.speed = 0.0f;
    value.template get<"modulate2">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate2">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate2">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate2">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate2">().modulation.sine.offset = 0.0f;
    value.template get<"modulate2">().modulation.sine.speed = 0.0f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().singularity_fade = 20.0f;
    value.template get<"source">().complexity = 0.0f;
    value.template get<"source">().pattern_freq = 5.436053f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.0f;
    return value;
  }
  // clang-format on
  // End generated params.

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, Iridescence<W, H>,
                                    IridescenceSpec>::choreography_hooks()
        .template extend<
            Iridescence,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &Iridescence::initial_params>>();
  }
};
