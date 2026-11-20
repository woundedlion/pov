/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/** @file EmberTracery.h
 * @brief Gold and jade meshes lead into a wandering lattice of embers.
 */

#include "core/platform/effect_identity.h"
#include "core/render/pullback/composed_tiling.h"

#include <array>
#include <string_view>

/** @brief Authored stages for EmberTracery. */
struct EmberTracerySpec : Pullback::Spec {
  static constexpr bool AUTHORED_TOPOLOGY = true;
  static constexpr auto TEMPORAL = Pullback::TemporalPolicy::phase();
  static constexpr bool HOLD_MODULATION_TOPOLOGY = true;

  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::RotationProvider<B, "rotate">>,
      Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
          Pullback::AuthoredDirectSurfaceParams, B, true, "displace">::Type>,
      Pullback::Stage::SampleSphere<typename Pullback::SourcePolicyFor<
          Pullback::Source::SphericalTilingSourceParams, B>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE, true>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate1", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE, true>>>>;
};

using EmberTraceryParams = Pullback::ParamsFor<EmberTracerySpec>;

/** @brief Gold and jade meshes lead into a wandering lattice of embers.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 */
template <int W, int H>
class EmberTracery : public Pullback::ComposedEffect<W, H, EmberTracery<W, H>,
                                                     EmberTracerySpec> {
public:
  using Params = EmberTraceryParams;

  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::EmberTracery::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "5fbd52fbf5251da79d88bcaea108f187513c095815a485224d859f4ce453a74f";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "e7b1aed8d33362269c066a30123aa64ed12061a9d46b103cf45703f34a6532d1";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 3> PRESET_IDS{
      "golden-mesh",
      "jade-weave",
      "ember-lattice"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  /** @brief Authored instance label of a native resource.
   * @param resource Native resource key.
   * @return Descriptor instance label. */
  static constexpr std::string_view authored_instance(std::string_view resource) {
    if (resource == "source") return "sample";
    if (resource == "color") return "colorize";
    return resource;
  }
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

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
    value.template get<"color">().palette_chroma = 0.7923977f;
    value.template get<"color">().palette_harmony = 5;
    value.template get<"color">().palette_hue = 0.0f;
    value.template get<"color">().palette_lightness = 0.62f;
    value.template get<"color">().palette_lightness_curve = 0;
    value.template get<"color">().palette_lightness_range = 0.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().palette_mode = 0;
    value.template get<"color">().palette_mutation = 0.0f;
    value.template get<"color">().palette_name = 17;
    value.template get<"color">().palette_name_b = 19;
    value.template get<"color">().palette_spread = 0.07f;
    value.template get<"color">().palette_type = 2;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"displace">().basis = 0;
    value.template get<"displace">().direction = 0.0f;
    value.template get<"displace">().scale = 1.0f;
    value.template get<"displace">().speed = 0.0f;
    value.template get<"displace">().strength = 0.0f;
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
    value.template get<"modulate">().modulation.mechanic.amount = 0.0f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate">().modulation.noise.scale = 1.0f;
    value.template get<"modulate">().modulation.noise.speed = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate">().modulation.sine.offset = 0.0f;
    value.template get<"modulate">().modulation.sine.speed = 0.0f;
    value.template get<"modulate1">().modulation.grid.axis_b.x = 1.0f;
    value.template get<"modulate1">().modulation.grid.axis_b.y = 0.0f;
    value.template get<"modulate1">().modulation.grid.axis_b.z = 0.0f;
    value.template get<"modulate1">().modulation.grid.axis.x = 0.0f;
    value.template get<"modulate1">().modulation.grid.axis.y = 1.0f;
    value.template get<"modulate1">().modulation.grid.axis.z = 0.0f;
    value.template get<"modulate1">().modulation.grid.frequency = 1.0f;
    value.template get<"modulate1">().modulation.grid.frequency_b = 1.0f;
    value.template get<"modulate1">().modulation.grid.offset = 0.0f;
    value.template get<"modulate1">().modulation.grid.offset_b = 0.0f;
    value.template get<"modulate1">().modulation.grid.speed = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis_b.x = 1.0f;
    value.template get<"modulate1">().modulation.interference.axis_b.y = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis_b.z = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis.x = 0.0f;
    value.template get<"modulate1">().modulation.interference.axis.y = 1.0f;
    value.template get<"modulate1">().modulation.interference.axis.z = 0.0f;
    value.template get<"modulate1">().modulation.interference.frequency = 1.0f;
    value.template get<"modulate1">().modulation.interference.frequency_b = 1.0f;
    value.template get<"modulate1">().modulation.interference.offset = 0.0f;
    value.template get<"modulate1">().modulation.interference.offset_b = 0.0f;
    value.template get<"modulate1">().modulation.interference.speed = 0.0f;
    value.template get<"modulate1">().modulation.mechanic.amount = 0.0f;
    value.template get<"modulate1">().modulation.enabled = true;
    value.template get<"modulate1">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate1">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate1">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate1">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate1">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate1">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate1">().modulation.noise.scale = 1.0f;
    value.template get<"modulate1">().modulation.noise.speed = 0.0f;
    value.template get<"modulate1">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate1">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate1">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate1">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate1">().modulation.sine.offset = 0.0f;
    value.template get<"modulate1">().modulation.sine.speed = 0.0f;
    value.template get<"rotate">().spin_speed = 0.0f;
    value.template get<"rotate">().wander = 0.0f;
    value.template get<"source">().line_softness = 0.5f;
    value.template get<"source">().line_thickness = 1.0f;
    value.template get<"source">().output = 0;
    value.template get<"source">().solid = 14;
    return value;
  }

  /**
   * @brief The preset at index in PRESET_IDS and how it departs.
   * @param index Preset number; an index past the last yields preset 0.
   * @return The preset's parameters and departure segue.
   */
  HS_FLASH_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.template get<"color">().palette_chroma = 1.0f;
      value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
      value.template get<"color">().palette_mode = 2;
      value.template get<"color">().palette_name = 9;
      value.template get<"color">().phase_oscillation_depth = 1.0f;
      value.template get<"color">().phase_oscillation_speed = 0.0070175435f;
      value.template get<"displace">().scale = 1.2462963f;
      value.template get<"displace">().speed = -0.0001f;
      value.template get<"displace">().strength = 0.18088235f;
      value.template get<"source">().line_softness = 1.398054f;
      value.template get<"source">().line_thickness = 0.0f;
      value.template get<"source">().solid = 52;
    }
    if (index == 2) {
      value.template get<"color">().palette_chroma = 1.0f;
      value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
      value.template get<"color">().palette_mode = 2;
      value.template get<"color">().palette_name = 10;
      value.template get<"color">().phase_oscillation_depth = 1.0f;
      value.template get<"color">().phase_oscillation_speed = 0.0070175435f;
      value.template get<"displace">().scale = 1.2462963f;
      value.template get<"displace">().speed = -0.0001f;
      value.template get<"displace">().strength = 0.18088235f;
      value.template get<"modulate">().modulation.mechanic.amount = 0.43235293f;
      value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::DIRECT;
      value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
      value.template get<"modulate">().modulation.source = Pattern::SourceKind::DIRECTIONAL_SINE;
      value.template get<"modulate">().modulation.noise.scale = 3.25564f;
      value.template get<"modulate">().modulation.noise.speed = 0.00012352942f;
      value.template get<"modulate">().modulation.sine.frequency = 0.3f;
      value.template get<"modulate">().modulation.sine.speed = 0.005882353f;
      value.template get<"modulate1">().modulation.mechanic.amount = -0.7294118f;
      value.template get<"modulate1">().modulation.evaluation = Pattern::Evaluation::DIRECT;
      value.template get<"modulate1">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
      value.template get<"rotate">().wander = 1.0f;
      value.template get<"source">().line_softness = 1.2350903f;
      value.template get<"source">().line_thickness = 0.0f;
      value.template get<"source">().solid = 50;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.

  /** @brief Binds the authored preset bank. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, EmberTracery<W, H>,
                                    EmberTracerySpec>::choreography_hooks()
        .template extend<
            EmberTracery,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &EmberTracery::initial_params>,
            Control::HookBinding<Control::ChoreographyHook::PRESET_ROW,
                                 &EmberTracery::preset>>();
  }
};
