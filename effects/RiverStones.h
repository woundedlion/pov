/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file RiverStones.h
 * @brief Ripple bursts across a hexagonal lattice of soft cells.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `RiverStones`. */
struct RiverStonesSpec : Pullback::Spec {
  static constexpr bool SOURCE_MORPH = true;
#if defined(ARDUINO_TEENSY40) || defined(ARDUINO_TEENSY41) ||                  \
    (HS_ENABLE_TEST_HOOKS && defined(HS_TEST_RIVER_STONES_HEX_TABLE))
  static constexpr bool HEX_MORPH_TABLE = true;
#endif
  static constexpr bool RIPPLE_RADIAL_TABLE = true;
  /// Samples animation at each column time.
  static constexpr auto TEMPORAL =
      Pullback::TemporalPolicy::fused(Pullback::PaletteCache::NONE);
  /// Generated palette harmony.
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  /// Spin and wander projection frame.
  static constexpr auto PROJECTION_MOTION =
      Pullback::ProjectionMotion::ANIMATED;

  /** @brief Camera, ripple, translating hex cells, palette and noise layers.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
          Pullback::SphereAnimation::BurstRippleParams, B, true,
          "displace1">::Type>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::EQUIRECTANGULAR, B>::Type>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::Warp::AuthoredAffineParams, B, "outer_warp", true,
          Pullback::AffineForm::TRANSLATES>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::FixedTessellationParams<
                                      Pattern::TessellationKind::HEXAGONAL,
                                      Pattern::TessellationFill::CELLS>,
                                  B>::Type,
                              Pullback::Weight::Projection,
                              Pullback::ProjectionCoverage::Weight>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::DirectOrDisplacementProvider<
              B, "modulate", Pattern::SourceKind::NOISE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate1", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate2", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE>>>>;
};
/// Parameter block derived from `RiverStonesSpec`.
using RiverStonesParams = Pullback::ParamsFor<RiverStonesSpec>;

/**
 * @brief Ripple bursts across a hexagonal lattice of soft cells.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class RiverStones : public Pullback::ComposedEffect<W, H, RiverStones<W, H>,
                                                    RiverStonesSpec> {
public:
  using Params = RiverStonesParams; ///< Live parameter struct.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::RiverStones::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "16da92705d798daf9c6602503e59214136e34148f5dae4e2f10a16cde25e482e";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "ff2fc61d3325d672e71958331b7d653041037841e378d3135e3508d5c205a31c";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "river-stones",
      "river-stones-2"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  /// Outer camera spin about +Y, in radians per frame.
  static constexpr float CAMERA_SPIN_RATE = 0.0f;
  // clang-format on
  // End generated identity.
  /// Params layout version; snapshots tagged with another are rejected.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 3;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 0.1f;
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
    value.template get<"displace1">().burst_rate = 0.015625f;
    value.template get<"displace1">().burst_size = 3;
    value.template get<"displace1">().center_azimuth = 0.0f;
    value.template get<"displace1">().center_polar = 1.5707964f;
    value.template get<"displace1">().decay = 0.0f;
    value.template get<"displace1">().period = 97.46764f;
    value.template get<"displace1">().random_center = 1;
    value.template get<"displace1">().strength = 0.15f;
    value.template get<"displace1">().thickness = 0.7f;
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
    value.template get<"modulate">().modulation.mechanic.amount = 0.25882354f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::DIRECT;
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
    value.template get<"modulate">().modulation.sine.speed = 0.052941177f;
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
    value.template get<"modulate1">().modulation.mechanic.amount = 0.4617647f;
    value.template get<"modulate1">().modulation.enabled = true;
    value.template get<"modulate1">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate1">().modulation.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
    value.template get<"modulate1">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate1">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate1">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate1">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate1">().modulation.noise.scale = 1.8775138f;
    value.template get<"modulate1">().modulation.noise.speed = 0.0001f;
    value.template get<"modulate1">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate1">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate1">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate1">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate1">().modulation.sine.offset = 0.0f;
    value.template get<"modulate1">().modulation.sine.speed = 0.0f;
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
    value.template get<"modulate2">().modulation.mechanic.amount = 0.6382353f;
    value.template get<"modulate2">().modulation.enabled = true;
    value.template get<"modulate2">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate2">().modulation.mechanic.kind = ColorModulation::Mechanic::OPACITY;
    value.template get<"modulate2">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate2">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate2">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate2">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate2">().modulation.noise.scale = 1.5627832f;
    value.template get<"modulate2">().modulation.noise.speed = 0.0001f;
    value.template get<"modulate2">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate2">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate2">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate2">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate2">().modulation.sine.offset = 0.0f;
    value.template get<"modulate2">().modulation.sine.speed = 0.0f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().cell_radius = 0.41925988f;
    value.template get<"source">().cell_scale = 2.5352101f;
    value.template get<"source">().cell_shape = 0.5292398f;
    value.template get<"source">().cell_softness = 0.10124387f;
    value.template get<"outer_warp">().lattice_period = 1.0f;
    value.template get<"outer_warp">().path_length = 1;
    value.template get<"outer_warp">().repeat_cells = 1.0f;
    value.template get<"outer_warp">().rotation_rate = 0.0f;
    value.template get<"outer_warp">().scale_x = 1.0f;
    value.template get<"outer_warp">().scale_y = 1.0f;
    value.template get<"outer_warp">().shear = 0.0f;
    value.template get<"outer_warp">().speed = 0.02f;
    value.template get<"outer_warp">().translation_x = 0.0f;
    value.template get<"outer_warp">().translation_y = 2.0f;
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
      value.template get<"modulate">().modulation.mechanic.amount = 2.0f;
      value.template get<"modulate">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
      value.template get<"source">().cell_scale = 4.774648f;
      value.template get<"source">().cell_softness = 0.32801026f;
      value.template get<"outer_warp">().path_length = 0;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, RiverStones<W, H>,
                                    RiverStonesSpec>::choreography_hooks()
        .template extend<
            RiverStones,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &RiverStones::initial_params>,
            Control::HookBinding<Control::ChoreographyHook::PRESET_ROW,
                                 &RiverStones::preset>>();
  }
};
