/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file AshCloud.h
 * @brief Composed pullback effect:
 *        a softly cut lattice curled across dodecahedral facets.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Pullback stage spec for `AshCloud`. */
struct AshCloudSpec : Pullback::Spec {
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr bool TEMPORAL_COLOR_PHASE = true;
  static constexpr bool PALETTE_LAB_CACHE = true;
  static constexpr bool CACHE_FULL_GAMUT_LUT = true;
  static constexpr bool RIPPLE_RADIAL_TABLE = true;
  static constexpr bool PERSISTENT_PREPARED_FRAME = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Stage chain: camera rotation, burst ripples, surface curl, dodecahedral
   *  kaleidoscope lens, projection, lattice source, value cutout, generated
   *  palette,
   *  noise and displacement hue modulation.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Displace<
          Pullback::Surface::BurstRipple<Pullback::SurfaceProvider<
              B, Pullback::SphereAnimation::BurstRippleParams, true,
              "displace">>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
              Pullback::SurfaceNoiseParams, B, true>::Type>,
          Pullback::Stage::Lens<Pullback::Lens::DodecahedralKaleidoscope>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::STEREOGRAPHIC, B>::Type>>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::LatticeSourceParams, B>::Type,
                              Pullback::Weight::Projection,
                              typename Pullback::CoveragePolicyFor<
                                  Pullback::ProjectionCoverageMode::WEIGHT, B,
                                  Pullback::CutoutValueParams>::Type>,
      Pullback::Stage::ApplyCoverage<Pullback::ValueCoverage::ValueCutout<
          Pullback::ValueProvider<B, Pullback::CutoutValueParams>>>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, Pullback::Color::BrightnessEnvelope::NONE,
                                  PALETTE_LAB_CACHE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::DIRECT, TEMPORAL_COLOR_PHASE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::DISPLACEMENT,
              Pattern::Evaluation::DIRECT, TEMPORAL_COLOR_PHASE>>>>;
};
/// Parameter block derived from `AshCloudSpec`.
using AshCloudParams = Pullback::ParamsFor<AshCloudSpec>;

/**
 * @brief A softly cut lattice curled across dodecahedral facets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AshCloud
    : public Pullback::ComposedEffect<W, H, AshCloud<W, H>, AshCloudSpec> {

  using Base = Pullback::ComposedEffect<W, H, AshCloud<W, H>, AshCloudSpec>;

public:
  /// Parameter block derived from `AshCloudSpec`.
  using Params = AshCloudParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::AshCloud::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "40605c23a59db06f66579a3b652c5ba45ab04e8f54c1e43f05e7453cf9bb80ff";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "38bddb9ffc08b44f637e3af83ecb9daaeafa83958523e4f94c8abc442c291e62";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "ash-cloud"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  /// Outer camera spin about +Y, in radians per frame.
  static constexpr float CAMERA_SPIN_RATE = 0.01975f;
  // clang-format on
  // End generated identity.
  /// Snapshot schema version; changes with the `Params` layout.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 5;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().mapping_frequency = 3.666f;
    value.template get<"color">().mapping_phase = 0.598f;
    value.template get<"color">().palette_chroma = 0.314f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.856f;
    value.template get<"color">().phase_oscillation_speed = 0.00448f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"value">().cutout_softness = 0.31635937f;
    value.template get<"value">().cutout_threshold = 0.365f;
    value.template get<"displace">().burst_rate = 0.015625f;
    value.template get<"displace">().center_azimuth = 0.0f;
    value.template get<"displace">().center_polar = 1.5707964f;
    value.template get<"displace">().decay = 0.1f;
    value.template get<"displace">().period = 80.0f;
    value.template get<"displace">().strength = 0.15f;
    value.template get<"displace">().thickness = 0.7f;
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
    value.template get<"modulate">().modulation.mechanic.amount = -1.3882353f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
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
    value.template get<"modulation">().modulation.mechanic.amount = -3.0823529f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulation">().modulation.noise.scale = 2.4029531f;
    value.template get<"modulation">().modulation.noise.speed = 0.00015f;
    value.template get<"projection">().spin_rate = 0.00175f;
    value.template get<"projection">().wander = 0.009f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().lattice_cell_scale = 0.98971874f;
    value.template get<"source">().lattice_radius = 0.37757313f;
    value.template get<"source">().lattice_shape_blend = 0.0f;
    value.template get<"source">().lattice_softness = 0.34563965f;
    value.template get<"surface">().scale = 5.7742186f;
    value.template get<"surface">().speed = 0.00040625f;
    value.template get<"surface">().strength = 0.043f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Evaluates the composed pixel pipeline from hot FLASH.
   * @param view Unit view direction.
   * @param frame Prepared column state.
   * @return Straight-alpha pixel color. */
  HS_HOT_FLASH_MEMBER static Color4 shade(const math::Vector &view,
                                          const typename Base::Frame &frame) {
    return Base::RenderPipeline::evaluate(view, frame.ctx, frame.prepared);
  }

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, AshCloud<W, H>,
                                    AshCloudSpec>::choreography_hooks()
        .template extend<
            AshCloud,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &AshCloud::initial_params>>();
  }
};
