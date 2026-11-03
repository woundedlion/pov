/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file WarpFest.h
 * @brief Ripple bursts and symmetric displacement fold a wave-sheared grid.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>
#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `WarpFest`. */
struct WarpFestSpec : Pullback::Spec {
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr bool AUTHORED_TOPOLOGY = true;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  /** @brief Authored sphere, projection, warp, source and color stages.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Rotate<
          Pullback::OuterCameraProvider<B, Pullback::AtlasProjectionParams>>,
      Pullback::Stage::Displace<
          Pullback::Surface::BurstRipple<Pullback::SurfaceProvider<
              B, Pullback::SphereAnimation::BurstRippleParams, true,
              "displace2">>>,
      Pullback::Stage::Displace<
          Pullback::Surface::Symmetric<Pullback::SurfaceProvider<
              B, Pullback::Surface::SymmetricParams, true, "displace">>>,
      Pullback::Stage::Lens<
          Pullback::Lens::ConfiguredKaleidoscope<Pullback::KeyedLensProvider<
              B, Pullback::Lens::KaleidoscopeParams, "lens">>>,
      Pullback::Stage::Project<Pullback::AuthoredProjection<
          Pullback::ProjectionKind::GNOMONIC_FOLDED,
          Pullback::ProjectionProvider<B, Pullback::AtlasProjectionParams>>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::AuthoredWaveShearParams, B, "outer_warp", true>::Type>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::MirrorParams, B, "inner_warp", true>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<
              Pullback::Source::AuthoredGridSampleParams, B>::Type,
          Pullback::Weight::Authored<Pullback::SourceProvider<
              B, Pullback::Source::AuthoredGridSampleParams>>,
          Pullback::ProjectionCoverage::Authored<Pullback::SourceProvider<
              B, Pullback::Source::AuthoredGridSampleParams>>>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE, false, true>>>>;
};
/// Parameter storage derived from `WarpFestSpec`.
using WarpFestParams = Pullback::ParamsFor<WarpFestSpec>;

/** @brief Ripple bursts and symmetric displacement fold a wave-sheared grid.
 * @tparam W Canvas width.
 * @tparam H Canvas height. */
template <int W, int H>
class WarpFest
    : public Pullback::ComposedEffect<W, H, WarpFest<W, H>, WarpFestSpec> {
public:
  using Params = WarpFestParams; ///< Live parameter record.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::WarpFest::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "ffa0c5273d090ee9aa9e1b97266245005c2756871aa24caebc8ff007254dc118";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "1e101518acd123ddb24c4e6490eba9d292ed92929253827ea0b30f571878381e";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "wave-mirror",
      "cup-hue"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  /// Snapshot parameter layout version.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  /**
   * @brief Parameters of the first preset; every preset varies from them.
   * @return The preset-0 `Params`.
   */
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().brightness_envelope = Pullback::Color::BrightnessEnvelope::NONE;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().palette_mode = 0;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"displace">().integrator = 0;
    value.template get<"displace">().kind = 0;
    value.template get<"displace">().speed = 0.0f;
    value.template get<"displace">().strength = 0.05f;
    value.template get<"displace">().symmetry = 0;
    value.template get<"displace2">().burst_rate = 0.0078125f;
    value.template get<"displace2">().burst_size = 3;
    value.template get<"displace2">().center_azimuth = 0.0f;
    value.template get<"displace2">().center_polar = 1.5707964f;
    value.template get<"displace2">().decay = 0.1f;
    value.template get<"displace2">().period = 80.0f;
    value.template get<"displace2">().random_center = 1;
    value.template get<"displace2">().strength = 0.15f;
    value.template get<"displace2">().thickness = 0.7f;
    value.template get<"lens">().symmetry = 3;
    value.template get<"modulation">().modulation.mechanic.amount = 0.721f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulation">().modulation.noise.scale = 1.0f;
    value.template get<"modulation">().modulation.noise.speed = 0.0f;
    value.template get<"projection">().frame = 0;
    value.template get<"projection">().hemisphere = 0;
    value.template get<"projection">().singularity_fade = 2.311f;
    value.template get<"source">().angle_rate = 0.027f;
    value.template get<"source">().complexity = 1.704f;
    value.template get<"source">().coverage_mode = 2;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 6.3287f;
    value.template get<"source">().pattern_mix = 0.0f;
    value.template get<"source">().speed = 0.04f;
    value.template get<"source">().weight_mode = 1;
    value.template get<"outer_warp">().field_angle = 2.2305307f;
    value.template get<"outer_warp">().frequency = 1.408f;
    value.template get<"outer_warp">().speed = -0.00325f;
    value.template get<"outer_warp">().strength = -0.176f;
    value.template get<"inner_warp">().cell_x = 1.0f;
    value.template get<"inner_warp">().cell_y = 1.0f;
    value.template get<"inner_warp">().offset_x = 0.0f;
    value.template get<"inner_warp">().offset_y = 0.0f;
    value.template get<"inner_warp">().rotation = 0.0f;
    value.template get<"inner_warp">().speed = 0.0f;
    return value;
  }

  /**
   * @brief The preset at index in PRESET_IDS and how it departs.
   * @param index Preset number; an index past the last yields preset 0.
   * @return The preset's parameters and departure segue.
   */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
      value.template get<"color">().palette_mode = 2;
      value.template get<"color">().opacity_high = 0.23684211f;
      value.template get<"color">().opacity_low = 0.9502924f;
      value.template get<"displace">().speed = 0.026470589f;
      value.template get<"displace">().strength = 0.24393544f;
      value.template get<"displace">().symmetry = 3;
      value.template get<"displace2">().burst_rate = 0.015625f;
      value.template get<"modulation">().modulation.mechanic.amount = 1.3411765f;
      value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::DIRECT;
      value.template get<"modulation">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
      value.template get<"modulation">().modulation.noise.scale = 1.9717969f;
      value.template get<"source">().secondary_rate = 0.0f;
      value.template get<"source">().pattern_freq = 9.9265785f;
      value.template get<"source">().speed = 0.07f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, WarpFest<W, H>,
                                    WarpFestSpec>::choreography_hooks()
        .template extend<
            WarpFest,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &WarpFest::initial_params>>();
  }
};
