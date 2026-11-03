/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file GlitchLace.h
 * @brief Twisted tangent fields lace a moving lattice across an S2 cube atlas.
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

/** @brief Stage spec for `GlitchLace`. */
struct GlitchLaceSpec : Pullback::Spec {
  static constexpr bool TEMPORAL_SAMPLING = true;
  static constexpr bool AUTHORED_TOPOLOGY = true;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  /** @brief Authored sphere, projection, warp, source and color stages.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Lens<
          Pullback::Lens::AnimatedTwist<Pullback::TwistProvider<B, "lens1">>,
          false>,
      Pullback::Stage::Rotate<Pullback::RotationProvider<B, "camera">>,
      Pullback::Stage::Displace<
          Pullback::Surface::RossbyHaurwitz<Pullback::SurfaceProvider<
              B, Pullback::Surface::RossbyHaurwitzParams, true, "displace">>>,
      Pullback::Stage::Displace<
          Pullback::Surface::Symmetric<Pullback::SurfaceProvider<
              B, Pullback::Surface::SymmetricParams, true, "displace1">>>,
      Pullback::Stage::Project<Pullback::AuthoredProjection<
          Pullback::ProjectionKind::S2,
          Pullback::ProjectionProvider<B, Pullback::AtlasProjectionParams>>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::ExplicitAffineParams, B, "outer_warp", true>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<
              Pullback::Source::AuthoredLatticeSampleParams, B>::Type,
          Pullback::Weight::Authored<Pullback::SourceProvider<
              B, Pullback::Source::AuthoredLatticeSampleParams>>,
          Pullback::ProjectionCoverage::Authored<Pullback::SourceProvider<
              B, Pullback::Source::AuthoredLatticeSampleParams>>>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulate", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE, false, true>>>>;
};
/// Parameter storage derived from `GlitchLaceSpec`.
using GlitchLaceParams = Pullback::ParamsFor<GlitchLaceSpec>;

/** @brief Twisted tangent fields lace a moving lattice across an S2 cube atlas.
 * @tparam W Canvas width.
 * @tparam H Canvas height. */
template <int W, int H>
class GlitchLace
    : public Pullback::ComposedEffect<W, H, GlitchLace<W, H>, GlitchLaceSpec> {
#if HS_ENABLE_TEST_HOOKS
  friend struct ::hs_test::ComposedWhiteBox<GlitchLace<W, H>>;
#endif
public:
  using Params = GlitchLaceParams; ///< Live parameter record.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::GlitchLace::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "127d39c31d03c927f811fa19a657afe2a004fde0aac3aa4f62442f93543a90ac";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "5255f1fd0d9203844567e350e794a7c9a6e0aef71169def9b1bdc0134b5d84ee";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "grid-worms"
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
    value.template get<"camera">().wander = 1.0f;
    value.template get<"color">().brightness_envelope = Pullback::Color::BrightnessEnvelope::NONE;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().palette_mode = 2;
    value.template get<"color">().palette_type = 0;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 0.9661458f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"displace">().drift = 0.081764705f;
    value.template get<"displace">().integrator = 0;
    value.template get<"displace">().omega = 1.3882353f;
    value.template get<"displace">().order = 3;
    value.template get<"displace">().strength = 0.09239978f;
    value.template get<"displace">().wave = 4.0f;
    value.template get<"displace1">().integrator = 0;
    value.template get<"displace1">().kind = 1;
    value.template get<"displace1">().strength = -1.0657501f;
    value.template get<"displace1">().symmetry = 3;
    value.template get<"lens1">().drift = 0.01632353f;
    value.template get<"lens1">().speed = 0.027205883f;
    value.template get<"lens1">().twist_angle = -10.259509f;
    value.template get<"lens1">().distribution = 0.999f;
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
    value.template get<"modulate">().modulation.mechanic.amount = 0.5882353f;
    value.template get<"modulate">().modulation.enabled = true;
    value.template get<"modulate">().modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
    value.template get<"modulate">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulate">().modulation.source = Pattern::SourceKind::NOISE;
    value.template get<"modulate">().modulation.noise.contrast = 0.0f;
    value.template get<"modulate">().modulation.noise.basis = Pattern::NoiseBasis::SIMPLEX;
    value.template get<"modulate">().modulation.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
    value.template get<"modulate">().modulation.noise.scale = 3.3737593f;
    value.template get<"modulate">().modulation.noise.speed = 0.00015882353f;
    value.template get<"modulate">().modulation.sine.axis.x = 0.0f;
    value.template get<"modulate">().modulation.sine.axis.y = 1.0f;
    value.template get<"modulate">().modulation.sine.axis.z = 0.0f;
    value.template get<"modulate">().modulation.sine.frequency = 1.0f;
    value.template get<"modulate">().modulation.sine.offset = 0.0f;
    value.template get<"modulate">().modulation.sine.speed = 0.0f;
    value.template get<"projection">().central_meridian = 0.0f;
    value.template get<"projection">().coordinate_scale = 1.0f;
    value.template get<"projection">().frame = 1;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"source">().coverage_mode = 3;
    value.template get<"source">().edge_width = 0.32748538f;
    value.template get<"source">().lattice_cell_scale = 2.420276f;
    value.template get<"source">().lattice_radius = 0.48583883f;
    value.template get<"source">().lattice_shape_blend = 1.0f;
    value.template get<"source">().lattice_softness = 0.22006065f;
    value.template get<"source">().weight_mode = 1;
    value.template get<"outer_warp">().lattice_period = 1.0f;
    value.template get<"outer_warp">().rotation_rate = 4.9895883f;
    value.template get<"outer_warp">().scale_x = 0.9565717f;
    value.template get<"outer_warp">().scale_y = 1.0f;
    value.template get<"outer_warp">().shear = 0.0f;
    value.template get<"outer_warp">().speed = 0.0036470587f;
    value.template get<"outer_warp">().translation_x = 4.0f;
    value.template get<"outer_warp">().translation_y = 0.25882354f;
    return value;
  }
  // clang-format on
  // End generated params.

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, GlitchLace<W, H>,
                                    GlitchLaceSpec>::choreography_hooks()
        .template extend<
            GlitchLace,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &GlitchLace::initial_params>>();
  }
};
