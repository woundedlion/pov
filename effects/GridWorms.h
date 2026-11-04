/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file GridWorms.h
 * @brief Composed pullback effect:
 *        an animated twist pulling a folded lattice into flowing worms.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `GridWorms`. */
struct GridWormsSpec : Pullback::Spec {
  /// Samples animation at each column time.
  static constexpr bool TEMPORAL_SAMPLING = true;
  /// Generated palette harmony.
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::ANALOGOUS;
  /// Projection-frame spin and wander.
  static constexpr bool ANIMATED_PROJECTION = true;
  /// Lens precedes camera rotation.
  static constexpr bool FACTORED_COLUMNS = false;
  /** @brief Twist, camera, folded projection, lattice, palette and sine hue.
   * @tparam B Effect binding. */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Lens<
          Pullback::Lens::AnimatedTwist<Pullback::TwistProvider<B, "lens1">>,
          false>,
      Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::GNOMONIC_FOLDED, B>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::SquareCellsParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              Pullback::ProjectionCoverageMode::WEIGHT, B, void>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::DIRECTIONAL_SINE,
              Pattern::Evaluation::DIRECT>>>>;
};
/// Parameter block derived from `GridWormsSpec`.
using GridWormsParams = Pullback::ParamsFor<GridWormsSpec>;

/**
 * @brief An animated axial twist pulling a folded lattice into flowing worms.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class GridWorms
    : public Pullback::ComposedEffect<W, H, GridWorms<W, H>, GridWormsSpec> {

public:
  using Params = GridWormsParams; ///< Live parameter struct.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::GridWorms::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "dab6a6a162b48430c41839d49a4c46d4efa053f0f67a1db427c4d6b9882264f1";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "8a1bb10bf172baf8d666885f61bbaa8e251cc78202976ca01c3aa8bdc3495273";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "grid-worms"
  };
  /// Frames each preset holds before advancing.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
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
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 1.0f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 0.9661458f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"lens1">().axis_lambda = 0.0f;
    value.template get<"lens1">().axis_phi = 0.0f;
    value.template get<"lens1">().drift = 0.01632353f;
    value.template get<"lens1">().speed = 0.027205883f;
    value.template get<"lens1">().twist_angle = -10.259509f;
    value.template get<"lens1">().distribution = 0.999f;
    value.template get<"modulation">().modulation.mechanic.amount = 0.11764706f;
    value.template get<"modulation">().modulation.enabled = true;
    value.template get<"modulation">().modulation.evaluation = Pattern::Evaluation::DIRECT;
    value.template get<"modulation">().modulation.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    value.template get<"modulation">().modulation.source = Pattern::SourceKind::DIRECTIONAL_SINE;
    value.template get<"modulation">().modulation.noise.scale = 0.14453939f;
    value.template get<"modulation">().modulation.noise.speed = -0.0000010416667f;
    value.template get<"modulation">().modulation.sine.frequency = 1.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 0.0f;
    value.template get<"projection">().singularity_fade = 1.0f;
    value.template get<"source">().cell_radius = 0.48583883f;
    value.template get<"source">().cell_scale = 2.420276f;
    value.template get<"source">().cell_shape = 1.0f;
    value.template get<"source">().cell_softness = 0.22006065f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, GridWorms<W, H>,
                                    GridWormsSpec>::choreography_hooks()
        .template extend<
            GridWorms,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &GridWorms::initial_params>>();
  }
};
