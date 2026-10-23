/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file KaleidoscopeHexSoft.h
 * @brief Composed pullback effect:
 *        a drifting twin wave reflected through a kaleidoscope.
 */

#include "core/platform/effect_identity.h"
#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

/** @brief Stage spec for `KaleidoscopeHexSoft`. */
struct KaleidoscopeHexSoftSpec : Pullback::Spec {
  /// Samples animation continuously across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = true;
  /// Whether each displayed column is followed by black.
  static constexpr bool STROBE_COLUMNS = false;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::TRIADIC;
  static constexpr bool ANIMATED_PROJECTION = true;
  /**
   * @brief Stage chain: camera rotation, kaleidoscope lens, stereographic projection, mirror warp, twin-wave source, generated palette.
   * @tparam B The effect's Binding.
   */
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<Pullback::Lens::Kaleidoscope>,
          Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
              Pullback::ProjectionKind::STEREOGRAPHIC, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::MirrorParams, B, "inner_warp", false>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::TwinWaveSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<
              Pullback::ProjectionCoverageMode::WEIGHT_SQUARED, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, Pullback::HueMode::NOISE,
                                  Pullback::Color::BrightnessEnvelope::NONE>>>>;
};
/// Parameter block derived from `KaleidoscopeHexSoftSpec`.
using KaleidoscopeHexSoftParams = Pullback::ParamsFor<KaleidoscopeHexSoftSpec>;

/**
 * @brief A drifting twin wave reflected through a kaleidoscope.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class KaleidoscopeHexSoft
    : public Pullback::ComposedEffect<W, H, KaleidoscopeHexSoft<W, H>,
                                      KaleidoscopeHexSoftSpec> {

public:
  using Params = KaleidoscopeHexSoftParams; ///< Live parameter struct.
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::KaleidoscopeHexSoft::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "376f024ceb9f880d96e96d495f397bef832f75d19768ee5ad0bd0c290f0ac30e";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "a54f426490227b6ad95b64d0a18bbd0bb65beaab6fb4053f37f85bb435a6e7ed";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "twin-wave"
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
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().hue_noise_scale = 2.2033439f;
    value.template get<"color">().hue_noise_speed = -0.00040800002f;
    value.template get<"color">().hue_shift_amount = 0.27f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.361f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 1.0f;
    value.template get<"projection">().singularity_fade = 4.971f;
    value.template get<"source">().angle_rate = 0.05f;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 4.9755f;
    value.template get<"source">().speed = 0.125f;
    value.template get<"inner_warp">().cell_x = 1.0f;
    value.template get<"inner_warp">().cell_y = 1.0f;
    value.template get<"inner_warp">().offset_x = 0.0f;
    value.template get<"inner_warp">().offset_y = 0.0f;
    value.template get<"inner_warp">().rotation = 0.0f;
    value.template get<"inner_warp">().speed = 0.0f;
    return value;
  }
  // clang-format on
  // End generated params.

public:
  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<
               W, H, KaleidoscopeHexSoft<W, H>,
               KaleidoscopeHexSoftSpec>::choreography_hooks()
        .template extend<
            KaleidoscopeHexSoft,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &KaleidoscopeHexSoft::initial_params>>();
  }
};
