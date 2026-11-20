/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/** @file JadeFiligree.h
 * @brief A golden mesh gives way to a flowing green filigree.
 */

#include "core/platform/effect_identity.h"
#include "core/render/pullback/composed_tiling.h"

#include <array>
#include <string_view>

/** @brief Authored stages for JadeFiligree. */
struct JadeFiligreeSpec : Pullback::Spec {
  static constexpr bool AUTHORED_TOPOLOGY = true;
  static constexpr auto TEMPORAL = Pullback::TemporalPolicy::phase();
  static constexpr auto PROJECTION_MOTION = Pullback::ProjectionMotion::STATIC;

  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Displace<typename Pullback::SurfacePolicyFor<
          Pullback::AuthoredDirectSurfaceParams, B, false, "displace">::Type>,
      Pullback::Stage::SampleSphere<typename Pullback::SourcePolicyFor<
          Pullback::Source::SphericalTilingSourceParams, B>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::AuthoredColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>>;
};

using JadeFiligreeParams = Pullback::ParamsFor<JadeFiligreeSpec>;

/** @brief A golden mesh gives way to a flowing green filigree.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 */
template <int W, int H>
class JadeFiligree : public Pullback::ComposedEffect<W, H, JadeFiligree<W, H>,
                                                     JadeFiligreeSpec> {
public:
  using Params = JadeFiligreeParams;

  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY = hs::effect_identities::JadeFiligree::VALUE;
  /// Canonical machine ID; equals the shader document's effect_id.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;
  /// SHA-256 hex of the shader document's parameter descriptor.
  static constexpr std::string_view DESCRIPTOR_DIGEST = "ff8577fff86d95b2e89b402fe7c860afd7b3688c4f5d251286c1430721396894";
  /// SHA-256 hex of the shader document's preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST = "8db8402081c4dc9a79bc9b5dc59c33b959c87a5cd630e64e997dd1b53c8f6bfa";
  /// Preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "golden-mesh",
      "jade-weave"
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
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.

  /** @brief Binds the authored preset bank. */
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<W, H, JadeFiligree<W, H>,
                                    JadeFiligreeSpec>::choreography_hooks()
        .template extend<
            JadeFiligree,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &JadeFiligree::initial_params>,
            Control::HookBinding<Control::ChoreographyHook::PRESET_ROW,
                                 &JadeFiligree::preset>>();
  }
};
