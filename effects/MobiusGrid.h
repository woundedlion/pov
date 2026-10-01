/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

struct MobiusGridSpec : Pullback::Spec {
  static constexpr Pullback::ProjectionKind PROJECTION =
      Pullback::ProjectionKind::STEREOGRAPHIC;
  static constexpr Pullback::TransferKind TRANSFER =
      Pullback::TransferKind::NONE;
  static constexpr Pullback::ProjectionCoverageMode COVERAGE =
      Pullback::ProjectionCoverageMode::WEIGHT_SQUARED;
  static constexpr Pullback::FieldCoverageKind FIELD_COVERAGE =
      Pullback::FieldCoverageKind::NONE;
  static constexpr PaletteHarmony HARMONY = PaletteHarmony::COMPLEMENTARY;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::PATH_LENGTH;
  static constexpr Pullback::Color::BrightnessEnvelope BRIGHTNESS =
      Pullback::Color::BrightnessEnvelope::CUP;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr Pullback::SurfacePlacement SURFACE_PLACEMENT =
      Pullback::SurfacePlacement::BEFORE_LENS;
  using LensPolicy = void;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<
              Pullback::Lens::Mobius<Pullback::LensProvider<B>>>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::MirrorParams, B, "inner_warp",
          HUE == Pullback::HueMode::PATH_LENGTH>::Type>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<Pullback::TwinWaveSourceParams,
                                             B>::Type,
          Pullback::Weight::Projection,
          typename Pullback::CoveragePolicyFor<COVERAGE, B, void>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
using MobiusGridParams = Pullback::ParamsFor<MobiusGridSpec>;

/**
 * @brief A continuously animated Mobius lens over a mirrored twin wave.
 * @details Presets store an authored lens snapshot. Once the timeline steps,
 * the circular animation drives b on the unit circle; authored b is not a
 * rendered-frame promise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class MobiusGrid
    : public Pullback::ComposedEffect<W, H, MobiusGrid<W, H>, MobiusGridSpec> {

public:
  using Params = MobiusGridParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "mobius-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "eef46eb35e8a9bf17957a06530bf4409f8cb4ffd0c120b28d85eb0c65dfd314e";
  static constexpr std::string_view PRESET_BANK_DIGEST = "eb6c7f90002c374d45ecf50f0e9ca9070fa7ce80c9d79119883697659bbba82c";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{
      "mobius-grid",
      "mobius-grid-2"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;
  static constexpr bool ANIMATED_MOBIUS = true;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.template get<"projection">().camera_wander = 1.0f;
    value.template get<"color">().brightness_bottom = 0.0f;
    value.template get<"color">().brightness_top = 1.0f;
    value.template get<"color">().hue_shift_amount = 0.312f;
    value.template get<"color">().mapping_frequency = 1.0f;
    value.template get<"color">().mapping_phase = 0.0f;
    value.template get<"color">().palette_chroma = 0.398f;
    value.template get<"color">().palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.template get<"color">().phase_oscillation_depth = 0.0f;
    value.template get<"color">().phase_oscillation_speed = 0.0f;
    value.template get<"color">().opacity_high = 1.0f;
    value.template get<"color">().opacity_low = 1.0f;
    value.template get<"lens">().mobius.a.im = 0.304f;
    value.template get<"lens">().mobius.a.re = -1.072f;
    value.template get<"lens">().mobius.b.im = 0.0f;
    value.template get<"lens">().mobius.b.re = 0.416f;
    value.template get<"lens">().mobius.c.im = 0.0f;
    value.template get<"lens">().mobius.c.re = 0.0f;
    value.template get<"lens">().mobius.d.im = 0.0f;
    value.template get<"lens">().mobius.d.re = 0.70710677f;
    value.template get<"projection">().spin_rate = 0.0f;
    value.template get<"projection">().wander = 1.0f;
    value.template get<"projection">().singularity_fade = 2.102f;
    value.template get<"source">().angle_rate = 0.027f;
    value.template get<"source">().secondary_rate = 0.8f;
    value.template get<"source">().pattern_freq = 10.158f;
    value.template get<"source">().speed = 0.245f;
    value.template get<"inner_warp">().cell_x = 1.0f;
    value.template get<"inner_warp">().cell_y = 1.0f;
    value.template get<"inner_warp">().offset_x = 0.0f;
    value.template get<"inner_warp">().offset_y = 0.0f;
    value.template get<"inner_warp">().rotation = 0.0f;
    value.template get<"inner_warp">().speed = 0.0f;
    return value;
  }

  /** @brief The preset at index in PRESET_IDS and how it departs. */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.template get<"inner_warp">().cell_x = 0.2791094f;
      value.template get<"inner_warp">().cell_y = 6.810328f;
      value.template get<"inner_warp">().speed = 0.005875f;
    }
    return {value, Segue::Preset::Lerp{480, math::ease_in_out_sin}};
  }
  // clang-format on
  // End generated params.

  HS_COLD_MEMBER void after_composed_init() {
    this->start_mobius_animation(1.0f, 160);
  }
};
