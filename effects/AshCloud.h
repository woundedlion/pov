/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using AshCloudParams =
    Pullback::Params<Pullback::LatticeSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::CutoutValueParams, Pullback::SurfaceNoiseParams>;
using AshCloudSpec = Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                                    Pullback::Lens::DodecahedralKaleidoscope,
                                    Pullback::TransferKind::NONE,
                                    Pullback::ProjectionCoverageMode::WEIGHT,
                                    Pullback::FieldCoverageKind::VALUE_CUTOUT>;

/**
 * @brief A softly cut lattice curled across dodecahedral facets.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AshCloud : public Pullback::ComposedEffect<
                     W, H, AshCloud<W, H>, AshCloudParams, AshCloudSpec,
                     PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
                     Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = AshCloudParams;
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "ash-cloud";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "6bde4043e907674b80f6b12bf66c873366359c9a7b1996e8804ae24da673518e";
  static constexpr std::string_view PRESET_BANK_DIGEST = "bb8b1dc6fe15beaf9a4c11730c902d9743a3035942e3b123e3760e25e50c214e";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "ash-cloud"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr float CAMERA_SPIN_RATE = 0.01975f;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.hue_noise_scale = 2.4029531f;
    value.color.hue_noise_speed = 0.00015f;
    value.color.hue_shift_amount = -3.248f;
    value.color.mapping_frequency = 3.666f;
    value.color.mapping_phase = 0.598f;
    value.color.palette_chroma = 0.314f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.phase_oscillation_depth = 0.856f;
    value.color.phase_oscillation_speed = 0.00448f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.value.cutout_softness = 0.31635937f;
    value.value.cutout_threshold = 0.365f;
    value.projection.spin_rate = 0.00175f;
    value.projection.wander = 0.009f;
    value.projection.singularity_fade = 1.0f;
    value.source.lattice_cell_scale = 0.98971874f;
    value.source.lattice_radius = 0.37757313f;
    value.source.lattice_shape_blend = 0.0f;
    value.source.lattice_softness = 0.34563965f;
    value.surface.scale = 5.7742186f;
    value.surface.speed = 0.00040625f;
    value.surface.strength = 0.043f;
    return value;
  }
  // clang-format on
  // End generated params.
};
