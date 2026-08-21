/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using CurlFacetsParams =
    Pullback::Params<Pullback::LatticeSourceParams, Pullback::NoWarpParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::CutoutValueParams, Pullback::SurfaceNoiseParams>;
using CurlFacetsSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION,
                   Pullback::FieldCoverageKind::VALUE_CUTOUT>;

template <int W, int H>
class CurlFacets : public Pullback::ComposedEffect<
                       W, H, CurlFacets<W, H>, CurlFacetsParams, CurlFacetsSpec,
                       PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
                       Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = CurlFacetsParams;
  static constexpr std::string_view EFFECT_ID = "curl-facets";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "c85fb6ca75e4ec2833ae441dd32497fc30353adcd2e47a802349b747885523be";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "21158c87c54f32be2cc89f4531d2619cbc325dc7261f110a03364dc652cfc9fe";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"curl-facets"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr float CAMERA_SPIN_RATE = 0.0197500009f;
  static constexpr int32_t SURFACE_NOISE_SEED = Pullback::EFFECT_NOISE_SEED;

  static constexpr Params initial_params() {
    Params value;
    value.source = {.lattice_cell_scale = 0.989718735f,
                    .lattice_shape_blend = 0.0f,
                    .lattice_softness = 0.345639646f,
                    .lattice_radius = 0.377573133f};
    value.projection.pole_fade = 1.0f;
    value.projection.spin_rate = 0.00175000005f;
    value.projection.wander = 0.00899999961f;
    value.projection.camera_wander = 1.0f;
    value.surface = {.scale = 5.77421856f,
                     .strength = 0.0430000015f,
                     .speed = 0.000406250008f};
    value.value = {.cutout_threshold = 0.36500001f,
                   .cutout_softness = 0.316359371f};
    value.color.hue_shift_amount = -3.24799991f;
    value.color.hue_noise_scale = 2.40295315f;
    value.color.hue_noise_speed = 0.000150000007f;
    value.color.palette_chroma = 0.31400001f;
    value.color.mapping_frequency = 3.66599989f;
    value.color.mapping_phase = 0.59799999f;
    value.color.phase_oscillation_depth = 0.856000006f;
    value.color.phase_oscillation_speed = 0.00448000012f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(CurlFacets)
