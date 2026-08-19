/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using VectorFacetsParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::VectorNoiseParams,
                     Pullback::MirrorParams>;
using VectorFacetsSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::LINEAR,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class VectorFacets
    : public Pullback::ComposedEffect<
          W, H, VectorFacets<W, H>, VectorFacetsParams, VectorFacetsSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::CUP, true> {

public:
  using Params = VectorFacetsParams;
  static constexpr std::string_view EFFECT_ID = "vector-facets";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "427d1301321002d21596847a73057cd9661c915d9bda361196e39d3aa398b922";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "fbafc54e0d5e87962f8763c86ca321acfc1c27a6ae4b0964786d2a9b1dc6314a";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"vector-mirror"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;
  static constexpr int32_t OUTER_NOISE_SEED = 1337;

  static HS_FLASH_MEMBER Color4
  shade(const Vector &view, const typename VectorFacets::FrameState &frame) {
    return VectorFacets::RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {4.9755f, 0.04f, 1.704f, 0.0f, 0.8f, 0.027f};
    value.projection.pole_fade = 2.311f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.strength = 0.138f;
    value.outer_warp.speed = -0.00005f;
    value.inner_warp.speed = 0.00327999983f;
    value.color.hue_shift_amount = 0.721f;
    value.color.palette_chroma = 1.0f;
    value.color.brightness_depth = 0.655f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/control/registry.h"
REGISTER_EFFECT(VectorFacets)
