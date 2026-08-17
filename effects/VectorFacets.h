/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "effects/fixed/FixedLookRuntime.h"

template <int W, int H>
class VectorFacets : public FixedLook::Runtime<
                         W, H, VectorFacets<W, H>,
                         FixedLook::Params<FixedLook::GridSourceParams,
                                           FixedLook::VectorNoiseParams,
                                           FixedLook::MirrorParams>,
                         PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::CUP, true> {
  using ParamsT =
      FixedLook::Params<FixedLook::GridSourceParams,
                        FixedLook::VectorNoiseParams, FixedLook::MirrorParams>;
  using Base =
      FixedLook::Runtime<W, H, VectorFacets<W, H>, ParamsT,
                         PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::CUP, true>;

public:
  using Params = ParamsT;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "vector-facets";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "97a8dd884898f8397325ee1b8fd3aa10e0cb1ad81a1a601d75f59d9709f1b161";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "240f9def2290da2810045910f102f68ae284368ae6fd8607736b6aed4b501d54";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"vector-mirror"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;
  static constexpr int32_t OUTER_NOISE_SEED = 1337;
  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding,
                                   FixedLook::OuterCameraProvider<Binding>>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity,
      Pullback::Lens::DodecahedralKaleidoscope, Pullback::Surface::Identity,
      Pullback::Projection::Gnomonic<
          FixedLook::ProjectionProvider<Binding>,
          Pullback::Projection::GnomonicHemisphere::FOLDED>>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding,
      Pullback::Warp::VectorNoise<FixedLook::WarpProvider<Binding, true>,
                                  NoiseBasis::SIMPLEX,
                                  Pullback::Warp::FlatEnvelope>,
      Pullback::Warp::MirrorTile<FixedLook::WarpProvider<Binding, false>>>;
  using SourceStage = Pullback::Stage::Source<
      Binding, Pullback::Source::Grid<FixedLook::SourceProvider<Binding>>>;
  using MaterialStage =
      Pullback::Stage::Material<Binding, Pullback::Weight::Projection,
                                Pullback::Transfer::Linear,
                                Pullback::Coverage::ProjectionSquared>;
  using ColorStage = Pullback::Stage::Color<
      Binding, Pullback::Color::GeneratedPalette<FixedLook::ColorProvider<
                   Binding, FixedLook::HueMode::NOISE,
                   Pullback::Color::BrightnessEnvelope::CUP>>>;
  using RenderPipeline =
      Pullback::Pipeline<Binding, OuterCameraStage, SurfaceStage,
                         PlanarWarpStage, SourceStage, MaterialStage,
                         ColorStage>;
  using FrameState = typename RenderPipeline::Frame;
  static HS_FLASH_MEMBER Color4 shade(const Vector &view,
                                      const FrameState &frame) {
    return RenderPipeline::shade(view, frame);
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

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(VectorFacets)
