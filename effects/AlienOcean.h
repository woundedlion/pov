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
class AlienOcean
    : public FixedLook::Runtime<
          W, H, AlienOcean<W, H>,
          FixedLook::Params<FixedLook::GridSourceParams,
                            FixedLook::MirrorParams, FixedLook::NoWarpParams,
                            FixedLook::NoLensParams,
                            FixedLook::EdgeValueParams>,
          PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::GridSourceParams, FixedLook::MirrorParams,
                        FixedLook::NoWarpParams, FixedLook::NoLensParams,
                        FixedLook::EdgeValueParams>;
  using Base =
      FixedLook::Runtime<W, H, AlienOcean<W, H>, ParamsT,
                         PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "alien-ocean";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "b6c0f04a02550a66af8c0f140822ddeccfbf110aa3c57bc6631397d0f7672a8d";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "85c8b3961f3d1de228ff43e90f9f0e957f453db2b30f2cd98222b305ab0e77f6";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"folded-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;

  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding,
                                   FixedLook::OuterCameraProvider<Binding>>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity, Pullback::Lens::Kaleidoscope,
      Pullback::Surface::Identity,
      Pullback::Projection::Gnomonic<
          FixedLook::ProjectionProvider<Binding>,
          Pullback::Projection::GnomonicHemisphere::FOLDED>>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding,
      Pullback::Warp::MirrorTile<FixedLook::WarpProvider<Binding, true>>,
      Pullback::Warp::Identity>;
  using SourceStage = Pullback::Stage::Source<
      Binding, Pullback::Source::Grid<FixedLook::SourceProvider<Binding>>>;
  using MaterialStage = Pullback::Stage::Material<
      Binding, Pullback::Weight::Projection, Pullback::Transfer::Linear,
      Pullback::Coverage::EdgeFade<FixedLook::ValueProvider<Binding>>>;
  using ColorStage = Pullback::Stage::Color<
      Binding, Pullback::Color::GeneratedPalette<FixedLook::ColorProvider<
                   Binding, FixedLook::HueMode::NOISE,
                   Pullback::Color::BrightnessEnvelope::NONE>>>;
  using RenderPipeline =
      Pullback::Pipeline<Binding, OuterCameraStage, SurfaceStage,
                         PlanarWarpStage, SourceStage, MaterialStage,
                         ColorStage>;
  using FrameState = typename RenderPipeline::Frame;

  static HS_FLASH_INLINE Color4 shade(const Vector &view,
                                      const FrameState &frame) {
    return RenderPipeline::shade(view, frame);
  }

  static constexpr Params initial_params() {
    Params value;
    value.source = {3.565f, 0.235f, 0.0f, 1.0f, 1.0f, 0.0f};
    value.projection.pole_fade = 1.4f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.rotation = 0.29530972f;
    value.outer_warp.cell_x = 5.381125f;
    value.outer_warp.offset_x = 1.344f;
    value.outer_warp.offset_y = -1.456f;
    value.value.edge_width = 0.5f;
    value.color.hue_shift_amount = 0.424f;
    value.color.hue_noise_scale = 2.2033439f;
    value.color.palette_chroma = 0.4f;
    return value;
  }

  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(AlienOcean)
