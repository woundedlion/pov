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
class MobiusGrid
    : public FixedLook::Runtime<
          W, H, MobiusGrid<W, H>,
          FixedLook::Params<FixedLook::TwinWaveSourceParams,
                            FixedLook::NoWarpParams, FixedLook::MirrorParams,
                            FixedLook::MobiusLensParams>,
          PaletteHarmony::COMPLEMENTARY, FixedLook::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::CUP> {
  using ParamsT =
      FixedLook::Params<FixedLook::TwinWaveSourceParams,
                        FixedLook::NoWarpParams, FixedLook::MirrorParams,
                        FixedLook::MobiusLensParams>;
  using Base = FixedLook::Runtime<W, H, MobiusGrid<W, H>, ParamsT,
                                  PaletteHarmony::COMPLEMENTARY,
                                  FixedLook::HueMode::PATH_LENGTH,
                                  Pullback::Color::BrightnessEnvelope::CUP>;

public:
  using Params = ParamsT;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "mobius-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "b7274aa41ad3c6f8ef3d75c1c978ba9f1ec03030144b02d06ba7964a0afe02e0";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "973d13d39e1203831d581ceada1de3dce6786552b4ae71f4db9987a137a4d801";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"mobius-grid",
                                                              "mobius-grid-2"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr bool ANIMATED_MOBIUS = true;
  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding,
                                   FixedLook::OuterCameraProvider<Binding>>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity,
      Pullback::Lens::Mobius<FixedLook::LensProvider<Binding>>,
      Pullback::Surface::Identity,
      Pullback::Projection::Stereographic<
          FixedLook::ProjectionProvider<Binding>>>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding, Pullback::Warp::Identity,
      Pullback::Warp::MirrorTile<
          FixedLook::WarpProvider<Binding, false, true>>>;
  using SourceStage = Pullback::Stage::Source<
      Binding, Pullback::Source::TwinWave<FixedLook::SourceProvider<Binding>>>;
  using MaterialStage =
      Pullback::Stage::Material<Binding, Pullback::Weight::Projection,
                                Pullback::Transfer::Linear,
                                Pullback::Coverage::ProjectionSquared>;
  using ColorStage = Pullback::Stage::Color<
      Binding, Pullback::Color::GeneratedPalette<FixedLook::ColorProvider<
                   Binding, FixedLook::HueMode::PATH_LENGTH,
                   Pullback::Color::BrightnessEnvelope::CUP>>>;
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
    value.source = {10.158f, 0.245f, 0.8f, 0.027f};
    value.projection.pole_fade = 2.102f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.lens.mobius = {-1.072f, 0.304f, 0.416f,      0.0f,
                         0.0f,    0.0f,   0.70710677f, 0.0f};
    value.color.hue_shift_amount = 0.312f;
    value.color.palette_chroma = 0.398f;
    return value;
  }
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1) {
      value.inner_warp.speed = 0.005875f;
      value.inner_warp.cell_x = 0.2791094f;
      value.inner_warp.cell_y = 6.810328f;
    }
    return value;
  }

public:
  HS_COLD_MEMBER void after_fixed_init() {
    this->start_mobius_animation(1.0f, 160);
  }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MobiusGrid)
