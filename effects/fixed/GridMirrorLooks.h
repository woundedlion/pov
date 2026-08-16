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
class GlitchGrid
    : public FixedLook::Runtime<
          W, H, GlitchGrid<W, H>,
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
      FixedLook::Runtime<W, H, GlitchGrid<W, H>, ParamsT,
                         PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "glitch-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "082daf2081c9c9f3a59c35d9d8c1faa4c531c23f1da4b6ff2e04ea9d66eed0a7";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "80b3716f38910adb47b3e7c15370d9132083dc94bc549556c1fb540dbf420320";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"folded-glitch"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;
  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding,
                                   FixedLook::OuterCameraProvider<Binding>>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity, Pullback::Lens::Glitch,
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
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

template <int W, int H>
class CosmicEyeball
    : public FixedLook::Runtime<
          W, H, CosmicEyeball<W, H>,
          FixedLook::Params<FixedLook::GridSourceParams,
                            FixedLook::MirrorParams, FixedLook::NoWarpParams,
                            FixedLook::NoLensParams,
                            FixedLook::EdgeValueParams>,
          PaletteHarmony::TRIADIC, FixedLook::HueMode::PATH_LENGTH,
          Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::GridSourceParams, FixedLook::MirrorParams,
                        FixedLook::NoWarpParams, FixedLook::NoLensParams,
                        FixedLook::EdgeValueParams>;
  using Base = FixedLook::Runtime<W, H, CosmicEyeball<W, H>, ParamsT,
                                  PaletteHarmony::TRIADIC,
                                  FixedLook::HueMode::PATH_LENGTH,
                                  Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "cosmic-eyeball";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "8cf940ce43ab8b2e1b583ecbd9e396ed1295b07339cefb135af82e6f2bc610de";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "a9e60a73e0e553691db4800ed847bb02acc5ef9f4617923c03b11f36a7425ff9";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"mirrored-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;
  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding,
                                   FixedLook::OuterCameraProvider<Binding>>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity, Pullback::Lens::Glitch,
      Pullback::Surface::Identity,
      Pullback::Projection::Stereographic<
          FixedLook::ProjectionProvider<Binding>>>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding,
      Pullback::Warp::MirrorTile<FixedLook::WarpProvider<Binding, true, true>>,
      Pullback::Warp::Identity>;
  using SourceStage = Pullback::Stage::Source<
      Binding, Pullback::Source::Grid<FixedLook::SourceProvider<Binding>>>;
  using MaterialStage = Pullback::Stage::Material<
      Binding, Pullback::Weight::Projection, Pullback::Transfer::Linear,
      Pullback::Coverage::EdgeFade<FixedLook::ValueProvider<Binding>>>;
  using ColorStage = Pullback::Stage::Color<
      Binding, Pullback::Color::GeneratedPalette<FixedLook::ColorProvider<
                   Binding, FixedLook::HueMode::PATH_LENGTH,
                   Pullback::Color::BrightnessEnvelope::NONE>>>;
  using RenderPipeline =
      Pullback::Pipeline<Binding, OuterCameraStage, SurfaceStage,
                         PlanarWarpStage, SourceStage, MaterialStage,
                         ColorStage>;
  static HS_FLASH_INLINE Color4 shade(const Vector &view,
                                      const FrameState &frame) {
    return RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {2.5477f, 0.235f, 1.854f, 0.0f, 1.0f, 0.0f};
    value.projection.pole_fade = 1.4f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.rotation = 0.295309722f;
    value.outer_warp.cell_x = 5.381125f;
    value.outer_warp.offset_x = 1.344f;
    value.outer_warp.offset_y = -1.456f;
    value.value.edge_width = 0.5f;
    value.color.hue_shift_amount = 2.048f;
    value.color.palette_chroma = 0.292f;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

template <int W, int H>
class EquatorGrid
    : public FixedLook::Runtime<
          W, H, EquatorGrid<W, H>,
          FixedLook::Params<FixedLook::GridSourceParams,
                            FixedLook::NoWarpParams, FixedLook::MirrorParams>,
          PaletteHarmony::ANALOGOUS, FixedLook::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::GridSourceParams, FixedLook::NoWarpParams,
                        FixedLook::MirrorParams>;
  using Base =
      FixedLook::Runtime<W, H, EquatorGrid<W, H>, ParamsT,
                         PaletteHarmony::ANALOGOUS, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "equator-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "d7422bfa1888ee66d90b8421cf02041cc32f72ca7cfd65a0e2ef5083d246cd03";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "278b43d19c17acc73e70c8f1879e9f5a6e9fc40e4f67355323adc926b32fe59b";
  static constexpr std::array<std::string_view, 3> PRESET_IDS{
      "double-map", "open-grid", "fine-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr bool USES_CENTRAL_MERIDIAN = true;
  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding,
                                   FixedLook::OuterCameraProvider<Binding>>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity,
      Pullback::Lens::DodecahedralKaleidoscope, Pullback::Surface::Identity,
      Pullback::Projection::Equirectangular<
          FixedLook::ProjectionProvider<Binding>>>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding, Pullback::Warp::Identity,
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
                   Pullback::Color::BrightnessEnvelope::NONE>>>;
  using RenderPipeline =
      Pullback::Pipeline<Binding, OuterCameraStage, SurfaceStage,
                         PlanarWarpStage, SourceStage, MaterialStage,
                         ColorStage>;
  static HS_FLASH_INLINE Color4 shade(const Vector &view,
                                      const FrameState &frame) {
    return RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {3.9407f, 0.0f, 3.0f, 1.0f, 0.8f, 0.0269999988f};
    value.projection.pole_fade = 2.14f;
    value.projection.wander = 0.165f;
    value.projection.camera_wander = 1.0f;
    value.inner_warp.speed = 0.00013f;
    value.inner_warp.cell_y = 0.997703135f;
    value.color.hue_shift_amount = 0.366f;
    value.color.hue_noise_scale = 1.47215629f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_frequency = 2.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1)
      value.color.mapping_frequency = 1.0f;
    else if (index == 2) {
      value.source.pattern_freq = 0.3985f;
      value.inner_warp.speed = 0.00058f;
      value.inner_warp.cell_y = 0.901890635f;
      value.color.mapping_frequency = 21.212f;
    }
    return value;
  }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(GlitchGrid)
REGISTER_EFFECT(CosmicEyeball)
REGISTER_EFFECT(EquatorGrid)
