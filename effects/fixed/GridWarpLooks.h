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
class SignalWeave
    : public FixedLook::Runtime<W, H, SignalWeave<W, H>,
                                FixedLook::Params<FixedLook::GridSourceParams,
                                                  FixedLook::WaveShearParams,
                                                  FixedLook::NoWarpParams>,
                                PaletteHarmony::TRIADIC,
                                FixedLook::HueMode::NOISE,
                                Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::GridSourceParams, FixedLook::WaveShearParams,
                        FixedLook::NoWarpParams>;
  using Base =
      FixedLook::Runtime<W, H, SignalWeave<W, H>, ParamsT,
                         PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "signal-weave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "2e37f7346be3aaa346492bb217b9f69696e25669863648000703d2f1792edb7c";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "fa21b459874c5b9f709b021bdfcfd91b983a0cbc4d72a396a6d9badca3ec7794";
  static constexpr uint16_t INITIAL_PRESET_DWELL_FRAMES = 120;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr std::array<std::string_view, 4> PRESET_IDS{
      "signal-weave", "signal-weave-2", "signal-weave-3", "signal-weave-4"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr bool ANIMATED_PROJECTION = true;
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
      Pullback::Warp::WaveShear<FixedLook::WarpProvider<Binding, true>>,
      Pullback::Warp::Identity>;
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
  static HS_HOT_FLASH_MEMBER Color4 shade(const Vector &view,
                                          const FrameState &frame) {
    return RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() { return make_params(0); }
  static constexpr Params preset_params(size_t index) {
    return make_params(index);
  }

public:
  HS_COLD_MEMBER void after_fixed_init() {
    this->hold_initial_preset(INITIAL_PRESET_DWELL_FRAMES);
  }

private:
  static constexpr Params make_params(size_t index) {
    constexpr float frequencies[] = {4.439f, 3.1447f, 7.5227f, 8.8162f};
    constexpr float complexities[] = {0.5f, 0.5f, 1.698f, 1.698f};
    constexpr float strengths[] = {0.5f, 2.72f, 0.0f, 1.376f};
    constexpr float speeds[] = {0.015625f, 0.00690625f, 0.00690625f,
                                0.00559375f};
    Params value;
    value.source = {
        frequencies[index], 0.245f, complexities[index], 0.0f, 0.0f, 0.0f};
    value.projection.pole_fade = 1.0f;
    value.projection.camera_wander = 0.8f;
    value.outer_warp.strength = strengths[index];
    value.outer_warp.speed = speeds[index];
    value.color.hue_shift_amount = 0.292f;
    value.color.hue_noise_scale = 0.6304219f;
    value.color.palette_chroma = 0.788f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
};

template <int W, int H>
class FacetWave
    : public FixedLook::Runtime<W, H, FacetWave<W, H>,
                                FixedLook::Params<FixedLook::GridSourceParams,
                                                  FixedLook::WaveShearParams,
                                                  FixedLook::MirrorParams>,
                                PaletteHarmony::TRIADIC,
                                FixedLook::HueMode::NOISE,
                                Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::GridSourceParams, FixedLook::WaveShearParams,
                        FixedLook::MirrorParams>;
  using Base =
      FixedLook::Runtime<W, H, FacetWave<W, H>, ParamsT,
                         PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "facet-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "ee241e340b3f0133d2acc8d6f7b2f6c6926e7132127349035c1d2321323d4ae3";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "8c9a7d04c494ead43710bda0ad354acbf02e4d91279efc6e21469e66057060b0";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"wave-mirror"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = false;
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
      Pullback::Warp::WaveShear<FixedLook::WarpProvider<Binding, true>>,
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
  static HS_HOT_FLASH_MEMBER Color4 shade(const Vector &view,
                                          const FrameState &frame) {
    return RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() {
    Params value;
    value.source = {6.3287f, 0.04f, 1.704f, 0.0f, 0.8f, 0.027f};
    value.projection.pole_fade = 2.311f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.strength = -0.176f;
    value.outer_warp.speed = -0.00325f;
    value.outer_warp.frequency = 1.408f;
    value.outer_warp.field_angle = 2.2305307f;
    value.color.hue_shift_amount = 0.721f;
    value.color.palette_chroma = 1.0f;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

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
  using FrameState = typename Base::Frame;
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
REGISTER_EFFECT(SignalWeave)
REGISTER_EFFECT(FacetWave)
REGISTER_EFFECT(VectorFacets)
