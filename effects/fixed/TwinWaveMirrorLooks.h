/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "effects/fixed/FixedLookRuntime.h"

#define HS_FIXED_TWIN_PIPELINE(LensPolicy, TrackPath)                          \
  using OuterCameraStage =                                                     \
      Pullback::Stage::OuterCamera<Binding,                                    \
                                   FixedLook::OuterCameraProvider<Binding>>;   \
  using SurfaceStage = Pullback::Stage::SurfaceProject<                        \
      Binding, Pullback::Surface::Identity, LensPolicy,                        \
      Pullback::Surface::Identity,                                             \
      Pullback::Projection::Stereographic<                                     \
          FixedLook::ProjectionProvider<Binding>>>;                            \
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<                         \
      Binding, Pullback::Warp::Identity,                                       \
      Pullback::Warp::MirrorTile<                                              \
          FixedLook::WarpProvider<Binding, false, TrackPath>>>;                \
  using SourceStage = Pullback::Stage::Source<                                 \
      Binding,                                                                 \
      Pullback::Source::TwinWave<FixedLook::SourceProvider<Binding>>>;         \
  using MaterialStage =                                                        \
      Pullback::Stage::Material<Binding, Pullback::Weight::Projection,         \
                                Pullback::Transfer::Linear,                    \
                                Pullback::Coverage::ProjectionSquared>

template <int W, int H>
class KaleidoWave
    : public FixedLook::Runtime<
          W, H, KaleidoWave<W, H>,
          FixedLook::Params<FixedLook::TwinWaveSourceParams,
                            FixedLook::NoWarpParams, FixedLook::MirrorParams>,
          PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::TwinWaveSourceParams,
                        FixedLook::NoWarpParams, FixedLook::MirrorParams>;
  using Base =
      FixedLook::Runtime<W, H, KaleidoWave<W, H>, ParamsT,
                         PaletteHarmony::TRIADIC, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "kaleido-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "605c55957f751d7baa43cd24ecf9d29ad57ff579657add5cea44e43b1214691e";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "3ea995834789f58dc5f71c3000fa7c94a3ff864c2e38c8ef8b2c3091dfe54a3f";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"twin-wave"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  HS_FIXED_TWIN_PIPELINE(Pullback::Lens::Kaleidoscope, false);
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
    value.source = {4.9755f, 0.125f, 0.8f, 0.05f};
    value.projection.pole_fade = 4.971f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.color.hue_shift_amount = 0.27f;
    value.color.hue_noise_scale = 2.2033439f;
    value.color.hue_noise_speed = -0.00040800002f;
    value.color.palette_chroma = 0.361f;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

template <int W, int H>
class HexWave
    : public FixedLook::Runtime<
          W, H, HexWave<W, H>,
          FixedLook::Params<FixedLook::TwinWaveSourceParams,
                            FixedLook::NoWarpParams, FixedLook::MirrorParams>,
          PaletteHarmony::ANALOGOUS, FixedLook::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::TwinWaveSourceParams,
                        FixedLook::NoWarpParams, FixedLook::MirrorParams>;
  using Base =
      FixedLook::Runtime<W, H, HexWave<W, H>, ParamsT,
                         PaletteHarmony::ANALOGOUS, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "hex-wave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "7c17e826920d96ebe397cd7939bbcc0a31bc431b6ff6545d99472a280ba5cbff";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "8d76e45aa0e6a6c76644991541bf513c531ebf423112fb36cdc93fc22516df87";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"hex-twin-wave"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  HS_FIXED_TWIN_PIPELINE(Pullback::Lens::HexagonalPrismKaleidoscope, false);
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
    value.source = {3.881f, 0.128598228f, 0.8f, 0.027f};
    value.projection.pole_fade = 4.971f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.color.hue_shift_amount = 0.226f;
    value.color.hue_noise_scale = 1.47215629f;
    value.color.hue_noise_speed = 0.000138f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_frequency = 1.341f;
    value.color.mapping_phase = -1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::BELL;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

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
  using FrameState = typename Base::Frame;
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
  HS_FIXED_TWIN_PIPELINE(
      Pullback::Lens::Mobius<FixedLook::LensProvider<Binding>>, true);
  using ColorStage = Pullback::Stage::Color<
      Binding, Pullback::Color::GeneratedPalette<FixedLook::ColorProvider<
                   Binding, FixedLook::HueMode::PATH_LENGTH,
                   Pullback::Color::BrightnessEnvelope::CUP>>>;
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

#undef HS_FIXED_TWIN_PIPELINE

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(KaleidoWave)
REGISTER_EFFECT(HexWave)
REGISTER_EFFECT(MobiusGrid)
