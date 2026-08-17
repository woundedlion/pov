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
class PrismLattice
    : public FixedLook::Runtime<
          W, H, PrismLattice<W, H>,
          FixedLook::Params<FixedLook::LatticeSourceParams,
                            FixedLook::PolarParams, FixedLook::WaveShearParams>,
          PaletteHarmony::ANALOGOUS, FixedLook::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {
  using ParamsT =
      FixedLook::Params<FixedLook::LatticeSourceParams, FixedLook::PolarParams,
                        FixedLook::WaveShearParams>;
  using Base =
      FixedLook::Runtime<W, H, PrismLattice<W, H>, ParamsT,
                         PaletteHarmony::ANALOGOUS, FixedLook::HueMode::NOISE,
                         Pullback::Color::BrightnessEnvelope::NONE>;

public:
  using Params = ParamsT;
  using FrameState = typename Base::Frame;
  using Binding = typename Base::PipelineBinding;
  static constexpr std::string_view EFFECT_ID = "prism-lattice";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "541b1c7919eaded5d18ce350d120634800aad9df577c21451e90cf5863b461e6";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "1804bcfac7caac6fb58e7d4129dd922442a19cd2395832e1ce03c9e41899d064";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"polar-wave"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding,
                                   FixedLook::OuterCameraProvider<Binding>>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity,
      Pullback::Lens::PentagonalPrismKaleidoscope, Pullback::Surface::Identity,
      Pullback::Projection::Stereographic<
          FixedLook::ProjectionProvider<Binding>>>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding,
      Pullback::Warp::PolarChart<FixedLook::WarpProvider<Binding, true>,
                                 Pullback::Warp::LinearPolar, 1>,
      Pullback::Warp::WaveShear<FixedLook::WarpProvider<Binding, false>>>;
  using SourceStage = Pullback::Stage::Source<
      Binding,
      Pullback::Source::PrimitiveLattice<FixedLook::SourceProvider<Binding>>>;
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
    value.source.lattice_cell_scale = 0.774140596f;
    value.source.lattice_shape_blend = 1.0f;
    value.source.lattice_softness = 0.377608389f;
    value.source.lattice_radius = 0.290762514f;
    value.projection.pole_fade = 2.273f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.speed = 0.000343749998f;
    value.inner_warp.speed = 0.000999999931f;
    value.inner_warp.strength = 0.0f;
    value.color.hue_shift_amount = 0.268000007f;
    value.color.hue_noise_scale = 2.0f;
    value.color.palette_chroma = 1.0f;
    value.color.mapping_phase = -0.165999994f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(PrismLattice)
