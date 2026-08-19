/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using PrismLatticeParams =
    Pullback::Params<Pullback::LatticeSourceParams, Pullback::PolarParams,
                     Pullback::WaveShearParams>;
using PrismLatticeSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::PentagonalPrismKaleidoscope,
                   Pullback::TransferKind::LINEAR,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class PrismLattice
    : public Pullback::ComposedEffect<
          W, H, PrismLattice<W, H>, PrismLatticeParams, PrismLatticeSpec,
          PaletteHarmony::ANALOGOUS, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = PrismLatticeParams;
  static constexpr std::string_view EFFECT_ID = "prism-lattice";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "d1c44c9ef778278e71a9624f1735d254e9744523b4f816dc878acbafeab89a63";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "1bb66847b71e59e0067d82149ae9028b615bc20b5cc4112e2efe29acad3eb756";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"polar-wave"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;

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

#include "core/control/registry.h"
REGISTER_EFFECT(PrismLattice)
