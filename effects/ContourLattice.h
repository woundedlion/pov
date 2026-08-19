/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using ContourLatticeParams =
    Pullback::Params<Pullback::LatticeSourceParams, Pullback::AffineParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::IsoValueParams>;
using ContourLatticeSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED, void,
                   Pullback::TransferKind::ISO_CONTOUR,
                   Pullback::CoverageKind::PROJECTION>;

template <int W, int H>
class ContourLattice
    : public Pullback::ComposedEffect<
          W, H, ContourLattice<W, H>, ContourLatticeParams, ContourLatticeSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = ContourLatticeParams;
  static constexpr std::string_view EFFECT_ID = "contour-lattice";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "ae0310aed51fdae9dec8112fa507eff4edbb12156482e315f456e78ac424bbbf";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "174131893c44f3d99b1d8c7f6fe42dabf0b83ca15fd226f04928cae2cb82373f";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"affine-contour"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;

  static constexpr Params initial_params() {
    Params value;
    value.source.lattice_cell_scale = 1.22925f;
    value.source.lattice_shape_blend = 1.0f;
    value.source.lattice_softness = 0.1608203f;
    value.source.lattice_radius = 0.332981884f;
    value.projection.spin_rate = 0.0208791979f;
    value.projection.wander = 0.00309175253f;
    value.projection.camera_wander = 1.0f;
    value.outer_warp.speed = 0.015625f;
    value.outer_warp.translation_x = 4.0f;
    value.outer_warp.translation_y = 4.0f;
    value.value.iso_level = 0.138f;
    value.value.iso_width = 0.227034181f;
    value.color.hue_shift_amount = 0.398f;
    value.color.hue_noise_scale = 0.8300313f;
    value.color.hue_noise_speed = 0.000212000014f;
    return value;
  }
  static constexpr Params preset_params(size_t) { return initial_params(); }
};

#include "core/control/registry.h"
REGISTER_EFFECT(ContourLattice)
