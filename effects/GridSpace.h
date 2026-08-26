/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using GridSpaceParams =
    Pullback::Params<Pullback::LatticeSourceParams, Pullback::AffineParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::IsoValueParams>;
using GridSpaceSpec = Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                                     void, Pullback::TransferKind::ISO_CONTOUR,
                                     Pullback::ProjectionCoverageMode::WEIGHT>;

/**
 * @brief An affine primitive lattice rendered as soft contours.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class GridSpace : public Pullback::ComposedEffect<
                      W, H, GridSpace<W, H>, GridSpaceParams, GridSpaceSpec,
                      PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
                      Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = GridSpaceParams;
  static constexpr std::string_view EFFECT_ID = "grid-space";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "738a1214551232ef6d4a87dffad982d1f3d2d3a43a0a192c8cd7601417ee2475";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "2229900cf32e7beedb2f68925a97a71598abc3b045fdbda07d8c0ad8cab9f32b";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"affine-contour"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

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
};

#include "core/control/registry.h"
REGISTER_EFFECT(GridSpace)
