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
  // Generated identity: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr std::string_view EFFECT_ID = "grid-space";
  static constexpr std::string_view DESCRIPTOR_DIGEST = "1e95a18ef294231c38dd3e793ee2c92849f7804ac68bbf61868672feb91f4528";
  static constexpr std::string_view PRESET_BANK_DIGEST = "992256905cd236a8dae29bbafc058efee46c7a98751abcf477ad1622ef0dd535";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{
      "affine-contour"
  };
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  // clang-format on
  // End generated identity.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;

  // Generated params: scripts/generate_composed_presets.mjs
  // clang-format off
  static constexpr Params initial_params() {
    Params value;
    value.projection.camera_wander = 1.0f;
    value.color.hue_noise_scale = 0.8300313f;
    value.color.hue_noise_speed = 0.00021200001f;
    value.color.hue_shift_amount = 0.398f;
    value.color.mapping_frequency = 1.0f;
    value.color.mapping_phase = 0.0f;
    value.color.palette_chroma = 0.62f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::LINEAR;
    value.color.phase_oscillation_depth = 0.0f;
    value.color.phase_oscillation_speed = 0.0f;
    value.color.opacity_high = 1.0f;
    value.color.opacity_low = 1.0f;
    value.projection.spin_rate = 0.020879198f;
    value.projection.wander = 0.0030917525f;
    value.projection.singularity_fade = 1.0f;
    value.source.lattice_cell_scale = 1.22925f;
    value.source.lattice_radius = 0.33298188f;
    value.source.lattice_shape_blend = 1.0f;
    value.source.lattice_softness = 0.1608203f;
    value.value.iso_level = 0.138f;
    value.value.iso_width = 0.22703418f;
    value.outer_warp.rotation_rate = 0.0f;
    value.outer_warp.scale_x = 1.0f;
    value.outer_warp.scale_y = 1.0f;
    value.outer_warp.shear = 0.0f;
    value.outer_warp.speed = 0.015625f;
    value.outer_warp.translation_x = 4.0f;
    value.outer_warp.translation_y = 4.0f;
    return value;
  }
  // clang-format on
  // End generated params.
};
