/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using AlienOceanParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::MirrorParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::EdgeValueParams>;
using AlienOceanSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::Kaleidoscope, Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::EDGE_FADE>;

/**
 * @brief A broad folded grid with slow mirrored drift.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AlienOcean
    : public Pullback::ComposedEffect<W, H, AlienOcean<W, H>, AlienOceanParams,
                                      AlienOceanSpec, PaletteHarmony::TRIADIC,
                                      Pullback::HueMode::NOISE,
                                      Pullback::Color::BrightnessEnvelope::NONE,
                                      /*AnimatedProjection=*/false> {

public:
  using Params = AlienOceanParams;
  static constexpr std::string_view EFFECT_ID = "alien-ocean";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "d1598bca5ca501586d6bd14a2767c7d8bfbcef4f447f9713543a4d1589360d16";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "0e1c9f3208591e7c5946c21a5c30a7d587c88b0f91c54888f2b4e16ff6679703";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"folded-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 3.565f,
                    .speed = 0.235f,
                    .complexity = 0.0f,
                    .pattern_mix = 1.0f,
                    .secondary_rate = 1.0f,
                    .angle_rate = 0.0f};
    value.projection.singularity_fade = 1.4f;
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
};

#include "core/control/registry.h"
REGISTER_EFFECT(AlienOcean)
