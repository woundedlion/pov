/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using AlienCoreParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::MirrorParams,
                     Pullback::NoWarpParams, Pullback::NoLensParams,
                     Pullback::EdgeValueParams>;
using AlienCoreSpec =
    Pullback::Spec<Pullback::ProjectionKind::GNOMONIC_FOLDED,
                   Pullback::Lens::Glitch, Pullback::TransferKind::NONE,
                   Pullback::ProjectionCoverageMode::EDGE_FADE>;

/**
 * @brief A mirrored grid folded by the glitch lens.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class AlienCore
    : public Pullback::ComposedEffect<W, H, AlienCore<W, H>, AlienCoreParams,
                                      AlienCoreSpec, PaletteHarmony::TRIADIC,
                                      Pullback::HueMode::NOISE,
                                      Pullback::Color::BrightnessEnvelope::NONE,
                                      /*AnimatedProjection=*/false> {

public:
  using Params = AlienCoreParams;
  static constexpr std::string_view EFFECT_ID = "alien-core";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "6e2bc1705f48b16d77b0a5c086dc5c187e485b2c65d753f66375114e9ef64013";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "8b76fda1895b0572c30c79e00baa135e3bbcc151c2c051f3c1ab0773be83c366";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"folded-glitch"};
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
    return value;
  }
};

#include "core/control/registry.h"
