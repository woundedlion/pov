/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

using SignalWeaveParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::WaveShearParams,
                     Pullback::NoWarpParams>;
using SignalWeaveSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::Glitch, Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

template <int W, int H>
class SignalWeave
    : public Pullback::ComposedEffect<
          W, H, SignalWeave<W, H>, SignalWeaveParams, SignalWeaveSpec,
          PaletteHarmony::TRIADIC, Pullback::HueMode::NOISE,
          Pullback::Color::BrightnessEnvelope::NONE> {

public:
  using Params = SignalWeaveParams;
  static constexpr std::string_view EFFECT_ID = "signal-weave";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "f12ffa86e951c04f5ebd810f4bf85d7cc9430ca1aad0525e1d637ebb0586a065";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "f8904a8ff85c3ad0b371e6dcaeea35d234decf146e1d799a53e98a9a7e712cb7";
  static constexpr uint16_t INITIAL_PRESET_DWELL_FRAMES = 16;
  // The 7 s Phantasm slot is 112 frames at 16 fps and the effect is rebuilt each
  // visit: the initial hold plus 3 crossfade+dwell cycles must fit, or the later
  // presets never render.
  static constexpr Segue::Lerp PRESET_SEGUE{6, ease_in_out_sin};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 26;
  static constexpr std::array<std::string_view, 4> PRESET_IDS{
      "signal-weave", "signal-weave-2", "signal-weave-3", "signal-weave-4"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr bool ANIMATED_PROJECTION = true;

  static HS_HOT_FLASH_MEMBER Color4
  shade(const Vector &view, const typename SignalWeave::FrameState &frame) {
    return SignalWeave::RenderPipeline::shade(view, frame);
  }
  static constexpr Params initial_params() { return make_params(0); }
  static constexpr Params preset_params(size_t index) {
    return make_params(index);
  }

public:
  HS_COLD_MEMBER void after_composed_init() {
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

#include "core/control/registry.h"
REGISTER_EFFECT(SignalWeave)
