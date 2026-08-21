/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/composed_effect.h"

namespace hs_test {
namespace facet_grid_tests {
struct FacetGridWhiteBox;
} // namespace facet_grid_tests
} // namespace hs_test

using FacetGridParams =
    Pullback::Params<Pullback::GridSourceParams, Pullback::NoWarpParams,
                     Pullback::MirrorParams>;
using FacetGridSpec =
    Pullback::Spec<Pullback::ProjectionKind::STEREOGRAPHIC,
                   Pullback::Lens::DodecahedralKaleidoscope,
                   Pullback::TransferKind::NONE,
                   Pullback::CoverageKind::PROJECTION_SQUARED>;

/**
 * @brief Mirrored grids folded through a dodecahedral stereographic lens.
 * @details Supplies the render pipeline and preset bank; Pullback::ComposedEffect
 * supplies parameter registration, preset choreography and the palette,
 * camera-walk and noise clocks. The kaleidoscopic fold lives in the surface
 * lens, so the surface itself is an identity and the mirror tiling happens in
 * the planar warp.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class FacetGrid : public Pullback::ComposedEffect<
                      W, H, FacetGrid<W, H>, FacetGridParams, FacetGridSpec,
                      PaletteHarmony::ANALOGOUS, Pullback::HueMode::NOISE,
                      Pullback::Color::BrightnessEnvelope::NONE> {
  friend struct ::hs_test::facet_grid_tests::FacetGridWhiteBox;

public:
  using Params = FacetGridParams;
  /// Registry identity, stable across class renames.
  static constexpr std::string_view EFFECT_ID = "facet-grid";
  /// SHA-256 of the canonicalized descriptor of
  /// `patterns/facet_grid.shader.json`; the browser editor matches it to
  /// recognize an imported document as this composed effect.
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "5c161697e1e62c963d9063a933ce016bf1b94a7a9444ca3f4178379b1bd293cc";
  /// SHA-256 of that document's canonicalized preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "d9f9c3ebee94b1144d9f66512a76168f9ed90a30f556a474d2d7f0ee763da5e8";
  /// Immutable preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 4> PRESET_IDS{
      "coupled-grid", "direct-grid", "double-map", "stretched-grid"};
  /// Bumped whenever the `Params` layout changes, rejecting stale snapshots.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 3;
  // The 19 s Phantasm slot is 304 frames at 16 fps and the effect is rebuilt
  // each visit: 4 dwells plus 3 crossfades must fit, or the later presets never
  // render.
  static constexpr Segue::Lerp PRESET_SEGUE{16, ease_in_out_sin};
  /// Frames a preset holds before the runtime begins the next transition.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 64;
  static constexpr bool ANIMATED_PROJECTION = true;

  /// Params the effect starts on, and the base every preset varies from.
  static constexpr Params initial_params() {
    Params value;
    value.source = {.pattern_freq = 2.82629991f,
                    .speed = 0.0f,
                    .complexity = 0.513f,
                    .pattern_mix = 0.0f,
                    .secondary_rate = 0.8f,
                    .angle_rate = 0.0269999988f};
    value.projection.pole_fade = 3.432f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.inner_warp.speed = 0.00013f;
    value.inner_warp.cell_y = 0.997703135f;
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.hue_shift_amount = 0.366f;
    value.color.hue_noise_scale = 1.47215629f;
    return value;
  }

  /**
   * @brief Params for the preset at @p index in PRESET_IDS.
   * @details `coupled-grid` is the initial preset; the other three raise the grid
   * complexity and fully mix in the secondary pattern, then vary its frequency,
   * projection wander, palette mapping frequency and — for `stretched-grid` —
   * the mirror-tile cell and rotation.
   */
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 0)
      return value;
    value.source.complexity = 3.0f;
    value.source.pattern_mix = 1.0f;
    if (index == 2) {
      value.source.pattern_freq = 3.9407f;
      value.projection.wander = 0.165f;
      value.color.mapping_frequency = 2.0f;
    } else if (index == 3) {
      value.source.pattern_freq = 2.9059f;
      value.projection.wander = 0.165f;
      value.inner_warp.speed = 0.0027299998f;
      value.inner_warp.rotation = 3.455752f;
      value.inner_warp.cell_x = 0.22321875f;
      value.inner_warp.cell_y = 5.085703f;
      value.color.mapping_frequency = 1.558f;
    }
    return value;
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(FacetGrid)
