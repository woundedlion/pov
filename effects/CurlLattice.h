/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "core/render/pullback/look.h"

namespace hs_test {
namespace curl_lattice_tests {
struct CurlLatticeWhiteBox;
} // namespace curl_lattice_tests
} // namespace hs_test

/**
 * @brief Fixed folded-sinusoidal lattice displaced by sphere-space curl noise.
 * @details Supplies the render pipeline and preset bank; Looks::Composed
 * supplies parameter registration, preset choreography and the palette,
 * camera-walk and noise clocks. The lattice source is read through a folded
 * sinusoidal projection, so the surface stage carries the curl displacement and
 * the warp stage is an identity.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
using CurlLatticeParams =
    Looks::Params<Looks::LatticeSourceParams, Looks::NoWarpParams,
                  Looks::NoWarpParams, Looks::NoLensParams,
                  Looks::LinearValueParams, Looks::SurfaceNoiseParams>;
using CurlLatticeSpec =
    Looks::Spec<Looks::ProjectionKind::FOLDED_SINUSOIDAL,
                Pullback::Lens::Identity, Looks::TransferKind::LINEAR,
                Looks::CoverageKind::PROJECTION>;

template <int W, int H>
class CurlLattice
    : public Looks::Composed<W, H, CurlLattice<W, H>, CurlLatticeParams,
                             CurlLatticeSpec, PaletteHarmony::TRIADIC,
                             Looks::HueMode::NOISE,
                             Pullback::Color::BrightnessEnvelope::CUP> {
  friend struct ::hs_test::curl_lattice_tests::CurlLatticeWhiteBox;

public:
  using Params = CurlLatticeParams;
  /// Registry identity, stable across class renames.
  static constexpr std::string_view EFFECT_ID = "curl-lattice";
  /// SHA-256 of the canonicalized descriptor of
  /// `patterns/curl_lattice.shader.json`; the browser editor matches it to
  /// recognize an imported document as this fixed effect.
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "504e5dc75bbd656d36b94b2752c0b6e3166ce80221f9b63d63127967323c96f8";
  /// SHA-256 of that document's canonicalized preset bank.
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "9b4b0152cd159b90bd3663505f53c8b6bdbcce846a37594135eff1e46aa8eaa3";
  /// Immutable preset identities, indexed by preset number.
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"open-curl",
                                                              "dense-curl"};
  /// Bumped whenever the `Params` layout changes, rejecting stale snapshots.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 4;
  /// Frames a preset holds before the runtime begins the next transition.
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr bool ANIMATED_PROJECTION = true;
  static constexpr bool USES_CENTRAL_MERIDIAN = true;
  static constexpr int32_t SURFACE_NOISE_SEED = 1337;

  /// Params the effect starts on, and the base every preset varies from.
  static constexpr Params initial_params() {
    Params value;
    value.source = {0.710265636f, 1.0f, 0.455532223f, 0.290762514f};
    value.projection.pole_fade = 20.0f;
    value.projection.wander = 1.0f;
    value.projection.camera_wander = 1.0f;
    value.surface = {1.78815627f, 0.0759999976f, 0.0f};
    value.color.palette_chroma = 1.0f;
    value.color.palette_mapping = Pullback::Color::PaletteMapping::CUP;
    value.color.mapping_phase = -0.165999994f;
    value.color.hue_shift_amount = 0.268000007f;
    value.color.hue_noise_scale = 2.0f;
    return value;
  }

  /**
   * @brief Params for the preset at @p index in PRESET_IDS.
   * @details `open-curl` is the initial look; `dense-curl` folds the lattice
   * through a shorter surface-noise wavelength.
   */
  static constexpr Params preset_params(size_t index) {
    Params value = initial_params();
    if (index == 1)
      value.surface.scale = 3.29720306f;
    return value;
  }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(CurlLattice)
