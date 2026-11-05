/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file authored_source.h
 * @brief Source policies retaining authored topology controls. */

#include "core/render/pullback/source.h"

namespace Pullback::Source {

/** @brief Projected noise controls with an authored octave basis. */
struct AuthoredProjectedNoiseSourceParams : ProjectedNoiseSourceParams {
  uint8_t basis = 0; ///< Noise-basis catalog index.
  /// Catalog values of the basis selector.
  static constexpr auto &BASIS_IDS = math::NOISE_BASIS_IDS;
  /// Authored topology selectors.
  static constexpr auto TOPOLOGY = std::array{
      TopologyField<AuthoredProjectedNoiseSourceParams>{
          "basis", &AuthoredProjectedNoiseSourceParams::basis, BASIS_IDS, 0},
  };
};
/** @brief Projected noise with a runtime-selected octave basis.
 * @tparam State Keyed source provider. */
template <typename State>
struct ConfiguredProjectedNoise
    : ProjectedNoise<State, math::NoiseBasis::SIMPLEX> {
  using Base =
      ProjectedNoise<State,
                     math::NoiseBasis::SIMPLEX>; ///< Static policy contract.
  using FrameState = typename Base::FrameState;  ///< Bound frame.
  using Prepared = typename Base::Prepared;      ///< Noise-loop offset.
  /** @brief Prepares and validates the selected noise basis.
   * @param frame Current frame.
   * @return Noise-loop offset. */
  static Prepared prepare(const FrameState &frame) {
    HS_CHECK(State::params(frame).basis <= 2, "invalid projected-noise basis");
    return Base::prepare(frame);
  }
  /** @brief Samples the authored noise basis.
   * @param input Projected coordinate.
   * @param frame Current frame.
   * @param prepared Noise-loop offset.
   * @return Signed source value. */
  static float sample(const PlaneSample &input, const FrameState &frame,
                      const Prepared &prepared) {
    if (fabsf(input.coords.re) > 0x1p20f || fabsf(input.coords.im) > 0x1p20f)
      return 0.0f;
    return Pattern::noise_contour(
        State::noise(frame),
        static_cast<math::NoiseBasis>(State::params(frame).basis),
        math::noise_projected_coordinate(input.coords,
                                         State::noise_scale(frame), prepared),
        State::noise_contrast(frame));
  }
};

} // namespace Pullback::Source
