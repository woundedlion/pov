/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/lens.h"

/** @file authored_lens.h
 * @brief Parameterized glitch and kaleidoscope policies. */

namespace Pullback::Lens {

/** @brief Angular multiplication controls. */
struct GlitchParams {
  float polar_multiplier = 2.0f;   ///< Polar-angle multiplier.
  float azimuth_multiplier = 3.0f; ///< Azimuth-angle multiplier.
  /// Polar and azimuth angle controls.
  static constexpr auto FIELDS = std::array{
      Field<GlitchParams>{"polar-multiplier", &GlitchParams::polar_multiplier,
                          "Polar Multiplier", -12.0f, 12.0f, FieldCurve::LERP},
      Field<GlitchParams>{
          "azimuth-multiplier", &GlitchParams::azimuth_multiplier,
          "Azimuth Multiplier", -12.0f, 12.0f, FieldCurve::LERP}};
};

/** @brief Authored kaleidoscope symmetry. */
struct KaleidoscopeParams {
  uint8_t symmetry = 0; ///< Shared kaleidoscope symmetry index.
  /// The kaleidoscope has no scalar controls.
  static constexpr std::array<Field<KaleidoscopeParams>, 0> FIELDS{};
  /// Catalog ids in symmetry-index order.
  static constexpr const char *SYMMETRY_IDS[] = {
      "azimuthal",        "tetrahedral",      "octahedral",
      "dodecahedral",     "triangular-prism", "square-prism",
      "pentagonal-prism", "hexagonal-prism",  "octagonal-prism"};
  /// Symmetry selector.
  static constexpr auto TOPOLOGY = std::array{TopologyField<KaleidoscopeParams>{
      "symmetry", &KaleidoscopeParams::symmetry, SYMMETRY_IDS, 0}};
};

/** @brief Folds a direction into an authored symmetry chamber.
 * @param input Unit direction.
 * @param symmetry Shared kaleidoscope symmetry index.
 * @return Folded direction. */
inline math::Vector configured_kaleidoscope(const math::Vector &input,
                                            uint8_t symmetry) {
  switch (symmetry) {
  case 1:
    return TetrahedralKaleidoscope::apply(input, 0);
  case 2:
    return OctahedralKaleidoscope::apply(input, 0);
  case 3:
    return lenses::dodecahedral_kaleidoscope_lens(input);
  case 4:
    return TriangularPrismKaleidoscope::apply(input, 0);
  case 5:
    return SquarePrismKaleidoscope::apply(input, 0);
  case 6:
    return PentagonalPrismKaleidoscope::apply(input, 0);
  case 7:
    return HexagonalPrismKaleidoscope::apply(input, 0);
  case 8:
    return OctagonalPrismKaleidoscope::apply(input, 0);
  default:
    return lenses::kaleidoscope_lens(input);
  }
}

/** @brief Glitch lens with authored angular multipliers.
 * @tparam State Keyed lens provider. */
template <typename State> struct ConfiguredGlitch : ApproximationDefaults {
  /** @brief Applies the authored angular multipliers.
   * @param input Unit direction.
   * @param frame Bound frame.
   * @return Folded direction. */
  static math::Vector apply(const math::Vector &input,
                            const typename State::FrameState &frame) {
    const auto &params = State::params(frame);
    if (params.polar_multiplier == 2.0f && params.azimuth_multiplier == 3.0f)
      return lenses::glitch_lens(input);
    return math::Vector::from_spherical(
        atan2f(input.z, input.x) * params.azimuth_multiplier,
        acosf(hs::clamp(input.y, -1.0f, 1.0f)) * params.polar_multiplier);
  }
};

/** @brief Kaleidoscope lens with authored symmetry.
 * @tparam State Keyed lens provider. */
template <typename State>
struct ConfiguredKaleidoscope : ApproximationDefaults {
  static constexpr bool STROBE_COLUMNS =
      true; ///< Black column strobes enabled.
  /** @brief Folds the direction into the selected chamber.
   * @param input Unit direction.
   * @param frame Bound frame.
   * @return Folded direction. */
  static math::Vector apply(const math::Vector &input,
                            const typename State::FrameState &frame) {
    return configured_kaleidoscope(input, State::params(frame).symmetry);
  }
};

} // namespace Pullback::Lens
