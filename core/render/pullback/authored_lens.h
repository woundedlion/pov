/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/lens.h"

/** @file authored_lens.h
 * @brief Glitch and kaleidoscope parameter families and their parameterized
 * policies. */

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

/** @brief Symmetry group of the kaleidoscope lens. */
enum class KaleidoscopeSymmetry : uint8_t {
  AZIMUTHAL = 0,        ///< Azimuthal wedge fold.
  TETRAHEDRAL = 1,      ///< Tetrahedral group.
  OCTAHEDRAL = 2,       ///< Octahedral group.
  DODECAHEDRAL = 3,     ///< Icosahedral/dodecahedral group.
  TRIANGULAR_PRISM = 4, ///< 3-fold prism group.
  SQUARE_PRISM = 5,     ///< 4-fold prism group.
  PENTAGONAL_PRISM = 6, ///< 5-fold prism group.
  HEXAGONAL_PRISM = 7,  ///< 6-fold prism group.
  OCTAGONAL_PRISM = 8   ///< 8-fold prism group.
};

/** @brief Evaluation of a fixed-symmetry kaleidoscope fold. */
enum class KaleidoscopeFold : uint8_t {
  DIRECT,                   ///< Folds any direction.
  TABLE,                    ///< Dodecahedral, by validated matrix lookup.
  FROM_DODECAHEDRAL_CHAMBER ///< Tetrahedral, from the dodecahedral chamber.
};

/** @brief Glitch family with the fixed default angular multipliers. */
struct FixedGlitchParams {};

/** @brief Kaleidoscope family with a compile-time symmetry and fold.
 * @tparam Symmetry Symmetry group.
 * @tparam Fold Fold evaluation. */
template <KaleidoscopeSymmetry Symmetry,
          KaleidoscopeFold Fold = KaleidoscopeFold::DIRECT>
  requires(Fold == KaleidoscopeFold::DIRECT ||
           (Fold == KaleidoscopeFold::TABLE &&
            Symmetry == KaleidoscopeSymmetry::DODECAHEDRAL) ||
           (Fold == KaleidoscopeFold::FROM_DODECAHEDRAL_CHAMBER &&
            Symmetry == KaleidoscopeSymmetry::TETRAHEDRAL))
struct FixedKaleidoscopeParams {
  static constexpr KaleidoscopeSymmetry SYMMETRY = Symmetry; ///< Symmetry.
  static constexpr KaleidoscopeFold FOLD = Fold;             ///< Fold.
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
  static_assert(std::size(SYMMETRY_IDS) ==
                static_cast<size_t>(KaleidoscopeSymmetry::OCTAGONAL_PRISM) + 1);
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
  switch (static_cast<KaleidoscopeSymmetry>(symmetry)) {
  case KaleidoscopeSymmetry::TETRAHEDRAL:
    return TetrahedralKaleidoscope::apply(input, 0);
  case KaleidoscopeSymmetry::OCTAHEDRAL:
    return OctahedralKaleidoscope::apply(input, 0);
  case KaleidoscopeSymmetry::DODECAHEDRAL:
    return lenses::dodecahedral_kaleidoscope_lens(input);
  case KaleidoscopeSymmetry::TRIANGULAR_PRISM:
    return TriangularPrismKaleidoscope::apply(input, 0);
  case KaleidoscopeSymmetry::SQUARE_PRISM:
    return SquarePrismKaleidoscope::apply(input, 0);
  case KaleidoscopeSymmetry::PENTAGONAL_PRISM:
    return PentagonalPrismKaleidoscope::apply(input, 0);
  case KaleidoscopeSymmetry::HEXAGONAL_PRISM:
    return HexagonalPrismKaleidoscope::apply(input, 0);
  case KaleidoscopeSymmetry::OCTAGONAL_PRISM:
    return OctagonalPrismKaleidoscope::apply(input, 0);
  case KaleidoscopeSymmetry::AZIMUTHAL:
    return lenses::kaleidoscope_lens(input);
  }
  HS_CHECK(false, "configured_kaleidoscope: invalid symmetry");
  return input;
}

/** @brief Folds a direction by authored angular multipliers.
 * @param input Unit direction.
 * @param params Angular multipliers.
 * @return Folded direction. */
inline math::Vector configured_glitch(const math::Vector &input,
                                      const GlitchParams &params) {
  if (params.polar_multiplier == 2.0f && params.azimuth_multiplier == 3.0f)
    return lenses::glitch_lens(input);
  return math::Vector::from_spherical(
      atan2f(input.z, input.x) * params.azimuth_multiplier,
      acosf(hs::clamp(input.y, -1.0f, 1.0f)) * params.polar_multiplier);
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
    return configured_glitch(input, State::params(frame));
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
