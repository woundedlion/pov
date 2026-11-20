/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spherical_tiling_source.h
 * @brief Shader Tiling parameter schema backed by the solids registry. */

#include "core/mesh/solids.h"
#include "core/render/pullback/fields.h"

namespace Pullback::Source {

/** @brief Scalar output of a spherical tiling. */
enum class TilingOutput : uint8_t { EDGES, FACES };

namespace detail {

/// Display names in flattened solid registry order.
inline constexpr auto TILING_SOLID_NAMES = [] {
  constexpr size_t LONGEST = [] {
    size_t length = 0;
    for (auto registry : Solids::all_registries())
      for (const auto &entry : registry)
        length = std::max(length, std::string_view(entry.name).size());
    return length;
  }();
  std::array<std::array<char, 2 * LONGEST + 1>, Solids::NUM_ENTRIES> names{};
  size_t index = 0;
  for (auto registry : Solids::all_registries()) {
    for (const auto &entry : registry) {
      auto &name = names[index++];
      size_t position = 0;
      for (char c : std::string_view(entry.name)) {
        if (c >= 'A' && c <= 'Z') {
          name[position++] = '-';
          c = static_cast<char>(c - 'A' + 'a');
        }
        name[position++] = c == '_' ? '-' : c;
      }
    }
  }
  return names;
}();

/// Stable solid ids in flattened registry order.
inline constexpr auto TILING_SOLID_IDS = [] {
  struct Values {
    const char *ids[Solids::NUM_ENTRIES];
  };
  Values values{};
  for (size_t index = 0; index < TILING_SOLID_NAMES.size(); ++index)
    values.ids[index] = TILING_SOLID_NAMES[index].data();
  return values;
}();

} // namespace detail

/** @brief Registry solid and scalar controls for Tiling. */
struct SphericalTilingSourceParams {
  static constexpr const char *OPERATOR_ID =
      "sample.spherical-tiling.v1"; ///< Catalog identity.
  static constexpr uint8_t DEFAULT_SOLID =
      Solids::SEED_TRUNCATED_ICOSAHEDRON; ///< Default registry entry.

  float line_thickness = 1.0f; ///< Full edge width in degrees.
  float line_softness =
      0.5f; ///< Fringe width in degrees for Edges or the deepest face.
  uint8_t solid = DEFAULT_SOLID; ///< Registry entry ordinal.
  uint8_t output = static_cast<uint8_t>(TilingOutput::EDGES); ///< Scalar mode.

  /// Scalar output choices.
  static constexpr const char *OUTPUT_IDS[] = {"edges", "faces"};
  static_assert(std::size(OUTPUT_IDS) ==
                static_cast<size_t>(TilingOutput::FACES) + 1);

  /// Width and softness controls.
  static constexpr auto FIELDS = std::array{
      Field<SphericalTilingSourceParams>{
          "line-thickness", &SphericalTilingSourceParams::line_thickness,
          "Edge Thickness (deg)", 0.0f, 10.0f, FieldCurve::LERP},
      Field<SphericalTilingSourceParams>{
          "line-softness", &SphericalTilingSourceParams::line_softness,
          "Edge Softness (deg)", 0.05f, 10.0f, FieldCurve::LOG_POSITIVE}};

  /// Registry solid and scalar output selectors.
  static constexpr auto TOPOLOGY =
      std::array{TopologyField<SphericalTilingSourceParams>{
                     "solid", &SphericalTilingSourceParams::solid,
                     detail::TILING_SOLID_IDS.ids, DEFAULT_SOLID},
                 TopologyField<SphericalTilingSourceParams>{
                     "output", &SphericalTilingSourceParams::output, OUTPUT_IDS,
                     static_cast<uint8_t>(TilingOutput::EDGES)}};
};

} // namespace Pullback::Source
