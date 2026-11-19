/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file planar_source.h
 * @brief Planar source selector and kind families. */

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string_view>

#include "core/render/pullback/fields.h"

namespace Pullback::Source {

/** @brief Planar source kinds, in catalog order. */
enum class PlanarKind : uint8_t {
  GRID,                          ///< Coupled grid.
  TWIN_WAVE,                     ///< Two-wave interference.
  BANDS,                         ///< Parallel bands.
  SPIRAL,                        ///< Stationary spiral.
  TESSELLATION_TRIANGULAR_CELLS, ///< Triangular cell primitives.
  TESSELLATION_TRIANGULAR_EDGES, ///< Triangular edge bands.
  TESSELLATION_SQUARE_CELLS,     ///< Square cell primitives.
  TESSELLATION_SQUARE_EDGES,     ///< Square edge bands.
  TESSELLATION_HEXAGONAL_CELLS,  ///< Hexagonal cell primitives.
  TESSELLATION_HEXAGONAL_EDGES,  ///< Hexagonal edge bands.
  PROJECTED_NOISE_SIMPLEX,       ///< Projected simplex noise.
  PROJECTED_NOISE_FBM3,          ///< Projected fractional Brownian noise.
  PROJECTED_NOISE_RIDGED3,       ///< Projected ridged noise.
};

/** @brief Catalog ids in PlanarKind order. */
inline constexpr const char *PLANAR_KIND_IDS[] = {
    "grid",
    "twin-wave",
    "bands",
    "spiral",
    "tessellation-triangular-cells",
    "tessellation-triangular-edges",
    "tessellation-square-cells",
    "tessellation-square-edges",
    "tessellation-hexagonal-cells",
    "tessellation-hexagonal-edges",
    "projected-noise-simplex",
    "projected-noise-fbm3",
    "projected-noise-ridged3"};

/** @brief Number of catalog source kinds. */
inline constexpr size_t PLANAR_KIND_COUNT = std::size(PLANAR_KIND_IDS);
static_assert(PLANAR_KIND_COUNT ==
              static_cast<size_t>(PlanarKind::PROJECTED_NOISE_RIDGED3) + 1);
static_assert(PLANAR_KIND_COUNT <= 16);

/** @brief Planar parameter families, in catalog prefix order. */
enum class PlanarFamily : uint8_t {
  GRID,            ///< Grid controls.
  TWIN_WAVE,       ///< Twin-wave controls.
  BANDS,           ///< Bands controls.
  SPIRAL,          ///< Spiral controls.
  TESSELLATION,    ///< Shared controls across tessellation variants.
  PROJECTED_NOISE, ///< Shared controls across noise bases.
};

/** @brief Catalog field prefixes in PlanarFamily order. */
inline constexpr std::string_view PLANAR_FAMILY_IDS[] = {
    "grid", "twin-wave", "bands", "spiral", "tessellation", "projected-noise"};

/** @brief Catalog ids of fields shared across all planar source kinds. */
inline constexpr std::string_view PLANAR_COMMON_IDS[] = {
    "source", "edge-width", "weight-mode", "coverage-mode"};

/** @brief Whether family prefixes and common ids are unambiguous.
 * @return True when none of the labels prefixes another. */
consteval bool planar_prefixes_unambiguous() {
  std::array<std::string_view,
             std::size(PLANAR_FAMILY_IDS) + std::size(PLANAR_COMMON_IDS)>
      labels{};
  size_t next = 0;
  for (const auto family : PLANAR_FAMILY_IDS)
    labels[next++] = family;
  for (const auto common : PLANAR_COMMON_IDS)
    labels[next++] = common;
  for (size_t first = 0; first < labels.size(); ++first)
    for (size_t second = first + 1; second < labels.size(); ++second)
      if (labels[first].starts_with(labels[second]) ||
          labels[second].starts_with(labels[first]))
        return false;
  return true;
}
static_assert(planar_prefixes_unambiguous(), "ambiguous planar source prefix");

/** @brief Planar kinds using tessellation cell fields. */
inline constexpr uint16_t PLANAR_TESSELLATION_CELLS =
    live_values(PlanarKind::TESSELLATION_TRIANGULAR_CELLS,
                PlanarKind::TESSELLATION_SQUARE_CELLS,
                PlanarKind::TESSELLATION_HEXAGONAL_CELLS);

/** @brief Planar kinds using tessellation edge fields. */
inline constexpr uint16_t PLANAR_TESSELLATION_EDGES =
    live_values(PlanarKind::TESSELLATION_TRIANGULAR_EDGES,
                PlanarKind::TESSELLATION_SQUARE_EDGES,
                PlanarKind::TESSELLATION_HEXAGONAL_EDGES);

/** @brief One live-kind mask per planar family. */
inline constexpr uint16_t PLANAR_FAMILY_MASKS[] = {
    live_values(PlanarKind::GRID),
    live_values(PlanarKind::TWIN_WAVE),
    live_values(PlanarKind::BANDS),
    live_values(PlanarKind::SPIRAL),
    PLANAR_TESSELLATION_CELLS | PLANAR_TESSELLATION_EDGES,
    live_values(PlanarKind::PROJECTED_NOISE_SIMPLEX,
                PlanarKind::PROJECTED_NOISE_FBM3,
                PlanarKind::PROJECTED_NOISE_RIDGED3)};
static_assert(std::size(PLANAR_FAMILY_IDS) == std::size(PLANAR_FAMILY_MASKS));

/** @brief Family of a valid catalog kind.
 * @param kind Selected source kind.
 * @return Its parameter family. */
constexpr PlanarFamily planar_family(PlanarKind kind) {
  const auto index = static_cast<uint8_t>(kind);
  if (index <= static_cast<uint8_t>(PlanarKind::SPIRAL))
    return static_cast<PlanarFamily>(index);
  if (index < static_cast<uint8_t>(PlanarKind::PROJECTED_NOISE_SIMPLEX))
    return PlanarFamily::TESSELLATION;
  return PlanarFamily::PROJECTED_NOISE;
}

/** @brief Source kinds belonging to the family named by a field id.
 * @param id Family-qualified field id, or a common field id.
 * @return The family kind mask, or zero for a common field. */
constexpr uint16_t planar_field_sources(std::string_view id) {
  size_t length = 0;
  uint16_t sources = 0;
  for (size_t index = 0; index < std::size(PLANAR_FAMILY_IDS); ++index) {
    const auto prefix = PLANAR_FAMILY_IDS[index];
    if (prefix.size() > length && id.size() > prefix.size() &&
        id.starts_with(prefix) && id[prefix.size()] == '-') {
      length = prefix.size();
      sources = PLANAR_FAMILY_MASKS[index];
    }
  }
  return sources;
}

} // namespace Pullback::Source
