/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/pattern/sources.h"

/** @file hex_morph_table.h
 * @brief Triangular interpolation of fixed hex-cell fields during preset morphs.
 */
namespace Pullback::Source {

/** @brief Periodic triangular samples in coordinates normal to hex cell faces. */
struct HexMorphTable {
  /// Intervals per lattice coordinate.
  static constexpr int STEPS = 32;
  /// Quantized field samples indexed by tile and lattice coordinates.
  uint16_t values[3][STEPS + 1][STEPS + 1]{};
  /** @var HexMorphTable::softness
   * @brief Cell softness of the resident profile. */
  /** @var HexMorphTable::shape
   * @brief Cell shape blend of the resident profile. */
  /// Cell radius of the resident profile, or negative before baking.
  float radius = -1.0f, softness = 0.0f, shape = 0.0f;

  /** @brief Tests whether the shape controls match the resident profile.
   * @param p Requested tessellation profile.
   * @return Whether the pre-baked field can serve the profile. */
  bool matches(const Pattern::TessellationCellParams &p) const {
    return p.cell_radius == radius && p.cell_softness == softness &&
           p.cell_shape == shape;
  }

  /** @brief Bakes periodic field samples for the authored cell shape.
   * @param p Tessellation profile to sample. */
  HS_SIZE_FLASH_MEMBER void bake(const Pattern::TessellationCellParams &p) {
    radius = p.cell_radius;
    softness = p.cell_softness;
    shape = p.cell_shape;
    auto sample_params = p;
    sample_params.cell_scale = 1.0f;
    for (int tile = 0; tile < 3; ++tile)
      for (int i = 0; i <= STEPS; ++i)
        for (int j = 0; j <= STEPS; ++j) {
          const float a = tile + static_cast<float>(i) / STEPS;
          const float b = static_cast<float>(j) / STEPS;
          const Pattern::PlanePosition position{
              {0.5f * (a + b), 0.5f * Pattern::INV_SQRT_3 * (a - b)}};
          const float field = Pattern::hexagonal_cells(position, sample_params);
          values[tile][i][j] =
              static_cast<uint16_t>((1.0f + field) * 32767.5f + 0.5f);
        }
  }

  /** @brief Quantizes a lattice coordinate, reducing large coordinates periodically.
   * @param x Lattice coordinate in periods.
   * @return Signed fixed-point table coordinate. */
  static HS_HOT_INLINE int fixed_coordinate(float x) {
    constexpr float factor = 256.0f * STEPS;
    if (fabsf(x) >= 2147483648.0f / factor)
      x = fmodf(x, 3.0f);
    return static_cast<int>(floorf(x * factor));
  }

  /** @brief Interpolates the wrapped field in the triangular sample grid.
   * @param coord Projected coordinate.
   * @param scale Cell frequency in the projected plane.
   * @return Interpolated signed field value. */
  HS_HOT_INLINE float sample(math::Complex coord, float scale) const {
    const float y = Pattern::SQRT_3 * coord.im;
    const int a = fixed_coordinate((coord.re + y) * scale);
    const int b = fixed_coordinate((coord.re - y) * scale);
    const unsigned tile = (2147483646u + (a >> 13) - (b >> 13)) % 3u;
    const unsigned i = (a >> 8) & (STEPS - 1);
    const unsigned j = (b >> 8) & (STEPS - 1);
    const int da = a & 255, db = b & 255;
    const auto &grid = values[tile];
    const int x = grid[i + 1][j], edge_y = grid[i][j + 1];
    int value;
    if (da + db <= 256) {
      const int base = grid[i][j];
      value = 256 * base + da * (x - base) + db * (edge_y - base);
    } else {
      const int base = grid[i + 1][j + 1];
      value =
          256 * base + (256 - da) * (edge_y - base) + (256 - db) * (x - base);
    }
    return static_cast<float>(value) / 256.0f;
  }
};

/// Immutable source tables borrowed by one transition, in departure order.
struct HexMorphPair {
  /// Table at the outgoing endpoint.
  const HexMorphTable *from = nullptr;
  /// Table at the incoming endpoint.
  const HexMorphTable *to = nullptr;
};

/// Two authored cell profiles baked at effect initialization.
struct HexMorphTables {
  /// Storage for the two supported profiles.
  HexMorphTable tables[2];
  /// Tables selected for the active transition.
  HexMorphPair active{};

  /** @brief Fills the owned tables at effect initialization.
   * @param a First authored profile.
   * @param b Second authored profile. */
  void prewarm(const Pattern::TessellationCellParams &a,
               const Pattern::TessellationCellParams &b) {
    tables[0].bake(a);
    tables[1].bake(b);
  }

  /// Returns an empty pair when either endpoint needs the analytic source.
  /// @param a Outgoing tessellation profile.
  /// @param b Incoming tessellation profile.
  /// @return Borrowed tables in endpoint order, or an empty pair.
  HS_SIZE_FLASH_MEMBER HexMorphPair
  select(const Pattern::TessellationCellParams &a,
         const Pattern::TessellationCellParams &b) const {
    const auto find = [&](const auto &p) -> const HexMorphTable * {
      for (const auto &table : tables)
        if (table.matches(p))
          return &table;
      return nullptr;
    };
    const auto *from = find(a);
    const auto *to = find(b);
    return from && to ? HexMorphPair{from, to} : HexMorphPair{};
  }
};

} // namespace Pullback::Source
