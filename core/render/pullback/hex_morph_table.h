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

/** @brief A triangular interpolation grid on a hex cell's 30-degree sector. */
struct HexMorphTable {
  static constexpr int STEPS = 60; ///< Intervals along each sector axis.
  uint16_t values[STEPS + 1][STEPS + 1]{};
  float radius = -1.0f;
  float softness = 0.0f;
  float shape = 0.0f;

  /// Scale is applied at sampling time. The other three controls define the grid.
  bool matches(const Pattern::TessellationCellParams &p) const {
    return p.cell_radius == radius && p.cell_softness == softness &&
           p.cell_shape == shape;
  }

  /// Samples the analytic hex primitive at each sector grid vertex.
  HS_SIZE_FLASH_MEMBER void bake(const Pattern::TessellationCellParams &p) {
    radius = p.cell_radius;
    softness = p.cell_softness;
    shape = p.cell_shape;
    for (int i = 0; i <= STEPS; ++i)
      for (int j = 0; j <= STEPS; ++j) {
        const float x = (0.5f / STEPS) * i;
        const float y = (0.5f * Pattern::INV_SQRT_3 / STEPS) * j;
        const float value =
            Pattern::cell_primitive<Pattern::TessellationKind::HEXAGONAL>(x, y,
                                                                          p);
        values[i][j] = static_cast<uint16_t>((1.0f + value) * 32767.5f + 0.5f);
      }
  }

  /// Interpolates the field encoded on [0, 65535] at a point on the working plane.
  HS_HOT_INLINE float sample(math::Complex coord, float scale) const {
    const float u = coord.re * scale;
    const float r = (coord.im * scale) * Pattern::INV_SQRT_3;
    float x = fabsf(u - nearbyintf(u));
    float row = fabsf(r - nearbyintf(r));
    if (x + 3.0f * row > 1.0f) {
      x = 0.5f - x;
      row = 0.5f - row;
    }
    const float y = Pattern::SQRT_3 * row;
    const float fold = fminf(0.5f * x - Pattern::HALF_SQRT_3 * y, 0.0f);
    const float s = (2.0f * STEPS) * (x - fold);
    const float t =
        (2.0f * Pattern::SQRT_3 * STEPS) * fabsf(y + Pattern::SQRT_3 * fold);
    const int i = std::min(static_cast<int>(s), STEPS - 1);
    const int j = std::min(static_cast<int>(t), STEPS - 1);
    const float a = values[i][j], b = values[i + 1][j + 1];
    const float ds = s - i, dt = t - j;
    if (ds >= dt) {
      const float mid = values[i + 1][j];
      return a + ds * (mid - a) + dt * (b - mid);
    }
    const float mid = values[i][j + 1];
    return a + dt * (mid - a) + ds * (b - mid);
  }
};

/// Immutable source tables borrowed by one transition, in departure order.
struct HexMorphPair {
  const HexMorphTable *from = nullptr;
  const HexMorphTable *to = nullptr;
};

/// Two authored cell profiles baked at effect initialization.
struct HexMorphTables {
  HexMorphTable tables[2];
  HexMorphPair active{};

  void prewarm(const Pattern::TessellationCellParams &a,
               const Pattern::TessellationCellParams &b) {
    tables[0].bake(a);
    tables[1].bake(b);
  }

  /// Returns an empty pair when either endpoint needs the analytic source.
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
