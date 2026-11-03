/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette_selection.h
 * @brief Workbench binding to the shared authored palette baker. */

#include "core/color/authored_palette.h"
#include "core/render/pullback/operators.h"

/** @brief Bakes the color operator's selected palette source. */
class ShaderPaletteResources : public AuthoredPaletteResources {
public:
  /** @brief Resolves the current source controls.
   * @param arena Persistent table allocation arena.
   * @param params Current color operator parameters.
   * @return Owned table, valid until the next source change. */
  const BakedPalette &
  prepare(Arena &arena,
          const Pullback::Interp::Op::GeneratedPaletteParams &params) {
    return AuthoredPaletteResources::prepare(arena, params,
                                             params.palette_chroma);
  }
};
