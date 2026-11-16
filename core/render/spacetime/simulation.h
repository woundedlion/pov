/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file simulation.h
 * @brief Canvas adapters for sweep-time simulation. */

#include "core/render/spacetime/surface.h"
#include "core/render/canvas.h"

namespace Spacetime {

/**
 * @brief Advances a simulation and shades it in sweep slabs, each slab after
 *        the substeps that end its share of the frame.
 * @tparam SUBSTEPS Substeps per frame.
 * @tparam SLABS Slab count, a divisor of `SUBSTEPS`.
 * @param canvas Target canvas; each slab shades under a TimeWindowScope.
 * @param advance Callable (int substeps) advancing the simulation.
 * @param shade Callable () shading the current state into the open window.
 * @details With temporal sampling disabled, runs every substep and then one
 * frame-wide shade.
 */
template <int SUBSTEPS, int SLABS, typename Advance, typename Shade>
HS_HOT_INLINE void simulate_in_slabs(Canvas &canvas, Advance &&advance,
                                     Shade &&shade) {
  static_assert(SLABS > 0 && SUBSTEPS % SLABS == 0);
  if (flat_surface()) {
    advance(SUBSTEPS);
    shade();
    return;
  }
  for (int slab = 0; slab < SLABS; ++slab) {
    advance(SUBSTEPS / SLABS);
    const TimeSlab columns = time_slab(slab, SLABS, canvas.width());
    TimeWindowScope window(canvas, columns.begin, columns.end);
    shade();
  }
}

} // namespace Spacetime
