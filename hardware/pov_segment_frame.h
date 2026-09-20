/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/canvas.h"

namespace pov {

/** @brief Retains the previous image outside this frame's segment half. */
inline void preserve_segment_half(Canvas &canvas) {
  HS_PROFILE(pov_preserve_half);
  const ClipRegion &clip = canvas.clip();
  const int X0 = clip.x_start == 0 ? clip.x_end : 0;
  const int X1 = clip.x_start == 0 ? canvas.width() : clip.x_start;
  Pixel *const dest = canvas.data();
  const Pixel *const source = canvas.prev_data();
  for (int y = clip.y_start; y < clip.y_end; ++y) {
    const int OFFSET = y * canvas.width() + X0;
    memcpy(dest + OFFSET, source + OFFSET, (X1 - X0) * sizeof(Pixel));
  }
}

} // namespace pov
