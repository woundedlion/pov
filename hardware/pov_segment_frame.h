/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file pov_segment_frame.h
 * @brief Preserves the inactive arm half after the pre-clear hook selects the
 *        segment clip for the newly acquired draw buffer.
 */
#pragma once

#include "core/render/canvas.h"

#include <cstring>

namespace pov {

/**
 * @brief Retains the previous image outside this frame's segment half.
 * @pre The display clip is one non-wrapping horizontal half, [0, width/2) or
 * [width/2, width), as returned by segment_clip() for an even canvas width.
 * Only rows within the display clip are preserved.
 */
inline void preserve_segment_half(Canvas &canvas) {
  HS_PROFILE(pov_preserve_half);
  const ClipRegion &clip = canvas.clip();
  const int X0 = clip.x_start == 0 ? clip.x_end : 0;
  const int X1 = clip.x_start == 0 ? canvas.width() : clip.x_start;
  Pixel *const dest = canvas.data();
  const Pixel *const source = canvas.prev_data();
  for (int y = clip.y_start; y < clip.y_end; ++y) {
    const int OFFSET = y * canvas.width() + X0;
    std::memcpy(dest + OFFSET, source + OFFSET, (X1 - X0) * sizeof(Pixel));
  }
}

} // namespace pov
