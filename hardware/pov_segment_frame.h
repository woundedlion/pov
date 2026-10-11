/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file pov_segment_frame.h
 * @brief Copies the inactive arm half from the displayed buffer after drawing,
 *        before the frame is queued (buffer-complete hook).
 */
#pragma once

#include "core/render/canvas.h"
#include "hardware/pov_segment_map.h"

#include <cstring>

namespace pov {

/**
 * @brief Selects the arm-A half for the window after the open one.
 * @param window_left True when the current window sweeps arm A's left half.
 * @return True when the queued frame will sweep arm A's left half.
 */
constexpr bool upcoming_window_arm_a_left(bool window_left) {
  return !window_left;
}

/**
 * @brief Retains the previous image outside this frame's segment half.
 * @pre The display clip is one non-wrapping horizontal half, [0, width/2) or
 * [width/2, width), as returned by segment_clip() for an even canvas width.
 * Only rows within the display clip are preserved.
 * @param canvas Canvas whose other half is copied from its previous buffer.
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

/**
 * @brief Draws the initial full-width band before installing segment hooks.
 * @param effect Newly constructed effect.
 * @param initial_clip Segment bounds for the initial arm-A-left window.
 * @param prepare Hook that selects the upcoming window after buffer acquisition.
 */
inline void prime_segment_frame(Effect &effect, const SegmentClip &initial_clip,
                                Effect::BufferReadyHook prepare) {
  const bool CLIPPED =
      segment_clip_applies(effect.needs_full_frame(), effect.persists_pixels());
  if (CLIPPED)
    effect.set_clip(initial_clip.y0, initial_clip.y1, initial_clip.x0,
                    initial_clip.x1);
  effect.set_clip_x(0, effect.width());
  effect.draw_frame();
  effect.set_buffer_ready_hook(prepare);
  if (CLIPPED)
    effect.set_buffer_complete_hook(preserve_segment_half);
}

} // namespace pov
