/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file pov_segment_frame.h
 * @brief Segment-frame priming, scheduling, opposite-half retention and scanout.
 */
#pragma once

#include "core/render/canvas.h"
#include "hardware/pov_handoff.h"
#include "hardware/pov_segment_map.h"

#include <cstring>

namespace pov {

/**
 * @brief Whether @p effect renders only its segment_clip() rectangle.
 * @param effect Effect whose needs_full_frame() and persists_pixels() gate the
 *        clip.
 * @return False when the effect reads cross-segment or prior-frame state.
 */
inline bool segment_clip_applies(const Effect &effect) {
  return segment_clip_applies(effect.needs_full_frame(),
                              effect.persists_pixels());
}

/**
 * @brief Sets the next logical frame and its display clip.
 * @param effect Effect with a newly acquired back buffer.
 * @param segment Physical arm and row mapping.
 * @param leds Total LEDs across the two arms.
 * @param segments Number of boards.
 * @param window Open display window.
 */
inline void prepare_segment_frame(Effect &effect, const SegmentMap &segment,
                                  int leds, int segments,
                                  DisplayWindow window) {
  const uint32_t open = window.index + 1;
  const uint32_t target =
      effect.can_skip_frames() ? open + 1 : effect.scheduled_frame() + 1;
  effect.set_scheduled_frame(target);
  if (!segment_clip_applies(effect))
    return;
  const bool target_left = window.arm_a_left != ((target - open) % 2 != 0);
  const SegmentClip clip =
      segment_clip(segment, target_left, leds, segments, effect.width());
  effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
}

/**
 * @brief Publishes a completed segment frame when its window is due or past.
 * @param effect Live effect with a completed or in-progress frame.
 * @param handoff Display window state.
 */
inline void publish_segment_frame(Effect &effect,
                                  const EffectHandoff<Effect> &handoff) {
  effect.advance_display_if([&] {
    const uint32_t target = effect.scheduled_frame();
    return target != 0 && target <= handoff.current_window().index + 1;
  });
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
  const bool CLIPPED = segment_clip_applies(effect);
  if (CLIPPED)
    effect.set_clip(initial_clip.y0, initial_clip.y1, initial_clip.x0,
                    initial_clip.x1);
  effect.set_clip_x(0, effect.width());
  effect.set_scheduled_frame(1);
  effect.draw_frame();
  effect.set_buffer_ready_hook(prepare);
  if (CLIPPED)
    effect.set_buffer_complete_hook(preserve_segment_half);
}

} // namespace pov
