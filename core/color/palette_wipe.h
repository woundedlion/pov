/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette_wipe.h
 * @brief Palette snapshot transitions and their rebake window.
 */

#include "color/generative_palette.h"
#include "color/baked_palette.h"

/**
 * @brief Steps a palette cross-fade rebake for one frame.
 * @tparam Source Palette type exposing Color4 get(float) const.
 * @param wipe_pending Set when a wipe was just armed; consumes the arming frame.
 * @param wipe_frames_remaining Frames left to rebake; decremented per step.
 * @param baked LUT rebaked from the wipe-mutated source while the wipe runs.
 * @param source Palette mutated in place by the in-flight ColorWipe.
 * @details A ColorWipe is armed mid-step and first steps next frame, so the
 * arming frame is skipped; each later frame rebakes the LUT the shader samples.
 */
template <typename Source>
inline void step_wipe_rebake(bool &wipe_pending, int &wipe_frames_remaining,
                             BakedPalette &baked, const Source &source) {
  if (wipe_pending) {
    wipe_pending = false;
  } else if (wipe_frames_remaining > 0) {
    baked.rebake(source);
    --wipe_frames_remaining;
  }
}

/**
 * @brief One palette cross-fade: the snapshots a ColorWipe interpolates and
 *        the rebake window that follows it.
 * @details The snapshots must outlive the wipe — a ColorWipe holds references
 * to them — so a caller whose rollover can fire mid-wipe gates on in_flight().
 */
struct PaletteWipe {
  GenerativePalette::Snapshot start;  /**< Current wipe's start. */
  GenerativePalette::Snapshot target; /**< Current wipe's target. */
  int frames_remaining =
      0; /**< Frames left to rebake the LUT for the in-flight wipe. */
  bool pending = false; /**< Armed this frame; it first steps next frame. */

  /** @brief True while an armed wipe still has rebake frames left. */
  bool in_flight() const { return frames_remaining > 0; }

  /**
   * @brief Captures the cross-fade endpoints and opens the rebake window.
   * @param from Palette the wipe starts from.
   * @param to Palette the wipe lands on.
   * @param frames Cross-fade duration, in frames.
   * @note Arms the state only; the caller schedules the ColorWipe over
   *       `start` and `target`.
   */
  void arm(const GenerativePalette &from, const GenerativePalette::Snapshot &to,
           int frames) {
    start = from.snapshot();
    target = to;
    frames_remaining = frames;
    pending = true;
  }

  /**
   * @brief Steps one frame of the rebake window.
   * @param baked LUT rebaked from the wipe-mutated source.
   * @param source Palette mutated in place by the in-flight ColorWipe.
   */
  template <typename Source>
  void step(BakedPalette &baked, const Source &source) {
    step_wipe_rebake(pending, frames_remaining, baked, source);
  }
};
