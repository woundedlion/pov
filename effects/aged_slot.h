/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/platform.h"

/**
 * @brief Shared age-driven lifecycle for recyclable ring slots.
 * @details RingShower and Thrusters both advance a slot by one `age` per frame,
 *          grow its radius linearly from an age+1 step, and recycle the slot once
 *          it outlives its bound. The growth window and life bound differ per
 *          effect (RingShower draws a random per-slot life; Thrusters uses a
 *          compile-time constant), so each supplies them at the call rather than
 *          storing a `life` here.
 */
struct AgedSlot {
  int age = 0; /**< Frames elapsed since (re)spawn. */

  /**
   * @brief Linear growth radius using the age+1 endpoint convention.
   * @param max_radius Peak radius reached at the end of the growth window.
   * @param grow_frames Frames over which the radius grows from a non-zero first
   *        step to max_radius, after which it holds.
   * @return Radius in [0, max_radius]: non-zero on the first drawn frame and
   *         exactly max_radius once age + 1 == grow_frames.
   */
  float radius_at(float max_radius, int grow_frames) const {
    float t = hs::clamp(static_cast<float>(age + 1) / grow_frames, 0.0f, 1.0f);
    return max_radius * t;
  }

  /**
   * @brief Whether the slot has reached the end of its life and is recyclable.
   * @param life Total visible frames for this slot.
   */
  bool expired(int life) const { return age >= life; }
};
