/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file temporal_history.h
 * @brief Retained frame endpoints, periodic clock arcs, and frame sampling
 *        shared by every column-time renderer.
 */

#include "core/math/3dmath.h"
#include "core/platform/attributes.h"
#include "core/render/pullback/color.h"
#include <cmath>
#include <concepts>
#include <type_traits>

namespace Pullback {

/** @brief Prepared shortest arc of a periodic clock. */
struct ClockInterval {
  float from;  ///< Previous clock value.
  float delta; ///< Shortest signed advance.

  /** @brief Prepares a periodic advance.
   * @param a Previous clock value.
   * @param b Current clock value.
   * @param period Clock period. */
  HS_FLASH_INLINE ClockInterval(float a, float b, float period)
      : from(a), delta(prepare_delta(a, b, period)) {}

  /** @brief Prepares a shortest periodic advance for external storage.
   * @param a Previous clock value.
   * @param b Current clock value.
   * @param period Clock period.
   * @return Signed clock advance. */
  HS_FLASH_INLINE static float prepare_delta(float a, float b, float period) {
    return std::remainder(b - a, period);
  }

  /** @brief Applies a prepared advance.
   * @param from Previous clock value.
   * @param delta Signed clock advance.
   * @param fraction Frame fraction in [0, 1].
   * @return Interpolated clock value. */
  HS_HOT_INLINE static float apply(float from, float delta, float fraction) {
    return from + delta * fraction;
  }

  /** @brief Evaluates the prepared clock.
   * @param fraction Frame fraction in [0, 1].
   * @return Clock value, equivalent modulo its period. */
  HS_HOT_INLINE float sample(float fraction) const {
    return apply(from, delta, fraction);
  }
};

/** @brief Interpolates the shortest arc of a periodic clock.
 * @param from Previous clock value.
 * @param to Current clock value.
 * @param fraction Fraction of one frame, in [0, 1].
 * @param period Clock period; the frame advance must be less than half of it.
 * @return Interpolated clock, equivalent modulo the period. */
inline float sample_clock(float from, float to, float fraction, float period) {
  return ClockInterval(from, to, period).sample(fraction);
}

/** @brief Samples an effect frame that interpolates itself.
 * @details The frame's sampler sees fraction zero, where a reset part of the
 * history takes the current endpoint rather than @p from.
 * @tparam Frame Frame with `static Frame sample(from, to, fraction)` for
 * fractions in [0, 1).
 * @param from Previous frame.
 * @param to Current frame.
 * @param fraction Fraction of one frame, in [0, 1].
 * @return @p to at 1 or above, else the frame's sample. */
template <typename Frame>
  requires requires(const Frame &frame, float fraction) {
    { Frame::sample(frame, frame, fraction) } -> std::same_as<Frame>;
  }
HS_FLASH_INLINE Frame sample_frame(const Frame &from, const Frame &to,
                                   float fraction) {
  if (fraction >= 1.0f)
    return to;
  return Frame::sample(from, to, fraction);
}

/** @brief Whether a frame borrows palette and hue tables its history must own.
 * @tparam Frame Frame context type. */
template <typename Frame>
concept BorrowsColorResources = requires(const Frame &frame) {
  frame.palette;
  frame.temporal;
};

/** @brief Optional retained endpoint storage. */
template <typename Frame, Color::HueMode Hue = Color::HueMode::NONE,
          bool Enabled = true, bool FloatHue = false>
struct TemporalHistory {};

/** @brief Retained previous endpoint of a frame that owns or shares
 *  immutable resources. */
template <typename Frame>
  requires(!BorrowsColorResources<Frame>)
struct TemporalHistory<Frame, Color::HueMode::NONE, true, false> {
  static_assert(std::is_trivially_copyable_v<Frame>);
  Frame previous{};         ///< Previous endpoint.
  bool initialized = false; ///< Whether a previous endpoint exists.

  /** @brief Retains the current endpoint before the next animation step.
   * @param frame Current frame. */
  void capture(const Frame &frame) {
    previous = frame;
    initialized = true;
  }
};

} // namespace Pullback
