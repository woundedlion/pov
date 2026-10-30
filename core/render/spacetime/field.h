/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file field.h
 * @brief Surface-time shading and retained frame intervals.
 */

#include "core/render/spacetime/surface.h"
#include "core/render/pullback/color.h"
#include <cmath>
#include <concepts>
#include <type_traits>
#include <utility>

namespace Spacetime {

/**
 * @brief Whether a scan shader reads the sweep time of the pixel it shades.
 * @tparam ShaderT Shader type.
 * @tparam Args The scan's shader arguments.
 * @details A shader declares it by its signature: one invocable with a
 * trailing `float u` after the scan's arguments receives, in a column-time
 * scan, the sweep fraction of the shaded column, in [0, 1), and evaluates its
 * time-varying state there. The fixed-meridian scans call the plain form.
 */
template <typename ShaderT, typename... Args>
inline constexpr bool reads_time =
    std::is_invocable_v<ShaderT &, Args..., float>;

/**
 * @brief Runs a shader at a frame fraction when it reads one.
 * @tparam ShaderT Shader type.
 * @tparam Args The shader's arguments.
 * @param shader Shader to run.
 * @param t Frame fraction in [0, 1].
 * @param args The shader's arguments.
 */
template <typename ShaderT, typename... Args>
__attribute__((always_inline)) inline void
shade_at_time(ShaderT &shader, float t, Args &&...args) {
  if constexpr (reads_time<ShaderT, Args...>)
    shader(std::forward<Args>(args)..., t);
  else
    shader(std::forward<Args>(args)...);
}

/**
 * @brief Runs a scan shader at a column's sweep time when it reads one.
 * @tparam W Canvas width in columns.
 * @tparam ShaderT Shader type.
 * @tparam Args The scan's shader arguments.
 * @param shader Shader to run.
 * @param x Shaded canvas column.
 * @param args The scan's shader arguments.
 */
template <int W, typename ShaderT, typename... Args>
__attribute__((always_inline)) inline void
shade_at_column(ShaderT &shader, int x, Args &&...args) {
  shade_at_time(shader, SamplingSurface<W>::time(x),
                std::forward<Args>(args)...);
}

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

/** @brief Value endpoint retained between draws.
 * @tparam T Trivially copyable endpoint without mutable borrowed resources. */
template <typename T> struct Retained {
  static_assert(std::is_trivially_copyable_v<T>);
  T previous{};             ///< Last committed endpoint.
  bool initialized = false; ///< Whether an endpoint has been captured.

  /** @brief Previous endpoint, or the current endpoint before priming and on a flat surface.
   * @param end Current endpoint.
   * @return An endpoint value independent of the argument's lifetime. */
  HS_HOT_INLINE constexpr T start(const T &end) const {
    return initialized && !flat_surface() ? previous : end;
  }

  /** @brief Commits an endpoint for the next draw.
   * @param end Endpoint used by this draw. */
  HS_HOT_INLINE constexpr void capture(const T &end) {
    previous = end;
    initialized = true;
  }

  /** @brief Discards continuity before the next draw. */
  HS_HOT_INLINE constexpr void reset() { initialized = false; }
};

/** @brief Optional retained endpoint storage. */
template <typename Frame, bool Enabled = true, bool Hue = false>
struct FrameHistory {};

/** @brief Retained previous endpoint without mutable borrowed resources.
 * @tparam Frame Frame context.
 * @tparam Hue Whether palette hue resources exist. */
template <typename Frame, bool Hue>
  requires(!BorrowsColorResources<Frame>)
struct FrameHistory<Frame, true, Hue> : Retained<Frame> {};

} // namespace Spacetime
