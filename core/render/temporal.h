/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file temporal.h
 * @brief Column-time sampling contract shared by every renderer family: the
 * fidelity toggle, the sweep fraction of a column, and capture-time solving
 * for plotted points.
 */

#include "core/animation/orientation.h"
#include "core/engine/concepts.h"
#include "core/math/3dmath.h"
#include "core/platform/attributes.h"
#include <algorithm>
#include <array>
#include <cmath>

namespace Temporal {

#if defined(__EMSCRIPTEN__) || HS_ENABLE_TEST_HOOKS
namespace Detail {
inline bool enabled = true; ///< Simulator sweep-time sampling selection.
}

/** @brief Selects sweep-time sampling for simulator presentation.
 * @param enabled Whether columns use their individual emission times. */
inline void set_temporal_sampling_enabled(bool enabled) {
  Detail::enabled = enabled;
}

/** @brief Whether columns use their individual emission times.
 * @return The simulator's temporal sampling mode. */
inline bool temporal_sampling_enabled() { return Detail::enabled; }
#else
/** @brief Whether columns use their individual emission times.
 * @return True for physical display rendering. */
constexpr bool temporal_sampling_enabled() { return true; }
#endif

/** @brief Column time within a two-arm half-turn.
 * @param column Absolute canvas column.
 * @param width Even canvas width.
 * @return Fraction of one effect frame, in [0, 1). */
constexpr float sweep_fraction(int column, int width) {
  return static_cast<float>(column % (width / 2)) / (width / 2);
}

/** @brief Largest azimuth travel, in half-turns per frame, solved by
 * capture(); faster or polar motion is emitted once at the current endpoint. */
inline constexpr float MAX_CAPTURE_DRIFT = 0.5f;

/** @brief Fixed-point convergence of a capture fraction, in half-turns. */
inline constexpr float CAPTURE_TOLERANCE = 1e-3f;

/** @brief Fixed-point refinements after the linear-motion estimate. */
inline constexpr int CAPTURE_ITERATIONS = 3;

/** @brief Frame fractions at which a display arm passes a moving point. */
struct Captures {
  static constexpr int CAPACITY = 2;      ///< Most arm passes in one frame.
  std::array<float, CAPACITY> fraction{}; ///< Frame fraction per pass.
  std::array<math::Vector, CAPACITY> position{}; ///< Point at that fraction.
  int count = 0;                                 ///< Live passes.
};

/** @brief Azimuth of a point in half-turns, so that a column's sweep
 * fraction is this value modulo one.
 * @param p Point on the sphere.
 * @return Azimuth divided by pi, in [-1, 1]. */
__attribute__((always_inline)) inline float half_turns(const math::Vector &p) {
  return math::fast_atan2(p.z, p.x) / math::PI_F;
}

/**
 * @brief Solves `u = (theta(u) / pi) mod 1` for a point moving over one frame.
 * @param position The point at frame fraction u in [0, 1]; unit length.
 * @return Zero, one or two passes with u in [0, 1). Motion faster than
 * `MAX_CAPTURE_DRIFT`, including azimuth swings near a pole, returns one pass
 * at u = 1.
 */
HS_FLASH_MEMBER inline Captures
capture(FunctionRef<math::Vector(float)> position) {
  Captures out;
  const math::Vector end = position(1.0f);
  const float s1 = half_turns(end);
  const float s0 = s1 + std::remainder(half_turns(position(0.0f)) - s1, 2.0f);
  const float drift = s1 - s0;
  if (!(std::fabs(drift) < MAX_CAPTURE_DRIFT)) {
    out.fraction[0] = 1.0f;
    out.position[0] = end;
    out.count = 1;
    return out;
  }
  const int first = static_cast<int>(std::floor(std::min(s0, s1))) - 1;
  const int last = static_cast<int>(std::floor(std::max(s0, s1))) + 1;
  for (int n = first; n <= last && out.count < Captures::CAPACITY; ++n) {
    float u = (s0 - n) / (1.0f - drift);
    if (u < -CAPTURE_TOLERANCE || u >= 1.0f + CAPTURE_TOLERANCE)
      continue;
    u = std::clamp(u, 0.0f, 1.0f);
    math::Vector p;
    float residual = 1.0f;
    for (int i = 0; i < CAPTURE_ITERATIONS && residual > CAPTURE_TOLERANCE;
         ++i) {
      p = position(u);
      const float next =
          s1 + std::remainder(half_turns(p) - s1, 2.0f) - static_cast<float>(n);
      residual = std::fabs(next - u);
      if (residual > CAPTURE_TOLERANCE)
        u = std::clamp(next, 0.0f, 1.0f);
    }
    if (residual <= CAPTURE_TOLERANCE && u < 1.0f) {
      out.fraction[out.count] = u;
      out.position[out.count] = p;
      ++out.count;
    }
  }
  return out;
}

/**
 * @brief Orientation along a frame's history at a frame fraction.
 * @tparam CAP Orientation capacity.
 * @param orientation History whose first entry is the previous frame's end
 *        pose and whose last entry is the current pose.
 * @param fraction Frame fraction in [0, 1].
 * @return Piecewise slerp between the neighbouring history entries.
 */
template <int CAP>
math::Quaternion sample_history(const math::Orientation<CAP> &orientation,
                                float fraction) {
  const int last = orientation.length() - 1;
  if (last == 0)
    return orientation.get();
  const float s = fraction * last;
  const int i = std::min(static_cast<int>(s), last - 1);
  return math::slerp(orientation.get(i), orientation.get(i + 1), s - i);
}

} // namespace Temporal
