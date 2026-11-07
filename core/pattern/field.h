/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file field.h
 * @brief Typed, immutable scalar-field adapters. */

#include <cassert>
#include <cmath>
#include <utility>
#include "core/pattern/types.h"
#include "core/math/periodic.h"

namespace Pattern {
/** @brief Checks a finite source parameter against an admitted interval.
 * @param value Parameter value.
 * @param low Inclusive lower bound.
 * @param high Inclusive upper bound.
 * @return Whether the parameter is admitted. */
inline bool parameter_in_range(float value, float low, float high) {
  return finite(value) && value >= low && value <= high;
}

/** @brief Advances a wrapped clock by a turn increment.
 * @param phase Previous phase in turns.
 * @param turns Increment in turns.
 * @return Updated phase in [0,1). */
HS_HOT_INLINE inline float advance_turn_phase(float phase, float turns) {
  return math::wrap_t(math::wrap_t(phase) + std::fmod(turns, 1.0f));
}

/** @brief Clamps finite numerical overshoot at the scalar boundary.
 * @param value Raw source sample.
 * @return Sample in [-1, 1]. */
HS_HOT_INLINE inline float bounded_scalar(float value) {
  assert(finite(value));
  return hs::clamp(value, -1.0f, 1.0f);
}

/** @brief Direct sphere evaluation of a source policy.
 * @tparam Source Immutable source with a sphere coordinate. */
template <class Source> struct DirectSphereField {
  using Coordinate =
      typename Source::Coordinate;            ///< Admitted sphere coordinate.
  using Prepared = typename Source::Prepared; ///< Immutable source inputs.
  static_assert(std::is_same_v<Coordinate, SphereDirection>);
  /** @brief Prepares the source from its declared inputs.
   * @param inputs Pure source preparation inputs.
   * @return Complete prepared source. */
  template <class... Inputs> static Prepared prepare(Inputs &&...inputs) {
    return Source::prepare(std::forward<Inputs>(inputs)...);
  }
  /** @brief Samples a prepared source.
   * @param coordinate Admitted direction.
   * @param prepared Immutable source inputs.
   * @return Scalar in [-1, 1]. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return bounded_scalar(Source::sample(coordinate, prepared));
  }
};
} // namespace Pattern
