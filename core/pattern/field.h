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

namespace Pattern {
/** @brief Clamps finite numerical overshoot at the scalar boundary.
 * @param value Raw source sample.
 * @return Sample in [-1, 1]. */
inline float bounded_scalar(float value) {
  assert(std::isfinite(value));
  return hs::clamp(value, -1.0f, 1.0f);
}

/** @brief Direct sphere evaluation of a source policy.
 * @tparam Source Immutable source with a sphere coordinate. */
template <class Source> struct DirectSphereField {
  using Coordinate = typename Source::Coordinate;
  using Prepared = typename Source::Prepared;
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
