/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include "core/render/pullback/modulation.h"

/** @file modulation_temporal.h
 * @brief Immutable modulation endpoint equality and completed-color blending. */
namespace Pullback::Modulation {
/** @brief Compares all source inputs independently of palette and mechanic.
 * @param a Previous endpoint.
 * @param b Current endpoint.
 * @return Whether both evaluate the same scalar field. */
inline bool same_field(const Prepared &a, const Prepared &b) {
  const auto identity = [](const Prepared &value) {
    const Config &c = value.config;
    return FieldIdentity{c.source, c.version, c.evaluation,
                         c.noise,  c.sine,    c.interference,
                         c.grid,   c.seed,    value.state.phase};
  };
  return identity(a) == identity(b);
}
/** @brief Interpolates completed straight-alpha colors with RGB16 rounding.
 * @param a Previous completed endpoint.
 * @param b Current completed endpoint.
 * @param fraction Original floating-point frame fraction.
 * @param weight Q16 RGB interpolation weight.
 * @return Independently interpolated RGB and alpha. */
inline Color4 interpolate_colors(const Color4 &a, const Color4 &b,
                                 float fraction, uint16_t weight) {
  return Color4(a.color.lerp16(b.color, weight),
                hs::lerp(a.alpha, b.alpha, fraction));
}
} // namespace Pullback::Modulation
