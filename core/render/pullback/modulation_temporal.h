/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include "core/render/pullback/modulation.h"

/** @file modulation_temporal.h
 * @brief Immutable modulation endpoint retention, equality and completed-color blending. */
namespace Pullback::Modulation {
/** @brief Retains an endpoint and rebinds its active cube to owned storage.
 * @param source Current endpoint, or null.
 * @param destination Retained endpoint storage.
 * @param table Retained byte-cube storage.
 * @return Retained endpoint, or null when the source is absent. */
inline const Prepared *retain(const Prepared *source, Prepared &destination,
                              std::span<int8_t> table) {
  if (source == nullptr)
    return nullptr;
  destination = *source;
  if (destination.cube.active) {
    assert(table.size() == Pattern::SphereLutView::SIZE);
    assert(destination.cube.linear_data == nullptr);
    std::copy_n(destination.cube.data, table.size(), table.data());
    destination.cube.data = table.data();
  }
  return &destination;
}
/** @brief Compares all source inputs independently of palette and mechanic.
 * @param a Previous endpoint.
 * @param b Current endpoint.
 * @return Whether both evaluate the same scalar field. */
HS_HOT_FLASH_MEMBER inline bool same_field(const Prepared &a,
                                           const Prepared &b) {
  const auto &start = a.config;
  const auto &end = b.config;
  if (start.source != end.source || start.version != end.version ||
      start.evaluation != end.evaluation || a.state.phase != b.state.phase)
    return false;
  switch (start.source) {
  case Pattern::SourceKind::DISPLACEMENT:
    return true;
  case Pattern::SourceKind::NOISE:
    return start.seed == end.seed && start.noise == end.noise;
  case Pattern::SourceKind::DIRECTIONAL_SINE:
    return start.sine == end.sine;
  case Pattern::SourceKind::SPHERE_INTERFERENCE:
    return start.interference == end.interference;
  case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
    return start.grid == end.grid;
  }
  return false;
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
