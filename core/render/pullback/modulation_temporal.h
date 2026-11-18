/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include "core/render/pullback/modulation.h"
#include <cstring>

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
    HS_CHECK(table.size() == Pattern::SphereLutView::SIZE,
             "modulation retain: table is not one sphere cube");
    HS_CHECK(destination.cube.linear_data == nullptr,
             "modulation retain: a linear cube cannot be retained");
    std::memcpy(table.data(), destination.cube.data, table.size_bytes());
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
  return Detail::visit_field_inputs(
      a.config, a.state.phase, b.config, b.state.phase,
      [](const auto &left, const auto &right) { return left == right; });
}
/** @brief Proves direct scalar equality for distinct immutable endpoints.
 * @param a First endpoint.
 * @param b Second endpoint.
 * @return Whether the same direct evaluation can serve both endpoints. */
HS_SIZE_FLASH_MEMBER inline bool same_direct_sample(const Prepared &a,
                                                    const Prepared &b) {
  if (&a == &b)
    return true;
  if (!a.evaluate || a.evaluate != b.evaluate ||
      a.source.index() != b.source.index())
    return false;
  const size_t slot = a.source.index();
  if (slot >= Detail::DIRECT_SAMPLERS.size() ||
      a.evaluate != Detail::DIRECT_SAMPLERS[slot]) {
    // Custom samplers may read any part of their prepared endpoint.
    return std::memcmp(&a, &b, sizeof(a)) == 0;
  }
  if (!same_field(a, b))
    return false;
  if (slot == 0) {
    const auto &x = std::get<0>(a.source);
    const auto &y = std::get<0>(b.source);
    return x.scale == y.scale && x.contrast == y.contrast &&
           x.basis == y.basis && x.variant == y.variant &&
           x.noise.same_noise_single_3d(y.noise) &&
           std::memcmp(&x.loop_offset, &y.loop_offset, sizeof(x.loop_offset)) ==
               0;
  }
  return [&]<size_t... I>(std::index_sequence<I...>) {
    return ((slot == I + 1 &&
             std::memcmp(&std::get<I + 1>(a.source), &std::get<I + 1>(b.source),
                         sizeof(std::get<I + 1>(a.source))) == 0) ||
            ...);
  }(std::make_index_sequence<std::variant_size_v<Source> - 1>{});
}

/** @brief Checks when color interpolation can use a single prepared field.
 * @param a Retained scalar source.
 * @param b Current scalar source.
 * @return Whether both endpoint samples are equal. */
HS_SIZE_FLASH_MEMBER inline bool same_endpoint_sample(const Prepared &a,
                                                      const Prepared &b) {
  if (a.evaluate != sample_cube || b.evaluate != sample_cube)
    return same_direct_sample(a, b);
  if (!a.cube.active || !b.cube.active || !same_field(a, b))
    return false;
  if (a.cube.linear_data || b.cube.linear_data)
    return a.cube.linear_data && b.cube.linear_data &&
           (a.cube.linear_data == b.cube.linear_data ||
            std::memcmp(a.cube.linear_data, b.cube.linear_data,
                        Pattern::SphereLutView::SIZE * sizeof(float)) == 0);
  return a.cube.data && b.cube.data &&
         (a.cube.data == b.cube.data ||
          std::memcmp(a.cube.data, b.cube.data, Pattern::SphereLutView::SIZE) ==
              0);
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
