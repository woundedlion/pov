/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file palette.h
 * @brief Runtime palette interface and source traits.
 */

#include "color/pixel.h"

/**
 * @brief Abstract base for all palettes.
 * @details Uniform color-lookup interface via a single vtable pointer.
 */
class Palette {
public:
  /**
   * @brief Samples the palette at a coordinate.
   * @param t Lookup coordinate, conventionally in [0, 1].
   * @return The color at t.
   */
  virtual Color4 get(float t) const = 0;
  /**
   * @brief Virtual destructor for polymorphic deletion.
   */
  virtual ~Palette() = default;
};

/**
 * @brief Reports whether a palette source wraps its lookup coordinate.
 * @tparam T Source or composition type.
 * @return T::WRAPS_COORDINATE, or false for a source that declares no marker.
 */
template <typename T> constexpr bool palette_wraps_coordinate() {
  if constexpr (requires { T::WRAPS_COORDINATE; })
    return T::WRAPS_COORDINATE;
  else
    return false;
}
