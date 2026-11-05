/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file noise_basis.h
 * @brief Dependency-free octave noise formula selection. */

#include <cstdint>

namespace math {

/** @brief Octave structure of one scalar sample. */
enum class NoiseBasis : uint8_t {
  /** One generator sample. */
  SIMPLEX,
  /** Three octaves at 1x, 2x, 4x, weighted 1, 1/2, 1/4. */
  FBM3,
  /** The same three octaves folded through 1 - |n|. */
  RIDGED3
};

/** @brief Serialized noise-basis ids in enumerator order. */
inline constexpr const char *NOISE_BASIS_IDS[] = {"simplex", "fbm3", "ridged3"};
static_assert(static_cast<uint8_t>(NoiseBasis::SIMPLEX) == 0 &&
              static_cast<uint8_t>(NoiseBasis::FBM3) == 1 &&
              static_cast<uint8_t>(NoiseBasis::RIDGED3) == 2);

} // namespace math
