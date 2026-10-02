/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <string_view>
#include "vendor/FastNoiseLite.h"

/**
 * @file runtime_seeds.h
 * @brief Deterministic noise and random-walk seeds shared by pullback runtimes.
 */

namespace Pullback {

/** @brief FNV-1a over bytes, optionally continuing an existing hash. */
constexpr uint32_t fnv1a(std::string_view text, uint32_t hash = 2166136261u) {
  for (char c : text)
    hash = (hash ^ static_cast<uint8_t>(c)) * 16777619u;
  return hash;
}

inline constexpr int32_t EFFECT_NOISE_SEED = 1337;
inline constexpr int32_t CAMERA_WALK_SEED = 1337;
inline constexpr int32_t PROJECTION_WALK_SEED = 7331;
inline constexpr int32_t HUE_NOISE_SEED = 6047;

inline void init_effect_noise(FastNoiseLite &noise, int32_t seed) {
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(seed);
  noise.SetFrequency(1.0f);
}

} // namespace Pullback
