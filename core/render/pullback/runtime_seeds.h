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

/** @brief FNV-1a over bytes, optionally continuing an existing hash.
 *  @param text Bytes to hash.
 *  @param hash Running hash; defaults to the FNV-1a offset basis.
 *  @return The updated hash. */
constexpr uint32_t fnv1a(std::string_view text, uint32_t hash = 2166136261u) {
  for (char c : text)
    hash = (hash ^ static_cast<uint8_t>(c)) * 16777619u;
  return hash;
}

/// Default seed for effect noise fields.
inline constexpr int32_t EFFECT_NOISE_SEED = 1337;
/// Seed for the camera orientation random walk.
inline constexpr int32_t CAMERA_WALK_SEED = 1337;
/// Seed for the projection-parameter random walk.
inline constexpr int32_t PROJECTION_WALK_SEED = 7331;
/// Seed for hue-modulation noise.
inline constexpr int32_t HUE_NOISE_SEED = 6047;

/** @brief Configures @p noise as unit-frequency OpenSimplex2 with @p seed.
 * @param noise Generator to configure in place.
 * @param seed Noise seed. */
inline void init_effect_noise(FastNoiseLite &noise, int32_t seed) {
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetSeed(seed);
  noise.SetFrequency(1.0f);
}

} // namespace Pullback
