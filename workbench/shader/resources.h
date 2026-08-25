/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_SHADER_WORKBENCH

/**
 * @file resources.h
 * @brief Noise-resource identity: the field key each stage of a configuration needs, and whether the union of two configurations fits the resident bank.
 */

#include "core/render/pullback/runtime_seeds.h"
#include "workbench/shader/config.h"
#include "workbench/shader/frame_state.h"

namespace Workbench {

HS_COLD_MEMBER inline constexpr NoiseFieldKey
warp_resource_key(const WarpStageSpec &spec) {
  return {NoiseDomain::PROJECTED_2D,
          spec.basis,
          spec.seed,
          spec.kind == WarpStageKind::CURL_FLOW
              ? NoiseChannelLayout::CURL_V1
              : (spec.basis == NoiseBasis::SIMPLEX
                     ? NoiseChannelLayout::DIRECT_VECTOR_V2
                     : NoiseChannelLayout::DIRECT_V1),
          1,
          1,
          static_cast<uint8_t>(spec.kind == WarpStageKind::CURL_FLOW ? 1 : 0),
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr NoiseFieldKey
source_resource_key(const Config &config) {
  return {config.slots.function == Function::NOISE_CONTOUR_SPHERE
              ? NoiseDomain::SPHERE_3D
              : NoiseDomain::PROJECTED_2D,
          config.params.source.noise_basis,
          config.params.source.noise_seed,
          NoiseChannelLayout::SCALAR_V1,
          1,
          1,
          0,
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr NoiseFieldKey
surface_noise_resource_key(const Config &config) {
  return {NoiseDomain::SPHERE_3D,
          config.params.surface_noise.basis,
          config.params.surface_noise.seed,
          config.slots.surface_noise == SurfaceNoise::CURL
              ? (config.params.surface_noise.basis == NoiseBasis::SIMPLEX
                     ? NoiseChannelLayout::CURL_ANALYTIC_V2
                     : NoiseChannelLayout::CURL_V1)
              : (config.params.surface_noise.basis == NoiseBasis::SIMPLEX
                     ? NoiseChannelLayout::DIRECT_VECTOR_V2
                     : NoiseChannelLayout::DIRECT_V1),
          1,
          1,
          static_cast<uint8_t>(
              config.slots.surface_noise == SurfaceNoise::CURL ? 1 : 0),
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr NoiseFieldKey color_noise_resource_key() {
  return {NoiseDomain::SPHERE_3D,
          NoiseBasis::SIMPLEX,
          Pullback::HUE_NOISE_SEED,
          NoiseChannelLayout::SCALAR_V1,
          1,
          1,
          0,
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr bool
append_resource_key(const NoiseFieldKey &key,
                    std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> &keys,
                    size_t &count) {
  for (size_t index = 0; index < count; ++index)
    if (keys[index] == key)
      return true;
  if (count == keys.size())
    return false;
  keys[count++] = key;
  return true;
}

HS_COLD_MEMBER inline constexpr bool append_config_resource_keys(
    const Config &config, std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> &keys,
    size_t &count) {
  if (warp_uses_noise(config.slots.warp_program.outer.kind) &&
      !append_resource_key(warp_resource_key(config.slots.warp_program.outer),
                           keys, count))
    return false;
  if (warp_uses_noise(config.slots.warp_program.inner.kind) &&
      !append_resource_key(warp_resource_key(config.slots.warp_program.inner),
                           keys, count))
    return false;
  if (is_noise_contour(config.slots.function) &&
      !append_resource_key(source_resource_key(config), keys, count))
    return false;
  if (config.slots.surface_noise != SurfaceNoise::NONE &&
      !append_resource_key(surface_noise_resource_key(config), keys, count))
    return false;
  if (config.slots.hue_shift == HueShiftMode::NOISE &&
      config.params.color.hue_shift_amount != 0.0f &&
      !append_resource_key(color_noise_resource_key(), keys, count))
    return false;
  return true;
}

HS_COLD_MEMBER inline constexpr bool resource_union_fits(const Config &from,
                                                         const Config &to) {
  std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> keys{};
  size_t count = 0;
  return append_config_resource_keys(from, keys, count) &&
         append_config_resource_keys(to, keys, count);
}

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
