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
 * @brief Noise-resource keys for each stage and whether two configurations'
 *        union fits the resident bank.
 */

#include "core/render/pullback/runtime_seeds.h"
#include "workbench/shader/config.h"
#include "workbench/shader/frame_state.h"

namespace Workbench {

HS_COLD_MEMBER inline constexpr math::NoiseFieldKey
warp_resource_key(const WarpStageSpec &spec) {
  return {math::NoiseDomain::PROJECTED_2D,
          spec.basis,
          spec.seed,
          spec.kind == WarpStageKind::CURL_FLOW
              ? math::NoiseChannelLayout::CURL_V1
              : (spec.basis == math::NoiseBasis::SIMPLEX
                     ? math::NoiseChannelLayout::DIRECT_VECTOR_V2
                     : math::NoiseChannelLayout::DIRECT_V1),
          1,
          1,
          static_cast<uint8_t>(spec.kind == WarpStageKind::CURL_FLOW ? 1 : 0),
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr math::NoiseFieldKey
source_resource_key(const Config &config) {
  return {config.slots.function == Function::NOISE_CONTOUR_SPHERE
              ? math::NoiseDomain::SPHERE_3D
              : math::NoiseDomain::PROJECTED_2D,
          config.params.source.noise_basis,
          config.params.source.noise_seed,
          math::NoiseChannelLayout::SCALAR_V1,
          1,
          1,
          0,
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr math::NoiseFieldKey
surface_noise_resource_key(const Config &config) {
  return {math::NoiseDomain::SPHERE_3D,
          config.params.surface_noise.basis,
          config.params.surface_noise.seed,
          config.slots.surface_noise == SurfaceNoise::CURL
              ? (config.params.surface_noise.basis == math::NoiseBasis::SIMPLEX
                     ? math::NoiseChannelLayout::CURL_ANALYTIC_V2
                     : math::NoiseChannelLayout::CURL_V1)
              : (config.params.surface_noise.basis == math::NoiseBasis::SIMPLEX
                     ? math::NoiseChannelLayout::DIRECT_VECTOR_V2
                     : math::NoiseChannelLayout::DIRECT_V1),
          1,
          1,
          static_cast<uint8_t>(
              config.slots.surface_noise == SurfaceNoise::CURL ? 1 : 0),
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr math::NoiseFieldKey color_noise_resource_key() {
  return {math::NoiseDomain::SPHERE_3D,
          math::NoiseBasis::SIMPLEX,
          Pullback::HUE_NOISE_SEED,
          math::NoiseChannelLayout::SCALAR_V1,
          1,
          1,
          0,
          FastNoiseLite::NoiseType_OpenSimplex2,
          1.0f};
}

HS_COLD_MEMBER inline constexpr bool
append_resource_key(const math::NoiseFieldKey &key,
                    std::array<math::NoiseFieldKey, MAX_NOISE_RESOURCES> &keys,
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
    const Config &config,
    std::array<math::NoiseFieldKey, MAX_NOISE_RESOURCES> &keys, size_t &count) {
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

static_assert(
    [] {
      Config from{};
      from.slots.function = Function::NOISE_CONTOUR;
      from.slots.warp_program.outer.kind = WarpStageKind::VECTOR_NOISE;
      from.slots.warp_program.inner.kind = WarpStageKind::CURL_FLOW;
      from.slots.surface_noise = SurfaceNoise::DIRECT;
      from.slots.hue_shift = HueShiftMode::NOISE;
      from.params.color.hue_shift_amount = 1.0f;
      std::array<math::NoiseFieldKey, MAX_NOISE_RESOURCES> keys{};
      size_t count = 0;
      return append_config_resource_keys(from, keys, count) &&
             count == MAX_NOISE_RESOURCES;
    }(),
    "MAX_NOISE_RESOURCES holds every noise key in one parameter topology");

} // namespace Workbench

#endif // HS_ENABLE_SHADER_WORKBENCH
