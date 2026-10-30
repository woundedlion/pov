/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_animation.h
 * @brief Ripple burst scheduling. */

#include "core/render/pullback/surface.h"
#include "core/render/pullback/runtime_seeds.h"

namespace Pullback::SphereAnimation {

/** @brief Four-wave ripple burst controls. */
struct BurstRippleParams : Surface::PeriodicRippleParams {
  float burst_rate = 1.0f / 128.0f; ///< Bursts per frame.
  uint8_t burst_size = 3;           ///< Ripples per burst minus one.
  uint8_t random_center = 1; ///< One random center per burst when enabled.
  static constexpr auto FIELDS = concat_fields<BurstRippleParams>(
      Surface::PeriodicRippleParams::FIELDS,
      std::array{Field<BurstRippleParams>{
          "burst-rate", &BurstRippleParams::burst_rate, "Burst Rate", 0.0f,
          1.0f / 64.0f, FieldCurve::LERP}});
};

/** @brief One scheduled wave in a ripple burst. */
struct RippleWave {
  float age = 1000.0f; ///< Frames since launch; negative before launch.
  uint32_t burst = 0;  ///< Burst identity shared by its waves.
};
/** @brief Ripple burst scheduler state. */
struct RippleState {
  float phase = 0;          ///< Fraction of the interval until the next burst.
  uint32_t seed = 0;        ///< Instance seed for burst origins.
  uint32_t burst_count = 0; ///< Bursts launched.
  uint32_t frames = 0;      ///< Frames advanced.
  std::array<RippleWave, 16> waves{}; ///< Scheduled ripple waves.
};

/** @brief Advances the scheduler by one frame.
 * @tparam Params Ripple controls.
 * @param state Scheduler state.
 * @param params Current controls. */
template <typename Params>
inline void advance_ripple(RippleState &state, const Params &params) {
  ++state.frames;
  for (auto &wave : state.waves)
    wave.age = fminf(wave.age + 1.0f, 1000.0f);
  state.phase += params.burst_rate;
  const bool FIRST = state.burst_count == 0;
  if (!FIRST && state.phase < 1.0f)
    return;
  if (!FIRST)
    state.phase -= 1.0f;
  const float AGE = FIRST ? 1.0f : state.phase / params.burst_rate;
  const size_t GROUP = (state.burst_count % 4) * 4;
  for (uint8_t i = 0; i <= params.burst_size; ++i)
    state.waves[GROUP + i] = {AGE - 16.0f * i, state.burst_count};
  ++state.burst_count;
}

/** @brief Deterministic uniform center for a burst.
 * @param seed Instance seed.
 * @param burst Burst identity.
 * @return Unit center direction. */
inline math::Vector ripple_center(uint32_t seed, uint32_t burst) {
  uint32_t hash = seed ^ (burst * 0x9e3779b9U);
  hash = (hash ^ (hash >> 16)) * 0x7feb352dU;
  hash = (hash ^ (hash >> 15)) * 0x846ca68bU;
  hash ^= hash >> 16;
  const float Y = 2.0f * static_cast<float>(hash & 0xffffU) / 65535.0f - 1.0f;
  const float AZIMUTH =
      static_cast<float>(hash >> 16) * (math::TWO_PI_F / 65536.0f);
  const float RADIUS = sqrtf(fmaxf(0.0f, 1.0f - Y * Y));
  return {RADIUS * cosf(AZIMUTH), Y, RADIUS * sinf(AZIMUTH)};
}

/** @brief Net rotation and absolute travel along one radial ray. */
struct RadialTransfer {
  float angle = 0.0f; ///< Signed rotation in radians.
  float path =
      0.0f; ///< Absolute angular travel; negative marks a pole crossing.
};

/** @brief Endpoint transfer curves for one burst center. */
struct RadialGroup {
  static constexpr size_t SAMPLES = 257; ///< Uniform samples over [0, pi].
  math::Vector center;                   ///< Shared burst center.
  float cos_threshold_min =
      -1.0f;                      ///< Innermost support bound across endpoints.
  float cos_threshold_max = 1.0f; ///< Outermost support bound across endpoints.
  std::array<RadialTransfer, SAMPLES> previous; ///< One-frame-old wave ages.
  std::array<RadialTransfer, SAMPLES> current;  ///< Current wave ages.
};

/** @brief Persistent radial curves baked once for a frame. */
struct RadialTables {
  std::array<RadialGroup, 4> groups; ///< Burst groups in launch order.
  uint8_t count = 0;                 ///< Active group count.
  bool enabled = false; ///< Whether the parameters admit radial interpolation.
};

/** @brief Active ripple kernels in launch order. */
struct PreparedRipples {
  std::array<Surface::PreparedRipple, 16> waves; ///< Active kernels.
  uint8_t count = 0;                             ///< Active wave count.
};

/** @brief Applies prepared waves in launch order.
 * @param input Unit direction.
 * @param prepared Prepared wave kernels.
 * @param track_path Whether to accumulate angular travel.
 * @return Displaced direction and travel. */
HS_HOT_ITCM_MEMBER inline SurfaceResult
sequential_ripples(const math::Vector &input, const PreparedRipples &prepared,
                   bool track_path) {
  SurfaceResult result{input, 0.0f};
  for (uint8_t index = 0; index < prepared.count; ++index) {
    const math::Vector START = result.sphere;
    result.sphere = ripple_transform(START, prepared.waves[index].ripple);
    if (track_path) {
      const math::Vector CROSS = math::cross(START, result.sphere);
      result.path_length += atan2f(sqrtf(math::dot(CROSS, CROSS)),
                                   math::dot(START, result.sphere));
    }
  }
  return result;
}

/** @brief Groups consecutive waves sharing a center into one radial basis.
 * @param input Unit direction.
 * @param prepared Prepared wave kernels in launch order.
 * @param track_path Whether to accumulate angular travel.
 * @return Displaced direction and travel. */
HS_HOT_ITCM_MEMBER inline SurfaceResult
grouped_ripples(const math::Vector &input, const PreparedRipples &prepared,
                bool track_path) {
  SurfaceResult result{input, 0.0f};
  uint8_t index = 0;
  while (index < prepared.count) {
    const math::Vector START = result.sphere;
    const math::Vector CENTER = prepared.waves[index].ripple.center;
    const float DISTANCE_COS = math::dot(CENTER, START);
    float distance = math::fast_acos(hs::clamp(DISTANCE_COS, -1.0f, 1.0f));
    float total_angle = 0.0f;
    do {
      const auto &wave = prepared.waves[index++].ripple;
      const float HALF_WIDTH = wave.half_width();
      if (wave.amplitude <= 0.001f ||
          distance < wave.phase - 2.0f * HALF_WIDTH ||
          distance > wave.phase + 2.0f * HALF_WIDTH)
        continue;
      const float T = (distance - wave.phase) / HALF_WIDTH * 2.0f;
      const float THETA =
          wave.amplitude * (1.0f - T * T) *
          math::fast_expf(-0.5f * T * T - wave.decay * distance);
      if (THETA == 0.0f)
        continue;
      distance += THETA;
      if (!(distance > 0.01f && distance < math::PI_F - 0.01f))
        return sequential_ripples(input, prepared, track_path);
      total_angle += THETA;
      if (track_path)
        result.path_length += fabsf(THETA);
    } while (index < prepared.count &&
             prepared.waves[index].ripple.center.x == CENTER.x &&
             prepared.waves[index].ripple.center.y == CENTER.y &&
             prepared.waves[index].ripple.center.z == CENTER.z);
    if (total_angle != 0.0f) {
      const math::Vector CROSS = math::cross(CENTER, START);
      const float CROSS_SQ = math::dot(CROSS, CROSS);
      if (!(CROSS_SQ >= 1e-4f) ||
          !(fabsf(math::dot(START, START) - 1.0f) <= 1e-4f))
        return sequential_ripples(input, prepared, track_path);
      const math::Vector AXIS = CROSS * (1.0f / sqrtf(CROSS_SQ));
      if (fabsf(total_angle) <= RIPPLE_SMALL_ANGLE_MAX) {
        const float H = total_angle * 0.5f, H2 = H * H;
        const float S = H * (1.0f - H2 * (1.0f / 6.0f));
        const float C = 1.0f - H2 * (0.5f - H2 * (1.0f / 24.0f));
        result.sphere = math::rotate(START, math::Quaternion(C, S * AXIS));
      } else
        result.sphere =
            math::rotate(START, math::make_rotation(AXIS, total_angle));
    }
  }
  return result;
}

/** @brief Evaluates a same-center wave sequence in radial coordinates.
 * @param distance Initial angular distance.
 * @param waves Prepared waves in launch order.
 * @param count Number of waves.
 * @return Net rotation and travel, or a negative path at a pole crossing. */
inline RadialTransfer
radial_transfer(float distance,
                const std::array<Surface::PreparedRipple, 4> &waves,
                uint8_t count) {
  RadialTransfer result;
  for (uint8_t i = 0; i < count; ++i) {
    const auto &wave = waves[i].ripple;
    const float HALF_WIDTH = wave.half_width();
    if (wave.amplitude <= 0.001f || distance < wave.phase - 2.0f * HALF_WIDTH ||
        distance > wave.phase + 2.0f * HALF_WIDTH)
      continue;
    const float T = (distance - wave.phase) / HALF_WIDTH * 2.0f;
    const float THETA = wave.amplitude * (1.0f - T * T) *
                        math::fast_expf(-0.5f * T * T - wave.decay * distance);
    if (THETA == 0.0f)
      continue;
    distance += THETA;
    if (!(distance > 0.01f && distance < math::PI_F - 0.01f))
      return {0.0f, -1.0f};
    result.angle += THETA;
    result.path += fabsf(THETA);
  }
  return result;
}

/** @brief Bakes both temporal endpoints of each burst's radial transfer.
 * @tparam Params Ripple controls.
 * @param params Current controls.
 * @param state Current scheduler state.
 * @param tables Persistent destination. */
template <typename Params>
inline void prepare_radial_tables(const Params &params,
                                  const RippleState &state,
                                  RadialTables &tables) {
  tables.count = 0;
  tables.enabled = params.strength >= 0.0f && params.strength <= 0.15f &&
                   params.thickness >= 0.7f && params.thickness <= 1.4f &&
                   params.period >= 60.0f && params.period <= 143.0f &&
                   params.decay >= 0.0f && params.decay <= 5.0f;
  if (!tables.enabled)
    return;
  for (uint32_t lag = 4; lag > 0; --lag) {
    const uint32_t BURST = state.burst_count - lag;
    std::array<Surface::PreparedRipple, 4> previous, current;
    uint8_t previous_count = 0, current_count = 0;
    for (size_t i = 0; i < 4; ++i) {
      const auto &wave = state.waves[(BURST % 4) * 4 + i];
      if (wave.burst != BURST)
        continue;
      if (wave.age >= 0.0f && wave.age < params.period)
        current[current_count++] =
            Surface::prepare_ripple(params, wave.age / params.period);
      const float PREVIOUS_AGE = wave.age - 1.0f;
      if (PREVIOUS_AGE >= 0.0f && PREVIOUS_AGE < params.period)
        previous[previous_count++] =
            Surface::prepare_ripple(params, PREVIOUS_AGE / params.period);
    }
    if (!previous_count && !current_count)
      continue;
    auto &group = tables.groups[tables.count];
    group.cos_threshold_min = -1.0f;
    group.cos_threshold_max = 1.0f;
    bool active = false;
    const auto include_support = [&](const auto &waves, uint8_t count) {
      for (uint8_t i = 0; i < count; ++i) {
        const auto &wave = waves[i].ripple;
        if (wave.amplitude <= 0.001f)
          continue;
        active = true;
        group.cos_threshold_min =
            fmaxf(group.cos_threshold_min, wave.cos_threshold_min);
        group.cos_threshold_max =
            fminf(group.cos_threshold_max, wave.cos_threshold_max);
      }
    };
    include_support(previous, previous_count);
    include_support(current, current_count);
    if (!active)
      continue;
    ++tables.count;
    group.center = params.random_center
                       ? ripple_center(state.seed, BURST)
                       : math::Vector::from_spherical(params.center_azimuth,
                                                      params.center_polar);
    for (size_t i = 0; i < RadialGroup::SAMPLES; ++i) {
      const float DISTANCE = math::PI_F * i / (RadialGroup::SAMPLES - 1);
      group.previous[i] = radial_transfer(DISTANCE, previous, previous_count);
      group.current[i] = radial_transfer(DISTANCE, current, current_count);
    }
  }
}

/** @brief Applies interpolated radial curves with grouped-kernel fallback.
 * @param input Unit direction.
 * @param prepared Current column's exact wave kernels.
 * @param track_path Whether to accumulate angular travel.
 * @param tables Optional current-owned transfer curves.
 * @param fraction Current endpoint weight.
 * @return Displaced direction and travel. */
HS_HOT_ITCM_MEMBER inline SurfaceResult
radial_ripples(const math::Vector &input, const PreparedRipples &prepared,
               bool track_path, const RadialTables *tables, float fraction) {
  if (!tables || !tables->enabled)
    return grouped_ripples(input, prepared, track_path);
  SurfaceResult result{input, 0.0f};
  for (uint8_t group_index = 0; group_index < tables->count; ++group_index) {
    const auto &group = tables->groups[group_index];
    const math::Vector START = result.sphere;
    const float DISTANCE_COS = math::dot(group.center, START);
    if (DISTANCE_COS > group.cos_threshold_min ||
        DISTANCE_COS < group.cos_threshold_max)
      continue;
    const float DISTANCE =
        math::fast_acos(hs::clamp(DISTANCE_COS, -1.0f, 1.0f));
    const float POSITION = DISTANCE * ((RadialGroup::SAMPLES - 1) / math::PI_F);
    if (!(POSITION >= 0.0f && POSITION <= RadialGroup::SAMPLES - 1))
      return grouped_ripples(input, prepared, track_path);
    const size_t LOW_INDEX =
        std::min(static_cast<size_t>(POSITION), RadialGroup::SAMPLES - 2);
    const float WEIGHT = POSITION - LOW_INDEX;
    const auto &a = group.previous[LOW_INDEX];
    const auto &b = group.previous[LOW_INDEX + 1];
    const auto &c = group.current[LOW_INDEX];
    const auto &d = group.current[LOW_INDEX + 1];
    if (a.path < 0.0f || b.path < 0.0f || c.path < 0.0f || d.path < 0.0f)
      return grouped_ripples(input, prepared, track_path);
    const float ANGLE = hs::lerp(hs::lerp(a.angle, b.angle, WEIGHT),
                                 hs::lerp(c.angle, d.angle, WEIGHT), fraction);
    if (track_path)
      result.path_length +=
          hs::lerp(hs::lerp(a.path, b.path, WEIGHT),
                   hs::lerp(c.path, d.path, WEIGHT), fraction);
    if (ANGLE == 0.0f)
      continue;
    const math::Vector CROSS = math::cross(group.center, START);
    const float CROSS_SQ = math::dot(CROSS, CROSS);
    if (!(CROSS_SQ >= 1e-4f) ||
        !(fabsf(math::dot(START, START) - 1.0f) <= 1e-4f))
      return grouped_ripples(input, prepared, track_path);
    const math::Vector AXIS = CROSS * (1.0f / sqrtf(CROSS_SQ));
    if (fabsf(ANGLE) <= RIPPLE_SMALL_ANGLE_MAX) {
      const float H = ANGLE * 0.5f, H2 = H * H;
      const float S = H * (1.0f - H2 * (1.0f / 6.0f));
      const float C = 1.0f - H2 * (0.5f - H2 * (1.0f / 24.0f));
      result.sphere = math::rotate(START, math::Quaternion(C, S * AXIS));
    } else
      result.sphere = math::rotate(START, math::make_rotation(AXIS, ANGLE));
  }
  return result;
}

/** @brief Prepares the active wave kernels.
 * @tparam Params Ripple controls.
 * @param params Current controls.
 * @param state Scheduler state.
 * @return Active kernels in launch order. */
template <typename Params>
inline PreparedRipples prepare_ripples(const Params &params,
                                       const RippleState &state) {
  PreparedRipples result;
  for (uint32_t lag = 4; lag > 0; --lag) {
    const uint32_t BURST = state.burst_count - lag;
    for (size_t i = 0; i < 4; ++i) {
      const auto &wave = state.waves[(BURST % 4) * 4 + i];
      if (wave.age < 0.0f || wave.age >= params.period || wave.burst != BURST)
        continue;
      auto prepared = Surface::prepare_ripple(params, wave.age / params.period);
      if (params.random_center)
        prepared.ripple.center = ripple_center(state.seed, BURST);
      result.waves[result.count++] = prepared;
    }
  }
  return result;
}

} // namespace Pullback::SphereAnimation
