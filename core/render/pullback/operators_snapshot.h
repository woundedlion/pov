/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file operators_snapshot.h
 * @brief Chain-operator runtime state snapshot codecs.
 */

#include "platform/build_features.h"
#if HS_ENABLE_CHAIN_INTERPRETER

#include "render/pullback/operators.h"

namespace Pullback::Interp {

namespace Detail {
inline bool snapshot_phase(float value) {
  return std::isfinite(value) && value >= 0.0f && value < 1.0f;
}
inline bool snapshot_angle(float value) {
  return std::isfinite(value) && fabsf(value) < math::TWO_PI_F;
}
inline bool snapshot_unit(const math::Vector &value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z) &&
         fabsf(math::dot(value, value) - 1.0f) < 1e-4f;
}
inline bool snapshot_unit(const math::Quaternion &value) {
  return std::isfinite(value.r) && std::isfinite(value.v.x) &&
         std::isfinite(value.v.y) && std::isfinite(value.v.z) &&
         fabsf(math::dot(value, value) - 1.0f) < 1e-4f;
}
} // namespace Detail

template <> struct RuntimeStateCodec<Op::SpatialWalkState> {
  static RuntimeSnapshot capture(const Op::SpatialWalkState &state) {
    return SpatialWalkSnapshot{state.noise_seed, state.walk_time,
                               state.position,   state.direction,
                               state.wander,     state.angular_velocity,
                               state.spin_phase};
  }
  static bool restore(Op::SpatialWalkState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<SpatialWalkSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_unit(value->position) ||
        !Detail::snapshot_unit(value->direction) ||
        fabsf(math::dot(value->position, value->direction)) > 1e-3f ||
        !Detail::snapshot_unit(value->wander) ||
        !std::isfinite(value->angular_velocity) ||
        !Detail::snapshot_angle(value->spin_phase))
      return false;
    Op::init_walk(state, value->noise_seed);
    state.walk_time = value->walk_time;
    state.position = value->position;
    state.direction = value->direction;
    state.wander = value->wander;
    state.angular_velocity = value->angular_velocity;
    state.spin_phase = value->spin_phase;
    return true;
  }
};

template <> struct RuntimeStateCodec<Op::SourceClockState> {
  static RuntimeSnapshot capture(const Op::SourceClockState &state) {
    return SourceClockSnapshot{state.primary, state.secondary, state.angle};
  }
  static bool restore(Op::SourceClockState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<SourceClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_angle(value->primary) ||
        !Detail::snapshot_angle(value->secondary) ||
        !Detail::snapshot_angle(value->angle))
      return false;
    state = {value->primary, value->secondary, value->angle};
    return true;
  }
};

template <> struct RuntimeStateCodec<Op::NoisePhaseState> {
  static RuntimeSnapshot capture(const Op::NoisePhaseState &state) {
    return NoiseClockSnapshot{state.phase, state.noise_seed};
  }
  static bool restore(Op::NoisePhaseState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<NoiseClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_phase(value->phase))
      return false;
    state.noise_seed = value->noise_seed;
    init_effect_noise(state.noise, state.noise_seed);
    state.phase = value->phase;
    return true;
  }
};

template <> struct RuntimeStateCodec<Op::WarpPhaseState> {
  static RuntimeSnapshot capture(const Op::WarpPhaseState &state) {
    return PhaseClockSnapshot{state.phase};
  }
  static bool restore(Op::WarpPhaseState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<PhaseClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_phase(value->phase))
      return false;
    state.phase = value->phase;
    return true;
  }
};

template <> struct RuntimeStateCodec<Op::RipplePhaseState> {
  static RuntimeSnapshot capture(const Op::RipplePhaseState &state) {
    return RippleClockSnapshot{state.phase};
  }
  static bool restore(Op::RipplePhaseState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<RippleClockSnapshot>(&snapshot);
    if (value == nullptr || !std::isfinite(value->phase) || value->phase < 0.0f)
      return false;
    state.phase = value->phase;
    return true;
  }
};

template <> struct RuntimeStateCodec<Op::AffineClockState> {
  static RuntimeSnapshot capture(const Op::AffineClockState &state) {
    return AffineClockSnapshot{state.phase, state.rotation};
  }
  static bool restore(Op::AffineClockState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<AffineClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_phase(value->phase) ||
        !Detail::snapshot_angle(value->rotation))
      return false;
    state = {value->phase, value->rotation};
    return true;
  }
};

template <> struct RuntimeStateCodec<Op::ColorClockState> {
  static RuntimeSnapshot capture(const Op::ColorClockState &state) {
    return ColorClockSnapshot{state.oscillation_phase, state.hue_noise_phase,
                              state.hue_noise_seed};
  }
  static bool restore(Op::ColorClockState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<ColorClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_phase(value->oscillation_phase) ||
        !Detail::snapshot_phase(value->hue_noise_phase))
      return false;
    state.oscillation_phase = value->oscillation_phase;
    state.hue_noise_phase = value->hue_noise_phase;
    state.hue_noise_seed = value->hue_noise_seed;
    return true;
  }
};

template <> struct RuntimeStateCodec<Op::SphericalRingsState> {
  static RuntimeSnapshot capture(const Op::SphericalRingsState &state) {
    return SphericalRingsSnapshot{
        std::get<SpatialWalkSnapshot>(
            RuntimeStateCodec<Op::SpatialWalkState>::capture(state.walk)),
        state.phase};
  }
  static bool restore(Op::SphericalRingsState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<SphericalRingsSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_angle(value->phase) ||
        !RuntimeStateCodec<Op::SpatialWalkState>::restore(state.walk,
                                                          value->walk))
      return false;
    state.phase = value->phase;
    return true;
  }
};

} // namespace Pullback::Interp

#endif
