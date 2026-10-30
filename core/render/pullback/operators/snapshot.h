/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file snapshot.h
 * @brief Chain-operator runtime state snapshot codecs.
 */

#include "core/platform/build_features.h"
#if HS_ENABLE_CHAIN_INTERPRETER

#include "core/render/pullback/operators.h"

namespace Pullback::Interp {

namespace Detail {
/**
 * @brief Whether @p value is a finite normalized phase in [0, 1).
 * @param value Phase to check.
 * @return True when admissible.
 */
inline bool snapshot_phase(float value) {
  return Pattern::finite(value) && value >= 0.0f && value < 1.0f;
}
/**
 * @brief Whether @p value is a finite angle with magnitude below 2π.
 * @param value Angle in radians.
 * @return True when admissible.
 */
inline bool snapshot_angle(float value) {
  return Pattern::finite(value) && fabsf(value) < math::TWO_PI_F;
}
/**
 * @brief Whether @p value is a finite unit vector.
 * @param value Vector to check.
 * @return True when admissible.
 */
inline bool snapshot_unit(const math::Vector &value) {
  return Pattern::finite(value.x) && Pattern::finite(value.y) &&
         Pattern::finite(value.z) &&
         fabsf(math::dot(value, value) - 1.0f) < 1e-4f;
}
/**
 * @brief Whether @p value is a finite unit quaternion.
 * @param value Quaternion to check.
 * @return True when admissible.
 */
inline bool snapshot_unit(const math::Quaternion &value) {
  return Pattern::finite(value.r) && Pattern::finite(value.v.x) &&
         Pattern::finite(value.v.y) && Pattern::finite(value.v.z) &&
         fabsf(math::dot(value, value) - 1.0f) < 1e-4f;
}
} // namespace Detail

/** @brief Snapshot codec of Op::SpatialWalkState. */
template <> struct RuntimeStateCodec<Op::SpatialWalkState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A SpatialWalkSnapshot.
   */
  static RuntimeSnapshot capture(const Op::SpatialWalkState &state) {
    return SpatialWalkSnapshot{state.noise_seed, state.walk_time,
                               state.position,   state.direction,
                               state.wander,     state.angular_velocity,
                               state.spin_phase};
  }
  /**
   * @brief Restores @p state from a SpatialWalkSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
  static bool restore(Op::SpatialWalkState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<SpatialWalkSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_unit(value->position) ||
        !Detail::snapshot_unit(value->direction) ||
        fabsf(math::dot(value->position, value->direction)) > 1e-3f ||
        !Detail::snapshot_unit(value->wander) ||
        !Pattern::finite(value->angular_velocity) ||
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

/** @brief Snapshot codec of Op::SourceClockState. */
template <> struct RuntimeStateCodec<Op::SourceClockState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A SourceClockSnapshot.
   */
  static RuntimeSnapshot capture(const Op::SourceClockState &state) {
    return SourceClockSnapshot{state.primary, state.secondary, state.angle};
  }
  /**
   * @brief Restores @p state from a SourceClockSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
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

/** @brief Snapshot codec of Op::NoisePhaseState. */
template <> struct RuntimeStateCodec<Op::NoisePhaseState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A NoiseClockSnapshot.
   */
  static RuntimeSnapshot capture(const Op::NoisePhaseState &state) {
    return NoiseClockSnapshot{state.phase, state.noise_seed};
  }
  /**
   * @brief Restores @p state from a NoiseClockSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
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

/** @brief Snapshot codec of Op::WarpPhaseState. */
template <> struct RuntimeStateCodec<Op::WarpPhaseState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A PhaseClockSnapshot.
   */
  static RuntimeSnapshot capture(const Op::WarpPhaseState &state) {
    return PhaseClockSnapshot{state.phase};
  }
  /**
   * @brief Restores @p state from a PhaseClockSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
  static bool restore(Op::WarpPhaseState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<PhaseClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_phase(value->phase))
      return false;
    state.phase = value->phase;
    return true;
  }
};

/** @brief Snapshot codec of Op::RipplePhaseState. */
template <> struct RuntimeStateCodec<Op::RipplePhaseState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A RippleClockSnapshot.
   */
  static RuntimeSnapshot capture(const Op::RipplePhaseState &state) {
    return static_cast<const RippleClockSnapshot &>(state);
  }
  /**
   * @brief Restores @p state from a RippleClockSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
  static bool restore(Op::RipplePhaseState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<RippleClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_phase(value->phase))
      return false;
    for (const auto &wave : value->waves)
      if (!Pattern::finite(wave.age) || wave.age < -49.0f || wave.age > 1000.0f)
        return false;
    static_cast<RippleClockSnapshot &>(state) = *value;
    return true;
  }
};

/** @brief Snapshot codec of Op::AffineClockState. */
template <> struct RuntimeStateCodec<Op::AffineClockState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A AffineClockSnapshot.
   */
  static RuntimeSnapshot capture(const Op::AffineClockState &state) {
    return AffineClockSnapshot{state.phase, state.rotation};
  }
  /**
   * @brief Restores @p state from a AffineClockSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
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

/** @brief Snapshot codec of Op::ColorClockState. */
template <> struct RuntimeStateCodec<Op::ColorClockState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A ColorClockSnapshot.
   */
  static RuntimeSnapshot capture(const Op::ColorClockState &state) {
    return ColorClockSnapshot{state.oscillation_phase};
  }
  /**
   * @brief Restores @p state from a ColorClockSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
  static bool restore(Op::ColorClockState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<ColorClockSnapshot>(&snapshot);
    if (value == nullptr || !Detail::snapshot_phase(value->oscillation_phase))
      return false;
    state.oscillation_phase = value->oscillation_phase;
    return true;
  }
};

/** @brief Snapshot codec of one modulation clock. */
template <> struct RuntimeStateCodec<Op::ModulationClockState> {
  /** @brief Captures the source clock.
   * @param state Clock state.
   * @return Captured clock. */
  static RuntimeSnapshot capture(const Op::ModulationClockState &state) {
    return ModulationClockSnapshot{state.modulation};
  }
  /** @brief Restores a validated source clock.
   * @param state State to update.
   * @param snapshot Captured clock.
   * @return Whether the clock is valid. */
  static bool restore(Op::ModulationClockState &state,
                      const RuntimeSnapshot &snapshot) {
    const auto *value = std::get_if<ModulationClockSnapshot>(&snapshot);
    if (!value ||
        !Detail::snapshot_phase(value->modulation.phase.phase_turns) ||
        (!value->modulation.initialized &&
         value->modulation.phase.phase_turns != 0.0f) ||
        value->modulation.version != 1 ||
        static_cast<unsigned>(value->modulation.source) >
            static_cast<unsigned>(Pattern::SourceKind::DISPLACEMENT))
      return false;
    state.modulation = value->modulation;
    return true;
  }
};

/** @brief Snapshot codec of Op::SphericalRingsState. */
template <> struct RuntimeStateCodec<Op::SphericalRingsState> {
  /**
   * @brief Captures @p state.
   * @param state State to capture.
   * @return A SphericalRingsSnapshot.
   */
  static RuntimeSnapshot capture(const Op::SphericalRingsState &state) {
    return SphericalRingsSnapshot{
        std::get<SpatialWalkSnapshot>(
            RuntimeStateCodec<Op::SpatialWalkState>::capture(state.walk)),
        state.phase};
  }
  /**
   * @brief Restores @p state from a SphericalRingsSnapshot.
   * @param state State to overwrite.
   * @param snapshot Captured snapshot.
   * @return False, leaving @p state untouched, on a wrong kind or
   *         an out-of-range value.
   */
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
