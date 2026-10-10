/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file runtime_snapshot.h
 * @brief Captured operator runtime state and restore codecs. */

#include <variant>
#include "render/pullback/runtime_seeds.h"
#include "math/3dmath.h"

namespace Pullback::Interp {

/** @brief Captured random-walk camera state of a spatial-walk operator. */
struct SpatialWalkSnapshot {
  int32_t noise_seed = 0;           ///< Seed of the walk noise field.
  uint32_t walk_time = 0;           ///< Walk steps taken, in frames.
  math::Vector position = math::UP; ///< Unit walker position on the sphere.
  /** Unit heading, perpendicular to `position`. */
  math::Vector direction = math::perpendicular_axis(math::UP);
  math::Quaternion wander;    ///< Unit accumulated wander rotation.
  float angular_velocity = 0; ///< Heading pivot rate, radians per frame.
  float spin_phase = 0;       ///< Spin angle in radians, |value| < 2 pi.
};
/** @brief Captured source clocks, each an angle in radians, |value| < 2 pi. */
struct SourceClockSnapshot {
  float primary = 0;   ///< Primary phase clock.
  float secondary = 0; ///< Secondary phase clock.
  float angle = 0;     ///< Pattern rotation angle.
};
/** @brief Captured clock and seed of a noise-driven operator. */
struct NoiseClockSnapshot {
  float phase = 0;        ///< Loop phase in [0, 1).
  int32_t noise_seed = 0; ///< Seed the noise field is rebuilt from.
};
/** @brief Captured warp phase clock. */
struct PhaseClockSnapshot {
  float phase = 0; ///< Loop phase in [0, 1).
};
/** @brief Captured ripple clock. */
struct RippleClockSnapshot {
  float phase = 0; ///< Frames into the current ripple period, non-negative.
};
/** @brief Captured affine warp clock. */
struct AffineClockSnapshot {
  float phase = 0; ///< Loop phase in [0, 1).
  float rotation =
      0; ///< Accumulated frame rotation in radians, |value| < 2 pi.
};
/** @brief Captured colour clocks of a generated-palette operator. */
struct ColorClockSnapshot {
  float oscillation_phase = 0;             ///< Mapping oscillation, [0, 1).
  float hue_noise_phase = 0;               ///< Hue-noise loop phase, [0, 1).
  int32_t hue_noise_seed = HUE_NOISE_SEED; ///< Seed of the hue-noise field.
};
/** @brief Captured state of the spherical-rings sampler. */
struct SphericalRingsSnapshot {
  SpatialWalkSnapshot walk; ///< Ring-axis walk.
  float phase = 0;          ///< Ring phase in radians, |value| < 2 pi.
};

/** @brief Any operator's captured runtime state; `std::monostate` for
    stateless operators. */
using RuntimeSnapshot =
    std::variant<std::monostate, SpatialWalkSnapshot, SourceClockSnapshot,
                 NoiseClockSnapshot, PhaseClockSnapshot, RippleClockSnapshot,
                 AffineClockSnapshot, ColorClockSnapshot,
                 SphericalRingsSnapshot>;

/**
 * @brief Versioned wire name of the snapshot's alternative.
 * @param snapshot Snapshot to name.
 * @return Static kind string, "none" for `std::monostate`.
 */
inline const char *runtime_snapshot_kind(const RuntimeSnapshot &snapshot) {
  constexpr const char *KINDS[] = {"none",
                                   "spatial-walk-v2",
                                   "source-clock-v1",
                                   "noise-clock-v1",
                                   "phase-clock-v1",
                                   "ripple-clock-v1",
                                   "affine-clock-v1",
                                   "color-clock-v1",
                                   "spherical-rings-v2"};
  static_assert(std::size(KINDS) == std::variant_size_v<RuntimeSnapshot>);
  return KINDS[snapshot.index()];
}

/** @brief Stateful operators specialize this codec; the primary captures monostate. */
template <typename State> struct RuntimeStateCodec {
  /** @brief Captures nothing.
      @return An empty (`std::monostate`) snapshot. */
  static RuntimeSnapshot capture(const State &) { return {}; }
  /** @brief Accepts only an empty snapshot; the state is untouched.
      @param snapshot Snapshot to restore from.
      @return Whether @p snapshot holds `std::monostate`. */
  static bool restore(State &, const RuntimeSnapshot &snapshot) {
    return std::holds_alternative<std::monostate>(snapshot);
  }
};

} // namespace Pullback::Interp
