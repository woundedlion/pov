/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <variant>
#include "math/3dmath.h"

namespace Pullback::Interp {

struct SpatialWalkSnapshot {
  int32_t noise_seed = 0;
  uint32_t walk_time = 0;
  math::Vector position = math::UP;
  math::Vector direction = math::perpendicular_axis(math::UP);
  math::Quaternion wander;
  math::Quaternion raw_orientation;
  float angular_velocity = 0;
  float spin_phase = 0;
  bool legacy = false;
};
struct SourceClockSnapshot {
  float primary = 0;
  float secondary = 0;
  float angle = 0;
};
struct NoiseClockSnapshot {
  float phase = 0;
  int32_t noise_seed = 0;
};
struct PhaseClockSnapshot {
  float phase = 0;
};
struct RippleClockSnapshot {
  float phase = 0;
};
struct AffineClockSnapshot {
  float phase = 0;
  float rotation = 0;
};
struct ColorClockSnapshot {
  float oscillation_phase = 0;
  float hue_noise_phase = 0;
  int32_t hue_noise_seed = 6047;
};
struct SphericalRingsSnapshot {
  SpatialWalkSnapshot walk;
  float phase = 0;
};

using RuntimeSnapshot =
    std::variant<std::monostate, SpatialWalkSnapshot, SourceClockSnapshot,
                 NoiseClockSnapshot, PhaseClockSnapshot, RippleClockSnapshot,
                 AffineClockSnapshot, ColorClockSnapshot,
                 SphericalRingsSnapshot>;

inline const char *runtime_snapshot_kind(const RuntimeSnapshot &snapshot) {
  constexpr const char *KINDS[] = {"none",
                                   "spatial-walk-v1",
                                   "source-clock-v1",
                                   "noise-clock-v1",
                                   "phase-clock-v1",
                                   "ripple-clock-v1",
                                   "affine-clock-v1",
                                   "color-clock-v1",
                                   "spherical-rings-v1"};
  static_assert(std::size(KINDS) == std::variant_size_v<RuntimeSnapshot>);
  return KINDS[snapshot.index()];
}

template <typename State> struct RuntimeStateCodec {
  static RuntimeSnapshot capture(const State &) { return {}; }
  static bool restore(State &, const RuntimeSnapshot &snapshot) {
    return std::holds_alternative<std::monostate>(snapshot);
  }
};

} // namespace Pullback::Interp
