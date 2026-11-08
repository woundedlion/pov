/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file types.h
 * @brief Passive pattern coordinates, recipes and source state. */

#include <cmath>
#include <cstdint>

#include "core/math/3dmath.h"
#include "core/math/noise_basis.h"

namespace Pattern {

using hs::finite;

/** @brief One-dimensional pattern coordinate. */
struct ScalarCoordinate {
  float value; ///< Scalar source coordinate.
};
/** @brief Plane point without a projection binding. */
struct PlanePosition {
  math::Complex value; ///< Plane position.
};
/** @brief Volumetric point whose radius is retained. */
struct Position3D {
  math::Vector value; ///< Volumetric position.
};

/** @brief Unit direction admitted from a finite, nondegenerate vector. */
class SphereDirection {
  math::Vector direction = math::Vector(0.0f, 1.0f, 0.0f);
  explicit SphereDirection(math::Vector direction) : direction(direction) {}

public:
  SphereDirection() = default;
  /** @brief Tests the squared-length admission boundary.
   * @param vector Candidate direction.
   * @return Whether its computed squared length is finite and normalizable. */
  static bool valid(const math::Vector &vector) {
    const float squared = math::dot(vector, vector);
    return finite(squared) && squared >= math::EPS_NORMALIZE_SQ;
  }
  /** @brief Admits and normalizes a direction.
   * @param vector Candidate direction.
   * @return Unit direction. */
  static SphereDirection from(const math::Vector &vector) {
    HS_CHECK(valid(vector), "pattern: invalid sphere direction");
    return SphereDirection(vector.normalized());
  }
  /** @brief Admits a unit vector without normalization.
   * @param vector Unit direction.
   * @return Direction with unchanged components. */
  static SphereDirection from_unit(const math::Vector &vector) {
    HS_AUDIT_CHECK(valid(vector));
    HS_AUDIT_CHECK(std::fabs(math::dot(vector, vector) - 1.0f) <=
                   math::EPS_UNIT_VEC_SQ);
    return SphereDirection(vector);
  }
  /** @brief Normalizes an admitted finite, nondegenerate vector.
   * @param vector Normalizable direction.
   * @return Unit direction. */
  static SphereDirection admit(const math::Vector &vector) {
    HS_AUDIT_CHECK(valid(vector));
    const float M2 =
        vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
    const float M = sqrtf(M2);
    return SphereDirection(
        math::Vector(vector.x / M, vector.y / M, vector.z / M));
  }
  /** @brief Reads the admitted unit direction.
   * @return Unit vector. */
  const math::Vector &value() const { return direction; }
};

/** @brief Serializable built-in source selection. */
enum class SourceKind : uint8_t {
  NOISE,
  DIRECTIONAL_SINE,
  SPHERE_INTERFERENCE,
  SPHERE_SMOOTH_GRID,
  DISPLACEMENT
};
/** @brief Explicit sphere evaluation strategy. */
enum class Evaluation : uint8_t { DIRECT, BAKED_CUBE };
/** @brief Noise octave formula. */
using NoiseBasis = math::NoiseBasis;
/** @brief Configured generator sampling or explicit octave/contrast sampling. */
enum class NoiseVariant : uint8_t { CONFIGURED_SINGLE, OCTAVE_CONTOUR };
/** @brief Maximum noise-loop speed exposed by modulation controls, in turns per update. */
inline constexpr float MAX_NOISE_SPEED = 0.001f;
/** @brief Maximum spatial noise scale admitted by a baked sphere cube. */
inline constexpr float MAX_BAKED_NOISE_SCALE = 8.0f;
/// Minimum scale of a direct noise source.
inline constexpr float MIN_NOISE_SCALE = 1.0f / 64.0f;
/// Maximum scale of a direct noise source.
inline constexpr float MAX_NOISE_SCALE = 64.0f;
/// Maximum contour contrast of a noise source.
inline constexpr float MAX_NOISE_CONTRAST = 8.0f;
/** @brief Maximum analytic speed exposed by modulation controls, in turns per update. */
inline constexpr float MAX_ANALYTIC_SPEED = 1.0f;
/** @brief Source phase in turns, wrapped after each engine update. */
struct PhaseState {
  float phase_turns = 0.0f; ///< Wrapped loop phase in turns.
  constexpr bool operator==(const PhaseState &) const = default;
};
/** @brief Source clocks in radians, each wrapped at two pi. */
struct RadianState {
  float primary = 0.0f;   ///< Primary source phase.
  float secondary = 0.0f; ///< Rate-coupled secondary phase.
  float angle = 0.0f;     ///< Independent source rotation.
  constexpr bool operator==(const RadianState &) const = default;
};
/** @brief Empty passive state of a stationary source. */
struct EmptyState {
  constexpr bool operator==(const EmptyState &) const = default;
};
/** @brief Parameters of an axial sine. */
struct SineParams {
  math::Vector axis = math::Vector(0.0f, 1.0f, 0.0f); ///< Wave axis.
  float frequency = 1.0f; ///< Cycles per unit coordinate.
  float offset = 0.0f;    ///< Phase offset in turns.
  float speed = 0.0f;     ///< Phase advance in turns per update.
  /** @brief Compares the complete wave recipe.
   * @param other Recipe to compare.
   * @return Whether every control matches. */
  constexpr bool operator==(const SineParams &other) const {
    return axis.x == other.axis.x && axis.y == other.axis.y &&
           axis.z == other.axis.z && frequency == other.frequency &&
           offset == other.offset && speed == other.speed;
  }
};
/** @brief Parameters of two axial sines sharing one loop phase. */
struct InterferenceParams : SineParams {
  math::Vector axis_b = math::Vector(1.0f, 0.0f, 0.0f); ///< Second wave axis.
  float frequency_b = 1.0f; ///< Second wave cycles per unit coordinate.
  float offset_b = 0.0f;    ///< Second wave phase offset in turns.
  /** @brief Compares both wave recipes.
   * @param other Recipe to compare.
   * @return Whether every control matches. */
  constexpr bool operator==(const InterferenceParams &other) const {
    return static_cast<const SineParams &>(*this) ==
               static_cast<const SineParams &>(other) &&
           axis_b.x == other.axis_b.x && axis_b.y == other.axis_b.y &&
           axis_b.z == other.axis_b.z && frequency_b == other.frequency_b &&
           offset_b == other.offset_b;
  }
};
/** @brief Parameters of a smooth sine grid with independently oriented axes. */
struct GridParams : InterferenceParams {};
/** @brief Parameters of a sphere noise source. */
struct NoiseParams {
  float scale = 1.0f;                     ///< Noise-space coordinate scale.
  float speed = 0.0f;                     ///< Loop advance in turns per update.
  NoiseBasis basis = NoiseBasis::SIMPLEX; ///< Octave formula.
  float contrast = 0.0f;                  ///< Nonnegative contour sharpening.
  NoiseVariant variant = NoiseVariant::CONFIGURED_SINGLE; ///< Sampling policy.
  constexpr bool operator==(const NoiseParams &) const = default;
};

/** @brief Passive parameters for `PlaneCoupledGrid`. */
struct GridSourceParams {
  float pattern_freq = 1.0f; /**< Scale applied to the warped plane
                                   coordinates before the grid is sampled. */
  float speed = 0.0f;        /**< Per-frame advance of the coupling phase. */
  float complexity = 0.0f;   /**< Amount of cross-axis coupling folded into the
                                   grid coordinates. */
  float pattern_mix = 0.0f;  /**< Blend from the coupled pattern at 0 to the
                                   direct sine product at 1. */
};

/** @brief Passive parameters for `PlaneInterference`. */
struct TwinWaveSourceParams {
  float pattern_freq = 1.0f;   /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;          /**< Per-frame advance of the primary phase. */
  float secondary_rate = 0.0f; /**< Secondary phase rate, as a multiple of
                                    `speed`. */
  float angle_rate = 0.0f; /**< Per-frame advance of the angle between the two
                                waves. */
};

/** @brief Passive parameters for `PlaneBands`. */
struct BandsSourceParams {
  float pattern_freq = 1.0f; ///< Radians per plane-coordinate unit.
  float angle = 0.0f; ///< Band angle counterclockwise from horizontal, radians.
};

/** @brief Coordinate chart consumed by the spiral source. */
enum class SpiralCoordinateSystem : uint8_t { PLANE, EQUIRECTANGULAR };

/** @brief Passive parameters for `PlaneSpiral`. */
struct SpiralSourceParams {
  float pattern_freq = 1.0f; /**< Plane-coordinate scale before sampling. */
  float arms = 3.0f; /**< Spiral arms, rounded to whole; 0 draws concentric
                          rings. */
  float center_azimuth = 0.0f; ///< Spherical winding-center azimuth, radians.
  float center_polar = 0.0f; ///< Winding-center colatitude from north, radians.
  uint8_t coordinate_system = static_cast<uint8_t>(
      SpiralCoordinateSystem::
          PLANE); ///< Cartesian plane or equirectangular chart.
};

/** @brief Passive parameters for `PlaneNoise`. */
struct NoiseSourceParams {
  float noise_scale = 1.0f;    /**< Spatial scale of the sampled field. */
  float noise_contrast = 0.0f; /**< Contour sharpening applied to the sample. */
  float noise_time_rate = 0.0f; /**< Per-frame advance of the noise time
                                     coordinate. */
};

/** @brief Cell shape of a `PlaneTessellation`. */
enum class TessellationKind : uint8_t {
  TRIANGULAR = 0,
  SQUARE = 1,
  HEXAGONAL = 2
};

/** @brief What a `PlaneTessellation` draws in each cell. */
enum class TessellationFill : uint8_t {
  CELLS = 0, ///< A primitive centred on each cell.
  EDGES = 1  ///< A band along the cell edges.
};

/** @brief Cell-primitive controls of a tessellation.
 * @details Lengths are in cell units: one cell inscribed diameter, which is
 * one lattice period for squares and hexagons and period / sqrt(3) for
 * triangles. */
struct TessellationCellParams {
  float cell_scale = 1.0f;     /**< Lattice periods per plane unit. */
  float cell_shape = 0.0f;     /**< Primitive, from a circle at 0 to the cell
                                    polygon at 1. */
  float cell_softness = 0.05f; /**< Half-width of the ramp across the
                                    primitive's boundary. */
  float cell_radius = 0.25f;   /**< Primitive inradius. */
};

/** @brief Edge-band controls of a tessellation, in cell units. */
struct TessellationEdgeParams {
  float cell_scale = 1.0f;      /**< Lattice periods per plane unit. */
  float line_thickness = 0.04f; /**< Edge half-width. */
  float line_softness = 0.02f;  /**< Width of the antialiased edge. */
};

/** @brief Passive parameters for `PlaneTessellation`: the controls of both
    fills. */
struct TessellationSourceParams : TessellationCellParams {
  float line_thickness = 0.04f; /**< Edge half-width in cell units. */
  float line_softness = 0.02f;  /**< Width of the antialiased edge. */
};

} // namespace Pattern
