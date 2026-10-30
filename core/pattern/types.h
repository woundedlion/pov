/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file types.h
 * @brief Passive pattern coordinates, recipes and source state. */

#include <cmath>
#include <bit>
#include <cstdint>

#include "core/math/3dmath.h"
#include "core/math/noise_basis.h"

namespace Pattern {

/** @brief Tests IEEE-754 finiteness independently of floating-point mode.
 * @param value Candidate scalar.
 * @return Whether the scalar has a finite exponent. */
inline bool finite(float value) {
  return (std::bit_cast<uint32_t>(value) & 0x7f800000u) != 0x7f800000u;
}

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
    assert(valid(vector));
    assert(std::fabs(math::dot(vector, vector) - 1.0f) <=
           math::EPS_UNIT_VEC_SQ);
    return SphereDirection(vector);
  }
  /** @brief Normalizes an admitted finite, nondegenerate vector.
   * @param vector Normalizable direction.
   * @return Unit direction. */
  static SphereDirection admit(const math::Vector &vector) {
    assert(valid(vector));
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
  SPHERE_SMOOTH_GRID
};
/** @brief Explicit sphere evaluation strategy. */
enum class Evaluation : uint8_t { DIRECT, BAKED_CUBE };
/** @brief Noise octave formula. */
using NoiseBasis = math::NoiseBasis;
/** @brief Configured generator sampling or explicit octave/contrast sampling. */
enum class NoiseVariant : uint8_t { CONFIGURED_SINGLE, OCTAVE_CONTOUR };
/** @brief Maximum noise-loop speed exposed by modulation controls, in turns per update. */
inline constexpr float MAX_NOISE_SPEED = 0.001f;
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

/** @brief Passive parameters for `PlaneRings`. */
struct RingsSourceParams {
  float pattern_freq = 1.0f; /**< Radial pattern frequency. */
  float speed = 0.0f;        /**< Phase advance per frame, in radians. */
};

/** @brief Passive parameters for `PlaneCoupledGrid`. */
struct GridSourceParams {
  float pattern_freq = 1.0f; /**< Scale applied to the warped plane
                                   coordinates before the grid is sampled. */
  float speed = 0.0f;        /**< Per-frame advance of the primary phase. */
  float complexity = 0.0f;   /**< Amount of cross-axis coupling folded into the
                                   grid coordinates. */
  float pattern_mix = 0.0f;  /**< Blend from the coupled pattern at 0 to the
                                   direct sine product at 1. */
  float secondary_rate = 0.0f; /**< Secondary phase rate, as a multiple of
                                    `speed`. */
  float angle_rate = 0.0f;     /**< Per-frame advance of the source rotation. */
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

/** @brief Passive parameters for `PlaneSpiral`. */
struct SpiralSourceParams {
  float pattern_freq = 1.0f; /**< Plane-coordinate scale before sampling. */
  float speed = 0.0f;        /**< Per-frame advance of the primary phase. */
  float angle_rate = 0.0f;   /**< Per-frame advance of the spiral rotation. */
};

/** @brief Passive parameters for `PlaneNoise`. */
struct NoiseSourceParams {
  float noise_scale = 1.0f;    /**< Spatial scale of the sampled field. */
  float noise_contrast = 0.0f; /**< Contour sharpening applied to the sample. */
  float noise_time_rate = 0.0f; /**< Per-frame advance of the noise time
                                     coordinate. */
};

/** @brief Passive parameters for the lattice kernel. */
struct LatticeSourceParams {
  float lattice_cell_scale = 1.0f;  /**< Lattice cells per plane unit. */
  float lattice_shape_blend = 0.0f; /**< Cell primitive, from a circle at 0 to a
                                         square at 1. */
  float lattice_softness = 0.05f;   /**< Half-width of the ramp across the
                                         primitive's boundary. */
  float lattice_radius = 0.25f;     /**< Primitive radius in cell units. */
};

/** @brief Passive parameters for `SphereRings`. */
struct SphericalRingsSourceParams {
  float ring_count = 6.0f;      /**< Number of bands from pole to pole. */
  float ring_thickness = 0.08f; /**< Band half-width, in radians. */
  float ring_softness = 0.02f;  /**< Angular width of the antialiased edge. */
  float speed = 0.0f;           /**< Per-frame advance of the band phase. */
  float spin_rate = 0.0f;       /**< Per-frame rotation of the band axis. */
  float wander = 0.0f;          /**< Fraction of the random walk applied. */
};

/** @brief Passive parameters for `PlaneEscapeFractal`. */
struct FractalSourceParams {
  float scale = 0.5f;      /**< Plane-coordinate scale before iteration. */
  float iterations = 8.0f; /**< Escape iterations in [2,16]. */
  float julia_mix = 0.0f;  /**< Blend from Mandelbrot at 0 to Julia at 1. */
  float julia_re = -0.8f;  /**< Julia seed real component. */
  float julia_im = 0.156f; /**< Julia seed imaginary component. */
  float contours = 4.0f;   /**< Exterior contour cycles across the orbit. */
  float speed = 0.0f;      /**< Per-frame rotation of the Julia seed. */
  float angle_rate = 0.0f; /**< Per-frame rotation of the source plane. */
};

/** @brief Passive parameters for `PlaneTessellation`. */
struct TessellationSourceParams {
  float cell_scale = 1.0f;      /**< Tessellation cells per plane unit. */
  float line_thickness = 0.04f; /**< Edge half-width in cell units. */
  float line_softness = 0.02f;  /**< Width of the antialiased edge. */
  float angle_rate = 0.0f;      /**< Per-frame rotation of the tessellation. */
};

/** @brief Cell shape of a `Tessellation` source. */
enum class TessellationKind : uint8_t {
  TRIANGULAR = 0,
  SQUARE = 1,
  HEXAGONAL = 2
};

} // namespace Pattern
