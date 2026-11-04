/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>

#include "core/math/noise_field.h"

#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"
#include "core/render/pullback/translation_lattice.h"
#include "core/math/3dmath.h"
#include "core/math/projection_patterns.h"
#include <iterator>
#include <limits>

/**
 * @file warp.h
 * @brief Planar warp policies.
 */

namespace Pullback {

namespace Warp {

/** @brief Continuous parameters of the planar curl flow. */
struct CurlFlowParams {
  float speed = 0.0f;    /**< Loop phase advance per frame. */
  float strength = 0.0f; /**< Signed curl displacement strength. */
  float scale = 1.0f;    /**< Spatial frequency of the curl field. */
  /** @brief Per-parameter id, member, label, range, curve and gate. */
  static constexpr auto FIELDS = std::array{
      Field<CurlFlowParams>{"speed", &CurlFlowParams::speed, nullptr, -0.02f,
                            0.02f, FieldCurve::LERP},
      Field<CurlFlowParams>{"strength", &CurlFlowParams::strength,
                            "Warp Strength", -1.0f / 32.0f, 1.0f / 32.0f,
                            FieldCurve::LERP},
      Field<CurlFlowParams>{"scale", &CurlFlowParams::scale, "Warp Scale",
                            1.0f / 64.0f, 4.0f, FieldCurve::LOG_POSITIVE},
  };
};
static_assert(field_ids_unique<CurlFlowParams>());
static_assert(field_defaults_in_range<CurlFlowParams>());

/** @brief Largest angular harmonic a `PolarChart` accepts. */
inline constexpr uint8_t MAX_POLAR_HARMONIC = 16;

/** @brief Fixed envelope policy: unit amplitude everywhere. */
struct FlatEnvelope {};
/** @brief Fixed envelope policy: amplitude follows the projection weight. */
struct ProjectionWeightEnvelope {};
/** @brief Fixed envelope policy: amplitude fades toward the projection edge. */
struct EdgeFadeEnvelope {};

/** @brief Envelope shaping a warp's amplitude across the projected domain. */
enum class Envelope : uint8_t {
  FLAT = 0,
  PROJECTION_WEIGHT = 1,
  EDGE_FADE = 2
};

/** @brief Serialized `Envelope` ids, indexed by enumerator value. */
inline constexpr const char *ENVELOPE_IDS[] = {"flat", "projection-weight",
                                               "edge-fade"};
static_assert(std::size(ENVELOPE_IDS) ==
              static_cast<size_t>(Envelope::EDGE_FADE) + 1);

/** @brief Activation relation of the fade band width, which only the edge-fade
    envelope reads. */
inline constexpr TopologyGate ENVELOPE_EDGE_FADE_GATE{
    "envelope", live_values(Envelope::EDGE_FADE)};
/** @brief Curl-flow integrator taking one Euler step. */
struct Euler1 {
  static constexpr uint8_t INTERVALS = 1; ///< Integration intervals.
};
/** @brief Curl-flow integrator taking two midpoint steps. */
struct Midpoint2 {
  static constexpr uint8_t INTERVALS = 2; ///< Integration intervals.
};
/** @brief Curl-flow integrator taking four midpoint steps. */
struct Midpoint4 {
  static constexpr uint8_t INTERVALS = 4; ///< Integration intervals.
};
/** @brief `PolarChart` mode mapping radius linearly. */
struct LinearPolar {};
/** @brief `PolarChart` mode mapping radius logarithmically. */
struct LogarithmicPolar {};

/**
 * @brief Warp parameters for the mirrored tiling
 *        (Pullback::Warp::MirrorTile).
 */
struct MirrorParams {
  float speed = 0.0f;    /**< Per-frame advance of the slot's phase. */
  float rotation = 0.0f; /**< Rotation of the fold lattice, in radians. */
  float cell_x = 1.0f;   /**< Mirror cell width in plane units. */
  float cell_y = 1.0f;   /**< Mirror cell height in plane units. */
  float offset_x = 0.0f; /**< Pre-fold translation along x; scrolls with the
                              slot's phase. */
  float offset_y = 0.0f; /**< Pre-fold translation along y; does not scroll. */

  /** @brief Per-parameter id, member, label, range, curve and gate. */
  static constexpr auto FIELDS = std::array{
      Field<MirrorParams>{"speed", &MirrorParams::speed, nullptr, -0.02f, 0.02f,
                          FieldCurve::LERP},
      Field<MirrorParams>{"rotation", &MirrorParams::rotation,
                          "Mirror Rotation", 0.0f, math::TWO_PI_F,
                          FieldCurve::SHORTEST_PERIODIC},
      Field<MirrorParams>{"cell-x", &MirrorParams::cell_x, "Mirror Cell X",
                          1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<MirrorParams>{"cell-y", &MirrorParams::cell_y, "Mirror Cell Y",
                          1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<MirrorParams>{"offset-x", &MirrorParams::offset_x,
                          "Mirror Offset X", -8.0f, 8.0f, FieldCurve::LERP},
      Field<MirrorParams>{"offset-y", &MirrorParams::offset_y,
                          "Mirror Offset Y", -8.0f, 8.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<MirrorParams>());
static_assert(field_defaults_in_range<MirrorParams>());

/** @brief Warp parameters for the sine shear (Pullback::Warp::WaveShear). */
struct WaveShearParams {
  float speed = 0.0f;       /**< Per-frame advance of the slot's phase. */
  float strength = 0.0f;    /**< Shear amplitude; 0 skips the stage. */
  float frequency = 1.0f;   /**< Spatial frequency along the field axis. */
  float field_angle = 0.0f; /**< Field axis direction, in radians. */
  float edge_width = 0.1f;  /**< Fade band width, read only under an
                                 EdgeFadeEnvelope; 0 is a hard cut. */

  /** @brief Per-parameter id, member, label, range, curve and gate. */
  static constexpr auto FIELDS = std::array{
      Field<WaveShearParams>{"speed", &WaveShearParams::speed, nullptr, -0.02f,
                             0.02f, FieldCurve::LERP},
      Field<WaveShearParams>{"strength", &WaveShearParams::strength,
                             "Warp Strength", -30.0f, 30.0f, FieldCurve::LERP},
      Field<WaveShearParams>{"frequency", &WaveShearParams::frequency,
                             "Warp Frequency", 0.01f, 32.0f,
                             FieldCurve::LOG_POSITIVE},
      Field<WaveShearParams>{"field-angle", &WaveShearParams::field_angle,
                             "Warp Field Angle", 0.0f, math::TWO_PI_F,
                             FieldCurve::SHORTEST_PERIODIC},
      edge_width_field(&WaveShearParams::edge_width, nullptr,
                       ENVELOPE_EDGE_FADE_GATE),
  };
};
static_assert(field_ids_unique<WaveShearParams>());
static_assert(field_defaults_in_range<WaveShearParams>());

/**
 * @brief Warp parameters for the noise-vector displacement
 *        (Pullback::Warp::VectorNoise).
 */
struct VectorNoiseParams {
  float speed = 0.0f;        /**< Per-frame advance of the slot's phase, which
                                  walks the noise loop. */
  float strength = 0.0f;     /**< Displacement amplitude; 0 skips the stage. */
  float scale = 1.0f;        /**< Spatial scale of the sampled field. */
  float vector_angle = 0.0f; /**< Rotation applied to the sampled vector, in
                                  radians. */
  float edge_width = 0.1f;   /**< Fade band width, read only under an
                                  EdgeFadeEnvelope; 0 is a hard cut. */

  /** @brief Per-parameter id, member, label, range, curve and gate. */
  static constexpr auto FIELDS = std::array{
      Field<VectorNoiseParams>{"speed", &VectorNoiseParams::speed, nullptr,
                               -0.02f, 0.02f, FieldCurve::LERP},
      Field<VectorNoiseParams>{"strength", &VectorNoiseParams::strength,
                               "Warp Strength", -30.0f, 30.0f,
                               FieldCurve::LERP},
      Field<VectorNoiseParams>{"scale", &VectorNoiseParams::scale, "Warp Scale",
                               1.0f / 64.0f, 64.0f, FieldCurve::LOG_POSITIVE},
      Field<VectorNoiseParams>{"vector-angle", &VectorNoiseParams::vector_angle,
                               "Warp Vector Angle", 0.0f, math::TWO_PI_F,
                               FieldCurve::SHORTEST_PERIODIC},
      edge_width_field(&VectorNoiseParams::edge_width, nullptr,
                       ENVELOPE_EDGE_FADE_GATE),
  };
};
static_assert(field_ids_unique<VectorNoiseParams>());
static_assert(field_defaults_in_range<VectorNoiseParams>());

/**
 * @brief Warp parameters for the affine frame change
 *        (Pullback::Warp::AffineFrame).
 * @details Translation windings count whole basis vectors of the next stage's
 * TranslationLattice, or of a square `lattice_period` lattice when that stage
 * reports none, times the whole repeat; the phase wrap is then seamless.
 */
struct AffineParams {
  float speed = 0.0f;         /**< Per-frame advance of the slot's phase. */
  float rotation_rate = 0.0f; /**< Frame rotation per frame, scaled by speed
                               and accumulated per affine warp instance. */
  float translation_x = 0.0f; /**< Repeats along the lattice's first basis
                                   vector per phase turn, rounded to whole. */
  float translation_y = 0.0f; /**< Repeats along the lattice's second basis
                                   vector per phase turn, rounded to whole. */
  float scale_x = 1.0f; /**< Scale along x, oscillated over the phase cycle. */
  float scale_y = 1.0f; /**< Scale along y, oscillated over the phase cycle. */
  float shear = 0.0f;   /**< Shear, oscillated over the phase cycle. */
  float repeat_cells = 1.0f;   /**< Lattice cells per repeat, rounded to
                                    whole. */
  float lattice_period = 1.0f; /**< Square lattice period, in plane units,
                                    used when the next stage reports no
                                    lattice. */

  /** @brief Per-parameter id, member, label, range, curve and gate. */
  static constexpr auto FIELDS = std::array{
      Field<AffineParams>{"speed", &AffineParams::speed, nullptr, -0.02f, 0.02f,
                          FieldCurve::LERP},
      Field<AffineParams>{"rotation-rate", &AffineParams::rotation_rate,
                          "Affine Rotation Rate", -4.0f * math::TWO_PI_F,
                          4.0f * math::TWO_PI_F, FieldCurve::LERP,
                          FieldGate::AFFINE_ROTATES},
      Field<AffineParams>{"translation-x", &AffineParams::translation_x,
                          "Repeat Windings A", -8.0f, 8.0f, FieldCurve::SNAP,
                          FieldGate::AFFINE_TRANSLATES},
      Field<AffineParams>{"translation-y", &AffineParams::translation_y,
                          "Repeat Windings B", -8.0f, 8.0f, FieldCurve::SNAP,
                          FieldGate::AFFINE_TRANSLATES},
      Field<AffineParams>{"scale-x", &AffineParams::scale_x, "Affine Scale X",
                          1.0f / 64.0f, 64.0f, FieldCurve::LOG_POSITIVE,
                          FieldGate::AFFINE_SCALES_SHEARS},
      Field<AffineParams>{"scale-y", &AffineParams::scale_y, "Affine Scale Y",
                          1.0f / 64.0f, 64.0f, FieldCurve::LOG_POSITIVE,
                          FieldGate::AFFINE_SCALES_SHEARS},
      Field<AffineParams>{"shear", &AffineParams::shear, "Affine Shear", -4.0f,
                          4.0f, FieldCurve::LERP,
                          FieldGate::AFFINE_SCALES_SHEARS},
      Field<AffineParams>{"repeat-cells", &AffineParams::repeat_cells,
                          "Repeat Cells", 1.0f, 16.0f, FieldCurve::SNAP,
                          FieldGate::AFFINE_REPEAT},
      Field<AffineParams>{"lattice-period", &AffineParams::lattice_period,
                          "Fallback Lattice Period", 1.0f / 64.0f, 100.0f,
                          FieldCurve::LOG_POSITIVE,
                          FieldGate::AFFINE_FALLBACK_LATTICE},
  };
};
static_assert(field_ids_unique<AffineParams>());
static_assert(field_defaults_in_range<AffineParams>());

/** @brief Parts of the affine frame change a stage evaluates, as bit flags;
    an omitted part is compiled out and behaves as its identity. */
enum class AffineForm : uint8_t {
  IDENTITY = 0,      ///< No part: the stage returns its input.
  ROTATES = 1,       ///< Frame rotation by the accumulated `rotation_rate`.
  TRANSLATES = 2,    ///< Phase-driven translation along the lattice.
  SCALES_SHEARS = 4, ///< Phase-oscillated scale and shear.
  GENERAL = ROTATES | TRANSLATES | SCALES_SHEARS, ///< Every part.
};

/** @brief Union of two affine forms.
    @param a First form.
    @param b Second form.
    @return The parts of either. */
constexpr AffineForm operator|(AffineForm a, AffineForm b) {
  return static_cast<AffineForm>(static_cast<uint8_t>(a) |
                                 static_cast<uint8_t>(b));
}

/** @brief Whether @p form evaluates every part of @p part.
    @param form Stage form.
    @param part Part or parts to test.
    @return True when @p form includes @p part. */
constexpr bool includes(AffineForm form, AffineForm part) {
  return (static_cast<uint8_t>(form) & static_cast<uint8_t>(part)) ==
         static_cast<uint8_t>(part);
}

/** @brief Advances the affine frame rotation by one frame.
    @param rotation Accumulated rotation in radians, wrapped into
    (-2pi, 2pi).
    @param params Speed and rotation rate. */
HS_COLD_MEMBER inline void advance_affine_rotation(float &rotation,
                                                   const AffineParams &params) {
  rotation =
      fmodf(rotation + params.speed * params.rotation_rate, math::TWO_PI_F);
}

/** @brief Warp parameters for the polar chart (Pullback::Warp::PolarChart). */
struct PolarParams {
  float speed = 0.0f;         /**< Per-frame advance of the slot's phase, which
                                   offsets the angular coordinate. */
  float radial_scale = 1.0f;  /**< Scale applied to the radial coordinate. */
  float radial_phase = 0.0f;  /**< Offset added to the radial coordinate. */
  float angular_phase = 0.0f; /**< Offset added to the angular coordinate. */

  /** @brief Per-parameter id, member, label, range, curve and gate. */
  static constexpr auto FIELDS = std::array{
      Field<PolarParams>{"speed", &PolarParams::speed, nullptr, -0.02f, 0.02f,
                         FieldCurve::LERP},
      Field<PolarParams>{"radial-scale", &PolarParams::radial_scale,
                         "Polar Radial Scale", 1.0f / 64.0f, 64.0f,
                         FieldCurve::LOG_POSITIVE},
      Field<PolarParams>{"radial-phase", &PolarParams::radial_phase,
                         "Polar Radial Phase", -math::TWO_PI_F, math::TWO_PI_F,
                         FieldCurve::LERP},
      Field<PolarParams>{"angular-phase", &PolarParams::angular_phase,
                         "Polar Angular Phase", -math::TWO_PI_F, math::TWO_PI_F,
                         FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<PolarParams>());
static_assert(field_defaults_in_range<PolarParams>());

/** @brief Warp parameters for the orbiting vortex (Pullback::Warp::Vortex). */
struct VortexParams {
  float speed = 0.0f;    /**< Per-frame advance of the slot's phase, which walks
                              the center around its orbit. */
  float center_x = 0.0f; /**< Orbit center along x, in plane units. */
  float center_y = 0.0f; /**< Orbit center along y, in plane units. */
  float radius = 1.0f;   /**< Radius at which the twist falls to half. */
  float turns = 0.0f;    /**< Twist at the center, in turns. */
  float center_orbit_radius = 0.0f; /**< Radius the center orbits over the
                                         phase cycle; 0 pins the center. */

  /** @brief Per-parameter id, member, label, range, curve and gate. */
  static constexpr auto FIELDS = std::array{
      Field<VortexParams>{"speed", &VortexParams::speed, nullptr, -0.02f, 0.02f,
                          FieldCurve::LERP},
      Field<VortexParams>{"center-x", &VortexParams::center_x,
                          "Vortex Center X", -4.0f, 4.0f, FieldCurve::LERP},
      Field<VortexParams>{"center-y", &VortexParams::center_y,
                          "Vortex Center Y", -4.0f, 4.0f, FieldCurve::LERP},
      Field<VortexParams>{"radius", &VortexParams::radius, "Vortex Radius",
                          1.0f / 64.0f, 8.0f, FieldCurve::LOG_POSITIVE},
      Field<VortexParams>{"turns", &VortexParams::turns, "Vortex Turns", -4.0f,
                          4.0f, FieldCurve::LERP},
      Field<VortexParams>{"center-orbit-radius",
                          &VortexParams::center_orbit_radius,
                          "Vortex Center Orbit", 0.0f, 4.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<VortexParams>());
static_assert(field_defaults_in_range<VortexParams>());

/** @brief Reduction an affine step applies to its displacement before
    measuring its path length. */
enum class RepeatReduction : uint8_t {
  UNWRAPPED,   ///< The raw displacement.
  RECTANGULAR, ///< Centered remainder along each axis.
  GENERAL,     ///< Closest vector of the reduced repeat basis.
};

/** @brief Repeat lattice an affine path length is measured modulo. */
struct PreparedRepeat {
  math::Complex a;           /**< Reduced basis vector; (period x, 0) under
                                  RECTANGULAR. */
  math::Complex b;           /**< Reduced basis vector; (0, period y) under
                                  RECTANGULAR. */
  math::Complex dual_a;      /**< Inverse-basis row: the `a` coordinate of d
                                  is dual_a.re * d.re + dual_a.im * d.im. */
  math::Complex dual_b;      /**< Inverse-basis row of the `b` coordinate. */
  RepeatReduction reduction; /**< Reduction applied per sample. */
};

/** @brief Affine warp coefficients, with the phase oscillation applied. */
struct PreparedAffine {
  float translation_x;   /**< Translation along x, in plane units. */
  float translation_y;   /**< Translation along y, in plane units. */
  float inverse_scale_x; /**< Reciprocal x scale at this frame's phase. */
  float inverse_scale_y; /**< Reciprocal y scale at this frame's phase. */
  float shear;           /**< Shear at this frame's phase. */
  PreparedRepeat repeat; /**< Path-length repeat lattice. */
};

/** @brief Mirror warp offsets, with the phase scroll already folded in. */
struct PreparedMirror {
  float offset_x; /**< Pre-fold translation along x. */
  float offset_y; /**< Pre-fold translation along y. */
};

/** @brief This frame's point on the noise field's closed loop. */
struct PreparedNoiseLoop {
  math::Vector
      offset; /**< Lattice offset the plane coordinate is taken against. */
};

/** @brief Rotation-only slot state, for families with no transform. */
struct PreparedRotation {
  float rotation_cos; /**< Cosine of the slot's rotation angle. */
  float rotation_sin; /**< Sine of the slot's rotation angle. */
};

/** @brief Vortex warp coefficients, with the phase orbit applied. */
struct PreparedVortex {
  float center_x;        /**< Orbited center along x. */
  float center_y;        /**< Orbited center along y. */
  float radius_sq;       /**< Square of the half-twist radius. */
  float angle_numerator; /**< Twist at the center, in radians. */
};

/** @brief Mirror slot state: the rotation pair plus the fold offsets. */
struct PreparedMirrorSlot {
  float rotation_cos; ///< Cosine of the fold-lattice rotation.
  float rotation_sin; ///< Sine of the fold-lattice rotation.
  struct {
    PreparedMirror mirror;
  } transform; ///< Fold offsets.
};

/** @brief Affine slot state: the rotation pair plus the frame coefficients. */
struct PreparedAffineSlot {
  float rotation_cos; ///< Cosine of the frame rotation.
  float rotation_sin; ///< Sine of the frame rotation.
  struct {
    PreparedAffine affine;
  } transform; ///< Frame coefficients.
};

/** @brief Vector-noise slot state: the rotation pair plus the loop point. */
struct PreparedVectorNoiseSlot {
  float rotation_cos; ///< Cosine of the vector rotation.
  float rotation_sin; ///< Sine of the vector rotation.
  struct {
    PreparedNoiseLoop noise_loop;
  } transform; ///< Noise-loop point.
};

/** @brief Vortex slot state: the vortex coefficients; no slot rotation. */
struct PreparedVortexSlot {
  struct {
    PreparedVortex vortex;
  } transform; ///< Vortex coefficients.
};

/**
 * @brief Prepares the wave-shear field rotation.
 * @param warp The slot's parameters.
 * @return Field-axis rotation.
 */
HS_FLASH_INLINE inline PreparedRotation prepare(const WaveShearParams &warp,
                                                float) {
  return {cosf(warp.field_angle), sinf(warp.field_angle)};
}

/**
 * @brief Mirror overload; the phase scrolls the x offset.
 * @param warp The slot's parameters.
 * @param phase The slot's phase clock, in cells of x scroll.
 * @return Fold-lattice rotation and wrapped offsets.
 */
HS_FLASH_INLINE inline PreparedMirrorSlot prepare(const MirrorParams &warp,
                                                  float phase) {
  PreparedMirrorSlot prepared{cosf(warp.rotation), sinf(warp.rotation), {}};
  prepared.transform.mirror = {
      math::wrap_t(warp.offset_x / warp.cell_x + phase) * warp.cell_x,
      math::wrap_t(warp.offset_y / warp.cell_y) * warp.cell_y};
  return prepared;
}

/**
 * @brief Vector-noise overload; the phase walks the noise loop.
 * @param warp The slot's parameters.
 * @param phase The slot's phase clock.
 * @return Vector rotation and noise-loop point.
 */
HS_FLASH_INLINE inline PreparedVectorNoiseSlot
prepare(const VectorNoiseParams &warp, float phase) {
  PreparedVectorNoiseSlot prepared{
      cosf(warp.vector_angle), sinf(warp.vector_angle), {}};
  prepared.transform.noise_loop = {math::noise_projected_loop_offset(phase)};
  return prepared;
}

/**
 * @brief Vortex overload; the phase orbits the center.
 * @param warp The slot's parameters.
 * @param phase The slot's phase clock, in orbits.
 * @return Orbited center, squared radius and twist in radians.
 */
HS_FLASH_INLINE inline PreparedVortexSlot prepare(const VortexParams &warp,
                                                  float phase) {
  const float orbit = math::TWO_PI_F * math::wrap_t(phase);
  return {{{warp.center_x + warp.center_orbit_radius * cosf(orbit),
            warp.center_y + warp.center_orbit_radius * sinf(orbit),
            warp.radius * warp.radius, math::TWO_PI_F * warp.turns}}};
}

/** @brief Translation basis an affine warp winds along.
    @param warp The slot's parameters.
    @param lattice The next stage's lattice; the zero basis selects the square
    `lattice_period` lattice.
    @return The basis in plane units. */
HS_FLASH_INLINE inline TranslationLattice
affine_basis(const AffineParams &warp, const TranslationLattice &lattice) {
  return lattice.present() ? lattice : square_lattice(warp.lattice_period);
}

/** @brief Translation of one phase turn: whole windings of the whole repeat.
    @param warp The slot's parameters.
    @param basis Translation basis from affine_basis().
    @return The translation in plane units. */
HS_FLASH_INLINE inline math::Complex
affine_turn_translation(const AffineParams &warp,
                        const TranslationLattice &basis) {
  const float repeat = std::round(warp.repeat_cells);
  const float windings_a = std::round(warp.translation_x);
  const float windings_b = std::round(warp.translation_y);
  return {
      windings_a * (repeat * basis.a.re) + windings_b * (repeat * basis.b.re),
      windings_a * (repeat * basis.a.im) + windings_b * (repeat * basis.b.im)};
}

/** @brief Reduces the repeat lattice an affine path length is measured
    modulo, with Lagrange reduction.
    @param a First repeat basis vector.
    @param b Second repeat basis vector, not parallel to @p a.
    @return The reduced basis, its inverse rows and the per-sample reduction. */
HS_COLD_MEMBER inline PreparedRepeat prepare_repeat(math::Complex a,
                                                    math::Complex b) {
  if (b.squared_magnitude() < a.squared_magnitude())
    std::swap(a, b);
  for (;;) {
    const float mu =
        std::round((a.re * b.re + a.im * b.im) / a.squared_magnitude());
    if (mu == 0.0f)
      break;
    b = {b.re - mu * a.re, b.im - mu * a.im};
    if (!(b.squared_magnitude() < a.squared_magnitude()))
      break;
    std::swap(a, b);
  }
  if (a.re == 0.0f && b.im == 0.0f)
    std::swap(a, b);
  const bool rectangular = a.im == 0.0f && b.re == 0.0f;
  if (rectangular) {
    a.re = fabsf(a.re);
    b.im = fabsf(b.im);
  }
  const float inverse_det = 1.0f / (a.re * b.im - a.im * b.re);
  return {a,
          b,
          {b.im * inverse_det, -b.re * inverse_det},
          {-a.im * inverse_det, a.re * inverse_det},
          rectangular ? RepeatReduction::RECTANGULAR
                      : RepeatReduction::GENERAL};
}

/**
 * @brief Affine overload; translation winds along the next stage's lattice.
 * @tparam Form Parts to prepare; an omitted part gets its identity
 *         coefficients.
 * @tparam Repeat Whether to prepare the path-length repeat lattice.
 * @param warp The slot's parameters.
 * @param phase The slot's phase clock.
 * @param frame_rotation Accumulated frame rotation for the slot.
 * @param lattice The next stage's TranslationLattice; the zero basis selects
 *        the square `lattice_period` lattice and unwrapped path length.
 * @return Frame rotation and affine coefficients.
 */
template <AffineForm Form = AffineForm::GENERAL, bool Repeat = true>
HS_FLASH_INLINE inline PreparedAffineSlot
prepare(const AffineParams &warp, float phase, float frame_rotation,
        const TranslationLattice &lattice) {
  PreparedAffineSlot prepared{1.0f, 0.0f, {{0.0f, 0.0f, 1.0f, 1.0f, 0.0f, {}}}};
  if constexpr (includes(Form, AffineForm::ROTATES)) {
    prepared.rotation_cos = cosf(frame_rotation);
    prepared.rotation_sin = sinf(frame_rotation);
  }
  if constexpr (includes(Form, AffineForm::SCALES_SHEARS)) {
    const float cycle_cos =
        warp.scale_x == 1.0f && warp.scale_y == 1.0f && warp.shear == 0.0f
            ? 1.0f
            : cosf(math::TWO_PI_F * math::wrap_t(phase));
    auto &affine = prepared.transform.affine;
    affine.inverse_scale_x =
        warp.scale_x == 1.0f ? 1.0f : 1.0f / powf(warp.scale_x, cycle_cos);
    affine.inverse_scale_y =
        warp.scale_y == 1.0f ? 1.0f : 1.0f / powf(warp.scale_y, cycle_cos);
    affine.shear = warp.shear * cycle_cos;
  }
  if constexpr (includes(Form, AffineForm::TRANSLATES) || Repeat) {
    const TranslationLattice basis = affine_basis(warp, lattice);
    if constexpr (includes(Form, AffineForm::TRANSLATES)) {
      const math::Complex turn = affine_turn_translation(warp, basis);
      prepared.transform.affine.translation_x = math::wrap_t(phase) * turn.re;
      prepared.transform.affine.translation_y = math::wrap_t(phase) * turn.im;
    }
    if (Repeat && lattice.present()) {
      const float repeat = std::round(warp.repeat_cells);
      prepared.transform.affine.repeat =
          prepare_repeat({repeat * basis.a.re, repeat * basis.a.im},
                         {repeat * basis.b.re, repeat * basis.b.im});
    }
  }
  return prepared;
}

/** @brief Bound on an affine warp's output magnitude over every phase.
    @param warp The slot's parameters.
    @param input_bound Bound on the input plane coordinates.
    @param lattice The next stage's TranslationLattice.
    @return Diagonal plus shear norm of the frame at its oscillation extreme,
    times @p input_bound, plus a full turn of translation. */
HS_COLD_MEMBER inline float plane_bound(const AffineParams &warp,
                                        float input_bound,
                                        const TranslationLattice &lattice) {
  const float x_gain = fmaxf(warp.scale_x, 1.0f / warp.scale_x);
  const float y_gain = fmaxf(warp.scale_y, 1.0f / warp.scale_y);
  const math::Complex turn =
      affine_turn_translation(warp, affine_basis(warp, lattice));
  return (fmaxf(x_gain, y_gain) + fabsf(warp.shear) * y_gain) * input_bound +
         std::hypot(turn.re, turn.im);
}

/** @brief A provider whose `prepare` builds the slot's per-frame state. */
template <typename State, typename Binding>
concept PreparedProvider = Detail::ProviderFor<State, Binding> &&
                           requires(const typename Binding::FrameState &frame) {
                             State::prepare(frame);
                           };

/** @brief A `PreparedProvider` that also hands out its parameter block. */
template <typename State, typename Binding>
concept ParamsPreparedProvider =
    PreparedProvider<State, Binding> && Detail::ParamsProvider<State, Binding>;

/** @brief Length of a stage delta, or zero when @p required is false.
    @param delta Stage displacement.
    @param required Whether to measure it.
    @return Length of @p delta, or 0. */
__attribute__((always_inline)) inline float
displacement(const math::Complex &delta, bool required) {
  if (!required)
    return 0.0f;
  const float squared = delta.squared_magnitude();
  if (squared > std::numeric_limits<float>::max())
    return std::hypot(delta.re, delta.im);
  return sqrtf(squared);
}

/**
 * @brief Packages a closed-form warp output with its displacement.
 * @param input Coordinate before the warp.
 * @param output Coordinate after the warp.
 * @param path_length_required Measure the step's path length; else 0.
 * @return @p output and the input-to-output distance.
 */
__attribute__((always_inline)) inline WarpStepResult
finish_closed_form(const math::Complex &input, const math::Complex &output,
                   bool path_length_required) {
  const math::Complex delta(output.re - input.re, output.im - input.im);
  return {output, displacement(delta, path_length_required)};
}

/**
 * @brief Envelope factor chosen at run time.
 * @param provenance Projection provenance of the sample.
 * @param edge_width Fade band width, read only by EDGE_FADE.
 * @param mode Envelope to evaluate.
 * @return 1 for FLAT, else the projection weight or edge fade.
 */
__attribute__((always_inline)) inline float
envelope(const ProjectionProvenance &provenance, float edge_width,
         Envelope mode) {
  if (mode == Envelope::PROJECTION_WEIGHT)
    return provenance.value_weight;
  if (mode == Envelope::EDGE_FADE)
    return Detail::edge_fade(provenance, edge_width);
  return 1.0f;
}

/**
 * @brief Envelope factor chosen at compile time.
 * @tparam EnvelopePolicy FlatEnvelope, ProjectionWeightEnvelope or
 *         EdgeFadeEnvelope.
 * @tparam Params Warp parameters; `edge_width` is read by EdgeFadeEnvelope.
 * @param provenance Projection provenance of the sample.
 * @param params The slot's parameters.
 * @return 1 for FlatEnvelope, else the projection weight or edge fade.
 */
template <typename EnvelopePolicy, typename Params>
__attribute__((always_inline)) inline float
fixed_envelope(const ProjectionProvenance &provenance, const Params &params) {
  if constexpr (std::is_same_v<EnvelopePolicy, ProjectionWeightEnvelope>)
    return provenance.value_weight;
  else if constexpr (std::is_same_v<EnvelopePolicy, EdgeFadeEnvelope>)
    return Detail::edge_fade(provenance, params.edge_width);
  else
    return 1.0f;
}

/** @brief Shortest representative of @p delta modulo a GENERAL repeat.
    @param delta Stage displacement.
    @param repeat Reduced repeat basis.
    @return @p delta minus the closest repeat lattice vector, taken among the
    corners of the basis cell containing @p delta. */
HS_HOT_FLASH_MEMBER inline math::Complex
closest_repeat_offset(const math::Complex &delta,
                      const PreparedRepeat &repeat) {
  const float u =
      floorf(repeat.dual_a.re * delta.re + repeat.dual_a.im * delta.im);
  const float v =
      floorf(repeat.dual_b.re * delta.re + repeat.dual_b.im * delta.im);
  const math::Complex base{delta.re - u * repeat.a.re - v * repeat.b.re,
                           delta.im - u * repeat.a.im - v * repeat.b.im};
  math::Complex best = base;
  for (const math::Complex &corner :
       {base - repeat.a, base - repeat.b, base - repeat.a - repeat.b})
    if (corner.squared_magnitude() < best.squared_magnitude())
      best = corner;
  return best;
}

/** @brief Displacement an affine step measures its path length on, under a
    reduction fixed at compile time.
    @tparam Reduction RECTANGULAR for an axis-aligned prepared repeat;
    GENERAL accepts any.
    @param delta Stage displacement.
    @param repeat Prepared repeat lattice.
    @return @p delta reduced modulo @p repeat. */
template <RepeatReduction Reduction>
__attribute__((always_inline)) inline math::Complex
repeat_offset(const math::Complex &delta, const PreparedRepeat &repeat) {
  if constexpr (Reduction == RepeatReduction::RECTANGULAR)
    return {math::centered_remainder(delta.re, repeat.a.re),
            math::centered_remainder(delta.im, repeat.b.im)};
  else if constexpr (Reduction == RepeatReduction::GENERAL)
    return closest_repeat_offset(delta, repeat);
  else
    return delta;
}

/** @brief Displacement an affine step measures its path length on.
    @param delta Stage displacement.
    @param repeat Prepared repeat lattice.
    @return @p delta reduced as `repeat.reduction` selects. */
__attribute__((always_inline)) inline math::Complex
repeat_offset(const math::Complex &delta, const PreparedRepeat &repeat) {
  if (repeat.reduction == RepeatReduction::RECTANGULAR) {
    const auto reduce_axis = [](float value, float period) {
      if (fabsf(value) < period * 0x1p31f)
        return math::centered_remainder(value, period);
      return std::remainder(value, period);
    };
    return {reduce_axis(delta.re, repeat.a.re),
            reduce_axis(delta.im, repeat.b.im)};
  }
  if (repeat.reduction == RepeatReduction::GENERAL)
    return repeat_offset<RepeatReduction::GENERAL>(delta, repeat);
  return delta;
}

/**
 * @brief Rotates into the affine frame, then undoes its scale and shear and
 *        subtracts its translation.
 * @tparam Form Parts to evaluate; an omitted part is skipped as identity.
 * @tparam Prepared PreparedAffineSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param prepared Frame rotation and affine coefficients.
 * @return Warped coordinate.
 */
template <AffineForm Form = AffineForm::GENERAL, typename Prepared>
__attribute__((always_inline)) inline math::Complex
affine_coords(const math::Complex &input, const Prepared &prepared) {
  const auto &affine = prepared.transform.affine;
  math::Complex output = input;
  if constexpr (includes(Form, AffineForm::ROTATES)) {
    const float c = prepared.rotation_cos;
    const float s = prepared.rotation_sin;
    output = {c * input.re + s * input.im, -s * input.re + c * input.im};
  }
  if constexpr (includes(Form, AffineForm::SCALES_SHEARS)) {
    const float y = output.im * affine.inverse_scale_y;
    output = {output.re * affine.inverse_scale_x - affine.shear * y, y};
  }
  if constexpr (includes(Form, AffineForm::TRANSLATES))
    output = {output.re - affine.translation_x,
              output.im - affine.translation_y};
  return output;
}

/**
 * @brief Affine frame change whose path length is measured under a reduction
 *        fixed at compile time.
 * @tparam Reduction RECTANGULAR for an axis-aligned prepared repeat;
 *         GENERAL accepts any.
 * @tparam Form Parts to evaluate; an omitted part is skipped as identity.
 * @tparam Prepared PreparedAffineSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param prepared Frame rotation and affine coefficients.
 * @param path_length_required Measure the step's path length modulo the
 *        prepared repeat lattice; else 0.
 * @return Warped coordinate and path length.
 */
template <RepeatReduction Reduction, AffineForm Form = AffineForm::GENERAL,
          typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
affine_frame(const math::Complex &input, const Prepared &prepared,
             bool path_length_required) {
  const math::Complex output = affine_coords<Form>(input, prepared);
  if (!path_length_required)
    return {output, 0.0f};
  return {output,
          displacement(repeat_offset<Reduction>(
                           output - input, prepared.transform.affine.repeat),
                       true)};
}

/**
 * @brief Affine frame change whose path length is measured under the prepared
 *        repeat's reduction.
 * @tparam Form Parts to evaluate; an omitted part is skipped as identity.
 * @tparam Prepared PreparedAffineSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param prepared Frame rotation and affine coefficients.
 * @param path_length_required Measure the step's path length modulo the
 *        prepared repeat lattice; else 0.
 * @return Warped coordinate and path length.
 */
template <AffineForm Form = AffineForm::GENERAL, typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
affine_frame(const math::Complex &input, const Prepared &prepared,
             bool path_length_required) {
  const math::Complex output = affine_coords<Form>(input, prepared);
  if (!path_length_required)
    return {output, 0.0f};
  return {output, displacement(repeat_offset(output - input,
                                             prepared.transform.affine.repeat),
                               true)};
}

/**
 * @brief Sine shear displacing across the field axis.
 * @tparam Params WaveShearParams-shaped parameter block.
 * @tparam Prepared PreparedRotation-shaped slot state.
 * @param input Plane coordinate.
 * @param params Strength and frequency; zero strength returns @p input.
 * @param phase The slot's phase clock, in cycles.
 * @param amplitude Shear amplitude, the strength times the envelope.
 * @param prepared Field-axis rotation.
 * @param path_length_required Measure the step's path length; else 0.
 * @return Sheared coordinate; path length is the shear offset's magnitude.
 */
template <typename Params, typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
wave_shear(const math::Complex &input, const Params &params, float phase,
           float amplitude, const Prepared &prepared,
           bool path_length_required) {
  if (params.strength == 0.0f)
    return {input, 0.0f};
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const float angle =
      params.frequency * (c * input.re + s * input.im) + math::TWO_PI_F * phase;
  // fast_sinf's reduction loses the low bits past this bound, and the plane
  // coordinate reaches STEREO_INF at the projection pole.
  const float offset =
      amplitude *
      math::fast_sinf(hs::clamp(angle, -projections::STEREO_PATTERN_ARG_LIMIT,
                                projections::STEREO_PATTERN_ARG_LIMIT));
  const math::Complex delta(-s * offset, c * offset);
  return {{input.re + delta.re, input.im + delta.im},
          path_length_required ? fabsf(offset) : 0.0f};
}

/**
 * @brief Folds the rotated plane into mirrored cells.
 * @tparam Params MirrorParams-shaped parameter block.
 * @tparam Prepared PreparedMirrorSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param params Cell width and height.
 * @param prepared Fold-lattice rotation and offsets.
 * @return Folded coordinate, rotated back to the input orientation.
 */
template <typename Params, typename Prepared>
__attribute__((always_inline)) inline math::Complex
mirror_tile_coords(const math::Complex &input, const Params &params,
                   const Prepared &prepared) {
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const float offset_x = prepared.transform.mirror.offset_x;
  const float offset_y = prepared.transform.mirror.offset_y;
  const float x = c * input.re + s * input.im + offset_x;
  const float y = -s * input.re + c * input.im + offset_y;
  const float folded_x =
      params.cell_x *
      (1.0f - 2.0f * fabsf(math::wrap_t(x / params.cell_x) - 0.5f));
  const float folded_y =
      params.cell_y *
      (1.0f - 2.0f * fabsf(math::wrap_t(y / params.cell_y) - 0.5f));
  return {c * folded_x - s * folded_y, s * folded_x + c * folded_y};
}

/**
 * @brief `mirror_tile_coords` as a warp step.
 * @tparam Params MirrorParams-shaped parameter block.
 * @tparam Prepared PreparedMirrorSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param params Cell width and height.
 * @param prepared Fold-lattice rotation and offsets.
 * @param path_length_required Measure the step's path length; else 0.
 * @return Folded coordinate and path length.
 */
template <typename Params, typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
mirror_tile(const math::Complex &input, const Params &params,
            const Prepared &prepared, bool path_length_required) {
  return finish_closed_form(input, mirror_tile_coords(input, params, prepared),
                            path_length_required);
}

/**
 * @brief Twists about the vortex center; the twist halves at the vortex
 *        radius.
 * @tparam Prepared PreparedVortexSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param prepared Center, squared radius and center twist.
 * @param path_length_required Measure the step's path length; else 0.
 * @return Twisted coordinate and path length.
 */
template <typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
vortex(const math::Complex &input, const Prepared &prepared,
       bool path_length_required) {
  const auto &vortex = prepared.transform.vortex;
  const float x = input.re - vortex.center_x;
  const float y = input.im - vortex.center_y;
  const float r_sq = x * x + y * y;
  const float angle = vortex.angle_numerator / (1.0f + r_sq / vortex.radius_sq);
  const float c = math::fast_cosf(angle);
  const float s = math::fast_sinf(angle);
  return finish_closed_form(
      input, {vortex.center_x + c * x - s * y, vortex.center_y + s * x + c * y},
      path_length_required);
}

/** @brief Bound on each curl-vector component. */
inline constexpr float CURL_VECTOR_COMPONENT_MAX = 4.0f;

/**
 * @brief Curl (-dn/dy, dn/dx) of the noise field by central differences.
 * @param input Plane coordinate.
 * @param noise Noise generator.
 * @param basis Octave basis.
 * @param scale Spatial scale of the sampled field.
 * @param loop_offset This frame's noise-loop point.
 * @return Curl vector, each component clamped to
 *         +-`CURL_VECTOR_COMPONENT_MAX`.
 */
HS_O3_FN inline math::Complex curl_vector(const math::Complex &input,
                                          const FastNoiseLite &noise,
                                          math::NoiseBasis basis, float scale,
                                          const math::Vector &loop_offset) {
  const math::Vector q =
      math::noise_projected_coordinate(input, scale, loop_offset);
  const float dx =
      (math::sample_noise_octaves(
           noise, basis,
           q + math::Vector(math::NOISE_STENCIL_RADIUS, 0.0f, 0.0f)) -
       math::sample_noise_octaves(
           noise, basis,
           q - math::Vector(math::NOISE_STENCIL_RADIUS, 0.0f, 0.0f))) /
      (2.0f * math::NOISE_STENCIL_RADIUS);
  const float dy =
      (math::sample_noise_octaves(
           noise, basis,
           q + math::Vector(0.0f, math::NOISE_STENCIL_RADIUS, 0.0f)) -
       math::sample_noise_octaves(
           noise, basis,
           q - math::Vector(0.0f, math::NOISE_STENCIL_RADIUS, 0.0f))) /
      (2.0f * math::NOISE_STENCIL_RADIUS);
  return {hs::clamp(-dy, -CURL_VECTOR_COMPONENT_MAX, CURL_VECTOR_COMPONENT_MAX),
          hs::clamp(dx, -CURL_VECTOR_COMPONENT_MAX, CURL_VECTOR_COMPONENT_MAX)};
}

/**
 * @brief Advects along the noise curl field.
 * @param input Plane coordinate.
 * @param noise Noise generator.
 * @param basis Octave basis.
 * @param intervals Steps to split @p distance into: 1 takes one Euler step,
 *        more take midpoint steps.
 * @param scale Spatial scale of the sampled field.
 * @param distance Signed advection distance; 0 returns @p input.
 * @param loop_offset This frame's noise-loop point.
 * @param path_length_required Measure the step's path length; else 0.
 * @return Advected coordinate and the summed step lengths.
 */
HS_O3_FN inline WarpStepResult
curl_flow(const math::Complex &input, const FastNoiseLite &noise,
          math::NoiseBasis basis, uint8_t intervals, float scale,
          float distance, const math::Vector &loop_offset,
          bool path_length_required) {
  if (distance == 0.0f)
    return {input, 0.0f};
  if (intervals == 1) {
    const math::Complex direction =
        curl_vector(input, noise, basis, scale, loop_offset);
    const math::Complex delta(distance * direction.re, distance * direction.im);
    return {{input.re + delta.re, input.im + delta.im},
            displacement(delta, path_length_required)};
  }
  math::Complex output = input;
  float path_length = 0.0f;
  const float step = distance / intervals;
  for (uint8_t index = 0; index < intervals; ++index) {
    const math::Complex first =
        curl_vector(output, noise, basis, scale, loop_offset);
    const math::Complex midpoint(output.re + 0.5f * step * first.re,
                                 output.im + 0.5f * step * first.im);
    const math::Complex direction =
        curl_vector(midpoint, noise, basis, scale, loop_offset);
    const math::Complex delta(step * direction.re, step * direction.im);
    output = {output.re + delta.re, output.im + delta.im};
    path_length += displacement(delta, path_length_required);
  }
  return {output, path_length};
}

/** @brief A chart change, not a step: the angular output is in radians, so
    the stage contributes no plane-unit path length.
    @tparam Params PolarParams-shaped parameter block.
    @param input Plane coordinate.
    @param params Radial scale and phases.
    @param phase The slot's phase clock, in cycles.
    @param logarithmic Use log radius instead of radius.
    @param harmonic Angular multiplier.
    @return (radial, angular) chart coordinate with zero path length. */
template <typename Params>
__attribute__((always_inline)) inline WarpStepResult
polar_chart(const math::Complex &input, const Params &params, float phase,
            bool logarithmic, uint8_t harmonic) {
  const float radius = input.magnitude();
  const float radial =
      logarithmic ? logf(fmaxf(radius, 1.0f / 4096.0f)) : radius;
  const math::Complex output(params.radial_scale * radial + params.radial_phase,
                             static_cast<float>(harmonic) *
                                     math::fast_atan2(input.im, input.re) +
                                 params.angular_phase + math::TWO_PI_F * phase);
  return {output, 0.0f};
}

/**
 * @brief Displaces by the rotated noise vector for a fixed basis.
 * @tparam BasisV Octave basis.
 * @tparam Params VectorNoiseParams-shaped parameter block.
 * @tparam Prepared PreparedVectorNoiseSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param params Strength and scale; zero strength returns @p input.
 * @param amplitude Displacement amplitude, the strength times the envelope.
 * @param noise Noise generator.
 * @param prepared Vector rotation and noise-loop point.
 * @param path_length_required Measure the step's path length; else 0.
 * @return Displaced coordinate and path length.
 */
template <math::NoiseBasis BasisV, typename Params, typename Prepared>
HS_HOT_FLASH_MEMBER inline WarpStepResult
vector_noise_fixed(const math::Complex &input, const Params &params,
                   float amplitude, const FastNoiseLite &noise,
                   const Prepared &prepared, bool path_length_required) {
  if (params.strength == 0.0f)
    return {input, 0.0f};
  const math::Vector q = math::noise_projected_coordinate(
      input, params.scale, prepared.transform.noise_loop.offset);
  float nx;
  float ny;
  if constexpr (BasisV == math::NoiseBasis::SIMPLEX) {
    const math::Vector field = math::sample_simplex_vector(noise, q);
    nx = field.x;
    ny = field.y;
  } else {
    nx = math::sample_noise_vector_channel(noise, BasisV, q, 0);
    ny = math::sample_noise_vector_channel(noise, BasisV, q, 1);
  }
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const math::Complex delta(amplitude * (c * nx - s * ny),
                            amplitude * (s * nx + c * ny));
  return {{input.re + delta.re, input.im + delta.im},
          displacement(delta, path_length_required)};
}

/**
 * @brief `vector_noise_fixed` dispatched on a run-time basis.
 * @tparam Params VectorNoiseParams-shaped parameter block.
 * @tparam Prepared PreparedVectorNoiseSlot-shaped slot state.
 * @param input Plane coordinate.
 * @param params Strength and scale; zero strength returns @p input.
 * @param amplitude Displacement amplitude, the strength times the envelope.
 * @param noise Noise generator.
 * @param basis Octave basis.
 * @param prepared Vector rotation and noise-loop point.
 * @param path_length_required Measure the step's path length; else 0.
 * @return Displaced coordinate and path length.
 */
template <typename Params, typename Prepared>
HS_HOT_FLASH_MEMBER inline WarpStepResult
vector_noise(const math::Complex &input, const Params &params, float amplitude,
             const FastNoiseLite &noise, math::NoiseBasis basis,
             const Prepared &prepared, bool path_length_required) {
  switch (basis) {
  case math::NoiseBasis::SIMPLEX:
    return vector_noise_fixed<math::NoiseBasis::SIMPLEX>(
        input, params, amplitude, noise, prepared, path_length_required);
  case math::NoiseBasis::FBM3:
    return vector_noise_fixed<math::NoiseBasis::FBM3>(
        input, params, amplitude, noise, prepared, path_length_required);
  case math::NoiseBasis::RIDGED3:
    return vector_noise_fixed<math::NoiseBasis::RIDGED3>(
        input, params, amplitude, noise, prepared, path_length_required);
  }
  return {input, 0.0f};
}

/**
 * @brief AffineFrame warp policy.
 * @tparam State Provider with Binding and FrameState types and
 * prepare(frame), path_length_required(frame) accessors.
 */
template <typename State> struct AffineFrame : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the affine slot
   *        state and the path-length flag.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      PreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::prepare(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_sin } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.affine.inverse_scale_x
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.affine.inverse_scale_y
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.affine.shear
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.affine.translation_x
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.affine.translation_y
        } -> std::convertible_to<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  /// Slot state `State::prepare` resolves.
  using Prepared = StatePrepared<State>;

  /// Parts `State::affine_form` declares; GENERAL when it declares none.
  static constexpr AffineForm FORM = [] {
    if constexpr (requires { State::affine_form(); })
      return State::affine_form();
    else
      return AffineForm::GENERAL;
  }();

  /**
   * @brief Applies the affine frame change.
   * @param input Plane coordinate.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Warped coordinate and the step's path length.
   */
  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    if constexpr (requires { State::repeat_reduction(); })
      return affine_frame<State::repeat_reduction(), FORM>(
          input, prepared, State::path_length_required(frame));
    else
      return affine_frame<FORM>(input, prepared,
                                State::path_length_required(frame));
  }
};

/**
 * @brief WaveShear warp policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), phase(frame), prepare(frame),
 * path_length_required(frame) accessors.
 */
template <typename State, typename EnvelopePolicy = FlatEnvelope>
struct WaveShear : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the shear
   *        parameters, phase, rotation and path-length flag, and `edge_width`
   *        under EdgeFadeEnvelope.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).strength } -> std::convertible_to<float>;
        { State::params(frame).frequency } -> std::convertible_to<float>;
        { State::phase(frame) } -> std::same_as<float>;
        { State::prepare(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_sin } -> std::convertible_to<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } &&
      (!std::is_same_v<EnvelopePolicy, EdgeFadeEnvelope> ||
       requires(const typename CandidateBinding::FrameState &frame) {
         { State::params(frame).edge_width } -> std::convertible_to<float>;
       });

  /// Slot state `State::prepare` resolves.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Applies the sine shear.
   * @param input Plane coordinate.
   * @param provenance Projection provenance the envelope reads.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Warped coordinate and the step's path length.
   */
  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const Prepared &prepared) {
    const auto &params = State::params(frame);
    const float amplitude =
        params.strength * fixed_envelope<EnvelopePolicy>(provenance, params);
    return wave_shear(input, params, State::phase(frame), amplitude, prepared,
                      State::path_length_required(frame));
  }
};

/**
 * @brief Vortex warp policy.
 * @tparam State Provider with Binding and FrameState types and
 * prepare(frame), path_length_required(frame) accessors.
 */
template <typename State> struct Vortex : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the vortex slot
   *        state and the path-length flag.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      PreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        {
          State::prepare(frame).transform.vortex.center_x
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.vortex.center_y
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.vortex.angle_numerator
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.vortex.radius_sq
        } -> std::convertible_to<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  /// Slot state `State::prepare` resolves.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Applies the vortex twist.
   * @param input Plane coordinate.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Warped coordinate and the step's path length.
   */
  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    return vortex(input, prepared, State::path_length_required(frame));
  }
};

/**
 * @brief MirrorTile warp policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), prepare(frame), path_length_required(frame) accessors.
 */
template <typename State> struct MirrorTile : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the cell size,
   *        mirror slot state and path-length flag.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).cell_x } -> std::convertible_to<float>;
        { State::params(frame).cell_y } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_sin } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.mirror.offset_x
        } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.mirror.offset_y
        } -> std::convertible_to<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  /// Slot state `State::prepare` resolves.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Applies the mirrored tiling.
   * @param input Plane coordinate.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Warped coordinate and the step's path length.
   */
  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    return mirror_tile(input, State::params(frame), prepared,
                       State::path_length_required(frame));
  }
};

/**
 * @brief PolarChart warp policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), phase(frame) accessors.
 */
template <typename State, typename PolarMode, uint8_t Harmonic>
struct PolarChart : ApproximationDefaults {
  static_assert(Harmonic >= 1 && Harmonic <= MAX_POLAR_HARMONIC);
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the polar
   *        parameters and phase.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).radial_scale } -> std::convertible_to<float>;
        { State::params(frame).radial_phase } -> std::convertible_to<float>;
        { State::params(frame).angular_phase } -> std::convertible_to<float>;
        { State::phase(frame) } -> std::same_as<float>;
      };

  /**
   * @brief Maps the coordinate to its polar chart; path length is 0.
   * @param input Plane coordinate.
   * @param frame Frame state.
   * @return Warped coordinate and the step's path length.
   */
  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame) {
    return polar_chart(input, State::params(frame), State::phase(frame),
                       std::is_same_v<PolarMode, LogarithmicPolar>, Harmonic);
  }
};

/**
 * @brief VectorNoise warp policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), noise(frame), prepare(frame),
 * path_length_required(frame) accessors.
 */
template <typename State, math::NoiseBasis BasisV, typename EnvelopePolicy>
struct VectorNoise : PrepareFromState<State> {
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the noise
   *        parameters, generator, slot state and path-length flag, and
   *        `edge_width` under EdgeFadeEnvelope.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).strength } -> std::convertible_to<float>;
        { State::params(frame).scale } -> std::convertible_to<float>;
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::prepare(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_sin } -> std::convertible_to<float>;
        {
          State::prepare(frame).transform.noise_loop.offset
        } -> std::convertible_to<math::Vector>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } &&
      (!std::is_same_v<EnvelopePolicy, EdgeFadeEnvelope> ||
       requires(const typename CandidateBinding::FrameState &frame) {
         { State::params(frame).edge_width } -> std::convertible_to<float>;
       });

  /// Slot state `State::prepare` resolves.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Applies the noise-vector displacement.
   * @param input Plane coordinate.
   * @param provenance Projection provenance the envelope reads.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Warped coordinate and the step's path length.
   */
  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const Prepared &prepared) {
    const auto &params = State::params(frame);
    return vector_noise_fixed<BasisV>(
        input, params,
        params.strength * fixed_envelope<EnvelopePolicy>(provenance, params),
        State::noise(frame), prepared, State::path_length_required(frame));
  }
};

/**
 * @brief CurlFlow warp policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), phase(frame), noise(frame), path_length_required(frame)
 * accessors.
 */
template <typename State, math::NoiseBasis BasisV, typename IntegratorPolicy,
          typename EnvelopePolicy = FlatEnvelope>
struct CurlFlow : ApproximationDefaults {
  static_assert(IntegratorPolicy::INTERVALS == 1 ||
                IntegratorPolicy::INTERVALS == 2 ||
                IntegratorPolicy::INTERVALS == 4);
  using Binding = typename State::Binding; ///< Binding of the provider `State`.
  using FrameState = typename State::FrameState; ///< Frame state `State` reads.

  /**
   * @brief Whether `State` under @p CandidateBinding supplies the curl
   *        parameters, phase, generator and path-length flag, and `edge_width`
   *        under EdgeFadeEnvelope.
   * @tparam CandidateBinding Binding being checked.
   */
  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).strength } -> std::convertible_to<float>;
        { State::params(frame).scale } -> std::convertible_to<float>;
        { State::phase(frame) } -> std::same_as<float>;
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } &&
      (!std::is_same_v<EnvelopePolicy, EdgeFadeEnvelope> ||
       requires(const typename CandidateBinding::FrameState &frame) {
         { State::params(frame).edge_width } -> std::convertible_to<float>;
       });

  /** @brief This frame's point on the plane domain's time loop. */
  using Prepared = math::Vector;

  /**
   * @brief Resolves the noise-loop point for the frame's phase.
   * @param frame Frame state.
   * @return Offset the plane coordinate is taken against.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return math::noise_projected_loop_offset(State::phase(frame));
  }

  /**
   * @brief Advects along the curl field.
   * @param input Plane coordinate.
   * @param provenance Projection provenance the envelope reads.
   * @param frame Frame state.
   * @param prepared This frame's `prepare` result.
   * @return Warped coordinate and the step's path length.
   */
  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const Prepared &prepared) {
    const auto &params = State::params(frame);
    const float amplitude =
        params.strength * fixed_envelope<EnvelopePolicy>(provenance, params);
    return curl_flow(input, State::noise(frame), BasisV,
                     IntegratorPolicy::INTERVALS, params.scale, amplitude,
                     prepared, State::path_length_required(frame));
  }
};

} // namespace Warp

} // namespace Pullback
