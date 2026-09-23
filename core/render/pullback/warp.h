/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"
#include "render/pullback/material.h"
#include "math/3dmath.h"
#include "math/projection_patterns.h"
#include <iterator>
#include <limits>

/**
 * @file warp.h
 * @brief Planar warp policies.
 */

namespace Pullback {

namespace Warp {

inline constexpr uint8_t MAX_POLAR_HARMONIC = 16;

struct FlatEnvelope {};
struct ProjectionWeightEnvelope {};
struct EdgeFadeEnvelope {};

/** @brief Envelope shaping a warp's amplitude across the projected domain. */
enum class Envelope : uint8_t {
  FLAT = 0,
  PROJECTION_WEIGHT = 1,
  EDGE_FADE = 2
};

inline constexpr const char *ENVELOPE_IDS[] = {"flat", "projection-weight",
                                               "edge-fade"};
static_assert(std::size(ENVELOPE_IDS) ==
              static_cast<size_t>(Envelope::EDGE_FADE) + 1);

/** @brief Activation relation of the fade band width, which only the edge-fade
    envelope reads. */
inline constexpr TopologyGate ENVELOPE_EDGE_FADE_GATE{
    "envelope", live_values(Envelope::EDGE_FADE)};
struct Euler1 {
  static constexpr uint8_t INTERVALS = 1;
};
struct Midpoint2 {
  static constexpr uint8_t INTERVALS = 2;
};
struct Midpoint4 {
  static constexpr uint8_t INTERVALS = 4;
};
struct LinearPolar {};
struct LogarithmicPolar {};

/**
 * @brief Warp slot placeholder for a slot whose warp policy is an identity.
 * @details `speed` still advances that slot's phase clock, so an effect can drive
 * a phase it exposes no warp for.
 */
struct NoWarpParams {
  float speed = 0.0f; /**< Per-frame advance of the slot's phase. */

  static constexpr auto FIELDS = std::array{
      Field<NoWarpParams>{"speed", &NoWarpParams::speed, nullptr, -0.02f, 0.02f,
                          FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<NoWarpParams>());
static_assert(field_defaults_in_range<NoWarpParams>());

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
 * @details Translation is scaled by the plane units per lattice cell that
 * Warp::prepare receives: the composed path reads that from its
 * LatticeSourceParams source, the chain path fixes it at 1 so translation is
 * in plane units. Only whole windings scroll seamlessly; a fractional
 * translation jumps when the phase wraps.
 */
struct AffineParams {
  float speed = 0.0f;         /**< Per-frame advance of the slot's phase. */
  float rotation_rate = 0.0f; /**< Frame rotation rate; read only in the outer
                                   slot. */
  float translation_x = 0.0f; /**< Translation along x per phase turn. */
  float translation_y = 0.0f; /**< Translation along y per phase turn. */
  float scale_x = 1.0f; /**< Scale along x, oscillated over the phase cycle. */
  float scale_y = 1.0f; /**< Scale along y, oscillated over the phase cycle. */
  float shear = 0.0f;   /**< Shear, oscillated over the phase cycle. */

  static constexpr auto FIELDS = std::array{
      Field<AffineParams>{"speed", &AffineParams::speed, nullptr, -0.02f, 0.02f,
                          FieldCurve::LERP},
      Field<AffineParams>{"rotation-rate", &AffineParams::rotation_rate,
                          "Affine Rotation Rate", -math::TWO_PI_F,
                          math::TWO_PI_F, FieldCurve::LERP},
      Field<AffineParams>{"translation-x", &AffineParams::translation_x,
                          "Affine Translation X", -4.0f, 4.0f,
                          FieldCurve::LERP},
      Field<AffineParams>{"translation-y", &AffineParams::translation_y,
                          "Affine Translation Y", -4.0f, 4.0f,
                          FieldCurve::LERP},
      Field<AffineParams>{"scale-x", &AffineParams::scale_x, "Affine Scale X",
                          1.0f / 64.0f, 64.0f, FieldCurve::LOG_POSITIVE},
      Field<AffineParams>{"scale-y", &AffineParams::scale_y, "Affine Scale Y",
                          1.0f / 64.0f, 64.0f, FieldCurve::LOG_POSITIVE},
      Field<AffineParams>{"shear", &AffineParams::shear, "Affine Shear", -4.0f,
                          4.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<AffineParams>());
static_assert(field_defaults_in_range<AffineParams>());

/** @brief Warp parameters for the polar chart (Pullback::Warp::PolarChart). */
struct PolarParams {
  float speed = 0.0f;         /**< Per-frame advance of the slot's phase, which
                                   offsets the angular coordinate. */
  float radial_scale = 1.0f;  /**< Scale applied to the radial coordinate. */
  float radial_phase = 0.0f;  /**< Offset added to the radial coordinate. */
  float angular_phase = 0.0f; /**< Offset added to the angular coordinate. */

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

/** @brief Affine warp coefficients, with the phase oscillation applied. */
struct PreparedAffine {
  float translation_x; /**< Translation along x, in plane units. */
  float translation_y; /**< Translation along y, in plane units. */
  float scale_x;       /**< Scale along x at this frame's phase. */
  float scale_y;       /**< Scale along y at this frame's phase. */
  float shear;         /**< Shear at this frame's phase. */
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
  float rotation_cos;
  float rotation_sin;
  struct {
    PreparedMirror mirror;
  } transform;
};

/** @brief Affine slot state: the rotation pair plus the frame coefficients. */
struct PreparedAffineSlot {
  float rotation_cos;
  float rotation_sin;
  struct {
    PreparedAffine affine;
  } transform;
};

/** @brief Vector-noise slot state: the rotation pair plus the loop point. */
struct PreparedVectorNoiseSlot {
  float rotation_cos;
  float rotation_sin;
  struct {
    PreparedNoiseLoop noise_loop;
  } transform;
};

/** @brief Vortex slot state: the vortex coefficients alone, the kernel spins
    the plane about the vortex center rather than a slot rotation. */
struct PreparedVortexSlot {
  struct {
    PreparedVortex vortex;
  } transform;
};

/**
 * @brief Resolves one warp slot's per-frame rotation and transform.
 * @details One overload per parameter family, each returning the slot type its
 * warp policy reads. Every overload takes the slot's parameters and its phase
 * clock; only the affine family rotates with the frame, so only its overload
 * takes an accumulated rotation.
 * @param warp The slot's parameters.
 */
HS_FLASH_INLINE inline PreparedRotation prepare(const WaveShearParams &warp,
                                                float) {
  return {cosf(warp.field_angle), sinf(warp.field_angle)};
}

HS_FLASH_INLINE inline PreparedMirrorSlot prepare(const MirrorParams &warp,
                                                  float phase) {
  PreparedMirrorSlot prepared{cosf(warp.rotation), sinf(warp.rotation), {}};
  prepared.transform.mirror = {
      math::wrap_t(warp.offset_x / warp.cell_x + phase) * warp.cell_x,
      math::wrap_t(warp.offset_y / warp.cell_y) * warp.cell_y};
  return prepared;
}

HS_FLASH_INLINE inline PreparedVectorNoiseSlot
prepare(const VectorNoiseParams &warp, float phase) {
  PreparedVectorNoiseSlot prepared{
      cosf(warp.vector_angle), sinf(warp.vector_angle), {}};
  prepared.transform.noise_loop = {noise_projected_loop_offset(phase)};
  return prepared;
}

HS_FLASH_INLINE inline PreparedVortexSlot prepare(const VortexParams &warp,
                                                  float phase) {
  const float orbit = math::TWO_PI_F * math::wrap_t(phase);
  return {{{warp.center_x + warp.center_orbit_radius * cosf(orbit),
            warp.center_y + warp.center_orbit_radius * sinf(orbit),
            warp.radius * warp.radius, math::TWO_PI_F * warp.turns}}};
}

/**
 * @brief Affine overload; translation is scaled from lattice cells to plane
 *        units by @p lattice_period.
 * @param warp The slot's parameters.
 * @param phase The slot's phase clock.
 * @param frame_rotation Accumulated frame rotation for the slot.
 * @param lattice_period Plane units per lattice cell.
 */
HS_FLASH_INLINE inline PreparedAffineSlot prepare(const AffineParams &warp,
                                                  float phase,
                                                  float frame_rotation,
                                                  float lattice_period) {
  PreparedAffineSlot prepared{cosf(frame_rotation), sinf(frame_rotation), {}};
  const float cycle_cos = cosf(math::TWO_PI_F * math::wrap_t(phase));
  prepared.transform.affine = {
      math::wrap_t(phase) * warp.translation_x * lattice_period,
      math::wrap_t(phase) * warp.translation_y * lattice_period,
      powf(warp.scale_x, cycle_cos), powf(warp.scale_y, cycle_cos),
      warp.shear * cycle_cos};
  return prepared;
}

template <typename State, typename Binding>
concept PreparedProvider = Detail::ProviderFor<State, Binding> &&
                           requires(const typename Binding::FrameState &frame) {
                             State::prepare(frame);
                           };

template <typename State, typename Binding>
concept ParamsPreparedProvider =
    PreparedProvider<State, Binding> && Detail::ParamsProvider<State, Binding>;

/** @brief Length of a stage delta, or zero when @p required is false. */
__attribute__((always_inline)) inline float
displacement(const math::Complex &delta, bool required) {
  if (!required)
    return 0.0f;
  const float SQUARED = delta.squared_magnitude();
  if (SQUARED > std::numeric_limits<float>::max())
    return std::hypot(delta.re, delta.im);
  return sqrtf(SQUARED);
}

__attribute__((always_inline)) inline WarpStepResult
finish_closed_form(const math::Complex &input, const math::Complex &output,
                   bool path_length_required) {
  const math::Complex delta(output.re - input.re, output.im - input.im);
  return {output, displacement(delta, path_length_required)};
}

__attribute__((always_inline)) inline float
envelope(const ProjectionProvenance &provenance, float edge_width,
         Envelope mode) {
  if (mode == Envelope::PROJECTION_WEIGHT)
    return provenance.value_weight;
  if (mode == Envelope::EDGE_FADE)
    return ProjectionCoverage::edge_fade(provenance, edge_width);
  return 1.0f;
}

template <typename Envelope, typename Params>
__attribute__((always_inline)) inline float
fixed_envelope(const ProjectionProvenance &provenance, const Params &params) {
  if constexpr (std::is_same_v<Envelope, ProjectionWeightEnvelope>)
    return provenance.value_weight;
  else if constexpr (std::is_same_v<Envelope, EdgeFadeEnvelope>)
    return ProjectionCoverage::edge_fade(provenance, params.edge_width);
  else
    return 1.0f;
}

template <typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
affine_frame(const math::Complex &input, const Prepared &prepared,
             bool path_length_required) {
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const auto &affine = prepared.transform.affine;
  const float rx = c * input.re + s * input.im;
  const float ry = -s * input.re + c * input.im;
  return finish_closed_form(input,
                            {rx / affine.scale_x -
                                 affine.shear * ry / affine.scale_y -
                                 affine.translation_x,
                             ry / affine.scale_y - affine.translation_y},
                            path_length_required);
}

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

template <typename Params, typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
mirror_tile(const math::Complex &input, const Params &params,
            const Prepared &prepared, bool path_length_required) {
  return finish_closed_form(input, mirror_tile_coords(input, params, prepared),
                            path_length_required);
}

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

inline constexpr float CURL_VECTOR_COMPONENT_MAX = 4.0f;

HS_FLASH_INLINE inline math::Complex
curl_vector(const math::Complex &input, const FastNoiseLite &noise,
            ::NoiseBasis basis, float scale, const math::Vector &loop_offset) {
  const math::Vector q = noise_projected_coordinate(input, scale, loop_offset);
  const float dx =
      (sample_noise_octaves(
           noise, basis, q + math::Vector(NOISE_STENCIL_RADIUS, 0.0f, 0.0f)) -
       sample_noise_octaves(
           noise, basis, q - math::Vector(NOISE_STENCIL_RADIUS, 0.0f, 0.0f))) /
      (2.0f * NOISE_STENCIL_RADIUS);
  const float dy =
      (sample_noise_octaves(
           noise, basis, q + math::Vector(0.0f, NOISE_STENCIL_RADIUS, 0.0f)) -
       sample_noise_octaves(
           noise, basis, q - math::Vector(0.0f, NOISE_STENCIL_RADIUS, 0.0f))) /
      (2.0f * NOISE_STENCIL_RADIUS);
  return {hs::clamp(-dy, -CURL_VECTOR_COMPONENT_MAX, CURL_VECTOR_COMPONENT_MAX),
          hs::clamp(dx, -CURL_VECTOR_COMPONENT_MAX, CURL_VECTOR_COMPONENT_MAX)};
}

HS_FLASH_INLINE inline WarpStepResult
curl_flow(const math::Complex &input, const FastNoiseLite &noise,
          ::NoiseBasis basis, uint8_t intervals, float scale, float distance,
          const math::Vector &loop_offset, bool path_length_required) {
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
    the stage contributes no plane-unit path length. */
template <typename Params>
__attribute__((always_inline)) inline WarpStepResult
polar_chart(const math::Complex &input, const Params &params, float phase,
            bool logarithmic, uint8_t harmonic) {
  const float radius = input.magnitude();
  const float radial =
      logarithmic ? logf(std::max(radius, 1.0f / 4096.0f)) : radius;
  const math::Complex output(params.radial_scale * radial + params.radial_phase,
                             static_cast<float>(harmonic) *
                                     math::fast_atan2(input.im, input.re) +
                                 params.angular_phase + math::TWO_PI_F * phase);
  return {output, 0.0f};
}

template <::NoiseBasis BasisV, typename Params, typename Prepared>
HS_FLASH_MEMBER inline WarpStepResult
vector_noise_fixed(const math::Complex &input, const Params &params,
                   float amplitude, const FastNoiseLite &noise,
                   const Prepared &prepared, bool path_length_required) {
  if (params.strength == 0.0f)
    return {input, 0.0f};
  const math::Vector q = noise_projected_coordinate(
      input, params.scale, prepared.transform.noise_loop.offset);
  float nx;
  float ny;
  if constexpr (BasisV == ::NoiseBasis::SIMPLEX) {
    const math::Vector field = sample_simplex_vector(noise, q);
    nx = field.x;
    ny = field.y;
  } else {
    nx = sample_noise_vector_channel(noise, BasisV, q, 0);
    ny = sample_noise_vector_channel(noise, BasisV, q, 1);
  }
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const math::Complex delta(amplitude * (c * nx - s * ny),
                            amplitude * (s * nx + c * ny));
  return {{input.re + delta.re, input.im + delta.im},
          displacement(delta, path_length_required)};
}

template <typename Params, typename Prepared>
HS_FLASH_MEMBER inline WarpStepResult
vector_noise(const math::Complex &input, const Params &params, float amplitude,
             const FastNoiseLite &noise, ::NoiseBasis basis,
             const Prepared &prepared, bool path_length_required) {
  switch (basis) {
  case ::NoiseBasis::SIMPLEX:
    return vector_noise_fixed<::NoiseBasis::SIMPLEX>(
        input, params, amplitude, noise, prepared, path_length_required);
  case ::NoiseBasis::FBM3:
    return vector_noise_fixed<::NoiseBasis::FBM3>(
        input, params, amplitude, noise, prepared, path_length_required);
  case ::NoiseBasis::RIDGED3:
    return vector_noise_fixed<::NoiseBasis::RIDGED3>(
        input, params, amplitude, noise, prepared, path_length_required);
  }
  return {input, 0.0f};
}

template <typename State> struct AffineFrame : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      PreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::prepare(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_sin } -> std::convertible_to<float>;
        State::prepare(frame).transform.affine;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    return affine_frame(input, prepared, State::path_length_required(frame));
  }
};

template <typename State, typename Envelope = FlatEnvelope>
struct WaveShear : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

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
      (!std::is_same_v<Envelope, EdgeFadeEnvelope> ||
       requires(const typename CandidateBinding::FrameState &frame) {
         { State::params(frame).edge_width } -> std::convertible_to<float>;
       });

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const Prepared &prepared) {
    const auto &params = State::params(frame);
    const float amplitude =
        params.strength * fixed_envelope<Envelope>(provenance, params);
    return wave_shear(input, params, State::phase(frame), amplitude, prepared,
                      State::path_length_required(frame));
  }
};

template <typename State> struct Vortex : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      PreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        State::prepare(frame).transform.vortex;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    return vortex(input, prepared, State::path_length_required(frame));
  }
};

template <typename State> struct MirrorTile : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).cell_x } -> std::convertible_to<float>;
        { State::params(frame).cell_y } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_sin } -> std::convertible_to<float>;
        State::prepare(frame).transform.mirror;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    using Instrumentation = typename Binding::Instrumentation;
    const auto start = Instrumentation::mark();
    const WarpStepResult result =
        mirror_tile(input, State::params(frame), prepared,
                    State::path_length_required(frame));
    Instrumentation::template span<ProfileEvent::MIRROR_TILE>(start);
    return result;
  }
};

template <typename State, typename PolarMode, uint8_t Harmonic>
struct PolarChart : ApproximationDefaults {
  static_assert(Harmonic >= 1 && Harmonic <= MAX_POLAR_HARMONIC);
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      Detail::ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).radial_scale } -> std::convertible_to<float>;
        { State::params(frame).radial_phase } -> std::convertible_to<float>;
        { State::params(frame).angular_phase } -> std::convertible_to<float>;
        { State::phase(frame) } -> std::same_as<float>;
      };

  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &,
        const FrameState &frame) {
    return polar_chart(input, State::params(frame), State::phase(frame),
                       std::is_same_v<PolarMode, LogarithmicPolar>, Harmonic);
  }
};

template <typename State, ::NoiseBasis BasisV, typename Envelope>
struct VectorNoise : ApproximationDefaults {
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsPreparedProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).strength } -> std::convertible_to<float>;
        { State::params(frame).scale } -> std::convertible_to<float>;
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::prepare(frame).rotation_cos } -> std::convertible_to<float>;
        { State::prepare(frame).rotation_sin } -> std::convertible_to<float>;
        State::prepare(frame).transform.noise_loop;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } &&
      (!std::is_same_v<Envelope, EdgeFadeEnvelope> ||
       requires(const typename CandidateBinding::FrameState &frame) {
         { State::params(frame).edge_width } -> std::convertible_to<float>;
       });

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const Prepared &prepared) {
    const auto &params = State::params(frame);
    return vector_noise_fixed<BasisV>(
        input, params,
        params.strength * fixed_envelope<Envelope>(provenance, params),
        State::noise(frame), prepared, State::path_length_required(frame));
  }
};

template <typename State, ::NoiseBasis BasisV, typename IntegratorPolicy,
          typename Envelope = FlatEnvelope>
struct CurlFlow : ApproximationDefaults {
  static_assert(IntegratorPolicy::INTERVALS == 1 ||
                IntegratorPolicy::INTERVALS == 2 ||
                IntegratorPolicy::INTERVALS == 4);
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

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
      (!std::is_same_v<Envelope, EdgeFadeEnvelope> ||
       requires(const typename CandidateBinding::FrameState &frame) {
         { State::params(frame).edge_width } -> std::convertible_to<float>;
       });

  /** @brief This frame's point on the plane domain's time loop. */
  using Prepared = math::Vector;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return noise_projected_loop_offset(State::phase(frame));
  }

  __attribute__((always_inline)) static WarpStepResult
  apply(const math::Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const Prepared &prepared) {
    const auto &params = State::params(frame);
    const float amplitude =
        params.strength * fixed_envelope<Envelope>(provenance, params);
    return curl_flow(input, State::noise(frame), BasisV,
                     IntegratorPolicy::INTERVALS, params.scale, amplitude,
                     prepared, State::path_length_required(frame));
  }
};

} // namespace Warp

} // namespace Pullback
