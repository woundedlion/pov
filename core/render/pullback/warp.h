/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

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
 * @brief Warp slot placeholder for a look whose warp policy is an identity.
 * @details `speed` still advances that slot's phase clock, so a look can drive
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
                          "Mirror Rotation", 0.0f, TWO_PI_F,
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

/** @brief Warp parameters for the sine shear (Pullback::Warp::WaveShear). */
struct WaveShearParams {
  float speed = 0.0f;       /**< Per-frame advance of the slot's phase. */
  float strength = 0.0f;    /**< Shear amplitude; 0 skips the stage. */
  float frequency = 1.0f;   /**< Spatial frequency along the field axis. */
  float field_angle = 0.0f; /**< Field axis direction, in radians. */
  float edge_width = 0.1f;  /**< Fade band width, read only under an
                                 EdgeFadeEnvelope. */

  static constexpr auto FIELDS = std::array{
      Field<WaveShearParams>{"speed", &WaveShearParams::speed, nullptr, -0.02f,
                             0.02f, FieldCurve::LERP},
      Field<WaveShearParams>{"strength", &WaveShearParams::strength,
                             "Warp Strength", -30.0f, 30.0f, FieldCurve::LERP},
      Field<WaveShearParams>{"frequency", &WaveShearParams::frequency,
                             "Warp Frequency", 0.01f, 32.0f,
                             FieldCurve::LOG_POSITIVE},
      Field<WaveShearParams>{"field-angle", &WaveShearParams::field_angle,
                             "Warp Field Angle", 0.0f, TWO_PI_F,
                             FieldCurve::SHORTEST_PERIODIC},
      Field<WaveShearParams>{"edge-width", &WaveShearParams::edge_width,
                             nullptr, 0.0f, 1.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<WaveShearParams>());

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
                                  EdgeFadeEnvelope. */

  static constexpr auto FIELDS = std::array{
      Field<VectorNoiseParams>{"speed", &VectorNoiseParams::speed, nullptr,
                               -0.02f, 0.02f, FieldCurve::LERP},
      Field<VectorNoiseParams>{"strength", &VectorNoiseParams::strength,
                               "Warp Strength", -30.0f, 30.0f,
                               FieldCurve::LERP},
      Field<VectorNoiseParams>{"scale", &VectorNoiseParams::scale, "Warp Scale",
                               1.0f / 64.0f, 64.0f, FieldCurve::LOG_POSITIVE},
      Field<VectorNoiseParams>{"vector-angle", &VectorNoiseParams::vector_angle,
                               "Warp Vector Angle", 0.0f, TWO_PI_F,
                               FieldCurve::SHORTEST_PERIODIC},
      Field<VectorNoiseParams>{"edge-width", &VectorNoiseParams::edge_width,
                               nullptr, 0.0f, 1.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<VectorNoiseParams>());

/**
 * @brief Warp parameters for the affine frame change
 *        (Pullback::Warp::AffineFrame).
 * @details Translation is expressed in lattice cells, so this family is only
 * valid alongside a LatticeSourceParams source.
 */
struct AffineParams {
  float speed = 0.0f;         /**< Per-frame advance of the slot's phase. */
  float rotation_rate = 0.0f; /**< Frame rotation rate; read only in the outer
                                   slot. */
  float translation_x = 0.0f; /**< Translation along x, in lattice cells per
                                   phase turn. */
  float translation_y = 0.0f; /**< Translation along y, in lattice cells per
                                   phase turn. */
  float scale_x = 1.0f; /**< Scale along x, oscillated over the phase cycle. */
  float scale_y = 1.0f; /**< Scale along y, oscillated over the phase cycle. */
  float shear = 0.0f;   /**< Shear, oscillated over the phase cycle. */

  static constexpr auto FIELDS = std::array{
      Field<AffineParams>{"speed", &AffineParams::speed, nullptr, -0.02f, 0.02f,
                          FieldCurve::LERP},
      Field<AffineParams>{"rotation-rate", &AffineParams::rotation_rate,
                          "Affine Rotation Rate", -TWO_PI_F, TWO_PI_F,
                          FieldCurve::LERP},
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
                         "Polar Radial Phase", -TWO_PI_F, TWO_PI_F,
                         FieldCurve::SHORTEST_PERIODIC},
      Field<PolarParams>{"angular-phase", &PolarParams::angular_phase,
                         "Polar Angular Phase", -TWO_PI_F, TWO_PI_F,
                         FieldCurve::SHORTEST_PERIODIC},
  };
};
static_assert(field_ids_unique<PolarParams>());

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
  float diagonal; /**< Offset added to both planar noise coordinates. */
  float z;        /**< Third noise coordinate. */
};

/** @brief Rotation-only slot state, for families with no transform. */
struct PreparedRotation {
  float rotation_cos; /**< Cosine of the slot's rotation angle. */
  float rotation_sin; /**< Sine of the slot's rotation angle. */
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
      wrap_t(warp.offset_x / warp.cell_x + phase) * warp.cell_x,
      wrap_t(warp.offset_y / warp.cell_y) * warp.cell_y};
  return prepared;
}

HS_FLASH_INLINE inline PreparedVectorNoiseSlot
prepare(const VectorNoiseParams &warp, float phase) {
  PreparedVectorNoiseSlot prepared{
      cosf(warp.vector_angle), sinf(warp.vector_angle), {}};
  const float angle = TWO_PI_F * wrap_t(phase);
  prepared.transform.noise_loop = {NOISE_LOOP_RADIUS * sinf(angle) *
                                       0.7071067811865475f,
                                   NOISE_LOOP_RADIUS * cosf(angle)};
  return prepared;
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
  const float cycle_cos = cosf(TWO_PI_F * wrap_t(phase));
  prepared.transform.affine = {
      wrap_t(phase) * warp.translation_x * lattice_period,
      wrap_t(phase) * warp.translation_y * lattice_period,
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
concept ParamsProvider =
    Detail::ProviderFor<State, Binding> &&
    requires(const typename Binding::FrameState &frame) {
      State::params(frame);
    } &&
    std::is_lvalue_reference_v<decltype(State::params(
        std::declval<const typename Binding::FrameState &>()))> &&
    std::is_const_v<std::remove_reference_t<decltype(State::params(
        std::declval<const typename Binding::FrameState &>()))>>;

template <typename State, typename Binding>
concept ParamsPreparedProvider =
    PreparedProvider<State, Binding> && ParamsProvider<State, Binding>;

/** @brief Length of a stage delta, or zero when @p required is false. */
__attribute__((always_inline)) inline float displacement(const Complex &delta,
                                                         bool required) {
  return required ? sqrtf(delta.re * delta.re + delta.im * delta.im) : 0.0f;
}

__attribute__((always_inline)) inline WarpStepResult
finish_closed_form(const Complex &input, const Complex &output,
                   bool path_length_required) {
  const Complex delta(output.re - input.re, output.im - input.im);
  return {output, delta, displacement(delta, path_length_required)};
}

__attribute__((always_inline)) inline float
envelope(const ProjectionProvenance &provenance, float edge_width,
         bool projection_weight, bool edge_fade) {
  if (projection_weight)
    return provenance.value_weight;
  if (edge_fade)
    return cubic_kernel(provenance.fade_edge_distance / edge_width);
  return 1.0f;
}

template <typename Envelope, typename Params>
__attribute__((always_inline)) inline float
fixed_envelope(const ProjectionProvenance &provenance, const Params &params) {
  if constexpr (std::is_same_v<Envelope, ProjectionWeightEnvelope>)
    return provenance.value_weight;
  else if constexpr (std::is_same_v<Envelope, EdgeFadeEnvelope>)
    return cubic_kernel(provenance.fade_edge_distance / params.edge_width);
  else
    return 1.0f;
}

template <typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
affine_frame(const Complex &input, const Prepared &prepared,
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
wave_shear(const Complex &input, const Params &params, float phase,
           float amplitude, const Prepared &prepared,
           bool path_length_required) {
  if (params.strength == 0.0f)
    return {input, Complex(), 0.0f};
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const float angle =
      params.frequency * (c * input.re + s * input.im) + TWO_PI_F * phase;
  const float offset = amplitude * sinf(angle);
  const Complex delta(-s * offset, c * offset);
  return {{input.re + delta.re, input.im + delta.im},
          delta,
          path_length_required ? fabsf(offset) : 0.0f};
}

template <typename Params, typename Prepared>
__attribute__((always_inline)) inline Complex
mirror_tile_coords(const Complex &input, const Params &params,
                   const Prepared &prepared) {
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const float offset_x = prepared.transform.mirror.offset_x;
  const float offset_y = prepared.transform.mirror.offset_y;
  const float x = c * input.re + s * input.im + offset_x;
  const float y = -s * input.re + c * input.im + offset_y;
  const float folded_x =
      params.cell_x * (1.0f - 2.0f * fabsf(wrap_t(x / params.cell_x) - 0.5f));
  const float folded_y =
      params.cell_y * (1.0f - 2.0f * fabsf(wrap_t(y / params.cell_y) - 0.5f));
  return {c * folded_x - s * folded_y, s * folded_x + c * folded_y};
}

template <typename Params, typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
mirror_tile(const Complex &input, const Params &params,
            const Prepared &prepared, bool path_length_required) {
  return finish_closed_form(input, mirror_tile_coords(input, params, prepared),
                            path_length_required);
}

template <typename Prepared>
__attribute__((always_inline)) inline WarpStepResult
vortex(const Complex &input, const Prepared &prepared,
       bool path_length_required) {
  const auto &vortex = prepared.transform.vortex;
  const float x = input.re - vortex.center_x;
  const float y = input.im - vortex.center_y;
  const float r_sq = x * x + y * y;
  const float angle = vortex.angle_numerator / (1.0f + r_sq / vortex.radius_sq);
  const float c = cosf(angle);
  const float s = sinf(angle);
  return finish_closed_form(
      input, {vortex.center_x + c * x - s * y, vortex.center_y + s * x + c * y},
      path_length_required);
}

inline constexpr float CURL_VECTOR_COMPONENT_MAX = 4.0f;

HS_FLASH_INLINE inline Complex curl_vector(const Complex &input,
                                           const FastNoiseLite &noise,
                                           ::NoiseBasis basis, float scale,
                                           float phase) {
  const Vector q = noise_projected_coordinate(input, scale, phase);
  const float dx =
      (sample_noise_octaves(noise, basis,
                            q + Vector(NOISE_STENCIL_RADIUS, 0.0f, 0.0f)) -
       sample_noise_octaves(noise, basis,
                            q - Vector(NOISE_STENCIL_RADIUS, 0.0f, 0.0f))) /
      (2.0f * NOISE_STENCIL_RADIUS);
  const float dy =
      (sample_noise_octaves(noise, basis,
                            q + Vector(0.0f, NOISE_STENCIL_RADIUS, 0.0f)) -
       sample_noise_octaves(noise, basis,
                            q - Vector(0.0f, NOISE_STENCIL_RADIUS, 0.0f))) /
      (2.0f * NOISE_STENCIL_RADIUS);
  return {hs::clamp(-dy, -CURL_VECTOR_COMPONENT_MAX, CURL_VECTOR_COMPONENT_MAX),
          hs::clamp(dx, -CURL_VECTOR_COMPONENT_MAX, CURL_VECTOR_COMPONENT_MAX)};
}

HS_FLASH_INLINE inline WarpStepResult
curl_flow(const Complex &input, const FastNoiseLite &noise, ::NoiseBasis basis,
          uint8_t intervals, float scale, float distance, float phase,
          bool path_length_required) {
  if (distance == 0.0f)
    return {input, Complex(), 0.0f};
  if (intervals == 1) {
    const Complex direction = curl_vector(input, noise, basis, scale, phase);
    const Complex delta(distance * direction.re, distance * direction.im);
    return {{input.re + delta.re, input.im + delta.im},
            delta,
            displacement(delta, path_length_required)};
  }
  Complex output = input;
  Complex net_delta;
  float path_length = 0.0f;
  const float step = distance / intervals;
  for (uint8_t index = 0; index < intervals; ++index) {
    const Complex first = curl_vector(output, noise, basis, scale, phase);
    const Complex midpoint(output.re + 0.5f * step * first.re,
                           output.im + 0.5f * step * first.im);
    const Complex direction = curl_vector(midpoint, noise, basis, scale, phase);
    const Complex delta(step * direction.re, step * direction.im);
    output = {output.re + delta.re, output.im + delta.im};
    net_delta = {net_delta.re + delta.re, net_delta.im + delta.im};
    path_length += displacement(delta, path_length_required);
  }
  return {output, net_delta, path_length};
}

template <typename Params>
__attribute__((always_inline)) inline WarpStepResult
polar_chart(const Complex &input, const Params &params, float phase,
            bool logarithmic, uint8_t harmonic, bool path_length_required) {
  const float radius = sqrtf(input.re * input.re + input.im * input.im);
  const float radial =
      logarithmic ? logf(std::max(radius, 1.0f / 4096.0f)) : radius;
  return finish_closed_form(
      input,
      {params.radial_scale * radial + params.radial_phase,
       static_cast<float>(harmonic) * fast_atan2(input.im, input.re) +
           params.angular_phase + TWO_PI_F * phase},
      path_length_required);
}

template <::NoiseBasis BasisV, typename Params, typename Prepared>
HS_FLASH_MEMBER inline WarpStepResult
vector_noise_fixed(const Complex &input, const Params &params, float amplitude,
                   const FastNoiseLite &noise, const Prepared &prepared,
                   bool path_length_required) {
  if (params.strength == 0.0f)
    return {input, Complex(), 0.0f};
  const auto &loop = prepared.transform.noise_loop;
  const Vector q(params.scale * input.re + loop.diagonal,
                 params.scale * input.im + loop.diagonal, loop.z);
  float nx;
  float ny;
  if constexpr (BasisV == ::NoiseBasis::SIMPLEX) {
    const Vector field = sample_simplex_vector(noise, q);
    nx = field.x;
    ny = field.y;
  } else {
    nx = sample_noise_vector_channel(noise, BasisV, q, 0);
    ny = sample_noise_vector_channel(noise, BasisV, q, 1);
  }
  const float c = prepared.rotation_cos;
  const float s = prepared.rotation_sin;
  const Complex delta(amplitude * (c * nx - s * ny),
                      amplitude * (s * nx + c * ny));
  return {{input.re + delta.re, input.im + delta.im},
          delta,
          displacement(delta, path_length_required)};
}

template <typename Params, typename Prepared>
HS_FLASH_MEMBER inline WarpStepResult
vector_noise(const Complex &input, const Params &params, float amplitude,
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
  return {input, Complex(), 0.0f};
}

template <typename State> struct AffineFrame : ExactPolicy {
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
  apply(const Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    return affine_frame(input, prepared, State::path_length_required(frame));
  }
};

template <typename State, typename Envelope = FlatEnvelope>
struct WaveShear : ExactPolicy {
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
  apply(const Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame, const Prepared &prepared) {
    const auto &params = State::params(frame);
    const float amplitude =
        params.strength * fixed_envelope<Envelope>(provenance, params);
    return wave_shear(input, params, State::phase(frame), amplitude, prepared,
                      State::path_length_required(frame));
  }
};

template <typename State> struct Vortex : ExactPolicy {
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
  apply(const Complex &input, const ProjectionProvenance &,
        const FrameState &frame, const Prepared &prepared) {
    return vortex(input, prepared, State::path_length_required(frame));
  }
};

template <typename State> struct MirrorTile : ExactPolicy {
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
  apply(const Complex &input, const ProjectionProvenance &,
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
struct PolarChart : ExactPolicy {
  static_assert(Harmonic >= 1 && Harmonic <= MAX_POLAR_HARMONIC);
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsProvider<State, CandidateBinding> &&
      requires(const typename CandidateBinding::FrameState &frame) {
        { State::params(frame).radial_scale } -> std::convertible_to<float>;
        { State::params(frame).radial_phase } -> std::convertible_to<float>;
        { State::params(frame).angular_phase } -> std::convertible_to<float>;
        { State::phase(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  __attribute__((always_inline)) static WarpStepResult
  apply(const Complex &input, const ProjectionProvenance &,
        const FrameState &frame) {
    return polar_chart(input, State::params(frame), State::phase(frame),
                       std::is_same_v<PolarMode, LogarithmicPolar>, Harmonic,
                       State::path_length_required(frame));
  }
};

template <typename State, ::NoiseBasis BasisV, typename Envelope>
struct VectorNoise : ExactPolicy {
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
  apply(const Complex &input, const ProjectionProvenance &provenance,
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
struct CurlFlow : ExactPolicy {
  static_assert(IntegratorPolicy::INTERVALS == 1 ||
                IntegratorPolicy::INTERVALS == 2 ||
                IntegratorPolicy::INTERVALS == 4);
  using Binding = typename State::Binding;
  using FrameState = typename State::FrameState;

  template <typename CandidateBinding>
  static constexpr bool PROVIDER_VALID =
      ParamsProvider<State, CandidateBinding> &&
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

  __attribute__((always_inline)) static WarpStepResult
  apply(const Complex &input, const ProjectionProvenance &provenance,
        const FrameState &frame) {
    const auto &params = State::params(frame);
    const float amplitude =
        params.strength * fixed_envelope<Envelope>(provenance, params);
    return curl_flow(input, State::noise(frame), BasisV,
                     IntegratorPolicy::INTERVALS, params.scale, amplitude,
                     State::phase(frame), State::path_length_required(frame));
  }
};

} // namespace Warp

} // namespace Pullback
