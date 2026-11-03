/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>

#include "core/math/noise_field.h"
#include "core/math/tangent_field.h"

#include <algorithm>

#include "core/animation/transformer.h"
#include "core/render/pullback/contract.h"
#include "core/render/pullback/fields.h"

/**
 * @file surface.h
 * @brief Surface-displacement policies.
 */

namespace Pullback {

namespace Surface {

/**
 * @brief Surface parameters for the sphere-space curl displacement
 *        (Pullback::Surface::CurlNoise).
 */
struct SurfaceNoiseParams {
  float scale = 1.0f;    /**< Spatial scale of the displacement field. */
  float strength = 0.0f; /**< Displacement distance; 0 skips the stage. */
  float speed = 0.0f;    /**< Per-frame advance of the field's loop phase. */

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      Pullback::Field<SurfaceNoiseParams>{"scale", &SurfaceNoiseParams::scale,
                                          "Surface Noise Scale", 1.0f / 64.0f,
                                          64.0f, FieldCurve::LOG_POSITIVE},
      Pullback::Field<SurfaceNoiseParams>{
          "strength", &SurfaceNoiseParams::strength, "Surface Noise Strength",
          -0.5f, 0.5f, FieldCurve::LERP},
      Pullback::Field<SurfaceNoiseParams>{"speed", &SurfaceNoiseParams::speed,
                                          "Surface Noise Speed", -1.0f / 64.0f,
                                          1.0f / 64.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<SurfaceNoiseParams>());
static_assert(field_defaults_in_range<SurfaceNoiseParams>());

/**
 * @brief Surface parameters for the direction-steered displacement
 *        (Pullback::Surface::DirectNoise).
 */
struct DirectSurfaceParams {
  float scale = 1.0f;     /**< Spatial scale of the displacement field. */
  float strength = 0.0f;  /**< Displacement distance; 0 skips the stage. */
  float speed = 0.0f;     /**< Per-frame advance of the field's loop phase. */
  float direction = 0.0f; /**< Tangent steering, in turns. */

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      Pullback::Field<DirectSurfaceParams>{"scale", &DirectSurfaceParams::scale,
                                           "Surface Noise Scale", 1.0f / 64.0f,
                                           64.0f, FieldCurve::LOG_POSITIVE},
      Pullback::Field<DirectSurfaceParams>{
          "strength", &DirectSurfaceParams::strength, "Surface Noise Strength",
          0.0f, 0.5f, FieldCurve::LERP},
      Pullback::Field<DirectSurfaceParams>{"speed", &DirectSurfaceParams::speed,
                                           "Surface Noise Speed", -1.0f / 64.0f,
                                           1.0f / 64.0f, FieldCurve::LERP},
      Pullback::Field<DirectSurfaceParams>{
          "direction", &DirectSurfaceParams::direction,
          "Surface Noise Direction", 0.0f, 1.0f, FieldCurve::SHORTEST_TURN},
  };
};
static_assert(field_ids_unique<DirectSurfaceParams>());
static_assert(field_defaults_in_range<DirectSurfaceParams>());

/** @brief Parameters for a periodically expanding spherical ripple. */
struct PeriodicRippleParams {
  float period = 80.0f;        /**< Frames per ripple cycle. */
  float strength = 0.15f;      /**< Peak angular displacement. */
  float decay = 0.1f;          /**< Attenuation with distance from center. */
  float thickness = 0.7f;      /**< Angular width of the wavelet. */
  float center_azimuth = 0.0f; /**< Center azimuth, in radians. */
  float center_polar =
      0.5f * math::PI_F; /**< Center polar angle, in radians. */

  /// Parameter registry: id, member, label, range and curve per field.
  static constexpr auto FIELDS = std::array{
      Pullback::Field<PeriodicRippleParams>{
          "period", &PeriodicRippleParams::period, "Ripple Period", 30.0f,
          143.0f, FieldCurve::LERP},
      Pullback::Field<PeriodicRippleParams>{
          "strength", &PeriodicRippleParams::strength, "Ripple Strength", 0.0f,
          0.15f, FieldCurve::LERP},
      Pullback::Field<PeriodicRippleParams>{
          "decay", &PeriodicRippleParams::decay, "Ripple Decay", 0.0f, 5.0f,
          FieldCurve::LERP},
      Pullback::Field<PeriodicRippleParams>{
          "thickness", &PeriodicRippleParams::thickness, "Ripple Thickness",
          0.7f, 1.4f, FieldCurve::LOG_POSITIVE},
      Pullback::Field<PeriodicRippleParams>{
          "center-azimuth", &PeriodicRippleParams::center_azimuth,
          "Ripple Center Azimuth", 0.0f, math::TWO_PI_F,
          FieldCurve::SHORTEST_PERIODIC},
      Pullback::Field<PeriodicRippleParams>{
          "center-polar", &PeriodicRippleParams::center_polar,
          "Ripple Center Polar", 0.0f, math::PI_F, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<PeriodicRippleParams>());
static_assert(field_defaults_in_range<PeriodicRippleParams>());

/** @brief Advances the periodic ripple phase by one frame.
 *  @param phase Ripple clock in frames, wrapped to [0, period).
 *  @param params Ripple parameters supplying the period. */
inline void advance_ripple_phase(float &phase,
                                 const PeriodicRippleParams &params) {
  phase = fmodf(phase + 1.0f, params.period);
}

/** @brief Converts the ripple clock in frames to a normalized cycle.
 *  @param phase Ripple clock in frames.
 *  @param params Ripple parameters supplying the period.
 *  @return @p phase over the period. */
inline float ripple_cycle(float phase, const PeriodicRippleParams &params) {
  return phase / params.period;
}

/** @brief This frame's point on the displacement field's closed loop. */
struct PreparedLoop {
  math::Vector loop_offset;           ///< Noise-space offset of the loop point.
  bool unit_default_gradient = false; ///< Generator transform was qualified.
};

/** @brief Resolves this frame's loop point from the loop phase.
 *  @param phase Position on the time loop, in turns.
 *  @return The prepared loop point. */
HS_FLASH_INLINE inline PreparedLoop prepare(float phase) {
  return {math::noise_sphere_loop_offset(phase)};
}

/** @brief Loop point plus the steering frame the direct displacement reads. */
struct PreparedDirect {
  math::Vector loop_offset; ///< Noise-space offset of the loop point.
  float direction_cos;      ///< cos of the steering angle.
  float direction_sin;      ///< sin of the steering angle.
};

/** @brief Prepared parameters for one frame of a periodic ripple. */
struct PreparedRipple {
  Animation::RippleWave ripple; ///< Synced ripple with the cycle envelope.
};

/** @brief Resolves the configured fixed ripple center.
 * @param params Ripple controls.
 * @return Unit center direction. */
HS_FLASH_INLINE inline math::Vector
prepare_ripple_center(const PeriodicRippleParams &params) {
  return math::Vector::from_spherical(params.center_azimuth,
                                      params.center_polar);
}

/** @brief Resolves a seamless ripple cycle around a supplied center.
 *  @param params Ripple parameters.
 *  @param cycle Normalized cycle; wrapped to [0, 1).
 *  @param center Prepared center direction.
 *  @return The synced ripple, its amplitude enveloped over the cycle. */
HS_FLASH_INLINE inline PreparedRipple
prepare_ripple(const PeriodicRippleParams &params, float cycle,
               const math::Vector &center) {
  Animation::RippleWave ripple;
  ripple.center = center;
  ripple.amplitude = params.strength;
  ripple.decay = params.decay;
  ripple.thickness = params.thickness;
  const float progress = math::wrap_t(cycle);
  ripple.phase = progress * math::PI_F;
  const float attack = fminf(progress * 10.0f, 1.0f);
  ripple.amplitude *= attack * attack * (1.0f - progress);
  ripple.sync();
  return {ripple};
}

/** @brief Resolves a seamless ripple cycle around the configured center.
 * @param params Ripple parameters.
 * @param cycle Normalized cycle; wrapped to [0, 1).
 * @return Synced ripple with its cycle envelope. */
HS_FLASH_INLINE inline PreparedRipple
prepare_ripple(const PeriodicRippleParams &params, float cycle) {
  return prepare_ripple(params, cycle, prepare_ripple_center(params));
}

/** @brief Applies a prepared periodic ripple and reports its angular travel.
 *  @param input Unit sphere direction.
 *  @param prepared This frame's ripple.
 *  @param path_length_required Whether to compute the angular travel.
 *  @return Displaced direction and travel in radians; travel is 0 when not
 *          required or when the point did not move. */
__attribute__((always_inline)) inline SurfaceResult
periodic_ripple(const math::Vector &input, const PreparedRipple &prepared,
                bool path_length_required) {
  const math::Vector displaced = ripple_transform(input, prepared.ripple);
  if (!path_length_required ||
      (displaced.x == input.x && displaced.y == input.y &&
       displaced.z == input.z))
    return {displaced, 0.0f};
  return {displaced,
          math::fast_acos(hs::clamp(math::dot(input, displaced), -1.0f, 1.0f))};
}

/** @brief Resolves the loop point and steering frame for DirectNoise.
 *  @param phase Position on the time loop, in turns.
 *  @param direction Steering angle, in turns.
 *  @return The prepared loop point and steering cos/sin. */
HS_FLASH_INLINE inline PreparedDirect prepare_direct(float phase,
                                                     float direction) {
  const float angle = math::TWO_PI_F * direction;
  return {prepare(phase).loop_offset, cosf(angle), sinf(angle)};
}

/** @brief Integration scheme of a curl-noise displacement. */
enum class Integrator : uint8_t {
  EULER,      ///< One step along the field at the start point.
  MIDPOINT,   ///< One midpoint step.
  MIDPOINT_2X ///< Two half-length midpoint steps.
};

/** @brief Integrator policy tag selecting Integrator::EULER. */
struct Euler {
  static constexpr Integrator VALUE = Integrator::EULER; ///< Selected scheme.
};

/** @brief Integrator policy tag selecting Integrator::MIDPOINT. */
struct Midpoint {
  static constexpr Integrator VALUE =
      Integrator::MIDPOINT; ///< Selected scheme.
};

/** @brief Integrator policy tag selecting Integrator::MIDPOINT_2X. */
struct Midpoint2 {
  static constexpr Integrator VALUE =
      Integrator::MIDPOINT_2X; ///< Selected scheme.
};

/**
 * @brief Arc length of a tangent step, if required.
 * @param step Tangent step; its length is the arc in radians.
 * @param required Whether to compute the length.
 * @return |step| when `required`, else 0.
 */
__attribute__((always_inline)) inline float
path_length(const math::Vector &step, bool required) {
  return required ? sqrtf(math::dot(step, step)) : 0.0f;
}

/**
 * @brief Moves a sphere point along a tangent step.
 * @param input Unit start point.
 * @param step Tangent at `input`; length at most half a radian.
 * @param path_length_required Whether to report the step's arc length.
 * @return The point reached and, when required, the arc length.
 */
__attribute__((always_inline)) inline SurfaceResult
finish_step(const math::Vector &input, const math::Vector &step,
            bool path_length_required) {
  return {math::sphere_exp_map_half_radian(input, step),
          path_length(step, path_length_required)};
}

/** @brief Takes one fixed step through a tangent field.
 * @tparam F Prepared tangent field.
 * @tparam IntegratorPolicy Fixed stepping scheme.
 * @param input Unit start point.
 * @param field Prepared velocity and bounds.
 * @param strength Signed displacement scale.
 * @param path_length_required Whether to report angular travel.
 * @return Displaced direction and accumulated arc length. */
template <math::tangent_field::TangentField F, typename IntegratorPolicy>
__attribute__((always_inline)) HS_O3_FN inline SurfaceResult
step_field(const math::Vector &input, const F &field, float strength,
           bool path_length_required) {
  if (strength == 0.0f)
    return {input, 0.0f};
  const bool SMALL = fabsf(strength) * field.speed_bound() <= 0.5f;
  const auto exp_map = [SMALL](const math::Vector &p, const math::Vector &v) {
    return SMALL ? math::sphere_exp_map_half_radian(p, v)
                 : math::sphere_exp_map(p, v);
  };
  const auto finish = [&](const math::Vector &p, const math::Vector &step) {
    return SurfaceResult{exp_map(p, step),
                         path_length(step, path_length_required)};
  };
  if constexpr (IntegratorPolicy::VALUE == Integrator::EULER)
    return finish(input, strength * field.velocity(input));
  else {
    const auto midpoint = [&](const math::Vector &p, float distance) {
      const math::Vector half_step = (0.5f * distance) * field.velocity(p);
      const math::Vector middle = exp_map(p, half_step);
      math::Vector velocity = field.velocity(middle);
      if (SMALL)
        return finish(p,
                      distance * math::parallel_transport(middle, p, velocity));
      const float angle = half_step.length();
      if (angle != 0.0f) {
        const math::Vector axis = math::cross(p, half_step) / angle;
        const float cosine = cosf(angle);
        const float sine = sinf(angle);
        velocity = cosine * velocity - sine * math::cross(axis, velocity) +
                   (1.0f - cosine) * math::dot(axis, velocity) * axis;
      }
      return finish(p, distance * velocity);
    };
    if constexpr (IntegratorPolicy::VALUE == Integrator::MIDPOINT)
      return midpoint(input, strength);
    else {
      const SurfaceResult first = midpoint(input, 0.5f * strength);
      SurfaceResult second = midpoint(first.sphere, 0.5f * strength);
      second.path_length += first.path_length;
      return second;
    }
  }
}

/** @brief Fixed-step surface policy over a prepared tangent field.
 * @tparam State Provider of prepare, strength and path_length_required.
 * @tparam IntegratorPolicy Fixed stepping scheme. */
template <typename State, typename IntegratorPolicy>
struct Field : PrepareFromState<State> {
  using FrameState = typename State::FrameState; ///< Provider frame.
  using Prepared = StatePrepared<State>;         ///< Prepared tangent field.

  /** @brief Checks the field provider contract.
   * @tparam Binding Chain binding. */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      math::tangent_field::TangentField<Prepared> &&
      requires(const typename Binding::FrameState &frame) {
        { State::strength(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
        { IntegratorPolicy::VALUE } -> std::convertible_to<Integrator>;
      };

  /** @brief Displaces one sphere point.
   * @param input Unit direction.
   * @param frame Provider frame.
   * @param prepared Prepared field.
   * @return Direction and angular travel. */
  HS_O3_FN static SurfaceResult apply(const math::Vector &input,
                                      const FrameState &frame,
                                      const Prepared &prepared) {
    return step_field<Prepared, IntegratorPolicy>(
        input, prepared, State::strength(frame),
        State::path_length_required(frame));
  }
};

/** @brief Prepared direct-noise tangent field. */
struct DirectNoiseField {
  const FastNoiseLite &noise; ///< Prepared generator.
  math::NoiseBasis basis;     ///< Octave structure.
  float scale;                ///< Spatial frequency.
  math::Vector loop_offset;   ///< Closed-loop position.
  float direction_cos;        ///< Steering cosine.
  float direction_sin;        ///< Steering sine.
  /** @brief Samples the steered tangent.
   * @param p Unit direction.
   * @return Bounded tangent velocity. */
  __attribute__((always_inline)) math::Vector
  velocity(const math::Vector &p) const {
    return math::sample_direct_tangent(
        noise, basis, math::noise_sphere_coordinate(p, scale, loop_offset), p,
        direction_cos, direction_sin);
  }
  /** @brief Uniform speed bound. @return One. */
  constexpr float speed_bound() const { return 1.0f; }
  /** @brief Covariant derivative bound. @return Unbounded noise derivative. */
  constexpr float lipschitz_bound() const { return INFINITY; }
};

/** @brief Prepared runtime-topology curl-noise tangent field. */
struct DynamicCurlNoiseField {
  const FastNoiseLite &noise; ///< Prepared generator.
  math::NoiseBasis basis;     ///< Octave structure.
  math::TangentLimit limit;   ///< Velocity saturation.
  float scale;                ///< Spatial frequency.
  math::Vector loop_offset;   ///< Closed-loop position.
  /** @brief Samples the bounded curl tangent.
   * @param p Unit direction.
   * @return Bounded tangent velocity. */
  __attribute__((always_inline)) math::Vector
  velocity(const math::Vector &p) const {
    return math::sample_curl_tangent(
        noise, basis, math::noise_sphere_coordinate(p, scale, loop_offset), p,
        limit);
  }
  /** @brief Uniform speed bound. @return One. */
  constexpr float speed_bound() const { return 1.0f; }
  /** @brief Covariant derivative bound. @return Unbounded noise derivative. */
  constexpr float lipschitz_bound() const { return INFINITY; }
};

/** @brief Prepared curl field retaining qualified gradient specialization.
 * @tparam Basis Octave structure.
 * @tparam Limit Velocity saturation.
 * @tparam InlineUnitGradient Whether to inline the qualified gradient. */
template <math::NoiseBasis Basis, math::TangentLimit Limit,
          bool InlineUnitGradient>
struct CurlNoiseField {
  const FastNoiseLite &noise; ///< Prepared generator.
  float scale;                ///< Spatial frequency.
  math::Vector loop_offset;   ///< Closed-loop position.
  bool unit_default_gradient; ///< Qualified generator transform.
  /** @brief Samples the qualified or general curl tangent.
   * @param p Unit direction.
   * @return Bounded tangent velocity. */
  __attribute__((always_inline)) math::Vector
  velocity(const math::Vector &p) const {
    const math::Vector coordinate =
        math::noise_sphere_coordinate(p, scale, loop_offset);
    if constexpr (Basis == math::NoiseBasis::SIMPLEX)
      if (unit_default_gradient)
        return math::sample_unit_simplex_curl_tangent<InlineUnitGradient>(
            noise, coordinate, p, Limit);
    return math::sample_curl_tangent(noise, Basis, coordinate, p, Limit);
  }
  /** @brief Uniform speed bound. @return One. */
  constexpr float speed_bound() const { return 1.0f; }
  /** @brief Covariant derivative bound. @return Unbounded noise derivative. */
  constexpr float lipschitz_bound() const { return INFINITY; }
};

/**
 * @brief Displaces a sphere point along a steered noise tangent.
 * @param input Unit point on the sphere.
 * @param noise Prepared noise generator.
 * @param basis Octave structure to apply.
 * @param scale Spatial scale of the noise field.
 * @param loop_offset Noise-space offset of this frame's loop point.
 * @param strength Displacement distance; 0 returns `input` unchanged.
 * @param direction_cos cos of the steering angle.
 * @param direction_sin sin of the steering angle.
 * @param path_length_required Whether to report the arc length.
 * @return Displaced point and its arc length when required.
 */
__attribute__((always_inline)) inline SurfaceResult
direct_noise(const math::Vector &input, const FastNoiseLite &noise,
             math::NoiseBasis basis, float scale,
             const math::Vector &loop_offset, float strength,
             float direction_cos, float direction_sin,
             bool path_length_required) {
  const DirectNoiseField field{noise,       basis,         scale,
                               loop_offset, direction_cos, direction_sin};
  return step_field<DirectNoiseField, Euler>(input, field, strength,
                                             path_length_required);
}

/**
 * @brief Curl-noise tangent at a sphere point.
 * @param input Unit point on the sphere.
 * @param noise Prepared noise generator.
 * @param basis Octave structure to apply.
 * @param limit Length limit applied to the tangent.
 * @param scale Spatial scale of the noise field.
 * @param loop_offset Noise-space offset of this frame's loop point.
 * @return Tangent at `input` of length at most 1.
 */
__attribute__((always_inline)) inline math::Vector
curl_field(const math::Vector &input, const FastNoiseLite &noise,
           math::NoiseBasis basis, math::TangentLimit limit, float scale,
           const math::Vector &loop_offset) {
  const math::Vector q =
      math::noise_sphere_coordinate(input, scale, loop_offset);
  return math::sample_curl_tangent(noise, basis, q, input, limit);
}

/**
 * @brief One midpoint-rule step through the curl field.
 * @param input Unit start point.
 * @param noise Prepared noise generator.
 * @param basis Octave structure to apply.
 * @param limit Length limit applied to the tangent.
 * @param scale Spatial scale of the noise field.
 * @param loop_offset Noise-space offset of this frame's loop point.
 * @param distance Step length scale; the midpoint field is transported back
 * to `input`.
 * @param path_length_required Whether to report the arc length.
 * @return The point reached and, when required, the arc length.
 */
__attribute__((always_inline)) inline SurfaceResult
curl_midpoint_step(const math::Vector &input, const FastNoiseLite &noise,
                   math::NoiseBasis basis, math::TangentLimit limit,
                   float scale, const math::Vector &loop_offset, float distance,
                   bool path_length_required) {
  const DynamicCurlNoiseField field{noise, basis, limit, scale, loop_offset};
  return step_field<DynamicCurlNoiseField, Midpoint>(input, field, distance,
                                                     path_length_required);
}

/**
 * @brief Displaces a sphere point through the curl field.
 * @param input Unit point on the sphere.
 * @param noise Prepared noise generator.
 * @param basis Octave structure to apply.
 * @param integrator Integration scheme.
 * @param limit Length limit applied to the tangent.
 * @param scale Spatial scale of the noise field.
 * @param loop_offset Noise-space offset of this frame's loop point.
 * @param strength Displacement distance; 0 returns `input` unchanged.
 * @param path_length_required Whether to report the total arc length.
 * @return Displaced point and its arc length when required.
 */
HS_O3_FN inline SurfaceResult
curl_noise(const math::Vector &input, const FastNoiseLite &noise,
           math::NoiseBasis basis, Integrator integrator,
           math::TangentLimit limit, float scale,
           const math::Vector &loop_offset, float strength,
           bool path_length_required) {
  const DynamicCurlNoiseField field{noise, basis, limit, scale, loop_offset};
  switch (integrator) {
  case Integrator::EULER:
    return step_field<DynamicCurlNoiseField, Euler>(input, field, strength,
                                                    path_length_required);
  case Integrator::MIDPOINT:
    return step_field<DynamicCurlNoiseField, Midpoint>(input, field, strength,
                                                       path_length_required);
  default:
    return step_field<DynamicCurlNoiseField, Midpoint2>(input, field, strength,
                                                        path_length_required);
  }
}

/**
 * @brief DirectNoise surface policy.
 * @tparam State Provider with Binding and FrameState types and
 * noise(frame), scale(frame), strength(frame), prepare(frame),
 * path_length_required(frame) accessors.
 */
template <typename State, math::NoiseBasis Basis>
struct DirectNoise : PrepareFromState<State> {
  using FrameState = typename State::FrameState; ///< The provider's frame.

  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::strength(frame) } -> std::same_as<float>;
        {
          State::prepare(frame).loop_offset
        } -> std::convertible_to<math::Vector>;
        { State::prepare(frame).direction_cos } -> std::convertible_to<float>;
        { State::prepare(frame).direction_sin } -> std::convertible_to<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  /// The provider's prepared type, e.g. PreparedDirect.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Displaces one sphere point.
   * @param input Unit point on the sphere.
   * @param frame Current frame state.
   * @param prepared This frame's prepared state.
   * @return Displaced point and its path length when required.
   */
  __attribute__((always_inline)) static SurfaceResult
  apply(const math::Vector &input, const FrameState &frame,
        const Prepared &prepared) {
    return direct_noise(input, State::noise(frame), Basis, State::scale(frame),
                        prepared.loop_offset, State::strength(frame),
                        prepared.direction_cos, prepared.direction_sin,
                        State::path_length_required(frame));
  }
};

/**
 * @brief CurlNoise surface policy.
 * @tparam State Provider with Binding and FrameState types and
 * noise(frame), scale(frame), strength(frame), prepare(frame),
 * path_length_required(frame) accessors.
 */
template <typename State, math::NoiseBasis Basis, typename IntegratorPolicy,
          math::TangentLimit Limit>
struct CurlNoise : PrepareFromState<State> {
  using FrameState = typename State::FrameState; ///< The provider's frame.
  static constexpr bool OWNED_GRADIENT_PREPARE = [] {
    if constexpr (requires { typename State::GradientPrepareOwner; })
      return std::is_same_v<typename State::GradientPrepareOwner, State>;
    else
      return false;
  }(); ///< Whether the provider owns generator qualification.
  static constexpr bool INLINE_UNIT_GRADIENT = [] {
    if constexpr (requires { IntegratorPolicy::INLINE_UNIT_GRADIENT; })
      return IntegratorPolicy::INLINE_UNIT_GRADIENT;
    else
      return false;
  }(); ///< Whether the qualified gradient core is emitted at the caller.

  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::strength(frame) } -> std::same_as<float>;
        {
          State::prepare(frame).loop_offset
        } -> std::convertible_to<math::Vector>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } && requires {
        { IntegratorPolicy::VALUE } -> std::convertible_to<Integrator>;
      };

  /// The provider's prepared type, e.g. PreparedLoop.
  using Prepared = StatePrepared<State>;

  /**
   * @brief Displaces one sphere point.
   * @param input Unit point on the sphere.
   * @param frame Current frame state.
   * @param prepared This frame's prepared state.
   * @return Displaced point and its path length when required.
   */
  __attribute__((always_inline)) HS_O3_FN static SurfaceResult
  apply(const math::Vector &input, const FrameState &frame,
        const Prepared &prepared) {
    bool qualified = false;
    if constexpr (Basis == math::NoiseBasis::SIMPLEX &&
                  OWNED_GRADIENT_PREPARE &&
                  requires { prepared.unit_default_gradient; })
      qualified = prepared.unit_default_gradient;
    using NoiseField = CurlNoiseField<Basis, Limit, INLINE_UNIT_GRADIENT>;
    const NoiseField field{State::noise(frame), State::scale(frame),
                           prepared.loop_offset, qualified};
    return step_field<NoiseField, IntegratorPolicy>(
        input, field, State::strength(frame),
        State::path_length_required(frame));
  }
};

/**
 * @brief PeriodicRipple surface policy.
 * @tparam State Provider with Binding and FrameState types and
 * params(frame), phase(frame), path_length_required(frame) accessors.
 * phase(frame) returns a normalized cycle, as produced by ripple_cycle().
 */
template <typename State> struct PeriodicRipple : ApproximationDefaults {
  using FrameState = typename State::FrameState; ///< The provider's frame.

  /**
   * @brief Whether `State` is a provider for `Binding` with every accessor
   * this policy reads.
   * @tparam Binding Chain binding to check against.
   */
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::params(frame) } -> std::same_as<const PeriodicRippleParams &>;
        { State::phase(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  using Prepared = PreparedRipple; ///< This frame's resolved ripple.

  /**
   * @brief Resolves this frame's ripple.
   * @param frame Current frame state.
   * @return prepare_ripple() of the frame's params and cycle.
   */
  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return prepare_ripple(State::params(frame), State::phase(frame));
  }

  /**
   * @brief Displaces one sphere point.
   * @param input Unit point on the sphere.
   * @param frame Current frame state.
   * @param prepared This frame's prepared state.
   * @return Displaced point and its path length when required.
   */
  HS_O3_FN static SurfaceResult apply(const math::Vector &input,
                                      const FrameState &frame,
                                      const Prepared &prepared) {
    return periodic_ripple(input, prepared, State::path_length_required(frame));
  }
};

} // namespace Surface

} // namespace Pullback
