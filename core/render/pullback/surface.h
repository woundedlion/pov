/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>

#include "animation/transformer.h"
#include "render/pullback/contract.h"
#include "render/pullback/fields.h"

/**
 * @file surface.h
 * @brief Surface-displacement policies.
 */

namespace Pullback {

namespace Surface {

/** @brief Surface stage placeholder for an effect with no displacement. */
struct NoSurfaceParams {
  static constexpr std::array<Field<NoSurfaceParams>, 0> FIELDS{};
};
static_assert(field_ids_unique<NoSurfaceParams>());

/**
 * @brief Surface parameters for the sphere-space curl displacement
 *        (Pullback::Surface::CurlNoise).
 */
struct SurfaceNoiseParams {
  float scale = 1.0f;    /**< Spatial scale of the displacement field. */
  float strength = 0.0f; /**< Displacement distance; 0 skips the stage. */
  float speed = 0.0f;    /**< Per-frame advance of the field's loop phase. */

  static constexpr auto FIELDS = std::array{
      Field<SurfaceNoiseParams>{"scale", &SurfaceNoiseParams::scale,
                                "Surface Noise Scale", 1.0f / 64.0f, 64.0f,
                                FieldCurve::LERP},
      Field<SurfaceNoiseParams>{"strength", &SurfaceNoiseParams::strength,
                                "Surface Noise Strength", -0.5f, 0.5f,
                                FieldCurve::LERP},
      Field<SurfaceNoiseParams>{"speed", &SurfaceNoiseParams::speed,
                                "Surface Noise Speed", -1.0f / 64.0f,
                                1.0f / 64.0f, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<SurfaceNoiseParams>());

/**
 * @brief Surface parameters for the direction-steered displacement
 *        (Pullback::Surface::DirectNoise).
 */
struct DirectSurfaceParams {
  float scale = 1.0f;     /**< Spatial scale of the displacement field. */
  float strength = 0.0f;  /**< Displacement distance; 0 skips the stage. */
  float speed = 0.0f;     /**< Per-frame advance of the field's loop phase. */
  float direction = 0.0f; /**< Tangent steering, in turns. */

  static constexpr auto FIELDS = std::array{
      Field<DirectSurfaceParams>{"scale", &DirectSurfaceParams::scale,
                                 "Surface Noise Scale", 1.0f / 64.0f, 64.0f,
                                 FieldCurve::LERP},
      Field<DirectSurfaceParams>{"strength", &DirectSurfaceParams::strength,
                                 "Surface Noise Strength", 0.0f, 0.5f,
                                 FieldCurve::LERP},
      Field<DirectSurfaceParams>{"speed", &DirectSurfaceParams::speed,
                                 "Surface Noise Speed", -1.0f / 64.0f,
                                 1.0f / 64.0f, FieldCurve::LERP},
      Field<DirectSurfaceParams>{"direction", &DirectSurfaceParams::direction,
                                 "Surface Noise Direction", 0.0f, 1.0f,
                                 FieldCurve::SHORTEST_TURN},
  };
};
static_assert(field_ids_unique<DirectSurfaceParams>());

/** @brief Parameters for a periodically expanding spherical ripple. */
struct PeriodicRippleParams {
  float period = 80.0f;        /**< Frames per ripple cycle. */
  float strength = 0.15f;      /**< Peak angular displacement. */
  float decay = 0.1f;          /**< Attenuation with distance from center. */
  float thickness = 0.7f;      /**< Angular width of the wavelet. */
  float center_azimuth = 0.0f; /**< Center azimuth, in radians. */
  float center_polar = 0.5f * PI_F; /**< Center polar angle, in radians. */

  static constexpr auto FIELDS = std::array{
      Field<PeriodicRippleParams>{"period", &PeriodicRippleParams::period,
                                  "Ripple Period", 30.0f, 143.0f,
                                  FieldCurve::LERP},
      Field<PeriodicRippleParams>{"strength", &PeriodicRippleParams::strength,
                                  "Ripple Strength", 0.0f, 0.15f,
                                  FieldCurve::LERP},
      Field<PeriodicRippleParams>{"decay", &PeriodicRippleParams::decay,
                                  "Ripple Decay", 0.0f, 5.0f, FieldCurve::LERP},
      Field<PeriodicRippleParams>{"thickness", &PeriodicRippleParams::thickness,
                                  "Ripple Thickness", 0.7f, 1.4f,
                                  FieldCurve::LOG_POSITIVE},
      Field<PeriodicRippleParams>{"center-azimuth",
                                  &PeriodicRippleParams::center_azimuth,
                                  "Ripple Center Azimuth", 0.0f, TWO_PI_F,
                                  FieldCurve::SHORTEST_PERIODIC},
      Field<PeriodicRippleParams>{
          "center-polar", &PeriodicRippleParams::center_polar,
          "Ripple Center Polar", 0.0f, PI_F, FieldCurve::LERP},
  };
};
static_assert(field_ids_unique<PeriodicRippleParams>());

/** @brief This frame's point on the displacement field's closed loop. */
struct PreparedLoop {
  Vector loop_offset;
};

/** @brief Resolves this frame's loop point from the loop phase. */
HS_FLASH_INLINE inline PreparedLoop prepare(float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  return {Vector(NOISE_LOOP_RADIUS * cosf(angle),
                 NOISE_LOOP_RADIUS * sinf(angle), 0.0f)};
}

/** @brief Loop point plus the steering frame the direct displacement reads. */
struct PreparedDirect {
  Vector loop_offset;
  float direction_cos;
  float direction_sin;
};

/** @brief Prepared parameters for one frame of a periodic ripple. */
struct PreparedRipple {
  Animation::RippleParams ripple;
};

/** @brief Resolves a seamless ripple cycle for ripple_transform(). */
HS_FLASH_INLINE inline PreparedRipple
prepare_ripple(const PeriodicRippleParams &params, float cycle) {
  Animation::RippleParams ripple;
  ripple.center =
      Vector::from_spherical(params.center_azimuth, params.center_polar);
  ripple.amplitude = params.strength;
  ripple.decay = params.decay;
  ripple.thickness = params.thickness;
  const float progress = wrap_t(cycle);
  ripple.phase = progress * PI_F;
  const float attack = std::min(progress * 10.0f, 1.0f);
  ripple.amplitude *= attack * attack * (1.0f - progress);
  ripple.sync();
  return {ripple};
}

/** @brief Applies a prepared periodic ripple and reports its angular travel. */
__attribute__((always_inline)) inline SurfaceResult
periodic_ripple(const Vector &input, const PreparedRipple &prepared,
                bool path_length_required) {
  const Vector displaced = ripple_transform(input, prepared.ripple);
  if (!path_length_required ||
      (displaced.x == input.x && displaced.y == input.y &&
       displaced.z == input.z))
    return {displaced, 0.0f};
  return {displaced, fast_acos(hs::clamp(dot(input, displaced), -1.0f, 1.0f))};
}

/** @brief Resolves the loop point and steering frame for DirectNoise. */
HS_FLASH_INLINE inline PreparedDirect prepare_direct(float phase,
                                                     float direction) {
  const float angle = TWO_PI_F * direction;
  return {prepare(phase).loop_offset, cosf(angle), sinf(angle)};
}

enum class Integrator : uint8_t { EULER, MIDPOINT, MIDPOINT_2X };

struct Euler {
  static constexpr Integrator VALUE = Integrator::EULER;
};

struct Midpoint {
  static constexpr Integrator VALUE = Integrator::MIDPOINT;
};

struct Midpoint2x {
  static constexpr Integrator VALUE = Integrator::MIDPOINT_2X;
};

__attribute__((always_inline)) inline float path_length(const Vector &step,
                                                        bool required) {
  return required ? sqrtf(dot(step, step)) : 0.0f;
}

__attribute__((always_inline)) inline SurfaceResult
finish_step(const Vector &input, const Vector &step,
            bool path_length_required) {
  return {sphere_exp_map_half_radian(input, step),
          path_length(step, path_length_required)};
}

__attribute__((always_inline)) inline SurfaceResult
direct_noise(const Vector &input, const FastNoiseLite &noise,
             ::NoiseBasis basis, float scale, const Vector &loop_offset,
             float strength, float direction_cos, float direction_sin,
             bool path_length_required) {
  if (strength == 0.0f)
    return {input, 0.0f};
  const Vector q = noise_sphere_coordinate(input, scale, loop_offset);
  const Vector tangent = sample_direct_tangent(noise, basis, q, input,
                                               direction_cos, direction_sin);
  return finish_step(input, strength * tangent, path_length_required);
}

__attribute__((always_inline)) inline Vector
curl_field(const Vector &input, const FastNoiseLite &noise, ::NoiseBasis basis,
           float scale, const Vector &loop_offset) {
  const Vector q = noise_sphere_coordinate(input, scale, loop_offset);
  return sample_curl_tangent(noise, basis, q, input);
}

__attribute__((always_inline)) inline SurfaceResult
curl_midpoint_step(const Vector &input, const FastNoiseLite &noise,
                   ::NoiseBasis basis, float scale, const Vector &loop_offset,
                   float distance, bool path_length_required) {
  const Vector first = curl_field(input, noise, basis, scale, loop_offset);
  const Vector midpoint =
      sphere_exp_map_half_radian(input, 0.5f * distance * first);
  const Vector midpoint_field =
      curl_field(midpoint, noise, basis, scale, loop_offset);
  return finish_step(
      input, distance * parallel_transport(midpoint, input, midpoint_field),
      path_length_required);
}

HS_FLASH_INLINE inline SurfaceResult
curl_noise(const Vector &input, const FastNoiseLite &noise, ::NoiseBasis basis,
           Integrator integrator, float scale, const Vector &loop_offset,
           float strength, bool path_length_required) {
  if (strength == 0.0f)
    return {input, 0.0f};
  if (integrator == Integrator::EULER)
    return finish_step(
        input, strength * curl_field(input, noise, basis, scale, loop_offset),
        path_length_required);
  if (integrator == Integrator::MIDPOINT)
    return curl_midpoint_step(input, noise, basis, scale, loop_offset, strength,
                              path_length_required);
  const SurfaceResult first =
      curl_midpoint_step(input, noise, basis, scale, loop_offset,
                         0.5f * strength, path_length_required);
  SurfaceResult second =
      curl_midpoint_step(first.sphere, noise, basis, scale, loop_offset,
                         0.5f * strength, path_length_required);
  second.path_length += first.path_length;
  return second;
}

template <typename State, ::NoiseBasis Basis>
struct DirectNoise : ApproximationDefaults {
  using FrameState = typename State::FrameState;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::strength(frame) } -> std::same_as<float>;
        { State::prepare(frame).loop_offset } -> std::convertible_to<Vector>;
        { State::prepare(frame).direction_cos } -> std::convertible_to<float>;
        { State::prepare(frame).direction_sin } -> std::convertible_to<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  __attribute__((always_inline)) static SurfaceResult
  apply(const Vector &input, const FrameState &frame,
        const Prepared &prepared) {
    return direct_noise(input, State::noise(frame), Basis, State::scale(frame),
                        prepared.loop_offset, State::strength(frame),
                        prepared.direction_cos, prepared.direction_sin,
                        State::path_length_required(frame));
  }
};

template <typename State, ::NoiseBasis Basis, typename IntegratorPolicy>
struct CurlNoise : ApproximationDefaults {
  using FrameState = typename State::FrameState;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::strength(frame) } -> std::same_as<float>;
        { State::prepare(frame).loop_offset } -> std::convertible_to<Vector>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } && requires {
        { IntegratorPolicy::VALUE } -> std::convertible_to<Integrator>;
      };

  using Prepared = std::remove_cvref_t<decltype(State::prepare(
      std::declval<const FrameState &>()))>;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return State::prepare(frame);
  }

  HS_FLASH_INLINE static SurfaceResult apply(const Vector &input,
                                             const FrameState &frame,
                                             const Prepared &prepared) {
    return curl_noise(input, State::noise(frame), Basis,
                      IntegratorPolicy::VALUE, State::scale(frame),
                      prepared.loop_offset, State::strength(frame),
                      State::path_length_required(frame));
  }
};

template <typename State> struct PeriodicRipple : ApproximationDefaults {
  using FrameState = typename State::FrameState;

  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::params(frame) } -> std::same_as<const PeriodicRippleParams &>;
        { State::phase(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  using Prepared = PreparedRipple;

  HS_FLASH_INLINE static Prepared prepare(const FrameState &frame) {
    return prepare_ripple(State::params(frame), State::phase(frame));
  }

  HS_FLASH_INLINE static SurfaceResult apply(const Vector &input,
                                             const FrameState &frame,
                                             const Prepared &prepared) {
    return periodic_ripple(input, prepared, State::path_length_required(frame));
  }
};

} // namespace Surface

} // namespace Pullback
