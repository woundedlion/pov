/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"

/**
 * @file surface.h
 * @brief Surface-displacement policies.
 */

namespace Pullback {

namespace Surface {

/** @brief Surface stage placeholder for a look that carries no displacement. */
struct NoSurfaceParams {};

/**
 * @brief Surface parameters for the sphere-space noise displacements
 *        (Pullback::Surface::CurlNoise and DirectNoise).
 */
struct SurfaceNoiseParams {
  float scale = 1.0f;    /**< Spatial scale of the displacement field. */
  float strength = 0.0f; /**< Displacement distance; 0 skips the stage. */
  float speed = 0.0f;    /**< Per-frame advance of the field's loop phase. */
};

/** @brief Surface stage per-frame state; empty unless the look displaces. */
template <typename SurfaceT> struct PreparedSurface {};
template <> struct PreparedSurface<SurfaceNoiseParams> {
  const FastNoiseLite *noise; /**< The runtime's surface noise field. */
  Vector loop_offset;         /**< This frame's point on the field's loop. */
};

/** @brief Binds the displacement field and this frame's loop point. */
HS_FLASH_INLINE inline PreparedSurface<SurfaceNoiseParams>
prepare(const FastNoiseLite &noise, float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  return {&noise, Vector(NOISE_LOOP_RADIUS * cosf(angle),
                         NOISE_LOOP_RADIUS * sinf(angle), 0.0f)};
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

struct Identity : ExactPolicy {
  template <typename FrameState>
  __attribute__((always_inline)) static SurfaceResult
  apply(const Vector &input, const FrameState &) {
    return {input, 0.0f};
  }
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
      input, distance * transport_tangent(midpoint, input, midpoint_field),
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

template <typename State, ::NoiseBasis Basis> struct DirectNoise : ExactPolicy {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::loop_offset(frame) } -> std::same_as<const Vector &>;
        { State::strength(frame) } -> std::same_as<float>;
        { State::direction_cos(frame) } -> std::same_as<float>;
        { State::direction_sin(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      };

  template <typename FrameState>
  __attribute__((always_inline)) static SurfaceResult
  apply(const Vector &input, const FrameState &frame) {
    return direct_noise(input, State::noise(frame), Basis, State::scale(frame),
                        State::loop_offset(frame), State::strength(frame),
                        State::direction_cos(frame),
                        State::direction_sin(frame),
                        State::path_length_required(frame));
  }
};

template <typename State, ::NoiseBasis Basis, typename IntegratorPolicy>
struct CurlNoise : ExactPolicy {
  template <typename Binding>
  static constexpr bool PROVIDER_VALID =
      Detail::ProviderFor<State, Binding> &&
      requires(const typename Binding::FrameState &frame) {
        { State::noise(frame) } -> std::same_as<const FastNoiseLite &>;
        { State::scale(frame) } -> std::same_as<float>;
        { State::loop_offset(frame) } -> std::same_as<const Vector &>;
        { State::strength(frame) } -> std::same_as<float>;
        { State::path_length_required(frame) } -> std::same_as<bool>;
      } && requires {
        { IntegratorPolicy::VALUE } -> std::convertible_to<Integrator>;
      };

  template <typename FrameState>
  HS_FLASH_INLINE static SurfaceResult apply(const Vector &input,
                                             const FrameState &frame) {
    return curl_noise(input, State::noise(frame), Basis,
                      IntegratorPolicy::VALUE, State::scale(frame),
                      State::loop_offset(frame), State::strength(frame),
                      State::path_length_required(frame));
  }
};

} // namespace Surface

} // namespace Pullback
