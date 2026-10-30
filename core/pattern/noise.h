/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file noise.h
 * @brief Explicit generator and octave noise policies. */

#include "core/math/noise_field.h"
#include "core/pattern/field.h"

namespace Pattern {

/**
 * @brief Octave noise sample with contrast sharpening.
 * @param noise Noise generator.
 * @param basis Octave basis.
 * @param coordinate Noise-space sample point.
 * @param contrast Nonnegative sharpening; 0 leaves the sample unchanged.
 * @return Sharpened sample in [-1, 1].
 */
HS_O3_FN inline float noise_contour(const FastNoiseLite &noise,
                                    math::NoiseBasis basis,
                                    const math::Vector &coordinate,
                                    float contrast) {
  const float sample = hs::clamp(
      math::sample_noise_octaves(noise, basis, coordinate), -1.0f, 1.0f);
  return sample * (1.0f + contrast) / (1.0f + contrast * fabsf(sample));
}

/** @brief Sphere noise with explicit formula selection and an owned generator. */
struct SphereNoise {
  using Coordinate = SphereDirection; ///< Admitted unit direction.
  using Params = NoiseParams;         ///< Noise recipe controls.
  using State = PhaseState;           ///< Noise loop phase in turns.
  /** @brief Immutable generator and prepared sample controls. */
  struct Prepared {
    FastNoiseLite noise;      ///< Owned configured generator.
    math::Vector loop_offset; ///< Prepared noise loop point.
    float scale;              ///< Noise-space coordinate scale.
    float contrast;           ///< Nonnegative contour sharpening.
    NoiseBasis basis;         ///< Octave formula.
    NoiseVariant variant;     ///< Sampling policy.
  };
  /** @brief Checks the recipe's finite, admitted controls.
   * @param params Source controls.
   * @return Whether preparation admits the recipe. */
  static bool valid(const Params &params) {
    return finite(params.scale) && params.scale >= 1.0f / 64.0f &&
           params.scale <= 64.0f && finite(params.speed) &&
           finite(params.contrast) && params.contrast >= 0.0f &&
           params.contrast <= 8.0f && params.basis <= NoiseBasis::RIDGED3 &&
           params.variant <= NoiseVariant::OCTAVE_CONTOUR &&
           (params.variant != NoiseVariant::CONFIGURED_SINGLE ||
            (params.basis == NoiseBasis::SIMPLEX && params.contrast == 0.0f));
  }
  /** @brief Initializes the loop clock.
   * @param params Source controls.
   * @return Zero loop phase. */
  static State init_state(const Params &params) {
    (void)params;
    return {};
  }
  /** @brief Advances one update in turns.
   * @param state Previous phase.
   * @param params Source controls.
   * @return Wrapped loop phase. */
  static State advance(const State &state, const Params &params) {
    HS_CHECK(valid(params) && finite(state.phase_turns),
             "pattern noise: invalid state");
    return {math::wrap_t(math::wrap_t(state.phase_turns) +
                         std::fmod(params.speed, 1.0f))};
  }
  /** @brief Copies the configured generator and resolves the loop offset.
   * @param params Source controls.
   * @param state Explicit loop phase.
   * @param noise Configured generator, copied for endpoint retention.
   * @return Complete immutable sample inputs. */
  static Prepared prepare(const Params &params, const State &state,
                          const FastNoiseLite &noise) {
    HS_CHECK(valid(params) && finite(state.phase_turns),
             "pattern noise: invalid parameters");
    return {noise,        math::noise_sphere_loop_offset(state.phase_turns),
            params.scale, params.contrast,
            params.basis, params.variant};
  }
  /** @brief Samples the selected configured or octave formula.
   * @param coordinate Admitted unit direction.
   * @param prepared Complete immutable generator and controls.
   * @return Signed scalar. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    const math::Vector q = math::noise_sphere_coordinate(
        coordinate.value(), prepared.scale, prepared.loop_offset);
    if (prepared.variant == NoiseVariant::CONFIGURED_SINGLE)
      return prepared.noise.GetNoiseSingle(q.x, q.y, q.z);
    return noise_contour(prepared.noise, prepared.basis, q, prepared.contrast);
  }
};

/** @brief Advances a noise loop by its declared turn rate.
 * @param phase Previous loop phase in turns.
 * @param params Noise source controls.
 * @return Wrapped next phase. */
inline float advance_noise_phase(float phase, const NoiseSourceParams &params) {
  return math::wrap_t(phase + params.noise_time_rate);
}

/** @brief Pure projected-plane octave noise with a turn-valued loop clock. */
struct PlaneNoise {
  using Coordinate = PlanePosition; ///< Consumer-supplied plane coordinate.
  using Params = NoiseSourceParams; ///< Canonical contour controls.
  using State = PhaseState;         ///< Noise loop phase in turns.
  /** @brief Complete immutable generator, loop point and shaping inputs. */
  struct Prepared {
    FastNoiseLite noise;      ///< Owned configured generator.
    math::Vector loop_offset; ///< Prepared projected loop point.
    Params params;            ///< Spatial scale, contrast and loop rate.
    math::NoiseBasis basis;   ///< Explicit octave formula.
  };
  /** @brief Checks the contour controls.
   * @param params Source controls.
   * @return Whether the controls are admitted. */
  static bool valid(const Params &params) {
    return finite(params.noise_scale) && params.noise_scale >= 1.0f / 64.0f &&
           params.noise_scale <= 64.0f && finite(params.noise_contrast) &&
           params.noise_contrast >= 0.0f && params.noise_contrast <= 8.0f &&
           finite(params.noise_time_rate) &&
           fabsf(params.noise_time_rate) <= 1.0f / 64.0f;
  }
  /** @brief Initializes the noise-loop phase.
   * @param params Source controls.
   * @return Zero phase in turns. */
  static State init_state(const Params &params) {
    (void)params;
    return {};
  }
  /** @brief Advances the loop by the source's declared turn rate.
   * @param state Previous loop phase.
   * @param params Source controls.
   * @return Wrapped next phase. */
  static State advance(const State &state, const Params &params) {
    HS_CHECK(valid(params) && finite(state.phase_turns),
             "pattern plane noise: invalid advance inputs");
    return {advance_noise_phase(state.phase_turns, params)};
  }
  /** @brief Copies the generator and resolves projected loop construction.
   * @param params Source controls.
   * @param state Explicit noise phase.
   * @param generator Configured generator retained by value.
   * @param basis Explicit octave formula.
   * @return Complete immutable source inputs. */
  static Prepared prepare(const Params &params, const State &state,
                          const FastNoiseLite &generator,
                          math::NoiseBasis basis) {
    HS_CHECK(valid(params) && finite(state.phase_turns) &&
                 basis <= math::NoiseBasis::RIDGED3,
             "pattern plane noise: invalid preparation inputs");
    return {generator, math::noise_projected_loop_offset(state.phase_turns),
            params, basis};
  }
  /** @brief Evaluates the canonical octave contour.
   * @param coordinate Plane position.
   * @param prepared Immutable generator and controls.
   * @return Signed contour sample. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return noise_contour(prepared.noise, prepared.basis,
                         math::noise_projected_coordinate(
                             coordinate.value, prepared.params.noise_scale,
                             prepared.loop_offset),
                         prepared.params.noise_contrast);
  }
};

} // namespace Pattern
