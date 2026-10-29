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
  using Coordinate = SphereDirection;
  using Params = NoiseParams;
  using State = PhaseState;
  /** @brief Immutable generator and prepared sample controls. */
  struct Prepared {
    FastNoiseLite noise;
    math::Vector loop_offset;
    float scale;
    float contrast;
    NoiseBasis basis;
    NoiseVariant variant;
  };
  /** @brief Checks the recipe's finite, admitted controls.
   * @param params Source controls.
   * @return Whether preparation admits the recipe. */
  static bool valid(const Params &params) {
    return std::isfinite(params.scale) && params.scale >= 1.0f / 64.0f &&
           params.scale <= 64.0f && std::isfinite(params.speed) &&
           std::isfinite(params.contrast) && params.contrast >= 0.0f &&
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
    HS_CHECK(valid(params) && std::isfinite(state.phase_turns),
             "pattern noise: invalid state");
    return {math::wrap_t(math::wrap_t(state.phase_turns) +
                         math::wrap_t(params.speed))};
  }
  /** @brief Copies the configured generator and resolves the loop offset.
   * @param params Source controls.
   * @param state Explicit loop phase.
   * @param noise Configured generator, copied for endpoint retention.
   * @return Complete immutable sample inputs. */
  static Prepared prepare(const Params &params, const State &state,
                          const FastNoiseLite &noise) {
    HS_CHECK(valid(params) && std::isfinite(state.phase_turns),
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
    return noise_contour(prepared.noise,
                         static_cast<math::NoiseBasis>(prepared.basis), q,
                         prepared.contrast);
  }
};

} // namespace Pattern
