/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file chamber_noise_cache.h
 * @brief Retained noise samples over a dodecahedral chamber chart. */

#include "core/render/pullback/modulation.h"

namespace Pullback::Modulation {

/** @brief One immutable endpoint sampled on a rectangular chamber chart. */
struct ChamberNoisePatch {
  static constexpr float QUALIFIED_SCALE =
      2.4029531f;                                    ///< Qualified noise scale.
  static constexpr float QUALIFIED_SPEED = 0.00015f; ///< Qualified loop speed.
  static constexpr float MAX_CLOCK_STEP =
      0.00021f; ///< Largest cached clock step.
  static constexpr float MAX_PROJECTION_STEP_SQ =
      4e-6f; ///< Quaternion step bound.

  static constexpr int U_STEPS = 16;             ///< Chart samples in x/y.
  static constexpr int V_STEPS = 24;             ///< Chart samples in z/y.
  static constexpr float U_MIN = -0.002f;        ///< Chart's lower x/y bound.
  static constexpr float U_MAX = 0.385f;         ///< Chart's upper x/y bound.
  static constexpr float V_MIN = -0.002f;        ///< Chart's lower z/y bound.
  static constexpr float V_MAX = 0.622f;         ///< Chart's upper z/y bound.
  std::array<float, U_STEPS * V_STEPS> values{}; ///< Row-major scalar samples.
  Pattern::NoiseParams recipe{};                 ///< Baked source recipe.
  math::Quaternion projection; ///< Endpoint projection orientation.
  float phase = 0.0f;          ///< Endpoint source clock.
  int32_t seed = 0;            ///< Endpoint generator seed.
  bool active = false;         ///< Whether samples have been published.

  /** @brief Tests the source recipe qualified for this approximation.
   * @param config Authored direct-noise controls.
   * @return Whether this source admits chamber caching. */
  static bool eligible(const Config &config) {
    return config.enabled && config.version == 1 &&
           config.seed == HUE_NOISE_SEED &&
           config.source == Pattern::SourceKind::NOISE &&
           config.evaluation == Pattern::Evaluation::DIRECT &&
           config.noise.variant == Pattern::NoiseVariant::CONFIGURED_SINGLE &&
           config.noise.basis == Pattern::NoiseBasis::SIMPLEX &&
           config.noise.contrast == 0.0f &&
           config.noise.scale == QUALIFIED_SCALE &&
           config.noise.speed == QUALIFIED_SPEED;
  }

  /** @brief Tests clock seams and antipodal projection steps.
   * @param a_phase Previous clock in [0, 1).
   * @param a Previous projection.
   * @param b_phase Current clock in [0, 1).
   * @param b Current projection.
   * @return Whether interpolation stays inside the cache step envelope. */
  static bool continuous(float a_phase, const math::Quaternion &a,
                         float b_phase, const math::Quaternion &b) {
    float delta = b_phase - a_phase;
    if (delta > 0.5f)
      delta -= 1.0f;
    if (delta < -0.5f)
      delta += 1.0f;
    const auto difference = a.v - b.v;
    const auto opposite = a.v + b.v;
    const float distance =
        (a.r - b.r) * (a.r - b.r) + math::dot(difference, difference);
    const float antipodal =
        (a.r + b.r) * (a.r + b.r) + math::dot(opposite, opposite);
    return fabsf(delta) <= MAX_CLOCK_STEP &&
           std::min(distance, antipodal) <= MAX_PROJECTION_STEP_SQ;
  }

  /** @brief Exact identity of the source and endpoint projection.
   * @param source Immutable direct source.
   * @param orientation Expected projection.
   * @return Whether this patch still describes that endpoint. */
  bool matches(const Prepared &source,
               const math::Quaternion &orientation) const {
    return active && eligible(source.config) && recipe == source.config.noise &&
           seed == source.config.seed &&
           phase == source.state.phase.phase_turns &&
           projection.r == orientation.r && projection.v.x == orientation.v.x &&
           projection.v.y == orientation.v.y &&
           projection.v.z == orientation.v.z;
  }

  /** @brief Bakes an endpoint without changing its direct source policy.
   * @param source Immutable prepared noise.
   * @param orientation Endpoint projection orientation. */
  HS_FLASH_MEMBER void bake(const Prepared &source,
                            const math::Quaternion &orientation) {
    const auto &noise = std::get<0>(source.source);
    for (int y = 0; y < V_STEPS; ++y) {
      const float V = V_MIN + (V_MAX - V_MIN) * y / (V_STEPS - 1);
      for (int x = 0; x < U_STEPS; ++x) {
        const float U = U_MIN + (U_MAX - U_MIN) * x / (U_STEPS - 1);
        const math::Vector chamber = math::Vector(U, 1.0f, V).normalized();
        const auto coordinate =
            Pattern::SphereDirection::admit(math::rotate(chamber, orientation));
        values[y * U_STEPS + x] = Pattern::bounded_scalar(
            Pattern::SphereNoise::sample(coordinate, noise));
      }
    }
    recipe = source.config.noise;
    seed = source.config.seed;
    phase = source.state.phase.phase_turns;
    projection = orientation;
    active = true;
  }

  /** @brief Bilinear scalar at admitted chart coordinates.
   * @param x Lower horizontal sample index.
   * @param y Lower vertical sample index.
   * @param u Horizontal interpolation weight.
   * @param v Vertical interpolation weight.
   * @return Baked scalar. */
  HS_HOT_INLINE float sample(int x, int y, float u, float v) const {
    const float A =
        hs::lerp(values[y * U_STEPS + x], values[y * U_STEPS + x + 1], u);
    const float B = hs::lerp(values[(y + 1) * U_STEPS + x],
                             values[(y + 1) * U_STEPS + x + 1], u);
    return hs::lerp(A, B, v);
  }
};

/** @brief A borrowed patch with its expected immutable orientation. */
struct ChamberNoiseView {
  const ChamberNoisePatch *patch = nullptr; ///< Retained runtime-owned samples.
  math::Quaternion projection;              ///< Expected endpoint orientation.
};

/** @brief Alternating endpoint storage for the chamber noise approximation. */
struct ChamberNoiseCache {
  std::array<ChamberNoisePatch, 2>
      patches{};                        ///< Retained and current endpoints.
  uint8_t slot = 0;                     ///< Current endpoint slot.
  math::Quaternion observed_projection; ///< Last observed endpoint orientation.
  float observed_phase = 0.0f;          ///< Last observed endpoint clock.
  bool observed = false; ///< Whether the preceding observation was eligible.
  bool accepted = false; ///< Whether that observation admitted a bake.

  /** @brief Reuses a matching endpoint or replaces the noncurrent slot.
   * @param source Immutable direct noise source.
   * @param orientation Endpoint projection orientation.
   * @return Current published samples, or an inactive view. */
  ChamberNoiseView prepare(const Prepared &source,
                           const math::Quaternion &orientation) {
    const float PHASE = source.state.phase.phase_turns;
    if (!ChamberNoisePatch::eligible(source.config) ||
        !(PHASE >= 0.0f && PHASE < 1.0f)) {
      observed = false;
      return {};
    }
    if (patches[slot].matches(source, orientation))
      return {&patches[slot], orientation};
    const bool REPEATED = observed && observed_phase == PHASE &&
                          observed_projection.r == orientation.r &&
                          observed_projection.v.x == orientation.v.x &&
                          observed_projection.v.y == orientation.v.y &&
                          observed_projection.v.z == orientation.v.z;
    if (REPEATED && !accepted)
      return {};
    accepted = !observed || ChamberNoisePatch::continuous(observed_phase,
                                                          observed_projection,
                                                          PHASE, orientation);
    observed_projection = orientation;
    observed_phase = PHASE;
    observed = true;
    if (!accepted)
      return {};
    slot ^= 1;
    patches[slot].bake(source, orientation);
    return {&patches[slot], orientation};
  }
};

/** @brief Samples two endpoint patches at one column's chamber coordinate.
 * @param sphere Projection-local sphere direction.
 * @param to_chamber Column's inverse projection.
 * @param previous Previous endpoint patch.
 * @param current Current endpoint patch.
 * @param fraction Column fraction in [0, 1].
 * @param scalar Receives the interpolated scalar on success.
 * @return Whether the direction lies inside the retained chart. */
HS_HOT_INLINE inline bool sample_chamber_noise(
    const math::Vector &sphere, const math::Quaternion &to_chamber,
    const ChamberNoisePatch &previous, const ChamberNoisePatch &current,
    float fraction, float &scalar) {
  const math::Vector chamber = math::rotate(sphere, to_chamber);
  if (!(chamber.y > 0.0f))
    return false;
  const float INVERSE = 1.0f / chamber.y;
  const float U = chamber.x * INVERSE;
  const float V = chamber.z * INVERSE;
  if (!(U >= ChamberNoisePatch::U_MIN && U <= ChamberNoisePatch::U_MAX &&
        V >= ChamberNoisePatch::V_MIN && V <= ChamberNoisePatch::V_MAX))
    return false;
  const float X = (U - ChamberNoisePatch::U_MIN) *
                  ((ChamberNoisePatch::U_STEPS - 1) /
                   (ChamberNoisePatch::U_MAX - ChamberNoisePatch::U_MIN));
  const float Y = (V - ChamberNoisePatch::V_MIN) *
                  ((ChamberNoisePatch::V_STEPS - 1) /
                   (ChamberNoisePatch::V_MAX - ChamberNoisePatch::V_MIN));
  const int IX = std::min(static_cast<int>(X), ChamberNoisePatch::U_STEPS - 2);
  const int IY = std::min(static_cast<int>(Y), ChamberNoisePatch::V_STEPS - 2);
  scalar = hs::lerp(previous.sample(IX, IY, X - IX, Y - IY),
                    current.sample(IX, IY, X - IX, Y - IY), fraction);
  return true;
}

} // namespace Pullback::Modulation
