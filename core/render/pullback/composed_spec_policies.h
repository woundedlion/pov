/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>

/** @file composed_spec_policies.h
 * @brief Typed runtime policies a composed effect's `Spec` selects.
 */

namespace Pullback {

/** @brief Column-time color sampling of a composed effect. */
enum class TemporalColor : uint8_t {
  OFF,   ///< One frame state per draw; no column-time sampling.
  PHASE, ///< Column-time sampling; eligible color sources sample one phase
         ///< per column after palette interpolation.
  FUSED, ///< Column-time sampling; an eligible final temporal color layer is
         ///< fused.
};

/** @brief Interval palette caches of a column-time sampled effect. */
enum class PaletteCache : uint8_t {
  NONE,             ///< No interval caches.
  LAB,              ///< Interval OKLab palette bases.
  LAB_RGB,          ///< OKLab bases and final hue colors over them.
  LAB_RGB_MIDPOINT, ///< OKLab bases and one admitted final-hue midpoint row.
};

/** @brief Direct-noise evaluation of phase-sampled color sources. */
enum class PhaseNoise : uint8_t {
  DIRECT,        ///< Evaluates direct noise at every sample.
  CHAMBER_CACHE, ///< Caches direct noise over the dodecahedral chamber at
                 ///< temporal endpoints.
};

/**
 * @brief Column-time sampling, color fusion and interval cache policy of a
 *        composed effect.
 * @details Built only by off(), phase() and fused(), so interval caches exist
 * only with column-time sampling and the chamber noise cache only with phase
 * sampling.
 */
class TemporalPolicy {
public:
  /** @brief One frame state per draw.
   * @return The policy without column-time sampling or caches. */
  static consteval TemporalPolicy off() {
    return {TemporalColor::OFF, PaletteCache::NONE, PhaseNoise::DIRECT};
  }
  /** @brief Column-time sampling with per-column color phases.
   * @param cache Interval palette caches.
   * @param noise Direct-noise evaluation of phase-sampled sources.
   * @return The phase-sampling policy. */
  static consteval TemporalPolicy
  phase(PaletteCache cache = PaletteCache::LAB_RGB,
        PhaseNoise noise = PhaseNoise::DIRECT) {
    return {TemporalColor::PHASE, cache, noise};
  }
  /** @brief Column-time sampling with a fused final color layer.
   * @param cache Interval palette caches.
   * @return The fused-color policy. */
  static consteval TemporalPolicy
  fused(PaletteCache cache = PaletteCache::LAB) {
    return {TemporalColor::FUSED, cache, PhaseNoise::DIRECT};
  }

  /** @brief Whether columns sample consecutive frame states.
   * @return True unless the policy is off(). */
  constexpr bool sampling() const { return color_mode != TemporalColor::OFF; }
  /** @brief Column-time color sampling.
   * @return The color mode. */
  constexpr TemporalColor color() const { return color_mode; }
  /** @brief Interval palette caches.
   * @return The cache selection. */
  constexpr PaletteCache cache() const { return palette_cache; }
  /** @brief Direct-noise evaluation of phase-sampled sources.
   * @return The noise evaluation. */
  constexpr PhaseNoise noise() const { return phase_noise; }

private:
  constexpr TemporalPolicy(TemporalColor color, PaletteCache cache,
                           PhaseNoise noise)
      : color_mode(color), palette_cache(cache), phase_noise(noise) {}

  TemporalColor color_mode;   ///< Column-time color sampling.
  PaletteCache palette_cache; ///< Interval palette caches.
  PhaseNoise phase_noise;     ///< Direct-noise evaluation.
};

/** @brief Motion of a composed effect's projection frame. */
enum class ProjectionMotion : uint8_t {
  STATIC,          ///< Identity projection frame.
  ANIMATED,        ///< Random-walk wander and Y spin on the identity base.
  ANIMATED_TILTED, ///< Wander and spin on projection_frame_base(true);
                   ///< authored as the `spin-wander-tilted` frame.
};

/** @brief Immutable dodecahedral tables copied into persistent storage. */
struct DodecahedralCache {
  bool matrices = false; ///< Copies the reflection matrices.
  bool faces = false;    ///< Copies the cube-face classifiers.
};

} // namespace Pullback
