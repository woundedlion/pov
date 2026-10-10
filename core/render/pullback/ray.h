/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"

/** @file ray.h
 * @brief Spherical ray renderer adapter for pullback stage pipelines. */

namespace Pullback {

/** @brief Frame-prepared spherical ray rendering without a new carrier.
 * @tparam Renderer Supplies prepare(frame) and shade(direction, frame,
 * prepared), returning Color4 from shade(). */
template <typename Renderer>
struct RayStage : Stage::Contract<RayStage<Renderer>, SphereSample, Color4> {
  using Policies = std::tuple<>; ///< No stage policies.

  /** @brief Forwards to `Renderer::prepare`.
   * @tparam Binding Pipeline binding supplying `FrameState`.
   * @param frame Per-frame state.
   * @return The renderer's prepared state. */
  template <typename Binding>
  static auto prepare(const typename Binding::FrameState &frame) {
    return Renderer::prepare(frame);
  }

  /** @brief Shades the sample's sphere direction via `Renderer::shade`.
   * @tparam Binding Pipeline binding supplying `FrameState`.
   * @tparam Prepared The renderer's prepared state type.
   * @param input Sphere sample; only `dir` is read.
   * @param frame Per-frame state.
   * @param prepared Result of `prepare` for this frame.
   * @return Shaded colour. */
  template <typename Binding, typename Prepared>
  __attribute__((always_inline)) static Color4
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const Prepared &prepared) {
    return Renderer::shade(input.dir, frame, prepared);
  }
};

} // namespace Pullback
