/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/pullback/contract.h"

namespace Pullback {

/** @brief Frame-prepared spherical ray rendering without a new carrier. */
template <typename Renderer>
struct RayStage : Stage::Contract<RayStage<Renderer>, SphereSample, Color4> {
  using Policies = std::tuple<>;

  template <typename Binding>
  static auto prepare(const typename Binding::FrameState &frame) {
    return Renderer::prepare(frame);
  }

  template <typename Binding, typename Prepared>
  __attribute__((always_inline)) static Color4
  run(const SphereSample &input, const typename Binding::FrameState &frame,
      const Prepared &prepared) {
    return Renderer::shade(input.dir, frame, prepared);
  }
};

} // namespace Pullback
