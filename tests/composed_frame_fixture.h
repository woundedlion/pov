/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/pullback/composed_effect.h"

namespace hs_test {

/** @brief Composed frame preparation with controlled parameters and clocks. */
struct ComposedFrameWhiteBox {
  template <typename FX>
  static void set_params(FX &effect, const typename FX::Params &params) {
    effect.adopt_params(params);
  }

  template <typename FX> static void advance(FX &effect) {
    effect.advance_runtime();
    effect.update_spatial_frames();
  }

  template <typename FX> static typename FX::FrameState frame(FX &effect) {
    effect.update_spatial_frames();
    return effect.prepare_frame();
  }
};

} // namespace hs_test
