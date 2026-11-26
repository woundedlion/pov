/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/filter.h.
 */
#pragma once

/** @file
 * @brief Shared filter test fixtures. */
#include <functional>

#include "tests/vec_test_util.h"

#include "core/animation/orientation.h"
#include "core/math/mobius.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstddef>
#include <span>
#include <type_traits>

#include "core/render/filter.h"
#include "core/render/filter/feedback_warp_cache.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/render/canvas.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/temporal_oracle.h"
#include <vector>

namespace hs_test::filter_tests {

template <typename P>
concept RawFramePlotter = requires(P &pipeline, Canvas &canvas) {
  pipeline.plot(canvas, 0, 0, Pixel{}, 0.0f, 1.0f);
};

template <typename P>
concept RawTimedPlotter = requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_at(canvas, math::X_AXIS, Pixel{}, 0.0f, 1.0f, 0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_transformed_at(canvas, math::X_AXIS, Pixel{}, 0.0f, 1.0f, 0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_projected_at(canvas, math::PixelCoords{}, Pixel{}, 0.0f, 1.0f,
                             0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_column_at(canvas, math::PixelCoords{}, 0, Pixel{}, 0.0f, 1.0f,
                          0.5f);
} || requires(P &pipeline, Canvas &canvas) {
  pipeline.plot_screen_at(canvas, 0.0f, 0.0f, Pixel{}, 0.0f, 1.0f, 0.5f);
};

template <typename P>
concept TerminalFlusher =
    requires(P &pipeline, Canvas &canvas) { pipeline.flush(canvas, 1.0f); };

template <typename P>
concept ReplacementFrameStarter = requires(P &pipeline, Canvas &canvas) {
  pipeline.begin_frame(canvas, 1.0f);
};
} // namespace hs_test::filter_tests
