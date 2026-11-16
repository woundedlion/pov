/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/scan.h using framebuffer readback, recording
 * callbacks, and pixel sinks.
 */
#pragma once

/** @file
 * @brief Shared scan test fixtures. */
#include <limits>
#include <functional>
#include <optional>
#include <array>
#include <cstdio>

#include "core/render/scan.h"
#include "core/render/scan/mesh_motion.h"
#include "core/render/sdf/volume.h"
#include "core/render/plot.h"
#include "core/render/filter.h"
#include "core/render/canvas.h"
#include "core/render/spacetime.h"
#include "core/math/geometry.h"
#include "tests/mesh_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/prepared_face_test_access.h"
#include "tests/volume_reference.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/vec_test_util.h"

#include <cfloat>
#include <cstring>
#include <memory>
#include <vector>

namespace hs_test::scan_tests {

template <int N, int LUT_N>
inline void build_ring_knots(float (&knots)[N][LUT_N + 1]) {
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < LUT_N; ++k) {
      const float t = 2.0f * math::PI_F * k / LUT_N;
      knots[i][k] = 0.06f * sinf((i + 2) * t) + 0.03f * cosf(3.0f * t + i);
    }
    knots[i][LUT_N] = knots[i][0];
  }
}

/**
 * @brief Pins Render::pole_lod_aggressiveness for a scope and restores it on exit.
 */
struct ScopedPoleLod {
  float saved; /**< Value in force at construction. */
  /**
   * @brief Pins the knob.
   * @param value Aggressiveness held for the scope.
   */
  explicit ScopedPoleLod(float value) : saved(Render::pole_lod_aggressiveness) {
    Render::pole_lod_aggressiveness = value;
  }
  ~ScopedPoleLod() { Render::pole_lod_aggressiveness = saved; }
  ScopedPoleLod(const ScopedPoleLod &) = delete;
  ScopedPoleLod &operator=(const ScopedPoleLod &) = delete;
};

/** @brief Contiguous ring storage that destroys every constructed element. */
template <size_t N> class DistortedRingArray {
public:
  DistortedRingArray() = default;
  DistortedRingArray(const DistortedRingArray &) = delete;
  DistortedRingArray &operator=(const DistortedRingArray &) = delete;
  ~DistortedRingArray() {
    while (count)
      data()[--count].~DistortedRing();
  }
  template <typename... Args> void emplace(size_t i, Args &&...args) {
    HS_CHECK(i == count && count < N, "ring fixture construction order");
    new (storage + i * sizeof(SDF::DistortedRing))
        SDF::DistortedRing(std::forward<Args>(args)...);
    ++count;
  }
  SDF::DistortedRing *data() {
    return reinterpret_cast<SDF::DistortedRing *>(storage);
  }

private:
  alignas(
      SDF::DistortedRing) unsigned char storage[N * sizeof(SDF::DistortedRing)];
  size_t count = 0;
};

/**
 * @brief Renders a draw into a fresh stub frame.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param draw Callable `(Pipeline<W, H> &, Canvas &)`.
 * @return Displayed frame.
 */
template <int W, int H, typename DrawFn>
inline std::vector<Pixel> render_stub_frame(DrawFn &&draw) {
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(fx);
    draw(pipeline, canvas);
  }
  fx.advance_display();
  std::vector<Pixel> out;
  capture_frame<W, H>(fx, out);
  return out;
}

/** @brief Applies the common partial scan clip without a margin. */
inline void apply_partial_clip(StubEffect &effect) {
  effect.set_clip(9, 53, 17, 81);
  effect.set_margin(0);
}
} // namespace hs_test::scan_tests
