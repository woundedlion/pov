/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Event tracing, stream merging, and verified supersample filtering.
 */
#pragma once

/** @file
 * @brief Shared ray events test fixtures. */

#include "core/render/ray/shade.h"
#include "core/render/sdf/lattice.h"
#include "core/render/sdf/lattice_trace.h"
#include "tests/test_harness.h"
#include "tests/test_fixture.h"
#include "tests/pixel_test_util.h"

namespace hs_test::ray_event_tests {

struct Streams {
  static constexpr size_t STREAM_COUNT = 5;
  std::array<Raycast::Contribution, STREAM_COUNT> heads{};
  std::array<bool, STREAM_COUNT> live{true, true, true, true, true};
  bool stalled = false;
  bool active(size_t i) const { return live[i]; }
  float distance(size_t i) const { return heads[i].t; }
  Raycast::Contribution candidate(size_t i) const { return heads[i]; }
  void advance(size_t i) { live[i] = stalled; }
};

struct SingleGroupStreams : Streams {
  static constexpr size_t GROUP_CAPACITY = 1;
};
} // namespace hs_test::ray_event_tests
