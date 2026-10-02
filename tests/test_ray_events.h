/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Event tracing, stream merging, and verified supersample filtering.
 */
#pragma once

#include "core/render/ray/shade.h"
#include "tests/test_cellular_wire.h"
#include "tests/test_lattice_trace.h"
#include "tests/test_harness.h"
#include "tests/test_fixture.h"

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

/** @brief Pins single group capacity. */
inline void test_single_group_capacity() {
  SingleGroupStreams streams;
  for (size_t i = 0; i < streams.STREAM_COUNT; ++i) {
    streams.heads[i].t = 1;
    streams.heads[i].coverage = .1f * static_cast<float>(i + 1);
  }
  int count = 0;
  auto result =
      Raycast::trace_events(streams, {0, 2}, {}, [&](const auto &hit) {
        HS_EXPECT_EQ(hit.coverage, .5f);
        ++count;
        return true;
      });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_EQ(result.counters.candidates, 5);
  streams.live.fill(true);
  streams.heads[1].merge_identity = 1;
  count = 0;
  result = Raycast::trace_events(streams, {0, 2}, {}, [&](const auto &hit) {
    HS_EXPECT_EQ(hit.coverage, .1f);
    ++count;
    return true;
  });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(result.counters.candidates, 2);
  HS_EXPECT_EQ(count, 1);
  streams.live.fill(true);
  streams.heads[1].merge_identity = 0;
  Raycast::TraceLimits limits;
  limits.max_candidates = 3;
  count = 0;
  result = Raycast::trace_events(streams, {0, 2}, limits, [&](const auto &hit) {
    HS_EXPECT_NEAR(hit.coverage, .3f, 1e-7f);
    ++count;
    return true;
  });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(count, 1);
  streams.live.fill(true);
  streams.heads[1].coverage = NAN;
  count = 0;
  result = Raycast::trace_events(streams, {0, 2}, {}, [&](const auto &hit) {
    HS_EXPECT_EQ(hit.coverage, .1f);
    ++count;
    return true;
  });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  HS_EXPECT_EQ(count, 1);
}

/** @brief Pins failure status survives flush. */
inline void test_failure_status_survives_flush() {
  SingleGroupStreams streams;
  streams.live.fill(false);
  streams.live[0] = streams.live[1] = true;
  streams.heads[0].t = streams.heads[1].t = 1;
  streams.heads[1].coverage = NAN;
  auto copy = streams;
  int count = 0;
  auto result = Raycast::trace_events(copy, {0, 2}, {}, [&](const auto &) {
    ++count;
    return false;
  });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  HS_EXPECT_EQ(count, 1);
  HS_EXPECT_EQ(result.counters.layers, 1);
  copy = streams;
  Raycast::TraceLimits limits;
  limits.max_layers = 0;
  result = Raycast::trace_events(copy, {0, 2}, limits,
                                 [](const auto &) { return false; });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  HS_EXPECT_EQ(result.counters.layers, 0);

  struct InvalidAfterAdvance : SingleGroupStreams {
    void advance(size_t i) {
      live[i] = false;
      heads[1].t = NAN;
      live[1] = true;
    }
  } invalid;
  invalid.live.fill(false);
  invalid.live[0] = true;
  invalid.heads[0].t = 1;
  result = Raycast::trace_events(invalid, {0, 2}, {},
                                 [](const auto &) { return false; });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  HS_EXPECT_EQ(result.counters.layers, 1);
}

/** @brief Event tracing orders, merges, bounds, and rejects invalid streams. */
inline void test_event_stream_contracts() {
  Streams streams;
  for (size_t i = 0; i < streams.STREAM_COUNT; ++i) {
    streams.heads[i].t = static_cast<float>(5 - i);
    streams.heads[i].merge_identity = i;
  }
  float previous = 0;
  int count = 0;
  auto result =
      Raycast::trace_events(streams, {0, 6}, {}, [&](const auto &hit) {
        HS_EXPECT_GT(hit.t, previous);
        previous = hit.t;
        ++count;
        return true;
      });
  HS_EXPECT_EQ(count, 5);
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(result.counters.candidates, 5);
  streams.live.fill(true);
  streams.heads[0].t = streams.heads[1].t = streams.heads[2].t = 1;
  streams.heads[0].merge_identity = streams.heads[1].merge_identity = 0;
  streams.heads[0].coverage = .25f;
  streams.heads[1].coverage = .75f;
  streams.heads[2].material = 1;
  count = 0;
  result = Raycast::trace_events(streams, {0, 6}, {}, [&](const auto &hit) {
    if (count == 0)
      HS_EXPECT_EQ(hit.coverage, .75f);
    ++count;
    return true;
  });
  HS_EXPECT_EQ(count, 4);
  streams.live.fill(true);
  Raycast::TraceLimits limits;
  limits.max_candidates = 2;
  count = 0;
  result = Raycast::trace_events(streams, {0, 6}, limits, [&](const auto &) {
    ++count;
    return true;
  });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(result.counters.candidates, 2);
  HS_EXPECT_EQ(count, 1);
  streams.live.fill(true);
  streams.stalled = true;
  result = Raycast::trace_events(streams, {0, 6}, {},
                                 [](const auto &) { return true; });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  streams.stalled = false;
  streams.heads[0].coverage = std::numeric_limits<float>::quiet_NaN();
  streams.live.fill(true);
  result = Raycast::trace_events(streams, {0, 6}, {},
                                 [](const auto &) { return true; });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::INVALID_QUERY);
  streams.heads[0].coverage = 1;
  streams.live.fill(true);
  result = Raycast::trace_events(streams, {0, 6}, {},
                                 [](const auto &) { return false; });
  HS_EXPECT_EQ(result.status, Raycast::TraceStatus::SATURATED);
  streams.live.fill(true);
  result = Raycast::trace_events(streams, {2, 3}, {}, [](const auto &hit) {
    HS_EXPECT_TRUE(hit.t == 2 || hit.t == 3);
    return true;
  });
  HS_EXPECT_EQ(result.counters.layers, 1);
}

/** @brief Verified filtering preserves aggregate counters and unresolved status. */
inline void test_verified_filter_contracts() {
  const std::array<math::Vector, 4> directions{math::X_AXIS, math::Y_AXIS,
                                               math::Z_AXIS, -math::X_AXIS};
  int sample = 0;
  const auto filtered = Raycast::verified_filter(directions, [&](const auto &) {
    Raycast::ShadedTrace trace;
    trace.trace.counters.queries = 2;
    trace.trace.counters.candidates = 3;
    trace.trace.counters.layers = 1;
    trace.trace.has_surface = sample < 2;
    trace.trace.contribution.verified = sample < 2;
    trace.trace.status = sample++ == 3 ? Raycast::TraceStatus::UNRESOLVED
                                       : Raycast::TraceStatus::SURFACE;
    trace.color = {Pixel(100, 200, 300), 1};
    return trace;
  });
  HS_EXPECT_EQ(filtered.color.alpha, .5f);
  HS_EXPECT_EQ(filtered.color.color.r, 100);
  HS_EXPECT_EQ(filtered.trace.counters.queries, 8);
  HS_EXPECT_EQ(filtered.trace.counters.candidates, 12);
  HS_EXPECT_EQ(filtered.trace.counters.layers, 4);
  HS_EXPECT_EQ(filtered.trace.status, Raycast::TraceStatus::UNRESOLVED);
}

inline int run_ray_event_tests() {
  hs_test::ModuleFixture fixture("ray_events");
  cellular_wire_tests::run_cellular_wire_cases();
  lattice_trace_tests::run_lattice_trace_cases();
  test_single_group_capacity();
  test_failure_status_survives_flush();
  test_event_stream_contracts();
  test_verified_filter_contracts();
  return fixture.result();
}

} // namespace hs_test::ray_event_tests
