/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/sdf/cellular_wire.h"
#include "tests/test_harness.h"

namespace hs_test::cellular_wire_tests {

inline void test_geometry() {
  using namespace SDF::CellularWire;
  const Geometry DIAMOND(Kind::DIAMOND);
  const Geometry HEXAGONAL(Kind::HEXAGONAL);
  const Geometry RHOMBIC(Kind::RHOMBIC);
  HS_EXPECT_EQ(DIAMOND.count, 16);
  HS_EXPECT_EQ(HEXAGONAL.count, 10);
  HS_EXPECT_EQ(RHOMBIC.count, 32);
  struct ExpectedGeometry {
    Kind kind;
    float length;
  };
  const ExpectedGeometry CASES[] = {{Kind::DIAMOND, sqrtf(3.f) * .25f},
                                    {Kind::HEXAGONAL, 1.f},
                                    {Kind::RHOMBIC, sqrtf(3.f) * .25f}};
  for (const auto &expected : CASES) {
    const Geometry geometry(expected.kind);
    const float LENGTH = expected.length;
    for (int i = 0; i < geometry.count; ++i) {
      const auto D = geometry.edges[i].b - geometry.edges[i].a;
      HS_EXPECT_NEAR(sqrtf(math::dot(D, D)), LENGTH, 1e-6f);
    }
  }
  const Edge EDGE{{0, 0, 0}, {0, 0, 1}};
  float t;
  HS_EXPECT_NEAR(closest({{-2, .05f, .5f}, math::X_AXIS, {0, 4}}, EDGE, t),
                 .05f, 1e-6f);
  HS_EXPECT_NEAR(t, 2.f, 1e-6f);
  HS_EXPECT_NEAR(closest({{.05f, 0, -2}, math::Z_AXIS, {0, 4}}, EDGE, t), .05f,
                 1e-6f);
  HS_EXPECT_NEAR(t, 2.f, 1e-6f);
  HS_EXPECT_NEAR(closest({{0, 0, -2}, -math::Z_AXIS, {0, 4}}, EDGE, t), 2.f,
                 1e-6f);
}

inline void test_box_slabs() {
  using namespace SDF::CellularWire;
  const math::Vector LO(1, -1, -1);
  const math::Vector HI(2, 1, 1);
  const BoxRay FORWARD({{}, math::X_AXIS, {0, 4}});
  HS_EXPECT_TRUE(box_overlap(FORWARD, LO, HI, 0, 4));
  HS_EXPECT_FALSE(box_overlap(FORWARD, LO, HI, 0, .5f));
  HS_EXPECT_FALSE(box_overlap(FORWARD, LO, HI, 3, 4));
  const BoxRay REVERSE({{3, 0, 0}, -math::X_AXIS, {0, 4}});
  HS_EXPECT_TRUE(box_overlap(REVERSE, LO, HI, 0, 4));
  const BoxRay PARALLEL_OUTSIDE({{0, 2, 0}, math::X_AXIS, {0, 4}});
  HS_EXPECT_FALSE(box_overlap(PARALLEL_OUTSIDE, LO, HI, 0, 4));
  const BoxRay GRAZING({{0, 1, 0}, math::X_AXIS, {0, 4}});
  HS_EXPECT_TRUE(box_overlap(GRAZING, LO, HI, 0, 4));
}

struct Palette {
  Color4 get(float) const { return {{40000, 30000, 20000}, 1}; }
};

inline void test_traversal() {
  using namespace SDF::CellularWire;
  alignas(Pixel)
      std::array<uint8_t, BakedPalette::required_arena_bytes() +
                              sizeof(HitStorage) + alignof(HitStorage)>
          buffer;
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage palette;
  palette.bake(arena, Palette{});
  HitStorage &storage = *arena.allocate_n<HitStorage>(1);
  Raycast::Appearance appearance;
  appearance.palette = &palette.view();
  appearance.inv_far = .1f;
  Raycast::TraceLimits limits;
  limits.max_candidates = 10000;
  Raycast::PreparedCamera camera;
  camera.interval = {0, 2};
  Geometry bounded_geometry;
  bounded_geometry.period = {10, 10, 10};
  bounded_geometry.lower = {0, 0, 0};
  bounded_geometry.upper = {2, 2, 2};
  bounded_geometry.count = 1;
  bounded_geometry.edges[0] = {{1, 1, .5f}, {1.1f, 1, .5f}};
  camera.center = {1.05f, 1, 0, 0};
  const auto INSIDE = shade(bounded_geometry, 1, .04f, camera, {1, 0}, limits,
                            appearance, storage, math::Z_AXIS);
  HS_EXPECT_EQ(INSIDE.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_GT(INSIDE.color.alpha, 0.f);
  bounded_geometry.edges[0] = {{1, 1, 1.2f}, {1.1f, 1, 1.2f}};
  const auto OUTSIDE = shade(bounded_geometry, 1, .04f, camera, {1, 0}, limits,
                             appearance, storage, math::Z_AXIS);
  HS_EXPECT_EQ(OUTSIDE.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(OUTSIDE.color.alpha, 0.f);
  const auto ZERO_ANGLE = shade(bounded_geometry, 1, .5f, camera, {0, 0},
                                limits, appearance, storage, math::Z_AXIS);
  HS_EXPECT_EQ(ZERO_ANGLE.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(ZERO_ANGLE.trace.counters.candidates, 0);
  bounded_geometry.count = 2;
  bounded_geometry.edges[0] = {{1, 1, .8f}, {1.1f, 1, .8f}};
  bounded_geometry.edges[1] = {{1, 1, .2f}, {1.1f, 1, .2f}};
  auto ordered_appearance = appearance;
  ordered_appearance.near_inv_span = .01f;
  const auto ORDERED = shade(bounded_geometry, 1, .04f, camera, {}, limits,
                             ordered_appearance, storage, math::Z_AXIS);
  HS_EXPECT_EQ(ORDERED.trace.status, Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(ORDERED.trace.counters.layers, 2);
  HS_EXPECT_GT(ORDERED.color.alpha, 0.f);
  HS_EXPECT_NEAR(storage.hits[0].t, .2f, 1e-6f);
  auto partial_limits = limits;
  partial_limits.max_candidates = 1;
  const auto PARTIAL =
      shade(bounded_geometry, 1, .04f, camera, {}, partial_limits,
            ordered_appearance, storage, math::Z_AXIS);
  HS_EXPECT_EQ(PARTIAL.trace.status, Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_EQ(PARTIAL.trace.counters.candidates, 1);
  HS_EXPECT_EQ(PARTIAL.trace.counters.layers, 0);
  HS_EXPECT_EQ(PARTIAL.color.alpha, 0.f);

  for (size_t i = 0; i < storage.hits.size(); ++i)
    storage.hits[i] = {10.f + static_cast<float>(i), 1, 0};
  retain_hit(storage, false, {1.f, 1, 0});
  retain_hit(storage, true, {100.f, 1, 0});
  std::sort(storage.hits.begin(), storage.hits.end(),
            [](const Hit &a, const Hit &b) { return a.t < b.t; });
  HS_EXPECT_EQ(storage.hits.front().t, 1.f);
  HS_EXPECT_EQ(storage.hits.back().t,
               10.f + static_cast<float>(storage.hits.size()) - 2.f);
  Geometry layered_geometry;
  layered_geometry.period = {1, 10, 10};
  layered_geometry.lower = {-1, 1, -.7f};
  layered_geometry.upper = {2, 1, 1.1f};
  layered_geometry.count = static_cast<int>(layered_geometry.edges.size());
  for (auto &edge : layered_geometry.edges)
    edge = {{-1, 1, -.7f}, {2, 1, 1.1f}};
  camera.center = {.5f, 1, 0, 0};
  camera.interval = {0, 1};
  const auto NEAREST_OVERFLOW =
      shade(layered_geometry, 1, .04f, camera, {}, limits, ordered_appearance,
            storage, math::Z_AXIS);
  HS_EXPECT_EQ(NEAREST_OVERFLOW.trace.status,
               Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_GT(NEAREST_OVERFLOW.trace.counters.candidates,
               static_cast<int>(storage.hits.size()));
  HS_EXPECT_EQ(NEAREST_OVERFLOW.trace.counters.layers, 1);
  HS_EXPECT_NEAR(NEAREST_OVERFLOW.trace.contribution.t, .2f, 1e-6f);
  HS_EXPECT_GT(NEAREST_OVERFLOW.color.alpha, 0.f);
  layered_geometry.count = 1;
  const auto COMPLETE_LAYERS =
      shade(layered_geometry, 1, .04f, camera, {}, limits, ordered_appearance,
            storage, math::Z_AXIS);
  HS_EXPECT_EQ(COMPLETE_LAYERS.trace.status,
               Raycast::TraceStatus::RANGE_COMPLETE);
  HS_EXPECT_EQ(COMPLETE_LAYERS.trace.counters.layers, 2);
  HS_EXPECT_NEAR(storage.hits[0].t, .2f, 1e-6f);
  HS_EXPECT_NEAR(storage.hits[1].t, .8f, 1e-6f);
  camera.interval = {0, 2};
  const Geometry WIDE(Kind::RHOMBIC);
  camera.center = {0, 0, -.5f, 0};
  auto overflow_appearance = appearance;
  overflow_appearance.near_inv_span = .01f;
  auto overflow_limits = limits;
  overflow_limits.max_layers = 10000;
  const auto HIT_OVERFLOW = shade(WIDE, 1, .45f, camera, {.0396f, 0},
                                  overflow_limits, overflow_appearance, storage,
                                  math::Vector(-11, -14, 1).normalized());
  HS_EXPECT_EQ(HIT_OVERFLOW.trace.status,
               Raycast::TraceStatus::BUDGET_EXHAUSTED);
  HS_EXPECT_GT(HIT_OVERFLOW.color.alpha, 0.f);
  HS_EXPECT_LT(HIT_OVERFLOW.trace.counters.steps, overflow_limits.max_steps);
  HS_EXPECT_LT(HIT_OVERFLOW.trace.counters.candidates,
               overflow_limits.max_candidates);
  HS_EXPECT_LT(HIT_OVERFLOW.trace.counters.layers, overflow_limits.max_layers);
  struct InvalidQuery {
    float cell_size;
    float wire_radius;
    Raycast::Footprint footprint;
    bool palette;
  };
  const InvalidQuery INVALID[] = {{0, .04f, {}, true},
                                  {-1, .04f, {}, true},
                                  {NAN, .04f, {}, true},
                                  {1, 0, {}, true},
                                  {1, -1, {}, true},
                                  {1, NAN, {}, true},
                                  {1, .04f, {-.002f, 0}, true},
                                  {1, .04f, {NAN, 0}, true},
                                  {1, .04f, {INFINITY, 0}, true},
                                  {1, .04f, {.002f, -1}, true},
                                  {1, .04f, {.002f, NAN}, true},
                                  {1, .04f, {.002f, INFINITY}, true},
                                  {1, .04f, {}, false}};
  for (const auto &query : INVALID) {
    auto invalid_appearance = appearance;
    if (!query.palette)
      invalid_appearance.palette = nullptr;
    HS_EXPECT_EQ(shade(WIDE, query.cell_size, query.wire_radius, camera,
                       query.footprint, limits, invalid_appearance, storage,
                       math::Z_AXIS)
                     .trace.status,
                 Raycast::TraceStatus::INVALID_QUERY);
  }
  for (auto kind : {Kind::DIAMOND, Kind::HEXAGONAL, Kind::RHOMBIC}) {
    const Geometry GEOMETRY(kind);
    const auto MID = (GEOMETRY.edges[0].a + GEOMETRY.edges[0].b) * .5f;
    camera.center = {MID.x, MID.y, MID.z - 1, 0};
    const auto FIRST = shade(GEOMETRY, 1, .04f, camera, {.002f, 0}, limits,
                             appearance, storage, math::Z_AXIS);
    HS_EXPECT_TRUE(FIRST.trace.status == Raycast::TraceStatus::RANGE_COMPLETE ||
                   FIRST.trace.status == Raycast::TraceStatus::SATURATED);
    HS_EXPECT_GT(FIRST.trace.counters.layers, 0);
    HS_EXPECT_GT(FIRST.color.alpha, 0.f);
    for (int axis = 0; axis < 3; ++axis) {
      auto shifted = camera;
      const float PERIOD[] = {GEOMETRY.period.x, GEOMETRY.period.y,
                              GEOMETRY.period.z};
      shifted.center[axis] += 2 * PERIOD[axis];
      const auto TRANSLATED = shade(GEOMETRY, 1, .04f, shifted, {.002f, 0},
                                    limits, appearance, storage, math::Z_AXIS);
      HS_EXPECT_NEAR(FIRST.color.alpha, TRANSLATED.color.alpha, 1e-5f);
    }
    auto bounded = limits;
    bounded.max_layers = 0;
    const auto NO_LAYERS = shade(GEOMETRY, 1, .04f, camera, {.002f, 0}, bounded,
                                 appearance, storage, math::Z_AXIS);
    HS_EXPECT_EQ(NO_LAYERS.trace.status,
                 Raycast::TraceStatus::BUDGET_EXHAUSTED);
    HS_EXPECT_EQ(NO_LAYERS.trace.counters.layers, 0);
    HS_EXPECT_EQ(NO_LAYERS.color.alpha, 0.0f);
    bounded = limits;
    bounded.max_steps = 0;
    const auto EXHAUSTED = shade(GEOMETRY, 1, .04f, camera, {}, bounded,
                                 appearance, storage, math::Z_AXIS);
    HS_EXPECT_EQ(EXHAUSTED.trace.status,
                 Raycast::TraceStatus::BUDGET_EXHAUSTED);
    bounded = limits;
    bounded.max_candidates = 0;
    const auto CANDIDATE_LIMIT = shade(GEOMETRY, 1, .04f, camera, {}, bounded,
                                       appearance, storage, math::Z_AXIS);
    HS_EXPECT_EQ(CANDIDATE_LIMIT.trace.status,
                 Raycast::TraceStatus::BUDGET_EXHAUSTED);
    auto invalid = camera;
    invalid.domain = Raycast::SamplingDomain::SLICE_4D;
    HS_EXPECT_EQ(shade(GEOMETRY, 1, .04f, invalid, {}, limits, appearance,
                       storage, math::Z_AXIS)
                     .trace.status,
                 Raycast::TraceStatus::INVALID_QUERY);
  }
}

inline void run_cellular_wire_cases() {
  test_geometry();
  test_box_slabs();
  test_traversal();
}

} // namespace hs_test::cellular_wire_tests
