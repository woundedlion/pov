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
  for (const auto &geometry : {DIAMOND, HEXAGONAL, RHOMBIC}) {
    const float LENGTH = geometry.count == 16   ? sqrtf(3.f) * .25f
                         : geometry.count == 32 ? sqrtf(3.f) * .25f
                                                : 1.f;
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
  alignas(Pixel) std::array<uint8_t, BakedPalette::required_arena_bytes()>
      buffer;
  Arena arena(buffer.data(), buffer.size());
  BakedPaletteStorage palette;
  palette.bake(arena, Palette{});
  Raycast::Appearance appearance;
  appearance.palette = &palette.view();
  appearance.inv_far = .1f;
  Raycast::TraceLimits limits;
  limits.max_candidates = 10000;
  Raycast::PreparedCamera camera;
  camera.interval = {0, 2};
  for (auto kind : {Kind::DIAMOND, Kind::HEXAGONAL, Kind::RHOMBIC}) {
    const Geometry GEOMETRY(kind);
    const auto MID = (GEOMETRY.edges[0].a + GEOMETRY.edges[0].b) * .5f;
    camera.center = {MID.x, MID.y, MID.z - 1, 0};
    const auto FIRST = shade(GEOMETRY, 1, .04f, camera, {.002f, 0}, limits,
                             appearance, math::Z_AXIS);
    HS_EXPECT_TRUE(FIRST.trace.has_surface);
    HS_EXPECT_FALSE(FIRST.trace.contribution.verified);
    HS_EXPECT_GT(FIRST.color.alpha, 0.f);
    for (int axis = 0; axis < 3; ++axis) {
      auto shifted = camera;
      const float PERIOD[] = {GEOMETRY.period.x, GEOMETRY.period.y,
                              GEOMETRY.period.z};
      shifted.center[axis] += 2 * PERIOD[axis];
      const auto TRANSLATED = shade(GEOMETRY, 1, .04f, shifted, {.002f, 0},
                                    limits, appearance, math::Z_AXIS);
      HS_EXPECT_NEAR(FIRST.color.alpha, TRANSLATED.color.alpha, 1e-5f);
    }
    auto bounded = limits;
    bounded.max_steps = 0;
    const auto EXHAUSTED =
        shade(GEOMETRY, 1, .04f, camera, {}, bounded, appearance, math::Z_AXIS);
    HS_EXPECT_EQ(EXHAUSTED.trace.status,
                 Raycast::TraceStatus::BUDGET_EXHAUSTED);
    bounded = limits;
    bounded.max_candidates = 0;
    const auto CANDIDATE_LIMIT =
        shade(GEOMETRY, 1, .04f, camera, {}, bounded, appearance, math::Z_AXIS);
    HS_EXPECT_EQ(CANDIDATE_LIMIT.trace.status,
                 Raycast::TraceStatus::BUDGET_EXHAUSTED);
    auto invalid = camera;
    invalid.domain = Raycast::SamplingDomain::SLICE_4D;
    HS_EXPECT_EQ(
        shade(GEOMETRY, 1, .04f, invalid, {}, limits, appearance, math::Z_AXIS)
            .trace.status,
        Raycast::TraceStatus::INVALID_QUERY);
  }
}

inline void run() {
  test_geometry();
  test_box_slabs();
  test_traversal();
}

} // namespace hs_test::cellular_wire_tests
