/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "effects/HyperLattice.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace hyper_lattice_tests {

namespace HL = HyperLatticeDetail;

inline void test_periodic_distance() {
  HS_EXPECT_EQ(HL::periodic_distance(0.0f), 0.0f);
  HS_EXPECT_EQ(HL::periodic_distance(1.0f), 0.0f);
  HS_EXPECT_NEAR(HL::periodic_distance(-1.25f), 0.25f, 1e-6f);
  HS_EXPECT_NEAR(HL::periodic_distance(3.5f), 0.5f, 1e-6f);
}

inline void test_edge_metrics() {
  const HL::Vec4 point{{0.0f, 0.02f, 0.4f, 0.1f}};
  const HL::EdgeMetric cubic = HL::edge_metric_3d(point, 0);
  HS_EXPECT_NEAR(cubic.distance_sq, 0.0004f, 1e-7f);
  HS_EXPECT_EQ(cubic.free_axis, uint8_t(2));

  const HL::EdgeMetric hyper = HL::edge_metric_4d(point, 0);
  HS_EXPECT_NEAR(hyper.distance_sq, 0.0104f, 1e-6f);
  HS_EXPECT_EQ(hyper.free_axis, uint8_t(2));
}

inline void test_so4_rotation() {
  HL::FrameState frame{};
  frame.params.dimension = 1.0f;
  frame.params.far_cells = 8.0f;
  frame.rotation_phase[3] = 0.5f * PI_F;
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  const HL::Vec4 rotated =
      prepared.world_to_lattice.apply({{1.0f, 0.0f, 0.0f, 0.0f}});
  HS_EXPECT_NEAR(rotated[0], 0.0f, 2e-4f);
  HS_EXPECT_NEAR(rotated[3], 1.0f, 2e-4f);
  float norm_sq = 0.0f;
  for (int axis = 0; axis < HL::DIMENSIONS; ++axis)
    norm_sq += rotated[axis] * rotated[axis];
  HS_EXPECT_NEAR(norm_sq, 1.0f, 2e-4f);

  frame.params.dimension = 0.0f;
  const HL::Vec4 cubic = HL::prepare_trace(frame).world_to_lattice.apply(
      {{1.0f, 0.0f, 0.0f, 0.0f}});
  HS_EXPECT_EQ(cubic[3], 0.0f);
}

inline void test_reflection_convention() {
  const HL::Vec4 center =
      HL::reflected_direction(X_AXIS, HL::ReflectionMode::CHROME);
  HS_EXPECT_EQ(center[0], 1.0f);
  HS_EXPECT_EQ(center[1], 0.0f);
  HS_EXPECT_EQ(center[2], 0.0f);

  const HL::Vec4 rim =
      HL::reflected_direction(Y_AXIS, HL::ReflectionMode::CHROME);
  HS_EXPECT_EQ(rim[0], -1.0f);
  HS_EXPECT_EQ(rim[1], 0.0f);
  const HL::Vec4 radial =
      HL::reflected_direction(Y_AXIS, HL::ReflectionMode::RADIAL);
  HS_EXPECT_EQ(radial[1], 1.0f);
}

inline void test_next_plane_is_strict() {
  HL::FrameState frame{};
  frame.params.far_cells = 8.0f;
  frame.origin = {{0.0f, 0.25f, 0.5f, 0.75f}};
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  HS_EXPECT_EQ(prepared.positive_offset[0], 1.0f);
  HS_EXPECT_EQ(prepared.negative_offset[0], 1.0f);
  HS_EXPECT_EQ(prepared.positive_offset[1], 0.75f);
  HS_EXPECT_EQ(prepared.negative_offset[1], 0.25f);
}

inline void test_trace_carrier() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.origin = {{0.25f, 0.0f, 0.31f, 0.43f}};
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  const Pullback::SphereSample input{X_AXIS, 0.25f};
  const Pullback::FieldSample output =
      HL::TraceStage::Bind<HL::Binding>::run(input, frame, prepared);
  HS_EXPECT_GT(output.coverage, 0.0f);
  HS_EXPECT_LE(output.coverage, 1.0f);
  HS_EXPECT_GE(output.value, 0.0f);
  HS_EXPECT_LE(output.value, 1.0f);
  HS_EXPECT_GT(output.path_length, input.path_length);
  HS_EXPECT_EQ(output.sphere.x, input.dir.x);
}

inline void test_hyperplane_event() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(3);
  frame.params.reflection = HL::ReflectionMode::RADIAL;
  frame.origin = {{0.0f, 0.0f, 0.31f, 0.25f}};
  frame.rotation_phase[3] = 0.5f * PI_F;
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  const HL::TraceHit hit = HL::trace(X_AXIS, prepared);
  HS_EXPECT_GT(hit.coverage, 0.0f);
  HS_EXPECT_EQ(hit.free_axis, uint8_t(2));
  HS_EXPECT_NEAR(hit.distance, 0.75f, 3e-4f);
}

inline void test_presets_and_pipeline() {
  using Effect = HyperLattice<96, 20>;
  static_assert(HL::RenderPipeline::STAGE_COUNT == 2);
  static_assert(HL::RenderPipeline::Validation::MONOTONE);
  static_assert(HL::RenderPipeline::Validation::ENTRY);
  static_assert(HL::RenderPipeline::Validation::EXIT);
  for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index)
    HS_EXPECT_TRUE(Effect::valid_params(Effect::preset_params(index)));
}

inline int run_hyper_lattice_tests() {
  hs_test::ModuleFixture fixture("hyper_lattice");
  test_periodic_distance();
  test_edge_metrics();
  test_so4_rotation();
  test_reflection_convention();
  test_next_plane_is_strict();
  test_trace_carrier();
  test_hyperplane_event();
  test_presets_and_pipeline();
  return fixture.result();
}

} // namespace hyper_lattice_tests
} // namespace hs_test
