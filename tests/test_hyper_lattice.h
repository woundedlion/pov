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

struct HyperLatticeWhiteBox {
  using Effect = HyperLattice<96, 20>;

  static HL::Vec4 origin(const Effect &effect) { return effect.origin; }
  static std::array<float, 6> rotation_phase(const Effect &effect) {
    return effect.rotation_phase;
  }
  static HL::Params &params(Effect &effect) { return effect.params; }
  static void step_depth_palette(Effect &effect) {
    effect.depth_palette.step();
  }
  static bool depth_palette_fading(const Effect &effect) {
    return effect.depth_palette.fading();
  }
  static Pixel depth_color(const Effect &effect, float amount) {
    return effect.depth_palette.palette().get(amount).color;
  }
  static Pixel axis_color(const Effect &effect, float amount) {
    return effect.axis_palette.get(amount).color;
  }
  static const BakedPalette *depth_palette(const Effect &effect) {
    return &effect.depth_palette.palette();
  }
  static const BakedPalette *axis_palette(const Effect &effect) {
    return &effect.axis_palette;
  }
};

inline void test_periodic_distance() {
  HS_EXPECT_EQ(HL::periodic_distance(0.0f), 0.0f);
  HS_EXPECT_EQ(HL::periodic_distance(1.0f), 0.0f);
  HS_EXPECT_NEAR(HL::periodic_distance(-1.25f), 0.25f, 1e-6f);
  HS_EXPECT_EQ(HL::periodic_distance(-0.5f), 0.5f);
  HS_EXPECT_EQ(HL::periodic_distance(1.5f), 0.5f);
  HS_EXPECT_EQ(HL::periodic_distance(2.5f), 0.5f);
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
  frame.params.mode = HL::LatticeMode::FOUR_D_SLICE;
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

  frame.params.mode = HL::LatticeMode::THREE_D;
  const HL::Vec4 cubic = HL::prepare_trace(frame).world_to_lattice.apply(
      {{1.0f, 0.0f, 0.0f, 0.0f}});
  HS_EXPECT_EQ(cubic[3], 0.0f);
}

inline void test_reflection_convention() {
  const HL::Vec4 center =
      HL::reflected_direction(X_AXIS, HL::ReflectionMode::CHROME, 1.0f);
  HS_EXPECT_NEAR(center[0], 1.0f, 2e-4f);
  HS_EXPECT_EQ(center[1], 0.0f);
  HS_EXPECT_EQ(center[2], 0.0f);

  const HL::Vec4 rim =
      HL::reflected_direction(Y_AXIS, HL::ReflectionMode::CHROME, 1.0f);
  HS_EXPECT_NEAR(rim[0], -0.5547f, 2e-3f);
  HS_EXPECT_NEAR(rim[1], 0.83205f, 2e-3f);
  const HL::Vec4 radial =
      HL::reflected_direction(Y_AXIS, HL::ReflectionMode::RADIAL, 1.0f);
  HS_EXPECT_EQ(radial[1], 1.0f);

  const HL::Vec4 open =
      HL::reflected_direction(Y_AXIS, HL::ReflectionMode::CHROME, 0.0f);
  HS_EXPECT_NEAR(open[0], 0.0f, 2e-4f);
  HS_EXPECT_NEAR(open[1], 1.0f, 2e-4f);
  const HL::Vec4 back = HL::reflected_direction(
      Vector(-1.0f, 0.0f, 0.0f), HL::ReflectionMode::CHROME, 1.0f);
  HS_EXPECT_NEAR(back[0], -1.0f, 2e-4f);
}

inline void test_resolution_aware_wire_coverage() {
  constexpr float LOW_RES = HL::pixel_half_angle<96, 20>();
  constexpr float HIGH_RES = HL::pixel_half_angle<288, 144>();
  static_assert(LOW_RES > HIGH_RES);

  constexpr float RADIUS = 0.1f;
  constexpr float HALF_WIDTH = 0.02f;
  constexpr float OFFSET = 0.01f;
  const float inside = HL::wire_coverage((RADIUS - OFFSET) * (RADIUS - OFFSET),
                                         RADIUS, HALF_WIDTH);
  const float boundary = HL::wire_coverage(RADIUS * RADIUS, RADIUS, HALF_WIDTH);
  const float outside = HL::wire_coverage((RADIUS + OFFSET) * (RADIUS + OFFSET),
                                          RADIUS, HALF_WIDTH);
  HS_EXPECT_NEAR(boundary, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(inside + outside, 1.0f, 1e-5f);
}

inline void test_near_field_fade() {
  constexpr float RADIUS = 0.1f;
  constexpr float NEAR_START = 1.5f * RADIUS;
  constexpr float NEAR_INV_SPAN = 1.0f / (2.5f * RADIUS);
  HS_EXPECT_EQ(HL::near_field_coverage(0.15f, NEAR_START, NEAR_INV_SPAN), 0.0f);
  HS_EXPECT_GT(HL::near_field_coverage(0.275f, NEAR_START, NEAR_INV_SPAN),
               0.0f);
  HS_EXPECT_LT(HL::near_field_coverage(0.275f, NEAR_START, NEAR_INV_SPAN),
               1.0f);
  HS_EXPECT_EQ(HL::near_field_coverage(0.4f, NEAR_START, NEAR_INV_SPAN), 1.0f);
}

inline void test_far_shell_fade() {
  HS_EXPECT_EQ(HL::shell_horizon_coverage(0, 2, 0.75f, 1.0f), 1.0f);
  HS_EXPECT_EQ(HL::shell_horizon_coverage(1, 2, 1.0f, 1.0f), 1.0f);
  HS_EXPECT_NEAR(HL::shell_horizon_coverage(1, 2, 1.5f, 1.0f), 0.5f, 1e-6f);
  HS_EXPECT_EQ(HL::shell_horizon_coverage(1, 2, 2.0f, 1.0f), 0.0f);
}

inline void test_pause_does_not_stop_motion() {
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);

  const HL::Vec4 origin_before = HyperLatticeWhiteBox::origin(effect);
  const auto rotation_before = HyperLatticeWhiteBox::rotation_phase(effect);
  effect.draw_frame();
  effect.advance_display();
  const HL::Vec4 origin_after = HyperLatticeWhiteBox::origin(effect);
  const auto rotation_after = HyperLatticeWhiteBox::rotation_phase(effect);
  HS_EXPECT_NE(origin_after[0], origin_before[0]);
  HS_EXPECT_NE(rotation_after[0], rotation_before[0]);

  HyperLatticeWhiteBox::params(effect).speed = 0.0f;
  const HL::Vec4 stopped_before = HyperLatticeWhiteBox::origin(effect);
  effect.draw_frame();
  effect.advance_display();
  const HL::Vec4 stopped_after = HyperLatticeWhiteBox::origin(effect);
  for (int axis = 0; axis < HL::DIMENSIONS; ++axis)
    HS_EXPECT_EQ(stopped_after[axis], stopped_before[axis]);
}

inline void test_depth_palette_mutates_slowly_while_paused() {
  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  effect.setAnimationsPaused(true);
  const Pixel initial = HyperLatticeWhiteBox::depth_color(effect, 0.5f);
  const Pixel axis = HyperLatticeWhiteBox::axis_color(effect, 0.5f);
  HyperLatticeWhiteBox::step_depth_palette(effect);
  HS_EXPECT_TRUE(HyperLatticeWhiteBox::depth_palette_fading(effect));
  HS_EXPECT_EQ(HyperLatticeWhiteBox::depth_color(effect, 0.5f), initial);

  for (int frame = 0; frame < 480; ++frame)
    HyperLatticeWhiteBox::step_depth_palette(effect);
  HS_EXPECT_NE(HyperLatticeWhiteBox::depth_color(effect, 0.5f), initial);
  HS_EXPECT_EQ(HyperLatticeWhiteBox::axis_color(effect, 0.5f), axis);
}

inline void test_next_plane_is_strict() {
  HS_EXPECT_EQ(HL::next_plane_offset(0.0f, true), 1.0f);
  HS_EXPECT_EQ(HL::next_plane_offset(0.0f, false), 1.0f);
  HS_EXPECT_EQ(HL::next_plane_offset(0.25f, true), 0.75f);
  HS_EXPECT_EQ(HL::next_plane_offset(0.25f, false), 0.25f);
}

inline void test_trace_layers_are_front_to_back() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.origin = {{0.25f, 0.0f, 0.31f, 0.43f}};
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  float previous_distance = 0.0f;
  int layers = 0;
  HL::trace_layers(X_AXIS, prepared, [&](const HL::TraceHit &hit) {
    HS_EXPECT_GT(hit.coverage, 0.0f);
    HS_EXPECT_LE(hit.coverage, 1.0f);
    HS_EXPECT_GT(hit.distance, previous_distance);
    previous_distance = hit.distance;
    ++layers;
    return true;
  });
  HS_EXPECT_EQ(layers, 2);
}

inline void test_layer_composite_reveals_background() {
  HL::LayerComposite composite;
  composite.add(Pixel(65535, 0, 0), 0.5f);
  composite.add(Pixel(0, 65535, 0), 1.0f);
  const Color4 result = composite.finish();
  HS_EXPECT_EQ(result.alpha, 1.0f);
  HS_EXPECT_NEAR(result.color.r, 32768, 1);
  HS_EXPECT_NEAR(result.color.g, 32768, 1);
  HS_EXPECT_EQ(result.color.b, 0);
}

inline void test_surface_origin_parallax() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.params.reflection = HL::ReflectionMode::RADIAL;
  frame.params.sphere_radius = 0.0f;
  frame.origin = {{0.25f, 0.0f, 0.31f, 0.43f}};
  const HL::TraceHit centered = HL::trace(X_AXIS, HL::prepare_trace(frame));
  frame.params.sphere_radius = 0.4f;
  const HL::TraceHit surfaced = HL::trace(X_AXIS, HL::prepare_trace(frame));
  HS_EXPECT_NEAR(centered.distance, 0.75f, 1e-6f);
  HS_EXPECT_NEAR(surfaced.distance, 0.35f, 1e-6f);
}

inline void test_hyperplane_event() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(1);
  frame.params.reflection = HL::ReflectionMode::RADIAL;
  frame.params.sphere_radius = 0.4f;
  frame.origin = {{0.0f, 0.0f, 0.31f, 0.25f}};
  frame.rotation_phase[3] = 0.5f * PI_F;
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  const HL::TraceHit hit = HL::trace(X_AXIS, prepared);
  HS_EXPECT_GT(hit.coverage, 0.0f);
  HS_EXPECT_EQ(hit.free_axis, uint8_t(2));
  HS_EXPECT_NEAR(hit.distance, 0.35f, 3e-4f);
}

inline void test_coincident_planes_form_one_layer() {
  HL::FrameState frame{};
  frame.params = HyperLattice<96, 20>::preset_params(0);
  frame.params.reflection = HL::ReflectionMode::RADIAL;
  frame.params.sphere_radius = 0.0f;
  frame.origin = {{0.25f, 0.25f, 0.31f, 0.43f}};
  const HL::PreparedTrace prepared = HL::prepare_trace(frame);
  constexpr float INV_SQRT_TWO = 0.707106781f;
  float previous_distance = 0.0f;
  int layers = 0;
  HL::trace_layers(Vector(INV_SQRT_TWO, INV_SQRT_TWO, 0.0f), prepared,
                   [&](const HL::TraceHit &hit) {
                     HS_EXPECT_GT(hit.distance, previous_distance);
                     previous_distance = hit.distance;
                     ++layers;
                     return true;
                   });
  HS_EXPECT_EQ(layers, 2);
}

inline void test_render_signature() {
  static constexpr Vector DIRECTIONS[] = {
      {1.0f, 0.0f, 0.0f},
      {-1.0f, 0.0f, 0.0f},
      {0.0f, 1.0f, 0.0f},
      {0.0f, 0.0f, 1.0f},
      {0.577350269f, 0.577350269f, 0.577350269f},
      {-0.707106781f, 0.707106781f, 0.0f},
      {0.301511345f, -0.904534034f, 0.301511345f},
      {0.267261242f, 0.534522484f, 0.801783726f},
      {-0.408248290f, 0.816496581f, -0.408248290f},
      {0.923879533f, 0.382683432f, 0.0f},
      {-0.270598050f, 0.653281482f, 0.707106781f},
      {0.5f, -0.5f, 0.707106781f},
  };
  static constexpr HL::Vec4 ORIGINS[] = {
      {{0.17f, 0.31f, 0.43f, 0.59f}},
      {{0.91f, 0.07f, 0.73f, 0.37f}},
  };
  static constexpr std::array<float, 6> ROTATIONS[] = {
      {0.0f, 0.3f, 0.7f, 0.0f, 0.0f, 0.0f},
      {1.1f, 2.3f, 0.4f, 0.0f, 0.0f, 0.0f},
  };

  HyperLatticeWhiteBox::Effect effect;
  effect.init();
  uint64_t signature = hs_test::FNV1A64_BASIS;
  for (size_t preset = 0; preset < std::size(ORIGINS); ++preset) {
    const HL::FrameState frame{
        HyperLatticeWhiteBox::Effect::preset_params(preset),
        ORIGINS[preset],
        ROTATIONS[preset],
        HL::pixel_half_angle<288, 144>(),
        HyperLatticeWhiteBox::depth_palette(effect),
        HyperLatticeWhiteBox::axis_palette(effect),
    };
    const HL::PreparedTrace prepared = HL::prepare_trace(frame);
    for (size_t sample = 0; sample < std::size(DIRECTIONS); ++sample) {
      const Color4 color =
          HL::shade({DIRECTIONS[sample], 0.0f}, frame, prepared);
      signature = hs_test::fnv1a64_channel(signature, color.color.r);
      signature = hs_test::fnv1a64_channel(signature, color.color.g);
      signature = hs_test::fnv1a64_channel(signature, color.color.b);
      signature = hs_test::fnv1a64_channel(signature, frac_to_q16(color.alpha));
    }
  }
  HS_EXPECT_EQ(signature, uint64_t(3173600656855287585ull));
}

inline void test_presets_and_pipeline() {
  using Effect = HyperLattice<96, 20>;
  static_assert(HL::RenderPipeline::STAGE_COUNT == 1);
  static_assert(HL::RenderPipeline::Validation::MONOTONE);
  static_assert(HL::RenderPipeline::Validation::ENTRY);
  static_assert(HL::RenderPipeline::Validation::EXIT);
  for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index)
    HS_EXPECT_TRUE(Effect::valid_params(Effect::preset_params(index)));
  static_assert(Effect::PRESET_IDS.size() == 2);
  static_assert(Effect::PRESET_IDS[0] == "cubic-flight");
  static_assert(Effect::PRESET_IDS[1] == "hypercube-flight");

  constexpr HL::Params preset1 = Effect::preset_params(0);
  static_assert(preset1.mode == HL::LatticeMode::THREE_D);
  static_assert(preset1.sphere_radius == 0.4f);
  static_assert(preset1.wire_radius == 0.055f);
  static_assert(preset1.softness == 0.08f);
  static_assert(preset1.far_cells == 4.198f);
  static_assert(preset1.aa_strength == 1.0f);
  static_assert(preset1.speed == 0.05f);
  static_assert(preset1.spin_3d == 0.015f);
  static_assert(preset1.spin_4d == 0.0f);
  static_assert(preset1.chrome_warp == 0.65f);
  static_assert(preset1.reflection == HL::ReflectionMode::RADIAL);
  static_assert(preset1.color == HL::ColorMode::DEPTH);
  static_assert(preset1.shells == HL::ShellCount::TWO);

  constexpr HL::Params preset2 = Effect::preset_params(1);
  static_assert(preset2.mode == HL::LatticeMode::FOUR_D_SLICE);
  static_assert(preset2.sphere_radius == 0.0f);
  static_assert(preset2.wire_radius == 0.03546f);
  static_assert(preset2.softness == 0.029612f);
  static_assert(preset2.far_cells == 7.264f);
  static_assert(preset2.aa_strength == 1.0f);
  static_assert(preset2.speed == 0.03f);
  static_assert(preset2.spin_3d == 0.01089f);
  static_assert(preset2.spin_4d == 0.015f);
  static_assert(preset2.chrome_warp == 0.65f);
  static_assert(preset2.reflection == HL::ReflectionMode::CHROME);
  static_assert(preset2.color == HL::ColorMode::DEPTH);
  static_assert(preset2.shells == HL::ShellCount::THREE);
}

inline void test_dimension_dropdown_and_mode_lerp() {
  using Effect = HyperLatticeWhiteBox::Effect;
  Effect effect;
  effect.init();
  const ParamDef *dimension = effect.getParameters().find("Dimension");
  HS_EXPECT_TRUE(dimension != nullptr);
  HS_EXPECT_TRUE(dimension->is_enum());
  HS_EXPECT_EQ(dimension->option_count, 3);
  HS_EXPECT_EQ(std::string_view(dimension->options[0]), std::string_view("3D"));
  HS_EXPECT_EQ(std::string_view(dimension->options[1]),
               std::string_view("Dimensional Rift"));
  HS_EXPECT_EQ(std::string_view(dimension->options[2]),
               std::string_view("4D Slice"));
  HS_EXPECT_EQ(std::string_view(dimension->export_options[2]),
               std::string_view("LatticeMode::FOUR_D_SLICE"));
  HS_EXPECT_TRUE(effect.updateParameter("Dimension", 2.0f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(HyperLatticeWhiteBox::params(effect).mode ==
                 HL::LatticeMode::FOUR_D_SLICE);

  HL::Params start;
  HL::Params target;
  target.mode = HL::LatticeMode::FOUR_D_SLICE;
  HL::Params blended;
  blended.lerp(start, target, 0.49f);
  HS_EXPECT_TRUE(blended.mode == HL::LatticeMode::THREE_D);
  blended.lerp(start, target, 0.5f);
  HS_EXPECT_TRUE(blended.mode == HL::LatticeMode::FOUR_D_SLICE);
}

inline int run_hyper_lattice_tests() {
  hs_test::ModuleFixture fixture("hyper_lattice");
  test_periodic_distance();
  test_edge_metrics();
  test_so4_rotation();
  test_reflection_convention();
  test_resolution_aware_wire_coverage();
  test_near_field_fade();
  test_far_shell_fade();
  test_pause_does_not_stop_motion();
  test_depth_palette_mutates_slowly_while_paused();
  test_next_plane_is_strict();
  test_trace_layers_are_front_to_back();
  test_layer_composite_reveals_background();
  test_surface_origin_parallax();
  test_hyperplane_event();
  test_coincident_planes_form_one_layer();
  test_render_signature();
  test_presets_and_pipeline();
  test_dimension_dropdown_and_mode_lerp();
  return fixture.result();
}

} // namespace hyper_lattice_tests
} // namespace hs_test
