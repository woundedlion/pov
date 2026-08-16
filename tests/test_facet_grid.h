/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <bit>
#include <limits>

#include "effects/FacetGrid.h"
#include "tests/test_shaderball.h"

namespace hs_test {
namespace facet_grid_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;
using ShaderBallWB = shaderball_tests::ShaderBallWhiteBox;

struct FacetGridWhiteBox {
  using FX = FacetGrid<SMALL_W, SMALL_H>;
  using FrameState = FX::FrameState;
  using Params = FX::Params;

  static constexpr const auto &presets() { return FX::PRESETS; }
  static const Params &params(const FX &effect) { return effect.params; }
  static bool transition_active(const FX &effect) {
    return effect.transition.active;
  }
  static uint16_t transition_evaluation(const FX &effect) {
    return effect.transition.evaluation;
  }
  static void set_transition_evaluation(FX &effect, uint16_t evaluation) {
    effect.transition.evaluation = evaluation;
  }
  static bool begin_automatic_transition(FX &effect) {
    return effect.advancePreset();
  }
  static void prepare_transition_value(FX &effect) {
    effect.prepare_transition_value();
  }
  static void finish_transition_evaluation(FX &effect) {
    effect.finish_transition_evaluation();
  }

  template <typename ReferenceFrame>
  static FrameState from_reference(const ReferenceFrame &frame) {
    Params params;
    params.source = {
        frame.params.source.pattern_freq,   frame.params.source.speed,
        frame.params.source.complexity,     frame.params.source.pattern_mix,
        frame.params.source.secondary_rate, frame.params.source.angle_rate};
    params.projection = {
        frame.params.projection.pole_fade, frame.params.projection.spin_rate,
        frame.params.projection.wander, frame.params.outer_camera.wander};
    params.mirror = {
        frame.params.warp.inner.speed,    frame.params.warp.inner.rotation,
        frame.params.warp.inner.cell_x,   frame.params.warp.inner.cell_y,
        frame.params.warp.inner.offset_x, frame.params.warp.inner.offset_y};
    params.color = {frame.params.color.palette_chroma,
                    frame.params.color.mapping_frequency,
                    frame.params.color.mapping_phase,
                    frame.params.color.phase_oscillation_depth,
                    frame.params.color.phase_oscillation_speed,
                    frame.params.color.value_opacity_low,
                    frame.params.color.value_opacity_high,
                    frame.params.color.hue_shift_amount,
                    frame.params.color.hue_noise_scale,
                    frame.params.color.hue_noise_speed};
    typename FX::PreparedMirror mirror{};
    mirror.transform.mirror.offset_x =
        frame.prepared_warp.inner.transform.mirror.offset_x;
    mirror.transform.mirror.offset_y =
        frame.prepared_warp.inner.transform.mirror.offset_y;
    mirror.rotation_cos = frame.prepared_warp.inner.rotation_cos;
    mirror.rotation_sin = frame.prepared_warp.inner.rotation_sin;
    return {frame.transforms.projection_conj,
            frame.transforms.outer_conj,
            {frame.prepared_source.primary, frame.prepared_source.secondary,
             frame.prepared_source.angle_cos, frame.prepared_source.angle_sin},
            mirror,
            frame.resources.generated_palette,
            frame.prepared_hue_rotation.lut,
            frame.prepared_hue_noise.lut,
            params,
            frame.clocks.palette_oscillation_phase};
  }

  static Color4 shade(const Vector &view, const FrameState &frame) {
    return FX::RenderPipeline::shade(view, frame);
  }
};

inline void expect_color_exact(const Color4 &actual, const Color4 &expected) {
  HS_EXPECT_EQ(actual.color.r, expected.color.r);
  HS_EXPECT_EQ(actual.color.g, expected.color.g);
  HS_EXPECT_EQ(actual.color.b, expected.color.b);
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual.alpha),
               std::bit_cast<uint32_t>(expected.alpha));
}

void test_facet_grid_identity_and_presets() {
  using WB = FacetGridWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "facet-grid");
  HS_EXPECT_EQ(WB::presets().size(), size_t{3});
  HS_EXPECT_EQ(sizeof(WB::Params), 26 * sizeof(float));
  HS_EXPECT_TRUE(sizeof(WB::FrameState) < sizeof(ShaderBallWB::FrameState));
  HS_EXPECT_TRUE(FX::PRESET_IDS[0] == "coupled-grid");
  HS_EXPECT_TRUE(FX::PRESET_IDS[1] == "direct-grid");
  HS_EXPECT_TRUE(FX::PRESET_IDS[2] == "double-map");
  HS_EXPECT_NEAR(WB::presets()[0].params.source.complexity, 0.513f, 0.0f);
  HS_EXPECT_NEAR(WB::presets()[1].params.source.complexity, 3.0f, 0.0f);
  HS_EXPECT_NEAR(WB::presets()[2].params.color.mapping_frequency, 2.0f, 0.0f);
  HS_EXPECT_EQ(FX::TRANSITION_DURATION, uint16_t{480});

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Pattern Freq",
      "Speed",
      "Source Angle Speed",
      "Complexity",
      "Pattern Mix",
      "Drift",
      "Pole Fade",
      "Projection Spin Speed",
      "Projection Wander",
      "Camera Wander",
      "Planar Warp 2 Speed",
      "Planar Warp 2 Rotation",
      "Planar Warp 2 Cell X",
      "Planar Warp 2 Cell Y",
      "Planar Warp 2 Offset X",
      "Planar Warp 2 Offset Y",
      "Palette Chroma",
      "Mapping Frequency",
      "Mapping Phase",
      "Phase Oscillation Depth",
      "Phase Oscillation Speed",
      "Value Opacity Low",
      "Value Opacity High",
      "Hue Shift Amount",
      "Hue Noise Scale",
      "Hue Noise Speed",
  };
  HS_EXPECT_EQ(effect.getParameters().size(), std::size(CONTROL_NAMES));
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

void test_facet_grid_transition_contract() {
  using WB = FacetGridWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{3});
  HS_EXPECT_TRUE(WB::begin_automatic_transition(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 WB::presets()[0].params.source.complexity, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});

  effect.setAnimationsPaused(true);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});
  effect.setAnimationsPaused(false);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION / 2);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(
      WB::params(effect).source.complexity,
      FixedPipeline::linear(WB::presets()[0].params.source.complexity,
                            WB::presets()[1].params.source.complexity, 0.5f),
      1e-6f);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 WB::presets()[1].params.source.complexity, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

void test_facet_grid_parameter_serialization() {
  using WB = FacetGridWhiteBox;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  auto snapshot = effect.serialize_parameters();
  snapshot.params.source.pattern_freq = 4.0f;
  snapshot.params.mirror.cell_y = 0.75f;
  snapshot.params.color.hue_noise_speed = 0.0005f;
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
  HS_EXPECT_NEAR(WB::params(effect).source.pattern_freq, 4.0f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).mirror.cell_y, 0.75f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).color.hue_noise_speed, 0.0005f, 0.0f);

  snapshot.schema_version += 1;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.schema_version = WB::FX::PARAMETER_SCHEMA_VERSION;
  snapshot.params.mirror.cell_x = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.params.mirror.cell_x = 9.0f;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
}

void test_facet_grid_shaderball_equivalence() {
  using WB = FacetGridWhiteBox;
  reset_effect_globals();
  ShaderBallWB::SB shaderball;
  shaderball.init();

  for (size_t preset : {size_t{11}, size_t{13}, size_t{14}}) {
    const ShaderBallWB::FrameState reference =
        ShaderBallWB::preset_frame(shaderball, preset);
    const WB::FrameState compiled = WB::from_reference(reference);
    for (int latitude_step = -9; latitude_step <= 9; ++latitude_step) {
      const float latitude = latitude_step * (0.5f * PI_F / 9.0f);
      const float radius = cosf(latitude);
      for (int longitude_step = 0; longitude_step < 37; ++longitude_step) {
        const float longitude = longitude_step * (TWO_PI_F / 37.0f);
        const Vector view(radius * cosf(longitude), sinf(latitude),
                          radius * sinf(longitude));
        expect_color_exact(WB::shade(view, compiled),
                           ShaderBallWB::stereographic_dodecahedral_grid_shade(
                               view, reference));
      }
    }
  }
}

inline int run_facet_grid_tests() {
  ModuleFixture fixture("facet_grid");
  test_facet_grid_identity_and_presets();
  test_facet_grid_transition_contract();
  test_facet_grid_parameter_serialization();
  test_facet_grid_shaderball_equivalence();
  return fixture.result();
}

} // namespace facet_grid_tests
} // namespace hs_test
