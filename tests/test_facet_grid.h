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

  static constexpr size_t PARAM_CAPACITY = FX::PARAM_CAPACITY;

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

  using Ctx = FixedLook::FrameState<Params>;

  template <typename ReferenceFrame>
  static Ctx from_reference(const ReferenceFrame &frame) {
    Params params;
    params.source = {
        frame.params.source.pattern_freq,   frame.params.source.speed,
        frame.params.source.complexity,     frame.params.source.pattern_mix,
        frame.params.source.secondary_rate, frame.params.source.angle_rate};
    params.projection = {
        frame.params.projection.pole_fade, frame.params.projection.spin_rate,
        frame.params.projection.wander, frame.params.outer_camera.wander,
        frame.params.projection.central_meridian};
    params.inner_warp = {
        frame.params.warp.inner.speed,    frame.params.warp.inner.rotation,
        frame.params.warp.inner.cell_x,   frame.params.warp.inner.cell_y,
        frame.params.warp.inner.offset_x, frame.params.warp.inner.offset_y};
    params.color = {frame.params.color.hue_shift_amount,
                    frame.params.color.hue_noise_scale,
                    frame.params.color.hue_noise_speed,
                    frame.params.color.palette_chroma,
                    frame.params.color.mapping_frequency,
                    frame.params.color.mapping_phase,
                    frame.params.color.phase_oscillation_depth,
                    frame.params.color.phase_oscillation_speed,
                    frame.params.color.brightness_depth,
                    frame.params.color.value_opacity_low,
                    frame.params.color.value_opacity_high,
                    static_cast<Pullback::Color::PaletteMapping>(
                        frame.slots.palette_mapping)};
    return {frame.transforms.projection_conj,
            frame.transforms.outer_conj,
            nullptr,
            nullptr,
            nullptr,
            frame.resources.generated_palette,
            frame.prepared_hue_rotation.lut,
            frame.prepared_hue_noise.lut,
            params,
            Pullback::Color::PaletteMappingWeights::single(
                static_cast<Pullback::Color::PaletteMapping>(
                    frame.slots.palette_mapping)),
            frame.clocks.source_primary,
            frame.clocks.source_secondary,
            frame.clocks.source_angle,
            0.0f,
            frame.clocks.warp_inner_phase,
            0.0f,
            0.0f,
            0.0f,
            frame.clocks.palette_oscillation_phase};
  }

  static Color4 shade(const Vector &view, const Ctx &frame) {
    return FX::RenderPipeline::shade(view, FX::RenderPipeline::prepare(frame));
  }
};

inline void expect_color_exact(const Color4 &actual, const Color4 &expected) {
  HS_EXPECT_EQ(actual.color.r, expected.color.r);
  HS_EXPECT_EQ(actual.color.g, expected.color.g);
  HS_EXPECT_EQ(actual.color.b, expected.color.b);
  HS_EXPECT_EQ(std::bit_cast<uint32_t>(actual.alpha),
               std::bit_cast<uint32_t>(expected.alpha));
}

inline void test_facet_grid_identity_and_presets() {
  using WB = FacetGridWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "facet-grid");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{4});
  HS_EXPECT_EQ(sizeof(WB::Params), 31 * sizeof(float));
  HS_EXPECT_TRUE(sizeof(WB::FrameState) < sizeof(ShaderBallWB::FrameState));
  HS_EXPECT_TRUE(FX::PRESET_IDS[0] == "coupled-grid");
  HS_EXPECT_TRUE(FX::PRESET_IDS[1] == "direct-grid");
  HS_EXPECT_TRUE(FX::PRESET_IDS[2] == "double-map");
  HS_EXPECT_TRUE(FX::PRESET_IDS[3] == "stretched-grid");
  HS_EXPECT_NEAR(FX::preset_params(0).source.complexity, 0.513f, 0.0f);
  HS_EXPECT_NEAR(FX::preset_params(1).source.complexity, 3.0f, 0.0f);
  HS_EXPECT_NEAR(FX::preset_params(2).color.mapping_frequency, 2.0f, 0.0f);
  const auto values = [](const WB::Params &params) {
    return std::array<float, 26>{
        params.source.pattern_freq,
        params.source.speed,
        params.source.angle_rate,
        params.source.complexity,
        params.source.pattern_mix,
        params.source.secondary_rate,
        params.projection.pole_fade,
        params.projection.spin_rate,
        params.projection.wander,
        params.projection.camera_wander,
        params.inner_warp.speed,
        params.inner_warp.rotation,
        params.inner_warp.cell_x,
        params.inner_warp.cell_y,
        params.inner_warp.offset_x,
        params.inner_warp.offset_y,
        params.color.palette_chroma,
        params.color.mapping_frequency,
        params.color.mapping_phase,
        params.color.phase_oscillation_depth,
        params.color.phase_oscillation_speed,
        params.color.opacity_low,
        params.color.opacity_high,
        params.color.hue_shift_amount,
        params.color.hue_noise_scale,
        params.color.hue_noise_speed,
    };
  };
  constexpr std::array<float, 26> STRETCHED_EXPECTED{
      2.9059f,     0.0f,      0.026999999f, 3.0f, 1.0f,          0.8f,
      3.432f,      0.0f,      0.165f,       1.0f, 0.0027299998f, 3.455752f,
      0.22321875f, 5.085703f, 0.0f,         0.0f, 1.0f,          1.558f,
      0.0f,        0.0f,      0.0f,         1.0f, 1.0f,          0.366f,
      1.4721563f,  0.0f,
  };
  const std::array<float, 26> stretched_actual = values(FX::preset_params(3));
  for (size_t index = 0; index < STRETCHED_EXPECTED.size(); ++index) {
    HS_EXPECT_EQ(std::bit_cast<uint32_t>(stretched_actual[index]),
                 std::bit_cast<uint32_t>(STRETCHED_EXPECTED[index]));
  }
  HS_EXPECT_EQ(FX::TRANSITION_DURATION, uint16_t{480});
  HS_EXPECT_TRUE(FixedLook::valid(FX::preset_params(3)));

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Pattern Freq",
      "Speed",
      "Complexity",
      "Pattern Mix",
      "Drift",
      "Source Angle Speed",
      "Pole Fade",
      "Projection Spin Speed",
      "Projection Wander",
      "Camera Wander",
      "Planar Warp 2 Speed",
      "Mirror Rotation",
      "Mirror Cell X",
      "Mirror Cell Y",
      "Mirror Offset X",
      "Mirror Offset Y",
      "Palette Chroma",
      "Palette Mapping",
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
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

inline void test_facet_grid_transition_contract() {
  using WB = FacetGridWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{4});
  HS_EXPECT_TRUE(WB::begin_automatic_transition(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset_params(0).source.complexity, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});

  effect.setAnimationsPaused(true);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{2});
  effect.setAnimationsPaused(false);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION / 2);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FixedPipeline::linear(FX::preset_params(0).source.complexity,
                                       FX::preset_params(1).source.complexity,
                                       0.5f),
                 1e-6f);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset_params(1).source.complexity, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

inline void test_facet_grid_parameter_serialization() {
  using WB = FacetGridWhiteBox;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  auto snapshot = effect.serialize_parameters();
  snapshot.params.source.pattern_freq = 4.0f;
  snapshot.params.inner_warp.cell_y = 0.75f;
  snapshot.params.color.hue_noise_speed = 0.0005f;
  snapshot.params.color.palette_mapping =
      Pullback::Color::PaletteMapping::REVERSE;
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
  HS_EXPECT_NEAR(WB::params(effect).source.pattern_freq, 4.0f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).inner_warp.cell_y, 0.75f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).color.hue_noise_speed, 0.0005f, 0.0f);
  HS_EXPECT_EQ(WB::params(effect).color.palette_mapping,
               Pullback::Color::PaletteMapping::REVERSE);

  snapshot.schema_version += 1;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.schema_version = WB::FX::PARAMETER_SCHEMA_VERSION;
  snapshot.params.inner_warp.cell_x = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.params.inner_warp.cell_x = 9.0f;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
}

inline void test_facet_grid_shaderball_equivalence() {
  using WB = FacetGridWhiteBox;
  reset_effect_globals();
  ShaderBallWB::SB shaderball;
  shaderball.init();

  for (size_t preset : {size_t{11}, size_t{13}, size_t{14}}) {
    const ShaderBallWB::FrameState reference =
        ShaderBallWB::preset_frame(shaderball, preset);
    const WB::Ctx compiled = WB::from_reference(reference);
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
