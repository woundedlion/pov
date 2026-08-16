/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <bit>
#include <limits>

#include "effects/CurlLattice.h"
#include "tests/test_shaderball.h"

namespace hs_test {
namespace curl_lattice_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;
using ShaderBallWB = shaderball_tests::ShaderBallWhiteBox;

struct CurlLatticeWhiteBox {
  using FX = CurlLattice<SMALL_W, SMALL_H>;
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
    params.lattice = {frame.params.source.lattice_cell_scale,
                      frame.params.source.lattice_shape_blend,
                      frame.params.source.lattice_softness,
                      frame.params.source.lattice_radius};
    params.projection = {frame.params.projection.pole_fade,
                         frame.params.projection.central_meridian,
                         frame.params.projection.spin_rate,
                         frame.params.projection.wander,
                         frame.params.outer_camera.wander};
    params.surface_noise = {frame.params.surface_noise.scale,
                            frame.params.surface_noise.strength,
                            frame.params.surface_noise.rate};
    params.color = {frame.params.color.palette_chroma,
                    frame.params.color.mapping_frequency,
                    frame.params.color.mapping_phase,
                    frame.params.color.phase_oscillation_depth,
                    frame.params.color.phase_oscillation_speed,
                    frame.params.color.brightness_depth,
                    frame.params.color.value_opacity_low,
                    frame.params.color.value_opacity_high,
                    frame.params.color.hue_shift_amount,
                    frame.params.color.hue_noise_scale,
                    frame.params.color.hue_noise_speed};
    return {frame.transforms.projection_conj,
            frame.transforms.outer_conj,
            frame.prepared_surface_noise.loop_offset,
            frame.resources.surface_noise,
            frame.resources.generated_palette,
            frame.prepared_hue_rotation.lut,
            frame.prepared_hue_noise.lut,
            params,
            Pullback::Color::PaletteMappingWeights::single(
                static_cast<Pullback::Color::PaletteMapping>(
                    frame.slots.palette_mapping)),
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

void test_curl_lattice_identity_and_presets() {
  using WB = CurlLatticeWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "curl-lattice");
  HS_EXPECT_EQ(WB::presets().size(), size_t{2});
  HS_EXPECT_EQ(sizeof(WB::Params), 24 * sizeof(float));
  HS_EXPECT_TRUE(sizeof(WB::FrameState) < sizeof(ShaderBallWB::FrameState));
  HS_EXPECT_TRUE(FX::PRESET_IDS[0] == "open-curl");
  HS_EXPECT_TRUE(FX::PRESET_IDS[1] == "dense-curl");
  HS_EXPECT_TRUE(WB::presets()[0].id == "open-curl");
  HS_EXPECT_TRUE(WB::presets()[1].id == "dense-curl");
  HS_EXPECT_NEAR(WB::presets()[0].params.surface_noise.scale, 1.78815627f,
                 0.0f);
  HS_EXPECT_NEAR(WB::presets()[1].params.surface_noise.scale, 3.29720306f,
                 0.0f);
  HS_EXPECT_EQ(FX::TRANSITION_DURATION, uint16_t{480});

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Lattice Cell Scale",
      "Lattice Shape",
      "Lattice Softness",
      "Lattice Radius",
      "Pole Fade",
      "Central Meridian",
      "Projection Spin Speed",
      "Projection Wander",
      "Camera Wander",
      "Surface Noise Scale",
      "Surface Noise Strength",
      "Surface Noise Speed",
      "Palette Chroma",
      "Palette Mapping",
      "Mapping Frequency",
      "Mapping Phase",
      "Phase Oscillation Depth",
      "Phase Oscillation Speed",
      "Brightness Depth",
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

void test_curl_lattice_transition_contract() {
  using WB = CurlLatticeWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{2});
  HS_EXPECT_TRUE(WB::begin_automatic_transition(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).surface_noise.scale,
                 WB::presets()[0].params.surface_noise.scale, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});

  effect.setAnimationsPaused(true);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});
  effect.setAnimationsPaused(false);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION / 2);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(
      WB::params(effect).surface_noise.scale,
      FixedPipeline::linear(WB::presets()[0].params.surface_noise.scale,
                            WB::presets()[1].params.surface_noise.scale, 0.5f),
      1e-6f);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).surface_noise.scale,
                 WB::presets()[1].params.surface_noise.scale, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

void test_curl_lattice_parameter_serialization() {
  using WB = CurlLatticeWhiteBox;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  auto snapshot = effect.serialize_parameters();
  snapshot.params.surface_noise.scale = 2.5f;
  snapshot.params.lattice.lattice_radius = 0.4f;
  snapshot.params.color.hue_noise_speed = 0.0005f;
  snapshot.params.color.palette_mapping = Pullback::Color::PaletteMapping::BELL;
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
  HS_EXPECT_NEAR(WB::params(effect).surface_noise.scale, 2.5f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).lattice.lattice_radius, 0.4f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).color.hue_noise_speed, 0.0005f, 0.0f);
  HS_EXPECT_EQ(WB::params(effect).color.palette_mapping,
               Pullback::Color::PaletteMapping::BELL);

  snapshot.schema_version += 1;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.schema_version = WB::FX::PARAMETER_SCHEMA_VERSION;
  snapshot.params.surface_noise.scale = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.params.surface_noise.scale = 65.0f;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
}

void test_curl_lattice_shaderball_equivalence() {
  using WB = CurlLatticeWhiteBox;
  reset_effect_globals();
  ShaderBallWB::SB shaderball;
  shaderball.init();

  for (size_t preset : {size_t{7}, size_t{8}}) {
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
        expect_color_exact(
            WB::shade(view, compiled),
            ShaderBallWB::sinusoidal_curl_shade(view, reference));
      }
    }
  }
}

inline int run_curl_lattice_tests() {
  ModuleFixture fixture("curl_lattice");
  test_curl_lattice_identity_and_presets();
  test_curl_lattice_transition_contract();
  test_curl_lattice_parameter_serialization();
  test_curl_lattice_shaderball_equivalence();
  return fixture.result();
}

} // namespace curl_lattice_tests
} // namespace hs_test
