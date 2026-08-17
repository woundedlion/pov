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

  template <typename ReferenceFrame>
  static FrameState from_reference(const ReferenceFrame &frame) {
    Params params;
    params.source = {frame.params.source.lattice_cell_scale,
                     frame.params.source.lattice_shape_blend,
                     frame.params.source.lattice_softness,
                     frame.params.source.lattice_radius};
    params.projection = {
        frame.params.projection.pole_fade, frame.params.projection.spin_rate,
        frame.params.projection.wander, frame.params.outer_camera.wander,
        frame.params.projection.central_meridian};
    params.surface = {frame.params.surface_noise.scale,
                      frame.params.surface_noise.strength,
                      frame.params.surface_noise.rate};
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
            {},
            {},
            {},
            nullptr,
            nullptr,
            frame.resources.generated_palette,
            frame.prepared_hue_rotation.lut,
            frame.prepared_hue_noise.lut,
            params,
            Pullback::Color::PaletteMappingWeights::single(
                static_cast<Pullback::Color::PaletteMapping>(
                    frame.slots.palette_mapping)),
            0.0f,
            0.0f,
            0.0f,
            frame.clocks.palette_oscillation_phase,
            {frame.resources.surface_noise,
             frame.prepared_surface_noise.loop_offset}};
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

inline void test_curl_lattice_identity_and_presets() {
  using WB = CurlLatticeWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "curl-lattice");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{2});
  HS_EXPECT_EQ(sizeof(WB::Params), 27 * sizeof(float));
  HS_EXPECT_TRUE(sizeof(WB::FrameState) < sizeof(ShaderBallWB::FrameState));
  HS_EXPECT_TRUE(FX::PRESET_IDS[0] == "open-curl");
  HS_EXPECT_TRUE(FX::PRESET_IDS[1] == "dense-curl");
  HS_EXPECT_NEAR(FX::preset_params(0).surface.scale, 1.78815627f, 0.0f);
  HS_EXPECT_NEAR(FX::preset_params(1).surface.scale, 3.29720306f, 0.0f);
  HS_EXPECT_EQ(FX::TRANSITION_DURATION, uint16_t{480});
  const auto values = [](const WB::Params &params) {
    return std::array<float, 23>{
        params.source.lattice_cell_scale,
        params.source.lattice_shape_blend,
        params.source.lattice_softness,
        params.source.lattice_radius,
        params.projection.pole_fade,
        params.projection.central_meridian,
        params.projection.spin_rate,
        params.projection.wander,
        params.projection.camera_wander,
        params.surface.scale,
        params.surface.strength,
        params.surface.speed,
        params.color.palette_chroma,
        params.color.mapping_frequency,
        params.color.mapping_phase,
        params.color.phase_oscillation_depth,
        params.color.phase_oscillation_speed,
        params.color.brightness_depth,
        params.color.opacity_low,
        params.color.opacity_high,
        params.color.hue_shift_amount,
        params.color.hue_noise_scale,
        params.color.hue_noise_speed,
    };
  };
  constexpr std::array<float, 23> DENSE_EXPECTED{
      0.710265636f, 1.0f, 0.455532223f,  0.290762514f, 20.0f,         0.0f,
      0.0f,         1.0f, 1.0f,          3.29720306f,  0.0759999976f, 0.0f,
      1.0f,         1.0f, -0.165999994f, 0.0f,         0.0f,          1.0f,
      1.0f,         1.0f, 0.268000007f,  2.0f,         0.0f,
  };
  const std::array<float, 23> dense_actual = values(FX::preset_params(1));
  for (size_t index = 0; index < DENSE_EXPECTED.size(); ++index) {
    HS_EXPECT_EQ(std::bit_cast<uint32_t>(dense_actual[index]),
                 std::bit_cast<uint32_t>(DENSE_EXPECTED[index]));
  }
  HS_EXPECT_TRUE(FixedLook::valid(FX::preset_params(1)));

  // The runtime rebuilds the hue-rotation LUT on the same predicate the
  // colorizer gates its view on, so both read dead at a zero shift amount.
  FixedLook::ColorParams shift;
  shift.hue_shift_amount = 0.0f;
  HS_EXPECT_FALSE(
      (FixedLook::hue_rotation_active<FixedLook::HueMode::NOISE>(shift)));
  shift.hue_shift_amount = 0.25f;
  HS_EXPECT_TRUE(
      (FixedLook::hue_rotation_active<FixedLook::HueMode::NOISE>(shift)));
  HS_EXPECT_TRUE(FX::preset_params(0).color.hue_shift_amount != 0.0f);
  HS_EXPECT_TRUE(FX::preset_params(1).color.hue_shift_amount != 0.0f);

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Lattice Cell Scale",
      "Lattice Shape",
      "Lattice Softness",
      "Lattice Radius",
      "Pole Fade",
      "Projection Spin Speed",
      "Projection Wander",
      "Camera Wander",
      "Central Meridian",
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
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

inline void test_curl_lattice_transition_contract() {
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
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(0).surface.scale, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});

  effect.setAnimationsPaused(true);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});
  effect.setAnimationsPaused(false);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION / 2);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FixedPipeline::linear(FX::preset_params(0).surface.scale,
                                       FX::preset_params(1).surface.scale,
                                       0.5f),
                 1e-6f);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(1).surface.scale, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

inline void test_curl_lattice_parameter_serialization() {
  using WB = CurlLatticeWhiteBox;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  auto snapshot = effect.serialize_parameters();
  snapshot.params.surface.scale = 2.5f;
  snapshot.params.source.lattice_radius = 0.4f;
  snapshot.params.color.hue_noise_speed = 0.0005f;
  snapshot.params.color.palette_mapping = Pullback::Color::PaletteMapping::BELL;
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
  HS_EXPECT_NEAR(WB::params(effect).surface.scale, 2.5f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).source.lattice_radius, 0.4f, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).color.hue_noise_speed, 0.0005f, 0.0f);
  HS_EXPECT_EQ(WB::params(effect).color.palette_mapping,
               Pullback::Color::PaletteMapping::BELL);

  snapshot.schema_version += 1;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.schema_version = WB::FX::PARAMETER_SCHEMA_VERSION;
  snapshot.params.surface.scale = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.params.surface.scale = 65.0f;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
}

inline void test_curl_lattice_shaderball_equivalence() {
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
