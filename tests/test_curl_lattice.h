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
    return {frame.transforms.projection_conj,
            frame.transforms.outer_conj,
            frame.prepared_surface_noise.loop_offset,
            frame.resources.surface_noise,
            frame.resources.generated_palette,
            frame.prepared_hue_rotation.lut,
            frame.prepared_hue_noise.lut,
            {},
            frame.params.surface_noise.scale};
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
  HS_EXPECT_EQ(sizeof(WB::Params), size_t{4});
  HS_EXPECT_TRUE(sizeof(WB::FrameState) < sizeof(ShaderBallWB::FrameState));
  HS_EXPECT_TRUE(FX::PRESET_IDS[0] == "open-curl");
  HS_EXPECT_TRUE(FX::PRESET_IDS[1] == "dense-curl");
  HS_EXPECT_TRUE(WB::presets()[0].id == "open-curl");
  HS_EXPECT_TRUE(WB::presets()[1].id == "dense-curl");
  HS_EXPECT_NEAR(WB::presets()[0].params.surface_noise_scale, 1.78815627f,
                 0.0f);
  HS_EXPECT_NEAR(WB::presets()[1].params.surface_noise_scale, 3.29720306f,
                 0.0f);
  HS_EXPECT_EQ(FX::TRANSITION_DURATION, uint16_t{480});
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
  HS_EXPECT_NEAR(WB::params(effect).surface_noise_scale,
                 WB::presets()[0].params.surface_noise_scale, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});

  effect.setAnimationsPaused(true);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_EQ(WB::transition_evaluation(effect), uint16_t{1});
  effect.setAnimationsPaused(false);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION / 2);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(
      WB::params(effect).surface_noise_scale,
      FixedPipeline::linear(WB::presets()[0].params.surface_noise_scale,
                            WB::presets()[1].params.surface_noise_scale, 0.5f),
      1e-6f);

  WB::set_transition_evaluation(effect, FX::TRANSITION_DURATION);
  WB::prepare_transition_value(effect);
  HS_EXPECT_NEAR(WB::params(effect).surface_noise_scale,
                 WB::presets()[1].params.surface_noise_scale, 0.0f);
  WB::finish_transition_evaluation(effect);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

void test_curl_lattice_parameter_serialization() {
  using WB = CurlLatticeWhiteBox;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  auto snapshot = effect.serialize_parameters();
  snapshot.surface_noise_scale = 2.5f;
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
  HS_EXPECT_NEAR(WB::params(effect).surface_noise_scale, 2.5f, 0.0f);

  snapshot.schema_version += 1;
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.schema_version = WB::FX::PARAMETER_SCHEMA_VERSION;
  snapshot.surface_noise_scale = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(effect.restore_parameters(snapshot));
  snapshot.surface_noise_scale = 65.0f;
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
        expect_color_exact(WB::shade(view, compiled),
                           ShaderBallWB::sinusoidal_curl_shade(view,
                                                               reference));
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
