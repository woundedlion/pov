/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <limits>

#include "core/math/interpolate.h"
#include "effects/KaleidoscopeSmooth.h"
#include "tests/pixel_test_util.h"
#include "tests/test_shader_workbench.h"

namespace hs_test {
namespace kaleidoscope_smooth_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;
using ShaderWorkbenchWB = shader_workbench_tests::ShaderWorkbenchWhiteBox;

struct KaleidoscopeSmoothWhiteBox {
  using FX = KaleidoscopeSmooth<SMALL_W, SMALL_H>;
  using Frame = FX::Frame;
  using Params = FX::Params;

  static constexpr size_t PARAM_CAPACITY = FX::PARAM_CAPACITY;

  static const Params &params(const FX &effect) { return effect.params; }
  static bool transition_active(const FX &effect) {
    return effect.transition.active;
  }
  static bool advance_preset(FX &effect) { return effect.advancePreset(); }
  static void drive_transition(FX &effect, float progress) {
    effect.run_blend(progress);
  }

  using Ctx = Pullback::FrameState<Params>;

  template <typename ReferenceFrame>
  static Ctx from_reference(const ReferenceFrame &frame) {
    Params params;
    params.source = {
        frame.params.source.pattern_freq,   frame.params.source.speed,
        frame.params.source.complexity,     frame.params.source.pattern_mix,
        frame.params.source.secondary_rate, frame.params.source.angle_rate};
    params.projection = {frame.params.projection.singularity_fade,
                         frame.params.projection.spin_rate,
                         frame.params.projection.wander,
                         frame.params.outer_camera.wander,
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
                    frame.params.color.brightness_bottom,
                    frame.params.color.brightness_top,
                    frame.params.color.opacity_low,
                    frame.params.color.opacity_high,
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

inline void test_kaleidoscope_smooth_identity_and_presets() {
  using WB = KaleidoscopeSmoothWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "kaleidoscope-smooth");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{4});
  HS_EXPECT_EQ(sizeof(WB::Params), 32 * sizeof(float));
  HS_EXPECT_TRUE(sizeof(WB::Frame) < sizeof(ShaderWorkbenchWB::FrameState));

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
      "Singularity Fade",
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
      "Opacity at Value 0",
      "Opacity at Value 1",
      "Hue Shift Amount",
      "Hue Noise Scale",
      "Hue Noise Speed",
  };
  HS_EXPECT_EQ(effect.getParameters().size(), std::size(CONTROL_NAMES));
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

inline void test_kaleidoscope_smooth_transition_contract() {
  using WB = KaleidoscopeSmoothWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{4});
  HS_EXPECT_TRUE(WB::advance_preset(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::drive_transition(effect, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset_params(0).source.complexity, 0.0f);

  WB::drive_transition(effect, 0.25f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 interp::linear(FX::preset_params(0).source.complexity,
                                FX::preset_params(1).source.complexity, 0.25f),
                 1e-6f);

  WB::drive_transition(effect, 0.5f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 interp::linear(FX::preset_params(0).source.complexity,
                                FX::preset_params(1).source.complexity, 0.5f),
                 1e-6f);

  WB::drive_transition(effect, 1.0f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset_params(1).source.complexity, 0.0f);
  HS_EXPECT_TRUE(WB::transition_active(effect));

  for (uint16_t frame = 4; frame < FX::TRANSITION_DURATION; ++frame)
    WB::drive_transition(effect, 0.5f);
  HS_EXPECT_NEAR(WB::params(effect).source.complexity,
                 FX::preset_params(1).source.complexity, 0.0f);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

inline void test_kaleidoscope_smooth_parameter_serialization() {
  using WB = KaleidoscopeSmoothWhiteBox;
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

inline void test_kaleidoscope_smooth_shader_workbench_equivalence() {
  using WB = KaleidoscopeSmoothWhiteBox;
  reset_effect_globals();
  ShaderWorkbenchWB::SB shader_workbench;
  shader_workbench.init();

  for (size_t preset : {size_t{11}, size_t{13}, size_t{14}}) {
    const ShaderWorkbenchWB::FrameState reference =
        ShaderWorkbenchWB::preset_frame(shader_workbench, preset);
    const WB::Ctx compiled = WB::from_reference(reference);
    for (int latitude_step = -9; latitude_step <= 9; ++latitude_step) {
      const float latitude = latitude_step * (0.5f * PI_F / 9.0f);
      const float radius = cosf(latitude);
      for (int longitude_step = 0; longitude_step < 37; ++longitude_step) {
        const float longitude = longitude_step * (TWO_PI_F / 37.0f);
        const Vector view(radius * cosf(longitude), sinf(latitude),
                          radius * sinf(longitude));
        expect_color_exact(
            WB::shade(view, compiled),
            ShaderWorkbenchWB::stereographic_dodecahedral_grid_shade(
                view, reference));
      }
    }
  }
}

inline int run_kaleidoscope_smooth_tests() {
  ModuleFixture fixture("kaleidoscope_smooth");
  test_kaleidoscope_smooth_identity_and_presets();
  test_kaleidoscope_smooth_transition_contract();
  test_kaleidoscope_smooth_parameter_serialization();
  test_kaleidoscope_smooth_shader_workbench_equivalence();
  return fixture.result();
}

} // namespace kaleidoscope_smooth_tests
} // namespace hs_test
