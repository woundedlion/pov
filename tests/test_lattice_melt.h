/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <bit>
#include <limits>

#include "core/math/interpolate.h"
#include "effects/LatticeMelt.h"
#include "tests/test_shader_workbench.h"

namespace hs_test {
namespace lattice_melt_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;
using ShaderWorkbenchWB = shader_workbench_tests::ShaderWorkbenchWhiteBox;

struct LatticeMeltWhiteBox {
  using FX = LatticeMelt<SMALL_W, SMALL_H>;
  using FrameState = FX::FrameState;
  using Params = FX::Params;

  static constexpr size_t PARAM_CAPACITY = FX::PARAM_CAPACITY;

  static const Params &params(const FX &effect) { return effect.params; }
  static bool transition_active(const FX &effect) {
    return effect.transition.active;
  }
  static bool begin_automatic_transition(FX &effect) {
    return effect.advancePreset();
  }
  static void tick_choreography(FX &effect) {
    effect.begin_automatic_transition();
  }
  static void saturate_timeline(FX &effect, float &sink) {
    while (Timeline::remaining() > 0)
      effect.timeline.add(0,
                          Animation::Transition(sink, 1.0f, 10, ease_linear));
  }
  static void clear_timeline(FX &effect) { effect.timeline.clear(); }
  static void drive_transition(FX &effect, float progress) {
    effect.preset_blend.lerp(effect.preset_blend, effect.preset_blend,
                             progress);
  }

  using Ctx = Pullback::FrameState<Params>;

  template <typename ReferenceFrame>
  static Ctx from_reference(const ReferenceFrame &frame) {
    Params params;
    params.source = {frame.params.source.lattice_cell_scale,
                     frame.params.source.lattice_shape_blend,
                     frame.params.source.lattice_softness,
                     frame.params.source.lattice_radius};
    params.projection = {frame.params.projection.singularity_fade,
                         frame.params.projection.spin_rate,
                         frame.params.projection.wander,
                         frame.params.outer_camera.wander,
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
                    frame.params.color.opacity_low,
                    frame.params.color.opacity_high,
                    static_cast<Pullback::Color::PaletteMapping>(
                        frame.slots.palette_mapping)};
    return {frame.transforms.projection_conj,
            frame.transforms.outer_conj,
            nullptr,
            nullptr,
            frame.resources.surface_noise,
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
            0.0f,
            0.0f,
            0.0f,
            0.0f,
            frame.clocks.surface_noise_time,
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

inline void test_lattice_melt_identity_and_presets() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  HS_EXPECT_TRUE(FX::EFFECT_ID == "lattice-melt");
  HS_EXPECT_EQ(FX::PRESET_IDS.size(), size_t{2});
  HS_EXPECT_EQ(sizeof(WB::Params), 27 * sizeof(float));
  HS_EXPECT_TRUE(sizeof(WB::FrameState) <
                 sizeof(ShaderWorkbenchWB::FrameState));
  for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index) {
    HS_EXPECT_FALSE(FX::PRESET_IDS[index].empty());
    HS_EXPECT_TRUE(Pullback::valid(FX::preset_params(index)));
    for (size_t earlier = 0; earlier < index; ++earlier)
      HS_EXPECT_TRUE(FX::PRESET_IDS[index] != FX::PRESET_IDS[earlier]);
  }
  HS_EXPECT_GT(FX::TRANSITION_DURATION, uint16_t{0});
  HS_EXPECT_GT(FX::PRESET_DWELL_FRAMES, uint16_t{0});

  // The runtime rebuilds the hue-rotation LUT on the same predicate the
  // colorizer gates its view on, so both read dead at a zero shift amount.
  Pullback::ColorParams shift;
  shift.hue_shift_amount = 0.0f;
  HS_EXPECT_FALSE(
      (Pullback::hue_rotation_active<Pullback::HueMode::NOISE>(shift)));
  shift.hue_shift_amount = 0.25f;
  HS_EXPECT_TRUE(
      (Pullback::hue_rotation_active<Pullback::HueMode::NOISE>(shift)));
  HS_EXPECT_TRUE(FX::preset_params(0).color.hue_shift_amount != 0.0f);
  HS_EXPECT_TRUE(FX::preset_params(1).color.hue_shift_amount != 0.0f);

  reset_effect_globals();
  FX effect;
  effect.init();
  static constexpr const char *CONTROL_NAMES[] = {
      "Lattice Cell Scale",      "Lattice Shape",
      "Lattice Softness",        "Lattice Radius",
      "Singularity Fade",        "Projection Spin Speed",
      "Projection Wander",       "Camera Wander",
      "Central Meridian",        "Surface Noise Scale",
      "Surface Noise Strength",  "Surface Noise Speed",
      "Palette Chroma",          "Palette Mapping",
      "Mapping Frequency",       "Mapping Phase",
      "Phase Oscillation Depth", "Phase Oscillation Speed",
      "Brightness Depth",        "Value Opacity Low",
      "Value Opacity High",      "Hue Shift Amount",
      "Hue Noise Scale",         "Hue Noise Speed",
  };
  HS_EXPECT_EQ(effect.getParameters().size(), std::size(CONTROL_NAMES));
  HS_EXPECT_EQ(effect.getParameters().capacity(), WB::PARAM_CAPACITY);
  for (const char *name : CONTROL_NAMES)
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
}

inline void test_lattice_melt_transition_contract() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{2});
  HS_EXPECT_TRUE(WB::begin_automatic_transition(effect));
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::drive_transition(effect, 0.0f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(0).surface.scale, 0.0f);

  effect.setAnimationsPaused(true);
  WB::drive_transition(effect, 0.25f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 interp::linear(FX::preset_params(0).surface.scale,
                                FX::preset_params(1).surface.scale, 0.25f),
                 1e-6f);
  effect.setAnimationsPaused(false);

  WB::drive_transition(effect, 0.5f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 interp::linear(FX::preset_params(0).surface.scale,
                                FX::preset_params(1).surface.scale, 0.5f),
                 1e-6f);

  WB::drive_transition(effect, 1.0f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(1).surface.scale, 0.0f);
  HS_EXPECT_TRUE(WB::transition_active(effect));

  for (uint16_t frame = 4; frame < FX::TRANSITION_DURATION; ++frame)
    WB::drive_transition(effect, 0.5f);
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(1).surface.scale, 0.0f);
  HS_EXPECT_FALSE(WB::transition_active(effect));
}

inline void test_lattice_melt_full_timeline_retries_transition() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  float sink = 0.0f;
  WB::saturate_timeline(effect, sink);
  const uint32_t dropped_before = Timeline::dropped_events();

  for (uint16_t f = 0; f < FX::PRESET_DWELL_FRAMES; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_EQ(Timeline::dropped_events(), dropped_before + 1);
  HS_EXPECT_FALSE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});

  WB::clear_timeline(effect);
  WB::tick_choreography(effect);
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
}

inline void test_lattice_melt_overshoot_finishes_on_frame_count() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_TRUE(WB::begin_automatic_transition(effect));

  bool saw_overshoot = false;
  for (uint16_t frame = 1; frame < FX::TRANSITION_DURATION; ++frame) {
    const float progress =
        ease_out_elastic(static_cast<float>(frame) / FX::TRANSITION_DURATION);
    WB::drive_transition(effect, progress);
    saw_overshoot |= progress > 1.0f;
    HS_EXPECT_TRUE(WB::transition_active(effect));
  }
  HS_EXPECT_TRUE(saw_overshoot);

  WB::drive_transition(effect, ease_out_elastic(1.0f));
  HS_EXPECT_FALSE(WB::transition_active(effect));
  HS_EXPECT_NEAR(WB::params(effect).surface.scale,
                 FX::preset_params(1).surface.scale, 0.0f);
}

/**
 * @brief Pins the dwell restart a manual parameter write owes the choreography.
 */
inline void test_lattice_melt_manual_write_restarts_dwell() {
  using WB = LatticeMeltWhiteBox;
  using FX = WB::FX;
  reset_effect_globals();
  FX effect;
  effect.init();

  for (uint16_t f = 0; f < FX::PRESET_DWELL_FRAMES; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_TRUE(WB::transition_active(effect));
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

  WB::drive_transition(effect, 0.5f);
  HS_EXPECT_TRUE(effect.updateParameter("Surface Noise Scale", 1.0f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_FALSE(WB::transition_active(effect));

  effect.setAnimationsPaused(false);
  for (int f = 0; f <= FX::PRESET_SEGUE.frames; ++f)
    WB::tick_choreography(effect);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
}

inline void test_lattice_melt_parameter_serialization() {
  using WB = LatticeMeltWhiteBox;
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

inline void test_lattice_melt_shader_workbench_equivalence() {
  using WB = LatticeMeltWhiteBox;
  reset_effect_globals();
  ShaderWorkbenchWB::SB shader_workbench;
  shader_workbench.init();

  for (size_t preset : {size_t{7}, size_t{8}}) {
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
            ShaderWorkbenchWB::sinusoidal_curl_shade(view, reference));
      }
    }
  }
}

inline int run_lattice_melt_tests() {
  ModuleFixture fixture("lattice_melt");
  test_lattice_melt_identity_and_presets();
  test_lattice_melt_transition_contract();
  test_lattice_melt_full_timeline_retries_transition();
  test_lattice_melt_overshoot_finishes_on_frame_count();
  test_lattice_melt_manual_write_restarts_dwell();
  test_lattice_melt_parameter_serialization();
  test_lattice_melt_shader_workbench_equivalence();
  return fixture.result();
}

} // namespace lattice_melt_tests
} // namespace hs_test
