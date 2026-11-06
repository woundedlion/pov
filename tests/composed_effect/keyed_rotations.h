/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file keyed_rotations.h
 * @brief Repeated authored rotations and column-time parity. */

#include "core/render/pullback/composed_effect.h"
#include "core/render/pullback/operators/table.h"
#include "tests/composed_frame_fixture.h"

namespace hs_test::keyed_rotation_tests {

/** @brief Two independent camera rotations over a planar source. */
struct RotationProbeSpec : Pullback::Spec {
  static constexpr bool ANIMATED_PROJECTION = false;
  static constexpr bool TEMPORAL_SAMPLING = true;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::RotationProvider<B, "first">>,
      Pullback::Stage::Rotate<Pullback::RotationProvider<B, "second">>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::STEREOGRAPHIC, B>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
          Pullback::SpiralSourceParams, B>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>>;
};

/** @brief Rotation controls that change through actual preset selection. */
struct RotationProbe
    : Pullback::ComposedEffect<96, 20, RotationProbe, RotationProbeSpec,
                               Control::HookOwnership::INHERITED> {
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"slow", "fast"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr std::array<PresetEntry<Params>, 2> PRESETS = [] {
    Params slow{};
    slow.template get<"first">() = {.wander = 0.2f, .spin_speed = 0.01f};
    slow.template get<"second">() = {.wander = 0.7f, .spin_speed = 0.03f};
    Params fast = slow;
    fast.template get<"first">() = {.wander = 0.8f, .spin_speed = 0.045f};
    fast.template get<"second">() = {.wander = 0.1f, .spin_speed = 0.0f};
    return std::array{PresetEntry<Params>{slow, Segue::Preset::Snap{}},
                      PresetEntry<Params>{fast, Segue::Preset::Snap{}}};
  }();
};

/** @brief Repeated keyed rotations match interpreter walks across preset changes. */
inline void test_keyed_rotation_runtime_parity() {
  using namespace Pullback;
  using Model = Interp::Op::Rotate;
  using Binding = RotationProbe::Binding;
  using First = Stage::Rotate<RotationProvider<Binding, "first">>;
  using Second = Stage::Rotate<RotationProvider<Binding, "second">>;
  reset_globals();
  RotationProbe effect;
  effect.init();
  HS_EXPECT_TRUE(effect.selectPreset(0));
  HS_EXPECT_TRUE(effect.getParameters().find("first.Spin Speed") != nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("second.Spin Speed") != nullptr);
  static_assert(!RotationProbe::Metadata::OUTER_CAMERA);
  HS_EXPECT_TRUE(effect.getParameters().find("Camera Wander") == nullptr);
  Model::State first{}, second{};
  Model::init(first,
              {"first", Model::ID, Interp::instance_hash("first", Model::ID)});
  Model::init(second, {"second", Model::ID,
                       Interp::instance_hash("second", Model::ID)});
  size_t preset = 0;
  for (int step = 0; step < 48; ++step) {
    const auto before = ComposedFrameWhiteBox::frame(effect);
    const auto previous_first = first;
    const auto previous_second = second;
    if (step == 12 || step == 29) {
      preset = step == 12 ? 1 : 0;
      HS_EXPECT_TRUE(effect.selectPreset(preset));
    }
    const auto params = RotationProbe::test_preset_params(preset);
    const auto &a = params.template get<"first">();
    const auto &b = params.template get<"second">();
    Model::Params first_params{.wander = a.wander, .spin_speed = a.spin_speed};
    Model::Params second_params{.wander = b.wander, .spin_speed = b.spin_speed};
    Model::advance(first, first_params);
    Model::advance(second, second_params);
    ComposedFrameWhiteBox::advance(effect);
    const auto after = ComposedFrameWhiteBox::frame(effect);
    HS_EXPECT_NEAR(after.resources.template get<"first">().spin_phase,
                   first.spin_phase, 0.0f);
    HS_EXPECT_NEAR(after.resources.template get<"second">().spin_phase,
                   second.spin_phase, 0.0f);
    const auto first_interval =
        Interp::Op::temporal_interval(previous_first, first);
    const auto second_interval =
        Interp::Op::temporal_interval(previous_second, second);
    FrameInterval interval(before, after);
    for (const float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
      const auto column = interval.at(fraction);
      auto first_state = first;
      auto second_state = second;
      Interp::Op::temporal_sample(first_state, first_interval, first_params,
                                  fraction);
      Interp::Op::temporal_sample(second_state, second_interval, second_params,
                                  fraction);
      const SphereSample input{math::Vector(0.2f, 0.7f, -0.3f).normalized(),
                               0.17f};
      const auto first_expected =
          Model::run(input, {}, first_params,
                     Model::prepare({}, first_params, first_state));
      const auto expected =
          Model::run(first_expected, {}, second_params,
                     Model::prepare({}, second_params, second_state));
      const auto actual = Second::template run<Binding>(
          First::template run<Binding>(
              input, column, First::template prepare<Binding>(column)),
          column, Second::template prepare<Binding>(column));
      HS_EXPECT_NEAR(actual.dir.x, expected.dir.x, 1e-6f);
      HS_EXPECT_NEAR(actual.dir.y, expected.dir.y, 1e-6f);
      HS_EXPECT_NEAR(actual.dir.z, expected.dir.z, 1e-6f);
      HS_EXPECT_NEAR(actual.path_length, 0.17f, 0.0f);
    }
  }
  const auto final = ComposedFrameWhiteBox::frame(effect);
  HS_EXPECT_TRUE(final.resources.template get<"first">().spin_phase !=
                 final.resources.template get<"second">().spin_phase);
}

} // namespace hs_test::keyed_rotation_tests
