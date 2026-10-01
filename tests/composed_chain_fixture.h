/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "tests/test_composed_effect.h"
#include "tests/test_shader_chain.h"

namespace hs_test::composed_chain_tests {

template <typename T>
T &mutable_state(Pullback::Interp::ChainProgram &program, size_t index) {
  return *static_cast<T *>(const_cast<void *>(program.state_block(index)));
}

template <typename FX> void verify_export(size_t preset) {
  using namespace composed_effect_tests;
  using namespace shader_chain_tests;
  HS_CONTEXT(FX::EFFECT_ID.data());
  HS_CONTEXT("preset", preset);
  std::string name(FX::EFFECT_ID);
  std::replace(name.begin(), name.end(), '-', '_');
  const std::string text = read_document(
      std::string(HS_PROMOTED_PATTERNS_DIR "/") + name + ".shader.json");
  JsonParser parser{text};
  const auto document = parser.parse_value();
  HS_EXPECT_FALSE(parser.failed);
  if (parser.failed)
    return;
  const auto &entries = document.find("descriptor")->find("chain")->items;
  const auto &presets = document.find("preset_bank")->find("presets")->items;
  const JsonValue *values = nullptr;
  for (const auto &entry : presets)
    if (entry.find("preset_id")->text == FX::PRESET_IDS[preset])
      values = entry.find("values");
  HS_EXPECT_TRUE(values != nullptr);
  if (!values)
    return;
  std::vector<In::ChainEntryRequest> requests;
  for (const auto &entry : entries)
    requests.push_back({entry.find("label")->text.c_str(),
                        entry.find("operator")->text.c_str()});
  auto fixture = std::make_unique<ProgramFixture>();
  HS_EXPECT_EQ(fixture->program.compile(requests).code, In::ChainStatus::OK);
  auto &program = fixture->program;
  for (const auto &op : program.ops())
    for (uint16_t index = 0; index < op.op->schema_count; ++index) {
      const auto &field = op.op->schema[index];
      const auto key = std::string(op.instance) + "." + field.id;
      const auto *value = values->find(key);
      if (!value)
        continue;
      void *address = op.op->runtime.param_address(
          program.param_block(&op - program.ops().data()), index);
      if (field.enum_count) {
        uint8_t selected = 0;
        while (selected < field.enum_count &&
               value->text != field.enum_ids[selected])
          ++selected;
        HS_EXPECT_LT(selected, field.enum_count);
        *static_cast<uint8_t *>(address) = selected;
      } else {
        *static_cast<float *>(address) = static_cast<float>(value->number);
      }
    }

  reset_globals();
  FX effect;
  effect.init();
  ComposedFrameWhiteBox::set_params(
      effect, effects_tests::preset_params_or_initial<FX>(preset));
  for (int frame = 0; frame < 3; ++frame)
    ComposedFrameWhiteBox::advance(effect);
  const auto own = ComposedFrameWhiteBox::frame(effect);
  for (size_t index = 0; index < program.ops().size(); ++index) {
    const auto &op = program.ops()[index];
    const std::string_view id = op.op->operator_id;
    if (id == "sphere.rotate.v2" || id.starts_with("project.")) {
      auto &state = mutable_state<In::Op::SpatialWalkState>(program, index);
      state.wander = (id == "sphere.rotate.v2" ? own.outer_conjugate
                                               : own.projection_conjugate)
                         .conjugate();
      state.spin_phase = 0.0f;
    } else if (id == "sphere.displace.curl.v2") {
      if constexpr (!std::is_void_v<typename FX::Params::surface_type>) {
        auto &state = mutable_state<In::Op::NoisePhaseState>(program, index);
        if constexpr (requires { own.surface_phase; }) {
          state.phase = own.surface_phase;
          state.noise = *own.surface_noise;
        } else {
          const auto &resource = own.resources.template get<"surface">();
          state.phase = resource.phase;
          state.noise = *resource.noise;
        }
      }
    } else if (id == "sample.grid.v2") {
      auto &state = mutable_state<In::Op::SourceClockState>(program, index);
      if constexpr (requires { own.source_primary; })
        state = {own.source_primary, own.source_secondary, own.source_angle};
      else {
        const auto &resource = own.resources.template get<"source">();
        state = {resource.primary, resource.secondary, resource.angle};
      }
    } else if (id == "warp.mirror-tile.v2") {
      if constexpr (!std::is_void_v<typename FX::Params::inner_warp_type>) {
        auto &state = mutable_state<In::Op::WarpPhaseState>(program, index);
        if constexpr (requires { own.inner_phase; })
          state.phase = own.inner_phase;
        else
          state.phase = own.resources.template get<"inner_warp">().phase;
      }
    } else if (id == "colorize.generated-palette.v3") {
      mutable_state<In::Op::ColorClockState>(program, index).oscillation_phase =
          own.palette_oscillation_phase;
    }
  }
  In::FrameContext context;
  context.palettes.fill(own.palette);
  context.hue_rotation_lut = own.hue_rotation_lut;
  context.hue_noise_lut = own.hue_noise_lut;
  program.prepare(context);
  const auto prepared = FX::RenderPipeline::prepare(own);
  size_t visible = 0;
  for (const auto &view : sweep_views()) {
    const auto expected = FX::shade(view, prepared);
    const auto actual = program.evaluate(view, context);
    expect_color_within(actual, expected, 1, 1e-6f);
    visible += expected.alpha > 0.0f &&
               (expected.color.r || expected.color.g || expected.color.b);
  }
  HS_EXPECT_GT(visible, 0u);
}

} // namespace hs_test::composed_chain_tests
