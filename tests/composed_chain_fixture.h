/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Document-built chain and composed frame parity fixtures.
 */
#pragma once

#include "tests/test_composed_effect.h"
#include "tests/test_shader_chain.h"

namespace hs_test::composed_chain_tests {

template <typename T>
T &mutable_state(Pullback::Interp::ChainProgram &program, size_t index) {
  const auto &layout = program.ops()[index].op->runtime.state;
  void *state = const_cast<void *>(program.state_block(index));
  const bool MATCHES = layout.size == sizeof(T) && layout.align == alignof(T) &&
                       reinterpret_cast<uintptr_t>(state) % alignof(T) == 0;
  HS_EXPECT_TRUE(MATCHES);
  HS_CHECK(MATCHES, "composed parity fixture state type mismatch");
  return *static_cast<T *>(state);
}

template <typename FX>
void verify_export(size_t preset, int channel_tolerance,
                   float alpha_tolerance) {
  using namespace composed_effect_tests;
  using namespace shader_chain_tests;
  HS_CONTEXT(FX::EFFECT_ID.data());
  HS_CONTEXT("preset", preset);
  const auto loaded = load_promoted_document<FX>();
  if (!loaded)
    return;
  const JsonValue &document = *loaded;
  const JsonValue *chain = document_chain(document);
  const JsonValue *presets = document_presets(document);
  HS_EXPECT_TRUE(chain != nullptr);
  HS_EXPECT_TRUE(presets != nullptr);
  if (!chain || !presets)
    return;
  const JsonValue *values = preset_values(*presets, FX::PRESET_IDS[preset]);
  HS_EXPECT_TRUE(values != nullptr);
  if (!values)
    return;
  std::vector<In::ChainEntryRequest> requests;
  for (const auto &entry : chain->items) {
    const auto *label = entry.find("label");
    const auto *op = entry.find("operator");
    HS_EXPECT_TRUE(label != nullptr);
    HS_EXPECT_TRUE(op != nullptr);
    if (!label || !op)
      return;
    requests.push_back({label->text.c_str(), op->text.c_str()});
  }
  auto fixture = std::make_unique<ProgramFixture>();
  const auto compiled = fixture->program.compile(requests).code;
  HS_EXPECT_EQ(compiled, In::ChainStatus::OK);
  if (compiled != In::ChainStatus::OK)
    return;
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
        if (selected == field.enum_count)
          return;
        *static_cast<uint8_t *>(address) = selected;
      } else {
        *static_cast<float *>(address) = static_cast<float>(value->number);
      }
    }

  reset_globals();
  FX effect;
  effect.init();
  ComposedFrameWhiteBox::set_params(effect,
                                    effects_tests::preset_params<FX>(preset));
  for (int frame = 0; frame < 3; ++frame)
    ComposedFrameWhiteBox::advance(effect);
  auto own = ComposedFrameWhiteBox::frame(effect);
  if constexpr (std::is_same_v<FX, GridWorms<SMALL_W, SMALL_H>>)
    own.outer_conjugate =
        math::make_rotation(math::Vector(0.3f, 0.7f, -0.2f).normalized(), 0.8f)
            .conjugate();
  for (size_t index = 0; index < program.ops().size(); ++index) {
    const auto &op = program.ops()[index];
    const std::string_view id = op.op->operator_id;
    if (id == "sphere.rotate.v2" || id.starts_with("project.")) {
      auto &state = mutable_state<In::Op::SpatialWalkState>(program, index);
      state.wander = (id == "sphere.rotate.v2" ? own.outer_conjugate
                                               : own.projection_conjugate)
                         .conjugate();
      state.spin_phase = 0.0f;
      if (id == "sphere.rotate.v2")
        own.params.visit([&]<typename Resource>(const auto &) {
          if constexpr (Resource::KIND == Pullback::ResourceKind::ROTATION)
            if (Resource::KEY.view() == op.instance) {
              const auto &clock = own.resources.template get<Resource::KEY>();
              state.wander = clock.wander;
              state.spin_phase = clock.spin_phase;
            }
        });
    } else if (id == "sphere.lens.twist.v2") {
      own.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::TWIST)
          if (Resource::KEY.view() == op.instance) {
            const auto &phase = own.resources.template get<Resource::KEY>();
            mutable_state<In::Op::TwistOrbitState>(program, index) = {
                phase.phi, phase.lambda};
          }
      });
    } else if (id == "sphere.displace.ripple.v2") {
      own.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::SURFACE &&
                      std::is_same_v<
                          typename Resource::Family,
                          Pullback::SphereAnimation::BurstRippleParams>)
          if (Resource::KEY.view() == op.instance)
            static_cast<In::RippleClockSnapshot &>(
                mutable_state<In::Op::RipplePhaseState>(program, index)) =
                own.resources.template get<Resource::KEY>();
      });
    } else if (id == "sphere.displace.symmetric.v2" ||
               id == "sphere.displace.rossby-haurwitz.v1" ||
               id == "sphere.displace.potential.v1") {
      own.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::SURFACE)
          if (Resource::KEY.view() == op.instance)
            mutable_state<In::Op::FieldPhaseState>(program, index).phase =
                own.resources.template get<Resource::KEY>().phase;
      });
    } else if (id == "warp.affine.v4" || id == "warp.wave-shear.v2" ||
               id == "warp.vortex.v2" || id == "warp.mirror-tile.v2") {
      const auto key = std::string_view(op.instance) == "warp2" ? "inner_warp"
                       : std::string_view(op.instance) == "warp1" ||
                               std::string_view(op.instance) == "warp"
                           ? "outer_warp"
                           : op.instance;
      own.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::WARP) {
          if (Resource::KEY.view() != key)
            return;
          const auto &clock = own.resources.template get<Resource::KEY>();
          if (id == "warp.affine.v4") {
            if constexpr (requires { clock.rotation; })
              mutable_state<In::Op::AffineClockState>(program, index) = {
                  clock.phase, clock.rotation};
          } else
            mutable_state<In::Op::WarpPhaseState>(program, index).phase =
                clock.phase;
        }
      });
    } else if (id == "sphere.displace.curl.v2" ||
               id == "sphere.displace.direct.v2") {
      if constexpr (!std::is_void_v<
                        typename FX::Params::template Family<"surface">>) {
        auto &state = mutable_state<In::Op::NoisePhaseState>(program, index);
        const auto &resource = own.resources.template get<"surface">();
        state.phase = resource.phase;
        state.noise = *resource.noise;
      }
    } else if (id == "sample.grid.v4") {
      auto &state = mutable_state<In::Op::SourceClockState>(program, index);
      const auto &resource = own.resources.template get<"source">();
      state = {resource.primary, resource.secondary, resource.angle};
    } else if (id == "colorize.generated-palette.v4") {
      mutable_state<In::Op::ColorClockState>(program, index).oscillation_phase =
          own.palette_oscillation_phase;
    }
  }
  In::FrameContext context;
  context.palettes.fill(own.palette);
  context.hue_rotation_lut = own.hue_rotation_lut;
  context.hue_rotation_gamut = own.hue_rotation_gamut;
  std::vector<Pullback::Modulation::Endpoint> endpoints;
  own.params.visit([&]<typename Resource>(const auto &) {
    if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
      for (size_t index = 0; index < program.ops().size(); ++index)
        if (std::string_view(program.ops()[index].instance) ==
            Resource::KEY.view()) {
          const auto key =
              mutable_state<In::Op::ModulationClockState>(program, index).key;
          endpoints.push_back(
              {key, own.resources.template get<Resource::KEY>().modulation});
        }
    }
  });
  context.modulations = endpoints;
  program.prepare(context);
  const auto prepared = FX::RenderPipeline::prepare(own);
  size_t visible = 0;
  for (const auto &view : sweep_views()) {
    const auto expected = FX::shade(view, prepared);
    const auto actual = program.evaluate(view, context);
    expect_color_within(actual, expected, channel_tolerance, alpha_tolerance);
    visible += expected.alpha > 0.0f &&
               (expected.color.r || expected.color.g || expected.color.b);
  }
  HS_EXPECT_GT(visible, 0u);
  if constexpr (std::is_same_v<FX, GridWorms<SMALL_W, SMALL_H>>) {
    auto blocks = std::make_unique<TemporalBlocks>();
    program.bind_temporal_storage(blocks->history, blocks->column);
    program.prepare(context);
    program.capture_history();
    auto &twist = mutable_state<In::Op::TwistOrbitState>(program, 0);
    const auto &twist_params =
        *reinterpret_cast<const In::Op::TwistChainParams *>(
            program.param_block(0));
    In::Op::LensTwist::advance(twist, twist_params);
    ComposedFrameWhiteBox::advance(effect);
    auto next = ComposedFrameWhiteBox::frame(effect);
    next.outer_conjugate = own.outer_conjugate;
    const auto &clock = next.resources.template get<"lens1">();
    HS_EXPECT_EQ(twist.phi, clock.phi);
    HS_EXPECT_EQ(twist.lambda, clock.lambda);
    HS_EXPECT_TRUE(program.history_continuous());
    for (const float fraction : {0.0f, 0.125f, 0.5f, 0.875f, 1.0f}) {
      const auto column = Pullback::FrameInterval(own, next).at(fraction);
      const auto frame = FX::RenderPipeline::prepare(column);
      program.prepare_column(context, context, fraction);
      for (const auto &view : sweep_views()) {
        const auto expected = FX::shade(view, frame);
        const auto actual =
            program.evaluate_column(view, context, context, fraction);
        expect_color_within(actual, expected, channel_tolerance,
                            alpha_tolerance);
      }
    }
  }
}

} // namespace hs_test::composed_chain_tests
