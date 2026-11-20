/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Document-built chain and composed frame parity fixtures.
 */
#pragma once

#include "tests/test_composed_effect.h"
#include "tests/shader_chain_fixture.h"

namespace hs_test::composed_chain_tests {

template <typename FX, typename BeforeCompare, typename AfterStatic>
void verify_export(size_t preset, int channel_tolerance, float alpha_tolerance,
                   BeforeCompare before_compare, AfterStatic after_static) {
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
  const auto slots = document_slots<FX>(*chain);
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
  using Selector = typename FX::Params::template Family<"source">;
  if constexpr (requires { typename Selector::PlanarWeights; })
    HS_EXPECT_TRUE(effect.selectPreset(preset));
  else
    ComposedFrameWhiteBox::set_params(effect,
                                      effects_tests::preset_params<FX>(preset));
  for (int frame = 0; frame < 3; ++frame)
    ComposedFrameWhiteBox::advance(effect);
  auto own = ComposedFrameWhiteBox::frame(effect);
  before_compare(own);
  for (size_t index = 0; index < program.ops().size(); ++index) {
    const auto &op = program.ops()[index];
    const std::string_view id = op.op->operator_id;
    if (id == "sphere.rotate.v2" || id.starts_with("project.")) {
      auto &state = shader_chain_tests::state_as<In::Op::SpatialWalkState>(
          program, index);
      state.wander = (id == "sphere.rotate.v2" ? own.outer_conjugate
                                               : own.projection_conjugate)
                         .conjugate();
      // The composed conjugate already carries the tilted base.
      if (id.starts_with("project.")) {
        auto &frame = shader_chain_tests::param_as<In::Op::ProjectChainParams>(
                          program, index)
                          .frame;
        if (frame ==
            static_cast<uint8_t>(In::Op::ProjectionFrame::SPIN_WANDER_TILTED))
          frame = static_cast<uint8_t>(In::Op::ProjectionFrame::SPIN_WANDER);
      }
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
            shader_chain_tests::state_as<In::Op::TwistOrbitState>(
                program, index) = {phase.phi, phase.lambda};
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
                shader_chain_tests::state_as<In::Op::RipplePhaseState>(program,
                                                                       index)) =
                own.resources.template get<Resource::KEY>();
      });
    } else if (id == "sphere.displace.symmetric.v2" ||
               id == "sphere.displace.rossby-haurwitz.v1" ||
               id == "sphere.displace.potential.v1") {
      own.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::SURFACE)
          if (Resource::KEY.view() == op.instance)
            shader_chain_tests::state_as<In::Op::FieldPhaseState>(program,
                                                                  index)
                .phase = own.resources.template get<Resource::KEY>().phase;
      });
    } else if (id == "warp.affine.v4" || id == "warp.wave-shear.v2" ||
               id == "warp.vortex.v2" || id == "warp.mirror-tile.v2" ||
               id == "warp.vector-noise.v2") {
      const auto slot =
          std::find_if(slots.begin(), slots.end(), [&](const auto &entry) {
            return entry.label == op.instance;
          });
      HS_EXPECT_TRUE(slot != slots.end());
      if (slot == slots.end())
        return;
      const auto key = document_resource_key(*slot);
      own.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::WARP) {
          if (Resource::KEY.view() != key)
            return;
          const auto &clock = own.resources.template get<Resource::KEY>();
          if (id == "warp.affine.v4") {
            if constexpr (requires { clock.rotation; })
              shader_chain_tests::state_as<In::Op::AffineClockState>(
                  program, index) = {clock.phase, clock.rotation};
          } else if (id == "warp.vector-noise.v2") {
            if constexpr (requires { clock.noise; }) {
              auto &state =
                  shader_chain_tests::state_as<In::Op::NoisePhaseState>(program,
                                                                        index);
              state.phase = clock.phase;
              state.noise = *clock.noise;
            }
          } else
            shader_chain_tests::state_as<In::Op::WarpPhaseState>(program, index)
                .phase = clock.phase;
        }
      });
    } else if (id == "sphere.displace.curl.v2" ||
               id == "sphere.displace.direct.v2") {
      own.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::SURFACE &&
                      (std::is_same_v<typename Resource::Family,
                                      Pullback::AuthoredDirectSurfaceParams> ||
                       Resource::KEY.view() == "surface"))
          if (Resource::KEY.view() == "surface" ||
              Resource::KEY.view() == op.instance) {
            auto &state = shader_chain_tests::state_as<In::Op::NoisePhaseState>(
                program, index);
            const auto &clock = own.resources.template get<Resource::KEY>();
            state.phase = clock.phase;
            state.noise = *clock.noise;
          }
      });
    } else if (id == "sample.planar.v1") {
      auto &state = shader_chain_tests::state_as<In::Op::PlanarSourceState>(
          program, index);
      if constexpr (FX::Params::template HAS<"source_grid">) {
        const auto &grid = own.resources.template get<"source_grid">();
        state.grid = {grid.primary, grid.secondary, grid.angle};
      }
      if constexpr (FX::Params::template HAS<"source_twin_wave">) {
        const auto &wave = own.resources.template get<"source_twin_wave">();
        state.twin_wave = {wave.primary, wave.secondary, wave.angle};
      }
      if constexpr (FX::Params::template HAS<"source_projected_noise">) {
        const auto &noise =
            own.resources.template get<"source_projected_noise">();
        state.noise.phase = noise.noise_time;
        state.noise.noise = *noise.noise;
      }
    } else if (id == "sample.grid.v4" || id == "sample.twin-wave.v3") {
      const auto &resource = own.resources.template get<"source">();
      if constexpr (requires { resource.primary; }) {
        auto &state = shader_chain_tests::state_as<In::Op::SourceClockState>(
            program, index);
        state = {resource.primary, resource.secondary, resource.angle};
      }
    } else if (id == "colorize.generated-palette.v4") {
      shader_chain_tests::state_as<In::Op::ColorClockState>(program, index)
          .oscillation_phase = own.palette_oscillation_phase;
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
              shader_chain_tests::state_as<In::Op::ModulationClockState>(
                  program, index)
                  .key;
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
  after_static(program, context, own, effect, channel_tolerance,
               alpha_tolerance);
}

} // namespace hs_test::composed_chain_tests
