/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Chain-interpreter core tests.
 *
 * Catalog regen: the golden at tests/data/shader_chain_catalog.json is
 * rewritten by the native build target:
 *   cmake --build --preset tests --target regenerate_shader_chain_catalog
 *
 * The golden's block sizes and alignments are the native ABI; the wasm32
 * catalog (regenerate_engine_catalog) is a separate file and is not
 * interchangeable with it.
 */
#pragma once

#include "core/math/mobius.h"
#include <bit>
#include <cstdio>
#include <cstring>
#include <memory>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "core/render/pullback/catalog_export.h"
#include "core/render/pullback/interpreter.h"
#include "effects/AlienCore.h"
#include "tests/composed_frame_fixture.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "workbench/shader/chain_host.h"

namespace hs_test {
namespace shader_chain_tests {

namespace PB = Pullback;
namespace In = Pullback::Interp;

// Admission-rule constants mirrored by scripts/shader_workbench.mjs.
static_assert(PB::Lens::MobiusLensParams::MOBIUS_MIN_DET_SQ == 1e-6f);
static_assert(PB::Warp::CURL_VECTOR_COMPONENT_MAX == 4.0f);
static_assert(In::Op::curl_intervals(In::Op::CurlIntegrator::EULER1) == 1);
static_assert(In::Op::curl_intervals(In::Op::CurlIntegrator::MIDPOINT2) == 2);
static_assert(In::Op::curl_intervals(In::Op::CurlIntegrator::MIDPOINT4) == 4);

#include "tests/shader_chain/fixtures.h"
#include "tests/shader_chain/effect_fixture.h"
#include "tests/shader_chain/slice_mirror.h"
#include "tests/shader_chain/lifecycle.h"
#include "tests/shader_chain/program_cases.h"
#include "tests/shader_chain/lifecycle_cases.h"
#include "tests/shader_chain/sphere_ops.h"
#include "tests/shader_chain/warp_ops.h"
#include "tests/shader_chain/projection_ops.h"
#include "tests/shader_chain/field_ops.h"
#include "tests/shader_chain/sample_ops.h"
#include "tests/shader_chain/palette_ops.h"
#include "tests/shader_chain/effect_cases.h"
/** @brief Layer caches retain independent fields across color-only edits. */
inline void test_shader_chain_modulation_bake_cache() {
  static_assert(In::instance_hash("mod-56428", "color.modulate.v1") ==
                In::instance_hash("mod-207280", "color.modulate.v1"));
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  WB::modulation_params(effect, "mod-56428");
  WB::modulation_params(effect, "mod-207280");
  auto &first = WB::modulation_params(effect, "mod-56428").modulation;
  auto &second = WB::modulation_params(effect, "mod-207280").modulation;
  first.mechanic.amount = .5f;
  second.mechanic.amount = .25f;
  ++second.seed;
  auto ctx = WB::frame_context(effect);
  HS_EXPECT_EQ(ctx.modulations.size(), 2u);
  const auto key_a = WB::modulation_clocks(effect, "mod-56428").key;
  const auto key_b = WB::modulation_clocks(effect, "mod-207280").key;
  auto &cache_a = WB::layer_resources(effect, key_a);
  auto &cache_b = WB::layer_resources(effect, key_b);
  HS_EXPECT_TRUE(&cache_a != &cache_b);
  HS_EXPECT_TRUE(ctx.modulations[0].prepared != ctx.modulations[1].prepared);
  const auto A = cache_a.modulation.bake_generation;
  const auto B = cache_b.modulation.bake_generation;
  HS_EXPECT_TRUE(cache_a.modulation.prepared.cube.active);
  HS_EXPECT_TRUE(cache_b.modulation.prepared.cube.active);
  HS_EXPECT_TRUE(cache_a.modulation.table.data() !=
                 cache_b.modulation.table.data());
  first.mechanic.amount = .75f;
  second.mechanic.kind = ColorModulation::Mechanic::OPACITY;
  WB::color_params(effect).palette_chroma = .1f;
  WB::advance_without_render(effect);
  WB::frame_context(effect);
  HS_EXPECT_EQ(cache_a.modulation.bake_generation, A);
  HS_EXPECT_EQ(cache_b.modulation.bake_generation, B);
  first.noise.scale = 3.0f;
  WB::frame_context(effect);
  HS_EXPECT_EQ(cache_a.modulation.bake_generation, A + 1);
  HS_EXPECT_EQ(cache_b.modulation.bake_generation, B);
  ++second.seed;
  WB::frame_context(effect);
  HS_EXPECT_EQ(cache_a.modulation.bake_generation, A + 1);
  HS_EXPECT_EQ(cache_b.modulation.bake_generation, B + 1);
  effect.draw_frame();
  const auto &previous = WB::previous_context(effect);
  HS_EXPECT_EQ(previous.modulations.size(), 2u);
  for (size_t index = 0; index < previous.modulations.size(); ++index) {
    const auto *retained = previous.modulations[index].prepared;
    HS_EXPECT_TRUE(retained != nullptr && retained->cube.active);
    HS_EXPECT_TRUE(retained->cube.data !=
                   ctx.modulations[index].prepared->cube.data);
  }
  effect.advance_display();
}

/** @brief Reorder preserves each layer's phase and disabled clocks advance. */
inline void test_shader_chain_modulation_clock_lifecycle() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  WB::modulation_params(effect, "first");
  WB::modulation_params(effect, "second");
  auto &first = WB::modulation_params(effect, "first").modulation;
  auto &second = WB::modulation_params(effect, "second").modulation;
  first.source = second.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  first.evaluation = second.evaluation = Pattern::Evaluation::DIRECT;
  first.sine.speed = .125f;
  second.sine.speed = .25f;
  WB::advance_without_render(effect);
  HS_EXPECT_NEAR(
      WB::modulation_clocks(effect, "first").modulation.phase.phase_turns,
      .125f, 0);
  HS_EXPECT_NEAR(
      WB::modulation_clocks(effect, "second").modulation.phase.phase_turns,
      .25f, 0);
  first.enabled = false;
  WB::advance_without_render(effect);
  HS_EXPECT_NEAR(
      WB::modulation_clocks(effect, "first").modulation.phase.phase_turns, .25f,
      0);
  HS_EXPECT_NEAR(
      WB::modulation_clocks(effect, "second").modulation.phase.phase_turns, .5f,
      0);
  const auto saved_a = WB::modulation_params(effect, "first");
  const auto saved_b = WB::modulation_params(effect, "second");
  const In::ChainEntryRequest reordered[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"},
      {"second", "color.modulate.v1"},
      {"first", "color.modulate.v1"}};
  HS_EXPECT_EQ(effect.set_chain(reordered).code, In::ChainStatus::OK);
  WB::modulation_params(effect, "first") = saved_a;
  WB::modulation_params(effect, "second") = saved_b;
  WB::advance_without_render(effect);
  HS_EXPECT_NEAR(
      WB::modulation_clocks(effect, "first").modulation.phase.phase_turns,
      .375f, 0);
  HS_EXPECT_NEAR(
      WB::modulation_clocks(effect, "second").modulation.phase.phase_turns,
      .75f, 0);
  effect.setAnimationsPaused(true);
  auto &paused = WB::modulation_params(effect, "first").modulation;
  paused.source = Pattern::SourceKind::DISPLACEMENT;
  WB::frame_context(effect);
  HS_EXPECT_EQ(WB::modulation_clocks(effect, "first").modulation.source,
               Pattern::SourceKind::DISPLACEMENT);
  paused.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  WB::frame_context(effect);
  HS_EXPECT_NEAR(
      WB::modulation_clocks(effect, "first").modulation.phase.phase_turns, 0.0f,
      0);
}

/** @brief Source and mechanic combinations use the selected layer's scalar. */
inline void test_shader_chain_modulation_combinations() {
  using Kind = ColorModulation::Mechanic;
  const PB::FieldSample FIELD{.413f, .6f, math::X_AXIS, 4.5f};
  const auto ctx = shared_resources().context();
  In::Op::GeneratedPaletteParams params;
  const auto palette = In::Op::prepare_generated_palette(ctx, params, {});
  for (auto source :
       {Pattern::SourceKind::NOISE, Pattern::SourceKind::DIRECTIONAL_SINE,
        Pattern::SourceKind::SPHERE_INTERFERENCE,
        Pattern::SourceKind::SPHERE_SMOOTH_GRID,
        Pattern::SourceKind::DISPLACEMENT}) {
    PB::Modulation::Prepared spatial;
    spatial.config.source = source;
    spatial.evaluate = +[](const PB::Modulation::Prepared &,
                           const Pattern::SphereDirection &) { return .5f; };
    for (auto mechanic :
         {Kind::HUE_ROTATION, Kind::POSITIVE_LIGHTNESS_LIFT, Kind::OPACITY}) {
      HS_CONTEXT("source / mechanic", int(source), int(mechanic));
      PB::Color::ModulationState layer;
      layer.current = &spatial;
      layer.source = source;
      layer.mechanic = {mechanic, .4f};
      layer.enabled = true;
      auto input = PB::Color::apply_generated_palette(FIELD, palette);
      Color4 expected = input.color;
      const float SCALAR =
          source == Pattern::SourceKind::DISPLACEMENT ? 4.5f : .5f;
      if (mechanic == Kind::HUE_ROTATION)
        expected.color = sample_hue_rotation_lut(
            palette.hue_rotation, FIELD.value, math::wrap_t(.4f * SCALAR));
      else
        expected = ColorModulation::apply(expected, SCALAR, layer.mechanic);
      const auto actual =
          PB::Color::finish(PB::Color::apply_modulation(input, layer));
      HS_EXPECT_TRUE(color4_identical(actual, expected));
      layer.enabled = false;
      HS_EXPECT_TRUE(color4_identical(
          PB::Color::finish(PB::Color::apply_modulation(input, layer)),
          input.color));
      layer.enabled = true;
      layer.mechanic.amount = 0;
      HS_EXPECT_TRUE(color4_identical(
          PB::Color::finish(PB::Color::apply_modulation(input, layer)),
          input.color));
    }
  }
}

/** @brief Temporal columns finish one blend of deferred color endpoints. */
inline void test_shader_chain_modulation_temporal_endpoints() {
  constexpr size_t BLOCKS = In::CHAIN_ARENA_BYTES / sizeof(std::max_align_t);
  std::vector<std::max_align_t> history(BLOCKS), column(BLOCKS);
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  program.bind_temporal_storage(reinterpret_cast<uint8_t *>(history.data()),
                                reinterpret_cast<uint8_t *>(column.data()));
  const In::ChainEntryRequest CHAIN[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"},
      {"lift", "color.modulate.v1"},
      {"hue", "color.modulate.v1"}};
  HS_EXPECT_EQ(program.compile(CHAIN).code, In::ChainStatus::OK);
  auto &lift = param_as<In::Op::ModulationParams>(program, 4).modulation;
  auto &hue = param_as<In::Op::ModulationParams>(program, 5).modulation;
  lift.evaluation = hue.evaluation = Pattern::Evaluation::DIRECT;
  lift.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, .6f};
  hue.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, .2f};
  std::array<PB::Modulation::Prepared, 2> previous, current;
  std::array<PB::Modulation::Endpoint, 2> old_endpoints, new_endpoints;
  for (size_t index = 0; index < 2; ++index) {
    auto &old = previous[index];
    old.config = index == 0 ? lift : hue;
    old.state.phase.phase_turns = .25f;
    old.evaluate = +[](const PB::Modulation::Prepared &p,
                       const Pattern::SphereDirection &) {
      return 2 * p.state.phase.phase_turns - 1;
    };
    current[index] = old;
    current[index].state.phase.phase_turns = .75f;
    const auto key =
        state_as<In::Op::ModulationClockState>(program, index + 4).key;
    old_endpoints[index] = {key, &old};
    new_endpoints[index] = {key, &current[index]};
  }
  auto old_ctx = shared_resources().context();
  old_ctx.modulations = old_endpoints;
  old_ctx.hue_rotation_lut = nullptr;
  old_ctx.palettes[0] = old_ctx.palettes[1];
  auto ctx = shared_resources().context();
  ctx.modulations = new_endpoints;
  ctx.hue_rotation_lut = nullptr;
  program.prepare(old_ctx);
  program.capture_history();
  for (float fraction : {0.0f, .375f, 1.0f}) {
    program.prepare_column(ctx, old_ctx, fraction);
    const auto VIEW = math::Vector(.3f, .7f, .2f).normalized();
    const auto field = run_entry<PB::FieldSample>(
        program, ctx, 2, projected_input(program, ctx, VIEW));
    const auto BASE = old_ctx.palettes[0]
                          ->get(field.value)
                          .lerp(ctx.palettes[0]->get(field.value), fraction);
    HueRotateBase base = make_hue_rotate_base(BASE);
    base.lab.L += (1.0f - base.lab.L) * (.3f * fraction);
    auto expected =
        hue_rotate_lut_gamut(base, math::wrap_t(-.1f + .2f * fraction));
    expected.alpha *= field.coverage;
    const auto actual = program.evaluate_column(VIEW, ctx, old_ctx, fraction);
    HS_EXPECT_NEAR(actual.color.r, expected.color.r, 1.0f);
    HS_EXPECT_NEAR(actual.color.g, expected.color.g, 1.0f);
    HS_EXPECT_NEAR(actual.color.b, expected.color.b, 1.0f);
    HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
  }
}

/** @brief Each modulation resolves its own endpoint regardless of view order. */
inline void test_shader_chain_modulation_endpoint_identity() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest CHAIN[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"},
      {"first", "color.modulate.v1"},
      {"second", "color.modulate.v1"}};
  HS_EXPECT_EQ(program.compile(CHAIN).code, In::ChainStatus::OK);
  std::array<PB::Modulation::Prepared, 2> prepared;
  std::array<PB::Modulation::Endpoint, 2> endpoints;
  for (size_t index = 0; index < 2; ++index) {
    auto &config =
        param_as<In::Op::ModulationParams>(program, index + 4).modulation;
    config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
    config.evaluation = Pattern::Evaluation::DIRECT;
    config.mechanic = {ColorModulation::Mechanic::OPACITY,
                       index == 0 ? 1.0f : 0.5f};
    prepared[index].config = config;
    prepared[index].state.phase.phase_turns = index == 0 ? -.5f : .5f;
    prepared[index].evaluate = +[](const PB::Modulation::Prepared &source,
                                   const Pattern::SphereDirection &) {
      return source.state.phase.phase_turns;
    };
    const auto key =
        state_as<In::Op::ModulationClockState>(program, index + 4).key;
    endpoints[1 - index] = {key, &prepared[index]};
  }
  auto context = shared_resources().context();
  context.modulations = endpoints;
  program.prepare(context);
  const auto VIEW = math::Vector(.3f, .7f, .2f).normalized();
  const auto FIELD = run_entry<PB::FieldSample>(
      program, context, 2, projected_input(program, context, VIEW));
  const float EXPECTED = FIELD.coverage * .25f * .875f;
  HS_EXPECT_GT(EXPECTED, 0.0f);
  HS_EXPECT_NEAR(program.evaluate(VIEW, context).alpha, EXPECTED, 1e-7f);
}

// The wire spellings are the JS contract.
inline const char *expected_chain_status_name(In::ChainStatus status) {
  switch (status) {
  case In::ChainStatus::OK:
    return "OK";
  case In::ChainStatus::NOT_CHAIN_EFFECT:
    return "NOT_CHAIN_EFFECT";
  case In::ChainStatus::MALFORMED_PAYLOAD:
    return "MALFORMED_PAYLOAD";
  case In::ChainStatus::EMPTY:
    return "EMPTY";
  case In::ChainStatus::TOO_LONG:
    return "TOO_LONG";
  case In::ChainStatus::UNKNOWN_OPERATOR:
    return "UNKNOWN_OPERATOR";
  case In::ChainStatus::DUPLICATE_INSTANCE:
    return "DUPLICATE_INSTANCE";
  case In::ChainStatus::MALFORMED_INSTANCE:
    return "MALFORMED_INSTANCE";
  case In::ChainStatus::ENTRY_FAMILY:
    return "ENTRY_FAMILY";
  case In::ChainStatus::EXIT_FAMILY:
    return "EXIT_FAMILY";
  case In::ChainStatus::CARRIER_MISMATCH:
    return "CARRIER_MISMATCH";
  case In::ChainStatus::ARENA_OVERFLOW:
    return "ARENA_OVERFLOW";
  case In::ChainStatus::PARAM_OVERFLOW:
    return "PARAM_OVERFLOW";
  case In::ChainStatus::MIGRATE_FAILED:
    return "MIGRATE_FAILED";
  }
  return nullptr;
}

inline void test_shader_chain_status_names() {
  for (unsigned raw = 0; raw <= UINT8_MAX; ++raw) {
    const auto status = static_cast<In::ChainStatus>(raw);
    const char *const expected = expected_chain_status_name(status);
    HS_EXPECT_STREQ(In::chain_status_name(status),
                    expected ? expected : "UNKNOWN");
  }
}

inline void test_shader_chain_snapshot_roundtrip() {
  using WB = ShaderChainWhiteBox;
  ChainSnapshot saved;
  std::array<std::array<Pixel, 96 * 20>, 5> frames;
  reset_globals();
  {
    WB::FX effect;
    effect.init();
    const In::ChainEntryRequest topology[] = {
        {"camera", "sphere.rotate.v2"},
        {"project", "project.stereographic.v2"},
        {"noise", "warp.vector-noise.v2"},
        {"outer", "warp.affine.v3"},
        {"sample", "sample.grid.v3"},
        {"colorize", "colorize.generated-palette.v4"},
        {"modulation", "color.modulate.v1"}};
    HS_EXPECT_EQ(effect.set_chain(topology).code, In::ChainStatus::OK);
    auto &program = WB::program(effect);
    param_as<In::Op::RotateChainParams>(program, 0).wander = 0.8f;
    param_as<In::Op::RotateChainParams>(program, 0).spin_rate = 0.001f;
    param_as<In::Op::GridSampleParams>(program, 4).speed = 0.02f;
    auto &layer = WB::modulation_params(effect).modulation;
    auto &color = WB::color_params(effect);
    layer.mechanic.amount = 0.7f;
    layer.noise.scale = 1.7f;
    color.palette_chroma = 0.43f;
    color.phase_oscillation_speed = 0.004f;
    layer.noise.speed = -0.0007f;
    for (int frame = 0; frame < 713; ++frame)
      WB::advance_without_render(effect);
    effect.setAnimationsPaused(true);
    saved = effect.snapshot();
    HS_EXPECT_TRUE(GeneratedPaletteBank::valid_snapshot(*saved.palette_bank));
    for (size_t index = 0; index < program.ops().size(); ++index) {
      const auto &op = program.ops()[index];
      std::vector<std::max_align_t> storage(
          (op.op->runtime.state.size + sizeof(std::max_align_t) - 1) /
          sizeof(std::max_align_t));
      op.op->runtime.init(storage.data(),
                          {op.instance, op.op->operator_id, op.stable_hash});
      const auto captured =
          op.op->runtime.capture_state(program.state_block(index));
      HS_EXPECT_TRUE(op.op->runtime.restore_state(storage.data(), captured));
      op.op->runtime.destroy(storage.data());
    }
    HS_EXPECT_TRUE(saved.palette_bank->cycles[0].next_sequence > 2);
    HS_EXPECT_TRUE(saved.palette_bank->cycles[1].display_dirty);
    for (size_t frame = 0; frame < frames.size(); ++frame) {
      effect.draw_frame();
      effect.advance_display();
      int lit = 0;
      for (int y = 0; y < 20; ++y)
        for (int x = 0; x < 96; ++x) {
          const auto pixel = effect.get_pixel(x, y);
          frames[frame][y * 96 + x] = pixel;
          lit += pixel.r != 0 || pixel.g != 0 || pixel.b != 0;
        }
      HS_EXPECT_GT(lit, 200);
    }
  }
  reset_globals();
  {
    WB::FX effect;
    effect.init();
    HS_EXPECT_EQ(effect.restore_snapshot(saved),
                 ChainSnapshotRestoreResult::APPLIED);
    HS_EXPECT_TRUE(effect.animations_paused());
    for (size_t frame = 0; frame < frames.size(); ++frame) {
      effect.draw_frame();
      effect.advance_display();
      for (int y = 0; y < 20; ++y)
        for (int x = 0; x < 96; ++x) {
          const auto expected = frames[frame][y * 96 + x];
          const auto actual = effect.get_pixel(x, y);
          HS_EXPECT_EQ(actual.r, expected.r);
          HS_EXPECT_EQ(actual.g, expected.g);
          HS_EXPECT_EQ(actual.b, expected.b);
        }
    }
  }
}

inline void test_shader_chain_snapshot_refusals() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  WB::modulation_params(effect);
  const auto live = effect.snapshot();
  // Candidates derive from a state that differs from the live one in every
  // section, so any section committed before a refusal is observable.
  param_as<In::Op::RotateChainParams>(WB::program(effect), 0).wander = 0.6f;
  WB::modulation_params(effect).modulation.mechanic.amount = .75f;
  WB::modulation_params(effect).modulation.noise.scale = 2.0f;
  for (int frame = 0; frame < 64; ++frame)
    WB::advance_without_render(effect);
  const auto saved = effect.snapshot();
  HS_EXPECT_EQ(effect.restore_snapshot(live),
               ChainSnapshotRestoreResult::APPLIED);
  const auto walk_time = [](const ChainSnapshot &snapshot) {
    if (!snapshot.runtime || snapshot.runtime->empty())
      return std::numeric_limits<float>::quiet_NaN();
    const auto *walk =
        std::get_if<In::SpatialWalkSnapshot>(&(*snapshot.runtime)[0].state);
    return walk ? walk->walk_time : std::numeric_limits<float>::quiet_NaN();
  };
  HS_EXPECT_FALSE(
      float_identical(saved.parameters[0].value, live.parameters[0].value));
  HS_EXPECT_NE(saved.palette_bank->cycles[0].frame,
               live.palette_bank->cycles[0].frame);
  HS_EXPECT_FALSE(float_identical(walk_time(saved), walk_time(live)));
  auto companion = saved;
  companion.modulation[0].recipe.mechanic.amount = 3.0f;
  HS_EXPECT_EQ(effect.restore_snapshot(companion),
               ChainSnapshotRestoreResult::APPLIED);
  HS_EXPECT_EQ(effect.snapshot().modulation[0].recipe.mechanic.amount, 3.0f);
  HS_EXPECT_EQ(effect.restore_snapshot(live),
               ChainSnapshotRestoreResult::APPLIED);
  const auto same_modulation = [](const ChainSnapshot &a,
                                  const ChainSnapshot &b) {
    HS_EXPECT_SIZE_OR_RETURN(a.modulation, b.modulation.size());
    for (size_t index = 0; index < a.modulation.size(); ++index) {
      HS_EXPECT_EQ(a.modulation[index].instance, b.modulation[index].instance);
      In::Op::ModulationParams left, right;
      left.modulation = a.modulation[index].recipe;
      right.modulation = b.modulation[index].recipe;
      HS_EXPECT_EQ(left.modulation.seed, right.modulation.seed);
      HS_EXPECT_EQ(left.modulation.version, right.modulation.version);
      for (size_t field = 0; field < left.EXTRA_SCHEMA.size(); ++field) {
        if (left.EXTRA_SCHEMA[field].enum_count)
          HS_EXPECT_EQ(
              *static_cast<const uint8_t *>(left.extra_address(field)),
              *static_cast<const uint8_t *>(right.extra_address(field)));
        else
          HS_EXPECT_TRUE(float_identical(
              *static_cast<const float *>(left.extra_address(field)),
              *static_cast<const float *>(right.extra_address(field))));
      }
    }
  };
  const auto generation = effect.getParameterSchemaGeneration();
  const auto refused = [&](ChainSnapshot candidate,
                           ChainSnapshotRestoreResult expected) {
    HS_EXPECT_EQ(effect.restore_snapshot(candidate), expected);
    HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), generation);
    const auto after = effect.snapshot();
    same_modulation(after, live);
    HS_EXPECT_SIZE_OR_RETURN(after.parameters, live.parameters.size());
    for (size_t index = 0; index < live.parameters.size(); ++index) {
      HS_EXPECT_TRUE(after.parameters[index].name ==
                     live.parameters[index].name);
      HS_EXPECT_TRUE(float_identical(after.parameters[index].value,
                                     live.parameters[index].value));
    }
    HS_EXPECT_EQ(after.palette_bank->cycles[0].frame,
                 live.palette_bank->cycles[0].frame);
    HS_EXPECT_TRUE(float_identical(walk_time(after), walk_time(live)));
  };
  HS_EXPECT_SIZE_OR_RETURN(saved.modulation, 1u);
  HS_EXPECT_NE(saved.modulation[0].recipe.mechanic.amount,
               live.modulation[0].recipe.mechanic.amount);
  auto candidate = saved;
  candidate.modulation[0].recipe.mechanic.amount = 10.0f;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.modulation[0].recipe.noise.speed = 100.0f;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.modulation[0].recipe.sine.offset = 1000.0f;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.modulation[0].recipe.sine.axis.x = 100.0f;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.modulation[0].recipe.noise.variant =
      static_cast<Pattern::NoiseVariant>(7);
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.parameters.push_back(
      {saved.modulation[0].instance + ".modulation-amount", .75f});
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.modulation.push_back({"unknown", saved.modulation[0].recipe});
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.modulation.clear();
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.schema_version = 1;
  refused(candidate, ChainSnapshotRestoreResult::UNSUPPORTED_VERSION);
  candidate = saved;
  candidate.parameters[0].value = std::numeric_limits<float>::quiet_NaN();
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.parameters.push_back(candidate.parameters[0]);
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.chain[0].operator_id = "invalid";
  refused(candidate, ChainSnapshotRestoreResult::INVALID_CHAIN);
  candidate = saved;
  candidate.runtime->pop_back();
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  std::get<In::SpatialWalkSnapshot>((*candidate.runtime)[0].state).wander =
      math::Quaternion(0, 0, 0, 0);
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.palette_bank->cycles[0].frame = GeneratedPaletteBank::FADE_FRAMES;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.chain.clear();
  refused(candidate, ChainSnapshotRestoreResult::INVALID_LENGTH);
  candidate = saved;
  candidate.chain.resize(In::MAX_CHAIN_OPS + 1, saved.chain[0]);
  refused(candidate, ChainSnapshotRestoreResult::INVALID_LENGTH);
  candidate = saved;
  candidate.parameters.resize(In::MAX_CHAIN_PARAMS + 1, saved.parameters[0]);
  refused(candidate, ChainSnapshotRestoreResult::INVALID_LENGTH);
  candidate = saved;
  candidate.runtime->resize(In::MAX_CHAIN_OPS + 1, (*saved.runtime)[0]);
  refused(candidate, ChainSnapshotRestoreResult::INVALID_LENGTH);
  candidate = saved;
  candidate.parameters.push_back({"camera.nope", 0.0f});
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.parameters.clear();
  candidate.parameters.push_back({"camera.wander", 2.0f});
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  for (float value : {4.0f, 0.5f}) {
    candidate = saved;
    candidate.parameters.clear();
    candidate.parameters.push_back({"sample.coverage-mode", value});
    refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  }
  candidate = saved;
  candidate.chain.insert(candidate.chain.begin() + 1,
                         {"lens", In::Op::LensMobius::ID});
  candidate.runtime.reset();
  candidate.parameters.clear();
  candidate.parameters.push_back({"lens.mobius-a-re", 0.0f});
  candidate.parameters.push_back({"lens.mobius-d-re", 0.0f});
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.chain[1].operator_id = In::Op::ProjectFoldedSinusoidal::ID;
  candidate.runtime.reset();
  candidate.parameters.clear();
  candidate.parameters.push_back({"sample.coverage-mode", 3.0f});
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.runtime->push_back((*candidate.runtime)[0]);
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  candidate.runtime->push_back({"unknown", (*candidate.runtime)[0].state});
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  (*candidate.runtime)[0].state = In::SourceClockSnapshot{};
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  std::get<In::SpatialWalkSnapshot>((*candidate.runtime)[0].state).direction =
      std::get<In::SpatialWalkSnapshot>((*candidate.runtime)[0].state).position;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  std::get<In::SpatialWalkSnapshot>((*candidate.runtime)[0].state).spin_phase =
      7.0f;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  for (auto &runtime : *candidate.runtime)
    if (auto *clock = std::get_if<In::SourceClockSnapshot>(&runtime.state))
      clock->angle = 7.0f;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  candidate = saved;
  for (auto &runtime : *candidate.runtime)
    if (auto *clock = std::get_if<In::ModulationClockSnapshot>(&runtime.state))
      clock->modulation.phase.phase_turns = 1.0f;
  refused(candidate, ChainSnapshotRestoreResult::INVALID_VALUE);
  In::Op::NoisePhaseState noise_state;
  HS_EXPECT_FALSE(In::RuntimeStateCodec<In::Op::NoisePhaseState>::restore(
      noise_state, In::NoiseClockSnapshot{1.0f, 1337}));
  In::Op::WarpPhaseState warp_state;
  HS_EXPECT_FALSE(In::RuntimeStateCodec<In::Op::WarpPhaseState>::restore(
      warp_state, In::PhaseClockSnapshot{1.0f}));
  In::Op::RipplePhaseState ripple_state;
  HS_EXPECT_FALSE(In::RuntimeStateCodec<In::Op::RipplePhaseState>::restore(
      ripple_state, In::RippleClockSnapshot{-1.0f}));
  In::Op::AffineClockState affine_state;
  HS_EXPECT_FALSE(In::RuntimeStateCodec<In::Op::AffineClockState>::restore(
      affine_state, In::AffineClockSnapshot{1.0f, 0.0f}));
  HS_EXPECT_FALSE(In::RuntimeStateCodec<In::Op::AffineClockState>::restore(
      affine_state, In::AffineClockSnapshot{0.0f, 7.0f}));
  In::Op::SphericalRingsState rings_state;
  const auto walk =
      std::get<In::SpatialWalkSnapshot>((*saved.runtime)[0].state);
  HS_EXPECT_FALSE(In::RuntimeStateCodec<In::Op::SphericalRingsState>::restore(
      rings_state, In::SphericalRingsSnapshot{walk, 7.0f}));
  candidate = saved;
  candidate.runtime.reset();
  candidate.palette_bank.reset();
  HS_EXPECT_EQ(effect.restore_snapshot(candidate),
               ChainSnapshotRestoreResult::APPLIED);
}

/** @brief Snapshot recipes and typed clocks retain matching source identities. */
inline void test_shader_chain_modulation_snapshot_identity() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  auto &controls = WB::modulation_params(effect);
  controls.modulation.noise.speed = .000125f;
  WB::advance_without_render(effect);
  const ChainSnapshot saved = effect.snapshot();
  const auto clock =
      [](const ChainSnapshot &snapshot) -> const In::ModulationClockSnapshot & {
    for (const auto &runtime : *snapshot.runtime)
      if (const auto *value =
              std::get_if<In::ModulationClockSnapshot>(&runtime.state))
        return *value;
    HS_CHECK(false, "missing test color clock");
    return std::get<In::ModulationClockSnapshot>(
        snapshot.runtime->back().state);
  };
  HS_EXPECT_EQ(clock(saved).modulation.phase.phase_turns, .000125f);
  auto mismatched = saved;
  for (auto &runtime : *mismatched.runtime)
    if (auto *value = std::get_if<In::ModulationClockSnapshot>(&runtime.state))
      value->modulation.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  HS_EXPECT_EQ(effect.restore_snapshot(mismatched),
               ChainSnapshotRestoreResult::INVALID_VALUE);
  HS_EXPECT_EQ(effect.restore_snapshot(saved),
               ChainSnapshotRestoreResult::APPLIED);
  auto uninitialized = saved;
  for (auto &runtime : *uninitialized.runtime)
    if (auto *value = std::get_if<In::ModulationClockSnapshot>(&runtime.state))
      value->modulation.initialized = false;
  HS_EXPECT_EQ(effect.restore_snapshot(uninitialized),
               ChainSnapshotRestoreResult::INVALID_VALUE);
  HS_EXPECT_EQ(effect.restore_snapshot(saved),
               ChainSnapshotRestoreResult::APPLIED);
  auto &changed = WB::modulation_params(effect).modulation;
  changed.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  changed.evaluation = Pattern::Evaluation::DIRECT;
  const ChainSnapshot switched = effect.snapshot();
  HS_EXPECT_EQ(clock(switched).modulation.source,
               Pattern::SourceKind::DIRECTIONAL_SINE);
  HS_EXPECT_EQ(clock(switched).modulation.phase.phase_turns, 0.0f);
  HS_EXPECT_TRUE(clock(switched).modulation.initialized);
  HS_EXPECT_EQ(effect.restore_snapshot(switched),
               ChainSnapshotRestoreResult::APPLIED);
}

inline int run_shader_chain_tests() {
  ModuleFixture fixture("shader_chain");
  HS_RUN(test_shader_chain_table_integrity);
  HS_RUN(test_shader_chain_state_restore_values);
  HS_RUN(test_shader_chain_table_behavior);
  HS_RUN(test_shader_chain_schema_and_field_ids);
  HS_RUN(test_shader_chain_schema_defaults_match_params);
  HS_RUN(test_shader_chain_instance_id_wellformed);
  HS_RUN(test_shader_chain_slot_and_hash_contract);
  HS_RUN(test_shader_chain_catalog_golden);
  HS_RUN(test_shader_chain_catalog_shape);
  HS_RUN(test_shader_chain_default_chain_renders);
  HS_RUN(test_shader_chain_param_address_channel);
  HS_RUN(test_shader_chain_parity_rotate_project);
  HS_RUN(test_shader_chain_parity_displace_curl);
  HS_RUN(test_shader_chain_parity_displace_direct);
  HS_RUN(test_shader_chain_parity_displace_ripple);
  HS_RUN(test_shader_chain_ripple_bursts);
  HS_RUN(test_shader_chain_twist_animation);
  HS_RUN(test_shader_chain_parity_lens_ops);
  HS_RUN(test_shader_chain_lens_controls);
  HS_RUN(test_shader_chain_parity_project_ops);
  HS_RUN(test_shader_chain_peirce_square_dispatch);
  HS_RUN(test_shader_chain_projection_frame_policy);
  HS_RUN(test_shader_chain_parity_project_hemispheres);
  HS_RUN(test_shader_chain_parity_field_ops);
  HS_RUN(test_shader_chain_parity_warp_affine_mirror);
  HS_RUN(test_shader_chain_parity_warp_wave_shear);
  HS_RUN(test_shader_chain_parity_warp_vortex);
  HS_RUN(test_shader_chain_parity_warp_vector_noise);
  HS_RUN(test_shader_chain_noise_instances_decorrelate);
  HS_RUN(test_shader_chain_parity_warp_polar_chart);
  HS_RUN(test_shader_chain_parity_warp_curl_flow);
  HS_RUN(test_shader_chain_parity_sample_variants);
  HS_RUN(test_shader_chain_parity_sample_twin_wave);
  HS_RUN(test_shader_chain_parity_sample_rings);
  HS_RUN(test_shader_chain_parity_sample_spherical_rings);
  HS_RUN(test_shader_chain_parity_sample_spiral);
  HS_RUN(test_shader_chain_parity_sample_lattice);
  HS_RUN(test_shader_chain_parity_sample_fractal);
  HS_RUN(test_shader_chain_parity_sample_tessellation);
  HS_RUN(test_shader_chain_large_finite_path_length);
  HS_RUN(test_shader_chain_composed_noise_domain);
  HS_RUN(test_shader_chain_noise_operator_domain);
  HS_RUN(test_shader_chain_parity_sample_projected_noise);
  HS_RUN(test_shader_chain_parity_sample_spherical_noise);
  HS_RUN(test_shader_chain_parity_colorize_variants);
  HS_RUN(test_shader_chain_palette_selection);
  HS_RUN(test_shader_chain_authored_palette_history);
  HS_RUN(test_shader_chain_combined_color_modulation);
  HS_RUN(test_shader_chain_combined_color_admission);
  HS_RUN(test_shader_chain_refusal_shape);
  HS_RUN(test_shader_chain_refusal_budget_overflows);
  HS_RUN(test_shader_chain_program_lifetime);
  HS_RUN(test_shader_chain_refusal_migrate_failed);
  HS_RUN(test_shader_chain_column_states_outlive_evaluation);
  HS_RUN(test_shader_chain_state_identity_migration);
  HS_RUN(test_shader_chain_state_continuity_slice);
  HS_RUN(test_shader_chain_operator_state_migration);
  HS_RUN(test_shader_chain_determinism);
  HS_RUN(test_shader_chain_param_names_and_budget);
  HS_RUN(test_shader_chain_composed_frame_parity);
  HS_RUN(test_shader_chain_kaleidoscope_strobe);
  HS_RUN(test_shader_chain_instantaneous_pixel_golden);
  HS_RUN(test_shader_chain_columns_follow_sweep_time);
  HS_RUN(test_shader_chain_segmented_matches_single_engine);
  HS_RUN(test_shader_chain_discrete_changes_reset_history);
  HS_RUN(test_shader_chain_fidelity_toggle_keeps_history);
  HS_RUN(test_shader_chain_palettes_blend_across_the_frame);
  HS_RUN(test_shader_chain_column_params_interpolate);
  HS_RUN(test_shader_chain_operator_states_sample_clocks);
  HS_RUN(test_shader_chain_effect_registers_params);
  HS_RUN(test_shader_chain_parameter_admission);
  HS_RUN(test_shader_chain_edge_distance_admission);
  HS_RUN(test_shader_chain_projection_plane_bound);
  HS_RUN(test_shader_chain_warp_plane_bound);
  HS_RUN(test_shader_chain_plane_growth_admission);
  HS_RUN(test_shader_chain_effect_rebind_generation);
  HS_RUN(test_shader_chain_effect_refusal_keeps_schema);
  HS_RUN(test_shader_chain_pause_semantics);
  HS_RUN(test_shader_chain_modulation_bake_cache);
  HS_RUN(test_shader_chain_modulation_clock_lifecycle);
  HS_RUN(test_shader_chain_modulation_combinations);
  HS_RUN(test_shader_chain_modulation_temporal_endpoints);
  HS_RUN(test_shader_chain_modulation_endpoint_identity);
  HS_RUN(test_shader_chain_status_names);
  HS_RUN(test_shader_chain_snapshot_roundtrip);
  HS_RUN(test_shader_chain_snapshot_refusals);
  HS_RUN(test_shader_chain_modulation_snapshot_identity);
  return fixture.result();
}

} // namespace shader_chain_tests
} // namespace hs_test
