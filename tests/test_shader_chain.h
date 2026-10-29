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
#include "tests/shader_chain/effect_cases.h"
/** @brief Spatial cache ignores palette/mechanic changes and tracks source inputs. */
inline void test_shader_chain_modulation_bake_cache() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  auto &color = WB::color_params(effect);
  color.modulation.mechanic.amount = .5f;
  auto ctx = WB::frame_context(effect);
  HS_EXPECT_TRUE(ctx.modulation != nullptr);
  HS_EXPECT_TRUE(ctx.modulation->cube.active);
  auto *lut = WB::pattern_lut(effect);
  constexpr int8_t POISON = -128;
  const auto poison = [&] {
    lut[0] = POISON;
    lut[Pattern::SphereLutView::SIZE - 1] = POISON;
  };
  const auto baked = [&](bool expected) {
    HS_EXPECT_EQ(lut[0] != POISON, expected);
    HS_EXPECT_EQ(lut[Pattern::SphereLutView::SIZE - 1] != POISON, expected);
  };
  poison();
  WB::frame_context(effect);
  baked(false);
  color.modulation.mechanic.amount = .25f;
  WB::frame_context(effect);
  baked(false);
  color.modulation.mechanic.kind = ColorModulation::Mechanic::OPACITY;
  WB::frame_context(effect);
  baked(false);
  const Pixel PALETTE = WB::palette_color(effect, .37f);
  color.palette_chroma = .1f;
  WB::advance_without_render(effect);
  WB::frame_context(effect);
  baked(false);
  HS_EXPECT_NE(WB::palette_color(effect, .37f), PALETTE);
  color.modulation.noise.scale = 3;
  WB::frame_context(effect);
  baked(true);
  poison();
  ++color.modulation.seed;
  WB::frame_context(effect);
  baked(true);
  poison();
  color.modulation.noise.speed = .001f;
  const float previous = WB::color_clocks(effect).modulation.phase.phase_turns;
  WB::program(effect).advance();
  HS_EXPECT_NE(WB::color_clocks(effect).modulation.phase.phase_turns, previous);
  WB::frame_context(effect);
  baked(true);
  poison();
  color.modulation.mechanic.amount = 0;
  ctx = WB::frame_context(effect);
  baked(false);
  HS_EXPECT_TRUE(ctx.hue_rotation_lut == nullptr);
}

/** @brief Source clocks retain same-source identity and reset on recipe changes. */
inline void test_shader_chain_modulation_clock_lifecycle() {
  namespace Mod = PB::Modulation;
  ColorModulation::Config config;
  config.enabled = false;
  config.noise.speed = -.125f;
  auto state = Mod::init(config);
  Mod::advance(state, config);
  HS_EXPECT_EQ(state.phase.phase_turns, .875f);
  config.mechanic.kind = ColorModulation::Mechanic::OPACITY;
  Mod::synchronize(state, config);
  HS_EXPECT_EQ(state.phase.phase_turns, .875f);
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.noise.scale = 3;
  Mod::synchronize(state, config);
  HS_EXPECT_EQ(state.phase.phase_turns, .875f);
  config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  config.sine.speed = .25f;
  Mod::synchronize(state, config);
  HS_EXPECT_EQ(state.phase.phase_turns, 0.0f);
  Mod::advance(state, config);
  HS_EXPECT_EQ(state.phase.phase_turns, .25f);
  config.sine.frequency = -1;
  HS_EXPECT_FALSE(Mod::valid(config));
  config.sine.frequency = 1;
  config.version = 2;
  HS_EXPECT_FALSE(Mod::valid(config));
  config.version = 1;
  config.source = static_cast<Pattern::SourceKind>(99);
  HS_EXPECT_FALSE(Mod::valid(config));
}
/** @brief Every source and mechanic compose through the operator's color binding. */
inline void test_shader_chain_modulation_combinations() {
  using Kind = ColorModulation::Mechanic;
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  HS_EXPECT_EQ(program.compile(DEFAULT_CHAIN).code, In::ChainStatus::OK);
  auto &params = param_as<In::Op::GeneratedPaletteParams>(program, 3);
  const PB::FieldSample SAMPLE{.37f, .6f, math::Vector(.2f, .7f, .4f), 0};
  const auto COORDINATE = Pattern::SphereDirection::from(SAMPLE.sphere);
  for (auto source :
       {Pattern::SourceKind::NOISE, Pattern::SourceKind::DIRECTIONAL_SINE,
        Pattern::SourceKind::SPHERE_INTERFERENCE,
        Pattern::SourceKind::SPHERE_SMOOTH_GRID})
    for (auto kind :
         {Kind::HUE_ROTATION, Kind::POSITIVE_LIGHTNESS_LIFT, Kind::OPACITY}) {
      params.modulation = {};
      params.modulation.source = source;
      params.modulation.evaluation = Pattern::Evaluation::DIRECT;
      params.modulation.mechanic = {kind, .7f};
      params.modulation.sine.offset = .13f;
      params.modulation.interference.offset = .13f;
      params.modulation.grid.offset = .13f;
      PB::Modulation::Resources<> resources;
      const auto clock = PB::Modulation::init(params.modulation);
      auto ctx = shared_resources().context();
      ctx.modulation = &resources.prepare(params.modulation, clock);
      ctx.hue_rotation_lut = nullptr;
      program.prepare(ctx);
      const Color4 ACTUAL = run_entry<Color4>(program, ctx, 3, SAMPLE);
      const float scalar = PB::Modulation::sample(*ctx.modulation, COORDINATE);
      Color4 expected = ctx.palettes[0]->get(SAMPLE.value);
      expected =
          ColorModulation::apply(expected, scalar, params.modulation.mechanic);
      expected.alpha *= SAMPLE.coverage;
      HS_EXPECT_EQ(ACTUAL.color, expected.color);
      HS_EXPECT_EQ(ACTUAL.alpha, expected.alpha);
      if (kind == Kind::OPACITY)
        HS_EXPECT_EQ(ACTUAL.color, ctx.palettes[0]->get(SAMPLE.value).color);
    }
}

inline int modulation_probe_calls = 0;
/** @brief Workbench columns blend completed alpha and hold current topology at zero. */
inline void test_shader_chain_modulation_temporal_endpoints() {
  using Kind = ColorModulation::Mechanic;
  for (int topology = 0; topology < 5; ++topology) {
    constexpr size_t BLOCKS =
        (In::CHAIN_ARENA_BYTES + sizeof(std::max_align_t) - 1) /
        sizeof(std::max_align_t);
    std::vector<std::max_align_t> history(BLOCKS), column(BLOCKS);
    auto fixture = std::make_unique<ProgramFixture>();
    auto &program = fixture->program;
    program.bind_temporal_storage(reinterpret_cast<uint8_t *>(history.data()),
                                  reinterpret_cast<uint8_t *>(column.data()));
    HS_EXPECT_EQ(program.compile(DEFAULT_CHAIN).code, In::ChainStatus::OK);
    auto &params = param_as<In::Op::GeneratedPaletteParams>(program, 3);
    params.modulation.evaluation = Pattern::Evaluation::DIRECT;
    params.modulation.mechanic = {Kind::OPACITY, .2f};
    PB::Modulation::Prepared previous;
    previous.config = params.modulation;
    previous.state.phase.phase_turns = .25f;
    previous.evaluate = +[](const PB::Modulation::Prepared &p,
                            const Pattern::SphereDirection &) {
      ++modulation_probe_calls;
      return 2 * p.state.phase.phase_turns - 1;
    };
    auto previous_ctx = shared_resources().context();
    previous_ctx.modulation = &previous;
    previous_ctx.hue_rotation_lut = nullptr;
    previous_ctx.palettes[0] = previous_ctx.palettes[1];
    program.prepare(previous_ctx);
    program.capture_history();
    params.modulation.mechanic.amount = .8f;
    switch (topology) {
    case 1:
      params.modulation.source = Pattern::SourceKind::DIRECTIONAL_SINE;
      break;
    case 2:
      params.modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
      break;
    case 3:
      params.modulation.mechanic.kind = Kind::HUE_ROTATION;
      break;
    case 4:
      params.modulation.enabled = false;
      break;
    }
    auto current = previous;
    current.config = params.modulation;
    current.state.phase.phase_turns = .75f;
    auto ctx = shared_resources().context();
    ctx.modulation = &current;
    ctx.hue_rotation_lut = nullptr;
    HS_EXPECT_TRUE(program.history_continuous());
    for (float fraction : {0.0f, .375f, 1.0f}) {
      program.prepare_column(ctx, previous_ctx, fraction);
      const auto &prepared =
          *reinterpret_cast<const PB::Color::GeneratedPaletteState *>(
              program.prepared_block(3));
      const float amount = topology == 0 ? hs::lerp(.2f, .8f, fraction) : .8f;
      HS_EXPECT_EQ(prepared.mechanic.amount, amount);
      const math::Vector VIEW = math::Vector(.3f, .7f, .2f).normalized();
      const auto sphere = run_entry<PB::SphereSample>(
          program, ctx, 0, PB::SphereSample{VIEW, 0});
      const auto plane = run_entry<PB::PlaneSample>(program, ctx, 1, sphere);
      const auto field = run_entry<PB::FieldSample>(program, ctx, 2, plane);
      modulation_probe_calls = 0;
      const Color4 actual =
          program.evaluate_column(VIEW, ctx, previous_ctx, fraction);
      if (topology == 4)
        HS_EXPECT_EQ(modulation_probe_calls, 0);
      Color4 before = previous_ctx.palettes[0]->get(field.value),
             after = ctx.palettes[0]->get(field.value);
      if (topology == 3) {
        before = hue_rotate_lut_gamut(before, .5f * amount);
        after = hue_rotate_lut_gamut(after, .5f * amount);
      } else if (topology != 4) {
        before.alpha *= hs::lerp(1.0f, topology == 0 ? .25f : .75f, amount);
        after.alpha *= hs::lerp(1.0f, .75f, amount);
      }
      before.alpha *= field.coverage;
      after.alpha *= field.coverage;
      const auto weight = frac_to_q16(fraction);
      HS_EXPECT_EQ(actual.color, before.color.lerp16(after.color, weight));
      HS_EXPECT_EQ(actual.alpha, hs::lerp(before.alpha, after.alpha, fraction));
    }
    if (topology == 0) {
      previous = current;
      previous_ctx.modulation = &previous;
      program.prepare(ctx);
      program.capture_history();
      program.prepare_column(ctx, previous_ctx, .375f);
      modulation_probe_calls = 0;
      (void)program.evaluate_column(math::Y_AXIS, ctx, previous_ctx, .375f);
      HS_EXPECT_EQ(modulation_probe_calls, 1);
    }
  }
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
        {"colorize", "colorize.generated-palette.v3"}};
    HS_EXPECT_EQ(effect.set_chain(topology).code, In::ChainStatus::OK);
    auto &program = WB::program(effect);
    param_as<In::Op::RotateChainParams>(program, 0).wander = 0.8f;
    param_as<In::Op::RotateChainParams>(program, 0).spin_rate = 0.001f;
    param_as<In::Op::GridSampleParams>(program, 4).speed = 0.02f;
    auto &color = WB::color_params(effect);
    color.modulation.mechanic.amount = 0.7f;
    color.modulation.noise.scale = 1.7f;
    color.palette_chroma = 0.43f;
    color.phase_oscillation_speed = 0.004f;
    color.modulation.noise.speed = -0.0007f;
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
  const auto live = effect.snapshot();
  // Candidates derive from a state that differs from the live one in every
  // section, so any section committed before a refusal is observable.
  param_as<In::Op::RotateChainParams>(WB::program(effect), 0).wander = 0.6f;
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
  const auto generation = effect.getParameterSchemaGeneration();
  const auto refused = [&](ChainSnapshot candidate,
                           ChainSnapshotRestoreResult expected) {
    HS_EXPECT_EQ(effect.restore_snapshot(candidate), expected);
    HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), generation);
    const auto after = effect.snapshot();
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
  auto candidate = saved;
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
    if (auto *clock = std::get_if<In::ColorClockSnapshot>(&runtime.state))
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
  auto &controls = WB::color_params(effect);
  controls.modulation.noise.speed = .125f;
  WB::advance_without_render(effect);
  const ChainSnapshot saved = effect.snapshot();
  const auto clock =
      [](const ChainSnapshot &snapshot) -> const In::ColorClockSnapshot & {
    for (const auto &runtime : *snapshot.runtime)
      if (const auto *value =
              std::get_if<In::ColorClockSnapshot>(&runtime.state))
        return *value;
    HS_CHECK(false, "missing test color clock");
    return std::get<In::ColorClockSnapshot>(snapshot.runtime->back().state);
  };
  HS_EXPECT_EQ(clock(saved).modulation.phase.phase_turns, .125f);
  auto mismatched = saved;
  for (auto &runtime : *mismatched.runtime)
    if (auto *value = std::get_if<In::ColorClockSnapshot>(&runtime.state))
      value->modulation.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  HS_EXPECT_EQ(effect.restore_snapshot(mismatched),
               ChainSnapshotRestoreResult::INVALID_VALUE);
  HS_EXPECT_EQ(effect.restore_snapshot(saved),
               ChainSnapshotRestoreResult::APPLIED);
  auto uninitialized = saved;
  for (auto &runtime : *uninitialized.runtime)
    if (auto *value = std::get_if<In::ColorClockSnapshot>(&runtime.state))
      value->modulation.initialized = false;
  HS_EXPECT_EQ(effect.restore_snapshot(uninitialized),
               ChainSnapshotRestoreResult::INVALID_VALUE);
  HS_EXPECT_EQ(effect.restore_snapshot(saved),
               ChainSnapshotRestoreResult::APPLIED);
  auto &changed = WB::color_params(effect).modulation;
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
  HS_RUN(test_shader_chain_status_names);
  HS_RUN(test_shader_chain_snapshot_roundtrip);
  HS_RUN(test_shader_chain_snapshot_refusals);
  HS_RUN(test_shader_chain_modulation_snapshot_identity);
  return fixture.result();
}

} // namespace shader_chain_tests
} // namespace hs_test
