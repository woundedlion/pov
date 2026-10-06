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

#include "math/mobius.h"
#include <cstdio>
#include <cstring>
#include <memory>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>

#include "core/render/pullback/catalog_export.h"
#include "core/render/pullback/interpreter.h"
#include "effects/AlienCore.h"
#include "tests/composed_frame_fixture.h"
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
/** @brief The hue tables are bound only on the modes that read them, and the
    hue-noise bake is skipped while both of its inputs hold still. */
inline void test_shader_chain_hue_lut_bake_cache() {
  using WB = ShaderChainWhiteBox;
  using HueShiftMode = In::Op::HueShiftMode;
  // Outside the quantizer's [-127, 127] range.
  constexpr int8_t POISON = -128;
  constexpr size_t LAST = PB::Color::HueNoiseLutView::SIZE - 1;
  reset_globals();
  WB::FX effect;
  effect.init();
  In::Op::GeneratedPaletteParams &color = WB::color_params(effect);
  int8_t *lut = WB::hue_noise_lut(effect);
  const auto poison = [lut] {
    lut[0] = POISON;
    lut[LAST] = POISON;
  };
  const auto expect_baked = [lut](bool baked) {
    HS_EXPECT_EQ(lut[0] != POISON, baked);
    HS_EXPECT_EQ(lut[LAST] != POISON, baked);
  };

  // Amount 0 disables both tables whatever the mode asks for.
  color.hue_mode = static_cast<uint8_t>(HueShiftMode::NOISE);
  color.hue_shift_amount = 0.0f;
  poison();
  In::FrameContext ctx = WB::frame_context(effect);
  HS_EXPECT_TRUE(ctx.hue_rotation_lut == nullptr);
  HS_EXPECT_TRUE(ctx.hue_noise_lut == nullptr);
  expect_baked(false);

  // NONE reads neither table, so a leftover amount buys no bake.
  color.hue_mode = static_cast<uint8_t>(HueShiftMode::NONE);
  color.hue_shift_amount = 0.5f;
  poison();
  ctx = WB::frame_context(effect);
  HS_EXPECT_TRUE(ctx.hue_rotation_lut == nullptr);
  HS_EXPECT_TRUE(ctx.hue_noise_lut == nullptr);
  expect_baked(false);

  // PATH_LENGTH reads the rotation table only.
  color.hue_mode = static_cast<uint8_t>(HueShiftMode::PATH_LENGTH);
  poison();
  ctx = WB::frame_context(effect);
  HS_EXPECT_TRUE(ctx.hue_rotation_lut != nullptr);
  HS_EXPECT_TRUE(ctx.hue_noise_lut == nullptr);
  expect_baked(false);

  // NOISE binds and bakes both.
  color.hue_mode = static_cast<uint8_t>(HueShiftMode::NOISE);
  color.hue_noise_scale = 1.5f;
  ctx = WB::frame_context(effect);
  HS_EXPECT_TRUE(ctx.hue_rotation_lut != nullptr);
  HS_EXPECT_TRUE(ctx.hue_noise_lut != nullptr);
  expect_baked(true);

  // Held inputs skip the bake.
  poison();
  ctx = WB::frame_context(effect);
  expect_baked(false);

  // A scale change re-bakes the whole table.
  color.hue_noise_scale = 3.0f;
  ctx = WB::frame_context(effect);
  expect_baked(true);

  // So does the loop phase moving under a held scale.
  poison();
  color.hue_noise_speed = 0.001f;
  const float held = WB::color_clocks(effect).hue_noise_phase;
  WB::program(effect).advance();
  HS_EXPECT_NE(WB::color_clocks(effect).hue_noise_phase, held);
  ctx = WB::frame_context(effect);
  expect_baked(true);
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
    color.hue_shift_amount = 0.7f;
    color.hue_noise_scale = 1.7f;
    color.palette_chroma = 0.43f;
    color.phase_oscillation_speed = 0.004f;
    color.hue_noise_speed = -0.0007f;
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
  reset_globals();
  ShaderChainWhiteBox::FX effect;
  effect.init();
  const auto saved = effect.snapshot();
  const auto generation = effect.getParameterSchemaGeneration();
  const auto refused = [&](ChainSnapshot candidate,
                           ChainSnapshotRestoreResult expected) {
    HS_EXPECT_EQ(effect.restore_snapshot(candidate), expected);
    HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), generation);
    const auto after = effect.snapshot();
    HS_EXPECT_EQ(after.parameters.size(), saved.parameters.size());
    for (size_t index = 0; index < saved.parameters.size(); ++index) {
      HS_EXPECT_TRUE(after.parameters[index].name ==
                     saved.parameters[index].name);
      HS_EXPECT_TRUE(float_identical(after.parameters[index].value,
                                     saved.parameters[index].value));
    }
    HS_EXPECT_EQ(after.palette_bank->cycles[0].frame,
                 saved.palette_bank->cycles[0].frame);
    HS_EXPECT_EQ(
        std::get<In::SpatialWalkSnapshot>((*after.runtime)[0].state).walk_time,
        std::get<In::SpatialWalkSnapshot>((*saved.runtime)[0].state).walk_time);
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
      clock->hue_noise_phase = 1.0f;
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

inline int run_shader_chain_tests() {
  ModuleFixture fixture("shader_chain");
  test_shader_chain_table_integrity();
  test_shader_chain_table_behavior();
  test_shader_chain_schema_and_field_ids();
  test_shader_chain_instance_id_wellformed();
  test_shader_chain_slot_and_hash_contract();
  test_shader_chain_catalog_golden();
  test_shader_chain_catalog_shape();
  test_shader_chain_default_chain_renders();
  test_shader_chain_param_address_channel();
  test_shader_chain_parity_rotate_project();
  test_shader_chain_parity_displace_curl();
  test_shader_chain_parity_displace_direct();
  test_shader_chain_parity_displace_ripple();
  test_shader_chain_parity_lens_ops();
  test_shader_chain_parity_project_ops();
  test_shader_chain_projection_frame_policy();
  test_shader_chain_parity_project_hemispheres();
  test_shader_chain_parity_field_ops();
  test_shader_chain_parity_warp_affine_mirror();
  test_shader_chain_parity_warp_wave_shear();
  test_shader_chain_parity_warp_vortex();
  test_shader_chain_parity_warp_vector_noise();
  test_shader_chain_noise_instances_decorrelate();
  test_shader_chain_parity_warp_polar_chart();
  test_shader_chain_parity_warp_curl_flow();
  test_shader_chain_parity_sample_variants();
  test_shader_chain_parity_sample_twin_wave();
  test_shader_chain_parity_sample_rings();
  test_shader_chain_parity_sample_spherical_rings();
  test_shader_chain_parity_sample_spiral();
  test_shader_chain_parity_sample_lattice();
  test_shader_chain_parity_sample_fractal();
  test_shader_chain_parity_sample_tessellation();
  test_shader_chain_large_finite_path_length();
  test_shader_chain_composed_noise_domain();
  test_shader_chain_noise_operator_domain();
  test_shader_chain_parity_sample_projected_noise();
  test_shader_chain_parity_sample_spherical_noise();
  test_shader_chain_parity_colorize_variants();
  test_shader_chain_refusal_shape();
  test_shader_chain_refusal_budget_overflows();
  test_shader_chain_program_lifetime();
  test_shader_chain_refusal_migrate_failed();
  test_shader_chain_state_identity_migration();
  test_shader_chain_state_continuity_slice();
  test_shader_chain_operator_state_migration();
  test_shader_chain_determinism();
  test_shader_chain_param_names_and_budget();
  test_shader_chain_composed_frame_parity();
  test_shader_chain_effect_registers_params();
  test_shader_chain_parameter_admission();
  test_shader_chain_edge_distance_admission();
  test_shader_chain_effect_rebind_generation();
  test_shader_chain_effect_refusal_keeps_schema();
  test_shader_chain_pause_semantics();
  test_shader_chain_hue_lut_bake_cache();
  test_shader_chain_status_names();
  test_shader_chain_snapshot_roundtrip();
  test_shader_chain_snapshot_refusals();
  return fixture.result();
}

} // namespace shader_chain_tests
} // namespace hs_test
