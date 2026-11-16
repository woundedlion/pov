/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "tests/shader_chain_fixture.h"
#include "tests/shader_chain/default_chain_mirror.h"

namespace hs_test::shader_chain_tests {

/** @brief Distribution controls produce monotone lighting probabilities. */
inline void test_shader_chain_probability_response() {
  using Op = In::Op::CoverageProbability;
  Op::Params params;
  const In::FrameContext ctx;
  Op::State state;
  auto prepared = Op::prepare(ctx, params, state);
  HS_EXPECT_NEAR(Op::probability(.25f, params, prepared), .25f, 1e-6f);
  params.field_min = .2f;
  params.field_max = .8f;
  params.probability_min = .1f;
  params.probability_max = .9f;
  HS_EXPECT_NEAR(Op::probability(.35f, params, prepared), .3f, 1e-6f);
  HS_EXPECT_NEAR(Op::probability(0, params, prepared), .1f, 1e-6f);
  HS_EXPECT_NEAR(Op::probability(1, params, prepared), .9f, 1e-6f);
  std::swap(params.field_min, params.field_max);
  std::swap(params.probability_min, params.probability_max);
  HS_EXPECT_NEAR(Op::probability(.35f, params, prepared), .3f, 1e-6f);
  params.field_min = params.field_max = .5f;
  HS_EXPECT_NEAR(Op::probability(.49f, params, prepared), .1f, 1e-6f);
  HS_EXPECT_NEAR(Op::probability(.5f, params, prepared), .9f, 1e-6f);

  params = {};
  params.distribution = 1;
  params.rate = 2;
  prepared = Op::prepare(ctx, params, state);
  HS_EXPECT_NEAR(Op::probability(.5f, params, prepared), .632120559f, 1e-6f);
  params.rate = 0;
  HS_EXPECT_EQ(Op::probability(1, params, prepared), 0);
  params.distribution = 2;
  params.deviation = .1f;
  prepared = Op::prepare(ctx, params, state);
  HS_EXPECT_NEAR(Op::probability(.3f, params, prepared), .022749858f, 1e-6f);
  HS_EXPECT_NEAR(Op::probability(.5f, params, prepared), .5f, 1e-6f);
  params.mean = .3f;
  prepared = Op::prepare(ctx, params, state);
  HS_EXPECT_GT(Op::probability(.3f, params, prepared), .49f);

  for (uint8_t distribution = 0; distribution < 3; ++distribution)
    for (float deviation : {.01f, .2f, 1.0f})
      for (float mean : {0.0f, .5f, 1.0f}) {
        params.distribution = distribution;
        params.deviation = deviation;
        params.mean = mean;
        params.rate = 16;
        prepared = Op::prepare(ctx, params, state);
        float previous = 0;
        for (int step = 0; step <= 100; ++step) {
          const float probability =
              Op::probability(step * .01f, params, prepared);
          HS_EXPECT_TRUE(std::isfinite(probability));
          HS_EXPECT_GE(probability, previous);
          HS_EXPECT_LE(probability, 1);
          previous = probability;
        }
      }
}

/** @brief Masks match event rates, reroll per frame and survive migration. */
inline void test_shader_chain_probability_mask() {
  using Op = In::Op::CoverageProbability;
  const In::FrameContext ctx;
  const In::InstanceId ID{"sparkles", Op::ID,
                          In::instance_hash("sparkles", Op::ID)};
  Op::State state;
  Op::init(state, ID);
  for (uint8_t distribution = 0; distribution < 3; ++distribution) {
    Op::Params params;
    params.distribution = distribution;
    params.rate = 2;
    const float EXPECTED = distribution == 1 ? .632120559f : .5f;
    const auto prepared = Op::prepare(ctx, params, state);
    Op::advance(state, params);
    const auto next = Op::prepare(ctx, params, state);
    auto seeded_params = params;
    seeded_params.seed = 1;
    const auto seeded = Op::prepare(ctx, seeded_params, state);
    size_t lit = 0;
    size_t changed = 0;
    size_t stream_changed = 0;
    constexpr size_t POINTS = 8192;
    for (size_t index = 0; index < POINTS; ++index) {
      const float ANGLE = math::TWO_PI_F * float(index) / float(POINTS);
      const PB::FieldSample INPUT{.5f, .6f, {cosf(ANGLE), sinf(ANGLE), 0}, 3};
      const auto output = Op::run(INPUT, ctx, params, prepared);
      HS_EXPECT_TRUE(output.coverage == 0 || output.coverage == INPUT.coverage);
      HS_EXPECT_EQ(output.value, INPUT.value);
      HS_EXPECT_EQ(output.path_length, INPUT.path_length);
      HS_EXPECT_EQ(output.sphere.x, INPUT.sphere.x);
      HS_EXPECT_EQ(output.sphere.y, INPUT.sphere.y);
      HS_EXPECT_EQ(output.sphere.z, INPUT.sphere.z);
      HS_EXPECT_TRUE(
          field_identical(output, Op::run(INPUT, ctx, params, prepared)));
      lit += output.coverage > 0;
      const auto following = Op::run(INPUT, ctx, params, next);
      changed += output.coverage != following.coverage;
      stream_changed += following.coverage !=
                        Op::run(INPUT, ctx, seeded_params, seeded).coverage;
      auto blocked = INPUT;
      blocked.coverage = 0;
      HS_EXPECT_EQ(Op::run(blocked, ctx, params, prepared).coverage, 0);
    }
    HS_EXPECT_NEAR(float(lit) / float(POINTS), EXPECTED, .025f);
    HS_EXPECT_GT(changed, POINTS / 4);
    HS_EXPECT_GT(stream_changed, POINTS / 4);
    Op::State migrated;
    HS_EXPECT_EQ(Op::migrate(migrated, state, ID), In::Status::OK);
    HS_EXPECT_EQ(migrated.frame, state.frame);
    HS_EXPECT_EQ(migrated.seed, state.seed);
    const auto saved = In::RuntimeStateCodec<Op::State>::capture(state);
    HS_EXPECT_TRUE(std::holds_alternative<In::ProbabilityClockSnapshot>(saved));
    Op::advance(state, params);
    HS_EXPECT_TRUE(In::RuntimeStateCodec<Op::State>::restore(state, saved));
    HS_EXPECT_EQ(state.frame, migrated.frame);
    HS_EXPECT_EQ(state.seed, migrated.seed);
    const auto interval = In::Op::temporal_interval(Op::State{}, state);
    for (float fraction : {0.0f, .5f, 1.0f}) {
      auto column = state;
      In::Op::temporal_sample(column, interval, params, fraction);
      HS_EXPECT_EQ(column.frame, state.frame);
      HS_EXPECT_EQ(column.seed, state.seed);
    }
  }
  Op::Params params;
  const auto prepared = Op::prepare(ctx, params, state);
  HS_EXPECT_EQ(Op::run({0, 1, math::X_AXIS, 0}, ctx, params, prepared).coverage,
               0);
  HS_EXPECT_EQ(Op::run({1, 1, math::X_AXIS, 0}, ctx, params, prepared).coverage,
               1);
  state.frame = UINT32_MAX;
  Op::advance(state, params);
  HS_EXPECT_EQ(state.frame, 0u);
}

// --- field-op mirrors -----------------------------------------------------

/** Mirror frame for the FIELD endomorphism batch. */
struct FieldMirrorFrame {
  In::Op::IsoContourChainParams iso;
  In::Op::SmoothBandsChainParams bands;
  In::Op::ValueCutoutChainParams cutout;
};

struct FieldMirrorBinding {
  using FrameState = FieldMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct IsoContourFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float iso_level(const FieldMirrorFrame &frame) {
    return frame.iso.iso_level;
  }
  static float iso_width(const FieldMirrorFrame &frame) {
    return frame.iso.iso_width;
  }
};

struct SmoothBandsFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float band_count(const FieldMirrorFrame &frame) {
    return frame.bands.band_count;
  }
  static float band_phase(const FieldMirrorFrame &frame) {
    return frame.bands.band_phase;
  }
};

struct ValueCutoutFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float cutout_threshold(const FieldMirrorFrame &frame) {
    return frame.cutout.cutout_threshold;
  }
  static float cutout_softness(const FieldMirrorFrame &frame) {
    return frame.cutout.cutout_softness;
  }
};

/** Compiles a chain with one FIELD endomorphism at entry 3 and applies the
    value set to every block. */
template <typename OpParams>
inline bool arm_field_op_chain(In::ChainProgram &program, const char *op_id,
                               int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v4"},
      {"op", op_id},
      {"colorize", "colorize.generated-palette.v4"},
  };
  return arm_chain<In::Op::RotateChainParams, In::Op::ProjectChainParams,
                   In::Op::GridSampleParams, OpParams,
                   In::Op::GeneratedPaletteParams>(program, chain, frames, set);
}

/** The FIELD op's input carrier for @p view: the erased camera, projection
    and sample entries run over the seed. */
inline PB::FieldSample field_input(In::ChainProgram &program,
                                   const In::FrameContext &ctx,
                                   const math::Vector &view) {
  return run_entry<PB::FieldSample>(program, ctx, 2,
                                    projected_input(program, ctx, view));
}

inline FieldMirrorFrame field_mirror(In::ChainProgram &program) {
  FieldMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[3].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::TransferIsoContour::ID)
    mirror.iso = param_as<In::Op::IsoContourChainParams>(program, 3);
  else if (id == In::Op::TransferSmoothBands::ID)
    mirror.bands = param_as<In::Op::SmoothBandsChainParams>(program, 3);
  else if (id == In::Op::CoverageValueCutout::ID)
    mirror.cutout = param_as<In::Op::ValueCutoutChainParams>(program, 3);
  return mirror;
}

/** Erased-vs-bound parity of the FIELD endomorphism at entry 3. */
template <typename BoundStage>
inline void expect_field_op_parity(In::ChainProgram &program,
                                   const In::FrameContext &ctx) {
  program.prepare(ctx);
  expect_op_parity<BoundStage>(
      program, ctx, 3, field_mirror(program),
      [&](const math::Vector &view) { return field_input(program, ctx, view); },
      field_identical);
}

template <typename BoundStage, typename OpParams>
inline void run_field_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  if (!arm_field_op_chain<OpParams>(program, op_id, 4, set))
    return;
  const In::FrameContext ctx = shared_resources().context();
  expect_field_op_parity<BoundStage>(program, ctx);
  program.clear();
}

inline void test_shader_chain_parity_field_ops() {
  for_each_value_set([&](ValueSet set) {
    run_field_parity<typename PB::Stage::Transfer<PB::Transfer::Ridge>::
                         template Bind<FieldMirrorBinding>,
                     PB::Transfer::NoValueParams>(In::Op::TransferRidge::ID,
                                                  set);
    run_field_parity<
        typename PB::Stage::Transfer<PB::Transfer::IsoContour<
            IsoContourFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::IsoContourChainParams>(In::Op::TransferIsoContour::ID, set);
    run_field_parity<
        typename PB::Stage::Transfer<PB::Transfer::SmoothBands<
            SmoothBandsFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::SmoothBandsChainParams>(In::Op::TransferSmoothBands::ID, set);
    run_field_parity<
        typename PB::Stage::ApplyCoverage<PB::ValueCoverage::ValueCutout<
            ValueCutoutFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::ValueCutoutChainParams>(In::Op::CoverageValueCutout::ID, set);
  });
}

} // namespace hs_test::shader_chain_tests
