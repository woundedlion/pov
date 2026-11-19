/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "tests/shader_chain_fixture.h"
#include "tests/shader_chain/default_chain_mirror.h"

namespace hs_test::shader_chain_tests {

// --- sample-op mirrors ----------------------------------------------------

/** Mirror frame for the sample source batch. */
struct SampleMirrorFrame {
  const FastNoiseLite *noise = nullptr;
  In::Op::TwinWaveSampleParams twin_wave;
  In::Op::BandsSampleParams bands;
  In::Op::SpiralSampleParams spiral;
  In::Op::TessellationSampleParams tessellation;
  In::Op::ProjectedNoiseSampleParams projected;
  float primary = 0.0f;
  float secondary = 0.0f;
  float angle = 0.0f;
  float noise_time = 0.0f;
  float edge_width = 0.1f;
};

struct SampleMirrorBinding {
  using FrameState = SampleMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct ClockedSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static Pattern::PreparedSource prepare(const FrameState &frame) {
    return Pattern::prepare_phases(frame.primary, frame.secondary, frame.angle);
  }
};

struct TwinWaveSampleMirror : ClockedSampleMirror {
  static const In::Op::TwinWaveSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.twin_wave;
  }
};

struct BandsSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const In::Op::BandsSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.bands;
  }
};

struct SpiralSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const In::Op::SpiralSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.spiral;
  }
};

struct TessellationSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const In::Op::TessellationSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.tessellation;
  }
};

struct ProjectedNoiseSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const FastNoiseLite &noise(const SampleMirrorFrame &frame) {
    return *frame.noise;
  }
  static float noise_scale(const SampleMirrorFrame &frame) {
    return frame.projected.noise_scale;
  }
  static float noise_time(const SampleMirrorFrame &frame) {
    return frame.noise_time;
  }
  static float noise_contrast(const SampleMirrorFrame &frame) {
    return frame.projected.noise_contrast;
  }
};

/** Edge-width provider for the mirror's edge-fade coverage policy. */
struct SampleEdgeMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static float edge_width(const SampleMirrorFrame &frame) {
    return frame.edge_width;
  }
};

/** Compiles a chain with one sample crossing at entry 2 and applies the value
    set to every block. */
template <typename OpParams>
inline bool arm_sample_op_chain(In::ChainProgram &program, const char *op_id,
                                int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"source", op_id},
      {"colorize", "colorize.generated-palette.v4"},
  };
  return arm_chain<In::Op::RotateChainParams, In::Op::ProjectChainParams,
                   OpParams, In::Op::GeneratedPaletteParams>(program, chain,
                                                             frames, set);
}

inline SampleMirrorFrame sample_mirror(In::ChainProgram &program,
                                       size_t source_index = 2) {
  SampleMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[source_index].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::SampleProjectedNoise::ID) {
    const auto &state =
        state_as<In::Op::NoisePhaseState>(program, source_index);
    mirror.noise = &state.noise;
    mirror.noise_time = state.phase;
    mirror.projected =
        param_as<In::Op::ProjectedNoiseSampleParams>(program, source_index);
    mirror.edge_width = mirror.projected.edge_width;
  } else if (id == In::Op::SampleTessellation::ID) {
    mirror.tessellation =
        param_as<In::Op::TessellationSampleParams>(program, source_index);
    mirror.edge_width = mirror.tessellation.edge_width;
  } else if (id == In::Op::SampleBands::ID) {
    mirror.bands = param_as<In::Op::BandsSampleParams>(program, source_index);
    mirror.edge_width = mirror.bands.edge_width;
  } else if (id == In::Op::SampleSpiralV3::ID) {
    mirror.spiral = param_as<In::Op::SpiralSampleParams>(program, source_index);
    mirror.edge_width = mirror.spiral.edge_width;
  } else {
    const auto &state =
        state_as<In::Op::SourceClockState>(program, source_index);
    mirror.primary = state.primary;
    mirror.secondary = state.secondary;
    mirror.angle = state.angle;
    if (id == In::Op::SampleTwinWaveV3::ID) {
      mirror.twin_wave =
          param_as<In::Op::TwinWaveSampleParams>(program, source_index);
      mirror.edge_width = mirror.twin_wave.edge_width;
    }
  }
  return mirror;
}

/** Erased-vs-bound parity of the sample crossing at entry 2. */
template <typename BoundStage>
inline void expect_sample_op_parity(In::ChainProgram &program,
                                    const In::FrameContext &ctx) {
  program.prepare(ctx);
  expect_op_parity<BoundStage>(
      program, ctx, 2, sample_mirror(program),
      [&](const math::Vector &view) {
        return projected_input(program, ctx, view);
      },
      field_identical);
}

/** Address of a topology enum8 inside entry @p index's param block, resolved
    through the erased param-address channel. */
inline uint8_t *topology_byte(In::ChainProgram &program, size_t index,
                              const char *field_id) {
  const In::OperatorDescriptor &op = *program.ops()[index].op;
  for (uint16_t field = 0; field < op.schema_count; ++field)
    if (std::string_view(op.schema[field].id) == field_id)
      return static_cast<uint8_t *>(
          op.runtime.param_address(program.param_block(index), field));
  return nullptr;
}

/** @brief Dispatches a coverage mode to its typed policy. */
template <typename EdgeP, typename Fn>
inline void with_coverage_policy(In::Op::ProjectionCoverageMode coverage,
                                 Fn &&fn) {
  switch (coverage) {
  case In::Op::ProjectionCoverageMode::NONE:
    fn.template operator()<PB::ProjectionCoverage::None>();
    break;
  case In::Op::ProjectionCoverageMode::WEIGHT:
    fn.template operator()<PB::ProjectionCoverage::Weight>();
    break;
  case In::Op::ProjectionCoverageMode::WEIGHT_SQUARED:
    fn.template operator()<PB::ProjectionCoverage::WeightSquared>();
    break;
  case In::Op::ProjectionCoverageMode::EDGE_FADE:
    fn.template operator()<PB::ProjectionCoverage::EdgeFade<EdgeP>>();
    break;
  case In::Op::ProjectionCoverageMode::EDGE_FADE_SQUARED:
    fn.template operator()<PB::ProjectionCoverage::EdgeFadeSquared<EdgeP>>();
    break;
  }
}

template <typename SourceP, typename WeightP>
inline void run_sample_coverage_cases(In::ChainProgram &program,
                                      const In::FrameContext &ctx,
                                      uint8_t coverage) {
  with_coverage_policy<SampleEdgeMirror>(
      static_cast<In::Op::ProjectionCoverageMode>(coverage),
      [&]<typename Policy>() {
        expect_sample_op_parity<typename PB::Stage::Sample<
            SourceP, WeightP, Policy>::template Bind<SampleMirrorBinding>>(
            program, ctx);
      });
}

/** Runs the full weight-by-coverage topology matrix of the sample crossing at
    entry 2 against @p SourceP. */
template <typename SourceP>
inline void run_sample_op_matrix(In::ChainProgram &program,
                                 const In::FrameContext &ctx) {
  uint8_t *weight = topology_byte(program, 2, "weight-mode");
  uint8_t *coverage = topology_byte(program, 2, "coverage-mode");
  HS_EXPECT_TRUE(weight != nullptr && coverage != nullptr);
  if (!weight || !coverage)
    return;
  for (uint8_t w = 0; w < std::size(In::Op::WEIGHT_MODE_IDS); ++w)
    for (uint8_t c = 0; c < std::size(In::Op::COVERAGE_MODE_IDS); ++c) {
      HS_CONTEXT("weight / coverage", w, c);
      *weight = w;
      *coverage = c;
      if (static_cast<In::Op::WeightMode>(w) == In::Op::WeightMode::NONE)
        run_sample_coverage_cases<SourceP, PB::Weight::None>(program, ctx, c);
      else
        run_sample_coverage_cases<SourceP, PB::Weight::Projection>(program, ctx,
                                                                   c);
    }
}

inline void test_shader_chain_parity_sample_twin_wave() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::TwinWaveSampleParams>(
            program, In::Op::SampleTwinWaveV3::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::TwinWave<TwinWaveSampleMirror>>(program,
                                                                     ctx);
  });
}

/** @brief Band source parity across signal and coverage selectors. */
inline void test_shader_chain_parity_sample_bands() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::BandsSampleParams>(
            program, In::Op::SampleBands::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Bands<BandsSampleMirror>>(program, ctx);
  });
}

/** @brief Bands column angles match direct midpoint evaluation across views. */
inline void test_shader_chain_bands_column_angle() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  if (!arm_sample_op_chain<In::Op::BandsSampleParams>(
          program, In::Op::SampleBands::ID, 4, ValueSet::DEFAULTS))
    return;
  const In::FrameContext ctx = shared_resources().context();
  auto blocks = std::make_unique<TemporalBlocks>();
  program.bind_temporal_storage(blocks->history, blocks->column);
  auto &params = param_as<In::Op::BandsSampleParams>(program, 2);
  for (const auto angles : {std::pair{5.5f, 5.6f}, std::pair{3.0f, 3.3f}}) {
    params.angle = angles.first;
    HS_EXPECT_TRUE(PB::Fields::valid(params));
    program.prepare(ctx);
    program.capture_history();
    params.angle = angles.second;
    program.prepare(ctx);
    program.prepare_column(ctx, ctx, 0.5f);
    const auto views = sweep_views();
    std::vector<Color4> actual;
    for (const auto &view : views)
      actual.push_back(program.evaluate_column(view, ctx, ctx, 0.5f));
    params.angle = (angles.first + angles.second) * 0.5f;
    program.prepare(ctx);
    for (size_t index = 0; index < views.size(); ++index) {
      const auto expected = program.evaluate(views[index], ctx);
      HS_EXPECT_EQ(actual[index].color, expected.color);
      HS_EXPECT_NEAR(actual[index].alpha, expected.alpha, 1e-6f);
    }
  }
  program.clear();
}

/** @brief Equirectangular bands at angle zero depend only on latitude. */
inline void test_shader_chain_bands_latitudes() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  const In::ChainEntryRequest CHAIN[] = {
      {"project", "project.equirectangular.v2"},
      {"source", "sample.bands.v1"},
      {"colorize", "colorize.generated-palette.v4"},
  };
  if (!arm_chain<In::Op::MeridianProjectChainParams, In::Op::BandsSampleParams,
                 In::Op::GeneratedPaletteParams>(program, CHAIN, 0,
                                                 ValueSet::DEFAULTS))
    return;
  auto &bands = param_as<In::Op::BandsSampleParams>(program, 1);
  bands.pattern_freq = 4.0f;
  bands.weight_mode = static_cast<uint8_t>(In::Op::WeightMode::NONE);
  bands.coverage_mode =
      static_cast<uint8_t>(In::Op::ProjectionCoverageMode::NONE);
  const auto ctx = shared_resources().context();
  program.prepare(ctx);
  for (const float latitude : {-0.9f, -0.3f, 0.0f, 0.4f, 1.1f}) {
    const float expected = 0.5f + 0.5f * std::sin(4.0f * latitude);
    for (const float longitude : {-2.7f, -1.2f, 0.0f, 0.8f, 2.9f}) {
      const math::Vector view{std::cos(latitude) * std::cos(longitude),
                              std::sin(latitude),
                              std::cos(latitude) * std::sin(longitude)};
      const auto plane = run_entry<PB::PlaneSample>(
          program, ctx, 0, PB::SphereSample{view, 0.0f});
      const auto actual = run_entry<PB::FieldSample>(program, ctx, 1, plane);
      HS_EXPECT_NEAR(actual.value, expected, 0.004f);
      HS_EXPECT_EQ(actual.coverage, 1.0f);
    }
  }
}

inline void test_shader_chain_parity_sample_spiral() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::SpiralSampleParams>(
            program, In::Op::SampleSpiralV3::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Spiral<SpiralSampleMirror>>(program, ctx);
  });
}

/** @brief Equirectangular spiral controls move a spherical winding center. */
inline void test_shader_chain_spiral_spherical_center() {
  const auto ctx = shared_resources().context();
  In::Op::SpiralSampleParams params;
  params.coordinate_system = 1;
  params.pattern_freq = 2.0f;
  params.arms = 3.0f;
  params.weight_mode = static_cast<uint8_t>(In::Op::WeightMode::NONE);
  params.coverage_mode =
      static_cast<uint8_t>(In::Op::ProjectionCoverageMode::NONE);
  PB::PlaneSample input{};
  input.coords = {0.4f, 0.7f};
  auto prepared = In::Op::SampleSpiralV3::prepare(ctx, params, {});
  const float expected =
      std::sin(2.0f * (math::PI_F / 2.0f - 0.7f) - params.arms * 0.4f);
  HS_EXPECT_NEAR(
      In::Op::SampleSpiralV3::run(input, ctx, params, prepared).value,
      0.5f + 0.5f * expected, 0.008f);
  params.arms = 0.0f;
  for (const float azimuth : {0.0f, 0.9f, 4.1f})
    for (const float polar : {0.0f, 0.8f, 2.2f}) {
      params.center_azimuth = azimuth;
      params.center_polar = polar;
      prepared = In::Op::SampleSpiralV3::prepare(ctx, params, {});
      const float dot =
          std::cos(0.7f) * std::sin(polar) * std::cos(0.4f - azimuth) +
          std::sin(0.7f) * std::cos(polar);
      const float reference = 0.5f + 0.5f * std::sin(2.0f * std::acos(dot));
      HS_EXPECT_NEAR(
          In::Op::SampleSpiralV3::run(input, ctx, params, prepared).value,
          reference, 0.002f);
    }
  params.center_azimuth = 0.0f;
  params.center_polar = math::PI_F / 2.0f;
  prepared = In::Op::SampleSpiralV3::prepare(ctx, params, {});
  input.coords = {};
  HS_EXPECT_NEAR(
      In::Op::SampleSpiralV3::run(input, ctx, params, prepared).value, 0.5f,
      0.001f);
  input.coords = {0.3f, 0.0f};
  HS_EXPECT_NEAR(
      In::Op::SampleSpiralV3::run(input, ctx, params, prepared).value,
      0.5f + 0.5f * std::sin(0.6f), 0.002f);
}

/** @brief Spherical spirals join both longitude boundaries for moved centers. */
inline void test_shader_chain_spiral_longitude_wrap() {
  const auto ctx = shared_resources().context();
  In::Op::SpiralSampleParams params;
  params.coordinate_system = 1;
  params.weight_mode = static_cast<uint8_t>(In::Op::WeightMode::NONE);
  params.coverage_mode =
      static_cast<uint8_t>(In::Op::ProjectionCoverageMode::NONE);
  for (const float arms : {0.0f, 1.0f, 3.0f, 12.0f})
    for (const float frequency : {0.1f, 2.0f, 20.0f})
      for (const float polar : {0.0f, 0.8f, math::PI_F}) {
        params.arms = arms;
        params.pattern_freq = frequency;
        params.center_azimuth = 0.9f;
        params.center_polar = polar;
        const auto prepared = In::Op::SampleSpiralV3::prepare(ctx, params, {});
        for (const float latitude :
             {-math::PI_F / 2.0f, -0.7f, 0.0f, 0.9f, math::PI_F / 2.0f}) {
          PB::PlaneSample left{}, right{};
          left.coords = {0.0f, latitude};
          right.coords = {math::TWO_PI_F, latitude};
          HS_EXPECT_EQ(
              In::Op::SampleSpiralV3::run(left, ctx, params, prepared).value,
              In::Op::SampleSpiralV3::run(right, ctx, params, prepared).value);
          left.coords.re = -math::PI_F;
          right.coords.re = math::PI_F;
          HS_EXPECT_EQ(
              In::Op::SampleSpiralV3::run(left, ctx, params, prepared).value,
              In::Op::SampleSpiralV3::run(right, ctx, params, prepared).value);
        }
      }
}

/** @brief Noise and speed operators reduce restored turn clocks before adding. */
inline void test_shader_chain_turn_clock_reduction() {
  for (const float phase : {-16777216.0f, 16777216.0f}) {
    In::Op::NoisePhaseState noise_state;
    noise_state.phase = phase;
    In::Op::ProjectedNoiseSampleParams noise_params;
    noise_params.noise_time_rate = 0.25f;
    In::Op::SampleProjectedNoise::advance(noise_state, noise_params);
    HS_EXPECT_EQ(noise_state.phase, 0.25f);
    In::Op::WarpPhaseState speed_state;
    speed_state.phase = phase;
    In::Op::WaveShearWarpParams speed_params;
    speed_params.speed = -0.25f;
    In::Op::WarpWaveShear::advance(speed_state, speed_params);
    HS_EXPECT_EQ(speed_state.phase, 0.75f);
  }
}

/**
 * @brief Noise operators return finite values for extreme planar coordinates.
 */
inline void test_shader_chain_noise_operator_domain() {
  const auto ctx = shared_resources().context();
  In::Op::NoisePhaseState state;
  In::Op::init_noise_phase(state, {"noise", "sample.projected-noise.v2", 1337});
  for (const auto basis : {math::NoiseBasis::SIMPLEX, math::NoiseBasis::FBM3,
                           math::NoiseBasis::RIDGED3}) {
    In::Op::ProjectedNoiseSampleParams sample_params;
    sample_params.noise_scale = 64.0f;
    sample_params.basis = static_cast<uint8_t>(basis);
    const auto sample_prepared =
        In::Op::SampleProjectedNoise::prepare(ctx, sample_params, state);
    In::Op::VectorNoiseWarpParams vector_params;
    vector_params.scale = 64.0f;
    vector_params.strength = 30.0f;
    vector_params.basis = static_cast<uint8_t>(basis);
    const auto vector_prepared =
        In::Op::WarpVectorNoise::prepare(ctx, vector_params, state);
    In::Op::CurlFlowWarpParams curl_params;
    curl_params.scale = 4.0f;
    curl_params.strength = 1.0f / 32.0f;
    curl_params.basis = static_cast<uint8_t>(basis);
    for (int signs = 0; signs < 4; ++signs) {
      PB::PlaneSample input{};
      input.coords = {(signs & 1) ? 0x1p20f : -0x1p20f,
                      (signs & 2) ? 0x1p20f : -0x1p20f};
      const auto sample = In::Op::SampleProjectedNoise::run(
          input, ctx, sample_params, sample_prepared);
      const auto coordinate = math::noise_projected_coordinate(
          input.coords, sample_params.noise_scale, sample_prepared.loop_offset);
      const float raw = Pattern::noise_contour(state.noise, basis, coordinate,
                                               sample_params.noise_contrast);
      const auto expected =
          In::Op::finish_sample(input, raw, sample_params, ctx);
      HS_EXPECT_EQ(sample.value, expected.value);
      HS_EXPECT_TRUE(std::isfinite(sample.value));
      const auto vector = In::Op::WarpVectorNoise::run(
          input, ctx, vector_params, vector_prepared);
      HS_EXPECT_TRUE(std::isfinite(vector.coords.re));
      HS_EXPECT_TRUE(std::isfinite(vector.coords.im));
      for (uint8_t integrator = 0; integrator < 3; ++integrator) {
        curl_params.integrator = integrator;
        const auto prepared =
            In::Op::WarpCurlFlow::prepare(ctx, curl_params, state);
        const auto curl =
            In::Op::WarpCurlFlow::run(input, ctx, curl_params, prepared);
        HS_EXPECT_TRUE(std::isfinite(curl.coords.re));
        HS_EXPECT_TRUE(std::isfinite(curl.coords.im));
      }
    }
    for (const float coordinate :
         {1e11f, -1e11f, std::numeric_limits<float>::infinity(),
          -std::numeric_limits<float>::infinity(),
          std::numeric_limits<float>::quiet_NaN()}) {
      for (int axis = 0; axis < 2; ++axis) {
        PB::PlaneSample input{};
        input.coords = {axis == 0 ? coordinate : 1.0f,
                        axis == 1 ? coordinate : 2.0f};
        const auto sample = In::Op::SampleProjectedNoise::run(
            input, ctx, sample_params, sample_prepared);
        HS_EXPECT_EQ(sample.value, 0.5f);
        const auto vector = In::Op::WarpVectorNoise::run(
            input, ctx, vector_params, vector_prepared);
        HS_EXPECT_TRUE(float_identical(vector.coords.re, input.coords.re));
        HS_EXPECT_TRUE(float_identical(vector.coords.im, input.coords.im));
        HS_EXPECT_EQ(vector.path_length, input.path_length);
        for (uint8_t integrator = 0; integrator < 3; ++integrator) {
          curl_params.integrator = integrator;
          const auto prepared =
              In::Op::WarpCurlFlow::prepare(ctx, curl_params, state);
          const auto curl =
              In::Op::WarpCurlFlow::run(input, ctx, curl_params, prepared);
          HS_EXPECT_TRUE(float_identical(curl.coords.re, input.coords.re));
          HS_EXPECT_TRUE(float_identical(curl.coords.im, input.coords.im));
          HS_EXPECT_EQ(curl.path_length, input.path_length);
        }
      }
    }
  }
}

/**
 * @brief Out-of-domain composed noise samples resolve to the neutral field
 * value 0.5.
 */
inline void test_shader_chain_composed_noise_domain() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest CHAIN[] = {
      {"project", In::Op::ProjectStereographic::ID},
      {"affine-a", In::Op::WarpAffineV4::ID},
      {"affine-b", In::Op::WarpAffineV4::ID},
      {"affine-c", In::Op::WarpAffineV4::ID},
      {"affine-d", In::Op::WarpAffineV4::ID},
      {"affine-e", In::Op::WarpAffineV4::ID},
      {"affine-f", In::Op::WarpAffineV4::ID},
      {"sample", In::Op::SampleProjectedNoise::ID},
      {"color", In::Op::ColorizeGeneratedPaletteV4::ID},
  };
  if (!expect_compiled(program.compile(CHAIN)))
    return;
  auto &projection = param_as<In::Op::ProjectChainParams>(program, 0);
  projection.frame = static_cast<uint8_t>(In::Op::ProjectionFrame::IDENTITY);
  for (int index = 1; index <= 6; ++index) {
    auto &params = param_as<In::Op::AffineWarpParams>(program, index);
    params.scale_x = params.scale_y = 1.0f / 64.0f;
    HS_EXPECT_TRUE(PB::Fields::valid(params));
  }
  const auto ctx = shared_resources().context();
  const PB::SphereSample sphere{math::Vector(1, 1, 0).normalized(), 0.0f};
  for (const auto basis : {math::NoiseBasis::SIMPLEX, math::NoiseBasis::FBM3,
                           math::NoiseBasis::RIDGED3}) {
    param_as<In::Op::ProjectedNoiseSampleParams>(program, 7).basis =
        static_cast<uint8_t>(basis);
    program.prepare(ctx);
    const PB::PlaneSample plane = run_to_plane(program, ctx, sphere, 7);
    HS_EXPECT_GT(plane.coords.re, 1e11f);
    const auto field = run_entry<PB::FieldSample>(program, ctx, 7, plane);
    HS_EXPECT_EQ(field.value, 0.5f);
    const Color4 actual = program.evaluate(sphere.dir, ctx);
    HS_EXPECT_TRUE(std::isfinite(actual.alpha));
  }
}

/**
 * @brief Lattice reduction keeps path length finite at large finite
 * coordinates.
 */
inline void test_shader_chain_large_finite_path_length() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  std::array<In::ChainEntryRequest, 15> chain{};
  std::array<std::string, 11> ids;
  chain[0] = {"project", In::Op::ProjectStereographic::ID};
  for (int index = 0; index < 11; ++index) {
    ids[index] = "warp" + std::to_string(index);
    chain[index + 1] = {ids[index], In::Op::WarpAffineV4::ID};
  }
  chain[12] = {"sample", In::Op::SampleGridV4::ID};
  chain[13] = {"color", In::Op::ColorizeGeneratedPaletteV4::ID};
  chain[14] = {"displacement", In::Op::ModulateColorV1::ID};
  if (!expect_compiled(program.compile(chain)))
    return;
  for (int index = 1; index <= 11; ++index) {
    auto &params = param_as<In::Op::AffineWarpParams>(program, index);
    params.scale_x = params.scale_y = 1.0f / 64.0f;
    HS_EXPECT_TRUE(PB::Fields::valid(params));
  }
  auto &modulation = param_as<In::Op::ModulationParams>(program, 14).modulation;
  modulation.source = Pattern::SourceKind::DISPLACEMENT;
  modulation.evaluation = Pattern::Evaluation::DIRECT;
  modulation.mechanic.amount = 1e-18f;
  auto ctx = shared_resources().context();
  std::array<Pixel, PB::Color::HueRotationLutView::SIZE> hue;
  for (size_t index = 0; index < hue.size(); ++index)
    hue[index] = Pixel(static_cast<uint16_t>((index % 16) * 4000), 0, 0);
  ctx.hue_rotation_lut = hue.data();
  program.prepare(ctx);
  const PB::SphereSample sphere{math::Vector(1, 0, 1).normalized(), 0.0f};
  const auto lattice = In::Op::SampleGridV4::translation_lattice(
      param_as<In::Op::GridSampleParams>(program, 12));
  PB::PlaneSample plane{};
  double expected_path = 0.0;
  for (int index = 0; index < 12; ++index) {
    const PB::PlaneSample next =
        index == 0 ? run_entry<PB::PlaneSample>(program, ctx, 0, sphere)
                   : run_entry<PB::PlaneSample>(program, ctx, index, plane);
    if (index != 0) {
      const math::Complex delta = next.coords - plane.coords;
      if (index == 11) {
        HS_EXPECT_GT(std::abs(delta.re / lattice.a.re), 0x1p31f);
        expected_path += std::hypot(std::remainder(delta.re, lattice.a.re),
                                    std::remainder(delta.im, lattice.b.im));
      } else {
        expected_path += std::hypot(static_cast<double>(delta.re),
                                    static_cast<double>(delta.im));
      }
    }
    plane = next;
  }
  HS_EXPECT_TRUE(std::isfinite(plane.path_length));
  HS_EXPECT_NEAR(plane.path_length, expected_path, expected_path * 1e-6);
  plane.path_length = static_cast<float>(expected_path);
  const auto expected_field =
      run_entry<PB::FieldSample>(program, ctx, 12, plane);
  const auto base =
      run_entry<PB::ColorSample>(program, ctx, 13, expected_field);
  const Color4 expected = run_entry<Color4>(program, ctx, 14, base);
  const Color4 actual = program.evaluate(sphere.dir, ctx);
  HS_EXPECT_GT(expected.color.r, 1000);
  HS_EXPECT_EQ(expected.color.g, 0);
  HS_EXPECT_EQ(expected.color.b, 0);
  HS_EXPECT_NEAR(actual.color.r, expected.color.r, 2);
  HS_EXPECT_EQ(actual.color.g, 0);
  HS_EXPECT_EQ(actual.color.b, 0);
}

/** Erased-vs-bound parity of every tessellation kind and fill. */
template <Pattern::TessellationKind Kind, Pattern::TessellationFill Fill>
inline void run_tessellation_parity(In::ChainProgram &program,
                                    const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::TessellationSampleParams>(program, 2);
  params.kind = static_cast<uint8_t>(Kind);
  params.fill = static_cast<uint8_t>(Fill);
  HS_CONTEXT("kind / fill", static_cast<int>(Kind), static_cast<int>(Fill));
  run_sample_op_matrix<
      PB::Source::Tessellation<TessellationSampleMirror, Kind, Fill>>(program,
                                                                      ctx);
}

inline void test_shader_chain_parity_sample_tessellation() {
  using K = Pattern::TessellationKind;
  using F = Pattern::TessellationFill;
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::TessellationSampleParams>(
            program, In::Op::SampleTessellation::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_tessellation_parity<K::TRIANGULAR, F::CELLS>(program, ctx);
    run_tessellation_parity<K::SQUARE, F::CELLS>(program, ctx);
    run_tessellation_parity<K::HEXAGONAL, F::CELLS>(program, ctx);
    run_tessellation_parity<K::TRIANGULAR, F::EDGES>(program, ctx);
    run_tessellation_parity<K::SQUARE, F::EDGES>(program, ctx);
    run_tessellation_parity<K::HEXAGONAL, F::EDGES>(program, ctx);
  });
}

inline void test_shader_chain_parity_sample_projected_noise() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::ProjectedNoiseSampleParams>(
            program, In::Op::SampleProjectedNoise::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    uint8_t *basis = topology_byte(program, 2, "basis");
    HS_EXPECT_TRUE(basis != nullptr);
    if (!basis)
      return;
    *basis = static_cast<uint8_t>(math::NoiseBasis::SIMPLEX);
    run_sample_op_matrix<PB::Source::ProjectedNoise<ProjectedNoiseSampleMirror,
                                                    math::NoiseBasis::SIMPLEX>>(
        program, ctx);
    *basis = static_cast<uint8_t>(math::NoiseBasis::FBM3);
    run_sample_op_matrix<PB::Source::ProjectedNoise<ProjectedNoiseSampleMirror,
                                                    math::NoiseBasis::FBM3>>(
        program, ctx);
    *basis = static_cast<uint8_t>(math::NoiseBasis::RIDGED3);
    run_sample_op_matrix<PB::Source::ProjectedNoise<ProjectedNoiseSampleMirror,
                                                    math::NoiseBasis::RIDGED3>>(
        program, ctx);
  });
}

template <typename WeightP>
inline void run_grid_coverage_cases(In::ChainProgram &program,
                                    const In::FrameContext &ctx,
                                    In::Op::ProjectionCoverageMode coverage) {
  constexpr auto ENV = PB::Color::BrightnessEnvelope::NONE;
  with_coverage_policy<MirrorValue>(coverage, [&]<typename Policy>() {
    expect_parity<WeightP, Policy, ENV>(program, ctx);
  });
}

inline void run_sample_variant(In::ChainProgram &program,
                               const In::FrameContext &ctx,
                               In::Op::WeightMode weight,
                               In::Op::ProjectionCoverageMode coverage) {
  auto &params = param_as<In::Op::GridSampleParams>(program, 2);
  params.weight_mode = static_cast<uint8_t>(weight);
  params.coverage_mode = static_cast<uint8_t>(coverage);
  if (weight == In::Op::WeightMode::NONE)
    run_grid_coverage_cases<PB::Weight::None>(program, ctx, coverage);
  else
    run_grid_coverage_cases<PB::Weight::Projection>(program, ctx, coverage);
}

inline void test_shader_chain_parity_sample_variants() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_default_chain(program, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    for (const In::Op::WeightMode weight :
         {In::Op::WeightMode::NONE, In::Op::WeightMode::PROJECTION})
      for (uint8_t c = 0; c < std::size(In::Op::COVERAGE_MODE_IDS); ++c) {
        const auto coverage = static_cast<In::Op::ProjectionCoverageMode>(c);
        HS_CONTEXT("weight / coverage", int(weight), int(coverage));
        run_sample_variant(program, ctx, weight, coverage);
      }
  });
}

/** @brief Palette mapping and envelopes match the typed color stage. */
template <In::Op::EnvelopeMode EnvelopeV>
void run_colorize_variant(In::ChainProgram &program,
                          const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::GeneratedPaletteParams>(program, 3);
  params.envelope_mode = static_cast<uint8_t>(EnvelopeV);
  for (uint8_t palette_mode = 0;
       palette_mode < std::size(In::Op::PALETTE_MODE_IDS); ++palette_mode) {
    params.palette_mode = palette_mode;
    for (uint8_t mapping = 0; mapping < std::size(In::Op::PALETTE_MAPPING_IDS);
         ++mapping) {
      params.mapping_mode = mapping;
      expect_parity<PB::Weight::Projection, PB::ProjectionCoverage::Weight,
                    EnvelopeV>(program, ctx);
      using Bound = PB::Stage::Colorize<PB::Color::GeneratedPalette<
          MirrorColor<EnvelopeV>>>::template Bind<MirrorBinding>;
      const PB::FieldSample SAMPLE{.4f, 1.0f, math::Y_AXIS, .37f};
      const MirrorFrame mirror = mirror_from(program, ctx);
      const Color4 expected = PB::finish_color_sample(
          Bound::run(SAMPLE, mirror, Bound::prepare(mirror)));
      HS_EXPECT_TRUE(color4_identical(
          run_entry<Color4>(program, ctx, 3, SAMPLE), expected));
    }
  }
}

template <size_t... Envelopes>
void run_colorize_envelopes(In::ChainProgram &program,
                            const In::FrameContext &ctx,
                            std::index_sequence<Envelopes...>) {
  (run_colorize_variant<static_cast<In::Op::EnvelopeMode>(Envelopes)>(program,
                                                                      ctx),
   ...);
}

inline void test_shader_chain_parity_colorize_variants() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_default_chain(program, 4, set))
      return;
    run_colorize_envelopes(
        program, shared_resources().context(),
        std::make_index_sequence<
            static_cast<size_t>(In::Op::EnvelopeMode::DESCENDING) + 1>{});
  });
}

/** @brief Adjacent displacement and spatial hue layers add and cancel. */
inline void test_shader_chain_combined_color_modulation() {
  const PB::FieldSample FIELD{.413f, .6f, math::Y_AXIS, 2.0f};
  for (bool table : {false, true}) {
    auto ctx = shared_resources().context();
    if (!table)
      ctx.hue_rotation_lut = nullptr;
    const auto palette = In::Op::prepare_generated_palette(
        ctx, In::Op::GeneratedPaletteParams{}, {});
    const auto base = PB::Color::apply_generated_palette(FIELD, palette);
    PB::Color::ModulationState path;
    path.enabled = true;
    path.source = Pattern::SourceKind::DISPLACEMENT;
    path.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, .125f};
    PB::Modulation::Prepared source;
    source.evaluate = +[](const PB::Modulation::Prepared &p,
                          const Pattern::SphereDirection &) {
      return p.state.phase.phase_turns;
    };
    PB::Color::ModulationState spatial;
    spatial.enabled = true;
    spatial.current = &source;
    spatial.source = Pattern::SourceKind::NOISE;
    for (float scalar : {-1.0f, 0.0f, 1.0f}) {
      source.state.phase.phase_turns = scalar;
      for (float amount : {0.0f, .125f, -.25f, .75f}) {
        spatial.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, amount};
        const auto actual = PB::Color::finish(PB::Color::apply_modulation(
            PB::Color::apply_modulation(base, path), spatial));
        const float rotation = math::wrap_t(.25f + scalar * amount);
        Color4 expected = base.color;
        if (rotation != 0)
          expected.color =
              table
                  ? sample_hue_rotation_lut(palette.hue_rotation, FIELD.value,
                                            rotation)
                  : ColorModulation::HueRotation::apply(expected, 1, {rotation})
                        .color;
        HS_EXPECT_TRUE(color4_identical(actual, expected));
      }
    }
  }
}

/** @brief Displacement is an ordinary independently enabled modulation source. */
inline void test_shader_chain_combined_color_admission() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  WB::ensure_modulation(effect, "displacement");
  WB::ensure_modulation(effect, "modulation");
  HS_EXPECT_EQ(effect.updateParameter("displacement.modulation-evaluation", 0),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displacement.modulation-source", 4),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displacement.modulation-amount", .5f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("modulation.modulation-amount", .25f),
               ParamSetResult::APPLIED);
  for (const auto *name :
       {"displacement.modulation-enabled", "modulation.modulation-enabled"}) {
    HS_EXPECT_EQ(effect.updateParameter(name, 0), ParamSetResult::APPLIED);
    HS_EXPECT_EQ(effect.updateParameter(name, 1), ParamSetResult::APPLIED);
  }
  HS_EXPECT_EQ(effect.updateParameter("colorize.hue-shift-mode", 1),
               ParamSetResult::UNKNOWN_PARAM);
}

// --- source motion equivalences --------------------------------------------

/** @brief The grid's reported square lattice translates its field onto
    itself, and half a period does not. */
inline void test_shader_chain_grid_lattice_symmetry() {
  In::Op::GridSampleParams params;
  params.pattern_freq = 3.7f;
  params.complexity = 1.4f;
  params.pattern_mix = 0.3f;
  const auto lattice = In::Op::SampleGridV4::translation_lattice(params);
  HS_EXPECT_NEAR(lattice.a.re, math::TWO_PI_F / 3.7f, 1e-7f);
  HS_EXPECT_EQ(lattice.a.im, 0.0f);
  HS_EXPECT_EQ(lattice.b.re, 0.0f);
  HS_EXPECT_NEAR(lattice.b.im, math::TWO_PI_F / 3.7f, 1e-7f);
  const auto prepared = Pattern::prepare_phases(0.9f, 0.0f, 0.0f);
  const auto field = [&](math::Complex point) {
    return Pattern::grid(
        Pattern::PlanePosition{
            projections::stereo_pattern_args(point, params.pattern_freq)},
        params, prepared);
  };
  float half_shift = 0.0f;
  for (const math::Complex point :
       {math::Complex{0.13f, -0.41f}, math::Complex{-1.7f, 0.6f},
        math::Complex{2.9f, 3.3f}, math::Complex{-0.05f, 0.02f}}) {
    for (const math::Complex basis : {lattice.a, lattice.b}) {
      HS_EXPECT_NEAR(field(point + basis), field(point), 2e-4f);
      HS_EXPECT_NEAR(field(point - basis), field(point), 2e-4f);
    }
    half_shift = std::max(
        half_shift,
        std::fabs(field(point + math::Complex{0.5f * lattice.a.re, 0.0f}) -
                  field(point)));
  }
  HS_EXPECT_GT(half_shift, 0.1f);
}

/** @brief An affine warp winding (-5, 4) grid periods per phase turn and
    turning at the source angle, before the stationary-direct grid with
    phase p + s, reproduces the drifting rotated coupled grid. */
inline void test_shader_chain_grid_motion_equivalence() {
  constexpr float FREQ = 4.4f, SPEED = 0.04f, DRIFT = 0.8f, ANGLE_RATE = 0.027f;
  const auto coupled_grid = [&](math::Complex point, float complexity,
                                float mix, float p, float s, float angle) {
    const auto scaled = projections::stereo_pattern_args(point, FREQ);
    const float c = std::cos(angle), sn = std::sin(angle);
    const float x = scaled.re * c + scaled.im * sn;
    const float y = -scaled.re * sn + scaled.im * c;
    const float coupled =
        math::fast_sinf(x + p + complexity * math::fast_sinf(y + p)) *
        math::fast_cosf(y - s + complexity * math::fast_cosf(x - s));
    const float direct = math::fast_sinf(x + p) * math::fast_cosf(y - s);
    return coupled + (direct - coupled) * mix;
  };
  In::Op::GridSampleParams grid;
  grid.pattern_freq = FREQ;
  grid.speed = SPEED * (1.0f + DRIFT);
  PB::Warp::AffineParams affine;
  affine.translation_x = -5.0f;
  affine.translation_y = 4.0f;
  affine.speed = SPEED / (5.0f * math::TWO_PI_F);
  affine.rotation_rate = ANGLE_RATE / affine.speed;
  const auto lattice = In::Op::SampleGridV4::translation_lattice(grid);
  float worst = 0.0f;
  for (const float complexity : {0.0f, 1.7f})
    for (const float mix : {0.0f, 0.4f, 1.0f})
      for (const int frame : {1, 7, 60, 333, 784}) {
        grid.complexity = complexity;
        grid.pattern_mix = mix;
        const double steps = frame;
        const float p = float(std::fmod(steps * SPEED, 2.0 * std::numbers::pi));
        const float s =
            float(std::fmod(steps * SPEED * DRIFT, 2.0 * std::numbers::pi));
        const float angle =
            float(std::fmod(steps * ANGLE_RATE, 2.0 * std::numbers::pi));
        const float phase =
            float(steps * affine.speed - std::floor(steps * affine.speed));
        const auto prepared = PB::Warp::prepare(affine, phase, angle, lattice);
        const auto source = Pattern::prepare_phases(
            float(std::fmod(steps * grid.speed, 2.0 * std::numbers::pi)), 0.0f,
            0.0f);
        for (const math::Complex point :
             {math::Complex{0.21f, -0.37f}, math::Complex{-1.3f, 0.8f},
              math::Complex{2.4f, 1.9f}}) {
          const float actual = Pattern::grid(
              Pattern::PlanePosition{projections::stereo_pattern_args(
                  PB::Warp::affine_coords(point, prepared), FREQ)},
              grid, source);
          worst = std::max(worst,
                           std::fabs(actual - coupled_grid(point, complexity,
                                                           mix, p, s, angle)));
        }
      }
  HS_EXPECT_LT(worst, 1e-3f);
}

/** @brief An affine rotation of -(a + p / 3) before the stationary
    three-armed spiral reproduces sin(r - 3 (azimuth + a) - p). */
inline void test_shader_chain_spiral_rotation_equivalence() {
  PB::Warp::AffineParams affine;
  float worst = 0.0f;
  for (const float p : {0.0f, 0.4f, 2.9f})
    for (const float a : {0.0f, 0.3f, 5.1f}) {
      const auto prepared = PB::Warp::prepare(affine, 0.0f, -(a + p / 3.0f),
                                              PB::TranslationLattice{});
      for (const math::Complex point :
           {math::Complex{0.6f, 0.25f}, math::Complex{-0.4f, 0.9f},
            math::Complex{-0.7f, -0.2f}, math::Complex{0.1f, -0.8f}}) {
        const float expected =
            std::sin(std::hypot(point.re, point.im) -
                     3.0f * (std::atan2(point.im, point.re) + a) - p);
        const float actual = Pattern::spiral(
            Pattern::PlanePosition{PB::Warp::affine_coords(point, prepared)},
            3.0f);
        worst = std::max(worst, std::fabs(actual - expected));
      }
    }
  HS_EXPECT_LT(worst, 0.02f);
}

/** @brief Zero spiral arms draws the concentric rings sin(r). */
inline void test_shader_chain_spiral_zero_arms_rings() {
  for (const math::Complex point :
       {math::Complex{}, math::Complex{0.6f, 0.25f}, math::Complex{-3.0f, 4.0f},
        math::Complex{0.0f, -7.5f}})
    HS_EXPECT_EQ(Pattern::spiral(Pattern::PlanePosition{point}, 0.0f),
                 math::fast_sinf(point.magnitude()));
  In::Op::SpiralSampleParams params;
  params.arms = 0.4f;
  HS_EXPECT_EQ(In::Op::SampleSpiralV3::prepare({}, params, {}).arms, 0.0f);
}

/** @brief Every planar source kind matches its legacy single-kind field. */
inline void test_shader_chain_planar_legacy_endpoints() {
  namespace Op = In::Op;
  using Kind = PB::Source::PlanarKind;
  using Planar = Op::SamplePlanar;
  Planar::Params planar;
  Planar::State state;
  Planar::init(state,
               {"sample", Planar::ID, In::instance_hash("sample", Planar::ID)});
  const In::FrameContext ctx;
  planar.grid.pattern_freq = 3.0f;
  planar.twin_wave.speed = .25f;
  planar.bands.angle = .75f;
  planar.spiral.coordinate_system = 1;
  planar.noise.noise_time_rate = .025f;
  planar.noise.noise_contrast = .4f;
  for (int frame = 0; frame < 5; ++frame)
    Planar::advance(state, planar);

  const auto compare = [&]<typename Legacy>(Kind kind,
                                            const typename Legacy::Params &p,
                                            const typename Legacy::State &old) {
    planar.source = static_cast<uint8_t>(kind);
    const auto a = Planar::prepare(ctx, planar, state);
    const auto b = Legacy::prepare(ctx, p, old);
    HS_EXPECT_EQ(a.count, 1u);
    for (const math::Complex xy : {math::Complex{-.25f, .5f},
                                   {.375f, .87f},
                                   {-.92f, -.1f},
                                   {.75f, -.6f}}) {
      PB::PlaneSample input{xy, {}, math::UP, .3f};
      input.provenance.value_weight = .75f;
      const auto actual = Planar::run(input, ctx, planar, a);
      const auto expected = Legacy::run(input, ctx, p, b);
      HS_EXPECT_EQ(actual.value, expected.value);
      HS_EXPECT_EQ(actual.coverage, expected.coverage);
      HS_EXPECT_EQ(actual.path_length, expected.path_length);
    }
  };
  Op::GridSampleParams grid;
  static_cast<PB::Source::GridSourceParams &>(grid) = planar.grid;
  compare.template operator()<Op::SampleGridV4>(Kind::GRID, grid, state.grid);
  Op::TwinWaveSampleParams twin_wave;
  static_cast<PB::Source::AuthoredTwinWaveSourceParams &>(twin_wave) =
      planar.twin_wave;
  compare.template operator()<Op::SampleTwinWaveV3>(Kind::TWIN_WAVE, twin_wave,
                                                    state.twin_wave);
  Op::BandsSampleParams bands;
  static_cast<PB::Source::BandsSourceParams &>(bands) = planar.bands;
  compare.template operator()<Op::SampleBands>(Kind::BANDS, bands, {});
  Op::SpiralSampleParams spiral;
  static_cast<PB::Source::SpiralSourceParams &>(spiral) = planar.spiral;
  compare.template operator()<Op::SampleSpiralV3>(Kind::SPIRAL, spiral, {});
  Op::TessellationSampleParams tessellation;
  static_cast<PB::Source::TessellationSourceParams &>(tessellation) =
      planar.tessellation;
  uint8_t option = static_cast<uint8_t>(Kind::TESSELLATION_TRIANGULAR_CELLS);
  for (uint8_t shape = 0; shape < 3; ++shape)
    for (uint8_t fill = 0; fill < 2; ++fill) {
      tessellation.kind = shape;
      tessellation.fill = fill;
      compare.template operator()<Op::SampleTessellation>(
          static_cast<Kind>(option++), tessellation, {});
    }

  Op::ProjectedNoiseSampleParams noise;
  static_cast<PB::Source::NoiseSourceParams &>(noise) = planar.noise;
  Op::SampleProjectedNoise::State legacy_noise;
  Op::SampleProjectedNoise::init(
      legacy_noise,
      {"sample", Op::SampleProjectedNoise::ID,
       In::instance_hash("sample", Op::SampleProjectedNoise::ID)});
  for (int frame = 0; frame < 5; ++frame)
    Op::SampleProjectedNoise::advance(legacy_noise, noise);
  HS_EXPECT_EQ(state.noise.noise_seed, legacy_noise.noise_seed);
  for (uint8_t basis = 0; basis < 3; ++basis) {
    noise.basis = basis;
    compare.template operator()<Op::SampleProjectedNoise>(
        static_cast<Kind>(option++), noise, legacy_noise);
  }
}

/** @brief New and departing source fields use their live endpoint values. */
inline void test_shader_chain_planar_column_weights_and_fields() {
  namespace Op = In::Op;
  using Kind = PB::Source::PlanarKind;
  using Planar = Op::SamplePlanar;
  constexpr uint8_t SPIRAL = static_cast<uint8_t>(Kind::SPIRAL);
  constexpr uint8_t CELLS =
      static_cast<uint8_t>(Kind::TESSELLATION_SQUARE_CELLS);
  constexpr uint8_t EDGES =
      static_cast<uint8_t>(Kind::TESSELLATION_SQUARE_EDGES);

  const auto fixture = std::make_unique<ProgramFixture>();
  const auto temporal = std::make_unique<TemporalBlocks>();
  auto &program = fixture->program;
  program.bind_temporal_storage(temporal->history, temporal->column);
  const In::ChainEntryRequest CHAIN[] = {
      {"project", Op::ProjectStereographic::ID},
      {"sample", Planar::ID},
      {"color", Op::ColorizeGeneratedPaletteV4::ID},
  };
  if (!expect_compiled(program.compile(CHAIN)))
    return;
  const auto ctx = shared_resources().context();
  auto &current = param_as<Planar::Params>(program, 1);
  current.source = SPIRAL;
  current.spiral.arms = 3.0f;
  current.spiral.coordinate_system = 1;
  current.tessellation.cell_radius = .15f;
  program.prepare(ctx);
  program.capture_history();

  const auto from = current;
  current.set_enum_blend(In::EnumBlend{SPIRAL, CELLS, .4f});
  current.tessellation.cell_radius = .4f;
  HS_EXPECT_TRUE(program.history_continuous());
  PB::PlaneSample input{{.375f, .87f}, {}, math::UP, 0.0f};
  input.provenance.value_weight = 1.0f;
  const auto &state = state_as<Planar::State>(program, 1);
  const auto endpoint = [&](const auto &params) {
    return Planar::run(input, ctx, params, Planar::prepare(ctx, params, state))
        .value;
  };
  auto cells = current;
  cells.source = CELLS;
  cells.set_enum_blend(std::nullopt);
  const auto evaluate_column = [&](const auto &a, const auto &b,
                                   float fraction) {
    alignas(Planar::Params) uint8_t column[sizeof(Planar::Params)];
    program.ops()[1].op->runtime.sample_params(column, &a, &b, fraction);
    return *std::launder(reinterpret_cast<const Planar::Params *>(column));
  };
  const auto first = evaluate_column(from, current, .5f);
  HS_EXPECT_EQ(first.tessellation.cell_radius, .4f);
  HS_EXPECT_EQ(first.weights().values[SPIRAL], .8f);
  HS_EXPECT_EQ(first.weights().values[CELLS], .2f);
  HS_EXPECT_EQ(Planar::prepare(ctx, first, state).count, 2u);
  HS_EXPECT_NEAR(endpoint(first), .8f * endpoint(from) + .2f * endpoint(cells),
                 1e-6f);

  program.prepare(ctx);
  program.capture_history();
  const auto previous = current;
  current = cells;
  current.spiral.arms = 9.0f;
  current.spiral.coordinate_system = 0;
  HS_EXPECT_TRUE(program.history_continuous());
  const auto last = evaluate_column(previous, current, .5f);
  HS_EXPECT_EQ(last.spiral.arms, 3.0f);
  HS_EXPECT_EQ(last.spiral.coordinate_system, 1u);
  HS_EXPECT_NEAR(last.weights().values[SPIRAL], .3f, 1e-6f);
  HS_EXPECT_NEAR(last.weights().values[CELLS], .7f, 1e-6f);

  current.source = EDGES;
  HS_EXPECT_FALSE(program.history_continuous());
  auto fill_to = current;
  fill_to.source = CELLS;
  fill_to.set_enum_blend(In::EnumBlend{CELLS, EDGES, .5f});
  fill_to.tessellation.cell_radius = .2f;
  fill_to.tessellation.line_thickness = .125f;
  const auto fill = evaluate_column(cells, fill_to, .5f);
  HS_EXPECT_EQ(fill.tessellation.line_thickness, .125f);
  HS_EXPECT_NEAR(fill.tessellation.cell_radius, .3f, 1e-6f);
}

template <typename Source> void check_tessellation_endpoint_fields() {
  namespace Op = In::Op;
  using Kind = Pattern::TessellationKind;
  using Fill = Pattern::TessellationFill;
  static int samples = 0;
  auto operators = In::OPERATOR_TABLE;
  for (auto &entry : operators)
    if (std::string_view(entry.operator_id) == Source::ID)
      entry.runtime.sample_source = [](const PB::PlaneSample &input,
                                       const void *params,
                                       const uint8_t *prepared) {
        ++samples;
        return Source::sample(
            input, *static_cast<const typename Source::Params *>(params),
            *reinterpret_cast<const typename Source::Prepared *>(prepared));
      };
  const auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, operators);
  const auto temporal = std::make_unique<TemporalBlocks>();
  auto &program = fixture->program;
  program.bind_temporal_storage(temporal->history, temporal->column);
  const In::ChainEntryRequest chain[] = {
      {"project", Op::ProjectEquirectangular::ID},
      {"sample", Source::ID},
      {"color", Op::ColorizeGeneratedPaletteV4::ID},
  };
  if (!expect_compiled(program.compile(chain)))
    return;
  const auto ctx = shared_resources().context();
  auto &current = param_as<typename Source::Params>(program, 1);
  const auto cells = [](auto &params) -> auto & {
    if constexpr (std::is_same_v<Source, Op::SamplePlanar>) {
      params.source = static_cast<uint8_t>(
          PB::Source::PlanarKind::TESSELLATION_HEXAGONAL_CELLS);
      return params.tessellation;
    } else {
      params.kind = static_cast<uint8_t>(Kind::HEXAGONAL);
      params.fill = static_cast<uint8_t>(Fill::CELLS);
      return params;
    }
  };
  auto from = current;
  auto to = current;
  cells(from).cell_scale = 2.5352101f;
  cells(from).cell_softness = .10124387f;
  cells(to).cell_scale = 4.774648f;
  cells(to).cell_softness = .32801026f;
  current = from;
  program.prepare(ctx);
  program.capture_history();

  current = PB::Fields::interpolate(from, to, .2f);
  program.set_source_morph(1, &from, &to, {1, 4, true, false});
  HS_EXPECT_TRUE(program.history_continuous());
  program.prepare(ctx);
  program.capture_history();
  current = PB::Fields::interpolate(from, to, .5f);
  program.set_source_morph(1, &from, &to, {2, 4, true, false});
  HS_EXPECT_TRUE(program.history_continuous());
  constexpr float COLUMN = .5f;
  const float mix = math::ease_in_out_sin(.375f);
  program.prepare_column(ctx, ctx, COLUMN);

  bool distinct = false;
  for (const auto &view : sweep_views()) {
    const auto plane = run_entry<PB::PlaneSample>(program, ctx, 0,
                                                  PB::SphereSample{view, 0.0f});
    const auto sample = [&](auto &params) {
      return Pattern::tessellation<Kind::HEXAGONAL, Fill::CELLS>(
          Pattern::PlanePosition{plane.coords},
          static_cast<const Pattern::TessellationCellParams &>(cells(params)));
    };
    const float raw = hs::lerp(sample(from), sample(to), mix);
    distinct |= std::abs(sample(current) - raw) > .1f;
    const auto field = Op::finish_sample(plane, raw, current, ctx);
    const auto expected = PB::finish_color_sample(
        run_entry<PB::ColorSample>(program, ctx, 2, field));
    samples = 0;
    const auto actual = program.evaluate_column(view, ctx, ctx, COLUMN);
    HS_EXPECT_EQ(samples, 2);
    HS_EXPECT_EQ(actual.color, expected.color);
    HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
  }
  HS_EXPECT_TRUE(distinct);

  program.prepare(ctx);
  program.capture_history();
  program.set_source_morph(1, &from, &to, {4, 4, true, false});
  HS_EXPECT_FALSE(program.history_continuous());
  program.set_source_morph(1, &from, &to, {3, 4, true, true});
  HS_EXPECT_FALSE(program.history_continuous());

  program.set_source_morph(1, &from, &from, {2, 4, true, false});
  program.prepare(ctx);
  samples = 0;
  program.evaluate(math::UP, ctx);
  HS_EXPECT_EQ(samples, 1);

  program.set_source_morph(1, &from, &to, {4, 4, true, false});
  program.prepare(ctx);
  samples = 0;
  program.evaluate(math::UP, ctx);
  HS_EXPECT_EQ(samples, 1);
}

/** @brief Both tessellation operators sample authored controls at column time. */
inline void test_shader_chain_same_kind_endpoint_fields() {
  check_tessellation_endpoint_fields<In::Op::SampleTessellation>();
  check_tessellation_endpoint_fields<In::Op::SamplePlanar>();
}

/** @brief The two source fields close their own affine repeats across a phase wrap. */
inline void test_shader_chain_affine_morph_phase_wrap() {
  namespace Op = In::Op;
  using Source = Op::SampleTessellation;
  static std::array<PB::PlaneSample, 2> sampled;
  static size_t count;
  auto operators = In::OPERATOR_TABLE;
  for (auto &entry : operators)
    if (std::string_view(entry.operator_id) == Source::ID)
      entry.runtime.sample_source = [](const PB::PlaneSample &input,
                                       const void *params,
                                       const uint8_t *prepared) {
        sampled.at(count++) = input;
        return Source::sample(
            input, *static_cast<const Source::Params *>(params),
            *reinterpret_cast<const Source::Prepared *>(prepared));
      };
  const auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, operators);
  const auto blocks = std::make_unique<TemporalBlocks>();
  auto &program = fixture->program;
  program.bind_temporal_storage(blocks->history, blocks->column);
  const In::ChainEntryRequest chain[] = {
      {"project", Op::ProjectEquirectangular::ID},
      {"warp", Op::WarpAffineV4::ID},
      {"source", Source::ID},
      {"color", Op::ColorizeGeneratedPaletteV4::ID}};
  if (!expect_compiled(program.compile(chain)))
    return;
  const auto ctx = shared_resources().context();
  auto &warp = param_as<Op::AffineWarpParams>(program, 1);
  warp.translation_x = 2.0f;
  warp.translation_y = -1.0f;
  warp.scale_x = 1.3f;
  warp.shear = .2f;
  warp.path_length = 1;
  auto &live_source = param_as<Source::Params>(program, 2);
  live_source.kind = static_cast<uint8_t>(Pattern::TessellationKind::HEXAGONAL);
  auto from = live_source, to = live_source;
  from.cell_scale = 2.5352101f;
  to.cell_scale = 4.774648f;
  live_source = PB::Fields::interpolate(from, to, .5f);
  program.set_source_morph(2, &from, &to, {2, 4, true, false});
  auto &clock = state_as<Op::AffineClockState>(program, 1);
  constexpr float before = .99998f, after = .00002f;
  clock = {before, .37f};
  program.prepare(ctx);
  program.capture_history();
  clock.phase = after;
  const auto interval = Op::temporal_interval({before, clock.rotation}, clock);
  const auto view = math::Vector(.3f, .7f, .4f).normalized();
  const auto plane =
      run_entry<PB::PlaneSample>(program, ctx, 0, PB::SphereSample{view, .0f});
  const auto check = [&](float phase) {
    const Op::AffineClockState at{phase, clock.rotation};
    const auto warped = [&](const Source::Params &source) {
      return Op::WarpAffineV4::run(
          plane, ctx, warp,
          Op::WarpAffineV4::prepare(ctx, warp, at,
                                    Source::translation_lattice(source)));
    };
    const auto live = warped(live_source);
    HS_EXPECT_EQ(count, size_t{2});
    for (size_t i = 0; i < 2; ++i) {
      const auto expected = warped(i == 0 ? from : to);
      HS_EXPECT_NEAR(sampled[i].coords.re, expected.coords.re, 2e-6f);
      HS_EXPECT_NEAR(sampled[i].coords.im, expected.coords.im, 2e-6f);
      HS_EXPECT_NEAR(sampled[i].path_length, live.path_length, 2e-6f);
    }
  };
  for (const float t : {0.0f, .25f, .75f, 1.0f}) {
    program.prepare_column(ctx, ctx, t);
    count = 0;
    program.evaluate_column(view, ctx, ctx, t);
    check(interval.phase.sample(t));
  }
  clock.phase = before;
  program.prepare(ctx);
  count = 0;
  program.evaluate(view, ctx);
  check(before);
}

/** @brief Planar preview admission checks translating and path-tracking affines. */
inline void test_shader_chain_planar_lattice_admission() {
  namespace Op = In::Op;
  using Kind = PB::Source::PlanarKind;
  Op::PlanarSampleParams from, to;
  to.source = static_cast<uint8_t>(Kind::BANDS);
  Op::AffineWarpParams affine;
  affine.path_length = 0;
  auto &bands = to.bands;
  bands.pattern_freq = 2.0f;
  HS_EXPECT_TRUE(Op::admits_planar_transition(from, to, affine, affine));
  affine.path_length = 1;
  HS_EXPECT_FALSE(Op::admits_planar_transition(from, to, affine, affine));
  affine.path_length = 0;
  affine.translation_x = 1.0f;
  HS_EXPECT_FALSE(Op::admits_planar_transition(from, to, affine, affine));
  bands.pattern_freq = 1.0f;
  HS_EXPECT_TRUE(Op::admits_planar_transition(from, to, affine, affine));
  bands.angle = math::PI_F / 2;
  HS_EXPECT_FALSE(Op::admits_planar_transition(from, to, affine, affine));
  to.source = static_cast<uint8_t>(Kind::SPIRAL);
  affine.lattice_period = math::TWO_PI_F;
  HS_EXPECT_TRUE(Op::admits_planar_transition(from, to, affine, affine));
}

/** @brief Non-pausing preview writes preserve clocks and normal writes snap. */
inline void test_shader_chain_planar_preview_writes_and_clock_migration() {
  namespace Op = In::Op;
  using WB = ShaderChainWhiteBox;
  using Kind = PB::Source::PlanarKind;
  reset_globals();
  WB::FX effect;
  effect.init();
  effect.setAnimationsPaused(false);
  auto &program = WB::program(effect);
  auto &grid = state_as<Op::SourceClockState>(program, 2);
  grid.primary = .875f;
  const In::ChainEntryRequest CHAIN[] = {
      {"camera", Op::Rotate::ID},
      {"project", Op::ProjectStereographic::ID},
      {"sample", Op::SamplePlanar::ID},
      {"colorize", Op::ColorizeGeneratedPaletteV4::ID},
  };
  if (!expect_compiled(effect.set_chain(CHAIN)))
    return;
  const auto &planar = state_as<Op::PlanarSourceState>(program, 2);
  HS_EXPECT_EQ(planar.grid.primary, .875f);
  const auto capture = effect.snapshot();
  HS_EXPECT_EQ(capture.schema_version, 6u);
  HS_EXPECT_EQ(effect.restore_snapshot(capture),
               ChainSnapshotRestoreResult::APPLIED);
  HS_EXPECT_EQ(state_as<Op::PlanarSourceState>(program, 2).grid.primary, .875f);

  constexpr uint8_t GRID = static_cast<uint8_t>(Kind::GRID);
  constexpr uint8_t SPIRAL = static_cast<uint8_t>(Kind::SPIRAL);
  const ShaderChainParameterWrite morph[]{
      {"sample.source", 0, In::EnumBlend{GRID, SPIRAL, .5f}},
      {"colorize.palette-mapping", 0, In::EnumBlend{0, 2, .25f}},
  };
  HS_EXPECT_EQ(effect.update_parameters(morph, true), ParamSetResult::APPLIED);
  HS_EXPECT_FALSE(effect.animations_paused());
  const auto &source = param_as<Op::PlanarSampleParams>(program, 2);
  HS_EXPECT_EQ(source.weights().values[SPIRAL], .5f);
  HS_EXPECT_EQ(WB::color_params(effect).mapping_blend.weights.values[2], .25f);

  program.prepare(shared_resources().context());
  program.capture_history();
  const ShaderChainParameterWrite snap[]{
      {"sample.source", static_cast<float>(Kind::BANDS)},
  };
  HS_EXPECT_EQ(effect.update_parameters(snap), ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_FALSE(program.history_continuous());
  HS_EXPECT_FALSE(source.source_blend.active);
  HS_EXPECT_FALSE(WB::color_params(effect).mapping_blend.active);
  HS_EXPECT_EQ(source.weights().exact, static_cast<uint8_t>(Kind::BANDS));
}

} // namespace hs_test::shader_chain_tests
