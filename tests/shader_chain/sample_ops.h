/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_shader_chain.h.

// --- sample-op mirrors ----------------------------------------------------

/** Mirror frame for the sample source batch. */
struct SampleMirrorFrame {
  const FastNoiseLite *noise = nullptr;
  In::Op::TwinWaveSampleParams twin_wave;
  In::Op::RingsSampleParams rings;
  In::Op::SphericalRingsSampleParams spherical_rings;
  In::Op::SpiralSampleParams spiral;
  In::Op::LatticeSampleParams lattice;
  In::Op::FractalSampleParams fractal;
  In::Op::TessellationSampleParams tessellation;
  In::Op::ProjectedNoiseSampleParams projected;
  In::Op::SphericalNoiseSampleParams spherical;
  math::Vector ring_axis = math::Y_AXIS;
  float ring_phase = 0.0f;
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

template <typename Derived> struct ClockedSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static PB::Source::PreparedSource prepare(const FrameState &frame) {
    return PB::Source::prepare(frame.primary, frame.secondary, frame.angle);
  }
};

struct TwinWaveSampleMirror : ClockedSampleMirror<TwinWaveSampleMirror> {
  static const In::Op::TwinWaveSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.twin_wave;
  }
};

struct RingsSampleMirror : ClockedSampleMirror<RingsSampleMirror> {
  static const In::Op::RingsSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.rings;
  }
};

struct SphericalRingsSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const In::Op::SphericalRingsSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.spherical_rings;
  }
  static PB::Source::PreparedSphericalRings
  prepare(const SampleMirrorFrame &frame) {
    return {frame.ring_axis, frame.ring_phase};
  }
};

struct SpiralSampleMirror : ClockedSampleMirror<SpiralSampleMirror> {
  static const In::Op::SpiralSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.spiral;
  }
};

struct LatticeSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const In::Op::LatticeSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.lattice;
  }
};

struct FractalSampleMirror : ClockedSampleMirror<FractalSampleMirror> {
  static const In::Op::FractalSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.fractal;
  }
};

struct TessellationSampleMirror
    : ClockedSampleMirror<TessellationSampleMirror> {
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

struct SphericalNoiseSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const FastNoiseLite &noise(const SampleMirrorFrame &frame) {
    return *frame.noise;
  }
  static float noise_scale(const SampleMirrorFrame &frame) {
    return frame.spherical.noise_scale;
  }
  static float noise_time(const SampleMirrorFrame &frame) {
    return frame.noise_time;
  }
  static float noise_contrast(const SampleMirrorFrame &frame) {
    return frame.spherical.noise_contrast;
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
inline void arm_sample_op_chain(In::ChainProgram &program, const char *op_id,
                                int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"source", op_id},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 1), set);
  apply_value_set(param_as<OpParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 3), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

template <typename OpParams>
inline void arm_spherical_sample_op_chain(In::ChainProgram &program,
                                          const char *op_id, int frames,
                                          ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"source", op_id},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<OpParams>(program, 1), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 2), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

inline SampleMirrorFrame sample_mirror(In::ChainProgram &program,
                                       size_t source_index = 2) {
  SampleMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[source_index].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::SampleProjectedNoise::ID ||
      id == In::Op::SampleSphericalNoise::ID) {
    const auto &state =
        state_as<In::Op::NoisePhaseState>(program, source_index);
    mirror.noise = &state.noise;
    mirror.noise_time = state.phase;
    if (id == In::Op::SampleProjectedNoise::ID) {
      mirror.projected =
          param_as<In::Op::ProjectedNoiseSampleParams>(program, source_index);
      mirror.edge_width = mirror.projected.edge_width;
    } else {
      mirror.spherical =
          param_as<In::Op::SphericalNoiseSampleParams>(program, source_index);
    }
  } else if (id == In::Op::SampleLattice::ID) {
    mirror.lattice =
        param_as<In::Op::LatticeSampleParams>(program, source_index);
    mirror.edge_width = mirror.lattice.edge_width;
  } else if (id == In::Op::SampleSphericalRings::ID) {
    const auto &state =
        state_as<In::Op::SphericalRingsState>(program, source_index);
    mirror.spherical_rings =
        param_as<In::Op::SphericalRingsSampleParams>(program, source_index);
    const math::Quaternion orientation =
        math::make_rotation(math::X_AXIS, state.walk.spin_phase) *
        state.walk.wander;
    mirror.ring_axis = math::rotate(math::Y_AXIS, orientation);
    mirror.ring_phase = state.phase;
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
    } else if (id == In::Op::SampleRings::ID) {
      mirror.rings = param_as<In::Op::RingsSampleParams>(program, source_index);
      mirror.edge_width = mirror.rings.edge_width;
    } else if (id == In::Op::SampleSpiral::ID) {
      mirror.spiral =
          param_as<In::Op::SpiralSampleParams>(program, source_index);
      mirror.edge_width = mirror.spiral.edge_width;
    } else if (id == In::Op::SampleFractal::ID) {
      mirror.fractal =
          param_as<In::Op::FractalSampleParams>(program, source_index);
      mirror.edge_width = mirror.fractal.edge_width;
    } else if (id == In::Op::SampleTessellation::ID) {
      mirror.tessellation =
          param_as<In::Op::TessellationSampleParams>(program, source_index);
      mirror.edge_width = mirror.tessellation.edge_width;
    }
  }
  return mirror;
}

template <typename SourceP>
inline void expect_spherical_sample_op_parity(In::ChainProgram &program,
                                              const In::FrameContext &ctx) {
  using BoundStage = typename PB::Stage::SampleSphere<SourceP>::template Bind<
      SampleMirrorBinding>;
  program.prepare(ctx);
  const SampleMirrorFrame mirror = sample_mirror(program, 1);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &rotate_op = *program.ops()[0].op;
  const In::OperatorDescriptor &sample_op = *program.ops()[1].op;
  int view_index = 0;
  for (const math::Vector &view : sweep_views()) {
    HS_CONTEXT("view", view_index++);
    const PB::SphereSample seed{view, 0.0f};
    alignas(In::SLOT_ALIGN) uint8_t rotated_out[In::SLOT_SIZE];
    rotate_op.runtime.run(&seed, rotated_out, ctx, program.param_block(0),
                          program.prepared_block(0));
    const auto &input =
        *std::launder(reinterpret_cast<PB::SphereSample *>(rotated_out));
    alignas(In::SLOT_ALIGN) uint8_t sampled_out[In::SLOT_SIZE];
    sample_op.runtime.run(&input, sampled_out, ctx, program.param_block(1),
                          program.prepared_block(1));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::FieldSample *>(sampled_out));
    const PB::FieldSample reference = BoundStage::run(input, mirror, prepared);
    HS_EXPECT_TRUE(field_identical(erased, reference));
  }
}

/** Erased-vs-bound parity of the sample crossing at entry 2. */
template <typename BoundStage>
inline void expect_sample_op_parity(In::ChainProgram &program,
                                    const In::FrameContext &ctx) {
  program.prepare(ctx);
  const SampleMirrorFrame mirror = sample_mirror(program);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[2].op;
  int view_index = 0;
  for (const math::Vector &view : sweep_views()) {
    HS_CONTEXT("view", view_index++);
    const PB::PlaneSample input = warp_input(program, ctx, view);
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&input, out, ctx, program.param_block(2),
                   program.prepared_block(2));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::FieldSample *>(out));
    const PB::FieldSample reference = BoundStage::run(input, mirror, prepared);
    HS_EXPECT_TRUE(field_identical(erased, reference));
  }
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

template <typename SourceP, typename WeightP>
inline void run_sample_coverage_cases(In::ChainProgram &program,
                                      const In::FrameContext &ctx,
                                      uint8_t coverage) {
  using EdgeP = PB::ProjectionCoverage::EdgeFade<SampleEdgeMirror>;
  switch (static_cast<In::Op::ProjectionCoverageMode>(coverage)) {
  case In::Op::ProjectionCoverageMode::NONE:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP,
        PB::ProjectionCoverage::None>::template Bind<SampleMirrorBinding>>(
        program, ctx);
    break;
  case In::Op::ProjectionCoverageMode::WEIGHT:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP,
        PB::ProjectionCoverage::Weight>::template Bind<SampleMirrorBinding>>(
        program, ctx);
    break;
  case In::Op::ProjectionCoverageMode::WEIGHT_SQUARED:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP, PB::ProjectionCoverage::WeightSquared>::
                                template Bind<SampleMirrorBinding>>(program,
                                                                    ctx);
    break;
  case In::Op::ProjectionCoverageMode::EDGE_FADE:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP, EdgeP>::template Bind<SampleMirrorBinding>>(program,
                                                                      ctx);
    break;
  }
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
  for (uint8_t w = 0; w < 2; ++w)
    for (uint8_t c = 0; c < 4; ++c) {
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
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::TwinWaveSampleParams>(
        program, In::Op::SampleTwinWaveV3::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::TwinWave<TwinWaveSampleMirror>>(program,
                                                                     ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_rings() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::RingsSampleParams>(
        program, In::Op::SampleRings::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Rings<RingsSampleMirror>>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_spherical_rings() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_spherical_sample_op_chain<In::Op::SphericalRingsSampleParams>(
        program, In::Op::SampleSphericalRings::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    expect_spherical_sample_op_parity<
        PB::Source::SphericalRings<SphericalRingsSampleMirror>>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_spiral() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::SpiralSampleParams>(
        program, In::Op::SampleSpiral::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Spiral<SpiralSampleMirror>>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_lattice() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::LatticeSampleParams>(
        program, In::Op::SampleLattice::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::PrimitiveLattice<LatticeSampleMirror>>(
        program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_fractal() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::FractalSampleParams>(
        program, In::Op::SampleFractal::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::EscapeFractal<FractalSampleMirror>>(
        program, ctx);
    program.clear();
  }
}

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
      const float raw = PB::Source::noise_contour(
          state.noise, basis, coordinate, sample_params.noise_contrast);
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

inline void test_shader_chain_composed_noise_domain() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest chain[] = {
      {"project", In::Op::ProjectStereographic::ID},
      {"affine-a", In::Op::WarpAffineV3::ID},
      {"affine-b", In::Op::WarpAffineV3::ID},
      {"affine-c", In::Op::WarpAffineV3::ID},
      {"affine-d", In::Op::WarpAffineV3::ID},
      {"affine-e", In::Op::WarpAffineV3::ID},
      {"affine-f", In::Op::WarpAffineV3::ID},
      {"sample", In::Op::SampleProjectedNoise::ID},
      {"color", In::Op::ColorizeGeneratedPaletteV3::ID},
  };
  HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
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
    PB::PlaneSample plane{};
    for (int index = 0; index < 7; ++index) {
      PB::PlaneSample next{};
      program.ops()[index].op->runtime.run(
          index == 0 ? static_cast<const void *>(&sphere)
                     : static_cast<const void *>(&plane),
          &next, ctx, program.param_block(index),
          program.prepared_block(index));
      plane = next;
    }
    HS_EXPECT_GT(plane.coords.re, 1e11f);
    PB::FieldSample field{};
    program.ops()[7].op->runtime.run(
        &plane, &field, ctx, program.param_block(7), program.prepared_block(7));
    HS_EXPECT_EQ(field.value, 0.5f);
    const Color4 actual = program.evaluate(sphere.dir, ctx);
    HS_EXPECT_TRUE(std::isfinite(actual.alpha));
  }
}

inline void test_shader_chain_large_finite_path_length() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  std::array<In::ChainEntryRequest, 14> chain{};
  std::array<std::string, 11> ids;
  chain[0] = {"project", In::Op::ProjectStereographic::ID};
  for (int index = 0; index < 11; ++index) {
    ids[index] = "warp" + std::to_string(index);
    chain[index + 1] = {ids[index], In::Op::WarpAffineV3::ID};
  }
  chain[12] = {"sample", In::Op::SampleGridV3::ID};
  chain[13] = {"color", In::Op::ColorizeGeneratedPaletteV3::ID};
  HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
  for (int index = 1; index <= 11; ++index) {
    auto &params = param_as<In::Op::AffineWarpParams>(program, index);
    params.scale_x = params.scale_y = 1.0f / 64.0f;
    HS_EXPECT_TRUE(PB::Fields::valid(params));
  }
  auto &color = param_as<In::Op::GeneratedPaletteParams>(program, 13);
  color.hue_mode = static_cast<uint8_t>(PB::Color::HueMode::PATH_LENGTH);
  color.hue_shift_amount = 1e-20f;
  auto ctx = shared_resources().context();
  ctx.projection_base = math::Quaternion();
  std::array<Pixel, PB::Color::HueRotationLutView::SIZE> hue;
  for (size_t index = 0; index < hue.size(); ++index)
    hue[index] = Pixel(static_cast<uint16_t>((index % 16) * 4000), 0, 0);
  ctx.hue_rotation_lut = hue.data();
  program.prepare(ctx);
  const PB::SphereSample sphere{math::Vector(1, 0, 1).normalized(), 0.0f};
  PB::PlaneSample plane{};
  double expected_path = 0.0;
  for (int index = 0; index < 12; ++index) {
    PB::PlaneSample next{};
    program.ops()[index].op->runtime.run(
        index == 0 ? static_cast<const void *>(&sphere)
                   : static_cast<const void *>(&plane),
        &next, ctx, program.param_block(index), program.prepared_block(index));
    if (index != 0)
      expected_path +=
          std::hypot(static_cast<double>(next.coords.re) - plane.coords.re,
                     static_cast<double>(next.coords.im) - plane.coords.im);
    plane = next;
  }
  HS_EXPECT_TRUE(std::isfinite(plane.path_length));
  HS_EXPECT_NEAR(plane.path_length, expected_path, expected_path * 1e-6);
  plane.path_length = static_cast<float>(expected_path);
  PB::FieldSample expected_field{};
  program.ops()[12].op->runtime.run(&plane, &expected_field, ctx,
                                    program.param_block(12),
                                    program.prepared_block(12));
  Color4 expected;
  program.ops()[13].op->runtime.run(&expected_field, &expected, ctx,
                                    program.param_block(13),
                                    program.prepared_block(13));
  const Color4 actual = program.evaluate(sphere.dir, ctx);
  HS_EXPECT_GT(expected.color.r, 1000);
  HS_EXPECT_NEAR(actual.color.r, expected.color.r, 2);
}

inline void test_shader_chain_parity_sample_tessellation() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::TessellationSampleParams>(
        program, In::Op::SampleTessellation::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    auto &params = param_as<In::Op::TessellationSampleParams>(program, 2);
    params.kind =
        static_cast<uint8_t>(PB::Source::TessellationKind::TRIANGULAR);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, PB::Source::TessellationKind::TRIANGULAR>>(
        program, ctx);
    params.kind = static_cast<uint8_t>(PB::Source::TessellationKind::SQUARE);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, PB::Source::TessellationKind::SQUARE>>(
        program, ctx);
    params.kind = static_cast<uint8_t>(PB::Source::TessellationKind::HEXAGONAL);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, PB::Source::TessellationKind::HEXAGONAL>>(
        program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_projected_noise() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::ProjectedNoiseSampleParams>(
        program, In::Op::SampleProjectedNoise::ID, 4, set);
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
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_spherical_noise() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_spherical_sample_op_chain<In::Op::SphericalNoiseSampleParams>(
        program, In::Op::SampleSphericalNoise::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    expect_spherical_sample_op_parity<PB::Source::SphericalNoise<
        SphericalNoiseSampleMirror, math::NoiseBasis::SIMPLEX>>(program, ctx);
    auto &params = param_as<In::Op::SphericalNoiseSampleParams>(program, 1);
    params.basis = static_cast<uint8_t>(math::NoiseBasis::FBM3);
    expect_spherical_sample_op_parity<PB::Source::SphericalNoise<
        SphericalNoiseSampleMirror, math::NoiseBasis::FBM3>>(program, ctx);
    params.basis = static_cast<uint8_t>(math::NoiseBasis::RIDGED3);
    expect_spherical_sample_op_parity<PB::Source::SphericalNoise<
        SphericalNoiseSampleMirror, math::NoiseBasis::RIDGED3>>(program, ctx);
    program.clear();
  }
}

inline void run_sample_variant(In::ChainProgram &program,
                               const In::FrameContext &ctx,
                               In::Op::WeightMode weight,
                               In::Op::ProjectionCoverageMode coverage) {
  auto &params = param_as<In::Op::GridSampleParams>(program, 2);
  params.weight_mode = static_cast<uint8_t>(weight);
  params.coverage_mode = static_cast<uint8_t>(coverage);
  constexpr auto HUE = PB::Color::HueMode::NOISE;
  constexpr auto ENV = PB::Color::BrightnessEnvelope::NONE;
  if (weight == In::Op::WeightMode::NONE) {
    switch (coverage) {
    case In::Op::ProjectionCoverageMode::NONE:
      expect_parity<PB::Weight::None, PB::ProjectionCoverage::None, HUE, ENV>(
          program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT:
      expect_parity<PB::Weight::None, PB::ProjectionCoverage::Weight, HUE, ENV>(
          program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT_SQUARED:
      expect_parity<PB::Weight::None, PB::ProjectionCoverage::WeightSquared,
                    HUE, ENV>(program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::EDGE_FADE:
      expect_parity<PB::Weight::None,
                    PB::ProjectionCoverage::EdgeFade<MirrorValue>, HUE, ENV>(
          program, ctx);
      break;
    }
  } else {
    switch (coverage) {
    case In::Op::ProjectionCoverageMode::NONE:
      expect_parity<PB::Weight::Projection, PB::ProjectionCoverage::None, HUE,
                    ENV>(program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT:
      expect_parity<PB::Weight::Projection, PB::ProjectionCoverage::Weight, HUE,
                    ENV>(program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT_SQUARED:
      expect_parity<PB::Weight::Projection,
                    PB::ProjectionCoverage::WeightSquared, HUE, ENV>(program,
                                                                     ctx);
      break;
    case In::Op::ProjectionCoverageMode::EDGE_FADE:
      expect_parity<PB::Weight::Projection,
                    PB::ProjectionCoverage::EdgeFade<MirrorValue>, HUE, ENV>(
          program, ctx);
      break;
    }
  }
}

inline void test_shader_chain_parity_sample_variants() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_default_chain(program, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    for (const In::Op::WeightMode weight :
         {In::Op::WeightMode::NONE, In::Op::WeightMode::PROJECTION})
      for (const In::Op::ProjectionCoverageMode coverage :
           {In::Op::ProjectionCoverageMode::NONE,
            In::Op::ProjectionCoverageMode::WEIGHT,
            In::Op::ProjectionCoverageMode::WEIGHT_SQUARED,
            In::Op::ProjectionCoverageMode::EDGE_FADE})
        run_sample_variant(program, ctx, weight, coverage);
    program.clear();
  }
}

template <In::Op::HueShiftMode HueV, In::Op::EnvelopeMode EnvelopeV>
void run_colorize_variant(In::ChainProgram &program,
                          const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::GeneratedPaletteParams>(program, 3);
  params.hue_mode = static_cast<uint8_t>(HueV);
  params.envelope_mode = static_cast<uint8_t>(EnvelopeV);
  for (uint8_t palette_mode = 0; palette_mode < 3; ++palette_mode) {
    params.palette_mode = palette_mode;
    for (uint8_t mapping = 0; mapping < 4; ++mapping) {
      params.mapping_mode = mapping;
      expect_parity<PB::Weight::Projection, PB::ProjectionCoverage::Weight,
                    HueV, EnvelopeV>(program, ctx);
    }
  }
}

/** Runs one hue mode against every brightness envelope. */
template <In::Op::HueShiftMode HueV, size_t... Envelopes>
void run_colorize_envelopes(In::ChainProgram &program,
                            const In::FrameContext &ctx,
                            std::index_sequence<Envelopes...>) {
  (run_colorize_variant<HueV, static_cast<In::Op::EnvelopeMode>(Envelopes)>(
       program, ctx),
   ...);
}

inline constexpr auto COLORIZE_ENVELOPES = std::make_index_sequence<
    static_cast<size_t>(In::Op::EnvelopeMode::DESCENDING) + 1>{};

inline void test_shader_chain_parity_colorize_variants() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_default_chain(program, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_colorize_envelopes<In::Op::HueShiftMode::NONE>(program, ctx,
                                                       COLORIZE_ENVELOPES);
    run_colorize_envelopes<In::Op::HueShiftMode::NOISE>(program, ctx,
                                                        COLORIZE_ENVELOPES);
    run_colorize_envelopes<In::Op::HueShiftMode::PATH_LENGTH>(
        program, ctx, COLORIZE_ENVELOPES);
    program.clear();
  }
}
