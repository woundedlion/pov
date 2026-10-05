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

/** Renders the sweep into @p out for a byte-identity comparison. */
inline void snapshot_render(In::ChainProgram &program,
                            const In::FrameContext &ctx,
                            std::array<Color4, 14> &out) {
  const auto views = sweep_views();
  for (size_t index = 0; index < views.size(); ++index)
    out[index] = program.evaluate(views[index], ctx);
}

inline void expect_refusal(In::ChainProgram &program,
                           std::span<const In::ChainEntryRequest> request,
                           In::ChainStatus expected, int16_t expected_index,
                           const In::FrameContext &ctx,
                           const std::array<Color4, 14> &baseline) {
  const In::ChainRefusal refusal = program.compile(request);
  HS_EXPECT_EQ(static_cast<int>(refusal.code), static_cast<int>(expected));
  HS_EXPECT_EQ(refusal.entry_index, expected_index);
  HS_EXPECT_EQ(program.ops().size(), 4u);
  std::array<Color4, 14> after;
  snapshot_render(program, ctx, after);
  for (size_t index = 0; index < after.size(); ++index)
    HS_EXPECT_TRUE(color4_identical(after[index], baseline[index]));
}

inline void test_shader_chain_refusal_shape() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 2, ValueSet::MAXIMUMS);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  std::array<Color4, 14> baseline;
  snapshot_render(program, ctx, baseline);

  expect_refusal(program, {}, In::ChainStatus::EMPTY, -1, ctx, baseline);

  std::array<In::ChainEntryRequest, In::MAX_CHAIN_OPS + 1> too_long;
  std::array<std::string, In::MAX_CHAIN_OPS + 1> labels;
  for (size_t index = 0; index < too_long.size(); ++index) {
    labels[index] = "cam" + std::to_string(index);
    too_long[index] = {labels[index], "sphere.rotate.v2"};
  }
  expect_refusal(program, too_long, In::ChainStatus::TOO_LONG, -1, ctx,
                 baseline);

  const In::ChainEntryRequest unknown[] = {
      {"camera", "sphere.rotate.v2"},
      {"warp", "warp.unknown.v2"},
  };
  expect_refusal(program, unknown, In::ChainStatus::UNKNOWN_OPERATOR, 1, ctx,
                 baseline);

  const In::ChainEntryRequest duplicate[] = {
      {"camera", "sphere.rotate.v2"},
      {"camera", "project.stereographic.v2"},
  };
  expect_refusal(program, duplicate, In::ChainStatus::DUPLICATE_INSTANCE, 1,
                 ctx, baseline);

  const In::ChainEntryRequest malformed[] = {
      {"camera", "sphere.rotate.v2"},
      {"Bad.Label", "project.stereographic.v2"},
  };
  expect_refusal(program, malformed, In::ChainStatus::MALFORMED_INSTANCE, 1,
                 ctx, baseline);

  const In::ChainEntryRequest bad_entry[] = {
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  expect_refusal(program, bad_entry, In::ChainStatus::ENTRY_FAMILY, 0, ctx,
                 baseline);

  const In::ChainEntryRequest bad_exit[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
  };
  expect_refusal(program, bad_exit, In::ChainStatus::EXIT_FAMILY, 2, ctx,
                 baseline);

  const In::ChainEntryRequest mismatch[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"camera2", "sphere.rotate.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  expect_refusal(program, mismatch, In::ChainStatus::CARRIER_MISMATCH, 2, ctx,
                 baseline);

  // State continuity: the shape refusals above left every clock untouched.
  const auto &source = state_as<In::Op::SourceClockState>(program, 2);
  const float speed = param_as<In::Op::GridSampleParams>(program, 2).speed;
  HS_EXPECT_GT(speed, 0.0f);
  HS_EXPECT_EQ(source.primary,
               fmodf(fmodf(speed, math::TWO_PI_F) + speed, math::TWO_PI_F));
  program.clear();
}

inline void test_shader_chain_refusal_budget_overflows() {
  auto oversized_table = In::OPERATOR_TABLE;
  oversized_table[0].runtime.param.size = 0xfffffffcu;
  auto oversized =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, oversized_table);
  HS_EXPECT_EQ(oversized->program.compile(DEFAULT_CHAIN).code,
               In::ChainStatus::ARENA_OVERFLOW);
  HS_EXPECT_FALSE(oversized->program.compiled());

  // Exact-fit boundary: capacity == used commits, capacity - 1 refuses.
  auto measured = std::make_unique<ProgramFixture>();
  arm_default_chain(measured->program, 0, ValueSet::DEFAULTS);
  const size_t needed = measured->program.used_bytes();
  measured->program.clear();

  auto exact = std::make_unique<ProgramFixture>(needed);
  const In::ChainRefusal fits = exact->program.compile(
      std::span<const In::ChainEntryRequest>(DEFAULT_CHAIN));
  HS_EXPECT_EQ(static_cast<int>(fits.code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(exact->program.used_bytes(), needed);

  auto small = std::make_unique<ProgramFixture>(needed - 1);
  const In::ChainRefusal overflow = small->program.compile(
      std::span<const In::ChainEntryRequest>(DEFAULT_CHAIN));
  HS_EXPECT_EQ(static_cast<int>(overflow.code),
               static_cast<int>(In::ChainStatus::ARENA_OVERFLOW));
  HS_EXPECT_EQ(overflow.entry_index, -1);
  HS_EXPECT_FALSE(small->program.compiled());

  // A committed program survives a later over-budget recompile.
  const In::ChainEntryRequest longer[] = {
      {"camera", "sphere.rotate.v2"},
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal rejected = exact->program.compile(longer);
  HS_EXPECT_EQ(static_cast<int>(rejected.code),
               static_cast<int>(In::ChainStatus::ARENA_OVERFLOW));
  HS_EXPECT_EQ(exact->program.ops().size(), 4u);
  const In::FrameContext ctx = shared_resources().context();
  exact->program.prepare(ctx);
  const Color4 color =
      exact->program.evaluate(math::Vector(1, 1, 1).normalized(), ctx);
  HS_EXPECT_TRUE(std::isfinite(color.alpha));
  exact->program.clear();

  // Schema-field budget: one fat operator overflows MAX_CHAIN_PARAMS.
  CountLifecycle::reset();
  auto fat =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  const In::ChainEntryRequest fat_chain[] = {
      {"fat", "test.fat.v2"},
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal params_overflow = fat->program.compile(fat_chain);
  HS_EXPECT_EQ(static_cast<int>(params_overflow.code),
               static_cast<int>(In::ChainStatus::PARAM_OVERFLOW));
  HS_EXPECT_EQ(params_overflow.entry_index, -1);
  HS_EXPECT_FALSE(fat->program.compiled());
  // Refused before layout: no lifecycle callback ran.
  HS_EXPECT_EQ(CountLifecycle::inits, 0);

  // Both caps admit an exact fit: MAX_CHAIN_OPS entries and MAX_CHAIN_PARAMS
  // schema fields compile.
  auto longest = std::make_unique<ProgramFixture>();
  std::array<In::ChainEntryRequest, In::MAX_CHAIN_OPS> at_cap;
  std::array<std::string, In::MAX_CHAIN_OPS> at_cap_labels;
  constexpr size_t TAIL = std::size(DEFAULT_CHAIN) - 1;
  for (size_t index = 0; index + TAIL < at_cap.size(); ++index) {
    at_cap_labels[index] = "cam" + std::to_string(index);
    at_cap[index] = {at_cap_labels[index], "sphere.rotate.v2"};
  }
  for (size_t index = 0; index < TAIL; ++index)
    at_cap[at_cap.size() - TAIL + index] = DEFAULT_CHAIN[1 + index];
  const In::ChainRefusal ops_fit = longest->program.compile(at_cap);
  HS_EXPECT_EQ(static_cast<int>(ops_fit.code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(longest->program.ops().size(), In::MAX_CHAIN_OPS);
  longest->program.clear();

  auto widest =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  const In::ChainEntryRequest exact_fit_chain[] = {
      {"filler", "test.exact-fit.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  size_t field_total = 0;
  for (const In::ChainEntryRequest &entry : exact_fit_chain)
    for (const In::OperatorDescriptor &op : extended_table())
      if (entry.operator_id == op.operator_id)
        field_total += op.schema_count;
  HS_EXPECT_EQ(field_total, In::MAX_CHAIN_PARAMS);
  const In::ChainRefusal params_fit = widest->program.compile(exact_fit_chain);
  HS_EXPECT_EQ(static_cast<int>(params_fit.code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_TRUE(widest->program.compiled());
  widest->program.clear();
}

inline void test_shader_chain_refusal_migrate_failed() {
  CountLifecycle::reset();
  auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  In::ChainProgram &program = fixture->program;
  const In::ChainEntryRequest chain[] = {
      {"counter", "test.count-a.v2"},
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal first = program.compile(chain);
  HS_EXPECT_EQ(static_cast<int>(first.code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(CountLifecycle::inits, 1);
  program.advance();
  program.advance();
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  std::array<Color4, 14> baseline;
  snapshot_render(program, ctx, baseline);

  CountLifecycle::fail_migrate = true;
  const int destroys_before = CountLifecycle::destroys;
  const In::ChainRefusal failed = program.compile(chain);
  HS_EXPECT_EQ(static_cast<int>(failed.code),
               static_cast<int>(In::ChainStatus::MIGRATE_FAILED));
  HS_EXPECT_EQ(failed.entry_index, 0);
  // The candidate tore down as a unit and the live program is untouched: the
  // failed dst was rolled back unconstructed, no live state was destroyed,
  // and the accumulated phase survives.
  HS_EXPECT_EQ(CountLifecycle::destroys, destroys_before);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  std::array<Color4, 14> after;
  snapshot_render(program, ctx, after);
  for (size_t index = 0; index < after.size(); ++index)
    HS_EXPECT_TRUE(color4_identical(after[index], baseline[index]));

  // A failing migrate deeper in the chain destroys the candidate states
  // constructed before it — and only those. "fresh" is a new pair (init at
  // entry 0), "counter" a surviving pair whose migrate fails at entry 1.
  const In::ChainEntryRequest deep_chain[] = {
      {"fresh", "test.count-b.v2"},
      {"counter", "test.count-a.v2"},
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const int inits_before = CountLifecycle::inits;
  const int migrates_before = CountLifecycle::migrates;
  const int deep_destroys_before = CountLifecycle::destroys;
  const In::ChainRefusal deep = program.compile(deep_chain);
  HS_EXPECT_EQ(static_cast<int>(deep.code),
               static_cast<int>(In::ChainStatus::MIGRATE_FAILED));
  HS_EXPECT_EQ(deep.entry_index, 1);
  // The fresh candidate was constructed, then torn down with the candidate
  // arena; the live program keeps its state.
  HS_EXPECT_EQ(CountLifecycle::inits, inits_before + 1);
  HS_EXPECT_EQ(CountLifecycle::migrates, migrates_before + 1);
  HS_EXPECT_EQ(CountLifecycle::destroys, deep_destroys_before + 1);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  std::array<Color4, 14> after_deep;
  snapshot_render(program, ctx, after_deep);
  for (size_t index = 0; index < after_deep.size(); ++index)
    HS_EXPECT_TRUE(color4_identical(after_deep[index], baseline[index]));
  CountLifecycle::fail_migrate = false;
  program.clear();
}

inline void test_shader_chain_state_identity_migration() {
  CountLifecycle::reset();
  auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  In::ChainProgram &program = fixture->program;
  const In::ChainEntryRequest first[] = {
      {"a", "test.count-a.v2"},
      {"b", "test.count-a.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(static_cast<int>(program.compile(first).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(CountLifecycle::inits, 2);
  program.advance();
  program.advance();
  program.advance();

  // Surviving pair migrates and keeps its phase; the removed instance is
  // destroyed; the new label gets a fresh init.
  const In::ChainEntryRequest second[] = {
      {"a", "test.count-a.v2"},
      {"c", "test.count-a.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const int destroys_before = CountLifecycle::destroys;
  HS_EXPECT_EQ(static_cast<int>(program.compile(second).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(CountLifecycle::migrates, 1);
  HS_EXPECT_EQ(CountLifecycle::inits, 3);
  // Commit destroyed the loser arena: old "a" and old "b".
  HS_EXPECT_EQ(CountLifecycle::destroys, destroys_before + 2);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 3.0f);
  HS_EXPECT_EQ(state_as<CountingState>(program, 1).accumulator, 0.0f);

  // Same label, different operator: a fresh init, never a migration.
  const In::ChainEntryRequest third[] = {
      {"a", "test.count-b.v2"},
      {"c", "test.count-a.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const int migrates_before = CountLifecycle::migrates;
  HS_EXPECT_EQ(static_cast<int>(program.compile(third).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 0.0f);
  // "a" re-inits under the new operator; "c" migrates.
  HS_EXPECT_EQ(CountLifecycle::inits, 4);
  HS_EXPECT_EQ(CountLifecycle::migrates, migrates_before + 1);
  program.clear();
  // Every live construction (init or successful migrate clone) was destroyed.
  HS_EXPECT_EQ(CountLifecycle::destroys,
               CountLifecycle::inits + CountLifecycle::migrates);
}

inline void test_shader_chain_state_continuity_slice() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 6, ValueSet::MAXIMUMS);
  const auto camera_before = state_as<In::Op::SpatialWalkState>(program, 0);
  const auto source_before = state_as<In::Op::SourceClockState>(program, 2);
  HS_EXPECT_GT(source_before.primary, 0.0f);

  const In::ChainEntryRequest edited[] = {
      {"camera", "sphere.rotate.v2"},
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(static_cast<int>(program.compile(edited).code),
               static_cast<int>(In::ChainStatus::OK));
  // An unchanged pair keeps its accumulated walk and clocks across an edit
  // elsewhere; the new instance starts fresh.
  const auto &camera_after = state_as<In::Op::SpatialWalkState>(program, 0);
  HS_EXPECT_EQ(std::memcmp(&camera_after.wander, &camera_before.wander,
                           sizeof(math::Quaternion)),
               0);
  HS_EXPECT_EQ(camera_after.spin_phase, camera_before.spin_phase);
  HS_EXPECT_EQ(camera_after.walk_time, camera_before.walk_time);
  const auto &fresh = state_as<In::Op::SpatialWalkState>(program, 1);
  HS_EXPECT_EQ(fresh.spin_phase, 0.0f);
  HS_EXPECT_EQ(fresh.walk_time, uint32_t{0});
  const auto &source_after = state_as<In::Op::SourceClockState>(program, 3);
  HS_EXPECT_EQ(source_after.primary, source_before.primary);
  HS_EXPECT_EQ(source_after.secondary, source_before.secondary);
  HS_EXPECT_EQ(source_after.angle, source_before.angle);
  program.clear();
}

inline void test_shader_chain_operator_state_migration() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest first[] = {
      {"ripple", "sphere.displace.ripple.v2"},
      {"project", "project.stereographic.v2"},
      {"wave", "warp.wave-shear.v2"},
      {"affine", "warp.affine.v3"},
      {"sample", "sample.projected-noise.v2"},
      {"colorize", "colorize.generated-palette.v3"}};
  const auto status = program.compile(first).code;
  HS_EXPECT_EQ(static_cast<int>(status), static_cast<int>(In::ChainStatus::OK));
  if (status != In::ChainStatus::OK)
    return;
  const_cast<In::Op::RipplePhaseState &>(
      state_as<In::Op::RipplePhaseState>(program, 0))
      .phase = 0.37f;
  const_cast<In::Op::WarpPhaseState &>(
      state_as<In::Op::WarpPhaseState>(program, 2))
      .phase = 0.59f;
  auto &affine = const_cast<In::Op::AffineClockState &>(
      state_as<In::Op::AffineClockState>(program, 3));
  affine.phase = 0.73f;
  affine.rotation = 1.19f;
  auto &noise = const_cast<In::Op::NoisePhaseState &>(
      state_as<In::Op::NoisePhaseState>(program, 4));
  noise.phase = 1.47f;
  noise.noise.SetSeed(923);
  const float sample = noise.noise.GetNoise(0.3f, 0.7f, 1.2f);
  const In::ChainEntryRequest edited[] = {{"camera", "sphere.rotate.v2"},
                                          first[0],
                                          first[1],
                                          first[2],
                                          first[3],
                                          first[4],
                                          first[5]};
  const auto edited_status = program.compile(edited).code;
  HS_EXPECT_EQ(static_cast<int>(edited_status),
               static_cast<int>(In::ChainStatus::OK));
  if (edited_status != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(state_as<In::Op::RipplePhaseState>(program, 1).phase, 0.37f);
  HS_EXPECT_EQ(state_as<In::Op::WarpPhaseState>(program, 3).phase, 0.59f);
  HS_EXPECT_EQ(state_as<In::Op::AffineClockState>(program, 4).phase, 0.73f);
  HS_EXPECT_EQ(state_as<In::Op::AffineClockState>(program, 4).rotation, 1.19f);
  HS_EXPECT_EQ(state_as<In::Op::NoisePhaseState>(program, 5).phase, 1.47f);
  HS_EXPECT_EQ(state_as<In::Op::NoisePhaseState>(program, 5)
                   .noise.GetNoise(0.3f, 0.7f, 1.2f),
               sample);
  const In::ChainEntryRequest rings[] = {
      {"rings", "sample.spherical-rings.v3"},
      {"colorize", "colorize.generated-palette.v3"}};
  const auto rings_status = program.compile(rings).code;
  HS_EXPECT_EQ(static_cast<int>(rings_status),
               static_cast<int>(In::ChainStatus::OK));
  if (rings_status != In::ChainStatus::OK)
    return;
  auto &ring = const_cast<In::Op::SphericalRingsState &>(
      state_as<In::Op::SphericalRingsState>(program, 0));
  ring.phase = 0.41f;
  ring.walk.spin_phase = 0.83f;
  ring.walk.walk_time = 167;
  const In::ChainEntryRequest edited_rings[] = {
      {"camera", "sphere.rotate.v2"}, rings[0], rings[1]};
  const auto final_status = program.compile(edited_rings).code;
  HS_EXPECT_EQ(static_cast<int>(final_status),
               static_cast<int>(In::ChainStatus::OK));
  if (final_status != In::ChainStatus::OK)
    return;
  const auto &after = state_as<In::Op::SphericalRingsState>(program, 1);
  HS_EXPECT_EQ(after.phase, 0.41f);
  HS_EXPECT_EQ(after.walk.spin_phase, 0.83f);
  HS_EXPECT_EQ(after.walk.walk_time, uint32_t{167});
}

inline void test_shader_chain_determinism() {
  auto first = std::make_unique<ProgramFixture>();
  auto second = std::make_unique<ProgramFixture>();
  arm_default_chain(first->program, 5, ValueSet::MAXIMUMS);
  arm_default_chain(second->program, 5, ValueSet::MAXIMUMS);
  FastNoiseLite authored_walk_noise;
  authored_walk_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  authored_walk_noise.SetSeed(
      static_cast<int32_t>(In::instance_hash("camera", "sphere.rotate.v2")));
  authored_walk_noise.SetFrequency(In::Op::WALK_OPTIONS.noise_scale);
  const auto &seeded_walk =
      state_as<In::Op::SpatialWalkState>(first->program, 0);
  HS_EXPECT_EQ(seeded_walk.walk_noise.GetNoise(0.25f, -0.5f, 0.75f),
               authored_walk_noise.GetNoise(0.25f, -0.5f, 0.75f));
  const In::FrameContext ctx = shared_resources().context();
  first->program.prepare(ctx);
  second->program.prepare(ctx);
  for (const math::Vector &view : sweep_views()) {
    const Color4 a = first->program.evaluate(view, ctx);
    const Color4 b = second->program.evaluate(view, ctx);
    HS_EXPECT_TRUE(color4_identical(a, b));
  }
  // Instance labels participate in stateful-resource identity.
  auto relabeled = std::make_unique<ProgramFixture>();
  const In::ChainEntryRequest renamed[] = {
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(static_cast<int>(relabeled->program.compile(renamed).code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(relabeled->program, 0),
                  ValueSet::MAXIMUMS);
  for (int frame = 0; frame < 5; ++frame)
    relabeled->program.advance();
  const auto &walk_a = state_as<In::Op::SpatialWalkState>(first->program, 0);
  const auto &walk_b =
      state_as<In::Op::SpatialWalkState>(relabeled->program, 0);
  HS_EXPECT_NE(
      std::memcmp(&walk_a.wander, &walk_b.wander, sizeof(math::Quaternion)), 0);
  first->program.clear();
  second->program.clear();
  relabeled->program.clear();
}

inline void test_shader_chain_param_names_and_budget() {
  static_assert(In::PER_PARAM_NAME_BYTES ==
                In::MAX_INSTANCE_ID + 1 + In::MAX_FIELD_ID + 1);
  static_assert(In::operator_schema_ids_fit_names());
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 0, ValueSet::DEFAULTS);

  // The names are exactly "{instance}.{field-id}", per schema entry.
  const auto ops = program.ops();
  for (size_t index = 0; index < ops.size(); ++index) {
    const In::OperatorDescriptor &op = *ops[index].op;
    for (uint16_t field = 0; field < op.schema_count; ++field) {
      std::string expected{ops[index].instance};
      expected += '.';
      expected += op.schema[field].id;
      HS_EXPECT_TRUE(expected == program.param_name(index, field));
    }
  }
  HS_EXPECT_TRUE(std::string_view(program.param_name(0, 0)) == "camera.wander");
  HS_EXPECT_TRUE(std::string_view(program.param_name(3, 0)) ==
                 "colorize.hue-shift-amount");

  // The committed footprint equals the accounted layout: cataloged blocks,
  // the per-op overhead, and the fixed per-field name reservation.
  size_t expected_bytes = 0;
  const auto align_to = [](size_t offset, size_t alignment) {
    return (offset + alignment - 1) & ~(alignment - 1);
  };
  for (const In::ChainProgram::ChainOp &op : ops) {
    const In::OperatorRuntime &runtime = op.op->runtime;
    expected_bytes = align_to(expected_bytes, runtime.param.align);
    expected_bytes += runtime.param.size;
    expected_bytes = align_to(expected_bytes, runtime.prepared.align);
    expected_bytes += runtime.prepared.size;
    expected_bytes = align_to(expected_bytes, runtime.state.align);
    expected_bytes += runtime.state.size;
    expected_bytes += In::PER_OP_OVERHEAD_BYTES;
    expected_bytes +=
        static_cast<size_t>(op.op->schema_count) * In::PER_PARAM_NAME_BYTES;
  }
  HS_EXPECT_EQ(program.used_bytes(), expected_bytes);
  program.clear();
}

/** Reaches the effect's committed program, palette state, and hue bakes. */
struct ShaderChainWhiteBox {
  using FX = ShaderChain<96, 20>;

  static const In::ChainProgram &program(const FX &effect) {
    return effect.program;
  }
  static In::ChainProgram &program(FX &effect) { return effect.program; }
  static void advance_without_render(FX &effect) { effect.advance_clocks(); }
  static Pixel palette_color(const FX &effect, float value) {
    return effect.generated_palettes.palette(In::Op::PaletteMode::TRIADIC)
        .get(value)
        .color;
  }

  /** The committed colorize instance's parameter block. */
  static In::Op::GeneratedPaletteParams &color_params(FX &effect) {
    return *reinterpret_cast<In::Op::GeneratedPaletteParams *>(
        effect.program.param_block(static_cast<size_t>(effect.colorize.index)));
  }
  /** The committed colorize instance's phase clocks. */
  static const In::Op::ColorClockState &color_clocks(const FX &effect) {
    return *static_cast<const In::Op::ColorClockState *>(
        effect.program.state_block(static_cast<size_t>(effect.colorize.index)));
  }
  /** The frame snapshot draw_frame() hands the program, hue bakes included. */
  static In::FrameContext frame_context(FX &effect) {
    return effect.make_frame_context(effect.colorize);
  }
  /** The resident hue-noise table; writable so a probe can poison it and see
      whether the next frame re-baked. */
  static int8_t *hue_noise_lut(FX &effect) {
    return effect.resources->hue_noise_lut.data();
  }
  /** The inputs the resident hue-noise table was baked from. */
  static float baked_noise_scale(const FX &effect) {
    return effect.resources->hue_noise_bake.scale;
  }
  static float baked_noise_phase(const FX &effect) {
    return effect.resources->hue_noise_bake.phase;
  }
};

inline void test_shader_chain_composed_frame_parity() {
  using FX = AlienCore<96, 20>;
  reset_globals();
  auto params = FX::initial_params();
  params.template get<"projection">().camera_wander = 0.0f;
  std::array<FX::FrameState, 5> references;
  {
    FX composed;
    composed.init();
    ComposedFrameWhiteBox::set_params(composed, params);
    for (size_t frame = 0; frame < references.size(); ++frame) {
      pin_frame_clock(static_cast<int>(frame) + 1);
      ComposedFrameWhiteBox::advance(composed);
      references[frame] = ComposedFrameWhiteBox::frame(composed);
    }
  }

  reset_globals();
  ShaderChainWhiteBox::FX chain;
  chain.init();
  const In::ChainEntryRequest topology[] = {
      {"camera", "sphere.rotate.v2"},
      {"lens", "sphere.lens.glitch.v2"},
      {"project", "project.gnomonic.v2"},
      {"warp", "warp.mirror-tile.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(static_cast<int>(chain.set_chain(topology).code),
               static_cast<int>(In::ChainStatus::OK));
  auto &program = ShaderChainWhiteBox::program(chain);
  param_as<In::Op::RotateChainParams>(program, 0).wander =
      params.template get<"projection">().camera_wander;
  auto &projection = param_as<In::Op::GnomonicChainParams>(program, 2);
  projection.frame = static_cast<uint8_t>(In::Op::ProjectionFrame::IDENTITY);
  projection.singularity_fade =
      params.template get<"projection">().singularity_fade;
  static_cast<PB::MirrorParams &>(param_as<In::Op::MirrorWarpParams>(
      program, 3)) = params.template get<"outer_warp">();
  auto &source = param_as<In::Op::GridSampleParams>(program, 4);
  static_cast<PB::GridSourceParams &>(source) = params.template get<"source">();
  source.edge_width = params.template get<"value">().edge_width;
  source.coverage_mode =
      static_cast<uint8_t>(PB::ProjectionCoverageMode::EDGE_FADE);
  auto &color = ShaderChainWhiteBox::color_params(chain);
  static_cast<PB::Color::ColorControls &>(color) =
      params.template get<"color">();
  color.mapping_mode =
      static_cast<uint8_t>(params.template get<"color">().palette_mapping);

  size_t visible = 0;
  for (int frame = 1; frame <= 5; ++frame) {
    HS_CONTEXT("frame", frame);
    pin_frame_clock(frame);
    chain.draw_frame();
    chain.advance_display();
    const In::FrameContext ctx = ShaderChainWhiteBox::frame_context(chain);
    auto reference = references[static_cast<size_t>(frame - 1)];
    reference.palette = ctx.palettes[color.palette_mode];
    reference.hue_rotation_lut = ctx.hue_rotation_lut;
    reference.hue_noise_lut = ctx.hue_noise_lut;
    const auto prepared = FX::RenderPipeline::prepare(reference);
    for (const math::Vector &view : sweep_views()) {
      const Color4 expected = FX::shade(view, prepared);
      const Color4 actual = program.evaluate(view, ctx);
      HS_EXPECT_NEAR(actual.color.r, expected.color.r, 1);
      HS_EXPECT_NEAR(actual.color.g, expected.color.g, 1);
      HS_EXPECT_NEAR(actual.color.b, expected.color.b, 1);
      HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
      visible += expected.alpha > 0.0f &&
                 (expected.color.r != 0 || expected.color.g != 0 ||
                  expected.color.b != 0);
    }
  }
  HS_EXPECT_GT(visible, 0u);
}

inline void test_shader_chain_effect_registers_params() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  size_t expected = 0;
  for (const In::ChainEntryRequest &entry : DEFAULT_CHAIN)
    expected += In::find_operator(entry.operator_id)->schema_count;
  const ParamList &params = effect.getParameters();
  HS_EXPECT_EQ(params.size(), expected);
  HS_EXPECT_TRUE(params.find("camera.wander") != nullptr);
  HS_EXPECT_TRUE(params.find("project.singularity-fade") != nullptr);
  HS_EXPECT_TRUE(params.find("sample.pattern-freq") != nullptr);
  HS_EXPECT_TRUE(params.find("colorize.palette-chroma") != nullptr);
  const ParamDef *coverage = params.find("sample.coverage-mode");
  HS_EXPECT_TRUE(coverage != nullptr);
  if (!coverage)
    return;
  HS_EXPECT_TRUE(coverage->is_enum());
  HS_EXPECT_EQ(coverage->option_count, 4);
  HS_EXPECT_TRUE(std::string_view(coverage->options[3]) == "edge-fade");
  HS_EXPECT_EQ(
      static_cast<int>(effect.updateParameter("sample.coverage-mode", 3.0f)),
      static_cast<int>(ParamSetResult::APPLIED));
  HS_EXPECT_EQ(param_as<In::Op::GridSampleParams>(
                   ShaderChainWhiteBox::program(effect), 2)
                   .coverage_mode,
               static_cast<uint8_t>(In::Op::ProjectionCoverageMode::EDGE_FADE));

  uint64_t lit = 0;
  for (int frame = 0; frame < 4; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
}

inline void test_shader_chain_parameter_admission() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"lens", "sphere.lens.mobius.v2"},
      {"project", "project.stereographic.v2"},
      {"warp", "warp.curl-flow.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(static_cast<int>(effect.set_chain(chain).code),
               static_cast<int>(In::ChainStatus::OK));
  const ParamDef *strength_param = effect.getParameters().find("warp.strength");
  HS_EXPECT_TRUE(strength_param != nullptr);
  if (!strength_param)
    return;
  const ParamDef &strength = *strength_param;
  HS_EXPECT_EQ(effect.updateParameter("warp.strength", 30.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.parameter_warning("warp.strength") == nullptr);
  HS_EXPECT_EQ(strength.get_requested(), strength.max);
  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_EQ(effect.accepted_parameter_value(strength), strength.max);
  effect.draw_frame();
  effect.advance_display();
  effect.updateParameter("warp.strength", 0.0f);
  HS_EXPECT_TRUE(effect.parameter_warning("warp.strength") == nullptr);
  HS_EXPECT_EQ(effect.accepted_parameter_value(strength), 0.0f);
  const ParamDef *d_re_param = effect.getParameters().find("lens.mobius-d-re");
  HS_EXPECT_TRUE(d_re_param != nullptr);
  if (!d_re_param)
    return;
  const ParamDef &d_re = *d_re_param;
  const float accepted_d_re = effect.accepted_parameter_value(d_re);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-a-re", 2.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-d-re", 0.0f),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-d-re") != nullptr);
  HS_EXPECT_EQ(d_re.get_requested(), accepted_d_re);
  HS_EXPECT_EQ(effect.accepted_parameter_value(d_re), accepted_d_re);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-a-re", 3.0f),
               ParamSetResult::APPLIED);
  const float saved_a =
      effect.getParameters().find("lens.mobius-a-re")->get_requested();
  const float saved_d = d_re.get_requested();
  HS_EXPECT_EQ(static_cast<int>(effect.set_chain(chain).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-a-re", saved_a),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-d-re", saved_d),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.getParameters().find("lens.mobius-a-re")->get_requested(),
               saved_a);
  HS_EXPECT_EQ(effect.getParameters().find("lens.mobius-d-re")->get_requested(),
               saved_d);
  effect.draw_frame();
  effect.advance_display();
  effect.updateParameter("lens.mobius-d-re", 1.0f);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-d-re") == nullptr);
  const ShaderChainParameterWrite valid_batch[] = {{"lens.mobius-a-re", 0.0f},
                                                   {"lens.mobius-b-re", 1.0f},
                                                   {"lens.mobius-c-re", 1.0f},
                                                   {"lens.mobius-d-re", 0.0f}};
  HS_EXPECT_EQ(effect.update_parameters(valid_batch), ParamSetResult::APPLIED);
  for (const auto &write : valid_batch)
    HS_EXPECT_EQ(effect.getParameters().find(write.name)->get_requested(),
                 write.value);
  const ShaderChainParameterWrite invalid_batch[] = {
      {"lens.mobius-b-re", 0.0f}, {"lens.mobius-c-re", 0.0f}};
  HS_EXPECT_EQ(effect.update_parameters(invalid_batch),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-b-re") != nullptr);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-a-re") == nullptr);
  const ShaderChainParameterWrite unknown_batch[] = {{"lens.mobius-a-re", 1.0f},
                                                     {"missing", 0.0f}};
  HS_EXPECT_EQ(effect.update_parameters(unknown_batch),
               ParamSetResult::UNKNOWN_PARAM);
  const ShaderChainParameterWrite nonfinite_batch[] = {
      {"lens.mobius-a-re", 1.0f}, {"lens.mobius-c-re", NAN}};
  HS_EXPECT_EQ(effect.update_parameters(nonfinite_batch),
               ParamSetResult::NON_FINITE);
  for (const auto &write : valid_batch)
    HS_EXPECT_EQ(effect.getParameters().find(write.name)->get_requested(),
                 write.value);
  effect.draw_frame();
  effect.advance_display();
}

inline void test_shader_chain_edge_distance_admission() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"project", "project.folded-sinusoidal.v2"},
      {"warp", "warp.wave-shear.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  HS_EXPECT_EQ(effect.updateParameter("sample.coverage-mode", 3),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning("sample.coverage-mode") != nullptr);
  const ParamDef *envelope = effect.getParameters().find("warp.envelope");
  const ParamDef *coverage_mode =
      effect.getParameters().find("sample.coverage-mode");
  HS_EXPECT_TRUE(envelope != nullptr && coverage_mode != nullptr);
  if (!envelope || !coverage_mode)
    return;
  const float ACCEPTED_ENVELOPE = envelope->get_requested();
  const ShaderChainParameterWrite writes[] = {{"warp.envelope", 2}};
  HS_EXPECT_EQ(effect.update_parameters(writes), ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.getParameters().find("warp.envelope")->get_requested(),
               ACCEPTED_ENVELOPE);
  HS_EXPECT_TRUE(effect.parameter_warning("warp.envelope") != nullptr);
  HS_EXPECT_EQ(
      effect.getParameters().find("sample.coverage-mode")->get_requested(),
      1.0f);
}

inline void test_shader_chain_effect_rebind_generation() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const uint32_t initial = effect.getParameterSchemaGeneration();
  const In::ChainEntryRequest edited[] = {
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal committed = effect.set_chain(edited);
  HS_EXPECT_EQ(static_cast<int>(committed.code),
               static_cast<int>(In::ChainStatus::OK));
  // set_chain rebinds before returning: fresh generation, fresh names.
  HS_EXPECT_NE(effect.getParameterSchemaGeneration(), initial);
  HS_EXPECT_TRUE(effect.getParameters().find("camera2.wander") != nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("camera.wander") == nullptr);
  effect.draw_frame();
  effect.advance_display();
}

inline void test_shader_chain_effect_refusal_keeps_schema() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const uint32_t committed = effect.getParameterSchemaGeneration();
  const size_t param_count = effect.getParameters().size();
  const In::ChainEntryRequest unknown[] = {
      {"camera", "sphere.rotate.v2"},
      {"warp", "warp.unknown.v2"},
  };
  const In::ChainRefusal refusal = effect.set_chain(unknown);
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::UNKNOWN_OPERATOR));
  HS_EXPECT_EQ(refusal.entry_index, 1);
  // Transactional at the effect layer too: definitions, generation, and the
  // committed program all survive, and the effect still renders.
  HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), committed);
  HS_EXPECT_EQ(effect.getParameters().size(), param_count);
  HS_EXPECT_TRUE(effect.getParameters().find("camera.wander") != nullptr);
  uint64_t lit = 0;
  for (int frame = 0; frame < 2; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
}

/** @brief Pause gates preset selection only: chain clocks and the visible
    palette keep advancing. */
inline void test_shader_chain_pause_semantics() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(static_cast<int>(effect.updateParameter("sample.speed", 0.05f)),
               static_cast<int>(ParamSetResult::APPLIED));
  HS_EXPECT_EQ(static_cast<int>(effect.updateParameter("camera.wander", 1.0f)),
               static_cast<int>(ParamSetResult::APPLIED));
  effect.setAnimationsPaused(true);
  const In::Op::SourceClockState source_before =
      state_as<In::Op::SourceClockState>(WB::program(effect), 2);
  const math::Quaternion wander_before =
      state_as<In::Op::SpatialWalkState>(WB::program(effect), 0).wander;
  const Pixel color_before = WB::palette_color(effect, 0.25f);

  for (int frame = 0; frame < 60; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }

  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_NE(
      state_as<In::Op::SourceClockState>(WB::program(effect), 2).primary,
      source_before.primary);
  HS_EXPECT_TRUE(
      state_as<In::Op::SpatialWalkState>(WB::program(effect), 0).wander !=
      wander_before);
  const Pixel color_after = WB::palette_color(effect, 0.25f);
  HS_EXPECT_TRUE(color_after.r != color_before.r ||
                 color_after.g != color_before.g ||
                 color_after.b != color_before.b);
  uint64_t lit = 0;
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
}
