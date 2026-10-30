/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

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

struct ClockedSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static Pattern::PreparedSource prepare(const FrameState &frame) {
    return Pattern::prepare(frame.primary, frame.secondary, frame.angle);
  }
};

struct TwinWaveSampleMirror : ClockedSampleMirror {
  static const In::Op::TwinWaveSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.twin_wave;
  }
};

struct RingsSampleMirror : ClockedSampleMirror {
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
  static Pattern::PreparedSphericalRings
  prepare(const SampleMirrorFrame &frame) {
    return {frame.ring_axis, frame.ring_phase};
  }
};

struct SpiralSampleMirror : ClockedSampleMirror {
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

struct FractalSampleMirror : ClockedSampleMirror {
  static const In::Op::FractalSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.fractal;
  }
};

struct TessellationSampleMirror : ClockedSampleMirror {
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

template <typename OpParams>
inline bool arm_spherical_sample_op_chain(In::ChainProgram &program,
                                          const char *op_id, int frames,
                                          ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"source", op_id},
      {"colorize", "colorize.generated-palette.v4"},
  };
  return arm_chain<In::Op::RotateChainParams, OpParams,
                   In::Op::GeneratedPaletteParams>(program, chain, frames, set);
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
  expect_op_parity<BoundStage>(
      program, ctx, 1, sample_mirror(program, 1),
      [&](const math::Vector &view) {
        return run_entry<PB::SphereSample>(program, ctx, 0,
                                           PB::SphereSample{view, 0.0f});
      },
      field_identical);
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

inline void test_shader_chain_parity_sample_rings() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::RingsSampleParams>(
            program, In::Op::SampleRings::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Rings<RingsSampleMirror>>(program, ctx);
  });
}

inline void test_shader_chain_parity_sample_spherical_rings() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_spherical_sample_op_chain<In::Op::SphericalRingsSampleParams>(
            program, In::Op::SampleSphericalRings::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    expect_spherical_sample_op_parity<
        PB::Source::SphericalRings<SphericalRingsSampleMirror>>(program, ctx);
  });
}

inline void test_shader_chain_parity_sample_spiral() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::SpiralSampleParams>(
            program, In::Op::SampleSpiral::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Spiral<SpiralSampleMirror>>(program, ctx);
  });
}

inline void test_shader_chain_parity_sample_lattice() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::LatticeSampleParams>(
            program, In::Op::SampleLattice::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::PrimitiveLattice<LatticeSampleMirror>>(
        program, ctx);
  });
}

inline void test_shader_chain_parity_sample_fractal() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::FractalSampleParams>(
            program, In::Op::SampleFractal::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::EscapeFractal<FractalSampleMirror>>(
        program, ctx);
  });
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
      {"color", In::Op::ColorizeGeneratedPaletteV4::ID},
  };
  const auto compile_status1 = program.compile(chain).code;
  HS_EXPECT_EQ(compile_status1, In::ChainStatus::OK);
  if (compile_status1 != In::ChainStatus::OK)
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

inline void test_shader_chain_large_finite_path_length() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  std::array<In::ChainEntryRequest, 15> chain{};
  std::array<std::string, 11> ids;
  chain[0] = {"project", In::Op::ProjectStereographic::ID};
  for (int index = 0; index < 11; ++index) {
    ids[index] = "warp" + std::to_string(index);
    chain[index + 1] = {ids[index], In::Op::WarpAffineV3::ID};
  }
  chain[12] = {"sample", In::Op::SampleGridV3::ID};
  chain[13] = {"color", In::Op::ColorizeGeneratedPaletteV4::ID};
  chain[14] = {"displacement", In::Op::ModulateColorV1::ID};
  const auto compile_status2 = program.compile(chain).code;
  HS_EXPECT_EQ(compile_status2, In::ChainStatus::OK);
  if (compile_status2 != In::ChainStatus::OK)
    return;
  for (int index = 1; index <= 11; ++index) {
    auto &params = param_as<In::Op::AffineWarpParams>(program, index);
    params.scale_x = params.scale_y = 1.0f / 64.0f;
    HS_EXPECT_TRUE(PB::Fields::valid(params));
  }
  auto &modulation = param_as<In::Op::ModulationParams>(program, 14).modulation;
  modulation.source = Pattern::SourceKind::DISPLACEMENT;
  modulation.evaluation = Pattern::Evaluation::DIRECT;
  modulation.mechanic.amount = 1e-20f;
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
    const PB::PlaneSample next =
        index == 0 ? run_entry<PB::PlaneSample>(program, ctx, 0, sphere)
                   : run_entry<PB::PlaneSample>(program, ctx, index, plane);
    if (index != 0)
      expected_path +=
          std::hypot(static_cast<double>(next.coords.re) - plane.coords.re,
                     static_cast<double>(next.coords.im) - plane.coords.im);
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
  HS_EXPECT_EQ(actual.color.b, 0);
}

inline void test_shader_chain_parity_sample_tessellation() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sample_op_chain<In::Op::TessellationSampleParams>(
            program, In::Op::SampleTessellation::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    auto &params = param_as<In::Op::TessellationSampleParams>(program, 2);
    params.kind = static_cast<uint8_t>(Pattern::TessellationKind::TRIANGULAR);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, Pattern::TessellationKind::TRIANGULAR>>(
        program, ctx);
    params.kind = static_cast<uint8_t>(Pattern::TessellationKind::SQUARE);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, Pattern::TessellationKind::SQUARE>>(program,
                                                                      ctx);
    params.kind = static_cast<uint8_t>(Pattern::TessellationKind::HEXAGONAL);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, Pattern::TessellationKind::HEXAGONAL>>(
        program, ctx);
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

inline void test_shader_chain_parity_sample_spherical_noise() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_spherical_sample_op_chain<In::Op::SphericalNoiseSampleParams>(
            program, In::Op::SampleSphericalNoise::ID, 4, set))
      return;
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
  });
}

inline void run_sample_variant(In::ChainProgram &program,
                               const In::FrameContext &ctx,
                               In::Op::WeightMode weight,
                               In::Op::ProjectionCoverageMode coverage) {
  auto &params = param_as<In::Op::GridSampleParams>(program, 2);
  params.weight_mode = static_cast<uint8_t>(weight);
  params.coverage_mode = static_cast<uint8_t>(coverage);
  constexpr auto ENV = PB::Color::BrightnessEnvelope::NONE;
  if (weight == In::Op::WeightMode::NONE) {
    switch (coverage) {
    case In::Op::ProjectionCoverageMode::NONE:
      expect_parity<PB::Weight::None, PB::ProjectionCoverage::None, ENV>(
          program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT:
      expect_parity<PB::Weight::None, PB::ProjectionCoverage::Weight, ENV>(
          program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT_SQUARED:
      expect_parity<PB::Weight::None, PB::ProjectionCoverage::WeightSquared,
                    ENV>(program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::EDGE_FADE:
      expect_parity<PB::Weight::None,
                    PB::ProjectionCoverage::EdgeFade<MirrorValue>, ENV>(program,
                                                                        ctx);
      break;
    }
  } else {
    switch (coverage) {
    case In::Op::ProjectionCoverageMode::NONE:
      expect_parity<PB::Weight::Projection, PB::ProjectionCoverage::None, ENV>(
          program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT:
      expect_parity<PB::Weight::Projection, PB::ProjectionCoverage::Weight,
                    ENV>(program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::WEIGHT_SQUARED:
      expect_parity<PB::Weight::Projection,
                    PB::ProjectionCoverage::WeightSquared, ENV>(program, ctx);
      break;
    case In::Op::ProjectionCoverageMode::EDGE_FADE:
      expect_parity<PB::Weight::Projection,
                    PB::ProjectionCoverage::EdgeFade<MirrorValue>, ENV>(program,
                                                                        ctx);
      break;
    }
  }
}

inline void test_shader_chain_parity_sample_variants() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_default_chain(program, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    for (const In::Op::WeightMode weight :
         {In::Op::WeightMode::NONE, In::Op::WeightMode::PROJECTION})
      for (const In::Op::ProjectionCoverageMode coverage :
           {In::Op::ProjectionCoverageMode::NONE,
            In::Op::ProjectionCoverageMode::WEIGHT,
            In::Op::ProjectionCoverageMode::WEIGHT_SQUARED,
            In::Op::ProjectionCoverageMode::EDGE_FADE}) {
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
  for (uint8_t palette_mode = 0; palette_mode < 3; ++palette_mode) {
    params.palette_mode = palette_mode;
    for (uint8_t mapping = 0; mapping < 4; ++mapping) {
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
        const float ROTATION = math::wrap_t(.25f + scalar * amount);
        Color4 expected = base.color;
        if (ROTATION != 0)
          expected.color =
              table
                  ? sample_hue_rotation_lut(palette.hue_rotation, FIELD.value,
                                            ROTATION)
                  : ColorModulation::HueRotation::apply(expected, 1, {ROTATION})
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
  WB::modulation_params(effect, "displacement");
  WB::modulation_params(effect, "modulation");
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
