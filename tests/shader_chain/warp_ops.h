/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- warp-op mirrors ------------------------------------------------------

/** Mirror frame for the planar warp batch. */
struct WarpMirrorFrame {
  const FastNoiseLite *noise = nullptr;
  In::Op::AffineWarpParams affine;
  In::Op::WaveShearWarpParams wave_shear;
  In::Op::VortexWarpParams vortex;
  In::Op::VectorNoiseWarpParams vector_noise;
  In::Op::MirrorWarpParams mirror;
  In::Op::PolarChartParams polar;
  In::Op::CurlFlowWarpParams curl;
  float phase = 0.0f;
  float rotation = 0.0f;
};

struct WarpMirrorBinding {
  using FrameState = WarpMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct AffineMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::AffineWarpParams &params(const FrameState &frame) {
    return frame.affine;
  }
  static PB::Warp::PreparedAffineSlot prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.affine, frame.phase, frame.rotation,
                             frame.affine.lattice_period);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct WaveShearMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::WaveShearWarpParams &params(const FrameState &frame) {
    return frame.wave_shear;
  }
  static PB::Warp::PreparedRotation prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.wave_shear, frame.phase);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct VectorNoiseMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::VectorNoiseWarpParams &params(const FrameState &frame) {
    return frame.vector_noise;
  }
  static PB::Warp::PreparedVectorNoiseSlot prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.vector_noise, frame.phase);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct VortexMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static PB::Warp::PreparedVortexSlot prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.vortex, frame.phase);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct MirrorTileMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::MirrorWarpParams &params(const FrameState &frame) {
    return frame.mirror;
  }
  static PB::Warp::PreparedMirrorSlot prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.mirror, frame.phase);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct PolarChartMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::PolarChartParams &params(const FrameState &frame) {
    return frame.polar;
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct CurlFlowMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::CurlFlowWarpParams &params(const FrameState &frame) {
    return frame.curl;
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static bool path_length_required(const FrameState &) { return true; }
};

/** Compiles a chain with one planar warp at entry 2 and applies the value set
    to every block. */
template <typename OpParams>
inline void arm_warp_op_chain(In::ChainProgram &program, const char *op_id,
                              int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"warp", op_id},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 1), set);
  apply_value_set(param_as<OpParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 3), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 4), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

/** Erased-vs-bound parity of the warp at entry 2. */
template <typename BoundStage>
inline void expect_warp_op_parity(In::ChainProgram &program,
                                  const In::FrameContext &ctx,
                                  const WarpMirrorFrame &mirror) {
  program.prepare(ctx);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[2].op;
  int view_index = 0;
  for (const math::Vector &view : sweep_views()) {
    HS_CONTEXT("view", view_index++);
    const PB::PlaneSample input = projected_input(program, ctx, view);
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&input, out, ctx, program.param_block(2),
                   program.prepared_block(2));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::PlaneSample *>(out));
    const PB::PlaneSample reference = BoundStage::run(input, mirror, prepared);
    HS_EXPECT_TRUE(plane_identical(erased, reference));
  }
}

inline WarpMirrorFrame warp_mirror(In::ChainProgram &program) {
  WarpMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[2].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::WarpAffineV3::ID) {
    const auto &state = state_as<In::Op::AffineClockState>(program, 2);
    mirror.phase = state.phase;
    mirror.rotation = state.rotation;
    mirror.affine = param_as<In::Op::AffineWarpParams>(program, 2);
  } else if (id == In::Op::WarpVectorNoise::ID ||
             id == In::Op::WarpCurlFlow::ID) {
    const auto &state = state_as<In::Op::NoisePhaseState>(program, 2);
    mirror.noise = &state.noise;
    mirror.phase = state.phase;
    if (id == In::Op::WarpVectorNoise::ID)
      mirror.vector_noise = param_as<In::Op::VectorNoiseWarpParams>(program, 2);
    else
      mirror.curl = param_as<In::Op::CurlFlowWarpParams>(program, 2);
  } else {
    mirror.phase = state_as<In::Op::WarpPhaseState>(program, 2).phase;
    if (id == In::Op::WarpWaveShear::ID)
      mirror.wave_shear = param_as<In::Op::WaveShearWarpParams>(program, 2);
    else if (id == In::Op::WarpVortex::ID)
      mirror.vortex = param_as<In::Op::VortexWarpParams>(program, 2);
    else if (id == In::Op::WarpMirrorTile::ID)
      mirror.mirror = param_as<In::Op::MirrorWarpParams>(program, 2);
    else if (id == In::Op::WarpPolarChart::ID)
      mirror.polar = param_as<In::Op::PolarChartParams>(program, 2);
  }
  return mirror;
}

inline void test_shader_chain_parity_warp_affine_mirror() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    const In::FrameContext ctx = shared_resources().context();
    arm_warp_op_chain<In::Op::AffineWarpParams>(
        program, In::Op::WarpAffineV3::ID, 4, set);
    using BoundAffine =
        typename PB::Stage::Warp<PB::Warp::AffineFrame<AffineMirrorProvider>>::
            template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<BoundAffine>(program, ctx, warp_mirror(program));
    program.clear();

    arm_warp_op_chain<In::Op::MirrorWarpParams>(
        program, In::Op::WarpMirrorTile::ID, 4, set);
    using BoundMirror = typename PB::Stage::Warp<PB::Warp::MirrorTile<
        MirrorTileMirrorProvider>>::template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<BoundMirror>(program, ctx, warp_mirror(program));
    program.clear();
  }
}

template <typename Envelope>
inline void run_wave_shear_variant(In::ChainProgram &program,
                                   const In::FrameContext &ctx,
                                   In::Op::WarpEnvelope envelope) {
  param_as<In::Op::WaveShearWarpParams>(program, 2).envelope =
      static_cast<uint8_t>(envelope);
  using Bound = typename PB::Stage::Warp<PB::Warp::WaveShear<
      WaveShearMirrorProvider, Envelope>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

inline void test_shader_chain_parity_warp_wave_shear() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::WaveShearWarpParams>(
        program, In::Op::WarpWaveShear::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_wave_shear_variant<PB::Warp::FlatEnvelope>(program, ctx,
                                                   In::Op::WarpEnvelope::FLAT);
    run_wave_shear_variant<PB::Warp::ProjectionWeightEnvelope>(
        program, ctx, In::Op::WarpEnvelope::PROJECTION_WEIGHT);
    run_wave_shear_variant<PB::Warp::EdgeFadeEnvelope>(
        program, ctx, In::Op::WarpEnvelope::EDGE_FADE);
    program.clear();
  }
}

inline void test_shader_chain_parity_warp_vortex() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::VortexWarpParams>(program, In::Op::WarpVortex::ID,
                                                4, set);
    const In::FrameContext ctx = shared_resources().context();
    using Bound =
        typename PB::Stage::Warp<PB::Warp::Vortex<VortexMirrorProvider>>::
            template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
    program.clear();
  }
}

template <math::NoiseBasis Basis, typename Envelope>
inline void run_vector_noise_variant(In::ChainProgram &program,
                                     const In::FrameContext &ctx,
                                     In::Op::WarpEnvelope envelope) {
  auto &params = param_as<In::Op::VectorNoiseWarpParams>(program, 2);
  params.basis = static_cast<uint8_t>(Basis);
  params.envelope = static_cast<uint8_t>(envelope);
  using Bound = typename PB::Stage::Warp<
      PB::Warp::VectorNoise<VectorNoiseMirrorProvider, Basis,
                            Envelope>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

template <math::NoiseBasis Basis>
inline void run_vector_noise_basis(In::ChainProgram &program,
                                   const In::FrameContext &ctx) {
  run_vector_noise_variant<Basis, PB::Warp::FlatEnvelope>(
      program, ctx, In::Op::WarpEnvelope::FLAT);
  run_vector_noise_variant<Basis, PB::Warp::ProjectionWeightEnvelope>(
      program, ctx, In::Op::WarpEnvelope::PROJECTION_WEIGHT);
  run_vector_noise_variant<Basis, PB::Warp::EdgeFadeEnvelope>(
      program, ctx, In::Op::WarpEnvelope::EDGE_FADE);
}

inline void test_shader_chain_parity_warp_vector_noise() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::VectorNoiseWarpParams>(
        program, In::Op::WarpVectorNoise::ID, 4, set);
    FastNoiseLite authored_noise;
    authored_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    authored_noise.SetSeed(static_cast<int32_t>(
        In::instance_hash("warp", In::Op::WarpVectorNoise::ID)));
    authored_noise.SetFrequency(1.0f);
    const auto &state = state_as<In::Op::NoisePhaseState>(program, 2);
    HS_EXPECT_EQ(state.noise.GetNoise(0.25f, -0.5f, 0.75f),
                 authored_noise.GetNoise(0.25f, -0.5f, 0.75f));
    const In::FrameContext ctx = shared_resources().context();
    run_vector_noise_basis<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_vector_noise_basis<math::NoiseBasis::FBM3>(program, ctx);
    run_vector_noise_basis<math::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}

/** Two instances of one noise operator own decorrelated fields, each seeded
    from its own instance identity. */
inline void test_shader_chain_noise_instances_decorrelate() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"warp-a", In::Op::WarpVectorNoise::ID},
      {"warp-b", In::Op::WarpVectorNoise::ID},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  const auto &first = state_as<In::Op::NoisePhaseState>(program, 2);
  const auto &second = state_as<In::Op::NoisePhaseState>(program, 3);
  bool differs = false;
  for (const math::Vector &view : sweep_views())
    if (first.noise.GetNoise(view.x, view.y, view.z) !=
        second.noise.GetNoise(view.x, view.y, view.z))
      differs = true;
  HS_EXPECT_TRUE(differs);
  // Deterministic: the seed is the instance's stable hash.
  FastNoiseLite authored;
  authored.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  authored.SetSeed(static_cast<int32_t>(
      In::instance_hash("warp-b", In::Op::WarpVectorNoise::ID)));
  authored.SetFrequency(1.0f);
  HS_EXPECT_EQ(second.noise.GetNoise(0.25f, -0.5f, 0.75f),
               authored.GetNoise(0.25f, -0.5f, 0.75f));
  program.clear();
}

template <typename Mode, uint8_t Harmonic>
inline void run_polar_variant(In::ChainProgram &program,
                              const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::PolarChartParams>(program, 2);
  params.mode =
      static_cast<uint8_t>(std::is_same_v<Mode, PB::Warp::LogarithmicPolar>
                               ? In::Op::PolarMode::LOGARITHMIC
                               : In::Op::PolarMode::LINEAR);
  params.harmonic = Harmonic - 1;
  using Bound = typename PB::Stage::Warp<
      PB::Warp::PolarChart<PolarChartMirrorProvider, Mode,
                           Harmonic>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

template <typename Mode>
inline void run_polar_harmonics(In::ChainProgram &program,
                                const In::FrameContext &ctx) {
  [&]<size_t... I>(std::index_sequence<I...>) {
    (run_polar_variant<Mode, static_cast<uint8_t>(I + 1)>(program, ctx), ...);
  }(std::make_index_sequence<PB::Warp::MAX_POLAR_HARMONIC>{});
}

inline void test_shader_chain_parity_warp_polar_chart() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::PolarChartParams>(
        program, In::Op::WarpPolarChart::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_polar_harmonics<PB::Warp::LinearPolar>(program, ctx);
    run_polar_harmonics<PB::Warp::LogarithmicPolar>(program, ctx);
    program.clear();
  }
}

template <math::NoiseBasis Basis, typename Integrator>
inline void run_curl_flow_variant(In::ChainProgram &program,
                                  const In::FrameContext &ctx,
                                  In::Op::CurlIntegrator integrator) {
  auto &params = param_as<In::Op::CurlFlowWarpParams>(program, 2);
  params.basis = static_cast<uint8_t>(Basis);
  params.integrator = static_cast<uint8_t>(integrator);
  using Bound = typename PB::Stage::Warp<
      PB::Warp::CurlFlow<CurlFlowMirrorProvider, Basis,
                         Integrator>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

template <math::NoiseBasis Basis>
inline void run_curl_flow_basis(In::ChainProgram &program,
                                const In::FrameContext &ctx) {
  run_curl_flow_variant<Basis, PB::Warp::Euler1>(
      program, ctx, In::Op::CurlIntegrator::EULER1);
  run_curl_flow_variant<Basis, PB::Warp::Midpoint2>(
      program, ctx, In::Op::CurlIntegrator::MIDPOINT2);
  run_curl_flow_variant<Basis, PB::Warp::Midpoint4>(
      program, ctx, In::Op::CurlIntegrator::MIDPOINT4);
}

inline void test_shader_chain_parity_warp_curl_flow() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::CurlFlowWarpParams>(
        program, In::Op::WarpCurlFlow::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_curl_flow_basis<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_curl_flow_basis<math::NoiseBasis::FBM3>(program, ctx);
    run_curl_flow_basis<math::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}
