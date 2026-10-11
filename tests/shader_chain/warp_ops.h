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
  arm_chain<In::Op::RotateChainParams, In::Op::ProjectChainParams, OpParams,
            In::Op::GridSampleParams, In::Op::GeneratedPaletteParams>(
      program, chain, frames, set);
}

/** Erased-vs-bound parity of the warp at entry 2. */
template <typename BoundStage>
inline void expect_warp_op_parity(In::ChainProgram &program,
                                  const In::FrameContext &ctx,
                                  const WarpMirrorFrame &mirror) {
  program.prepare(ctx);
  expect_op_parity<BoundStage>(
      program, ctx, 2, mirror,
      [&](const math::Vector &view) {
        return projected_input(program, ctx, view);
      },
      plane_identical);
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
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
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
  });
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
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    arm_warp_op_chain<In::Op::WaveShearWarpParams>(
        program, In::Op::WarpWaveShear::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_wave_shear_variant<PB::Warp::FlatEnvelope>(program, ctx,
                                                   In::Op::WarpEnvelope::FLAT);
    run_wave_shear_variant<PB::Warp::ProjectionWeightEnvelope>(
        program, ctx, In::Op::WarpEnvelope::PROJECTION_WEIGHT);
    run_wave_shear_variant<PB::Warp::EdgeFadeEnvelope>(
        program, ctx, In::Op::WarpEnvelope::EDGE_FADE);
  });
}

inline void test_shader_chain_parity_warp_vortex() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    arm_warp_op_chain<In::Op::VortexWarpParams>(program, In::Op::WarpVortex::ID,
                                                4, set);
    const In::FrameContext ctx = shared_resources().context();
    using Bound =
        typename PB::Stage::Warp<PB::Warp::Vortex<VortexMirrorProvider>>::
            template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
  });
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
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    arm_warp_op_chain<In::Op::VectorNoiseWarpParams>(
        program, In::Op::WarpVectorNoise::ID, 4, set);
    const FastNoiseLite authored_noise =
        make_oracle_noise("warp", In::Op::WarpVectorNoise::ID, 1.0f);
    const auto &state = state_as<In::Op::NoisePhaseState>(program, 2);
    HS_EXPECT_EQ(state.noise.GetNoise(0.25f, -0.5f, 0.75f),
                 authored_noise.GetNoise(0.25f, -0.5f, 0.75f));
    const In::FrameContext ctx = shared_resources().context();
    run_vector_noise_basis<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_vector_noise_basis<math::NoiseBasis::FBM3>(program, ctx);
    run_vector_noise_basis<math::NoiseBasis::RIDGED3>(program, ctx);
  });
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
  const FastNoiseLite authored =
      make_oracle_noise("warp-b", In::Op::WarpVectorNoise::ID, 1.0f);
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
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    arm_warp_op_chain<In::Op::PolarChartParams>(
        program, In::Op::WarpPolarChart::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_polar_harmonics<PB::Warp::LinearPolar>(program, ctx);
    run_polar_harmonics<PB::Warp::LogarithmicPolar>(program, ctx);
  });
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
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    arm_warp_op_chain<In::Op::CurlFlowWarpParams>(
        program, In::Op::WarpCurlFlow::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_curl_flow_basis<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_curl_flow_basis<math::NoiseBasis::FBM3>(program, ctx);
    run_curl_flow_basis<math::NoiseBasis::RIDGED3>(program, ctx);
  });
}

/** @brief Sweeps one warp's numeric parameter corners against its plane bound.
 * @tparam Model Plane-to-plane operator model. */
template <typename Model> inline void expect_warp_plane_bound() {
  using Params = typename Model::Params;
  typename Model::State state{};
  if constexpr (requires { Model::init(state, In::InstanceId{}); })
    Model::init(state,
                {"bound", Model::ID, In::instance_hash("bound", Model::ID)});
  const In::FrameContext ctx = shared_resources().context();
  const auto probe = [&](const Params &params) {
    for (float phase : {0.0f, 0.125f, 0.375f, 0.75f}) {
      state.phase = phase;
      const int rotations = std::is_same_v<Model, In::Op::WarpAffineV3> ? 8 : 1;
      for (int rotation = 0; rotation < rotations; ++rotation) {
        if constexpr (std::is_same_v<Model, In::Op::WarpAffineV3>)
          state.rotation = rotation * math::PI_F / 4.0f;
        const auto prepared = Model::prepare(ctx, params, state);
        for (float input_bound : {0.0f, 0.25f, 1.0f, 4.0f, 32.0f}) {
          const float bound = Model::plane_bound(params, input_bound);
          for (float radius_fraction : {0.0f, 0.5f, 1.0f}) {
            for (int direction = 0; direction < 32; ++direction) {
              const double angle = 2.0 * std::acos(-1.0) * direction / 32;
              PB::PlaneSample input{};
              input.coords = {input_bound * radius_fraction *
                                  static_cast<float>(std::cos(angle)),
                              input_bound * radius_fraction *
                                  static_cast<float>(std::sin(angle))};
              input.provenance.value_weight = 1.0f;
              input.provenance.fade_edge_distance = 2.0f;
              input.sphere = math::UP;
              const auto output = Model::run(input, ctx, params, prepared);
              HS_EXPECT_LE(std::hypot(output.coords.re, output.coords.im),
                           bound + 5e-5f * std::max(1.0f, bound));
            }
          }
        }
      }
    }
  };
  HS_CONTEXT(Model::ID);
  probe(Params{});
  if constexpr (std::is_same_v<Model, In::Op::WarpVortex>) {
    Params displaced{};
    displaced.center_x = 1.0f;
    displaced.radius = 8.0f;
    displaced.turns = 0.5f;
    probe(displaced);
  }
  constexpr size_t CORNERS = size_t{1} << Params::FIELDS.size();
  for (size_t corner = 0; corner < CORNERS; ++corner) {
    Params params{};
    for (size_t field = 0; field < Params::FIELDS.size(); ++field) {
      const auto &descriptor = Params::FIELDS[field];
      params.*descriptor.member =
          corner & (size_t{1} << field) ? descriptor.max : descriptor.min;
    }
    if constexpr (std::is_same_v<Model, In::Op::WarpVectorNoise>) {
      for (uint8_t basis = 0; basis < 3; ++basis)
        for (uint8_t envelope = 0; envelope < 3; ++envelope) {
          params.basis = basis;
          params.envelope = envelope;
          probe(params);
        }
    } else if constexpr (std::is_same_v<Model, In::Op::WarpPolarChart>) {
      for (uint8_t mode = 0; mode < 2; ++mode)
        for (uint8_t harmonic : {uint8_t{0}, uint8_t{15}}) {
          params.mode = mode;
          params.harmonic = harmonic;
          probe(params);
        }
    } else if constexpr (std::is_same_v<Model, In::Op::WarpCurlFlow>) {
      for (uint8_t basis = 0; basis < 3; ++basis)
        for (uint8_t integrator = 0; integrator < 3; ++integrator) {
          params.basis = basis;
          params.integrator = integrator;
          probe(params);
        }
    } else {
      probe(params);
    }
  }
}

/** @brief Warp bounds contain outputs at parameter extremes. */
inline void test_shader_chain_warp_plane_bound() {
  expect_warp_plane_bound<In::Op::WarpAffineV3>();
  expect_warp_plane_bound<In::Op::WarpWaveShear>();
  expect_warp_plane_bound<In::Op::WarpVortex>();
  expect_warp_plane_bound<In::Op::WarpVectorNoise>();
  expect_warp_plane_bound<In::Op::WarpMirrorTile>();
  expect_warp_plane_bound<In::Op::WarpPolarChart>();
  expect_warp_plane_bound<In::Op::WarpCurlFlow>();
}
