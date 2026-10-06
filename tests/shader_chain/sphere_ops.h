/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- sphere-op mirrors ----------------------------------------------------

/** Mirror frame for the sphere endomorphism batch; populated from the erased
    program's own blocks so both sides read identical raw state. */
struct SphereOpMirrorFrame {
  const FastNoiseLite *noise = nullptr;
  In::Op::CurlDisplaceParams curl;
  In::Op::DirectDisplaceParams direct;
  PB::Surface::PeriodicRippleParams ripple;
  math::MobiusParams mobius;
  float phase = 0.0f;
};

struct SphereOpMirrorBinding {
  using FrameState = SphereOpMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct CurlMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static float scale(const FrameState &frame) { return frame.curl.scale; }
  static float strength(const FrameState &frame) { return frame.curl.strength; }
  static PB::Surface::PreparedLoop prepare(const FrameState &frame) {
    return PB::Surface::prepare(frame.phase);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct DirectMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static float scale(const FrameState &frame) { return frame.direct.scale; }
  static float strength(const FrameState &frame) {
    return frame.direct.strength;
  }
  static PB::Surface::PreparedDirect prepare(const FrameState &frame) {
    return PB::Surface::prepare_direct(frame.phase, frame.direct.direction);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct RippleMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const PB::Surface::PeriodicRippleParams &
  params(const FrameState &frame) {
    return frame.ripple;
  }
  static float phase(const FrameState &frame) {
    return PB::Surface::ripple_cycle(frame.phase, frame.ripple);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct MobiusMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const math::MobiusParams &params(const FrameState &frame) {
    return frame.mobius;
  }
};

/** Compiles a chain with one sphere endomorphism at entry 1 and applies the
    value set to every block. */
template <typename OpParams>
inline void arm_sphere_op_chain(In::ChainProgram &program, const char *op_id,
                                int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"op", op_id},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<OpParams>(program, 1), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 3), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 4), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

/** Erased-vs-bound parity of the sphere endomorphism at entry 1: identical
    seeds, bit-identical outputs across the sweep. */
template <typename BoundStage>
inline void expect_sphere_op_parity(In::ChainProgram &program,
                                    const In::FrameContext &ctx,
                                    const SphereOpMirrorFrame &mirror) {
  program.prepare(ctx);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  int view_index = 0;
  for (const math::Vector &view : sweep_views()) {
    HS_CONTEXT("view", view_index++);
    const PB::SphereSample seed{view, 0.25f};
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&seed, out, ctx, program.param_block(1),
                   program.prepared_block(1));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::SphereSample *>(out));
    const PB::SphereSample reference = BoundStage::run(seed, mirror, prepared);
    HS_EXPECT_TRUE(sphere_identical(erased, reference));
  }
}

inline SphereOpMirrorFrame sphere_op_mirror(In::ChainProgram &program) {
  SphereOpMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::DisplaceCurl::ID || id == In::Op::DisplaceDirect::ID) {
    const auto &state = state_as<In::Op::NoisePhaseState>(program, 1);
    mirror.noise = &state.noise;
    mirror.phase = state.phase;
  }
  if (id == In::Op::DisplaceRipple::ID) {
    const auto &state = state_as<In::Op::RipplePhaseState>(program, 1);
    mirror.ripple = param_as<PB::Surface::PeriodicRippleParams>(program, 1);
    mirror.phase = state.phase;
  }
  if (id == In::Op::DisplaceCurl::ID)
    mirror.curl = param_as<In::Op::CurlDisplaceParams>(program, 1);
  if (id == In::Op::DisplaceDirect::ID)
    mirror.direct = param_as<In::Op::DirectDisplaceParams>(program, 1);
  if (id == In::Op::LensMobius::ID) {
    const auto &params = param_as<In::Op::MobiusChainParams>(program, 1);
    mirror.mobius =
        math::MobiusParams{params.a_re, params.a_im, params.b_re, params.b_im,
                           params.c_re, params.c_im, params.d_re, params.d_im};
  }
  return mirror;
}

template <math::NoiseBasis Basis, typename Integrator>
inline void run_curl_displace_variant(In::ChainProgram &program,
                                      const In::FrameContext &ctx) {
  using Bound = typename PB::Stage::Displace<
      PB::Surface::CurlNoise<CurlMirrorProvider, Basis,
                             Integrator>>::template Bind<SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
}

template <math::NoiseBasis Basis>
inline void run_curl_displace_basis(In::ChainProgram &program,
                                    const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::CurlDisplaceParams>(program, 1);
  params.basis = static_cast<uint8_t>(Basis);
  params.integrator = static_cast<uint8_t>(PB::Surface::Integrator::EULER);
  run_curl_displace_variant<Basis, PB::Surface::Euler>(program, ctx);
  params.integrator = static_cast<uint8_t>(PB::Surface::Integrator::MIDPOINT);
  run_curl_displace_variant<Basis, PB::Surface::Midpoint>(program, ctx);
  params.integrator =
      static_cast<uint8_t>(PB::Surface::Integrator::MIDPOINT_2X);
  run_curl_displace_variant<Basis, PB::Surface::Midpoint2>(program, ctx);
}

inline void test_shader_chain_parity_displace_curl() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sphere_op_chain<In::Op::CurlDisplaceParams>(
        program, In::Op::DisplaceCurl::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_curl_displace_basis<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_curl_displace_basis<math::NoiseBasis::FBM3>(program, ctx);
    run_curl_displace_basis<math::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}

template <math::NoiseBasis Basis>
inline void run_direct_displace_variant(In::ChainProgram &program,
                                        const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::DirectDisplaceParams>(program, 1);
  params.basis = static_cast<uint8_t>(Basis);
  using Bound = typename PB::Stage::Displace<PB::Surface::DirectNoise<
      DirectMirrorProvider, Basis>>::template Bind<SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
}

inline void test_shader_chain_parity_displace_direct() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sphere_op_chain<In::Op::DirectDisplaceParams>(
        program, In::Op::DisplaceDirect::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_direct_displace_variant<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_direct_displace_variant<math::NoiseBasis::FBM3>(program, ctx);
    run_direct_displace_variant<math::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_displace_ripple() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sphere_op_chain<PB::Surface::PeriodicRippleParams>(
        program, In::Op::DisplaceRipple::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    using Bound = typename PB::Stage::Displace<PB::Surface::PeriodicRipple<
        RippleMirrorProvider>>::template Bind<SphereOpMirrorBinding>;
    expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
    program.clear();
  }
}

template <typename LensPolicy, typename OpParams>
inline void run_lens_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_sphere_op_chain<OpParams>(program, op_id, 3, set);
  if constexpr (std::is_same_v<OpParams, In::Op::MobiusChainParams>) {
    if (set != ValueSet::DEFAULTS)
      param_as<OpParams>(program, 1).d_re *= -1.0f;
  }
  const In::FrameContext ctx = shared_resources().context();
  using Bound = typename PB::Stage::Lens<LensPolicy>::template Bind<
      SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
  program.clear();
}

template <typename Policy>
inline void run_kaleidoscope_variant(In::ChainProgram &program,
                                     const In::FrameContext &ctx,
                                     In::Op::KaleidoscopeSymmetry symmetry) {
  param_as<In::Op::KaleidoscopeChainParams>(program, 1).symmetry =
      static_cast<uint8_t>(symmetry);
  using Bound =
      typename PB::Stage::Lens<Policy>::template Bind<SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
}

inline void test_shader_chain_parity_lens_ops() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    run_lens_parity<PB::Lens::Glitch, PB::Lens::NoLensParams>(
        In::Op::LensGlitch::ID, set);
    run_lens_parity<PB::Lens::Mobius<MobiusMirrorProvider>,
                    In::Op::MobiusChainParams>(In::Op::LensMobius::ID, set);
  }

  run_lens_parity<PB::Lens::Twist, In::Op::TwistChainParams>(
      In::Op::LensTwist::ID, ValueSet::DEFAULTS);
  for (const float rate : {-12.0f, 0.0f, 12.0f}) {
    const PB::SphereSample input{math::Vector(0.6f, 0.5f, 0.6244998f), 0.0f};
    In::Op::TwistChainParams params;
    params.twist_rate = rate;
    const auto result =
        In::Op::LensTwist::run(input, shared_resources().context(), params, {});
    const float angle = rate * input.dir.y;
    HS_EXPECT_NEAR(result.dir.x,
                   input.dir.x * cosf(angle) - input.dir.z * sinf(angle),
                   2e-3f);
    HS_EXPECT_EQ(result.dir.y, input.dir.y);
    HS_EXPECT_NEAR(result.dir.z,
                   input.dir.x * sinf(angle) + input.dir.z * cosf(angle),
                   2e-3f);
  }

  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_sphere_op_chain<In::Op::KaleidoscopeChainParams>(
      program, In::Op::LensKaleidoscope::ID, 3, ValueSet::DEFAULTS);
  const In::FrameContext ctx = shared_resources().context();
  using Symmetry = In::Op::KaleidoscopeSymmetry;
  run_kaleidoscope_variant<PB::Lens::Kaleidoscope>(program, ctx,
                                                   Symmetry::AZIMUTHAL);
  run_kaleidoscope_variant<PB::Lens::TetrahedralKaleidoscope>(
      program, ctx, Symmetry::TETRAHEDRAL);
  run_kaleidoscope_variant<PB::Lens::OctahedralKaleidoscope>(
      program, ctx, Symmetry::OCTAHEDRAL);
  run_kaleidoscope_variant<PB::Lens::DodecahedralKaleidoscope>(
      program, ctx, Symmetry::DODECAHEDRAL);
  run_kaleidoscope_variant<PB::Lens::TriangularPrismKaleidoscope>(
      program, ctx, Symmetry::TRIANGULAR_PRISM);
  run_kaleidoscope_variant<PB::Lens::SquarePrismKaleidoscope>(
      program, ctx, Symmetry::SQUARE_PRISM);
  run_kaleidoscope_variant<PB::Lens::PentagonalPrismKaleidoscope>(
      program, ctx, Symmetry::PENTAGONAL_PRISM);
  run_kaleidoscope_variant<PB::Lens::HexagonalPrismKaleidoscope>(
      program, ctx, Symmetry::HEXAGONAL_PRISM);
  run_kaleidoscope_variant<PB::Lens::OctagonalPrismKaleidoscope>(
      program, ctx, Symmetry::OCTAGONAL_PRISM);
  program.clear();
}
