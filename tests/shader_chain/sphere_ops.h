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
  PB::SphereAnimation::PreparedRipples ripples;
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
  static const PB::SphereAnimation::PreparedRipples &
  prepare(const FrameState &frame) {
    return frame.ripples;
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
inline bool arm_sphere_op_chain(In::ChainProgram &program, const char *op_id,
                                int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"op", op_id},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"},
  };
  return arm_chain<In::Op::RotateChainParams, OpParams,
                   In::Op::ProjectChainParams, In::Op::GridSampleParams,
                   In::Op::GeneratedPaletteParams>(program, chain, frames, set);
}

/** Erased-vs-bound parity of the sphere endomorphism at entry 1: identical
    seeds, bit-identical outputs across the sweep. */
template <typename BoundStage>
inline void expect_sphere_op_parity(In::ChainProgram &program,
                                    const In::FrameContext &ctx,
                                    const SphereOpMirrorFrame &mirror) {
  program.prepare(ctx);
  expect_op_parity<BoundStage>(
      program, ctx, 1, mirror,
      [](const math::Vector &view) { return PB::SphereSample{view, 0.25f}; },
      sphere_identical);
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
    mirror.ripples = In::Op::DisplaceRipple::prepare(
        shared_resources().context(),
        param_as<In::Op::RippleChainParams>(program, 1), state);
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

template <math::NoiseBasis Basis, typename Integrator, math::TangentLimit Limit>
inline void run_curl_displace_variant(In::ChainProgram &program,
                                      const In::FrameContext &ctx) {
  using Bound = typename PB::Stage::Displace<
      PB::Surface::CurlNoise<CurlMirrorProvider, Basis, Integrator,
                             Limit>>::template Bind<SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
}

template <math::NoiseBasis Basis, math::TangentLimit Limit>
inline void run_curl_displace_integrators(In::ChainProgram &program,
                                          const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::CurlDisplaceParams>(program, 1);
  params.basis = static_cast<uint8_t>(Basis);
  params.limit = static_cast<uint8_t>(Limit);
  params.integrator = static_cast<uint8_t>(PB::Surface::Integrator::EULER);
  run_curl_displace_variant<Basis, PB::Surface::Euler, Limit>(program, ctx);
  params.integrator = static_cast<uint8_t>(PB::Surface::Integrator::MIDPOINT);
  run_curl_displace_variant<Basis, PB::Surface::Midpoint, Limit>(program, ctx);
  params.integrator =
      static_cast<uint8_t>(PB::Surface::Integrator::MIDPOINT_2X);
  run_curl_displace_variant<Basis, PB::Surface::Midpoint2, Limit>(program, ctx);
}

template <math::NoiseBasis Basis>
inline void run_curl_displace_basis(In::ChainProgram &program,
                                    const In::FrameContext &ctx) {
  run_curl_displace_integrators<Basis, math::TangentLimit::SMOOTH>(program,
                                                                   ctx);
  run_curl_displace_integrators<Basis, math::TangentLimit::CLAMP>(program, ctx);
}

inline void test_shader_chain_parity_displace_curl() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sphere_op_chain<In::Op::CurlDisplaceParams>(
            program, In::Op::DisplaceCurl::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_curl_displace_basis<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_curl_displace_basis<math::NoiseBasis::FBM3>(program, ctx);
    run_curl_displace_basis<math::NoiseBasis::RIDGED3>(program, ctx);
  });
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
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sphere_op_chain<In::Op::DirectDisplaceParams>(
            program, In::Op::DisplaceDirect::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_direct_displace_variant<math::NoiseBasis::SIMPLEX>(program, ctx);
    run_direct_displace_variant<math::NoiseBasis::FBM3>(program, ctx);
    run_direct_displace_variant<math::NoiseBasis::RIDGED3>(program, ctx);
  });
}

inline void test_shader_chain_parity_displace_ripple() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_sphere_op_chain<In::Op::RippleChainParams>(
            program, In::Op::DisplaceRipple::ID, 50, set))
      return;
    param_as<In::Op::RippleChainParams>(program, 1).random_center = 0;
    const In::FrameContext ctx = shared_resources().context();
    using Bound = typename PB::Stage::Displace<PB::Surface::BurstRipple<
        RippleMirrorProvider>>::template Bind<SphereOpMirrorBinding>;
    expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
  });
}

/** @brief Ripple bursts retain their stagger, cadence and shared centers. */
inline void test_shader_chain_ripple_bursts() {
  using Ripple = In::Op::DisplaceRipple;
  Ripple::Params params;
  params.burst_rate = 1.0f / 64.0f;
  Ripple::State state;
  Ripple::init(state,
               {"ripple", Ripple::ID, In::instance_hash("ripple", Ripple::ID)});
  const auto ctx = shared_resources().context();
  for (int frame = 0; frame < 17; ++frame)
    Ripple::advance(state, params);
  auto prepared = Ripple::prepare(ctx, params, state);
  HS_EXPECT_EQ(prepared.count, 2);
  HS_EXPECT_NEAR(prepared.waves[0].ripple.phase, 17.0f / 80.0f * math::PI_F,
                 1e-6f);
  HS_EXPECT_NEAR(prepared.waves[1].ripple.phase, 1.0f / 80.0f * math::PI_F,
                 1e-6f);
  const auto center = prepared.waves[0].ripple.center;
  HS_EXPECT_EQ(prepared.waves[1].ripple.center.x, center.x);
  HS_EXPECT_EQ(prepared.waves[1].ripple.center.y, center.y);
  HS_EXPECT_EQ(prepared.waves[1].ripple.center.z, center.z);
  const auto snapshot = In::RuntimeStateCodec<Ripple::State>::capture(state);
  Ripple::State restored;
  HS_EXPECT_TRUE(
      In::RuntimeStateCodec<Ripple::State>::restore(restored, snapshot));
  for (int frame = 17; frame < 63; ++frame) {
    Ripple::advance(state, params);
    Ripple::advance(restored, params);
  }
  HS_EXPECT_EQ(state.burst_count, 1U);
  Ripple::advance(state, params);
  Ripple::advance(restored, params);
  HS_EXPECT_EQ(state.burst_count, 2U);
  HS_EXPECT_TRUE(
      snapshot_words(In::RuntimeStateCodec<Ripple::State>::capture(state)) ==
      snapshot_words(In::RuntimeStateCodec<Ripple::State>::capture(restored)));
  prepared = Ripple::prepare(ctx, params, state);
  HS_EXPECT_EQ(prepared.count, 5);
  const auto next_center = prepared.waves[4].ripple.center;
  HS_EXPECT_TRUE((next_center - center).magnitude() > 0.01f);
  params.random_center = 0;
  params.center_azimuth = 0.7f;
  params.center_polar = 1.2f;
  prepared = Ripple::prepare(ctx, params, state);
  const auto manual = math::Vector::from_spherical(0.7f, 1.2f);
  for (uint8_t i = 0; i < prepared.count; ++i) {
    HS_EXPECT_EQ(prepared.waves[i].ripple.center.x, manual.x);
    HS_EXPECT_EQ(prepared.waves[i].ripple.center.y, manual.y);
    HS_EXPECT_EQ(prepared.waves[i].ripple.center.z, manual.z);
  }
  params.burst_rate = 0;
  for (int frame = 0; frame < 200; ++frame)
    Ripple::advance(state, params);
  HS_EXPECT_EQ(state.burst_count, 2U);
  HS_EXPECT_EQ(Ripple::prepare(ctx, params, state).count, 0);
  Ripple::State single;
  params.burst_size = 0;
  for (int frame = 0; frame < 50; ++frame)
    Ripple::advance(single, params);
  HS_EXPECT_EQ(Ripple::prepare(ctx, params, single).count, 1);
}

template <typename LensPolicy, typename OpParams>
inline void run_lens_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  if (!arm_sphere_op_chain<OpParams>(program, op_id, 3, set))
    return;
  if constexpr (std::is_same_v<OpParams, In::Op::MobiusChainParams>) {
    if (set != ValueSet::DEFAULTS)
      param_as<OpParams>(program, 1).d_re *= -1.0f;
  }
  if constexpr (std::is_same_v<OpParams, In::Op::TwistChainParams>)
    const_cast<In::Op::TwistOrbitState &>(
        state_as<In::Op::TwistOrbitState>(program, 1)) = {};
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
  run_lens_parity<PB::Lens::Glitch, In::Op::GlitchChainParams>(
      In::Op::LensGlitch::ID, ValueSet::DEFAULTS);
  for_each_value_set([&](ValueSet set) {
    run_lens_parity<PB::Lens::Mobius<MobiusMirrorProvider>,
                    In::Op::MobiusChainParams>(In::Op::LensMobius::ID, set);
  });

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
  if (!arm_sphere_op_chain<In::Op::KaleidoscopeChainParams>(
          program, In::Op::LensKaleidoscope::ID, 3, ValueSet::DEFAULTS))
    return;
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

/** @brief Lens controls alter direction while preserving path length. */
inline void test_shader_chain_lens_controls() {
  const auto *twist = require_operator(In::Op::LensTwist::ID);
  const auto *glitch = require_operator(In::Op::LensGlitch::ID);
  if (!twist || !glitch)
    return;
  const auto ctx = shared_resources().context();
  const PB::SphereSample input{math::Vector(0.6f, 0.5f, sqrtf(0.39f)), 0.25f};
  for (const float phi : {0.4f, math::PI_F * 0.5f, 2.2f}) {
    for (const float lambda : {0.0f, 0.7f, math::PI_F * 0.5f}) {
      In::Op::TwistChainParams params;
      In::Op::LensTwist::State state;
      state.phi = phi;
      state.lambda = lambda;
      const auto prepared = In::Op::LensTwist::prepare(ctx, params, state);
      const auto result = In::Op::LensTwist::run(input, ctx, params, prepared);
      const double ax = sin(phi) * cos(lambda);
      const double ay = cos(phi);
      const double az = sin(phi) * sin(lambda);
      const double h = input.dir.x * ax + input.dir.y * ay + input.dir.z * az;
      const double c = cos(params.twist_rate * h);
      const double s = sin(params.twist_rate * h);
      HS_EXPECT_NEAR(result.dir.x,
                     input.dir.x * c +
                         (input.dir.y * az - input.dir.z * ay) * s +
                         ax * h * (1 - c),
                     2e-3);
      HS_EXPECT_NEAR(result.dir.y,
                     input.dir.y * c +
                         (input.dir.z * ax - input.dir.x * az) * s +
                         ay * h * (1 - c),
                     2e-3);
      HS_EXPECT_NEAR(result.dir.z,
                     input.dir.z * c +
                         (input.dir.x * ay - input.dir.y * ax) * s +
                         az * h * (1 - c),
                     2e-3);
      HS_EXPECT_EQ(result.path_length, input.path_length);
    }
  }
  for (const float polar : {-12.0f, 0.0f, 1.0f, 2.0f, 12.0f}) {
    for (const float azimuth : {-12.0f, 0.0f, 1.0f, 3.0f, 12.0f}) {
      In::Op::GlitchChainParams params;
      *static_cast<float *>(glitch->runtime.param_address(&params, 0)) = polar;
      *static_cast<float *>(glitch->runtime.param_address(&params, 1)) =
          azimuth;
      const auto result = In::Op::LensGlitch::run(input, ctx, params, {});
      const double p = acos(input.dir.y) * polar;
      const double a = atan2(input.dir.z, input.dir.x) * azimuth;
      HS_EXPECT_NEAR(result.dir.x, sin(p) * cos(a), 3e-6);
      HS_EXPECT_NEAR(result.dir.y, cos(p), 3e-6);
      HS_EXPECT_NEAR(result.dir.z, sin(p) * sin(a), 3e-6);
      HS_EXPECT_EQ(result.path_length, input.path_length);
    }
  }
}

/** @brief Twist origin angles circle independently and interpolate across wraps. */
inline void test_shader_chain_twist_animation() {
  using Twist = In::Op::LensTwist;
  const auto ctx = shared_resources().context();
  const auto expect_axis = [&](const Twist::State &state, double phi,
                               double lambda) {
    const auto axis = Twist::prepare(ctx, {}, state).axis;
    HS_EXPECT_NEAR(axis.x, sin(phi) * cos(lambda), 1e-4);
    HS_EXPECT_NEAR(axis.y, cos(phi), 1e-4);
    HS_EXPECT_NEAR(axis.z, sin(phi) * sin(lambda), 1e-4);
  };
  for (const float speed : {0.0f, 0.01f, 0.04f}) {
    for (const float drift : {0.0f, 0.013f}) {
      Twist::Params params;
      params.speed = speed;
      params.drift = drift;
      Twist::State state{0.8f, 0.3f};
      for (int frame = 1; frame <= 700; ++frame) {
        Twist::advance(state, params);
        expect_axis(state, 0.8f + double(frame) * speed,
                    0.3f + double(frame) * drift);
      }
    }
  }
  Twist::Params params;
  params.speed = 0.04f;
  params.drift = 0.013f;
  const Twist::State start{math::TWO_PI_F - 0.02f, math::TWO_PI_F - 0.006f};
  auto end = start;
  Twist::advance(end, params);
  HS_EXPECT_TRUE(end.phi < 0.03f);
  HS_EXPECT_TRUE(end.lambda < 0.01f);
  const auto interval = In::Op::temporal_interval(start, end);
  for (const float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
    auto column = end;
    In::Op::temporal_sample(column, interval, params, fraction);
    expect_axis(column, double(start.phi) + params.speed * fraction,
                double(start.lambda) + params.drift * fraction);
  }
  const auto snapshot = In::RuntimeStateCodec<Twist::State>::capture(end);
  Twist::State restored;
  HS_EXPECT_TRUE(
      In::RuntimeStateCodec<Twist::State>::restore(restored, snapshot));
  Twist::advance(end, params);
  Twist::advance(restored, params);
  HS_EXPECT_EQ(restored.phi, end.phi);
  HS_EXPECT_EQ(restored.lambda, end.lambda);
  auto invalid = std::get<In::SourceClockSnapshot>(snapshot);
  invalid.secondary = math::TWO_PI_F;
  HS_EXPECT_FALSE(
      In::RuntimeStateCodec<Twist::State>::restore(restored, invalid));
  HS_EXPECT_EQ(restored.phi, end.phi);
  HS_EXPECT_EQ(restored.lambda, end.lambda);
}
