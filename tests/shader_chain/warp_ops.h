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
  PB::TranslationLattice lattice;
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
                             frame.lattice);
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
inline bool arm_warp_op_chain(In::ChainProgram &program, const char *op_id,
                              int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"warp", op_id},
      {"sample", "sample.grid.v4"},
      {"colorize", "colorize.generated-palette.v4"},
  };
  return arm_chain<In::Op::RotateChainParams, In::Op::ProjectChainParams,
                   OpParams, In::Op::GridSampleParams,
                   In::Op::GeneratedPaletteParams>(program, chain, frames, set);
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
  if (id == In::Op::WarpAffineV4::ID) {
    const auto &state = state_as<In::Op::AffineClockState>(program, 2);
    mirror.phase = state.phase;
    mirror.rotation = state.rotation;
    mirror.affine = param_as<In::Op::AffineWarpParams>(program, 2);
    mirror.lattice = In::Op::SampleGridV4::translation_lattice(
        param_as<In::Op::GridSampleParams>(program, 3));
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
    if (!arm_warp_op_chain<In::Op::AffineWarpParams>(
            program, In::Op::WarpAffineV4::ID, 4, set))
      return;
    using BoundAffine =
        typename PB::Stage::Warp<PB::Warp::AffineFrame<AffineMirrorProvider>>::
            template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<BoundAffine>(program, ctx, warp_mirror(program));
    program.clear();

    if (!arm_warp_op_chain<In::Op::MirrorWarpParams>(
            program, In::Op::WarpMirrorTile::ID, 4, set))
      return;
    using BoundMirror = typename PB::Stage::Warp<PB::Warp::MirrorTile<
        MirrorTileMirrorProvider>>::template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<BoundMirror>(program, ctx, warp_mirror(program));
  });
}

/** @brief A model's prepared block, passing @p lattice when it reads one.
 * @tparam Model Operator model. */
template <typename Model>
inline typename Model::Prepared prepare_with_lattice(
    const In::FrameContext &ctx, const typename Model::Params &params,
    const typename Model::State &state, const PB::TranslationLattice &lattice) {
  if constexpr (requires { Model::prepare(ctx, params, state, lattice); })
    return Model::prepare(ctx, params, state, lattice);
  else
    return Model::prepare(ctx, params, state);
}

/** @brief A model's plane bound, passing @p lattice when it reads one.
 * @tparam Model Operator model. */
template <typename Model>
inline float plane_bound_with_lattice(const typename Model::Params &params,
                                      float input_bound,
                                      const PB::TranslationLattice &lattice) {
  if constexpr (requires { Model::plane_bound(params, input_bound, lattice); })
    return Model::plane_bound(params, input_bound, lattice);
  else
    return Model::plane_bound(params, input_bound);
}

/** @brief Hexagonal lattice of the given side. */
inline PB::TranslationLattice hex_lattice(float side) {
  return {{side, 0.0f}, {0.5f * side, 0.5f * std::sqrt(3.0f) * side}};
}

/** @brief Winding resets over each tessellation kind preserve the sampled
    field and displacement hue. */
inline void test_shader_chain_affine_lattice_wrap() {
  using Model = In::Op::WarpAffineV4;
  using Sample = In::Op::SampleTessellation;
  const In::FrameContext ctx = shared_resources().context();
  for (const auto kind :
       {Pattern::TessellationKind::SQUARE, Pattern::TessellationKind::HEXAGONAL,
        Pattern::TessellationKind::TRIANGULAR}) {
    HS_CONTEXT("kind", static_cast<int>(kind));
    Sample::Params source;
    source.kind = static_cast<uint8_t>(kind);
    source.cell_scale = 4.0f;
    source.cell_shape = 1.0f;
    source.cell_softness = 0.2f;
    source.weight_mode = 0;
    source.coverage_mode = 0;
    const Sample::Prepared sample{};
    const PB::TranslationLattice downstream =
        Sample::translation_lattice(source);
    const bool square = kind == Pattern::TessellationKind::SQUARE;
    HS_EXPECT_NEAR(downstream.a.re, 0.25f, 0.0f);
    HS_EXPECT_NEAR(downstream.a.im, 0.0f, 0.0f);
    HS_EXPECT_NEAR(downstream.b.re, square ? 0.0f : 0.125f, 0.0f);
    HS_EXPECT_NEAR(downstream.b.im, square ? 0.25f : 0.2165063509f, 1e-7f);
    const float longest_reduced =
        square ? std::sqrt(2.0f) * 0.125f : 0.25f / std::sqrt(3.0f);
    PB::Color::ModulationState modulation;
    modulation.enabled = true;
    modulation.source = Pattern::SourceKind::DISPLACEMENT;
    modulation.mechanic.amount = 2.0f;
    for (const float direction : {-1.0f, 1.0f}) {
      Model::Params params;
      params.lattice_period = 3.0f;
      params.translation_x = direction * 1.7f;
      params.translation_y = direction * -2.2f;
      params.scale_x = 1.3f;
      params.scale_y = 0.7f;
      params.shear = 0.3f;
      const float wa = direction * 2.0f;
      const float wb = direction * -2.0f;
      const math::Complex shift{downstream.a.re * wa + downstream.b.re * wb,
                                downstream.a.im * wa + downstream.b.im * wb};
      const Model::State before{0.999999f, 0.73f};
      const Model::State after{0.0f, 0.73f};
      const auto from = Model::prepare(ctx, params, before, downstream);
      const auto to = Model::prepare(ctx, params, after, downstream);
      bool varied = false;
      float first_hue = 0.0f;
      for (int index = 0; index < 21; ++index) {
        PB::PlaneSample input{};
        input.coords = {0.071f * index - 0.61f, 0.039f * index - 0.17f};
        input.sphere = math::UP;
        input.path_length = 0.125f;
        const auto warped_from = Model::run(input, ctx, params, from);
        const auto warped_to = Model::run(input, ctx, params, to);
        HS_EXPECT_NEAR(warped_to.coords.re - warped_from.coords.re, shift.re,
                       2e-6f);
        HS_EXPECT_NEAR(warped_to.coords.im - warped_from.coords.im, shift.im,
                       2e-6f);
        PB::ColorSample color_from{};
        color_from.field = Sample::run(warped_from, ctx, source, sample);
        PB::ColorSample color_to{};
        color_to.field = Sample::run(warped_to, ctx, source, sample);
        HS_EXPECT_NEAR(color_from.field.value, color_to.field.value, 5e-5f);
        HS_EXPECT_NEAR(color_from.field.path_length, color_to.field.path_length,
                       3e-6f);
        HS_EXPECT_GE(color_to.field.path_length, input.path_length);
        HS_EXPECT_LE(color_to.field.path_length,
                     input.path_length + longest_reduced + 1e-6f);
        PB::Color::apply_modulation_in_place(color_from, modulation);
        PB::Color::apply_modulation_in_place(color_to, modulation);
        HS_EXPECT_NEAR(color_from.pending_hue, color_to.pending_hue, 6e-6f);
        if (index == 0)
          first_hue = color_to.pending_hue;
        else if (fabsf(color_to.pending_hue - first_hue) > 0.01f)
          varied = true;
      }
      HS_EXPECT_TRUE(varied);
    }
  }
}

/** @brief Square-lattice displacement remains continuous across its half-cell
    seam. */
inline void test_shader_chain_affine_lattice_distance() {
  using Model = In::Op::WarpAffineV4;
  const In::FrameContext ctx = shared_resources().context();
  const PB::TranslationLattice lattice = PB::square_lattice(0.4f);
  Model::Params params;
  params.translation_x = 1.0f;
  PB::PlaneSample input{};
  input.sphere = math::UP;
  input.path_length = 0.125f;
  for (float phase : {0.499999f, 0.5f, 0.500001f}) {
    const auto prepared = Model::prepare(ctx, params, {phase, 0.0f}, lattice);
    const auto output = Model::run(input, ctx, params, prepared);
    HS_EXPECT_NEAR(output.path_length, 0.325f, 1e-6f);
    HS_EXPECT_NEAR(output.coords.re, -0.4f * phase, 1e-7f);
  }
  params.translation_x = 1.51f;
  params.translation_y = -2.51f;
  const auto prepared = Model::prepare(ctx, params, {0.125f, 0.0f}, lattice);
  const auto output = Model::run(input, ctx, params, prepared);
  HS_EXPECT_NEAR(output.coords.re, -0.1f, 1e-7f);
  HS_EXPECT_NEAR(output.coords.im, 0.15f, 1e-7f);
  HS_EXPECT_NEAR(output.path_length, 0.125 + std::sqrt(0.1 * 0.1 + 0.15 * 0.15),
                 1e-7f);
  const auto end = Model::prepare(ctx, params, {0.999999f, 0.0f}, lattice);
  const auto end_output = Model::run(input, ctx, params, end);
  const float bound = Model::plane_bound(params, 0.0f, lattice);
  HS_EXPECT_NEAR(bound, 0.4 * std::sqrt(13.0), 2e-6f);
  HS_EXPECT_LE(std::hypot(end_output.coords.re, end_output.coords.im), bound);
  for (const auto &field : Model::Params::FIELDS)
    if (std::string_view(field.id) == "translation-x" ||
        std::string_view(field.id) == "translation-y")
      HS_EXPECT_EQ(field.curve, PB::FieldCurve::SNAP);
}

/** @brief Without a downstream lattice, translation winds the fallback period
    and the path length is the unwrapped displacement. */
inline void test_shader_chain_affine_fallback_period() {
  using Model = In::Op::WarpAffineV4;
  const In::FrameContext ctx = shared_resources().context();
  Model::Params params;
  params.lattice_period = 0.4f;
  params.repeat_cells = 2.0f;
  params.translation_x = 1.0f;
  params.translation_y = -1.0f;
  PB::PlaneSample input{};
  input.sphere = math::UP;
  const auto prepared = Model::prepare(ctx, params, {0.9f, 0.0f}, {});
  const auto output = Model::run(input, ctx, params, prepared);
  HS_EXPECT_NEAR(output.coords.re, -0.72f, 1e-6f);
  HS_EXPECT_NEAR(output.coords.im, 0.72f, 1e-6f);
  HS_EXPECT_NEAR(output.path_length, 0.72f * std::sqrt(2.0f), 1e-6f);
  HS_EXPECT_NEAR(Model::plane_bound(params, 0.0f, {}), 0.8f * std::sqrt(2.0f),
                 1e-6f);
  const auto lattice =
      Model::prepare(ctx, params, {0.9f, 0.0f}, PB::square_lattice(0.4f));
  HS_EXPECT_NEAR(Model::run(input, ctx, params, lattice).path_length,
                 0.08f * std::sqrt(2.0f), 1e-6f);
}

/** @brief A hexagonal lattice winds along its basis and measures path length
    to the closest repeat lattice vector. */
inline void test_shader_chain_affine_hex_lattice() {
  using Model = In::Op::WarpAffineV4;
  const In::FrameContext ctx = shared_resources().context();
  const PB::TranslationLattice lattice = hex_lattice(0.3f);
  Model::Params params;
  params.translation_x = 1.0f;
  params.translation_y = 2.0f;
  params.repeat_cells = 2.0f;
  PB::PlaneSample input{};
  input.sphere = math::UP;
  const auto turn = Model::prepare(ctx, params, {0.5f, 0.0f}, lattice);
  HS_EXPECT_NEAR(turn.affine.translation_x, 0.5f * 2.0f * (0.3f + 2.0f * 0.15f),
                 1e-6f);
  HS_EXPECT_NEAR(turn.affine.translation_y,
                 0.5f * 2.0f * 2.0f * 0.15f * std::sqrt(3.0f), 1e-6f);
  HS_EXPECT_EQ(turn.affine.repeat.reduction,
               PB::Warp::RepeatReduction::GENERAL);
  const PB::TranslationLattice repeat{{0.6f, 0.0f},
                                      {0.3f, 0.3f * std::sqrt(3.0f)}};
  for (const float rotation : {0.0f, 0.7f, 2.9f}) {
    for (const float phase : {0.0f, 0.13f, 0.5f, 0.77f, 0.999f}) {
      const auto prepared =
          Model::prepare(ctx, params, {phase, rotation}, lattice);
      for (int index = 0; index < 37; ++index) {
        input.coords = {0.173f * index - 3.1f, 0.211f * (index % 11) - 1.2f};
        const auto output = Model::run(input, ctx, params, prepared);
        const double dx = output.coords.re - input.coords.re;
        const double dy = output.coords.im - input.coords.im;
        double shortest = std::numeric_limits<double>::infinity();
        for (int i = -40; i <= 40; ++i)
          for (int j = -40; j <= 40; ++j)
            shortest = std::min(
                shortest, std::hypot(dx - i * repeat.a.re - j * repeat.b.re,
                                     dy - i * repeat.a.im - j * repeat.b.im));
        HS_EXPECT_NEAR(output.path_length, shortest, 2e-5);
      }
    }
  }
}

/** @brief Repeat reduction finds the closest vector for unreduced bases; the
    general reduction also holds for axis-aligned ones. */
inline void test_shader_chain_affine_repeat_reduction() {
  const std::array<std::pair<math::Complex, math::Complex>, 4> bases{{
      {{0.25f, 0.0f}, {0.0f, 0.25f}},
      {{0.3f, 0.0f}, {1.35f, 0.26f}},
      {{0.2f, 0.1f}, {-0.7f, 0.45f}},
      {{0.0f, 0.5f}, {0.3f, 0.0f}},
  }};
  for (const auto &[a, b] : bases) {
    const PB::Warp::PreparedRepeat repeat = PB::Warp::prepare_repeat(a, b);
    for (int index = 0; index < 61; ++index) {
      const math::Complex delta{0.37f * index - 9.0f,
                                0.29f * (index % 13) - 1.9f};
      const math::Complex reduced = PB::Warp::repeat_offset(delta, repeat);
      const math::Complex general =
          PB::Warp::repeat_offset<PB::Warp::RepeatReduction::GENERAL>(delta,
                                                                      repeat);
      double shortest = std::numeric_limits<double>::infinity();
      for (int i = -80; i <= 80; ++i)
        for (int j = -80; j <= 80; ++j)
          shortest = std::min(
              shortest, std::hypot(double{delta.re} - i * a.re - j * b.re,
                                   double{delta.im} - i * a.im - j * b.im));
      HS_EXPECT_NEAR(std::hypot(reduced.re, reduced.im), shortest, 2e-5);
      HS_EXPECT_NEAR(std::hypot(general.re, general.im), shortest, 2e-5);
    }
  }
  HS_EXPECT_EQ(PB::Warp::prepare_repeat({0.0f, 0.5f}, {0.3f, 0.0f}).reduction,
               PB::Warp::RepeatReduction::RECTANGULAR);
}

/** @brief A multi-cell repeat scales the winding reset and the hue modulus. */
inline void test_shader_chain_affine_repeat_cells() {
  using Model = In::Op::WarpAffineV4;
  const In::FrameContext ctx = shared_resources().context();
  In::Op::SampleTessellation::Params lattice;
  lattice.cell_scale = 4.0f;
  lattice.cell_shape = 1.0f;
  lattice.cell_softness = 0.2f;
  lattice.weight_mode = 0;
  lattice.coverage_mode = 0;
  const In::Op::SampleTessellation::Prepared sample{};
  const PB::TranslationLattice downstream =
      In::Op::SampleTessellation::translation_lattice(lattice);
  Model::Params params;
  params.repeat_cells = 2.6f;
  params.translation_x = 2.0f;
  params.translation_y = -1.0f;
  params.rotation_rate = 0.4f;
  const auto from = Model::prepare(ctx, params, {0.999999f, 0.6f}, downstream);
  const auto to = Model::prepare(ctx, params, {0.0f, 0.6f}, downstream);
  HS_EXPECT_NEAR(to.affine.repeat.a.re, 0.75f, 1e-7f);
  float longest = 0.0f;
  for (int index = 0; index < 41; ++index) {
    PB::PlaneSample input{};
    input.coords = {0.083f * index - 1.7f, 0.051f * index - 0.9f};
    input.sphere = math::UP;
    const auto warped_from = Model::run(input, ctx, params, from);
    const auto warped_to = Model::run(input, ctx, params, to);
    HS_EXPECT_NEAR(warped_to.coords.re - warped_from.coords.re, 1.5f, 2e-5f);
    HS_EXPECT_NEAR(warped_to.coords.im - warped_from.coords.im, -0.75f, 2e-5f);
    HS_EXPECT_NEAR(
        In::Op::SampleTessellation::run(warped_from, ctx, lattice, sample)
            .value,
        In::Op::SampleTessellation::run(warped_to, ctx, lattice, sample).value,
        2e-4f);
    HS_EXPECT_NEAR(warped_from.path_length, warped_to.path_length, 2e-5f);
    HS_EXPECT_LE(warped_to.path_length, std::sqrt(2.0f) * 0.375f + 1e-6f);
    longest = std::max(longest, warped_to.path_length);
  }
  HS_EXPECT_GT(longest, std::sqrt(2.0f) * 0.125f);
  for (const auto &field : Model::Params::FIELDS)
    if (std::string_view(field.id) == "repeat-cells")
      HS_EXPECT_EQ(field.curve, PB::FieldCurve::SNAP);
}

/** @brief Both travel directions retain a continuous affine column interval. */
inline void test_shader_chain_affine_interval() {
  using Model = In::Op::WarpAffineV4;
  const In::FrameContext ctx = shared_resources().context();
  const PB::TranslationLattice lattice = PB::square_lattice(0.25f);
  for (const float speed : {-0.002f, 0.002f}) {
    Model::Params params;
    params.speed = speed;
    params.rotation_rate = 0.7f;
    params.translation_x = 2.0f;
    params.translation_y = -3.0f;
    const Model::State before{speed > 0.0f ? 0.999f : 0.001f, 0.71f};
    auto after = before;
    Model::advance(after, params);
    HS_EXPECT_NEAR(after.phase, speed > 0.0f ? 0.001f : 0.999f, 1e-7f);
    const auto interval = In::Op::temporal_interval(before, after);
    PB::PlaneSample input{};
    input.coords = {0.13f, -0.21f};
    input.sphere = math::UP;
    float previous_length = 0.0f;
    for (int index = 0; index <= 8; ++index) {
      const float fraction = index / 8.0f;
      auto state = after;
      In::Op::temporal_sample(state, interval, params, fraction);
      const float expected_phase = before.phase + speed * fraction;
      HS_EXPECT_NEAR(state.phase, expected_phase, 2e-7f);
      HS_EXPECT_NEAR(state.rotation,
                     before.rotation + speed * params.rotation_rate * fraction,
                     1e-6f);
      const auto prepared = Model::prepare(ctx, params, state, lattice);
      const auto output = Model::run(input, ctx, params, prepared);
      if (index > 0)
        HS_EXPECT_LT(fabsf(output.path_length - previous_length), 0.001f);
      previous_length = output.path_length;
    }
  }
}

/** @brief The interpreter hands each affine entry the lattice of the entry
    immediately after it. */
inline void test_shader_chain_affine_reads_next_lattice() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest chain[] = {
      {"project", In::Op::ProjectStereographic::ID},
      {"far", In::Op::WarpAffineV4::ID},
      {"near", In::Op::WarpAffineV4::ID},
      {"sample", In::Op::SampleTessellation::ID},
      {"color", In::Op::ColorizeGeneratedPaletteV4::ID},
  };
  HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
  HS_EXPECT_TRUE(program.ops()[3].op->runtime.translation_lattice != nullptr);
  HS_EXPECT_TRUE(program.ops()[2].op->runtime.translation_lattice == nullptr);
  param_as<In::Op::SampleTessellation::Params>(program, 3).cell_scale = 4.0f;
  for (const size_t index : {size_t{1}, size_t{2}}) {
    auto &params = param_as<In::Op::AffineWarpParams>(program, index);
    params.lattice_period = 3.0f;
    params.translation_x = 1.0f;
    state_as<In::Op::AffineClockState>(program, index).phase = 0.5f;
  }
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  const auto translation = [&](size_t index) {
    return prepared_as<PB::Warp::PreparedAffineSlot>(program, index)
        .affine.translation_x;
  };
  HS_EXPECT_NEAR(translation(1), 1.5f, 1e-7f);
  HS_EXPECT_NEAR(translation(2), 0.125f, 1e-7f);
  const int overflow = In::plane_bound_overflow(
      program.ops(), [&](size_t index) { return program.param_block(index); });
  HS_EXPECT_EQ(overflow, -1);
  HS_EXPECT_NEAR(
      program.ops()[2].op->runtime.plane_bound(
          program.param_block(2), 0.0f,
          In::ChainProgram::downstream_lattice(
              program.ops(), 2,
              [&](size_t index) { return program.param_block(index); })),
      0.25f, 1e-7f);
}

/** @brief A column prepares each affine entry from the next entry's
    column-sampled params. */
inline void test_shader_chain_affine_column_lattice() {
  auto blocks = std::make_unique<TemporalBlocks>();
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  program.bind_temporal_storage(blocks->history, blocks->column);
  const In::ChainEntryRequest chain[] = {
      {"project", In::Op::ProjectStereographic::ID},
      {"warp", In::Op::WarpAffineV4::ID},
      {"sample", In::Op::SampleTessellation::ID},
      {"color", In::Op::ColorizeGeneratedPaletteV4::ID},
  };
  HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
  auto &params = param_as<In::Op::AffineWarpParams>(program, 1);
  params.translation_x = 1.0f;
  state_as<In::Op::AffineClockState>(program, 1).phase = 0.5f;
  auto &source = param_as<In::Op::SampleTessellation::Params>(program, 2);
  source.cell_scale = 2.0f;
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  program.capture_history();
  source.cell_scale = 8.0f;
  program.prepare(ctx);
  HS_EXPECT_TRUE(program.history_continuous());
  program.prepare_column(ctx, ctx, 0.5f);
  const auto &prepared = prepared_as<PB::Warp::PreparedAffineSlot>(program, 1);
  HS_EXPECT_NEAR(prepared.affine.translation_x, 0.125f, 1e-6f);
  program.prepare(ctx);
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
    if (!arm_warp_op_chain<In::Op::WaveShearWarpParams>(
            program, In::Op::WarpWaveShear::ID, 4, set))
      return;
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
    if (!arm_warp_op_chain<In::Op::VortexWarpParams>(
            program, In::Op::WarpVortex::ID, 4, set))
      return;
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
    if (!arm_warp_op_chain<In::Op::VectorNoiseWarpParams>(
            program, In::Op::WarpVectorNoise::ID, 4, set))
      return;
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
      {"sample", "sample.grid.v4"},
      {"colorize", "colorize.generated-palette.v4"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  if (refusal.code != In::ChainStatus::OK)
    return;
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
    if (!arm_warp_op_chain<In::Op::PolarChartParams>(
            program, In::Op::WarpPolarChart::ID, 4, set))
      return;
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
    if (!arm_warp_op_chain<In::Op::CurlFlowWarpParams>(
            program, In::Op::WarpCurlFlow::ID, 4, set))
      return;
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
  PB::TranslationLattice lattice{};
  const auto probe = [&](const Params &params) {
    for (float phase : {0.0f, 0.125f, 0.375f, 0.75f}) {
      state.phase = phase;
      const int rotations = std::is_same_v<Model, In::Op::WarpAffineV4> ? 8 : 1;
      for (int rotation = 0; rotation < rotations; ++rotation) {
        if constexpr (std::is_same_v<Model, In::Op::WarpAffineV4>)
          state.rotation = rotation * math::PI_F / 4.0f;
        const auto prepared =
            prepare_with_lattice<Model>(ctx, params, state, lattice);
        for (float input_bound : {0.0f, 0.25f, 1.0f, 4.0f, 32.0f}) {
          const float bound =
              plane_bound_with_lattice<Model>(params, input_bound, lattice);
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
    } else if constexpr (std::is_same_v<Model, In::Op::WarpAffineV4>) {
      for (const PB::TranslationLattice candidate :
           {PB::TranslationLattice{}, hex_lattice(0.3f)}) {
        lattice = candidate;
        probe(params);
      }
      lattice = {};
    } else {
      probe(params);
    }
  }
}

/** @brief Warp bounds contain outputs at parameter extremes. */
inline void test_shader_chain_warp_plane_bound() {
  expect_warp_plane_bound<In::Op::WarpAffineV4>();
  expect_warp_plane_bound<In::Op::WarpWaveShear>();
  expect_warp_plane_bound<In::Op::WarpVortex>();
  expect_warp_plane_bound<In::Op::WarpVectorNoise>();
  expect_warp_plane_bound<In::Op::WarpMirrorTile>();
  expect_warp_plane_bound<In::Op::WarpPolarChart>();
  expect_warp_plane_bound<In::Op::WarpCurlFlow>();
}
