/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file spherical_sources.h
 * @brief Spherical source carriers, composed policies and instance clocks. */

/** @brief Composed frame for one spherical source resource.
 * @tparam Family Source parameter family. */
template <typename Family> struct SphericalSourceFrame {
  using Resource =
      PB::ParameterResource<"source", Family, PB::ResourceKind::SOURCE>;
  using Params = PB::ComposedDetail::ParameterSet<
      PB::ComposedDetail::ResourceList<Resource>>;
  using Frame = PB::FrameIntervalState<Params>;
  using Binding = PB::Binding<Frame>;
  using Policy = typename PB::SourcePolicyFor<Family, Binding>::Type;
  using Stage = PB::Stage::SampleSphere<Policy>;
  static_assert(
      std::is_same_v<typename PB::ComposedDetail::PolicyResources<Policy>::Type,
                     PB::ComposedDetail::ResourceList<Resource>>);
};

/** @brief Spherical noise adapters preserve the canonical source and carrier. */
inline void test_shader_chain_spherical_noise_policy() {
  using Model = In::Op::SampleSphericalNoise;
  using Fixture = SphericalSourceFrame<Model::Params>;
  using Binding = Fixture::Binding;
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"source">();
  params.noise_scale = 3.5f;
  params.noise_contrast = 1.7f;
  params.noise_time_rate = -1.0f / 128.0f;
  Model::State state;
  Model::init(state,
              {"source", Model::ID, In::instance_hash("source", Model::ID)});
  state.phase = 0.37f;
  auto &clock = frame.resources.template get<"source">();
  clock.noise = &state.noise;
  clock.noise_time = state.phase;
  const auto ORACLE_NOISE = make_oracle_noise("source", Model::ID, 1.0f);
  for (uint8_t basis = 0; basis < 3; ++basis) {
    params.basis = basis;
    const auto PREPARED = Model::prepare({}, params, state);
    const auto TYPED = Fixture::Stage::template prepare<Binding>(frame);
    const auto ORACLE = Pattern::SphereNoise::prepare(
        {3.5f, -1.0f / 128.0f, static_cast<math::NoiseBasis>(basis), 1.7f,
         Pattern::NoiseVariant::OCTAVE_CONTOUR},
        {0.37f}, ORACLE_NOISE);
    for (const auto &view : sweep_views()) {
      const float UNIT_RAW = Pattern::SphereNoise::sample(
          Pattern::SphereDirection::from(view), ORACLE);
      for (float scale : {1.0f, std::sqrt(1.003f)}) {
        const PB::SphereSample INPUT{scale * view, 0.375f};
        const auto ACTUAL = Model::run(INPUT, {}, params, PREPARED);
        const auto COMPOSED =
            Fixture::Stage::template run<Binding>(INPUT, frame, TYPED);
        const float RAW = Pattern::SphereNoise::sample(
            Pattern::SphereDirection::from(INPUT.dir), ORACLE);
        HS_EXPECT_NEAR(ACTUAL.value, 0.5f + 0.5f * RAW, 1.0e-7f);
        HS_EXPECT_NEAR(ACTUAL.value, 0.5f + 0.5f * UNIT_RAW, 2.0e-5f);
        HS_EXPECT_TRUE(field_identical(ACTUAL, COMPOSED));
        HS_EXPECT_EQ(ACTUAL.coverage, 1.0f);
        HS_EXPECT_EQ(ACTUAL.path_length, 0.375f);
        HS_EXPECT_EQ(ACTUAL.sphere.x, INPUT.dir.x);
        HS_EXPECT_EQ(ACTUAL.sphere.y, INPUT.dir.y);
        HS_EXPECT_EQ(ACTUAL.sphere.z, INPUT.dir.z);
      }
    }
  }
}

/** @brief Spherical cell selectors agree between typed and runtime crossings. */
inline void test_shader_chain_spherical_cells_policy() {
  using Model = In::Op::SampleSphericalCells;
  using Fixture = SphericalSourceFrame<Model::Params>;
  using Binding = Fixture::Binding;
  typename Fixture::Frame frame{};
  auto &params = frame.params.template get<"source">();
  params.count = 17.0f;
  params.seed = 391.0f;
  for (uint8_t layout = 0;
       layout <= static_cast<uint8_t>(Pattern::SphereCellLayout::SCATTERED);
       ++layout) {
    params.layout = layout;
    for (uint8_t output = 0;
         output <= static_cast<uint8_t>(Pattern::SphereCellOutput::CELL_VALUE);
         ++output) {
      params.output = output;
      const auto PREPARED = Model::prepare({}, params, {});
      const auto TYPED = Fixture::Stage::template prepare<Binding>(frame);
      for (const auto &view : sweep_views()) {
        const float UNIT_RAW = Pattern::SphereCells::sample(
            Pattern::SphereDirection::from(view), PREPARED);
        for (float scale : {1.0f, std::sqrt(1.003f)}) {
          const PB::SphereSample INPUT{scale * view, 1.25f};
          const auto ACTUAL = Model::run(INPUT, {}, params, PREPARED);
          const auto COMPOSED =
              Fixture::Stage::template run<Binding>(INPUT, frame, TYPED);
          HS_EXPECT_TRUE(field_identical(ACTUAL, COMPOSED));
          const float RAW = Pattern::SphereCells::sample(
              Pattern::SphereDirection::from(INPUT.dir), PREPARED);
          HS_EXPECT_NEAR(ACTUAL.value, 0.5f + 0.5f * RAW, 1.0e-7f);
          if (output !=
              static_cast<uint8_t>(Pattern::SphereCellOutput::CELL_VALUE))
            HS_EXPECT_NEAR(ACTUAL.value, 0.5f + 0.5f * UNIT_RAW, 3.0e-6f);
          HS_EXPECT_TRUE(ACTUAL.value >= 0.0f && ACTUAL.value <= 1.0f);
          HS_EXPECT_EQ(ACTUAL.coverage, 1.0f);
          HS_EXPECT_EQ(ACTUAL.path_length, 1.25f);
          HS_EXPECT_EQ(ACTUAL.sphere.x, INPUT.dir.x);
          HS_EXPECT_EQ(ACTUAL.sphere.y, INPUT.dir.y);
          HS_EXPECT_EQ(ACTUAL.sphere.z, INPUT.dir.z);
        }
      }
    }
  }
  params.layout = static_cast<uint8_t>(Pattern::SphereCellLayout::OCTAHEDRON);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::CENTERS);
  const PB::SphereSample INPUT{math::X_AXIS, 2.0f};
  const auto PREPARED = Model::prepare({}, params, {});
  HS_EXPECT_NEAR(Model::run(INPUT, {}, params, PREPARED).value, 1.0f, 1.0e-6f);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::BORDERS);
  params.line_thickness = 0.1f;
  params.line_softness = 0.2f;
  const float ANGLE =
      std::numbers::pi_v<float> / 4.0f -
      0.2f / std::sqrt(6.0f / (4.0f * std::numbers::pi_v<float>));
  const math::Vector BORDER{std::cos(ANGLE), std::sin(ANGLE), 0.0f};
  const auto BORDER_PREPARED = Model::prepare({}, params, {});
  const auto BORDER_TYPED = Fixture::Stage::template prepare<Binding>(frame);
  for (float scale : {1.0f, std::sqrt(1.003f)}) {
    const PB::SphereSample DRIFTED{scale * BORDER, 2.0f};
    const auto ACTUAL = Model::run(DRIFTED, {}, params, BORDER_PREPARED);
    const auto COMPOSED =
        Fixture::Stage::template run<Binding>(DRIFTED, frame, BORDER_TYPED);
    HS_EXPECT_NEAR(ACTUAL.value, 0.5f, 2.0e-6f);
    HS_EXPECT_NEAR(COMPOSED.value, 0.5f, 2.0e-6f);
    HS_EXPECT_EQ(ACTUAL.sphere.x, DRIFTED.dir.x);
    HS_EXPECT_EQ(COMPOSED.sphere.x, DRIFTED.dir.x);
    HS_EXPECT_EQ(ACTUAL.path_length, 2.0f);
    HS_EXPECT_EQ(COMPOSED.path_length, 2.0f);
  }
}

/** @brief A spherical source reaches cutout and palette without a projection. */
inline void test_shader_chain_spherical_source_pipeline() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest CHAIN[] = {
      {"source", In::Op::SampleSphericalCells::ID},
      {"cutout", In::Op::CoverageValueCutout::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  const auto RESULT = program.compile(CHAIN);
  HS_EXPECT_EQ(RESULT.code, In::ChainStatus::OK);
  if (RESULT.code != In::ChainStatus::OK)
    return;
  auto &params = param_as<In::Op::SampleSphericalCells::Params>(program, 0);
  params.layout = static_cast<uint8_t>(Pattern::SphereCellLayout::OCTAHEDRON);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::CENTERS);
  auto &cutout = param_as<In::Op::ValueCutoutChainParams>(program, 1);
  cutout.cutout_threshold = 0.75f;
  cutout.cutout_softness = 0.0f;
  const auto CTX = shared_resources().context();
  program.prepare(CTX);
  const PB::SphereSample CENTER{math::X_AXIS, 0.25f};
  const PB::SphereSample EDGE{math::Vector(1, 1, 0).normalized(), 0.25f};
  const auto CENTER_FIELD = run_entry<PB::FieldSample>(program, CTX, 0, CENTER);
  const auto EDGE_FIELD = run_entry<PB::FieldSample>(program, CTX, 0, EDGE);
  HS_EXPECT_NEAR(CENTER_FIELD.value, 1.0f, 1.0e-6f);
  HS_EXPECT_TRUE(EDGE_FIELD.value < 0.75f);
  const auto VISIBLE =
      run_entry<PB::FieldSample>(program, CTX, 1, CENTER_FIELD);
  const auto HIDDEN = run_entry<PB::FieldSample>(program, CTX, 1, EDGE_FIELD);
  HS_EXPECT_EQ(VISIBLE.coverage, 1.0f);
  HS_EXPECT_EQ(HIDDEN.coverage, 0.0f);
  HS_EXPECT_EQ(HIDDEN.path_length, 0.25f);
  HS_EXPECT_EQ(HIDDEN.sphere.x, EDGE.dir.x);
  HS_EXPECT_EQ(HIDDEN.sphere.y, EDGE.dir.y);
  const auto COLOR = program.evaluate(CENTER.dir, CTX);
  const auto CUT = program.evaluate(EDGE.dir, CTX);
  HS_EXPECT_EQ(COLOR.alpha, 1.0f);
  HS_EXPECT_EQ(CUT.alpha, 0.0f);
  HS_EXPECT_TRUE(carrier_live(COLOR));
}

/** @brief Spherical noise keeps its instance clock through edits and restore. */
inline void test_shader_chain_spherical_noise_lifecycle() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  using Model = In::Op::SampleSphericalNoise;
  const In::ChainEntryRequest CHAIN[] = {
      {"source", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  auto result = program.compile(CHAIN);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  param_as<Model::Params>(program, 0).noise_time_rate = 1.0f / 128.0f;
  for (int frame = 0; frame < 3; ++frame)
    program.advance();
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 3.0f / 128.0f);
  const auto CTX = shared_resources().context();
  program.prepare(CTX);
  SweepColors before;
  snapshot_render(program, CTX, before);
  const auto &runtime = program.ops()[0].op->runtime;
  const auto SAVED = runtime.capture_state(program.state_block(0));
  program.advance();
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 4.0f / 128.0f);
  HS_EXPECT_TRUE(runtime.restore_state(program.state_block(0), SAVED));
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 3.0f / 128.0f);
  program.prepare(CTX);
  SweepColors restored;
  snapshot_render(program, CTX, restored);
  expect_sweep_identical(restored, before);
  const auto SEED = state_as<Model::State>(program, 0).noise_seed;
  const In::ChainEntryRequest EDITED[] = {
      {"camera", "sphere.rotate.v2"},
      {"source", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  result = program.compile(EDITED);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(state_as<Model::State>(program, 1).phase, 3.0f / 128.0f);
  HS_EXPECT_EQ(state_as<Model::State>(program, 1).noise_seed, SEED);
  const In::ChainEntryRequest RENAMED[] = {
      {"other", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  result = program.compile(RENAMED);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 0.0f);
  HS_EXPECT_TRUE(state_as<Model::State>(program, 0).noise_seed != SEED);
}
