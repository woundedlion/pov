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
  const auto oracle_noise = make_oracle_noise("source", Model::ID, 1.0f);
  for (uint8_t basis = 0; basis < 3; ++basis) {
    params.basis = basis;
    const auto prepared = Model::prepare({}, params, state);
    const auto typed = Fixture::Stage::template prepare<Binding>(frame);
    const auto oracle = Pattern::SphereNoise::prepare(
        {3.5f, -1.0f / 128.0f, static_cast<math::NoiseBasis>(basis), 1.7f,
         Pattern::NoiseVariant::OCTAVE_CONTOUR},
        {0.37f}, oracle_noise);
    for (const auto &view : sweep_views()) {
      const float unit_raw = Pattern::SphereNoise::sample(
          Pattern::SphereDirection::from(view), oracle);
      for (float scale : {1.0f, std::sqrt(1.003f)}) {
        const PB::SphereSample input{scale * view, 0.375f};
        const auto actual = Model::run(input, {}, params, prepared);
        const auto composed =
            Fixture::Stage::template run<Binding>(input, frame, typed);
        const float raw = Pattern::SphereNoise::sample(
            Pattern::SphereDirection::from(input.dir), oracle);
        HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * raw, 1.0e-7f);
        HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * unit_raw, 2.0e-5f);
        HS_EXPECT_TRUE(field_identical(actual, composed));
        HS_EXPECT_EQ(actual.coverage, 1.0f);
        HS_EXPECT_EQ(actual.path_length, 0.375f);
        HS_EXPECT_EQ(actual.sphere.x, input.dir.x);
        HS_EXPECT_EQ(actual.sphere.y, input.dir.y);
        HS_EXPECT_EQ(actual.sphere.z, input.dir.z);
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
      const auto prepared = Model::prepare({}, params, {});
      const auto typed = Fixture::Stage::template prepare<Binding>(frame);
      for (const auto &view : sweep_views()) {
        const float unit_raw = Pattern::SphereCells::sample(
            Pattern::SphereDirection::from(view), prepared);
        for (float scale : {1.0f, std::sqrt(1.003f)}) {
          const PB::SphereSample input{scale * view, 1.25f};
          const auto actual = Model::run(input, {}, params, prepared);
          const auto composed =
              Fixture::Stage::template run<Binding>(input, frame, typed);
          HS_EXPECT_TRUE(field_identical(actual, composed));
          const float raw = Pattern::SphereCells::sample(
              Pattern::SphereDirection::from(input.dir), prepared);
          HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * raw, 1.0e-7f);
          if (output !=
              static_cast<uint8_t>(Pattern::SphereCellOutput::CELL_VALUE))
            HS_EXPECT_NEAR(actual.value, 0.5f + 0.5f * unit_raw, 3.0e-6f);
          HS_EXPECT_TRUE(actual.value >= 0.0f && actual.value <= 1.0f);
          HS_EXPECT_EQ(actual.coverage, 1.0f);
          HS_EXPECT_EQ(actual.path_length, 1.25f);
          HS_EXPECT_EQ(actual.sphere.x, input.dir.x);
          HS_EXPECT_EQ(actual.sphere.y, input.dir.y);
          HS_EXPECT_EQ(actual.sphere.z, input.dir.z);
        }
      }
    }
  }
  params.layout = static_cast<uint8_t>(Pattern::SphereCellLayout::OCTAHEDRON);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::CENTERS);
  const PB::SphereSample INPUT{math::X_AXIS, 2.0f};
  const auto prepared = Model::prepare({}, params, {});
  HS_EXPECT_NEAR(Model::run(INPUT, {}, params, prepared).value, 1.0f, 1.0e-6f);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::BORDERS);
  params.line_thickness = 0.1f;
  params.line_softness = 0.2f;
  const float angle =
      std::numbers::pi_v<float> / 4.0f -
      0.2f / std::sqrt(6.0f / (4.0f * std::numbers::pi_v<float>));
  const math::Vector border{std::cos(angle), std::sin(angle), 0.0f};
  const auto border_prepared = Model::prepare({}, params, {});
  const auto border_typed = Fixture::Stage::template prepare<Binding>(frame);
  for (float scale : {1.0f, std::sqrt(1.003f)}) {
    const PB::SphereSample drifted{scale * border, 2.0f};
    const auto actual = Model::run(drifted, {}, params, border_prepared);
    const auto composed =
        Fixture::Stage::template run<Binding>(drifted, frame, border_typed);
    HS_EXPECT_NEAR(actual.value, 0.5f, 2.0e-6f);
    HS_EXPECT_NEAR(composed.value, 0.5f, 2.0e-6f);
    HS_EXPECT_EQ(actual.sphere.x, drifted.dir.x);
    HS_EXPECT_EQ(composed.sphere.x, drifted.dir.x);
    HS_EXPECT_EQ(actual.path_length, 2.0f);
    HS_EXPECT_EQ(composed.path_length, 2.0f);
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
  const auto result = program.compile(CHAIN);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  auto &params = param_as<In::Op::SampleSphericalCells::Params>(program, 0);
  params.layout = static_cast<uint8_t>(Pattern::SphereCellLayout::OCTAHEDRON);
  params.output = static_cast<uint8_t>(Pattern::SphereCellOutput::CENTERS);
  auto &cutout = param_as<In::Op::ValueCutoutChainParams>(program, 1);
  cutout.cutout_threshold = 0.75f;
  cutout.cutout_softness = 0.0f;
  const auto ctx = shared_resources().context();
  program.prepare(ctx);
  const PB::SphereSample CENTER{math::X_AXIS, 0.25f};
  const PB::SphereSample edge{math::Vector(1, 1, 0).normalized(), 0.25f};
  const auto center_field = run_entry<PB::FieldSample>(program, ctx, 0, CENTER);
  const auto edge_field = run_entry<PB::FieldSample>(program, ctx, 0, edge);
  HS_EXPECT_NEAR(center_field.value, 1.0f, 1.0e-6f);
  HS_EXPECT_TRUE(edge_field.value < 0.75f);
  const auto visible =
      run_entry<PB::FieldSample>(program, ctx, 1, center_field);
  const auto hidden = run_entry<PB::FieldSample>(program, ctx, 1, edge_field);
  HS_EXPECT_EQ(visible.coverage, 1.0f);
  HS_EXPECT_EQ(hidden.coverage, 0.0f);
  HS_EXPECT_EQ(hidden.path_length, 0.25f);
  HS_EXPECT_EQ(hidden.sphere.x, edge.dir.x);
  HS_EXPECT_EQ(hidden.sphere.y, edge.dir.y);
  const auto color = program.evaluate(CENTER.dir, ctx);
  const auto cut = program.evaluate(edge.dir, ctx);
  HS_EXPECT_EQ(color.alpha, 1.0f);
  HS_EXPECT_EQ(cut.alpha, 0.0f);
  HS_EXPECT_TRUE(carrier_live(color));
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
  const auto ctx = shared_resources().context();
  program.prepare(ctx);
  SweepColors before;
  snapshot_render(program, ctx, before);
  const auto &runtime = program.ops()[0].op->runtime;
  const auto saved = runtime.capture_state(program.state_block(0));
  program.advance();
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 4.0f / 128.0f);
  HS_EXPECT_TRUE(runtime.restore_state(program.state_block(0), saved));
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 3.0f / 128.0f);
  program.prepare(ctx);
  SweepColors restored;
  snapshot_render(program, ctx, restored);
  expect_sweep_identical(restored, before);
  const auto seed = state_as<Model::State>(program, 0).noise_seed;
  const In::ChainEntryRequest EDITED[] = {
      {"camera", "sphere.rotate.v2"},
      {"source", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  result = program.compile(EDITED);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(state_as<Model::State>(program, 1).phase, 3.0f / 128.0f);
  HS_EXPECT_EQ(state_as<Model::State>(program, 1).noise_seed, seed);
  const In::ChainEntryRequest RENAMED[] = {
      {"other", Model::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  result = program.compile(RENAMED);
  HS_EXPECT_EQ(result.code, In::ChainStatus::OK);
  if (result.code != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(state_as<Model::State>(program, 0).phase, 0.0f);
  HS_EXPECT_TRUE(state_as<Model::State>(program, 0).noise_seed != seed);
}
