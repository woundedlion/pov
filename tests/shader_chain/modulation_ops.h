/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file modulation_ops.h
 * @brief Ordered color modulation and retained endpoint checks. */

/** @brief Repeated field stages own separate caches and survive reordering. */
inline void test_shader_chain_modulation_instances() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  const In::ChainEntryRequest chain[] = {
      {"sample", "sample.spherical-noise.v3"},
      {"colorize", "colorize.generated-palette.v3"},
      {"first", "modulate.field.v1"},
      {"second", "modulate.field.v1"},
  };
  HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
  auto &palette = param_as<In::Op::GeneratedPaletteParams>(program, 1);
  palette.modulation.enabled = false;
  auto &first = param_as<In::Op::FieldModulationParams>(program, 2).modulation;
  auto &second = param_as<In::Op::FieldModulationParams>(program, 3).modulation;
  first.seed = 1234;
  first.noise.speed = .00025f;
  first.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, .6f};
  second.seed = 4321;
  second.noise.speed = -.0005f;
  second.mechanic = {ColorModulation::Mechanic::HUE_ROTATION, .4f};
  program.advance();
  const auto ctx = shared_resources().context();
  program.prepare(ctx);
  const auto &state_a = state_as<In::Op::FieldModulationState>(program, 2);
  const auto &state_b = state_as<In::Op::FieldModulationState>(program, 3);
  HS_EXPECT_NE(state_a.resources.get(), state_b.resources.get());
  HS_EXPECT_NE(state_a.resources->prepared.cube.data,
               state_b.resources->prepared.cube.data);
  HS_EXPECT_EQ(state_a.clock.modulation.phase.phase_turns, .00025f);
  HS_EXPECT_EQ(state_b.clock.modulation.phase.phase_turns, .9995f);
  PB::Modulation::Resources<> oracle_a, oracle_b;
  const auto &source_a = oracle_a.prepare(first, state_a.clock.modulation);
  const auto &source_b = oracle_b.prepare(second, state_b.clock.modulation);
  size_t order_changes = 0;
  for (const auto &view : sweep_views()) {
    const auto field =
        run_entry<PB::FieldSample>(program, ctx, 0, PB::SphereSample{view, 0});
    const auto direction = Pattern::SphereDirection::admit(field.sphere);
    const Color4 base = ctx.palettes[0]->get(field.value);
    const float scalar_a = PB::Modulation::sample(source_a, direction);
    const float scalar_b = PB::Modulation::sample(source_b, direction);
    const Color4 expected = ColorModulation::apply(
        ColorModulation::apply(base, scalar_a, first.mechanic), scalar_b,
        second.mechanic);
    const auto actual = program.evaluate(view, ctx);
    HS_EXPECT_TRUE(color4_identical(actual, expected));
    const auto reversed = ColorModulation::apply(
        ColorModulation::apply(base, scalar_b, second.mechanic), scalar_a,
        first.mechanic);
    order_changes += !color4_identical(expected, reversed);
  }
  HS_EXPECT_GT(order_changes, size_t{0});
  const auto clock_a = state_a.clock;
  const auto clock_b = state_b.clock;
  const auto first_config = first;
  const auto second_config = second;
  const In::ChainEntryRequest reordered[] = {chain[0], chain[1], chain[3],
                                             chain[2]};
  const auto reordered_status =
      program.compile(reordered, [&](auto ops, auto *base) {
        reinterpret_cast<In::Op::GeneratedPaletteParams *>(base +
                                                           ops[1].param_offset)
            ->modulation.enabled = false;
        reinterpret_cast<In::Op::FieldModulationParams *>(base +
                                                          ops[2].param_offset)
            ->modulation = second_config;
        reinterpret_cast<In::Op::FieldModulationParams *>(base +
                                                          ops[3].param_offset)
            ->modulation = first_config;
        return true;
      });
  HS_EXPECT_EQ(reordered_status.code, In::ChainStatus::OK);
  HS_EXPECT_EQ(
      state_as<In::Op::FieldModulationState>(program, 2).clock.modulation,
      clock_b.modulation);
  HS_EXPECT_EQ(
      state_as<In::Op::FieldModulationState>(program, 3).clock.modulation,
      clock_a.modulation);
  HS_EXPECT_EQ(param_as<In::Op::FieldModulationParams>(program, 2).modulation,
               second_config);
  HS_EXPECT_EQ(param_as<In::Op::FieldModulationParams>(program, 3).modulation,
               first_config);
  param_as<In::Op::FieldModulationParams>(program, 2).modulation.enabled =
      false;
  param_as<In::Op::FieldModulationParams>(program, 3).modulation.enabled =
      false;
  program.prepare(ctx);
  for (const auto &view : sweep_views()) {
    const auto field =
        run_entry<PB::FieldSample>(program, ctx, 0, PB::SphereSample{view, 0});
    HS_EXPECT_TRUE(color4_identical(program.evaluate(view, ctx),
                                    ctx.palettes[0]->get(field.value)));
  }
}

/** @brief Columns blend completed nonlinear color tails once. */
inline void test_shader_chain_modulation_columns() {
  constexpr size_t BLOCKS = In::CHAIN_ARENA_BYTES / sizeof(std::max_align_t);
  std::vector<std::max_align_t> history(BLOCKS), column(BLOCKS);
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  program.bind_temporal_storage(reinterpret_cast<uint8_t *>(history.data()),
                                reinterpret_cast<uint8_t *>(column.data()));
  const In::ChainEntryRequest chain[] = {
      {"displace", "sphere.displace.curl.v2"},
      {"sample", "sample.spherical-noise.v3"},
      {"colorize", "colorize.generated-palette.v3"},
      {"hue", "modulate.hue-shift.v1"},
      {"first", "modulate.field.v1"},
      {"second", "modulate.field.v1"},
  };
  HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
  param_as<In::Op::CurlDisplaceParams>(program, 0).strength = .4f;
  auto &palette = param_as<In::Op::GeneratedPaletteParams>(program, 2);
  palette.modulation.enabled = false;
  param_as<In::Op::HueShiftModulationParams>(program, 3).hue_shift_amount = .7f;
  auto &a = param_as<In::Op::FieldModulationParams>(program, 4).modulation;
  auto &b = param_as<In::Op::FieldModulationParams>(program, 5).modulation;
  for (auto *config : {&a, &b}) {
    config->source = Pattern::SourceKind::DIRECTIONAL_SINE;
    config->evaluation = Pattern::Evaluation::DIRECT;
    config->sine.frequency = .4f;
  }
  a.sine.speed = .125f;
  a.sine.offset = .13f;
  a.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, .4f};
  b.sine.speed = -.0625f;
  b.sine.axis = math::X_AXIS;
  b.mechanic = {ColorModulation::Mechanic::OPACITY, .8f};
  program.advance();
  auto previous_ctx = shared_resources().context();
  previous_ctx.hue_rotation_lut = nullptr;
  previous_ctx.palettes[0] = previous_ctx.palettes[1];
  program.prepare(previous_ctx);
  program.capture_history();
  const auto *old_a_owner =
      state_as<In::Op::FieldModulationState>(program, 4).resources.get();
  program.advance();
  const auto ctx = shared_resources().context();
  HS_EXPECT_TRUE(program.history_continuous());
  const auto scalar = [](const Pattern::SineParams &params,
                         const math::Vector &direction, float phase) {
    return math::sinf_0_2pi(
        math::TWO_PI_F *
        math::wrap_t(params.frequency *
                         math::dot(params.axis.normalized(), direction) +
                     math::wrap_t(phase + params.offset)));
  };
  size_t nonzero_path = 0;
  for (float fraction : {0.0f, .375f, 1.0f}) {
    program.prepare_column(ctx, previous_ctx, fraction);
    const auto &current_a =
        *reinterpret_cast<const In::Op::ModulateField::Prepared *>(
            program.prepared_block(4));
    HS_EXPECT_EQ(current_a.source, &old_a_owner->prepared);
    for (const auto &view : sweep_views()) {
      const auto sphere = run_entry<PB::SphereSample>(
          program, ctx, 0, PB::SphereSample{view, 0});
      const auto field = run_entry<PB::FieldSample>(program, ctx, 1, sphere);
      const auto direction =
          Pattern::SphereDirection::admit(field.sphere).value();
      nonzero_path += field.path_length > 0;
      const auto endpoint = [&](Color4 input, float a_phase, float b_phase) {
        input = ColorModulation::HueRotation::apply(input, field.path_length,
                                                    {.7f});
        input = ColorModulation::apply(
            input, scalar(a.sine, direction, a_phase), a.mechanic);
        return ColorModulation::apply(input, scalar(b.sine, direction, b_phase),
                                      b.mechanic);
      };
      const auto before =
          endpoint(previous_ctx.palettes[0]->get(field.value), .125f, .9375f);
      const auto after =
          endpoint(ctx.palettes[0]->get(field.value), .25f, .875f);
      const auto expected = before.lerp(after, fraction);
      const auto actual =
          program.evaluate_column(view, ctx, previous_ctx, fraction);
      HS_EXPECT_NEAR(static_cast<float>(actual.color.r),
                     static_cast<float>(expected.color.r), 2.0f);
      HS_EXPECT_NEAR(static_cast<float>(actual.color.g),
                     static_cast<float>(expected.color.g), 2.0f);
      HS_EXPECT_NEAR(static_cast<float>(actual.color.b),
                     static_cast<float>(expected.color.b), 2.0f);
      HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
    }
  }
  HS_EXPECT_GT(nonzero_path, size_t{0});
}

/** @brief Snapshots preserve each field recipe and its independent source clock. */
inline void test_shader_chain_modulation_stage_snapshots() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"sample", "sample.spherical-noise.v3"},
      {"colorize", "colorize.generated-palette.v3"},
      {"first", "modulate.field.v1"},
      {"hue", "modulate.hue-shift.v1"},
      {"second", "modulate.field.v1"},
  };
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  auto &program = WB::program(effect);
  auto &a = param_as<In::Op::FieldModulationParams>(program, 2).modulation;
  auto &b = param_as<In::Op::FieldModulationParams>(program, 4).modulation;
  a.noise.speed = .00025f;
  a.mechanic.amount = .4f;
  b.source = Pattern::SourceKind::SPHERE_SMOOTH_GRID;
  b.evaluation = Pattern::Evaluation::DIRECT;
  b.grid.speed = .125f;
  b.mechanic = {ColorModulation::Mechanic::OPACITY, .7f};
  for (int frame = 0; frame < 3; ++frame)
    WB::advance_without_render(effect);
  const auto saved = effect.snapshot();
  HS_EXPECT_SIZE_OR_RETURN(saved.modulation, 3);
  HS_EXPECT_EQ(saved.modulation[1].instance, std::string("first"));
  HS_EXPECT_EQ(saved.modulation[2].instance, std::string("second"));
  const auto recipe_a = a;
  const auto recipe_b = b;
  const auto state_a = state_as<In::Op::FieldModulationState>(program, 2).clock;
  const auto state_b = state_as<In::Op::FieldModulationState>(program, 4).clock;
  a.noise.scale = 3;
  b.grid.frequency = 3;
  WB::advance_without_render(effect);
  HS_EXPECT_EQ(effect.restore_snapshot(saved),
               ChainSnapshotRestoreResult::APPLIED);
  HS_EXPECT_EQ(param_as<In::Op::FieldModulationParams>(program, 2).modulation,
               recipe_a);
  HS_EXPECT_EQ(param_as<In::Op::FieldModulationParams>(program, 4).modulation,
               recipe_b);
  HS_EXPECT_EQ(
      state_as<In::Op::FieldModulationState>(program, 2).clock.modulation,
      state_a.modulation);
  HS_EXPECT_EQ(
      state_as<In::Op::FieldModulationState>(program, 4).clock.modulation,
      state_b.modulation);
  auto invalid = saved;
  for (auto &entry : *invalid.runtime)
    if (entry.instance == "second")
      std::get<In::ColorClockSnapshot>(entry.state).modulation.source =
          Pattern::SourceKind::NOISE;
  HS_EXPECT_EQ(effect.restore_snapshot(invalid),
               ChainSnapshotRestoreResult::INVALID_VALUE);
}
