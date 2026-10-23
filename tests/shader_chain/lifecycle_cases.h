/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Transactional refusals, migration, determinism and budgets.

using SweepColors =
    std::array<Color4, std::tuple_size_v<decltype(sweep_views())>>;

/** Renders the sweep into @p out for a byte-identity comparison. */
inline void snapshot_render(In::ChainProgram &program,
                            const In::FrameContext &ctx, SweepColors &out) {
  const auto views = sweep_views();
  for (size_t index = 0; index < views.size(); ++index)
    out[index] = program.evaluate(views[index], ctx);
}

/** @brief Expects @p actual to match @p expected byte for byte. */
inline void expect_sweep_identical(const SweepColors &actual,
                                   const SweepColors &expected) {
  for (size_t index = 0; index < actual.size(); ++index) {
    HS_CONTEXT("view", static_cast<int>(index));
    HS_EXPECT_TRUE(color4_identical(actual[index], expected[index]));
  }
}

/** @brief Expects at least one swept color to be live. */
inline void expect_sweep_live(const SweepColors &colors) {
  HS_EXPECT_TRUE(std::any_of(colors.begin(), colors.end(),
                             [](const Color4 &c) { return carrier_live(c); }));
}

inline void expect_refusal(In::ChainProgram &program,
                           std::span<const In::ChainEntryRequest> request,
                           In::ChainStatus expected, int16_t expected_index,
                           const In::FrameContext &ctx,
                           const SweepColors &baseline) {
  const In::ChainRefusal refusal = program.compile(request);
  HS_EXPECT_EQ(static_cast<int>(refusal.code), static_cast<int>(expected));
  HS_EXPECT_EQ(refusal.entry_index, expected_index);
  HS_EXPECT_SIZE_OR_RETURN(program.ops(), 4);
  SweepColors after;
  snapshot_render(program, ctx, after);
  expect_sweep_identical(after, baseline);
}

inline void test_shader_chain_refusal_shape() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  if (!arm_default_chain(program, 2, ValueSet::MAXIMUMS))
    return;
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  SweepColors baseline;
  snapshot_render(program, ctx, baseline);
  expect_sweep_live(baseline);
  const auto camera_before = state_as<In::Op::SpatialWalkState>(program, 0);
  const auto source_before = state_as<In::Op::SourceClockState>(program, 2);
  const auto color_before = state_as<In::Op::ColorClockState>(program, 3);
  HS_EXPECT_GT(source_before.primary, 0.0f);

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

  // State continuity: the shape refusals left every clock untouched.
  const auto &camera = state_as<In::Op::SpatialWalkState>(program, 0);
  HS_EXPECT_EQ(std::memcmp(&camera.position, &camera_before.position,
                           sizeof(math::Vector)),
               0);
  HS_EXPECT_EQ(std::memcmp(&camera.direction, &camera_before.direction,
                           sizeof(math::Vector)),
               0);
  HS_EXPECT_EQ(std::memcmp(&camera.wander, &camera_before.wander,
                           sizeof(math::Quaternion)),
               0);
  HS_EXPECT_EQ(camera.angular_velocity, camera_before.angular_velocity);
  HS_EXPECT_EQ(camera.spin_phase, camera_before.spin_phase);
  HS_EXPECT_EQ(camera.walk_time, camera_before.walk_time);
  HS_EXPECT_EQ(camera.noise_seed, camera_before.noise_seed);
  const auto &source = state_as<In::Op::SourceClockState>(program, 2);
  HS_EXPECT_EQ(source.primary, source_before.primary);
  HS_EXPECT_EQ(source.secondary, source_before.secondary);
  HS_EXPECT_EQ(source.angle, source_before.angle);
  const auto &color = state_as<In::Op::ColorClockState>(program, 3);
  HS_EXPECT_EQ(color.hue_noise_seed, color_before.hue_noise_seed);
  HS_EXPECT_EQ(color.oscillation_phase, color_before.oscillation_phase);
  HS_EXPECT_EQ(color.hue_noise_phase, color_before.hue_noise_phase);
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
  if (!arm_default_chain(measured->program, 0, ValueSet::DEFAULTS))
    return;
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
  const In::FrameContext ctx = shared_resources().context();
  exact->program.prepare(ctx);
  SweepColors baseline;
  snapshot_render(exact->program, ctx, baseline);
  expect_sweep_live(baseline);
  const In::ChainRefusal rejected = exact->program.compile(longer);
  HS_EXPECT_EQ(static_cast<int>(rejected.code),
               static_cast<int>(In::ChainStatus::ARENA_OVERFLOW));
  HS_EXPECT_EQ(exact->program.ops().size(), 4u);
  SweepColors after;
  snapshot_render(exact->program, ctx, after);
  expect_sweep_identical(after, baseline);
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
  SweepColors baseline;
  snapshot_render(program, ctx, baseline);
  expect_sweep_live(baseline);

  CountLifecycle::fail_migrate = true;
  const int destroys_before = CountLifecycle::destroys;
  const In::ChainRefusal failed = program.compile(chain);
  HS_EXPECT_EQ(static_cast<int>(failed.code),
               static_cast<int>(In::ChainStatus::MIGRATE_FAILED));
  HS_EXPECT_EQ(failed.entry_index, 0);
  // The live program is untouched: no state destroyed, phase and render kept.
  HS_EXPECT_EQ(CountLifecycle::destroys, destroys_before);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  SweepColors after;
  snapshot_render(program, ctx, after);
  expect_sweep_identical(after, baseline);

  // A failing migrate deeper in the chain destroys only the candidate states
  // constructed before it: "fresh" inits at entry 0, "counter" fails at 1.
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
  // The fresh candidate was constructed then torn down; the live program
  // keeps its state.
  HS_EXPECT_EQ(CountLifecycle::inits, inits_before + 1);
  HS_EXPECT_EQ(CountLifecycle::migrates, migrates_before + 1);
  HS_EXPECT_EQ(CountLifecycle::destroys, deep_destroys_before + 1);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  SweepColors after_deep;
  snapshot_render(program, ctx, after_deep);
  expect_sweep_identical(after_deep, baseline);
  CountLifecycle::fail_migrate = false;
  program.clear();
}

/** @brief Sampled column states outlive column evaluation and are destroyed
    on reuse, history capture, recompile, clear and program destruction. */
inline void test_shader_chain_column_states_outlive_evaluation() {
  struct TemporalBlocks {
    alignas(std::max_align_t) uint8_t history[In::CHAIN_ARENA_BYTES];
    alignas(std::max_align_t) uint8_t column[In::CHAIN_ARENA_BYTES];
  };
  LifetimeRegistry::reset();
  auto blocks = std::make_unique<TemporalBlocks>();
  auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, lifetime_table());
  In::ChainProgram &program = fixture->program;
  program.bind_temporal_storage(blocks->history, blocks->column);
  const In::ChainEntryRequest chain[] = {
      {"life", "test.lifetime.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::FrameContext ctx = shared_resources().context();
  const auto views = sweep_views();
  const auto render_column = [&](float fraction) {
    program.prepare_column(ctx, ctx, fraction);
    for (const math::Vector &view : views)
      program.evaluate_column(view, ctx, ctx, fraction);
  };
  const auto arm = [&] {
    HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
    program.prepare(ctx);
    program.capture_history();
    program.advance();
    program.prepare(ctx);
    HS_EXPECT_TRUE(program.history_continuous());
  };

  arm();
  HS_EXPECT_EQ(LifetimeRegistry::live.size(), 2u);
  render_column(0.25f);
  HS_EXPECT_EQ(LifetimeRegistry::live.size(), 3u);
  render_column(0.75f);
  HS_EXPECT_EQ(LifetimeRegistry::live.size(), 3u);
  program.prepare(ctx);
  program.capture_history();
  HS_EXPECT_EQ(LifetimeRegistry::live.size(), 2u);

  program.advance();
  program.prepare(ctx);
  render_column(0.5f);
  HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
  HS_EXPECT_EQ(LifetimeRegistry::live.size(), 1u);

  arm();
  render_column(0.5f);
  program.clear();
  HS_EXPECT_EQ(LifetimeRegistry::live.size(), 0u);

  arm();
  render_column(0.5f);
  fixture.reset();
  HS_EXPECT_EQ(LifetimeRegistry::live.size(), 0u);
  HS_EXPECT_EQ(LifetimeRegistry::dead_reads, 0);
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
  const auto compile_status1 = program.compile(first).code;
  HS_EXPECT_EQ(compile_status1, In::ChainStatus::OK);
  if (compile_status1 != In::ChainStatus::OK)
    return;
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
  const auto compile_status2 = program.compile(second).code;
  HS_EXPECT_EQ(compile_status2, In::ChainStatus::OK);
  if (compile_status2 != In::ChainStatus::OK)
    return;
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
  const auto compile_status3 = program.compile(third).code;
  HS_EXPECT_EQ(compile_status3, In::ChainStatus::OK);
  if (compile_status3 != In::ChainStatus::OK)
    return;
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
  if (!arm_default_chain(program, 6, ValueSet::MAXIMUMS))
    return;
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
  const auto compile_status4 = program.compile(edited).code;
  HS_EXPECT_EQ(compile_status4, In::ChainStatus::OK);
  if (compile_status4 != In::ChainStatus::OK)
    return;
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
  if (!arm_default_chain(first->program, 5, ValueSet::MAXIMUMS))
    return;
  if (!arm_default_chain(second->program, 5, ValueSet::MAXIMUMS))
    return;
  const FastNoiseLite authored_walk_noise = make_oracle_noise(
      "camera", "sphere.rotate.v2", In::Op::WALK_OPTIONS.noise_scale);
  const auto &seeded_walk =
      state_as<In::Op::SpatialWalkState>(first->program, 0);
  HS_EXPECT_EQ(seeded_walk.walk_noise.GetNoise(0.25f, -0.5f, 0.75f),
               authored_walk_noise.GetNoise(0.25f, -0.5f, 0.75f));
  const In::FrameContext ctx = shared_resources().context();
  first->program.prepare(ctx);
  second->program.prepare(ctx);
  SweepColors first_colors;
  SweepColors second_colors;
  snapshot_render(first->program, ctx, first_colors);
  snapshot_render(second->program, ctx, second_colors);
  expect_sweep_live(first_colors);
  expect_sweep_identical(second_colors, first_colors);
  // Instance labels participate in stateful-resource identity.
  auto relabeled = std::make_unique<ProgramFixture>();
  const In::ChainEntryRequest renamed[] = {
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const auto compile_status5 = relabeled->program.compile(renamed).code;
  HS_EXPECT_EQ(compile_status5, In::ChainStatus::OK);
  if (compile_status5 != In::ChainStatus::OK)
    return;
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
  if (!arm_default_chain(program, 0, ValueSet::DEFAULTS))
    return;

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

  // The committed footprint equals the accounted layout: operator blocks,
  // the per-op overhead, and the fixed per-field name reservation.
  size_t expected_bytes = 0;
  const auto align_to = [](size_t offset, size_t alignment) {
    return (offset + alignment - 1) & ~(alignment - 1);
  };
  for (const In::ChainProgram::ChainOp &op : ops) {
    const In::OperatorRuntime &runtime = op.op->runtime;
    expected_bytes = align_to(expected_bytes, runtime.param.align);
    expected_bytes += runtime.param.size;
    expected_bytes = align_to(expected_bytes, std::max(runtime.prepared.align,
                                                       runtime.interval.align));
    expected_bytes += std::max(runtime.prepared.size, runtime.interval.size);
    expected_bytes = align_to(expected_bytes, runtime.state.align);
    expected_bytes += runtime.state.size;
    expected_bytes += In::PER_OP_OVERHEAD_BYTES;
    expected_bytes +=
        static_cast<size_t>(op.op->schema_count) * In::PER_PARAM_NAME_BYTES;
  }
  HS_EXPECT_EQ(program.used_bytes(), expected_bytes);
  program.clear();
}
