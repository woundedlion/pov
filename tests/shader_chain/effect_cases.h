/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Composed frame parity and ShaderChain effect lifecycle.

inline void test_shader_chain_composed_frame_parity() {
  using FX = AlienCore<96, 20>;
  reset_globals();
  auto params = FX::initial_params();
  params.template get<"projection">().camera_wander = 0.0f;
  std::array<FX::FrameState, 5> references;
  {
    FX composed;
    composed.init();
    ComposedFrameWhiteBox::set_params(composed, params);
    for (size_t frame = 0; frame < references.size(); ++frame) {
      pin_frame_clock(static_cast<int>(frame) + 1);
      ComposedFrameWhiteBox::advance(composed);
      references[frame] = ComposedFrameWhiteBox::frame(composed);
    }
  }

  reset_globals();
  ShaderChainWhiteBox::FX chain;
  chain.init();
  const In::ChainEntryRequest topology[] = {
      {"camera", "sphere.rotate.v2"},
      {"lens", "sphere.lens.glitch.v2"},
      {"project", "project.gnomonic.v2"},
      {"warp", "warp.mirror-tile.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const auto compile_status1 = chain.set_chain(topology).code;
  HS_EXPECT_EQ(compile_status1, In::ChainStatus::OK);
  if (compile_status1 != In::ChainStatus::OK)
    return;
  auto &program = ShaderChainWhiteBox::program(chain);
  param_as<In::Op::RotateChainParams>(program, 0).wander =
      params.template get<"projection">().camera_wander;
  auto &projection = param_as<In::Op::GnomonicChainParams>(program, 2);
  projection.frame = static_cast<uint8_t>(In::Op::ProjectionFrame::IDENTITY);
  projection.singularity_fade =
      params.template get<"projection">().singularity_fade;
  static_cast<PB::MirrorParams &>(param_as<In::Op::MirrorWarpParams>(
      program, 3)) = params.template get<"outer_warp">();
  auto &source = param_as<In::Op::GridSampleParams>(program, 4);
  static_cast<PB::GridSourceParams &>(source) = params.template get<"source">();
  source.edge_width = params.template get<"value">().edge_width;
  source.coverage_mode =
      static_cast<uint8_t>(PB::ProjectionCoverageMode::EDGE_FADE);
  auto &color = ShaderChainWhiteBox::color_params(chain);
  static_cast<PB::Color::ColorControls &>(color) =
      params.template get<"color">();
  color.mapping_mode =
      static_cast<uint8_t>(params.template get<"color">().palette_mapping);

  size_t visible = 0;
  for (int frame = 1; frame <= 5; ++frame) {
    HS_CONTEXT("frame", frame);
    pin_frame_clock(frame);
    chain.draw_frame();
    chain.advance_display();
    const In::FrameContext ctx = ShaderChainWhiteBox::frame_context(chain);
    auto reference = references[static_cast<size_t>(frame - 1)];
    reference.palette = ctx.palettes[color.palette_mode];
    reference.hue_rotation_lut = ctx.hue_rotation_lut;
    reference.hue_noise_lut = ctx.hue_noise_lut;
    const auto prepared = FX::RenderPipeline::prepare(reference);
    for (const math::Vector &view : sweep_views()) {
      const Color4 expected = FX::shade(view, prepared);
      const Color4 actual = program.evaluate(view, ctx);
      expect_color_within(actual, expected, 1, 1e-6f);
      visible += expected.alpha > 0.0f &&
                 (expected.color.r != 0 || expected.color.g != 0 ||
                  expected.color.b != 0);
    }
  }
  HS_EXPECT_GT(visible, 0u);
}

inline void test_shader_chain_effect_registers_params() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  size_t expected = 0;
  for (const In::ChainEntryRequest &entry : DEFAULT_CHAIN) {
    const auto *op = require_operator(entry.operator_id);
    if (!op)
      return;
    expected += op->schema_count;
  }
  const ParamList &params = effect.getParameters();
  HS_EXPECT_EQ(params.size(), expected);
  HS_EXPECT_TRUE(params.find("camera.wander") != nullptr);
  HS_EXPECT_TRUE(params.find("project.singularity-fade") != nullptr);
  HS_EXPECT_TRUE(params.find("sample.pattern-freq") != nullptr);
  HS_EXPECT_TRUE(params.find("colorize.palette-chroma") != nullptr);
  const ParamDef *coverage = params.find("sample.coverage-mode");
  HS_EXPECT_TRUE(coverage != nullptr);
  if (!coverage)
    return;
  HS_EXPECT_TRUE(coverage->is_enum());
  HS_EXPECT_EQ(coverage->option_count, 4);
  HS_EXPECT_TRUE(std::string_view(coverage->options[3]) == "edge-fade");
  HS_EXPECT_EQ(
      static_cast<int>(effect.updateParameter("sample.coverage-mode", 3.0f)),
      static_cast<int>(ParamSetResult::APPLIED));
  HS_EXPECT_EQ(param_as<In::Op::GridSampleParams>(
                   ShaderChainWhiteBox::program(effect), 2)
                   .coverage_mode,
               static_cast<uint8_t>(In::Op::ProjectionCoverageMode::EDGE_FADE));

  hs_test::run_frames(effect, 4);
  HS_EXPECT_GT((frame_energy<96, 20>(effect)), 0u);
}

inline void test_shader_chain_parameter_admission() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"lens", "sphere.lens.mobius.v2"},
      {"project", "project.stereographic.v2"},
      {"warp", "warp.curl-flow.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const auto compile_status2 = effect.set_chain(chain).code;
  HS_EXPECT_EQ(compile_status2, In::ChainStatus::OK);
  if (compile_status2 != In::ChainStatus::OK)
    return;
  const ParamDef *strength_param = effect.getParameters().find("warp.strength");
  HS_EXPECT_TRUE(strength_param != nullptr);
  if (!strength_param)
    return;
  const ParamDef &strength = *strength_param;
  HS_EXPECT_EQ(effect.updateParameter("warp.strength", 30.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.parameter_warning("warp.strength") == nullptr);
  HS_EXPECT_EQ(strength.get_requested(), strength.max);
  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_EQ(effect.accepted_parameter_value(strength), strength.max);
  effect.draw_frame();
  effect.advance_display();
  effect.updateParameter("warp.strength", 0.0f);
  HS_EXPECT_TRUE(effect.parameter_warning("warp.strength") == nullptr);
  HS_EXPECT_EQ(effect.accepted_parameter_value(strength), 0.0f);
  const ParamDef *d_re_param = effect.getParameters().find("lens.mobius-d-re");
  HS_EXPECT_TRUE(d_re_param != nullptr);
  if (!d_re_param)
    return;
  const ParamDef &d_re = *d_re_param;
  const float accepted_d_re = effect.accepted_parameter_value(d_re);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-a-re", 2.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-d-re", 0.0f),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-d-re") != nullptr);
  HS_EXPECT_EQ(d_re.get_requested(), accepted_d_re);
  HS_EXPECT_EQ(effect.accepted_parameter_value(d_re), accepted_d_re);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-a-re", 3.0f),
               ParamSetResult::APPLIED);
  const float saved_a = requested_param(effect, "lens.mobius-a-re");
  const float saved_d = d_re.get_requested();
  const auto compile_status3 = effect.set_chain(chain).code;
  HS_EXPECT_EQ(compile_status3, In::ChainStatus::OK);
  if (compile_status3 != In::ChainStatus::OK)
    return;
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-a-re", saved_a),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-d-re", saved_d),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(requested_param(effect, "lens.mobius-a-re"), saved_a);
  HS_EXPECT_EQ(requested_param(effect, "lens.mobius-d-re"), saved_d);
  effect.draw_frame();
  effect.advance_display();
  effect.updateParameter("lens.mobius-d-re", 1.0f);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-d-re") == nullptr);
  const ShaderChainParameterWrite valid_batch[] = {{"lens.mobius-a-re", 0.0f},
                                                   {"lens.mobius-b-re", 1.0f},
                                                   {"lens.mobius-c-re", 1.0f},
                                                   {"lens.mobius-d-re", 0.0f}};
  HS_EXPECT_EQ(effect.update_parameters(valid_batch), ParamSetResult::APPLIED);
  for (const auto &write : valid_batch)
    HS_EXPECT_EQ(requested_param(effect, write.name), write.value);
  const ShaderChainParameterWrite invalid_batch[] = {
      {"lens.mobius-b-re", 0.0f}, {"lens.mobius-c-re", 0.0f}};
  HS_EXPECT_EQ(effect.update_parameters(invalid_batch),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-b-re") != nullptr);
  HS_EXPECT_TRUE(effect.parameter_warning("lens.mobius-a-re") == nullptr);
  const ShaderChainParameterWrite unknown_batch[] = {{"lens.mobius-a-re", 1.0f},
                                                     {"missing", 0.0f}};
  HS_EXPECT_EQ(effect.update_parameters(unknown_batch),
               ParamSetResult::UNKNOWN_PARAM);
  const ShaderChainParameterWrite nonfinite_batch[] = {
      {"lens.mobius-a-re", 1.0f}, {"lens.mobius-c-re", NAN}};
  HS_EXPECT_EQ(effect.update_parameters(nonfinite_batch),
               ParamSetResult::NON_FINITE);
  for (const auto &write : valid_batch)
    HS_EXPECT_EQ(requested_param(effect, write.name), write.value);
  effect.draw_frame();
  effect.advance_display();
}

inline void test_shader_chain_edge_distance_admission() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"project", "project.folded-sinusoidal.v2"},
      {"warp", "warp.wave-shear.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  HS_EXPECT_EQ(effect.updateParameter("sample.coverage-mode", 3),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning("sample.coverage-mode") != nullptr);
  const ParamDef *envelope = effect.getParameters().find("warp.envelope");
  const ParamDef *coverage_mode =
      effect.getParameters().find("sample.coverage-mode");
  HS_EXPECT_TRUE(envelope != nullptr && coverage_mode != nullptr);
  if (!envelope || !coverage_mode)
    return;
  const float ACCEPTED_ENVELOPE = envelope->get_requested();
  const ShaderChainParameterWrite writes[] = {{"warp.envelope", 2}};
  HS_EXPECT_EQ(effect.update_parameters(writes), ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(requested_param(effect, "warp.envelope"), ACCEPTED_ENVELOPE);
  HS_EXPECT_TRUE(effect.parameter_warning("warp.envelope") != nullptr);
  HS_EXPECT_EQ(requested_param(effect, "sample.coverage-mode"), 1.0f);
}

/** @brief Projections stay within their declared plane bound. */
inline void test_shader_chain_projection_plane_bound() {
  std::vector<math::Vector> views;
  for (const math::Vector &view : sweep_views())
    views.push_back(view);
  for (const math::Vector &view :
       {math::Vector(1e-4f, 1, 0), math::Vector(0, 1, 1e-4f),
        math::Vector(1, 1e-6f, 0), math::Vector(0, -1e-6f, 1)})
    views.push_back(view.normalized());
  const float NEAR_POLE_Y = std::nextafter(1.0f, 0.0f);
  views.emplace_back(sqrtf(1.0f + PB::DIRECTION_NORM_SQUARED_TOLERANCE -
                           NEAR_POLE_Y * NEAR_POLE_Y),
                     NEAR_POLE_Y, 0.0f);
  constexpr int SWEEP = 4096;
  for (int index = 0; index < SWEEP; ++index) {
    const float y = 1.0f - 2.0f * (index + 0.5f) / SWEEP;
    const float radius = sqrtf(fmaxf(0.0f, 1.0f - y * y));
    const float angle = 2.39996323f * static_cast<float>(index);
    views.push_back(
        math::Vector(radius * cosf(angle), y, radius * sinf(angle)));
  }
  size_t projections = 0;
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    if (op.input != In::CarrierId::SPHERE || op.output != In::CarrierId::PLANE)
      continue;
    HS_CONTEXT(op.operator_id);
    ++projections;
    auto fixture = std::make_unique<ProgramFixture>();
    auto &program = fixture->program;
    const In::ChainEntryRequest chain[] = {
        {"project", op.operator_id},
        {"sample", In::Op::SampleGridV3::ID},
        {"color", In::Op::ColorizeGeneratedPaletteV3::ID},
    };
    const auto compile_status4 = program.compile(chain).code;
    HS_EXPECT_EQ(compile_status4, In::ChainStatus::OK);
    if (compile_status4 != In::ChainStatus::OK)
      return;
    param_as<In::Op::ProjectChainParams>(program, 0).frame =
        static_cast<uint8_t>(In::Op::ProjectionFrame::IDENTITY);
    const auto ctx = shared_resources().context();
    program.prepare(ctx);
    HS_EXPECT_TRUE(op.runtime.plane_bound != nullptr);
    if (op.runtime.plane_bound == nullptr)
      continue;
    const float bound = op.runtime.plane_bound(program.param_block(0), 0.0f);
    for (const math::Vector &view : views) {
      const PB::PlaneSample plane = run_entry<PB::PlaneSample>(
          program, ctx, 0, PB::SphereSample{view, 0.0f});
      HS_EXPECT_LE(std::hypot(plane.coords.re, plane.coords.im), bound);
    }
  }
  HS_EXPECT_GT(projections, 0u);
}

/** @brief A chain whose combined plane growth could overflow is refused on
    every admission path; ordinary extreme compositions still admit. */
inline void test_shader_chain_plane_growth_admission() {
  using WB = ShaderChainWhiteBox;
  constexpr int WARPS = 23;
  constexpr float MIN_SCALE = 1.0f / 64.0f;
  std::array<std::string, WARPS> ids;
  std::array<std::string, WARPS> scale_x;
  std::array<std::string, WARPS> scale_y;
  std::array<In::ChainEntryRequest, WARPS + 3> chain{};
  chain[0] = {"project", In::Op::ProjectStereographic::ID};
  for (int index = 0; index < WARPS; ++index) {
    ids[index] = "warp" + std::to_string(index);
    scale_x[index] = ids[index] + ".scale-x";
    scale_y[index] = ids[index] + ".scale-y";
    chain[index + 1] = {ids[index], In::Op::WarpAffineV3::ID};
  }
  chain[WARPS + 1] = {"sample", In::Op::SampleGridV3::ID};
  chain[WARPS + 2] = {"colorize", In::Op::ColorizeGeneratedPaletteV3::ID};

  reset_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  std::vector<ShaderChainParameterWrite> all;
  for (int index = 0; index < WARPS; ++index) {
    all.push_back({scale_x[index].c_str(), MIN_SCALE});
    all.push_back({scale_y[index].c_str(), MIN_SCALE});
  }
  HS_EXPECT_EQ(effect.update_parameters(all), ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning(scale_x[0].c_str()) ==
                 In::PLANE_GROWTH_WARNING);
  HS_EXPECT_EQ(requested_param(effect, scale_x[0].c_str()), 1.0f);

  int admitted = 0;
  for (; admitted < WARPS; ++admitted) {
    const ShaderChainParameterWrite step[] = {
        {scale_x[admitted].c_str(), MIN_SCALE},
        {scale_y[admitted].c_str(), MIN_SCALE}};
    if (effect.update_parameters(step) != ParamSetResult::APPLIED)
      break;
  }
  HS_EXPECT_EQ(admitted, 4);
  HS_EXPECT_TRUE(effect.parameter_warning(scale_x[4].c_str()) ==
                 In::PLANE_GROWTH_WARNING);
  HS_EXPECT_EQ(effect.updateParameter(scale_y[4].c_str(), MIN_SCALE),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning(scale_y[4].c_str()) ==
                 In::PLANE_GROWTH_WARNING);
  HS_EXPECT_EQ(requested_param(effect, scale_y[4].c_str()), 1.0f);
  HS_EXPECT_EQ(effect.updateParameter(scale_x[4].c_str(), 1.5f),
               ParamSetResult::APPLIED);

  const ChainSnapshot live = effect.snapshot();
  ChainSnapshot overflowing = live;
  for (auto &parameter : overflowing.parameters)
    if (parameter.name == scale_y[4])
      parameter.value = MIN_SCALE;
  HS_EXPECT_EQ(effect.restore_snapshot(overflowing),
               ChainSnapshotRestoreResult::INVALID_VALUE);
  HS_EXPECT_EQ(requested_param(effect, scale_y[4].c_str()), 1.0f);
  HS_EXPECT_EQ(effect.restore_snapshot(live),
               ChainSnapshotRestoreResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  auto &program = WB::program(effect);
  const In::FrameContext ctx = WB::frame_context(effect);
  for (const math::Vector &view : sweep_views()) {
    const PB::PlaneSample plane =
        run_to_plane(program, ctx, PB::SphereSample{view, 0.0f}, WARPS + 1);
    HS_EXPECT_LE(std::hypot(plane.coords.re, plane.coords.im),
                 In::MAX_PLANE_BOUND);
    HS_EXPECT_TRUE(std::isfinite(plane.path_length));
  }

  const In::ChainEntryRequest ordinary[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.gnomonic.v2"},
      {"warp", "warp.affine.v3"},
      {"polar", "warp.polar-chart.v2"},
      {"sample", "sample.lattice.v2"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(effect.set_chain(ordinary).code, In::ChainStatus::OK);
  const ShaderChainParameterWrite extreme[] = {
      {"warp.scale-x", MIN_SCALE},    {"warp.scale-y", 64.0f},
      {"warp.shear", -4.0f},          {"warp.translation-x", 4.0f},
      {"warp.translation-y", -4.0f},  {"warp.lattice-period", 100.0f},
      {"polar.radial-scale", 64.0f},  {"polar.radial-phase", 6.0f},
      {"polar.angular-phase", -6.0f}, {"polar.harmonic", 15.0f}};
  HS_EXPECT_EQ(effect.update_parameters(extreme), ParamSetResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
}

inline void test_shader_chain_effect_rebind_generation() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const uint32_t initial = effect.getParameterSchemaGeneration();
  const In::ChainEntryRequest edited[] = {
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal committed = effect.set_chain(edited);
  HS_EXPECT_EQ(static_cast<int>(committed.code),
               static_cast<int>(In::ChainStatus::OK));
  // set_chain rebinds before returning: fresh generation, fresh names.
  HS_EXPECT_NE(effect.getParameterSchemaGeneration(), initial);
  HS_EXPECT_TRUE(effect.getParameters().find("camera2.wander") != nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("camera.wander") == nullptr);
  effect.draw_frame();
  effect.advance_display();
}

inline void test_shader_chain_effect_refusal_keeps_schema() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const uint32_t committed = effect.getParameterSchemaGeneration();
  const size_t param_count = effect.getParameters().size();
  const In::ChainEntryRequest unknown[] = {
      {"camera", "sphere.rotate.v2"},
      {"warp", "warp.unknown.v2"},
  };
  const In::ChainRefusal refusal = effect.set_chain(unknown);
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::UNKNOWN_OPERATOR));
  HS_EXPECT_EQ(refusal.entry_index, 1);
  // Definitions, generation, and the committed program survive, and the
  // effect still renders.
  HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), committed);
  HS_EXPECT_EQ(effect.getParameters().size(), param_count);
  HS_EXPECT_TRUE(effect.getParameters().find("camera.wander") != nullptr);
  hs_test::run_frames(effect, 2);
  HS_EXPECT_GT((frame_energy<96, 20>(effect)), 0u);
}

/** @brief Pause gates preset selection only: chain clocks and the visible
    palette keep advancing. */
inline void test_shader_chain_pause_semantics() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(static_cast<int>(effect.updateParameter("sample.speed", 0.05f)),
               static_cast<int>(ParamSetResult::APPLIED));
  HS_EXPECT_EQ(static_cast<int>(effect.updateParameter("camera.wander", 1.0f)),
               static_cast<int>(ParamSetResult::APPLIED));
  effect.setAnimationsPaused(true);
  const In::Op::SourceClockState source_before =
      state_as<In::Op::SourceClockState>(WB::program(effect), 2);
  const math::Quaternion wander_before =
      state_as<In::Op::SpatialWalkState>(WB::program(effect), 0).wander;
  const Pixel color_before = WB::palette_color(effect, 0.25f);

  hs_test::run_frames(effect, 60);

  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_NE(
      state_as<In::Op::SourceClockState>(WB::program(effect), 2).primary,
      source_before.primary);
  HS_EXPECT_TRUE(
      state_as<In::Op::SpatialWalkState>(WB::program(effect), 0).wander !=
      wander_before);
  const Pixel color_after = WB::palette_color(effect, 0.25f);
  HS_EXPECT_TRUE(color_after != color_before);
  HS_EXPECT_GT((frame_energy<96, 20>(effect)), 0u);
}

/** @brief ShaderChain with temporal sampling disabled matches its committed
 *  render over two programs, clips and a hue-noise colorize. */
inline void test_shader_chain_instantaneous_pixel_golden() {
  using WB = ShaderChainWhiteBox;
  constexpr int W = 96;
  constexpr int H = 20;
  reset_globals();
  const ScopedTemporalSampling instantaneous(false);
  HS_EXPECT_TRUE(Spacetime::flat_surface());
  WB::FX chain;
  chain.init();
  const In::ChainEntryRequest topology[] = {
      {"camera", "sphere.rotate.v2"},
      {"lens", "sphere.lens.glitch.v2"},
      {"project", "project.gnomonic.v2"},
      {"warp", "warp.mirror-tile.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  uint64_t hash = HASH_SEED64;
  for (int frame = 0; frame < 100; ++frame) {
    pin_frame_clock(frame);
    if (frame == 20)
      chain.set_clip(2, H - 2, 40, 60);
    if (frame == 30) {
      chain.set_clip(0, H, 0, 8);
      chain.set_margin(4);
    }
    if (frame == 40) {
      chain.set_margin(0);
      chain.set_clip(0, H, 0, W);
    }
    if (frame == 60) {
      HS_EXPECT_EQ(chain.set_chain(topology).code, In::ChainStatus::OK);
      WB::color_params(chain).hue_shift_amount = 0.3f;
    }
    chain.draw_frame();
    chain.advance_display();
    hash = frame_hash<W, H>(chain, hash);
  }
  hs::clear_mock_time();
  expect_pixel_golden("ShaderChain 96x20", hash, 0xbfbf4393aff25990ull);
}

/** @brief View direction the scan shades for a 96x20 pixel.
 * @param x Column.
 * @param y Row.
 * @return Unit view direction. */
inline math::Vector chain_scan_view(int x, int y) {
  const float sin_phi = math::TrigLUT<96, 20>::sin_phi[y];
  return math::Vector(sin_phi * math::TrigLUT<96, 20>::cos_theta(x),
                      math::TrigLUT<96, 20>::cos_phi[y],
                      sin_phi * math::TrigLUT<96, 20>::sin_theta[x]);
}

/** Camera spin of the column-time ShaderChain fixtures, radians per frame. */
inline constexpr float CHAIN_SPIN = 0.05f;

/** @brief Spins the default chain's camera over a fine grid.
 * @param chain Initialised effect. */
inline void spin_default_chain(ShaderChainWhiteBox::FX &chain) {
  In::ChainProgram &program = ShaderChainWhiteBox::program(chain);
  param_as<In::Op::RotateChainParams>(program, 0).spin_rate = CHAIN_SPIN;
  param_as<In::Op::GridSampleParams>(program, 2).pattern_freq = 8.0f;
}

/** @brief Pixel the drawn chain shades with its camera at one spin angle and
 *  the current colour resources.
 * @param chain Effect after a draw.
 * @param view Unit view direction.
 * @param angle Camera spin, radians.
 * @return Premultiplied pixel. */
inline Pixel chain_pixel_at_spin(ShaderChainWhiteBox::FX &chain,
                                 const math::Vector &view, float angle) {
  using WB = ShaderChainWhiteBox;
  In::ChainProgram &program = WB::program(chain);
  const In::FrameContext ctx = WB::frame_context(chain);
  auto *camera = const_cast<In::Op::Rotate::Prepared *>(
      reinterpret_cast<const In::Op::Rotate::Prepared *>(
          program.prepared_block(0)));
  const In::Op::Rotate::Prepared saved = *camera;
  camera->conjugate = math::make_rotation(math::Y_AXIS, angle).conjugate();
  const Color4 color = program.evaluate(view, ctx);
  *camera = saved;
  return color.color * color.alpha;
}

/** @brief Largest channel gap between a drawn column and the chain shaded at
 *  one camera spin.
 * @param chain Effect after a draw.
 * @param x Column.
 * @param angle Camera spin, radians.
 * @return Largest RGB channel gap over the column. */
inline int chain_column_gap(ShaderChainWhiteBox::FX &chain, int x,
                            float angle) {
  int worst = 0;
  for (int y = 0; y < 20; ++y)
    worst = std::max<int>(
        worst, max_channel_gap(
                   chain.get_pixel(x, y),
                   chain_pixel_at_spin(chain, chain_scan_view(x, y), angle)));
  return worst;
}

/** Channel gap allowed between a drawn column and the spin evaluator; covers
 *  the generated palettes' per-frame fade. */
inline constexpr int CHAIN_TIME_TOLERANCE = 128;

/** @brief ShaderChain columns follow their sweep time through the camera spin
 *  wrap, matching an independent spin evaluator across the frame seam. */
inline void test_shader_chain_columns_follow_sweep_time() {
  reset_globals();
  ShaderChainWhiteBox::FX chain;
  chain.init();
  spin_default_chain(chain);
  int worst = 0;
  int held = 0;
  for (int draw = 0; draw < 130; ++draw) {
    pin_frame_clock(draw);
    chain.draw_frame();
    chain.advance_display();
    if (draw > 2 && draw < 123)
      continue;
    HS_CONTEXT("draw", draw);
    for (int x = 0; x < 96; ++x) {
      const float time =
          draw == 0 ? 1.0f : draw + Spacetime::SamplingSurface<96>::time(x);
      worst = std::max(worst, chain_column_gap(chain, x, CHAIN_SPIN * time));
      held =
          std::max(held, chain_column_gap(chain, x, CHAIN_SPIN * (draw + 1)));
    }
  }
  hs::clear_mock_time();
  std::printf("  ShaderChain sweep time: worst gap %d, held-endpoint gap %d\n",
              worst, held);
  HS_EXPECT_LE(worst, CHAIN_TIME_TOLERANCE);
  HS_EXPECT_GT(held, 16 * CHAIN_TIME_TOLERANCE);
}

/** @brief Segmented ShaderChain engines, including clips crossing both arm
 *  boundaries, draw the same pixels as one full-frame engine. */
inline void test_shader_chain_segmented_matches_single_engine() {
  using WB = ShaderChainWhiteBox;
  struct Segment {
    int x0, x1, margin;
  };
  constexpr std::array<Segment, 4> SEGMENTS{
      {{0, 40, 0}, {40, 60, 0}, {60, 96, 0}, {0, 8, 4}}};
  constexpr int DRAWS = 12;
  const auto render = [](const Segment &segment, auto &&inspect) {
    reset_globals();
    WB::FX engine;
    engine.init();
    spin_default_chain(engine);
    engine.set_clip(0, 20, segment.x0, segment.x1);
    engine.set_margin(segment.margin);
    for (int draw = 0; draw < DRAWS; ++draw) {
      pin_frame_clock(draw);
      engine.draw_frame();
      engine.advance_display();
      inspect(engine, draw);
    }
    hs::clear_mock_time();
  };
  std::vector<std::vector<Pixel>> single(DRAWS);
  render({0, 96, 0}, [&](const WB::FX &engine, int draw) {
    capture_frame<96, 20>(engine, single[draw]);
  });
  int mismatches = 0;
  int compared = 0;
  for (const Segment &segment : SEGMENTS)
    render(segment, [&](const WB::FX &engine, int draw) {
      for (int column = segment.x0 - segment.margin;
           column < segment.x1 + segment.margin; ++column) {
        const int x = (column + 96) % 96;
        for (int y = 0; y < 20; ++y) {
          ++compared;
          mismatches += engine.get_pixel(x, y) != single[draw][y * 96 + x];
        }
      }
    });
  HS_EXPECT_EQ(compared, DRAWS * (96 + 16) * 20);
  HS_EXPECT_EQ(mismatches, 0);
}

/** @brief A topology change, a snapshot restore and the first frame draw the
 *  current endpoint in every column; an unchanged frame does not. */
inline void test_shader_chain_discrete_changes_reset_history() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX chain;
  chain.init();
  spin_default_chain(chain);
  const auto current_gap = [&] {
    const In::FrameContext ctx = WB::frame_context(chain);
    int worst = 0;
    for (int x = 0; x < 96; ++x)
      for (int y = 0; y < 20; ++y) {
        const Color4 color =
            WB::program(chain).evaluate(chain_scan_view(x, y), ctx);
        worst =
            std::max<int>(worst, max_channel_gap(chain.get_pixel(x, y),
                                                 color.color * color.alpha));
      }
    return worst;
  };
  int frame = 0;
  const auto draw = [&] {
    pin_frame_clock(frame++);
    chain.draw_frame();
    chain.advance_display();
  };
  draw();
  HS_EXPECT_EQ(current_gap(), 0);
  draw();
  HS_EXPECT_GT(current_gap(), 0);
  WB::color_params(chain).mapping_mode =
      static_cast<uint8_t>(PB::Color::PaletteMapping::CUP);
  draw();
  HS_EXPECT_EQ(current_gap(), 0);
  draw();
  const ChainSnapshot snapshot = chain.snapshot();
  draw();
  draw();
  HS_EXPECT_GT(current_gap(), 0);
  HS_EXPECT_EQ(chain.restore_snapshot(snapshot),
               ChainSnapshotRestoreResult::APPLIED);
  draw();
  HS_EXPECT_EQ(current_gap(), 0);
  hs::clear_mock_time();
}

/** @brief Disabling temporal sampling draws the current endpoint while the
 *  history advances, so re-enabling samples from the last drawn frame. */
inline void test_shader_chain_fidelity_toggle_keeps_history() {
  reset_globals();
  ShaderChainWhiteBox::FX chain;
  chain.init();
  spin_default_chain(chain);
  int worst = 0;
  int draw = 0;
  for (const bool enabled : {true, true, false, false, true}) {
    HS_CONTEXT("draw", draw);
    Spacetime::set_flat_surface(!enabled);
    pin_frame_clock(draw);
    chain.draw_frame();
    chain.advance_display();
    for (int x = 0; x < 96; ++x) {
      const float time = draw == 0 || !enabled
                             ? draw + 1.0f
                             : draw + Spacetime::SamplingSurface<96>::time(x);
      worst = std::max(worst, chain_column_gap(chain, x, CHAIN_SPIN * time));
    }
    ++draw;
  }
  Spacetime::set_flat_surface(false);
  hs::clear_mock_time();
  HS_EXPECT_LE(worst, CHAIN_TIME_TOLERANCE);
}

/** @brief Columns blend the colours of the retained and current palettes at
 *  their sweep time, and the retained palettes are owned copies. */
inline void test_shader_chain_palettes_blend_across_the_frame() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX chain;
  chain.init();
  for (int draw = 0; draw < 2; ++draw) {
    pin_frame_clock(draw);
    chain.draw_frame();
    chain.advance_display();
  }
  const In::FrameContext before = WB::frame_context(chain);
  std::vector<Pixel> previous(96 * 20);
  for (int x = 0; x < 96; ++x)
    for (int y = 0; y < 20; ++y) {
      const Color4 color =
          WB::program(chain).evaluate(chain_scan_view(x, y), before);
      previous[y * 96 + x] = color.color * color.alpha;
    }
  const In::FrameContext &retained = WB::previous_context(chain);
  for (size_t index = 0; index < before.palettes.size(); ++index) {
    HS_EXPECT_TRUE(retained.palettes[index] != before.palettes[index]);
    HS_EXPECT_EQ(retained.palettes[index]->get(0.37f).color,
                 before.palettes[index]->get(0.37f).color);
  }
  WB::color_params(chain).palette_chroma = 0.05f;
  pin_frame_clock(2);
  chain.draw_frame();
  chain.advance_display();
  const In::FrameContext after = WB::frame_context(chain);
  int worst = 0;
  int changed = 0;
  for (int x = 0; x < 96; ++x) {
    const float fraction = Spacetime::SamplingSurface<96>::time(x);
    for (int y = 0; y < 20; ++y) {
      const Color4 color =
          WB::program(chain).evaluate(chain_scan_view(x, y), after);
      const Pixel current = color.color * color.alpha;
      const Pixel &old = previous[y * 96 + x];
      changed += max_channel_gap(old, current) > 1024;
      const auto blend = [&](uint16_t a, uint16_t b) {
        return static_cast<uint16_t>(a + (b - a) * fraction + 0.5f);
      };
      const Pixel expected(blend(old.r, current.r), blend(old.g, current.g),
                           blend(old.b, current.b));
      worst = std::max<int>(worst,
                            max_channel_gap(chain.get_pixel(x, y), expected));
    }
  }
  hs::clear_mock_time();
  std::printf("  ShaderChain palette blend: worst gap %d, changed %d\n", worst,
              changed);
  HS_EXPECT_GT(changed, 0);
  HS_EXPECT_LE(worst, CHAIN_TIME_TOLERANCE);
}

/** @brief Column params follow their transition curves between the retained
 *  and current param blocks. */
inline void test_shader_chain_column_params_interpolate() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX chain;
  chain.init();
  pin_frame_clock(0);
  chain.draw_frame();
  chain.advance_display();
  WB::color_params(chain).mapping_phase = 0.4f;
  pin_frame_clock(1);
  chain.draw_frame();
  chain.advance_display();
  In::ChainProgram &program = WB::program(chain);
  const In::FrameContext ctx = WB::frame_context(chain);
  HS_EXPECT_TRUE(program.history_continuous());
  program.capture_history();
  WB::color_params(chain).mapping_phase = 0.8f;
  program.prepare_column(ctx, ctx, 0.25f);
  const auto &prepared =
      *reinterpret_cast<const PB::Color::GeneratedPaletteState *>(
          program.prepared_block(3));
  HS_EXPECT_NEAR(prepared.mapping_offset, 0.5f, 1e-6f);
  program.prepare(ctx);
}

/** @brief Every shipped operator with instance state samples its clocks, each
 *  along the shortest arc of its period. */
inline void test_shader_chain_operator_states_sample_clocks() {
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    HS_CONTEXT(op.operator_id);
    HS_EXPECT_TRUE(op.runtime.sample_params != nullptr);
    HS_EXPECT_TRUE(op.runtime.state.size == sizeof(In::EmptyState) ||
                   op.runtime.sample_state != nullptr);
  }
  const float TWO_PI = math::TWO_PI_F;
  In::Op::SpatialWalkState walk_from, walk_to, walk;
  walk_from.spin_phase = TWO_PI - 0.1f;
  walk_to.spin_phase = 0.1f;
  walk_from.wander = math::make_rotation(math::X_AXIS, 0.2f);
  walk_to.wander = -math::make_rotation(math::X_AXIS, 0.4f);
  walk = walk_to;
  In::Op::temporal_sample(walk, In::Op::temporal_interval(walk_from, walk_to),
                          In::Op::RotateChainParams{}, 0.25f);
  HS_EXPECT_NEAR(std::remainder(walk.spin_phase, TWO_PI), -0.05f, 1e-5f);
  HS_EXPECT_NEAR(std::fabs(math::dot(walk.wander,
                                     math::make_rotation(math::X_AXIS, 0.25f))),
                 1.0f, 1e-6f);
  In::Op::SourceClockState source_from{TWO_PI - 0.2f, 0.2f, 3.0f};
  In::Op::SourceClockState source_to{0.2f, TWO_PI - 0.2f, 3.4f};
  In::Op::SourceClockState source = source_to;
  In::Op::temporal_sample(
      source, In::Op::temporal_interval(source_from, source_to), 0, 0.5f);
  HS_EXPECT_NEAR(std::remainder(source.primary, TWO_PI), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(std::remainder(source.secondary, TWO_PI), 0.0f, 1e-5f);
  HS_EXPECT_NEAR(source.angle, 3.2f, 1e-5f);
  In::Op::WarpPhaseState warp_from{0.9f}, warp_to{0.1f}, warp = warp_to;
  In::Op::temporal_sample(warp, In::Op::temporal_interval(warp_from, warp_to),
                          0, 0.5f);
  HS_EXPECT_NEAR(std::remainder(warp.phase, 1.0f), 0.0f, 1e-6f);
  In::Op::AffineClockState affine_from{0.95f, TWO_PI - 0.1f};
  In::Op::AffineClockState affine_to{0.05f, 0.1f};
  In::Op::AffineClockState affine = affine_to;
  In::Op::temporal_sample(
      affine, In::Op::temporal_interval(affine_from, affine_to), 0, 0.5f);
  HS_EXPECT_NEAR(std::remainder(affine.phase, 1.0f), 0.0f, 1e-6f);
  HS_EXPECT_NEAR(std::remainder(affine.rotation, TWO_PI), 0.0f, 1e-5f);
  In::Op::RipplePhaseState ripple_from{29.0f}, ripple_to{0.0f},
      ripple = ripple_to;
  PB::Surface::PeriodicRippleParams ripple_params;
  ripple_params.period = 30.0f;
  In::Op::temporal_sample(ripple,
                          In::Op::temporal_interval(ripple_from, ripple_to),
                          ripple_params, 0.5f);
  HS_EXPECT_NEAR(ripple.phase, 29.5f, 1e-5f);
  In::Op::ColorClockState color_from, color_to;
  color_from.oscillation_phase = 0.8f;
  color_to.oscillation_phase = 0.2f;
  color_from.hue_noise_phase = 0.3f;
  color_to.hue_noise_phase = 0.4f;
  In::Op::ColorClockState color = color_to;
  In::Op::temporal_sample(
      color, In::Op::temporal_interval(color_from, color_to), 0, 0.5f);
  HS_EXPECT_NEAR(std::remainder(color.oscillation_phase, 1.0f), 0.0f, 1e-6f);
  HS_EXPECT_EQ(color.hue_noise_phase, 0.4f);
}
