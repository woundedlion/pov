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
  HS_EXPECT_EQ(static_cast<int>(chain.set_chain(topology).code),
               static_cast<int>(In::ChainStatus::OK));
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
      HS_EXPECT_NEAR(actual.color.r, expected.color.r, 1);
      HS_EXPECT_NEAR(actual.color.g, expected.color.g, 1);
      HS_EXPECT_NEAR(actual.color.b, expected.color.b, 1);
      HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
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
  for (const In::ChainEntryRequest &entry : DEFAULT_CHAIN)
    expected += In::find_operator(entry.operator_id)->schema_count;
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

  uint64_t lit = 0;
  for (int frame = 0; frame < 4; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
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
  HS_EXPECT_EQ(static_cast<int>(effect.set_chain(chain).code),
               static_cast<int>(In::ChainStatus::OK));
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
  const float saved_a =
      effect.getParameters().find("lens.mobius-a-re")->get_requested();
  const float saved_d = d_re.get_requested();
  HS_EXPECT_EQ(static_cast<int>(effect.set_chain(chain).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-a-re", saved_a),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("lens.mobius-d-re", saved_d),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.getParameters().find("lens.mobius-a-re")->get_requested(),
               saved_a);
  HS_EXPECT_EQ(effect.getParameters().find("lens.mobius-d-re")->get_requested(),
               saved_d);
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
    HS_EXPECT_EQ(effect.getParameters().find(write.name)->get_requested(),
                 write.value);
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
    HS_EXPECT_EQ(effect.getParameters().find(write.name)->get_requested(),
                 write.value);
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
  HS_EXPECT_EQ(effect.getParameters().find("warp.envelope")->get_requested(),
               ACCEPTED_ENVELOPE);
  HS_EXPECT_TRUE(effect.parameter_warning("warp.envelope") != nullptr);
  HS_EXPECT_EQ(
      effect.getParameters().find("sample.coverage-mode")->get_requested(),
      1.0f);
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
    HS_EXPECT_EQ(program.compile(chain).code, In::ChainStatus::OK);
    param_as<In::Op::ProjectChainParams>(program, 0).frame =
        static_cast<uint8_t>(In::Op::ProjectionFrame::IDENTITY);
    const auto ctx = shared_resources().context();
    program.prepare(ctx);
    HS_EXPECT_TRUE(op.runtime.plane_bound != nullptr);
    if (op.runtime.plane_bound == nullptr)
      continue;
    const float bound = op.runtime.plane_bound(program.param_block(0), 0.0f);
    for (const math::Vector &view : views) {
      const PB::SphereSample sphere{view, 0.0f};
      PB::PlaneSample plane{};
      op.runtime.run(&sphere, &plane, ctx, program.param_block(0),
                     program.prepared_block(0));
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
  HS_EXPECT_EQ(effect.getParameters().find(scale_x[0].c_str())->get_requested(),
               1.0f);

  int admitted = 0;
  for (; admitted < WARPS; ++admitted) {
    const ShaderChainParameterWrite step[] = {
        {scale_x[admitted].c_str(), MIN_SCALE},
        {scale_y[admitted].c_str(), MIN_SCALE}};
    if (effect.update_parameters(step) != ParamSetResult::APPLIED)
      break;
  }
  HS_EXPECT_EQ(admitted, 5);
  HS_EXPECT_TRUE(effect.parameter_warning(scale_x[5].c_str()) ==
                 In::PLANE_GROWTH_WARNING);
  HS_EXPECT_EQ(effect.updateParameter(scale_y[5].c_str(), MIN_SCALE),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(effect.parameter_warning(scale_y[5].c_str()) ==
                 In::PLANE_GROWTH_WARNING);
  HS_EXPECT_EQ(effect.getParameters().find(scale_y[5].c_str())->get_requested(),
               1.0f);
  HS_EXPECT_EQ(effect.updateParameter(scale_x[5].c_str(), 1.5f),
               ParamSetResult::APPLIED);

  const ChainSnapshot live = effect.snapshot();
  ChainSnapshot overflowing = live;
  for (auto &parameter : overflowing.parameters)
    if (parameter.name == scale_y[5])
      parameter.value = MIN_SCALE;
  HS_EXPECT_EQ(effect.restore_snapshot(overflowing),
               ChainSnapshotRestoreResult::INVALID_VALUE);
  HS_EXPECT_EQ(effect.getParameters().find(scale_y[5].c_str())->get_requested(),
               1.0f);
  HS_EXPECT_EQ(effect.restore_snapshot(live),
               ChainSnapshotRestoreResult::APPLIED);
  effect.draw_frame();
  effect.advance_display();
  auto &program = WB::program(effect);
  const In::FrameContext ctx = WB::frame_context(effect);
  for (const math::Vector &view : sweep_views()) {
    const PB::SphereSample sphere{view, 0.0f};
    PB::PlaneSample plane{};
    for (int index = 0; index <= WARPS; ++index) {
      PB::PlaneSample next{};
      program.ops()[index].op->runtime.run(
          index == 0 ? static_cast<const void *>(&sphere)
                     : static_cast<const void *>(&plane),
          &next, ctx, program.param_block(index),
          program.prepared_block(index));
      plane = next;
    }
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
  // Transactional at the effect layer too: definitions, generation, and the
  // committed program all survive, and the effect still renders.
  HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), committed);
  HS_EXPECT_EQ(effect.getParameters().size(), param_count);
  HS_EXPECT_TRUE(effect.getParameters().find("camera.wander") != nullptr);
  uint64_t lit = 0;
  for (int frame = 0; frame < 2; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
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

  for (int frame = 0; frame < 60; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }

  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_NE(
      state_as<In::Op::SourceClockState>(WB::program(effect), 2).primary,
      source_before.primary);
  HS_EXPECT_TRUE(
      state_as<In::Op::SpatialWalkState>(WB::program(effect), 0).wander !=
      wander_before);
  const Pixel color_after = WB::palette_color(effect, 0.25f);
  HS_EXPECT_TRUE(color_after.r != color_before.r ||
                 color_after.g != color_before.g ||
                 color_after.b != color_before.b);
  uint64_t lit = 0;
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
}
