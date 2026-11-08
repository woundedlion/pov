/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- cases ----------------------------------------------------------------

inline void test_shader_chain_program_lifetime() {
  static_assert(!std::is_copy_constructible_v<In::ChainProgram>);
  static_assert(!std::is_copy_assignable_v<In::ChainProgram>);

  CountLifecycle::reset();
  auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  const In::ChainEntryRequest chain[] = {
      {"counter", "test.count-a.v2"},
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v4"},
      {"colorize", "colorize.generated-palette.v4"},
  };
  if (!expect_compiled(fixture->program.compile(chain)))
    return;
  HS_EXPECT_EQ(CountLifecycle::inits, 1);
  HS_EXPECT_EQ(CountLifecycle::destroys, 0);

  fixture.reset();
  HS_EXPECT_EQ(CountLifecycle::destroys, 1);
}

inline void test_shader_chain_table_behavior() {
  const auto ctx = shared_resources().context();
  for (const auto &op : In::OPERATOR_TABLE) {
    HS_CONTEXT(op.operator_id);
    size_t lit = 0;
    auto fixture = std::make_unique<ProgramFixture>();
    std::array<In::ChainEntryRequest, 5> chain{};
    size_t count = 0;
    const auto append = [&](const char *instance, const char *id) {
      chain[count++] = {instance, id};
    };
    if (op.input != In::CarrierId::SPHERE)
      append("project", "project.stereographic.v2");
    if (op.input == In::CarrierId::FIELD || op.input == In::CarrierId::COLOR)
      append("sample", "sample.grid.v4");
    if (op.input == In::CarrierId::COLOR)
      append("colorize", "colorize.generated-palette.v4");
    append("subject", op.operator_id);
    if (op.output == In::CarrierId::SPHERE)
      append("project", "project.stereographic.v2");
    if (op.output == In::CarrierId::SPHERE || op.output == In::CarrierId::PLANE)
      append("sample", "sample.grid.v4");
    if (op.output != In::CarrierId::COLOR)
      append("colorize", "colorize.generated-palette.v4");
    const auto status = fixture->program.compile(
        std::span<const In::ChainEntryRequest>(chain.data(), count));
    if (!expect_compiled(status))
      continue;
    for (int frame = 0; frame < 3; ++frame) {
      fixture->program.advance();
      fixture->program.prepare(ctx);
      for (const auto &view : sweep_views()) {
        const auto color = fixture->program.evaluate(view, ctx);
        lit += color.alpha > 0 && !is_black(color.color);
        const auto repeated = fixture->program.evaluate(view, ctx);
        HS_EXPECT_EQ(repeated.color, color.color);
        HS_EXPECT_EQ(repeated.alpha, color.alpha);
        HS_EXPECT_TRUE(std::isfinite(color.alpha));
        HS_EXPECT_GE(color.alpha, 0.0f);
        HS_EXPECT_LE(color.alpha, 1.0f);
      }
    }
    HS_EXPECT_GT(lit, size_t{0});
  }
}

/**
 * @brief Verifies each stateful operator's snapshot restores every captured
 *        field into a differently seeded instance.
 */
inline void test_shader_chain_state_restore_values() {
  constexpr int FRAMES = 37;
  const auto blocks = [](size_t bytes) {
    return std::vector<std::max_align_t>(
        (bytes + sizeof(std::max_align_t) - 1) / sizeof(std::max_align_t));
  };
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    if (op.runtime.state.size == sizeof(In::EmptyState))
      continue;
    HS_CONTEXT(op.operator_id);
    auto params = blocks(op.runtime.param.size);
    op.runtime.construct_params(params.data());
    for (uint16_t index = 0; index < op.schema_count; ++index) {
      const In::ParamFieldInfo &field = op.schema[index];
      if (field.enum_count == 0)
        *static_cast<float *>(op.runtime.param_address(params.data(), index)) =
            field.min + 0.37f * (field.max - field.min);
    }
    if (std::string_view(op.operator_id) == In::Op::ModulateColorV1::ID) {
      static_cast<In::Op::ModulationParams *>(
          static_cast<void *>(params.data()))
          ->modulation.noise.variant = Pattern::NoiseVariant::OCTAVE_CONTOUR;
      static_cast<In::Op::ModulationParams *>(
          static_cast<void *>(params.data()))
          ->modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
      static_cast<In::Op::ModulationParams *>(
          static_cast<void *>(params.data()))
          ->modulation.noise.scale = 4.0f;
    }
    const auto *param_bytes = reinterpret_cast<const uint8_t *>(params.data());
    auto advanced = blocks(op.runtime.state.size);
    auto fresh = blocks(op.runtime.state.size);
    auto restored = blocks(op.runtime.state.size);
    op.runtime.init(advanced.data(), {"a", op.operator_id, 1337});
    op.runtime.init(fresh.data(), {"a", op.operator_id, 1337});
    op.runtime.init(restored.data(), {"b", op.operator_id, 4242});
    for (int frame = 0; frame < FRAMES; ++frame)
      op.runtime.advance(advanced.data(), param_bytes);
    const auto a = op.runtime.capture_state(advanced.data());
    HS_EXPECT_TRUE(snapshot_words(a) !=
                   snapshot_words(op.runtime.capture_state(fresh.data())));
    HS_EXPECT_TRUE(op.runtime.restore_state(restored.data(), a));
    HS_EXPECT_TRUE(snapshot_words(op.runtime.capture_state(restored.data())) ==
                   snapshot_words(a));
    op.runtime.destroy(advanced.data());
    op.runtime.destroy(fresh.data());
    op.runtime.destroy(restored.data());
  }
}

inline void test_shader_chain_table_integrity() {
  static_assert(In::operator_ids_unique());
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    HS_EXPECT_TRUE(op.operator_id != nullptr && op.display_name != nullptr);
    HS_EXPECT_LE(static_cast<int>(op.input), static_cast<int>(op.output));
    HS_EXPECT_TRUE(op.runtime.construct_params != nullptr);
    HS_EXPECT_TRUE(op.runtime.init != nullptr);
    HS_EXPECT_TRUE(op.runtime.migrate != nullptr);
    HS_EXPECT_TRUE(op.runtime.destroy != nullptr);
    HS_EXPECT_TRUE(op.runtime.advance != nullptr);
    HS_EXPECT_TRUE(op.runtime.prepare != nullptr);
    HS_EXPECT_TRUE(op.runtime.run != nullptr);
    HS_EXPECT_TRUE(op.runtime.param_address != nullptr);
    HS_EXPECT_TRUE(op.runtime.capture_state != nullptr);
    HS_EXPECT_TRUE(op.runtime.restore_state != nullptr);
    HS_EXPECT_GT(op.runtime.param.size, 0u);
    HS_EXPECT_GT(op.runtime.state.size, 0u);
    std::vector<std::max_align_t> state(
        (op.runtime.state.size + sizeof(std::max_align_t) - 1) /
        sizeof(std::max_align_t));
    op.runtime.init(state.data(), {"snapshot", op.operator_id, 1337});
    const auto captured = op.runtime.capture_state(state.data());
    HS_EXPECT_EQ(std::holds_alternative<std::monostate>(captured),
                 op.runtime.state.size == sizeof(In::EmptyState));
    HS_EXPECT_TRUE(op.runtime.restore_state(state.data(), captured));
    op.runtime.destroy(state.data());
    for (uint16_t index = 0; index < op.schema_count; ++index) {
      const In::ParamFieldInfo &field = op.schema[index];
      HS_EXPECT_TRUE(field.id != nullptr);
      if (field.topology)
        HS_EXPECT_GE(field.enum_count, 2);
      else
        HS_EXPECT_LE(field.min, field.max);
    }
    const bool wants_oracle = op.approximate;
    HS_EXPECT_EQ(op.oracle != PB::ApproximationOracleId::NONE, wants_oracle);
    HS_EXPECT_EQ(op.metric_count > 0, wants_oracle);
  }
  for (const auto *id :
       {"project.peirce.v2", "project.peirce-square-fast.v2",
        "project.bonne.v2", "project.airocean.v2", "warp.affine.v2",
        "warp.affine.v3", "warp.affine-periodic.v1", "sample.grid.v2",
        "sample.twin-wave.v2", "colorize.generated-palette.v2"})
    HS_EXPECT_TRUE(In::find_operator(id) == nullptr);
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v4") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v1") == nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.displace.curl.v2") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.displace.ripple.v2") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.lens.kaleidoscope.v2") != nullptr);
  const auto *rotate_operator = require_operator("sphere.rotate.v2");
  if (!rotate_operator)
    return;
  HS_EXPECT_FALSE(rotate_operator->approximate);
  const auto *peirce_operator = require_operator("project.peirce.v3");
  if (!peirce_operator)
    return;
  HS_EXPECT_FALSE(peirce_operator->approximate);
  const auto *colorize_found =
      require_operator("colorize.generated-palette.v4");
  if (!colorize_found)
    return;
  const In::OperatorDescriptor &colorize = *colorize_found;
  HS_EXPECT_FALSE(colorize.approximate);
  const auto *modulate = require_operator("color.modulate.v1");
  if (!modulate)
    return;
  HS_EXPECT_TRUE(modulate->approximate);
  HS_EXPECT_EQ(
      static_cast<int>(modulate->oracle),
      static_cast<int>(PB::ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS));
  const auto *peirce_fast_found =
      require_operator("project.peirce-square-fast.v3");
  if (!peirce_fast_found)
    return;
  const In::OperatorDescriptor &peirce_fast = *peirce_fast_found;
  HS_EXPECT_TRUE(peirce_fast.approximate);
  HS_EXPECT_EQ(peirce_fast.oracle,
               PB::ApproximationOracleId::PEIRCE_FAST_SQUARE);
}

/** @brief Pins invariants operator schema fields and identities. */
inline void test_shader_chain_schema_invariants() {
  static_assert(std::is_same_v<In::Op::ProjectionCoverageMode,
                               PB::ProjectionCoverageMode>);
  static_assert(!In::topology_defaults_match<MismatchedTopologyDefaults>());
  static_assert(PB::field_ids_unique<In::Op::RotateChainParams>());
  static_assert(PB::field_ids_unique<In::Op::ProjectChainParams>());
  static_assert(PB::field_ids_unique<In::Op::GridSampleParams>());
  static_assert(PB::field_ids_unique<In::Op::GeneratedPaletteParams>());
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::Rotate>));
  static_assert(
      In::schema_ids_unique(In::SCHEMA<In::Op::ProjectStereographic>));
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::SampleGridV4>));
  static_assert(
      In::schema_ids_unique(In::SCHEMA<In::Op::ColorizeGeneratedPaletteV4>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::SampleGridV4>));
  static_assert(
      In::topology_wellformed(In::SCHEMA<In::Op::ColorizeGeneratedPaletteV4>));

  static_assert(PB::field_ids_unique<In::Op::CurlDisplaceParams>());
  static_assert(PB::field_ids_unique<In::Op::DirectDisplaceParams>());
  static_assert(PB::field_ids_unique<PB::Surface::PeriodicRippleParams>());
  static_assert(PB::field_ids_unique<In::Op::MobiusChainParams>());
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::DisplaceCurl>));
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::DisplaceRipple>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::DisplaceCurl>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::LensKaleidoscope>));
}

/** @brief Pins color operator schema fields and identities. */
inline void test_shader_chain_schema_color() {
  const auto *colorize_found =
      require_operator("colorize.generated-palette.v4");
  if (!colorize_found)
    return;
  const In::OperatorDescriptor &colorize = *colorize_found;
  constexpr size_t COLOR_FIELDS = In::Op::GeneratedPaletteParams::FIELDS.size();
  HS_EXPECT_SIZE_OR_RETURN(colorize.schema_span(), COLOR_FIELDS + 8);
  HS_EXPECT_STREQ(colorize.schema[COLOR_FIELDS].id, "palette-mode");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS].enum_count, 3);
  HS_EXPECT_STREQ(colorize.schema[COLOR_FIELDS + 1].id, "palette-mapping");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 1].enum_count, 4);
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 1].enum_def, 2);
  HS_EXPECT_STREQ(colorize.schema[COLOR_FIELDS + 2].id, "brightness-envelope");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 2].enum_count, 5);
  const auto *modulate = require_operator("color.modulate.v1");
  HS_EXPECT_TRUE(modulate != nullptr);
  if (modulate) {
    HS_EXPECT_EQ(modulate->input, In::CarrierId::COLOR);
    HS_EXPECT_EQ(modulate->output, In::CarrierId::COLOR);
    HS_EXPECT_EQ(modulate->schema_count,
                 In::Op::ModulationParams::EXTRA_SCHEMA.size());
    const auto *source = schema_field(*modulate, "modulation-source");
    HS_EXPECT_TRUE(source != nullptr);
    if (source) {
      HS_EXPECT_EQ(source->enum_count, 5);
      HS_EXPECT_STREQ(source->enum_ids[4], "displacement");
    }
  }
}

/** @brief Pins sphere operator schema fields and identities. */
inline void test_shader_chain_schema_sphere() {
  const auto *rotate_found = require_operator("sphere.rotate.v2");
  if (!rotate_found)
    return;
  const In::OperatorDescriptor &rotate = *rotate_found;
  HS_EXPECT_SIZE_OR_RETURN(rotate.schema_span(), 2);
  HS_EXPECT_STREQ(rotate.schema[0].id, "wander");
  HS_EXPECT_STREQ(rotate.schema[1].id, "spin-speed");

  const auto *glitch_operator = require_operator("sphere.lens.glitch.v2");
  if (!glitch_operator)
    return;
  HS_EXPECT_SIZE_OR_RETURN(glitch_operator->schema_span(), 2);
  HS_EXPECT_STREQ(glitch_operator->schema[0].id, "polar-multiplier");
  HS_EXPECT_STREQ(glitch_operator->schema[1].id, "azimuth-multiplier");
  const auto *twist_found = require_operator("sphere.lens.twist.v2");
  if (!twist_found)
    return;
  const In::OperatorDescriptor &twist = *twist_found;
  HS_EXPECT_SIZE_OR_RETURN(twist.schema_span(), 6);
  HS_EXPECT_STREQ(twist.schema[0].id, "twist-angle");
  HS_EXPECT_STREQ(twist.schema[1].id, "twist-distribution");
  HS_EXPECT_EQ(twist.schema[1].min, 0.0f);
  HS_EXPECT_EQ(twist.schema[1].max, 1.0f);
  HS_EXPECT_STREQ(twist.schema[2].id, "axis-phi");
  HS_EXPECT_EQ(twist.schema[2].def, 0.0f);
  HS_EXPECT_EQ(twist.schema[2].min, 0.0f);
  HS_EXPECT_EQ(twist.schema[2].max, math::PI_F);
  HS_EXPECT_STREQ(twist.schema[3].id, "axis-lambda");
  HS_EXPECT_EQ(twist.schema[3].def, 0.0f);
  HS_EXPECT_EQ(twist.schema[3].min, -math::PI_F);
  HS_EXPECT_EQ(twist.schema[3].max, math::PI_F);
  HS_EXPECT_STREQ(twist.schema[4].id, "speed");
  HS_EXPECT_EQ(twist.schema[4].min, 0.0f);
  HS_EXPECT_STREQ(twist.schema[5].id, "drift");
  HS_EXPECT_EQ(twist.schema[5].min, 0.0f);
  const auto *kaleidoscope_found =
      require_operator("sphere.lens.kaleidoscope.v2");
  if (!kaleidoscope_found)
    return;
  const In::OperatorDescriptor &kaleidoscope = *kaleidoscope_found;
  HS_EXPECT_SIZE_OR_RETURN(kaleidoscope.schema_span(), 1);
  HS_EXPECT_TRUE(kaleidoscope.schema[0].topology);
  HS_EXPECT_EQ(kaleidoscope.schema[0].enum_count, 9);
  HS_EXPECT_STREQ(kaleidoscope.schema[0].enum_ids[0], "azimuthal");
  HS_EXPECT_STREQ(kaleidoscope.schema[0].enum_ids[8], "octagonal-prism");
  const auto *mobius_found = require_operator("sphere.lens.mobius.v3");
  if (!mobius_found)
    return;
  const In::OperatorDescriptor &mobius = *mobius_found;
  HS_EXPECT_SIZE_OR_RETURN(mobius.schema_span(), 6);
  HS_EXPECT_STREQ(mobius.schema[0].id, "spin-x");
  HS_EXPECT_STREQ(mobius.schema[5].id, "push-z");
  HS_EXPECT_EQ(mobius.schema[0].min, -math::TWO_PI_F);
  HS_EXPECT_EQ(mobius.schema[0].max, math::TWO_PI_F);
  HS_EXPECT_EQ(mobius.schema[0].def, 0.0f);
  const auto *curl_found = require_operator("sphere.displace.curl.v2");
  if (!curl_found)
    return;
  const In::OperatorDescriptor &curl = *curl_found;
  constexpr size_t CURL_FIELDS = In::Op::CurlDisplaceParams::FIELDS.size();
  HS_EXPECT_SIZE_OR_RETURN(curl.schema_span(), CURL_FIELDS + 3);
  HS_EXPECT_STREQ(curl.schema[CURL_FIELDS].id, "basis");
  HS_EXPECT_EQ(curl.schema[CURL_FIELDS].enum_count, 3);
  HS_EXPECT_STREQ(curl.schema[CURL_FIELDS + 1].id, "integrator");
  HS_EXPECT_EQ(curl.schema[CURL_FIELDS + 1].enum_count, 3);
  HS_EXPECT_STREQ(curl.schema[CURL_FIELDS + 1].enum_ids[2], "midpoint-2x");
  HS_EXPECT_STREQ(curl.schema[CURL_FIELDS + 2].id, "limit");
  HS_EXPECT_EQ(curl.schema[CURL_FIELDS + 2].enum_count, 2);
  HS_EXPECT_STREQ(curl.schema[CURL_FIELDS + 2].enum_ids[0], "smooth");
  const auto *direct_found = require_operator("sphere.displace.direct.v2");
  if (!direct_found)
    return;
  const In::OperatorDescriptor &direct = *direct_found;
  constexpr size_t DIRECT_FIELDS = In::Op::DirectDisplaceParams::FIELDS.size();
  HS_EXPECT_SIZE_OR_RETURN(direct.schema_span(), DIRECT_FIELDS + 1);
  HS_EXPECT_STREQ(direct.schema[DIRECT_FIELDS].id, "basis");
  HS_EXPECT_EQ(direct.schema[DIRECT_FIELDS].enum_count, 3);
  const auto *ripple_found = require_operator("sphere.displace.ripple.v2");
  if (!ripple_found)
    return;
  const In::OperatorDescriptor &ripple = *ripple_found;
  HS_EXPECT_EQ(ripple.schema_count,
               In::Op::RippleChainParams::FIELDS.size() + 2);
  HS_EXPECT_STREQ(ripple.schema[0].id, "period");
  HS_EXPECT_STREQ(ripple.schema[5].id, "center-polar");
}

/** @brief Pins warp operator schema fields and identities. */
inline void test_shader_chain_schema_warp() {
  // Warp batch: every op is PLANE->PLANE with "speed" first; the polar chart
  // carries the harmonic list; curl-flow exposes basis and integrator.
  for (const char *id :
       {"warp.affine.v4", "warp.vortex.v2", "warp.wave-shear.v2",
        "warp.vector-noise.v2", "warp.mirror-tile.v2", "warp.polar-chart.v2",
        "warp.curl-flow.v2"}) {
    const auto *warp_found = require_operator(id);
    if (!warp_found)
      return;
    const In::OperatorDescriptor &warp = *warp_found;
    HS_EXPECT_EQ(warp.input, In::CarrierId::PLANE);
    HS_EXPECT_EQ(warp.output, In::CarrierId::PLANE);
    HS_EXPECT_STREQ(warp.schema[0].id, "speed");
  }
  const auto *polar_found = require_operator("warp.polar-chart.v2");
  if (!polar_found)
    return;
  const In::OperatorDescriptor &polar = *polar_found;
  constexpr size_t POLAR_FIELDS = In::Op::PolarChartParams::FIELDS.size();
  HS_EXPECT_SIZE_OR_RETURN(polar.schema_span(), POLAR_FIELDS + 2);
  HS_EXPECT_STREQ(polar.schema[POLAR_FIELDS].id, "mode");
  HS_EXPECT_EQ(polar.schema[POLAR_FIELDS].enum_count, 2);
  const In::ParamFieldInfo &harmonic = polar.schema[POLAR_FIELDS + 1];
  HS_EXPECT_STREQ(harmonic.id, "harmonic");
  HS_EXPECT_EQ(harmonic.enum_count, PB::Warp::MAX_POLAR_HARMONIC);
  HS_EXPECT_STREQ(harmonic.enum_ids[0], "h1");
  HS_EXPECT_STREQ(harmonic.enum_ids[15], "h16");
  const auto *shear_found = require_operator("warp.wave-shear.v2");
  if (!shear_found)
    return;
  const In::OperatorDescriptor &shear = *shear_found;
  const In::ParamFieldInfo &shear_envelope =
      shear.schema[shear.schema_count - 1];
  HS_EXPECT_STREQ(shear_envelope.id, "envelope");
  HS_EXPECT_EQ(shear_envelope.enum_count, 3);
  HS_EXPECT_STREQ(shear_envelope.enum_ids[1], "projection-weight");
  const auto *vector_noise_found = require_operator("warp.vector-noise.v2");
  if (!vector_noise_found)
    return;
  const In::OperatorDescriptor &vector_noise = *vector_noise_found;
  constexpr size_t VECTOR_NOISE_FIELDS =
      In::Op::VectorNoiseWarpParams::FIELDS.size();
  HS_EXPECT_SIZE_OR_RETURN(vector_noise.schema_span(), VECTOR_NOISE_FIELDS + 2);
  HS_EXPECT_STREQ(vector_noise.schema[VECTOR_NOISE_FIELDS].id, "basis");
  HS_EXPECT_EQ(vector_noise.schema[VECTOR_NOISE_FIELDS].enum_count, 3);
  HS_EXPECT_STREQ(vector_noise.schema[VECTOR_NOISE_FIELDS + 1].id, "envelope");
  HS_EXPECT_EQ(vector_noise.schema[VECTOR_NOISE_FIELDS + 1].enum_count, 3);
  const auto *curl_flow_found = require_operator("warp.curl-flow.v2");
  if (!curl_flow_found)
    return;
  const In::OperatorDescriptor &curl_flow = *curl_flow_found;
  constexpr size_t CURL_FLOW_FIELDS = In::Op::CurlFlowWarpParams::FIELDS.size();
  HS_EXPECT_SIZE_OR_RETURN(curl_flow.schema_span(), CURL_FLOW_FIELDS + 2);
  HS_EXPECT_STREQ(curl_flow.schema[CURL_FLOW_FIELDS].id, "basis");
  HS_EXPECT_EQ(curl_flow.schema[CURL_FLOW_FIELDS].enum_count, 3);
  const In::ParamFieldInfo &curl_integrator =
      curl_flow.schema[CURL_FLOW_FIELDS + 1];
  HS_EXPECT_STREQ(curl_integrator.id, "integrator");
  HS_EXPECT_EQ(curl_integrator.enum_count, 3);
  HS_EXPECT_STREQ(curl_integrator.enum_ids[2], "midpoint-4");
  const auto *vortex_found = require_operator("warp.vortex.v2");
  if (!vortex_found)
    return;
  const In::OperatorDescriptor &vortex = *vortex_found;
  HS_EXPECT_SIZE_OR_RETURN(vortex.schema_span(),
                           In::Op::VortexWarpParams::FIELDS.size());
  HS_EXPECT_STREQ(vortex.schema[1].id, "center-x");
}

/** @brief Pins sample operator schema fields and identities. */
inline void test_shader_chain_schema_sample() {
  // Schema order is the family table then the topology enum8s; defaults come
  // from the default-constructed family.
  const auto *sample_found = require_operator("sample.grid.v4");
  if (!sample_found)
    return;
  const In::OperatorDescriptor &sample = *sample_found;
  constexpr size_t GRID_FIELDS = In::Op::GridSampleParams::FIELDS.size();
  HS_EXPECT_SIZE_OR_RETURN(sample.schema_span(), GRID_FIELDS + 2);
  HS_EXPECT_STREQ(sample.schema[0].id, "pattern-freq");
  HS_EXPECT_STREQ(sample.schema[GRID_FIELDS - 1].id, "edge-width");
  const In::ParamFieldInfo &weight = sample.schema[GRID_FIELDS];
  HS_EXPECT_STREQ(weight.id, "weight-mode");
  HS_EXPECT_TRUE(weight.topology);
  HS_EXPECT_EQ(weight.enum_count, 2);
  HS_EXPECT_EQ(weight.enum_def, 1);
  HS_EXPECT_STREQ(weight.enum_ids[1], "projection");
  const In::ParamFieldInfo &coverage = sample.schema[GRID_FIELDS + 1];
  HS_EXPECT_EQ(coverage.enum_count, 5);
  HS_EXPECT_STREQ(coverage.enum_ids[1], "weight");
  HS_EXPECT_STREQ(coverage.enum_ids[2], "weight-squared");
  HS_EXPECT_STREQ(coverage.enum_ids[3], "edge-fade");
  HS_EXPECT_STREQ(coverage.enum_ids[4], "edge-fade-squared");

  // Projected sample batch: every crossing carries the union edge-width field
  // and the weight/coverage topology pair; only projected-noise adds a basis.
  for (const char *id :
       {"sample.grid.v4", "sample.twin-wave.v3", "sample.spiral.v3",
        "sample.bands.v1", "sample.tessellation.v3",
        "sample.projected-noise.v2"}) {
    const auto *crossing_found = require_operator(id);
    if (!crossing_found)
      return;
    const In::OperatorDescriptor &crossing = *crossing_found;
    HS_EXPECT_EQ(crossing.input, In::CarrierId::PLANE);
    HS_EXPECT_EQ(crossing.output, In::CarrierId::FIELD);
    bool has_edge_width = false;
    bool has_weight = false;
    bool has_coverage = false;
    bool has_basis = false;
    for (uint16_t field = 0; field < crossing.schema_count; ++field) {
      const std::string_view field_id{crossing.schema[field].id};
      has_edge_width |= field_id == "edge-width";
      has_weight |= field_id == "weight-mode";
      has_coverage |= field_id == "coverage-mode";
      has_basis |= field_id == "basis";
    }
    HS_EXPECT_TRUE(has_edge_width && has_weight && has_coverage);
    HS_EXPECT_EQ(has_basis,
                 std::string_view(id) == "sample.projected-noise.v2");
  }
  const In::OperatorDescriptor &tessellation =
      *In::find_operator("sample.tessellation.v3");
  const In::ParamFieldInfo *tessellation_kind =
      schema_field(tessellation, "kind");
  HS_EXPECT_TRUE(tessellation_kind != nullptr);
  if (tessellation_kind != nullptr)
    HS_EXPECT_EQ(tessellation_kind->enum_count, 3);
  const In::ParamFieldInfo *tessellation_fill =
      schema_field(tessellation, "fill");
  HS_EXPECT_TRUE(tessellation_fill != nullptr);
  if (tessellation_fill != nullptr)
    HS_EXPECT_EQ(tessellation_fill->enum_count, 2);
  for (const char *field : {"cell-shape", "cell-softness", "cell-radius",
                            "line-thickness", "line-softness"}) {
    const In::ParamFieldInfo *gated = schema_field(tessellation, field);
    HS_EXPECT_TRUE(gated != nullptr);
    if (gated)
      HS_EXPECT_STREQ(gated->gated_by, "fill");
  }
  const auto *cell_scale = schema_field(tessellation, "cell-scale");
  HS_EXPECT_TRUE(cell_scale != nullptr);
  if (cell_scale)
    HS_EXPECT_EQ(cell_scale->gated_by, nullptr);
  const In::ParamFieldInfo *projected_basis =
      schema_field(*In::find_operator("sample.projected-noise.v2"), "basis");
  HS_EXPECT_TRUE(projected_basis != nullptr);
  if (projected_basis != nullptr)
    HS_EXPECT_EQ(projected_basis->enum_count, 3);
}

/** @brief Pins projection operator schema fields and identities. */
inline void test_shader_chain_schema_projection() {
  // Projection batch: the meridian-consuming projections extend the shared
  // family with central-meridian; stereographic, gnomonic, and the fast square
  // Peirce do not.
  for (const char *id :
       {"project.stereographic.v2", "project.folded-sinusoidal.v2",
        "project.equirectangular.v2", "project.gnomonic.v2",
        "project.peirce.v3", "project.peirce-square-fast.v3",
        "project.bonne.v3", "project.airocean.v3", "project.tetrahedral.v1",
        "project.cubic.v1", "project.octahedral.v1", "project.dodecahedral.v1",
        "project.icosahedral.v1"}) {
    const auto *projection_found = require_operator(id);
    if (!projection_found)
      return;
    const In::OperatorDescriptor &projection = *projection_found;
    HS_EXPECT_EQ(projection.input, In::CarrierId::SPHERE);
    HS_EXPECT_EQ(projection.output, In::CarrierId::PLANE);
    bool has_meridian = false;
    const auto *frame_policy = schema_field(projection, "frame");
    HS_EXPECT_TRUE(frame_policy != nullptr);
    if (frame_policy != nullptr) {
      HS_EXPECT_TRUE(frame_policy->topology);
      HS_EXPECT_EQ(frame_policy->enum_count, 3);
      HS_EXPECT_EQ(frame_policy->enum_def, 1);
      HS_EXPECT_STREQ(frame_policy->enum_ids[0], "identity");
      HS_EXPECT_STREQ(frame_policy->enum_ids[1], "spin-wander");
      HS_EXPECT_STREQ(frame_policy->enum_ids[2], "spin-wander-tilted");
    }
    for (uint16_t field = 0; field < projection.schema_count; ++field)
      has_meridian |=
          std::string_view(projection.schema[field].id) == "central-meridian";
    const std::string_view id_view{id};
    HS_EXPECT_EQ(has_meridian, id_view != "project.stereographic.v2" &&
                                   id_view != "project.gnomonic.v2" &&
                                   id_view != "project.peirce-square-fast.v3");
  }
  const auto *gnomonic_found = require_operator("project.gnomonic.v2");
  if (!gnomonic_found)
    return;
  const In::OperatorDescriptor &gnomonic = *gnomonic_found;
  const In::ParamFieldInfo &gnomonic_hemisphere =
      gnomonic.schema[gnomonic.schema_count - 1];
  HS_EXPECT_STREQ(gnomonic_hemisphere.id, "hemisphere");
  HS_EXPECT_EQ(gnomonic_hemisphere.enum_count, 3);
  HS_EXPECT_STREQ(gnomonic_hemisphere.enum_ids[0], "folded");
  const auto *bonne_found = require_operator("project.bonne.v3");
  if (!bonne_found)
    return;
  const In::OperatorDescriptor &bonne = *bonne_found;
  const In::ParamFieldInfo &bonne_hemisphere =
      bonne.schema[bonne.schema_count - 1];
  HS_EXPECT_STREQ(bonne_hemisphere.id, "hemisphere");
  HS_EXPECT_EQ(bonne_hemisphere.enum_count, 2);
  HS_EXPECT_STREQ(bonne_hemisphere.enum_ids[0], "north");
}

/** @brief Pins field operator schema fields and identities. */
inline void test_shader_chain_schema_field() {
  // Field batch: ridge is schema-free; iso-contour, smooth-bands and the value
  // cutout carry their own families.
  const auto *ridge_operator = require_operator("field.transfer.ridge.v2");
  if (!ridge_operator)
    return;
  HS_EXPECT_EQ(ridge_operator->schema_count, 0);
  const auto *iso_found = require_operator("field.transfer.iso-contour.v2");
  if (!iso_found)
    return;
  const In::OperatorDescriptor &iso = *iso_found;
  HS_EXPECT_SIZE_OR_RETURN(iso.schema_span(), 2);
  HS_EXPECT_STREQ(iso.schema[0].id, "iso-level");
  HS_EXPECT_STREQ(iso.schema[1].id, "iso-width");
  const auto *bands_found = require_operator("field.transfer.smooth-bands.v2");
  if (!bands_found)
    return;
  const In::OperatorDescriptor &bands = *bands_found;
  HS_EXPECT_SIZE_OR_RETURN(bands.schema_span(), 2);
  HS_EXPECT_STREQ(bands.schema[0].id, "band-count");
  HS_EXPECT_STREQ(bands.schema[1].id, "band-phase");
  const auto *cutout_found = require_operator("field.coverage.value-cutout.v2");
  if (!cutout_found)
    return;
  const In::OperatorDescriptor &cutout = *cutout_found;
  HS_EXPECT_SIZE_OR_RETURN(cutout.schema_span(), 2);
  HS_EXPECT_STREQ(cutout.schema[0].id, "cutout-threshold");
  HS_EXPECT_STREQ(cutout.schema[1].id, "cutout-softness");
}

inline void test_shader_chain_instance_id_wellformed() {
  HS_EXPECT_TRUE(In::instance_id_wellformed("camera"));
  HS_EXPECT_TRUE(In::instance_id_wellformed("warp1"));
  HS_EXPECT_TRUE(In::instance_id_wellformed("a"));
  HS_EXPECT_TRUE(In::instance_id_wellformed("wave-shear-2"));
  HS_EXPECT_FALSE(In::instance_id_wellformed(""));
  HS_EXPECT_FALSE(In::instance_id_wellformed("Camera"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("1cam"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam.era"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam--x"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam-"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("-cam"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam era"));
  const std::string at_cap(In::MAX_INSTANCE_ID, 'a');
  HS_EXPECT_TRUE(In::instance_id_wellformed(at_cap));
  HS_EXPECT_FALSE(In::instance_id_wellformed(at_cap + "a"));
}

inline void test_shader_chain_slot_and_hash_contract() {
  static_assert(In::SLOT_SIZE == sizeof(PB::ColorSample));
  static_assert(In::SLOT_ALIGN == alignof(PB::ColorSample));
  static_assert(In::SLOT_SIZE >= sizeof(Color4));
  HS_EXPECT_LE(In::SLOT_SIZE, sizeof(void *) == 8 ? 96u : 76u);
  HS_EXPECT_EQ(In::SLOT_ALIGN, alignof(void *));
  HS_EXPECT_EQ(In::carrier_id_of<PB::SphereSample>(), In::CarrierId::SPHERE);
  HS_EXPECT_EQ(In::carrier_id_of<PB::ColorSample>(), In::CarrierId::COLOR);
  constexpr uint32_t HASH = In::instance_hash("camera", "sphere.rotate.v2");
  static_assert(HASH == 0x220a459au);
  HS_EXPECT_NE(HASH, In::instance_hash("camera2", "sphere.rotate.v2"));
  HS_EXPECT_NE(HASH, In::instance_hash("camera", "project.stereographic.v2"));
  // The separator keeps ("ab", "c") and ("a", "bc") distinct.
  HS_EXPECT_NE(In::instance_hash("ab", "c"), In::instance_hash("a", "bc"));
}

/** @brief Parses catalog JSON text; a parse failure fails the calling test.
 * @param text Catalog JSON.
 * @return The parsed root. */
inline JsonValue parse_catalog(const std::string &text) {
  JsonParser parser{text};
  JsonValue root = parser.parse_value();
  parser.skip_space();
  HS_EXPECT_FALSE(parser.failed);
  HS_EXPECT_EQ(parser.position, text.size());
  return root;
}

/** @brief The building ABI's catalog, exported in-process and parsed.
 * @return The parsed root. */
inline JsonValue native_catalog() {
  std::string text;
  In::append_catalog_json(text);
  return parse_catalog(text);
}

/** @brief A required member of a parsed catalog object.
 * @param object Object to read.
 * @param key Member name.
 * @return The member, or a null value after failing the calling test. */
inline const JsonValue &catalog_member(const JsonValue &object,
                                       const char *key) {
  static const JsonValue MISSING;
  const JsonValue *value = object.find(key);
  HS_EXPECT(value != nullptr, key);
  return value != nullptr ? *value : MISSING;
}

/** @brief A block member as an unsigned size.
 * @param block Parsed `{"size", "align"}` object.
 * @param key `"size"` or `"align"`.
 * @return The member's value. */
inline size_t catalog_block_field(const JsonValue &block, const char *key) {
  return static_cast<size_t>(catalog_member(block, key).number);
}

/** @brief Whether two parsed JSON values are equal, member order included.
 * @param a First value.
 * @param b Second value.
 * @return True when kind, scalar and every child agree. */
inline bool json_equal(const JsonValue &a, const JsonValue &b) {
  if (a.kind != b.kind || a.boolean != b.boolean || a.number != b.number ||
      a.text != b.text || a.items.size() != b.items.size() ||
      a.member_keys != b.member_keys)
    return false;
  for (size_t index = 0; index < a.items.size(); ++index)
    if (!json_equal(a.items[index], b.items[index]))
      return false;
  for (size_t index = 0; index < a.member_values.size(); ++index)
    if (!json_equal(a.member_values[index], b.member_values[index]))
      return false;
  return true;
}

/** @brief A parsed operator entry without its ABI-dependent `blocks` member.
 * @param op Parsed operator entry.
 * @return The entry's copy minus `blocks`. */
inline JsonValue without_blocks(const JsonValue &op) {
  JsonValue copy;
  copy.kind = op.kind;
  for (size_t index = 0; index < op.member_keys.size(); ++index) {
    if (op.member_keys[index] == "blocks")
      continue;
    copy.member_keys.push_back(op.member_keys[index]);
    copy.member_values.push_back(op.member_values[index]);
  }
  return copy;
}

/** @brief Cataloged runtime block names. */
inline constexpr const char *CATALOG_BLOCKS[] = {"param", "prepared", "state",
                                                 "interval"};

/** @brief Every cataloged operator id is unique and ends in `.v<N>`. */
inline void test_shader_chain_catalog_ids_versioned_and_unique() {
  const JsonValue root = native_catalog();
  const auto &operators = catalog_member(root, "operators").items;
  HS_EXPECT_EQ(operators.size(), In::OPERATOR_TABLE.size());
  std::vector<std::string> seen;
  for (const JsonValue &op : operators) {
    const std::string &id = catalog_member(op, "id").text;
    HS_CONTEXT(id.c_str());
    const size_t mark = id.rfind(".v");
    HS_EXPECT_TRUE(mark != std::string::npos && mark > 0);
    if (mark == std::string::npos || mark == 0)
      continue;
    const std::string_view version = std::string_view(id).substr(mark + 2);
    HS_EXPECT_FALSE(version.empty());
    HS_EXPECT_TRUE(version.empty() || version.front() != '0');
    HS_EXPECT_TRUE(std::all_of(version.begin(), version.end(),
                               [](char c) { return c >= '0' && c <= '9'; }));
    HS_EXPECT_TRUE(std::find(seen.begin(), seen.end(), id) == seen.end());
    seen.push_back(id);
  }
}

/** @brief Every numeric default lies in [min, max] and every topology default
 * is one of its values. */
inline void test_shader_chain_catalog_defaults_within_domain() {
  const JsonValue root = native_catalog();
  size_t numeric = 0;
  size_t topology = 0;
  for (const JsonValue &op : catalog_member(root, "operators").items) {
    HS_CONTEXT(catalog_member(op, "id").text.c_str());
    for (const JsonValue &param : catalog_member(op, "params").items) {
      HS_CONTEXT(catalog_member(param, "id").text.c_str());
      const JsonValue &def = catalog_member(param, "default");
      if (param.find("topology") != nullptr) {
        ++topology;
        const auto &values = catalog_member(param, "values").items;
        HS_EXPECT_FALSE(values.empty());
        HS_EXPECT_TRUE(def.kind == JsonValue::Kind::STRING);
        HS_EXPECT_TRUE(std::any_of(
            values.begin(), values.end(),
            [&def](const JsonValue &value) { return value.text == def.text; }));
      } else {
        ++numeric;
        HS_EXPECT_TRUE(def.kind == JsonValue::Kind::NUMBER);
        HS_EXPECT_LE(catalog_member(param, "min").number, def.number);
        HS_EXPECT_LE(def.number, catalog_member(param, "max").number);
      }
    }
  }
  HS_EXPECT_GT(numeric, 0u);
  HS_EXPECT_GT(topology, 0u);
}

/** @brief Whether a table entry carries its model's id and block layouts.
 * @tparam Model Operator model.
 * @param op Table entry.
 * @return True when id, param, prepared and state layouts match. */
template <typename Model>
consteval bool entry_matches_model(const In::OperatorDescriptor &op) {
  const auto matches = [](const In::BlockLayout &block, size_t size,
                          size_t align) {
    return block.size == size && block.align == align;
  };
  using Params = typename Model::Params;
  using Prepared = typename Model::Prepared;
  using State = typename Model::State;
  return std::string_view(op.operator_id) == Model::ID &&
         matches(op.runtime.param, sizeof(Params), alignof(Params)) &&
         matches(op.runtime.prepared, sizeof(Prepared), alignof(Prepared)) &&
         matches(op.runtime.state, sizeof(State), alignof(State));
}

/** @brief Whether every table entry matches the model at its list position.
 * @tparam Models Operator models.
 * @return True when every entry matches. */
template <typename... Models>
consteval bool table_matches_models(PB::Detail::TypeList<Models...>) {
  size_t index = 0;
  return (entry_matches_model<Models>(In::OPERATOR_TABLE[index++]) && ...);
}

static_assert(table_matches_models(In::OperatorModels{}));

/** @brief Cataloged block layouts equal the operator table's. */
inline void test_shader_chain_catalog_blocks_match_layouts() {
  const JsonValue root = native_catalog();
  const auto &operators = catalog_member(root, "operators").items;
  HS_EXPECT_SIZE_OR_RETURN(operators, In::OPERATOR_TABLE.size());
  for (size_t index = 0; index < operators.size(); ++index) {
    const In::OperatorDescriptor &expected = In::OPERATOR_TABLE[index];
    HS_CONTEXT(expected.operator_id);
    const JsonValue &blocks = catalog_member(operators[index], "blocks");
    const In::BlockLayout layouts[] = {
        expected.runtime.param, expected.runtime.prepared,
        expected.runtime.state, expected.runtime.interval};
    for (size_t kind = 0; kind < std::size(CATALOG_BLOCKS); ++kind) {
      HS_CONTEXT(CATALOG_BLOCKS[kind]);
      const JsonValue &block = catalog_member(blocks, CATALOG_BLOCKS[kind]);
      HS_EXPECT_EQ(catalog_block_field(block, "size"), layouts[kind].size);
      HS_EXPECT_EQ(catalog_block_field(block, "align"), layouts[kind].align);
    }
  }
}

/** @brief A full chain's fixed overhead and every operator's single-entry
 * footprint fit the cataloged arena. */
inline void test_shader_chain_catalog_fits_arena() {
  const JsonValue root = native_catalog();
  const JsonValue &budgets = catalog_member(root, "budgets");
  const auto budget = [&budgets](const char *key) {
    return static_cast<size_t>(catalog_member(budgets, key).number);
  };
  const size_t arena = budget("arena_bytes");
  const size_t per_op = budget("per_op_overhead_bytes");
  const size_t per_param = budget("per_param_name_bytes");
  HS_EXPECT_LE(budget("max_chain_ops") * per_op +
                   budget("max_params") * per_param,
               arena);
  for (const JsonValue &op : catalog_member(root, "operators").items) {
    HS_CONTEXT(catalog_member(op, "id").text.c_str());
    const JsonValue &blocks = catalog_member(op, "blocks");
    const JsonValue &interval = catalog_member(blocks, "interval");
    size_t cursor = 0;
    for (const char *kind : {"param", "prepared", "state"}) {
      const JsonValue &block = catalog_member(blocks, kind);
      size_t size = catalog_block_field(block, "size");
      size_t align = catalog_block_field(block, "align");
      if (std::string_view(kind) == "prepared") {
        size = std::max(size, catalog_block_field(interval, "size"));
        align = std::max(align, catalog_block_field(interval, "align"));
      }
      HS_EXPECT_GT(align, 0u);
      if (align == 0)
        break;
      cursor = (cursor + align - 1) / align * align + size;
    }
    cursor += per_op + catalog_member(op, "params").items.size() * per_param;
    HS_EXPECT_LE(cursor, arena);
  }
}

/** @brief The native catalog equals the wasm32 one outside block layouts,
 * where only pointer-bearing blocks widen. */
inline void test_shader_chain_catalog_matches_wasm_outside_blocks() {
  const std::string wasm_text = read_text_file(HS_ENGINE_CATALOG_PATH);
  HS_EXPECT_FALSE(wasm_text.empty());
  const JsonValue wasm = parse_catalog(wasm_text);
  const JsonValue native = native_catalog();
  for (const char *key : {"catalog_version", "budgets", "carriers"}) {
    HS_CONTEXT(key);
    HS_EXPECT_TRUE(
        json_equal(catalog_member(native, key), catalog_member(wasm, key)));
  }
  const auto &native_ops = catalog_member(native, "operators").items;
  const auto &wasm_ops = catalog_member(wasm, "operators").items;
  HS_EXPECT_SIZE_OR_RETURN(wasm_ops, native_ops.size());
  size_t widest = 0;
  for (size_t index = 0; index < native_ops.size(); ++index) {
    const JsonValue &host_op = native_ops[index];
    const JsonValue &wasm_op = wasm_ops[index];
    HS_CONTEXT(catalog_member(host_op, "id").text.c_str());
    HS_EXPECT_TRUE(
        json_equal(without_blocks(host_op), without_blocks(wasm_op)));
    size_t widening = 0;
    for (const char *kind : CATALOG_BLOCKS) {
      HS_CONTEXT(kind);
      const JsonValue &host_block =
          catalog_member(catalog_member(host_op, "blocks"), kind);
      const JsonValue &wasm_block =
          catalog_member(catalog_member(wasm_op, "blocks"), kind);
      const size_t host_size = catalog_block_field(host_block, "size");
      const size_t wasm_size = catalog_block_field(wasm_block, "size");
      const size_t host_align = catalog_block_field(host_block, "align");
      const size_t wasm_align = catalog_block_field(wasm_block, "align");
      if (std::string_view(kind) == "param")
        HS_EXPECT_TRUE(json_equal(host_block, wasm_block));
      HS_EXPECT_LE(wasm_size, host_size);
      if (host_size == wasm_size && host_align == wasm_align)
        continue;
      // wasm32 pointers are 4-byte aligned.
      HS_EXPECT_EQ(host_align, alignof(void *));
      HS_EXPECT_EQ(wasm_align, 4u);
      if (host_size >= wasm_size)
        widening += host_size - wasm_size;
    }
    widest = std::max(widest, widening);
  }
  const JsonValue &budgets = catalog_member(wasm, "budgets");
  HS_EXPECT_LT(
      widest *
          static_cast<size_t>(catalog_member(budgets, "max_chain_ops").number),
      static_cast<size_t>(catalog_member(budgets, "arena_bytes").number));
}

inline void test_shader_chain_catalog_shape() {
  std::string catalog;
  In::append_catalog_json(catalog);
  const auto contains = [&catalog](const char *needle) {
    return catalog.find(needle) != std::string::npos;
  };
  HS_EXPECT_TRUE(contains("\"catalog_version\":2"));
  const std::string budgets =
      "\"budgets\":{\"max_chain_ops\":" + std::to_string(In::MAX_CHAIN_OPS) +
      ",\"arena_bytes\":" + std::to_string(In::CHAIN_ARENA_BYTES) +
      ",\"max_params\":" + std::to_string(In::MAX_CHAIN_PARAMS) +
      ",\"max_instance_id_length\":" + std::to_string(In::MAX_INSTANCE_ID) +
      ",\"per_op_overhead_bytes\":" +
      std::to_string(In::PER_OP_OVERHEAD_BYTES) +
      ",\"per_param_name_bytes\":" + std::to_string(In::PER_PARAM_NAME_BYTES) +
      "}";
  HS_EXPECT_TRUE(contains(budgets.c_str()));
  HS_EXPECT_TRUE(
      contains("\"carriers\":[\"sphere\",\"plane\",\"field\",\"color\"]"));
  HS_EXPECT_TRUE(contains("\"id\":\"sphere.rotate.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"project.stereographic.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"warp.vortex.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.grid.v4\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.tessellation.v3\""));
  HS_EXPECT_TRUE(contains("\"id\":\"colorize.generated-palette.v4\""));
  HS_EXPECT_TRUE(contains("\"id\":\"weight-mode\",\"topology\":true,"
                          "\"values\":[\"none\",\"projection\"],"
                          "\"default\":\"projection\""));
  HS_EXPECT_TRUE(contains("\"values\":[\"none\",\"weight\",\"weight-squared\","
                          "\"edge-fade\",\"edge-fade-squared\"]"));
  HS_EXPECT_TRUE(contains("\"approximate\":true"));
  // The cataloged block layouts are the real ABI figures.
  std::string sample_blocks = "\"blocks\":{\"param\":{\"size\":";
  sample_blocks += std::to_string(sizeof(In::Op::GridSampleParams));
  HS_EXPECT_TRUE(contains(sample_blocks.c_str()));
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    HS_CONTEXT(op.operator_id);
    HS_EXPECT_TRUE(op.display_name != nullptr && *op.display_name != '\0');
    if (op.display_name == nullptr)
      continue;
    const std::string expected = std::string("\"id\":\"") + op.operator_id +
                                 "\",\"name\":\"" + op.display_name + "\"";
    HS_EXPECT_TRUE(contains(expected.c_str()));
  }
}

inline void test_shader_chain_default_chain_renders() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  HS_EXPECT_FALSE(program.compiled());
  if (!arm_default_chain(program, 3, ValueSet::DEFAULTS))
    return;
  HS_EXPECT_TRUE(program.compiled());
  HS_EXPECT_SIZE_OR_RETURN(program.ops(), 4);
  HS_EXPECT_STREQ(program.ops()[0].instance, "camera");
  HS_EXPECT_STREQ(program.ops()[3].instance, "colorize");
  HS_EXPECT_EQ(program.ops()[0].stable_hash,
               In::instance_hash("camera", "sphere.rotate.v2"));
  HS_EXPECT_GT(program.used_bytes(), 0u);
  HS_EXPECT_LE(program.used_bytes(), In::CHAIN_ARENA_BYTES);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  bool any_opaque = false;
  for (const math::Vector &view : sweep_views()) {
    const Color4 color = program.evaluate(view, ctx);
    HS_EXPECT_TRUE(std::isfinite(color.alpha));
    HS_EXPECT_GE(color.alpha, 0.0f);
    if (color.alpha > 0.0f &&
        (color.color.r | color.color.g | color.color.b) != 0)
      any_opaque = true;
  }
  HS_EXPECT_TRUE(any_opaque);
  program.clear();
  HS_EXPECT_FALSE(program.compiled());
}

inline void test_shader_chain_param_address_channel() {
  const auto verify_crossing_address = []<typename Model>() {
    typename Model::Params params;
    const auto before = params;
    const auto *descriptor_found = require_operator(Model::ID);
    if (!descriptor_found)
      return;
    const auto &descriptor = *descriptor_found;
    constexpr size_t FIELDS = Model::Params::FIELDS.size();
    auto *coverage = static_cast<uint8_t *>(
        descriptor.runtime.param_address(&params, FIELDS + 1));
    *coverage = static_cast<uint8_t>(In::Op::ProjectionCoverageMode::NONE);
    HS_EXPECT_EQ(params.coverage_mode, 0);
    for (const auto &field : Model::Params::FIELDS)
      HS_EXPECT_EQ(params.*field.member, before.*field.member);
  };
  verify_crossing_address.template operator()<In::Op::SampleGridV4>();
  verify_crossing_address.template operator()<In::Op::SampleTwinWaveV3>();
  verify_crossing_address.template operator()<In::Op::SampleSpiralV3>();
  verify_crossing_address.template operator()<In::Op::SampleBands>();
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  if (!arm_default_chain(program, 1, ValueSet::DEFAULTS))
    return;
  const In::OperatorDescriptor &sample_op = *program.ops()[2].op;
  uint8_t *block = program.param_block(2);
  auto &typed = *reinterpret_cast<In::Op::GridSampleParams *>(block);
  // Every schema index resolves inside the block.
  for (uint16_t index = 0; index < sample_op.schema_count; ++index) {
    const uint8_t *address = static_cast<const uint8_t *>(
        sample_op.runtime.param_address(block, index));
    HS_EXPECT_GE(address, block);
    HS_EXPECT_LT(address, block + sample_op.runtime.param.size);
  }
  auto *freq = static_cast<float *>(sample_op.runtime.param_address(block, 0));
  *freq = 5.5f;
  HS_EXPECT_EQ(typed.pattern_freq, 5.5f);
  constexpr size_t GRID_FIELDS = In::Op::GridSampleParams::FIELDS.size();
  auto *coverage_mode = static_cast<uint8_t *>(
      sample_op.runtime.param_address(block, GRID_FIELDS + 1));
  *coverage_mode =
      static_cast<uint8_t>(In::Op::ProjectionCoverageMode::EDGE_FADE);
  HS_EXPECT_EQ(typed.coverage_mode,
               static_cast<uint8_t>(In::Op::ProjectionCoverageMode::EDGE_FADE));
  // The write is visible to the render: identical view, different coverage.
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  const math::Vector view = math::Y_AXIS;
  const Color4 faded = program.evaluate(view, ctx);
  *coverage_mode = static_cast<uint8_t>(In::Op::ProjectionCoverageMode::NONE);
  program.prepare(ctx);
  const Color4 full = program.evaluate(view, ctx);
  HS_EXPECT_LT(faded.alpha, full.alpha);
  HS_EXPECT_LT(faded.alpha, 0.5f);
  program.clear();
}

inline void test_shader_chain_parity_rotate_project() {
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_default_chain(program, 5, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    program.prepare(ctx);
    // The erased prepared conjugates equal the same composition recomputed
    // from the state blocks.
    const MirrorFrame mirror = mirror_from(program, ctx);
    const auto &camera_prepared =
        prepared_as<In::Op::Rotate::Prepared>(program, 0);
    const auto &project_prepared =
        prepared_as<In::Op::ProjectStereographic::Prepared>(program, 1);
    HS_EXPECT_EQ(std::memcmp(&camera_prepared.conjugate,
                             &mirror.camera_conjugate,
                             sizeof(math::Quaternion)),
                 0);
    HS_EXPECT_EQ(std::memcmp(&project_prepared.conjugate,
                             &mirror.projection_conjugate,
                             sizeof(math::Quaternion)),
                 0);
    // Op-level kernel parity of the erased camera and projection.
    using BoundRotate = PB::Stage::Rotate<MirrorCamera>::Bind<MirrorBinding>;
    using BoundProject = PB::Stage::Project<
        PB::Projection::Stereographic<MirrorProjection>>::Bind<MirrorBinding>;
    expect_op_parity<BoundRotate>(
        program, ctx, 0, mirror,
        [](const math::Vector &view) { return PB::SphereSample{view, 0.0f}; },
        sphere_identical);
    expect_op_parity<BoundProject>(
        program, ctx, 1, mirror,
        [&](const math::Vector &view) {
          return run_entry<PB::SphereSample>(program, ctx, 0,
                                             PB::SphereSample{view, 0.0f});
        },
        plane_identical);
  });
}

/** @brief Erased parameter defaults match every registered schema field. */
inline void test_shader_chain_schema_defaults_match_params() {
  for (const auto &op : In::OPERATOR_TABLE) {
    HS_CONTEXT(op.operator_id);
    std::vector<std::max_align_t> storage(
        (op.runtime.param.size + sizeof(std::max_align_t) - 1) /
        sizeof(std::max_align_t));
    op.runtime.construct_params(storage.data());
    const auto BEGIN = reinterpret_cast<uintptr_t>(storage.data());
    const auto END = BEGIN + op.runtime.param.size;
    for (uint16_t index = 0; index < op.schema_count; ++index) {
      const auto &field = op.schema[index];
      HS_CONTEXT(field.id);
      const auto *address = op.runtime.param_address(storage.data(), index);
      const auto VALUE = reinterpret_cast<uintptr_t>(address);
      const size_t BYTES = field.enum_count ? sizeof(uint8_t) : sizeof(float);
      HS_EXPECT_GE(VALUE, BEGIN);
      HS_EXPECT_LE(VALUE, END);
      HS_EXPECT_LE(BYTES, END - VALUE);
      if (VALUE < BEGIN || VALUE > END || BYTES > END - VALUE)
        continue;
      if (field.enum_count)
        HS_EXPECT_EQ(*static_cast<const uint8_t *>(address), field.enum_def);
      else
        HS_EXPECT_TRUE(
            float_identical(*static_cast<const float *>(address), field.def));
    }
  }
}
