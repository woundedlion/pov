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
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  HS_EXPECT_EQ(static_cast<int>(fixture->program.compile(chain).code),
               static_cast<int>(In::ChainStatus::OK));
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
    if (op.input == In::CarrierId::FIELD)
      append("sample", "sample.grid.v3");
    append("subject", op.operator_id);
    if (op.output == In::CarrierId::SPHERE)
      append("project", "project.stereographic.v2");
    if (op.output == In::CarrierId::SPHERE || op.output == In::CarrierId::PLANE)
      append("sample", "sample.grid.v3");
    if (op.output != In::CarrierId::COLOR)
      append("colorize", "colorize.generated-palette.v3");
    const auto status = fixture->program.compile(
        std::span<const In::ChainEntryRequest>(chain.data(), count));
    HS_EXPECT_EQ(static_cast<int>(status.code),
                 static_cast<int>(In::ChainStatus::OK));
    if (status.code != In::ChainStatus::OK)
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

inline void test_shader_chain_table_integrity() {
  static_assert(In::operator_ids_unique());
  static_assert(In::operator_table_monotone());
  HS_EXPECT_EQ(In::OPERATOR_TABLE.size(), 38u);
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
        "sample.grid.v2", "sample.twin-wave.v2",
        "colorize.generated-palette.v2"})
    HS_EXPECT_TRUE(In::find_operator(id) == nullptr);
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v3") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v1") == nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.displace.curl.v2") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.displace.ripple.v2") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.lens.kaleidoscope.v2") != nullptr);
  HS_EXPECT_FALSE(In::find_operator("sphere.rotate.v2")->approximate);
  HS_EXPECT_FALSE(In::find_operator("project.peirce.v3")->approximate);
  const In::OperatorDescriptor &colorize =
      *In::find_operator("colorize.generated-palette.v3");
  HS_EXPECT_TRUE(colorize.approximate);
  HS_EXPECT_EQ(
      static_cast<int>(colorize.oracle),
      static_cast<int>(PB::ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS));
  const In::OperatorDescriptor &peirce_fast =
      *In::find_operator("project.peirce-square-fast.v3");
  HS_EXPECT_TRUE(peirce_fast.approximate);
  HS_EXPECT_EQ(static_cast<int>(peirce_fast.oracle),
               static_cast<int>(PB::ApproximationOracleId::PEIRCE_FAST_SQUARE));
  HS_EXPECT_EQ(peirce_fast.metric_count, 3);
  size_t approximate_count = 0;
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE)
    approximate_count += op.approximate ? 1 : 0;
  HS_EXPECT_EQ(approximate_count, 2u);
}

/** Schema field of @p op with id @p field_id, or null. */
inline const In::ParamFieldInfo *schema_field(const In::OperatorDescriptor &op,
                                              std::string_view field_id) {
  for (const In::ParamFieldInfo &field : op.schema_span())
    if (std::string_view(field.id) == field_id)
      return &field;
  return nullptr;
}

inline void test_shader_chain_schema_and_field_ids() {
  static_assert(std::is_same_v<In::Op::HueShiftMode, PB::Color::HueMode>);
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
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::SampleGridV3>));
  static_assert(
      In::schema_ids_unique(In::SCHEMA<In::Op::ColorizeGeneratedPaletteV3>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::SampleGridV3>));
  static_assert(
      In::topology_wellformed(In::SCHEMA<In::Op::ColorizeGeneratedPaletteV3>));

  static_assert(PB::field_ids_unique<In::Op::CurlDisplaceParams>());
  static_assert(PB::field_ids_unique<In::Op::DirectDisplaceParams>());
  static_assert(PB::field_ids_unique<PB::Surface::PeriodicRippleParams>());
  static_assert(PB::field_ids_unique<In::Op::MobiusChainParams>());
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::DisplaceCurl>));
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::DisplaceRipple>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::DisplaceCurl>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::LensKaleidoscope>));

  // Schema order is the family table then the topology enum8s; defaults come
  // from the default-constructed family.
  const In::OperatorDescriptor &sample = *In::find_operator("sample.grid.v3");
  constexpr size_t GRID_FIELDS = In::Op::GridSampleParams::FIELDS.size();
  HS_EXPECT_EQ(sample.schema_count, GRID_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(sample.schema[0].id) == "pattern-freq");
  HS_EXPECT_EQ(sample.schema[0].min, 0.01f);
  HS_EXPECT_EQ(sample.schema[0].max, 64.0f);
  HS_EXPECT_TRUE(std::string_view(sample.schema[GRID_FIELDS - 1].id) ==
                 "edge-width");
  const In::ParamFieldInfo &weight = sample.schema[GRID_FIELDS];
  HS_EXPECT_TRUE(std::string_view(weight.id) == "weight-mode");
  HS_EXPECT_TRUE(weight.topology);
  HS_EXPECT_EQ(weight.enum_count, 2);
  HS_EXPECT_EQ(weight.enum_def, 1);
  HS_EXPECT_TRUE(std::string_view(weight.enum_ids[1]) == "projection");
  const In::ParamFieldInfo &coverage = sample.schema[GRID_FIELDS + 1];
  HS_EXPECT_EQ(coverage.enum_count, 4);
  HS_EXPECT_TRUE(std::string_view(coverage.enum_ids[1]) == "weight");
  HS_EXPECT_TRUE(std::string_view(coverage.enum_ids[2]) == "weight-squared");
  HS_EXPECT_TRUE(std::string_view(coverage.enum_ids[3]) == "edge-fade");

  const In::OperatorDescriptor &colorize =
      *In::find_operator("colorize.generated-palette.v3");
  constexpr size_t COLOR_FIELDS = PB::Color::ColorParams::FIELDS.size();
  HS_EXPECT_EQ(colorize.schema_count, COLOR_FIELDS + 4);
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS].id) ==
                 "palette-mode");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS].enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS + 1].id) ==
                 "palette-mapping");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 1].enum_count, 4);
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 1].enum_def, 2);
  const In::ParamFieldInfo &hue = colorize.schema[COLOR_FIELDS + 2];
  HS_EXPECT_EQ(hue.enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(hue.enum_ids[0]) == "none");
  HS_EXPECT_TRUE(std::string_view(hue.enum_ids[1]) == "noise");
  HS_EXPECT_TRUE(std::string_view(hue.enum_ids[2]) == "path-length");
  HS_EXPECT_EQ(hue.enum_def, static_cast<uint8_t>(PB::Color::HueMode::NOISE));
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS + 3].id) ==
                 "brightness-envelope");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 3].enum_count, 5);

  const In::OperatorDescriptor &rotate = *In::find_operator("sphere.rotate.v2");
  HS_EXPECT_EQ(rotate.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(rotate.schema[0].id) == "wander");
  HS_EXPECT_TRUE(std::string_view(rotate.schema[1].id) == "spin-speed");
  HS_EXPECT_EQ(rotate.schema[1].max, 0.05f);

  HS_EXPECT_EQ(In::find_operator("sphere.lens.glitch.v2")->schema_count, 0);
  const In::OperatorDescriptor &twist =
      *In::find_operator("sphere.lens.twist.v2");
  HS_EXPECT_EQ(twist.schema_count, 1);
  HS_EXPECT_TRUE(std::string_view(twist.schema[0].id) == "twist-rate");
  HS_EXPECT_EQ(twist.schema[0].def, lenses::TWIST_RATE);
  const In::OperatorDescriptor &kaleidoscope =
      *In::find_operator("sphere.lens.kaleidoscope.v2");
  HS_EXPECT_EQ(kaleidoscope.schema_count, 1);
  HS_EXPECT_TRUE(kaleidoscope.schema[0].topology);
  HS_EXPECT_EQ(kaleidoscope.schema[0].enum_count, 9);
  HS_EXPECT_TRUE(std::string_view(kaleidoscope.schema[0].enum_ids[0]) ==
                 "azimuthal");
  HS_EXPECT_TRUE(std::string_view(kaleidoscope.schema[0].enum_ids[8]) ==
                 "octagonal-prism");
  const In::OperatorDescriptor &mobius =
      *In::find_operator("sphere.lens.mobius.v2");
  HS_EXPECT_EQ(mobius.schema_count, 8);
  HS_EXPECT_TRUE(std::string_view(mobius.schema[0].id) == "mobius-a-re");
  HS_EXPECT_TRUE(std::string_view(mobius.schema[7].id) == "mobius-d-im");
  HS_EXPECT_EQ(mobius.schema[0].min, -4.0f);
  HS_EXPECT_EQ(mobius.schema[0].max, 4.0f);
  HS_EXPECT_EQ(mobius.schema[0].def, 0.7071067811865475f);
  const In::OperatorDescriptor &curl =
      *In::find_operator("sphere.displace.curl.v2");
  constexpr size_t CURL_FIELDS = In::Op::CurlDisplaceParams::FIELDS.size();
  HS_EXPECT_EQ(curl.schema_count, CURL_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(curl.schema[CURL_FIELDS].id) == "basis");
  HS_EXPECT_EQ(curl.schema[CURL_FIELDS].enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(curl.schema[CURL_FIELDS + 1].id) ==
                 "integrator");
  HS_EXPECT_EQ(curl.schema[CURL_FIELDS + 1].enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(curl.schema[CURL_FIELDS + 1].enum_ids[2]) ==
                 "midpoint-2x");
  const In::OperatorDescriptor &direct =
      *In::find_operator("sphere.displace.direct.v2");
  constexpr size_t DIRECT_FIELDS = In::Op::DirectDisplaceParams::FIELDS.size();
  HS_EXPECT_EQ(direct.schema_count, DIRECT_FIELDS + 1);
  HS_EXPECT_TRUE(std::string_view(direct.schema[DIRECT_FIELDS].id) == "basis");
  HS_EXPECT_EQ(direct.schema[DIRECT_FIELDS].enum_count, 3);
  const In::OperatorDescriptor &ripple =
      *In::find_operator("sphere.displace.ripple.v2");
  HS_EXPECT_EQ(ripple.schema_count,
               PB::Surface::PeriodicRippleParams::FIELDS.size());
  HS_EXPECT_TRUE(std::string_view(ripple.schema[0].id) == "period");
  HS_EXPECT_EQ(ripple.schema[0].min, 30.0f);
  HS_EXPECT_EQ(ripple.schema[0].max, 143.0f);
  HS_EXPECT_EQ(ripple.schema[0].def, 80.0f);
  HS_EXPECT_EQ(ripple.schema[1].max, 0.15f);
  HS_EXPECT_EQ(ripple.schema[1].def, 0.15f);
  HS_EXPECT_EQ(ripple.schema[2].max, 5.0f);
  HS_EXPECT_EQ(ripple.schema[2].def, 0.1f);
  HS_EXPECT_EQ(ripple.schema[3].min, 0.7f);
  HS_EXPECT_EQ(ripple.schema[3].def, 0.7f);
  HS_EXPECT_TRUE(std::string_view(ripple.schema[5].id) == "center-polar");
  In::Op::RipplePhaseState ripple_state;
  PB::Surface::PeriodicRippleParams ripple_params;
  for (int frame = 0; frame < 40; ++frame)
    In::Op::DisplaceRipple::advance(ripple_state, ripple_params);
  HS_EXPECT_EQ(ripple_state.phase, 40.0f);
  for (int frame = 0; frame < 40; ++frame)
    In::Op::DisplaceRipple::advance(ripple_state, ripple_params);
  HS_EXPECT_EQ(ripple_state.phase, 0.0f);

  // Warp batch: every op is PLANE->PLANE with "speed" first; the polar chart
  // carries the full sixteen-harmonic list; curl-flow exposes basis and integrator.
  for (const char *id :
       {"warp.affine.v3", "warp.vortex.v2", "warp.wave-shear.v2",
        "warp.vector-noise.v2", "warp.mirror-tile.v2", "warp.polar-chart.v2",
        "warp.curl-flow.v2"}) {
    const In::OperatorDescriptor &warp = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(warp.input),
                 static_cast<int>(In::CarrierId::PLANE));
    HS_EXPECT_EQ(static_cast<int>(warp.output),
                 static_cast<int>(In::CarrierId::PLANE));
    HS_EXPECT_TRUE(std::string_view(warp.schema[0].id) == "speed");
  }
  const In::OperatorDescriptor &polar =
      *In::find_operator("warp.polar-chart.v2");
  constexpr size_t POLAR_FIELDS = In::Op::PolarChartParams::FIELDS.size();
  HS_EXPECT_EQ(polar.schema_count, POLAR_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(polar.schema[POLAR_FIELDS].id) == "mode");
  HS_EXPECT_EQ(polar.schema[POLAR_FIELDS].enum_count, 2);
  const In::ParamFieldInfo &harmonic = polar.schema[POLAR_FIELDS + 1];
  HS_EXPECT_TRUE(std::string_view(harmonic.id) == "harmonic");
  HS_EXPECT_EQ(harmonic.enum_count, PB::Warp::MAX_POLAR_HARMONIC);
  HS_EXPECT_TRUE(std::string_view(harmonic.enum_ids[0]) == "h1");
  HS_EXPECT_TRUE(std::string_view(harmonic.enum_ids[15]) == "h16");
  const In::OperatorDescriptor &shear =
      *In::find_operator("warp.wave-shear.v2");
  const In::ParamFieldInfo &shear_envelope =
      shear.schema[shear.schema_count - 1];
  HS_EXPECT_TRUE(std::string_view(shear_envelope.id) == "envelope");
  HS_EXPECT_EQ(shear_envelope.enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(shear_envelope.enum_ids[1]) ==
                 "projection-weight");
  const In::OperatorDescriptor &vector_noise =
      *In::find_operator("warp.vector-noise.v2");
  constexpr size_t VECTOR_NOISE_FIELDS =
      In::Op::VectorNoiseWarpParams::FIELDS.size();
  HS_EXPECT_EQ(vector_noise.schema_count, VECTOR_NOISE_FIELDS + 2);
  HS_EXPECT_TRUE(
      std::string_view(vector_noise.schema[VECTOR_NOISE_FIELDS].id) == "basis");
  HS_EXPECT_EQ(vector_noise.schema[VECTOR_NOISE_FIELDS].enum_count, 3);
  HS_EXPECT_TRUE(
      std::string_view(vector_noise.schema[VECTOR_NOISE_FIELDS + 1].id) ==
      "envelope");
  HS_EXPECT_EQ(vector_noise.schema[VECTOR_NOISE_FIELDS + 1].enum_count, 3);
  const In::OperatorDescriptor &curl_flow =
      *In::find_operator("warp.curl-flow.v2");
  constexpr size_t CURL_FLOW_FIELDS = In::Op::CurlFlowWarpParams::FIELDS.size();
  HS_EXPECT_EQ(curl_flow.schema_count, CURL_FLOW_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(curl_flow.schema[CURL_FLOW_FIELDS].id) ==
                 "basis");
  HS_EXPECT_EQ(curl_flow.schema[CURL_FLOW_FIELDS].enum_count, 3);
  const In::ParamFieldInfo &curl_integrator =
      curl_flow.schema[CURL_FLOW_FIELDS + 1];
  HS_EXPECT_TRUE(std::string_view(curl_integrator.id) == "integrator");
  HS_EXPECT_EQ(curl_integrator.enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(curl_integrator.enum_ids[2]) == "midpoint-4");
  const In::OperatorDescriptor &vortex = *In::find_operator("warp.vortex.v2");
  HS_EXPECT_EQ(vortex.schema_count, In::Op::VortexWarpParams::FIELDS.size());
  HS_EXPECT_TRUE(std::string_view(vortex.schema[1].id) == "center-x");

  // Projected sample batch: every crossing carries the union edge-width field
  // and the weight/coverage topology pair; only projected-noise adds a basis.
  for (const char *id :
       {"sample.grid.v3", "sample.twin-wave.v3", "sample.rings.v2",
        "sample.spiral.v2", "sample.lattice.v2", "sample.fractal.v2",
        "sample.tessellation.v2", "sample.projected-noise.v2"}) {
    const In::OperatorDescriptor &crossing = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(crossing.input),
                 static_cast<int>(In::CarrierId::PLANE));
    HS_EXPECT_EQ(static_cast<int>(crossing.output),
                 static_cast<int>(In::CarrierId::FIELD));
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
  const In::ParamFieldInfo *tessellation_kind =
      schema_field(*In::find_operator("sample.tessellation.v2"), "kind");
  HS_EXPECT_TRUE(tessellation_kind != nullptr);
  if (tessellation_kind != nullptr)
    HS_EXPECT_EQ(tessellation_kind->enum_count, 3);
  const In::ParamFieldInfo *projected_basis =
      schema_field(*In::find_operator("sample.projected-noise.v2"), "basis");
  HS_EXPECT_TRUE(projected_basis != nullptr);
  if (projected_basis != nullptr)
    HS_EXPECT_EQ(projected_basis->enum_count, 3);

  for (const char *id :
       {"sample.spherical-rings.v3", "sample.spherical-noise.v3"}) {
    const In::OperatorDescriptor &crossing = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(crossing.input),
                 static_cast<int>(In::CarrierId::SPHERE));
    HS_EXPECT_EQ(static_cast<int>(crossing.output),
                 static_cast<int>(In::CarrierId::FIELD));
    for (uint16_t field = 0; field < crossing.schema_count; ++field) {
      const std::string_view field_id{crossing.schema[field].id};
      HS_EXPECT_FALSE(field_id == "edge-width" || field_id == "weight-mode" ||
                      field_id == "coverage-mode");
    }
  }

  // Projection batch: the meridian-consuming projections extend the shared
  // family with central-meridian; stereographic, gnomonic, and the fast square
  // Peirce do not.
  for (const char *id :
       {"project.stereographic.v2", "project.folded-sinusoidal.v2",
        "project.equirectangular.v2", "project.gnomonic.v2",
        "project.peirce.v3", "project.peirce-square-fast.v3",
        "project.bonne.v3", "project.airocean.v3"}) {
    const In::OperatorDescriptor &projection = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(projection.input),
                 static_cast<int>(In::CarrierId::SPHERE));
    HS_EXPECT_EQ(static_cast<int>(projection.output),
                 static_cast<int>(In::CarrierId::PLANE));
    bool has_meridian = false;
    const In::ParamFieldInfo *frame_policy = nullptr;
    for (const In::ParamFieldInfo &field : projection.schema_span())
      if (std::string_view(field.id) == "frame")
        frame_policy = &field;
    HS_EXPECT_TRUE(frame_policy != nullptr);
    if (frame_policy != nullptr) {
      HS_EXPECT_TRUE(frame_policy->topology);
      HS_EXPECT_EQ(frame_policy->enum_count, 2);
      HS_EXPECT_EQ(frame_policy->enum_def, 1);
      HS_EXPECT_TRUE(std::string_view(frame_policy->enum_ids[0]) == "identity");
      HS_EXPECT_TRUE(std::string_view(frame_policy->enum_ids[1]) ==
                     "spin-wander");
    }
    for (uint16_t field = 0; field < projection.schema_count; ++field)
      has_meridian |=
          std::string_view(projection.schema[field].id) == "central-meridian";
    const std::string_view id_view{id};
    HS_EXPECT_EQ(has_meridian, id_view != "project.stereographic.v2" &&
                                   id_view != "project.gnomonic.v2" &&
                                   id_view != "project.peirce-square-fast.v3");
  }
  const In::OperatorDescriptor &gnomonic =
      *In::find_operator("project.gnomonic.v2");
  const In::ParamFieldInfo &gnomonic_hemisphere =
      gnomonic.schema[gnomonic.schema_count - 1];
  HS_EXPECT_TRUE(std::string_view(gnomonic_hemisphere.id) == "hemisphere");
  HS_EXPECT_EQ(gnomonic_hemisphere.enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(gnomonic_hemisphere.enum_ids[0]) == "folded");
  const In::OperatorDescriptor &bonne = *In::find_operator("project.bonne.v3");
  const In::ParamFieldInfo &bonne_hemisphere =
      bonne.schema[bonne.schema_count - 1];
  HS_EXPECT_TRUE(std::string_view(bonne_hemisphere.id) == "hemisphere");
  HS_EXPECT_EQ(bonne_hemisphere.enum_count, 2);
  HS_EXPECT_TRUE(std::string_view(bonne_hemisphere.enum_ids[0]) == "north");

  // Field batch: ridge is schema-free; iso-contour, smooth-bands and the value
  // cutout carry their own families.
  HS_EXPECT_EQ(In::find_operator("field.transfer.ridge.v2")->schema_count, 0);
  const In::OperatorDescriptor &iso =
      *In::find_operator("field.transfer.iso-contour.v2");
  HS_EXPECT_EQ(iso.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(iso.schema[0].id) == "iso-level");
  HS_EXPECT_TRUE(std::string_view(iso.schema[1].id) == "iso-width");
  const In::OperatorDescriptor &bands =
      *In::find_operator("field.transfer.smooth-bands.v2");
  HS_EXPECT_EQ(bands.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(bands.schema[0].id) == "band-count");
  HS_EXPECT_EQ(bands.schema[0].def, 4.0f);
  HS_EXPECT_TRUE(std::string_view(bands.schema[1].id) == "band-phase");
  const In::OperatorDescriptor &cutout =
      *In::find_operator("field.coverage.value-cutout.v2");
  HS_EXPECT_EQ(cutout.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(cutout.schema[0].id) == "cutout-threshold");
  HS_EXPECT_TRUE(std::string_view(cutout.schema[1].id) == "cutout-softness");
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
  static_assert(In::SLOT_SIZE == sizeof(PB::PlaneSample));
  static_assert(In::SLOT_ALIGN == alignof(PB::PlaneSample));
  static_assert(In::SLOT_SIZE >= sizeof(Color4));
  HS_EXPECT_EQ(In::SLOT_SIZE, 44u);
  HS_EXPECT_EQ(In::SLOT_ALIGN, 4u);
  HS_EXPECT_EQ(static_cast<int>(In::carrier_id_of<PB::SphereSample>()),
               static_cast<int>(In::CarrierId::SPHERE));
  HS_EXPECT_EQ(static_cast<int>(In::carrier_id_of<Color4>()),
               static_cast<int>(In::CarrierId::COLOR));
  constexpr uint32_t HASH = In::instance_hash("camera", "sphere.rotate.v2");
  static_assert(HASH == In::instance_hash("camera", "sphere.rotate.v2"));
  HS_EXPECT_NE(HASH, In::instance_hash("camera2", "sphere.rotate.v2"));
  HS_EXPECT_NE(HASH, In::instance_hash("camera", "project.stereographic.v2"));
  // The separator keeps ("ab", "c") and ("a", "bc") distinct.
  HS_EXPECT_NE(In::instance_hash("ab", "c"), In::instance_hash("a", "bc"));
}

inline std::string read_file(const char *path) {
  std::string content;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  std::FILE *file = std::fopen(path, "rb");
#pragma clang diagnostic pop
  if (file == nullptr)
    return content;
  char buffer[4096];
  size_t bytes;
  while ((bytes = std::fread(buffer, 1, sizeof(buffer), file)) > 0)
    content.append(buffer, bytes);
  std::fclose(file);
  return content;
}

inline void test_shader_chain_catalog_golden() {
  std::string catalog;
  In::append_catalog_json(catalog);
  catalog += '\n';
  const std::string golden = read_file(HS_SHADER_CHAIN_CATALOG_PATH);
  HS_EXPECT_FALSE(golden.empty());
  HS_EXPECT_TRUE(catalog == golden);
  if (catalog != golden)
    std::printf("shader_chain catalog drift: %zu generated vs %zu golden "
                "bytes; rebuild the regenerate_shader_chain_catalog target"
                "\n",
                catalog.size(), golden.size());
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
  HS_EXPECT_TRUE(contains("\"id\":\"sample.grid.v3\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.spherical-rings.v3\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.fractal.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.tessellation.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"colorize.generated-palette.v3\""));
  HS_EXPECT_TRUE(contains("\"id\":\"weight-mode\",\"topology\":true,"
                          "\"values\":[\"none\",\"projection\"],"
                          "\"default\":\"projection\""));
  HS_EXPECT_TRUE(contains("\"values\":[\"none\",\"weight\",\"weight-squared\","
                          "\"edge-fade\"]"));
  HS_EXPECT_TRUE(contains("\"approximate\":true"));
  // The cataloged block layouts are the real ABI figures.
  std::string sample_blocks = "\"blocks\":{\"param\":{\"size\":";
  sample_blocks += std::to_string(sizeof(In::Op::GridSampleParams));
  HS_EXPECT_TRUE(contains(sample_blocks.c_str()));
}

inline void test_shader_chain_default_chain_renders() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  HS_EXPECT_FALSE(program.compiled());
  arm_default_chain(program, 3, ValueSet::DEFAULTS);
  HS_EXPECT_TRUE(program.compiled());
  HS_EXPECT_EQ(program.ops().size(), 4u);
  HS_EXPECT_TRUE(std::string_view(program.ops()[0].instance) == "camera");
  HS_EXPECT_TRUE(std::string_view(program.ops()[3].instance) == "colorize");
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
  const auto verify_extended_crossing = []<typename Model>() {
    typename Model::Params params;
    const auto before = params;
    const auto &descriptor = *In::find_operator(Model::ID);
    constexpr size_t FIELDS = Model::Params::FIELDS.size();
    auto *coverage = static_cast<uint8_t *>(
        descriptor.runtime.param_address(&params, FIELDS + 1));
    *coverage = static_cast<uint8_t>(In::Op::ProjectionCoverageMode::NONE);
    HS_EXPECT_EQ(params.coverage_mode, 0);
    for (const auto &field : Model::Params::FIELDS)
      HS_EXPECT_EQ(params.*field.member, before.*field.member);
    const auto *drift = std::find_if(
        descriptor.schema, descriptor.schema + descriptor.schema_count,
        [](const auto &field) {
          return std::string_view(field.id) == "drift";
        });
    HS_EXPECT_TRUE(drift != descriptor.schema + descriptor.schema_count);
    if (drift != descriptor.schema + descriptor.schema_count)
      HS_EXPECT_EQ(drift->max, 2.0f);
  };
  verify_extended_crossing.template operator()<In::Op::SampleGridV3>();
  verify_extended_crossing.template operator()<In::Op::SampleTwinWaveV3>();
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 1, ValueSet::DEFAULTS);
  const In::OperatorDescriptor &sample_op = *program.ops()[2].op;
  uint8_t *block = program.param_block(2);
  auto &typed = *reinterpret_cast<In::Op::GridSampleParams *>(block);
  // Every schema index resolves inside the block, floats first.
  for (uint16_t index = 0; index < sample_op.schema_count; ++index) {
    const uint8_t *address = static_cast<const uint8_t *>(
        sample_op.runtime.param_address(block, index));
    HS_EXPECT_TRUE(address >= block &&
                   address < block + sample_op.runtime.param.size);
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
  const math::Vector view = math::Vector(0, 0, -1);
  const Color4 faded = program.evaluate(view, ctx);
  *coverage_mode = static_cast<uint8_t>(In::Op::ProjectionCoverageMode::NONE);
  program.prepare(ctx);
  const Color4 full = program.evaluate(view, ctx);
  HS_EXPECT_LT(faded.alpha, full.alpha);
  HS_EXPECT_LT(faded.alpha, 0.5f);
  program.clear();
}

inline void test_shader_chain_parity_rotate_project() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_default_chain(program, 5, set);
    const In::FrameContext ctx = shared_resources().context();
    program.prepare(ctx);
    // The erased prepared conjugates equal the same composition recomputed
    // from the state blocks.
    const MirrorFrame mirror = mirror_from(program, ctx);
    const auto &camera_prepared =
        *reinterpret_cast<const In::Op::Rotate::Prepared *>(
            program.prepared_block(0));
    const auto &project_prepared =
        *reinterpret_cast<const In::Op::ProjectStereographic::Prepared *>(
            program.prepared_block(1));
    HS_EXPECT_EQ(std::memcmp(&camera_prepared.conjugate,
                             &mirror.camera_conjugate,
                             sizeof(math::Quaternion)),
                 0);
    HS_EXPECT_EQ(std::memcmp(&project_prepared.conjugate,
                             &mirror.projection_conjugate,
                             sizeof(math::Quaternion)),
                 0);
    // Op-level kernel parity for the two sphere-family ops.
    using BoundRotate = PB::Stage::Rotate<MirrorCamera>::Bind<MirrorBinding>;
    using BoundProject = PB::Stage::Project<
        PB::Projection::Stereographic<MirrorProjection>>::Bind<MirrorBinding>;
    const In::OperatorDescriptor &rotate_op = *program.ops()[0].op;
    const In::OperatorDescriptor &project_op = *program.ops()[1].op;
    for (const math::Vector &view : sweep_views()) {
      const PB::SphereSample seed{view, 0.0f};
      alignas(In::SLOT_ALIGN) uint8_t erased_out[In::SLOT_SIZE];
      rotate_op.runtime.run(&seed, erased_out, ctx, program.param_block(0),
                            program.prepared_block(0));
      const auto &erased_rotated =
          *std::launder(reinterpret_cast<PB::SphereSample *>(erased_out));
      const PB::SphereSample reference_rotated =
          BoundRotate::run(seed, mirror, {});
      HS_EXPECT_TRUE(sphere_identical(erased_rotated, reference_rotated));
      alignas(In::SLOT_ALIGN) uint8_t projected_out[In::SLOT_SIZE];
      project_op.runtime.run(&erased_rotated, projected_out, ctx,
                             program.param_block(1), program.prepared_block(1));
      const auto &erased_projected =
          *std::launder(reinterpret_cast<PB::PlaneSample *>(projected_out));
      const PB::PlaneSample reference_projected =
          BoundProject::run(reference_rotated, mirror, {});
      HS_EXPECT_TRUE(plane_identical(erased_projected, reference_projected));
    }
    program.clear();
  }
}
