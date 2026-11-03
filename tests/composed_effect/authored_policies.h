/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file authored_policies.h
 * @brief Authored displacement, projection and modulation bindings. */

template <typename Family>
using AuthoredSurfaceParams = Pullback::ComposedDetail::ParameterSet<
    Pullback::ComposedDetail::ResourceList<Pullback::ParameterResource<
        "surface", Family, Pullback::ResourceKind::SURFACE>>>;

template <typename Family>
using AuthoredSurfaceFrame =
    Pullback::FrameIntervalState<AuthoredSurfaceParams<Family>>;

template <typename Family>
using AuthoredSurfaceProvider =
    Pullback::SurfaceProvider<Pullback::Binding<AuthoredSurfaceFrame<Family>>,
                              Family, true>;

/** @brief Composed tangent policies retain authored phase and integration settings. */
inline void test_authored_tangent_bindings_match_operators() {
  using namespace Pullback;
  using namespace Surface;
  const auto INPUT = math::Vector(0.31f, 0.57f, -0.76f).normalized();
  for (uint8_t integrator = 0; integrator < 3; ++integrator) {
    for (float phase : {-2.3f, 0.0f, 1.7f}) {
      AuthoredSurfaceFrame<RossbyHaurwitzParams> frame{};
      auto &params = frame.params.template get<"surface">();
      params = {.strength = -0.17f,
                .omega = 0.3f,
                .wave = 0.7f,
                .order = 4,
                .drift = 0.03f,
                .integrator = integrator};
      frame.resources.template get<"surface">().phase = phase;
      using Policy =
          RossbyHaurwitz<AuthoredSurfaceProvider<RossbyHaurwitzParams>>;
      const auto ACTUAL = Policy::apply(INPUT, frame, Policy::prepare(frame));
      Interp::Op::RossbyHaurwitzParams reference;
      reference.strength = params.strength;
      reference.omega = params.omega;
      reference.wave = params.wave;
      reference.order = params.order;
      const auto EXPECTED = Interp::Op::DisplaceRossbyHaurwitz::step(
          INPUT, Interp::Op::DisplaceRossbyHaurwitz::field(reference, phase),
          reference.strength, integrator);
      HS_EXPECT_NEAR(ACTUAL.sphere.x, EXPECTED.sphere.x, 1e-6f);
      HS_EXPECT_NEAR(ACTUAL.sphere.y, EXPECTED.sphere.y, 1e-6f);
      HS_EXPECT_NEAR(ACTUAL.sphere.z, EXPECTED.sphere.z, 1e-6f);
      HS_EXPECT_NEAR(ACTUAL.path_length, EXPECTED.path_length, 1e-6f);

      for (uint8_t symmetry = 0; symmetry < 9; ++symmetry)
        for (uint8_t kind = 0; kind < 2; ++kind) {
          AuthoredSurfaceFrame<SymmetricParams> symmetric{};
          auto &controls = symmetric.params.template get<"surface">();
          controls = {.strength = 0.23f,
                      .speed = -0.02f,
                      .integrator = integrator,
                      .symmetry = symmetry,
                      .kind = kind};
          symmetric.resources.template get<"surface">().phase = phase;
          using Symmetry = Symmetric<AuthoredSurfaceProvider<SymmetricParams>>;
          const auto OUTPUT =
              Symmetry::apply(INPUT, symmetric, Symmetry::prepare(symmetric));
          Interp::Op::SymmetricDisplaceParams baseline;
          baseline.symmetry = symmetry;
          baseline.kind = kind;
          const auto REFERENCE = Interp::Op::DisplaceSymmetric::step(
              INPUT, Interp::Op::DisplaceSymmetric::field(baseline, phase),
              0.23f, integrator);
          HS_EXPECT_NEAR(OUTPUT.sphere.x, REFERENCE.sphere.x, 1e-6f);
          HS_EXPECT_NEAR(OUTPUT.sphere.y, REFERENCE.sphere.y, 1e-6f);
          HS_EXPECT_NEAR(OUTPUT.sphere.z, REFERENCE.sphere.z, 1e-6f);
          HS_EXPECT_NEAR(OUTPUT.path_length, REFERENCE.path_length, 1e-6f);
        }
    }
  }
  AuthoredSurfaceFrame<PotentialParams> potential{};
  auto &params = potential.params.template get<"surface">();
  params.orientation_azimuth = 0.47f;
  params.orientation_polar = 1.24f;
  params.terms[0] = {0, 1, 3, 2, 0.8f, -0.3f};
  params.terms[1] = {2, 0, 1, 0, -0.2f, 0};
  potential.resources.template get<"surface">().phase = -0.73f;
  using Policy = Potential<AuthoredSurfaceProvider<PotentialParams>>;
  Interp::Op::PotentialDisplaceParams reference;
  reference.orientation_azimuth = params.orientation_azimuth;
  reference.orientation_polar = params.orientation_polar;
  for (size_t i = 0; i < params.terms.size(); ++i) {
    const auto &term = params.terms[i];
    reference.terms[i] = {term.type, term.kind,      term.l,
                          term.m,    term.amplitude, term.phase};
  }
  const auto ACTUAL =
      Policy::apply(INPUT, potential, Policy::prepare(potential));
  const auto EXPECTED = Interp::Op::DisplacePotential::step(
      INPUT, Interp::Op::DisplacePotential::field(reference, -0.73f),
      params.strength, 0);
  HS_EXPECT_NEAR(ACTUAL.sphere.x, EXPECTED.sphere.x, 1e-6f);
  HS_EXPECT_NEAR(ACTUAL.sphere.y, EXPECTED.sphere.y, 1e-6f);
  HS_EXPECT_NEAR(ACTUAL.sphere.z, EXPECTED.sphere.z, 1e-6f);
  HS_EXPECT_NEAR(ACTUAL.path_length, EXPECTED.path_length, 1e-6f);
}

/** @brief Tangent fields interpolate drift across the radian wrap seam. */
inline void test_authored_tangent_interval_uses_radians() {
  using Frame = AuthoredSurfaceFrame<Pullback::SymmetricParams>;
  Frame from{}, to{};
  from.resources.template get<"surface">().phase = math::TWO_PI_F - 0.04f;
  to.resources.template get<"surface">().phase = 0.06f;
  const Pullback::FrameInterval<Frame> INTERVAL(from, to);
  const auto MIDDLE = INTERVAL.at(0.5f);
  HS_EXPECT_NEAR(MIDDLE.resources.template get<"surface">().phase,
                 math::TWO_PI_F + 0.01f, 1e-6f);
}

/** @brief Potential interpolation keeps signed angles and snaps harmonic topology. */
inline void test_authored_potential_interpolation_and_admission() {
  Pullback::PotentialParams from, to;
  from.orientation_azimuth = 3.0f;
  to.orientation_azimuth = -3.0f;
  from.terms[0].phase = 3.0f;
  to.terms[0].phase = -3.0f;
  from.terms[0].amplitude = 0.4f;
  to.terms[0].amplitude = 0.9f;
  const auto MIDDLE = Pullback::interpolate(from, to, 0.5f);
  HS_EXPECT_TRUE(Pullback::valid(MIDDLE));
  HS_EXPECT_NEAR(fabsf(MIDDLE.orientation_azimuth), math::PI_F, 1e-6f);
  HS_EXPECT_NEAR(fabsf(MIDDLE.terms[0].phase), math::PI_F, 1e-6f);
  HS_EXPECT_EQ(MIDDLE.terms[0].amplitude, 0.4f);
  to.terms[0].m = 2;
  HS_EXPECT_FALSE(Pullback::valid(to));
  to.terms[0].l = 3;
  HS_EXPECT_TRUE(Pullback::valid(to));
  HS_EXPECT_EQ(Pullback::interpolate(from, to, 0.5f).terms[0].m, uint8_t(0));
  HS_EXPECT_EQ(Pullback::interpolate(from, to, 1.0f).terms[0].m, uint8_t(2));
}

/** @brief Authored noise policies honor every basis, stepper and limit. */
inline void test_authored_surface_noise_topology() {
  using namespace Pullback;
  FastNoiseLite noise;
  init_effect_noise(noise, 137);
  const auto INPUT = math::Vector(-0.37f, 0.58f, 0.72f).normalized();
  for (uint8_t basis = 0; basis < 3; ++basis) {
    AuthoredSurfaceFrame<AuthoredDirectSurfaceParams> direct{};
    auto &p = direct.params.template get<"surface">();
    p.scale = 2.3f;
    p.strength = 0.11f;
    p.direction = 0.31f;
    p.basis = basis;
    auto &resource = direct.resources.template get<"surface">();
    resource.phase = 0.27f;
    resource.noise = &noise;
    using Directed = Surface::ConfiguredDirect<
        AuthoredSurfaceProvider<AuthoredDirectSurfaceParams>>;
    const auto PREPARED = Directed::prepare(direct);
    const auto ACTUAL = Directed::apply(INPUT, direct, PREPARED);
    const auto EXPECTED = Surface::direct_noise(
        INPUT, noise, static_cast<math::NoiseBasis>(basis), p.scale,
        PREPARED.loop_offset, p.strength, PREPARED.direction_cos,
        PREPARED.direction_sin, true);
    HS_EXPECT_NEAR(ACTUAL.sphere.x, EXPECTED.sphere.x, 1e-6f);
    HS_EXPECT_NEAR(ACTUAL.sphere.y, EXPECTED.sphere.y, 1e-6f);
    HS_EXPECT_NEAR(ACTUAL.path_length, EXPECTED.path_length, 1e-6f);
    for (uint8_t integrator = 0; integrator < 3; ++integrator)
      for (uint8_t limit = 0; limit < 2; ++limit) {
        AuthoredSurfaceFrame<AuthoredSurfaceNoiseParams> curl{};
        auto &params = curl.params.template get<"surface">();
        params.scale = 2.1f;
        params.strength = -0.19f;
        params.basis = basis;
        params.integrator = integrator;
        params.limit = limit;
        curl.resources.template get<"surface">().phase = 0.27f;
        curl.resources.template get<"surface">().noise = &noise;
        using Curled = Surface::ConfiguredCurl<
            AuthoredSurfaceProvider<AuthoredSurfaceNoiseParams>>;
        const auto PREPARE = Curled::prepare(curl);
        const auto OUTPUT = Curled::apply(INPUT, curl, PREPARE);
        const auto REFERENCE = Surface::curl_noise(
            INPUT, noise, static_cast<math::NoiseBasis>(basis),
            static_cast<Surface::Integrator>(integrator),
            static_cast<math::TangentLimit>(limit), params.scale,
            PREPARE.loop_offset, params.strength, true);
        HS_EXPECT_NEAR(OUTPUT.sphere.x, REFERENCE.sphere.x, 1e-6f);
        HS_EXPECT_NEAR(OUTPUT.sphere.y, REFERENCE.sphere.y, 1e-6f);
        HS_EXPECT_NEAR(OUTPUT.sphere.z, REFERENCE.sphere.z, 1e-6f);
        HS_EXPECT_NEAR(OUTPUT.path_length, REFERENCE.path_length, 1e-6f);
      }
  }
}

/** @brief Flexible modulation preserves spatial storage across source changes. */
inline void test_authored_dynamic_modulation_storage() {
  using namespace Pullback;
  using Family = Modulation::Params<Pattern::SourceKind::NOISE,
                                    Pattern::Evaluation::BAKED_CUBE, true>;
  using Resource = ParameterResource<"layer", Family, ResourceKind::MODULATION>;
  ComposedDetail::ResourceNoise<Resource> storage;
  Family params;
  params.modulation.source = Pattern::SourceKind::DISPLACEMENT;
  params.modulation.evaluation = Pattern::Evaluation::DIRECT;
  HS_EXPECT_TRUE(valid(params));
  const auto &DISPLACEMENT =
      storage.modulation_resources.prepare(params.modulation, {});
  HS_EXPECT_FALSE(DISPLACEMENT.cube.active);
  params.modulation.source = Pattern::SourceKind::NOISE;
  params.modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
  HS_EXPECT_TRUE(valid(params));
  const auto &NOISE =
      storage.modulation_resources.prepare(params.modulation, {});
  HS_EXPECT_TRUE(NOISE.cube.active);
  HS_EXPECT_TRUE(NOISE.evaluate != nullptr);
  params.modulation.noise.scale = -1.0f;
  HS_EXPECT_FALSE(valid(params));
  using Minimal = Modulation::Params<Pattern::SourceKind::DISPLACEMENT,
                                     Pattern::Evaluation::DIRECT, true>;
  using MinimalResource =
      ParameterResource<"minimal", Minimal, ResourceKind::MODULATION>;
  ComposedDetail::ResourceNoise<MinimalResource> flexible;
  ComposedDetail::ResourceModulationHistory<MinimalResource> retained;
  Minimal switches;
  HS_EXPECT_TRUE(valid(switches));
  switches.modulation.source = Pattern::SourceKind::NOISE;
  switches.modulation.evaluation = Pattern::Evaluation::BAKED_CUBE;
  HS_EXPECT_TRUE(valid(switches));
  const auto &BAKED =
      flexible.modulation_resources.prepare(switches.modulation, {});
  HS_EXPECT_TRUE(BAKED.cube.active);
  const auto *COPY =
      Modulation::retain(&BAKED, retained.modulation, retained.table);
  HS_EXPECT_TRUE(COPY->cube.data == retained.table.data());
  HS_EXPECT_TRUE(COPY->cube.data != BAKED.cube.data);
  HS_EXPECT_EQ(retained.table.size(), size_t(Pattern::SphereLutView::SIZE));
}

template <Pullback::ProjectionKind Kind, typename Operator>
inline void verify_authored_projection() {
  using namespace Pullback;
  using Params = ComposedDetail::ParameterSet<
      ComposedDetail::ResourceList<ParameterResource<
          "projection", AtlasProjectionParams, ResourceKind::PROJECTION>>>;
  using Frame = FrameState<Params>;
  using Provider = ProjectionProvider<Binding<Frame>, AtlasProjectionParams>;
  using Policy = AuthoredProjection<Kind, Provider>;
  Frame frame{};
  auto &params = frame.params.template get<"projection">();
  typename Operator::Params chain;
  params.singularity_fade = chain.singularity_fade = 2.7f;
  if constexpr (requires { chain.central_meridian; })
    params.central_meridian = chain.central_meridian = 0.41f;
  if constexpr (requires { chain.coordinate_scale; })
    params.coordinate_scale = chain.coordinate_scale = 1.7f;
  if constexpr (requires { chain.layout_scroll; })
    params.layout_scroll = chain.layout_scroll = -0.37f;
  if constexpr (requires { chain.standard_parallel; })
    params.standard_parallel = chain.standard_parallel = 0.63f;
  params.frame = chain.frame = 0;
  const int LAYOUTS = Kind == ProjectionKind::PEIRCE     ? 4
                      : Kind == ProjectionKind::AIROCEAN ? 2
                                                         : 1;
  const int HEMISPHERES = Kind == ProjectionKind::GNOMONIC_FOLDED ? 3
                          : Kind == ProjectionKind::BONNE         ? 2
                                                                  : 1;
  for (int layout = 0; layout < LAYOUTS; ++layout)
    for (int hemisphere = 0; hemisphere < HEMISPHERES; ++hemisphere) {
      params.layout = static_cast<uint8_t>(layout);
      params.hemisphere = static_cast<uint8_t>(hemisphere);
      if constexpr (requires { chain.layout; })
        chain.layout = params.layout;
      if constexpr (requires { chain.hemisphere; })
        chain.hemisphere = params.hemisphere;
      const auto PREPARED = Policy::prepare(frame);
      const auto REFERENCE = Operator::prepare({}, chain, {});
      for (const auto &view :
           {math::Vector(0.31f, 0.57f, -0.76f).normalized(),
            math::Vector(-0.58f, -0.69f, 0.42f).normalized()}) {
        const auto ACTUAL = Policy::project(view, frame, PREPARED);
        const auto EXPECTED = [&] {
          if constexpr (requires { Operator::project(view, chain, REFERENCE); })
            return Operator::project(view, chain, REFERENCE);
          else
            return Operator::project(view, chain);
        }();
        HS_EXPECT_NEAR(ACTUAL.coords.re, EXPECTED.coords.re, 1e-6f);
        HS_EXPECT_NEAR(ACTUAL.coords.im, EXPECTED.coords.im, 1e-6f);
        HS_EXPECT_NEAR(ACTUAL.provenance.value_weight,
                       EXPECTED.provenance.value_weight, 1e-6f);
        HS_EXPECT_NEAR(ACTUAL.provenance.fade_edge_distance,
                       EXPECTED.provenance.fade_edge_distance, 1e-6f);
        HS_EXPECT_EQ(ACTUAL.provenance.region_id,
                     EXPECTED.provenance.region_id);
        HS_EXPECT_EQ(ACTUAL.provenance.component_id,
                     EXPECTED.provenance.component_id);
        HS_EXPECT_EQ(ACTUAL.provenance.boundary_flags,
                     EXPECTED.provenance.boundary_flags);
      }
    }
  frame.projection_conjugate = math::make_rotation(math::Y_AXIS, 0.7f);
  const auto VIEW = math::Vector(0.3f, 0.2f, 0.8f).normalized();
  const auto IDENTITY = math::rotate(VIEW, Policy::frame_conjugate(frame));
  HS_EXPECT_NEAR(IDENTITY.x, VIEW.x, 1e-6f);
  HS_EXPECT_NEAR(IDENTITY.z, VIEW.z, 1e-6f);
  params.frame = 1;
  const auto ANIMATED = math::rotate(VIEW, Policy::frame_conjugate(frame));
  HS_EXPECT_GT(fabsf(ANIMATED.x - VIEW.x), 0.1f);
}

/** @brief Authored projection families preserve every catalog projection layout. */
inline void test_authored_projection_bindings_match_operators() {
  using namespace Pullback;
  using namespace Interp::Op;
  verify_authored_projection<ProjectionKind::STEREOGRAPHIC,
                             ProjectStereographic>();
  verify_authored_projection<ProjectionKind::FOLDED_SINUSOIDAL,
                             ProjectFoldedSinusoidal>();
  verify_authored_projection<ProjectionKind::EQUIRECTANGULAR,
                             ProjectEquirectangular>();
  verify_authored_projection<ProjectionKind::GNOMONIC_FOLDED,
                             ProjectGnomonic>();
  verify_authored_projection<ProjectionKind::PEIRCE, ProjectPeirceV3>();
  verify_authored_projection<ProjectionKind::PEIRCE_SQUARE_FAST,
                             ProjectPeirceSquareFastV3>();
  verify_authored_projection<ProjectionKind::BONNE, ProjectBonneV3>();
  verify_authored_projection<ProjectionKind::AIROCEAN, ProjectAiroceanV3>();
  verify_authored_projection<ProjectionKind::S2, ProjectS2>();
}

/** @brief Keyed lenses preserve every symmetry and authored angle multiplier. */
inline void test_authored_lens_bindings_match_operators() {
  using namespace Pullback;
  using Params = ComposedDetail::ParameterSet<ComposedDetail::ResourceList<
      ParameterResource<"glitch", Lens::GlitchParams, ResourceKind::LENS>,
      ParameterResource<"symmetry", Lens::KaleidoscopeParams,
                        ResourceKind::LENS>>>;
  using Frame = FrameState<Params>;
  using Glitch = Lens::ConfiguredGlitch<
      KeyedLensProvider<Binding<Frame>, Lens::GlitchParams, "glitch">>;
  using Kaleidoscope = Lens::ConfiguredKaleidoscope<
      KeyedLensProvider<Binding<Frame>, Lens::KaleidoscopeParams, "symmetry">>;
  Frame frame{};
  const auto INPUT = math::Vector(0.39f, -0.57f, 0.72f).normalized();
  for (uint8_t symmetry = 0; symmetry < 9; ++symmetry) {
    frame.params.template get<"symmetry">().symmetry = symmetry;
    const auto ACTUAL = Kaleidoscope::apply(INPUT, frame);
    const auto EXPECTED = Interp::Op::kaleidoscope_lens(INPUT, symmetry);
    HS_EXPECT_NEAR(ACTUAL.x, EXPECTED.x, 1e-6f);
    HS_EXPECT_NEAR(ACTUAL.y, EXPECTED.y, 1e-6f);
    HS_EXPECT_NEAR(ACTUAL.z, EXPECTED.z, 1e-6f);
  }
  for (float polar : {2.0f, -1.5f})
    for (float azimuth : {3.0f, 4.25f}) {
      auto &params = frame.params.template get<"glitch">();
      params.polar_multiplier = polar;
      params.azimuth_multiplier = azimuth;
      const auto ACTUAL = Glitch::apply(INPUT, frame);
      Interp::Op::GlitchChainParams reference;
      reference.polar_multiplier = polar;
      reference.azimuth_multiplier = azimuth;
      const auto EXPECTED =
          Interp::Op::LensGlitch::run({INPUT, 0}, {}, reference, {}).dir;
      HS_EXPECT_NEAR(ACTUAL.x, EXPECTED.x, 1e-6f);
      HS_EXPECT_NEAR(ACTUAL.y, EXPECTED.y, 1e-6f);
      HS_EXPECT_NEAR(ACTUAL.z, EXPECTED.z, 1e-6f);
    }
}

/** @brief Composed fixture with a flexible modulation recipe. */
struct DynamicModulationControlsSpec : Pullback::Spec {
  static constexpr bool ANIMATED_PROJECTION = false;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Displace<Pullback::Surface::PeriodicRipple<
          Pullback::SurfaceProvider<B, Pullback::PeriodicRippleParams, true>>>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::STEREOGRAPHIC, B>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
          Pullback::GridSourceParams, B>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>,
      Pullback::Stage::Modulate<
          Pullback::Color::Modulate<Pullback::ModulationProvider<
              B, "modulation", Pattern::SourceKind::NOISE,
              Pattern::Evaluation::BAKED_CUBE, false, true>>>>;
};

/** @brief Native controls and snapshot adoption for flexible modulation. */
struct DynamicModulationControlsProbe
    : Pullback::ComposedEffect<8, 8, DynamicModulationControlsProbe,
                               DynamicModulationControlsSpec,
                               Control::HookOwnership::INHERITED> {
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr auto PRESETS = [] {
    std::array<PresetEntry<Params>, 2> rows{};
    rows[0].params.template get<"modulation">().modulation.mechanic.amount =
        0.5f;
    auto &config = rows[1].params.template get<"modulation">().modulation;
    config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
    config.evaluation = Pattern::Evaluation::DIRECT;
    return rows;
  }();
  /** @brief Live recipe controls.
   * @return Mutable canonical configuration. */
  auto &configuration() {
    return params.template get<"modulation">().modulation;
  }
};

/** @brief Dynamic modulation registers canonical controls and refreshes admission. */
inline void test_authored_dynamic_modulation_controls() {
  using namespace Pullback;
  hs_test::reset_globals();
  DynamicModulationControlsProbe effect;
  effect.init();
  for (size_t index = 0; index < Modulation::Controls::EXTRA_SCHEMA.size();
       ++index) {
    const auto &field = Modulation::Controls::EXTRA_SCHEMA[index];
    const auto *parameter = effect.getParameters().find(field.id);
    HS_EXPECT_TRUE(parameter != nullptr);
    if (!parameter)
      continue;
    HS_EXPECT_TRUE(parameter->target == Modulation::Controls::address(
                                            effect.configuration(), index));
    HS_EXPECT_EQ(parameter->option_count, int(field.enum_count));
    HS_EXPECT_EQ(parameter->readonly,
                 !Modulation::Controls::active(effect.configuration(), index));
    HS_EXPECT_TRUE(parameter->animated);
  }
  HS_EXPECT_EQ(effect.updateParameter("modulation-evaluation", 0),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("modulation-source", 1),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.getParameters().find("pattern-scale")->readonly);
  HS_EXPECT_FALSE(effect.getParameters().find("sine-frequency")->readonly);
  HS_EXPECT_EQ(effect.updateParameter("sine-frequency", 2.3f),
               ParamSetResult::APPLIED);
  HS_EXPECT_NEAR(effect.configuration().sine.frequency, 2.3f, 1e-6f);
  HS_EXPECT_EQ(effect.updateParameter("pattern-scale", 2.3f),
               ParamSetResult::READONLY);
  HS_EXPECT_EQ(effect.updateParameter("modulation-source", 0),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(effect.getParameters().find("pattern-noise-basis")->readonly);
  HS_EXPECT_EQ(effect.updateParameter("pattern-noise-variant", 1),
               ParamSetResult::APPLIED);
  HS_EXPECT_FALSE(effect.getParameters().find("pattern-noise-basis")->readonly);
  HS_EXPECT_EQ(effect.updateParameter("pattern-noise-basis", 2),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("pattern-contrast", 0.7f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.configuration().noise.basis,
               Pattern::NoiseBasis::RIDGED3);
  HS_EXPECT_EQ(effect.updateParameter("modulation-evaluation", 1),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("modulation-source", 1),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.configuration().source, Pattern::SourceKind::NOISE);
  auto snapshot = effect.serialize_parameters();
  auto &recipe = snapshot.params.template get<"modulation">().modulation;
  recipe.source = Pattern::SourceKind::SPHERE_INTERFERENCE;
  recipe.evaluation = Pattern::Evaluation::DIRECT;
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
  HS_EXPECT_TRUE(effect.getParameters().find("pattern-scale")->readonly);
  HS_EXPECT_FALSE(
      effect.getParameters().find("interference-frequency-b")->readonly);
  HS_EXPECT_EQ(effect.updateParameter("interference-frequency-b", 4.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.configuration().interference.frequency_b, 4.0f);
  HS_EXPECT_TRUE(
      effect.edit_document_modulation("modulation", [](auto &config) {
        config.source = Pattern::SourceKind::DISPLACEMENT;
        return ParamSetResult::APPLIED;
      }) == ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(
      effect.getParameters().find("interference-frequency-b")->readonly);
}
