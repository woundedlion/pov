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

/** @brief Parameter interpolation preserves both alias directions and recipe topology. */
inline void test_authored_interpolation_destination_aliases() {
  using namespace Pullback;
  using Mod = Modulation::Params<Pattern::SourceKind::NOISE,
                                 Pattern::Evaluation::BAKED_CUBE, true>;
  using Params = ComposedDetail::ParameterSet<ComposedDetail::ResourceList<
      ParameterResource<"a", Mod, ResourceKind::MODULATION>,
      ParameterResource<"b", Mod, ResourceKind::MODULATION>,
      ParameterResource<"surface", Surface::PotentialParams,
                        ResourceKind::SURFACE>>>;
  Params from{}, to{};
  from.template get<"a">().modulation.noise.scale = 1.0f;
  to.template get<"a">().modulation.noise.scale = 16.0f;
  from.template get<"a">().modulation.mechanic.amount = 0.0f;
  to.template get<"a">().modulation.mechanic.amount = 2.0f;
  to.template get<"b">().modulation.source = Pattern::SourceKind::DISPLACEMENT;
  to.template get<"b">().modulation.evaluation = Pattern::Evaluation::DIRECT;
  from.template get<"surface">().orientation_azimuth = 3.0f;
  to.template get<"surface">().orientation_azimuth = -3.0f;
  from.template get<"surface">().terms[0].amplitude = 0.25f;
  to.template get<"surface">().terms[0].amplitude = 0.75f;
  for (const float progress : {0.0f, 0.5f, 1.0f}) {
    for (const int alias : {0, 1, 2}) {
      Params a = from, b = to, destination{};
      Params &out = alias == 1 ? a : alias == 2 ? b : destination;
      interpolate_into(out, a, b, progress);
      const float scales[] = {1.0f, 4.0f, 16.0f};
      HS_EXPECT_NEAR(out.template get<"a">().modulation.noise.scale,
                     scales[static_cast<int>(progress * 2)], 1e-6f);
      HS_EXPECT_NEAR(out.template get<"a">().modulation.mechanic.amount,
                     2.0f * progress, 0.0f);
      HS_EXPECT_EQ(out.template get<"b">().modulation.source,
                   Pattern::SourceKind::DISPLACEMENT);
      HS_EXPECT_EQ(out.template get<"b">().modulation.evaluation,
                   Pattern::Evaluation::DIRECT);
      const float azimuth = progress == 0.0f   ? 3.0f
                            : progress == 1.0f ? -3.0f
                                               : math::PI_F;
      HS_EXPECT_NEAR(
          remainderf(out.template get<"surface">().orientation_azimuth -
                         azimuth,
                     math::TWO_PI_F),
          0.0f, 1e-6f);
      HS_EXPECT_NEAR(out.template get<"surface">().terms[0].amplitude,
                     progress < 1.0f ? 0.25f : 0.75f, 0.0f);
    }
  }
}

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

/** @brief Chained polynomial displacements preserve approximate lens norms. */
inline void test_authored_tangent_fields_after_twist() {
  using namespace Pullback::Surface;
  const auto rossby = rossby_haurwitz_field(1.3882353f, 4.0f, 3, 0.0f);
  const auto symmetric = symmetric_field(3, 1, 0.0f);
  float worst_norm_error = 0.0f;
  float worst_tangent_residual = 0.0f;
  for (uint8_t integrator = 0; integrator < 3; ++integrator)
    for (int latitude = 1; latitude < 64; ++latitude)
      for (int longitude = 0; longitude < 128; ++longitude) {
        const float polar = latitude * math::PI_F / 64;
        const float azimuth = longitude * math::TWO_PI_F / 128;
        const math::Vector input(sinf(polar) * cosf(azimuth), cosf(polar),
                                 sinf(polar) * sinf(azimuth));
        const auto lensed = lenses::twist_lens(input, -10.259509f, 0.999f);
        const auto first =
            step_tangent_field(lensed, rossby, 0.09239978f, integrator, true);
        const auto second = step_tangent_field(first.sphere, symmetric,
                                               -1.0657501f, integrator, true);
        worst_norm_error =
            std::max(worst_norm_error,
                     fabsf(math::dot(second.sphere, second.sphere) - 1));
        worst_tangent_residual = std::max(
            worst_tangent_residual,
            fabsf(math::dot(first.sphere, symmetric.velocity(first.sphere))));
        HS_EXPECT_TRUE(std::isfinite(second.path_length));
      }
  HS_EXPECT_LT(worst_norm_error, Pullback::DIRECTION_NORM_SQUARED_TOLERANCE);
  HS_EXPECT_LT(worst_tangent_residual, 1e-6f);
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
    HS_EXPECT_NEAR(ACTUAL.sphere.z, EXPECTED.sphere.z, 1e-6f);
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
  if constexpr (requires { chain.area_correction; })
    params.area_correction = chain.area_correction = 0.58f;
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
  verify_authored_projection<ProjectionKind::TETRAHEDRAL, ProjectTetrahedral>();
  verify_authored_projection<ProjectionKind::CUBIC, ProjectCubic>();
  verify_authored_projection<ProjectionKind::OCTAHEDRAL, ProjectOctahedral>();
  verify_authored_projection<ProjectionKind::DODECAHEDRAL,
                             ProjectDodecahedral>();
  verify_authored_projection<ProjectionKind::ICOSAHEDRAL, ProjectIcosahedral>();
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
              Pattern::Evaluation::BAKED_CUBE, true, false>>>>;
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
  HS_EXPECT_NEAR(effect.configuration().noise.contrast, 0.7f, 1e-6f);
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
  HS_EXPECT_EQ(effect.updateParameter("modulation-source", 3),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("grid-frequency-b", 1.9f),
               ParamSetResult::APPLIED);
  HS_EXPECT_NEAR(effect.configuration().grid.frequency_b, 1.9f, 1e-6f);
  HS_EXPECT_EQ(effect.updateParameter("modulation-amount", 0.37f),
               ParamSetResult::APPLIED);
  HS_EXPECT_NEAR(effect.configuration().mechanic.amount, 0.37f, 1e-6f);
  HS_EXPECT_TRUE(
      effect.edit_document_modulation("modulation", [](auto &config) {
        config.source = Pattern::SourceKind::DISPLACEMENT;
        return ParamSetResult::APPLIED;
      }) == ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(
      effect.getParameters().find("interference-frequency-b")->readonly);
}

/** @brief Composed fixture for nested polynomial controls. */
struct PotentialControlsSpec : Pullback::Spec {
  static constexpr bool ANIMATED_PROJECTION = false;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Displace<
          Pullback::Surface::Potential<Pullback::SurfaceProvider<
              B, Pullback::PotentialParams, false, "displace">>>,
      Pullback::Stage::Project<typename Pullback::ProjectionPolicyFor<
          Pullback::ProjectionKind::STEREOGRAPHIC, B>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
          Pullback::GridSourceParams, B>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>>;
};

/** @brief Native registration and harmonic-domain admission probe. */
struct PotentialControlsProbe
    : Pullback::ComposedEffect<8, 8, PotentialControlsProbe,
                               PotentialControlsSpec,
                               Control::HookOwnership::INHERITED> {
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr std::array<PresetEntry<Params>, 1> PRESETS{};
  /** @brief Live potential controls.
   * @return Mutable polynomial parameter family. */
  auto &configuration() { return params.template get<"displace">(); }
};

/** @brief Nested term controls bind native selectors and reject invalid harmonics. */
inline void test_authored_potential_nested_controls() {
  hs_test::reset_globals();
  PotentialControlsProbe effect;
  effect.init();
  using Family = Pullback::PotentialParams;
  HS_EXPECT_EQ(Family::EXTRA_SCHEMA.size(), size_t(24));
  for (size_t index = 0; index < Family::EXTRA_SCHEMA.size(); ++index) {
    const auto &field = Family::EXTRA_SCHEMA[index];
    const auto *parameter = effect.getParameters().find(
        (std::string("displace.") + field.id).c_str());
    HS_EXPECT_TRUE(parameter != nullptr);
    if (!parameter)
      continue;
    HS_EXPECT_TRUE(parameter->target ==
                   effect.configuration().extra_address(index));
    HS_EXPECT_EQ(parameter->option_count, int(field.enum_count));
    HS_EXPECT_TRUE(parameter->animated);
  }
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-l", 3),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-m", 2),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-amplitude", 1.3f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-phase", 2.7f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.configuration().terms[0].l, uint8_t(3));
  HS_EXPECT_EQ(effect.configuration().terms[0].m, uint8_t(2));
  HS_EXPECT_EQ(effect.configuration().terms[0].amplitude, 1.3f);
  HS_EXPECT_EQ(effect.configuration().terms[0].phase, 2.7f);
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-l", 1),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.configuration().terms[0].l, uint8_t(3));
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-type", 2),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-m", 8),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-type", 0),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(effect.configuration().terms[0].type, uint8_t(2));
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-m", 2),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("displace.term-1-type", 0),
               ParamSetResult::APPLIED);
  for (size_t index = 1; index < 4; ++index) {
    const char *const AMPLITUDE[] = {"displace.term-2-amplitude",
                                     "displace.term-3-amplitude",
                                     "displace.term-4-amplitude"};
    const float VALUE = float(index) * -0.7f;
    HS_EXPECT_EQ(effect.updateParameter(AMPLITUDE[index - 1], VALUE),
                 ParamSetResult::APPLIED);
    HS_EXPECT_EQ(effect.configuration().terms[index].amplitude, VALUE);
  }
}

/** @brief Projected noise source fixture with independently varied resource spelling.
 * @tparam Key Source resource key.
 * @tparam Authored Whether to use authored node identities. */
template <Pullback::ResourceKey Key, bool Authored>
struct SourceSeedSpec : Pullback::Spec {
  static constexpr bool AUTHORED_TOPOLOGY = Authored;
  static constexpr bool ANIMATED_PROJECTION = false;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Rotate<
          Pullback::OuterCameraProvider<B, Pullback::AtlasProjectionParams>>,
      Pullback::Stage::Project<Pullback::AuthoredProjection<
          Pullback::ProjectionKind::STEREOGRAPHIC,
          Pullback::ProjectionProvider<B, Pullback::AtlasProjectionParams>>>,
      Pullback::Stage::Sample<Pullback::Source::ProjectedNoise<
          Pullback::SourceProvider<B, Pullback::ProjectedNoiseSourceParams,
                                   Key>,
          math::NoiseBasis::SIMPLEX>>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>>;
};

/** @brief Initializes source noise through a declared authored label.
 * @tparam Key Source resource key.
 * @tparam Authored Whether to use authored node identities. */
template <Pullback::ResourceKey Key, bool Authored>
struct SourceSeedProbe
    : Pullback::ComposedEffect<8, 8, SourceSeedProbe<Key, Authored>,
                               SourceSeedSpec<Key, Authored>,
                               Control::HookOwnership::INHERITED> {
  using Base = Pullback::ComposedEffect<8, 8, SourceSeedProbe<Key, Authored>,
                                        SourceSeedSpec<Key, Authored>,
                                        Control::HookOwnership::INHERITED>;
  using Params = typename Base::Params;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr std::array<PresetEntry<Params>, 1> PRESETS{};
  /** @brief Authored node label corresponding to a source resource.
   * @param resource Semantic resource key.
   * @return Declared source label or the unchanged resource. */
  static constexpr std::string_view
  authored_instance(std::string_view resource) {
    return resource == Key.view() ? "authored_source" : resource;
  }
};

template <Pullback::ResourceKey Key, bool Authored>
inline void verify_source_seed(uint32_t expected_seed) {
  hs_test::reset_globals();
  SourceSeedProbe<Key, Authored> effect;
  effect.init();
  const auto FRAME = effect.frame_for_test();
  const auto &ACTUAL = *FRAME.resources.template get<Key>().noise;
  FastNoiseLite expected;
  Pullback::init_effect_noise(expected, static_cast<int32_t>(expected_seed));
  for (const auto &point :
       {math::Vector(0.13f, -0.77f, 0.91f), math::Vector(-2.3f, 1.17f, 0.42f)})
    HS_EXPECT_NEAR(ACTUAL.GetNoiseSingle(point.x, point.y, point.z),
                   expected.GetNoiseSingle(point.x, point.y, point.z), 1e-6f);
}

/** @brief Projection frame seeded from an authored projection node. */
struct ProjectionSeedSpec : Pullback::Spec {
  static constexpr bool AUTHORED_TOPOLOGY = true;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Rotate<
          Pullback::OuterCameraProvider<B, Pullback::AtlasProjectionParams>>,
      Pullback::Stage::Project<Pullback::AuthoredProjection<
          Pullback::ProjectionKind::STEREOGRAPHIC,
          Pullback::ProjectionProvider<B, Pullback::AtlasProjectionParams>>>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
          Pullback::GridSourceParams, B>::Type>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>>;
};
/** @brief Projection walk seed binding probe. */
struct ProjectionSeedProbe
    : Pullback::ComposedEffect<8, 8, ProjectionSeedProbe, ProjectionSeedSpec,
                               Control::HookOwnership::INHERITED> {
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr std::array<PresetEntry<Params>, 1> PRESETS{};
  /** @brief Authored label of the projection node.
   * @param resource Semantic resource key.
   * @return Declared projection label or the unchanged resource. */
  static constexpr std::string_view
  authored_instance(std::string_view resource) {
    return resource == "projection" ? "authored_projection" : resource;
  }
};

/** @brief Authored noise seeds depend on node identity while fixed seeds stay stable. */
inline void test_authored_noise_seed_identity() {
  static_assert(
      Pullback::instance_hash("authored_source", "sample.projected-noise.v2") ==
      0xC4BCAA76u);
  static_assert(Pullback::instance_hash("authored_projection",
                                        "project.stereographic.v2") ==
                0x2058155Cu);
  verify_source_seed<"source", true>(0xC4BCAA76u);
  verify_source_seed<"noisy", true>(0xC4BCAA76u);
  verify_source_seed<"source", false>(1337u);
  verify_source_seed<"noisy", false>(0x9C442C8Cu);
  hs_test::reset_globals();
  ProjectionSeedProbe effect;
  effect.init();
  const auto &ACTUAL = hs_test::ComposedFrameWhiteBox::projection_noise(effect);
  FastNoiseLite expected;
  Pullback::init_effect_noise(expected, static_cast<int32_t>(0x2058155Cu));
  expected.SetFrequency(Animation::RandomWalk<8>::Options{}.noise_scale);
  for (const auto &point :
       {math::Vector(0.13f, -0.77f, 0.91f), math::Vector(-2.3f, 1.17f, 0.42f)})
    HS_EXPECT_NEAR(ACTUAL.GetNoiseSingle(point.x, point.y, point.z),
                   expected.GetNoiseSingle(point.x, point.y, point.z), 1e-6f);
}
