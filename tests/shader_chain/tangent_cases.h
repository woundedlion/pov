/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Fixed-step field admission, clocks and advisories. */

struct RotationTestField {
  math::Vector velocity(const math::Vector &p) const {
    return math::cross(math::Y_AXIS, p);
  }
  float speed_bound() const { return 1.0f; }
  float lipschitz_bound() const { return 1.0f; }
};

/** @brief Every field integrator preserves equatorial rotation and arc length. */
inline void test_tangent_step_integrators() {
  const RotationTestField field;
  const math::Vector input(1, 0, 0);
  for (const float strength :
       {-math::PI_F, -0.5001f, -0.5f, 0.0f, 0.5f, 0.5001f, math::PI_F}) {
    const auto a =
        PB::Surface::step_field<RotationTestField, PB::Surface::Euler>(
            input, field, strength, true);
    const auto b =
        PB::Surface::step_field<RotationTestField, PB::Surface::Midpoint>(
            input, field, strength, true);
    const auto c =
        PB::Surface::step_field<RotationTestField, PB::Surface::Midpoint2>(
            input, field, strength, true);
    const math::Vector expected(cosf(strength), 0, -sinf(strength));
    for (const auto &result : {a, b, c}) {
      HS_EXPECT_NEAR(result.sphere.x, expected.x, 2e-5f);
      HS_EXPECT_NEAR(result.sphere.y, expected.y, 2e-5f);
      HS_EXPECT_NEAR(result.sphere.z, expected.z, 2e-5f);
      HS_EXPECT_NEAR(result.path_length, fabsf(strength), 1e-5f);
      HS_EXPECT_NEAR(result.sphere.length(), 1.0f, 2e-5f);
    }
    const auto no_path =
        PB::Surface::step_field<RotationTestField, PB::Surface::Midpoint2>(
            input, field, strength, false);
    HS_EXPECT_EQ(no_path.path_length, 0.0f);
  }
}

/** @brief Mode domains and the uniform arc cap reject inadmissible recipes. */
inline void test_tangent_operator_admission() {
  In::Op::PotentialDisplaceParams p;
  p.terms[0].amplitude = 1;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) == nullptr);
  p.terms[0].m = 2;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  p.terms[0].m = 0;
  p.terms[0].l = 9;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  p.terms[0].l = 1.5f;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  p.terms[0].l = 1;
  p.terms[0].amplitude = 2;
  p.strength = math::PI_F;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  p.strength = -math::PI_F / 2;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) == nullptr);
  p.terms[1].amplitude = NAN;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  In::Op::RossbyHaurwitzParams r;
  for (float order : {0.0f, 1.5f, 8.0f, NAN}) {
    r.order = order;
    HS_EXPECT_TRUE(In::Op::DisplaceRossbyHaurwitz::validate(r) != nullptr);
  }
  In::Op::SymmetricDisplaceParams s;
  s.symmetry = 9;
  HS_EXPECT_TRUE(In::Op::DisplaceSymmetric::validate(s) != nullptr);
}

/** @brief Admitted field boundaries stay finite under every stepping scheme. */
inline void test_tangent_admitted_domain() {
  for (uint8_t symmetry = 0; symmetry < 9; ++symmetry)
    for (uint8_t kind = 0; kind < 2; ++kind)
      for (uint8_t integrator = 0; integrator < 3; ++integrator) {
        In::Op::SymmetricDisplaceParams params;
        params.symmetry = symmetry;
        params.kind = kind;
        params.integrator = integrator;
        params.strength = math::PI_F;
        HS_EXPECT_TRUE(In::Op::DisplaceSymmetric::validate(params) == nullptr);
        const auto field = In::Op::DisplaceSymmetric::field(params, 0);
        for (int index = 0; index < 24; ++index) {
          const auto p =
              math::Vector::from_spherical(index * 0.817f, index * 0.137f);
          const auto out = In::Op::DisplaceSymmetric::step(
              p, field, params.strength, integrator);
          HS_EXPECT_TRUE(std::isfinite(out.sphere.x) &&
                         std::isfinite(out.sphere.y) &&
                         std::isfinite(out.sphere.z));
          HS_EXPECT_NEAR(out.sphere.length(), 1.0f, 5e-4f);
          HS_EXPECT_LE(out.path_length, math::PI_F + 1e-4f);
        }
      }
}

/** @brief Mixed four-term recipes respect the admitted cap under each integrator. */
inline void test_tangent_mixed_boundaries() {
  In::Op::PotentialDisplaceParams params;
  params.orientation_azimuth = 0.7f;
  params.orientation_polar = 1.1f;
  for (size_t i = 0; i < params.terms.size(); ++i) {
    auto &term = params.terms[i];
    term.type = static_cast<uint8_t>(i);
    term.kind = static_cast<uint8_t>(i % 2);
    term.l = 8;
    term.m = 7;
    term.amplitude = i % 2 ? -4.0f : 4.0f;
  }
  for (const float sign : {-1.0f, 1.0f})
    for (uint8_t scheme = 0; scheme < 3; ++scheme) {
      params.strength = sign * math::PI_F / 16;
      params.integrator = scheme;
      HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(params) == nullptr);
      const auto field = In::Op::DisplacePotential::field(params, 0.9f);
      for (int i = 0; i < 24; ++i) {
        const auto p = math::Vector::from_spherical(0.817f * i, 0.137f * i);
        const auto result =
            In::Op::DisplacePotential::step(p, field, params.strength, scheme);
        HS_EXPECT_NEAR(result.sphere.length(), 1, 5e-4f);
        HS_EXPECT_LE(result.path_length, math::PI_F + 1e-4f);
      }
    }
}

/** @brief Field drift accumulates continuously and survives snapshot restore. */
inline void test_tangent_phase_clock() {
  In::Op::FieldPhaseState state;
  In::Op::RossbyHaurwitzParams params;
  params.drift = 0.04f;
  In::Op::DisplaceRossbyHaurwitz::advance(state, params);
  const auto old = state;
  params.drift = 0.01f;
  In::Op::DisplaceRossbyHaurwitz::advance(state, params);
  HS_EXPECT_NEAR(state.phase, 0.05f, 1e-7f);
  auto column = state;
  In::Op::temporal_sample(column, In::Op::temporal_interval(old, state), params,
                          0.5f);
  HS_EXPECT_NEAR(column.phase, 0.045f, 1e-7f);
  const auto saved =
      In::RuntimeStateCodec<In::Op::FieldPhaseState>::capture(state);
  state.phase = 0;
  HS_EXPECT_TRUE(
      In::RuntimeStateCodec<In::Op::FieldPhaseState>::restore(state, saved));
  HS_EXPECT_NEAR(state.phase, 0.05f, 1e-7f);
}

/** @brief The advisory follows the mix and permits strengths beyond its threshold. */
inline void test_tangent_dynamic_advisory() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"field", In::Op::DisplaceSymmetric::ID},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"}};
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  float first = 0;
  HS_EXPECT_TRUE(effect.parameter_fold_strength("field.strength", first));
  HS_EXPECT_GT(first, 0.0f);
  HS_EXPECT_LT(first, 1.0f);
  HS_EXPECT_EQ(effect.updateParameter("field.amplitude", 2.0f),
               ParamSetResult::APPLIED);
  float second = 0;
  HS_EXPECT_TRUE(effect.parameter_fold_strength("field.strength", second));
  HS_EXPECT_NEAR(second, first / 2, 1e-6f);
  const float beyond = 2 * second;
  HS_EXPECT_EQ(effect.updateParameter("field.strength", beyond),
               ParamSetResult::APPLIED);
  HS_EXPECT_NEAR(requested_param(effect, "field.strength"), beyond, 1e-6f);
  HS_EXPECT_TRUE(effect.parameter_warning("field.strength") == nullptr);
  HS_EXPECT_FALSE(effect.parameter_fold_strength("field.amplitude", second));
  HS_EXPECT_EQ(effect.updateParameter("field.strength", math::PI_F),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_NEAR(requested_param(effect, "field.strength"), beyond, 1e-6f);
}
