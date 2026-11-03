/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Fixed-step field admission, clocks and advisories. */

/** @brief Gated integer modes require an ungated selector and a partial mask. */
inline void test_tangent_integer_gates() {
  constexpr auto VALID = In::SCHEMA<In::Op::DisplacePotential>;
  HS_EXPECT_TRUE(In::gates_wellformed(VALID));
  constexpr auto UNKNOWN = [] {
    auto schema = In::SCHEMA<In::Op::DisplacePotential>;
    for (auto &field : schema)
      if (std::string_view(field.id) == "term-1-l")
        field.gated_by = "missing";
    return schema;
  }();
  HS_EXPECT_FALSE(In::gates_wellformed(UNKNOWN));
  constexpr auto CYCLE = [] {
    auto schema = In::SCHEMA<In::Op::DisplacePotential>;
    for (auto &field : schema)
      if (std::string_view(field.id) == "term-1-type") {
        field.gated_by = "term-1-l";
        field.gate_values = 1;
      }
    return schema;
  }();
  HS_EXPECT_FALSE(In::gates_wellformed(CYCLE));
  constexpr auto ALL_VALUES = [] {
    auto schema = In::SCHEMA<In::Op::DisplacePotential>;
    for (auto &field : schema)
      if (std::string_view(field.id) == "term-1-l")
        field.gate_values = 15;
    return schema;
  }();
  HS_EXPECT_FALSE(In::gates_wellformed(ALL_VALUES));
}

/** @brief Engine mode controls expose bounded integer dropdowns and their gates. */
inline void test_tangent_integer_controls() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"field", In::Op::DisplacePotential::ID},
      {"wave", In::Op::DisplaceRossbyHaurwitz::ID},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"}};
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  for (const char *id : {"field.term-1-l", "field.term-1-m"}) {
    const auto *parameter = require_param(effect, id);
    if (!parameter)
      continue;
    HS_EXPECT_TRUE(parameter->is_enum());
    HS_EXPECT_EQ(parameter->option_count, 9);
    HS_EXPECT_EQ(std::string_view(parameter->options[0]), "zero");
    HS_EXPECT_EQ(std::string_view(parameter->options[8]), "eight");
  }
  const auto *order = require_param(effect, "wave.order");
  if (order) {
    HS_EXPECT_EQ(order->option_count, 7);
    HS_EXPECT_EQ(std::string_view(order->options[0]), "one");
    HS_EXPECT_EQ(std::string_view(order->options[6]), "seven");
    HS_EXPECT_EQ(order->get_requested(), 3.0f);
  }
  const auto *op = require_operator(In::Op::DisplacePotential::ID);
  for (const auto &field : op->schema_span())
    if (std::string_view(field.id) == "term-1-l" ||
        std::string_view(field.id) == "term-1-m") {
      HS_EXPECT_TRUE(field.topology);
      HS_EXPECT_EQ(std::string_view(field.gated_by), "term-1-type");
      HS_EXPECT_EQ(field.gate_values, 1);
    }
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-l", 3.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-m", 4.0f),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_EQ(requested_param(effect, "field.term-1-m"), 0.0f);
  HS_EXPECT_EQ(effect.updateParameter("wave.order", 6.0f),
               ParamSetResult::APPLIED);
  const auto &program = ShaderChainWhiteBox::program(effect);
  const auto &params = *reinterpret_cast<const In::Op::RossbyHaurwitzParams *>(
      program.param_block(1));
  const auto field = In::Op::DisplaceRossbyHaurwitz::field(params, 0.0f);
  HS_EXPECT_EQ(field.terms[1].term.l, 8);
  HS_EXPECT_EQ(field.terms[1].term.m, 7);
}

struct RotationTestField {
  math::Vector velocity(const math::Vector &p) const {
    return math::cross(math::Y_AXIS, p);
  }
  float speed_bound() const { return 1.0f; }
  float lipschitz_bound() const { return 1.0f; }
};

/** @brief Symmetric sectoral fields use the positive real-power invariant. */
inline void test_symmetric_sectoral_phase() {
  constexpr uint8_t SYMMETRIES[] = {0, 4, 5, 6, 7, 8};
  constexpr int ORDERS[] = {6, 3, 4, 5, 6, 8};
  for (size_t i = 0; i < std::size(SYMMETRIES); ++i) {
    const float longitude = math::PI_F / (2 * ORDERS[i]);
    const math::Vector p(cosf(longitude), 0, sinf(longitude));
    const math::Vector azimuth(-p.z, 0, p.x);
    In::Op::SymmetricDisplaceParams params;
    params.symmetry = SYMMETRIES[i];
    for (uint8_t kind = 0; kind < 2; ++kind) {
      params.kind = kind;
      const auto field = In::Op::DisplaceSymmetric::field(params, 0);
      const math::Vector expected =
          (kind == 0 ? math::Y_AXIS : -azimuth) / 1.01f;
      const auto velocity = field.velocity(p);
      HS_EXPECT_NEAR(velocity.x, expected.x, 2e-5f);
      HS_EXPECT_NEAR(velocity.y, expected.y, 2e-5f);
      HS_EXPECT_NEAR(velocity.z, expected.z, 2e-5f);
    }
  }
}

/** @brief Symmetric lenses expose one gain spanning the legacy product domain. */
inline void test_symmetric_single_strength() {
  using Symmetric = In::Op::DisplaceSymmetric;
  HS_EXPECT_STREQ(Symmetric::ID, "sphere.displace.symmetric.v2");
  const auto *op = require_operator(Symmetric::ID);
  if (!op)
    return;
  HS_EXPECT_EQ(op->schema_count, 5);
  const In::ParamFieldInfo *strength = nullptr;
  for (const auto &field : op->schema_span()) {
    HS_EXPECT_FALSE(std::string_view(field.id) == "amplitude");
    if (std::string_view(field.id) == "strength")
      strength = &field;
  }
  HS_EXPECT_TRUE(strength != nullptr);
  if (!strength)
    return;
  HS_EXPECT_EQ(strength->min, -4.0f * math::PI_F);
  HS_EXPECT_EQ(strength->max, 4.0f * math::PI_F);
  HS_EXPECT_EQ(strength->def, 0.05f);
  for (const float value : {-4.0f * math::PI_F, 0.0f, 4.0f * math::PI_F}) {
    Symmetric::Params params;
    params.strength = value;
    HS_EXPECT_TRUE(Symmetric::validate(params) == nullptr);
    const auto field = Symmetric::field(params, 0.0f);
    HS_EXPECT_EQ(field.speed_bound(), 1.0f);
  }
  Symmetric::Params outside;
  outside.strength = std::nextafter(4.0f * math::PI_F, INFINITY);
  HS_EXPECT_TRUE(Symmetric::validate(outside) != nullptr);
  reset_globals();
  ShaderChainWhiteBox::FX effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"symmetric", Symmetric::ID},
      {"project", In::Op::ProjectStereographic::ID},
      {"sample", In::Op::SampleGridV3::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  HS_EXPECT_TRUE(effect.getParameters().find("symmetric.amplitude") == nullptr);
  HS_EXPECT_EQ(effect.updateParameter("symmetric.strength", 4.0f * math::PI_F),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(requested_param(effect, "symmetric.strength"),
               4.0f * math::PI_F);
  const In::ChainEntryRequest retired[] = {
      {"symmetric", "sphere.displace.symmetric.v1"},
      {"project", In::Op::ProjectStereographic::ID},
      {"sample", In::Op::SampleGridV3::ID},
      {"colorize", In::Op::ColorizeGeneratedPaletteV4::ID}};
  HS_EXPECT_EQ(effect.set_chain(retired).code,
               In::ChainStatus::UNKNOWN_OPERATOR);
}

/** @brief One Strength reproduces the legacy gain product for every invariant. */
inline void test_symmetric_legacy_gain_equivalence() {
  using Symmetric = In::Op::DisplaceSymmetric;
  using namespace math::tangent_field;
  struct Invariant {
    TermType type;
    uint8_t degree;
    float sign;
  };
  constexpr Invariant INVARIANTS[] = {
      {TermType::MODE, 6, 1},       {TermType::TETRAHEDRAL, 0, 1},
      {TermType::OCTAHEDRAL, 0, 1}, {TermType::DODECAHEDRAL, 0, 1},
      {TermType::MODE, 3, -1},      {TermType::MODE, 4, 1},
      {TermType::MODE, 5, -1},      {TermType::MODE, 6, 1},
      {TermType::MODE, 8, 1},
  };
  for (uint8_t symmetry = 0; symmetry < std::size(INVARIANTS); ++symmetry)
    for (uint8_t kind = 0; kind < 2; ++kind)
      for (uint8_t integrator = 0; integrator < 3; ++integrator)
        for (const float old_strength : {-math::PI_F, 0.25f, math::PI_F})
          for (const float old_amplitude : {-4.0f, -1.25f, 0.0f, 0.5f, 4.0f})
            for (const float phase : {-0.7f, 0.6f}) {
              const auto &invariant = INVARIANTS[symmetry];
              std::array<PotentialTerm, 4> terms{};
              terms[0] = {invariant.type,
                          static_cast<Kind>(kind),
                          invariant.degree,
                          invariant.degree,
                          invariant.sign * old_amplitude,
                          0};
              const PotentialField legacy(
                  terms, 0.0f, math::make_rotation(math::Y_AXIS, phase));
              Symmetric::Params params;
              params.symmetry = symmetry;
              params.kind = kind;
              params.integrator = integrator;
              params.strength = old_strength * old_amplitude;
              const auto field = Symmetric::field(params, phase);
              for (int index = 0; index < 3; ++index) {
                const auto p = math::Vector::from_spherical(
                    0.31f + index * 1.1f, 0.61f + index * 0.5f);
                const auto expected =
                    Symmetric::step(p, legacy, old_strength, integrator);
                const auto actual =
                    Symmetric::run({p, 0.2f}, {}, params, field);
                HS_EXPECT_NEAR(actual.dir.x, expected.sphere.x, 1e-3f);
                HS_EXPECT_NEAR(actual.dir.y, expected.sphere.y, 1e-3f);
                HS_EXPECT_NEAR(actual.dir.z, expected.sphere.z, 1e-3f);
                HS_EXPECT_NEAR(actual.path_length, 0.2f + expected.path_length,
                               1e-3f);
              }
            }
}

/** @brief Every invariant rotates its field and displacement bodily about y. */
inline void test_symmetric_frame_rotation() {
  using Symmetric = In::Op::DisplaceSymmetric;
  const auto rotate_y = [](const math::Vector &p, double angle) {
    return math::Vector(cos(angle) * p.x + sin(angle) * p.z, p.y,
                        -sin(angle) * p.x + cos(angle) * p.z);
  };
  const auto expect_vector = [](const math::Vector &actual,
                                const math::Vector &expected, float tolerance) {
    HS_EXPECT_NEAR(actual.x, expected.x, tolerance);
    HS_EXPECT_NEAR(actual.y, expected.y, tolerance);
    HS_EXPECT_NEAR(actual.z, expected.z, tolerance);
  };
  for (uint8_t symmetry = 0; symmetry < 9; ++symmetry)
    for (uint8_t kind = 0; kind < 2; ++kind) {
      Symmetric::Params params;
      params.symmetry = symmetry;
      params.kind = kind;
      params.strength = 0.2f;
      const auto fixed = Symmetric::field(params, 0.0f);
      float movement = 0.0f;
      for (const float phase : {-1.13f, 0.37f}) {
        const auto rotating = Symmetric::prepare({}, params, {phase});
        HS_EXPECT_EQ(rotating.speed_bound(), fixed.speed_bound());
        HS_EXPECT_EQ(rotating.lipschitz_bound(), fixed.lipschitz_bound());
        for (int sample = 0; sample < 8; ++sample) {
          const auto p = math::Vector::from_spherical(0.29f + sample * 0.59f,
                                                      0.41f + sample * 0.23f);
          const auto local = rotate_y(p, -phase);
          const auto expected = rotate_y(fixed.velocity(local), phase);
          const auto actual = rotating.velocity(p);
          expect_vector(actual, expected, 2e-5f);
          movement += (actual - fixed.velocity(p)).length();
          for (uint8_t integrator = 0; integrator < 3; ++integrator) {
            params.integrator = integrator;
            const auto stepped = Symmetric::step(local, fixed, params.strength,
                                                 params.integrator);
            const auto output =
                Symmetric::run({p, 0.25f}, {}, params, rotating);
            expect_vector(output.dir, rotate_y(stepped.sphere, phase), 3e-5f);
            HS_EXPECT_NEAR(output.path_length, 0.25f + stepped.path_length,
                           1e-5f);
          }
        }
      }
      HS_EXPECT_GT(movement, 0.01f);
    }
}

/** @brief Speed holds, reverses and continuously samples the instance rotation. */
inline void test_symmetric_rotation_clock() {
  using Symmetric = In::Op::DisplaceSymmetric;
  const auto *op = require_operator(Symmetric::ID);
  if (!op)
    return;
  HS_EXPECT_TRUE(op->runtime.sample_state != nullptr);
  HS_EXPECT_TRUE(op->runtime.prepare_interval != nullptr);
  const In::ParamFieldInfo *speed = nullptr;
  for (const auto &field : op->schema_span())
    if (std::string_view(field.id) == "speed")
      speed = &field;
  HS_EXPECT_TRUE(speed != nullptr);
  if (!speed)
    return;
  HS_EXPECT_STREQ(speed->display_name, "Speed");
  HS_EXPECT_EQ(speed->def, 0.0f);
  HS_EXPECT_EQ(speed->min, -0.1f);
  HS_EXPECT_EQ(speed->max, 0.1f);
  Symmetric::Params params;
  Symmetric::State state{0.7f};
  Symmetric::advance(state, params);
  HS_EXPECT_EQ(state.phase, 0.7f);
  params.speed = 0.04f;
  Symmetric::advance(state, params);
  HS_EXPECT_NEAR(state.phase, 0.74f, 1e-7f);
  params.speed = -0.01f;
  Symmetric::advance(state, params);
  HS_EXPECT_NEAR(state.phase, 0.73f, 1e-7f);
  for (const float rate : {-0.04f, 0.04f}) {
    params.speed = rate;
    const Symmetric::State from{rate > 0 ? math::TWO_PI_F - 0.02f
                                         : -math::TWO_PI_F + 0.02f};
    auto to = from;
    Symmetric::advance(to, params);
    HS_EXPECT_LT(fabsf(to.phase), 0.03f);
    const auto interval = In::Op::temporal_interval(from, to);
    for (const float fraction : {0.0f, 0.25f, 0.5f, 1.0f}) {
      auto column = to;
      In::Op::temporal_sample(column, interval, params, fraction);
      const double expected = double(from.phase) + rate * fraction;
      HS_EXPECT_NEAR(cos(column.phase), cos(expected), 1e-6);
      HS_EXPECT_NEAR(sin(column.phase), sin(expected), 1e-6);
    }
    const auto snapshot = In::RuntimeStateCodec<Symmetric::State>::capture(to);
    Symmetric::State restored;
    HS_EXPECT_TRUE(
        In::RuntimeStateCodec<Symmetric::State>::restore(restored, snapshot));
    HS_EXPECT_EQ(restored.phase, to.phase);
    Symmetric::advance(restored, params);
    Symmetric::advance(to, params);
    HS_EXPECT_EQ(restored.phase, to.phase);
  }
}

/** @brief Legacy snapshots restore zero rotation; new snapshots require clocks. */
inline void test_symmetric_legacy_snapshot() {
  using Symmetric = In::Op::DisplaceSymmetric;
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  const In::ChainEntryRequest chain[] = {
      {"symmetric", Symmetric::ID},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"}};
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  HS_EXPECT_EQ(effect.updateParameter("symmetric.speed", -0.04f),
               ParamSetResult::APPLIED);
  for (int frame = 0; frame < 7; ++frame)
    WB::advance_without_render(effect);
  const auto saved = effect.snapshot();
  const auto phase_of = [](const ChainSnapshot &snapshot) {
    for (const auto &entry : *snapshot.runtime)
      if (entry.instance == "symmetric")
        return std::get<In::SourceClockSnapshot>(entry.state).primary;
    return std::numeric_limits<float>::quiet_NaN();
  };
  HS_EXPECT_NEAR(phase_of(saved), -0.28f, 1e-6f);
  auto missing = saved;
  std::erase_if(*missing.runtime, [](const auto &entry) {
    return entry.instance == "symmetric";
  });
  const auto generation = effect.getParameterSchemaGeneration();
  HS_EXPECT_EQ(effect.restore_snapshot(missing),
               ChainSnapshotRestoreResult::INVALID_VALUE);
  HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), generation);
  HS_EXPECT_EQ(phase_of(effect.snapshot()), phase_of(saved));
  HS_EXPECT_EQ(requested_param(effect, "symmetric.speed"), -0.04f);
  auto legacy = missing;
  std::erase_if(legacy.parameters, [](const auto &parameter) {
    return parameter.name == "symmetric.speed";
  });
  HS_EXPECT_EQ(effect.restore_snapshot(legacy),
               ChainSnapshotRestoreResult::APPLIED);
  HS_EXPECT_EQ(phase_of(effect.snapshot()), 0.0f);
  HS_EXPECT_EQ(requested_param(effect, "symmetric.speed"), 0.0f);
  HS_EXPECT_EQ(effect.restore_snapshot(saved),
               ChainSnapshotRestoreResult::APPLIED);
  HS_EXPECT_EQ(phase_of(effect.snapshot()), phase_of(saved));
}

/** @brief Every field integrator preserves equatorial rotation and arc length. */
inline void test_tangent_step_integrators() {
  const RotationTestField field;
  const math::Vector input(1, 0, 0);
  for (const float strength :
       {-16.0f * math::PI_F, -4.0f * math::PI_F, -2.0f * math::PI_F,
        -math::PI_F, -0.5001f, -0.5f, 0.0f, 0.5f, 0.5001f, math::PI_F,
        2.0f * math::PI_F, 4.0f * math::PI_F, 16.0f * math::PI_F}) {
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

struct MeridianTestField {
  math::Vector velocity(const math::Vector &p) const {
    return math::Y_AXIS - p.y * p +
           (0.5f * (1.0f - p.x)) * (math::Z_AXIS - p.z * p);
  }
  float speed_bound() const { return 2.0f; }
  float lipschitz_bound() const { return 4.0f; }
};

/** @brief Midpoint transport follows its known geodesic through antipodes. */
inline void test_tangent_antipodal_midpoint() {
  const MeridianTestField field;
  for (const float strength :
       {-8.0f * math::PI_F, -2.0f * math::PI_F, 2.0f * math::PI_F - 0.001f,
        2.0f * math::PI_F, 2.0f * math::PI_F + 0.001f, 6.0f * math::PI_F}) {
    const double cosine = cos(0.5 * double(strength));
    const double z = 0.5 * (1.0 - cosine);
    const double arc = fabs(double(strength)) * sqrt(cosine * cosine + z * z);
    const double sinc = sin(arc) / arc;
    const auto result =
        PB::Surface::step_field<MeridianTestField, PB::Surface::Midpoint>(
            math::X_AXIS, field, strength, true);
    HS_EXPECT_NEAR(result.sphere.x, cos(arc), 5e-5f);
    HS_EXPECT_NEAR(result.sphere.y, sinc * strength * cosine, 5e-5f);
    HS_EXPECT_NEAR(result.sphere.z, sinc * strength * z, 5e-5f);
    HS_EXPECT_NEAR(result.path_length, arc, 5e-5f);
    HS_EXPECT_NEAR(result.sphere.length(), 1.0f, 5e-5f);
  }
  const RotationTestField rotation;
  for (const auto result :
       {PB::Surface::step_field<RotationTestField, PB::Surface::Midpoint>(
            math::Y_AXIS, rotation, 8 * math::PI_F, true),
        PB::Surface::step_field<RotationTestField, PB::Surface::Midpoint2>(
            math::Y_AXIS, rotation, 8 * math::PI_F, true)}) {
    HS_EXPECT_EQ(result.sphere.x, 0.0f);
    HS_EXPECT_EQ(result.sphere.y, 1.0f);
    HS_EXPECT_EQ(result.sphere.z, 0.0f);
    HS_EXPECT_EQ(result.path_length, 0.0f);
  }
}

/** @brief Scalar and harmonic domains reject inadmissible recipes. */
inline void test_tangent_operator_admission() {
  In::Op::PotentialDisplaceParams p;
  p.terms[0].amplitude = 1;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) == nullptr);
  p.terms[0].m = 2;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  p.terms[0].m = 0;
  p.terms[0].l = 9;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  p.terms[0].l = 1;
  p.terms[0].amplitude = 2;
  p.strength = math::PI_F;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) == nullptr);
  p.strength = -math::PI_F / 2;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) == nullptr);
  p.terms[1].amplitude = NAN;
  HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(p) != nullptr);
  In::Op::RossbyHaurwitzParams r;
  for (uint8_t order : {7, 8, 255}) {
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
        params.strength = 4.0f * math::PI_F;
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
          HS_EXPECT_LE(out.path_length, 4.0f * math::PI_F + 1e-4f);
        }
      }
}

/** @brief Mixed four-term recipes stay finite at the scalar range endpoints. */
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
      params.strength = sign * math::PI_F;
      params.integrator = scheme;
      HS_EXPECT_TRUE(In::Op::DisplacePotential::validate(params) == nullptr);
      const auto field = In::Op::DisplacePotential::field(params, 0.9f);
      for (int i = 0; i < 24; ++i) {
        const auto p = math::Vector::from_spherical(0.817f * i, 0.137f * i);
        const auto result =
            In::Op::DisplacePotential::step(p, field, params.strength, scheme);
        HS_EXPECT_NEAR(result.sphere.length(), 1, 5e-4f);
        HS_EXPECT_TRUE(std::isfinite(result.path_length));
        HS_EXPECT_LE(result.path_length, 16.0f * math::PI_F + 1e-4f);
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
      {"field", In::Op::DisplacePotential::ID},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v4"}};
  HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  const std::array replay{
      std::pair<std::string, float>{"field.strength", 0.1f},
      std::pair<std::string, float>{"field.term-1-amplitude", 1.0f},
      std::pair<std::string, float>{"field.term-1-l", 6.0f},
      std::pair<std::string, float>{"field.term-1-m", 6.0f}};
  effect.replay_parameter_writes(replay);
  HS_EXPECT_EQ(requested_param(effect, "field.strength"), 0.1f);
  HS_EXPECT_EQ(requested_param(effect, "field.term-1-amplitude"), 1.0f);
  float first = 0;
  HS_EXPECT_TRUE(effect.parameter_fold_strength("field.strength", first));
  HS_EXPECT_GT(first, 0.0f);
  HS_EXPECT_LT(first, 1.0f);
  HS_EXPECT_EQ(effect.updateParameter("field.term-1-amplitude", 2.0f),
               ParamSetResult::APPLIED);
  float second = 0;
  HS_EXPECT_TRUE(effect.parameter_fold_strength("field.strength", second));
  HS_EXPECT_NEAR(second, first / 2, 1e-6f);
  const float beyond = 2 * second;
  HS_EXPECT_EQ(effect.updateParameter("field.strength", beyond),
               ParamSetResult::APPLIED);
  HS_EXPECT_NEAR(requested_param(effect, "field.strength"), beyond, 1e-6f);
  HS_EXPECT_TRUE(effect.parameter_warning("field.strength") == nullptr);
  HS_EXPECT_FALSE(
      effect.parameter_fold_strength("field.term-1-amplitude", second));
  HS_EXPECT_EQ(effect.updateParameter("field.strength", math::PI_F),
               ParamSetResult::APPLIED);
  HS_EXPECT_NEAR(requested_param(effect, "field.strength"), math::PI_F, 1e-6f);
}
