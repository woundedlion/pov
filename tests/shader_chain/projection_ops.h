/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- projection-op mirrors ------------------------------------------------

/** Mirror frame for the projection batch. */
struct ProjMirrorFrame {
  math::Quaternion conjugate;
  In::Op::MeridianProjectChainParams meridian;
  In::Op::GnomonicChainParams gnomonic;
  In::Op::ProjectBonneV3::Params bonne;
  In::Op::ProjectPeirceV3::Params peirce;
  In::Op::ProjectPeirceSquareFastV3::Params peirce_fast;
  In::Op::ProjectAiroceanV3::Params airocean;
  In::Op::PolyhedralAtlasChainParams atlas;
};

struct ProjMirrorBinding {
  using FrameState = ProjMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct ProjMirrorBase {
  using Binding = ProjMirrorBinding;
  using FrameState = ProjMirrorFrame;
  static const math::Quaternion &conjugate(const ProjMirrorFrame &frame) {
    return frame.conjugate;
  }
};

struct MeridianProjMirror : ProjMirrorBase {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.meridian.central_meridian;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.meridian.singularity_fade;
  }
};

struct GnomonicProjMirror : ProjMirrorBase {
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.gnomonic.singularity_fade;
  }
};

struct PeirceProjMirror : ProjMirrorBase {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.peirce.central_meridian;
  }
  static float layout_scroll(const ProjMirrorFrame &frame) {
    return frame.peirce.layout_scroll;
  }
  static float coordinate_scale(const ProjMirrorFrame &frame) {
    return frame.peirce.coordinate_scale;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.peirce.singularity_fade;
  }
};

struct PeirceFastProjMirror : ProjMirrorBase {
  static constexpr bool ZERO_CENTRAL_MERIDIAN = true;
  static float coordinate_scale(const ProjMirrorFrame &frame) {
    return frame.peirce_fast.coordinate_scale;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.peirce_fast.singularity_fade;
  }
};

struct PeirceFastMovingMeridianMirror : PeirceFastProjMirror {
  static constexpr bool ZERO_CENTRAL_MERIDIAN = false;
};

struct BonneProjMirror : ProjMirrorBase {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.bonne.central_meridian;
  }
  static float standard_parallel(const ProjMirrorFrame &frame) {
    return frame.bonne.standard_parallel;
  }
  static float coordinate_scale(const ProjMirrorFrame &frame) {
    return frame.bonne.coordinate_scale;
  }
};

struct AiroceanProjMirror : ProjMirrorBase {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.airocean.central_meridian;
  }
  static float coordinate_scale(const ProjMirrorFrame &frame) {
    return frame.airocean.coordinate_scale;
  }
};

struct AtlasProjMirror : ProjMirrorBase {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.atlas.central_meridian;
  }
  static float coordinate_scale(const ProjMirrorFrame &frame) {
    return frame.atlas.coordinate_scale;
  }
  static float area_correction(const ProjMirrorFrame &frame) {
    return frame.atlas.area_correction;
  }
};

inline bool is_atlas_operator(std::string_view id) {
  return id == In::Op::ProjectTetrahedral::ID ||
         id == In::Op::ProjectCubic::ID ||
         id == In::Op::ProjectOctahedral::ID ||
         id == In::Op::ProjectDodecahedral::ID ||
         id == In::Op::ProjectIcosahedral::ID;
}

/** Compiles a chain with one projection at entry 1 and applies the value set
    to every block. */
template <typename OpParams>
inline bool arm_project_op_chain(In::ChainProgram &program, const char *op_id,
                                 int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", op_id},
      {"sample", "sample.grid.v4"},
      {"colorize", "colorize.generated-palette.v4"},
  };
  return arm_chain<In::Op::RotateChainParams, OpParams,
                   In::Op::GridSampleParams, In::Op::GeneratedPaletteParams>(
      program, chain, frames, set);
}

/** @brief Source edge fading reaches the center of a cubic atlas tile at
    every scale. */
inline void test_shader_chain_atlas_full_face_edge_fade() {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  if (!arm_project_op_chain<In::Op::PolyhedralAtlasChainParams>(
          program, In::Op::ProjectCubic::ID, 0, ValueSet::DEFAULTS))
    return;
  auto &projection = param_as<In::Op::PolyhedralAtlasChainParams>(program, 1);
  auto &source = param_as<In::Op::GridSampleParams>(program, 2);
  projection.frame = static_cast<uint8_t>(In::Op::ProjectionFrame::IDENTITY);
  source.coverage_mode =
      static_cast<uint8_t>(In::Op::ProjectionCoverageMode::EDGE_FADE);
  const auto *width_field = schema_field(*program.ops()[2].op, "edge-width");
  HS_EXPECT_TRUE(width_field != nullptr);
  if (!width_field)
    return;
  const auto ctx = shared_resources().context();
  const auto &atlas = projections::POLYHEDRAL_ATLAS<Solids::Cube>;
  const math::Vector face = atlas.selectors[0].normalized();
  math::Vector neighbor = atlas.selectors[1].normalized();
  for (const math::Vector &selector : atlas.selectors)
    if (fabsf(math::dot(selector.normalized(), face)) < 1e-4f)
      neighbor = selector.normalized();
  const math::Vector points[] = {face + neighbor, face + neighbor * 0.5f, face};
  const float COVERAGE[] = {0.0f, 0.5f, 1.0f};
  for (float scale : {0.25f, 1.0f, 4.0f}) {
    projection.coordinate_scale = scale;
    source.edge_width = 0.5f * scale;
    HS_EXPECT_TRUE(source.edge_width <= width_field->max);
    program.prepare(ctx);
    for (size_t index = 0; index < std::size(points); ++index) {
      const PB::SphereSample seed{points[index].normalized(), 0.0f};
      const auto plane = run_entry<PB::PlaneSample>(program, ctx, 1, seed);
      const auto field = run_entry<PB::FieldSample>(program, ctx, 2, plane);
      HS_EXPECT_NEAR(field.coverage, COVERAGE[index], 1e-6f);
    }
  }
  program.clear();
}

inline ProjMirrorFrame project_mirror(In::ChainProgram &program) {
  ProjMirrorFrame mirror;
  const auto &state = state_as<In::Op::SpatialWalkState>(program, 1);
  mirror.conjugate = (math::make_rotation(math::Y_AXIS, state.spin_phase) *
                      PB::projection_frame_base(false) * state.wander)
                         .conjugate();
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::ProjectGnomonic::ID) {
    mirror.gnomonic = param_as<In::Op::GnomonicChainParams>(program, 1);
  } else if (id == In::Op::ProjectBonneV3::ID) {
    mirror.bonne = param_as<In::Op::ProjectBonneV3::Params>(program, 1);
  } else if (id == In::Op::ProjectPeirceSquareFastV3::ID) {
    mirror.peirce_fast =
        param_as<In::Op::ProjectPeirceSquareFastV3::Params>(program, 1);
  } else if (id == In::Op::ProjectPeirceV3::ID) {
    mirror.peirce = param_as<In::Op::ProjectPeirceV3::Params>(program, 1);
  } else if (id == In::Op::ProjectAiroceanV3::ID) {
    mirror.airocean = param_as<In::Op::ProjectAiroceanV3::Params>(program, 1);
  } else if (is_atlas_operator(id)) {
    mirror.atlas = param_as<In::Op::PolyhedralAtlasChainParams>(program, 1);
  } else {
    mirror.meridian = param_as<In::Op::MeridianProjectChainParams>(program, 1);
  }
  return mirror;
}

template <typename Model> inline void expect_project_frame_policy() {
  HS_CONTEXT(Model::ID);
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  if (!arm_project_op_chain<typename Model::Params>(program, Model::ID, 0,
                                                    ValueSet::DEFAULTS))
    return;
  auto &params = param_as<typename Model::Params>(program, 1);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  const auto &prepared = prepared_as<typename Model::Prepared>(program, 1);
  const math::Quaternion base;
  HS_EXPECT_TRUE(prepared.conjugate == base);
  for (const auto &axis : {math::X_AXIS, math::Y_AXIS, math::Z_AXIS}) {
    const auto actual = math::rotate(axis, prepared.conjugate);
    HS_EXPECT_NEAR(actual.x, axis.x, 1e-6f);
    HS_EXPECT_NEAR(actual.y, axis.y, 1e-6f);
    HS_EXPECT_NEAR(actual.z, axis.z, 1e-6f);
  }
  params.spin_rate = 0.03f;
  params.wander = 1.0f;
  for (int frame = 0; frame < 4; ++frame)
    program.advance();
  const auto &state = state_as<In::Op::SpatialWalkState>(program, 1);
  const uint32_t walk_time = state.walk_time;
  const float spin = state.spin_phase;
  const math::Quaternion wander = state.wander;
  params.frame = static_cast<uint8_t>(In::Op::ProjectionFrame::IDENTITY);
  for (int frame = 0; frame < 3; ++frame)
    program.advance();
  HS_EXPECT_EQ(state.walk_time, walk_time);
  HS_EXPECT_EQ(state.spin_phase, spin);
  HS_EXPECT_EQ(std::memcmp(&state.wander, &wander, sizeof(math::Quaternion)),
               0);
  program.prepare(ctx);
  for (const math::Vector &view : sweep_views()) {
    const PB::SphereSample seed{view, 0.25f};
    const auto actual = run_entry<PB::PlaneSample>(program, ctx, 1, seed);
    const auto projected = [&] {
      if constexpr (requires { Model::project(view, params, prepared); })
        return Model::project(view, params, prepared);
      else
        return Model::project(view, params);
    }();
    const auto expected = PB::Kernel::project(seed, view, projected);
    HS_EXPECT_TRUE(plane_identical(actual, expected));
  }
  params.frame = static_cast<uint8_t>(In::Op::ProjectionFrame::SPIN_WANDER);
  program.advance();
  HS_EXPECT_EQ(state.walk_time, walk_time + 1);
  HS_EXPECT_NE(state.spin_phase, spin);
  program.prepare(ctx);
  const math::Quaternion expected =
      (math::make_rotation(math::Y_AXIS, state.spin_phase) *
       PB::projection_frame_base(false) * state.wander)
          .conjugate();
  HS_EXPECT_NE(std::memcmp(&expected, &base, sizeof(math::Quaternion)), 0);
  HS_EXPECT_EQ(
      std::memcmp(&prepared.conjugate, &expected, sizeof(math::Quaternion)), 0);
  params.frame =
      static_cast<uint8_t>(In::Op::ProjectionFrame::SPIN_WANDER_TILTED);
  program.advance();
  HS_EXPECT_EQ(state.walk_time, walk_time + 2);
  program.prepare(ctx);
  const math::Quaternion tilted =
      (math::make_rotation(math::Y_AXIS, state.spin_phase) *
       PB::projection_frame_base(true) * state.wander)
          .conjugate();
  HS_EXPECT_EQ(
      std::memcmp(&prepared.conjugate, &tilted, sizeof(math::Quaternion)), 0);
  const math::Vector pole =
      math::rotate(math::Vector(0, 0, -1), PB::projection_frame_base(true));
  HS_EXPECT_NEAR(pole.x, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(pole.y, -1.0f, 1e-6f);
  HS_EXPECT_NEAR(pole.z, 0.0f, 1e-6f);
  program.clear();
}

inline void test_shader_chain_projection_frame_policy() {
  expect_project_frame_policy<In::Op::ProjectStereographic>();
  expect_project_frame_policy<In::Op::ProjectFoldedSinusoidal>();
  expect_project_frame_policy<In::Op::ProjectEquirectangular>();
  expect_project_frame_policy<In::Op::ProjectGnomonic>();
  expect_project_frame_policy<In::Op::ProjectPeirceV3>();
  expect_project_frame_policy<In::Op::ProjectPeirceSquareFastV3>();
  expect_project_frame_policy<In::Op::ProjectBonneV3>();
  expect_project_frame_policy<In::Op::ProjectAiroceanV3>();
  expect_project_frame_policy<In::Op::ProjectTetrahedral>();
  expect_project_frame_policy<In::Op::ProjectCubic>();
  expect_project_frame_policy<In::Op::ProjectOctahedral>();
  expect_project_frame_policy<In::Op::ProjectDodecahedral>();
  expect_project_frame_policy<In::Op::ProjectIcosahedral>();
}

/** Erased-vs-bound parity of the projection at entry 1. */
template <typename BoundStage>
inline void expect_project_op_parity(In::ChainProgram &program,
                                     const In::FrameContext &ctx) {
  program.prepare(ctx);
  expect_op_parity<BoundStage>(
      program, ctx, 1, project_mirror(program),
      [](const math::Vector &view) { return PB::SphereSample{view, 0.25f}; },
      plane_identical);
}

template <typename BoundStage>
inline void expect_project_parameter_control(In::ChainProgram &program,
                                             const In::FrameContext &ctx,
                                             const ProjMirrorFrame &changed) {
  program.prepare(ctx);
  const ProjMirrorFrame mirror = project_mirror(program);
  const auto prepared = BoundStage::prepare(mirror);
  const auto changed_prepared = BoundStage::prepare(changed);
  int differences = 0;
  for (const auto &view : sweep_views()) {
    const PB::SphereSample seed{view, 0.25f};
    differences +=
        !plane_identical(BoundStage::run(seed, mirror, prepared),
                         BoundStage::run(seed, changed, changed_prepared));
  }
  HS_EXPECT_TRUE(differences > 0);
  expect_project_op_parity<BoundStage>(program, ctx);
}

template <typename Policy, typename OpParams>
inline void run_project_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  if (!arm_project_op_chain<OpParams>(program, op_id, 4, set))
    return;
  const In::FrameContext ctx = shared_resources().context();
  using Bound =
      typename PB::Stage::Project<Policy>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
  if constexpr (requires(OpParams p) { p.coordinate_scale; }) {
    auto &params = param_as<OpParams>(program, 1);
    params = OpParams{};
    params.coordinate_scale = 0.5f;
    auto changed = project_mirror(program);
    changed.peirce_fast.coordinate_scale = 1.0f;
    changed.airocean.coordinate_scale = 1.0f;
    changed.atlas.coordinate_scale = 1.0f;
    expect_project_parameter_control<Bound>(program, ctx, changed);
  }
  program.clear();
}

template <typename Solid> inline void run_atlas_parity(ValueSet set) {
  run_project_parity<PB::Projection::PolyhedralAtlas<AtlasProjMirror, Solid>,
                     In::Op::PolyhedralAtlasChainParams>(
      In::Op::ProjectPolyhedralAtlas<Solid>::ID, set);
}

template <projections::PeirceLayout Layout>
inline void run_peirce_variant(ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  if (!arm_project_op_chain<In::Op::ProjectPeirceV3::Params>(
          program, In::Op::ProjectPeirceV3::ID, 4, set))
    return;
  param_as<In::Op::ProjectPeirceV3::Params>(program, 1).layout =
      static_cast<uint8_t>(Layout);
  using Bound = typename PB::Stage::Project<
      PB::Projection::Peirce<PeirceProjMirror, static_cast<uint8_t>(Layout),
                             true>>::template Bind<ProjMirrorBinding>;
  const auto ctx = shared_resources().context();
  expect_project_op_parity<Bound>(program, ctx);
  auto &params = param_as<In::Op::ProjectPeirceV3::Params>(program, 1);
  params = {};
  params.layout = static_cast<uint8_t>(Layout);
  params.coordinate_scale = 0.5f;
  auto changed = project_mirror(program);
  changed.peirce.coordinate_scale = 1.0f;
  expect_project_parameter_control<Bound>(program, ctx, changed);
  if constexpr (Layout == projections::PeirceLayout::HORIZONTAL ||
                Layout == projections::PeirceLayout::VERTICAL) {
    params.layout_scroll = 0.375f;
    changed = project_mirror(program);
    changed.peirce.layout_scroll = 0.0f;
    expect_project_parameter_control<Bound>(program, ctx, changed);
  }
  program.clear();
}

/** PeirceSquare takes the fast kernel at a zero central meridian and the exact
    square Peirce kernel elsewhere. */
inline void test_shader_chain_peirce_square_dispatch() {
  using Square = typename PB::Stage::Project<PB::Projection::PeirceSquare<
      PeirceProjMirror>>::template Bind<ProjMirrorBinding>;
  using Fast = typename PB::Stage::Project<PB::Projection::PeirceFastSquare<
      PeirceFastProjMirror>>::template Bind<ProjMirrorBinding>;
  using Exact = typename PB::Stage::Project<PB::Projection::Peirce<
      PeirceProjMirror, static_cast<uint8_t>(projections::PeirceLayout::SQUARE),
      true>>::template Bind<ProjMirrorBinding>;
  ProjMirrorFrame frame{};
  frame.conjugate = math::make_rotation(math::Y_AXIS, 0.3f).conjugate();
  frame.peirce.layout_scroll = 0.25f;
  frame.peirce.coordinate_scale = 0.75f;
  frame.peirce.singularity_fade = 0.5f;
  frame.peirce_fast.coordinate_scale = frame.peirce.coordinate_scale;
  frame.peirce_fast.singularity_fade = frame.peirce.singularity_fade;
  for (const float meridian : {0.0f, 0.7f}) {
    HS_CONTEXT(meridian == 0.0f ? "zero meridian" : "moving meridian");
    frame.peirce.central_meridian = meridian;
    const auto square = Square::prepare(frame);
    const auto fast = Fast::prepare(frame);
    const auto exact = Exact::prepare(frame);
    int fast_matches = 0;
    int exact_matches = 0;
    for (const auto &view : sweep_views()) {
      const PB::SphereSample seed{view, 0.25f};
      const PB::PlaneSample got = Square::run(seed, frame, square);
      fast_matches += plane_identical(got, Fast::run(seed, frame, fast));
      exact_matches += plane_identical(got, Exact::run(seed, frame, exact));
    }
    const int views = static_cast<int>(sweep_views().size());
    if (meridian == 0.0f) {
      HS_EXPECT_EQ(fast_matches, views);
      HS_EXPECT_TRUE(exact_matches < views);
    } else {
      HS_EXPECT_EQ(exact_matches, views);
      HS_EXPECT_TRUE(fast_matches < views);
    }
  }
}

inline void test_shader_chain_parity_project_ops() {
  static_assert(PB::Projection::PeirceFastSquare<
                PeirceFastProjMirror>::PROVIDER_VALID<ProjMirrorBinding>);
  static_assert(
      !PB::Projection::PeirceFastSquare<
          PeirceFastMovingMeridianMirror>::PROVIDER_VALID<ProjMirrorBinding>);
  static_assert(!PB::Projection::PeirceFastSquare<
                PeirceProjMirror>::PROVIDER_VALID<ProjMirrorBinding>);
  for_each_value_set([&](ValueSet set) {
    run_project_parity<PB::Projection::FoldedSinusoidal<MeridianProjMirror>,
                       In::Op::ProjectFoldedSinusoidal::Params>(
        In::Op::ProjectFoldedSinusoidal::ID, set);
    run_project_parity<PB::Projection::Equirectangular<MeridianProjMirror>,
                       In::Op::MeridianProjectChainParams>(
        In::Op::ProjectEquirectangular::ID, set);
    run_peirce_variant<projections::PeirceLayout::DIAMOND>(set);
    run_peirce_variant<projections::PeirceLayout::SQUARE>(set);
    run_peirce_variant<projections::PeirceLayout::HORIZONTAL>(set);
    run_peirce_variant<projections::PeirceLayout::VERTICAL>(set);
    run_project_parity<PB::Projection::PeirceFastSquare<PeirceFastProjMirror>,
                       In::Op::ProjectPeirceSquareFastV3::Params>(
        In::Op::ProjectPeirceSquareFastV3::ID, set);
    run_project_parity<
        PB::Projection::Airocean<AiroceanProjMirror, false, true>,
        In::Op::ProjectAiroceanV3::Params>(In::Op::ProjectAiroceanV3::ID, set);
    run_atlas_parity<Solids::Tetrahedron>(set);
    run_atlas_parity<Solids::Cube>(set);
    run_atlas_parity<Solids::Octahedron>(set);
    run_atlas_parity<Solids::Dodecahedron>(set);
    run_atlas_parity<Solids::Icosahedron>(set);
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    if (!arm_project_op_chain<In::Op::ProjectAiroceanV3::Params>(
            program, In::Op::ProjectAiroceanV3::ID, 4, set))
      return;
    param_as<In::Op::ProjectAiroceanV3::Params>(program, 1).layout = 1;
    using Horizontal = typename PB::Stage::Project<PB::Projection::Airocean<
        AiroceanProjMirror, true, true>>::template Bind<ProjMirrorBinding>;
    expect_project_op_parity<Horizontal>(program, shared_resources().context());
    program.clear();
  });
}

template <PB::Projection::GnomonicHemisphere Hemisphere>
inline void run_gnomonic_variant(In::ChainProgram &program,
                                 const In::FrameContext &ctx) {
  param_as<In::Op::GnomonicChainParams>(program, 1).hemisphere =
      static_cast<uint8_t>(Hemisphere);
  using Bound = typename PB::Stage::Project<PB::Projection::Gnomonic<
      GnomonicProjMirror, Hemisphere>>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
}

template <bool North>
inline void run_bonne_variant(In::ChainProgram &program,
                              const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::ProjectBonneV3::Params>(program, 1);
  params.hemisphere = North ? 0 : 1;
  using Bound = typename PB::Stage::Project<PB::Projection::Bonne<
      BonneProjMirror, North>>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
  params = {};
  params.hemisphere = North ? 0 : 1;
  params.coordinate_scale = 0.5f;
  auto changed = project_mirror(program);
  changed.bonne.coordinate_scale = 1.0f;
  expect_project_parameter_control<Bound>(program, ctx, changed);
}

inline void test_shader_chain_parity_project_hemispheres() {
  using Hemisphere = PB::Projection::GnomonicHemisphere;
  for_each_value_set_program([&](In::ChainProgram &program, ValueSet set) {
    if (!arm_project_op_chain<In::Op::GnomonicChainParams>(
            program, In::Op::ProjectGnomonic::ID, 4, set))
      return;
    const In::FrameContext ctx = shared_resources().context();
    run_gnomonic_variant<Hemisphere::FOLDED>(program, ctx);
    run_gnomonic_variant<Hemisphere::FRONT>(program, ctx);
    run_gnomonic_variant<Hemisphere::BACK>(program, ctx);
    program.clear();

    if (!arm_project_op_chain<In::Op::ProjectBonneV3::Params>(
            program, In::Op::ProjectBonneV3::ID, 4, set))
      return;
    run_bonne_variant<true>(program, ctx);
    if (!arm_project_op_chain<In::Op::ProjectBonneV3::Params>(
            program, In::Op::ProjectBonneV3::ID, 4, set))
      return;
    run_bonne_variant<false>(program, ctx);
  });
}
