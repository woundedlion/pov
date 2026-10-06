/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_shader_chain.h.

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
};

struct ProjMirrorBinding {
  using FrameState = ProjMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

template <typename Derived> struct ProjMirrorBase {
  using Binding = ProjMirrorBinding;
  using FrameState = ProjMirrorFrame;
  static const math::Quaternion &conjugate(const ProjMirrorFrame &frame) {
    return frame.conjugate;
  }
};

struct MeridianProjMirror : ProjMirrorBase<MeridianProjMirror> {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.meridian.central_meridian;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.meridian.singularity_fade;
  }
};

struct GnomonicProjMirror : ProjMirrorBase<GnomonicProjMirror> {
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.gnomonic.singularity_fade;
  }
};

struct PeirceProjMirror : ProjMirrorBase<PeirceProjMirror> {
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

struct PeirceFastProjMirror : ProjMirrorBase<PeirceFastProjMirror> {
  static constexpr bool ZERO_CENTRAL_MERIDIAN = true;
  static float coordinate_scale(const ProjMirrorFrame &frame) {
    return frame.peirce_fast.coordinate_scale;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.peirce_fast.singularity_fade;
  }
};

struct BonneProjMirror : ProjMirrorBase<BonneProjMirror> {
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

struct AiroceanProjMirror : ProjMirrorBase<AiroceanProjMirror> {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.airocean.central_meridian;
  }
  static float coordinate_scale(const ProjMirrorFrame &frame) {
    return frame.airocean.coordinate_scale;
  }
};

/** Compiles a chain with one projection at entry 1 and applies the value set
    to every block. */
template <typename OpParams>
inline void arm_project_op_chain(In::ChainProgram &program, const char *op_id,
                                 int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", op_id},
      {"sample", "sample.grid.v3"},
      {"colorize", "colorize.generated-palette.v3"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<OpParams>(program, 1), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 3), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

inline ProjMirrorFrame project_mirror(In::ChainProgram &program,
                                      const In::FrameContext &ctx) {
  ProjMirrorFrame mirror;
  const auto &state = state_as<In::Op::SpatialWalkState>(program, 1);
  mirror.conjugate = (math::make_rotation(math::Y_AXIS, state.spin_phase) *
                      ctx.projection_base * state.wander)
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
  } else {
    mirror.meridian = param_as<In::Op::MeridianProjectChainParams>(program, 1);
  }
  return mirror;
}

template <typename Model> inline void expect_project_frame_policy() {
  HS_CONTEXT(Model::ID);
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_project_op_chain<typename Model::Params>(program, Model::ID, 0,
                                               ValueSet::DEFAULTS);
  auto &params = param_as<typename Model::Params>(program, 1);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  const auto &prepared = *reinterpret_cast<const typename Model::Prepared *>(
      program.prepared_block(1));
  const math::Quaternion base = ctx.projection_base.conjugate();
  HS_EXPECT_TRUE(prepared.conjugate == base);
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
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  for (const math::Vector &view : sweep_views()) {
    const PB::SphereSample seed{view, 0.25f};
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&seed, out, ctx, program.param_block(1),
                   program.prepared_block(1));
    const auto &actual =
        *std::launder(reinterpret_cast<PB::PlaneSample *>(out));
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
       ctx.projection_base * state.wander)
          .conjugate();
  HS_EXPECT_NE(std::memcmp(&expected, &base, sizeof(math::Quaternion)), 0);
  HS_EXPECT_EQ(
      std::memcmp(&prepared.conjugate, &expected, sizeof(math::Quaternion)), 0);
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
}

/** Erased-vs-bound parity of the projection at entry 1. */
template <typename BoundStage>
inline void expect_project_op_parity(In::ChainProgram &program,
                                     const In::FrameContext &ctx) {
  program.prepare(ctx);
  const ProjMirrorFrame mirror = project_mirror(program, ctx);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  int view_index = 0;
  for (const math::Vector &view : sweep_views()) {
    HS_CONTEXT("view", view_index++);
    const PB::SphereSample seed{view, 0.25f};
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&seed, out, ctx, program.param_block(1),
                   program.prepared_block(1));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::PlaneSample *>(out));
    const PB::PlaneSample reference = BoundStage::run(seed, mirror, prepared);
    HS_EXPECT_TRUE(plane_identical(erased, reference));
  }
}

template <typename BoundStage>
inline void expect_project_parameter_control(In::ChainProgram &program,
                                             const In::FrameContext &ctx,
                                             const ProjMirrorFrame &changed) {
  program.prepare(ctx);
  const ProjMirrorFrame mirror = project_mirror(program, ctx);
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
  arm_project_op_chain<OpParams>(program, op_id, 4, set);
  const In::FrameContext ctx = shared_resources().context();
  using Bound =
      typename PB::Stage::Project<Policy>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
  if constexpr (requires(OpParams p) { p.coordinate_scale; }) {
    auto &params = param_as<OpParams>(program, 1);
    params = OpParams{};
    params.coordinate_scale = 0.5f;
    auto changed = project_mirror(program, ctx);
    changed.peirce.coordinate_scale = 1.0f;
    changed.peirce_fast.coordinate_scale = 1.0f;
    changed.bonne.coordinate_scale = 1.0f;
    changed.airocean.coordinate_scale = 1.0f;
    expect_project_parameter_control<Bound>(program, ctx, changed);
  }
  program.clear();
}

template <projections::PeirceLayout Layout>
inline void run_peirce_variant(ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  auto &program = fixture->program;
  arm_project_op_chain<In::Op::ProjectPeirceV3::Params>(
      program, In::Op::ProjectPeirceV3::ID, 4, set);
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
  auto changed = project_mirror(program, ctx);
  changed.peirce.coordinate_scale = 1.0f;
  expect_project_parameter_control<Bound>(program, ctx, changed);
  if constexpr (Layout == projections::PeirceLayout::HORIZONTAL ||
                Layout == projections::PeirceLayout::VERTICAL) {
    params.layout_scroll = 0.375f;
    changed = project_mirror(program, ctx);
    changed.peirce.layout_scroll = 0.0f;
    expect_project_parameter_control<Bound>(program, ctx, changed);
  }
  program.clear();
}

inline void test_shader_chain_parity_project_ops() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    run_project_parity<PB::Projection::FoldedSinusoidal<MeridianProjMirror>,
                       In::Op::MeridianProjectChainParams>(
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
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_project_op_chain<In::Op::ProjectAiroceanV3::Params>(
        program, In::Op::ProjectAiroceanV3::ID, 4, set);
    param_as<In::Op::ProjectAiroceanV3::Params>(program, 1).layout = 1;
    using Horizontal = typename PB::Stage::Project<PB::Projection::Airocean<
        AiroceanProjMirror, true, true>>::template Bind<ProjMirrorBinding>;
    expect_project_op_parity<Horizontal>(program, shared_resources().context());
    program.clear();
  }
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
  param_as<In::Op::ProjectBonneV3::Params>(program, 1).hemisphere =
      North ? 0 : 1;
  using Bound = typename PB::Stage::Project<PB::Projection::Bonne<
      BonneProjMirror, North>>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
}

inline void test_shader_chain_parity_project_hemispheres() {
  using Hemisphere = PB::Projection::GnomonicHemisphere;
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    HS_CONTEXT(value_set_name(set));
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_project_op_chain<In::Op::GnomonicChainParams>(
        program, In::Op::ProjectGnomonic::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_gnomonic_variant<Hemisphere::FOLDED>(program, ctx);
    run_gnomonic_variant<Hemisphere::FRONT>(program, ctx);
    run_gnomonic_variant<Hemisphere::BACK>(program, ctx);
    program.clear();

    arm_project_op_chain<In::Op::ProjectBonneV3::Params>(
        program, In::Op::ProjectBonneV3::ID, 4, set);
    run_bonne_variant<true>(program, ctx);
    run_bonne_variant<false>(program, ctx);
    program.clear();
  }
}
