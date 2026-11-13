/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Pullback death cases.

inline void chain_invalid_layout(unsigned variant) {
  using namespace Pullback::Interp;
  ChainProgram program;
  OperatorDescriptor descriptor = OPERATOR_TABLE[0];
  alignas(std::max_align_t) uint8_t block_a[512], block_b[512];
  size_t capacity = sizeof(block_a);
  switch (variant) {
  case 0:
    descriptor.runtime.param.align = 0;
    break;
  case 1:
    descriptor.runtime.prepared.align = 3;
    break;
  case 2:
    descriptor.runtime.state.align = 2 * alignof(std::max_align_t);
    break;
  case 3:
    descriptor.runtime.param.size = 0;
    break;
  case 4:
    descriptor.runtime.state.size = 3;
    break;
  case 5:
    capacity = std::numeric_limits<uint32_t>::max();
    break;
  }
  program.bind_storage(block_a, block_b, capacity, {&descriptor, 1});
}

/** @brief Registered field curves retain their catalog spellings. */
inline void check_catalog_field_curve_companion() {
  using Curve = Control::FieldCurve;
  const std::pair<Curve, const char *> CURVES[] = {
      {Curve::LERP, "lerp"},
      {Curve::LOG_POSITIVE, "log-positive"},
      {Curve::SHORTEST_PERIODIC, "shortest-periodic"},
      {Curve::SHORTEST_TURN, "shortest-turn"},
      {Curve::SNAP, "snap"}};
  for (const auto &[curve, name] : CURVES)
    HS_EXPECT(std::strcmp(Pullback::Interp::Detail::curve_name(curve), name) ==
                  0,
              "registered field curve spelling");
}

/** @brief Catalog export rejects model-only midpoint curves. */
inline void case_catalog_curve_midpoint() {
  (void)Pullback::Interp::Detail::curve_name(Control::FieldCurve::MIDPOINT);
}

/** @brief Catalog export rejects model-only raw-linear curves. */
inline void case_catalog_curve_raw_linear() {
  (void)Pullback::Interp::Detail::curve_name(Control::FieldCurve::RAW_LINEAR);
}

/** @brief Catalog export rejects unknown field curves. */
inline void case_catalog_curve_unknown() {
  (void)Pullback::Interp::Detail::curve_name(
      static_cast<Control::FieldCurve>(opaque<uint8_t>(0xff)));
}

/** @brief Prepares a valid palette before rejecting one selector.
 * @tparam Field Palette selector member. */
template <auto Field> void prepare_palette_with_invalid() {
  using Op = Pullback::Interp::Op::ColorizeGeneratedPaletteV4;
  Op::Params params;
  Op::State state;
  BakedPalette palette;
  Pullback::Interp::FrameContext context;
  context.palettes.fill(&palette);
  (void)Op::prepare(context, params, state);
  params.*Field = opaque<uint8_t>(0xff);
  death_admitted();
  (void)Op::prepare(context, params, state);
}

/** @brief Death case: the generated-palette operator rejects an unknown
    brightness envelope. */
inline void case_pullback_operator_invalid_brightness_envelope() {
  prepare_palette_with_invalid<
      &Pullback::Interp::Op::GeneratedPaletteParams::envelope_mode>();
}

/**
 * @brief Death case: an operator table whose entry decreases carrier family
 *        rank must trap.
 * @details compile() only matches adjacent carriers.
 */
inline void case_chain_table_rank_decreases() {
  static Pullback::Interp::ChainProgram program;
  static Pullback::Interp::OperatorDescriptor descriptor{};
  descriptor.input = Pullback::Interp::CarrierId::COLOR;
  descriptor.output = Pullback::Interp::CarrierId::SPHERE;
  alignas(std::max_align_t) static uint8_t block_a[64];
  alignas(std::max_align_t) static uint8_t block_b[64];
  program.bind_storage(
      block_a, block_b, opaque<size_t>(sizeof(block_a)),
      std::span<const Pullback::Interp::OperatorDescriptor>(&descriptor, 1));
}

inline void case_chain_overlapping_storage() {
  Pullback::Interp::ChainProgram program;
  alignas(std::max_align_t) uint8_t block[128];
  program.bind_storage(block, block + alignof(std::max_align_t), 64);
}

inline void case_chain_identical_storage() {
  Pullback::Interp::ChainProgram program;
  alignas(std::max_align_t) uint8_t block[64];
  program.bind_storage(block, block, sizeof(block));
}

inline void case_chain_zero_alignment() { chain_invalid_layout(0); }

inline void case_chain_non_power_alignment() { chain_invalid_layout(1); }

inline void case_chain_overaligned_block() { chain_invalid_layout(2); }

inline void case_chain_zero_size() { chain_invalid_layout(3); }

inline void case_chain_misaligned_size() { chain_invalid_layout(4); }

inline void case_chain_capacity_overflow() { chain_invalid_layout(5); }

inline void case_pullback_project_nonunit_direction() {
  using namespace hs_test::pullback_tests;
  using Project =
      Pullback::Stage::Project<CountingProjectionPolicy>::Bind<TestBinding>;
  const TestFrame frame;
  const Pullback::SphereSample input{math::Vector(opaque(2.0f), 0.0f, 0.0f),
                                     0.0f};
  const auto result = Project::run(input, frame, Project::prepare(frame));
  (void)result;
}

/** @brief Traps when a factored column basis produces a nonunit direction. */
inline void case_pullback_column_nonunit_direction() {
  using FX = GridSpace<96, 20>;
  using Pipeline = FX::RenderPipeline;
  reset_globals();
  FX effect;
  effect.init();
  const auto frame = Pipeline::prepare(effect.frame_for_test());
  const Pipeline::ColumnBasis basis{math::Vector(opaque(2.0f), 0, 0),
                                    math::Vector(0, 1, 0)};
  const auto result = Pipeline::shade_column(1.0f, 0.0f, frame, basis);
  (void)result;
}

/** @brief Prepares an operator with an invalid enum field. */
template <typename Op, auto Field>
inline void prepare_operator_with(uint8_t value = 0xff) {
  typename Op::Params params;
  params.*Field = opaque(value);
  typename Op::State state;
  Pullback::Interp::FrameContext context{};
  (void)Op::prepare(context, params, state);
}

/** @brief Death case: a Sample operator rejects an unknown coverage mode. */
inline void case_pullback_operator_invalid_coverage_mode() {
  prepare_operator_with<
      Pullback::Interp::Op::SampleGridV4,
      &Pullback::Interp::Op::GridSampleParams::coverage_mode>();
}

/** @brief Death case: a Sample operator rejects an unknown weight mode. */
inline void case_pullback_operator_invalid_weight_mode() {
  prepare_operator_with<Pullback::Interp::Op::SampleGridV4,
                        &Pullback::Interp::Op::GridSampleParams::weight_mode>();
}

/** @brief Death case: a warp operator rejects an unknown envelope. */
inline void case_pullback_operator_invalid_warp_envelope() {
  prepare_operator_with<Pullback::Interp::Op::WarpWaveShear,
                        &Pullback::Interp::Op::WaveShearWarpParams::envelope>();
}

/** @brief Death case: the polar chart rejects an unknown polar mode. */
inline void case_pullback_operator_invalid_polar_mode() {
  prepare_operator_with<Pullback::Interp::Op::WarpPolarChart,
                        &Pullback::Interp::Op::PolarChartParams::mode>();
}

/** @brief Death case: the polar chart rejects an out-of-range harmonic. */
inline void case_pullback_operator_invalid_polar_harmonic() {
  prepare_operator_with<Pullback::Interp::Op::WarpPolarChart,
                        &Pullback::Interp::Op::PolarChartParams::harmonic>();
}

/** @brief Death case: the curl-flow operator rejects an unknown integrator. */
inline void case_pullback_operator_invalid_curl_integrator() {
  prepare_operator_with<
      Pullback::Interp::Op::WarpCurlFlow,
      &Pullback::Interp::Op::CurlFlowWarpParams::integrator>();
}

/** @brief Death case: the curl displacement rejects an unknown integrator. */
inline void case_pullback_operator_invalid_surface_integrator() {
  prepare_operator_with<
      Pullback::Interp::Op::DisplaceCurl,
      &Pullback::Interp::Op::CurlDisplaceParams::integrator>();
}

/** @brief Death case: the curl displacement rejects an unknown tangent limit. */
inline void case_pullback_operator_invalid_tangent_limit() {
  prepare_operator_with<Pullback::Interp::Op::DisplaceCurl,
                        &Pullback::Interp::Op::CurlDisplaceParams::limit>();
}

/** @brief Ripple advance rejects an oversized burst. */
inline void case_pullback_operator_invalid_ripple_burst() {
  using Ripple = Pullback::Interp::Op::DisplaceRipple;
  Ripple::Params params;
  Ripple::State state;
  params.burst_size = opaque<uint8_t>(4);
  Ripple::advance(state, params);
}

/** @brief Tangent displacement preparation rejects an unknown integrator. */
inline void case_pullback_operator_invalid_field_integrator() {
  using Field = Pullback::Interp::Op::DisplacePotential;
  Field::Params params;
  params.integrator = opaque<uint8_t>(0xff);
  (void)Field::prepare({}, params, {});
}

/** @brief Symmetric displacement rejects an unknown symmetry. */
inline void case_pullback_operator_invalid_field_symmetry() {
  using Field = Pullback::Interp::Op::DisplaceSymmetric;
  Field::Params params;
  params.symmetry = opaque<uint8_t>(0xff);
  (void)Field::prepare({}, params, {});
}

/** @brief Symmetric displacement rejects an unknown field kind. */
inline void case_pullback_operator_invalid_field_kind() {
  using Field = Pullback::Interp::Op::DisplaceSymmetric;
  Field::Params params;
  params.kind = opaque<uint8_t>(0xff);
  (void)Field::prepare({}, params, {});
}

/** @brief Rejects an unknown Bonne hemisphere. */
inline void case_pullback_operator_invalid_bonne_hemisphere() {
  prepare_operator_with<
      Pullback::Interp::Op::ProjectBonneV3,
      &Pullback::Interp::Op::ProjectBonneV3::Params::hemisphere>();
}

/** @brief Death case: the gnomonic projection rejects an unknown hemisphere. */
inline void case_pullback_operator_invalid_gnomonic_hemisphere() {
  prepare_operator_with<
      Pullback::Interp::Op::ProjectGnomonic,
      &Pullback::Interp::Op::GnomonicChainParams::hemisphere>();
}

inline void case_pullback_operator_invalid_airocean_layout() {
  prepare_operator_with<
      Pullback::Interp::Op::ProjectAiroceanV3,
      &Pullback::Interp::Op::ProjectAiroceanV3::Params::layout>();
}

inline void case_pullback_operator_invalid_peirce_layout() {
  prepare_operator_with<Pullback::Interp::Op::ProjectPeirceV3,
                        &Pullback::Interp::Op::ProjectPeirceV3::Params::layout>(
      std::size(Pullback::Interp::Op::PEIRCE_LAYOUT_IDS));
}

/** @brief Death case: a noise-driven operator rejects an unknown basis. */
inline void case_pullback_operator_invalid_noise_basis() {
  prepare_operator_with<Pullback::Interp::Op::DisplaceCurl,
                        &Pullback::Interp::Op::CurlDisplaceParams::basis>();
}

/** @brief Death case: the tessellation source rejects an unknown kind. */
inline void case_pullback_operator_invalid_tessellation_kind() {
  prepare_operator_with<
      Pullback::Interp::Op::SampleTessellation,
      &Pullback::Interp::Op::TessellationSampleParams::kind>();
}

/** @brief Death case: the tessellation source rejects an unknown fill. */
inline void case_pullback_operator_invalid_tessellation_fill() {
  prepare_operator_with<
      Pullback::Interp::Op::SampleTessellation,
      &Pullback::Interp::Op::TessellationSampleParams::fill>();
}

/** @brief Death case: the kaleidoscope lens rejects an unknown symmetry. */
inline void case_pullback_operator_invalid_kaleidoscope_symmetry() {
  prepare_operator_with<
      Pullback::Interp::Op::LensKaleidoscope,
      &Pullback::Interp::Op::KaleidoscopeChainParams::symmetry>();
}

/** @brief Authored kaleidoscope dispatch rejects an unknown symmetry. */
inline void case_authored_kaleidoscope_invalid_symmetry() {
  (void)Pullback::Lens::configured_kaleidoscope(math::Y_AXIS,
                                                opaque<uint8_t>(0xff));
}

/** @brief A modulation operator rejects an unknown source. */
inline void case_pullback_operator_invalid_modulation_source() {
  Pullback::Interp::Op::ModulationParams params;
  (void)Pullback::Modulation::init(params.modulation);
  death_admitted();
  params.modulation.source =
      static_cast<Pattern::SourceKind>(opaque<uint8_t>(0xff));
  (void)Pullback::Modulation::init(params.modulation);
}

/** @brief Death case: the generated-palette operator rejects an unknown palette
    mode. */
inline void case_pullback_operator_invalid_palette_mode() {
  prepare_palette_with_invalid<
      &Pullback::Interp::Op::GeneratedPaletteParams::palette_mode>();
}

/** @brief Death case: the generated-palette operator rejects an unknown palette
    mapping. */
inline void case_pullback_operator_invalid_palette_mapping() {
  prepare_palette_with_invalid<
      &Pullback::Interp::Op::GeneratedPaletteParams::mapping_mode>();
}

/** @brief A field transfer must preserve the normalized value domain. */
inline void case_field_transfer_outside_range() {
  Pullback::Kernel::transfer(Pullback::FieldSample{}, opaque(1.1f));
}

/** @brief Field coverage cannot grow after the sample crossing. */
inline void case_field_coverage_increases() {
  Pullback::Kernel::coverage(Pullback::FieldSample{}, opaque(1.1f));
}

/** @brief The Sample crossing rejects a coverage factor outside [0, 1]. */
inline void case_field_sample_coverage_outside_range() {
  Pullback::Kernel::sample(Pullback::PlaneSample{}, 0.0f, opaque(1.1f));
}

inline void case_projection_invalid_frame_advance() {
  using Op = Pullback::Interp::Op::ProjectStereographic;
  Op::Params params;
  params.frame = opaque<uint8_t>(255);
  Op::State state;
  Op::advance(state, params);
}

inline void case_projection_invalid_frame_prepare() {
  using Op = Pullback::Interp::Op::ProjectStereographic;
  Op::Params params;
  params.frame = opaque<uint8_t>(255);
  Op::State state;
  (void)Op::prepare(Pullback::Interp::FrameContext{}, params, state);
}

inline void case_sample_plane_nan() {
  const float nan_value = opaque(std::numeric_limits<float>::quiet_NaN());
  (void)Pullback::Kernel::sample(Pullback::PlaneSample{}, nan_value, 1.0f);
}

inline void case_sample_sphere_nan() {
  const float nan_value = opaque(std::numeric_limits<float>::quiet_NaN());
  (void)Pullback::Kernel::sample(Pullback::SphereSample{}, nan_value);
}

/** @brief Arena-sized blocks for column-time chain death cases. */
struct ChainTemporalBlocks {
  alignas(std::max_align_t)
      uint8_t a[Pullback::Interp::CHAIN_ARENA_BYTES]; ///< First program arena.
  alignas(std::max_align_t)
      uint8_t b[Pullback::Interp::CHAIN_ARENA_BYTES]; ///< Second program arena.
  alignas(std::max_align_t) uint8_t
      history[Pullback::Interp::CHAIN_ARENA_BYTES]; ///< Retained endpoint.
  alignas(std::max_align_t)
      uint8_t column[Pullback::Interp::CHAIN_ARENA_BYTES]; ///< Column image.
};

/** @brief Compiles a three-operator chain into static storage.
 * @param program Program to compile.
 * @param blocks Storage to bind; the temporal blocks are bound when
 *        @p temporal is true.
 * @param temporal Whether to bind column-time storage. */
inline void compile_temporal_chain(Pullback::Interp::ChainProgram &program,
                                   ChainTemporalBlocks &blocks, bool temporal) {
  namespace In = Pullback::Interp;
  program.bind_storage(blocks.a, blocks.b);
  if (temporal)
    program.bind_temporal_storage(blocks.history, blocks.column);
  const In::ChainEntryRequest REQUEST[] = {
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v4"},
      {"color", "colorize.generated-palette.v4"}};
  if (program.compile(REQUEST).code != In::ChainStatus::OK)
    std::exit(0);
}

/** @brief Column-time storage requires bound program storage. */
inline void case_chain_temporal_before_storage() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_temporal_storage(blocks.history, blocks.column);
}

/** @brief Column-time storage rejects a null block. */
inline void case_chain_temporal_null_block() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_storage(blocks.a, blocks.b);
  program.bind_temporal_storage(blocks.history, {});
}

/** @brief Column-time storage rejects a misaligned block. */
inline void case_chain_temporal_misaligned_block() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_storage(blocks.a, blocks.b);
  program.bind_temporal_storage(std::span(blocks.history).subspan(1),
                                blocks.column);
}

/** @brief Column-time storage rejects a short block. */
inline void case_chain_temporal_short_block() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_storage(blocks.a, blocks.b);
  program.bind_temporal_storage(std::span(blocks.history).first(1),
                                blocks.column);
}

/** @brief The temporal blocks must be disjoint from one another. */
inline void case_chain_temporal_overlap() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_storage(blocks.a, blocks.b);
  program.bind_temporal_storage(blocks.history, blocks.history);
}

/** @brief A temporal block must be disjoint from both program blocks. */
inline void case_chain_temporal_program_overlap() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_storage(blocks.a, blocks.b);
  program.bind_temporal_storage(blocks.history, blocks.b);
}

/** @brief Column-time storage binds once. */
inline void case_chain_temporal_rebinding() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_storage(blocks.a, blocks.b);
  program.bind_temporal_storage(blocks.history, blocks.column);
  program.bind_temporal_storage(blocks.history, blocks.column);
}

/** @brief Retaining an endpoint requires column-time storage. */
inline void case_chain_capture_without_temporal_storage() {
  static ChainTemporalBlocks blocks;
  static Pullback::Interp::ChainProgram program;
  compile_temporal_chain(program, blocks, false);
  program.capture_history();
}

/** @brief Retaining an endpoint requires a committed program. */
inline void case_chain_capture_without_program() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  program.bind_storage(blocks.a, blocks.b);
  program.bind_temporal_storage(blocks.history, blocks.column);
  program.capture_history();
}

/** @brief A column requires a retained endpoint of the committed program. */
inline void case_chain_column_without_history() {
  static ChainTemporalBlocks blocks;
  static Pullback::Interp::ChainProgram program;
  compile_temporal_chain(program, blocks, true);
  const Pullback::Interp::FrameContext CTX{};
  program.prepare_column(CTX, CTX, 0.5f);
}

/** @brief Lifecycle event invalidating sampled history. */
enum class ColumnReset { ADVANCE, PREPARE, CLEAR_HISTORY };

/** @brief Frame lifecycle changes invalidate cached column continuity. */
template <ColumnReset Reset> inline void chain_column_rechecks_history() {
  static ChainTemporalBlocks blocks;
  static Pullback::Interp::ChainProgram program;
  compile_temporal_chain(program, blocks, true);
  BakedPalette palette;
  Pullback::Interp::FrameContext ctx{};
  ctx.palettes.fill(&palette);
  program.prepare(ctx);
  program.capture_history();
  program.prepare_column(ctx, ctx, 0.5f);
  if constexpr (Reset == ColumnReset::CLEAR_HISTORY) {
    program.clear_history();
  } else {
    bool changed = false;
    const auto ops = program.ops();
    for (size_t index = 0; index < ops.size() && !changed; ++index) {
      const auto &op = ops[index];
      for (uint16_t field = 0; field < op.op->schema_count; ++field) {
        if (op.op->schema[field].enum_count > 1 &&
            (op.op->runtime.temporal_hold_field == nullptr ||
             !op.op->runtime.temporal_hold_field(field))) {
          auto *value = static_cast<uint8_t *>(
              op.op->runtime.param_address(program.param_block(index), field));
          *value = static_cast<uint8_t>((*value + 1) %
                                        op.op->schema[field].enum_count);
          changed = true;
          break;
        }
      }
    }
    if (!changed)
      std::exit(0);
    if constexpr (Reset == ColumnReset::ADVANCE)
      program.advance();
    else
      program.prepare(ctx);
  }
  program.prepare_column(ctx, ctx, 0.75f);
}

inline void case_chain_column_after_advance() {
  chain_column_rechecks_history<ColumnReset::ADVANCE>();
}

inline void case_chain_column_after_prepare() {
  chain_column_rechecks_history<ColumnReset::PREPARE>();
}

inline void case_chain_column_after_clear_history() {
  chain_column_rechecks_history<ColumnReset::CLEAR_HISTORY>();
}

/** @brief Readiness requires a committed chain. */
inline void case_chain_ready_without_program() {
  Pullback::Interp::ChainProgram program;
  program.check_ready();
}

/** @brief A committed chain requires preparation before evaluation. */
inline void case_chain_ready_before_prepare() {
  static ChainTemporalBlocks blocks;
  Pullback::Interp::ChainProgram program;
  compile_temporal_chain(program, blocks, false);
  program.check_ready();
}
