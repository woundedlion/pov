/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Pullback death cases.

inline void case_pullback_mobius_degenerate() {
  Pullback::Interp::Op::MobiusChainParams params;
  params.a_re = opaque(0.0f);
  params.d_re = opaque(0.0f);
  Pullback::Interp::Op::LensMobius::State state;
  Pullback::Interp::FrameContext context{};
  (void)Pullback::Interp::Op::LensMobius::prepare(context, params, state);
}

/** @brief Death case: the generated-palette operator rejects an unknown
    brightness envelope. */
inline void case_pullback_operator_invalid_brightness_envelope() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.envelope_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
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
  if (result.coords.re != 0.0f)
    std::printf("x");
}

/** @brief Death case: a Sample operator rejects an unknown coverage mode. */
inline void case_pullback_operator_invalid_coverage_mode() {
  Pullback::Interp::Op::GridSampleParams params;
  params.coverage_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::SourceClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::SampleGridV3::prepare(context, params, state)
          .primary != 0.0f)
    std::printf("x");
}

/** @brief Death case: a Sample operator rejects an unknown weight mode. */
inline void case_pullback_operator_invalid_weight_mode() {
  Pullback::Interp::Op::GridSampleParams params;
  params.weight_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::SourceClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::SampleGridV3::prepare(context, params, state)
          .primary != 0.0f)
    std::printf("x");
}

/** @brief Death case: a warp operator rejects an unknown envelope. */
inline void case_pullback_operator_invalid_warp_envelope() {
  Pullback::Interp::Op::WaveShearWarpParams params;
  params.envelope = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::WarpPhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpWaveShear::prepare(context, params, state)
          .phase != 0.0f)
    std::printf("x");
}

/** @brief Death case: the polar chart rejects an unknown polar mode. */
inline void case_pullback_operator_invalid_polar_mode() {
  Pullback::Interp::Op::PolarChartParams params;
  params.mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::WarpPhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpPolarChart::prepare(context, params, state)
          .phase != 0.0f)
    std::printf("x");
}

/** @brief Death case: the polar chart rejects an out-of-range harmonic. */
inline void case_pullback_operator_invalid_polar_harmonic() {
  Pullback::Interp::Op::PolarChartParams params;
  params.harmonic = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::WarpPhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpPolarChart::prepare(context, params, state)
          .phase != 0.0f)
    std::printf("x");
}

/** @brief Death case: the curl-flow operator rejects an unknown integrator. */
inline void case_pullback_operator_invalid_curl_integrator() {
  Pullback::Interp::Op::CurlFlowWarpParams params;
  params.integrator = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::NoisePhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::WarpCurlFlow::prepare(context, params, state)
          .intervals != 0)
    std::printf("x");
}

/** @brief Death case: the curl displacement rejects an unknown integrator. */
inline void case_pullback_operator_invalid_surface_integrator() {
  Pullback::Interp::Op::CurlDisplaceParams params;
  params.integrator = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::NoisePhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::DisplaceCurl::prepare(context, params, state)
          .noise != nullptr)
    std::printf("x");
}

/** @brief Rejects an unknown Bonne hemisphere. */
inline void case_pullback_operator_invalid_bonne_hemisphere() {
  Pullback::Interp::Op::ProjectBonneV3::Params params;
  params.hemisphere = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ProjectBonneV3::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectBonneV3::prepare(context, params, state);
}

/** @brief Death case: the gnomonic projection rejects an unknown hemisphere. */
inline void case_pullback_operator_invalid_gnomonic_hemisphere() {
  Pullback::Interp::Op::GnomonicChainParams params;
  params.hemisphere = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ProjectGnomonic::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectGnomonic::prepare(context, params, state);
}

inline void case_pullback_operator_invalid_airocean_layout() {
  Pullback::Interp::Op::ProjectAiroceanV3::Params params;
  params.layout = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ProjectAiroceanV3::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectAiroceanV3::prepare(context, params, state);
}

inline void case_pullback_operator_invalid_peirce_layout() {
  Pullback::Interp::Op::ProjectPeirceV3::Params params;
  params.layout =
      opaque<uint8_t>(std::size(Pullback::Interp::Op::PEIRCE_LAYOUT_IDS));
  Pullback::Interp::Op::ProjectPeirceV3::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::ProjectPeirceV3::prepare(context, params, state);
}

/** @brief Death case: a noise-driven operator rejects an unknown basis. */
inline void case_pullback_operator_invalid_noise_basis() {
  Pullback::Interp::Op::CurlDisplaceParams params;
  params.basis = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::NoisePhaseState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::DisplaceCurl::prepare(context, params, state)
          .noise != nullptr)
    std::printf("x");
}

/** @brief Death case: the tessellation source rejects an unknown kind. */
inline void case_pullback_operator_invalid_tessellation_kind() {
  Pullback::Interp::Op::TessellationSampleParams params;
  params.kind = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::SourceClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::SampleTessellation::prepare(context, params, state)
          .primary != 0.0f)
    std::printf("x");
}

/** @brief Death case: the kaleidoscope lens rejects an unknown symmetry. */
inline void case_pullback_operator_invalid_kaleidoscope_symmetry() {
  Pullback::Interp::Op::KaleidoscopeChainParams params;
  params.symmetry = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::LensKaleidoscope::State state;
  Pullback::Interp::FrameContext context{};
  Pullback::Interp::Op::LensKaleidoscope::prepare(context, params, state);
}

/** @brief Death case: the generated-palette operator rejects an unknown hue mode. */
inline void case_pullback_operator_invalid_hue_mode() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.hue_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
}

/** @brief Death case: the generated-palette operator rejects an unknown palette
    mode. */
inline void case_pullback_operator_invalid_palette_mode() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.palette_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
}

/** @brief Death case: the generated-palette operator rejects an unknown palette
    mapping. */
inline void case_pullback_operator_invalid_palette_mapping() {
  Pullback::Interp::Op::GeneratedPaletteParams params;
  params.mapping_mode = opaque<uint8_t>(0xff);
  Pullback::Interp::Op::ColorClockState state;
  Pullback::Interp::FrameContext context{};
  if (Pullback::Interp::Op::ColorizeGeneratedPaletteV3::prepare(context, params,
                                                                state)
          .palette != nullptr)
    std::printf("x");
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
  const float NAN_VALUE = opaque(std::numeric_limits<float>::quiet_NaN());
  (void)Pullback::Kernel::sample(Pullback::PlaneSample{}, NAN_VALUE, 1.0f);
}

inline void case_sample_sphere_nan() {
  const float NAN_VALUE = opaque(std::numeric_limits<float>::quiet_NaN());
  (void)Pullback::Kernel::sample(Pullback::SphereSample{}, NAN_VALUE);
}
