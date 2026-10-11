/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- field-op mirrors -----------------------------------------------------

/** Mirror frame for the FIELD endomorphism batch. */
struct FieldMirrorFrame {
  In::Op::IsoContourChainParams iso;
  In::Op::SmoothBandsChainParams bands;
  In::Op::ValueCutoutChainParams cutout;
};

struct FieldMirrorBinding {
  using FrameState = FieldMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct IsoContourFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float iso_level(const FieldMirrorFrame &frame) {
    return frame.iso.iso_level;
  }
  static float iso_width(const FieldMirrorFrame &frame) {
    return frame.iso.iso_width;
  }
};

struct SmoothBandsFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float band_count(const FieldMirrorFrame &frame) {
    return frame.bands.band_count;
  }
  static float band_phase(const FieldMirrorFrame &frame) {
    return frame.bands.band_phase;
  }
};

struct ValueCutoutFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float cutout_threshold(const FieldMirrorFrame &frame) {
    return frame.cutout.cutout_threshold;
  }
  static float cutout_softness(const FieldMirrorFrame &frame) {
    return frame.cutout.cutout_softness;
  }
};

/** Compiles a chain with one FIELD endomorphism at entry 3 and applies the
    value set to every block. */
template <typename OpParams>
inline void arm_field_op_chain(In::ChainProgram &program, const char *op_id,
                               int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v3"},
      {"op", op_id},
      {"colorize", "colorize.generated-palette.v3"},
  };
  arm_chain<In::Op::RotateChainParams, In::Op::ProjectChainParams,
            In::Op::GridSampleParams, OpParams, In::Op::GeneratedPaletteParams>(
      program, chain, frames, set);
}

/** The FIELD op's input carrier for @p view: the erased camera, projection
    and sample entries run over the seed. */
inline PB::FieldSample field_input(In::ChainProgram &program,
                                   const In::FrameContext &ctx,
                                   const math::Vector &view) {
  return run_entry<PB::FieldSample>(program, ctx, 2,
                                    projected_input(program, ctx, view));
}

inline FieldMirrorFrame field_mirror(In::ChainProgram &program) {
  FieldMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[3].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::TransferIsoContour::ID)
    mirror.iso = param_as<In::Op::IsoContourChainParams>(program, 3);
  else if (id == In::Op::TransferSmoothBands::ID)
    mirror.bands = param_as<In::Op::SmoothBandsChainParams>(program, 3);
  else if (id == In::Op::CoverageValueCutout::ID)
    mirror.cutout = param_as<In::Op::ValueCutoutChainParams>(program, 3);
  return mirror;
}

/** Erased-vs-bound parity of the FIELD endomorphism at entry 3. */
template <typename BoundStage>
inline void expect_field_op_parity(In::ChainProgram &program,
                                   const In::FrameContext &ctx) {
  program.prepare(ctx);
  expect_op_parity<BoundStage>(
      program, ctx, 3, field_mirror(program),
      [&](const math::Vector &view) { return field_input(program, ctx, view); },
      field_identical);
}

template <typename BoundStage, typename OpParams>
inline void run_field_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_field_op_chain<OpParams>(program, op_id, 4, set);
  const In::FrameContext ctx = shared_resources().context();
  expect_field_op_parity<BoundStage>(program, ctx);
  program.clear();
}

inline void test_shader_chain_parity_field_ops() {
  for_each_value_set([&](ValueSet set) {
    run_field_parity<typename PB::Stage::Transfer<PB::Transfer::Ridge>::
                         template Bind<FieldMirrorBinding>,
                     PB::Transfer::NoValueParams>(In::Op::TransferRidge::ID,
                                                  set);
    run_field_parity<
        typename PB::Stage::Transfer<PB::Transfer::IsoContour<
            IsoContourFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::IsoContourChainParams>(In::Op::TransferIsoContour::ID, set);
    run_field_parity<
        typename PB::Stage::Transfer<PB::Transfer::SmoothBands<
            SmoothBandsFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::SmoothBandsChainParams>(In::Op::TransferSmoothBands::ID, set);
    run_field_parity<
        typename PB::Stage::ApplyCoverage<PB::ValueCoverage::ValueCutout<
            ValueCutoutFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::ValueCutoutChainParams>(In::Op::CoverageValueCutout::ID, set);
  });
}
