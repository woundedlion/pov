/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- template mirror of the slice chain -----------------------------------

struct MirrorFrame {
  math::Quaternion camera_conjugate;
  math::Quaternion projection_conjugate;
  In::Op::ProjectChainParams projection;
  In::Op::GridSampleParams sample;
  In::Op::GeneratedPaletteParams color;
  float source_primary = 0.0f;
  float source_secondary = 0.0f;
  float source_angle = 0.0f;
  float oscillation_phase = 0.0f;
  const BakedPalette *palette = nullptr;
  const Pixel *hue_rotation_lut = nullptr;
  const int8_t *hue_noise_lut = nullptr;
};

struct MirrorBinding {
  using FrameState = MirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct MirrorCamera {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static const math::Quaternion &conjugate(const MirrorFrame &frame) {
    return frame.camera_conjugate;
  }
};

struct MirrorProjection {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static const math::Quaternion &conjugate(const MirrorFrame &frame) {
    return frame.projection_conjugate;
  }
  static float singularity_fade(const MirrorFrame &frame) {
    return frame.projection.singularity_fade;
  }
};

struct MirrorSource {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static const PB::Source::GridSourceParams &params(const MirrorFrame &frame) {
    return frame.sample;
  }
  static PB::Source::PreparedSource prepare(const MirrorFrame &frame) {
    return PB::Source::prepare(frame.source_primary, frame.source_secondary,
                               frame.source_angle);
  }
};

struct MirrorValue {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static float edge_width(const MirrorFrame &frame) {
    return frame.sample.edge_width;
  }
};

template <PB::Color::HueMode HueV, PB::Color::BrightnessEnvelope BrightV>
struct MirrorColor {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static PB::Color::PaletteMappingWeights
  mapping_weights(const MirrorFrame &frame) {
    return PB::Color::PaletteMappingWeights::single(
        static_cast<PB::Color::PaletteMapping>(frame.color.mapping_mode));
  }
  static float mapping_frequency(const MirrorFrame &frame) {
    return frame.color.mapping_frequency;
  }
  static float mapping_phase(const MirrorFrame &frame) {
    return frame.color.mapping_phase;
  }
  static float oscillation_depth(const MirrorFrame &frame) {
    return frame.color.phase_oscillation_depth;
  }
  static float oscillation_phase(const MirrorFrame &frame) {
    return frame.oscillation_phase;
  }
  static const BakedPalette &palette(const MirrorFrame &frame) {
    return *frame.palette;
  }
  static PB::Color::HueMode hue_mode(const MirrorFrame &) { return HueV; }
  static float hue_shift_amount(const MirrorFrame &frame) {
    return frame.color.hue_shift_amount;
  }
  static PB::Color::HueRotationLutView hue_rotation(const MirrorFrame &frame) {
    return {frame.hue_rotation_lut, frame.color.hue_shift_amount != 0.0f};
  }
  static PB::Color::HueNoiseLutView hue_noise(const MirrorFrame &frame) {
    return {frame.hue_noise_lut, HueV == PB::Color::HueMode::NOISE &&
                                     frame.color.hue_shift_amount != 0.0f};
  }
  static PB::Color::BrightnessEnvelope
  brightness_envelope(const MirrorFrame &) {
    return BrightV;
  }
  static float brightness_bottom(const MirrorFrame &frame) {
    return frame.color.brightness_bottom;
  }
  static float brightness_top(const MirrorFrame &frame) {
    return frame.color.brightness_top;
  }
  static float opacity_low(const MirrorFrame &frame) {
    return frame.color.opacity_low;
  }
  static float opacity_high(const MirrorFrame &frame) {
    return frame.color.opacity_high;
  }
};

template <typename WeightP, typename CoverageP, PB::Color::HueMode HueV,
          PB::Color::BrightnessEnvelope BrightV>
using MirrorPipeline = PB::Pipeline<
    MirrorBinding, PB::Stage::Rotate<MirrorCamera>,
    PB::Stage::Project<PB::Projection::Stereographic<MirrorProjection>>,
    PB::Stage::Sample<PB::Source::Grid<MirrorSource>, WeightP, CoverageP>,
    PB::Stage::Colorize<
        PB::Color::GeneratedPalette<MirrorColor<HueV, BrightV>>>>;

/** Mirror frame populated from the erased program's blocks: providers read
    the same raw state the erased prepare reads. @p color is entry 3's family
    in the v3 layout. */
inline MirrorFrame mirror_from(In::ChainProgram &program,
                               const In::FrameContext &ctx,
                               const In::Op::GeneratedPaletteParams &color) {
  MirrorFrame frame;
  const auto &camera = state_as<In::Op::SpatialWalkState>(program, 0);
  const auto &projection = state_as<In::Op::SpatialWalkState>(program, 1);
  const auto &source = state_as<In::Op::SourceClockState>(program, 2);
  const auto &clock = state_as<In::Op::ColorClockState>(program, 3);
  frame.camera_conjugate =
      (math::make_rotation(math::Y_AXIS, camera.spin_phase) * camera.wander)
          .conjugate();
  frame.projection_conjugate =
      (math::make_rotation(math::Y_AXIS, projection.spin_phase) *
       ctx.projection_base * projection.wander)
          .conjugate();
  frame.projection = param_as<In::Op::ProjectChainParams>(program, 1);
  frame.sample = param_as<In::Op::GridSampleParams>(program, 2);
  frame.color = color;
  frame.source_primary = source.primary;
  frame.source_secondary = source.secondary;
  frame.source_angle = source.angle;
  frame.oscillation_phase = clock.oscillation_phase;
  frame.palette = ctx.palettes[frame.color.palette_mode];
  frame.hue_rotation_lut = ctx.hue_rotation_lut;
  frame.hue_noise_lut = ctx.hue_noise_lut;
  return frame;
}

inline MirrorFrame mirror_from(In::ChainProgram &program,
                               const In::FrameContext &ctx) {
  return mirror_from(program, ctx,
                     param_as<In::Op::GeneratedPaletteParams>(program, 3));
}

template <typename Pipe>
void expect_frame_parity(In::ChainProgram &program, const In::FrameContext &ctx,
                         const MirrorFrame &mirror) {
  const typename Pipe::Frame reference_frame = Pipe::prepare(mirror);
  int view_index = 0;
  for (const math::Vector &view : sweep_views()) {
    HS_CONTEXT("view", view_index++);
    const Color4 erased = program.evaluate(view, ctx);
    const Color4 reference = Pipe::shade(view, reference_frame);
    HS_EXPECT_TRUE(color4_identical(erased, reference));
  }
}

template <typename WeightP, typename CoverageP, PB::Color::HueMode HueV,
          PB::Color::BrightnessEnvelope BrightV>
void expect_parity(In::ChainProgram &program, const In::FrameContext &ctx) {
  program.prepare(ctx);
  expect_frame_parity<MirrorPipeline<WeightP, CoverageP, HueV, BrightV>>(
      program, ctx, mirror_from(program, ctx));
}

/** Compiles the default chain and steps it @p frames times. */
inline void arm_default_chain(In::ChainProgram &program, int frames,
                              ValueSet set) {
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(DEFAULT_CHAIN));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 1), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 3), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

/** @brief Runs the erased camera and projection entries (0, 1) over the seed. */
inline PB::PlaneSample projected_input(In::ChainProgram &program,
                                       const In::FrameContext &ctx,
                                       const math::Vector &view) {
  const PB::SphereSample seed{view, 0.0f};
  alignas(In::SLOT_ALIGN) uint8_t rotated[In::SLOT_SIZE];
  program.ops()[0].op->runtime.run(&seed, rotated, ctx, program.param_block(0),
                                   program.prepared_block(0));
  alignas(In::SLOT_ALIGN) uint8_t projected[In::SLOT_SIZE];
  program.ops()[1].op->runtime.run(rotated, projected, ctx,
                                   program.param_block(1),
                                   program.prepared_block(1));
  return *std::launder(reinterpret_cast<PB::PlaneSample *>(projected));
}
