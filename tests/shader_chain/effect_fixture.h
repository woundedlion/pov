/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_shader_chain.h.

// ShaderChain effect white-box fixture.

/** Reaches the effect's committed program, palette state, and hue bakes. */
struct ShaderChainWhiteBox {
  using FX = ShaderChain<96, 20>;

  static const In::ChainProgram &program(const FX &effect) {
    return effect.program;
  }
  static In::ChainProgram &program(FX &effect) { return effect.program; }
  static void advance_without_render(FX &effect) { effect.advance_clocks(); }
  static Pixel palette_color(const FX &effect, float value) {
    return effect.generated_palettes.palette(In::Op::PaletteMode::TRIADIC)
        .get(value)
        .color;
  }

  /** The committed colorize instance's parameter block. */
  static In::Op::GeneratedPaletteParams &color_params(FX &effect) {
    return *reinterpret_cast<In::Op::GeneratedPaletteParams *>(
        effect.program.param_block(static_cast<size_t>(effect.colorize.index)));
  }
  /** The committed colorize instance's phase clocks. */
  static const In::Op::ColorClockState &color_clocks(const FX &effect) {
    return *static_cast<const In::Op::ColorClockState *>(
        effect.program.state_block(static_cast<size_t>(effect.colorize.index)));
  }
  /** The frame snapshot draw_frame() hands the program, hue bakes included. */
  static In::FrameContext frame_context(FX &effect) {
    return effect.make_frame_context(effect.colorize);
  }
  /** The resident hue-noise table; writable so a probe can poison it and see
      whether the next frame re-baked. */
  static int8_t *hue_noise_lut(FX &effect) {
    return effect.resources->hue_noise_lut.data();
  }
  /** The inputs the resident hue-noise table was baked from. */
  static float baked_noise_scale(const FX &effect) {
    return effect.resources->hue_noise_bake.scale;
  }
  static float baked_noise_phase(const FX &effect) {
    return effect.resources->hue_noise_bake.phase;
  }
};
