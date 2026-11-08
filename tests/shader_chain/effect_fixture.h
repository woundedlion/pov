/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ShaderChain effect white-box fixture.

/** Reaches the effect's committed program, palette state, and hue bakes. */
struct ShaderChainWhiteBox {
  using FX = ShaderChain<96, 20>;

  static const In::ChainProgram &program(const FX &effect) {
    return effect.program;
  }
  static In::ChainProgram &program(FX &effect) { return effect.program; }
  static int colorize_index(FX &effect) {
    return effect.find_colorize_tap().index;
  }
  static void advance_without_render(FX &effect) { effect.advance_clocks(); }
  /** The generated palette bank's current state. */
  static GeneratedPaletteBank::Snapshot palette_bank(const FX &effect) {
    return effect.generated_palettes.snapshot();
  }
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
  static size_t modulation_index(const FX &effect,
                                 std::string_view label = "modulation") {
    for (size_t index = 0; index < effect.program.ops().size(); ++index)
      if (std::string_view(effect.program.ops()[index].instance) == label)
        return index;
    HS_CHECK(false, "missing test modulation node");
    return 0;
  }
  /** @brief Appends a named modulation node if absent. */
  static void ensure_modulation(FX &effect, const char *label = "modulation") {
    for (const auto &op : effect.program.ops())
      if (std::string_view(op.instance) == label)
        return;
    std::vector<In::ChainEntryRequest> chain;
    for (const auto &op : effect.program.ops())
      chain.push_back({op.instance, op.op->operator_id});
    chain.push_back({label, In::Op::ModulateColorV1::ID});
    HS_EXPECT_EQ(effect.set_chain(chain).code, In::ChainStatus::OK);
  }
  /** @brief Existing node's parameter block; traps when absent. */
  static In::Op::ModulationParams &
  modulation_params(FX &effect, const char *label = "modulation") {
    return *reinterpret_cast<In::Op::ModulationParams *>(
        effect.program.param_block(modulation_index(effect, label)));
  }
  static const In::Op::ModulationClockState &
  modulation_clocks(const FX &effect, std::string_view label = "modulation") {
    return *static_cast<const In::Op::ModulationClockState *>(
        effect.program.state_block(modulation_index(effect, label)));
  }
  static auto &layer_resources(FX &effect, std::string_view key) {
    return effect.modulation_resources(key);
  }
  /** The retained previous frame context, its palettes and hue tables owned
      by the effect. */
  static const In::FrameContext &previous_context(const FX &effect) {
    return effect.resources->context.previous;
  }
  /** The frame snapshot draw_frame() hands the program, hue bakes included. */
  static In::FrameContext frame_context(FX &effect) {
    return effect.make_frame_context(effect.colorize);
  }
};
