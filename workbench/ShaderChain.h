/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <span>
#include <string_view>

/**
 * @file ShaderChain.h
 * @brief ShaderChain: the chain-interpreter preview effect. Renders an
 *        arbitrary compiled operator chain and exposes every chain parameter
 *        as "{instance}.{field-id}".
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "core/render/pullback/interpreter.h"
#include "core/render/pullback/runtime_seeds.h"

namespace hs_test {
namespace shader_chain_tests {
struct ShaderChainWhiteBox;
} // namespace shader_chain_tests
} // namespace hs_test

/**
 * @brief Stage-program interpreter effect over the pullback operator table.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Owns the two chain arenas, the shared color resources the
 * FrameContext borrows (three generated-palette cyclers, hue LUTs, gamut LUT),
 * and the dynamic parameter schema. Presets are the document layer's concern;
 * the effect registers chain parameters only.
 */
template <int W, int H> class ShaderChain : public Effect {
public:
  static constexpr std::string_view EFFECT_ID = "shader-chain";

  using ChainEntryRequest = Pullback::Interp::ChainEntryRequest;
  using ChainRefusal = Pullback::Interp::ChainRefusal;
  using ChainStatus = Pullback::Interp::ChainStatus;

  HS_COLD_MEMBER ShaderChain() : Effect(W, H, {.strobe = true}) {}

  /** @brief Allocates chain storage and color resources, compiles the default
      chain, and registers its parameters. */
  HS_COLD_MEMBER void init() override {
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    resources = persistent_arena.make<Resources>();
    resources->hue_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    resources->hue_noise.SetSeed(Pullback::HUE_NOISE_SEED);
    resources->hue_noise.SetFrequency(1.0f);
    auto *block_a = static_cast<uint8_t *>(persistent_arena.allocate(
        Pullback::Interp::CHAIN_ARENA_BYTES, alignof(std::max_align_t)));
    auto *block_b = static_cast<uint8_t *>(persistent_arena.allocate(
        Pullback::Interp::CHAIN_ARENA_BYTES, alignof(std::max_align_t)));
    program.bind_storage(block_a, block_b);

    init_gamut_lut(persistent_arena, GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS);
    generated_palettes.init(persistent_arena, DEFAULT_CHROMA, ease_in_out_sin);

    static constexpr ChainEntryRequest DEFAULT_CHAIN[] = {
        {"camera", "sphere.rotate.v2"},
        {"project", "project.stereographic.v2"},
        {"sample", "sample.grid.v2"},
        {"colorize", "colorize.generated-palette.v2"},
    };
    const ChainRefusal refusal =
        set_chain(std::span<const ChainEntryRequest>(DEFAULT_CHAIN));
    HS_CHECK(refusal.code == ChainStatus::OK,
             "ShaderChain: the default chain must compile");
  }

  /**
   * @brief Compiles a program shape transactionally.
   * @return The compile refusal; {OK, -1} on commit.
   * @details On commit the parameter definitions are rebuilt and the schema
   * generation bumped BEFORE returning, so preset values never apply against a
   * stale definition snapshot. A refusal leaves the previous program, its
   * registered definitions, and all live instance state untouched.
   */
  HS_COLD_MEMBER ChainRefusal
  set_chain(std::span<const ChainEntryRequest> request) {
    const ChainRefusal refusal = program.compile(request);
    if (refusal.code == ChainStatus::OK)
      rebind_chain_parameters();
    return refusal;
  }

  /** @brief Advances clocks and palettes, prepares the program, and renders
      one frame. */
  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    ++frame_index;
    program.advance();
    const ColorizeTap tap = colorize_tap();
    update_palette_chroma(tap.params != nullptr ? tap.params->palette_chroma
                                                : DEFAULT_CHROMA);
    step_generated_palettes(tap.params != nullptr ? tap.params->palette_mode
                                                  : uint8_t{0});
    const Pullback::Interp::FrameContext ctx = make_frame_context(tap);
    program.prepare(ctx);
    program.check_ready();
    const ChainShader shader{&program, &ctx};
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

private:
  friend struct ::hs_test::shader_chain_tests::ShaderChainWhiteBox;

  /** @brief Per-sample functor over the compiled program. */
  struct ChainShader {
    const Pullback::Interp::ChainProgram *program;
    const Pullback::Interp::FrameContext *ctx;

    HS_FLASH_MEMBER Color4 operator()(const Vector &view) const {
      return program->evaluate(view, *ctx);
    }
  };

  /** @brief Engine-owned shared resources the FrameContext borrows. */
  struct Resources {
    std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE>
        hue_rotation_lut{};
    std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut{};
    FastNoiseLite hue_noise;
    float baked_noise_scale = 0.0f; /**< 0 forces the first active bake. */
    float baked_noise_phase = 0.0f;
  };

  /** @brief The committed program's colorize entry, when present. */
  struct ColorizeTap {
    int index = -1;
    const Pullback::Interp::Op::GeneratedPaletteParams *params = nullptr;
  };

  ColorizeTap colorize_tap() {
    const auto ops = program.ops();
    for (size_t index = ops.size(); index-- > 0;)
      if (std::string_view(ops[index].op->operator_id) ==
          Pullback::Interp::Op::ColorizeGeneratedPalette::ID)
        return {static_cast<int>(index),
                reinterpret_cast<
                    const Pullback::Interp::Op::GeneratedPaletteParams *>(
                    program.param_block(index))};
    return {};
  }

  /**
   * @brief Rebuilds the registered parameter definitions from the committed
   *        program.
   * @details Name strings live in the winning arena (ChainProgram::
   * param_name), so this must run immediately after every successful compile:
   * the compile that produced the program reset the loser arena holding the
   * previously registered names.
   */
  HS_COLD_MEMBER void rebind_chain_parameters() {
    reset_parameters();
    const auto ops = program.ops();
    for (size_t index = 0; index < ops.size(); ++index) {
      const Pullback::Interp::OperatorDescriptor &op = *ops[index].op;
      uint8_t *block = program.param_block(index);
      for (uint16_t field = 0; field < op.schema_count; ++field) {
        const Pullback::Interp::ParamFieldInfo &info = op.schema[field];
        const char *name = program.param_name(index, field);
        void *address = op.runtime.param_address(block, field);
        if (info.enum_count > 0)
          register_animated_enum8_param(name, static_cast<uint8_t *>(address),
                                        info.enum_ids, info.enum_count);
        else
          register_animated_param(name, static_cast<float *>(address), info.min,
                                  info.max);
      }
    }
  }

  /** @brief Builds the per-frame snapshot, baking hue LUTs when active. */
  HS_FLASH_MEMBER Pullback::Interp::FrameContext
  make_frame_context(const ColorizeTap &tap) {
    Pullback::Interp::FrameContext ctx;
    ctx.frame = frame_index;
    ctx.time = static_cast<float>(frame_index) * FRAME_SECONDS;
    ctx.anims_paused = anims_paused;
    ctx.projection_base = make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
    using PaletteMode = Pullback::Interp::Op::PaletteMode;
    ctx.palettes = {&generated_palettes.palette(PaletteMode::TRIADIC),
                    &generated_palettes.palette(PaletteMode::COMPLEMENTARY),
                    &generated_palettes.palette(PaletteMode::ANALOGOUS)};
    if (tap.params == nullptr || tap.params->hue_shift_amount == 0.0f)
      return ctx;
    const size_t palette_index = tap.params->palette_mode < ctx.palettes.size()
                                     ? tap.params->palette_mode
                                     : 0;
    Pullback::Color::prepare_hue_rotation_lut(
        std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
            resources->hue_rotation_lut),
        *ctx.palettes[palette_index]);
    ctx.hue_rotation_lut = resources->hue_rotation_lut.data();
    if (static_cast<Pullback::Interp::Op::HueShiftMode>(tap.params->hue_mode) !=
        Pullback::Interp::Op::HueShiftMode::NOISE)
      return ctx;
    const auto &clocks =
        *static_cast<const Pullback::Interp::Op::ColorClockState *>(
            program.state_block(static_cast<size_t>(tap.index)));
    if (resources->baked_noise_scale != tap.params->hue_noise_scale ||
        resources->baked_noise_phase != clocks.hue_noise_phase) {
      Pullback::Color::prepare_hue_noise_lut(
          std::span<int8_t, Pullback::Color::HueNoiseLutView::SIZE>(
              resources->hue_noise_lut),
          resources->hue_noise, tap.params->hue_noise_scale,
          clocks.hue_noise_phase);
      resources->baked_noise_scale = tap.params->hue_noise_scale;
      resources->baked_noise_phase = clocks.hue_noise_phase;
    }
    ctx.hue_noise_lut = resources->hue_noise_lut.data();
    return ctx;
  }

  void step_generated_palettes(uint8_t visible_mode) {
    generated_palettes.step(
        static_cast<Pullback::Interp::Op::PaletteMode>(visible_mode));
  }

  HS_COLD_MEMBER void update_palette_chroma(float chroma) {
    generated_palettes.set_chroma(chroma);
  }

  static constexpr float DEFAULT_CHROMA =
      Pullback::Color::ColorParams{}.palette_chroma;
  /** Chain schema capacity; the effect registers no globals of its own. */
  static constexpr size_t PARAM_CAPACITY = Pullback::Interp::MAX_CHAIN_PARAMS;
  /** Nominal frame period behind FrameContext::time. The engine has no wall
      clock — every operator rate is per frame — so this only gives the context
      a monotonic seconds axis. */
  static constexpr float FRAME_SECONDS = 1.0f / 30.0f;

  Pullback::Interp::ChainProgram program;
  Resources *resources = nullptr;
  EffectPaletteRecipes::GeneratedPaletteBank generated_palettes;
  uint32_t frame_index = 0;
};

#include "core/control/registry.h"
REGISTER_EFFECT(ShaderChain)

#endif // HS_ENABLE_CHAIN_INTERPRETER
