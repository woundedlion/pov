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
 * @file chain_host.h
 * @brief ShaderChain: the chain-interpreter preview effect. Renders an
 *        arbitrary compiled operator chain and exposes every chain parameter
 *        as "{instance}.{field-id}".
 */

#include "core/color/palette_cycler.h"
#include "core/engine/engine.h"
#include "core/render/pullback/interpreter.h"
#include "core/render/pullback/runtime_seeds.h"
#include "chain_snapshot.h"

namespace hs_test {
namespace shader_chain_tests {
struct ShaderChainWhiteBox;
} // namespace shader_chain_tests
} // namespace hs_test

/** @brief One named value in an atomic chain parameter update. */
struct ShaderChainParameterWrite {
  const char *name;
  float value;
};

/**
 * @brief Stage-program interpreter effect over the pullback operator table.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Owns the two chain arenas, the shared color resources the
 * FrameContext borrows (three generated-palette cyclers and hue LUTs),
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
    use_parameter_storage(persistent_arena,
                          persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    resources = persistent_arena.make<Resources>();
    resources->hue_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    resources->hue_noise.SetSeed(Pullback::HUE_NOISE_SEED);
    resources->hue_noise.SetFrequency(1.0f);
    static_assert(
        Pullback::Interp::CHAIN_ARENA_BYTES % sizeof(std::max_align_t) == 0);
    auto *block_a = reinterpret_cast<uint8_t *>(
        persistent_arena.allocate_n<std::max_align_t>(
            Pullback::Interp::CHAIN_ARENA_BYTES / sizeof(std::max_align_t)));
    auto *block_b = reinterpret_cast<uint8_t *>(
        persistent_arena.allocate_n<std::max_align_t>(
            Pullback::Interp::CHAIN_ARENA_BYTES / sizeof(std::max_align_t)));
    program.bind_storage(block_a, block_b);

    generated_palettes.init(persistent_arena, DEFAULT_CHROMA,
                            math::ease_in_out_sin);

    static constexpr ChainEntryRequest DEFAULT_CHAIN[] = {
        {"camera", "sphere.rotate.v2"},
        {"project", "project.stereographic.v2"},
        {"sample", "sample.grid.v3"},
        {"colorize", "colorize.generated-palette.v3"},
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
    if (refusal.code == ChainStatus::OK) {
#if HS_ENABLE_PARAM_GUI_BRIDGE
      refused_name = nullptr;
      refusal_warning = nullptr;
#endif
      colorize = find_colorize_tap();
      rebind_chain_parameters();
    }
    return refusal;
  }

  /** @brief Borrows the compiled entries until the next program replacement. */
  std::span<const Pullback::Interp::ChainProgram::ChainOp> chain_ops() const {
    return program.ops();
  }

  ChainSnapshot snapshot() const {
    ChainSnapshot out;
    out.animations_paused = animations_paused();
    out.palette_bank = generated_palettes.snapshot();
    out.runtime.emplace();
    const auto ops = program.ops();
    for (size_t index = 0; index < ops.size(); ++index) {
      const auto &op = ops[index];
      out.chain.push_back({op.instance, op.op->operator_id});
      for (uint16_t field = 0; field < op.op->schema_count; ++field) {
        const auto &info = op.op->schema[field];
        const void *address = op.op->runtime.param_address(
            const_cast<uint8_t *>(program.param_block(index)), field);
        const float value =
            info.enum_count > 0
                ? static_cast<float>(*static_cast<const uint8_t *>(address))
                : *static_cast<const float *>(address);
        out.parameters.push_back({program.param_name(index, field), value});
      }
      auto state = op.op->runtime.capture_state(program.state_block(index));
      if (!std::holds_alternative<std::monostate>(state))
        out.runtime->push_back({op.instance, std::move(state)});
    }
    return out;
  }

  HS_COLD_MEMBER ChainSnapshotRestoreResult
  restore_snapshot(const ChainSnapshot &snapshot) {
    using Result = ChainSnapshotRestoreResult;
    using namespace Pullback::Interp;
    if (snapshot.schema_version != ChainSnapshot::SCHEMA_VERSION)
      return Result::UNSUPPORTED_VERSION;
    if (snapshot.chain.empty() || snapshot.chain.size() > MAX_CHAIN_OPS ||
        snapshot.parameters.size() > MAX_CHAIN_PARAMS ||
        (snapshot.runtime && snapshot.runtime->size() > MAX_CHAIN_OPS))
      return Result::INVALID_LENGTH;
    if (snapshot.palette_bank &&
        !GeneratedPaletteBank::valid_snapshot(*snapshot.palette_bank))
      return Result::INVALID_VALUE;
    ChainEntryRequest requests[MAX_CHAIN_OPS];
    for (size_t index = 0; index < snapshot.chain.size(); ++index)
      requests[index] = {snapshot.chain[index].instance,
                         snapshot.chain[index].operator_id};
    Result validation = Result::INVALID_VALUE;
    const auto refusal = program.compile(
        std::span<const ChainEntryRequest>(requests, snapshot.chain.size()),
        [&](std::span<const ChainProgram::ChainOp> ops, uint8_t *base) {
          for (size_t write_index = 0; write_index < snapshot.parameters.size();
               ++write_index) {
            const auto &write = snapshot.parameters[write_index];
            if (!std::isfinite(write.value))
              return false;
            for (size_t earlier = 0; earlier < write_index; ++earlier)
              if (snapshot.parameters[earlier].name == write.name)
                return false;
            bool found = false;
            for (const auto &op : ops)
              for (uint16_t field = 0; field < op.op->schema_count; ++field) {
                const auto &info = op.op->schema[field];
                if (write.name != std::string(op.instance) + "." + info.id)
                  continue;
                if (write.value < (info.enum_count > 0 ? 0.0f : info.min) ||
                    write.value > (info.enum_count > 0
                                       ? static_cast<float>(info.enum_count - 1)
                                       : info.max))
                  return false;
                void *address =
                    op.op->runtime.param_address(base + op.param_offset, field);
                if (info.enum_count > 0) {
                  if (std::floor(write.value) != write.value)
                    return false;
                  *static_cast<uint8_t *>(address) =
                      static_cast<uint8_t>(write.value);
                } else {
                  *static_cast<float *>(address) = write.value;
                }
                found = true;
              }
            if (!found)
              return false;
          }
          bool edge_distance_available = false;
          size_t runtime_count = 0;
          for (const auto &op : ops) {
            void *params = base + op.param_offset;
            if (op.op->runtime.validate(params))
              return false;
            if (op.op->input == CarrierId::SPHERE &&
                op.op->output == CarrierId::PLANE)
              edge_distance_available = op.op->edge_distance_available;
            if (!edge_distance_available)
              for (uint16_t field = 0; field < op.op->schema_count; ++field) {
                const auto &info = op.op->schema[field];
                if (info.enum_count == 0 ||
                    (std::strcmp(info.id, "coverage-mode") != 0 &&
                     std::strcmp(info.id, "envelope") != 0))
                  continue;
                const auto value = *static_cast<uint8_t *>(
                    op.op->runtime.param_address(params, field));
                if (value < info.enum_count &&
                    std::strcmp(info.enum_ids[value], "edge-fade") == 0)
                  return false;
              }
            if (!snapshot.runtime)
              continue;
            void *state = base + op.state_offset;
            const auto initial = op.op->runtime.capture_state(state);
            if (std::holds_alternative<std::monostate>(initial))
              continue;
            const ChainSnapshot::Runtime *selected = nullptr;
            for (const auto &entry : *snapshot.runtime)
              if (entry.instance == op.instance) {
                if (selected != nullptr)
                  return false;
                selected = &entry;
              }
            if (selected == nullptr ||
                !op.op->runtime.restore_state(state, selected->state))
              return false;
            ++runtime_count;
          }
          return !snapshot.runtime || runtime_count == snapshot.runtime->size();
        },
        false);
    if (refusal.code != ChainStatus::OK)
      return refusal.code == ChainStatus::MALFORMED_PAYLOAD
                 ? validation
                 : Result::INVALID_CHAIN;
    colorize = find_colorize_tap();
    rebind_chain_parameters();
    GeneratedPaletteBank::Snapshot palette;
    palette.chroma =
        colorize.palette_chroma ? *colorize.palette_chroma : DEFAULT_CHROMA;
    palette.hues.fill(GeneratedPaletteBank::HUE_STEP);
    generated_palettes.restore_snapshot(
        snapshot.palette_bank.value_or(palette));
    resources->hue_noise_bake = {};
    setAnimationsPaused(snapshot.animations_paused);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    refused_name = nullptr;
    refusal_warning = nullptr;
#endif
    return Result::APPLIED;
  }

  /** @brief Advances clocks and palettes, prepares the program, and renders
      one frame. */
  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    program.advance();
    const ColorizeTap &tap = colorize;
    update_palette_chroma(tap.palette_chroma != nullptr ? *tap.palette_chroma
                                                        : DEFAULT_CHROMA);
    step_generated_palettes(tap.palette_mode != nullptr ? *tap.palette_mode
                                                        : uint8_t{0});
    const Pullback::Interp::FrameContext ctx = make_frame_context(tap);
    program.prepare(ctx);
    program.check_ready();
    const ChainShader shader{&program, &ctx};
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Validates the final parameter state, then commits every write. */
  ParamSetResult
  update_parameters(std::span<const ShaderChainParameterWrite> writes) {
    if (writes.size() > Pullback::Interp::MAX_CHAIN_PARAMS)
      return ParamSetResult::TOO_LONG;
    alignas(std::max_align_t)
        uint8_t candidates[Pullback::Interp::MAX_CHAIN_OPS][PARAM_BYTES];
    const auto ops = program.ops();
    for (size_t index = 0; index < ops.size(); ++index) {
      const auto &runtime = ops[index].op->runtime;
      runtime.construct_params(candidates[index]);
      std::memcpy(candidates[index], program.param_block(index),
                  runtime.param.size);
    }
    bool animated = false;
    for (const auto &write : writes) {
      if (write.name == nullptr)
        return ParamSetResult::MALFORMED_PAYLOAD;
      const ParamDef *parameter = getParameters().find(write.name);
      if (parameter == nullptr)
        return ParamSetResult::UNKNOWN_PARAM;
      float value = write.value;
      const ParamSetResult result = parameter->normalize(value);
      if (result != ParamSetResult::APPLIED)
        return result;
      for (size_t index = 0; index < ops.size(); ++index)
        for (uint16_t field = 0; field < ops[index].op->schema_count; ++field)
          if (std::strcmp(write.name, program.param_name(index, field)) == 0) {
            ParamDef proposed = *parameter;
            proposed.target =
                ops[index].op->runtime.param_address(candidates[index], field);
            write_parameter_unchecked(proposed, value);
          }
      animated |= parameter->animated;
    }
    for (size_t index = 0; index < ops.size(); ++index)
      if (const char *warning = validate_parameters(index, candidates[index])) {
        refused_name = nullptr;
        for (const auto &write : writes) {
          for (uint16_t field = 0; field < ops[index].op->schema_count; ++field)
            if (std::strcmp(write.name, program.param_name(index, field)) ==
                0) {
              refused_name = program.param_name(index, field);
              break;
            }
          if (refused_name != nullptr)
            break;
        }
        refusal_warning = warning;
        return ParamSetResult::INADMISSIBLE;
      }
    for (size_t index = 0; index < ops.size(); ++index)
      std::memcpy(program.param_block(index), candidates[index],
                  ops[index].op->runtime.param.size);
    refused_name = nullptr;
    refusal_warning = nullptr;
    if (animated)
      setAnimationsPaused(true);
    parameter_written();
    return ParamSetResult::APPLIED;
  }

  const char *parameter_warning(const char *name) const override {
    return refused_name != nullptr && std::strcmp(name, refused_name) == 0
               ? refusal_warning
               : nullptr;
  }
#endif

private:
  friend struct ::hs_test::shader_chain_tests::ShaderChainWhiteBox;

  /** @brief Per-sample functor over the compiled program. */
  struct ChainShader {
    const Pullback::Interp::ChainProgram *program;
    const Pullback::Interp::FrameContext *ctx;

    HS_FLASH_MEMBER Color4 operator()(const math::Vector &view) const {
      return program->evaluate(view, *ctx);
    }
  };

  /** @brief Engine-owned shared resources the FrameContext borrows. */
  struct Resources {
    int32_t hue_noise_seed = Pullback::HUE_NOISE_SEED;
    std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE>
        hue_rotation_lut{};
    std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut{};
    FastNoiseLite hue_noise;
    Pullback::Color::HueNoiseBakeCache hue_noise_bake;
  };

  /** @brief The committed program's colorize entry, when present. */
  struct ColorizeTap {
    int index = -1;
    const float *palette_chroma = nullptr;
    const float *hue_shift_amount = nullptr;
    const float *hue_noise_scale = nullptr;
    const uint8_t *palette_mode = nullptr;
    const uint8_t *hue_mode = nullptr;
  };

  template <typename Params>
  static ColorizeTap make_colorize_tap(int index, const Params *params) {
    return {index,
            &params->palette_chroma,
            &params->hue_shift_amount,
            &params->hue_noise_scale,
            &params->palette_mode,
            &params->hue_mode};
  }

  HS_COLD_MEMBER ColorizeTap find_colorize_tap() {
    const auto ops = program.ops();
    for (size_t index = ops.size(); index-- > 0;)
      if (std::string_view(ops[index].op->operator_id) ==
          Pullback::Interp::Op::ColorizeGeneratedPaletteV3::ID)
        return make_colorize_tap(
            static_cast<int>(index),
            reinterpret_cast<
                const Pullback::Interp::Op::ColorizeGeneratedPaletteV3::Params
                    *>(program.param_block(index)));
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

#if HS_ENABLE_PARAM_GUI_BRIDGE
  const char *validate_parameters(size_t index, void *params) {
    const auto ops = program.ops();
    const auto &op = *ops[index].op;
    if (const char *warning = op.runtime.validate(params))
      return warning;
    bool edge_distance_available = false;
    for (size_t upstream = 0; upstream < index; ++upstream)
      if (ops[upstream].op->input == Pullback::Interp::CarrierId::SPHERE &&
          ops[upstream].op->output == Pullback::Interp::CarrierId::PLANE)
        edge_distance_available = ops[upstream].op->edge_distance_available;
    if (!edge_distance_available)
      for (uint16_t field = 0; field < op.schema_count; ++field) {
        const auto &info = op.schema[field];
        if (info.enum_count == 0 ||
            (std::strcmp(info.id, "coverage-mode") != 0 &&
             std::strcmp(info.id, "envelope") != 0))
          continue;
        const auto value =
            *static_cast<uint8_t *>(op.runtime.param_address(params, field));
        if (value < info.enum_count &&
            std::strcmp(info.enum_ids[value], "edge-fade") == 0)
          return "Edge-fade requires a projection with edge distance";
      }
    return nullptr;
  }

  static constexpr size_t PARAM_BYTES = [] {
    size_t largest = 0;
    for (const auto &op : Pullback::Interp::OPERATOR_TABLE)
      largest = std::max(largest, static_cast<size_t>(op.runtime.param.size));
    return (largest + alignof(std::max_align_t) - 1) /
           alignof(std::max_align_t) * alignof(std::max_align_t);
  }();

  bool parameter_write_admitted(const ParamDef &parameter,
                                float value) override {
    const auto ops = program.ops();
    for (size_t index = 0; index < ops.size(); ++index)
      for (uint16_t field = 0; field < ops[index].op->schema_count; ++field)
        if (std::strcmp(parameter.name, program.param_name(index, field)) ==
            0) {
          const auto &runtime = ops[index].op->runtime;
          alignas(std::max_align_t) uint8_t candidate[PARAM_BYTES];
          runtime.construct_params(candidate);
          std::memcpy(candidate, program.param_block(index),
                      runtime.param.size);
          ParamDef proposed = parameter;
          proposed.target = runtime.param_address(candidate, field);
          write_parameter_unchecked(proposed, value);
          refusal_warning = validate_parameters(index, candidate);
          refused_name = refusal_warning != nullptr ? parameter.name : nullptr;
          return refusal_warning == nullptr;
        }
    return true;
  }

  const char *refused_name = nullptr;
  const char *refusal_warning = nullptr;
#endif

  /** @brief Builds the per-frame snapshot, baking hue LUTs when active. */
  HS_FLASH_MEMBER Pullback::Interp::FrameContext
  make_frame_context(const ColorizeTap &tap) {
    Pullback::Interp::FrameContext ctx;
    ctx.projection_base = Pullback::projection_base_orientation();
    using PaletteMode = Pullback::Interp::Op::PaletteMode;
    ctx.palettes = {&generated_palettes.palette(PaletteMode::TRIADIC),
                    &generated_palettes.palette(PaletteMode::COMPLEMENTARY),
                    &generated_palettes.palette(PaletteMode::ANALOGOUS)};
    if (tap.hue_shift_amount == nullptr || *tap.hue_shift_amount == 0.0f)
      return ctx;
    // The colorize stage reads neither table under HueShiftMode::NONE, so a
    // leftover amount must not buy a palette resample.
    const auto hue_mode =
        static_cast<Pullback::Interp::Op::HueShiftMode>(*tap.hue_mode);
    if (hue_mode == Pullback::Interp::Op::HueShiftMode::NONE)
      return ctx;
    const size_t palette_index =
        *tap.palette_mode < ctx.palettes.size() ? *tap.palette_mode : 0;
    Pullback::Color::prepare_hue_rotation_lut(
        std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
            resources->hue_rotation_lut),
        *ctx.palettes[palette_index]);
    ctx.hue_rotation_lut = resources->hue_rotation_lut.data();
    if (hue_mode != Pullback::Interp::Op::HueShiftMode::NOISE)
      return ctx;
    const auto &clocks =
        *static_cast<const Pullback::Interp::Op::ColorClockState *>(
            program.state_block(static_cast<size_t>(tap.index)));
    if (resources->hue_noise_seed != clocks.hue_noise_seed) {
      resources->hue_noise.SetSeed(clocks.hue_noise_seed);
      resources->hue_noise_seed = clocks.hue_noise_seed;
      resources->hue_noise_bake = {};
    }
    resources->hue_noise_bake.refresh(
        resources->hue_noise_lut, resources->hue_noise, *tap.hue_noise_scale,
        clocks.hue_noise_phase);
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

  Pullback::Interp::ChainProgram program;
  /** Refreshed on every commit; the param block it points at stays put until
      the next compile. */
  ColorizeTap colorize;
  Resources *resources = nullptr;
  GeneratedPaletteBank generated_palettes;

  // Against the browser module's arena, not the build's: this effect never
  // reaches the device, and the host suite's arena is far too loose to catch a
  // widened chain arena before it traps at construction in the browser.
  static constexpr size_t FOOTPRINT_BYTES =
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(Resources) +
      alignof(Resources) +
      2 * (Pullback::Interp::CHAIN_ARENA_BYTES + alignof(std::max_align_t)) +
      GeneratedPaletteBank::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= WASM_PERSISTENT_BUDGET,
                "ShaderChain persistent footprint exceeds the browser "
                "module's default partition");
};

#include "core/control/registry.h"
#endif // HS_ENABLE_CHAIN_INTERPRETER
