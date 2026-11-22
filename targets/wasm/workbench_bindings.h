/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file workbench_bindings.h
 * @brief Typed WASM authoring adapters bound to an effect incarnation.
 */
#pragma once

#include "targets/wasm/ts_types.h"
#include <emscripten/bind.h>
#include "targets/wasm/effect_factory.h"
#include "targets/wasm/payload_clone.h"
#include "targets/wasm/param_marshal.h"
static_assert(HS_ENABLE_CHAIN_INTERPRETER);
#include "core/render/pullback/catalog_export.h"
#include "targets/wasm/chain_snapshot_codec.h"
#include <functional>
#include <memory>
#include <span>
#include <vector>

// Caller property access can re-enter embind, including delete().
/// True while any caller payload is being decoded.
static bool payload_decode_active = false;
/// RAII latch on `payload_decode_active`; traps on a re-entrant decode.
struct PayloadDecodeGuard {
  bool *decoding; ///< Optional per-handle flag mirrored for the guard's life.
  /**
   * @brief Enters a decode; traps if one is already active.
   * @param decoding Per-handle flag set for the guard's life, or null.
   */
  explicit PayloadDecodeGuard(bool *decoding = nullptr) : decoding(decoding) {
    HS_CHECK(!payload_decode_active,
             "re-entrant engine decode from a caller accessor");
    payload_decode_active = true;
    if (decoding)
      *decoding = true;
  }
  ~PayloadDecodeGuard() {
    if (decoding)
      *decoding = false;
    payload_decode_active = false;
  }
};

/** @brief False once the module aborted or was marked dead. */
// clang-format off
EM_JS(bool, workbench_module_alive, (), {
  return !ABORT && !Module['HS_MODULE_DEAD'];
});
// clang-format on

/** @brief WASM-only effect incarnation shared by engine and capability handles. */
struct WorkbenchBindingState {
  Effect *effect = nullptr; ///< Live effect, or null between effects.
  const FactoryEntry *entry = nullptr; ///< Factory entry `effect` came from.
  uint64_t generation = 0;             ///< Bumped on every effect teardown.
  int width = 0;                       ///< Active resolution width, in pixels.
  int height = 0;                      ///< Active resolution height, in pixels.
  bool alive = true;   ///< False once the owning engine is deleted.
  bool paused = false; ///< Engine-owned animation pause; seeds each new effect.
};

/** @brief Base for capability handles pinned to one effect incarnation. */
class WorkbenchBindings {
public:
  /// Reapplies captured state to a rebuilt effect; true on success.
  using RebuildRestore =
      std::function<bool(const std::shared_ptr<WorkbenchBindingState> &)>;
  /**
   * @brief Pins the handle to the state's current generation.
   * @param state Shared incarnation state owned by the engine.
   */
  explicit WorkbenchBindings(std::shared_ptr<WorkbenchBindingState> state)
      : state(std::move(state)), generation(this->state->generation) {}
  /**
   * @brief Whether the originating engine and effect incarnation remain live.
   * @return True while the handle may be used.
   */
  bool isValid() const {
    return workbench_module_alive() && state->alive && state->effect &&
           state->entry && state->generation == generation;
  }

protected:
  std::shared_ptr<WorkbenchBindingState> state; ///< Engine-shared incarnation.
  /**
   * @brief Runs a callback when the live effect has the requested factory
   *        type.
   * @tparam EffectT Effect class template to match.
   * @tparam Callback Callable taking `EffectT<W, H> &`.
   * @param callback Invoked with the live effect on a type match.
   * @return True when the callback ran.
   */
  template <template <int, int> class EffectT, typename Callback>
  bool with_effect(Callback &&callback) {
    if (!isValid())
      return false;
    bool invoked = false;
    hs_wasm::dispatch_resolution(
        state->width, state->height, [&]<int W, int H>() {
          if (state->entry->type_key == effect_type_key<EffectT<W, H>>()) {
            callback(static_cast<EffectT<W, H> &>(*state->effect));
            invoked = true;
          }
        });
    return invoked;
  }

private:
  const uint64_t generation;
};
/** @brief JS handle for program and snapshot access on a ShaderChain effect. */
class ShaderChainBindings : public WorkbenchBindings {
public:
  /// Inherits the incarnation-pinning constructor.
  /// @param state Shared incarnation state owned by the engine.
  using WorkbenchBindings::WorkbenchBindings;
  ~ShaderChainBindings() {
    HS_CHECK(!decoding,
             "delete() of a chain handle from a caller accessor during decode");
  }
  /**
   * @brief Whether the originating engine and effect incarnation remain live.
   * @return True while the handle may be used.
   */
  bool isValid() const { return WorkbenchBindings::isValid(); }
  /**
   * @brief Encodes the loaded chain's full state.
   * @return Snapshot JS object, or null if the handle is stale or not a chain.
   */
  SnapshotValue getSnapshot() {
    auto result = emscripten::val::null();
    with_effect<ShaderChain>([&](auto &chain) {
      result = hs_wasm::ChainSnapshotCodec::encode(chain.snapshot());
    });
    return result;
  }

  /**
   * @brief Decodes and applies a snapshot from getSnapshot().
   * @param caller_input Snapshot JS object; cloned before decoding.
   * @return APPLIED, NOT_SHADER_CHAIN if the handle is stale or the effect
   *         changed during decode, or the decode/restore refusal.
   */
  ChainSnapshotRestoreResult
  restoreSnapshot(const SnapshotInputValue &caller_input) {
    using Result = ChainSnapshotRestoreResult;
    const PayloadDecodeGuard guard(&decoding);
    if (!isValid())
      return Result::NOT_SHADER_CHAIN;
    const emscripten::val input = clone_payload(caller_input);
    if (!isValid())
      return Result::NOT_SHADER_CHAIN;
    ChainSnapshot snapshot;
    const auto result = hs_wasm::ChainSnapshotCodec::decode(input, snapshot);
    if (result != Result::APPLIED)
      return result;
    Result restored = Result::NOT_SHADER_CHAIN;
    with_effect<ShaderChain>(
        [&](auto &chain) { restored = chain.restore_snapshot(snapshot); });
    if (restored == Result::APPLIED) {
      hs_wasm::check_param_capacity(*state->effect);
      state->paused = state->effect->animations_paused();
    }
    return restored;
  }
  /**
   * @brief Captures the chain's snapshot for reapplication after a rebuild.
   * @return Restorer for the rebuilt effect, or empty if no chain is loaded.
   */
  RebuildRestore capture_rebuild_state() {
    RebuildRestore restore;
    with_effect<ShaderChain>([&](auto &chain) {
      restore = [snapshot = chain.snapshot()](const auto &state) {
        ShaderChainBindings bindings(state);
        bool restored = false;
        bindings.with_effect<ShaderChain>([&](auto &target) {
          restored = target.restore_snapshot(snapshot) ==
                     ChainSnapshotRestoreResult::APPLIED;
        });
        return restored;
      };
    });
    return restore;
  }

  /**
   * @brief Compiles a chain program shape on the loaded ShaderChain effect.
   * @param caller_entries JS array of {instance, operator} string pairs: the
   *        ordered program shape only.
   * @return JS object {code, status, entryIndex}: code is "APPLIED" on commit,
   *         otherwise the refusal name; status is the ChainStatus enum.
   *         entryIndex names the offending entry, -1 for a whole-chain refusal.
   * @details On APPLIED the parameter definitions are rebuilt before this
   * returns, so values can be set by "{instance}.{field-id}" immediately.
   * NOT_CHAIN_EFFECT also covers accessors that swap the loaded effect out
   * mid-decode. A refusal commits nothing, but side effects of caller
   * accessors are not rolled back.
   */
  ChainResultValue setShaderChain(const ChainEntriesValue &caller_entries) {
    const PayloadDecodeGuard decode_guard(&decoding);
    using Pullback::Interp::ChainStatus;
    if (!isValid())
      return chain_result(ChainStatus::NOT_CHAIN_EFFECT, -1);
    // Payload cloning can invoke getters that replace the addressed chain.
    const emscripten::val entries = clone_payload(caller_entries);
    if (!isValid())
      return chain_result(ChainStatus::NOT_CHAIN_EFFECT, -1);
    if (!hs_wasm::is_js_array(entries))
      return chain_result(ChainStatus::MALFORMED_PAYLOAD, -1);
    const size_t count = entries["length"].as<size_t>();
    // The length cap precedes per-entry decode, so an oversized payload is
    // TOO_LONG.
    if (count > Pullback::Interp::MAX_CHAIN_OPS)
      return chain_result(ChainStatus::TOO_LONG, -1);
    std::vector<std::string> instances(count);
    std::vector<std::string> operators(count);
    for (size_t index = 0; index < count; ++index) {
      const emscripten::val entry = entries[index];
      if (entry.isNull() || entry.isUndefined())
        return chain_result(ChainStatus::MALFORMED_PAYLOAD,
                            static_cast<int>(index));
      const emscripten::val instance = entry["instance"];
      const emscripten::val operator_id = entry["operator"];
      if (!instance.isString() || !operator_id.isString())
        return chain_result(ChainStatus::MALFORMED_PAYLOAD,
                            static_cast<int>(index));
      instances[index] = instance.as<std::string>();
      operators[index] = operator_id.as<std::string>();
    }
    std::vector<Pullback::Interp::ChainEntryRequest> request(count);
    for (size_t index = 0; index < count; ++index)
      request[index] = {instances[index], operators[index]};
    Pullback::Interp::ChainRefusal refusal{ChainStatus::NOT_CHAIN_EFFECT, -1};
    with_effect<ShaderChain>([&]<typename SC>(SC &chain) {
      refusal = chain.set_chain(
          std::span<const Pullback::Interp::ChainEntryRequest>(request));
    });
    if (refusal.code == ChainStatus::OK)
      hs_wasm::check_param_capacity(*state->effect);
    return chain_result(refusal.code, refusal.entry_index);
  }

  /**
   * @brief Atomically applies named writes after validating the final state.
   * @return APPLIED, or MALFORMED_PAYLOAD for invalid entry shapes, TOO_LONG
   * for capacity overflow, NO_EFFECT for a missing or replaced chain,
   * UNKNOWN_PARAM for an unknown name, READONLY for a protected parameter,
   * NON_FINITE for NaN/infinite values, or INADMISSIBLE for cross-field conflicts.
   * @param caller_entries JS array of `{name, value}` writes; cloned before
   *        decoding.
   */
  ParamSetResult
  setShaderChainParameters(const ParameterWritesValue &caller_entries) {
    return set_parameters(caller_entries, false);
  }

  /** @brief Applies a frame of preview values with animation still running.
   * @param caller_entries Named writes with numeric values or explicit blends.
   * @return APPLIED or the first refusal, with no partial writes. */
  ParamSetResult setShaderChainTransitionParameters(
      const TransitionParameterWritesValue &caller_entries) {
    return set_parameters(caller_entries, true);
  }

  /** @brief Applies an authored planar morph frame without pausing animation.
   * @param caller_entries Named preview parameter writes.
   * @param caller_morph Source endpoint values and the authored edge frame.
   * @return APPLIED or the first refusal, with no partial writes. */
  ParamSetResult setShaderChainSourceMorphParameters(
      const TransitionParameterWritesValue &caller_entries,
      const SourceMorphPreviewValue &caller_morph) {
    return set_parameters(caller_entries, true, &caller_morph);
  }

  /** @brief Applies planar and ordered modulation source transitions atomically.
   * @param caller_entries Named preview writes.
   * @param caller_morph Optional planar source transition.
   * @param caller_modulations Modulation transitions by instance.
   * @return Parameter admission status. */
  ParamSetResult setShaderChainMorphParameters(
      const TransitionParameterWritesValue &caller_entries,
      const NullableSourceMorphPreviewValue &caller_morph,
      const ModulationMorphsPreviewValue &caller_modulations) {
    return set_parameters(caller_entries, true, &caller_morph,
                          &caller_modulations);
  }

  /** @brief Applies a preview frame with a palette crossfade.
   * @param caller_entries Named preview writes.
   * @param caller_morph Optional planar source transition.
   * @param caller_modulations Modulation transitions by instance.
   * @param caller_palette Palette transition endpoints and weight.
   * @return Parameter admission status. */
  ParamSetResult setShaderChainPaletteMorphParameters(
      const TransitionParameterWritesValue &caller_entries,
      const NullableSourceMorphPreviewValue &caller_morph,
      const ModulationMorphsPreviewValue &caller_modulations,
      const PaletteMorphPreviewValue &caller_palette) {
    return set_parameters(caller_entries, true, &caller_morph,
                          &caller_modulations, &caller_palette);
  }

  /** @brief Checks effective float32 source lattices for a preset edge.
   * @param caller_program Ordered entries of the authored chain.
   * @param caller_from Outgoing preset parameters by qualified name.
   * @param caller_to Incoming preset parameters by qualified name.
   * @return True when the edge is compatible with every upstream affine. */
  static bool
  admitsPlanarSourceTransition(const ChainEntriesValue &caller_program,
                               const TransitionEndpointValue &caller_from,
                               const TransitionEndpointValue &caller_to) {
    const PayloadDecodeGuard decode_guard;
    const auto program = clone_payload(caller_program);
    const auto from = clone_payload(caller_from);
    const auto to = clone_payload(caller_to);
    if (!hs_wasm::is_js_array(program) || !hs_wasm::is_js_object(from) ||
        !hs_wasm::is_js_object(to))
      return false;
    const size_t count = program["length"].as<size_t>();
    if (count > Pullback::Interp::MAX_CHAIN_OPS)
      return false;
    namespace Op = Pullback::Interp::Op;
    std::vector<std::string> affines;
    for (size_t index = 0; index < count; ++index) {
      const auto entry = program[index];
      if (!hs_wasm::is_js_object(entry) || !entry["instance"].isString() ||
          !entry["operator"].isString())
        return false;
      const std::string instance = entry["instance"].as<std::string>();
      const std::string op = entry["operator"].as<std::string>();
      if (!Pullback::Interp::instance_id_wellformed(instance))
        return false;
      if (op == Op::WarpAffineV4::ID) {
        affines.push_back(instance);
        continue;
      }
      if (op != Op::SamplePlanar::ID)
        continue;
      Op::PlanarSampleParams source_from, source_to;
      if (!read_endpoint<Op::SamplePlanar>(from, instance, source_from) ||
          !read_endpoint<Op::SamplePlanar>(to, instance, source_to))
        return false;
      if (source_from.source == source_to.source)
        continue;
      for (const std::string &affine : affines) {
        Op::AffineWarpParams affine_from, affine_to;
        if (!read_endpoint<Op::WarpAffineV4>(from, affine, affine_from) ||
            !read_endpoint<Op::WarpAffineV4>(to, affine, affine_to))
          return false;
        if (!Op::admits_planar_transition(source_from, source_to, affine_from,
                                          affine_to))
          return false;
      }
    }
    return true;
  }

  /**
   * @brief Exports the chain-interpreter operator catalog.
   * @return The catalog JSON — budgets, carriers, and every operator-table
   *         entry. Block sizes are wasm32 ABI figures: a pointer-bearing
   *         operator's `prepared` block is narrower than on LP64.
   */
  static std::string getShaderChainCatalog() {
    std::string catalog;
    Pullback::Interp::append_catalog_json(catalog);
    return catalog;
  }

  /**
   * @brief Lists the loaded program shape.
   * @return JS array of {instance, operator} in chain order, or null if the
   *         handle is stale or not a chain.
   */
  ProgramValue getProgram() {
    emscripten::val result = emscripten::val::null();
    with_effect<ShaderChain>([&]<typename SC>(SC &chain) {
      result = emscripten::val::array();
      size_t index = 0;
      for (const auto &entry : chain.chain_ops()) {
        emscripten::val operation = emscripten::val::object();
        operation.set("instance", std::string(entry.instance));
        operation.set("operator", std::string(entry.op->operator_id));
        result.set(index++, operation);
      }
    });
    return result;
  }

  /** @brief Admitted scalar intervals and potential-field topology index bounds.
   * @return Parameter names and inclusive ranges, empty for a stale handle. */
  ParameterRangesValue getShaderChainParameterRanges() {
    emscripten::val result = emscripten::val::array();
    with_effect<ShaderChain>([&]<typename SC>(SC &chain) {
      size_t index = 0;
      chain.parameter_ranges(
          [&](const char *name, float minimum, float maximum) {
            emscripten::val range = emscripten::val::object();
            range.set("name", std::string(name));
            range.set("min", minimum);
            range.set("max", maximum);
            result.set(index++, range);
          });
    });
    return result;
  }

private:
  template <typename Model>
  static bool read_endpoint(const emscripten::val &values,
                            const std::string &instance,
                            typename Model::Params &params) {
    using namespace Pullback::Interp;
    for (uint16_t field = 0; field < SCHEMA<Model>.size(); ++field) {
      const auto &info = SCHEMA<Model>[field];
      const emscripten::val value = values[instance + "." + info.id];
      if (value.isUndefined())
        continue;
      if (!value.isNumber())
        return false;
      const double number = value.as<double>();
      if (!std::isfinite(number))
        return false;
      void *address =
          Detail::ErasedAdapter<Model>::param_address(&params, field);
      if (info.enum_count != 0) {
        if (number < 0.0 || number >= info.enum_count ||
            number != std::floor(number))
          return false;
        *static_cast<uint8_t *>(address) = static_cast<uint8_t>(number);
      } else {
        const float scalar = static_cast<float>(number);
        if (!(scalar >= info.min && scalar <= info.max))
          return false;
        *static_cast<float *>(address) = scalar;
      }
    }
    return true;
  }

  struct SourceMorphEndpoint {
    std::vector<std::string> names;
    std::vector<ShaderChainParameterWrite> writes;
  };

  static ParamSetResult decode_source_endpoint(const emscripten::val &value,
                                               SourceMorphEndpoint &endpoint) {
    if (!hs_wasm::is_js_object(value) || hs_wasm::is_js_array(value))
      return ParamSetResult::MALFORMED_PAYLOAD;
    const auto keys =
        emscripten::val::global("Object").call<emscripten::val>("keys", value);
    const size_t count = keys["length"].as<size_t>();
    if (count > Pullback::Interp::MAX_CHAIN_PARAMS)
      return ParamSetResult::TOO_LONG;
    endpoint.names.resize(count);
    endpoint.writes.resize(count);
    for (size_t index = 0; index < count; ++index) {
      auto &name = endpoint.names[index];
      name = keys[index].as<std::string>();
      if (name.find('\0') != std::string::npos)
        return ParamSetResult::UNKNOWN_PARAM;
      const auto scalar = value[name];
      if (!scalar.isNumber())
        return ParamSetResult::MALFORMED_PAYLOAD;
      endpoint.writes[index] = {name.c_str(), scalar.as<float>()};
    }
    return ParamSetResult::APPLIED;
  }

  struct DecodedColorMorph {
    std::string instance;
    SourceMorphEndpoint from, to;
  };

  template <typename Write>
  static ParamSetResult decode_color_morph(const emscripten::val &entry,
                                           DecodedColorMorph &decoded,
                                           Write &write) {
    if (!hs_wasm::is_js_object(entry) || !entry["instance"].isString() ||
        !entry["mix"].isNumber())
      return ParamSetResult::MALFORMED_PAYLOAD;
    const auto snap = entry["snap"];
    if (!snap.isUndefined() && snap.typeOf().as<std::string>() != "boolean")
      return ParamSetResult::MALFORMED_PAYLOAD;
    decoded.instance = entry["instance"].as<std::string>();
    if (!Pullback::Interp::instance_id_wellformed(decoded.instance))
      return ParamSetResult::UNKNOWN_PARAM;
    const auto from = decode_source_endpoint(entry["from"], decoded.from);
    if (from != ParamSetResult::APPLIED)
      return from;
    const auto to = decode_source_endpoint(entry["to"], decoded.to);
    if (to != ParamSetResult::APPLIED)
      return to;
    const double mix = entry["mix"].as<double>();
    if (!std::isfinite(mix))
      return ParamSetResult::NON_FINITE;
    if (mix < 0.0 || mix > 1.0)
      return ParamSetResult::INADMISSIBLE;
    write = {decoded.instance, decoded.from.writes, decoded.to.writes,
             static_cast<float>(mix), !snap.isUndefined() && snap.as<bool>()};
    return ParamSetResult::APPLIED;
  }

  ParamSetResult
  set_parameters(const emscripten::val &caller_entries, bool preview,
                 const emscripten::val *caller_morph = nullptr,
                 const emscripten::val *caller_modulations = nullptr,
                 const emscripten::val *caller_palette = nullptr) {
    const PayloadDecodeGuard decode_guard(&decoding);
    if (!isValid())
      return ParamSetResult::NO_EFFECT;
    const bool has_morph = caller_morph != nullptr &&
                           !(caller_modulations && caller_morph->isNull());
    const emscripten::val entries = clone_payload(caller_entries);
    const auto morph =
        has_morph ? clone_payload(*caller_morph) : emscripten::val::undefined();
    const auto modulations = caller_modulations == nullptr
                                 ? emscripten::val::undefined()
                                 : clone_payload(*caller_modulations);
    const auto palette = caller_palette == nullptr
                             ? emscripten::val::undefined()
                             : clone_payload(*caller_palette);
    if (!isValid())
      return ParamSetResult::NO_EFFECT;
    if (!hs_wasm::is_js_array(entries))
      return ParamSetResult::MALFORMED_PAYLOAD;
    const size_t count = entries["length"].as<size_t>();
    if (count > Pullback::Interp::MAX_CHAIN_PARAMS)
      return ParamSetResult::TOO_LONG;
    std::vector<std::string> names(count);
    std::vector<ShaderChainParameterWrite> writes(count);
    for (size_t index = 0; index < count; ++index) {
      const emscripten::val entry = entries[index];
      if (entry.isNull() || entry.isUndefined())
        return ParamSetResult::MALFORMED_PAYLOAD;
      const emscripten::val name = entry["name"];
      const emscripten::val value = entry["value"];
      if (!name.isString())
        return ParamSetResult::MALFORMED_PAYLOAD;
      names[index] = name.as<std::string>();
      if (names[index].find('\0') != std::string::npos)
        return ParamSetResult::UNKNOWN_PARAM;
      auto &write = writes[index];
      write.name = names[index].c_str();
      if (value.isNumber()) {
        write.value = value.as<float>();
      } else {
        if (!preview || !hs_wasm::is_js_object(value) ||
            !value["from"].isNumber() || !value["to"].isNumber() ||
            !value["mix"].isNumber())
          return ParamSetResult::MALFORMED_PAYLOAD;
        const double from = value["from"].as<double>();
        const double to = value["to"].as<double>();
        const double mix = value["mix"].as<double>();
        if (!std::isfinite(from) || !std::isfinite(to) || !std::isfinite(mix))
          return ParamSetResult::NON_FINITE;
        if (from < 0.0 || to < 0.0 || from > UINT8_MAX || to > UINT8_MAX ||
            from != std::floor(from) || to != std::floor(to) || mix < 0.0 ||
            mix > 1.0)
          return ParamSetResult::INADMISSIBLE;
        write.blend = Pullback::Interp::EnumBlend{static_cast<uint8_t>(from),
                                                  static_cast<uint8_t>(to),
                                                  static_cast<float>(mix)};
      }
    }
    std::string instance;
    SourceMorphEndpoint from, to;
    ShaderChainSourceMorphWrite source{};
    if (has_morph) {
      if (!hs_wasm::is_js_object(morph) || !morph["instance"].isString() ||
          !morph["evaluation"].isNumber() || !morph["duration"].isNumber() ||
          !morph["easing"].isString())
        return ParamSetResult::MALFORMED_PAYLOAD;
      const auto snap = morph["snap"];
      if (!snap.isUndefined() && snap.typeOf().as<std::string>() != "boolean")
        return ParamSetResult::MALFORMED_PAYLOAD;
      const double evaluation = morph["evaluation"].as<double>();
      const double duration = morph["duration"].as<double>();
      if (!std::isfinite(evaluation) || !std::isfinite(duration))
        return ParamSetResult::NON_FINITE;
      const std::string easing = morph["easing"].as<std::string>();
      if (duration < 1.0 || duration > UINT32_MAX || evaluation < 0.0 ||
          evaluation > duration || evaluation != std::floor(evaluation) ||
          duration != std::floor(duration) ||
          (easing != "LINEAR" && easing != "EASE_IN_OUT_SIN"))
        return ParamSetResult::INADMISSIBLE;
      instance = morph["instance"].as<std::string>();
      if (!Pullback::Interp::instance_id_wellformed(instance))
        return ParamSetResult::UNKNOWN_PARAM;
      const auto from_result = decode_source_endpoint(morph["from"], from);
      if (from_result != ParamSetResult::APPLIED)
        return from_result;
      const auto to_result = decode_source_endpoint(morph["to"], to);
      if (to_result != ParamSetResult::APPLIED)
        return to_result;
      source = {instance,
                from.writes,
                to.writes,
                {static_cast<uint32_t>(evaluation),
                 static_cast<uint32_t>(duration), easing == "EASE_IN_OUT_SIN",
                 !snap.isUndefined() && snap.as<bool>()}};
    }
    std::vector<DecodedColorMorph> layers;
    std::vector<ShaderChainModulationMorphWrite> layer_writes;
    if (caller_modulations != nullptr) {
      if (!hs_wasm::is_js_array(modulations))
        return ParamSetResult::MALFORMED_PAYLOAD;
      const size_t length = modulations["length"].as<size_t>();
      if (length > Pullback::Interp::MAX_CHAIN_OPS)
        return ParamSetResult::TOO_LONG;
      layers.resize(length);
      layer_writes.resize(length);
      for (size_t i = 0; i < length; ++i) {
        const auto result =
            decode_color_morph(modulations[i], layers[i], layer_writes[i]);
        if (result != ParamSetResult::APPLIED)
          return result;
      }
    }
    DecodedColorMorph palette_endpoints;
    ShaderChainPaletteMorphWrite palette_write{};
    if (caller_palette != nullptr) {
      const auto result =
          decode_color_morph(palette, palette_endpoints, palette_write);
      if (result != ParamSetResult::APPLIED)
        return result;
    }
    ParamSetResult result = ParamSetResult::NO_EFFECT;
    with_effect<ShaderChain>([&]<typename SC>(SC &chain) {
      result = chain.update_parameters(
          writes, preview, has_morph ? &source : nullptr, layer_writes,
          caller_palette ? &palette_write : nullptr);
    });
    if (result == ParamSetResult::APPLIED)
      state->paused = state->effect->animations_paused();
    return result;
  }

  bool decoding = false;
  /** @brief Builds the JavaScript result for a chain operation.
   * @param code Engine operation status.
   * @param entry_index Entry index to report.
   * @return Object with code, status and entryIndex fields.
   */
  static emscripten::val chain_result(Pullback::Interp::ChainStatus code,
                                      int entry_index) {
    emscripten::val result = emscripten::val::object();
    result.set("code", emscripten::val(std::string(
                           code == Pullback::Interp::ChainStatus::OK
                               ? "APPLIED"
                               : Pullback::Interp::chain_status_name(code))));
    result.set("status", emscripten::val(code));
    result.set("entryIndex", entry_index);
    return result;
  }
};

/**
 * @brief ShaderChain handle for the live effect.
 * @param state Engine-shared incarnation state.
 * @return New handle, or null if the module or engine is dead or the live
 *         effect is not a ShaderChain.
 */
inline std::shared_ptr<ShaderChainBindings> acquire_shader_chain_bindings(
    const std::shared_ptr<WorkbenchBindingState> &state) {
  if (!workbench_module_alive() || !state->alive || !state->effect ||
      !state->entry)
    return nullptr;
  std::shared_ptr<ShaderChainBindings> result;
  hs_wasm::dispatch_resolution(
      state->width, state->height, [&]<int W, int H>() {
        if (state->entry->type_key == effect_type_key<ShaderChain<W, H>>())
          result = std::make_shared<ShaderChainBindings>(state);
      });
  return result;
}

/** @brief Registers the workbench adapter enums and classes with embind. */
static void bind_workbench_adapters() {
  emscripten::enum_<ChainSnapshotRestoreResult>("ChainSnapshotRestoreResult")
      .value("APPLIED", ChainSnapshotRestoreResult::APPLIED)
      .value("NOT_SHADER_CHAIN", ChainSnapshotRestoreResult::NOT_SHADER_CHAIN)
      .value("UNSUPPORTED_VERSION",
             ChainSnapshotRestoreResult::UNSUPPORTED_VERSION)
      .value("INVALID_LENGTH", ChainSnapshotRestoreResult::INVALID_LENGTH)
      .value("INVALID_VALUE", ChainSnapshotRestoreResult::INVALID_VALUE)
      .value("INVALID_CHAIN", ChainSnapshotRestoreResult::INVALID_CHAIN);
  emscripten::class_<ShaderChainBindings>("ShaderChainBindings")
      .smart_ptr<std::shared_ptr<ShaderChainBindings>>(
          "ShaderChainBindingsHandle")
      .function("isValid", &ShaderChainBindings::isValid)
      .function("getSnapshot", &ShaderChainBindings::getSnapshot)
      .function("restoreSnapshot", &ShaderChainBindings::restoreSnapshot)
      .function("setShaderChain", &ShaderChainBindings::setShaderChain)
      .function("setShaderChainParameters",
                &ShaderChainBindings::setShaderChainParameters)
      .function("setShaderChainTransitionParameters",
                &ShaderChainBindings::setShaderChainTransitionParameters)
      .function("setShaderChainSourceMorphParameters",
                &ShaderChainBindings::setShaderChainSourceMorphParameters)
      .function("setShaderChainMorphParameters",
                &ShaderChainBindings::setShaderChainMorphParameters)
      .function("setShaderChainPaletteMorphParameters",
                &ShaderChainBindings::setShaderChainPaletteMorphParameters)
      .class_function("admitsPlanarSourceTransition",
                      &ShaderChainBindings::admitsPlanarSourceTransition)
      .function("getShaderChainParameterRanges",
                &ShaderChainBindings::getShaderChainParameterRanges)
      .function("getProgram", &ShaderChainBindings::getProgram)
      .class_function("getShaderChainCatalog",
                      &ShaderChainBindings::getShaderChainCatalog);
}
