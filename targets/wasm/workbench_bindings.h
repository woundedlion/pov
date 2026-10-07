/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file workbench_bindings.h
 * @brief Typed WASM authoring adapters bound to an effect incarnation.
 */
#pragma once
#include <emscripten/bind.h>
#include "targets/wasm/effect_factory.h"
#include "targets/wasm/payload_clone.h"
#include "targets/wasm/param_marshal.h"
#if HS_ENABLE_CHAIN_INTERPRETER
#include "core/render/pullback/catalog_export.h"
#include "targets/wasm/chain_snapshot_codec.h"
#endif
#include <functional>
#include <memory>
#include <span>
#include <vector>

#if HS_ENABLE_CHAIN_INTERPRETER
// Caller property access can re-enter embind, including delete().
static bool snapshot_decode_active = false;
struct SnapshotDecodeGuard {
  bool *decoding;
  explicit SnapshotDecodeGuard(bool *decoding = nullptr) : decoding(decoding) {
    HS_CHECK(!snapshot_decode_active,
             "re-entrant engine decode from a caller accessor");
    snapshot_decode_active = true;
    if (decoding)
      *decoding = true;
  }
  ~SnapshotDecodeGuard() {
    if (decoding)
      *decoding = false;
    snapshot_decode_active = false;
  }
};
#endif

// clang-format off
EM_JS(bool, workbench_module_alive, (), {
  return !ABORT && !Module['HS_MODULE_DEAD'];
});
// clang-format on

/** @brief WASM-only effect incarnation shared by engine and capability handles. */
struct WorkbenchBindingState {
  Effect *effect = nullptr;
  const FactoryEntry *entry = nullptr;
  uint64_t generation = 0;
  int width = 0;
  int height = 0;
  bool alive = true;
  bool paused = false;
};

class WorkbenchBindings {
public:
  using RebuildRestore =
      std::function<bool(const std::shared_ptr<WorkbenchBindingState> &)>;
  explicit WorkbenchBindings(std::shared_ptr<WorkbenchBindingState> state)
      : state(std::move(state)), generation(this->state->generation) {}
  /** @brief Whether the originating engine and effect incarnation remain live. */
  bool isValid() const {
    return workbench_module_alive() && state->alive && state->effect &&
           state->entry && state->generation == generation;
  }

protected:
  std::shared_ptr<WorkbenchBindingState> state;
#if HS_ENABLE_CHAIN_INTERPRETER
  /** @brief Runs a callback when the live effect has the requested factory type. */
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
#endif

  static bool is_array(const emscripten::val &value) {
    return emscripten::val::global("Array").call<bool>("isArray", value);
  }

private:
  const uint64_t generation;
};
#if HS_ENABLE_CHAIN_INTERPRETER
class ShaderChainBindings : public WorkbenchBindings {
public:
  using WorkbenchBindings::WorkbenchBindings;
  ~ShaderChainBindings() {
    HS_CHECK(!decoding,
             "delete() of a chain handle from a caller accessor during decode");
  }
  bool isValid() const { return WorkbenchBindings::isValid(); }
  emscripten::val getSnapshot() {
    auto result = emscripten::val::null();
    with_effect<ShaderChain>([&](auto &chain) {
      result = hs_wasm::ChainSnapshotCodec::encode(chain.snapshot());
    });
    return result;
  }

  ChainSnapshotRestoreResult
  restoreSnapshot(const emscripten::val &caller_input) {
    using Result = ChainSnapshotRestoreResult;
    const SnapshotDecodeGuard guard(&decoding);
    if (!isValid())
      return Result::NOT_SHADER_CHAIN;
    const uint64_t owner_generation = state->generation;
    const auto schema_generation =
        state->effect->getParameterSchemaGeneration();
    const Effect *const owner = state->effect;
    ChainSnapshot snapshot;
    const auto result = hs_wasm::ChainSnapshotCodec::decode(
        clone_payload(caller_input), snapshot);
    if (!isValid() || state->generation != owner_generation ||
        state->effect != owner ||
        state->effect->getParameterSchemaGeneration() != schema_generation)
      return Result::NOT_SHADER_CHAIN;
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
  emscripten::val setShaderChain(const emscripten::val &caller_entries) {
    const SnapshotDecodeGuard decode_guard(&decoding);
    using Pullback::Interp::ChainStatus;
    if (!with_effect<ShaderChain>([]<typename SC>(SC &) {}))
      return chain_result(ChainStatus::NOT_CHAIN_EFFECT, -1);
    // Payload cloning can invoke getters that replace the addressed chain.
    const uint64_t owner_generation = state->generation;
    const Effect *const owner = state->effect;
    const void *const owner_type_key = state->entry->type_key;
    const emscripten::val entries = clone_payload(caller_entries);
    if (!is_array(entries))
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
    if (state->generation != owner_generation || state->effect != owner ||
        state->entry->type_key != owner_type_key)
      return chain_result(ChainStatus::NOT_CHAIN_EFFECT, -1);
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
   */
  ParamSetResult
  setShaderChainParameters(const emscripten::val &caller_entries) {
    const SnapshotDecodeGuard decode_guard(&decoding);
    if (!with_effect<ShaderChain>([]<typename SC>(SC &) {}))
      return ParamSetResult::NO_EFFECT;
    const uint64_t owner_generation = state->generation;
    const Effect *const owner = state->effect;
    const auto schema_generation =
        state->effect->getParameterSchemaGeneration();
    const emscripten::val entries = clone_payload(caller_entries);
    if (!is_array(entries))
      return ParamSetResult::MALFORMED_PAYLOAD;
    const size_t count = entries["length"].as<size_t>();
    if (count > Pullback::Interp::MAX_CHAIN_PARAMS)
      return ParamSetResult::TOO_LONG;
    std::vector<std::string> names(count);
    std::vector<float> values(count);
    for (size_t index = 0; index < count; ++index) {
      const emscripten::val entry = entries[index];
      if (entry.isNull() || entry.isUndefined())
        return ParamSetResult::MALFORMED_PAYLOAD;
      const emscripten::val name = entry["name"];
      const emscripten::val value = entry["value"];
      if (!name.isString())
        return ParamSetResult::MALFORMED_PAYLOAD;
      if (!value.isNumber())
        return ParamSetResult::MALFORMED_PAYLOAD;
      names[index] = name.as<std::string>();
      values[index] = value.as<float>();
    }
    if (state->generation != owner_generation || state->effect != owner ||
        state->effect->getParameterSchemaGeneration() != schema_generation)
      return ParamSetResult::NO_EFFECT;
    std::vector<ShaderChainParameterWrite> writes(count);
    for (size_t index = 0; index < count; ++index)
      writes[index] = {names[index].c_str(), values[index]};
    ParamSetResult result = ParamSetResult::NO_EFFECT;
    with_effect<ShaderChain>([&]<typename SC>(SC &chain) {
      result = chain.update_parameters(writes);
    });
    if (result == ParamSetResult::APPLIED)
      state->paused = state->effect->animations_paused();
    return result;
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

  emscripten::val getProgram() {
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

private:
  bool decoding = false;
  /** @brief Result with ChainStatus, commit or refusal code, and entry index. */
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
#endif

template <template <int, int> class EffectT, typename Bindings>
std::shared_ptr<Bindings> acquire_workbench_bindings(
    const std::shared_ptr<WorkbenchBindingState> &state) {
  if (!workbench_module_alive() || !state->alive || !state->effect ||
      !state->entry)
    return nullptr;
  std::shared_ptr<Bindings> result;
  hs_wasm::dispatch_resolution(
      state->width, state->height, [&]<int W, int H>() {
        if (state->entry->type_key == effect_type_key<EffectT<W, H>>())
          result = std::make_shared<Bindings>(state);
      });
  return result;
}

#if HS_ENABLE_CHAIN_INTERPRETER
inline std::shared_ptr<ShaderChainBindings> acquire_shader_chain_bindings(
    const std::shared_ptr<WorkbenchBindingState> &state) {
  return acquire_workbench_bindings<ShaderChain, ShaderChainBindings>(state);
}
#endif

static void bind_workbench_adapters() {
#if HS_ENABLE_CHAIN_INTERPRETER
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
      .function("getProgram", &ShaderChainBindings::getProgram)
      .class_function("getShaderChainCatalog",
                      &ShaderChainBindings::getShaderChainCatalog);
#endif
}
