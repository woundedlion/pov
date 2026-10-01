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
#endif
#include <cmath>
#include <functional>
#include <memory>
#include <span>
#include <vector>

#if HS_ENABLE_SHADER_WORKBENCH
/**
 * @brief Outcome of restoreFullConfigSnapshot().
 * @details Exposed to JS as the Module.FullConfigRestoreResult embind enum;
 *          compare against its values, never by truthiness. Every value but
 *          APPLIED leaves the effect exactly as it was.
 */
enum class FullConfigRestoreResult : uint8_t {
  APPLIED,              /**< Snapshot installed. */
  NOT_SHADER_WORKBENCH, /**< The loaded effect has no full configuration. */
  UNSUPPORTED_VERSION,  /**< schemaVersion is not the current schema. */
  INVALID_LENGTH,       /**< Snapshot missing, or an array whose length is not
                            the field count. */
  INVALID_VALUE,        /**< A field or runtime value outside what its slot
                            admits. */
  INVALID_ACCEPTED,     /**< Fields each in range but a combination the effect
                            will not render. */
  INVALID_PENDING,      /**< Pending list is absent, is not a set of in-range
                            field indices, or does not match where accepted and
                            requested differ. Retry with an empty list. */
};
#endif // HS_ENABLE_SHADER_WORKBENCH

#if HS_ENABLE_SHADER_WORKBENCH || HS_ENABLE_CHAIN_INTERPRETER
// Caller property access can re-enter embind, including delete().
static bool snapshot_decode_active = false;
struct SnapshotDecodeGuard {
  SnapshotDecodeGuard() {
    HS_CHECK(!snapshot_decode_active,
             "re-entrant engine decode from a caller accessor");
    snapshot_decode_active = true;
  }
  ~SnapshotDecodeGuard() { snapshot_decode_active = false; }
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
#if HS_ENABLE_SHADER_WORKBENCH || HS_ENABLE_CHAIN_INTERPRETER
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

  /** @brief Why a snapshot array failed to decode. */
  enum class ArrayDecode : uint8_t {
    OK,
    BAD_LENGTH, /**< Not an array, or not the field count long. */
    BAD_VALUE,  /**< An element outside what its slot admits. */
  };

  static bool is_array(const emscripten::val &value) {
    return emscripten::val::global("Array").call<bool>("isArray", value);
  }

  static bool whole_uint32(double value) {
    return std::isfinite(value) && value >= 0.0 &&
           value <= static_cast<double>(UINT32_MAX) &&
           value == std::floor(value);
  }

  template <size_t N>
  static emscripten::val uint32_array(const std::array<uint32_t, N> &values) {
    emscripten::val output = emscripten::val::array();
    for (size_t index = 0; index < N; ++index)
      output.set(index, values[index]);
    return output;
  }

  template <size_t N>
  static emscripten::val float_array(const std::array<float, N> &values) {
    emscripten::val output = emscripten::val::array();
    for (size_t index = 0; index < N; ++index)
      output.set(index, values[index]);
    return output;
  }

  template <size_t N>
  static ArrayDecode decode_uint32_array(const emscripten::val &input,
                                         std::array<uint32_t, N> &output) {
    if (!is_array(input) || input["length"].as<size_t>() != N)
      return ArrayDecode::BAD_LENGTH;
    for (size_t index = 0; index < N; ++index) {
      const emscripten::val element = input[index];
      if (!element.isNumber())
        return ArrayDecode::BAD_VALUE;
      const double number = element.as<double>();
      if (!whole_uint32(number))
        return ArrayDecode::BAD_VALUE;
      output[index] = static_cast<uint32_t>(number);
    }
    return ArrayDecode::OK;
  }

  template <size_t N>
  static ArrayDecode decode_runtime(const emscripten::val &input,
                                    std::array<float, N> &output) {
    if (!is_array(input) || input["length"].as<size_t>() != N)
      return ArrayDecode::BAD_LENGTH;
    for (size_t index = 0; index < N; ++index) {
      const emscripten::val element = input[index];
      if (!element.isNumber())
        return ArrayDecode::BAD_VALUE;
      output[index] = element.as<float>();
      if (!std::isfinite(output[index]))
        return ArrayDecode::BAD_VALUE;
    }
    return ArrayDecode::OK;
  }

  void check_param_capacity() const {
    HS_CHECK(state->effect->getParameters().size() <=
                 hs_wasm::ParamStreams::CAPACITY,
             "workbench parameter capacity exceeded");
  }

private:
  const uint64_t generation;
};
#if HS_ENABLE_SHADER_WORKBENCH
class LegacyShaderBindings : public WorkbenchBindings {
public:
  using WorkbenchBindings::WorkbenchBindings;
  bool isValid() const { return WorkbenchBindings::isValid(); }
  RebuildRestore capture_rebuild_state() {
    RebuildRestore restore;
    with_effect<Shader>([&]<typename SB>(SB &shader) {
      const auto snapshot = std::make_shared<typename SB::FullConfigSnapshot>(
          shader.capture_full_config_snapshot());
      restore = [snapshot](const auto &state) {
        LegacyShaderBindings bindings(state);
        bool restored = false;
        bindings.with_effect<Shader>([&]<typename Target>(Target &target) {
          if constexpr (std::is_same_v<typename SB::FullConfigSnapshot,
                                       typename Target::FullConfigSnapshot>)
            restored = target.restore_full_config_snapshot(*snapshot) ==
                       Target::ConfigRestoreResult::APPLIED;
        });
        return restored;
      };
    });
    return restore;
  }
#if HS_ENABLE_SHADER_WORKBENCH
  /**
   * @brief Returns the current Shader workbench's versioned full-state snapshot.
   * @return JS object {schemaVersion, accepted, requested, pendingFieldIds,
   *         hasRuntime, runtime}, or null when the loaded effect is not
   *         Shader workbench.
   * @details `accepted` and `requested` are CONFIG_FIELD_COUNT-long arrays of
   *          uint32-encoded field values in ConfigFieldId order;
   *          getFullConfigFieldDefinitions() names the indices.
   *          `pendingFieldIds` lists the fields carrying an unresolved edit,
   *          which at the current schema are exactly the fields where the two
   *          arrays differ. `runtime` is the animation clock state and is
   *          meaningful only when `hasRuntime`. The whole configuration crosses
   *          as one object because Shader vets slots and params together:
   *          replaying the parameter stream entry by entry walks through
   *          combinations it refuses.
   */
  emscripten::val getFullConfigSnapshot() {
    emscripten::val output = emscripten::val::null();
    with_effect<Shader>([&]<typename SB>(SB &shader) {
      const typename SB::FullConfigSnapshot snapshot =
          shader.capture_full_config_snapshot();
      output = emscripten::val::object();
      output.set("schemaVersion", snapshot.schema_version);
      output.set("accepted", uint32_array(snapshot.accepted));
      output.set("requested", uint32_array(snapshot.requested));
      emscripten::val pending_ids = emscripten::val::array();
      size_t pending_count = 0;
      for (size_t index = 0; index < snapshot.pending.size(); ++index)
        if (snapshot.pending[index] != 0)
          pending_ids.set(pending_count++, index);
      output.set("pendingFieldIds", pending_ids);
      output.set("hasRuntime", snapshot.has_runtime);
      output.set("runtime", float_array(snapshot.runtime));
    });
    return output;
  }

  /**
   * @brief Atomically restores a current Shader workbench snapshot.
   * @param caller_input Object in getFullConfigSnapshot()'s shape.
   * @return APPLIED, or the reason the snapshot was refused.
   * @details Rejections leave the effect exactly as it was, so a failed restore
   *          needs no rollback. NOT_SHADER_WORKBENCH covers the loaded effect;
   *          UNSUPPORTED_VERSION a schemaVersion without a supported migration;
   *          INVALID_LENGTH a missing snapshot or an array whose length is not
   *          the field count; INVALID_VALUE a field or runtime value outside
   *          what its slot admits; INVALID_ACCEPTED fields each in range but a
   *          combination the effect will not render; INVALID_PENDING a pending
   *          list that is absent, that is not a set of in-range field indices,
   *          or that does not match where accepted and requested differ.
   *          NOT_SHADER_WORKBENCH also covers an input whose accessors swap
   *          the loaded effect out while the snapshot is being decoded.
   */
  FullConfigRestoreResult
  restoreFullConfigSnapshot(const emscripten::val &caller_input) {
    FullConfigRestoreResult result =
        FullConfigRestoreResult::NOT_SHADER_WORKBENCH;
    with_effect<Shader>([&]<typename SB>(SB &shader) {
      // Payload cloning can invoke getters that replace or delete the owner.
      const SnapshotDecodeGuard decode_guard;
      const uint64_t owner_generation = state->generation;
      const Effect *const owner = state->effect;
      const void *const owner_type_key = state->entry->type_key;
      const emscripten::val input = clone_payload(caller_input);
      if (input.isUndefined() || input.isNull()) {
        result = FullConfigRestoreResult::INVALID_LENGTH;
        return;
      }
      const emscripten::val schema_version = input["schemaVersion"];
      if (!schema_version.isNumber()) {
        result = FullConfigRestoreResult::UNSUPPORTED_VERSION;
        return;
      }
      const double schema_number = schema_version.as<double>();
      if (!whole_uint32(schema_number)) {
        result = FullConfigRestoreResult::UNSUPPORTED_VERSION;
        return;
      }
      const auto version = static_cast<uint32_t>(schema_number);
      if (!SB::config_version_supported(version)) {
        result = FullConfigRestoreResult::UNSUPPORTED_VERSION;
        return;
      }
      auto decode_and_restore = [&](auto &snapshot) {
        snapshot.schema_version = version;
        const emscripten::val has_runtime = input["hasRuntime"];
        if (!has_runtime.isTrue() && !has_runtime.isFalse()) {
          result = FullConfigRestoreResult::INVALID_VALUE;
          return;
        }
        snapshot.has_runtime = has_runtime.as<bool>();
        ArrayDecode decoded =
            decode_uint32_array(input["accepted"], snapshot.accepted);
        if (decoded == ArrayDecode::OK)
          decoded = decode_uint32_array(input["requested"], snapshot.requested);
        if (decoded == ArrayDecode::OK && snapshot.has_runtime)
          decoded = decode_runtime(input["runtime"], snapshot.runtime);
        if (decoded != ArrayDecode::OK) {
          result = decoded == ArrayDecode::BAD_LENGTH
                       ? FullConfigRestoreResult::INVALID_LENGTH
                       : FullConfigRestoreResult::INVALID_VALUE;
          return;
        }
        const emscripten::val pending_ids = input["pendingFieldIds"];
        if (!is_array(pending_ids)) {
          result = FullConfigRestoreResult::INVALID_PENDING;
          return;
        }
        const size_t pending_count = pending_ids["length"].as<size_t>();
        if (pending_count > snapshot.pending.size()) {
          result = FullConfigRestoreResult::INVALID_PENDING;
          return;
        }
        for (size_t index = 0; index < pending_count; ++index) {
          const emscripten::val field_id = pending_ids[index];
          if (!field_id.isNumber()) {
            result = FullConfigRestoreResult::INVALID_PENDING;
            return;
          }
          const double field_number = field_id.as<double>();
          if (!whole_uint32(field_number) ||
              field_number >= snapshot.pending.size()) {
            result = FullConfigRestoreResult::INVALID_PENDING;
            return;
          }
          uint8_t &pending =
              snapshot.pending[static_cast<size_t>(field_number)];
          if (pending != 0) {
            result = FullConfigRestoreResult::INVALID_PENDING;
            return;
          }
          pending = 1;
        }
        if (state->generation != owner_generation || state->effect != owner ||
            state->entry->type_key != owner_type_key) {
          result = FullConfigRestoreResult::NOT_SHADER_WORKBENCH;
          return;
        }
        result = map_restore_result<SB>(
            shader.restore_full_config_snapshot(snapshot));
        if (result == FullConfigRestoreResult::APPLIED)
          check_param_capacity();
      };
      if (version == SB::LEGACY_CONFIG_SCHEMA_VERSION) {
        typename SB::LegacyFullConfigSnapshot snapshot;
        decode_and_restore(snapshot);
      } else {
        typename SB::FullConfigSnapshot snapshot;
        decode_and_restore(snapshot);
      }
    });
    return result;
  }

  /**
   * @brief Returns stable Shader workbench field IDs and diagnostic names.
   * @return JS array of {id, name} in ConfigFieldId order, or null when the
   *         loaded effect is not Shader.
   * @details `id` is the index into a snapshot's accepted/requested arrays and
   *          the value pendingFieldIds carries; `name` is the field's dotted
   *          config path. A caller labels a field through this rather than a
   *          hardcoded index, which moves when the schema gains a field.
   */
  emscripten::val getFullConfigFieldDefinitions() {
    emscripten::val output = emscripten::val::null();
    with_effect<Shader>([&]<typename SB>(SB &) {
      output = emscripten::val::array();
      for (size_t index = 0; index < SB::CONFIG_FIELD_COUNT; ++index) {
        emscripten::val field = emscripten::val::object();
        field.set("id", index);
        field.set("name", SB::config_field_name(
                              static_cast<typename SB::ConfigFieldId>(index)));
        output.set(index, field);
      }
    });
    return output;
  }

#endif // HS_ENABLE_SHADER_WORKBENCH
private:
#if HS_ENABLE_SHADER_WORKBENCH
  template <typename SB>
  static FullConfigRestoreResult
  map_restore_result(typename SB::ConfigRestoreResult result) {
    switch (result) {
    case SB::ConfigRestoreResult::APPLIED:
      return FullConfigRestoreResult::APPLIED;
    case SB::ConfigRestoreResult::UNSUPPORTED_VERSION:
      return FullConfigRestoreResult::UNSUPPORTED_VERSION;
    case SB::ConfigRestoreResult::INVALID_VALUE:
      return FullConfigRestoreResult::INVALID_VALUE;
    case SB::ConfigRestoreResult::INVALID_ACCEPTED:
      return FullConfigRestoreResult::INVALID_ACCEPTED;
    case SB::ConfigRestoreResult::INVALID_PENDING:
      return FullConfigRestoreResult::INVALID_PENDING;
    }
    __builtin_unreachable();
  }
#endif // HS_ENABLE_SHADER_WORKBENCH
};
#endif
#if HS_ENABLE_CHAIN_INTERPRETER
class ShaderChainBindings : public WorkbenchBindings {
public:
  using WorkbenchBindings::WorkbenchBindings;
  bool isValid() const { return WorkbenchBindings::isValid(); }
  RebuildRestore capture_rebuild_state() {
    RebuildRestore restore;
    with_effect<ShaderChain>([&]<typename SC>(SC &chain) {
      std::vector<std::pair<std::string, std::string>> entries;
      for (const auto &entry : chain.chain_ops())
        entries.emplace_back(entry.instance, entry.op->operator_id);
      std::vector<std::pair<std::string, float>> parameters;
      chain.refresh_parameter_display();
      for (const auto &parameter : chain.getParameters())
        if (!parameter.readonly)
          parameters.emplace_back(parameter.name, parameter.get_requested());
      restore = [entries = std::move(entries),
                 parameters = std::move(parameters)](const auto &state) {
        ShaderChainBindings bindings(state);
        bool restored = false;
        bindings.with_effect<ShaderChain>([&]<typename Target>(Target &target) {
          std::vector<Pullback::Interp::ChainEntryRequest> requests;
          for (const auto &[instance, operation] : entries)
            requests.push_back({instance, operation});
          if (target.set_chain(requests).code !=
              Pullback::Interp::ChainStatus::OK)
            return;
          std::vector<ShaderChainParameterWrite> writes;
          for (const auto &[name, value] : parameters)
            writes.push_back({name.c_str(), value});
          restored =
              target.update_parameters(writes) == ParamSetResult::APPLIED;
        });
        return restored;
      };
    });
    return restore;
  }

#if HS_ENABLE_CHAIN_INTERPRETER
  /**
   * @brief Compiles a chain program shape on the loaded ShaderChain effect.
   * @param caller_entries JS array of {instance, operator} string pairs — the ordered
   *        program shape and nothing else. No values, no offsets, no family
   *        tags.
   * @return JS object {code, status, entryIndex}: code is "APPLIED" on commit,
   *         otherwise the refusal name; status is the ChainStatus enum.
   *         entryIndex names the offending entry, -1 for a whole-chain refusal.
   * @details Synchronous: on APPLIED the parameter definitions are already
   * rebuilt and the schema generation bumped before this returns, so the
   * caller applies preset values by "{instance}.{field-id}" name immediately
   * (apply order: setShaderChain -> values -> syncEffectGui -> invalidate).
   * The boundary rejects a non-array payload or a non-string entry field as
   * MALFORMED_PAYLOAD and never traps; NOT_CHAIN_EFFECT reports that the
   * loaded effect is not ShaderChain, and covers an input whose accessors swap
   * the loaded effect out while the entries are being decoded. A refusal
   * leaves the previous program, its parameter definitions, the generation,
   * and all instance state untouched.
   */
  emscripten::val setShaderChain(const emscripten::val &caller_entries) {
    const SnapshotDecodeGuard decode_guard;
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
    // The length cap precedes per-entry decode — compile()'s own shape order —
    // so an oversized payload is TOO_LONG, never an unbounded decode.
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
      check_param_capacity();
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
    const SnapshotDecodeGuard decode_guard;
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
   *         entry. Budgets, carriers, operator ids and parameter schemas
   *         match the native suite's golden pin; the block sizes do not.
   *         These are wasm32 ABI figures, so a pointer-bearing operator's
   *         `prepared` block is narrower here than in the LP64 golden.
   *         An editor budgets arena bytes against these figures, which
   *         are the ones this module's own runtime allocates from.
   */
  static std::string getShaderChainCatalog() {
    std::string catalog;
    Pullback::Interp::append_catalog_json(catalog);
    return catalog;
  }
#endif // HS_ENABLE_CHAIN_INTERPRETER

#if HS_ENABLE_CHAIN_INTERPRETER
  /** @brief Result with enum status, legacy string code, and entry index. */
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
#endif // HS_ENABLE_CHAIN_INTERPRETER
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

#if HS_ENABLE_SHADER_WORKBENCH
inline std::shared_ptr<LegacyShaderBindings> acquire_legacy_shader_bindings(
    const std::shared_ptr<WorkbenchBindingState> &state) {
  return acquire_workbench_bindings<Shader, LegacyShaderBindings>(state);
}
#endif
#if HS_ENABLE_CHAIN_INTERPRETER
inline std::shared_ptr<ShaderChainBindings> acquire_shader_chain_bindings(
    const std::shared_ptr<WorkbenchBindingState> &state) {
  return acquire_workbench_bindings<ShaderChain, ShaderChainBindings>(state);
}
#endif

static void bind_workbench_adapters() {
#if HS_ENABLE_SHADER_WORKBENCH
  emscripten::class_<LegacyShaderBindings>("LegacyShaderBindings")
      .smart_ptr<std::shared_ptr<LegacyShaderBindings>>(
          "LegacyShaderBindingsHandle")
      .function("isValid", &LegacyShaderBindings::isValid)
      .function("getFullConfigSnapshot",
                &LegacyShaderBindings::getFullConfigSnapshot)
      .function("restoreFullConfigSnapshot",
                &LegacyShaderBindings::restoreFullConfigSnapshot)
      .function("getFullConfigFieldDefinitions",
                &LegacyShaderBindings::getFullConfigFieldDefinitions);
#endif
#if HS_ENABLE_CHAIN_INTERPRETER
  emscripten::class_<ShaderChainBindings>("ShaderChainBindings")
      .smart_ptr<std::shared_ptr<ShaderChainBindings>>(
          "ShaderChainBindingsHandle")
      .function("isValid", &ShaderChainBindings::isValid)
      .function("setShaderChain", &ShaderChainBindings::setShaderChain)
      .function("setShaderChainParameters",
                &ShaderChainBindings::setShaderChainParameters)
      .function("getProgram", &ShaderChainBindings::getProgram)
      .class_function("getShaderChainCatalog",
                      &ShaderChainBindings::getShaderChainCatalog);
#endif
}
