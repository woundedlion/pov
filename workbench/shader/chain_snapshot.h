/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <optional>
#include <string>
#include <vector>
#include "core/color/palette_cycler.h"
#include "core/render/pullback/runtime_snapshot.h"

enum class ChainSnapshotRestoreResult {
  APPLIED,
  NOT_SHADER_CHAIN,
  UNSUPPORTED_VERSION,
  INVALID_LENGTH,
  INVALID_VALUE,
  INVALID_CHAIN
};

/** @brief Serializable ShaderChain state: program, values and runtime. */
struct ChainSnapshot {
  /** @brief One program entry. */
  struct Entry {
    std::string instance;    ///< Instance name, unique within the chain.
    std::string operator_id; ///< Operator-table ID.
  };
  /** @brief One parameter value. */
  struct Parameter {
    std::string name; ///< "{instance}.{field-id}" parameter name.
    float value;      ///< Value; enum fields store the option index.
  };
  /** @brief Captured per-instance runtime state. */
  struct Runtime {
    std::string instance;                    ///< Owning instance name.
    Pullback::Interp::RuntimeSnapshot state; ///< Operator runtime state.
  };
  static constexpr uint32_t SCHEMA_VERSION = 2; ///< Current snapshot schema.
  uint32_t schema_version = SCHEMA_VERSION;     ///< Schema this snapshot uses.
  std::vector<Entry> chain;          ///< Program entries, in chain order.
  std::vector<Parameter> parameters; ///< Values to restore, unique names.
  /// Runtime state for stateful instances; absent starts them fresh.
  std::optional<std::vector<Runtime>> runtime;
  /// Generated-palette bank state; absent resets the bank to its defaults.
  std::optional<GeneratedPaletteBank::Snapshot> palette_bank;
  bool animations_paused = false; ///< Effect animations-paused flag.
};
