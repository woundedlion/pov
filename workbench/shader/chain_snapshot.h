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

struct ChainSnapshot {
  struct Entry {
    std::string instance;
    std::string operator_id;
  };
  struct Parameter {
    std::string name;
    float value;
  };
  struct Runtime {
    std::string instance;
    Pullback::Interp::RuntimeSnapshot state;
  };
  static constexpr uint32_t SCHEMA_VERSION = 2;
  uint32_t schema_version = SCHEMA_VERSION;
  std::vector<Entry> chain;
  std::vector<Parameter> parameters;
  std::optional<std::vector<Runtime>> runtime;
  std::optional<GeneratedPaletteBank::Snapshot> palette_bank;
  bool animations_paused = false;
};
