/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file param_marshal.h
 * @brief Pure (no-Emscripten) parameter-marshaling layer for the WASM bridge.
 *
 * The JS frontend consumes two parallel streams from the engine: the parameter
 * *definitions* (getParameterDefinitions) and the per-frame *values*
 * (getParamValues). They MUST agree on order — values[i] describes
 * definitions[i] — or every slider mis-binds. Both are derived here from a
 * single pass over Effect::getParameters(), so the order is defined in exactly
 * one place. This layer carries no Emscripten dependency (wasm.cpp adds the
 * emscripten::val translation on top), so the contract is host-unit-testable
 * without an Emscripten toolchain — see tests/test_param_marshal.h.
 */
#pragma once

#include <vector>

#include "core/canvas.h" // Effect, Effect::ParamDef

namespace hs_wasm {

/**
 * @brief One parameter as the JS boundary sees it, in definition order.
 *
 * Mirrors what getParameterDefinitions() emits per entry. `is_bool` toggles
 * carry no meaningful range (the GUI keys off the type and ignores min/max for
 * them), matching wasm.cpp's deliberate omission of min/max for bool params.
 */
struct ParamView {
  const char *name;
  float value;
  float min;
  float max;
  bool is_bool;
  bool animated;
  bool readonly;
};

/**
 * @brief Snapshot an effect's parameters into `out`, in definition order.
 *
 * The single source of ordering for the definition stream. `out` is caller-
 * owned so a reused vector amortizes its allocation.
 */
inline void collect_param_views(const Effect &effect,
                                std::vector<ParamView> &out) {
  out.clear();
  for (const auto &def : effect.getParameters()) {
    out.push_back(ParamView{def.name, def.get(), def.min, def.max,
                            def.is_bool(), def.animated, def.readonly});
  }
}

/**
 * @brief Fill `out` with current parameter values, in the SAME order as
 *        collect_param_views() — this is the getParamValues() stream.
 *
 * Iterates the identical Effect::getParameters() sequence, so values[i] always
 * corresponds to view[i]. `out.clear()` retains capacity, so a caller that has
 * reserved MAX-params up front (wasm.cpp's `paramValues`) gets the same zero-
 * reallocation guarantee the per-frame memory-view contract depends on.
 */
inline void fill_param_values(const Effect &effect, std::vector<float> &out) {
  out.clear();
  for (const auto &def : effect.getParameters()) {
    out.push_back(def.get());
  }
}

} // namespace hs_wasm
