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
 * @details Mirrors what getParameterDefinitions() emits per entry. `is_bool`
 *          toggles carry no meaningful range (the GUI keys off the type and
 *          ignores min/max for them).
 */
struct ParamView {
  const char *name; /**< Parameter name, as exposed to the JS boundary. */
  float value;      /**< Current value in the parameter's native units. */
  float min;        /**< Inclusive lower bound; ignored when is_bool. */
  float max;        /**< Inclusive upper bound; ignored when is_bool. */
  bool is_bool;     /**< True if the parameter is a boolean toggle. */
  bool animated;    /**< True if the parameter is currently animated. */
  bool readonly;    /**< True if the parameter is read-only (not editable). */
};

/**
 * @brief Snapshot an effect's parameters into `out`, in definition order.
 * @param effect Effect whose getParameters() sequence defines the order.
 * @param out Destination vector, cleared then filled in definition order;
 *            caller-owned so a reused vector amortizes its allocation.
 * @details The single source of ordering for the definition stream.
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
 * @param effect Effect whose getParameters() values are read, in order.
 * @param out Destination vector, cleared (retaining capacity) then filled so
 *            that out[i] corresponds to collect_param_views()'s view[i].
 * @details Iterates the identical Effect::getParameters() sequence, so
 *          values[i] always corresponds to view[i]. `out.clear()` retains
 *          capacity, so a caller that has reserved MAX-params up front gets a
 *          zero-reallocation guarantee for the per-frame memory-view contract.
 */
inline void fill_param_values(const Effect &effect, std::vector<float> &out) {
  out.clear();
  for (const auto &def : effect.getParameters()) {
    out.push_back(def.get());
  }
}

} // namespace hs_wasm
