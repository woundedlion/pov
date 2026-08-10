/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file param_marshal.h
 * @brief Pure (no-Emscripten) parameter-marshaling layer for the WASM bridge.
 *
 * The JS frontend consumes two parallel streams from the engine: the parameter
 * *definitions* (getParameterDefinitions) and the per-frame *values*
 * (getParamValues). They MUST agree on order — values[i] describes
 * definitions[i] — or every slider mis-binds. This layer takes two independent
 * passes over Effect::getParameters() at different times. Effect implementations
 * expose a schema-generation token; the WASM bridge joins it with effect
 * replacement so dynamic rebinds invalidate stale positional bindings. This
 * layer carries no Emscripten dependency (wasm.cpp adds the
 * emscripten::val translation on top), so the contract is host-unit-testable
 * without an Emscripten toolchain — see tests/test_param_marshal.h.
 */
#pragma once

#include <vector>

#include "core/render/canvas.h" // Effect, Effect::ParamDef

namespace hs_wasm {

/** @brief Joins engine replacement and effect-local schema changes into one token. */
class ParamGenerationTracker {
public:
  void replace(uint32_t schema_generation) {
    observed_schema_generation = schema_generation;
    ++generation_value;
  }

  void observe(uint32_t schema_generation) {
    if (schema_generation == observed_schema_generation)
      return;
    observed_schema_generation = schema_generation;
    ++generation_value;
  }

  uint32_t generation() const { return generation_value; }

private:
  uint32_t observed_schema_generation = 0;
  uint32_t generation_value = 0;
};

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
  bool is_integer;  /**< True if the target stores whole numbers, so the GUI
                       steps by one; set for enums and plain integers alike. */
  bool animated;    /**< True if the parameter is currently animated. */
  bool readonly;    /**< True if the parameter is read-only (not editable). */
  bool preset;      /**< True if preset exports include the parameter. */
  const char *const *options; /**< Enum option labels, or null for a plain
                                 param; the value is the selected index. */
  int option_count; /**< Number of option labels; > 0 marks an enum. */
  const char *const *export_options; /**< C++ enum literals, or null. */
};

/**
 * @brief Snapshot an effect's parameters into `out`, in definition order.
 * @param effect Effect whose getParameters() sequence defines the order.
 * @param out Destination vector, cleared then filled in definition order;
 *            caller-owned so a reused vector amortizes its allocation.
 * @details Captures the order supplied by Effect's registered ParamList. Pair
 *          the snapshot with the effect's schema-generation token.
 */
inline void collect_param_views(const Effect &effect,
                                std::vector<ParamView> &out) {
  out.clear();
  for (const auto &def : effect.getParameters()) {
    out.push_back(ParamView{def.name, def.get(), def.min, def.max,
                            def.is_bool(), def.is_integer(), def.animated,
                            def.readonly, def.preset, def.options,
                            def.option_count, def.export_options});
  }
}

/**
 * @brief Fill `out` with current parameter values, in the SAME order as
 *        collect_param_views() — this is the getParamValues() stream.
 * @param effect Effect whose getParameters() values are read, in order.
 * @param out Destination vector, cleared (retaining capacity) then filled so
 *            that out[i] corresponds to collect_param_views()'s view[i].
 * @details Iterates Effect's registered ParamList independently of
 *          collect_param_views(). Callers reject the stream when its schema
 *          generation differs from the definition snapshot. The stream carries raw
 *          floats even for is_bool params (a bool is emitted as 0.0/1.0, not a
 *          JS boolean); consumers key off the definition type from
 *          collect_param_views(), not this value. `out.clear()` retains
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
