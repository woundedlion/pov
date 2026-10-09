/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file payload_clone.h
 * @brief Exception-contained structured clones of caller-owned JS payloads.
 */
#pragma once

#include <emscripten.h>
#include <emscripten/val.h>

// clang-format off
EM_JS(emscripten::EM_VAL, clone_payload_handle, (emscripten::EM_VAL input), {
  try {
    return Emval.toHandle(structuredClone(Emval.toValue(input)));
  } catch (error) {
    if (ABORT || Module['HS_MODULE_DEAD'] || error instanceof WebAssembly.RuntimeError) throw error;
    return Emval.toHandle(null);
  }
});
// clang-format on

/**
 * @brief Returns a `structuredClone` of @p input, so later reads cannot run
 *        caller getters or observe caller mutation.
 * @details An uncloneable value yields JS `null`, which callers must reject as
 * an invalid payload. The clone error is rethrown when the module has aborted
 * or is marked dead, or when it is a `WebAssembly.RuntimeError`.
 * @return The clone, or JS `null` when @p input cannot be cloned.
 */
inline emscripten::val clone_payload(const emscripten::val &input) {
  return emscripten::val::take_ownership(
      clone_payload_handle(input.as_handle()));
}
