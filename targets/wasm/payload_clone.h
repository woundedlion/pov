/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
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

inline emscripten::val clone_payload(const emscripten::val &input) {
  return emscripten::val::take_ownership(
      clone_payload_handle(input.as_handle()));
}
