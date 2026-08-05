/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#ifdef __EMSCRIPTEN__

#include "targets/wasm/engine_bindings.h"
#include "targets/wasm/mesh_ops_bindings.h"
#include "targets/wasm/palette_bindings.h"
#include "targets/wasm/math_exports.h"

/**
 * @brief Registers the HolosphereEngine and MeshOps bindings with Embind so
 *        JavaScript can construct and call them.
 */
EMSCRIPTEN_BINDINGS(holosphere_engine) {
  bind_engine();
  bind_mesh_ops();
  bind_palette_ops();
  bind_math_exports();
}

#endif
