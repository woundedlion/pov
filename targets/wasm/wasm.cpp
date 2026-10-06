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
 * @brief Registers HolosphereEngine, ShaderChainBindings, MeshOps, PaletteOps
 *        and the free math exports with Embind.
 */
EMSCRIPTEN_BINDINGS(holosphere_engine) {
  bind_engine();
  bind_mesh_ops();
  bind_palette_ops();
  bind_math_exports();
}

#endif
