// FastNoiseLite configuration for Holosphere.
//
// PATCH RECORD — FastNoiseLite.h is upstream 1.1.1 plus four in-tree edits,
// none of which fails to compile if a version bump drops it. They are pinned by
// tests/check_vendor_patches.cmake (CTest: unit_vendor_patches):
//   - the `#include "FastNoiseLite_config.h"` that pulls in this file
//   - three FASTNOISELITE_ONLY_OPENSIMPLEX2 guards (SetRotationType3D, both
//     GenNoiseSingle overloads)
//   - HS_O3_FN on SingleOpenSimplex2 (selective -O3, docs/selective_o3_spec.md)
//
// This macro hard-routes GenNoiseSingle straight to OpenSimplex2 (bypassing the
// per-noise-type switch) and skips the 3D warp-transform setup. It does not
// itself remove the Cellular/Perlin/Value/OpenSimplex2S, fractal, or domain-warp
// code — those definitions stay in the header and are dropped by the compiler as
// dead code once nothing references them, which is what yields the lean binary.
// Remove this define to re-enable the full FastNoiseLite feature set.
#define FASTNOISELITE_ONLY_OPENSIMPLEX2

#include "../engine/platform.h"
