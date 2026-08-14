/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// FastNoiseLite configuration for Holosphere. First-party: the only file in
// core/vendor/ that is not upstream.
//
// PATCH RECORD — FastNoiseLite.h is upstream 1.1.1 plus seven in-tree edits,
// none of which fails to compile if a version bump drops it. They are pinned by
// tests/check_vendor_patches.cmake (CTest: unit_vendor_patches):
//   - the `#include "FastNoiseLite_config.h"` that pulls in this file
//   - four FASTNOISELITE_ONLY_OPENSIMPLEX2 guards (SetRotationType3D, both
//     GenNoiseSingle overloads, vector-noise dispatch)
//   - HS_O3_FN on SingleOpenSimplex2 (selective -O3)
//   - raw-octave paths used by first-party basis and derivative policies
//   - a raw vector-noise path used by spherical tangent displacement
//   - an analytic raw OpenSimplex2 gradient path used by spherical curl
//
// HS_O3_FN comes from engine/platform.h below, so the vendored header is no
// longer a standalone drop-in: a version bump has to re-apply the patches
// against this file rather than swap FastNoiseLite.h in isolation.
#pragma once

// This macro hard-routes GenNoiseSingle straight to OpenSimplex2 (bypassing the
// per-noise-type switch) and skips the 3D warp-transform setup. It does not
// itself remove the Cellular/Perlin/Value/OpenSimplex2S, fractal, or domain-warp
// code — those definitions stay in the header and are dropped by the compiler as
// dead code once nothing references them, which is what yields the lean binary.
// Remove this define to re-enable the full FastNoiseLite feature set.
#define FASTNOISELITE_ONLY_OPENSIMPLEX2

#include "engine/platform.h"
