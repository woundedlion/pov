/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// FastNoiseLite configuration for Holosphere. First-party: the only file in
// core/vendor/ that is not upstream.
//
// PATCH RECORD — FastNoiseLite.h is upstream 1.1.1 plus six in-tree edits, and
// this file carries a seventh. Dropping the config include, guards, HS_O3_FN
// placement, or platform include can fail silently; dropping the raw-path
// additions fails to compile at their callers. All seven are pinned by
// tests/check_vendor_patches.cmake (CTest: unit_vendor_patches):
//   - the `#include "FastNoiseLite_config.h"` that pulls in this file
//   - three FASTNOISELITE_ONLY_OPENSIMPLEX2 guards (both
//     GenNoiseSingle overloads, vector-noise dispatch)
//   - selective-O3 placement on the scalar and vector OpenSimplex2 leaves
//   - raw-octave paths used by first-party basis and derivative policies
//   - a raw vector-noise path used by spherical tangent displacement
//   - an analytic raw OpenSimplex2 gradient path used by spherical curl
//   - the `#include "platform/platform.h"` below, which supplies HS_O3_FN
//
// HS_O3_FN comes from platform/platform.h below. Version bumps must preserve
// the configured patches.
#pragma once

// This macro hard-routes GenNoiseSingle straight to OpenSimplex2 (bypassing the
// per-noise-type switch). It does not
// itself remove the Cellular/Perlin/Value/OpenSimplex2S, fractal, or domain-warp
// code — those definitions stay in the header and are dropped by the compiler as
// dead code once nothing references them, which is what yields the lean binary.
// Remove this define to re-enable the full FastNoiseLite feature set.
#define FASTNOISELITE_ONLY_OPENSIMPLEX2

#include "platform/platform.h"
