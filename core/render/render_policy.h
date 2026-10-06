/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file render_policy.h
 * @brief Shared shape geometry and pole shading policy.
 */
#include "platform/build_features.h"

namespace Render {

/**
 * @brief Inner/outer radius ratio for star shapes (1/φ² ≈ 0.382).
 */
inline constexpr float STAR_INNER_RATIO = 0.382f;

/** @brief Longest column run a single shade may be splatted across. */
inline constexpr int POLE_LOD_MAX_RUN = 32;

/**
 * @brief Aggressiveness of near-pole azimuthal shading decimation; 0 disables.
 * @details Horizontal pixel pitch scales with sin(phi). The column run is
 * aggressiveness / sin(phi); the footprint depends on the display aspect,
 * LED angular size and per-column exposure. At 0 every run is one column.
 * Firmware fixes it at HS_POLE_LOD_DEFAULT; zero compiles out the decimated
 * walk there. A nonzero value requires hardware calibration.
 */
#ifndef HS_POLE_LOD_DEFAULT
#define HS_POLE_LOD_DEFAULT 0.0f
#endif
#ifdef ARDUINO
inline constexpr float pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
/** @brief Whether the decimated scan walk is compiled into this build. */
inline constexpr bool POLE_LOD_ENABLED = HS_POLE_LOD_DEFAULT > 0.0f;
#else
inline float pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
/**
 * @brief Whether the decimated scan walk is compiled into this build.
 * @details Always on host and WASM, where an aggressiveness of 0 disables it
 *          per scan.
 */
inline constexpr bool POLE_LOD_ENABLED = true;
#endif

} // namespace Render
