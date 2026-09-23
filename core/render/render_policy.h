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

/**
 * @brief Longest column run a single shade may be splatted across.
 * @details Bounds the near-pole run so one shade never covers a visually
 *          significant arc, and keeps a run short relative to the narrowest
 *          clip segment.
 */
inline constexpr int POLE_LOD_MAX_RUN = 32;

/**
 * @brief Aggressiveness of near-pole azimuthal shading decimation; 0 disables.
 * @details A row at colatitude phi has horizontal pixel pitch sin(phi) times
 *          the vertical, so 1/sin(phi) columns share one physical LED
 *          footprint and need only one shade between them. The column run is
 *          `aggressiveness / sin(phi)`, so 1.0 tracks that footprint exactly
 *          and smaller values stay inside it. At 0 every run is one column and
 *          the scan is bit-identical to an undecimated walk.
 *
 *          The true masking width depends on the LED's angular size and the
 *          per-column exposure, so this is a hardware-calibrated knob rather
 *          than a derived constant. Firmware has no setter, so the starting
 *          value comes from HS_POLE_LOD_DEFAULT.
 */
#ifndef HS_POLE_LOD_DEFAULT
#define HS_POLE_LOD_DEFAULT 0.0f
#endif
#ifdef ARDUINO
inline constexpr float pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
/**
 * @brief Whether the decimated scan walk is compiled into this build.
 * @details Firmware has no setter, so the aggressiveness is a constant and a
 *          build left at 0 drops the decimation path outright.
 */
inline constexpr bool POLE_LOD_ENABLED = HS_POLE_LOD_DEFAULT > 0.0f;
#else
inline float pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
/**
 * @brief Whether the decimated scan walk is compiled into this build.
 * @details Host and WASM builds can raise `pole_lod_aggressiveness` at runtime,
 *          so the path is always compiled in and an aggressiveness of 0
 *          disables it per scan instead.
 */
inline constexpr bool POLE_LOD_ENABLED = true;
#endif

} // namespace Render
