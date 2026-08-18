/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file constants.h
 * @brief Engine-wide compile-time constants: canvas bounds, star geometry,
 *        and pole-LOD tuning.
 */

#include "platform/platform.h" // CANVAS_W/CANVAS_H

/**
 * @brief Inner/outer radius ratio for star shapes (1/φ² ≈ 0.382).
 */
inline constexpr float STAR_INNER_RATIO = 0.382f;

/**
 * @brief Maximum horizontal resolution (width) for effects, from CANVAS_W.
 * @details Sizes the shared framebuffers (Effect::buffer_a/buffer_b) and bounds
 *          the Effect constructor, so a target that renders one small canvas
 *          reserves only what it draws (-DCANVAS_W/-DCANVAS_H per env). Builds
 *          that instantiate several resolutions — host tests, WASM, Phantasm —
 *          keep the 288x144 default, which bounds every one of them.
 */
inline constexpr int MAX_W = CANVAS_W;

/**
 * @brief Maximum vertical resolution (height) for effects, from CANVAS_H.
 */
inline constexpr int MAX_H = CANVAS_H;

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
inline constexpr bool POLE_LOD_ENABLED = HS_POLE_LOD_DEFAULT > 0.0f;
#else
inline float pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
inline constexpr bool POLE_LOD_ENABLED = true;
#endif
