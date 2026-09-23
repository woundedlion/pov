/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file constants.h
 * @brief Engine-wide compile-time canvas bounds.
 */

#include "platform/build_features.h" // CANVAS_W/CANVAS_H

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
