/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @brief Rotations Per Minute of the POV display.
 */
static constexpr unsigned int RPM = 480;

/**
 * @brief The number of physical LED pixels on the strip.
 */
static constexpr int NUM_PIXELS = 40;

/**
 * @brief Half the number of pixels (physical height).
 */
static constexpr int H = NUM_PIXELS / 2;

/**
 * @brief Virtual height used for mapping (H + buffer).
 */
static constexpr int H_VIRT = H + 3;

/**
 * @brief Maximum horizontal resolution (width) for effects.
 */
static constexpr int MAX_W = 96;

/**
 * @brief Macro to calculate the 1D array index from 2D coordinates (x, y).
 * @param x The horizontal coordinate (column).
 * @param y The vertical coordinate (row).
 * @return The 1D index.
 */
inline constexpr int XY(int x, int y) { return x * H + y; }
