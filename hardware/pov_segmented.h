/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_segmented.h
 * @brief Multi-Teensy segmented POV display driver for Phantasm.
 *
 * Phantasm uses 4 Teensys, each controlling a segment of the full LED run.
 * Each Teensy reads its hardware ID at boot to determine which Y-segment
 * it owns, renders the full effect canvas, but only pushes its segment
 * to the physical LEDs.
 *
 * TODO: Implement segmented driver.
 */
#pragma once

// TODO: Implement POVSegmented<TOTAL_PIXELS, NUM_SEGMENTS, RPM>
// - Hardware ID detection (GPIO pins or EEPROM)
// - Y-clipping: each Teensy only pushes its segment rows to DMA
// - Same DMA pipeline (dma_led.h) as single driver
// - Effects render full canvas; only the ISR is segment-aware
